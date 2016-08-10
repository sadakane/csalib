/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "typedef.h"
#include "mman.h"
#include "csa.h"
#include "lf_dna2.h"

//#ifndef _BITVEC_T_
//#define _BITVEC_T_
//typedef u64 bitvec_t;
//#endif
#define DD (sizeof(bitvec_t)*8)

#define logSIGMADNA 3
//#define SIGMADNA (1<<logSIGMADNA)
#define SIGMADNA 5

//#define SB 32
//#define MB 256
#define logLB 16
#define LB (1<<logLB)


static int blog(i64 x)
{
i64 l;
  l = -1;
  while (x>0) {
    x>>=1;
    l++;
  }
  return l;
}

static u64 getuint(uchar *s, i64 i, i64 w)
{
  u64 x;
  i64 j;
  s += i*w;
  x = 0;
  for (j=0; j<w; j++) {
    x += ((u64)(*s++)) << (j*8);
  }
  return x;
}

static void putuint(uchar *s, i64 i, i64 x, i64 w)
{
  i64 j;
  s += i*w;
  for (j=0; j<w; j++) {
    *s++ = x & 0xff;
    x >>= 8;
  }
}

static void writeint(int k,i64 x,FILE *f)
{
  int i;
  for (i=k-1; i>=0; i--) {
    fputc(x & 0xff,f);
    x >>= 8;
  }
}

static void writeuint(int k,u64 x,FILE *f)
{
  int i;
  for (i=k-1; i>=0; i--) {
    fputc(x & 0xff,f);
    x >>= 8;
  }
}



static int setbit(bitvec_t *B, i64 i,bitvec_t x)
{
  i64 j,l;
  j = i / DD;
  l = i % DD;
  if (x==0) B[j] &= (~(1L<<(DD-1-l)));
  else if (x==1) B[j] |= (1L<<(DD-1-l));
  else {
    printf("error setbit x=%ld\n",x);
    exit(1);
  }
  return x;
}

static int setbits(bitvec_t *B, i64 i, bitvec_t x, int w)
{
  i64 j;
  
  for (j=0; j<w; j++) {
    setbit(B,i+j,(x>>(w-1-j)) & 1);
  }
  return 0;
}


static int lf_dna_BW_sub(lf_dna2 *lf,i64 i)
{
  bitvec_t x;
  x = lf->BW[i/(DD/logSIGMADNA)];
  x = (x >> (DD-logSIGMADNA-(i % (DD/logSIGMADNA))*logSIGMADNA)) & ((1<<logSIGMADNA)-1);
  return x;
}


static bitvec_t popcount(bitvec_t x)
{
  bitvec_t r;

  r = x;
  r = ((r & 0xaaaaaaaaaaaaaaaa)>>1) + (r & 0x5555555555555555);
  r = ((r & 0xcccccccccccccccc)>>2) + (r & 0x3333333333333333);
  r = ((r>>4) + r) & 0x0f0f0f0f0f0f0f0f;

  r *= 0x0101010101010101;
  r >>= 64-8;

  return r;
}

#define MP (0x1249249249249249ULL*2) // 0010010010010010010010010010010010010010010010010010010010010010
static i64 lf_dna_rankc_sub(lf_dna2 *lf, i64 i, int c)
{
  i64 j,r;
  i64 i2;
  int mb,lb;
  bitvec_t x,m,*p;
  static bitvec_t masktbl[SIGMADNA] = {0, MP*1, MP*2, MP*3, MP*4};

  if (i >= lf->last) i--;
  if (i < 0) return 0;

  mb = lf->mb;
  lb = lf->lb;

  r = getuint(lf->RL,(i/lb)*SIGMADNA + c,lf->k);
  if (c != 0) { // ACGT
    r += lf->RM[(i/mb)*(SIGMADNA-1) + (c-1)];
  } else {
    r += ((i/lb)/mb)*mb;
    r -= lf->RM[(i/mb)*(SIGMADNA-1) + 0];
    r -= lf->RM[(i/mb)*(SIGMADNA-1) + 1];
    r -= lf->RM[(i/mb)*(SIGMADNA-1) + 2];
    r -= lf->RM[(i/mb)*(SIGMADNA-1) + 3];
  }
  i2 = (i/mb)*mb;
  p = &lf->BW[i2/(DD/logSIGMADNA)];

  for (j=0; j+(DD/logSIGMADNA)-1 <= (i-i2); j+=(DD/logSIGMADNA)) {
    x = (*p++) ^ masktbl[c];
    x = (x | (x>>1) | (x>>2)) & MP;
    r += (DD/logSIGMADNA) - POPCOUNT3(x);
  }
  x = (*p) ^ masktbl[c];
  x = (x | (x>>1) | (x>>2)) & MP;
  m = MP >> (((i-i2) - j + 1) * logSIGMADNA);
  x |= m;
  r += (DD/logSIGMADNA) - POPCOUNT3(x);

  return r;
}


static int convertchar(uchar t)
{
  int c = 0;
  switch (t) {
  case '$':
  case 'N':  c = 0;  break;
  case 'a':
  case 'A':  c = 1;  break;
  case 'c':
  case 'C':  c = 2;  break;
  case 'g':
  case 'G':  c = 3;  break;
  case 't':
  case 'T':  c = 4;  break;
  default:  printf("error char = %c [%02x]\n",t,t);
  }
  return c;
}

static void make_tbl(lf_dna2 *lf)
{
  i64 n;
  i64 i;
  i64 C[SIGMADNA];
  int c,mb, lb;
  bitvec_t *BW;
  int k;

  n = lf->n;
  k = lf->k;
  BW = lf->BW;
//  mb = 1 << lf->lmb;
  mb = lf->mb;
  lb = lf->lb;

  lf->RL = mymalloc(k*((n+lb-1)/lb+1)*SIGMADNA);
  lf->RM = mymalloc(sizeof(lf->RM[0])*((n+mb-1)/mb+1)*(SIGMADNA-1));

  for (i=0;i<((n+lb-1)/lb)*SIGMADNA;i++) {
    putuint(lf->RL,i,0,k);
  }
  for (i=0;i<((n+mb-1)/mb)*(SIGMADNA-1);i++) lf->RM[i] = 0;

  for (c = 0; c < SIGMADNA; c++) C[c] = 0;

  fprintf(stderr,"making tables...\n");
  for (i = 0; i < n; i++) {
    if (i % 1000000 == 0) {
      fprintf(stderr,"%ld\r",i/1000);
      fflush(stderr);
    }
    if (i % lb == 0) {
      for (c = 0; c < SIGMADNA; c++) {
        putuint(lf->RL,(i / lb)*SIGMADNA + c,C[c],k);
      }
    }
    if (i % mb == 0) {
      for (c = 1; c < SIGMADNA; c++) {
        lf->RM[(i / mb)*(SIGMADNA-1) + (c-1)]
            = (u16)(C[c] - getuint(lf->RL,(i / lb)*SIGMADNA + c,k));
      }
    }
    c = lf_dna_BW_sub(lf,i);
    C[c]++;
  }
  
}

void lf_dna2_options(CSA *csa, char *p)
{
  lf_dna2 *lf;
  csa->psi_struc = lf = mymalloc(sizeof(lf_dna2));
//  lf->l = 128;
  lf->mb = 126;
  if (p[0] == 0) goto end;
  p++;
  if (p[0] != ':') {
    sscanf(p,"%d",&lf->mb);
    if (lf->mb % (DD/logSIGMADNA)) {
      fprintf(stderr, "BW_DNA2 MB=%d must be a multiple of %d\n",lf->mb, (DD/logSIGMADNA));
      exit(1);
    }
//    printf("L = %d\n",lf->l);
    while (p[0] != ':' && p[0] != 0) p++;
  }
  if (p[0] == 0) goto end;
  p++;
end:
  fprintf(stderr, "BW_DNA2 MB=%d\n",lf->mb);
}

static int lf_dna_BW(CSA *csa, i64 i)
{
  lf_dna2 *lf;
  lf = (lf_dna2 *)csa->psi_struc;

  if (i == lf->last) return -1; // EOF
  if (i > lf->last) i--;
  
  return csa->AtoC[lf_dna_BW_sub(lf,i)];
}

static i64 lf_dna_rankc(CSA *csa, i64 i, int c)
{
  lf_dna2 *lf;

  lf = (lf_dna2 *)csa->psi_struc;
  return lf_dna_rankc_sub(lf,i,c);
}


i64 lf_dna2_makeindex(CSA *csa, char *fname)
{
  i64 last;
  FILE *in,*out;
  bitvec_t *BW;
  char *fbw, *flst, *fbw2, *fidx;
  int k,c,c2;
  i64 size;
  i64 n,i;
  lf_dna2 *lf;
  int lb,mb;
  
  lf = (lf_dna2 *)csa->psi_struc;

  size = 0;

  mb = lf->mb;
  lb = lf->lb = (LB/mb)*mb;
//  lf->lmb = blog(L);
//  L = 1 << lf->lmb;
//  lf->l = L;

  k = strlen(fname);
  fbw = mymalloc(k+5);
  flst = mymalloc(k+5);
  fbw2 = mymalloc(k+5);
  fidx = mymalloc(k+5);
  sprintf(fbw,"%s.bw",fname);
  sprintf(flst,"%s.lst",fname);
  sprintf(fbw2,"%s.bw3",fname);
  sprintf(fidx,"%s.bwe",fname);

  in = fopen(flst,"r");
  if (in == NULL) {
    perror("lf_dna_makeindex:");  exit(1);
  }
  fscanf(in,"%ld",&last);
  printf("last = %ld\n",last);
  lf->last = last;
  fclose(in);


  in = fopen(fbw,"rb");
  if (in == NULL) {
    printf("lf_dna_makeindex: cannot open %s\n",fbw);
    exit(1);
  }
  fseek(in,0,SEEK_END);
  n = ftell(in);
  fseek(in,0,SEEK_SET);
  printf("n=%ld\n",n);
  csa->n = lf->n = n;
  csa->sigma = 256;

  for (i=0; i<csa->sigma; i++) csa->C[i] = 0;

  BW = mymalloc(sizeof(*BW) * (n/(DD/logSIGMADNA)+1));
  for (i=0; i<(n/(DD/logSIGMADNA)+1); i++) BW[i] = 0;
  lf->BW = BW;

  fprintf(stderr,"packing...\n");
  for (i = 0; i < n; i++) {
    i64 i2,i3;
    if (i % 1000000 == 0) {
      fprintf(stderr,"%ld\r",i/1000);
      fflush(stderr);
    }
    c = fgetc(in);
    csa->C[c]++;
    c2 = convertchar(c);
    i2 = i / (DD/logSIGMADNA);
    i3 = i % (DD/logSIGMADNA);
    setbits(BW,i2*DD+i3*logSIGMADNA,c2,logSIGMADNA);
  }
  fclose(in);

  out = fopen(fbw2,"w");
  if (out == NULL) {
    printf("lf_dna2_makeindex: cannot open %s\n",fbw2);
  }
  for (i=0; i<n/(DD/logSIGMADNA)+1; i++) {
    writeuint(sizeof(*BW),BW[i],out);
    size += sizeof(*BW);
  }
  fclose(out);

  lf->k = k = (blog(n+1)+1+8-1)/8;

  make_tbl(lf);

  out = fopen(fidx,"w");
  if (out == NULL) {
    printf("lf_dna_makeindex: cannot open %s\n",fidx);
  }

  writeint(1,ID_LF,out);
  writeint(1,k,out); /* #bytes of integer */
  writeint(k,n,out);
  writeint(k,last,out);
  writeint(k,mb,out);
  size += 1+1+3*k;

  writeint(1,ID_BWT_DNA2,out);
  size += 1;

  for (i=0;i<((n+lb-1)/lb)*SIGMADNA;i++) {
    writeint(k,getuint(lf->RL,i,k),out);
    size += k;
  }
  for (i=0;i<((n+mb-1)/mb)*(SIGMADNA-1);i++) {
    writeint(sizeof(lf->RM[0]),lf->RM[i],out);
    size += sizeof(lf->RM[0]);
  }
  fclose(out);

  csa->BW = lf_dna_BW;
  csa->rankc = lf_dna_rankc;

  csa->LF = csa_LF;

  free(fbw);
  free(flst);
  free(fbw2);
  free(fidx);

  return size;
}


void lf_dna2_read(CSA *csa, char *fname)
{
  char *fbw, *fbwi, *fname2;
  int k,id;

  i64 psize1,psize2;
  i64 n;
  lf_dna2 *lf;
  uchar *p, *q;
  int lb, mb;
  
  csa->psi_struc = lf = (lf_dna2 *)mymalloc(sizeof(lf_dna2));

  k = strlen(fname);
  fname2 = mymalloc(k-4+1);
  strncpy(fname2,fname,k-4);
  fname2[k-4] = 0;
  k -= 5;

  fbw = mymalloc(k+5+1);
  fbwi = mymalloc(k+5+1);

  sprintf(fbwi,"%s.bwe",fname2);

//  printf("psi_read: read %s\n",fbwi);

  lf->mapidx = mymmap(fbwi);
  if (lf->mapidx->addr==NULL) {
    perror("psi1_read: mmap2\n");
    exit(1);
  }
  p = q = (uchar *)lf->mapidx->addr;
  psize1 = lf->mapidx->len;

  id = getuint(p,0,1);  p += 1;
  if (id != ID_LF) {
    printf("lf_dna2_read: id = %d is not supported.\n",id);
    exit(1);
  }
  lf->k = k = getuint(p,0,1);  p += 1;
  lf->n = n = getuint(p,0,k);  p += k;
  lf->last = getuint(p,0,k);  p += k;
  lf->mb = mb = getuint(p,0,k);  p += k;
  lf->lb = lb = (LB/mb)*mb;
#if 0
  lf->lmb = blog(l);
  if ((1 << lf->lmb) != l) {
    printf("L=%d must be a power of 2.\n",l);
    exit(1);
  }
#endif

  id = getuint(p,0,1);  p += 1;
  lf->id = id;
  csa->id = id; // new

//  printf("lf_dna_read: psi_id = %d L = %d\n",id,l);
  switch (id) {
    case ID_BWT_DNA2:
      fprintf(stderr, "#lf format = BWT_DNA2\n");
      sprintf(fbw,"%s.bw3",fname2);
      break;
    default:
      fprintf(stderr, "lf_dna2_read: ID %d is not supported.\n",id);
      break;
  }

  lf->RL = (uchar *)p;  p += ((n+lb-1)/lb)*SIGMADNA * k;
  lf->RM = (u16 *)p;  p += ((n+mb-1)/mb)*(SIGMADNA-1) * sizeof(lf->RM[0]);

//  printf("lf_dna_read: map %s\n",fbw);
  lf->mapbwt = mymmap(fbw);
  if (lf->mapbwt->addr==NULL) {
    perror("psi1_read: mmap1\n");
    exit(1);
  }
  lf->BW = (bitvec_t *)lf->mapbwt->addr;
  psize2 = lf->mapbwt->len;
//  printf("psize2 = %ld\n",psize2);

//  printf("lf_dna_read: psize1 = %ld psize2 = %ld\n",psize1,psize2);
  lf->psize = psize1 + psize2;

  free(fbw);
  free(fbwi);
  free(fname2);

// user-specific functions
  csa->BW = lf_dna_BW;
  csa->rankc = lf_dna_rankc;

// default functions
  csa->BW_LF = csa_BW_LF;
  csa->LF = csa_LF;
  csa->selectc = csa_selectc;
  csa->psi = csa_psi_by_rankc_naive;
  csa->psi_succ = csa_psi_succ_naive;
  csa->psi_pred = csa_psi_pred_naive;
  csa->lookup = csa_lookup_lf;
  csa->inverse = csa_inverse_lf;
  csa->text = csa_text_lf;
  csa->search = csa_search_lf;
  csa->searchsub = csa_searchsub_lf;

}

