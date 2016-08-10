/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "csa.h"

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


CSAFILE *csa_fdopen(CSA *csa, const char *mode)
{
  CSAFILE *csafile;
  csafile = mymalloc(sizeof(*csafile));
  csafile->csa = csa;
  csafile->mode = 0; //
  csafile->pos = 0;
  csafile->rank = csa->inverse(csa, csafile->pos);
  csafile->b = csafile->e = -1;
  csafile->bufsize = max(csa->D2, 1);
  csafile->errorno = 0;

  csafile->rankbuf = mymalloc((csafile->bufsize+1) * sizeof(rank_t));
  csafile->textbuf = mymalloc(csafile->bufsize * csa->k2);

  return csafile;
}

i64 csa_fread(void *ptr, i64 size, i64 nmemb, CSAFILE *csafile)
{
  i64 i, nmemb_read;
  pos_t s, t;
  
  s = csafile->pos;
  t = s + size * nmemb - 1;
  if (t >= csafile->csa->n) {
    nmemb_read = (csafile->csa->n - s) / size;
    t = s + size * nmemb_read - 1;
  } else {
    nmemb_read = nmemb;
  }
//  csafile->csa->text(ptr, csafile->csa, s, t); // rankを更新できないから使わない
  for (i=s; i<=t; i++) {
//    *(uchar *)ptr++ = csa_fgetc(csafile);
    putuint((uchar *)ptr, i-s, csa_fgetc(csafile), csafile->csa->k2);
  }
  return nmemb_read;
}

static void fgetc_sub_lf(CSAFILE *csafile)
{
  i64 t;
  rank_t rank;
  pos_t pos;
  int c;
  pos_t b, e;
  CSA *csa;
  
  csa = csafile->csa;

  if (csafile->bufsize == 1) {
    csafile->rankbuf[0] = csafile->rank;
    rank = csa->psi(csa, csafile->rank);
//    csafile->textbuf[0] = csa->T(csa, rank);
    putuint(csafile->textbuf, 0, csa->T(csa, rank), csafile->csa->k2);
    return;
  }

  t = csa->D2;
  b = csafile->b;  e = csafile->e;

  if (e+1 < csa->n) {
    rank = csa->inverse(csa, e+1);
    pos = e+1;
    csafile->rankbuf[pos-b] = rank;
    while (--pos >= b) {
      rank = csa->BW_LF(csa, rank, &c);
      csafile->rankbuf[pos-b] = rank;
//      csafile->textbuf[pos-b] = c;
      putuint(csafile->textbuf, pos-b, c, csafile->csa->k2);
    }
  } else {
    pos = csa->n;
    rank = 0;
    csafile->rankbuf[pos-b] = rank;
    while (--pos >= b) {
      rank = csa->BW_LF(csa, rank, &c);
      csafile->rankbuf[pos-b] = rank;
//      csafile->textbuf[pos-b] = c;
      putuint(csafile->textbuf, pos-b, c, csafile->csa->k2);
    }
    csafile->e = csa->n-1;
  }
}

static void fgetc_sub(CSAFILE *csafile)
{
  i64 i;
  pos_t pos;
  rank_t rank;
  CSA *csa;

  csa = csafile->csa;

  if (csa_type(csa) == CSA_PSI) {
    rank = csa->inverse(csa, csafile->b);
    for (i=0; i<csafile->bufsize; i++) {
      csafile->rankbuf[i] = rank;
//      csafile->textbuf[i] = csa->T(csa, rank);
      putuint(csafile->textbuf, i, csa->T(csa, rank), csafile->csa->k2);
      rank = csa->psi(csa, rank);
    }
    csafile->rankbuf[i] = rank;
  } else {
    fgetc_sub_lf(csafile);
  }
}

///////////////////////////////////////////
// 次の文字を返す
// 実行する前は，pos と rank は次に得られる文字の位置とランクを表している
// つまり c = T[SA[rank]] を返す
// pos は変化しない
///////////////////////////////////////////
int csa_fgetc_nd(CSAFILE *csafile)
{
  int c;
  pos_t pos;
  rank_t rank;

  i64 i;

  pos = csafile->pos;
  if (pos >= csafile->csa->n) return -1;

  if (pos < csafile->b || pos > csafile->e) { // out of decoded region
    csafile->b = (pos / csafile->bufsize) * csafile->bufsize;
    csafile->e = min(csafile->b + csafile->bufsize - 1, csafile->csa->n-1);
    fgetc_sub(csafile);
  }

//  c = csafile->textbuf[pos - csafile->b];
  c = getuint(csafile->textbuf, pos - csafile->b, csafile->csa->k2);
  return c;
}

///////////////////////////////////////////
// 次の文字を返す
// 実行する前は，pos と rank は次に得られる文字の位置とランクを表している
// つまり c = T[SA[rank]] を返す
// 実行後は pos は 1 増える
///////////////////////////////////////////
int csa_fgetc(CSAFILE *csafile)
{
  int c;
  c = csa_fgetc_nd(csafile);
  csafile->pos++;
  csafile->rank = csafile->rankbuf[csafile->pos - csafile->b];
  return c;
}


///////////////////////////////////////////
// 1つ前の文字(BW)を返す
// 実行する前は，pos と rank は次に得られる文字の1つ右の接尾辞の位置とランクを表している
// つまり c = T[SA[rank]-1] = T[pos-1] を返す．
// pos は 変化しない
///////////////////////////////////////////
int csa_fgetbw_nd(CSAFILE *csafile)
{
  int c;
  pos_t pos;
  rank_t rank;

  i64 i;

  pos = csafile->pos;
  if (pos-1 >= csafile->csa->n || pos-1 < 0) return -1;

  if (pos-1 < csafile->b || pos-1 > csafile->e) { // out of decoded region
    csafile->b = ((pos-1) / csafile->bufsize) * csafile->bufsize;
    csafile->e = min(csafile->b + csafile->bufsize - 1, csafile->csa->n-1);
    fgetc_sub(csafile);
  }

  pos--;
//  c = csafile->textbuf[pos - csafile->b];
  c = getuint(csafile->textbuf, pos - csafile->b, csafile->csa->k2);
  return c;
}

///////////////////////////////////////////
// 1つ前の文字(BW)を返す
// 実行する前は，pos と rank は次に得られる文字の1つ右の接尾辞の位置とランクを表している
// つまり c = T[SA[rank]-1] = T[pos-1] を返す．
// 実行後は pos は 1 減る
///////////////////////////////////////////
int csa_fgetbw(CSAFILE *csafile)
{
  int c;
  c = csa_fgetbw_nd(csafile);
  csafile->pos--;
  csafile->rank = csafile->rankbuf[csafile->pos - csafile->b];
  return c;
}

///////////////////////////////////////////
// ポインタの位置を設定する
///////////////////////////////////////////
int csa_fseek(CSAFILE *csafile, i64 offset, int whence)
{
  pos_t new_pos;
  switch (whence) {
    case SEEK_SET:
      new_pos = offset;
      break;
    case SEEK_CUR:
      new_pos = csafile->pos + offset;
      break;
    case SEEK_END:
      new_pos = csafile->csa->n + offset;
      break;
    default:
      csafile->errorno = EINVAL;
      return -1;
      break;
  }
  if (csafile->b <= new_pos && new_pos <= csafile->e) {
    csafile->pos = new_pos;
    csafile->rank = csafile->rankbuf[csafile->pos - csafile->b];
  } else {
    csafile->pos = new_pos;
    if (new_pos < csafile->csa->n) {
      csafile->rank = csafile->csa->inverse(csafile->csa, new_pos);
      csafile->b = csafile->e = new_pos - 1;
    } else {
      csafile->rank = -1;
    }
  }
  return 0;
}

///////////////////////////////////////////
// 現在のポインタの位置を返す
///////////////////////////////////////////
i64 csa_ftell(CSAFILE *csafile)
{
  return csafile->pos;
}

///////////////////////////////////////////
// ポインタの位置を先頭にする
///////////////////////////////////////////
void csa_rewind(CSAFILE *csafile)
{
  csa_fseek(csafile, 0, SEEK_SET);
  csafile->errorno = 0;
}

