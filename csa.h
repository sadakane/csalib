/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
/* compressed suffix arrays 

 */

#ifndef _CSA_H_
#define _CSA_H_

#define VERSION 2009071800
//#define VERSION 2010111000

#define CSA_PSI 1
#define CSA_BW 2

#define ID_HEADER 0x00
#define ID_SA     0x01
#define ID_ISA    0x02
#define ID_PSI    0x03
#define ID_LF     0x04
#define ID_HUFFMAN     0x05
#define ID_WAVELET 0x06
#define ID_HUFFMAN2     0x07
#define ID_DBLOCK     0x08
#define ID_DARRAY     0x09
#define ID_CRAM     0x0a
#define ID_CRAMBG     0x0b
#define ID_REPAIR  0x0c
#define ID_IREPAIR  0x0d
#define ID_DEBRUIJN 0x0e
#define ID_WAVELET2 0x0f
#define ID_DENSEARRAY 0x10
#define ID_WAVELET3 0x11
#define ID_BITVECTOR 0x12
#define ID_VLA 0x13
#define ID_VLSTREAM 0x14
#define ID_ONLINEBP 0x15
#define ID_VLSTREAMA 0x16
#define ID_BP 0x17

#define ID_ARRAY_FL 0x00
#define ID_DIFF_GAMMA 0x01
#define ID_DIFF_GAMMA_RL 0x02
#define ID_BWT_DNA 0x03
#define ID_BWT_WT 0x04
#define ID_BWT_WT_HUF 0x05
#define ID_SPARSE4 0x06
#define ID_BWT_WT_DENSE 0x07
#define ID_DIFF_GAMMA_SPARSE 0x08
#define ID_DIFF_GAMMA_RL_SPARSE 0x09
#define ID_BWT_WT_SPARSE4 0x0a
#define ID_DIFF_GAMMA_RR 0x0b
#define ID_BWT_WT_RR 0x0c
#define ID_BWT_HUF 0x0d
#define ID_BWT_BIT 0x0e
#define ID_BWT_DNA2 0x0f


#include "typedef.h"

#ifndef min
#define min(x,y) ((x)<(y)?(x):(y))
#endif
#ifndef max
#define max(x,y) ((x)>(y)?(x):(y))
#endif

typedef i64 rank_t, pos_t;

typedef struct csa {
  i64 n; /* length of the text */
  int k; /* number of bytes in an integer */
  int m; /* the number of distinct characters in the text */
  int D; /* interval between two SA values stored explicitly */
  int D2; /* interval between two inverse SA values stored explicitly */
  int sigma; /* alphabet size */
  int k2; /* the number of bytes in a character */
  i64 *C; /* frequency of characters */
  i64 *K; /* table of cumulative frequency */
  int *CtoA;
  int *AtoC; /* table of character codes */
  uchar *SA,*ISA;
  int id; /* format of psi/bw */
  void *psi_struc;
  void *mapidx;
  i64 psize, isize;

  pos_t (*lookup)(struct csa *csa, rank_t i);
  rank_t (*psi)(struct csa *csa,rank_t i);
  rank_t (*psi_succ)(struct csa *csa, rank_t pr, rank_t l, rank_t r);
  rank_t (*psi_pred)(struct csa *csa, rank_t pl, rank_t l, rank_t r);
  rank_t (*LF)(struct csa *csa,rank_t i);
  rank_t (*rankc)(struct csa *csa,rank_t i, int c);
  rank_t (*selectc)(struct csa *csa,rank_t i, int c);
  rank_t (*inverse)(struct csa *csa, pos_t suf);
  void (*text)(uchar *p,struct csa *csa,pos_t i,pos_t j);
  i64 (*substring)(uchar *p,struct csa *csa,rank_t r,i64 len);
  i64 (*substring_lf)(uchar *p,struct csa *csa,rank_t r,i64 len);
  int (*T)(struct csa *csa,rank_t i);
  int (*BW)(struct csa *csa,rank_t i);
  int (*BW_rank)(struct csa *csa,i64 i, rank_t *r);
  rank_t (*BW_LF)(struct csa *csa, rank_t i, int *bw);
  int (*head)(struct csa *csa,rank_t i);
  i64 (*search)(uchar *key,i64 keylen,struct csa *csa,rank_t *li,rank_t *ri);
  i64 (*searchsub)(int c, struct csa *csa, rank_t *ll, rank_t *rr);
  int (*child_l)(struct csa *csa, rank_t l, rank_t r, uchar *tail, rank_t *ll, rank_t *rr);
  int (*child_r)(struct csa *csa, i64 len, rank_t l, rank_t r, uchar *tail, rank_t *ll, rank_t *rr);

} CSA;


/* calculate SA[i] */
pos_t csa_lookup(CSA *csa, rank_t i);

/* calculate Psi[i] */
rank_t csa_psi(CSA *csa,rank_t i);

/* calculate SA^{-1}[i] */
rank_t csa_inverse(CSA *csa, pos_t suf);

/* decode T[i..j] into p */
void csa_text(unsigned char *p,CSA *csa, pos_t i, pos_t j);

/* decode T[SA[pos]..SA[pos]+len-1] into p */
i64 csa_substring(unsigned char *p,CSA *csa,rank_t r,i64 len);
i64 csa_substring_lf(uchar *p,CSA *csa,rank_t r,i64 len);
i64 csa_substring_lf_naive(uchar *p,CSA *csa,rank_t r,i64 len);

void csa_new_from_bwt(int argc, char *argv[]);
int csa_read(CSA *SA, int argc, char *argv[]);

i64 csa_search(uchar *key,i64 keylen,CSA *csa,rank_t *li,rank_t *ri);
i64 csa_search_lf(uchar *key, i64 keylen, CSA *csa, rank_t *ll, rank_t *rr);
rank_t csa_searchsub(int c, CSA *csa, rank_t *ll, rank_t *rr);
rank_t csa_searchsub_lf(int c, CSA *csa, rank_t *ll, rank_t *rr);
i64 csa_search_r(i64 keylen,int c, CSA *csa,rank_t *li,rank_t *ri);
int csa_left_diverse(CSA *csa, rank_t l, rank_t r);
int csa_right_diverse(CSA *csa, rank_t l, rank_t r, i64 length);
i64 csa_search_prefix(CSA *csa, uchar *pattern, i64 length, rank_t *s, rank_t *t);


int csa_child_l(CSA *csa, rank_t l, rank_t r, uchar *head, rank_t *ll, rank_t *rr);
int csa_child_r(CSA *csa, i64 len, rank_t l, rank_t r, uchar *tail, rank_t *ll, rank_t *rr);
//int csa_child_r0(CSA *csa, i64 len, rank_t l, rank_t r, uchar *tail, rank_t *ll, rank_t *rr);
int csa_child_rs(CSA *csa, uchar *substr, i64 len, uchar *tail, rank_t *ll, rank_t *rr);
pos_t csa_lookup_lf(CSA *csa, rank_t i);
rank_t csa_inverse_lf(CSA *csa, pos_t suf);
void csa_text_lf(uchar *p,CSA *csa, pos_t s, pos_t t);
rank_t csa_LF(CSA *csa, rank_t i);
rank_t csa_psi_by_rankc_naive(CSA *csa, rank_t i);
rank_t csa_selectc(CSA *csa, i64 i, int c);
rank_t csa_psi_pred_naive(CSA *csa, rank_t pr, rank_t l, rank_t r);
rank_t csa_psi_succ_naive(CSA *csa, rank_t pr, rank_t l, rank_t r);
int csa_BW_by_psi(CSA *csa, rank_t i);
rank_t csa_BW_LF(CSA *csa, rank_t i, int *bw);
rank_t csa_BW_LF_by_psi(CSA *csa, rank_t i, int *bw);
rank_t csa_LF_by_psi(CSA *csa, rank_t i);
int csa_BW_rank(CSA *csa,i64 i, rank_t *r);
int csa_T(CSA *csa,rank_t i);
int csa_head_rank(CSA *csa,rank_t i);
//void csa_error(void);
int csa_type(CSA *csa);

void *mymalloc(size_t n);
void myfree(void *p, size_t s);
//u64 getuint(uchar *s, i64 i, i64 w);

void bw_to_psi(FILE *out, CSA *csa, char *fbw, char *flst, int *k);


// stream
typedef struct {
  CSA *csa;
  pos_t b, e; // textbuf に格納されている部分文字列の先頭と最後の位置
  pos_t pos; // 次に出力する文字の位置
  rank_t rank; // 次に出力する文字のランク
  uchar *textbuf; // 復元した部分文字列
  rank_t *rankbuf; // textbufの各文字から始まる接尾辞のランク
  i64 bufsize; // textbufのバイト数
  u64 mode; // ファイルをオープンするモード（未使用）
  int errorno; // エラー番号
} CSAFILE;

CSAFILE *csa_fdopen(CSA *csa, const char *mode);
i64 csa_fread(void *ptr, i64 size, i64 nmemb, CSAFILE *csafile);
int csa_fgetc(CSAFILE *csafile);
int csa_fgetc_nd(CSAFILE *csafile);
int csa_fgetbw(CSAFILE *csafile);
int csa_fgetbw_nd(CSAFILE *csafile);
int csa_fseek(CSAFILE *csafile, i64 offset, int whence);
i64 csa_ftell(CSAFILE *csafile);
void csa_rewind(CSAFILE *csafile);

// unicode
typedef int unicode_t;

int utf8_head(int c);
int utf8_tail(int c);
int utf8_len(int c);
int unicode_len(unicode_t code);

unicode_t string_to_unicode(uchar **p);
int unicode_to_string(uchar **p, unicode_t code);
int valid_utf8(uchar *p, int len, uchar **b, uchar **e);

int csa_utf8_child_l(CSA *csa, rank_t l, rank_t r, unicode_t **head, rank_t **L, rank_t **R, int **len);
int csa_utf8_child_r(CSA *csa, int plen, rank_t l, rank_t r, unicode_t **head, rank_t **L, rank_t **R, int **len);

int csa_utf8_T_psi(CSA *csa, rank_t *x, unicode_t *code);
int csa_utf8_BW_LF(CSA *csa, rank_t *x, unicode_t *code);
unicode_t csa_fgetwc(CSAFILE *csafile);
unicode_t csa_fgetwbw(CSAFILE *csafile);

#endif // _CSA_H_
