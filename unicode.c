/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
//#include "unicode.h"
#include "csa.h"

#define MAXLEN 6

/////////////////////////////////////////////////////////////
// utf8の１バイト目なら 1 を返す
/////////////////////////////////////////////////////////////
int utf8_head(int c)
{
  return ((0 <= c && c <= 0x7f) || (0xc0 <= c && c <= 0xfd));
}

/////////////////////////////////////////////////////////////
// utf8の２バイト目以降なら 1 を返す
/////////////////////////////////////////////////////////////
int utf8_tail(int c)
{
  return (0x80 <= c && c <= 0xbf);
}

/////////////////////////////////////////////////////////////
// utf8 の1バイト目から符号長を求める
/////////////////////////////////////////////////////////////
int utf8_len(int c)
{
  int k;
  if (0 <= c && c <= 0x7f) k = 1;
  else if (0xc2 <= c && c <= 0xdf) k = 2;
  else if (0xe0 <= c && c <= 0xef) k = 3;
  else if (0xf0 <= c && c <= 0xf7) k = 4;
  else if (0xf8 <= c && c <= 0xfb) k = 5;
  else if (0xfc <= c && c <= 0xfd) k = 6;
  else k = -1;
  return k;
}

/////////////////////////////////////////////////////////////
// unicodeからそのutf8のバイト数を求める
/////////////////////////////////////////////////////////////
int unicode_len(unicode_t code)
{
  int len;
  if (code <= 0x7f) len = 1;
  else if (code <= 0x000007ff) len = 2;
  else if (code <= 0x0000ffff) len = 3;
  else if (code <= 0x0001ffff) len = 4;
  else if (code <= 0x03ffffff) len = 5;
  else if (code <= 0x7fffffff) len = 6;
  else len = -1;
  return len;
}

/////////////////////////////////////////////////////////////
// p で表されるutf8の文字列からunicodeを求める
/////////////////////////////////////////////////////////////
unicode_t string_to_unicode(uchar **p)
{
  unicode_t code;
  uchar *q;
  int i,len;

  q = *p;
  len = utf8_len(*q);
  if (len < 0) {
    code = -1;
    goto end;
  }
  if (len == 1) {
    code = *q++ & 0x7f;
  } else {
    code = *q++ & ((1<<(7-len))-1);
  }
  for (i=2; i<=len; i++) {
    code <<= 6;
    code |= *q++ & 0x3f;
  }

end:
  *p = q;
  return code;
}

/////////////////////////////////////////////////////////////
// unicodeをutf8の文字列に変換し， p で表される領域に書く
/////////////////////////////////////////////////////////////
int unicode_to_string(uchar **p, unicode_t code)
{
  int i,c,len;
  uchar *q;
  if (code <= 0x7f) len = 1;
  else if (code <= 0x000007ff) len = 2;
  else if (code <= 0x0000ffff) len = 3;
  else if (code <= 0x0001ffff) len = 4;
  else if (code <= 0x03ffffff) len = 5;
  else if (code <= 0x7fffffff) len = 6;
  else return -1;

  q = *p;
  if (len == 1) {
    *q++ = code;
    goto end;
  }
  c = (0xff << (8-len)) & 0xff;
  c |= (code >> (6*(len-1)));
  *q++ = c;
  for (i=2; i<=len; i++) {
    c = 0x80;
    c |= (code >> (6*(len-i))) & 0x3f;
    *q++ = c;
  }

end:
  *p = q;
  return 0;
}

/////////////////////////////////////////////////////////////
// p で始まる長さ len のバイト列の中で正しいutf8列を抜き出す
// b から e までの領域が答え
/////////////////////////////////////////////////////////////
int valid_utf8(uchar *p, int len, uchar **b, uchar **e)
{
  int c, k, i;

  *b = *e = NULL;
  while (len > 0) {
    c = *p;
    k = utf8_len(c);
    if (k > 0) break;
    p++;
    len--;
  }
  if (len == 0) return 0;
  *b = p;
  while (len > 0) {
    c = *p++;
    k = utf8_len(c);
    if (k <= 0) return 0;
    if (k > len) return 0;
    for (i=1; i<k; i++) {
      c = *p++;
      if (!utf8_tail(c)) return 0;
    }
    len -= k;
    *e = p-1;
  }
  return 1;
}

typedef struct {
  unicode_t head; // 文字コード
  rank_t l,r; // headから始まる接尾辞のランク
  int len; // 文字コードの元の(utf8での)バイト数
} unicode_s;

static int cmpset(const void *p, const void *q)
{
  unicode_s *pp, *qq;
  pp = (unicode_s *)p;  qq = (unicode_s *)q;
  return pp->head - qq->head;
}

static int csa_utf8_child_l_sub(CSA *csa, rank_t l, rank_t r, unicode_s **set_p, uchar *unicode_buf, int d, int num, int maxc)
{
  uchar head_tmp[csa->m];
  i64 L_tmp[csa->m], R_tmp[csa->m];
  int i,c,k;
  unicode_t code;
  uchar *p;
  unicode_s *set;

  set = *set_p;

  k = csa->child_l(csa,l,r,head_tmp, L_tmp, R_tmp);
  for (i=0; i<k; i++) {
    c = head_tmp[i];
    *unicode_buf = c;
    if (utf8_head(c)) {
      if (utf8_len(c) != d) continue; // utf8の途中で終わる場合は捨てる
      p = unicode_buf;
      code = string_to_unicode(&p);
      if (code < 0) {
        printf("invalid unicode\n");  exit(1);
      }
      if (num >= maxc) {// 格納する領域がない
        maxc = (maxc == 0) ? csa->sigma : maxc*2;
        set = realloc(set, maxc * sizeof(*set));
        if (set == NULL) {
          printf("realloc\n");  exit(1);
        }
      }
      set[num].head = code;
      set[num].l = L_tmp[i];
      set[num].r = R_tmp[i];
      set[num].len = d;
      num++;
    } else { // utf8の先頭でなければ再帰
     *set_p = set;
      num = csa_utf8_child_l_sub(csa, L_tmp[i], R_tmp[i], set_p, unicode_buf-1, d+1, num, maxc);
      set = *set_p;
    }
  }
  *set_p = set;
  return num;
}

/////////////////////////////////////////////////////////////
// ランク [l,r] のパタンの左に現れるUnicode文字集合を求める
/////////////////////////////////////////////////////////////
int csa_utf8_child_l(CSA *csa, rank_t l, rank_t r, unicode_t **head, rank_t **L, rank_t **R, int **len)
{
  uchar unicode_buf[MAXLEN];
  int i, num;
  unicode_s *set;

  set = NULL;
  num = csa_utf8_child_l_sub(csa, l, r, &set, &unicode_buf[MAXLEN-1], 1,0,0);

  qsort(set, num, sizeof(unicode_s), cmpset);

  *head = mymalloc(num * sizeof(**head));
  *L = mymalloc(num * sizeof(**L));
  *R = mymalloc(num * sizeof(**R));
  *len = mymalloc(num * sizeof(**len));

  for (i=0; i<num; i++) {
    (*head)[i] = set[i].head;
    (*L)[i] = set[i].l;
    (*R)[i] = set[i].r;
    (*len)[i] = set[i].len;
  }

#if 0
  printf("utf8_child_l num = %d\n", num);
  for (i=0; i<num; i++) {
    p = unicode_buf;
    unicode_to_string(&p, (*head)[i]);
    unicode_buf[(*len)[i]] = 0;
    printf("[%ld,%ld] len = %d char = %s\n", (*L)[i], (*R)[i], (*len)[i], unicode_buf);
  }
#endif

  free(set);
  return num;
}

static int csa_utf8_child_r_sub(CSA *csa, rank_t l, rank_t r, unicode_s **set_p, uchar *unicode_buf, int d, int clen, int num, int maxc, int plen)
{
  uchar head_tmp[csa->m];
  i64 L_tmp[csa->m], R_tmp[csa->m];
  int i,k;
  unicode_t code;
  uchar *p;
  unicode_s *set;

  set = *set_p;

  k = csa->child_r(csa,plen+d, l,r, head_tmp, L_tmp, R_tmp);
  if (d+1 == clen) { // 文字コードの最後のバイトならリストに追加
    for (i=0; i<k; i++) {
      unicode_buf[d] = head_tmp[i];
      p = unicode_buf;
      code = string_to_unicode(&p);
      if (code < 0) {
        printf("invalid unicode\n");  exit(1);
      }
      if (num >= maxc) {// 格納する領域がない
        maxc = (maxc == 0) ? csa->sigma : maxc*2;
        set = realloc(set, maxc * sizeof(*set));
        if (set == NULL) {
          printf("realloc\n");  exit(1);
        }
      }
      unicode_buf[d+1] = 0;
//      printf("%d head = %d (%s) [%ld,%ld] len = %d\n", num, code, unicode_buf, L_tmp[i], R_tmp[i], d+1);
      set[num].head = code;
      set[num].l = L_tmp[i];
      set[num].r = R_tmp[i];
      set[num].len = d+1;
      num++;
    }
  } else { // 再帰
    for (i=0; i<k; i++) {
      unicode_buf[d] = head_tmp[i];
      *set_p = set;
      num = csa_utf8_child_r_sub(csa, L_tmp[i], R_tmp[i], set_p, unicode_buf, d+1, clen, num, maxc, plen);
      set = *set_p;
    }
  }
  *set_p = set;
  return num;
}

/////////////////////////////////////////////////////////////
// 長さ plen バイト，ランク [l,r] のパタンの右に現れるUnicode文字集合を求める
/////////////////////////////////////////////////////////////
int csa_utf8_child_r(CSA *csa, int plen, rank_t l, rank_t r, unicode_t **head, rank_t **L, rank_t **R, int **len)
{
  uchar head_tmp[csa->m];
  i64 L_tmp[csa->m], R_tmp[csa->m];
  uchar unicode_buf[MAXLEN];
  int i, num, k, c, maxc, d;
  unicode_s *set;

  maxc = 0;
  k = csa->child_r(csa,plen, l,r,head_tmp, L_tmp, R_tmp);
  num = 0;  d = 0;
  set = NULL;
  for (i=0; i<k; i++) {
    c = head_tmp[i];
    unicode_buf[d] = c;
    if (utf8_head(c)) {
      if (utf8_len(c) == 1) {
        if (num >= maxc) {// 格納する領域がない
          maxc = (maxc == 0) ? csa->sigma : maxc*2;
          set = realloc(set, maxc * sizeof(*set));
          if (set == NULL) {
            printf("realloc\n");  exit(1);
          }
        }
        unicode_buf[d+1] = 0;
        //printf("%d head = %d (%s) [%ld,%ld] len = %d\n", num, c, unicode_buf, L_tmp[i], R_tmp[i], d+1);
        set[num].head = c;
        set[num].l = L_tmp[i];
        set[num].r = R_tmp[i];
        set[num].len = d+1;
        num++;
      } else {
        num = csa_utf8_child_r_sub(csa, L_tmp[i], R_tmp[i], &set, unicode_buf, d+1, utf8_len(c), num, maxc, plen);
      }
    }
  }

  qsort(set, num, sizeof(unicode_s), cmpset);

  if (num > 0) {
    *head = mymalloc(num * sizeof(**head));
    *L = mymalloc(num * sizeof(**L));
    *R = mymalloc(num * sizeof(**R));
    *len = mymalloc(num * sizeof(**len));

    for (i=0; i<num; i++) {
      (*head)[i] = set[i].head;
      (*L)[i] = set[i].l;
      (*R)[i] = set[i].r;
      (*len)[i] = set[i].len;
    }
  }

#if 0
  printf("utf8_child_r num = %d\n", num);
  for (i=0; i<num; i++) {
    p = unicode_buf;
    unicode_to_string(&p, (*head)[i]);
    unicode_buf[(*len)[i]] = 0;
    printf("[%ld,%ld] len = %d char = %s\n", (*L)[i], (*R)[i], (*len)[i], unicode_buf);
  }
#endif

  free(set);
  return num;
}


/////////////////////////////////////////////////////////////
// ランク x の接尾辞の先頭のutf8での1文字と，その文字を除いた接尾辞のランクを求める
// 入力: x
// 出力: code (文字のunicode), x (新しいランク)
// 関数の返り値: -1 (正しいutf8ではない), 0 (正常)
/////////////////////////////////////////////////////////////
int csa_utf8_T_psi(CSA *csa, rank_t *x, unicode_t *code)
{
  int i, c, len;
  rank_t xx;
  uchar unicode_buf[MAXLEN], *p;

  xx = *x;
  c = csa->T(csa, xx);
  len = utf8_len(c);
  if (len == -1) return -1;
  unicode_buf[0] = c;
  xx = csa->psi(csa, xx);

  for (i=1; i<len; i++) {
    c = csa->T(csa, xx);
    if (!utf8_tail(c)) return -1;
    unicode_buf[i] = c;
    xx = csa->psi(csa, xx);
  }
  p = &unicode_buf[0];
  *code = string_to_unicode(&p);
  *x = xx;

  return 0;
}

/////////////////////////////////////////////////////////////
// ランク x の接尾辞の直前のutf8での1文字と，その文字を加えた接尾辞のランクを求める
// 入力: x
// 出力: code (文字のunicode), x (新しいランク)
// 関数の返り値: -1 (正しいutf8ではない), 0 (正常)
/////////////////////////////////////////////////////////////
int csa_utf8_BW_LF(CSA *csa, rank_t *x, unicode_t *code)
{
  int i, c;
  rank_t xx;
  uchar unicode_buf[MAXLEN], *p;

  xx = *x;
  for (i=1; i<=MAXLEN; i++) {
    xx = csa->BW_LF(csa, xx, &c);
    if (c == -1) {
      if (i == 1) {
        *code = -1;  *x = 0;
        return 0; //?
      } else {
        return -1;
      }
    }
    unicode_buf[MAXLEN-i] = c;
    if (utf8_head(c)) break;
  }

  if (utf8_len(c) != i) return -1;

  p = &unicode_buf[MAXLEN-i];
  *code = string_to_unicode(&p);
  *x = xx;

  return 0;
}

/////////////////////////////////////////////////////////////
// 次のUnicode文字を返す
/////////////////////////////////////////////////////////////
unicode_t csa_fgetwc(CSAFILE *csafile)
{
  unicode_t code;
  int i, c, len;
  
  c = csa_fgetc(csafile);
  len = utf8_len(c);
  if (len == -1) {
    csafile->errorno = EILSEQ;
    return -1;
  }
  if (len == 1) {
    code = c & 0x7f;
  } else {
    code = c & ((1<<(7-len))-1);
  }

  for (i=2; i<=len; i++) {
    c = csa_fgetc(csafile);
    if (!utf8_tail(c)) {
      csafile->errorno = EILSEQ;
      return -1;
    }
    code <<= 6;
    code |= c & 0x3f;
  }

  return code;
}

/////////////////////////////////////////////////////////////
// 前のUnicode文字を返す
/////////////////////////////////////////////////////////////
unicode_t csa_fgetwbw(CSAFILE *csafile)
{
  unicode_t code;
  int i, c, len;
  uchar unicode_buf[MAXLEN], *p;

  for (i=1; i<=MAXLEN; i++) {
    c = csa_fgetbw(csafile);
    if (c == -1) {
      return -1;
    }
    unicode_buf[MAXLEN-i] = c;
    if (utf8_head(c)) break;
  }
  if (utf8_len(c) != i) return -1;

  p = &unicode_buf[MAXLEN-i];
  code = string_to_unicode(&p);

  return code;
}
