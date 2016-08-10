/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#ifndef _UNICODE_H
#define _UNICODE_H

#include "csa.h"

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

#endif // _UNICODE_H_
