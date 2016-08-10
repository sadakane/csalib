/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdlib.h>
#include <stdio.h>
#include "csa.h"
#include "unicode.h"


int main(int argc, char *argv[])
{
  i64 n;
  CSA SA;

  if (argc<3) {
    fprintf(stderr, "syntax: %s {indexfiles}\n", argv[0]);
    return 1;
  }

  csa_read(&SA,argc-1, argv+1);
  n = SA.n;

#if 0
{
  int i;
  rank_t x;
  unicode_t code;
  uchar buf[6], *p;
  x = SA.inverse(&SA, 0);
  for (i = 0; i < 1000; i++) {
    if (csa_utf8_T_psi(&SA, &x, &code) == -1) {
      printf("???\n");
    }
    p = &buf[0];
    unicode_to_string(&p, code);
    buf[unicode_len(code)] = 0;
    printf("%d code = %d (%s) rank = %ld\n", i, code, &buf[0], x);
  }
  for (i = 0; i < 1000; i++) {
    if (csa_utf8_BW_LF(&SA, &x, &code) == -1) {
      printf("???\n");
    }
    p = &buf[0];
    unicode_to_string(&p, code);
    buf[unicode_len(code)] = 0;
    printf("%d code = %d (%s) rank = %ld\n", i, code, &buf[0], x);
  }
}
#endif

#if 1
{
  int i;
  rank_t x;
  unicode_t code;
  uchar buf[6], *p;
  CSAFILE *csafile;

  csafile = csa_fdopen(&SA, NULL);
  for (i = 0; i < 1000; i++) {
    code = csa_fgetwc(csafile);
    x = csafile->rank;
    p = &buf[0];
    unicode_to_string(&p, code);
    buf[unicode_len(code)] = 0;
    printf("%d code = %d (%s) rank = %ld\n", i, code, &buf[0], x);
  }
  for (i = 0; i < 1000; i++) {
    code = csa_fgetwbw(csafile);
    x = csafile->rank;
    p = &buf[0];
    unicode_to_string(&p, code);
    buf[unicode_len(code)] = 0;
    printf("%d code = %d (%s) rank = %ld\n", i, code, &buf[0], x);
  }
}
#endif

  return 0;
}
