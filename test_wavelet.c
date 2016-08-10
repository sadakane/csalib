#include <stdlib.h>
#include <stdio.h>
#include "csa.h"
#include "wavelet.h"

int main(int argc, char *argv[])
{
  i64 i, j, j2, n, m, m2, pos;
  wavelet wt;

  wt_read(&wt, argv[1]);

  for (i=0; i<wt.n; i++) {
    printf("%ld c = %ld (%c)\n", i, wt.access(&wt, i), wt.access(&wt, i));
  }
}
