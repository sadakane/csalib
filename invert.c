/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <stdlib.h>

#define BUFSIZE (1<<10)

char buf[BUFSIZE];

int main(int argc, char *argv[])
{
  FILE *in, *out;
  off_t ofs;
  off_t i,j,n,m;

  in = fopen(argv[1],"r");
  out = fopen(argv[2],"w");

  fseeko(in, 0, SEEK_END);
  n = ftell(in);
  fseeko(in, 0, SEEK_SET);

  ofs = n;
  while (1) {
    printf("%ld \r", ofs);  fflush(stdout);
    m = BUFSIZE;  if (ofs < m) m = ofs;
    fseeko(in, ofs-m, SEEK_SET);
    j = fread(buf, 1, m, in);
    if (j != m) {printf("read %ld bytes m=%ld\n");  exit(1);}
    for (j=m-1; j>=0; j--) {
      fputc(buf[j], out);
    }
    if (ofs == 0) break;
    ofs -= m;
  }

  
  fclose(in);
  fclose(out);

}
