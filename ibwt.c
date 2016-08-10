/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

typedef long i64;

int main(int argc, char *argv[])
{
   FILE *f;
   i64 i,j,n,a,b;
   unsigned char *bw,*s,*t;
   int c;
   i64 *LF;
   i64 last;
   i64 C[256+1];

   if (argc<3) {
     printf("usage: %s bw last\n",argv[0]);
     return 0;
   }

   f=fopen(argv[1], "rb");
   fseek(f, 0L, SEEK_END);  n=ftell(f);  fseek(f, 0L, SEEK_SET);
   printf("n = %ld\n",n);

   bw = malloc(sizeof(*bw) * n);
   s = malloc(sizeof(*bw) * (n+1));
   t = malloc(sizeof(*bw) * (n+1));
   LF = malloc(sizeof(*LF) * (n+1));
   if (bw == NULL || s == NULL || t == NULL || LF == NULL) {
     printf("not enough memory.\n");  exit(1);
   }

   printf("reading...\n");
   for (i=0; i<n; i++) {
     if (i % 1000000 == 0) {
       fprintf(stderr,"%ld \r",i,1000);
       fflush(stderr);
     }
     bw[i] = getc(f);
   }
   fclose(f);

   f=fopen(argv[2], "rb");
   fscanf(f,"%ld",&last);
   printf("last = %ld\n",last);

   printf("counting...\n");
   for (i=0; i<=256; i++) C[i] = 0;
   for (i=0; i<n; i++) {
     if (i % 1000000 == 0) {
       fprintf(stderr,"%ld \r",i,1000);
       fflush(stderr);
     }
     C[bw[i]]++;
   }
   
   a = 1;
   for (c=0; c<256; c++) {
     b = C[c];
     C[c] = a;
     for (i=a; i<a+b; i++) s[i] = c;
     a += b;
   }

   printf("making LF...\n");
   for (i=0; i<=n; i++) {
     if (i % 1000000 == 0) {
       fprintf(stderr,"%ld \r",i,1000);
       fflush(stderr);
     }
     if (i == last) {
       LF[i] = 0;
     } else {
       j = i;  if (j > last) j--;
       c = bw[j];
       LF[i] = C[c];
       C[c]++;
     }
   }

   printf("inverting...\n");
   j = 0;
   for (i=n; i>=1; i--) {
     if (i % 1000000 == 0) {
       fprintf(stderr,"%ld \r",i,1000);
       fflush(stderr);
     }
     j = LF[j];
     t[i-1] = s[j];
   }

   printf("writing text...\n");
   f = fopen("output.txt","wb");
   for (i=0; i<n; i++) {
     if (i % 1000000 == 0) {
       fprintf(stderr,"%ld \r",i,1000);
       fflush(stderr);
     }
     fputc(t[i],f);
   }
   fclose(f);

   free(LF);
   free(bw);
}
