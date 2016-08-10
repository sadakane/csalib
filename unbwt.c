/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#if 1
 #include <sys/timeb.h>
#else
 #include <sys/time.h>
 #include <sys/resource.h>
#endif
#include "typedef.h"
#include "csa.h"
#include "psi1.h"
#include "lf_dna.h"

#ifndef min
#define min(x,y) ((x)<(y)?(x):(y))
#endif

#if 1
typedef struct timeb mytimestruct;
void mygettime(mytimestruct *t)
{
  ftime(t);
}
double mylaptime(mytimestruct *before,mytimestruct *after)
{
  double t;
  t = after->time - before->time;
  t += (double)(after->millitm - before->millitm)/1000;
  return t;
}
#else
typedef mytimestruct struct rusage;
void mygettime(mytimestruct *t)
{
  getrusage(RUSAGE_SELF,t);
}
double mylaptime(mytimestruct *before,mytimestruct *after)
{
  double t;
  t = after->ru_utime.tv_sec - before->ru_utime.tv_sec;
  t += (double)(after->ru_utime.tv_usec
		- before->ru_utime.tv_usec)/1000000;
  return t;
}
#endif

#define PAGE (1<<10)
uchar buf[PAGE];

int main(int argc, char *argv[])
{
  i64 i,n;
  CSA csa;
  mytimestruct before,after;
  double t;

   if (argc<2) {
      fprintf(stderr, "syntax: suftest file\n");
      return 1;
   }

   csa_read(&csa,argc,argv);
   n = csa.n;

   mygettime(&before);
   {
     int m;
     FILE *out;
     out = fopen("output.dec","w");
     i = 0;
     while (i < n) {
       if ((i/PAGE) % PAGE == 0) {
         fprintf(stderr,"%ld \r",i/PAGE);  fflush(stderr);
       }
       m = PAGE/csa.k2;
       if (i+m >= n) m = n-i;
       csa.text(buf,&csa,i,i+m-1);
       fwrite(buf,csa.k2,m,out);
       i += m;
     }
     fwrite(buf,1,0,out);
     fclose(out);
   }
   mygettime(&after);
   t = mylaptime(&before,&after);
   fprintf(stderr,"time %f sec\n",t);
   return 0;
}
