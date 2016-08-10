/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "csa.h"

int main(int argc, char *argv[])
{
   csa_new_from_bwt(argc, argv);
   printf("write ok\n");

   return 0;
}
