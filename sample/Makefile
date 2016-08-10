all: leftfrequent rightfrequent maximalfrequent proximity approxfrequent substringfreq unicode_sample

#CC = gcc43
CC = gcc

#CFLAGS = -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -msse4.2 -O9 -fPIC -I..
#CFLAGS = -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -O9 -fPIC -I..
CFLAGS = -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -g -fPIC -I..

leftfrequent: leftfrequent.c ../csa.h ../csa.a
	$(CC)  $(CFLAGS) -o leftfrequent leftfrequent.c ../csa.a
rightfrequent: rightfrequent.c ../csa.h ../csa.a
	$(CC)  $(CFLAGS) -o rightfrequent rightfrequent.c ../csa.a
maximalfrequent: maximalfrequent.c ../csa.h ../csa.a
	$(CC)  $(CFLAGS) -o maximalfrequent maximalfrequent.c ../csa.a
proximity: proximity.c ../csa.h ../csa.a
	$(CC)  $(CFLAGS) -o proximity proximity.c ../csa.a
approxfrequent: approxfrequent.c ../csa.h ../csa.a
	$(CC)  $(CFLAGS) -o approxfrequent approxfrequent.c ../csa.a
substringfreq: substringfreq.c ../csa.h ../csa.a
	$(CC)  $(CFLAGS) -o substringfreq substringfreq.c ../csa.a
unicode_sample: unicode_sample.c ../csa.h ../csa.a
	$(CC)  $(CFLAGS) -o unicode_sample unicode_sample.c ../csa.a
