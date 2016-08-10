#ifndef _MYMMAP_H_
#define _MYMMAP_H_


#if 0
#include <windows.h>
#else
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/mman.h>
#endif

#ifdef WIN32
#define PAGE_READONLY          0x02     
#define SECTION_MAP_READ    0x0004
#define FILE_MAP_READ       SECTION_MAP_READ
#endif

typedef struct {
  void *addr;
  size_t len;
#ifdef WIN32
  HANDLE h1,h2;
#else
  int fd;
#endif
} MMAP;

MMAP *mymmap (char *fname);
MMAP *mymmap_w (char *fname, size_t len);
MMAP *mymmap_anom (size_t len);
int mymunmap (MMAP *m);
void *mymremap(MMAP *m, size_t new_len);

#endif
