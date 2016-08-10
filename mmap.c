#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "typedef.h"
#include "mman.h"

#ifndef min
 #define min(x,y) (((x)<(y))?(x):(y))
#endif

#ifdef WIN32
MMAP *mymmap (char *fname)
{
void *base;
HANDLE fd,h;
size_t len;
MMAP *m;
  m = malloc(sizeof(*m));
  if (m==NULL) {perror("mymmap malloc");  exit(1);}
  fd = CreateFile(fname,GENERIC_READ,FILE_SHARE_READ,NULL,OPEN_EXISTING,
         FILE_ATTRIBUTE_NORMAL,NULL);
  if (fd==INVALID_HANDLE_VALUE) {
    printf("createfile\n");
    exit(1);
  }
  m->h1 = fd;
  len = GetFileSize(fd,0);
  len++;
  m->len = len;
  //h = CreateFileMapping (fd, NULL, PAGE_READONLY, 0, len, NULL);
  h = CreateFileMapping (fd, NULL, PAGE_READONLY, 0, 0, NULL);
  if (h==NULL) {
    printf("createfilemapping\n");
    exit(1);
  }
  m->h2 = h;
  //base = MapViewOfFile (h, FILE_MAP_READ, 0, 0, len);
  base = MapViewOfFile (h, FILE_MAP_READ, 0, 0, 0);
  if (base==NULL) {
    printf("mapviewoffile\n");
    return NULL;
  }
  m->addr = base;
  return m;
}

int mymunmap (MMAP *m)
{
 UnmapViewOfFile (m->addr);
 CloseHandle(m->h2);
 CloseHandle(m->h1);
 return 0;
}                

#else

MMAP *mymmap (char *fname)
{
int fd;
size_t len;
MMAP *m;
struct stat statbuf;
caddr_t base;
  m = malloc(sizeof(*m));
  if (m==NULL) {perror("mymmap malloc");  exit(1);}

  stat(fname,&statbuf);
  len = statbuf.st_size;
  fd = open(fname,O_RDONLY);
//  fd = open(fname,O_RDWR);
  if (fd == -1) {
    perror("open2\n");
    exit(1);
  }
  base = (void *)mmap(0,len,PROT_READ,MAP_SHARED,fd,0);
//  base = (void *)mmap(0,len,PROT_READ | PROT_WRITE,MAP_SHARED,fd,0);
  if (base==(caddr_t)-1) {
    perror("mmap1\n");
    exit(1);
  }
  m->addr = (void *)base;
  m->fd = fd;
  m->len = len;
  return m;
}

MMAP *mymmap_w (char *fname, size_t len)
{
  int fd;
  MMAP *m;
  struct stat statbuf;
  caddr_t base;
  size_t old_len;
  char c;

  m = malloc(sizeof(*m));
  if (m==NULL) {perror("mymmap_w malloc");  exit(1);}

  if (stat(fname,&statbuf) == 0) {
    old_len = statbuf.st_size;
  } else { // does not exist?
    old_len = 0;
  }

//  fd = open(fname,O_RDWR | O_CREAT | O_TRUNC);
  fd = open(fname,O_RDWR | O_CREAT);
  if (fd == -1) {
    perror("mymmap_w: open1\n");
    exit(1);
  }
  fchmod(fd, 0644);
#if 0
  if (lseek(fd, len-1, SEEK_SET) == -1) {
    perror("mymmap_w: lseek"); exit(1);
  }
  if (write(fd, &c, 1) != 1) {
    perror("mymmap_w: write"); exit(1);
  }
#else
  if (len > old_len) {
    uchar *buf;
    long l,bs,s;

    if (lseek(fd, old_len, SEEK_SET) == -1) {
      perror("mymmap_w: lseek"); exit(1);
    }
    bs = 1<<16;
    buf = malloc(bs);
    if (buf==NULL) {perror("mymmap_w: malloc buf");  exit(1);}
    for (l=0; l<bs; l++) buf[l] = 0;
    l = len - old_len;
    while (l > 0) {
      printf("write %ld \r",l);  fflush(stdout);
      s = min(l,bs);
      if (write(fd, buf, s) != s) {
        perror("mymmap_w: write"); exit(1);
      }
      l -= s;
    }
    free(buf);
  }
#endif
  lseek(fd, 0, SEEK_SET);
  
  base = (void *)mmap(0,len,PROT_READ | PROT_WRITE,MAP_SHARED,fd,0);
  if (base==(caddr_t)-1) {
    perror("mymmap_w\n");
    exit(1);
  }
  m->addr = (void *)base;
  m->fd = fd;
  m->len = len;
  return m;
}

void *mymremap(MMAP *m, size_t new_len)
{
  uchar *buf;
  long l,bs;
  caddr_t base;
struct stat statbuf;

  if (msync(m->addr, m->len, MS_ASYNC)) {
    perror("msync\n");
    exit(1);
  }
  if (munmap(m->addr,m->len)==-1) {
    perror("mymremap 1:");
  }

if (m->fd != -1) {
  if (new_len > m->len) {
    bs = 1<<16;
    buf = malloc(bs);
    if (buf==NULL) {perror("mymremap: malloc buf");  exit(1);}
    for (l=0; l<bs; l++) buf[l] = 0;

    if (lseek(m->fd, m->len, SEEK_SET) == -1) {
      perror("mymremap: lseek"); exit(1);
    }
    l = new_len - m->len;
    while (l > 0) {
      long s;
//      printf("write %ld \r",l);  fflush(stdout);
      s = min(l, bs);
      if (write(m->fd, buf, s) != s) {
        perror("mymremap: write"); exit(1);
      }
      l -= s;
    }
    free(buf);

  } else {
#if 0
    if (lseek(m->fd, new_len, SEEK_SET) == -1) {
      perror("mymremap: lseek2"); exit(1);
    }
    if (write(m->fd, &l, 0) != 0) {
      perror("mymremap: write2"); exit(1);
    }
#else
    if (ftruncate(m->fd, new_len)) {
      perror("mymremap: ftruncate"); exit(1);
    }
#endif
  }
}
  if (m->fd == -1) {
    base = (void *)mmap(0,new_len,PROT_READ | PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS,m->fd,0);
  } else {
    base = (void *)mmap(0,new_len,PROT_READ | PROT_WRITE,MAP_SHARED,m->fd,0);
  }
  if (base==(caddr_t)-1) {
    perror("mymremap: mmap\n");
    exit(1);
  }
  m->addr = (void *)base;
  m->len = new_len;

  return m->addr;
}

MMAP *mymmap_anom (size_t len)
{
int fd;
MMAP *m;
caddr_t base;
  m = malloc(sizeof(*m));
  if (m==NULL) {perror("mymmap malloc");  exit(1);}

  fd = -1;
  base = (void *)mmap(0,len,PROT_READ | PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS,fd,0);
  if (base==(caddr_t)-1) {
    perror("mymmap_anom\n");
    exit(1);
  }
  m->addr = (void *)base;
  m->fd = fd;
  m->len = len;
  return m;
}

int mymunmap (MMAP *m)
{
  if (munmap(m->addr,m->len)==-1) {
    perror("munmap 1:");
  }
  close(m->fd);
  return 0;
}                
#endif
