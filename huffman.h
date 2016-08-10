/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#ifndef _HUFFMAN_H_
#define _HUFFMAN_H_

#include "csa.h"

#define HUFTBLWIDTH 8
#define HUFTBLSIZ (1<<HUFTBLWIDTH)
#define Huffman_MAXLEN 32

typedef struct {
  int n;  // alphabet size
//  int *v;
  int *left, *right;
  byte *clen;
  u64 *code;
  short tbl[HUFTBLSIZ];
} Huffman;

typedef struct {
  int n;  // alphabet size
  int m;  // number of characters with positive frequency
//  int *v;
  int *left, *right;
  byte *clen;
  u64 *code;
//  short tbl[HUFTBLSIZ];
  int *tbl;
  int tbl_width;
} Huffman2;

void freeHuffman(Huffman *p);
void freeHuffman2(Huffman2 *p);
int DecodeHuffman(Huffman *h, u64 x);
int DecodeHuffman2(Huffman2 *h, u64 x);
int DecodeHuffman_tbl(Huffman *h, u64 x);
int DecodeHuffman2_tbl(Huffman2 *h, u64 x);
Huffman *MakeHuffmanTree(int n, double *freq);
Huffman2 *MakeHuffman2Tree(int n, double *freq);
Huffman *MakeHuffmanTree2(int n, double *freq);
//Huffman2 *MakeHuffman2Tree2(int n, double *freq, int tbl_width);
Huffman2 *MakeHuffman2Tree2(int n, i64 *freq, int tbl_width);
void Huffman_write(Huffman *h, FILE *out);
void Huffman2_write(Huffman2 *h, FILE *out);
Huffman *Huffman_read(uchar **map);
Huffman2 *Huffman2_read(uchar **map);
i64 Huffman2_usedmemory(Huffman2 *h);

#endif // _HUFFMAN_H_
