CSA library
http://researchmap.jp/sada/csalib/
http://code.google.com/p/csalib/

Kunihiko Sadakane
National Institute of Informatics (NII)
sada@nii.ac.jp

Usage:

mkcsa filename [-P[n]{:L}{:O}] [-I{:D}{:D2}] [-Cc{:x}]
  It constructs a compressed full-text index from the BWT.
  The program reads filename.bw and filename.lst, and outputs
  Psi/BW.  It outputs sampled SA and inverse SA if -I option is given.
  -P[n] specifies the compression method.

    -P1    psi encoded by the gamma function
    -P2    psi encoded by the gamma function and run-length codes
    -P3    BW for DNA
    -P4    BW using Huffman-shaped wavelet tree.  Bit-vectors are compressed.
    -P5    BW using Huffman-shaped wavelet tree.  Bit-vectors are further compressed.
    -P6    psi encoded by the sparse array
    -P7    BW using Huffman-shaped wavelet tree.  Bit-vectors are not compressed (dense array).
    -P10   BW using Huffman-shaped wavelet tree.  Bit-vectors are compressed with sparse array.
    -P11   psi encoded by lengths of zero-runs and one-runs
    -P12   BW using Huffman-shaped wavelet tree.  Bit-vectors are compressed with lengths of zero-runs and one-runs.

  If -I is given, the common data (frequencies of characters, sampled SA) are output.
  Output filename is filename.idx.  This file is compatible with any compression method.

  L, O, D, D2 are options.
     L specifies the sampling rate.
     D specifies the sampling rate of SA
     D2 specifies the sampling rate of inverse SA

  The option -C specifies the alphabet size of the input BWT.  Every consecutive c bytes of
  BWT (in little endian) form a character.  If no x is given, alphabet size is set to the maximum
  value for c bytes, for example, it is 256 for -C1, and 65536 for -C2.  If x is given together
  with c, the alphabet size is set to x.  Note that in this case each character of BWT must be
  smaller than x.

  Recomended parameters are
     -P1:256 (good balance of size and speed)
     -P4:512 (good compression)
     -P12:512 (good compression)
     -P3:512 (absolutely good for DNA)

  The input files filename.bw and filename.lst can be created from a file "filename"
  by dbwt or ssss commands.  They are available on the project Web page.

csa filename.idx filename.{ext}
  It tests the index.

unbwt filename.idx filename.{ext}
  It decodes the compressed file.  Output filename is "output.dec".

API:
int csa_read(CSA *SA, int argc, char *argv[]);
  read the index into memory.  The structure SA stores pointers to supported functions.

SA->search(unsigned char *key,i64 keylen,CSA *csa,rank_t *li,rank_t *ri)
  search the index for pattern "key" of length "keylen".
  The output is the range of the suffix array [li,ri]
SA->lookup(CSA *csa, rank_t i)
  return SA[i]
SA->inverse(CSA *csa, pos_t suf)
  return ISA[i] (inverse SA)
SA->text(uchar *p,CSA *csa,pos_t i,pos_t j)
  extract T[i,j] to "p"

There are many other functions not listed here.

Compile:
  The library is compiled with gcc4.3 to use the Intel SSE4.2 technology supported
  by Core i7 processor.  SSE4.2 has the popcount instruction, which improves the speed
  of succinct data structures.  The library does not check if the CPU has SSE4.2.
  In this case, the program will terminate with illegal instruction error.
  If the -msse4.2 option is not given, the library will use software popcount.

Changes:
2011-12-20: Fixed a bug in csa_new_from_bwt in csa.c.
2010-11-11: Supported multi-byte alphabet.
2010-11-09: Improved index construction time for BW.  Fixed bugs in csa_new_from_bwt in csa.c, and
            decode_enum2 and comparray_construct in comparray.c.
2010-09-18: Added csa_fgetc_nd, csa_fgetbw_nd in csafile.c.
2010-08-27: Added BW_LF, csa_BW_LF, csa_BW_by_psi and csa_type in csa.c.  Added unicode.c, csafile.c.
            Fixed a bug in csa_BW_by_psi.  Supported Ruby 1.9.2.  Added csa_child_rs in csa.c.
2010-08-10: Fixed a bug in psi1.c.  Improved the speed for successor/predecessor of Psi.
            Added csa_left_diverse and csa_right_diverse in csa.c.
            Added sample programs (left/right maximal frequent patterns).
2010-07-18: Added ruby binding and Pizza&Chile interface.
2010-07-17: First release.
