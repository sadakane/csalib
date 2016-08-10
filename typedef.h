/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#ifndef _TYPEDEF_H_
#define _TYPEDEF_H_

#ifndef i32
typedef int i32;
#endif
#ifndef u32
typedef unsigned int u32;
#endif
#ifndef i64
typedef long i64;
/*typedef long long i64;*/  /*  for 32-bit gcc */
#endif
#ifndef u64
typedef unsigned long u64;
#endif
#ifndef u16
typedef unsigned short u16;
#endif
#ifndef i16
typedef short i16;
#endif



#ifndef uchar
typedef unsigned char uchar;
#endif

#ifndef ushort
typedef unsigned short ushort;
#endif

#ifndef uint
typedef unsigned int uint;
#endif

#ifndef ulong
typedef unsigned long ulong;
#endif


#ifndef byte
typedef unsigned char byte;
#endif
#ifndef word
typedef unsigned short word;
#endif
#ifndef dword
typedef unsigned int dword;
#endif
#ifndef qword
typedef unsigned long qword;
#endif

#ifndef _BITVEC_T_
#define _BITVEC_T_
typedef u64 bitvec_t;
#endif


#endif // _TYPEDEF_H_
