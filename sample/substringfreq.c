/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "csa.h"
#include "mman.h"

#define FLAG_ONLYMAXIMAL 1
#define FLAG_ALLSUBSTR 2

////////////////////////////////////////////////////////////////
// key �̊e�����������csa���ł̕p�x�����߂�
////////////////////////////////////////////////////////////////
void substring_frequency(CSA *csa, uchar *key, i64 len, int flag)
{
  i64 pl, pr;
  i64 ps, pt; // ����������̍��[�ƉE�[
  i64 p;
  i64 qs, qt, ql, qr; // ����܂łɌ�����������������̈ʒu�ƃ����N
  i64 ol, or;

  ps = pt = 0;
  qs = -1;

  while (pt < len) {
//    if (pt % 10000 == 0) {fprintf(stderr,"%ld \r", pt); fflush(stderr);}
    pl = 0; pr = csa->n;
    for (p=pt; p>=ps; p--) { // key[ps..pt] ������
      ol = pl;  or = pr;
      if (csa->searchsub(key[p], csa, &pl, &pr) == -1) break; // ��v��������Ȃ�
      if (flag & FLAG_ALLSUBSTR) { // �S�Ă̕����������\��
        printf("found pos = [%ld,%ld] len = %ld freq = %ld\n", p, pt, pt-p+1, pr-pl+1);
      }
    }
    if (p >= ps) { // key[ps..pt] �Ƃ̊��S��v��������Ȃ�����
      if (qs >= 0 && (flag & FLAG_ONLYMAXIMAL)) { // �܂��\�����Ă��Ȃ��ɑ啔�������񂪂���
        printf("maximal pos = [%ld,%ld] len = %ld freq = %ld\n", qs, qt, qt-qs+1, qr-ql+1);
      }
      qs = p+1;  qt = pt;  ql = ol;  qr = or; // �����������̂�ۑ�
      ps = p+1; // ���[���~�X�}�b�`��1�E�ɂ���
    } else { // ��v����������
      qs = ps;  qt = pt;  ql = pl;  qr = pr; // �����������̂�ۑ�
    }

    if (!(flag & FLAG_ONLYMAXIMAL) && !(flag & FLAG_ALLSUBSTR)) {
      // �ɑ�ł͂Ȃ���������Ȃ����̂��\��
      printf("pos = [%d,%d] len = %d freq = %ld\n", qs, qt, qt-qs+1, qr-ql+1);
    }

    pt++; // �E�[��1�i�߂�
  }
  if (qs >= 0) { // �܂��\�����Ă��Ȃ��ɑ啔�������񂪂���
    printf("maximal pos = [%ld,%ld] len = %ld freq = %ld\n", qs, qt, qt-qs+1, qr-ql+1);
  }

}

int main(int argc, char *argv[])
{
  CSA SA;
  MMAP *map;
  int flag = 0;
  int ac;

  if (argc<4) {
    fprintf(stderr, "syntax: %s [-f<flag>] pattern-file {indexfiles}\n", argv[0]);
    return 1;
  }
  for (ac = 1; ac < argc; ac++) {
    if (argv[ac][0] != '-') break;
    if (toupper(argv[ac][1]) == 'F') flag = atoi(&argv[ac][2]);
  }

  map = mymmap(argv[ac]);
  csa_read(&SA,argc-(ac+1), &argv[ac+1]);

  substring_frequency(&SA, map->addr, map->len, flag);

  return 0;
}
