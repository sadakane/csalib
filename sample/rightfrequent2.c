/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdlib.h>
#include <stdio.h>
#include "csa.h"
#include "unicode.h"

////////////////////////////////////////////////////////////////
// �E�ɑ�ȕp�o�p�^����񋓂���
// �p�^���̒����� maxlen �ȉ�
// �e�L�X�g��Unicode (utf8) �Ɖ��肵�C�����̐؂�ڂ����������̂������߂�
////////////////////////////////////////////////////////////////
void frequent_rightmaximal(CSA *csa, i64 l, i64 r, i64 th, uchar *key, int len, int maxlen)
{
  int i,c,k;
  int maximal;
  i64 ll, rr;
  unicode_t *head;
  i64 *L, *R;
  int *len2;
  uchar *p;

  if (len >= maxlen) return; // �p�^�����̍ő�l�𒴂�����I��

  k = csa_utf8_child_r(csa,len,l,r,&head, &L, &R, &len2); // �p�^�� P �̉E���ɏo�����镶���W�������߂�

  maximal = 1;
  for (i=0; i<k; i++) {
    c = head[i];
    ll = L[i];  rr = R[i];
//    printf("%d %d [%ld,%ld]\n",i,c,ll,rr);
    if (rr - ll + 1 >= th) { // �p�^�� Pc �̕p�x��臒l�ȏ�Ȃ�ΒT���𑱂���
      p = &key[len];
      unicode_to_string(&p, c);
      frequent_rightmaximal(csa, ll, rr, th, key, len+len2[i], maxlen);
      maximal = 0; // Pc ���p�o�Ȃ̂ŁCP �͋ɑ�ł͂Ȃ�
//      if (valid_utf8(key,len+1)) maximal = 0;
    }
  }

  if (maximal) { // P ���ɑ�Ȃ�Ε\��
    key[len] = 0;
    printf("[%ld,%ld] key = %s\n",l,r,key);
  }
  if (k > 0) {free(head);  free(L);  free(R);  free(len2);}
}

#define MAXLEN 100
int main(int argc, char *argv[])
{
  i64 n;
  CSA SA;
  double t;
  i64 th;
  uchar key[MAXLEN+1];

  if (argc<3) {
    fprintf(stderr, "syntax: %s threshold {indexfiles}\n", argv[0]);
    return 1;
  }

  t = atof(argv[1]);

  csa_read(&SA,argc-2, argv+2);
  n = SA.n;
  if (t < 1.0) th = n * t; else th = t;
  printf("n = %ld threshold = %ld\n", n, th);

  key[MAXLEN] = 0;
  frequent_rightmaximal(&SA, 0, n, th, key, 0, MAXLEN);

  return 0;
}
