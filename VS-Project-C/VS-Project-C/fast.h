#ifndef FAST_H
#define FAST_H
#define  BARRIER   70
typedef struct { int x, y; } xy;
typedef unsigned char  byte;
/////�ǵ��⺯��////
///����ͼ������im,ͼ�����:xsize��ysize���Լ�������޺ͼ����Ľǵ�����numcorners,�����ؽǵ������
void  fast_corner_detect(const byte* im, int xsize, int ysize, int barrier, int *num, xy *ret);

#endif
