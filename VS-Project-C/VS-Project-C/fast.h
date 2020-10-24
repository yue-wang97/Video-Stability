#ifndef FAST_H
#define FAST_H
#define  BARRIER   70
typedef struct { int x, y; } xy;
typedef unsigned char  byte;
/////角点检测函数////
///传入图像数据im,图像宽、高:xsize和ysize，以及检测门限和检测出的角点数量numcorners,并返回角点的坐标
void  fast_corner_detect(const byte* im, int xsize, int ysize, int barrier, int *num, xy *ret);

#endif
