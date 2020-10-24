#ifndef FAST_H
#define FAST_H
#define  BARRIER   60
typedef struct { int x, y; } xy;
typedef unsigned char  byte;

xy*  fast_corner_detect(const byte* im, int xsize, int ysize, int barrier, int *numcorners);

#endif
