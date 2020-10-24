/************************************************************************
Demo software: Invariant keypoint matching.
Author: David Lowe

defs.h:
This file contains the headers for a sample program to read images and
keypoints, then perform simple keypoint matching.
*************************************************************************/

/* From the standard C libaray: */
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>

/*------------------------------ Macros  ---------------------------------*/

#define ABS(x)    (((x) > 0) ? (x) : (-(x)))
#define MAX(x,y)  (((x) > (y)) ? (x) : (y))
#define MIN(x,y)  (((x) < (y)) ? (x) : (y))
typedef unsigned int   size;
typedef __int64 int64;
typedef unsigned char uchar;


#define  CV_MALLOC_ALIGN    64
#define CV_MAT_MAGIC_VAL    0x42420000
#define IPL_DEPTH_SIGN 0x80000000
#define INT_MAX       2147483647
#define CV_CN_MAX     512
#define CV_CN_SHIFT   3
#define CV_MAT_CONT_FLAG_SHIFT  14
#define CV_MAT_CONT_FLAG        (1 << CV_MAT_CONT_FLAG_SHIFT)
#define CV_DEPTH_MAX  (1 << CV_CN_SHIFT)

#define CV_8U   0
#define IPL_DEPTH_8U     8
#define CV_32F  5

#define CV_MAT_DEPTH_MASK       (CV_DEPTH_MAX - 1)
#define CV_MAT_DEPTH(flags)     ((flags) & CV_MAT_DEPTH_MASK)
#define CV_MAKETYPE(depth,cn) (CV_MAT_DEPTH(depth) + (((cn)-1) << CV_CN_SHIFT))
#define CV_32FC1 CV_MAKETYPE(CV_32F,1)




#define CV_MAT_TYPE_MASK        (CV_DEPTH_MAX*CV_CN_MAX - 1)
#define CV_MAT_TYPE(flags)      ((flags) & CV_MAT_TYPE_MASK)
#define CV_MAT_CN_MASK          ((CV_CN_MAX - 1) << CV_CN_SHIFT)
#define CV_MAT_CN(flags)        ((((flags) & CV_MAT_CN_MASK) >> CV_CN_SHIFT) + 1)
#define CV_MAT_DEPTH_MASK       (CV_DEPTH_MAX - 1)
#define CV_MAT_DEPTH(flags)     ((flags) & CV_MAT_DEPTH_MASK)
#define CV_ELEM_SIZE(type) \
    (CV_MAT_CN(type) << ((((sizeof(size_t)/4+1)*16384|0x3a50) >> CV_MAT_DEPTH(type)*2) & 3))
/*---------------------------- Structures --------------------------------*/

typedef struct CvSize2
{
	int width;
	int height;
}
CvSize2;
CvSize2  cvSize2(int width, int height);

typedef struct CvMat
{
	int type;
	int step;
	int* refcount;
	int hdr_refcount;
	union {
		uchar* ptr;
		short* s;
		int* i;
		float* fl;
		double* db;
	} data;
	union {
		int rows;
		int height;
	};
	union {
		int cols;
		int width;
	};
} CvMat;

void  cvmSet2(CvMat* mat, int row, int col, double value);
CvMat* cvCreateMat2(int height, int width, int type);
double  cvmGet(const CvMat* mat, int row, int col);
void  myMatcopy(CvMat *A, CvMat *B);

typedef struct IplImage
{
	int  nSize;             /**< sizeof(IplImage) */
	int imageSize;
	int  nChannels;         /**< Most of OpenCV functions support 1,2,3 or 4 channels */	   
	int  depth;             
	int  width;             /**< Image width in pixels.                           */
	int  height;            /**< Image height in pixels.                          */	  
	uchar *imageData;              
	int  widthStep;         /**< Size of aligned image row in bytes.    */
}
IplImage;

IplImage * cvCreateImage2(CvSize2 size, int depth, int channels);
void InitializeImage(IplImage *image);
void   Array2Image(IplImage *img, uchar *a);
IplImage * cvCloneImage(IplImage* src);

