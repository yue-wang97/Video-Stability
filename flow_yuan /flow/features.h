#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "cv.h"
#include "cxcore.h"
#include "highgui.h"


using namespace std;
using namespace cv;

struct greaterThanPtr :
	public std::binary_function<const float *, const float *, bool>
{
	bool operator () (const float * a, const float * b) const
		// Ensure a fully deterministic result of the sort
	{
		return (*a > *b) ? true : (*a < *b) ? false : (a > b);
	}
};

void goodFeatures(uchar* image, IplImage * cornerMap1, int width, int height,int maxCorners, float qualityLevel, float minDistance);
void cornerMinEigenVal(uchar* src, CvMat * eigenv, int width, int height);
void box_filter(CvMat* img, CvMat* result); 
void Sobel66(uchar* _src, CvMat * _dst, int width, int height, int ddepth, int dx, int dy, int ksize, double scale, double delta, int borderType = BORDER_DEFAULT);
void getSobelKernels66(CvMat * _kx, CvMat* _ky, int dx, int dy, int _ksize, bool normalize, int ktype);
void calcMinEigenVal(CvMat * _cov1, CvMat * _cov2, CvMat * _cov3, CvMat * _dst);
void threshold11(CvMat* _src, CvMat* _dst, float thresh, float maxval, int type);
void minMaxLoc11(CvMat* src, float* _minval, float* _maxval, int* _minidx, int* _maxidx);
void dilation(uchar* data, int width, int height);
void quick_sort(float *s, int l, int r, Point2f * m);
void MakeImageBlock(IplImage *block, IplImage *image,
	int centPosI, int centPosJ);
