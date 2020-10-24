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

void goodFeatures(InputArray image, OutputArray corners,int maxCorners, double qualityLevel, double minDistance);
void cornerMinEigenVal(Mat& src, Mat& eigenv);
void box_filter(IplImage* img, IplImage* result);
void Sobel66(InputArray _src, OutputArray _dst, int ddepth, int dx, int dy, int ksize, double scale, double delta, int borderType = BORDER_DEFAULT);
void getSobelKernels66(OutputArray _kx, OutputArray _ky, int dx, int dy, int _ksize, bool normalize, int ktype);
void calcHarris(const Mat& _cov, Mat& _dst, double k);
double getThreshVal_Otsu_8u11(const Mat &_src);
void threshold11(InputArray _src, OutputArray _dst, double thresh, double maxval, int type);
void minMaxLoc11(const Mat& src, double* _minval, double* _maxval, int* _minidx, int* _maxidx);
void dilation(uchar* data, int width, int height);
