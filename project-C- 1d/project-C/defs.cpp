#include"defs.h"

///////////////////////////////////////////////////////CvMat/////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////cvCreateMat///////////////////////////////////////////////////////

uchar** alignPtr(uchar** ptr, int n )
{
	return (uchar**)(((size)ptr + n - 1) & -n);
}
void* fastMalloc(size size)
{
	uchar* udata = (uchar*)malloc(size + sizeof(void*) + CV_MALLOC_ALIGN);
	uchar** adata = alignPtr((uchar**)udata + 1, CV_MALLOC_ALIGN);
	adata[-1] = udata;
	return adata;
}
void * cvAlloc(size size)
{
	return fastMalloc(size);
}

void* cvAlignPtr(const void* ptr, int align )
{
	//CV_DbgAssert((align & (align - 1)) == 0);
	return (void*)(((size)ptr + align - 1) & ~(size)(align - 1));
}


/////////////////////////////////////////////////////cvCreateMat///////////////////////////////////////////////////////////


/////////////////////////////////////////////////Image///////////////////////////////////////////////
/////////////////////////////////////cvCreateImage////////////////////////////////////////////////////
IplImage2* cvInitImageHeader(IplImage2 * image, CvSize2 size, int depth,
	int channels)
{
	const char *colorModel, *channelSeq;
	memset(image, 0, sizeof(*image));
	image->nSize = sizeof(*image);
	image->width = size.width;
	image->height = size.height;
	image->nChannels = MAX(channels, 1);
	image->depth = depth;
	image->widthStep = (((image->width * image->nChannels *
		(image->depth & ~IPL_DEPTH_SIGN) + 7) / 8) + 4 - 1) & (~(4 - 1));
	const int imageSize_tmp = (int)image->widthStep*(int)image->height;
	image->imageSize = (int)imageSize_tmp;
	return image;
}

void cvCreateData_img(IplImage2* arr)
{
	IplImage2* img = (IplImage2*)arr;
	img->imageData = (char*)cvAlloc((size_t)img->imageSize);		
}

IplImage2 * cvCreateImageHeader(CvSize2 size, int depth, int channels)
{
	IplImage2 *img = 0;

	img = (IplImage2 *)cvAlloc(sizeof(*img));
	cvInitImageHeader(img, size, depth, channels);
	return img;
}

IplImage2 * cvCreateImage2(CvSize2 size, int depth, int channels)
{
	IplImage2 *img = cvCreateImageHeader(size, depth, channels);
	cvCreateData_img(img);

	return img;
}
/////////////////////////////////////cvCreateImage////////////////////////////////////////////////////

//////////////////////////////////////InitializeImage////////////////////////////////////////////////////
void InitializeImage(IplImage2 *image)
{
	int i, j, k;
	int height = image->height;
	int width = image->width;
	int channels = image->nChannels;
	int step = image->widthStep;
	uchar *imageData = (uchar *)image->imageData;
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			for (k = 0; k < channels; k++) {
				imageData[i*step + j*channels + k]
					= 0;
			}
		}
	}
}
//////////////////////////////////////InitializeImage////////////////////////////////////////////////////


///////////////////////////////////////////Array2Image//////////////////////////////////////////////////////
void   Array2Image(IplImage2 *img, char *a)
{
	int i, j, k;
	int height = img->height;
	int width = img->width;
	int channels = img->nChannels;
	int step = img->widthStep;
	char *imageData = img->imageData;
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			for (k = 0; k < channels; k++) {
				imageData[i*step + j*channels + k]
					= a[i*step + j*channels + k];
			}
		}
	}
}

///////////////////////////////////////////Array2Image//////////////////////////////////////////////////////

//////////////////////////////////////////cvCloneImage////////////////////////////////////////////////////
IplImage2* cvCloneImage(IplImage2* src)
{
	IplImage2* dst = 0;
	dst = (IplImage2*)cvAlloc(sizeof(*dst));
	memcpy(dst, src, sizeof(*src));
	if (src->imageData)
	{
		int size = src->imageSize;
		cvCreateData_img(dst);
		memcpy(dst->imageData, src->imageData, size);
	}
	return dst;
}
//////////////////////////////////////////cvCloneImage////////////////////////////////////////////////////



/////////////////////////////////////////////////////cvSize2/////////////////////////////////////////////////////////
CvSize2  cvSize2(int width, int height)
{
	CvSize2 s;

	s.width = width;
	s.height = height;

	return s;
}
/////////////////////////////////////////////////////cvSize2/////////////////////////////////////////////////////////



void SetH(float* h, int row, int col, float value)
{
	int type,step;
	type = CV_MAT_TYPE(CV_32FC1);
	step = CV_ELEM_SIZE(type);
	((float*)(void*)(h+step*row))[col] = (float)value;
}
float  GetH(float *h, int row, int col)
{
	int type, step;
	type = CV_MAT_TYPE(CV_32FC1);
	step = CV_ELEM_SIZE(type);
	return ((float*)(void*)(h+step*row))[col];
}
void copyH(float *a, float* b)
{
	int n = sizeof(a) / sizeof(a[0]);
	int i = 0;
	for (i; i < n; i++)
		b[i] = a[i];
}