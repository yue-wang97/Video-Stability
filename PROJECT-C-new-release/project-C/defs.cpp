#include"defs.h"

///////////////////////////////////////////////////////CvMat/////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////cvCreateMat///////////////////////////////////////////////////////

static void icvCheckHuge(CvMat2* arr)
{
	if ((int)arr->step*arr->rows > INT_MAX)
		arr->type &= ~CV_MAT_CONT_FLAG;
}

void* fastMalloc(size size)
{
	uchar* udata = (uchar*)malloc(size);
	return udata;
}
void * cvAlloc(size size)
{
	return fastMalloc(size);
}
CvMat2* cvCreateMatHeader(int rows, int cols, int type)
{
	type = CV_MAT_TYPE(type);
	int min_step = CV_ELEM_SIZE(type);
	min_step *= cols;
	CvMat2* arr = (CvMat2*)cvAlloc(sizeof(*arr));
	arr->step = min_step;
	arr->type = CV_MAT_MAGIC_VAL | type | CV_MAT_CONT_FLAG;
	arr->rows = rows;
	arr->cols = cols;
	arr->data.ptr = 0;
	icvCheckHuge(arr);
	return arr;
}

void cvCreateData(CvMat2* arr)
{
	size step, total_size;
	CvMat2* mat = (CvMat2*)arr;
	step = mat->step;
	if (mat->rows == 0 || mat->cols == 0)
		return;
	if (step == 0)
		step = CV_ELEM_SIZE(mat->type)*mat->cols;
	int _total_size = (int)step*mat->rows;
	total_size = (size)_total_size;
	mat->data.ptr = (uchar*)cvAlloc(total_size);
}

CvMat2* cvCreateMat2(int height, int width, int type)
{
	CvMat2* arr = cvCreateMatHeader(height, width, type);
	cvCreateData(arr);
	return arr;
}

/////////////////////////////////////////////////////cvCreateMat///////////////////////////////////////////////////////////

////////////////////////////////////////cvmSet///////////////////////////////////////////////////
void  cvmSet2(CvMat2* mat, int row, int col, float value)
{
	int type;
	type = CV_MAT_TYPE(mat->type);
	if (type == CV_32FC1)
		((float*)(void*)(mat->data.ptr + (size)mat->step*row))[col] = (float)value;
	else
	{
		((double*)(void*)(mat->data.ptr + (size)mat->step*row))[col] = value;
	}
}
////////////////////////////////////////cvmSet////////////////////////////////////////////////////

//////////////////////////////////////////cvmGet//////////////////////////////////////////////////
float  cvmGet(const CvMat2* mat, int row, int col)
{
	int type;
	type = CV_MAT_TYPE(mat->type);
	if (type == CV_32FC1)
		return ((float*)(void*)(mat->data.ptr + (size)mat->step*row))[col];
}
/////////////////////////////////////////cvmGet//////////////////////////////////////////////////////

////////////////////////////////////////myMatcopy//////////////////////////////////////////////////////
void  myMatcopy(CvMat2 *A, CvMat2 *B)
{
	int i = A->rows;
	int j = A->cols;
	int m = 0, n = 0;
	for (m = 0; m < i; m++)
		for (n = 0; n < j; n++)
			cvmSet2(B, m, n, cvmGet(A, m, n));
}
////////////////////////////////////////myMatcopy//////////////////////////////////////////////////////


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
	img->imageData = (char*)cvAlloc((size)img->imageSize);
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

//////////////////////////////////////////cvCloneImage////////////////////////////////////////////////////
void cvCloneImage(IplImage2* src,IplImage2* dst)
{
	int i, j;
	int height = src->height;
	int width = src->width;
	//int channels = src->nChannels;
	int step = src->widthStep;
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {		
			dst->imageData[i*step + j ]
					= src->imageData[i*step + j];
		}
	}
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



void cvReleaseImageHeader(IplImage2** image)
{
	if (*image)
	{
		free(*image);
		*image = 0;
	}
}

void cvReleaseData(IplImage2* img)
{
	free(img->imageData);
	img->imageData = 0;	
}

void ReleaseImage(IplImage2 ** image)
{
	if (*image)
	{
		IplImage2* img = *image;
		*image = 0;
		cvReleaseData(img);
		cvReleaseImageHeader(&img);
	}
}


void  cvDecRefData(CvMat2* mat)
{
	    free(mat->data.ptr);
		mat->data.ptr = NULL;
}

void ReleaseMat(CvMat2** array)
{
	if (*array)
	{
		CvMat2* arr = *array;
		*array = 0;
		cvDecRefData(arr);
		free(arr);
		arr = 0;
	}
}