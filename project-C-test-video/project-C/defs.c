#include"defs.h"

///////////////////////////////////////////////////////CvMat/////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////cvCreateMat///////////////////////////////////////////////////////

static void icvCheckHuge(CvMat* arr)
{
	if ((int64)arr->step*arr->rows > INT_MAX)
		arr->type &= ~CV_MAT_CONT_FLAG;
}
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


CvMat* cvCreateMatHeader(int rows, int cols, int type)
{
	type = CV_MAT_TYPE(type);

	
	int min_step = CV_ELEM_SIZE(type);
	
	min_step *= cols;

	CvMat* arr = (CvMat*)cvAlloc(sizeof(*arr));

	arr->step = min_step;
	arr->type = CV_MAT_MAGIC_VAL | type | CV_MAT_CONT_FLAG;
	arr->rows = rows;
	arr->cols = cols;
	arr->data.ptr = 0;
	arr->refcount = 0;
	arr->hdr_refcount = 1;

	icvCheckHuge(arr);
	return arr;
}

void cvCreateData(CvMat* arr)
{
	size_t step, total_size;
	CvMat* mat = (CvMat*)arr;
	step = mat->step;

	if (mat->rows == 0 || mat->cols == 0)
		return;

	if (step == 0)
		step = CV_ELEM_SIZE(mat->type)*mat->cols;

	int64 _total_size = (int64)step*mat->rows + sizeof(int) + CV_MALLOC_ALIGN;
	total_size = (size_t)_total_size;
	mat->refcount = (int*)cvAlloc((size_t)total_size);
	mat->data.ptr = (uchar*)cvAlignPtr(mat->refcount + 1, CV_MALLOC_ALIGN);
	*mat->refcount = 1;
}

CvMat* cvCreateMat2(int height, int width, int type)
{
	CvMat* arr = cvCreateMatHeader(height, width, type);
	cvCreateData(arr);

	return arr;
}

/////////////////////////////////////////////////////cvCreateMat///////////////////////////////////////////////////////////

////////////////////////////////////////cvmSet///////////////////////////////////////////////////
void  cvmSet2(CvMat* mat, int row, int col, double value)
{
	int type;
	type = CV_MAT_TYPE(mat->type);
	if (type == CV_32FC1)
		((float*)(void*)(mat->data.ptr + (size_t)mat->step*row))[col] = (float)value;
	else
	{
		((double*)(void*)(mat->data.ptr + (size_t)mat->step*row))[col] = value;
	}
}
////////////////////////////////////////cvmSet////////////////////////////////////////////////////

//////////////////////////////////////////cvmGet//////////////////////////////////////////////////
double  cvmGet(const CvMat* mat, int row, int col)
{
	int type;
	type = CV_MAT_TYPE(mat->type);
	if (type == CV_32FC1)
		return ((float*)(void*)(mat->data.ptr + (size_t)mat->step*row))[col];
}
/////////////////////////////////////////cvmGet//////////////////////////////////////////////////////

////////////////////////////////////////myMatcopy//////////////////////////////////////////////////////
void  myMatcopy(CvMat *A, CvMat *B)
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
IplImage* cvInitImageHeader(IplImage * image, CvSize2 size, int depth,
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
	const int64 imageSize_tmp = (int64)image->widthStep*(int64)image->height;
	image->imageSize = (int)imageSize_tmp;
	return image;
}

void cvCreateData_img(IplImage* arr)
{
	IplImage* img = (IplImage*)arr;
	img->imageData = (char*)cvAlloc((size_t)img->imageSize);		
}

IplImage * cvCreateImageHeader(CvSize2 size, int depth, int channels)
{
	IplImage *img = 0;

	img = (IplImage *)cvAlloc(sizeof(*img));
	cvInitImageHeader(img, size, depth, channels);
	return img;
}

IplImage * cvCreateImage2(CvSize2 size, int depth, int channels)
{
	IplImage *img = cvCreateImageHeader(size, depth, channels);
	cvCreateData_img(img);

	return img;
}
/////////////////////////////////////cvCreateImage////////////////////////////////////////////////////

//////////////////////////////////////InitializeImage////////////////////////////////////////////////////
void InitializeImage(IplImage *image)
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
void   Array2Image(IplImage *img, uchar *a)
{
	int i, j, k;
	int height = img->height;
	int width = img->width;
	int channels = img->nChannels;
	int step = img->widthStep;
	uchar *imageData = img->imageData;
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
IplImage* cvCloneImage(IplImage* src)
{
	IplImage* dst = 0;
	dst = (IplImage*)cvAlloc(sizeof(*dst));
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