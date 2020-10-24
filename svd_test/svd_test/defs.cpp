#include"defs.h"

///////////////////////////////////////////////////////CvMat/////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////cvCreateMat///////////////////////////////////////////////////////

static void icvCheckHuge(CvMat2* arr)
{
	if ((int64)arr->step*arr->rows > INT_MAX)
		arr->type &= ~CV_MAT_CONT_FLAG;
}
uchar1** alignPtr(uchar1** ptr, int n = (int)sizeof(uchar1*))
{
	return (uchar1**)(((size)ptr + n - 1) & -n);
}
void* fastMalloc(size size)
{
	uchar1* udata = (uchar1*)malloc(size + sizeof(void*) + CV_MALLOC_ALIGN);
	uchar1** adata = alignPtr((uchar1**)udata + 1, CV_MALLOC_ALIGN);
	adata[-1] = udata;
	return adata;
}
void * cvAlloc(size size)
{
	return fastMalloc(size);
}

void* cvAlignPtr(const void* ptr, int align = 32)
{
	//CV_DbgAssert((align & (align - 1)) == 0);
	return (void*)(((size)ptr + align - 1) & ~(size)(align - 1));
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
	arr->refcount = 0;
	arr->hdr_refcount = 1;

	icvCheckHuge(arr);
	return arr;
}

void cvCreateData(CvMat2* arr)
{
	size_t step, total_size;
	CvMat2* mat = (CvMat2*)arr;
	step = mat->step;

	if (mat->rows == 0 || mat->cols == 0)
		return;

	if (step == 0)
		step = CV_ELEM_SIZE(mat->type)*mat->cols;

	int64 _total_size = (int64)step*mat->rows + sizeof(int) + CV_MALLOC_ALIGN;
	total_size = (size_t)_total_size;
	mat->refcount = (int*)cvAlloc((size_t)total_size);
	mat->data.ptr = (uchar1*)cvAlignPtr(mat->refcount + 1, CV_MALLOC_ALIGN);
	*mat->refcount = 1;
}

CvMat2* cvCreateMat2(int height, int width, int type)
{
	CvMat2* arr = cvCreateMatHeader(height, width, type);
	cvCreateData(arr);

	return arr;
}

/////////////////////////////////////////////////////cvCreateMat///////////////////////////////////////////////////////////

////////////////////////////////////////cvmSet///////////////////////////////////////////////////
void  cvmSet2(CvMat2* mat, int row, int col, double value)
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
double  cvmGet(const CvMat2* mat, int row, int col)
{
	int type;
	type = CV_MAT_TYPE(mat->type);
	if (type == CV_32FC1)
		return ((float*)(void*)(mat->data.ptr + (size_t)mat->step*row))[col];
}
/////////////////////////////////////////cvmGet//////////////////////////////////////////////////////


//cvCreateImage
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
	const int64 imageSize_tmp = (int64)image->widthStep*(int64)image->height;
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


//cvSize2
CvSize2  cvSize2(int width, int height)
{
	CvSize2 s;

	s.width = width;
	s.height = height;

	return s;
}

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