//
// file : utility.h
//----------------------------------------
// this file contains utility functions to
// deal with general processes.
//

#define OFF 0
#define ON 1
#define EPS 0.5
#define IMPLEMENTATION 2
#define min(a, b) ((a <= b) ? a : b)
#define max(a, b) ((a >= b) ? a : b)
void Array2CvMat(float *arr, CvMat2 *cvArr, int row, int column);
void CvMat2Array(CvMat2 *cvArr, float *arr, int row, int column);
void MakeImageBlock(IplImage2 *block, IplImage2 *image,
					int centPosI, int centPosJ);
void TransformImage(IplImage2 *inImage, IplImage2 *outImage, int dx, int dy);
void InitializeimageData(uchar * img, int height, int width, int step);
void   Array2Image(IplImage2 *img, char *b);
void TransformImage1(IplImage2 *inImage, char * Data, int dx, int dy);
