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

void MakeImageBlock(IplImage *block, IplImage *image,
					int centPosI, int centPosJ);
void InitializeimageData(char * img);
void InitializeImage(IplImage *image);
void  Array2Image(IplImage *img, char *a);
void TransformImage1(IplImage *inImage, char * outData, int dx, int dy);
//void TransformImage(IplImage *inImage, IplImage *outImage, int dx, int dy);
