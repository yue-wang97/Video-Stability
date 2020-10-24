

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

void MakeImageBlock(IplImage2 *block, IplImage2 *image,
					int centPosI, int centPosJ);
void InitializeimageData(char * img);
void TransformImage1(IplImage2 *inImage, char * outData, int dx, int dy);
void TransformImage(IplImage2 *inImage, IplImage2 *outImage, int dx, int dy);
