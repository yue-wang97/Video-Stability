#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include<time.h>
#include <opencv2/opencv.hpp>  
#include <opencv2/highgui/highgui.hpp>  
#include <iostream>  
#include "cv.h"
#include "cxcore.h"
#include "highgui.h"
#include"utility.h"
//#include"defs.h"
#include"fast.h"
#include"correspEstimation.h"
#include"utility.h"
#include"homographyEst.h"
#include"filter.h"
#include"videosta.h"

extern int k;


void Result1(IplImage *inImage1, IplImage *inImage2,
	IplImage *cornerMap1, IplImage *cornerMap2);
void Result2(IplImage *inImage1, IplImage *inImage2, CorspMap *corspMap,
	char *outputName);

using namespace std;
using namespace cv;


int main() {
	clock_t start_time = clock();
	CvCapture* cap = cvCreateFileCapture("my02.avi");//72
	int	height = (int)cvGetCaptureProperty(cap, CV_CAP_PROP_FRAME_HEIGHT);
	int width = (int)cvGetCaptureProperty(cap, CV_CAP_PROP_FRAME_WIDTH);
	IplImage*srcimg = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 1);
	CvVideoWriter *outputVideo;
	outputVideo = cvCreateVideoWriter("camera.avi", CV_FOURCC('M', 'J', 'P', 'G'), 20, cvGetSize(srcimg), 0);
	int max_frames = cvGetCaptureProperty(cap, CV_CAP_PROP_FRAME_COUNT);
	Initializevideosta(height, width);
	while (k<5) {
		IplImage *inImage = 0;
		IplImage *Image1 = 0;
		Image1 = cvQueryFrame(cap);
		inImage = cvCreateImage(cvSize(Image1->width, Image1->height), IPL_DEPTH_8U, 1);
		cvCvtColor(Image1, inImage, CV_RGB2GRAY);
		videosta(inImage->imageData, height, width);
		//Array2Image(srcimg, src);
		cvWriteFrame(outputVideo, inImage);
		printf("Frame: %d / %d \n", k, max_frames);
	}
	cvReleaseVideoWriter(&outputVideo);
	clock_t end_time = clock();
	printf("Running time is:%f ms\n", (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000);
	getchar();
	return 0;
}






void MarkCornerPoints(IplImage *image, IplImage *cornerMap) {
	int i, j;
	uchar *cornerMapData = 0;
	int height = cornerMap->height;
	int width = cornerMap->width;
	int mapStep = cornerMap->widthStep;
	cornerMapData = (uchar *)cornerMap->imageData;
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			if (cornerMapData[i*mapStep + j] == true) {
				cvCircle(image, cvPoint(j, i), 2, cvScalar(255, 255, 255), 2);
			}
		}
	}
}
void CombineTwoImages(IplImage *image1, IplImage *image2,
	IplImage *outImage)
{
	int i, j, k;
	uchar *outImageData = 0, *image1Data = 0, *image2Data = 0;
	int height = image1->height;
	int width = image1->width;
	int step = image1->widthStep;
	int channels = image1->nChannels;
	int outWidth = outImage->width;
	int outHeight = outImage->height;
	int outStep = outImage->widthStep;
	if (outWidth == width * 2 && outHeight == height) {
	}
	else if (outWidth == width && outHeight == height * 2) {
	}
	else {
		printf("image combining error\n");
		exit(0);
	}
	outImageData = (uchar *)outImage->imageData;
	image1Data = (uchar *)image1->imageData;
	image2Data = (uchar *)image2->imageData;
	for (i = 0; i < outHeight; i++) {
		for (j = 0; j < outWidth; j++) {

			for (k = 0; k < channels; k++) {
				if (i < height && j < width) {
					outImageData[i*outStep + j*channels + k]
						= image1Data[i*step + j*channels + k];
				}
				else if ((i >= height && j < width)) {
					outImageData[i*outStep + j*channels + k]
						= image2Data[(i - height)*step + j*channels + k];
				}
				else if ((i < height && j >= width)) {
					outImageData[i*outStep + j*channels + k]
						= image2Data[i*step + (j - width)*channels + k];
				}
				else {
					printf("there is no i > height & j > width \n");
					exit(0);
				}
			}
		}
	}
}
void WriteImage(IplImage *image, char *imageName) {
	if (!cvSaveImage(imageName, image)) {
		printf("Could not save: %s\n", imageName);
	}
}


void Result1(IplImage *inImage1, IplImage *inImage2,
	IplImage *cornerMap1, IplImage *cornerMap2)
{
	IplImage *outImage1 = 0, *outImage2 = 0, *outImage3 = 0;
	int height, width, channels;
	// create the output images

	height = inImage1->height;
	width = inImage1->width;
	channels = inImage1->nChannels;
	outImage1 = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, channels);
	outImage2 = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, channels);
	// draw a circle on the corner point
	cvCopy(inImage1, outImage1);
	cvCopy(inImage2, outImage2);
	MarkCornerPoints(outImage1, cornerMap1);
	MarkCornerPoints(outImage2, cornerMap2);
	// create the output images : 2 images in the same output image
	outImage3 = cvCreateImage(cvSize(width * 2, height), IPL_DEPTH_8U, channels);
	CombineTwoImages(outImage1, outImage2, outImage3);
	// display the result image
	cvNamedWindow("output image", CV_WINDOW_AUTOSIZE);
	cvShowImage("output image", outImage3);
	cvWaitKey(0);
	cvDestroyWindow("output image");
	// write output image
	WriteImage(outImage1, "result_corner1.jpg");
	WriteImage(outImage2, "result_corner2.jpg");
	WriteImage(outImage3, "result_step1.jpg");
	CombineTwoImages(inImage1, inImage2, outImage3);
	WriteImage(outImage3, "inputs.jpg");
	cvReleaseImage(&outImage1);
	cvReleaseImage(&outImage2);
	cvReleaseImage(&outImage3);
}

void DrawCorrespLine(IplImage *domainImage, IplImage *rangeImage,
	IplImage *outImage, CorspMap *corspMap)
{
	int a, b, i;
	uchar *outImageData = 0, *domainImageData = 0, *rangeImageData = 0;
	CvPoint rangePos;
	CvPoint domainPos;
	int height = rangeImage->height;
	int width = rangeImage->width;
	if (height == outImage->height && width * 2 == outImage->width) {
		a = 1; b = 0;
	}
	else if (height * 2 == outImage->height && width == outImage->width) {
		a = 0; b = 1;
	}
	else {
		printf("Error\n");
		exit(0);
	}
	outImageData = (uchar *)outImage->imageData;
	domainImageData = (uchar *)domainImage->imageData;
	rangeImageData = (uchar *)rangeImage->imageData;
	// create output image that contains ref & test images
	CombineTwoImages(domainImage, rangeImage, outImage);
	// draw correspondence lines & corner points
	for (i = 0; i < corspMap->len; i++) {
		rangePos = cvPoint(corspMap->rangeImagePositionJ[i] + width * a,
			corspMap->rangeImagePositionI[i] + height * b);
		domainPos = cvPoint(corspMap->domainImagePositionJ[i],
			corspMap->domainImagePositionI[i]);
		cvCircle(outImage, domainPos, 2, cvScalar(255, 255, 255), 2);
		cvCircle(outImage, rangePos, 2, cvScalar(255, 255, 255), 2);
		cvLine(outImage, domainPos, rangePos, cvScalar(0, 0, 0), 1);
	}
}

void Result2(IplImage *inImage1, IplImage *inImage2, CorspMap *corspMap,
	char *outputName)
{
	IplImage *outImage = 0;
	int height, width, channels;
	// create the output images
	height = inImage1->height;
	width = inImage1->width;
	channels = inImage1->nChannels;
	// create the output images : 2 images in the same output image
	outImage = cvCreateImage(cvSize(width * 2, height), IPL_DEPTH_8U, channels);
	// draw correspondence lines on the resultant images
	DrawCorrespLine(inImage1, inImage2, outImage, corspMap);
	// display the result image
	cvNamedWindow("output image", CV_WINDOW_AUTOSIZE);
	cvShowImage("output image", outImage);
	cvWaitKey(0);
	cvDestroyWindow("output image");
	// write output image
	WriteImage(outImage, outputName);
	cvReleaseImage(&outImage);
}


/*int main() {
	clock_t start_time = clock();
	CvCapture* cap = cvCreateFileCapture("car1.avi");//72
	int	height = (int)cvGetCaptureProperty(cap, CV_CAP_PROP_FRAME_HEIGHT);
	int width = (int)cvGetCaptureProperty(cap, CV_CAP_PROP_FRAME_WIDTH);
	IplImage *cornerMap1 = 0, *cornerMap2 = 0;
	CorspMap corspMap, inlierMap;
	CvMat *H;
	MRF  mrf_x, mrf_y;
	InitializeMRFilter(2, 5, &mrf_x);
	InitializeMRFilter(2, 5, &mrf_y);
	xy *x_y1, *x_y2;
	int num1, num2;
	IplImage *inImagea = 0, *inImageb = 0;
	IplImage *Image1 = 0, *Image2 = 0;
	Image1 = cvQueryFrame(cap);
	inImagea = cvCreateImage(cvSize(Image1->width, Image1->height), IPL_DEPTH_8U, 1);
	cvCvtColor(Image1, inImagea, CV_RGB2GRAY);
//	Image2 = cvQueryFrame(cap);
	Image2 = Image1;
	inImageb = cvCreateImage(cvSize(Image2->width, Image2->height), IPL_DEPTH_8U, 1);
	cvCvtColor(Image2, inImageb, CV_RGB2GRAY);
	IplImage*srcimg= cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 1);
	//srcimg = Image1;
	CvVideoWriter *outputVideo;
	outputVideo = cvCreateVideoWriter("camera.avi", CV_FOURCC('M', 'J', 'P', 'G'), 20, cvGetSize(srcimg),0);
	int max_frames = cvGetCaptureProperty(cap, CV_CAP_PROP_FRAME_COUNT);
	int kk = 1;
	while (kk<201) {
		cornerMap1 = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 1);
		cornerMap2 = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 1);
		int Step = cornerMap1->widthStep;
		x_y1 = fast_corner_detect((uchar*)inImagea->imageData, width, height, BARRIER, &num1);	
		x_y2 = fast_corner_detect((uchar*)inImageb->imageData, width, height, BARRIER, &num2);
		uchar * cornerMap1Data, *cornerMap2Data;
		cornerMap1Data = (uchar*)cornerMap1->imageData;
		cornerMap2Data = (uchar*)cornerMap2->imageData;
		for (int i = 0; i < height; i++)
			for (int j = 0; j < width; j++) {
				cornerMap1Data[i*Step + j] = 0;
				cornerMap2Data[i*Step + j] = 0;
			}
		int p = 0, xx, yy;
		int q = 0;
		while (p <num1) {
			xx = x_y1[p].x;
			yy = x_y1[p++].y;
			cornerMap1Data[yy*Step + xx] = 1;
		}
		while (q <num2) {
			xx = x_y2[q].x;
			yy = x_y2[q++].y;
			cornerMap2Data[yy*Step + xx] = 1;
		}
		InitializeCorspMap(&corspMap);
		CorrespEstimation(inImagea, inImageb, cornerMap1, cornerMap2, &corspMap);
		H = cvCreateMat(3, 3, CV_32FC1);
		InitializeCorspMap(&inlierMap);
		RansacHomograhyEstimation(&corspMap, &inlierMap, H);
		HomograhyEstimation(&inlierMap, H);
		float s = cvmGet(H, 2, 2);
		float dx = cvmGet(H, 0, 2) / s;
		float dy = cvmGet(H, 1, 2) / s;
		if (fabs(dx) < 20) {
			MRFilter(dx, 1, &mrf_x);
			dx = dx + mrf_x.x_diff; //进行真正的运动矫正
		}
		else
		{
			dx = 0;
		}
		if (fabs(dy) < 20) {
			MRFilter(dy, 1, &mrf_y);
			dy = dy + mrf_y.x_diff;
		}
		else
		{
			dy = 0;
		}
		TransformImage(inImagea, srcimg, dx, dy);
		inImagea = cvCloneImage(inImageb);	
		Image1 = cvCloneImage(Image2);
		cvWriteFrame(outputVideo, srcimg);
		Image2 = cvQueryFrame(cap);
		cvCvtColor(Image2, inImageb, CV_RGB2GRAY);
		printf("Frame: %d / %d \n", kk, max_frames);
		kk++;
	}
	cvReleaseVideoWriter(&outputVideo);
	clock_t end_time = clock();
	printf("Running time is:%f ms\n", (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000);
	getchar();
	return 0;
}*/