
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <opencv2/opencv.hpp>  
#include <opencv2/highgui/highgui.hpp>  
#include <iostream>  
#include "cv.h"
#include "cxcore.h"
#include "highgui.h"
#include "utility.h"
#include "corner.h"
#include "correspEstimation.h"
#include "homographyEst.h"
#include"fast.h"
#include"goodfeatures.h"
using namespace cv;
void Result1(IplImage *inImage1, IplImage *inImage2,
			 IplImage *cornerMap1, IplImage *cornerMap2);
void Result2(IplImage *inImage1, IplImage *inImage2, CorspMap *corspMap,
			 char *outputName);
void Result3(IplImage *inImage1, IplImage *inImage2, CvMat *H,
			 char *outputName);
void Result33(IplImage *inImage1, IplImage *inImage2, int dx, int dy,
	char *outputName);
void Result444(IplImage *inImage1, IplImage *inImage2, IplImage *cornerMap1, IplImage *cornerMap2,
	 CorspMap *corspMap);
void MarkCornerPoints(IplImage *image, IplImage *cornerMap);
void DrawCorrespLine(IplImage *domainImage, IplImage *rangeImage,
					 IplImage *outImage, CorspMap *corspMap);
void Sharpen(IplImage *image, IplImage *result)
{
	//result.create(image.size(), image.type());
	int step = image->widthStep;
	int step1 = result->widthStep;
	//处理边界内部的像素点, 图像最外围的像素点应该额外处理
	for (int j = 1; j<image->height - 1; j++)
	{
		//前一行像素点
		uchar* previous = (uchar*)image->imageData + (j - 1)*step;
		//待处理的当前行
		uchar* current = (uchar*)image->imageData + j*step;
		//下一行
		uchar* next = (uchar*)image->imageData + (j + 1)*step;
		uchar * output = (uchar*)result->imageData + j*step1;
		for (int i = 1; i<(image->width - 1); i++)
		{
			//输出图像的遍历指针与当前行的指针同步递增, 以每行的每一个像素点的每一个通道值为一个递增量, 因为要考虑到图像的通道数
			int out1 = (5 * current[i * 3] - current[(i - 1) * 3] - current[(i + 1) * 3] - previous[i * 3] - next[i * 3]);
			int out2 = (5 * current[i * 3 + 1] - current[(i - 1) * 3 + 1] - current[(i + 1) * 3 + 1] - previous[i * 3 + 1] - next[i * 3 + 1]);
			int out3 = (5 * current[i * 3 + 2] - current[(i - 1) * 3 + 2] - current[(i + 1) * 3 + 2] - previous[i * 3 + 2] - next[i * 3 + 2]);
			if (out1 < 0)
				out1 = 0;
			else if (out1 >255)
				out1 = 255;//防止数据溢出
			output[i * 3] = (uchar)out1;

			if (out2 < 0)
				out2 = 0;
			else if (out2 >255)
				out2 = 255;
			output[i * 3 + 1] = (uchar)out2;

			if (out3 < 0)
				out3 = 0;
			else if (out3 >255)
				out3 = 255;
			output[i * 3] = (uchar)out3;

		}
	}
	//处理边界, 外围像素点设为 0
	for (int i = 0; i < image->width; i++)
	{
		image->imageData[i] = 0;
		(image->imageData + (image->height - 1)*image->widthStep)[i] = 0;
	}
	for (int i = 0; i < image->height; i++)
	{
		(image->imageData + i*image->widthStep)[0] = 0;
		(image->imageData + i*image->widthStep)[image->height - 1] = 0;
	}
	/*	result.row(0).setTo(Scalar(0, 0, 0));
	result.row(result.rows - 1).setTo(Scalar(0, 0, 0));
	result.col(0).setTo(Scalar(0, 0, 0));
	result.col(result.cols - 1).setTo(Scalar(0, 0, 0));*/
}

int main(int argc, char **argv) {
	CvCapture* cap = cvCreateFileCapture("test.mp4");//72
	int  height= (int)cvGetCaptureProperty(cap, CV_CAP_PROP_FRAME_HEIGHT);
	int width = (int)cvGetCaptureProperty(cap, CV_CAP_PROP_FRAME_WIDTH);
	IplImage *cornerMap1 = 0, *cornerMap2 = 0;
	CorspMap corspMap, inlierMap;
	CvMat *H;
/*	MRF  mrf_x, mrf_y;
	InitializeMRFilter(2, 5, &mrf_x);
	InitializeMRFilter(2, 5, &mrf_y);*/
	xy *x_y1, *x_y2;
	int num1, num2;
	IplImage *inImage1 = 0, *inImage2 = 0;
	IplImage *Image1 = 0, *Image2 = 0;
	Image1 = cvQueryFrame(cap);
	inImage1 = cvCreateImage(cvSize(Image1->width, Image1->height), IPL_DEPTH_8U, 1);
	cvCvtColor(Image1, inImage1, CV_RGB2GRAY);
	//	Image2 = cvQueryFrame(cap);
	Image2 = Image1;
	inImage2 = cvCreateImage(cvSize(Image2->width, Image2->height), IPL_DEPTH_8U, 1);
	cvCvtColor(Image2, inImage2, CV_RGB2GRAY);
	IplImage*srcimg= cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 1);
	//srcimg = Image1;
	
	CvVideoWriter *outputVideo;
	outputVideo = cvCreateVideoWriter("camera.avi", CV_FOURCC('M', 'J', 'P', 'G'), 20, cvGetSize(srcimg));
	int max_frames = cvGetCaptureProperty(cap, CV_CAP_PROP_FRAME_COUNT);
	int k = 1;
	while (k<10) {
	cornerMap1 = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 1);
	cornerMap2 = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 1);

	int Step = cornerMap1->widthStep;
/*	int threshold = 10000;	
	HarrisCornerDectect(inImage1, cornerMap1, threshold);
	HarrisCornerDectect(inImage2, cornerMap2, threshold);*/
	//图像锐化 
	IplImage *inImage11 = cvCreateImage(cvSize(inImage1->width, inImage1->height), IPL_DEPTH_8U, 1);
	IplImage *inImage22 = cvCreateImage(cvSize(inImage2->width, inImage2->height), IPL_DEPTH_8U, 1);
	Sharpen(inImage1, inImage11);
	Sharpen(inImage2, inImage22);

	x_y1 = fast_corner_detect((unsigned char*)inImage11->imageData, width, height, BARRIER, &num1);
	x_y2 = fast_corner_detect((unsigned char*)inImage22->imageData, width, height, BARRIER, &num2);
	uchar * cornerMap1Data, *cornerMap2Data;
	cornerMap1Data = (uchar *)cornerMap1->imageData;
	cornerMap2Data = (uchar *)cornerMap2->imageData;
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
	Result1(inImage1, inImage2, cornerMap1, cornerMap2);
	
	InitializeCorspMap(&corspMap); // initialize correspondence map
	// estimate correspondences
	// (inImage1 : domain image, inImage2 : range image)
	CorrespEstimation(inImage1, inImage2, cornerMap1, cornerMap2, &corspMap);
	Result2(inImage1, inImage2, &corspMap, "result_step2.jpg");
	printf("Number of Correspondences = %d\n", corspMap.len);
	
	H = cvCreateMat(3, 3, CV_32FC1);
	InitializeCorspMap(&inlierMap); // initialize inliers map
	RansacHomograhyEstimation(&corspMap, &inlierMap, H);
	// write and plot results
	Result2(inImage1, inImage2, &inlierMap, "result_step3_inliers.jpg");
	//Result3(inImage1, inImage2, H, "result_step3.jpg");
	printf("Number of inliers = %d\n", inlierMap.len);
	
	HomograhyEstimation(&inlierMap, H);
	Mat T(H->rows, H->cols, CV_32F, H->data.fl);
	cout << T << endl;
	float s = cvmGet(H, 2, 2);
	float dx = cvmGet(H, 0, 2) / s;
	float dy = cvmGet(H, 1, 2) / s;
	Result33(inImage1, inImage1, (int)dx, (int)dy, "result_step4.jpg");
//	Result3(inImage1, inImage2, H, "result_step4.jpg");
	// release the images and matrix
	//cvReleaseImage(&inImage1);
	//cvReleaseImage(&inImage2);
	//cvReleaseMat(&H);
	//4

//	Result444(inImage1,  inImage2, cornerMap1, cornerMap2, &corspMap);
//	TransformImage(Image1, srcimg,(int)dx,(int) dy);
	inImage1 = cvCloneImage(inImage2);
	Image1 = cvCloneImage(Image2);
	Image2 = cvQueryFrame(cap);
	cvCvtColor(Image2, inImage2, CV_RGB2GRAY);
	printf("Frame: %d / %d \n", k, max_frames);
	k++;
	}
	return 0;
}
//
// function : Result1
// usage : Result1(inImage1, inImage2, cornerMap1, cornerMap2);
// ----------------------------------------------------------------
// This function writes and plots the results of step 1
// (interest points).
//
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
//
// function : Result2
// usage : Result2(inImage1, inImage2, corspMap, name);
// -----------------------------------------------------
// This function writes and plots the results
// (putative correspondences).
//
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
//
// function : Result3
// usage : Result3(inImage1, inImage2, H, name);
// -----------------------------------------------------
// This function writes and plots the transformed image
// results.
//
void Result3(IplImage *inImage1, IplImage *inImage2, CvMat *H,
			 char *outputName)
{
	IplImage *outImage1 = 0;
	int height, width, channels;
	// get the input image data
	height = inImage1->height;
	width = inImage1->width;
	channels = inImage1->nChannels;
	outImage1 = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, channels);
	// transform inImage1 using H
	TransformImage(inImage1, outImage1, H);
	// create the output images : 2 images in the same output image
	IplImage * resultImage;
	resultImage = cvCreateImage(cvSize(width * 2, height), IPL_DEPTH_8U, channels);
	CombineTwoImages(outImage1, inImage2, resultImage);
	// display the result image
	cvNamedWindow("output image", CV_WINDOW_AUTOSIZE);
	cvShowImage("output image", resultImage);
	cvWaitKey(0);
	cvDestroyWindow("output image");
	// write output image
	WriteImage(resultImage, outputName);
	char name[80];
	sprintf_s(name, "transformed_%s", outputName);
	WriteImage(outImage1, name);
	cvReleaseImage(&resultImage);
	cvReleaseImage(&outImage1);
}



	void Result444(IplImage *inImage1, IplImage *inImage2, IplImage *cornerMap1, IplImage *cornerMap2,
		CorspMap *corspMap)
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

			int a=1,b=0,i, j,x1,y1,x2,y2;
			
			
			uchar *cornerMapData3 = 0, *cornerMapData4 = 0;
			
			int mapStep = cornerMap1->widthStep;
			cornerMapData3 = (uchar *)cornerMap1->imageData;
			cornerMapData4 = (uchar *)cornerMap2->imageData;

			for (i = 0; i < corspMap->len; i++) 
			{
				x1 = corspMap->domainImagePositionJ[i];
				y1 = corspMap->domainImagePositionI[i];
					cornerMapData3[y1*mapStep + x1] = false;
				
				x2 = corspMap->rangeImagePositionJ[i];
				y2 = corspMap->rangeImagePositionI[i] ;
					cornerMapData4[y2*mapStep + x2] = false;

		}


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
		WriteImage(outImage3, "result_step444.jpg");
		CombineTwoImages(inImage1, inImage2, outImage3);
		WriteImage(outImage3, "inputs.jpg");
		cvReleaseImage(&outImage1);
		cvReleaseImage(&outImage2);
		cvReleaseImage(&outImage3);
	}



//
// function : MarkCornerPoints
// usage : MarkCornerPoints(image, cornerMap);
// -------------------------------------------------
// This function draws marks in input image corresponding
// to the corner map.
//
void MarkCornerPoints(IplImage *image, IplImage *cornerMap) {
	int i, j;
	uchar *cornerMapData = 0;
	int height = cornerMap->height;
	int width = cornerMap->width;
	int mapStep = cornerMap->widthStep;
	cornerMapData = (uchar *)cornerMap->imageData;
	for(i = 0; i < height; i++){
		for(j = 0; j < width; j++){
			if(cornerMapData[i*mapStep + j] == true){
				cvCircle(image, cvPoint(j, i), 2, cvScalar(255,255 , 255), 2);
			}
		}
	}
}
//
// function : DrawCorrespLine
// usage : DrawCorrespLine(image1, image2, outImage, corspMap);
// -------------------------------------------------------------
// This function draws lines for correspondences of image 1 and 2.
//
void DrawCorrespLine(IplImage *domainImage, IplImage *rangeImage,
					 IplImage *outImage, CorspMap *corspMap)
{
	int a, b, i;
	uchar *outImageData = 0, *domainImageData = 0, *rangeImageData = 0;
	CvPoint rangePos;
	CvPoint domainPos;
	int height = rangeImage->height;
	int width = rangeImage->width;
	if(height == outImage->height && width * 2 == outImage->width){
		a = 1; b = 0;
	}else if(height * 2 == outImage->height && width == outImage->width){
		a = 0; b = 1;
	}else{
		printf("Error\n");
		exit(0);
	}
	outImageData = (uchar *)outImage->imageData;
	domainImageData = (uchar *)domainImage->imageData;
	rangeImageData = (uchar *)rangeImage->imageData;
	// create output image that contains ref & test images
	CombineTwoImages(domainImage, rangeImage, outImage);
	// draw correspondence lines & corner points
	for(i = 0; i < corspMap->len; i++){
		rangePos = cvPoint(corspMap->rangeImagePositionJ[i] + width * a,
			corspMap->rangeImagePositionI[i] + height * b);
		domainPos = cvPoint(corspMap->domainImagePositionJ[i],
			corspMap->domainImagePositionI[i]);
		cvCircle(outImage, domainPos, 2, cvScalar(255,255, 255), 2);
		cvCircle(outImage, rangePos, 2, cvScalar(255, 255, 255), 2);
		cvLine(outImage, domainPos, rangePos, cvScalar(0,0,0), 1);
	}
}
void Result33(IplImage *inImage1, IplImage *inImage2, int dx, int dy,
	char *outputName)
{
	IplImage *outImage1 = 0;
	int height, width, channels;
	// get the input image data
	height = inImage1->height;
	width = inImage1->width;
	channels = inImage1->nChannels;
	outImage1 = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, channels);
	// transform inImage1 using H
	TransformImage1(inImage1, outImage1, dx, dy);
	// create the output images : 2 images in the same output image
	IplImage * resultImage;
	resultImage = cvCreateImage(cvSize(width * 2, height), IPL_DEPTH_8U, channels);
	CombineTwoImages(outImage1, inImage2, resultImage);
	// display the result image
	cvNamedWindow("output image", CV_WINDOW_AUTOSIZE);
	cvShowImage("output image", resultImage);
	cvWaitKey(0);
	cvDestroyWindow("output image");
	// write output image
	WriteImage(resultImage, outputName);
	char name[80];
	sprintf_s(name, "transformed_%s", outputName);
	WriteImage(outImage1, name);
	cvReleaseImage(&resultImage);
	cvReleaseImage(&outImage1);
}
