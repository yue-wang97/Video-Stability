
/*#ifdef _DEBUG

#include <stdlib.h>
#include <crtdbg.h>

#ifndef DEBUG_NEW
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define new DEBUG_NEW
#endif

#endif*/
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
#include"videosta.h"
//#include <crtdbg.h> 

//extern int k;





int main()
{
	//EnableMemLeakCheck();

	clock_t start_time = clock();
	CvCapture* cap = cvCreateFileCapture("my02.avi");
	int	height = (int)cvGetCaptureProperty(cap, CV_CAP_PROP_FRAME_HEIGHT);
	int width = (int)cvGetCaptureProperty(cap, CV_CAP_PROP_FRAME_WIDTH);
	IplImage *img1 = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 1);
	CvVideoWriter *outputVideo;
	outputVideo = cvCreateVideoWriter("camera.avi", CV_FOURCC('M', 'J', 'P', 'G'), 20, cvGetSize(img1),0);
	int max_frames = cvGetCaptureProperty(cap, CV_CAP_PROP_FRAME_COUNT);
	int kk = 0;
	int * k = &kk;
	//初始化
	Initializevideosta(height, width);
	while (*k<5) {
		IplImage *Image1;
		Image1 = cvQueryFrame(cap);
		IplImage *inImage = cvCreateImage(cvSize(Image1->width, Image1->height), IPL_DEPTH_8U, 1);
		cvCvtColor(Image1, inImage, CV_RGB2GRAY);
		//调用videosta函数
		
		videosta(inImage->imageData, height, width,k);
//		videosta(pBuf, height, width, k);

		cvWriteFrame(outputVideo, inImage);
		cvReleaseImage(&inImage);
		printf("Frame: %d / %d \n",*k,max_frames);
	}
	cvReleaseVideoWriter(&outputVideo);
	cvReleaseImage(&img1);
	clock_t end_time = clock();
	printf( "Running time is:%f ms\n", (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000);
	sopt_videosta();
	getchar();
	
	return 0;
//	_CrtDumpMemoryLeaks();
	
}
