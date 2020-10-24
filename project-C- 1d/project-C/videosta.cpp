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
#include"H.h"
#include"filter.h"
#include"videosta.h"
int k ;
IplImage *inImage1=0;
IplImage *inImage2=0;
MRF  mrf_x;
MRF  mrf_y;

void Initializevideosta(int height,int width) {
	 k = 0;
	inImage1 = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 1);
	inImage2 = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 1);
	//MRF  mrf_x;
	//MRF  mrf_y;
	InitializeMRFilter(2, 5, &mrf_x);
	InitializeMRFilter(2, 5, &mrf_y);
}

void videosta(char *img, int height, int width)
{
	Array2Image(inImage2, img);
	if (k == 0)
	{
		inImage1 = cvCloneImage(inImage2);
	}
	IplImage *cornerMap1 = 0, *cornerMap2 = 0;
	CorspMap corspMap, inlierMap;
	float H[9];
	cornerMap1 = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 1);
	cornerMap2 = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 1);
	int Step = cornerMap1->widthStep;
	xy *x_y1, *x_y2;
	int num1, num2;
	x_y1 = fast_corner_detect((uchar*)inImage1->imageData, width, height, BARRIER, &num1);
	x_y2 = fast_corner_detect((uchar*)inImage2->imageData, width, height, BARRIER, &num2);
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
	CorrespEstimation(inImage1, inImage2, cornerMap1, cornerMap2, &corspMap);
	//H = cvCreateMat2(3, 3, CV_32FC1);
	InitializeCorspMap(&inlierMap); 
	RansacHomograhyEstimation(&corspMap, &inlierMap, H);
	HomograhyEstimation(&inlierMap, H);
	float s = GetH(H, 2, 2);
	float dx = GetH(H, 0, 2) / s;
	float dy = GetH(H, 1, 2) / s;
/*	if (fabs(dx) < 20) {
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
	}*/
	MRFilter(dx, 1, &mrf_x);
	dx = dx + mrf_x.x_diff; //进行真正的运动矫正
	MRFilter(dy, 1, &mrf_y);
	dy = dy + mrf_y.x_diff;
	TransformImage1(inImage1, img, dx, dy);
	inImage1 = cvCloneImage(inImage2);
	k++;
	//free(cornerMap1);
	//free(cornerMap2);
}

void  sopt_videosta()
{
	free(inImage1);
	free(inImage2);
}
