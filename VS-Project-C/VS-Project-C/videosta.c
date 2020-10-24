#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//#include "point.h"
#include "defs.h"
#include "utility.h"
#include "filter.h"
#include "svd.h"
#include "correspEstimation.h"
#include "homographyEst.h"
#include "fast.h"
#include "videosta.h"
//int k;
 IplImage2 *inImage1 = 0;
 IplImage2 *inImage2 = 0;
 MRF  mrf_x;
 MRF  mrf_y;

void Initializevideosta(int height, int width) {
//	k = 0;
	inImage1 = cvCreateImage2(cvSize2(width, height), IPL_DEPTH_8U, 1);//初始化参考帧
	inImage2 = cvCreateImage2(cvSize2(width, height), IPL_DEPTH_8U, 1);//初始化当前帧
	InitializeMRFilter(2, 5, &mrf_x);//初始化x方向滤波器
	InitializeMRFilter(2, 5, &mrf_y);//初始化y方向滤波器
}

void videosta(char *img, int height, int width,int *k)
{
	int i, j;
	//获取当前帧
	Array2Image(inImage2, img);//用IplImage2的结构进行操作
	if (*k == 0)//第一次调用时将参考帧和当前帧都设为第一帧
	{
		cvCloneImage(inImage2,inImage1);
	}
	IplImage2 *cornerMap1 = 0, *cornerMap2 = 0;
	CorspMap corspMap, inlierMap;//存匹配角点、内点
	CvMat2 *H;
	H = cvCreateMat2(3, 3, CV_32FC1);
	cornerMap1 = cvCreateImage2(cvSize2(width, height), IPL_DEPTH_8U, 1);//存参考帧的角点
	cornerMap2 = cvCreateImage2(cvSize2(width, height), IPL_DEPTH_8U, 1);//存当前帧的角点
	int Step = cornerMap1->widthStep;
	//xy *x_y1, *x_y2;
	xy x_y1[512];
	xy x_y2[512];
	int num1, num2;
	//fast角点检测
	fast_corner_detect((uchar*)inImage1->imageData, width, height, BARRIER, &num1,x_y1);//BARRIER为检测的门限，一般选用60
    fast_corner_detect((uchar*)inImage2->imageData, width, height, BARRIER, &num2,x_y2);
	uchar * cornerMap1Data, *cornerMap2Data;
	cornerMap1Data = (uchar*)cornerMap1->imageData;//参考帧的角点
	cornerMap2Data = (uchar*)cornerMap2->imageData;//当前帧的角点
	//初始化
	for (i = 0; i < height; i++)
		for (j = 0; j < width; j++) {
			cornerMap1Data[i*Step + j] = 0;
			cornerMap2Data[i*Step + j] = 0;
		}
	int p = 0, xx, yy;
	int q = 0;
	while (p <num1) {
		xx = x_y1[p].x;
		yy = x_y1[p++].y;
		cornerMap1Data[yy*Step + xx] = 1;//角点处设为true
	}
	while (q <num2) {
		xx = x_y2[q].x;
		yy = x_y2[q++].y;
		cornerMap2Data[yy*Step + xx] = 1;//角点处设为true
	}
	InitializeCorspMap(&corspMap);
	//粗匹配
	CorrespEstimation(inImage1, inImage2, cornerMap1, cornerMap2, &corspMap);
	InitializeCorspMap(&inlierMap);
	//淘汰错误匹配
	RansacHomograhyEstimation(&corspMap, &inlierMap, H);
	//计算H矩阵
	HomograhyEstimation(&inlierMap, H);
	//取出平移变量dx,dy
	float s = cvmGet(H, 2, 2);//缩放系数
	float dx = cvmGet(H, 0, 2) / s;//平移系数dx
	float dy = cvmGet(H, 1, 2) / s;//平移系数dy
	//对运动参数进行滤波
	MRFilter(dx, 1, &mrf_x);
	MRFilter(dy, 1, &mrf_y);
	dx = dx + mrf_x.x_diff; //进行真正的运动矫正
	dy = dy + mrf_y.x_diff;
	//对参考帧进行变换
	TransformImage1(inImage1, img, dx, dy);
	//当前帧为下一次的参考帧
	cvCloneImage(inImage2,inImage1);
	(*k)++;
	ReleaseImage(&cornerMap1);
	ReleaseImage(&cornerMap2);
	ReleaseMat(&H);
//	free(x_y1);
//	free(x_y2);
}

void  sopt_videosta()
{
	ReleaseImage(&inImage1);
	ReleaseImage(&inImage2);
	//free(&mrf_x);
	//free(&mrf_y);
	//free(&k);
}
