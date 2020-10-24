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
	inImage1 = cvCreateImage2(cvSize2(width, height), IPL_DEPTH_8U, 1);//��ʼ���ο�֡
	inImage2 = cvCreateImage2(cvSize2(width, height), IPL_DEPTH_8U, 1);//��ʼ����ǰ֡
	InitializeMRFilter(2, 5, &mrf_x);//��ʼ��x�����˲���
	InitializeMRFilter(2, 5, &mrf_y);//��ʼ��y�����˲���
}

void videosta(char *img, int height, int width,int *k)
{
	int i, j;
	//��ȡ��ǰ֡
	Array2Image(inImage2, img);//��IplImage2�Ľṹ���в���
	if (*k == 0)//��һ�ε���ʱ���ο�֡�͵�ǰ֡����Ϊ��һ֡
	{
		cvCloneImage(inImage2,inImage1);
	}
	IplImage2 *cornerMap1 = 0, *cornerMap2 = 0;
	CorspMap corspMap, inlierMap;//��ƥ��ǵ㡢�ڵ�
	CvMat2 *H;
	H = cvCreateMat2(3, 3, CV_32FC1);
	cornerMap1 = cvCreateImage2(cvSize2(width, height), IPL_DEPTH_8U, 1);//��ο�֡�Ľǵ�
	cornerMap2 = cvCreateImage2(cvSize2(width, height), IPL_DEPTH_8U, 1);//�浱ǰ֡�Ľǵ�
	int Step = cornerMap1->widthStep;
	//xy *x_y1, *x_y2;
	xy x_y1[512];
	xy x_y2[512];
	int num1, num2;
	//fast�ǵ���
	fast_corner_detect((uchar*)inImage1->imageData, width, height, BARRIER, &num1,x_y1);//BARRIERΪ�������ޣ�һ��ѡ��60
    fast_corner_detect((uchar*)inImage2->imageData, width, height, BARRIER, &num2,x_y2);
	uchar * cornerMap1Data, *cornerMap2Data;
	cornerMap1Data = (uchar*)cornerMap1->imageData;//�ο�֡�Ľǵ�
	cornerMap2Data = (uchar*)cornerMap2->imageData;//��ǰ֡�Ľǵ�
	//��ʼ��
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
		cornerMap1Data[yy*Step + xx] = 1;//�ǵ㴦��Ϊtrue
	}
	while (q <num2) {
		xx = x_y2[q].x;
		yy = x_y2[q++].y;
		cornerMap2Data[yy*Step + xx] = 1;//�ǵ㴦��Ϊtrue
	}
	InitializeCorspMap(&corspMap);
	//��ƥ��
	CorrespEstimation(inImage1, inImage2, cornerMap1, cornerMap2, &corspMap);
	InitializeCorspMap(&inlierMap);
	//��̭����ƥ��
	RansacHomograhyEstimation(&corspMap, &inlierMap, H);
	//����H����
	HomograhyEstimation(&inlierMap, H);
	//ȡ��ƽ�Ʊ���dx,dy
	float s = cvmGet(H, 2, 2);//����ϵ��
	float dx = cvmGet(H, 0, 2) / s;//ƽ��ϵ��dx
	float dy = cvmGet(H, 1, 2) / s;//ƽ��ϵ��dy
	//���˶����������˲�
	MRFilter(dx, 1, &mrf_x);
	MRFilter(dy, 1, &mrf_y);
	dx = dx + mrf_x.x_diff; //�����������˶�����
	dy = dy + mrf_y.x_diff;
	//�Բο�֡���б任
	TransformImage1(inImage1, img, dx, dy);
	//��ǰ֡Ϊ��һ�εĲο�֡
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
