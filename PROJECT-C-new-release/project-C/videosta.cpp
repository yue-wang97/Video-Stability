#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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

/*void Sharpen(IplImage2 *image, IplImage2 *result)
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

}*/
/*void box_filter(IplImage2* img, IplImage2* result) {
	//init part
	double s;
	int width = img->width, height = img->height;
	int m_w = 5, m_h = 5;//window_size
	int boxwidth = width - m_w, boxheight = height - m_h;
	int *sum = (int*)malloc(boxwidth *boxheight * sizeof(double));
	int *buff = (int*)malloc(width * sizeof(double));
	memset(sum, 0, boxwidth *boxheight * sizeof(int));
	memset(buff, 0, width * sizeof(int));
	//set buff:from 0 to 4 rows,per col
	int x, y, j;
	for (y = 0; y<m_h; y++) {
		for (x = 0; x<width; x++) {
			uchar pixel = (uchar)img->imageData[y *width + x];
			buff[x] += pixel;
			//printf("%d:%d\n", x, buff[x]);
		}
	}
	for (y = 0; y<height - m_h; y++) {
		int Xsum = 0;
		for (j = 0; j<m_w; j++) {
			Xsum += buff[j];//sum of pixel from (0,0) to (m_h,m_w) (also x = 0)
		}

		for (x = 0; x<boxwidth; x++) {
			if (x != 0) {
				Xsum = Xsum - buff[x - 1] + buff[m_w - 1 + x];//Xsum:sum of cols range from x to x+m_w ,rows range from 0 to 4
			}
			sum[y*boxwidth + x] = (float)Xsum;
		}

		for (x = 0; x<width; x++) {
			uchar pixel = (uchar)img->imageData[y *width + x];//img[y *width + x];    
			uchar pixel2 = (uchar)img->imageData[(y + m_h) *width + x];//img[(y+mheight) *width + x];  
			buff[x] = buff[x] - pixel + pixel2;
		}
	}
	//cout << buff << endl;
	//�������õ�ÿ����ĺͣ���������result
	for (y = 0; y<height - 5; y++) {
		for (x = 0; x<width; x++) {
			if (y>m_h / 2 && y<height - m_h / 2 && x>m_w / 2 && x<width - m_w / 2) {
				s = sum[(y - m_h / 2) *boxwidth + (x - m_h / 2)] / (m_h*m_w);
				result->imageData[y *width + x]=s;
			}
			else {
				s = -1;
				result->imageData[y *width + x] = s;
			}//end else
		}//end the first for
	}//end the second for
}*/
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
	//图像锐化 
	/*IplImage2 *inImage11=cvCreateImage2(cvSize2(width, height), IPL_DEPTH_8U, 1);
	IplImage2 *inImage22 = cvCreateImage2(cvSize2(width, height), IPL_DEPTH_8U, 1);
	Sharpen(inImage1, inImage11);
	Sharpen(inImage2, inImage22);
	//图像滤波
	box_filter(inImage11, inImage11);
	box_filter(inImage22, inImage22);*/

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
	cvCloneImage(inImage2, inImage1);
	 
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
