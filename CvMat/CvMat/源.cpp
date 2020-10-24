#include"defs.h"
#include"stdio.h"
#include <iostream>
#include<opencv2\opencv.hpp>
using namespace cv;
using namespace std;

/*void main()
{
	CvMat2 * H;
	H = cvCreateMat2(3, 3, CV_32FC1);
	CvMat2 * N;
	N = cvCreateMat2(3, 3, CV_32FC1);
	cvmSet2(H, 0,0,1);
	cvmSet2(H, 0, 1, 1); 
	cvmSet2(H, 0, 2, 1);
	cvmSet2(H, 1, 0, 2);
	cvmSet2(H, 1, 1, 2);
	cvmSet2(H, 1, 2, 2);
	cvmSet2(H, 2, 0, 3);
	cvmSet2(H, 2, 1, 3);
	cvmSet2(H, 2, 2, 3);
	printf("%lf    %lf     %lf\n", cvmGet(H, 0, 0), cvmGet(H, 0, 1), cvmGet(H, 0, 2));
	printf("%lf    %lf     %lf\n", cvmGet(H, 1, 0), cvmGet(H, 1, 1), cvmGet(H, 1, 2));
	printf("%lf    %lf     %lf\n", cvmGet(H, 2, 0), cvmGet(H, 2, 1), cvmGet(H, 2, 2));
	printf("\n");
	printf("%lf    %lf     %lf\n", cvmGet(N, 0, 0), cvmGet(N, 0, 1), cvmGet(N, 0, 2));
	printf("%lf    %lf     %lf\n", cvmGet(N, 1, 0), cvmGet(N, 1, 1), cvmGet(N, 1, 2));
	printf("%lf    %lf     %lf\n", cvmGet(N, 2, 0), cvmGet(N, 2, 1), cvmGet(N, 2, 2));
	myMatcopy(H, N);
	printf("\n");
	printf("%lf    %lf     %lf\n", cvmGet(N, 0, 0), cvmGet(N, 0, 1), cvmGet(N, 0, 2));
	printf("%lf    %lf     %lf\n", cvmGet(N, 1, 0), cvmGet(N, 1, 1), cvmGet(N, 1, 2));
	printf("%lf    %lf     %lf\n", cvmGet(N, 2, 0), cvmGet(N, 2, 1), cvmGet(N, 2, 2));
	getchar();*/
/*	IplImage2 *i,*m;
	i = cvCreateImage2(cvSize2(640, 512), IPL_DEPTH_8U, 1);
	m = cvCreateImage2(cvSize2(640, 512), IPL_DEPTH_8U, 1);
	int j = 0;
	int k = 0;
	for (j = 0; j<2; j++)
		for (k = 0; k<2; k++)
			*(i->imageData+k*i->widthStep + j)='1';
	for (j = 0; j<2; j++)
		for (k = 0; k<2; k++)
			printf("%c\n", *(i->imageData + k*i->widthStep + j));
	printf("********************\n");
	for(j=0;j<2;j++)
		for(k=0;k<2;k++)
			printf("%c\n", *(m->imageData+k*m->widthStep+j));
	printf("********************\n");
	m = cvCloneImage(i);
	for (j = 0; j<2; j++)
		for (k = 0; k<2; k++)
			printf("%c\n", *(m->imageData + k*m->widthStep + j));
	getchar();
	}*/

	// Sharpen.cpp : Defines the entry point for the console application.
	/*-----CODE FOR FUN---------------
	-------CREATED BY Dream_Whui------
	-------2015-3-9-------------------------*/
	//锐化处理,采用两种方法，1.手动编写处理相邻像素键的关系 2.核函数


void box_filter(CvMat* img, CvMat* result) {
	//init part
	double s;
	int width = img->cols, height = img->rows;
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
			uchar pixel = (uchar)cvmGet(img, y, x);
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
			uchar pixel = (uchar)cvmGet(img, y, x);//img[y *width + x];    
			uchar pixel2 = (uchar)cvmGet(img, y + m_h, x);//img[(y+mheight) *width + x];  
			buff[x] = buff[x] - pixel + pixel2;
		}
	}
	//cout << buff << endl;
	//�������õ�ÿ����ĺͣ���������result
	for (y = 0; y<height - 5; y++) {
		for (x = 0; x<width; x++) {
			if (y>m_h / 2 && y<height - m_h / 2 && x>m_w / 2 && x<width - m_w / 2) {
				s = sum[(y - m_h / 2) *boxwidth + (x - m_h / 2)] / (m_h*m_w);
				cvmSet(result, y, x, s);
			}
			else {
				s = -1;
				cvmSet(result, y, x, s);
			}//end else
		}//end the first for
	}//end the second for
}
void box_filter(IplImage2* img, IplImage2* result) {
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
				result->imageData[y *width + x] = s;
			}
			else {
				s = -1;
				result->imageData[y *width + x] = s;
			}//end else
		}//end the first for
	}//end the second for
}
void Sharpen(IplImage *image, IplImage *result)
{
	//result.create(image.size(), image.type());
	int step = image->widthStep;
	int step1 = result->widthStep;
	//处理边界内部的像素点, 图像最外围的像素点应该额外处理
	for (int j = 1; j<image->height - 1; j++)
	{
		//前一行像素点
		uchar* previous = (uchar*) image->imageData+(j - 1)*step;
		//待处理的当前行
		uchar* current = (uchar*)image->imageData+j*step;
		//下一行
		uchar* next = (uchar*)image->imageData + (j + 1)*step;
		uchar * output = (uchar*)result->imageData+j*step1;
		for (int i = 1; i<(image->width - 1); i++)
		{
	//输出图像的遍历指针与当前行的指针同步递增, 以每行的每一个像素点的每一个通道值为一个递增量, 因为要考虑到图像的通道数
			int out1 = 5 * current[i * 3] - current[(i - 1) * 3] - current[(i + 1) * 3] - previous[i * 3] - next[i * 3];
			int out2= 5 * current[i * 3 + 1] - current[(i - 1) * 3 + 1] - current[(i + 1) * 3 + 1] - previous[i * 3 + 1] - next[i * 3 + 1];
			int out3 = 5 * current[i * 3 + 2] - current[(i - 1) * 3 + 2] - current[(i + 1) * 3 + 2] - previous[i * 3 + 2] - next[i * 3 + 2];
			if (out1 < 0)
				out1 = 0;
			else if (out1 > 255)
				out1=255;//防止数据溢出
			output[i * 3] = (uchar) out1;
			
			if (out2 < 0)
				out2 = 0;
			else if (out2 > 255)
				out2 = 255;
			output[i * 3+1] = (uchar)out2;

			if (out3 < 0)
				out3 = 0;
			else if (out3 > 255)
				out3= 255;
			output[i * 3] = (uchar)out3;
			
		}
	}
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
void Sharpen0(IplImage *image, IplImage *result)
{
	//result.create(image.size(), image.type());
	int step = image->widthStep;
	int step1 = result->widthStep;
	for (int j = 1; j<image->height - 1; j++)
	{
		uchar* previous = (uchar*)image->imageData + (j - 1)*step;
		uchar* current = (uchar*)image->imageData + j*step;
		uchar* next = (uchar*)image->imageData + (j + 1)*step;

		uchar * output = (uchar*)result->imageData + j*step1;
		for (int i = 1; i<(image->width - 1); i++)
		{
			output[i * 3] = saturate_cast<uchar>(5 * current[i * 3] - current[(i - 1) * 3] - current[(i + 1) * 3] - previous[i * 3] - next[i * 3]);
			output[i * 3 + 1] = saturate_cast<uchar>(5 * current[i * 3 + 1] - current[(i - 1) * 3 + 1] - current[(i + 1) * 3 + 1] - previous[i * 3 + 1] - next[i * 3 + 1]);
			output[i * 3 + 2] = saturate_cast<uchar>(5 * current[i * 3 + 2] - current[(i - 1) * 3 + 2] - current[(i + 1) * 3 + 2] - previous[i * 3 + 2] - next[i * 3 + 2]);
		}
	}
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
/*void SharpenKernel(const Mat &image, Mat &result)
{
	Mat kernel(3, 3, CV_32F, Scalar(0));
	kernel.at<float>(1, 1) = 5;
	kernel.at<float>(0, 1) = -1;
	kernel.at<float>(2, 1) = -1;
	kernel.at<float>(1, 0) = -1;
	kernel.at<float>(1, 2) = -1;
	filter2D(image, result, image.depth(), kernel);
	result.row(0).setTo(Scalar(0, 0, 0));
	result.row(result.rows - 1).setTo(Scalar(0, 0, 0));
	result.col(0).setTo(Scalar(0, 0, 0));
	result.col(result.cols - 1).setTo(Scalar(0, 0, 0));
}*/

int main()
{
	IplImage* image = cvLoadImage("111.jpg",1);
	
	IplImage *result=cvCreateImage(cvSize(image->width, image->height), IPL_DEPTH_8U, 1);
	IplImage *img = cvCreateImage(cvSize(image->width, image->height), IPL_DEPTH_8U, 1);
	cvCvtColor(image, img, CV_RGB2GRAY);
	namedWindow("image");
	cvShowImage("image", img);
	
	//滤波
	CvMat *mat = cvCreateMat(img->height, img->width, CV_64FC1);    //注意height和width的顺序
	cvConvert(img, mat);
	box_filter(mat, mat);
	IplImage *img1 = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_8U, 1);
	cvConvert(mat, img1);
	namedWindow("image_222");
	cvShowImage("image_222", img1);

	//锐化
	Sharpen(img, result);
	namedWindow("image_process");
	cvShowImage("image_process", result);


	waitKey();
	return 0;
}