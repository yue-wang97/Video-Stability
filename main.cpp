#include<iostream>
#include <cv.h>    
#include <highgui.h>    
#include <stdio.h> 
using namespace std;
using namespace cv;

int Width=640;
int Height=480;
int W = Width / 4;
int H = Height / 4;
int horPrj[240] = { 0 };//Height-2*H
int verPrj[320] = { 0 };//Width-2*W

void init(uchar *pGrayData)
{
	int m=0, n=0, i, j; 
	for ( i = H; i <(Height-H); i++)
	{
		horPrj[m] = 0;
		for (j = W; j < Width-W; j++)
		{
			horPrj[m] += pGrayData[i * Width + j];
		}
		horPrj[m] /= (Width-2*W);
		m++;
	}
	for ( j = W; j < Width-W; j++)
	{
		verPrj[n] = 0;
		for (i = H; i < Height-H; i++)
		{
			verPrj[n] += pGrayData[i * Width + j];
		}
		verPrj[n] /= (Height-2*H);
		n++;
	}
}
int main(int argc, char* argv[])
{
	int i, j;
	long numFrame = 0;
	//投影数组
	int CurHorPrj[480] = { 0 };
	int CurVerPrj[640] = { 0 };
	int m, n;
	int DevX = 0, DevY = 0;
	IplImage * pFrame;
	CvCapture * pCapture = cvCaptureFromAVI("AAA.avi");
	CvVideoWriter * writer = cvCreateVideoWriter("out.avi", CV_FOURCC('M', 'J', 'P', 'G'), 20, cvSize(Width,Height),0);
	IplImage * pGrayImg = cvCreateImage(cvSize(Width, Height), 8, 1);
	uchar * pGrayData = (uchar *)pGrayImg->imageData;
	int nStep = pGrayImg->widthStep / (sizeof(uchar));//Width
	IplImage * pStableImg = cvCreateImage(cvSize(Width, Height), 8, 1);
	uchar * pStableData = (uchar *)pStableImg->imageData;
	int nCStep = pStableImg->widthStep / (sizeof(uchar));//Width
	cvZero(pStableImg);
	cvNamedWindow("Origin");
	cvNamedWindow("Stabilization");
	while (pFrame = cvQueryFrame(pCapture))
	{
		cvCvtColor(pFrame, pGrayImg, CV_BGR2GRAY);
		if (numFrame == 0)
		{
			init(pGrayData);
		}
		numFrame++;
		m = 0;
		n = 0;
		//当前帧投影
		for (i = 0; i < Height; i++)
		{
			CurHorPrj[m] = 0;
			for (j = W; j < Width-W; j++)
			{
				CurHorPrj[m] += pGrayData[i*nStep + j];

			}
			CurHorPrj[m] /= (Width-2*W);
			m++;
		}
		for (j = 0; j < Width; j++)
		{
			CurVerPrj[n] = 0;
			for (i = H; i < Height-H; i++)
			{
				CurVerPrj[n] += pGrayData[i*nStep + j];
			}
			CurVerPrj[n] /= (Height - 2*H);
			n++;
		}
		//相关运算
		long MinY = 1000000000;
		long SumY = 0;
		long MinX = 1000000000;
		long SumX = 0;
		for (i = H-50; i < H+50; i++)
		{
			SumY = 0;
			for (j = 0; j < 2*H; j++)
			{
				SumY += abs(horPrj[j] - CurHorPrj[j + i]);
			}
			if (SumY < MinY)
			{
				MinY = SumY;
				DevY = i - H;
			}
		}
		for (i = W-50; i < W+50; i++)
		{
			SumX = 0;
			for (j = 0; j < 2*W; j++)
			{
				SumX += abs(verPrj[j] - CurVerPrj[j + i]);
			}
			if (SumX < MinX)
			{
				MinX = SumX;
				DevX = i - W;
			}
		}
		cout << "Y " << DevY << "  X " << DevX << endl;
		//运动补偿
		for (i = 50; i < Height-50; i++)
		{
			for (j = 50; j < Width-50; j++)
			{
				pStableData[i*nCStep + j] = pGrayData[(i + DevY)*nCStep + (j + DevX)];
			}
		}
		cvWriteFrame(writer, pStableImg);
		cvShowImage("Origin", pFrame);
		cvShowImage("Stabilization", pStableImg);
		char c = cvWaitKey(50);
		if (c == 27)
			break;
	}
	cvReleaseVideoWriter(&writer);
	cvReleaseCapture(&pCapture);
	cvDestroyAllWindows();
	return 0;
}
