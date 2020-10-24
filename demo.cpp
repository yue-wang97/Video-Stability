#include <cv.h>    
#include <highgui.h>    
#include <stdio.h>    

using namespace std;

int horPrj[240] = {0};
int verPrj[320] = {0};

void init(uchar *pGrayData)
{
	int m,n,i,j;
	for (m = 0, i = 120; i < 360 ; i++)
	{
		horPrj[m] = 0;
		for ( j =  160 ; j < 480; j++)
		{
			horPrj[m] += pGrayData[ i*640 + j ] ; 					
		}
		horPrj[m] /= 320;
		m++;
	}

	for (n = 0, j = 160; j < 480 ; j++)
	{
		verPrj[n] = 0;
		for ( i =  120 ; i < 360; i++)
		{
			verPrj[n] += pGrayData[ i*640 + j ] ;
		}
		verPrj[n] /= 240;
		n++;
	}
}

int main(int argc, char* argv[])
{
	IplImage * pFrame;
	CvCapture * pCapture = cvCaptureFromAVI( "D:\\AAA.avi" );
	//CvCapture * pCapture = cvCaptureFromCAM(0);
	//	CvVideoWriter * writer = cvCreateVideoWriter("AAA.avi",-1,15,cvSize(640,480));

	//pFrame = cvQueryFrame( pCapture );

	IplImage * pGrayImg = cvCreateImage( cvSize(640,480), 8, 1 );
	//pGrayImg->origin = pFrame->origin;
	uchar * pGrayData = (uchar *)pGrayImg->imageData;
	int nStep = pGrayImg->widthStep/(sizeof(uchar));

	IplImage * pStableImg = cvCreateImage( cvSize(640,480), 8, 1 );
	//pStableImg->origin = pFrame->origin;
	uchar * pStableData = (uchar *)pStableImg->imageData;
	int nCStep = pStableImg->widthStep/(sizeof(uchar));
	cvZero(pStableImg);

	//IplImage * pHisHorImg = cvCreateImage( cvSize(256,240), 8, 3 );
	//pHisHorImg->origin = 1;
	//cvZero(pHisHorImg);

	//IplImage * pHisVerImg = cvCreateImage( cvSize(320,256), 8, 3 );
	//pHisVerImg->origin = 1;
	//cvZero(pHisVerImg);

	//IplImage * pCHisHorImg = cvCreateImage( cvSize(256,480), 8, 3 );
	//pCHisHorImg->origin = 1;
	//cvZero(pCHisHorImg);
	//
	//IplImage * pCHisVerImg = cvCreateImage( cvSize(640,256), 8, 3 );
	//pCHisVerImg->origin = 1;
	//cvZero(pCHisVerImg);


	cvNamedWindow("Origin");
	cvNamedWindow("Stabilization");
	//	cvNamedWindow("HisHor");
	//	cvNamedWindow("HisVer");
	//	cvNamedWindow("CHisHor");
	//	cvNamedWindow("CHisVer");

	//int nWidth = cvGetSize(pFrame).width;
	//int nHeight = cvGetSize(pFrame).height;
	int i,j;
	long numFrame = 0;

	//投影数组


	int CurHorPrj[480] = {0};
	int CurVerPrj[640] = {0};

	int m , n ;

	int DevX = 0, DevY = 0;

	//pFrame = cvQueryFrame( pCapture );
	//cvWaitKey(1000);
	//pFrame = cvQueryFrame( pCapture );

	//cvCvtColor( pFrame, pGrayImg, CV_BGR2GRAY );
	//cvSmooth( pGrayImg, pGrayImg, CV_MEDIAN) ;

	//参考帧投影
	// 	for ( i = 120; i < 360 ; i++)
	// 	{
	// 		temp = i*nStep + 160;
	// 		horPrj[m] = 0;
	// 		for ( j =  160 ; j < 480; j++)
	// 		{
	// 			horPrj[m] += pGrayData[ temp++ ] ; 					
	// 		}
	// 		horPrj[m] /= 320;
	// 		m++;
	// 	}
	// 
	// 	for ( j = 160; j < 480 ; j++)
	// 	{
	// 		temp = 120*nStep + j;
	// 		verPrj[n] = 0;
	// 		for ( i =  120 ; i < 360; i++)
	// 		{
	// 			verPrj[n] += pGrayData[ temp ] ;
	// 			temp += nStep;
	// 		}
	// 		verPrj[n] /= 240;
	// 		n++;
	// 	}
	// 
	// 	for ( i = 0; i < 240 ; i++)
	// 	{
	// 		cvLine(pHisHorImg, cvPoint(i,0), cvPoint(i,horPrj[i]), CV_RGB(255,100,100), 1, 8);
	// 	}
	// 	for ( i = 0; i < 320 ; i++)
	// 	{
	// 		cvLine(pHisVerImg, cvPoint(i,0), cvPoint(i,verPrj[i]), CV_RGB(100,255,100), 1, 8);
	// 	}




	//for ( i = 0; i < 240 ; i++)
	//{
	//	cvLine(pHisHorImg, cvPoint(0,i), cvPoint(horPrj[i],i), CV_RGB(255,100,100), 1, 8);
	//}
	//for ( i = 0; i < 320 ; i++)
	//{
	//	cvLine(pHisVerImg, cvPoint(i,0), cvPoint(i,verPrj[i]), CV_RGB(100,255,100), 1, 8);
	//}

	// 	cvShowImage("HisHor",pHisHorImg);
	//	cvShowImage("HisVer",pHisVerImg);

	while( pFrame = cvQueryFrame( pCapture ) )
	{
		cvCvtColor( pFrame, pGrayImg, CV_BGR2GRAY );
		if(numFrame == 0)
		{
			init(pGrayData);
		}
		//		cvWriteFrame(writer,pFrame);
		numFrame++;
		
		//cvSmooth( pGrayImg, pGrayImg, CV_MEDIAN) ;

		m = 0;
		n = 0;

		//当前帧投影
		for ( i = 0; i < 480 ; i++)
		{
			CurHorPrj[m] = 0;
			for ( j = 160 ; j < 480; j++)
			{
				CurHorPrj[m] += pGrayData[ i*nStep + j ] ;

			}
			CurHorPrj[m] /= 320;
			m++;
		}

		for ( j = 0; j < 640 ; j++)
		{
			CurVerPrj[n] = 0;
			for ( i = 120 ; i < 360; i++)
			{
				CurVerPrj[n] += pGrayData[ i*nStep + j ] ;
			}
			CurVerPrj[n] /= 240;
			n++;
		}

		//for ( i = 0; i < 480 ; i++)
		//{
		//	cvLine(pCHisHorImg, cvPoint(0,i), cvPoint(CurHorPrj[i],i), CV_RGB(255,100,100), 1, 8);
		//}
		//cvLine(pCHisHorImg, cvPoint(0,120), cvPoint(255,120), CV_RGB(255,100,100), 1, 8);

		//for ( i = 0; i < 640 ; i++)
		//{
		//	cvLine(pCHisVerImg, cvPoint(i,0), cvPoint(i,CurVerPrj[i]), CV_RGB(100,255,100), 1, 8);
		//}
		//cvLine(pCHisVerImg, cvPoint(160,0), cvPoint(160,255), CV_RGB(100,255,100), 1, 8);


		//相关运算
		long MinY = 1000000000;
		long SumY = 0;
		long MinX = 1000000000;
		long SumX = 0;

		for ( i = 20; i < 220; i++)
		{
			SumY = 0;
			for ( j = 0; j < 240; j++)
			{
				SumY += abs(horPrj[j] - CurHorPrj[j+i]);
			}
			if (SumY < MinY)
			{
				MinY = SumY;
				DevY = i - 120;
			}			  
		}

		for ( i = 60; i < 260; i++)
		{
			SumX = 0;
			for ( j = 0; j < 320; j++)
			{
				SumX += abs(verPrj[j] - CurVerPrj[j+i]);
			}
			if (SumX < MinX)
			{
				MinX = SumX;
				DevX = i - 160;
			}			  
		}

		cout<<"Y "<<DevY<<"  X "<<DevX<<endl;

		//运动补偿

		for ( i = 100; i < 380 ; i++)
		{
			for ( j =  100 ; j < 540; j++)
			{
				pStableData[i*nCStep + j] = pGrayData[(i+DevY)*nCStep + (j+DevX)] ; 
				//pStableData[i*nCStep + 3*j + 1] = pFrame->imageData[(i+DevY)*nCStep + 3*(j+DevX) + 1] ; 
				//pStableData[i*nCStep + 3*j + 2] = pFrame->imageData[(i+DevY)*nCStep + 3*(j+DevX) + 2] ; 
			}

		}



		cvShowImage("Origin",pFrame);

		cvShowImage("Stabilization",pStableImg);
		//		cvResetImageROI(pFrame);

		//		cvShowImage("CHisHor",pCHisHorImg);
		//		cvShowImage("CHisVer",pCHisVerImg);
		//cvZero(pCHisHorImg);
		//cvZero(pCHisVerImg);

		char c = cvWaitKey(50);
		if( c == 27 )
			break;
	}

	//	cvReleaseVideoWriter(&writer);
	cvReleaseCapture( &pCapture );
	cvDestroyAllWindows();

	return 0;
}
