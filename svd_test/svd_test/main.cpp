
#include <math.h>
#include <opencv2/opencv.hpp>  
#include <opencv2/highgui/highgui.hpp>  
#include <iostream>  
#include "cv.h"
#include "cxcore.h"
#include "highgui.h"
#include"filter.h"

using namespace std;
using namespace cv;

#include "svd.h"
#include "stdio.h"
void main()
{
	int i, j;
	double a[12] = { 1.0,1.0,-1.0,2.0,1.0,0.0,1.0,-1.0,0.0,-1.0,2.0,1.0 };
	double b[12] = { 1.0,1.0,-1.0,-1.0,2.0,1.0,0.0,2.0,1.0,-1.0,0.0,1.0 };
	double u[16], v[9], c[12], d[12];

	for (i = 0; i<16; i++)
		u[i] = 0;
	double eps = 0.000001;
	i = dluav(a, 4, 3, u, v, eps, 5);

	printf("\n");
	getchar();

}
void main()
{
	int i, j;
	/*double a[4][3]={{1.0,1.0,-1.0},{2.0,1.0,0.0},{1.0,-1.0,0.0},{-1.0,2.0,1.0}};
	double b[3][4]={{1.0,1.0,-1.0,-1.0},{2.0,1.0,0.0,2.0},{1.0,-1.0,0.0,1.0}};
	double u[4][4],v[3][3],c[4][3],d[3][4];*/
	double a[12] = { 1.0,1.0,-1.0,2.0,1.0,0.0,1.0,-1.0,0.0,-1.0,2.0,1.0 };
	double u[16], v[9], c[12], d[12];

	for (i = 0; i<16; i++)
		u[i] = 0;
	double eps = 0.000001;
	//i = dluav(a, 4, 3, u, v, eps, 5);
	CvMat *D, *V, *U, *A;

	A = cvCreateMat(4, 3, CV_32FC1);
	D = cvCreateMat(4, 3, CV_32FC1);
	U = cvCreateMat(4, 4, CV_32FC1);
	V = cvCreateMat(3, 3, CV_32FC1);
	Array2CvMat1(a, A, 4, 3);
	cvSVD(A, D, U, V, CV_SVD_U_T | CV_SVD_V_T);// A = U^T D V : opencv setting
	CvMatArray1(A, a, 4, 3);
	CvMatArray1(D, d, 4, 3);
	CvMatArray1(U, u, 4, 4);
	CvMatArray1(V, v, 3, 3);
	printf("\n");

	printf("\nMAT U Is:\n");
	for (i = 0; i <= 3; i++)
	{
		for (j = 0; j <= 3; j++)
			printf("%e ", u[i * 4 + j]);
		printf("\n");
	}
	printf("\n");

	printf("MAT D Is:\n");
	for (i = 0; i <= 3; i++)
	{
		for (j = 0; j <= 2; j++)
			printf("%e ", a[i * 3 + j]);
		printf("\n");
	}
	printf("\n");

	printf("MAT V IS:\n");
	for (i = 0; i <= 2; i++)
	{
		for (j = 0; j <= 2; j++)
			printf("%e ", v[i * 3 + j]);
		printf("\n");
	}
	printf("\n");

	printf("MAT A Is:\n");
	for (i = 0; i <= 3; i++)
	{
		for (j = 0; j <= 2; j++)
			printf("%e ", a[i * 3 + j]);
		printf("\n");
	}

	printf("\n\n");
	getchar();

}