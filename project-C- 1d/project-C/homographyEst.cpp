//
// file : homographyEst.cpp
//------------------------------------------------
// this file contains functions for homography
// estimation
//
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <opencv2/opencv.hpp>  
#include <opencv2/highgui/highgui.hpp>  
#include <iostream>  
#include "cv.h"
#include "cxcore.h"
#include "highgui.h"
//#include"defs.h"
#include "correspEstimation.h"
#include "utility.h"

#include "svd.h"
#include "H.h"
#include "homographyEst.h"


using namespace std;
using namespace cv;


void HomograhyEstimation(CorspMap *inlierMap, float *H)
{
	int i;
	int numOfCorresp = inlierMap->len;
//	float (*domainPositions)[2]=malloc(sizeof(int) *numOfCorresp* 2);
//	float (*rangePositions)[2] = malloc(sizeof(int) *numOfCorresp * 2);
	float(*domainPositions)[2] = new float[numOfCorresp][2];
	float(*rangePositions)[2] = new float[numOfCorresp][2];
	for(i = 0; i < numOfCorresp; i++){
		// I -> y, J -> x
		domainPositions[i][1] = inlierMap->domainImagePositionI[i];
		domainPositions[i][0] = inlierMap->domainImagePositionJ[i];
		rangePositions[i][1] = inlierMap->rangeImagePositionI[i];
		rangePositions[i][0] = inlierMap->rangeImagePositionJ[i];
	}
	// compute the homography
	ComputeHomography(domainPositions, rangePositions, numOfCorresp, H);
	free(domainPositions);
	free(rangePositions);
}
void RansacHomograhyEstimation(CorspMap *corspMap, CorspMap *inlierMap,
							   float *H)
{
	int numOfCorresp = 4;
//	float (*domainPositions)[2] = malloc(sizeof(int) *numOfCorresp * 2);
//	float (*rangePositions)[2] = malloc(sizeof(int) *numOfCorresp * 2);
	float(*domainPositions)[2] = new float[numOfCorresp][2];
	float(*rangePositions)[2] = new float[numOfCorresp][2];
	int i, sampleCount = 0;
	int pos;
	float p, e, outlierProb;
	int numOfInliers, totalNumOfPoints, maxNumOfInliers = 1;
	CorspMap tempInlierMap;
	float Htmp[9];
	p = 0.99;
	e = 0.5;
	float N = 1000;
	while(N > sampleCount){
		// pick 4 corresponding points
		for(i = 0; i < numOfCorresp; i++){
			pos = rand() % corspMap->len; // select random positions
			// I -> y, J -> x
			domainPositions[i][1] = corspMap->domainImagePositionI[pos];
			domainPositions[i][0] = corspMap->domainImagePositionJ[pos];
			rangePositions[i][1] = corspMap->rangeImagePositionI[pos];
			rangePositions[i][0] = corspMap->rangeImagePositionJ[pos];
		}
		// check whether samples are good or not.
		// if the selected samples are good, then do homography estimation
		// else reselect samples.
		if(IsGoodSample(domainPositions, numOfCorresp) &&
			IsGoodSample(rangePositions, numOfCorresp)){
			// compute the homography
			ComputeHomography(domainPositions, rangePositions, numOfCorresp, Htmp);
			// calculate the distance for each correspondences
			// compute the number of inliers
			InitializeCorspMap(&tempInlierMap);
			CalculateDistance(Htmp, corspMap, &tempInlierMap);
			// choose H with the largest number of inliears
			numOfInliers = tempInlierMap.len;
			if(numOfInliers >= maxNumOfInliers){
				maxNumOfInliers = numOfInliers;
				CopyCorspMap(inlierMap, &tempInlierMap);
			//	myMatcopy(Htmp, H);
				copyH(Htmp, H);
			}
			// adaptive algorithm for determining the number of RANSAC samples
			// textbook algorithm 4.6
			totalNumOfPoints = corspMap->len;
			outlierProb = 1 - ((float)maxNumOfInliers / (float)totalNumOfPoints);
			e = (e < outlierProb ? e : outlierProb);
			N = log(1 - p) / log(1 - pow((1 - e), numOfCorresp));
			sampleCount += 1;
		}
	}
	free(domainPositions);
	free(rangePositions);
}
//
// function : ComputeHomography
// usage : ComputeHomography(domainPositions, rangePositions,
// numOfCorresp, H);
// -----------------------------------------------------------
// This function calculate the homography, H, using the set of
// given pairs of correspondences.
// Before computing the homography, data normalization will be
// performed. Then, it solve Ah = 0 using SVD to get H.
//
void ComputeHomography(float domainPositions[][2], float rangePositions[][2],
					   int numOfCorresp, float *H)
{
	int column = 9, row;
	float x1, y1, w1, x2, y2, w2;
	int i, ii, jj;
	float h[9];
	if(numOfCorresp == 4){
		row = 3; // eq 4.1 : make a thin matrix to solve SVD in opencv
	}else if(numOfCorresp > 4){
		row = 2; //eq 4.3
	}else{
		printf("Need more correspondence points! for computing H.\n");
		exit(0);
	}
//	float a[row * column];
//float *meanB1 = new float[channels];
	float *a= (float*)malloc(sizeof(float)*(row * column));
	float *aa = (float*)malloc(sizeof(float)*numOfCorresp * row*column);
	float *v = (float*)malloc(sizeof(float)*column*column);
	float *u = (float*)malloc(sizeof(float)*numOfCorresp * row*numOfCorresp * row);
//	CvMat2 *A;
//	CvMat2 *T1, *T2;
//	CvMat2 *D, *V, *U;
//	CvMat2 *invT2, *temp;
	// normalization
//	T1 = cvCreateMat2(3, 3, CV_32FC1);
//	T2 = cvCreateMat2(3, 3, CV_32FC1);
	float T1[9], T2[9],invT2[9],temp[9];
	DataNormalization(numOfCorresp, domainPositions, T1);
	DataNormalization(numOfCorresp, rangePositions, T2);
	// set A
//	A = cvCreateMat2(numOfCorresp * row, column, CV_32FC1);
	for(i = 0; i < numOfCorresp; i++){
		x1 = domainPositions[i][0];
		y1 = domainPositions[i][1];
		w1 = 1;
		x2 = rangePositions[i][0];
		y2 = rangePositions[i][1];
		w2 = 1;
		// set Ai
		// [0, 0, 0, -w2*x1, -w2*y1, -w2*w1, y2*x1, y2*y1, y2*w1]
		// [w2*x1, w2*y1, w2*w1, 0, 0, 0, -x2*x1, -x2*y1, -x2*w1]
		a[0] = 0; a[1] = 0; a[2] = 0;
		a[3] = -w2*x1; a[4] = -w2*y1; a[5] = -w2*w1;
		a[6] = y2*x1; a[7] = y2*y1; a[8] = y2*w1;
		a[9] = w2*x1; a[10] = w2*y1; a[11] = w2*w1;
		a[12] = 0; a[13] = 0; a[14] = 0;
		a[15] = -x2*x1; a[16] = -x2*y1; a[17] = -x2*w1;
		if(row == 3){ // eq 4.1 : make a thin matrix to solve SVD in opencv
			a[18] = -y2*x1; a[19] = -y2*y1; a[20] = -y2*w1;
			a[21] = x2*x1; a[22] = x2*y1; a[23] = x2*w1;
			a[24] = 0; a[25] = 0; a[26] = 0;
		}
		// assemble Ai into a matrix A
		for(ii = 0; ii < row; ii++){
			for(jj = 0; jj < column; jj++){
				SetH(aa, (ii + i*row)*column, jj, a[ii*column + jj]);
			}
		}
	}
	// calculate H
//	Htmp = cvCreateMat2(3, 3, CV_32FC1);
//	D = cvCreateMat2(numOfCorresp*row, column, CV_32FC1);
//	U = cvCreateMat2(numOfCorresp*row, numOfCorresp*row, CV_32FC1);
//	V = cvCreateMat2(column, column, CV_32FC1);
	

	//	cvSVD(A, D, U, V, CV_SVD_U_T|CV_SVD_V_T);// A = U^T D V : opencv setting
	// take last column of V
	
//	CvMat2Array(A, aa, numOfCorresp * row, column);
//	CvMat2Array(U, u, numOfCorresp*row, numOfCorresp*row);
//	CvMat2Array(V, v, column, column);
	float eps = 0.000001;
	int i2;
	i2 = dluav(aa, numOfCorresp * row, column, u, v, eps, (max(numOfCorresp * row, column) + 1));
//	Array2CvMat(v, V, column, column);

//	for(i = 0; i < column; i++){
//		h[i] = -cvmGet(V, column-1, i);
//	}
	for (i = 0; i < column; i++) {
		h[i] = -GetH(v, (column - 1)*3, i);
	}
//	Array2CvMat(h, Htmp, 3, 3);
	// denormalization : H = invT2 * Htmp * T1 <- Htmp = T2 * H * invT1
//	invT2 = cvCreateMat2(3, 3, CV_32FC1);
//	temp = cvCreateMat2(3, 3, CV_32FC1);
	myInvert(T2, invT2);
	arymul1(invT2, h, temp);
	arymul1(temp, T1, H);
	free(a);
	free(aa);
	free(v);
	free(u);

	// release matrices
	/*cvReleaseMat(&T1); cvReleaseMat(&T2);
	cvReleaseMat(&A); cvReleaseMat(&Htmp);
	cvReleaseMat(&D); cvReleaseMat(&U); cvReleaseMat(&V);
	cvReleaseMat(&T1); cvReleaseMat(&temp);*/
}
//
// function : DataNormalization
// usage : DataNormalization(numOfx, x, T);
// ------------------------------------------------------
// This function normalizes x and returns the similarity
// transform, T.
// The centroid of x will be transformed into (0,0).
// The average distance of normalized x will be sqrt(2).
//
void DataNormalization(int numOfPositions, float positions[][2], float *t)
{
	int i;
	float sumI = 0, sumJ = 0, meanI = 0, meanJ = 0;
	float squareDist = 0, sumDist = 0, meanDist = 0;
	float scale = 0;
	float x, y, xx, yy, ww;
	// calculate the centroid
	for(i = 0; i < numOfPositions; i++){
		sumI += positions[i][0];
		sumJ += positions[i][1];
	}
	meanI = sumI / numOfPositions;
	meanJ = sumJ / numOfPositions;
	// calculate the mean distance
	for(i = 0; i < numOfPositions; i++){
		squareDist = pow(positions[i][0] - meanI, 2)
			+ pow(positions[i][1] - meanJ, 2);
		sumDist += pow(squareDist, 0.5);
	}
	meanDist = sumDist / numOfPositions;
	// set the similarity transform
	scale = pow(1, 0.5) / meanDist;
	t[0] = scale;
	t[1] = 0;
	t[2] = -scale * meanI;
	t[3] = 0;
	t[4] = scale;
	t[5] = -scale * meanJ;
	t[6] = 0;
	t[7] = 0;
	t[8] = 1;	
	// data normalization
	for(i = 0; i < numOfPositions; i++){
		x = positions[i][0];
		y = positions[i][1];
		xx = t[0] * x + t[1] * y + t[2];
		yy = t[3] * x + t[4] * y + t[5];
		ww = t[6] * x + t[7] * y + t[8];
		xx = xx / ww;
		yy = yy / ww;
		positions[i][0] = xx;
		positions[i][1] = yy;
	}
}
//
// function : CalculateDistance
// usage : CalculateDistance(H, corspMap, inlierMap);
// ---------------------------------------------------
// This function calculates distance of data using
// symmetric transfer error. Then, compute inliers
// that consist with H.
//
void CalculateDistance(float *H, CorspMap *corspMap, CorspMap *inlierMap)
{
	int i;
	int x1, y1, x2, y2;
	float x1Trans, y1Trans, w1Trans, x2Trans, y2Trans, w2Trans;
	float dist2x1AndInvHx2, dist2x2AndHx1, dist2Trans;
	float tSquare = T_SQUARE;
	//CvMat2 *invH = cvCreateMat2(3, 3, CV_32FC1);
	float invH[9];
	myInvert(H, invH);
	// use d^2_transfer as distance measure
	for(i = 0; i < corspMap->len; i++){
		// I -> y, J -> x
		x1 = corspMap->domainImagePositionJ[i];
		y1 = corspMap->domainImagePositionI[i];
		x2 = corspMap->rangeImagePositionJ[i];
		y2 = corspMap->rangeImagePositionI[i];
		// calculate x_trans = H * x
		x2Trans = GetH(H, 0, 0) * x1 + GetH(H, 0, 1) * y1
			+ GetH(H, 0, 2);
		y2Trans = GetH(H, 1, 0) * x1 + GetH(H, 1, 1) * y1
			+ GetH(H, 1, 2);
		w2Trans = GetH(H, 2, 0) * x1 + GetH(H, 2, 1) * y1
			+ GetH(H, 2, 2);
		x2Trans = x2Trans / w2Trans;
		y2Trans = y2Trans / w2Trans;
		// calculate x'_trans = H^(-1) * x'
		x1Trans = GetH(invH, 0, 0) * x2 + GetH(invH, 0, 1) * y2+ GetH(invH, 0, 2);
		y1Trans = GetH(invH, 1, 0) * x2 + GetH(invH, 1, 1) * y2
			+ GetH(invH, 1, 2);
		w1Trans = GetH(invH, 2, 0) * x2 + GetH(invH, 2, 1) * y2
			+ GetH(invH, 2, 2);
		x1Trans = x1Trans / w1Trans;
		y1Trans = y1Trans / w1Trans;
		// calculate the square distance (symmetric transfer error)
		dist2x1AndInvHx2 = pow(x1 - x1Trans, 2) + pow(y1 - y1Trans, 2);
		dist2x2AndHx1 = pow(x2 - x2Trans, 2) + pow(y2 - y2Trans, 2);
		dist2Trans = dist2x1AndInvHx2 + dist2x2AndHx1;
		if(dist2Trans < tSquare){
			UpdateCorrespMap(inlierMap, y1, x1, y2, x2);
		}
	}
	// release matrices
	//cvReleaseMat(&invH);
}
//
// function : IsGoodSample
// usage : r = IsGoodSample(points, numOfPoints)
// -------------------------------------------------
// This function checks colinearity of all given points.
//
int IsGoodSample(float points[][2], int numOfPoints)
{
	int r;
	int i, j, k;
//	CvMat2 *A = cvCreateMat2(3, 1, CV_32FC1);
//	CvMat2 *B = cvCreateMat2(3, 1, CV_32FC1);
//	CvMat2 *C = cvCreateMat2(3, 1, CV_32FC1);
	float A[3], B[3], C[3];
	i = 0;
	j = i + 1;
	k = j + 1;
	r = 0;
	// check colinearity recursively
	while(1){
		// set point vectors
		SetH(A, 0, 0, points[i][0]);
		SetH(A, 1, 0, points[i][1]);
		SetH(A, 2, 0, 1);
		SetH(B, 0, 0, points[j][0]);
		SetH(B, 1, 0, points[j][1]);
		SetH(B, 2, 0, 1);
		SetH(C, 0, 0, points[k][0]);
		SetH(C, 1, 0, points[k][1]);
		SetH(C, 2, 0, 1);
		// check linearity
		r = IsColinear(A, B, C) || r;
		// update point index
		if(k < numOfPoints - 1){
			k += 1;
		}else{
			if(j < numOfPoints - 2){
				j += 1;
				k = j + 1;
			}else{
				if(i < numOfPoints - 3){
					i += 1;
					j = i + 1;
					k = j + 1;
				}else{
					break;
				}
			}
		}
	}
	return(!r);
}
//
// function : IsColinear
// usage : r = IsColinear(A, B, C);
// --------------------------------------
// This function checks the colinearity of
// the given 3 points A, B, and C.
// If these are colinear, it returns false.
//
int IsColinear(float* A, float* B, float* C)
{
	float x1, x2, x3, y1, y2, y3;
	int type, step;
	type = CV_MAT_TYPE(CV_32FC1);
	step = CV_ELEM_SIZE(type);
	x1 = *A;
	y1 = *(float *)(A + step);
	x2 = *B;
	y2 = *(float *)(B + step);
	x3 = *C;
	y3 = *(float *)(C + step);
	if (fabs((x1 - x2) * y3 - (y2 - y1) * x3 + y1*x2 - x1*y2) < EPS)
		return 0;
	else
		return 1;
}

