//
// file : homographyEst.cpp
//------------------------------------------------
// this file contains functions for homography
// estimation
//
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include"defs.h"
#include "utility.h"
#include "correspEstimation.h"
#include "homographyEst.h"
#include"svd.h"


void HomograhyEstimation(CorspMap *inlierMap, float H[3][3])
{
	int i;
	int numOfCorresp = inlierMap->len;
//	float(*domainPositions)[2] = (float(*)[2])malloc(sizeof(float) *numOfCorresp * 2);
//	float(*rangePositions)[2] = (float(*)[2])malloc(sizeof(float) *numOfCorresp * 2);
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
							   float H[3][3])
{
	int numOfCorresp = 4;
	//int(*a)[2] = (int(*)[2])malloc(sizeof(int) * 3 * 2);
	float (*domainPositions)[2] = (float(*)[2])malloc(sizeof(float) *numOfCorresp * 2);
	float (*rangePositions)[2] = (float(*)[2])malloc(sizeof(float) *numOfCorresp * 2);
//	float(*domainPositions)[2] = new float[numOfCorresp][2];
//	float(*rangePositions)[2] = new float[numOfCorresp][2];
	int i, sampleCount = 0;
	int pos;
	float p, e, outlierProb;
	int numOfInliers, totalNumOfPoints, maxNumOfInliers = 1;
	CorspMap tempInlierMap;
	//CvMat2 *Htmp;
	//Htmp = cvCreateMat2(3, 3, CV_32FC1);
	float Htmp[3][3];
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
//				copy2(Htmp, H);
				int i,j;
				for (i = 0; i < 3; i++)
					for (j = 0; j < 3; i++) 
						H[i][j] = Htmp[i][j];		

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
					   int numOfCorresp, float H[3][3])
{
	int column = 9, row;
	float x1, y1, w1, x2, y2, w2;
	int i, ii, jj;
	
	if(numOfCorresp == 4){
		row = 3; // eq 4.1 : make a thin matrix to solve SVD in opencv
	}else if(numOfCorresp > 4){
		row = 2; //eq 4.3
	}else{
		printf("Need more correspondence points! for computing H.\n");
		exit(0);
	}
	//float a[row * column];
	//float *meanB1 = new float[channels];
	float *a= (float*)malloc(sizeof(float)*(row * column));
	/*CvMat2 *A, *Htmp;
	CvMat2 *T1, *T2;
	CvMat2 *V, *U;
	CvMat2 *invT2, *temp;*/
	// normalization
	float T1[3][3],T2[3][3];
	DataNormalization(numOfCorresp, domainPositions, T1);
	DataNormalization(numOfCorresp, rangePositions, T2);
	// set A
//	A = cvCreateMat2(numOfCorresp * row, column, CV_32FC1);
	float(*A)[9] = (float(*)[9])malloc(sizeof(float) * numOfCorresp * row * 9);
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
		a[3] =-w2*x1; a[4] = -w2*y1; a[5] = -w2*w1;
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
				A[ii + i*row][jj]= a[ii*column+jj];
			}
		}
	}
	// calculate H
	float Htmp[3][3];
//	U = cvCreateMat2(numOfCorresp*row, numOfCorresp*row, CV_32FC1);
/*	float **U;
	U = (float**)malloc(sizeof(float*)*numOfCorresp*row);
	int i;
	for (i = 0; i<numOfCorresp*row; i++)
		U[i] = (float*)malloc(sizeof(float)*numOfCorresp*row);*/
	float *u = (float*)malloc(sizeof(float)*numOfCorresp * row*numOfCorresp * row);
	float *v = (float*)malloc(sizeof(float)*column*column);
//	float V[9][9];
	//	cvSVD(A, D, U, V, CV_SVD_U_T|CV_SVD_V_T);// A = U^T D V : opencv setting
	// take last column of V
    float *aa = (float*)malloc(sizeof(float)*numOfCorresp * row*column);
/*	float *v = (float*)malloc(sizeof(float)*column*column);
	float *u = (float*)malloc(sizeof(float)*numOfCorresp * row*numOfCorresp * row);
	CvMat2Array(A, aa, numOfCorresp * row, column);
	CvMat2Array(U, u, numOfCorresp*row, numOfCorresp*row);
	CvMat2Array(V, v, column, column);*/
	float eps = 0.000001;
	int i2;
	int  j, k = 0;
	for (i = 0; i < numOfCorresp * row; i++)
		for (j = 0; j < column; j++)
			aa[k++] = A[i][j];
	i2 = dluav(aa, numOfCorresp * row, column,u, v, eps, (max(numOfCorresp * row, column) + 1));
	//Array2CvMat(v, V, column, column);

	/*for(i = 0; i < column; i++){
		h[i] = -V[column-1][i];
	}*/
	k = 0;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
		{ 
			Htmp[i][j] = -v[8*36+k];
			k++;
		}
	//Array2CvMat(h, Htmp, 3, 3);
	// denormalization : H = invT2 * Htmp * T1 <- Htmp = T2 * H * invT1
	float invT2[3][3], temp[3][3];
	myInvert(T2, invT2);
	arymul1(invT2, Htmp, temp);
	arymul1(temp, T1, H);
	free(a);
	free(A);
/*	for (i = 0; i<numOfCorresp*row; i++)
		free(U[i]);
	free(U);*/
	free(u);
	free(v);
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
void DataNormalization(int numOfPositions, float positions[][2], float T[3][3])
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
	float t[9] = {scale, 0, -scale * meanI,
		0, scale, -scale * meanJ,
		0, 0, 1};
	T[0][0] = scale;
	T[0][1] = 0;
	T[0][2] = -scale * meanI;
	T[1][0] = 0;
	T[1][1] = scale;
	T[1][2] = -scale * meanI;
	T[2][0] = 0;
	T[2][1] = 0;
	T[2][2] = 1;
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
void CalculateDistance(float H[3][3], CorspMap *corspMap, CorspMap *inlierMap)
{
	int i;
	int x1, y1, x2, y2;
	float x1Trans, y1Trans, w1Trans, x2Trans, y2Trans, w2Trans;
	float dist2x1AndInvHx2, dist2x2AndHx1, dist2Trans;
	float tSquare = T_SQUARE;
	float invH[3][3];
	myInvert(H, invH);
	// use d^2_transfer as distance measure
	for(i = 0; i < corspMap->len; i++){
		// I -> y, J -> x
		x1 = corspMap->domainImagePositionJ[i];
		y1 = corspMap->domainImagePositionI[i];
		x2 = corspMap->rangeImagePositionJ[i];
		y2 = corspMap->rangeImagePositionI[i];
		// calculate x_trans = H * x
		x2Trans = H[0][0] * x1 + H[0][1] * y1
			+ H[0][2];
		y2Trans = H[1][0] * x1 + H[1][1] * y1
			+ H[1][2];
		w2Trans = H[2][0] * x1 + H[2][1] * y1
			+ H[2][2];
		x2Trans = x2Trans / w2Trans;
		y2Trans = y2Trans / w2Trans;
		// calculate x'_trans = H^(-1) * x'
		x1Trans = invH[0][0] * x2 +invH[0][1] * y2+ invH[0][2];
		y1Trans = invH[1][0] * x2 + invH[1][1] * y2
			+ invH[1][2];
		w1Trans = invH[2][0] * x2 + invH[2][1] * y2
			+ invH[2][2];
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
	float A[3][1];
	float B[3][1];
	float C[3][1];
	i = 0;
	j = i + 1;
	k = j + 1;
	r = 0;
	// check colinearity recursively
	while(1){
		// set point vectors
		A[0][0]= points[i][0];
		A[1][0]= points[i][1];
		A[2][0]=1;
		B[0][0]= points[j][0];
		B[1][0]= points[j][1];
		B[2][0]= 1;
		C[0][0]= points[k][0];
		C[1][0]= points[k][1];
		C[2][0]= 1;
		// check linearity
		r = IsColinear(A, B, C) || r;
		// update point index
		if(k < numOfPoints - 1){
			k += 1;
		}
		else{
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
int IsColinear(float A[3][1], float B[3][1], float C[3][1])
{
	float x1, x2, x3, y1, y2, y3;
	x1 = A[0][0];
	y1 = A[1][0];
	x2 = B[0][0];
	y2 =B[1][0];
	x3 = C[0][0];
	y3 =C[1][0];
	if (fabs((x1 - x2) * y3 - (y2 - y1) * x3 + y1*x2 - x1*y2) < EPS)
		return 0;
	else
		return 1;
}

