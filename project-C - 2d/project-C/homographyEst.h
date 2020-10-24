//
// file : homographyEst.h
//------------------------------------------------
// this file contains functions for homography
// estimation
//
#define T_SQUARE 36 // t = sqrt(6)* sigma and set sigma = sqrt(6)
void HomograhyEstimation(CorspMap *inlierMap, float H[3][3]);
void RansacHomograhyEstimation(CorspMap *corspMap, CorspMap *inlierMap,
							   float H[3][3]);
void ComputeHomography(float domainPosiitons[][2], float rangePosiitons[][2],
					   int numOfCorresp, float H[3][3]);
void DataNormalization(int numOfPositions, float positions[][2], float T[3][3]);
void CalculateDistance(float H[3][3], CorspMap *corspMap, CorspMap *inlierMap);
int IsGoodSample(float points[][2], int numOfPoints);
int IsColinear(float A[3][1], float B[3][1], float C[3][1]);
