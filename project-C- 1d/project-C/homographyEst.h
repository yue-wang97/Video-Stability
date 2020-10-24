//
// file : homographyEst.h
//------------------------------------------------
// this file contains functions for homography
// estimation
//
#define T_SQUARE 36 // t = sqrt(6)* sigma and set sigma = sqrt(6)
void HomograhyEstimation(CorspMap *inlierMap, float* H);
void RansacHomograhyEstimation(CorspMap *corspMap, CorspMap *inlierMap,
							   float* H);
void ComputeHomography(float domainPosiitons[][2], float rangePosiitons[][2],
					   int numOfCorresp, float *H);
void DataNormalization(int numOfPositions, float positions[][2], float* t);
void CalculateDistance(float* H, CorspMap *corspMap, CorspMap *inlierMap);
int IsGoodSample(float points[][2], int numOfPoints);
int IsColinear(float* A, float* B, float* C);
