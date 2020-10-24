//
// file : homographyEst.h
//------------------------------------------------
// this file contains functions for homography
// estimation
//
#define T_SQUARE 36 // t = sqrt(6)* sigma and set sigma = sqrt(6)
void HomograhyEstimation(CorspMap *inlierMap, CvMat2 *H);
void RansacHomograhyEstimation(CorspMap *corspMap, CorspMap *inlierMap,
							   CvMat2 *H);
void ComputeHomography(float domainPosiitons[][2], float rangePosiitons[][2],
					   int numOfCorresp, CvMat2 *H);
void DataNormalization(int numOfPositions, float positions[][2], CvMat2 *T);
void CalculateDistance(CvMat2 *H, CorspMap *corspMap, CorspMap *inlierMap);
int IsGoodSample(float points[][2], int numOfPoints);
int IsColinear(CvMat2 *A, CvMat2 *B, CvMat2 *C);
