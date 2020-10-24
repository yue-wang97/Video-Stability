//
// file : correspEstimation.h
//-----------------------------
// this file contains functions for correspondence
// estimation
//
#define WINDOWSIZE_NCC 16
#define SEARCH_RANGE   128
#define MAX_POINT_SIZE 500
typedef struct{
	int len;
	int rangeImagePositionI[MAX_POINT_SIZE];
	int rangeImagePositionJ[MAX_POINT_SIZE];
	int domainImagePositionI[MAX_POINT_SIZE];
	int domainImagePositionJ[MAX_POINT_SIZE];
}CorspMap;
void CorrespEstimation(IplImage2 *domainImage, IplImage2 *rangeImage,
					   IplImage2 *domainCornerMap, IplImage2 *rangeCornerMap,
					   CorspMap *corspMap);
int FindCorrespodence(IplImage2 *rangeImage, IplImage2 *domainBlock,
					   IplImage2 *rangeCornerMap,
					   int cornerPosI, int cornerPosJ,
					   int *pCorspPosI, int *pCorspPosJ,
					   int searchRange);
float NCC(IplImage2 *block1, IplImage2 *block2);
void UpdateCorrespMap(CorspMap *corspMap, int domainPosI, int domainPosJ,
					  int rangePosI, int rangePosJ);
void InitializeCorspMap(CorspMap *corspMap);
void CopyCorspMap(CorspMap *dst, CorspMap *src);