//
// file : correspEstimation.h
//-----------------------------
// this file contains functions for correspondence
// estimation
//
#ifndef CORRESPESTIMATION_H_
#define CORRESPESTIMATION_H_

//#include "defs.h"
#define WINDOWSIZE_NCC 16
#define SEARCH_RANGE   128
#define MAX_POINT_SIZE 500//允许最多匹配500个点

typedef struct{
	int len;//匹配角点的个数
	int rangeImagePositionI[MAX_POINT_SIZE];//当前帧角点坐标y
	int rangeImagePositionJ[MAX_POINT_SIZE];//当前帧角点坐标x
	int domainImagePositionI[MAX_POINT_SIZE];//参考帧角点坐标y
	int domainImagePositionJ[MAX_POINT_SIZE];//参考帧角点坐标x
}CorspMap;//用来存储对应匹配的角点
void CorrespEstimation(IplImage2 *domainImage, IplImage2 *rangeImage,
					   IplImage2 *domainCornerMap, IplImage2 *rangeCornerMap,
					   CorspMap *corspMap);//对参考帧和当前帧的角点进行粗匹配
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
//void getCorspMap(CorspMap *corspMap, Point2 *point1_innext, Point2 *point2_temp, int len);

#endif