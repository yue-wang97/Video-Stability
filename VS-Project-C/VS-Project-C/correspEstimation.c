//
// file : correspEstimation.cpp
//-----------------------------
// this file contains functions for correspondence
// estimation
//
#pragma once
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//#include "point.h"
#include "defs.h"
#include "utility.h"
#include "correspEstimation.h"

////粗匹配函数///
///根据参考帧，当前帧及两帧检测出的角点，采用NCC（ 归一化互相关）的方法，通过计算模板图像和匹配图像
///的互相关值，来确定匹配的程度从而进行角点间的粗匹配
void CorrespEstimation(IplImage2 *domainImage, IplImage2 *rangeImage,
					   IplImage2 *domainCornerMap, IplImage2 *rangeCornerMap,
					   CorspMap *corspMap)
{
	int blockHeight = WINDOWSIZE_NCC;
	int blockWidth = WINDOWSIZE_NCC;
	IplImage2 *domainBlock = 0;
	uchar *cornerMapData = 0;
	int height, width, channels, mapStep;
	int i, j;
	int corspPosI = 0,corspPosJ = 0;
	int searchRange = SEARCH_RANGE;
	channels = domainImage->nChannels;
	mapStep = domainCornerMap->widthStep;
	height = domainImage->height;
	width = domainImage->width;
	// create the domain image block（模板图像）
	domainBlock = cvCreateImage2(cvSize2(blockWidth, blockHeight), IPL_DEPTH_8U, channels);
	InitializeImage(domainBlock);
	cornerMapData = (uchar *)domainCornerMap->imageData;
	for(i = 0; i < height; i++){
		for(j = 0; j < width; j++){
			if(cornerMapData[i*mapStep + j]){//参考帧i,j处有角点
				// make range image block
				MakeImageBlock(domainBlock, domainImage, i, j);
				// find the corresponding pixel position in range image
				// that will give the largest NCC in the search range
				if(FindCorrespodence(rangeImage, domainBlock, rangeCornerMap, i, j,
					&corspPosI, &corspPosJ, searchRange)){
					// store the corresponding position
					// & the domain image corner position
					UpdateCorrespMap(corspMap, i, j, corspPosI, corspPosJ);
				}
			}
		}
	}
	ReleaseImage(&domainBlock);
}
int FindCorrespodence(IplImage2 *rangeImage, IplImage2 *domainBlock,
					   IplImage2 *rangeCornerMap,
					   int cornerPosI, int cornerPosJ,
					   int *pCorspPosI, int *pCorspPosJ,
					   int searchRange)
{
	IplImage2 *rangeBlock = 0;
	int channels, blockHeight, blockWidth;
	int i, j;
	float value, maxNcc = -2.0;
	int r; // half searching range
	uchar *cornerMapData = 0;
	int mapStep = rangeCornerMap->widthStep;
	int height, width, iBegin, jBegin, iEnd, jEnd;
	height = rangeImage->height;
	width = rangeImage->width;
	// create the range image block（匹配图像，以cornerPosI，cornerPosJ（参考帧角点位置）为中心）
	blockHeight = domainBlock->height; // WINDOWSIZE_NCC
	blockWidth = domainBlock->width; // WINDOWSIZE_NCC
	channels = domainBlock->nChannels;
	rangeBlock = cvCreateImage2(cvSize2(blockWidth, blockHeight), IPL_DEPTH_8U, channels);
	InitializeImage(rangeBlock);
	r = searchRange / 2;
	cornerMapData = (uchar *)rangeCornerMap->imageData;
	// calculate NCC only in the search range
	// find the range image pixel that gives the maximum ncc
	// this pixel will be assigned as the correspondence
	// pixel of the domain image pixel
	//以cornerPosI，cornerPosJ（参考帧角点位置）为中心进行当前帧的角点的搜索
	iBegin = max(cornerPosI - r, 0);
	jBegin = max(cornerPosJ - r, 0);
	iEnd = min(cornerPosI + r + 1, height - 1);
	jEnd = min(cornerPosJ + r + 1, width - 1) ;
	for(i = iBegin; i < iEnd; i++){
		for(j = jBegin; j < jEnd; j++){
			// calculate NCC at only range image corner points
			if(cornerMapData[i*mapStep + j] ){//遍历该匹配图像的角点
				// make range image block
				MakeImageBlock(rangeBlock, rangeImage, i, j);
				// calculate ncc
				value = NCC(rangeBlock, domainBlock);
				// take the position that gives maximum NCC
				if(value >= maxNcc){
					maxNcc = value;//选择取最大value值时的角点（pCorspPosI,pCorspPosJ)与参考帧角点（cornerPosI,cornerPosJ）匹配
					*pCorspPosI = i;
					*pCorspPosJ = j;
				}
			}
		}
	}
	ReleaseImage(&rangeBlock);
	// if no corner in search range, do nothing
	if(maxNcc == -0.2){
		return 0;
	}else{
		return 1;
	}
}
//
// function : NCC
// usage : nccValue = NCC(block1, block2);
//---------------------------------------------
// This function returns ncc value between
// block1 and block2.
//
float NCC(IplImage2 *block1, IplImage2 *block2)
{
	int i, j, k;
	uchar *block1Data = 0, *block2Data = 0;
	int height = block1->height;
	int width = block1->width;
	int channels = block1->nChannels;
	int step = block1->widthStep;
	float *meanB1= (float*)malloc(sizeof(float )*channels);
	if ((meanB1 = (float *)malloc(sizeof(float)*channels)) != NULL)
	{
		printf("分配成功");
	}
	else {
		printf("分配失败");
		exit(0);

	}
	float *meanB2= (float*)malloc(sizeof(float)*channels);
	float *varB1= (float*)malloc(sizeof(float)*channels);
	float *varB2= (float*)malloc(sizeof(float)*channels);
	float *numerTerm= (float*)malloc(sizeof(float)*channels);
	float *denomTerm= (float*)malloc(sizeof(float)*channels);
	float ncc = 0;
	block1Data = (uchar *)block1->imageData;
	block2Data = (uchar *)block2->imageData;
	// initialize
	for(k = 0; k < channels; k++){
		meanB1[k] = 0;
		meanB2[k] = 0;
		varB1[k] = 0;
		varB2[k] = 0;
		numerTerm[k] = 0;
	}
	// calculate mean values
	for(i = 0; i < height; i++){
		for(j = 0; j < width; j++){
			for(k = 0; k < channels; k++){
				meanB1[k] += (float)block1Data[i*step + j*channels + k];
				meanB2[k] += (float)block2Data[i*step + j*channels + k];
			}
		}
	}
	for(k = 0; k < channels; k++){
		meanB1[k] = meanB1[k] / (height * width);
		meanB2[k] = meanB2[k] / (height * width);
	}
	for(i = 0; i < height; i++){
		for(j = 0; j < width; j++){
			for(k = 0; k < channels; k++){
				numerTerm[k] += ((float)block1Data[i*step + j*channels + k]
					- meanB1[k])
					* ((float)block2Data[i*step + j*channels + k]
					- meanB2[k]);
				varB1[k] += pow(((float)block1Data[i*step + j*channels + k]
					- meanB1[k]), 2);
				varB2[k] += pow(((float)block2Data[i*step + j*channels + k]
					- meanB2[k]), 2);
			}
		}
	}
	for(k = 0; k < channels; k++){
		denomTerm[k] = pow((varB1[k]*varB2[k]), (0.5));
		if(denomTerm[k] == 0){
			ncc += 0;
		}else{
			ncc += numerTerm[k] / denomTerm[k];
		}
	}	
	free(meanB1);
	free(meanB2);
	free(varB1);
	free(varB2);
	free(numerTerm);
	free(denomTerm);
	meanB1 = 0;
	meanB2 = 0;
	varB1 = 0;
	varB2 = 0;
	numerTerm = 0;
	denomTerm = 0;
	// we can calculate NCC for color image blocks as
	// the average of each color NCCs
	ncc = ncc / channels;
	return(ncc);
}
void UpdateCorrespMap(CorspMap *corspMap, int domainPosI, int domainPosJ,
					  int rangePosI, int rangePosJ)
{
	int len;
	len = corspMap->len;
	if(corspMap->len >= MAX_POINT_SIZE){
		printf("the length of corspMap > MAX_POINT_SIZE \n");
		len = MAX_POINT_SIZE - 1;
	}
	corspMap->rangeImagePositionI[len] = rangePosI;
	corspMap->rangeImagePositionJ[len] = rangePosJ;
	corspMap->domainImagePositionI[len] = domainPosI;
	corspMap->domainImagePositionJ[len] = domainPosJ;
	corspMap->len = len + 1;
}
void InitializeCorspMap(CorspMap *corspMap)
{
	int i;
	for( i = 0; i < MAX_POINT_SIZE; i++){
		corspMap->rangeImagePositionI[i] = 0;
		corspMap->rangeImagePositionJ[i] = 0;
		corspMap->domainImagePositionI[i] = 0;
		corspMap->domainImagePositionJ[i] = 0;
	}
	corspMap->len = 0;
}
void CopyCorspMap(CorspMap *dst, CorspMap *src)
{
	int i;
	for(i = 0; i < MAX_POINT_SIZE; i++){
		dst->rangeImagePositionI[i] = src->rangeImagePositionI[i];
		dst->rangeImagePositionJ[i] = src->rangeImagePositionJ[i];
		dst->domainImagePositionI[i] = src->domainImagePositionI[i];
		dst->domainImagePositionJ[i] = src->domainImagePositionJ[i];
	}
	dst->len = src->len;
}

/*void getCorspMap(CorspMap *corspMap, Point2 *point1_innext, Point2 *point2_temp, int len)
{
	int j = 0;
	for (int i = 0; i < len; i++)
	{
		if (point1_innext[i].x != -1 && point1_innext[i].y != -1 && point2_temp[i].x != -1 && point2_temp[i].y != -1)
		{
			corspMap->rangeImagePositionI[j] = point2_temp[i].y;
			corspMap->rangeImagePositionJ[j] = point2_temp[i].x;
			corspMap->domainImagePositionI[j] = point1_innext[i].y;
			corspMap->domainImagePositionJ[j] = point1_innext[i].x;
			if (j < MAX_POINT_SIZE - 2)
			{
				j++;
			}
			else
			{
				j++; //len最长是499
				break;
			}
		}
	}
	corspMap->len = j;  //粗匹配得到的点对数
}*/