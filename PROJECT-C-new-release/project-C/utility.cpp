//
// file : utility.cpp
//----------------------------------------
// this file contains utility functions to
// deal with general processes.
//
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "defs.h"
#include "utility.h"

void Array2CvMat(float *arr, CvMat2 *cvArr, int row, int column)
{
	int i, j;
	for(i = 0; i < row; i++){
		for(j = 0; j < column; j++){
			cvmSet2(cvArr, i, j, arr[i*column + j]);
		}
	}
}
void CvMat2Array(CvMat2 *cvArr, float *arr, int row, int column)
{
	int i, j;
	for(i = 0; i < row; i++){
		for(j = 0; j < column; j++){
			arr[i*column + j] = cvmGet(cvArr, i, j);
		}
	}
}


//
// function : MakeImageBlock
// usage : MakeImageBlock(block, image, centPosI, centPosJ);
// ------------------------------------------------------------
// This function copies a block region of the image into a block
// for example, if block size is 3 by 3 and the position of the block
// is i, j on the image. the resultant block will be 3 by 3 and
// the block will be copied by image(i-1, j-1) ... image(i+1, j+1).
//
void MakeImageBlock(IplImage2 *block, IplImage2 *image,
					int centPosI, int centPosJ)
{
	uchar *blockData = 0, *imageData = 0;
	int blockHeight, blockWidth, imageHeight, imageWidth;
	int blockStep, channels, imageStep;
	int i, j, k, posI, posJ;
	blockHeight = block->height;
	blockWidth = block->width;
	imageHeight = image->height;	
	imageWidth = image->width;
	channels = block->nChannels;
	blockStep = block->widthStep;
	imageStep = image->widthStep;
	blockData = (uchar *)block->imageData;
	imageData = (uchar *)image->imageData;
	for(i = 0; i < blockHeight; i++){
		for(j = 0; j < blockWidth; j++){
			for(k = 0; k < channels; k++){
				posI = centPosI + i - blockHeight / 2;
				posJ = centPosJ + j - blockWidth / 2;
				posI = min(max(posI, 0), imageHeight - 1);
				posJ = min(max(posJ, 0), imageWidth - 1);
				blockData[i*blockStep + j*channels + k]
					= imageData[posI*imageStep + posJ*channels + k];
			}
		}
	}
}


void InitializeimageData(uchar * img,int height,int width,int step) {
	int i, j, k;
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			for (k = 0; k < 1; k++) {
				img[i*step + j + k]
					= 0;
			}
		}
	}
}
void   Array2Image(IplImage2 *img, char *b)//将
{
	int i, j, k;
	int height = img->height;
	int width = img->width;
	int channels = img->nChannels;
	int step = img->widthStep;
	uchar *imageData = (uchar*)img->imageData;
	uchar* a = (uchar*)b;
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			for (k = 0; k < channels; k++) {
				imageData[i*step + j*channels + k]
					= a[i*step + j*channels + k];
			}
		}
	}
}

//
// function : TransformImage
// usage : TransformImage(inImage, Data,dx,dy);
// ---------------------------------------------
// This function transforms input image using dx,dy
//
void TransformImage1(IplImage2 *inImage, char * Data, int dx, int dy)//dx,dy为平移量
{
	uchar *inData;
	int height, width, step;
	height = inImage->height;
	width = inImage->width;
	step = inImage->widthStep;
	inData = (uchar*)inImage->imageData;
	uchar* outData = (uchar*)Data;
	InitializeimageData(outData,height,width,step);
	int i, j;
	if ((dx >= 0) && (dy >= 0))
	{
		for (i = 0; i<height - dy; i++)
			for (j = 0; j<width - dx; j++)
			{
				outData[i*step + j] = inData[(i + dy)*step + j + dx];
			}
	}
	else if ((dx >= 0) && (dy<0))
	{
		for (i = height - 1; i >= abs(dy); i--)
			for (j = 0; j<width - dx; j++)
			{
				outData[i*step + j] = inData[(i + dy)*step + j + dx];
			}
	}
	else if ((dx<0) && (dy >= 0))
	{
		for (i = 0; i<height - dy; i++)
			for (j = width - 1; j >= abs(dx); j--)
			{
				outData[i*step + j] = inData[(i + dy)*step + j + dx];
			}
	}
	else
	{
		for (i = height - 1; i >= abs(dy); i--)
			for (j = width - 1; j >= abs(dx); j--)
			{
				outData[i*step + j] = inData[(i + dy)*step + j + dx];
			}
	}
}

