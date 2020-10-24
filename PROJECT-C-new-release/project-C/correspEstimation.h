//
// file : correspEstimation.h
//-----------------------------
// this file contains functions for correspondence
// estimation
//9 256 500  18 128 300
#define WINDOWSIZE_NCC 9
#define SEARCH_RANGE   256
#define MAX_POINT_SIZE 500//�������ƥ��300����
typedef struct{
	int len;//ƥ��ǵ�ĸ���
	int rangeImagePositionI[MAX_POINT_SIZE];//��ǰ֡�ǵ�����y
	int rangeImagePositionJ[MAX_POINT_SIZE];//��ǰ֡�ǵ�����x
	int domainImagePositionI[MAX_POINT_SIZE];//�ο�֡�ǵ�����y
	int domainImagePositionJ[MAX_POINT_SIZE];//�ο�֡�ǵ�����x
}CorspMap;//�����洢��Ӧƥ��Ľǵ�
void CorrespEstimation(IplImage2 *domainImage, IplImage2 *rangeImage,
					   IplImage2 *domainCornerMap, IplImage2 *rangeCornerMap,
					   CorspMap *corspMap);//�Բο�֡�͵�ǰ֡�Ľǵ���д�ƥ��
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