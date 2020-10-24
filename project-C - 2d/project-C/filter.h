#pragma once
/**
** filter.h
**/

///////////////////////// First Order Lag Filter /////////////////////////
typedef struct {
	float k;		//alpha value

	unsigned char num;

	float x;	//estimated dx 
} FOLF;

int InitializeFOLFilter(float k, FOLF* pfol);
// kȡֵ0��1��kԽ���ͺ�ԽС���Թ۲����ݵ��˲�ԽС���Թ۲����ݵ��涯��Խ��

void FOLFilter(float x, FOLF* pfol);
// assert: pfol!=NULL

///////////////////////// Alpha Beta filter /////////////////////////
typedef struct {
	float alpha1;		//(1-alpha) value
	float beta;			//beta value

	unsigned char num;

	float dx;	//estimated dx 
	float v;	//estimated v 
} AlphaBeta;

int InitializeAlphaBeta(float alpha, float beta, AlphaBeta* pab);
// alpha��betaȡֵ0��1��alpha-beta�˲��ɿ�Ϊ���ٶ��Ͻ���һ���ͺ��˲���
// alphaԽ�󣬶Թ۲����ݵ��˲�ԽС���Թ۲����ݵ��涯��Խ��
// betaԽ�󣬶��ٶȵ��˲�ԽС���涯��Խ��

void AlphaBetaFilter(float dx, float dt, AlphaBeta* pab);
// assert: pab!=NULL && dt>0

int AlphaBetaPredict(AlphaBeta* pab, float dt, float* dx_predicted);
// assert: pab!=NULL

///////////////////////// Moving Average Speed Filter /////////////////////////
#define MASF_WINDOW_SIZE		10
typedef struct {
	float alpha1;			//(1-alpha) value
	
	unsigned char head, rear, num;
	float sum_v;
	float v0[MASF_WINDOW_SIZE];

	float dx;			//estimated dx 
	float v_average;		//average v 
} MASF;

int InitializeMASFilter(float alpha, MASF* pmas);
// alphaȡֵ0��1��MAS�˲��ɿ�Ϊ��alpha-beta�˲������ƶ�ƽ�����������ٶȣ�
// alphaԽ�󣬶Թ۲����ݵ��˲�ԽС���Թ۲����ݵ��涯��Խ��
// MASF_WINDOW_SIZEԽ�󣬶��ٶȵ��˲�Խ��

void MASFilter(float dx, float dt, MASF* pmas);
// assert: pmas!=NULL && dt>0

int MASFPredict(MASF* pmas, float dt, float* dx_predicted);
// assert: pmas!=NULL

///////////////////////// Moving Regression Filter /////////////////////////
#define MRF_WINDOW_SIZE		5
typedef struct {
	float t_max, dt_min;

	float x_current, t_current;

	unsigned char head, rear, num;
	float x[MRF_WINDOW_SIZE];
	float t[MRF_WINDOW_SIZE];

	float sum_xt, sum_t2, sum_x, sum_t;

	float x_diff;	//estimated x_diff
} MRF;

int InitializeMRFilter(float t_max, float dt_min, MRF* pmr);
// MRF_WINDOW_SIZEԽ�󣬶��ٶȵ��˲�Խ��

void MRFilter(float dx, float dt, MRF* pmr);
// assert: pmr!=NULL && dt>0
