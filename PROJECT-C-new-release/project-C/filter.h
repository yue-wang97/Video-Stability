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
// k取值0～1，k越大，滞后越小，对观测数据的滤波越小，对观测数据的随动性越大；

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
// alpha、beta取值0～1。alpha-beta滤波可看为在速度上进行一阶滞后滤波；
// alpha越大，对观测数据的滤波越小，对观测数据的随动性越大；
// beta越大，对速度的滤波越小，随动性越大。

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
// alpha取值0～1。MAS滤波可看为在alpha-beta滤波中用移动平均方法计算速度；
// alpha越大，对观测数据的滤波越小，对观测数据的随动性越大；
// MASF_WINDOW_SIZE越大，对速度的滤波越大。

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
// MRF_WINDOW_SIZE越大，对速度的滤波越大。

void MRFilter(float dx, float dt, MRF* pmr);
// assert: pmr!=NULL && dt>0
