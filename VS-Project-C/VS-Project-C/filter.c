#include "filter.h"

#define DX_PREDICT(dx0, v0, dt) ( (dx0) + (dt)*(v0) )
#define DX_PREDICTED_ERROR(dx_predicted, dx) ( (dx_predicted) - (dx) )
#define ALPHA_FILTER_1(alpha1, dx_predicted_error) ( (alpha1) * ( dx_predicted_error ) )
#define ALPHA_FILTER_2(alpha1, dx0, v0, dx, dt)    ( (alpha1) * ( DX_PREDICTED_ERROR( (DX_PREDICT(dx0, v0, dt)),  dx ) ) )

///////////////////////// First Order Lag Filter /////////////////////////
int InitializeFOLFilter(float k, FOLF* pfol)
{
	if (!pfol)
		return -1;

	pfol->k = k;
	pfol->x = 0.0f;

	pfol->num = 0;

	return 0;
}

// assert: pfol!=NULL
void FOLFilter(float x, FOLF* pfol)
{
	if (pfol->num == 0)
	{
		pfol->x = x;

		pfol->num = 1;
		return;
	}

	float k = pfol->k;
	float x_old = pfol->x;

	pfol->x = k * x + (1 - k) * x_old;

	return;
}
//////////////////////////////////////////////////

///////////////////////// Alpha Beta filter /////////////////////////
int InitializeAlphaBeta(float alpha, float beta, AlphaBeta* pab) 
{
	if (!pab)
		return -1;

	pab->alpha1 =1.0f - alpha;
	pab->beta = beta;

	pab->num = 0;

	pab->dx = 0.0f;
	pab->v = 0.0f;

	return 0;
}

// assert: pab!=NULL && dt>0
void AlphaBetaFilter(float dx, float dt, AlphaBeta* pab)
{
	if (pab->num == 0)
	{
		pab->dx = 0.0f;
		pab->v = dx / dt;

		pab->num = 1;
		return;
	}

	float dx0 = pab->dx;
	float v0 = pab->v;

	float dx1 = DX_PREDICT(dx0, v0, dt);					//pridicted value
	float err1 = DX_PREDICTED_ERROR(dx1, dx);			//prediction error 

	pab->dx = ALPHA_FILTER_1(pab->alpha1, err1);		//filtering
	pab->v = v0 - pab->beta * err1 / dt ;
	
	return;
}

// assert: pab!=NULL
int AlphaBetaPredict(AlphaBeta* pab, float dt, float* dx_predicted)
{
	if (pab->num == 0)
	{
		*dx_predicted = 0.0f;
		return -1;
	}

	*dx_predicted = DX_PREDICT(pab->dx, pab->v, dt);
	return 0;
}
//////////////////////////////////////////////////



///////////////////////// Moving Average Speed Filter /////////////////////////
int InitializeMASFilter(float alpha, MASF* pmas)
{
	if (!pmas)
		return -1;

	pmas->alpha1 = 1.0f - alpha;

	pmas->head = 0; pmas->rear = 0; pmas->num = 0;
	pmas->sum_v = 0.0f;
	for (int i = 0; i < MASF_WINDOW_SIZE; i++)
		pmas->v0[i] = 0.0f;

	pmas->dx = 0.0f;
	pmas->v_average = 0.0f;

	return 0;
}

// assert: pmas!=NULL && dt>0
void MASFilter(float dx, float dt, MASF* pmas)
{
	float v = dx / dt;

	if (pmas->num == 0)
	{
		pmas->head = 0; pmas->rear = 0; 	pmas->num = 1;
		pmas->sum_v = v;

		pmas->v0[0] = v;

		pmas->dx = 0.0f;
		pmas->v_average = v;

		return;
	}

	// Alpha filter 
	float dx_filtered = ALPHA_FILTER_2(pmas->alpha1, pmas->dx, pmas->v_average, dx, dt);

	// remove the old point if necessary, add the new point
	if (pmas->num == MASF_WINDOW_SIZE)
	{
		unsigned char rear = pmas->rear;
		pmas->sum_v -= pmas->v0[rear];
		pmas->rear = (rear+1) % MASF_WINDOW_SIZE;
	}
	else
		pmas->num++;

	//add new point
	pmas->head = (pmas->head+1) % MASF_WINDOW_SIZE;
	pmas->v0[pmas->head] = v;

	pmas->sum_v += v;

	// average speed
	pmas->dx = dx_filtered;
	pmas->v_average = pmas->sum_v / pmas->num;

	return;
}

// assert: pmas!=NULL
int MASFPredict(MASF* pmas, float dt, float* dx_predicted)
{
	if (pmas->num == 0)
	{
		*dx_predicted = 0.0f;
		return -1;
	}

	*dx_predicted = DX_PREDICT(pmas->dx, pmas->v_average, dt);
	return 0;
}

//////////////////////////////////////////////////

///////////////////////// Moving Regression Filter /////////////////////////
int InitializeMRFilter(float t_max, float dt_min, MRF* pmr)
{
	if (!pmr)
		return -1;

	pmr->t_max  = t_max;	 
	pmr->dt_min = dt_min;

	pmr->x_current = 0.0f; pmr->t_current = 0.0f;

	float *x = pmr->x;
	float *t = pmr->t;
	for (int i = 0; i < MRF_WINDOW_SIZE; i++)
	{
		x[i] = 0.0f;
		t[i] = 0.0f;
	}

	pmr->sum_xt = 0.0f;
	pmr->sum_t2 = 0.0f;
	pmr->sum_x = 0.0f;
	pmr->sum_t = 0.0f;

	pmr->x_diff = 0.0f;

	//由于是增量输入，第一点已经存入. As a result, num is always greater than 0
	pmr->head = 0; pmr->rear = 0;
	pmr->num = 1; 

	return 0;
}

// assert: pmr!=NULL
void MRFAdjust(MRF* pmr)
{
	float *t = pmr->t;
	float t0 = t[pmr->head];
	if ( t0 < pmr->t_max )
		return;

	float *x = pmr->x;
	float x0 = x[pmr->head];

	pmr->x_current -= x0;
	pmr->t_current -= t0;

	int num = pmr->num;
	unsigned char p = pmr->rear;
	int i;
	for( i=0; i<num; i++)
	{
		x[p] -= x0;
		t[p] -= t0;

		p = (p+1) % MRF_WINDOW_SIZE;
	}

	pmr->sum_xt -= ( t0*pmr->sum_x + x0*pmr->sum_t - num*x0*t0 );
	pmr->sum_t2 -= ( 2*t0*pmr->sum_t - num*t0*t0 );
	pmr->sum_x  -= num*x0;
	pmr->sum_t  -= num*t0;

	return;
}

int MRFAddPoint(MRF* pmr)
// assert: pmr!=NULL && dt>0
{
	float x_current = pmr->x_current;
	float t_current = pmr->t_current;

	float *x = pmr->x;
	float *t = pmr->t;

	unsigned char head = pmr->head;
	if ( (t_current - t[head]) < pmr->dt_min ) // the gap of t is too small, the new point is not added into the quene
		return -1;

	if (pmr->num == MRF_WINDOW_SIZE) // remove the last point
	{
		unsigned char rear = pmr->rear;

		float x0 = x[rear];
		float t0 = t[rear];

		pmr->sum_xt -= x0*t0;
		pmr->sum_t2 -= t0*t0;
		pmr->sum_x  -= x0;
		pmr->sum_t  -= t0;

		pmr->rear = (rear + 1) % MRF_WINDOW_SIZE;
	}
	else
		pmr->num++;

	//add the new point
	head = (head + 1) % MRF_WINDOW_SIZE;
	pmr->head = head;

	x[head] = x_current;
	t[head] = t_current;

	pmr->sum_xt += x_current*t_current;
	pmr->sum_t2 += t_current*t_current;
	pmr->sum_x  += x_current;
	pmr->sum_t  += t_current;

	return 0;
}

// assert: pmr!=NULL && dt>0
void MRFilter(float dx, float dt, MRF* pmr)
{
	float x_current = pmr->x_current + dx;
	pmr->x_current  = x_current;
	float t_current = pmr->t_current + dt;
	pmr->t_current  = t_current;;

	int num;
	float sum_xt, sum_t2, sum_x, sum_t;
	if (MRFAddPoint(pmr) ) 
	{
		// the point has not been added. add it here
		sum_xt = pmr->sum_xt + x_current*t_current;
		sum_t2 = pmr->sum_t2 + t_current*t_current;
		sum_x  = pmr->sum_x  + x_current;
		sum_t  = pmr->sum_t  + t_current;

		num = pmr->num + 1;
	}
	else
	{
		// the point has been added.
		sum_xt = pmr->sum_xt;
		sum_t2 = pmr->sum_t2;
		sum_x = pmr->sum_x;
		sum_t = pmr->sum_t;

		num = pmr->num;
	}
	
	// Regression
	float ave_x = sum_x / num;
	float ave_t = sum_t / num;

	float b = (sum_xt - num*ave_x*ave_t) / (sum_t2 - num*ave_t*ave_t);
	float a = ave_x - b*ave_t;

	// dx
	pmr->x_diff = (b*t_current + a) - x_current;

	// adjust
	MRFAdjust(pmr);

	return;
}

//////////////////////////////////////////////////