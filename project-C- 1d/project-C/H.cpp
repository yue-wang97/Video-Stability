#include"H.h"

void SetH(float* h, int row, int col, float value)
{
	int type, step;
	type = CV_MAT_TYPE(CV_32FC1);
	step = CV_ELEM_SIZE(type);
	((float*)(void*)(h + step*row))[col] = (float)value;
}
float  GetH(float *h, int row, int col)
{
	int type, step;
	type = CV_MAT_TYPE(CV_32FC1);
	step = CV_ELEM_SIZE(type);
	return ((float*)(void*)(h + step*row))[col];
}
void copyH(float *a, float* b)
{
	int n = sizeof(a) / sizeof(a[0]);
	int i = 0;
	for (i; i < n; i++)
		b[i] = a[i];
}