#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "cv.h"
#include "cxcore.h"
#include "highgui.h"


using namespace std;
using namespace cv;

struct greaterThanPtr :
	public std::binary_function<const float *, const float *, bool>
{
	bool operator () (const float * a, const float * b) const
		// Ensure a fully deterministic result of the sort
	{
		return (*a > *b) ? true : (*a < *b) ? false : (a > b);
	}
};

void goodFeatures(InputArray image, OutputArray corners,int maxCorners, double qualityLevel, double minDistance,
	InputArray mask = noArray(), int blockSize = 3,
	bool useHarrisDetector = false, double k = 0.04);
