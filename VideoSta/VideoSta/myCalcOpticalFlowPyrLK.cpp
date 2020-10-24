void myCalcOpticalFlowPyrLK(InputArray _prevImg, InputArray _nextImg,
	InputArray _prevPts, InputOutputArray _nextPts,
	OutputArray _status, OutputArray _err,
	Size winSize, int maxLevel,
	TermCriteria criteria,
	int flags, double minEigThreshold)
{
	Mat prevPtsMat = _prevPts.getMat();
	const int derivDepth = DataType<cv::detail::deriv_type>::depth;

	CV_Assert(maxLevel >= 0 && winSize.width > 2 && winSize.height > 2);

	int level = 0, i, npoints;
	CV_Assert((npoints = prevPtsMat.checkVector(2, CV_32F, true)) >= 0);

	if (npoints == 0)
	{
		_nextPts.release();
		_status.release();
		_err.release();
		return;
	}

	if (!(flags & OPTFLOW_USE_INITIAL_FLOW))
		_nextPts.create(prevPtsMat.size(), prevPtsMat.type(), -1, true);

	Mat nextPtsMat = _nextPts.getMat();
	CV_Assert(nextPtsMat.checkVector(2, CV_32F, true) == npoints);

	const Point2f* prevPts = (const Point2f*)prevPtsMat.data;
	Point2f* nextPts = (Point2f*)nextPtsMat.data;

	_status.create((int)npoints, 1, CV_8U, -1, true);
	Mat statusMat = _status.getMat(), errMat;
	CV_Assert(statusMat.isContinuous());
	uchar* status = statusMat.data;
	float* err = 0;

	for (i = 0; i < npoints; i++)
		status[i] = true;

	if (_err.needed())
	{
		_err.create((int)npoints, 1, CV_32F, -1, true);
		errMat = _err.getMat();
		CV_Assert(errMat.isContinuous());
		err = (float*)errMat.data;
	}

	vector<Mat> prevPyr, nextPyr;
	int levels1 = -1;
	int lvlStep1 = 1;
	int levels2 = -1;
	int lvlStep2 = 1;

	if (_prevImg.kind() == _InputArray::STD_VECTOR_MAT)
	{
		_prevImg.getMatVector(prevPyr);

		levels1 = int(prevPyr.size()) - 1;
		CV_Assert(levels1 >= 0);

		if (levels1 % 2 == 1 && prevPyr[0].channels() * 2 == prevPyr[1].channels() && prevPyr[1].depth() == derivDepth)
		{
			lvlStep1 = 2;
			levels1 /= 2;
		}

		// ensure that pyramid has reqired padding
		if (levels1 > 0)
		{
			Size fullSize;
			Point ofs;
			prevPyr[lvlStep1].locateROI(fullSize, ofs);
			CV_Assert(ofs.x >= winSize.width && ofs.y >= winSize.height
				&& ofs.x + prevPyr[lvlStep1].cols + winSize.width <= fullSize.width
				&& ofs.y + prevPyr[lvlStep1].rows + winSize.height <= fullSize.height);
		}

		if (levels1 < maxLevel)
			maxLevel = levels1;
	}

	if (_nextImg.kind() == _InputArray::STD_VECTOR_MAT)
	{
		_nextImg.getMatVector(nextPyr);

		levels2 = int(nextPyr.size()) - 1;
		CV_Assert(levels2 >= 0);

		if (levels2 % 2 == 1 && nextPyr[0].channels() * 2 == nextPyr[1].channels() && nextPyr[1].depth() == derivDepth)
		{
			lvlStep2 = 2;
			levels2 /= 2;
		}

		// ensure that pyramid has reqired padding
		if (levels2 > 0)
		{
			Size fullSize;
			Point ofs;
			nextPyr[lvlStep2].locateROI(fullSize, ofs);
			CV_Assert(ofs.x >= winSize.width && ofs.y >= winSize.height
				&& ofs.x + nextPyr[lvlStep2].cols + winSize.width <= fullSize.width
				&& ofs.y + nextPyr[lvlStep2].rows + winSize.height <= fullSize.height);
		}

		if (levels2 < maxLevel)
			maxLevel = levels2;
	}

	if (levels1 < 0)
		maxLevel = buildOpticalFlowPyramid(_prevImg, prevPyr, winSize, maxLevel, false);

	if (levels2 < 0)
		maxLevel = buildOpticalFlowPyramid(_nextImg, nextPyr, winSize, maxLevel, false);

	if ((criteria.type & TermCriteria::COUNT) == 0)
		criteria.maxCount = 30;
	else
		criteria.maxCount = std::min(std::max(criteria.maxCount, 0), 100);
	if ((criteria.type & TermCriteria::EPS) == 0)
		criteria.epsilon = 0.01;
	else
		criteria.epsilon = std::min(std::max(criteria.epsilon, 0.), 10.);
	criteria.epsilon *= criteria.epsilon;

	// dI/dx ~ Ix, dI/dy ~ Iy
	Mat derivIBuf;
	if (lvlStep1 == 1)
		derivIBuf.create(prevPyr[0].rows + winSize.height * 2, prevPyr[0].cols + winSize.width * 2, CV_MAKETYPE(derivDepth, prevPyr[0].channels() * 2));

	for (level = maxLevel; level >= 0; level--)
	{
		Mat derivI;
		if (lvlStep1 == 1)
		{
			Size imgSize = prevPyr[level * lvlStep1].size();
			Mat _derivI(imgSize.height + winSize.height * 2,
				imgSize.width + winSize.width * 2, derivIBuf.type(), derivIBuf.data);
			derivI = _derivI(Rect(winSize.width, winSize.height, imgSize.width, imgSize.height));
			calcSharrDeriv(prevPyr[level * lvlStep1], derivI);
			copyMakeBorder(derivI, _derivI, winSize.height, winSize.height, winSize.width, winSize.width, BORDER_CONSTANT | BORDER_ISOLATED);
		}
		else
			derivI = prevPyr[level * lvlStep1 + 1];

		CV_Assert(prevPyr[level * lvlStep1].size() == nextPyr[level * lvlStep2].size());
		CV_Assert(prevPyr[level * lvlStep1].type() == nextPyr[level * lvlStep2].type());

#ifdef HAVE_TEGRA_OPTIMIZATION
		typedef tegra::LKTrackerInvoker<cv::detail::LKTrackerInvoker> LKTrackerInvoker;
#else
		typedef cv::detail::LKTrackerInvoker LKTrackerInvoker;
#endif

		parallel_for_(Range(0, npoints), LKTrackerInvoker(prevPyr[level * lvlStep1], derivI,
			nextPyr[level * lvlStep2], prevPts, nextPts,
			status, err,
			winSize, criteria, level, maxLevel,
			flags, (float)minEigThreshold));
	}
}