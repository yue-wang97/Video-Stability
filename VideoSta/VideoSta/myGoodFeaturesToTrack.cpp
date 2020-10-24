void myGoodFeaturesToTrack(InputArray _image, OutputArray _corners,
	int maxCorners, double qualityLevel, double minDistance,
	InputArray _mask, int blockSize,
	bool useHarrisDetector, double harrisK)
{
	Mat image = _image.getMat(), mask = _mask.getMat();

	CV_Assert(qualityLevel > 0 && minDistance >= 0 && maxCorners >= 0);
	CV_Assert(mask.empty() || (mask.type() == CV_8UC1 && mask.size() == image.size()));

	Mat eig, tmp;
	if (useHarrisDetector)
		cornerHarris(image, eig, blockSize, 3, harrisK);
	else
		cornerMinEigenVal(image, eig, blockSize, 3);

	double maxVal = 0;
	minMaxLoc(eig, 0, &maxVal, 0, 0, mask);
	threshold(eig, eig, maxVal*qualityLevel, 0, THRESH_TOZERO);
	dilate(eig, tmp, Mat());

	Size imgsize = image.size();

	vector<const float*> tmpCorners;

	// collect list of pointers to features - put them into temporary image
	for (int y = 1; y < imgsize.height - 1; y++)
	{
		const float* eig_data = (const float*)eig.ptr(y);
		const float* tmp_data = (const float*)tmp.ptr(y);
		const uchar* mask_data = mask.data ? mask.ptr(y) : 0;

		for (int x = 1; x < imgsize.width - 1; x++)
		{
			float val = eig_data[x];
			if (val != 0 && val == tmp_data[x] && (!mask_data || mask_data[x]))
				tmpCorners.push_back(eig_data + x);
		}
	}

	sort(tmpCorners, greaterThanPtr<float>());
	vector<Point2f> corners;
	size_t i, j, total = tmpCorners.size(), ncorners = 0;

	if (minDistance >= 1)
	{
		// Partition the image into larger grids
		int w = image.cols;
		int h = image.rows;

		const int cell_size = cvRound(minDistance);
		const int grid_width = (w + cell_size - 1) / cell_size;
		const int grid_height = (h + cell_size - 1) / cell_size;

		std::vector<std::vector<Point2f> > grid(grid_width*grid_height);

		minDistance *= minDistance;

		for (i = 0; i < total; i++)
		{
			int ofs = (int)((const uchar*)tmpCorners[i] - eig.data);
			int y = (int)(ofs / eig.step);
			int x = (int)((ofs - y*eig.step) / sizeof(float));

			bool good = true;

			int x_cell = x / cell_size;
			int y_cell = y / cell_size;

			int x1 = x_cell - 1;
			int y1 = y_cell - 1;
			int x2 = x_cell + 1;
			int y2 = y_cell + 1;

			// boundary check
			x1 = std::max(0, x1);
			y1 = std::max(0, y1);
			x2 = std::min(grid_width - 1, x2);
			y2 = std::min(grid_height - 1, y2);

			for (int yy = y1; yy <= y2; yy++)
			{
				for (int xx = x1; xx <= x2; xx++)
				{
					vector <Point2f> &m = grid[yy*grid_width + xx];

					if (m.size())
					{
						for (j = 0; j < m.size(); j++)
						{
							float dx = x - m[j].x;
							float dy = y - m[j].y;

							if (dx*dx + dy*dy < minDistance)
							{
								good = false;
								goto break_out;
							}
						}
					}
				}
			}

		break_out:

			if (good)
			{
				// printf("%d: %d %d -> %d %d, %d, %d -- %d %d %d %d, %d %d, c=%d\n",
				//    i,x, y, x_cell, y_cell, (int)minDistance, cell_size,x1,y1,x2,y2, grid_width,grid_height,c);
				grid[y_cell*grid_width + x_cell].push_back(Point2f((float)x, (float)y));

				corners.push_back(Point2f((float)x, (float)y));
				++ncorners;

				if (maxCorners > 0 && (int)ncorners == maxCorners)
					break;
			}
		}
	}
	else
	{
		for (i = 0; i < total; i++)
		{
			int ofs = (int)((const uchar*)tmpCorners[i] - eig.data);
			int y = (int)(ofs / eig.step);
			int x = (int)((ofs - y*eig.step) / sizeof(float));

			corners.push_back(Point2f((float)x, (float)y));
			++ncorners;
			if (maxCorners > 0 && (int)ncorners == maxCorners)
				break;
		}
	}

	Mat(corners).convertTo(_corners, _corners.fixedType() ? _corners.type() : CV_32F);

	/*
	for( i = 0; i < total; i++ )
	{
	int ofs = (int)((const uchar*)tmpCorners[i] - eig.data);
	int y = (int)(ofs / eig.step);
	int x = (int)((ofs - y*eig.step)/sizeof(float));

	if( minDistance > 0 )
	{
	for( j = 0; j < ncorners; j++ )
	{
	float dx = x - corners[j].x;
	float dy = y - corners[j].y;
	if( dx*dx + dy*dy < minDistance )
	break;
	}
	if( j < ncorners )
	continue;
	}

	corners.push_back(Point2f((float)x, (float)y));
	++ncorners;
	if( maxCorners > 0 && (int)ncorners == maxCorners )
	break;
	}
	*/
}