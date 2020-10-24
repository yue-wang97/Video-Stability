#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "goodfeatures.h"

using namespace std;
using namespace cv;


void cornerMinEigenVal(Mat& src, Mat& eigenv)
{
	int depth = src.depth();
	double scale = (double)6;
	if (depth == CV_8U)
		scale *= 255.;
	scale = 1. / scale;
	Mat Dx, Dy;   //保存每个像素点的水平方向和垂直方向的一阶差分  

	Sobel66(src, Dx, CV_32F, 1, 0,3, scale, 0);
	Sobel66(src, Dy, CV_32F, 0, 1, 3, scale, 0);
	//getSobelKernels(src, Dx, 1, 0, 0,CV_32F);
	//getSobelKernels(src, Dy, 0, 1, 0, CV_32F);
	Size size = src.size();
	Mat cov(size, CV_32FC3);
	int i, j;
	for (i = 0; i < size.height; i++)
	{
		float* cov_data = (float*)(cov.data + i*cov.step);
		const float* dxdata = (const float*)(Dx.data + i*Dx.step);
		const float* dydata = (const float*)(Dy.data + i*Dy.step);
		for (j = 0; j < size.width; j++)
		{
			float dx = dxdata[j];
			float dy = dydata[j];
			cov_data[j * 3] = dx*dx;  //第一个通道存dx*dx,即M矩阵左上角的元素  
			cov_data[j * 3 + 1] = dx*dy;//第二个通道存dx*dy,即M矩阵左下角和右上角的元素  
			cov_data[j * 3 + 2] = dy*dy;//第三个通道存dy*dy,即M矩阵右下角的元素  
		}
	}
	/*boxFilter(cov, cov, cov.depth(), Size(3, 3), //计算邻域上的差分相关矩阵（block_size×block_size）  
		Point(-1, -1), false, borderType);*/
	IplImage cov_img = IplImage(cov);
	box_filter(&cov_img, &cov_img);
	calcHarris(cov, eigenv, 0.);
}

void calcHarris(const Mat& _cov, Mat& _dst, double k)
{
	int i, j;
	Size size = _cov.size();
	if (_cov.isContinuous() && _dst.isContinuous())
	{
		size.width *= size.height;
		size.height = 1;
	}

	for (i = 0; i < size.height; i++)
	{
		const float* cov = (const float*)(_cov.data + _cov.step*i);
		float* dst = (float*)(_dst.data + _dst.step*i);
	//	float dst[800];
		j = 0;
		for (; j < size.width; j++)
		{
			float a = cov[j * 3];
			float b = cov[j * 3 + 1];
			float c = cov[j * 3 + 2];
			dst[j] = (float)(a*c - b*b - k*(a + c)*(a + c));  //计算每个像素对应角点响应函数R  
		}
	}
}

void getSobelKernels66(OutputArray _kx, OutputArray _ky, int dx, int dy, int _ksize, bool normalize, int ktype)
{
	int i, j, ksizeX = _ksize, ksizeY = _ksize;
	if (ksizeX == 1 && dx > 0)
		ksizeX = 3;
	if (ksizeY == 1 && dy > 0)
		ksizeY = 3;

	CV_Assert(ktype == CV_32F || ktype == CV_64F);

	_kx.create(ksizeX, 1, ktype, -1, true);
	_ky.create(ksizeY, 1, ktype, -1, true);
	Mat kx = _kx.getMat();
	Mat ky = _ky.getMat();

	if (_ksize % 2 == 0 || _ksize > 31)
		CV_Error(CV_StsOutOfRange, "The kernel size must be odd and not larger than 31");
	vector<int> kerI(std::max(ksizeX, ksizeY) + 1);

	CV_Assert(dx >= 0 && dy >= 0 && dx + dy > 0);

	for (int k = 0; k < 2; k++)
	{
		Mat* kernel = k == 0 ? &kx : &ky;
		int order = k == 0 ? dx : dy;
		int ksize = k == 0 ? ksizeX : ksizeY;

		CV_Assert(ksize > order);

		if (ksize == 1)
			kerI[0] = 1;
		else if (ksize == 3)
		{
			if (order == 0)
				kerI[0] = 1, kerI[1] = 2, kerI[2] = 1;
			else if (order == 1)
				kerI[0] = -1, kerI[1] = 0, kerI[2] = 1;
			else
				kerI[0] = 1, kerI[1] = -2, kerI[2] = 1;
		}
		else
		{
			int oldval, newval;
			kerI[0] = 1;
			for (i = 0; i < ksize; i++)
				kerI[i + 1] = 0;

			for (i = 0; i < ksize - order - 1; i++)
			{
				oldval = kerI[0];
				for (j = 1; j <= ksize; j++)
				{
					newval = kerI[j] + kerI[j - 1];
					kerI[j - 1] = oldval;
					oldval = newval;
				}
			}

			for (i = 0; i < order; i++)
			{
				oldval = -kerI[0];
				for (j = 1; j <= ksize; j++)
				{
					newval = kerI[j - 1] - kerI[j];
					kerI[j - 1] = oldval;
					oldval = newval;
				}
			}
		}

		Mat temp(kernel->rows, kernel->cols, CV_32S, &kerI[0]);
		double scale = !normalize ? 1. : 1. / (1 << (ksize - order - 1));
		temp.convertTo(*kernel, ktype, scale);
	}
}

void Sobel66(InputArray _src, OutputArray _dst, int ddepth, int dx, int dy, int ksize, double scale, double delta, int borderType)
{
	Mat src = _src.getMat();
	if (ddepth < 0)
		ddepth = src.depth();
	_dst.create(src.size(), CV_MAKETYPE(ddepth, src.channels()));
	Mat dst = _dst.getMat();
	int ktype = std::max(CV_32F, std::max(ddepth, src.depth()));

	Mat kx, ky;
	getSobelKernels66(kx, ky, dx, dy, ksize, false, ktype);
	if (scale != 1)
	{
		// usually the smoothing part is the slowest to compute,
		// so try to scale it instead of the faster differenciating part
		if (dx == 0)
			kx *= scale;
		else
			ky *= scale;
	}
	sepFilter2D(src, dst, ddepth, kx, ky, Point(-1, -1), delta, borderType);
}

void box_filter(IplImage* img,IplImage* result){
     //init part
     CvScalar s;
     int width = img->width, height = img->height;
     int m_w = 5,m_h = 5;//window_size
     int boxwidth = width - m_w, boxheight = height - m_h;
     int *sum = (int*)malloc(boxwidth *boxheight*sizeof(double));
     int *buff= (int*)malloc(width*sizeof(double));
     memset(sum,0,boxwidth *boxheight*sizeof(int));
     memset(buff,0,width*sizeof(int));    
     //set buff:from 0 to 4 rows,per col
     int x,y,j;
     for(y=0; y<m_h; y++){
         for(x=0; x<width; x++){
             uchar pixel = CV_IMAGE_ELEM(img,uchar,y,x);
             buff[x] += pixel;
             //printf("%d:%d\n",x,buff[x]);
         }    
     }
     for(y=0; y<height - m_h;y++){
         int Xsum = 0;
 
         for(j=0; j<m_w; j++){
             Xsum += buff[j];//sum of pixel from (0,0) to (m_h,m_w) (also x = 0)
         }
 
         for(x=0; x<boxwidth; x++){
             if(x!=0){
                 Xsum = Xsum-buff[x-1]+buff[m_w-1+x];//Xsum:sum of cols range from x to x+m_w ,rows range from 0 to 4
             }
            sum[y*boxwidth + x] = (float) Xsum;        
         }
 
         for(x=0; x<width; x++){
             uchar pixel = CV_IMAGE_ELEM(img,uchar,y,x);//img[y *width + x];    
             uchar pixel2= CV_IMAGE_ELEM(img,uchar,y+m_h,x);//img[(y+mheight) *width + x];    
             buff[x] = buff[x] - pixel + pixel2;      
         }
     }
     //遍历，得到每个点的和，传给矩阵result
     for( y=0; y<height-5; y++){
         for( x=0; x<width; x++){
             if(y>m_h/2 && y<height - m_h/2 && x>m_w/2 && x<width - m_w/2){
                s.val[0] =  sum[(y - m_h/2) *boxwidth + (x - m_h/2)]/(m_h*m_w);
                 cvSet2D(result,y,x,s);
             }else{
                 s.val[0] = -1;
                cvSet2D(result,y,x,s);
            }//end else
         }//end the first for
     }//end the second for
 }

double getThreshVal_Otsu_8u11(const Mat &_src)
{
	Size size = _src.size();
	const int N = 256;
	int i, j, h[N] = { 0 };

	for (i = 0; i < size.height; i++) {
		const uchar* src = _src.ptr(i);
		j = 0;
		for (; j <= size.width - 4; j += 4) {
			int v0 = src[j], v1 = src[j + 1];
			h[v0]++; h[v1]++;
			v0 = src[j + 2]; v1 = src[j + 3];
			h[v0]++; h[v1]++;
		}
		for (; j < size.width; j++)
			h[src[j]]++;
	}

	double mu = 0, scale = 1. / (size.width*size.height);
	for (i = 0; i < N; i++)
		mu += i*(double)h[i];

	mu *= scale;
	double mu1 = 0, q1 = 0;
	double max_sigma = 0, max_val = 0;

	for (i = 0; i < N; i++) {
		double p_i, q2, mu2, sigma;

		p_i = h[i] * scale;
		mu1 *= q1;
		q1 += p_i;
		q2 = 1. - q1;

		if (std::min(q1, q2) < FLT_EPSILON || std::max(q1, q2) > 1. - FLT_EPSILON)
			continue;

		mu1 = (mu1 + i*p_i) / q1;
		mu2 = (mu - q1*mu1) / q2;
		sigma = q1*q2*(mu1 - mu2)*(mu1 - mu2);
		if (sigma > max_sigma) {
			max_sigma = sigma;
			max_val = i;
		}
	}

	return max_val;
}

void threshold11(InputArray _src, OutputArray _dst, double thresh, double maxval, int type)
{
	Mat src = _src.getMat();
	bool use_otsu = (type & THRESH_OTSU) != 0;
	type &= THRESH_MASK;

	if (use_otsu)
	{
		CV_Assert(src.type() == CV_8UC1);
		thresh = getThreshVal_Otsu_8u11(src);
	}

	_dst.create(src.size(), src.type());
	Mat dst = _dst.getMat();

	if (src.depth() ==3)
	{
		int ithresh = cvFloor(thresh);
		thresh = ithresh;
		int imaxval = cvRound(maxval);
		
		imaxval = saturate_cast<short>(imaxval);

		if (ithresh < SHRT_MIN || ithresh >= SHRT_MAX)
		{
			if (type == THRESH_BINARY || type == THRESH_BINARY_INV ||
				((type == THRESH_TRUNC || type == THRESH_TOZERO_INV) && ithresh < SHRT_MIN) ||
				(type == THRESH_TOZERO&& ithresh >= SHRT_MAX))
			{
				int v = type == THRESH_BINARY ? (ithresh >= SHRT_MAX ? 0 : imaxval) :
					type == THRESH_BINARY_INV ? (ithresh >= SHRT_MAX ? imaxval : 0) :
					/*type == THRESH_TRUNC ?imaxval :*/ 0;
				dst.setTo(v);
			}
			else
				src.copyTo(dst);
		}
		thresh = ithresh;
		maxval = imaxval;
	}
	
	/*parallel_for_(Range(0, dst.rows),
		ThresholdRunner(src, dst, thresh, maxval, type),
		dst.total() / (double)(1 << 16));*/
}

void minMaxLoc11(const Mat& src, double* _minval, double* _maxval, int* _minidx, int* _maxidx) {
	
	int i, j;
	Size size = src.size();
	float maxval = -FLT_MAX;
	const int  *maxidx = 0;
	for (i = 0; i < size.height; i++)//遍历所有像素点
	{
		const float* cov = (const float*)(src.data + src.step*i);
		j = 0;
		for (; j < size.width; j++)
		{
			float v = cov[j];
			if (v > maxval) {
				maxval = v;
			}
		}
		if (_maxval)
			*_maxval = maxval;
	}
	/*
	MatConstIterator it = src.begin;
	Size i, N = src.size(), d = src.hdr ? src.hdr->dims : 0; 
	int type = src.type();
	const int *minidx = 0, *maxidx = 0;
	if (type == CV_32F) {
		float minval = FLT_MAX, maxval = -FLT_MAX; 
		for (i = 0; i < N; i++, ++it) {
			float v = *(const float*)it.ptr; 
			if (v > maxval) {
				maxval = v;
				maxidx = it.node()->idx;
			}
		}  
		if (_maxval)
			*_maxval = maxval;
	}*/
	/*else if(type == CV_64F) {
		double minval = DBL_MAX, maxval = -DBL_MAX; for (i = 0; i < N; i++, ++it) {
			double v = *(const double*)it.ptr; if (v < minval) {
				minval = v;
				minidx = it.node()->idx;
			} if (v > maxval) {
				maxval = v;
				maxidx = it.node()->idx;
			}
		} if (_minval)
			*_minval = minval; if (_maxval)
			*_maxval = maxval;
	}*/
}

void dilation(uchar* data, int width, int height)
{
	int i, j, index, sum, flag;
	sum = height * width * sizeof(uchar);
	uchar *tmpdata = (uchar*)malloc(sum);
	memcpy((char*)tmpdata, (char*)data, sum);
	for (i = 1; i < height - 1; i++)
	{
		for (j = 1; j < width - 1; j++)
		{
			flag = 1;
			for (int m = i - 1; m < i + 2; m++)
			{
				for (int n = j - 1; n < j + 2; n++)
				{
					if (tmpdata[i * width + j] == 255
						|| tmpdata[m * width + n] == 255)
					{
						flag = 0;
						break;
					}
				}
				if (flag == 0)
				{
					break;
				}
			}
			if (flag == 0)
			{
				data[i * width + j] = 255;
			}
			else
			{
				data[i * width + j] = 0;
			}
		}
	}
	free(tmpdata);
}

void goodFeatures(InputArray _image, OutputArray _corners,int maxCorners, double qualityLevel, double minDistance)
{
	//如果需要对_image全图操作，则给_mask传入cv::Mat()，否则传入感兴趣区域
	Mat image = _image.getMat();
	
	Mat eig(image.rows,image.cols, CV_32F);
	Mat	tmp(image.rows, image.cols, CV_32F);   //eig存储每个像素协方差矩阵的最小特征值，tmp用来保存经膨胀后的eig

	//calcMinEigenVal(image, eig);
	cornerMinEigenVal(image, eig); //计算每个像素对应的协方差矩阵的最小特征值，保存在eig中 

	double maxVal = 0;
	minMaxLoc11(eig, 0, &maxVal, 0, 0);   //maxVal保存了eig的最大值
	threshold11(eig, eig, maxVal*qualityLevel, 0, THRESH_TOZERO);  //阈值设置为maxVal乘以qualityLevel，大于此阈值的保持不变，小于此阈值的都设为0

	eig.copyTo(tmp);															 //默认用3*3的核膨胀，膨胀之后，除了局部最大值点和原来相同，其它非局部最大值点被  
	dilation(tmp.data,tmp.size().width,tmp.size().height);

	//3*3邻域内的最大值点取代，如不理解，可看一下灰度图像的膨胀原理  
	//dilate(eig, tmp, Mat());  //tmp中保存了膨胀之后的eig

	Size imgsize = image.size();

	vector<const float*> tmpCorners;  //存放粗选出的角点地址

									  // collect list of pointers to features - put them into temporary image 
	for (int y = 1; y < imgsize.height - 1; y++)
	{
		const float* eig_data = (const float*)eig.ptr(y);  //获得eig第y行的首地址
		const float* tmp_data = (const float*)tmp.ptr(y);  //获得tmp第y行的首地址

		for (int x = 1; x < imgsize.width - 1; x++)
		{
			float val = eig_data[x];
			if (val != 0 && val == tmp_data[x] )  //val == tmp_data[x]说明这是局部极大值
				tmpCorners.push_back(eig_data + x);  //保存其位置
		}
	}

	//-----------此分割线以上是根据特征值粗选出的角点，我们称之为弱角点----------//
	//-----------此分割线以下还要根据minDistance进一步筛选角点，仍然能存活下来的我们称之为强角点----------//
	std::sort(tmpCorners.begin(), tmpCorners.end(), greaterThanPtr()); //按特征值降序排列，注意这一步很重要，后面的很多编程思路都是建立在这个降序排列的基础上
	vector<Point2f> corners;
	size_t i, j, total = tmpCorners.size(), ncorners = 0;

	//下面的程序有点稍微难理解，需要自己仔细想想
	if (minDistance >= 1)
	{
		// Partition the image into larger grids
		int w = image.cols;
		int h = image.rows;

		const int cell_size = cvRound(minDistance);   //向最近的整数取整

													  //这里根据cell_size构建了一个矩形窗口grid(虽然下面的grid定义的是vector<vector>，而并不是我们这里说的矩形窗口，但为了便于理解,还是将grid想象成一个grid_width * grid_height的矩形窗口比较好)，除以cell_size说明grid窗口里相差一个像素相当于_image里相差minDistance个像素，至于为什么加上cell_size - 1后面会讲
		const int grid_width = (w + cell_size - 1) / cell_size;
		const int grid_height = (h + cell_size - 1) / cell_size;

		std::vector<std::vector<Point2f> > grid(grid_width*grid_height);  //vector里面是vector，grid用来保存获得的强角点坐标

		minDistance *= minDistance;  //平方，方面后面计算，省的开根号

		for (i = 0; i < total; i++)     // 刚刚粗选的弱角点，都要到这里来接收新一轮的考验
		{
			int ofs = (int)((const uchar*)tmpCorners[i] - eig.data);  //tmpCorners中保存了角点的地址，eig.data返回eig内存块的首地址
			int y = (int)(ofs / eig.step);   //角点在原图像中的行
			int x = (int)((ofs - y*eig.step) / sizeof(float));  //在原图像中的列

			bool good = true;  //先认为当前角点能接收考验，即能被保留下来

			int x_cell = x / cell_size;  //x_cell，y_cell是角点（y,x）在grid中的对应坐标
			int y_cell = y / cell_size;

			int x1 = x_cell - 1;  // (y_cell，x_cell）的4邻域像素
			int y1 = y_cell - 1;  //现在知道为什么前面grid_width定义时要加上cell_size - 1了吧，这是为了使得（y,x）在grid中的4邻域像素都存在，也就是说(y_cell，x_cell）不会成为边界像素
			int x2 = x_cell + 1;
			int y2 = y_cell + 1;

			// boundary check，再次确认x1,y1,x2或y2不会超出grid边界
			x1 = std::max(0, x1);  //比较0和x1的大小
			y1 = std::max(0, y1);
			x2 = std::min(grid_width - 1, x2);
			y2 = std::min(grid_height - 1, y2);

			//记住grid中相差一个像素，相当于_image中相差了minDistance个像素
			for (int yy = y1; yy <= y2; yy++)  // 行
			{
				for (int xx = x1; xx <= x2; xx++)  //列
				{
					vector <Point2f> &m = grid[yy*grid_width + xx];  //引用

					if (m.size())  //如果(y_cell，x_cell)的4邻域像素，也就是(y,x)的minDistance邻域像素中已有被保留的强角点
					{
						for (j = 0; j < m.size(); j++)   //当前角点周围的强角点都拉出来跟当前角点比一比
						{
							float dx = x - m[j].x;
							float dy = y - m[j].y;
							//注意如果(y,x)的minDistance邻域像素中已有被保留的强角点，则说明该强角点是在(y,x)之前就被测试过的，又因为tmpCorners中已按照特征值降序排列（特征值越大说明角点越好），这说明先测试的一定是更好的角点，也就是已保存的强角点一定好于当前角点，所以这里只要比较距离，如果距离满足条件，可以立马扔掉当前测试的角点
							if (dx*dx + dy*dy < minDistance)
							{
								good = false;
								goto break_out;
							}
						}
					}
				}   // 列
			}    //行

		break_out:

			if (good)
			{
				// printf("%d: %d %d -> %d %d, %d, %d -- %d %d %d %d, %d %d, c=%d\n",
				//    i,x, y, x_cell, y_cell, (int)minDistance, cell_size,x1,y1,x2,y2, grid_width,grid_height,c);
				grid[y_cell*grid_width + x_cell].push_back(Point2f((float)x, (float)y));

				corners.push_back(Point2f((float)x, (float)y));
				++ncorners;

				if (maxCorners > 0 && (int)ncorners == maxCorners)  //由于前面已按降序排列，当ncorners超过maxCorners的时候跳出循环直接忽略tmpCorners中剩下的角点，反正剩下的角点越来越弱
					break;
			}
		}
	}
	else    //除了像素本身，没有哪个邻域像素能与当前像素满足minDistance < 1,因此直接保存粗选的角点
	{
		for (i = 0; i < total; i++)
		{
			int ofs = (int)((const uchar*)tmpCorners[i] - eig.data);
			int y = (int)(ofs / eig.step);   //粗选的角点在原图像中的行
			int x = (int)((ofs - y*eig.step) / sizeof(float));  //在图像中的列

			corners.push_back(Point2f((float)x, (float)y));
			++ncorners;
			if (maxCorners > 0 && (int)ncorners == maxCorners)
				break;
		}
	}

	Mat(corners).convertTo(_corners, _corners.fixedType() ? _corners.type() : CV_32F);
}