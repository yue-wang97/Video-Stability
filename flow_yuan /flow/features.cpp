#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "features.h"

using namespace std;
using namespace cv;


void cornerMinEigenVal(uchar* src, CvMat * eigenv,int width,int height)
{
	double scale = (double)6;
	scale *= 255.;
	scale = 1. / scale;
//	Mat Dx, Dy;   //保存每个像素点的水平方向和垂直方向的一阶差分  
				  //	Sobel(src, Dx, CV_32F, 1, 0, 3, scale, 0, BORDER_DEFAULT);
				  //	Sobel(src, Dy, CV_32F, 0, 1, 3, scale, 0, BORDER_DEFAULT);
	CvMat * Dx = cvCreateMat(height, width, CV_32FC1);
	CvMat * Dy = cvCreateMat(height, width, CV_32FC1);
	Sobel66(src, Dx, width,height, CV_32FC1, 1, 0, 3, scale, 0);
	Sobel66(src, Dy, width, height, CV_32FC1, 0, 1, 3, scale, 0);
	//getSobelKernels(src, Dx, 1, 0, 0,CV_32F);
	//getSobelKernels(src, Dy, 0, 1, 0, CV_32F);
//	Size size = height*width;
//	Mat cov(size, CV_32FC3);
	CvMat * cov1 = cvCreateMat(height, width, CV_32FC1);
	CvMat * cov2 = cvCreateMat(height, width, CV_32FC1);
	CvMat * cov3 = cvCreateMat(height, width, CV_32FC1);
	int i, j;
	for (i = 0; i < height; i++)
	{
		float* cov_data1 = (float*)(cov1->data.ptr + i*cov1->step);
		float* cov_data2 = (float*)(cov2->data.ptr + i*cov2->step);
		float* cov_data3 = (float*)(cov3->data.ptr + i*cov3->step);
		float* dxdata = (float*)(Dx->data.ptr + i*Dx->step);
		float* dydata = (float*)(Dy->data.ptr + i*Dy->step);
		for (j = 0; j <width; j++)
		{
			float dx = dxdata[j];
			float dy = dydata[j];
			cov_data1[j ] = dx*dx;  //第一个通道存dx*dx,即M矩阵左上角的元素  
			cov_data2[j ] = dx*dy;//第二个通道存dx*dy,即M矩阵左下角和右上角的元素  
			cov_data3[j ] = dy*dy;//第三个通道存dy*dy,即M矩阵右下角的元素  
		}
	}
//	boxFilter(cov, cov, cov.depth(), Size(3, 3), //计算邻域上的差分相关矩阵（block_size×block_size）  
//		Point(-1, -1), false, BORDER_DEFAULT);
	//IplImage cov_img = IplImage(cov);
	box_filter(cov1, cov1);
	box_filter(cov2, cov2);
	box_filter(cov3, cov3);
	//calcHarris(cov, eigenv, 0.);

	calcMinEigenVal(cov1,cov2,cov3,eigenv);
}
void calcMinEigenVal(CvMat * _cov1, CvMat * _cov2, CvMat * _cov3,CvMat * _dst)
{
	int i, j;
	int width = _cov1->cols;
	int height = _cov1->rows;
//	if (_cov.isContinuous() && _dst.isContinuous())
//	{
		width *= height;
		height = 1;
//	}

	for (i = 0; i < height; i++)//遍历所有像素点  
	{
		float* cov1 = (float*)(_cov1->data.ptr + _cov1->step*i);
		float* cov2 = (float*)(_cov2->data.ptr + _cov2->step*i);
		float* cov3 = (float*)(_cov3->data.ptr + _cov3->step*i);
		float* dst = (float*)(_dst->data.ptr + _dst->step*i);
		j = 0;
		for (; j < width; j++)
		{
			float a = cov1[3] * 0.5f;//cov[j*3]保存矩阵M左上角元素  
			float b = cov2[j];   //cov[j*3+1]保存左下角和右上角元素  
			float c = cov3[j] * 0.5f;//cov[j*3+2]右下角元素  
			dst[j] = (float)((a + c) - std::sqrt((a - c)*(a - c) + b*b));//求最小特征值，一元二次方程求根公式  
		}
	}
}



void getSobelKernels66(CvMat * _kx, CvMat* _ky, int dx, int dy, int _ksize, bool normalize, int ktype)
{
	int i, j, ksizeX = _ksize, ksizeY = _ksize;

	//CV_Assert(ktype == CV_32F || ktype == CV_64F);

//	_kx.create(ksizeX, 1, ktype, -1, true);
//	_ky.create(ksizeY, 1, ktype, -1, true);
//	Mat kx = _kx.getMat();
//	Mat ky = _ky.getMat();
//	CvMat * kx = cvCreateMat(ksizeX, 1, ktype);
//	CvMat * ky = cvCreateMat(ksizeX, 1, ktype);
//	kx = cvCloneMat(_kx);
//	ky = cvCloneMat(_ky);
//	if (_ksize % 2 == 0 || _ksize > 31)
//		CV_Error(CV_StsOutOfRange, "The kernel size must be odd and not larger than 31");
//	vector<int> kerI(std::max(ksizeX, ksizeY) + 1);
//	int * kerI = (int *)malloc((max(ksizeX, ksizeY) + 1) * sizeof(int));
//	CV_Assert(dx >= 0 && dy >= 0 && dx + dy > 0);

	for (int k = 0; k < 2; k++)
	{
		CvMat* kernel = k == 0 ? _kx : _ky;
		int order = k == 0 ? dx : dy;
		int ksize = k == 0 ? ksizeX : ksizeY;

	//	CV_Assert(ksize > order);

		 if (ksize == 3)
		{
			if (order == 0)
			//	kerI[0] = 1, kerI[1] = 2, kerI[2] = 1;
			{
				cvmSet(kernel, 0, 0, 1);
				cvmSet(kernel, 1, 0, 2);
				cvmSet(kernel, 2, 0, 1);

			}
			else if (order == 1)
			//	kerI[0] = -1, kerI[1] = 0, kerI[2] = 1;
			{
				cvmSet(kernel, 0, 0, -1);
				cvmSet(kernel, 1, 0, 0);
				cvmSet(kernel, 2, 0, 1);

			}
			else
			//	kerI[0] = 1, kerI[1] = -2, kerI[2] = 1;
			{
				cvmSet(kernel, 0, 0, 1);
				cvmSet(kernel, 1, 0, -2);
				cvmSet(kernel, 2, 0, 1);

			}
		}
		

	//	Mat temp(kernel->rows, kernel->cols, CV_32S, &kerI[0]);
	//	double scale = !normalize ? 1. : 1. / (1 << (ksize - order - 1));
	//	temp.convertTo(*kernel, ktype, scale);
	}
}

void conv2(uchar* x, CvMat* y, CvMat* z2)
{
	int i, j;
	int n, m;
	int N1 = z2->rows;
	int M1 = z2->cols;
	int N2 = y->rows;
	int M2 = y->cols;
	int step = z2->step;
	CvMat *z=cvCreateMat(N1 + N2 - 1,M1 + M2 - 1,y->type);
	for (i = 0; i<N1 + N2 - 1; i++)
		for (j = 0; j<M1 + M2 - 1; j++)
		{
			float temp = 0;
			for (m = 0; m<N1; m++)
				for (n = 0; n<M1; n++)
					if ((i - m) >= 0 && (i - m)<N2 && (j - n) >= 0 && (j - n)<M2)
						temp += x[m*step+n] * cvmGet(y,i - m,j - n);
			cvmSet(z,i,j,temp);
		}
	for (i = 0; i<N1; i++)
		for (j = 0; j<M1; j++)
		{
			cvmSet(z2,i,j,cvmGet(z,i + (N2 - 1) / 2,j + (M2 - 1) / 2));
		}
}
void pointmull(CvMat* a, float b)
{
	int row = a->rows;
	int col = a->cols;
	int i, j;
	for(i=0;i<row;i++)
		for(j=0;j<col;j++)
			cvmSet(a, i, j, cvmGet(a, i, j)*b);
}

void Sobel66(uchar* _src, CvMat *  _dst, int width,int height,int ddepth, int dx, int dy, int ksize, double scale, double delta, int borderType)
{	
	int ktype = CV_32FC1;
	CvMat * kx = cvCreateMat(3, 1, ktype);
	CvMat * ky = cvCreateMat(3, 1, ktype);
	getSobelKernels66(kx, ky, dx, dy, ksize, false, ktype);
	if (scale != 1)
	{
		// usually the smoothing part is the slowest to compute,
		// so try to scale it instead of the faster differenciating part
		if (dx == 0)
			//kx *= scale;
			pointmull(kx, (float)scale);
		else
			//ky *= scale;
			pointmull(ky, (float)scale);
	}
//	sepFilter2D(src, dst, ddepth, kx, ky, Point(-1, -1), delta, borderType);
	conv2(_src, kx,_dst);
	conv2(_dst->data.ptr, ky, _dst);
}

void box_filter(CvMat* img, CvMat* result) {
	//init part
	double s;
	int width = img->cols, height = img->rows;
	int m_w = 5, m_h = 5;//window_size
	int boxwidth = width - m_w, boxheight = height - m_h;
	int *sum = (int*)malloc(boxwidth *boxheight * sizeof(double));
	int *buff = (int*)malloc(width * sizeof(double));
	memset(sum, 0, boxwidth *boxheight * sizeof(int));
	memset(buff, 0, width * sizeof(int));
	//set buff:from 0 to 4 rows,per col
	int x, y, j;
	for (y = 0; y<m_h; y++) {
		for (x = 0; x<width; x++) {
			uchar pixel = (uchar)cvmGet(img, y,x);
			buff[x] += pixel;
			printf("%d:%d\n", x, buff[x]);
		}
	}
	for (y = 0; y<height - m_h; y++) {
		int Xsum = 0;
		for (j = 0; j<m_w; j++) {
			Xsum += buff[j];//sum of pixel from (0,0) to (m_h,m_w) (also x = 0)
		}

		for (x = 0; x<boxwidth; x++) {
			if (x != 0) {
				Xsum = Xsum - buff[x - 1] + buff[m_w - 1 + x];//Xsum:sum of cols range from x to x+m_w ,rows range from 0 to 4
			}
			sum[y*boxwidth + x] = (float)Xsum;
		}

		for (x = 0; x<width; x++) { 
			uchar pixel = (uchar)cvmGet(img, y, x);//img[y *width + x];    
			uchar pixel2 = (uchar)cvmGet(img, y + m_h, x);//img[(y+mheight) *width + x];  
			buff[x] = buff[x] - pixel + pixel2;
		}
	}
	cout << buff << endl;
	//�������õ�ÿ����ĺͣ���������result
	for (y = 0; y<height - 5; y++) {
		for (x = 0; x<width; x++) {
			if (y>m_h / 2 && y<height - m_h / 2 && x>m_w / 2 && x<width - m_w / 2) {
				s = sum[(y - m_h / 2) *boxwidth + (x - m_h / 2)] / (m_h*m_w);
				cvmSet(result, y, x, s);
			}
			else {
				s = -1;
				cvmSet(result, y, x, s);
			}//end else
		}//end the first for
	}//end the second for
}

void threshold11(CvMat* _src, CvMat* _dst, float thresh, float maxval, int type)
{
	int rows = _src->rows;
	int cols = _src->cols;
	int step = _src->step;
	int i, j;
	float tt;
	for (i = 0; i<rows; i++)
		for (j = 0; j < cols; j++)
		{
			tt = cvmGet(_src,i, j);
			if (tt < thresh)
				cvmSet(_dst,i, j,0);
			else
				cvmSet(_dst, i, j, tt);
		}
}


void minMaxLoc11(CvMat *src, float* _minval, float* _maxval, int* _minidx, int* _maxidx) {

	int i, j;
	int height = src->rows;
	int width = src->cols;
	float maxval = -FLT_MAX;
	int  *maxidx = 0;
	for (i = 0; i < height; i++)//遍历所有像素点
	{
		float* cov = (float*)(src->data.ptr + src->step*i);
		j = 0;
		for (; j < width; j++)
		{
			float v = cov[j];
			if (v > maxval) {
				maxval = v;
			}
		}
		if (_maxval)
			*_maxval = maxval;
	}
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


void goodFeatures(uchar* _image, IplImage* cornerMap1,int width,int height,int maxCorners, float qualityLevel, float minDistance)
{
	//如果需要对_image全图操作，则给_mask传入cv::Mat()，否则传入感兴趣区域
//	Mat image = _image.getMat();

//	Mat eig(image.rows, image.cols, CV_32F);
//	Mat	tmp(image.rows, image.cols, CV_32F);  //eig存储每个像素协方差矩阵的最小特征值，tmp用来保存经膨胀后的eig
	CvMat * eig = cvCreateMat(height, width, CV_32FC1);
	CvMat * tmp = cvCreateMat(height, width, CV_32FC1);
											   //calcMinEigenVal(image, eig);
	cornerMinEigenVal(_image, eig,width,height); //计算每个像素对应的协方差矩阵的最小特征值，保存在eig中 

	float maxVal = 0;
	minMaxLoc11(eig, 0, &maxVal, 0, 0);   //maxVal保存了eig的最大值
	threshold11(eig, eig, maxVal*qualityLevel, 0, THRESH_TOZERO);  //阈值设置为maxVal乘以qualityLevel，大于此阈值的保持不变，小于此阈值的都设为0

	tmp= cvCloneMat(eig);															 //默认用3*3的核膨胀，膨胀之后，除了局部最大值点和原来相同，其它非局部最大值点被  
	dilation(tmp->data.ptr,width,height);

	//3*3邻域内的最大值点取代，如不理解，可看一下灰度图像的膨胀原理  
	//dilate(eig, tmp, Mat());  //tmp中保存了膨胀之后的eig

//	Size imgsize = image.size();

//	vector<const float*> tmpCorners;  //存放粗选出的角点地址
	float tmpCorners[500];
	Point2f Corners[500];
	IplImage *img = cvCreateImage(CvSize(width, height), IPL_DEPTH_8U, 1);
	for (int i = 0; i < height; i++)
		for (int j = 0; j < width; j++)
			img->imageData[i*img->widthStep + j] = 0;
	int num = 0;
									  // collect list of pointers to features - put them into temporary image 
	for (int y = 1; y < height - 1; y++)
	{
		float* eig_data = (float*)(eig->data.ptr+eig->step*y);  //获得eig第y行的首地址
		float* tmp_data = (float*)(tmp->data.ptr + tmp->step*y);  //获得tmp第y行的首地址

		for (int x = 1; x < width - 1; x++)
		{
			float val = eig_data[x];
			if (val != 0 && val == tmp_data[x])//val == tmp_data[x]说明这是局部极大值
			{
				tmpCorners[num] = val;
				Corners[num].x = x;
				Corners[num].y = y;
				img->imageData[(y-1)*img->widthStep + x-1] = 1;
				num++;
				//	tmpCorners.push_back(eig_data + x);  //保存其位置
			}
		}
	}

	//-----------此分割线以上是根据特征值粗选出的角点，我们称之为弱角点----------//
	//-----------此分割线以下还要根据minDistance进一步筛选角点，仍然能存活下来的我们称之为强角点----------//
//	std::sort(tmpCorners.begin(), tmpCorners.end(), greaterThanPtr()); //按特征值降序排列，注意这一步很重要，后面的很多编程思路都是建立在这个降序排列的基础上
	quick_sort(tmpCorners, 0, num - 1, Corners);
	//Point2f corners[200];
//	size_t i, j, total = tmpCorners.size(), ncorners = 0;
	int i, j, total = num, ncorners = 0;
	//下面的程序有点稍微难理解，需要自己仔细想想
	if (minDistance >= 1)
	{
		// Partition the image into larger grids
		int w = width;
		int h = height;
		int cell_size = 2*cvRound(minDistance);   //向最近的整数取整

													  //这里根据cell_size构建了一个矩形窗口grid(虽然下面的grid定义的是vector<vector>，而并不是我们这里说的矩形窗口，但为了便于理解,还是将grid想象成一个grid_width * grid_height的矩形窗口比较好)，除以cell_size说明grid窗口里相差一个像素相当于_image里相差minDistance个像素，至于为什么加上cell_size - 1后面会讲
//	    int grid_width = (w + cell_size - 1) / cell_size;
//		int grid_height = (h + cell_size - 1) / cell_size;

		//std::vector<std::vector<Point2f> > grid(grid_width*grid_height);  //vector里面是vector，grid用来保存获得的强角点坐标
		IplImage *grid = cvCreateImage(cvSize(cell_size, cell_size), IPL_DEPTH_8U, 1);
		//InitializeImage(domainBlock);
		
		minDistance *= minDistance;  //平方，方面后面计算，省的开根号

		for (i = 0; i < total; i++)     // 刚刚粗选的弱角点，都要到这里来接收新一轮的考验
		{
	//		int ofs = (int)((const uchar*)tmpCorners[i] - eig->data.ptr);  //tmpCorners中保存了角点的地址，eig.data返回eig内存块的首地址
	//		int y = (int)(ofs / eig.step);   //角点在原图像中的行
	//		int x = (int)((ofs - y*eig.step) / sizeof(float));  //在原图像中的列
			int y = Corners[i].y;
			int x = Corners[i].x;

			bool good = true;  //先认为当前角点能接收考验，即能被保留下来
			MakeImageBlock(grid, img, y, x);
	/*		int x_cell = x / cell_size;  //x_cell，y_cell是角点（y,x）在grid中的对应坐标
			int y_cell = y / cell_size;

			int x1 = x_cell - 1;  // (y_cell，x_cell）的4邻域像素
			int y1 = y_cell - 1;  //现在知道为什么前面grid_width定义时要加上cell_size - 1了吧，这是为了使得（y,x）在grid中的4邻域像素都存在，也就是说(y_cell，x_cell）不会成为边界像素
			int x2 = x_cell + 1;
			int y2 = y_cell + 1;

			// boundary check，再次确认x1,y1,x2或y2不会超出grid边界
			x1 = std::max(0, y);  //比较0和x1的大小
			y1 = std::max(0, y1);
			x2 = std::min(grid_width - 1, x2);
			y2 = std::min(grid_height - 1, y2);*/
			int iBegin= std::max(y - cell_size/2, 0);
			int jBegin = std::max(x - cell_size / 2, 0);
			int iEnd = std::min(y + cell_size / 2 + 1, height - 1);
			int jEnd = std::min(x + cell_size / 2 + 1, width - 1);
			for (i = iBegin; i < iEnd; i++) {
				for (j = jBegin; j < jEnd; j++) {
					// calculate NCC at only range image corner points
					if (img->imageData[i*img->widthStep + j]) {
						float dx = x - j;
						float dy = y - i;
						if (dx*dx + dy*dy < minDistance)
						{
							good = false;
							goto break_out;
						}
					}
				}
			}
		break_out:
			//记住grid中相差一个像素，相当于_image中相差了minDistance个像素
	/*		for (int yy = y1; yy <= y2; yy++)  // 行
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

		break_out:*/

			if (good)
			{
				// printf("%d: %d %d -> %d %d, %d, %d -- %d %d %d %d, %d %d, c=%d\n",
				//    i,x, y, x_cell, y_cell, (int)minDistance, cell_size,x1,y1,x2,y2, grid_width,grid_height,c);
			//	grid[y_cell*grid_width + x_cell].push_back(Point2f((float)x, (float)y));

			//	corners.push_back(Point2f((float)x, (float)y));
			//	++ncorners;
			//	corners[ncorners].x = (float)x;
			//	corners[ncorners++].y = (float)y;
				ncorners++;
				cornerMap1->imageData[y*cornerMap1->widthStep + x] = 1;
				if (maxCorners > 0 && (int)ncorners == maxCorners)  //由于前面已按降序排列，当ncorners超过maxCorners的时候跳出循环直接忽略tmpCorners中剩下的角点，反正剩下的角点越来越弱
					break;
			}
		}
	}
	else    //除了像素本身，没有哪个邻域像素能与当前像素满足minDistance < 1,因此直接保存粗选的角点
	{
		for (i = 0; i < total; i++)
		{
	//		int ofs = (int)((const uchar*)tmpCorners[i] - eig.data);
	//		int y = (int)(ofs / eig.step);   //粗选的角点在原图像中的行
	//		int x = (int)((ofs - y*eig.step) / sizeof(float));  //在图像中的列
			int y = Corners[i].y;
			int x = Corners[i].x;
			ncorners++;
		//	corners[ncorners].x = (float)x;
		//	corners[ncorners++].y = (float)y;
			cornerMap1->imageData[y*cornerMap1->widthStep + x] = 1;
			if (maxCorners > 0 && (int)ncorners == maxCorners)
				break;
		}
	}

//	Mat(corners).convertTo(_corners, _corners.fixedType() ? _corners.type() : CV_32F);
}


void quick_sort(float *s, int l, int r,Point2f * m)
{
	int i, j;
	float x,x1,y1;
	if (l < r)
	{
		i = l;
		j = r;
		x = s[i];
		x1 = m[i].x;
		y1 = m[i].y;
		while (i < j)
		{
			while (i < j && s[j] > x)
				j--; /* 从右向左找第一个小于x的数 */
			if (i < j)
			{
				s[i] = s[j];
				m[i].x = m[j].x;
				m[i].y = m[j].y;
				i++;
			}
			while (i < j && s[i] < x)
				i++; /* 从左向右找第一个大于x的数 */
			if (i < j)
			{
				s[j] = s[i];
				m[j].x = m[i].x;
				m[j].y = m[i].y;
				j--;
			}
		}
		s[i] = x;
		m[i].x = x1;
		m[i].y = y1;
		quick_sort(s, l, i - 1,m); /* 递归调用 */
		quick_sort(s, i + 1, r,m);
	}
}


void MakeImageBlock(IplImage *block, IplImage *image,
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
	for (i = 0; i < blockHeight; i++) {
		for (j = 0; j < blockWidth; j++) {
			for (k = 0; k < channels; k++) {
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





