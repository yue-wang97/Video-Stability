//#include "stdafx.h"
#include <opencv2/opencv.hpp>  
#include <opencv2/highgui/highgui.hpp>  
#include <iostream>  
#include <cassert>  
#include <cmath>  
#include"goodfeatures.h"
#include <fstream>  


using namespace std;
using namespace cv;

// This video stablisation smooths the global trajectory using a sliding average window  
//const int SMOOTHING_RADIUS = 15; // In frames. The larger the more stable the video, but less reactive to sudden panning  
const int HORIZONTAL_BORDER_CROP = 20; // In pixels. Crops the border to reduce the black borders from stabilisation being too noticeable.  

									   // 1. Get previous to current frame transformation (dx, dy, da) for all frames  
									   // 2. Accumulate the transformations to get the image trajectory  
									   // 3. Smooth out the trajectory using an averaging window  
									   // 4. Generate new set of previous to current transform, such that the trajectory ends up being the same as the smoothed trajectory  
									   // 5. Apply the new transformation to the video  


void Result1(IplImage *inImage1, IplImage *inImage2,
	IplImage *cornerMap1, IplImage *cornerMap2);
struct Trajectory
{
	Trajectory() {}
	Trajectory(double _x, double _y, double _a) {
		x = _x;
		y = _y;
		a = _a;
	}
	// "+"  
	friend Trajectory operator+(const Trajectory &c1, const Trajectory  &c2) {
		return Trajectory(c1.x + c2.x, c1.y + c2.y, c1.a + c2.a);
	}
	//"-"  
	friend Trajectory operator-(const Trajectory &c1, const Trajectory  &c2) {
		return Trajectory(c1.x - c2.x, c1.y - c2.y, c1.a - c2.a);
	}
	//"*"  
	friend Trajectory operator*(const Trajectory &c1, const Trajectory  &c2) {
		return Trajectory(c1.x*c2.x, c1.y*c2.y, c1.a*c2.a);
	}
	//"/"  
	friend Trajectory operator/(const Trajectory &c1, const Trajectory  &c2) {
		return Trajectory(c1.x / c2.x, c1.y / c2.y, c1.a / c2.a);
	}
	//"="  
	Trajectory operator =(const Trajectory &rx) {
		x = rx.x;
		y = rx.y;
		a = rx.a;
		return Trajectory(x, y, a);
	}

	double x;
	double y;
	double a; // angle  
};
//  
int main(int argc, char **argv)
{

	// For further analysis  
	ofstream out_transform("prev_to_cur_transformation.txt");
	ofstream out_trajectory("trajectory.txt");
	ofstream out_smoothed_trajectory("smoothed_trajectory.txt");
	ofstream out_new_transform("new_prev_to_cur_transformation.txt");

	//	ofstream out_trajectory("trajectory.txt");

	VideoCapture cap("test.mp4");
	//VideoCapture cap("my02.avi");
	//VideoCapture cap("car1.avi");
	assert(cap.isOpened());

	Mat cur, cur_grey;
	Mat prev, prev_grey;

	cap >> prev;//get the first frame.ch  
	prev.copyTo(cur);
	cvtColor(prev, prev_grey, COLOR_BGR2GRAY);

	// Step 1 - Get previous to current frame transformation (dx, dy, da) for all frames  
	// Accumulated frame to frame transform  
	double a = 0;
	double x = 0;
	double y = 0;
	// Step 2 - Accumulate the transformations to get the image trajectory  
	vector <Trajectory> trajectory; // trajectory at all frames  
									//  
									// Step 3 - Smooth out the trajectory using an averaging window  
	vector <Trajectory> smoothed_trajectory; // trajectory at all frames  
	Trajectory X;//posteriori state estimate  
	Trajectory  X_;//priori estimate  
	Trajectory P;// posteriori estimate error covariance  
	Trajectory P_;// priori estimate error covariance  
	Trajectory K;//gain  
	Trajectory  z;//actual measurement  
	double pstd = 4e-3;//can be changed  4e-3
	double cstd = 0.25;//can be changed  0.25
					   /*他们被假设成高斯白噪声(White Gaussian Noise)，他们的covariance 分别是Q，R
					   （这里我们假设他们不随系统状态变化而变化）它们值的选取很重要，
					   一般要通过实验统计得到，它们分布代表了状态空间估计的误差和测量的误差。。
					   */
	Trajectory Q(pstd, pstd, 0);// process noise covariance  
	Trajectory R(cstd, cstd, 0);// measurement noise covariance   
								// Step 4 - Generate new set of previous to current transform, such that the trajectory ends up being the same as the smoothed trajectory   
								// Step 5 - Apply the new transformation to the video  
								//cap.set(CV_CAP_PROP_POS_FRAMES, 0);  
	Mat T(2, 3, CV_64F);
	//vert_border是图像稳像后大约要剪切的边缘大小，其是与HORIZONTAL_BORDER_CROP成比例的，其实可以随意
	int vert_border = HORIZONTAL_BORDER_CROP * prev.rows / prev.cols; // get the aspect ratio correct  
																	  //VideoWriter outputVideo("compare.avi", CV_FOURCC('M', 'P', '4', '2'), 20, cvSize(prev.rows, prev.cols), 1);
	IplImage*srcimg;
	srcimg = &IplImage(prev);
	CvVideoWriter*outputVideo;
	//这里的VideoWriter其实是可以用的，是自己定义图像cvSize(prev.rows, prev.cols)时值要与图像大小一致，否则制作的视频没用
	outputVideo = cvCreateVideoWriter("my.avi", CV_FOURCC('M', 'J', 'P', 'G'), 20, cvGetSize(srcimg));
	int k = 1;
	//获得图像总帧数
	int max_frames = cap.get(CV_CAP_PROP_FRAME_COUNT);
	Mat last_T;
	Mat prev_grey_, cur_grey_;
	int framecoun = 0;
	while (cur.data != NULL) {
		cvtColor(cur, cur_grey, COLOR_BGR2GRAY);
		framecoun++;
		// vector from prev to cur  
		vector <Point2f> prev_corner, cur_corner;
		vector <Point2f> prev_corner2, cur_corner2;
		vector <uchar> status;
		vector <float> err;
		IplImage  inImage1 = (IplImage)prev_grey;
		int height = inImage1.height;
		int width = inImage1.width;
		int Step = inImage1.widthStep;
		IplImage * cornerMap1 = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 1);
		uchar * cornerMap1Data;
		cornerMap1Data = (uchar *)cornerMap1->imageData;
		for (int i = 0; i < height; i++)
			for (int j = 0; j < width; j++) {
				cornerMap1Data[i*Step + j] = 0;
			}
		//要记得参数的意义，200代表检测角点的个数，0.01代表角点的质量，一般0.4，越大质量越好，30角点间的像素距离
	//	goodFeatures(prev_grey.data, cornerMap1, width, height, 200, 0.01, 30);
		goodFeatures(prev_grey, prev_corner, 200, 0.01, 30);
		int p = 0, xx, yy;
		while (p <prev_corner.size()) {
			xx = prev_corner[p].x;
			yy = prev_corner[p++].y;
			cornerMap1Data[yy*Step + xx] = 1;
		}
		Result1(&inImage1, &inImage1, cornerMap1, cornerMap1);
		calcOpticalFlowPyrLK(prev_grey, cur_grey, prev_corner, cur_corner, status, err);

		// weed out bad matches  
		for (size_t i = 0; i < status.size(); i++) {
			if (status[i]) {
				prev_corner2.push_back(prev_corner[i]);
				cur_corner2.push_back(cur_corner[i]);
			}
		}

		// translation + rotation only  其是使用RANsc算法，故有噪声点也能很好的选取代表点
		Mat T = estimateRigidTransform(prev_corner2, cur_corner2, false); // false = rigid transform, no scaling/shearing  
																		  // in rare cases no transform is found. We'll just use the last known good transform.  
																		  /*
																		  这里是为了防止程序崩掉，也可以用来在当检测的点个数达不到estimateRigidTransform要求的3个点时，
																		  进行使用上次的值，不至于崩掉
																		  */
		if (T.data == NULL) {
			last_T.copyTo(T);
		}
		T.copyTo(last_T);
		// decompose T  
		double dx = T.at<double>(0, 2);
		double dy = T.at<double>(1, 2);
		//	double a = atan2(T.at<double>(1, 0), T.at<double>(0, 0));
		//prev_to_cur_transform.push_back(TransformParam(dx, dy, da));  
		out_transform << k << " " << dx << " " << dy << " " << endl;
		// Accumulated frame to frame transform  
		x += dx;//这里是使用主运动来稳像
		y += dy;
		//	a += da;

		out_trajectory << k << " " << x << " " << y << " " << a << endl;

		z = Trajectory(x, y, a);
		//  
		if (k == 1) {
			// intial guesses  
			X = Trajectory(0, 0, 0); //Initial estimate,  set 0  
			P = Trajectory(1, 1, 1); //set error variance,set 1  
		}
		else
		{
			//time update（prediction）
			//记住A矩阵代表着状态矢量与测量矢量之间的关系，一般一个测量矢量就要有两个状态矢量分别为:X，delt_X 
			//一般X与delt_X有关，delt_X与X无关，即构造矩阵时X的那一行有两个1，delt_X的那一行只有一个1.
			X_ = X; //X_(k) = X(k-1);  
			P_ = P + Q; //P_(k) = P(k-1)+Q;  
						// measurement update（correction）  
			K = P_ / (P_ + R); //gain;K(k) = P_(k)/( P_(k)+R );  
			X = X_ + K*(z - X_); //z-X_ is residual,X(k) = X_(k)+K(k)*(z(k)-X_(k));   
			P = (Trajectory(1, 1, 1) - K)*P_; //P(k) = (1-K(k))*P_(k);  
		}
		//smoothed_trajectory.push_back(X);  
		out_smoothed_trajectory << k << " " << X.x << " " << X.y << " " << X.a << endl;
		//-  
		// target - current  diff_x是估计的主运动与现实的主运动之间的差值，这里把其当做偏差值
		double diff_x = X.x - x;//  
		double diff_y = X.y - y;
		//	double diff_a = X.a - a;

		dx = dx + diff_x; //进行真正的运动矫正
		dy = dy + diff_y;
		//	da = da + diff_a;

		out_new_transform << k << " " << dx << " " << dy << " " << endl;
		//  
		/*	T.at<double>(0, 0) = cos(da);
		T.at<double>(0, 1) = -sin(da);
		T.at<double>(1, 0) = sin(da);
		T.at<double>(1, 1) = cos(da);*/

		T.at<double>(0, 2) = dx;
		T.at<double>(1, 2) = dy;
		Mat cur2;

		warpAffine(prev, cur2, T, cur.size());
		//cur2 = cur2(Range(vert_border, cur2.rows - vert_border), Range(HORIZONTAL_BORDER_CROP, cur2.cols - HORIZONTAL_BORDER_CROP));

		// Resize cur2 back to cur size, for better side by side comparison  
		//resize(cur2, cur2, cur.size());
		srcimg = &IplImage(cur2);
		cvWriteFrame(outputVideo, srcimg);
		// Now draw the original and stablised side by side for coolness  
		/*		//把图像现实在同一个画布上
		Mat canvas = Mat::zeros(cur.rows, cur.cols * 2 + 10, cur.type());
		prev.copyTo(canvas(Range::all(), Range(0, cur2.cols)));
		cur2.copyTo(canvas(Range::all(), Range(cur2.cols + 10, cur2.cols * 2 + 10)));

		// If too big to fit on the screen, then scale it down by 2, hopefully it'll fit :)
		if (canvas.cols > 1920) {
		resize(canvas, canvas, Size(canvas.cols / 2, canvas.rows / 2));
		}
		imshow("before and after", canvas);
		waitKey(10);*/
		//  
		prev = cur.clone();//cur.copyTo(prev);  
		cur_grey.copyTo(prev_grey);
		cout << "Frame: " << k << "/" << max_frames << " - good optical flow: " << prev_corner2.size() << endl;
		k++;
		cap >> cur;
	}
	cvReleaseVideoWriter(&outputVideo);
	return 1;
}
void MarkCornerPoints(IplImage *image, IplImage *cornerMap) {
	int i, j;
	uchar *cornerMapData = 0;
	int height = cornerMap->height;
	int width = cornerMap->width;
	int mapStep = cornerMap->widthStep;
	cornerMapData = (uchar *)cornerMap->imageData;
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			if (cornerMapData[i*mapStep + j] == true) {
				cvCircle(image, cvPoint(j, i), 2, cvScalar(0, 255, 0), 2);
			}
		}
	}
}
void CombineTwoImages(IplImage *image1, IplImage *image2,
	IplImage *outImage)
{
	int i, j, k;
	uchar *outImageData = 0, *image1Data = 0, *image2Data = 0;
	int height = image1->height;
	int width = image1->width;
	int step = image1->widthStep;
	int channels = image1->nChannels;
	int outWidth = outImage->width;
	int outHeight = outImage->height;
	int outStep = outImage->widthStep;
	if (outWidth == width * 2 && outHeight == height) {
	}
	else if (outWidth == width && outHeight == height * 2) {
	}
	else {
		printf("image combining error\n");
		exit(0);
	}
	outImageData = (uchar *)outImage->imageData;
	image1Data = (uchar *)image1->imageData;
	image2Data = (uchar *)image2->imageData;
	for (i = 0; i < outHeight; i++) {
		for (j = 0; j < outWidth; j++) {

			for (k = 0; k < channels; k++) {
				if (i < height && j < width) {
					outImageData[i*outStep + j*channels + k]
						= image1Data[i*step + j*channels + k];
				}
				else if ((i >= height && j < width)) {
					outImageData[i*outStep + j*channels + k]
						= image2Data[(i - height)*step + j*channels + k];
				}
				else if ((i < height && j >= width)) {
					outImageData[i*outStep + j*channels + k]
						= image2Data[i*step + (j - width)*channels + k];
				}
				else {
					printf("there is no i > height & j > width \n");
					exit(0);
				}
			}
		}
	}
}
void WriteImage(IplImage *image, char *imageName) {
	if (!cvSaveImage(imageName, image)) {
		printf("Could not save: %s\n", imageName);
	}
}
//
void Result1(IplImage *inImage1, IplImage *inImage2,
	IplImage *cornerMap1, IplImage *cornerMap2)
{
	IplImage *outImage1 = 0, *outImage2 = 0, *outImage3 = 0;
	int height, width, channels;
	// create the output images

	height = inImage1->height;
	width = inImage1->width;
	channels = inImage1->nChannels;
	outImage1 = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, channels);
	outImage2 = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, channels);
	// draw a circle on the corner point
	cvCopy(inImage1, outImage1);
	cvCopy(inImage2, outImage2);
	MarkCornerPoints(outImage1, cornerMap1);
	MarkCornerPoints(outImage2, cornerMap2);
	// create the output images : 2 images in the same output image
	outImage3 = cvCreateImage(cvSize(width * 2, height), IPL_DEPTH_8U, channels);
	CombineTwoImages(outImage1, outImage2, outImage3);
	// display the result image
	cvNamedWindow("output image", CV_WINDOW_AUTOSIZE);
	cvShowImage("output image", outImage3);
	cvWaitKey(0);
	cvDestroyWindow("output image");
	// write output image
	WriteImage(outImage1, "result_corner1.jpg");
	WriteImage(outImage2, "result_corner2.jpg");
	WriteImage(outImage3, "result_step1.jpg");
	CombineTwoImages(inImage1, inImage2, outImage3);
	WriteImage(outImage3, "inputs.jpg");
	cvReleaseImage(&outImage1);
	cvReleaseImage(&outImage2);
	cvReleaseImage(&outImage3);
}