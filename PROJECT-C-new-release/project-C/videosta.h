////初始化函数////
///初始化参数，包括参考帧、当前帧及滤波器
void Initializevideosta(int height, int width);

////稳像函数////
///传入当前帧的数据img，图像的宽和高进行稳像，流程为对参考帧和当前帧进行角点检测，之后进行粗匹配，再采用Ransac算法淘汰错误匹配，
///求出H矩阵，从而利用x方向、y方向的一阶滤波器对平移参数dx,dy进行滤波，最后对参考帧进行平移变换完成稳像过程
void videosta(char *img, int height, int width,int *k);

////结束函数////
///释放相关内存
void  sopt_videosta();
