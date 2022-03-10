#include<opencv2\opencv.hpp>   
#include<opencv2\highgui\highgui.hpp>
#include"ImageQuilting.h" 
#include"time.h"
using namespace std;
using namespace cv;

int main()
{
	clock_t start, finish;
	start = clock();
	ImageQuilting thepic("./TestImg/5.jpg",20 ,20,20,20);
	thepic.ImageGenerate();
	finish = clock();
	double Total_time = (double)(finish - start);
	printf("\n函数运行时间：%0.3f毫秒 \n", Total_time);

}