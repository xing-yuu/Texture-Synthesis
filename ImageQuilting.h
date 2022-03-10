#pragma once
#include<opencv2\opencv.hpp>   
#include<opencv2\highgui\highgui.hpp>
#include<iostream>
#include<cstring>
#include<queue>
#include<vector>
#include"time.h"
using namespace std;
using namespace cv;
class ImageQuilting
{
	private:
		Mat *input_images;//输入的原始图像
		Mat *quilting_results;//生成的图像
		Mat *tem_rol;//求宽度的时候用的距离矩阵,横向的一行
		Mat *tem_coll;//求宽度的时候用的距离矩阵,纵向的左侧
		Mat *tem_colr;//求宽度的时候用的距离矩阵,纵向的右侧
		int OverlapX;//重叠部分竖条的宽度
		int OverlapY;//重叠部分横条的高度
		int block_x;//生成图像的长度
		int block_y;//生成图像的高度
		int select_x;//选取块的长度
		int select_y;//选取块的宽度
		int num_x;//最后结果横向生成多少块
		int num_y;//最后结果纵向生成多少块
		int Bestx, Besty;//最后选定的地方
		float e1[100][100];//路径规划用
		float e2[100][100];//路径规划用
		float e3[100][100];//路径规划用
	public:
		ImageQuilting(string path, int the_blockx, int the_blocky, int tnum_x, int tnum_y);//每一块的大小，横
		void ImageGenerate();
		void Getdistance(Mat* coll, int X, int Y);//第一行，只需要求左侧即可
		void Getdistance(Mat* rol, int X, int Y, bool flag);//第一列，只需要求上侧即可
		void Getdistance(Mat* col, Mat* rol, int X, int Y);//正常上侧加左侧
		void Getdistance(Mat* coll, Mat* colr, int X, int Y, bool flag);//第一行最后一个，左侧加右侧 
		void Getdistance(Mat* coll, Mat* rol, Mat* colr, int X, int Y);//正常行的最后一个，左侧加右侧加上侧
};

