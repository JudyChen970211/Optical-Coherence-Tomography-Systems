// matlab_transform.cpp : 定義主控台應用程式的進入點。
//

#include	"stdafx.h"
#include	<iostream>
#include	<stdio.h>
#include	<stdlib.h>
#include	<fstream>
#include    <vector>
#include	<math.h>
#include    <complex>
#include	<opencv2/opencv.hpp>  
#include	<opencv2/imgproc/imgproc.hpp>
#include	<cv.hpp>

#define MAX 2048
#define M_PI 3.1415926535897932384

using namespace cv;
using namespace std;
enum
{
	NO_AScan = 1000,
	NO_perAScan = 2048,
	No_ReSample = 2048,
	No_PostFFT = 1024,
	Letter_size = 12,
	Num = 1024000
};

// %------------------------------Variable--------------------------------------------

#define data_name = "1";
#define output_name = "outbmp";
#define loop = "data from MEDT";

//----------------------------------FFT-----------------------------------------------

int log2(int N)    /*function to calculate the log2(.) of int numbers*/
{
	int k = N, i = 0;
	while (k) {
		k >>= 1;
		i++;
	}
	return i - 1;
}

int check(int n)    //checking if the number of element is a power of 2
{
	return n > 0 && (n & (n - 1)) == 0;
}

int reverse(int N, int n)    //calculating revers number
{
	int j, p = 0;
	for (j = 1; j <= log2(N); j++) {
		if (n & (1 << (log2(N) - j)))
			p |= 1 << (j - 1);
	}
	return p;
}

void ordina(complex<float>* f1, int N) //using the reverse order in the array
{
	complex<float> f2[MAX];
	for (int i = 0; i < N; i++)
		f2[i] = f1[reverse(N, i)];
	for (int j = 0; j < N; j++)
		f1[j] = f2[j];
}

void transform(complex<float>* f, int N) //
{
	ordina(f, N);    //first: reverse order
	complex<float> *W;
	W = (complex<float> *)malloc(N / 2 * sizeof(complex<float>));
	W[1] = polar(1., -2. * M_PI / N);
	W[0] = 1;
	for (int i = 2; i < N / 2; i++)
		W[i] = pow(W[1], i);
	int n = 1;
	int a = N / 2;
	for (int j = 0; j < log2(N); j++) {
		for (int i = 0; i < N; i++) {
			if (!(i & n)) {
				complex<float> temp = f[i];
				complex<float> Temp = W[(i * a) % (n * a)] * f[i + n];
				f[i] = temp + Temp;
				f[i + n] = temp - Temp;
			}
		}
		n *= 2;
		a = a / 2;
	}
}

void FFT(complex<float>* f, int N, double d)
{
	transform(f, N);
	for (int i = 0; i < N; i++)
		f[i] *= d; //multiplying by step
}

void Forward_FFT( int now, vector<float> do_all, float D1_raw[2048], float table[2048][3], complex<float>  D2_Resamp[2048] ) {

	for (int j = now*NO_perAScan, k = 0; j < now*NO_perAScan + NO_perAScan && j < do_all.size(); j++, k++)
		D1_raw[k] = do_all.at(j);

	for (int j = 0; j < No_ReSample; j++) {
		int num = table[j][2] - 1;
		float temp1 = D1_raw[num] * table[j][0];
		num++;
		float temp2 = D1_raw[num] * table[j][1];
		D2_Resamp[j] = temp1 + temp2;
	} // for




} // void 

void Back_FFT( int now, float D4_absFFT[2048], complex<float>  D2_Resamp[2048], vector<float> & Img ) {

	for (int j = 0; j < MAX; j++)
		D4_absFFT[j] = abs(D2_Resamp[j]);

	D4_absFFT[0] = 1;
	D4_absFFT[1] = 1;
	float max_D4 = 0.0;
	// max(D4_absFFT(1, 1:No_PostFFT));
	for (int k = 0; k < 2047; k++)
		if (D4_absFFT[k] > max_D4) max_D4 = D4_absFFT[k];
	
	int nowBegin = No_PostFFT * now; // No_PostFFT*(ind_AScan - 1)
	int nowEnd = No_PostFFT*now + No_PostFFT; // No_PostFFT*(ind_AScan - 1) + No_PostFFT 
											//  Img(1, No_PostFFT*(ind_AScan - 1) + 1:No_PostFFT*(ind_AScan - 1) + No_PostFFT) = D4_absFFT(1, 1:No_PostFFT) / max(D4_absFFT(1, 1:No_PostFFT));

	for (int k = nowBegin, s = 0; k < nowEnd; k++, s++) {
		Img[k] = D4_absFFT[s] / max_D4;
		Img[k] = ((20 * log10(Img[k])) / 40) * 255 + 255;
		if (Img[k] == -INFINITY) Img[k] = 255.0;
	} // for
	cout << nowEnd << " = END " << endl;

} // 

//--------------------------------------------------------------------------------------------

int main()
{
	// %-------------- - Interpolation 整理系數 整數小數分離-------------------------------- -
	float table[2048][3];
	float a = 0.0;
	vector<float> do_all;
	vector <float> Img;
	// 開檔 // 

	fstream fin;
	fin.open("table2048.txt", ios::in);
	if (!fin) {
		cout << "無法讀入檔案\n";
		system("pause");
		return 0;
	} // 

	for (int i = 0; i < 2048; i++) {
		fin >> table[i][0];
		// table ++ // 
		table[i][0]++;
		// table ++ // 
		int interger;
		interger = (int)table[i][0] + 1;
		table[i][2] = table[i][0];         // 整數 // 
		table[i][0] = interger - table[i][0];
		table[i][1] = 1 - table[i][0]; // table[i][1] 為小數部分
		table[i][2] = table[i][2] - table[i][1];  // table[i][2] 為整數部分
	} // for
	fin.close();
	// 開檔 //

	fin.open("D0_all.txt", ios::in);
	for (int i = 0; fin >> a; i++)
		do_all.push_back(a);
	// D0_all = abs(D0_all);  //abs(x)：純量的絕對值或向量的長度 ????


	float D1_raw[2048] = { 0.0 };
	complex<float> D2_Resamp[4][2048];
	float D3_Hammin[2048] = { 0.0 };
	float D4_absFFT[2048] = { 0.0 };
	float D5_log[1024] = { 0.0 };
	for (int i = 0; i < Num; i++)
		Img.push_back(0.0);
	// ind_AScan = 1 : 1000
	for (int i = 0; i < 1000; ) {
		int s = i ;
		bool in = false;
		for (int k = 0; i % 4 != 0 || !in; i++, k++) {
			Forward_FFT(i, do_all, D1_raw, table, D2_Resamp[k]);
			if (i % 4 == 3 ) in = true;
		} // for

		
		for (int k = 0; k < 4; k++) {
			FFT(D2_Resamp[k], MAX, 1);
			for (int j = 0; j < MAX; j++)
				D4_absFFT[j] = abs(D2_Resamp[k][j]);
		} // for

		for ( int j = 0 ; j < 4 ; j++,s++ )
			Back_FFT(s, D4_absFFT,  D2_Resamp[j], Img);


	} // for

	/*
	%---------------------------------- - 為了show圖把1D變成2D---------------------------------- -
		Img_temp = zeros(No_PostFFT, No_AScan); % 1024x1000
		for i = 1 : No_AScan
			Img_temp(:, i) = (Img(1, (i - 1)*No_PostFFT + 1 : (i*No_PostFFT)))'; %'
			end
			Img = Img_temp;
	clear Img_temp;
	fprintf(1, ['原始影像最大值 = ', num2str(max(max(Img)), '%.2f'), '(S)\n']);
	%--------------------------------------------------------------------------
		*/

	
	Mat img, img2;
	
	float img_show[1000][1024];
	
	img2.create(1024, 1000, CV_8U);
	int numImgSize = 0;
	fstream fin2;

	for (int i = 0; i < img2.cols; i++) {
		for (int j = 0; j < img2.rows; j++, numImgSize++) {
			img_show[i][j] = Img[numImgSize];
		} // for 
	} //for


	for (int i = 0; i < img2.rows; i++) {
		for (int j = 0; j < img2.cols; j++)
			img2.at<uchar>(i, j) = img_show[j][i];
	} //for

	namedWindow("Image2", WINDOW_AUTOSIZE);
	namedWindow("Image", WINDOW_AUTOSIZE);

	resize(img2,img, Size(512,500),(0.0),(0.0), INTER_LINEAR);
	imshow("Image", img);
	imshow("Image2", img2);
	waitKey(0);
	system("pause");
} // main() 

