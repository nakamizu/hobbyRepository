/*パラメトリック曲線の全長を求めるプログラム*/

#include <stdlib.h>
#include <math.h>
#include <iostream>
using namespace std;
#define KEY_ESC 27
#define N 4

class basic {
	int i;
	double y0, y1, y2, y3;
	double temp[N];
public:
	basic();
	void makeBasicFunc(double t, double *p); /* 基底関数Bn(t)を入力tをもとに与えられた配列に入れる関数*/
};

basic::basic()
{
	y0 = 0.0;
	y1 = 0.0;
	y2 = 0.0;
	y3 = 0.0;
	temp[0] = 0.0;
	temp[1] = 0.0;
	temp[2] = 0.0;
	temp[3] = 0.0;
}

void basic::makeBasicFunc(double t, double *p) {

	y0 = (1.0 - t)*(1.0 - t)*(1.0 - t);
	y1 = 3.0*(1.0 - t)*(1.0 - t)*t;
	y2 = 3.0*(1.0 - t)*t*t;
	y3 = t * t*t;

	temp[0] = y0;
	temp[1] = y1;
	temp[2] = y2;
	temp[3] = y3;

	for (i = 0; i < N; i++) {
		*p = temp[i];
		++p;
	}
}

int main() {
	int i, j, k;
	double t, dx, dy;
	double d = 0.0;
	double prex = 1.0;
	double prey = 2.0;
	for (i = 1; i < 10000000; i++) {
		t = i / 10000000.0;

		basic ob;

		double a[N] = { 0.0, 0.0, 0.0, 0.0 };

		ob.makeBasicFunc(t, a);

		double b[][4] = { a[0], a[1], a[2], a[3] }; /*Bn(t)を二次元配列に置き換える*/

		double q[][2] = { { 1.0, 2.0 },{ 2.0, 3.0 },{ 3.0, 1.0 },{ 5.0, 1.0 } }; /*制御点Qn*/

		double p[][2] = { 0.0, 0.0 }; /*曲線の位置ベクトルの初期化*/

		for (j = 0; j < 2; j++) {
			for (k = 0; k < 4; k++) {
				p[0][j] += b[0][k] * q[k][j]; /*曲線の各位置ベクトルを求める*/
			}
		}

		dx = p[0][0] - prex; /* 前の点のx座標prexと現時点でのx座標の微小差dx */
		dy = p[0][1] - prey; /* 前の点のy座標preyと現時点でのy座標の微小差dy */
		d += sqrt(dx*dx + dy * dy); /* 前の点からの距離 */
		prex = p[0][0];
		prey = p[0][1];
	}
	cout << "曲線の長さは[" << d;
	cout << "]です。" << "\n";

	return 0;
}