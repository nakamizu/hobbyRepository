/*通過点とそのときのtの値から制御点を求め、その制御点からパラメトリック曲線を描画するプログラム*/

#include <stdlib.h>
#include <gl/glut.h>
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

basic::basic() /*コンストラクタ*/
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

void printString(float x, float y, const char* str, int length) {
	float z = -1.0f;
	glRasterPos3f(x, y, z);

	for (int i = 0; i < length; i++) {
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, str[i]);
	}
}

void display()
{
	basic ob1, ob2, ob3, ob4, ob5;
	int i, j, k, m, n, l;
	double t, y0, y1, y2, y3;
	float x, y;

	double a[N] = { 0.0, 0.0, 0.0, 0.0 };
	double b[N] = { 0.0, 0.0, 0.0, 0.0 };
	double c[N] = { 0.0, 0.0, 0.0, 0.0 };
	double d[N] = { 0.0, 0.0, 0.0, 0.0 };


	glClear(GL_COLOR_BUFFER_BIT);

	ob1.makeBasicFunc(0.0, a);
	ob2.makeBasicFunc(0.2, b);
	ob3.makeBasicFunc(0.6, c);
	ob4.makeBasicFunc(1.0, d);

	double e[][N * 2] = { { 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 },
	{ 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 },
	{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 },
	{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 } };

	double e_reverse[][N] = { { 0.0, 0.0, 0.0, 0.0 },
	{ 0.0, 0.0, 0.0, 0.0 },
	{ 0.0, 0.0, 0.0, 0.0 },
	{ 0.0, 0.0, 0.0, 0.0 } };

	for (m = 0; m < N; m++)
	{
		e[0][m] = a[m];
		e[1][m] = b[m];
		e[2][m] = c[m];
		e[3][m] = d[m];
	}

	double pivot, mul;

	//ここから逆行列を求める掃き出し法開始

	// 対角成分が1で正規化された階段行列を作る
	for (i = 0; i < N; ++i)
	{
		// 対角成分の選択、この値で行成分を正規化
		pivot = e[i][i];
		for (j = 0; j < N * 2; ++j)
		{
			e[i][j] = (1 / pivot) * e[i][j];
		}

		// 階段行列を作る為に、現在の行より下の行につ
		// i列目の成分が0になるような基本変形をする
		for (k = i + 1; k < N; ++k)
		{
			mul = e[k][i];
			for (n = i; n < N * 2; ++n)
			{
				e[k][n] = e[k][n] - mul * e[i][n];
			}
		}
	}

	// 下から上に向かって簡約化
	for (i = N - 1; i > 0; --i)
	{
		for (k = i - 1; k >= 0; --k)
		{
			mul = e[k][i];
			for (n = i; n < N * 2; ++n)
			{
				e[k][n] = e[k][n] - mul * e[i][n];
			}
		}
	}

	for (i = 0; i < N; i++) {
		for (j = N; j < 2 * N; j++) {
			e_reverse[i][j - N] = e[i][j];
		}
	}

	//掃き出し法終了

	double r[][2] = { { 1.0, 2.0 },{ 2.0, 3.0 },{ 3.0, 1.0 },{ 5.0, 1.0 } }; //通過点

	double q[][2] = { { 0.0, 0.0 },{ 0.0, 0.0 },{ 0.0, 0.0 },{ 0.0, 0.0 } }; //制御点を入れる二次元配列

	glPointSize(5.0);
	glBegin(GL_POINTS);
	glColor3f(1.0, 0.0, 0.0);
	for (i = 0; i < 4; i++) {
		cout << "( ";
		for (j = 0; j < 2; j++) {
			for (k = 0; k < 4; k++) {
				q[i][j] += e_reverse[i][k] * r[k][j]; //求めた基底関数の逆行列と通過点の行列の積から制御点を求める
			}
			cout << q[i][j] << " ";
		}
		cout << ")\n";
		glVertex2d(q[i][0], q[i][1]); /*制御点の描画*/
	}
	glEnd();

	glPointSize(1.0);

	glBegin(GL_LINE_STRIP);
	glColor3f(1.0, 1.0, 0.0);


	double bn[][N] = { 0.0, 0.0, 0.0, 0.0 };

	for (i = 0; i <= 100; i++) { //制御点が求まったのでこれをもとに曲線の描画
		basic ob5;

		t = i / 100.0;

		double temp[N] = { 0.0, 0.0, 0.0, 0.0 };

		ob5.makeBasicFunc(t, temp);

		double bn[][N] = { temp[0], temp[1], temp[2], temp[3] };

		double p[][2] = { 0.0, 0.0 };

		for (j = 0; j < 2; j++) {
			for (k = 0; k < 4; k++) {
				p[0][j] += bn[0][k] * q[k][j];
			}
		}
		glVertex2d(p[0][0], p[0][1]);
	}
	glEnd();

	glPointSize(1.0); /*座標軸の表示*/
	glBegin(GL_LINES);
	glColor3f(1.0, 1.0, 1.0);
	glVertex2f(0.0, -3.0);
	glVertex2f(5.0, -3.0);
	glVertex2f(0.0, 6.0);
	glVertex2f(0.0, -3.0);
	glVertex2f(5.0, 6.0);
	glVertex2f(5.0, -3.0);
	glVertex2f(5.0, 6.0);
	glVertex2f(0.0, 6.0);
	glEnd();

	glEnable(GL_LINE_STIPPLE); /*破線でグリッドの表示*/
	glLineStipple(1, 0xF0F0);

	glBegin(GL_LINES);
	for (x = 1.0; x < 5.0; x++) {
		glVertex2f(x, -3.0);
		glVertex2f(x, 6.0);
	}
	for (y = -2.0; y < 6.0; y++) {
		glVertex2f(0.0, y);
		glVertex2f(5.0, y);
	}
	glEnd();

	glDisable(GL_LINE_STIPPLE);

	glColor3f(1.0, 1.0, 1.0); /*座標値表示*/
	printString(0.00, -3.35, "0", 1);
	printString(1.00, -3.35, "1", 1);
	printString(2.00, -3.35, "2", 1);
	printString(3.00, -3.35, "3", 1);
	printString(4.00, -3.35, "4", 1);
	printString(5.00, -3.35, "5", 1);

	printString(-0.15, 6.00, "6", 1);
	printString(-0.15, 5.00, "5", 1);
	printString(-0.15, 4.00, "4", 1);
	printString(-0.15, 3.00, "3", 1);
	printString(-0.15, 2.00, "2", 1);
	printString(-0.15, 1.00, "1", 1);
	printString(-0.15, 0.00, "0", 1);
	printString(-0.25, -1.00, "-1", 2);
	printString(-0.25, -2.00, "-2", 2);
	printString(-0.25, -3.00, "-3", 3);

	glFlush();
}

void myKbd(unsigned char key, int x, int y)
{
	if (key == KEY_ESC) exit(0);
}


void myinit()
{
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(500, 400);
	glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
	glutCreateWindow("Graph of Parametric curve");

	glClearColor(0.0, 0.0, 0.0, 1.0);

	glutKeyboardFunc(myKbd);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluOrtho2D(-0.5, 5.5, -3.5, 6.5);
}

int main(int argc, char **argv)
{
	glutInit(&argc, argv);
	myinit();
	glutDisplayFunc(display);
	glutMainLoop();
	return 0;
}



