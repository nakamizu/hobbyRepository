/*パラメトリック曲線上に接線を引くプログラム*/

#include <stdlib.h>
#include <gl/glut.h>
#include <math.h>
#include <iostream>
using namespace std;
#define KEY_ESC 27
#define M_PI 3.1415926536
#define N 4
#define ARRAYSIZE 100

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

void printString(float x, float y, const char* str, int length) { /*指定された座標上に文字を表示する関数*/
	float z = -1.0f;
	glRasterPos3f(x, y, z);

	for (int i = 0; i < length; i++) {
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, str[i]);
	}
}

void glDrawArrowd(double x0, double y0, double z0,
	double x1, double y1, double z1) /*指定された座標に矢印を表示する関数*/
{
	GLUquadricObj *arrows[2];
	double x2, y2, z2, len, ang;

	x2 = x1 - x0; y2 = y1 - y0; z2 = z1 - z0;
	len = sqrt(x2*x2 + y2 * y2 + z2 * z2);
	if (len != 0.0) {
		ang = acos(z2*len / (sqrt(x2*x2 + y2 * y2 + z2 * z2)*len)) / M_PI * 180.0;

		glPushMatrix();
		glTranslated(x0, y0, z0);
		glRotated(ang, -y2 * len, x2*len, 0.0);
		arrows[0] = gluNewQuadric();
		gluQuadricDrawStyle(arrows[0], GLU_FILL);
		gluCylinder(arrows[0], 0.01, 0.01, (len - 0.05), 8, 8);
		glPushMatrix();
		glTranslated(0.0, 0.0, (len - 0.05));
		arrows[1] = gluNewQuadric();
		gluQuadricDrawStyle(arrows[1], GLU_FILL);
		gluCylinder(arrows[1], 0.05, 0.0f, 0.05, 8, 8);
		glPopMatrix();
		glPopMatrix();
	}
}

void display()
{
	int i = 1;
	int j, k, l;
	double t, dx, dy, dxdt, dydt, nx, ny;
	double prex = 1.0;
	double prey = 2.0;

	double p[ARRAYSIZE][2] = {};

	glClear(GL_COLOR_BUFFER_BIT);

	for (i = 1; i <= ARRAYSIZE; i++) {
		basic ob;

		double rest = i % (ARRAYSIZE / 10); /*剰余で0.1単位かどうかを判断*/

		t = i / (double)ARRAYSIZE;

		double a[N] = { 0.0, 0.0, 0.0, 0.0 };

		ob.makeBasicFunc(t, a);

		double b[][4] = { a[0], a[1], a[2], a[3] };

		double q[][2] = { { 1.0, 2.0 },{ 2.0, 3.0 },{ 3.0, 1.0 },{ 5.0, 1.0 } };

		for (j = 0; j < 2; j++) {
			for (k = 0; k < 4; k++) {
				p[i - 1][j] += b[0][k] * q[k][j];
			}
		}

		dx = p[i - 1][0] - prex;
		dy = p[i - 1][1] - prey;
		dxdt = dx / t;
		nx = dxdt + p[i - 1][0];
		dydt = dy / t;
		ny = dydt + p[i - 1][1];

		if (rest == 0.0) { /*tが0.1刻みの値になったときに接線ベクトルと位置ベクトルを求める&図示*/
			glColor3f(1.0, 0.0, 0.0);
			glDrawArrowd(p[i - 1][0], p[i - 1][1], 0.0, nx, ny, 0.0);
			glColor3f(0.0, 0.0, 1.0);
			glDrawArrowd(0.0, 0.0, 0.0, p[i - 1][0], p[i - 1][1], 0.0);
			cout << "t = " << t << "のとき\n";
			cout << "位置ベクトル : (" << p[i - 1][0] << ", " << p[i - 1][1] << ")\n";
			cout << "接線ベクトル : (" << dxdt << ", " << dydt << ") 大きさ:" << sqrt(dxdt*dxdt + dydt * dydt);
			cout << "\n\n";
		}

		prex = p[i - 1][0];
		prey = p[i - 1][1];
	}

	glColor3f(1.0, 1.0, 0.0);
	glBegin(GL_LINE_STRIP);
	for (l = 0; l < ARRAYSIZE; l++) {
		glVertex2d(p[l][0], p[l][1]);
	}
	glEnd();

	glBegin(GL_LINES);
	glColor3f(1.0, 1.0, 1.0);
	glVertex2f(0.0, 0.0);
	glVertex2f(5.0, 0.0);
	glVertex2f(0.0, 3.0);
	glVertex2f(0.0, 0.0);
	glVertex2f(5.0, 3.0);
	glVertex2f(5.0, 0.0);
	glVertex2f(5.0, 3.0);
	glVertex2f(0.0, 3.0);
	glEnd();

	printString(-0.20, -0.20, "0", 1);
	printString(1.00, -0.20, "1", 1);
	printString(2.00, -0.20, "2", 1);
	printString(3.00, -0.20, "3", 1);
	printString(4.00, -0.20, "4", 1);
	printString(5.00, -0.20, "5", 1);

	printString(-0.20, 3.0, "3", 1);
	printString(-0.20, 2.0, "2", 1);
	printString(-0.20, 1.0, "1", 1);

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
	gluOrtho2D(-0.30, 5.3, -0.30, 3.3);
}


int main(int argc, char **argv)
{
	glutInit(&argc, argv);
	myinit();
	glutDisplayFunc(display);
	glutMainLoop();
	return 0;
}