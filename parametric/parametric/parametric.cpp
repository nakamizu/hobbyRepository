#include <stdlib.h>
#include <gl/glut.h>
#define KEY_ESC 27
#define L 1
#define M 4
#define N 3

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

void printString(float x, float y, char* str, int length) { /* 指定した座標上に文字を表示する関数 */
	float z = -1.0f;
	glRasterPos3f(x, y, z);

	for (int i = 0; i < length; i++) {
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, str[i]);
	}
}

void display()
{
	int s = 0;
	int t = 0;
	int i, j, k = 0;
	double u = 1.0;
	double v = 1.0;

	double b[][M] = { { -1.0, 3.0, -3.0, 1.0 },{ 3.0, -6.0, 3.0, 0.0 },{ -3.0, 3.0, 0.0, 0.0 },{ 1.0, 0.0, 0.0, 0.0 } };
	double q[][M][M] = {
		{ { 0.0, 1.0, 2.0, 3.0 },{ 0.0, 1.0, 2.0, 3.0 },{ 0.0, 1.0, 2.0, 3.0 },{ 0.0, 1.0, 2.0, 3.0 } },
	{ { 0.0, 0.0, 0.0, 0.0 },{ 1.0, 1.0, 1.0, 1.0 },{ 2.0, 2.0, 2.0, 2.0 },{ 3.0, 3.0, 3.0, 3.0 } },
	{ { 0.0, 1.0, 1.0, 0.0 },{ 1.0, 2.0, 2.0, 2.0 },{ 1.0, 1.0, 2.0, 2.0 },{ 1.0, 2.0, 1.0, 2.0 } } };


	glClear(GL_COLOR_BUFFER_BIT);

	glBegin(GL_POINT);

	glColor3f(1.0, 1.0, 0.0);

	for (s = 0; s < 100; s++) {

		u = s / 100.0;

		for (t = 0; t < 100; t++) {

			v = t / 100.0;

			double temp1[][M] = { 0.0, 0.0, 0.0, 0.0 };
			double temp2[][L] = { { 0.0 },{ 0.0 },{ 0.0 },{ 0.0 } };
			double temp3[][L][M] = {
				{ 0.0, 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 0.0, 0.0 } };
			double temp4[N] = { 0.0, 0.0, 0.0 };

			double U[][M] = { u*u*u, u*u, u, 1.0 };
			double V[][L] = { { v*v*v },{ v*v },{ v },{ 1.0 } };

			for (int i = 0; i < L; i++) {
				for (int j = 0; j < M; j++) {
					for (int k = 0; k < M; k++) {
						temp1[i][j] += U[i][k] * b[k][j];
					}
				}
			}

			for (int i = 0; i < M; i++) {
				for (int j = 0; j < L; j++) {
					for (int k = 0; k < M; k++) {
						temp2[i][j] += b[i][k] * V[k][j];
					}
				}
			}

			for (int i = 0; i < L; i++) {
				for (int j = 0; j < M; j++) {
					for (int k = 0; k < N; k++) {
						for (int l = 0; l < M; l++) {
							temp3[k][i][j] += temp1[i][l] * q[k][l][j];
						}
					}
				}
			}

			for (int i = 0; i < N; i++) {
				for (int j = 0; j < M; j++) {
					temp4[i] += temp3[i][0][j] * temp2[j][0];
				}
			}
			glVertex3d(temp4[1], temp4[2], temp4[3]);
		}
	}
	glEnd();

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
	gluOrtho2D(0.5, 5.5, 0.5, 3.5);
}


int main(int argc, char **argv)
{
	glutInit(&argc, argv);
	myinit();
	glutDisplayFunc(display);
	glutMainLoop();
	return 0;
}
