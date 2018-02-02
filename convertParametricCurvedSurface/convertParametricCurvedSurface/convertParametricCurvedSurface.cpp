/*UV空間の図形をパラメトリック曲面に変換するプログラム*/

#include <stdlib.h>
#include <gl/glut.h>
#include <gl.h>
#include <iostream>
using namespace std;
#define KEY_ESC 27
#define L 1
#define M 4
#define N 3
#define PAI 3.14159265
#include <math.h>

void glDrawArrowd(double x0, double y0, double z0,
	double x1, double y1, double z1) {
	GLUquadricObj *arrows[2];
	double x2, y2, z2, len, ang;

	x2 = x1 - x0; y2 = y1 - y0; z2 = z1 - z0;
	len = sqrt(x2*x2 + y2 * y2 + z2 * z2);
	if (len != 0.0) {
		ang = acos(z2*len / (sqrt(x2*x2 + y2 * y2 + z2 * z2)*len)) / PAI * 180.0;

		glPushMatrix();
		glTranslated(x0, y0, z0);
		glRotated(ang, -y2 * len, x2*len, 0.0);
		arrows[0] = gluNewQuadric();
		gluQuadricDrawStyle(arrows[0], GLU_FILL);
		gluCylinder(arrows[0], len / 80, len / 80, len*0.9, 8, 8);
		glPushMatrix();
		glTranslated(0.0, 0.0, len*0.9);
		arrows[1] = gluNewQuadric();
		gluQuadricDrawStyle(arrows[1], GLU_FILL);
		gluCylinder(arrows[1], len / 30, 0.0f, len / 10, 8, 8);
		glPopMatrix();
		glPopMatrix();
	}
}

void printString(float x, float y, const char* str, int length) { /* 指定した座標上に文字を表示する関数 */
	float z = -1.0f;
	glRasterPos3f(x, y, z);

	for (int i = 0; i < length; i++) {
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, str[i]);
	}
}

void display()
{
	int i, j, k, m = 0;
	double t = 0.0;
	double circle[36][2] = {};
	double r = 0.2;
	int n = 36;

	glClear(GL_COLOR_BUFFER_BIT);

	glColor3f(0.0, 0.0, 0.0);

	glPushMatrix();
	int s = 0;
	int s2 = 0;
	double u = 1.0;
	double v = 1.0;
	double z = 0.0;

	double p[100][2] = {}; /* パラメトリック曲線P(t) */

	double b[][M] = { { -1.0, 3.0, -3.0, 1.0 },{ 3.0, -6.0, 3.0, 0.0 },{ -3.0, 3.0, 0.0, 0.0 },{ 1.0, 0.0, 0.0, 0.0 } };

	double q[][M][M] = {
		{ { 0.0, 1.0, 2.0, 3.0 },{ 0.0, 1.0, 2.0, 3.0 },{ 0.0, 1.0, 2.0, 3.0 },{ 0.0, 1.0, 2.0, 3.0 } },
	{ { 0.0, 0.0, 0.0, 0.0 },{ 1.0, 1.0, 1.0, 1.0 },{ 2.0, 2.0, 2.0, 2.0 },{ 3.0, 3.0, 3.0, 3.0 } },
	{ { 0.0, 1.0, 1.0, 0.0 },{ 1.0, 2.0, 2.0, 2.0 },{ 1.0, 1.0, 2.0, 2.0 },{ 1.0, 2.0, 1.0, 2.0 } } };

	gluLookAt(2.5, -2.3, 4.5, 1.5, 1.5, 1.568, 0.0, 0.0, 1.0);

	glPointSize(2.0);

	glBegin(GL_POINTS);

	glColor3f(0.0, 1.0, 1.0);

	for (s = 0; s < 100; s++) {

		u = s / 100.0;

		for (s2 = 0; s2 < 100; s2++) {

			v = s2 / 100.0;

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
			glColor3f(0.0, 1.0, 1.0);
			glVertex3d(temp4[0], temp4[1], temp4[2]);
		}
	}
	glEnd();

	for (i = 0; i <= 20; i++) {

		t = i / 20.0;

		double a[M] = { (1.0 - t)*(1.0 - t)*(1.0 - t), 3.0*(1.0 - t)*(1.0 - t)*t, 3.0*(1.0 - t)*t*t, t*t*t };

		double bn[][M] = { a[0], a[1], a[2], a[3] };

		double qn[][2] = { { 0.8, 0.8 },{ 0.6, 0.7 },{ 0.4, 0.9 },{ 0.2, 0.8 } };

		for (j = 0; j < 2; j++) {
			for (k = 0; k < M; k++) {
				p[i][j] += bn[0][k] * qn[k][j]; //Pn(t)を求める
			}
		}
	}

	glColor3f(0.0, 0.0, 0.0);
	glBegin(GL_POINTS);
	for (s = 0; s <= 20; s++) {

		u = p[s][1];
		v = p[s][0];

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
		glVertex3d(temp4[0], temp4[1], temp4[2]);
	}

	glEnd();

	glColor3f(0.0, 0.0, 0.0);  // 円の色(RGBA)
							   // 円周上の座標(x,y)を計算して円を描画
	for (i = 0; i < n; i++) {
		circle[i][0] = r * cos(2.0 * 3.14 * ((double)i / n)) + 0.5;
		circle[i][1] = r * sin(2.0 * 3.14 * ((double)i / n)) + 0.5;
	}

	glBegin(GL_POINTS);
	for (s = 0; s <= 36; s++) {

		u = circle[s][0];
		v = circle[s][1];

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
		glVertex3d(temp4[0], temp4[1], temp4[2]);
	}

	glEnd();

	//直線1
	glBegin(GL_POINTS);
	for (s = 20; s <= 80; s++) {

		u = 0.2;
		v = s / 100.0;

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
		glVertex3d(temp4[0], temp4[1], temp4[2]);
	}

	//直線2
	for (s = 20; s <= 80; s++) {

		u = s / 100.0;

		v = 0.2;

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
		glVertex3d(temp4[0], temp4[1], temp4[2]);
	}

	for (s = 20; s <= 80; s++) {

		u = s / 100.0;
		v = 0.8;

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
		glVertex3d(temp4[0], temp4[1], temp4[2]);
	}
	glEnd();

	glColor3f(0.0, 0.0, 1.0);
	glColor3f(0.0, 0.0, 1.0);
	glDrawArrowd(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
	printString(0.7, 0.2, "x", 1);
	glColor3f(0.0, 1.0, 0.0);
	glDrawArrowd(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	printString(-0.5, 1.0, "y", 1);
	glColor3f(1.0, 0.0, 0.0);
	glDrawArrowd(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
	printString(-1.0, 0.5, "z", 1);

	glPopMatrix();
	glFlush();
}

void myKbd(unsigned char key, int x, int y)
{
	if (key == KEY_ESC) exit(0);
}

void myInit()
{
	int width = 800, height = 800;
	float aspect = (float)width / (float)height;

	glutInitWindowPosition(0, 0);
	glutInitWindowSize(width, height);
	glutInitDisplayMode(GLUT_RGBA);
	glutCreateWindow("Graph of xyz");

	glClearColor(1.0, 1.0, 1.0, 0.0);

	glutKeyboardFunc(myKbd);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(75.0, aspect, 2.0, 7.0);
	glMatrixMode(GL_MODELVIEW);
}


int main(int argc, char **argv)
{
	glutInit(&argc, argv);
	myInit();
	glutDisplayFunc(display);
	glutMainLoop();
	return 0;
}