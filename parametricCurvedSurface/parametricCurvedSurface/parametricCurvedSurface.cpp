/*アイソライン付きでパラメトリック曲面を描画するプログラム*/

#include <stdlib.h>
#include <gl/glut.h>
#include <iostream>
using namespace std;
#define KEY_ESC 27
#define L 1
#define M 4
#define N 3
#define PAI 3.1415926536

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

void polarview(float distance, float twist, float elevation, float azimuth) {
	glTranslatef(0.0, 0.0, -distance);
	glRotatef(-twist, 0.0, 0.0, 1.0);
	glRotatef(-elevation, 1.0, 0.0, 0.0);
	glRotatef(-azimuth, 0.0, 1.0, 0.0);
}

void display()
{
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT);

	glPushMatrix();
	int s = 0;
	int t = 0;
	double u = 1.0;
	double v = 1.0;
	double x, y, z, x1, y1, z1, x2, y2, z2 = 0.0;
	double dxdt, dydt, dzdt, dxdt1, dydt1, dzdt1, vx, vy, vz = 0.0;
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

			if (s == 30 && t == 40) {
				x = temp4[0];
				y = temp4[1];
				z = temp4[2];
			}
			else if (s == 50 && t == 50) {
				cout << temp4[0] << " " << temp4[1] << " " << temp4[2] << "\n";
			}
			else if (s == 31 && t == 40) {
				x1 = temp4[0];
				y1 = temp4[1];
				z1 = temp4[2];
			}
			else if (s == 30 && t == 41) {
				x2 = temp4[0];
				y2 = temp4[1];
				z2 = temp4[2];
			}
			else if (s == 30) {
				glColor3f(0.0, 1.0, 0.0);
				glVertex3d(temp4[0], temp4[1], temp4[2]);
			}
			else if (t == 40) {
				glColor3f(1.0, 0.0, 0.0);
				glVertex3d(temp4[0], temp4[1], temp4[2]);
			}
			else {
				glColor3f(0.0, 1.0, 1.0);
				glVertex3d(temp4[0], temp4[1], temp4[2]);
			}
		}
	}
	glEnd();

	dxdt = (x1 - x) / 0.01;
	dydt = (y1 - y) / 0.01;
	dzdt = (z1 - z) / 0.01;

	cout << "V = 0.3の接線ベクトル\n";
	cout << dxdt << " " << dydt << " " << dzdt << "\n";

	dxdt1 = (x2 - x) / 0.01;
	dydt1 = (y2 - y) / 0.01;
	dzdt1 = (z2 - z) / 0.01;

	cout << "U = 0.4の接線ベクトル\n";
	cout << dxdt1 << " " << dydt1 << " " << dzdt1 << "\n";

	glBegin(GL_LINES);
	glColor3f(1.0, 0.0, 0.0);
	glVertex3d(x, y, z);
	glVertex3d(x + dxdt, y + dydt, z + dzdt);
	glEnd();

	glBegin(GL_LINES);
	glColor3f(0.0, 1.0, 0.0);
	glVertex3d(x, y, z);
	glVertex3d(x + dxdt1, y + dydt1, z + dzdt1);
	glEnd();

	glBegin(GL_LINES);
	glColor3f(0.0, 0.0, 0.0);
	glVertex3d(x, y, z);
	vx = dydt * dzdt1 - dzdt * dydt1;
	vy = dzdt * dxdt1 - dxdt * dzdt1;
	vz = dxdt * dydt1 - dydt * dxdt1;
	glVertex3d(x + vx / sqrt(vx*vx + vy * vy + vz * vz), y + vy / sqrt(vx*vx + vy * vy + vz * vz), z + vz / sqrt(vx*vx + vy * vy + vz * vz));
	glEnd();

	glColor3f(1.0, 1.0, 0.0);
	glPointSize(5.0);
	glBegin(GL_POINTS);
	glVertex3d(x, y, z);
	glEnd();

	glBegin(GL_LINES);
	glColor3f(0.0, 0.0, 1.0);

	glVertex3d(0.0, 0.0, 0.0);
	glVertex3d(1.0, 0.0, 1.0);
	glVertex3d(1.0, 0.0, 1.0);
	glVertex3d(2.0, 0.0, 1.0);
	glVertex3d(2.0, 0.0, 1.0);
	glVertex3d(3.0, 0.0, 0.0);

	glVertex3d(0.0, 1.0, 1.0);
	glVertex3d(1.0, 1.0, 2.0);
	glVertex3d(1.0, 1.0, 2.0);
	glVertex3d(2.0, 1.0, 2.0);
	glVertex3d(2.0, 1.0, 2.0);
	glVertex3d(3.0, 1.0, 2.0);

	glVertex3d(0.0, 2.0, 1.0);
	glVertex3d(1.0, 2.0, 1.0);
	glVertex3d(1.0, 2.0, 1.0);
	glVertex3d(2.0, 2.0, 2.0);
	glVertex3d(2.0, 2.0, 2.0);
	glVertex3d(3.0, 2.0, 2.0);

	glVertex3d(0.0, 3.0, 1.0);
	glVertex3d(1.0, 3.0, 2.0);
	glVertex3d(1.0, 3.0, 2.0);
	glVertex3d(2.0, 3.0, 1.0);
	glVertex3d(2.0, 3.0, 1.0);
	glVertex3d(3.0, 3.0, 2.0);

	glVertex3d(0.0, 0.0, 0.0);
	glVertex3d(0.0, 1.0, 1.0);
	glVertex3d(0.0, 1.0, 1.0);
	glVertex3d(0.0, 2.0, 1.0);
	glVertex3d(0.0, 2.0, 1.0);
	glVertex3d(0.0, 3.0, 1.0);

	glVertex3d(1.0, 0.0, 1.0);
	glVertex3d(1.0, 1.0, 2.0);
	glVertex3d(1.0, 1.0, 2.0);
	glVertex3d(1.0, 2.0, 1.0);
	glVertex3d(1.0, 2.0, 1.0);
	glVertex3d(1.0, 3.0, 2.0);

	glVertex3d(2.0, 0.0, 1.0);
	glVertex3d(2.0, 1.0, 2.0);
	glVertex3d(2.0, 1.0, 2.0);
	glVertex3d(2.0, 2.0, 2.0);
	glVertex3d(2.0, 2.0, 2.0);
	glVertex3d(2.0, 3.0, 1.0);

	glVertex3d(3.0, 0.0, 0.0);
	glVertex3d(3.0, 1.0, 2.0);
	glVertex3d(3.0, 1.0, 2.0);
	glVertex3d(3.0, 2.0, 2.0);
	glVertex3d(3.0, 2.0, 2.0);
	glVertex3d(3.0, 3.0, 2.0);
	glEnd();

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

void myInit(char *progname)
{
	int width = 800, height = 800;
	float aspect = (float)width / (float)height;

	glutInitWindowPosition(0, 0);
	glutInitWindowSize(width, height);
	glutInitDisplayMode(GLUT_RGBA);
	glutCreateWindow(progname);

	glClearColor(0.0, 0.0, 0.0, 1.0);
	glutKeyboardFunc(myKbd);

	glMatrixMode(GL_PROJECTION);

	glLoadIdentity();
	gluPerspective(75.0, aspect, 1.0, 10.0);
	glMatrixMode(GL_MODELVIEW);
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	myInit(argv[0]);
	glutDisplayFunc(display);
	glutMainLoop();
	return(0);
}