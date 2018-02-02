#include <stdlib.h>
#include <glut.h>
#include <gl.h>
#include <glu.h>
#include <math.h>
#include <iostream>
using namespace std;
#define KEY_ESC 27
#define L 1
#define M 4
#define N 3
#define PAI 3.14159265

void drawBase(); //パラメトリック曲面の描画を行う関数
void determineNearest(); //ある外部点から曲面上の最近点を求める関数
void drawAxis(); //xyz座標軸を描画する関数

//矢印の描画を行う関数
void glDrawArrowd(double x0, double y0, double z0,double x1, double y1, double z1) {
	GLUquadricObj *arrows[2];
	double x2, y2, z2, len, ang;

	x2 = x1 - x0; y2 = y1 - y0; z2 = z1 - z0;
	len = sqrt(x2*x2 + y2*y2 + z2*z2);
	if (len != 0.0) {
		ang = acos(z2*len / (sqrt(x2*x2 + y2*y2 + z2*z2)*len)) / PAI*180.0;

		glPushMatrix();
		glTranslated(x0, y0, z0);
		glRotated(ang, -y2*len, x2*len, 0.0);
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

//指定した座標上に文字を表示する関数 
void printString(float x, float y, char* str, int length) {
	float z = -1.0f;
	glRasterPos3f(x, y, z);

	for (int i = 0; i < length; i++) {
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, str[i]);
	}
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT);

	glPushMatrix();

	gluLookAt(2.5, -2.3, 4.5, 1.5, 1.5, 1.568, 0.0, 0.0, 1.0);

	//曲面の描画を行う関数
	drawBase();

	//最近点を求める&表示する関数
	determineNearest();

	//座標軸の描画を行う関数
	drawAxis();

	glPopMatrix();
	glFlush();
}

void myKbd(unsigned char key, int x, int y)
{
	if (key == KEY_ESC) exit(0);
}


void drawBase() {
	int s = 0;
	int t = 0;
	double u = 1.0;
	double v = 1.0;
	double b[][M] = { { -1.0, 3.0, -3.0, 1.0 },{ 3.0, -6.0, 3.0, 0.0 },{ -3.0, 3.0, 0.0, 0.0 },{ 1.0, 0.0, 0.0, 0.0 } };
	double q[][M][M] = {
		{ { 0.0, 1.0, 2.0, 3.0 },{ 0.0, 1.0, 2.0, 3.0 },{ 0.0, 1.0, 2.0, 3.0 },{ 0.0, 1.0, 2.0, 3.0 } },
		{ { 0.0, 0.0, 0.0, 0.0 },{ 1.0, 1.0, 1.0, 1.0 },{ 2.0, 2.0, 2.0, 2.0 },{ 3.0, 3.0, 3.0, 3.0 } },
		{ { 0.0, 1.0, 1.0, 0.0 },{ 1.0, 2.0, 2.0, 2.0 },{ 1.0, 1.0, 2.0, 2.0 },{ 1.0, 2.0, 1.0, 2.0 } } };

	glPointSize(2.0);

	glBegin(GL_POINTS);

	glColor3f(0.0, 1.0, 1.0);

	for (s = 0; s < 100; s++) {

		v = s / 100.0;

		for (t = 0; t < 100; t++) {

			u = t / 100.0;

			double temp1[][M] = { 0.0, 0.0, 0.0, 0.0 };
			double temp2[][L] = { { 0.0 },{ 0.0 },{ 0.0 },{ 0.0 } };
			double temp3[][L][M] = {
				{ 0.0, 0.0, 0.0, 0.0 },
				{ 0.0, 0.0, 0.0, 0.0 },
				{ 0.0, 0.0, 0.0, 0.0 } };
			double temp4[N] = { 0.0, 0.0, 0.0 };

			double U[][L] = { { u*u*u },{ u*u },{ u },{ 1.0 } };
			double V[][M] = { v*v*v, v*v, v, 1.0 };

			for (int i = 0; i < L; i++) {
				for (int j = 0; j < M; j++) {
					for (int k = 0; k < M; k++) {
						temp1[i][j] += V[i][k] * b[k][j];
					}
				}
			}

			for (int i = 0; i < M; i++) {
				for (int j = 0; j < L; j++) {
					for (int k = 0; k < M; k++) {
						temp2[i][j] += b[i][k] * U[k][j];
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
			if (s == 40) {
				glColor3f(0.0, 1.0, 0.0);
				glVertex3d(temp4[0], temp4[1], temp4[2]);
			}
			else if (t == 30) {
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
}

void determineNearest() {
	int s = 0;
	int t = 0;
	double u = 1.0;
	double v = 1.0;
	double dxdt, dydt, dzdt, dxdt2, dydt2, dzdt2 = 0.0;
	double b[][M] = { { -1.0, 3.0, -3.0, 1.0 },{ 3.0, -6.0, 3.0, 0.0 },{ -3.0, 3.0, 0.0, 0.0 },{ 1.0, 0.0, 0.0, 0.0 } };
	double q[][M][M] = {
		{ { 0.0, 1.0, 2.0, 3.0 },{ 0.0, 1.0, 2.0, 3.0 },{ 0.0, 1.0, 2.0, 3.0 },{ 0.0, 1.0, 2.0, 3.0 } },
		{ { 0.0, 0.0, 0.0, 0.0 },{ 1.0, 1.0, 1.0, 1.0 },{ 2.0, 2.0, 2.0, 2.0 },{ 3.0, 3.0, 3.0, 3.0 } },
		{ { 0.0, 1.0, 1.0, 0.0 },{ 1.0, 2.0, 2.0, 2.0 },{ 1.0, 1.0, 2.0, 2.0 },{ 1.0, 2.0, 1.0, 2.0 } } };

	int flag = 0;
	double delta_u, delta_v, inner_product, inner_product2 = 0.0;
	double length, length1, length2, cos1, cos2, degree1, degree2 = 0.0;
	double d_u, d_v, d_xy = 0.0;
	u = 0.4;
	v = 0.3;
	double past[N] = { 0.0, 0.0, 0.0 };

	for (s = 0; s < 20; s++) {

		d_u = u + 0.01;
		d_v = v + 0.01;

		double temp1[][M] = { 0.0, 0.0, 0.0, 0.0 };
		double temp2[][L] = { { 0.0 },{ 0.0 },{ 0.0 },{ 0.0 } };
		double temp3[][L][M] = {
			{ 0.0, 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 0.0, 0.0 } };
		double temp4[N] = { 0.0, 0.0, 0.0 };

		double u_temp1[][M] = { 0.0, 0.0, 0.0, 0.0 };
		double u_temp2[][L] = { { 0.0 },{ 0.0 },{ 0.0 },{ 0.0 } };
		double u_temp3[][L][M] = {
			{ 0.0, 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 0.0, 0.0 } };
		double u_temp4[N] = { 0.0, 0.0, 0.0 };

		double v_temp1[][M] = { 0.0, 0.0, 0.0, 0.0 };
		double v_temp2[][L] = { { 0.0 },{ 0.0 },{ 0.0 },{ 0.0 } };
		double v_temp3[][L][M] = {
			{ 0.0, 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 0.0, 0.0 } };
		double v_temp4[N] = { 0.0, 0.0, 0.0 };

		double U[][M] = { u*u*u, u*u, u, 1.0 };
		double V[][L] = { { v*v*v },{ v*v },{ v },{ 1.0 } };

		double d_U[][M] = { d_u*d_u*d_u, d_u*d_u, d_u, 1.0 };

		double d_V[][L] = { { d_v*d_v*d_v },{ d_v*d_v },{ d_v },{ 1.0 } };

		for (int i = 0; i < L; i++) {
			for (int j = 0; j < M; j++) {
				for (int k = 0; k < M; k++) {
					temp1[i][j] += U[i][k] * b[k][j];
					u_temp1[i][j] += d_U[i][k] * b[k][j];
					v_temp1[i][j] += U[i][k] * b[k][j];
				}
			}
		}

		for (int i = 0; i < M; i++) {
			for (int j = 0; j < L; j++) {
				for (int k = 0; k < M; k++) {
					temp2[i][j] += b[i][k] * V[k][j];
					u_temp2[i][j] += b[i][k] * V[k][j];
					v_temp2[i][j] += b[i][k] * d_V[k][j];
				}
			}
		}

		for (int i = 0; i < L; i++) {
			for (int j = 0; j < M; j++) {
				for (int k = 0; k < N; k++) {
					for (int l = 0; l < M; l++) {
						temp3[k][i][j] += temp1[i][l] * q[k][l][j];
						u_temp3[k][i][j] += u_temp1[i][l] * q[k][l][j];
						v_temp3[k][i][j] += v_temp1[i][l] * q[k][l][j];
					}
				}
			}
		}

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < M; j++) {
				temp4[i] += temp3[i][0][j] * temp2[j][0];
				u_temp4[i] += u_temp3[i][0][j] * u_temp2[j][0];
				v_temp4[i] += v_temp3[i][0][j] * v_temp2[j][0];
			}
		}

		if (flag == 1) {
			cout << "注目点\n";
			cout << "UV(" << u << ", " << v << ")\n";
			cout << "xyz(" << temp4[0] << ", " << temp4[1] << ", " << temp4[2] << ")\n\n";
			glBegin(GL_LINES);
			glColor3f(1.0, 1.0, 0.0);
			glVertex3d(1.0, 1.0, 3.0);
			glVertex3d(temp4[0], temp4[1], temp4[2]);
			glEnd();

			d_xy = sqrt((temp4[0] - past[0])*(temp4[0] - past[0]) + (temp4[1] - past[1])*(temp4[1] - past[1]) + (temp4[2] - past[2])*(temp4[2] - past[2]));

			cout << "UV空間での移動量  U : " << delta_u << "  V : " << delta_v << "\n";
			cout << "xyz空間での移動量 : " << d_xy << "\n\n";

			break;
		}

		length = sqrt((1.0 - temp4[0])*(1.0 - temp4[0]) + (1.0 - temp4[1])*(1.0 - temp4[1]) + (3.0 - temp4[2])*(3.0 - temp4[2]));

		dxdt = (u_temp4[0] - temp4[0]) / 0.01;
		dydt = (u_temp4[1] - temp4[1]) / 0.01;
		dzdt = (u_temp4[2] - temp4[2]) / 0.01;
		length1 = sqrt(dxdt*dxdt + dydt*dydt + dzdt*dzdt);

		dxdt2 = (v_temp4[0] - temp4[0]) / 0.01;
		dydt2 = (v_temp4[1] - temp4[1]) / 0.01;
		dzdt2 = (v_temp4[2] - temp4[2]) / 0.01;
		length2 = sqrt(dxdt2*dxdt2 + dydt2*dydt2 + dzdt2*dzdt2);

		inner_product = (1.0 - temp4[0])* dxdt + (1.0 - temp4[1])* dydt + (3.0 - temp4[2])* dzdt;
		inner_product2 = (1.0 - temp4[0])* dxdt2 + (1.0 - temp4[1])* dydt2 + (3.0 - temp4[2])* dzdt2;

		cos1 = inner_product / (length * length1);
		cos2 = inner_product2 / (length * length2);

		degree1 = acos(cos1) * 180.0 / PAI;
		degree2 = acos(cos2) * 180.0 / PAI;

		//sは止めたい場合によって適宜変更する

		//if(s == 14){
		//if(s == 16){
		if (s == 15) {
			cout << "試行回数:" << s + 1 << "\n\n";
			if (fabs(90.0 - degree1) < 0.5 && fabs(90.0 - degree2) < 0.5) {
				cout << "最近点\n";
				cout << temp4[0] << " " << temp4[1] << " " << temp4[2] << "\n\n";
			}
			cout << "角度\n";
			cout << degree1 << " " << degree2 << "\n\n";

			glBegin(GL_LINES);
			glColor3f(0.0, 0.0, 0.0);
			glVertex3d(1.0, 1.0, 3.0);
			glVertex3d(temp4[0], temp4[1], temp4[2]);
			glEnd();

			glBegin(GL_LINES);
			glColor3f(1.0, 0.0, 0.0);
			glVertex3d(temp4[0], temp4[1], temp4[2]);
			glVertex3d(temp4[0] + dxdt, temp4[1] + dydt, temp4[2] + dzdt);
			glEnd();

			glBegin(GL_LINES);
			glColor3f(0.0, 1.0, 0.0);
			glVertex3d(temp4[0], temp4[1], temp4[2]);
			glVertex3d(temp4[0] + dxdt2, temp4[1] + dydt2, temp4[2] + dzdt2);
			glEnd();

			glEnable(GL_BLEND);

			glBegin(GL_POLYGON);
			glColor4f(1.0, 0.0, 0.0, 0.6);
			glVertex3d(temp4[0] + (dxdt - dxdt2) / 3.0, temp4[1] + (dydt - dydt2) / 3.0, temp4[2] + (dzdt - dzdt2) / 3.0);
			glVertex3d(temp4[0] + (dxdt2 - dxdt) / 3.0, temp4[1] + (dydt2 - dydt) / 3.0, temp4[2] + (dzdt2 - dzdt) / 3.0);
			glVertex3d(temp4[0] - (dxdt + dxdt2) / 3.0, temp4[1] - (dydt + dydt2) / 3.0, temp4[2] - (dzdt + dzdt2) / 3.0);
			glEnd();

			glBegin(GL_POLYGON);
			glVertex3d(temp4[0] + (dxdt + dxdt2) / 3.0, temp4[1] + (dydt + dydt2) / 3.0, temp4[2] + (dzdt + dzdt2) / 3.0);
			glVertex3d(temp4[0] + (dxdt - dxdt2) / 3.0, temp4[1] + (dydt - dydt2) / 3.0, temp4[2] + (dzdt - dzdt2) / 3.0);
			glVertex3d(temp4[0] + (dxdt2 - dxdt) / 3.0, temp4[1] + (dydt2 - dydt) / 3.0, temp4[2] + (dzdt2 - dzdt) / 3.0);
			glEnd();

			glDisable(GL_BLEND);

			past[0] = temp4[0];
			past[1] = temp4[1];
			past[2] = temp4[2];

			flag = 1;
		}

		delta_u = 0.7 * inner_product / (temp4[0] * temp4[0] + temp4[1] * temp4[1] + temp4[2] * temp4[2]);
		delta_v = 0.7 * inner_product2 / (temp4[0] * temp4[0] + temp4[1] * temp4[1] + temp4[2] * temp4[2]);

		u = u + delta_u;
		v = v + delta_v;
	}
}

void drawAxis() {
	glColor3f(0.0, 0.0, 1.0);
	glDrawArrowd(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
	printString(0.7, 0.2, "x", 1);
	glColor3f(0.0, 1.0, 0.0);
	glDrawArrowd(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	printString(-0.5, 1.0, "y", 1);
	glColor3f(1.0, 0.0, 0.0);
	glDrawArrowd(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
	printString(-1.0, 0.5, "z", 1);
}

void myInit(char *progname)
{
	int width = 800, height = 800;
	float aspect = (float)width / (float)height;

	glutInitWindowPosition(0, 0);
	glutInitWindowSize(width, height);
	glutInitDisplayMode(GLUT_RGBA);
	glutCreateWindow(progname);

	glClearColor(1.0, 1.0, 1.0, 0.0);
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
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);
	glutMainLoop();
	return(0);
}
