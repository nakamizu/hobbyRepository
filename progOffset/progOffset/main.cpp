/*�Ȑ��̃I�t�Z�b�g��`�悷��v���O����*/

#include <stdlib.h>
#include <gl/glut.h>
#include <gl.h>

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

//�Ȑ��̕`��
int drawCurve(double t, double qn[M][N]) {
	try {
		int i = 0;
		int j = 0;
		double bn[4] = {};

		bn[0] = (1 - t) * (1 - t) * (1 - t);
		bn[1] = 3 * (1 - t) * (1 - t) * t;
		bn[2] = 3 * (1 - t) * t * t;
		bn[3] = t * t * t;

		double p[3] = { 0.0, 0.0, 0.0 };

		for (i = 0; i < 3; i++) {
			for (j = 0; j < 4; j++) {
				p[i] += bn[j] * qn[j][i];
			}
		}
		glVertex3d(p[0], p[1], p[2]);

		if (t == 0.7) {
			cout << "(" << p[0] << ", " << p[1] << ", " << p[2] << ")\n";
		}


		return 0;
	}
	catch (...) {
		return -1;
	}
}

//�ʒu�x�N�g���A�ڐ��x�N�g�������߂�
int curveVector(double t, double qn[M][N], double *pos, double *tan) {
	try {
		int i = 0;
		int j = 0;
		int k = 0;
		double u = 0.0;
		double bn[4] = {};
		double next_bn[4] = {};

		bn[0] = (1 - t) * (1 - t) * (1 - t);
		bn[1] = 3 * (1 - t) * (1 - t) * t;
		bn[2] = 3 * (1 - t) * t * t;
		bn[3] = t * t * t;

		u = t + 0.001;

		next_bn[0] = (1 - u) * (1 - u) * (1 - u);
		next_bn[1] = 3 * (1 - u) * (1 - u) * u;
		next_bn[2] = 3 * (1 - u) * u * u;
		next_bn[3] = u * u * u;

		double p[3] = { 0.0, 0.0, 0.0 };
		double next_p[3] = { 0.0, 0.0, 0.0 };

		for (i = 0; i < 3; i++) {
			for (j = 0; j < 4; j++) {
				p[i] += bn[j] * qn[j][i];
				next_p[i] += next_bn[j] * qn[j][i];
			}
		}

		for (k = 0; k < 3; k++) {
			*pos = p[k];
			*tan = (next_p[k] - p[k]) / 0.001;
			if (k < 2) {
				++tan;
				++pos;
			}
		}
		return 0;
	}
	catch (...) {
		return - 1;
	}
}

//�|���o���@�ɂ�萧��_��������
int curveControl(double sample[M], double *q, double p[M][N]) {
	try {
		int row;
		int i, j, k, m, n;

		double bn[][M] = { { 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0 } };


		for (row = 0; row < M; row++) {
			bn[row][0] = (1 - sample[row])*(1 - sample[row])*(1 - sample[row]);
			bn[row][1] = 3 * (1 - sample[row]) * (1 - sample[row]) * sample[row];
			bn[row][2] = 3 * (1 - sample[row]) * sample[row] * sample[row];
			bn[row][3] = sample[row] * sample[row] * sample[row];
		}

		double e[][M * 2] = { { 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 } };

		double e_reverse[][M] = { { 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0 } };

		for (m = 0; m < M; m++)
		{
			e[0][m] = bn[m][0];
			e[1][m] = bn[m][1];
			e[2][m] = bn[m][2];
			e[3][m] = bn[m][3];
		}

		double pivot, mul;

		//��������t�s������߂�|���o���@�J�n

		// �Ίp������1�Ő��K�����ꂽ�K�i�s������
		for (i = 0; i < M; ++i)
		{
			// �Ίp�����̑I���A���̒l�ōs�����𐳋K��
			pivot = e[i][i];
			for (j = 0; j < M * 2; ++j)
			{
				e[i][j] = (1 / pivot) * e[i][j];
			}
			// �K�i�s������ׂɁA���݂̍s��艺�̍s�ɂ�
			// i��ڂ̐�����0�ɂȂ�悤�Ȋ�{�ό`������
			for (k = i + 1; k < M; ++k)
			{
				mul = e[k][i];
				for (n = i; n < M * 2; ++n)
				{
					e[k][n] = e[k][n] - mul * e[i][n];
				}
			}
		}

		// �������Ɍ������ĊȖ�
		for (i = M - 1; i > 0; --i)
		{
			for (k = i - 1; k >= 0; --k)
			{
				mul = e[k][i];
				for (n = i; n < M * 2; ++n)
				{
					e[k][n] = e[k][n] - mul * e[i][n];
				}
			}
		}

		for (i = 0; i < M; i++) {
			for (j = M; j < 2 * M; j++) {
				e_reverse[i][j - M] = e[i][j];
			}
		}

		//�|���o���@�I��

		for (i = 0; i < M; i++) {
			for (j = 0; j < N; j++) {
				for (k = 0; k < M; k++) {
					*q += e_reverse[k][i] * p[k][j]; //���߂����֐��̋t�s��ƒʉߓ_�̍s��̐ς��琧��_�����߂�
				}

				++q;
			}
		}

		return 0;
	}
	catch (...) {
		return -1;
	}
}

void printString(float x, float y, const char* str, int length) { /* �w�肵�����W��ɕ�����\������֐� */
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
	int i;
	int j = 0;
	double s = 0;
	double d[M] = {0.0, 0.33, 0.66, 0.99};/*�T���v�����O�ʒu*/
	double bn[M] = {}; /*���֐�*/
	double point[][N] = { { 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 } }; /*�ʂ�_*/
	double q[M][N] = { { 1.0, 2.0, 0.0},{ 2.0, 3.0, 1.0},{ 3.0, 1.0, 2.0},{ 5.0, 1.0, 2.0} }; /*����_*/
	double next_q[][N] = { { 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 } };
	double position[N] = {}; /*�ʒu�x�N�g��*/
	double tangent[N] = {}; /*�ڐ��x�N�g��*/
	double unit_tangent[N] = {}; /*�P�ʐڐ��x�N�g��*/
	double inv_tan[N] = {};	/*�i�s�����Ƌt�����ł��鍶����*/
	double direction[N] = {}; /*�I�t�Z�b�g����*/
	double unit_direction[N] = {}; /*�I�t�Z�b�g�����P�ʃx�N�g��*/

	gluLookAt(-3.0, -2.2, 3.0, 1.3, 1.3, 1.568, 0.0, 0.0, 1.0);

	glBegin(GL_LINES);

	glColor3f(1.0, 0.0, 0.0);
	glVertex3d(q[0][0], q[0][1], q[0][2]);
	glVertex3d(q[1][0], q[1][1], q[1][2]);
	glVertex3d(q[1][0], q[1][1], q[1][2]);
	glVertex3d(q[2][0], q[2][1], q[2][2]);
	glVertex3d(q[2][0], q[2][1], q[2][2]);
	glVertex3d(q[3][0], q[3][1], q[3][2]);

	glEnd();

	glPointSize(2.0);

	glBegin(GL_POINTS);

	glColor3f(1.0, 0.0, 0.0);
	glVertex3d(q[0][0], q[0][1], q[0][2]);
	glVertex3d(q[1][0], q[1][1], q[1][2]);
	glVertex3d(q[2][0], q[2][1], q[2][2]);
	glVertex3d(q[3][0], q[3][1], q[3][2]);

	glEnd();


	glBegin(GL_POINTS);

	glColor3f(0.0, 1.0, 1.0);

	for (i = 0; i <= 100; i++) {

		s = i / 100.0;

		//�I���W�i���Ȑ��̕`��
		if (i == 70) {
			cout << "t = 0.7�̂Ƃ��̈ʒu ";
		}

		drawCurve(s, q);

	}

	glEnd();

	for (i = 0; i <= 100; i++){
		s = i / 100.0;

		if (i == 0 || i == 33 || i == 66 || i == 99) {
			
			//�ʒu�x�N�g���A�ڐ��x�N�g�������߂�.
			curveVector(s, q, position, tangent);

			unit_tangent[0] = tangent[0] / sqrt(tangent[0] * tangent[0] + tangent[1] * tangent[1] + tangent[2] * tangent[2]);
			unit_tangent[1] = tangent[1] / sqrt(tangent[0] * tangent[0] + tangent[1] * tangent[1] + tangent[2] * tangent[2]);
			unit_tangent[2] = tangent[2] / sqrt(tangent[0] * tangent[0] + tangent[1] * tangent[1] + tangent[2] * tangent[2]);

			//�ڐ��x�N�g���̕`��
			glBegin(GL_LINES);
			glColor3f(0.0, 0.0, 0.0);
			glVertex3d(position[0], position[1], position[2]);
			glVertex3d((position[0]+ unit_tangent[0]/5.0), (position[1] + unit_tangent[1]/5.0), (position[2] + unit_tangent[2]/5.0));
			glEnd();
			
			inv_tan[0] = -tangent[0];
			inv_tan[1] = -tangent[1];
			inv_tan[2] = -tangent[2];

			direction[0] = inv_tan[1] * 1.0;
			direction[1] = -1.0 * inv_tan[2];
			direction[2] = 0.0;

			unit_direction[0] = 0.5 * direction[0] / sqrt(direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2]);
			unit_direction[1] = 0.5 * direction[1] / sqrt(direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2]);
			unit_direction[2] = 0.5 * direction[2] / sqrt(direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2]);

			point[j][0] = position[0] + unit_direction[0];
			point[j][1] = position[1] + unit_direction[1];
			point[j][2] = position[2] + unit_direction[2];

			//�I�t�Z�b�g�����ɓ�����`��
			glBegin(GL_LINES);
			glColor3f(0.0, 1.0, 0.0);
			glVertex3d(position[0], position[1], position[2]);
			glVertex3d(point[j][0],point[j][1],point[j][2]);
			glEnd();

			j++;
		}
	}

	curveControl(d, *next_q, point);

	glBegin(GL_POINTS);
	glColor3f(0.0, 0.0, 0.0);
	for (i = 0; i <= 100; i++) {
		s = i / 100.0;
		if (i == 70) {
			cout << "t = 0.7�̂Ƃ��̈ړ���̈ʒu ";
		}
		drawCurve(s, next_q);
	}
	glEnd();

	glColor3f(0.0, 0.0, 1.0);
	glDrawArrowd(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
	printString(1.1, 0.5, "x", 1);
	glColor3f(0.0, 1.0, 0.0);
	glDrawArrowd(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	printString(0.5, 1.0, "y", 1);
	glColor3f(1.0, 0.0, 0.0);
	glDrawArrowd(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
	printString(2.0, 1.3, "z", 1);

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



