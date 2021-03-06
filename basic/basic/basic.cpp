/*基底関数Bn(t)を描画するプログラム*/

#include <stdlib.h>
#include <gl.h>
#include <gl/glut.h>
#define KEY_ESC 27

void printString(float x, float y, const char* str, int length) /*指定した座標上に文字を表示する関数*/
{
	float z = -1.0f;
	glRasterPos3f(x, y, z);

	for (int i = 0; i < length; i++) {
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, str[i]);
	}
}

void display()
{
	int  i;
	float x, y;

	glClear(GL_COLOR_BUFFER_BIT);

	glBegin(GL_LINE_STRIP); /*B0(t)の描画*/
	glColor3f(0.0, 0.0, 1.0);
	for (i = 0; i <= 100; i++) {
		x = i / 100.0;
		y = (1 - x)*(1 - x)*(1 - x);
		glVertex2f(x, y);
	}
	glEnd();

	glBegin(GL_LINE_STRIP); /*B1(t)の描画*/
	glColor3f(0.0, 1.0, 0.0);
	for (i = 0; i <= 100; i++) {
		x = i / 100.0;
		y = 3.0 * (1 - x) * (1 - x) * x;
		glVertex2f(x, y);
	}
	glEnd();

	glBegin(GL_LINE_STRIP); /*B2(t)の描画*/
	glColor3f(1.0, 0.0, 0.0);
	for (i = 0; i <= 100; i++) {
		x = i / 100.0;
		y = 3.0 * (1 - x) * x * x;
		glVertex2f(x, y);
	}
	glEnd();

	glBegin(GL_LINE_STRIP); /*B3(t)の描画*/
	glColor3f(1.0, 1.0, 1.0);
	for (i = 0; i <= 100; i++) {
		x = i / 100.0;
		y = x * x * x;
		glVertex2f(x, y);
	}
	glEnd();

	glBegin(GL_LINES); /* 座標軸の描画*/
	glColor3f(1.0, 1.0, 1.0);
	glVertex2f(0.0, 0.0);
	glVertex2f(1.0, 0.0);
	glVertex2f(0.0, 1.0);
	glVertex2f(0.0, 0.0);
	glVertex2f(1.0, 1.0);
	glVertex2f(1.0, 0.0);
	glVertex2f(1.0, 1.0);
	glVertex2f(0.0, 1.0);
	glEnd();

	/* 文字の表示*/
	printString(-0.05, -0.05, "0", 1);
	printString(0.99, -0.05, "1", 1);
	printString(-0.05, 0.99, "1", 1);

	glColor3f(0.0, 0.0, 1.0);
	printString(1.04, 0.5, "- B0(t)", 7);

	glColor3f(0.0, 1.0, 0.0);
	printString(1.04, 0.4, "- B1(t)", 7);

	glColor3f(1.0, 0.0, 0.0);
	printString(1.04, 0.3, "- B2(t)", 7);

	glColor3f(1.0, 1.0, 1.0);
	printString(1.04, 0.2, "- B3(t)", 7);

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
	glutCreateWindow("Graph of basis function");

	glClearColor(0.0, 0.0, 0.0, 1.0);

	glutKeyboardFunc(myKbd);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluOrtho2D(-0.1, 1.2, -0.1, 1.1);
}


int main(int argc, char **argv)
{
	glutInit(&argc, argv);
	myinit();
	glutDisplayFunc(display);
	glutMainLoop();
	return 0;
}
