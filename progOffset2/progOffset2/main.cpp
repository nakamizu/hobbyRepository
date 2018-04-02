/*曲面のオフセットを描画するプログラム*/

#include <stdlib.h>
#include <gl/glut.h>
#include <gl.h>

#include <iostream>
using namespace std;

#define KEY_ESC 27
#define L 16
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

//ガウスの消去法
int gauss(double a[L][L + 1], double *qn)
{
	try {
		int i, j, k, l, pivot;
		double x[L];
		double p, q, m, b[1][L + 1];

		for (i = 0; i < L; i++) {
			m = 0;
			pivot = i;

			for (l = i; l < L; l++) {
				/*i列の中で一番値が大きい行を選ぶ*/
				if (fabs(a[l][i]) > m) {
					m = fabs(a[l][i]);
					pivot = l;
				}
			}

			/*pivotがiと違えば、行の入れ替え*/
			if (pivot != i) {
				for (j = 0; j < L + 1; j++) {
					b[0][j] = a[i][j];
					a[i][j] = a[pivot][j];
					a[pivot][j] = b[0][j];
				}
			}
		}

		for (k = 0; k < L; k++) {
			p = a[k][k];              //対角要素を保存
									  /*対角要素は1になることがわかっているので直接代入*/
			a[k][k] = 1;

			for (j = k + 1; j < L + 1; j++) {
				a[k][j] /= p;
			}

			for (i = k + 1; i < L; i++) {
				q = a[i][k];

				for (j = k + 1; j < L + 1; j++) {
					a[i][j] -= q * a[k][j];
				}
				/*0となることがわかっているので直接代入*/
				a[i][k] = 0;
			}
		}

		/*解の計算*/
		for (i = L - 1; i >= 0; i--) {
			x[i] = a[i][L];
			for (j = L - 1; j > i; j--) {
				x[i] -= a[i][j] * x[j];
			}
		}

		for (i = 0; i < L; i++) {
			*qn = x[i];
			++qn;
		}

		return 0;
	}
	catch (...) {
		return -1;
	}
}

//曲面の描画
int drawCurve(double u, double v, double qn[L][N]) {
	try {
		double u_bn[4] = {};
		double v_bn[4] = {};

		u_bn[0] = (1 - u) * (1 - u) * (1 - u);
		u_bn[1] = 3 * (1 - u) * (1 - u) * u;
		u_bn[2] = 3 * (1 - u) * u * u;
		u_bn[3] = u * u * u;

		v_bn[0] = (1 - v) * (1 - v) * (1 - v);
		v_bn[1] = 3 * (1 - v) * (1 - v) * v;
		v_bn[2] = 3 * (1 - v) * v * v;
		v_bn[3] = v * v * v;

		double a[4][3] =
		{
			{ 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0 } };

		double p[3] = { 0.0, 0.0, 0.0 };


		for (int j = 0; j < N; j++) {
			for (int l = 0; l < M; l++) {
				a[0][j] += u_bn[l] * qn[l][j];
				a[1][j] += u_bn[l] * qn[l + 4][j];
				a[2][j] += u_bn[l] * qn[l + 8][j];
				a[3][j] += u_bn[l] * qn[l + 12][j];
			}
		}

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < M; j++) {
				p[i] += a[j][i] * v_bn[j];
			}
		}

		glVertex3d(p[0], p[1], p[2]);

		if (u == 0.7 && v == 0.7) {
			cout << "(" << p[0] << ", " << p[1] << ", " << p[2] << ")\n";
		}

		return 0;
	}
	catch (...) {
		return -1;
	}
}

//接線ベクトルPu、接線ベクトルPvを求める
int curveVector(double u, double v, double qn[M][M][N], double *pos, double *pu, double *pv) {
	try {
		int i = 0;
		int j = 0;
		int k = 0;
		double next_u = 0.0;
		double next_v = 0.0;
		double u_bn[4] = {};
		double v_bn[4] = {};
		double next_u_bn[4] = {};
		double next_v_bn[4] = {};

		u_bn[0] = (1 - u) * (1 - u) * (1 - u);
		u_bn[1] = 3 * (1 - u) * (1 - u) * u;
		u_bn[2] = 3 * (1 - u) * u * u;
		u_bn[3] = u * u * u;

		v_bn[0] = (1 - v) * (1 - v) * (1 - v);
		v_bn[1] = 3 * (1 - v) * (1 - v) * v;
		v_bn[2] = 3 * (1 - v) * v * v;
		v_bn[3] = v * v * v;

		next_u = u + 0.001;
		next_v = v + 0.001;

		next_u_bn[0] = (1 - next_u) * (1 - next_u) * (1 - next_u);
		next_u_bn[1] = 3 * (1 - next_u) * (1 - next_u) * next_u;
		next_u_bn[2] = 3 * (1 - next_u) * next_u * next_u;
		next_u_bn[3] = next_u * next_u * next_u;

		next_v_bn[0] = (1 - next_v) * (1 - next_v) * (1 - next_v);
		next_v_bn[1] = 3 * (1 - next_v) * (1 - next_v) * next_v;
		next_v_bn[2] = 3 * (1 - next_v) * next_v * next_v;
		next_v_bn[3] = next_v * next_v * next_v;

		double temp[M][M] = { { 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0 } };

		double next_temp[M][M] = { { 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0 } };

		double p[3] = { 0.0, 0.0, 0.0 };
		double next_pu[3] = { 0.0, 0.0, 0.0 };
		double next_pv[3] = { 0.0, 0.0, 0.0 };

		for (int j = 0; j < N; j++) {
			for (int k = 0; k < M; k++) {
				for (int l = 0; l < M; l++) {
					temp[k][j] += u_bn[l] * qn[l][k][j];
					next_temp[k][j] += next_u_bn[l] * qn[l][k][j];
				}
			}
		}

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < M; j++) {
				p[i] += temp[j][i] * v_bn[j];
				next_pu[i] += next_temp[j][i] * v_bn[j];
				next_pv[i] += temp[j][i] * next_v_bn[j];

			}
		}


		for (k = 0; k < 3; k++) {
			*pos = p[k];
			*pu = (next_pu[k] - p[k]) / 0.001;
			*pv = (next_pv[k] - p[k]) / 0.001;
			if (k < 2) {
				++pos;
				++pu;
				++pv;
			}
		}
		return 0;
	}
	catch (...) {
		return -1;
	}
}

//制御点を見つける
int curveControl(double u_sample[M], double v_sample[M], double *q, double p[L][N]) {
	try {
		int row;
		int i, j, k;
		int count = 0;

		double q1[L] = {};
		double q2[L] = {};
		double q3[L] = {};

		double answer[][N] = {
			{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 } };

		double bn1[][L + 1] = {
			{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } };

		double bn2[][L + 1] = {
			{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } };

		double bn3[][L + 1] = {
			{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } };

		double u_bn[][M] = {
			{ 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0 } };

		double v_bn[][M] = {
			{ 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0 } };

		for (row = 0; row < M; row++) {
			u_bn[row][0] = (1 - u_sample[row])*(1 - u_sample[row])*(1 - u_sample[row]);
			u_bn[row][1] = 3 * (1 - u_sample[row]) * (1 - u_sample[row]) * u_sample[row];
			u_bn[row][2] = 3 * (1 - u_sample[row]) * u_sample[row] * u_sample[row];
			u_bn[row][3] = u_sample[row] * u_sample[row] * u_sample[row];

			v_bn[row][0] = (1 - v_sample[row])*(1 - v_sample[row])*(1 - v_sample[row]);
			v_bn[row][1] = 3 * (1 - v_sample[row]) * (1 - v_sample[row]) * v_sample[row];
			v_bn[row][2] = 3 * (1 - v_sample[row]) * v_sample[row] * v_sample[row];
			v_bn[row][3] = v_sample[row] * v_sample[row] * v_sample[row];
		}

		for (i = 0; i < 4; i++) {
			for (j = 0; j < 4; j++) {
				for (k = 0; k < 4; k++) {
					bn1[count][k] = u_bn[i][k] * v_bn[j][0];
					bn1[count][k + 4] = u_bn[i][k] * v_bn[j][1];
					bn1[count][k + 8] = u_bn[i][k] * v_bn[j][2];
					bn1[count][k + 12] = u_bn[i][k] * v_bn[j][3];

					bn2[count][k] = u_bn[i][k] * v_bn[j][0];
					bn2[count][k + 4] = u_bn[i][k] * v_bn[j][1];
					bn2[count][k + 8] = u_bn[i][k] * v_bn[j][2];
					bn2[count][k + 12] = u_bn[i][k] * v_bn[j][3];

					bn3[count][k] = u_bn[i][k] * v_bn[j][0];
					bn3[count][k + 4] = u_bn[i][k] * v_bn[j][1];
					bn3[count][k + 8] = u_bn[i][k] * v_bn[j][2];
					bn3[count][k + 12] = u_bn[i][k] * v_bn[j][3];
				}

				count++;
			}
		}

		double b1[L] = {};
		double b2[L] = {};
		double b3[L] = {};

		for (i = 0; i < L; i++) {
			b1[i] = p[i][0];
			b2[i] = p[i][1];
			b3[i] = p[i][2];
		}

		for (i = 0; i < L; i++) {
			bn1[i][L] = b1[i];
			bn2[i][L] = b2[i];
			bn3[i][L] = b3[i];
		}

		gauss(bn1, q1);
		gauss(bn2, q2);
		gauss(bn3, q3);

		for (i = 0; i < L; i++) {
			answer[i][0] = q1[i];
			answer[i][1] = q2[i];
			answer[i][2] = q3[i];
		}

		for (i = 0; i < L; i++) {
			for (j = 0; j < N; j++) {
				*q = answer[i][j];
				++q;
			}
		}

		return 0;
	}
	catch (...) {
		return -1;
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
	int i = 0;
	int j = 0;
	int k = 0;
	int s, t;
	double a = 0.0;
	double b = 0.0;

	double u_d[M] = { 0.1, 0.33, 0.66, 0.99 };/*uのサンプリング位置*/
	double v_d[M] = { 0.1, 0.33, 0.66, 0.99 };/*vのサンプリング位置*/

	double point[][M][N] = {
		{ { 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 } },
	{ { 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 } },
	{ { 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 } },
	{ { 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 } } };

	double pn[][N] = {
		{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },
	{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },
	{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },
	{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 } };

	double q[][M][N] = {
		{ { 0.0, 0.0, 0.0 },{ 1.0, 0.0, 1.0 },{ 2.0, 0.0, 1.0 },{ 3.0, 0.0, 0.0 } },
	{ { 0.0, 1.0, 1.0 },{ 1.0, 1.0, 2.0 },{ 2.0, 1.0, 2.0 },{ 3.0, 1.0, 2.0 } },
	{ { 0.0, 2.0, 1.0 },{ 1.0, 2.0, 1.0 },{ 2.0, 2.0, 2.0 },{ 3.0, 2.0, 2.0 } },
	{ { 0.0, 3.0, 1.0 },{ 1.0, 3.0, 2.0 },{ 2.0, 3.0, 1.0 },{ 3.0, 3.0, 2.0 } } };

	double qn[][N] = {
		{ 0.0, 0.0, 0.0 },{ 1.0, 0.0, 1.0 },{ 2.0, 0.0, 1.0 },{ 3.0, 0.0, 0.0 },
	{ 0.0, 1.0, 1.0 },{ 1.0, 1.0, 2.0 },{ 2.0, 1.0, 2.0 },{ 3.0, 1.0, 2.0 },
	{ 0.0, 2.0, 1.0 },{ 1.0, 2.0, 1.0 },{ 2.0, 2.0, 2.0 },{ 3.0, 2.0, 2.0 },
	{ 0.0, 3.0, 1.0 },{ 1.0, 3.0, 2.0 },{ 2.0, 3.0, 1.0 },{ 3.0, 3.0, 2.0 } };

	double next_q[][N] = {
		{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },
	{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },
	{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },
	{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 } };


	double position[N] = {}; /*位置ベクトル*/
	double u_tangent[N] = {}; /*接線ベクトルPu*/
	double v_tangent[N] = {}; /*接線ベクトルPv*/
	double unit_tangent[N] = {}; /*単位接線ベクトル*/
	double inv_tan[N] = {};	/*進行方向と逆方向である左方向*/
	double direction[N] = {}; /*オフセット方向*/
	double unit_direction[N] = {}; /*オフセット方向単位ベクトル*/

	gluLookAt(-3.0, -2.2, 3.0, 1.3, 1.3, 1.568, 0.0, 0.0, 1.0);

	glEnable(GL_BLEND);

	//オリジナル曲面の描画
	glPointSize(2.0);
	glBegin(GL_POINTS);
	glColor4f(0.0, 1.0, 1.0, 1.0);

	for (s = 0; s < 100; s++) {

		a = s / 100.0;

		for (t = 0; t < 100; t++) {

			b = t / 100.0;

			if (s == 70 && t == 70) {
				cout << "u = 0.7, v = 0.7のときの位置 ";
			}

			drawCurve(a, b, qn);
		}
	}

	glEnd();

	for (s = 0; s < 100; s = s++) {

		if (s == 10 || s == 33 || s == 66 || s == 99) {

			a = s / 100.0;

			for (t = 0; t < 100; t = t++) {

				b = t / 100.0;

				if (t == 10 || t == 33 || t == 66 || t == 99) {
					//位置ベクトル、接線ベクトルを求める.
					curveVector(a, b, q, position, u_tangent, v_tangent);

					unit_tangent[0] = u_tangent[0] / sqrt(u_tangent[0] * u_tangent[0] + u_tangent[1] * u_tangent[1] + u_tangent[2] * u_tangent[2]);
					unit_tangent[1] = u_tangent[1] / sqrt(u_tangent[0] * u_tangent[0] + u_tangent[1] * u_tangent[1] + u_tangent[2] * u_tangent[2]);
					unit_tangent[2] = u_tangent[2] / sqrt(u_tangent[0] * u_tangent[0] + u_tangent[1] * u_tangent[1] + u_tangent[2] * u_tangent[2]);

					//接線ベクトルの描画
					glBegin(GL_LINES);
					glColor4f(0.0, 0.0, 1.0, 1.0);
					glVertex3d(position[0], position[1], position[2]);
					glVertex3d((position[0] + unit_tangent[0]), (position[1] + unit_tangent[1]), (position[2] + unit_tangent[2]));
					glEnd();

					direction[0] = v_tangent[1] * u_tangent[2] - v_tangent[2] * u_tangent[1];
					direction[1] = v_tangent[2] * u_tangent[0] - v_tangent[0] * u_tangent[2];
					direction[2] = v_tangent[0] * u_tangent[1] - v_tangent[1] * u_tangent[0];

					if (direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2] != 0) {
						unit_direction[0] = 0.5 * direction[0] / sqrt(direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2]);
						unit_direction[1] = 0.5 * direction[1] / sqrt(direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2]);
						unit_direction[2] = 0.5 * direction[2] / sqrt(direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2]);
					}
					else {
						unit_direction[0] = 0.0;
						unit_direction[1] = 0.0;
						unit_direction[2] = 0.0;
					}

					point[j][k][0] = position[0] + unit_direction[0];
					point[j][k][1] = position[1] + unit_direction[1];
					point[j][k][2] = position[2] + unit_direction[2];


					//オフセット方向に動線を描画
					glBegin(GL_LINES);
					glColor4f(0.0, 0.0, 0.0, 1.0);
					glVertex3d(position[0], position[1], position[2]);
					glVertex3d(point[j][k][0], point[j][k][1], point[j][k][2]);
					glEnd();

					j++;
				}
			}
			j = 0;
			k++;
		}
	}
	k = 0;


	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			pn[i][j] = point[i][0][j];
			pn[i + 4][j] = point[i][1][j];
			pn[i + 8][j] = point[i][2][j];
			pn[i + 12][j] = point[i][3][j];
		}
	}


	curveControl(u_d, v_d, *next_q, pn);

	//オフセット曲面の描画
	glBegin(GL_POINTS);
	glColor4f(1.0, 1.0, 0.0, 0.5);
	for (s = 0; s < 100; s++) {

		a = s / 100.0;

		for (t = 0; t < 100; t++) {

			b = t / 100.0;

			if (s == 70 && t == 70) {
				cout << "u = 0.7, v = 0.7のときの移動後の位置 ";
			}

			drawCurve(a, b, next_q);
		}
	}
	glEnd();

	glDisable(GL_BLEND);

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
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

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