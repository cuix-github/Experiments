#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <glew.h>
#include <glut.h>
#include "Helpers.h"

#pragma comment(lib, "glut32.lib")
#pragma comment(lib, "glew32.lib")

#define I2D(ni,i,j) ((i + (ni)*(j)))

GLuint gl_PBO, gl_Tex;

float *f0, *f1, *f2, *f3, *f4, *f5, *f6, *f7, *f8;
float *tmpf0, *tmpf1, *tmpf2, *tmpf3, *tmpf4, *tmpf5, *tmpf6, *tmpf7, *tmpf8;
float *cmap, *plotvar;
float dt;
int *solid;
int win_x, win_y;
unsigned int *cmap_rgba, *plot_rgba;

float *vel_u, *vel_v;
float *w;
float tau, faceq1, faceq2, faceq3;
float vxin, roout;
float width, height;
float streamline_length;
int N = 6;
int ni, nj;
int ncol;
int ipos_old, jpos_old, draw_solid_flag;

void display(void);
void resize(int w, int h);
void mouse(int button, int state, int x, int y);
void mouse_motion(int x, int y);
void shutdown(void);

void stream(void);
void collide(void);
void solid_BC(void);
void per_BC(void);
void in_BC(void);
void ex_BC_crude(void);
void apply_BCs(void);

unsigned int get_col(float min, float max, float val);

static void pre_display(void){
	glViewport(0, 0, win_x, win_y);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, 1.0, 0.0, 1.0);
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// Make the pixel looks round.
	glEnable(GL_POINT_SMOOTH);
}

static void draw_vector_field(float * u, float * v, float lineWidth, float r, float g, float b)
{
	int i, j;
	float x, y, h;
	int N = ni = nj;
	h = 1.0f / N;

	glColor3f(r, g, b);
	glLineWidth(lineWidth);

	glBegin(GL_LINES);
	for (j = 1; j <= N; j++) {
		x = (j - 0.5f)*h;
		for (i = 1; i <= N; i++) {
			y = (i - 0.5f)*h;
			glVertex2f(y, x);
			glVertex2f(y + u[I2D(N, i, j)] * streamline_length / N,
				x + v[I2D(N, i, j)] * streamline_length / N);
		}
	}
	glEnd();

	glBegin(GL_POINTS);
	glPointSize(1.0f);
	for (j = 1; j <= N; j++) {
		x = (j - 0.5f)*h;
		for (i = 1; i <= N; i++) {
			y = (i - 0.5f)*h;
			glVertex2f(x, y);
		}
	}
	glEnd();
}

static void post_display(void){
	glutSwapBuffers();
}

void shutdown(void){
	if (f0) free(f0);
	if (f1) free(f1);
	if (f2) free(f2);
	if (f3) free(f3);
	if (f4) free(f4);
	if (f5) free(f5);
	if (f6) free(f6);
	if (f7) free(f7);
	if (f8) free(f8);

	if (tmpf0) free(tmpf0);
	if (tmpf1) free(tmpf1);
	if (tmpf2) free(tmpf2);
	if (tmpf3) free(tmpf3);
	if (tmpf4) free(tmpf4);
	if (tmpf5) free(tmpf5);
	if (tmpf6) free(tmpf6);
	if (tmpf7) free(tmpf7);
	if (tmpf8) free(tmpf8);

	if (vel_u) free(vel_u);
	if (vel_v) free(vel_v);
	if (w) free(w);
	if (plotvar) free(plotvar);
	if (plot_rgba) free(plot_rgba);
	if (solid) free(solid);
}

void stream(void)
{
	int i, j, im1, ip1, jm1, jp1, i0;

	for (j = 1; j<=nj; j++) {
		jm1 = j - 1;
		jp1 = j + 1;
		if (j == 0) jm1 = 0;
		if (j == (nj - 1)) jp1 = nj - 1;
		for (i = 1; i<=ni; i++) {
			i0 = I2D(ni, i, j);
			im1 = i - 1;
			ip1 = i + 1;
			if (i == 0) im1 = 0;
			if (i == (ni - 1)) ip1 = ni - 1;
			tmpf1[i0] = f1[I2D(ni, im1, j)];
			tmpf2[i0] = f2[I2D(ni, i, jp1)];
			tmpf3[i0] = f3[I2D(ni, ip1, j)];
			tmpf4[i0] = f4[I2D(ni, i, jm1)];
			tmpf5[i0] = f5[I2D(ni, im1, jp1)];
			tmpf6[i0] = f6[I2D(ni, ip1, jp1)];
			tmpf7[i0] = f7[I2D(ni, ip1, jm1)];
			tmpf8[i0] = f8[I2D(ni, im1, jm1)];
		}
	}

	for (j = 1; j<=nj; j++) {
		for (i = 1; i<=ni; i++) {
			i0 = I2D(ni, i, j);
			f1[i0] = tmpf1[i0];
			f2[i0] = tmpf2[i0];
			f3[i0] = tmpf3[i0];
			f4[i0] = tmpf4[i0];
			f5[i0] = tmpf5[i0];
			f6[i0] = tmpf6[i0];
			f7[i0] = tmpf7[i0];
			f8[i0] = tmpf8[i0];
		}
	}
}

void collide(void)
{
	int i, j, i0;
	float ro, rovx, rovy, vx, vy, v_sq_term;
	float f0eq, f1eq, f2eq, f3eq, f4eq, f5eq, f6eq, f7eq, f8eq;
	float rtau, rtau1;
	float h = 1.f / ni;
	float c = h / dt;
	float c2 = pow(c, 2);

	rtau = 1.f / tau;
	rtau1 = 1.f - rtau;

	for (j = 1; j<=nj; j++) {
		for (i = 1; i<=ni; i++) {

			i0 = I2D(ni, i, j);

			ro = f0[i0] + f1[i0] + f2[i0] + f3[i0] + f4[i0] + f5[i0] + f6[i0] + f7[i0] + f8[i0];
			rovx = f1[i0] - f3[i0] + f5[i0] - f6[i0] - f7[i0] + f8[i0];
			rovy = f2[i0] - f4[i0] + f5[i0] + f6[i0] - f7[i0] - f8[i0];
			vx = rovx / ro;
			vy = rovy / ro;
			vel_u[i0] = vx;
			vel_v[i0] = vy;
		}
	}

	//computeVortConf(N, vel_u, vel_v, 0.01f, 0.2f);

	for (int j = 1; j <= nj; j++){
		for (int i = 1; i <= ni; i++){
			i0 = I2D(ni, i, j);

			ro = f0[i0] + f1[i0] + f2[i0] + f3[i0] + f4[i0] + f5[i0] + f6[i0] + f7[i0] + f8[i0];

			vx = vel_u[i0];
			vy = vel_v[i0];

			v_sq_term = 1.5f*(vx*vx + vy*vy);

			f0eq = ro * faceq1 * (1.f - v_sq_term);
			f1eq = ro * faceq2 * (1.f + 3.f*vx + 4.5f*vx*vx - v_sq_term);
			f2eq = ro * faceq2 * (1.f + 3.f*vy + 4.5f*vy*vy - v_sq_term);
			f3eq = ro * faceq2 * (1.f - 3.f*vx + 4.5f*vx*vx - v_sq_term);
			f4eq = ro * faceq2 * (1.f - 3.f*vy + 4.5f*vy*vy - v_sq_term);
			f5eq = ro * faceq3 * (1.f + 3.f*(vx + vy) + 4.5f*(vx + vy)*(vx + vy) - v_sq_term);
			f6eq = ro * faceq3 * (1.f + 3.f*(-vx + vy) + 4.5f*(-vx + vy)*(-vx + vy) - v_sq_term);
			f7eq = ro * faceq3 * (1.f + 3.f*(-vx - vy) + 4.5f*(-vx - vy)*(-vx - vy) - v_sq_term);
			f8eq = ro * faceq3 * (1.f + 3.f*(vx - vy) + 4.5f*(vx - vy)*(vx - vy) - v_sq_term);

			f0[i0] = rtau1 * f0[i0] + rtau * f0eq;
			f1[i0] = rtau1 * f1[i0] + rtau * f1eq;
			f2[i0] = rtau1 * f2[i0] + rtau * f2eq;
			f3[i0] = rtau1 * f3[i0] + rtau * f3eq;
			f4[i0] = rtau1 * f4[i0] + rtau * f4eq;
			f5[i0] = rtau1 * f5[i0] + rtau * f5eq;
			f6[i0] = rtau1 * f6[i0] + rtau * f6eq;
			f7[i0] = rtau1 * f7[i0] + rtau * f7eq;
			f8[i0] = rtau1 * f8[i0] + rtau * f8eq;
		}
	}
}

void solid_BC(void)
{
	int i, j, i0;
	float f1old, f2old, f3old, f4old, f5old, f6old, f7old, f8old;

	for (j = 1; j<=nj; j++){
		for (i = 1; i<=ni; i++){
			i0 = I2D(ni, i, j);
			if (solid[i0] == 0) {
				f1old = f1[i0];
				f2old = f2[i0];
				f3old = f3[i0];
				f4old = f4[i0];
				f5old = f5[i0];
				f6old = f6[i0];
				f7old = f7[i0];
				f8old = f8[i0];

				f1[i0] = f3old;
				f2[i0] = f4old;
				f3[i0] = f1old;
				f4[i0] = f2old;
				f5[i0] = f7old;
				f6[i0] = f8old;
				f7[i0] = f5old;
				f8[i0] = f6old;
			}
		}
	}
}

void per_BC(void)
{
	int i0, i1, i;

	for (i = 0; i<ni; i++){
		i0 = I2D(ni, i, 0);
		i1 = I2D(ni, i, nj - 1);
		f2[i0] = f2[i1];
		f5[i0] = f5[i1];
		f6[i0] = f6[i1];
		f4[i1] = f4[i0];
		f7[i1] = f7[i0];
		f8[i1] = f8[i0];
	}
}

void in_BC(void)
{
	int i0, j;
	float f4new, f7new, f8new, vx_term;
	float h = 1 / ni;
	float c = h / dt;
	float c2 = pow(c, 2);

	vx_term = 1.f + 3.f*vxin+ 4.5f*vxin*vxin;
	f4new = roout * faceq2 * vx_term;
	f7new = roout * faceq3 * vx_term;
	f8new = roout * faceq3 * vx_term;

	int emit_idx = I2D(ni, ni / 2, 10);

	for (j = nj / 2 - 2; j <= nj / 2 + 2; j++){
		i0 = I2D(ni, j, 20);
		f4[emit_idx] = f4new;
		f7[emit_idx] = f7new;
		f8[emit_idx] = f8new;
	}

}

void ex_BC_crude(void)
{
	int i0, i1, j;

	for (j = 0; j<nj; j++){
		i0 = I2D(ni, ni - 1, j);
		i1 = i0 - 1;
		f2[i0] = f2[i1];
		f5[i0] = f5[i1];
		f6[i0] = f6[i1];
	}
}

void apply_BCs(void)
{
	per_BC();

	solid_BC();

	in_BC();

	ex_BC_crude();
}

void display(void)
{
	//int i, j, ip1, jp1, i0, icol, i1, i2, i3, i4, isol;
	//float minvar, maxvar, frac;
	//
	//minvar = 0.0;
	//maxvar = 0.2;
	//
	//stream();
	//apply_BCs();
	//collide();
	//computeVortConf(ni, vel_u, vel_v, dt, 0.55f);
	//
	//for (j = 0; j<nj; j++){
	//	for (i = 0; i<ni; i++){
	//		i0 = I2D(ni, i, j);
	//		//plotvar[i0] = sqrt(pow(vel_u[i0], 2) + pow(vel_v[i0], 2));
	//		frac = (plotvar[i0] - minvar) / (maxvar - minvar);
	//		icol = frac*ncol;
	//		isol = (int)solid[i0];
	//		plot_rgba[i0] = isol*cmap_rgba[icol];
	//	}
	//}
	//
	////glBufferData(GL_PIXEL_UNPACK_BUFFER_ARB, ni*nj*sizeof(unsigned int),
	////	(void **)plot_rgba, GL_STREAM_COPY);
	//
	////glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, ni, nj, GL_RGBA, GL_UNSIGNED_BYTE, 0);
	//
	//glClear(GL_COLOR_BUFFER_BIT);
	//glBegin(GL_POINTS);
	//glPointSize(10.0f);
	//glColor3f(1.0f, 1.0f, 1.0f);
	////glTexCoord2f(0.0, 0.0);
	////glVertex3f(0.0, 0.0, 0.0);
	////glTexCoord2f(1.0, 0.0);
	////glVertex3f(win_x, 0.0, 0.0);
	////glTexCoord2f(1.0, 1.0);
	////glVertex3f(win_x, win_y, 0.0);
	////glTexCoord2f(0.0, 1.0);
	////glVertex3f(0.0, win_y, 0.0);
	//glVertex2f(0.5f, 0.5f);
	//glEnd();
	//glutSwapBuffers();
	pre_display();
	stream();
	apply_BCs();
	collide();
	
	draw_vector_field(vel_u, vel_v, 1.0f, 0.0f, 0.0f, 0.0f);
	glEnd();
	post_display();
}

void resize(int w, int h)
{
	//width = w;
	//height = h;
	//glViewport(0, 0, w, h);
	//glMatrixMode(GL_PROJECTION);
	//glLoadIdentity();
	//glOrtho(0., win_x, 0., win_y, -200., 200.);
	//glMatrixMode(GL_MODELVIEW);
	//glLoadIdentity();

	glutReshapeWindow(w, h);
	win_x = w;
	win_y = h;
}

void mouse(int button, int state, int x, int y)
{
	float xx, yy;

	if ((button == GLUT_LEFT_BUTTON) && (state == GLUT_DOWN)) {
		draw_solid_flag = 0;
		xx = x;
		yy = y;
		ipos_old = xx / width*ni;
		jpos_old = (height - yy) / height*nj;
	}

	if ((button == GLUT_RIGHT_BUTTON) && (state == GLUT_DOWN)) {
		draw_solid_flag = 1;
		xx = x;
		yy = y;
		ipos_old = xx / width*ni;
		jpos_old = (height - yy) / height*nj;
	}
}

void mouse_motion(int x, int y)
{
	float xx, yy, frac;
	int ipos, jpos, i, j, i1, i2, j1, j2, jlast, jnext;
	xx = x;
	yy = y;
	ipos = (int)(xx / width*(float)ni);
	jpos = (int)((height - yy) / height*(float)nj);

	if (ipos <= ipos_old){
		i1 = ipos;
		i2 = ipos_old;
		j1 = jpos;
		j2 = jpos_old;
	}
	else {
		i1 = ipos_old;
		i2 = ipos;
		j1 = jpos_old;
		j2 = jpos;
	}

	jlast = j1;

	for (i = i1; i <= i2; i++){
		if (i1 != i2) {
			frac = (float)(i - i1) / (float)(i2 - i1);
			jnext = (int)(frac*(j2 - j1)) + j1;
		}
		else {
			jnext = j2;
		}
		if (jnext >= jlast) {
			solid[I2D(ni, i, jlast)] = draw_solid_flag;
			for (j = jlast; j <= jnext; j++){
				solid[I2D(ni, i, j)] = draw_solid_flag;
			}
		}
		else {
			solid[I2D(ni, i, jlast)] = draw_solid_flag;
			for (j = jnext; j <= jlast; j++){
				solid[I2D(ni, i, j)] = draw_solid_flag;
			}
		}
		jlast = jnext;
	}


	ipos_old = ipos;
	jpos_old = jpos;
}
int main(int argc, char **argv)
{
	int array_size_2d, totpoints, i;
	float rcol, gcol, bcol;

	FILE *fp_col;

	N = 128;
	ni = N;
	nj = N;
	dt = 0.01f;
	vxin = 0.15f;
	roout = 1.0f;
	tau = 0.51f;
	win_x = 640;
	win_y = 960;
	streamline_length = 20.0f;

	printf("ni = %d\n", ni);
	printf("nj = %d\n", nj);
	printf("vxin = %f\n", vxin);
	printf("roout = %f\n", roout);
	printf("tau = %f\n", tau);


	totpoints = (ni + 2) * (nj + 2);
	array_size_2d = (ni + 2) * (nj + 2)*sizeof(float);
	int array_size_2d_omega = (ni + 4) * (nj + 4) * sizeof(float);

	f0 = (float*)malloc(array_size_2d);
	f1 = (float*)malloc(array_size_2d);
	f2 = (float*)malloc(array_size_2d);
	f3 = (float*)malloc(array_size_2d);
	f4 = (float*)malloc(array_size_2d);
	f5 = (float*)malloc(array_size_2d);
	f6 = (float*)malloc(array_size_2d);
	f7 = (float*)malloc(array_size_2d);
	f8 = (float*)malloc(array_size_2d);

	tmpf0 = (float*)malloc(array_size_2d);
	tmpf1 = (float*)malloc(array_size_2d);
	tmpf2 = (float*)malloc(array_size_2d);
	tmpf3 = (float*)malloc(array_size_2d);
	tmpf4 = (float*)malloc(array_size_2d);
	tmpf5 = (float*)malloc(array_size_2d);
	tmpf6 = (float*)malloc(array_size_2d);
	tmpf7 = (float*)malloc(array_size_2d);
	tmpf8 = (float*)malloc(array_size_2d);

	vel_u = (float*)malloc(array_size_2d);
	vel_v = (float*)malloc(array_size_2d);
	w = (float*)malloc(array_size_2d_omega);

	plotvar = (float*)malloc(array_size_2d);

	plot_rgba = (unsigned int*)malloc((ni + 2) * (nj + 2)*sizeof(unsigned int));

	solid = (int*)malloc((ni + 2) * (nj + 2)*sizeof(int));

	faceq1 = 4.f / 9.f;
	faceq2 = 1.f / 9.f;
	faceq3 = 1.f / 36.f;

	float h = 1 / ni;
	float c = h / dt;
	float c2 = pow(c, 2);

	int emit_idx = I2D(ni, ni / 2, 2);
	for (i = 0; i<totpoints; i++) {
		if (i == emit_idx){
			f0[i] = faceq1 * roout * (1.f - 1.5f*vxin*vxin);
			f1[i] = faceq2 * roout * (1.f + 3.f*vxin + 4.5f*vxin*vxin- 1.5f*vxin*vxin);
			f2[i] = faceq2 * roout * (1.f - 1.5f*vxin*vxin);
			f3[i] = faceq2 * roout * (1.f - 3.f*vxin + 4.5f*vxin*vxin - 1.5f*vxin*vxin);
			f4[i] = faceq2 * roout * (1.f - 1.5f*vxin*vxin);
			f5[i] = faceq3 * roout * (1.f + 3.f*vxin + 4.5f*vxin*vxin - 1.5f*vxin*vxin);
			f6[i] = faceq3 * roout * (1.f - 3.f*vxin + 4.5f*vxin*vxin - 1.5f*vxin*vxin);
			f7[i] = faceq3 * roout * (1.f - 3.f*vxin + 4.5f*vxin*vxin - 1.5f*vxin*vxin);
			f8[i] = faceq3 * roout * (1.f + 3.f*vxin + 4.5f*vxin*vxin - 1.5f*vxin*vxin);
			plotvar[i] = vxin;
			solid[i] = 1;
		}
		else{
			f0[i] = faceq1 * roout * (1.f);
			f1[i] = faceq2 * roout * (1.f);
			f2[i] = faceq2 * roout * (1.f);
			f3[i] = faceq2 * roout * (1.f);
			f4[i] = faceq2 * roout * (1.f);
			f5[i] = faceq3 * roout * (1.f);
			f6[i] = faceq3 * roout * (1.f);
			f7[i] = faceq3 * roout * (1.f);
			f8[i] = faceq3 * roout * (1.f);
			//plotvar[i] = vxin;
			solid[i] = 1;
		}
	}

	//fopen_s(&fp_col, "cmap.dat", "r");
	//if (fp_col == NULL) {
	//	printf("Error: can't open cmap.dat \n");
	//	return 1;
	//}
	//
	//fscanf_s(fp_col, "%d", &ncol);
	//cmap_rgba = (unsigned int *)malloc(ncol*sizeof(unsigned int));
	//
	//for (i = 0; i<ncol; i++){
	//	fscanf_s(fp_col, "%f%f%f", &rcol, &gcol, &bcol);
	//	cmap_rgba[i] = ((int)(255.0f) << 24) |
	//		((int)(bcol * 255.0f) << 16) |
	//		((int)(gcol * 255.0f) << 8) |
	//		((int)(rcol * 255.0f) << 0);
	//}
	//fclose(fp_col);

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(win_x, win_y);
	glutInitWindowPosition((glutGet(GLUT_SCREEN_WIDTH) - win_x) / 2,
		(glutGet(GLUT_SCREEN_HEIGHT) - win_y) / 2);
	glutCreateWindow("2D LB");

	//printf("Loading extensions: %s\n", glewGetErrorString(glewInit()));
	//if (!glewIsSupported(
	//	"GL_VERSION_2_0 "
	//	"GL_ARB_pixel_buffer_object "
	//	"GL_EXT_framebuffer_object "
	//	)){
	//	fprintf(stderr, "ERROR: Support for necessary OpenGL extensions missing.");
	//	fflush(stderr);
	//	return -1;
	//}

	glClearColor(0.0, 0.0, 0.0, 0.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, win_x, 0., win_y, -200.0, 200.0);

	//glEnable(GL_TEXTURE_2D);
	//glGenTextures(1, &gl_Tex);
	//glBindTexture(GL_TEXTURE_2D, gl_Tex);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, ni, nj, 0,
	//	GL_RGBA, GL_UNSIGNED_BYTE, NULL);
	//
	//glGenBuffers(1, &gl_PBO);
	//glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, gl_PBO);
	//printf("Buffer created.\n");

	printf("Starting GLUT main loop...\n");
	glutDisplayFunc(display);
	glutReshapeFunc(resize);
	glutIdleFunc(display);
	glutMouseFunc(mouse);
	glutMotionFunc(mouse_motion);
	glutMainLoop();
	shutdown();
	return 0;
}
