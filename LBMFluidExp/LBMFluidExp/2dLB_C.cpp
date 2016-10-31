#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <glew.h>
#include <glut.h>
#include "Helpers.h"

#pragma comment(lib, "glew32.lib")
#pragma comment(lib, "glew32s.lib")



#define I2D(ni,i,j) (((ni)*(j)) + i)

GLuint gl_PBO, gl_Tex;

int win_x = 960;
int win_y = 640;

float *f0,*f1,*f2,*f3,*f4,*f5,*f6,*f7,*f8;
float *tmpf0,*tmpf1,*tmpf2,*tmpf3,*tmpf4,*tmpf5,*tmpf6,*tmpf7,*tmpf8;
float *vel_u, *vel_v;
float *cmap,*plotvar;
int *solid;
unsigned int *cmap_rgba, *plot_rgba;  //rgba arrays for plotting

float tau,faceq1,faceq2,faceq3; 
float vxin, roout;
float width, height;
int ni,nj;
int ncol;
int streamline_length;
int ipos_old,jpos_old, draw_solid_flag;
int win_id;

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

static void pre_display(void)
{
	glViewport(0, 0, win_x, win_y);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, 1.0, 0.0, 1.0);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// Make the pixel looks round.
	// glEnable(GL_POINT_SMOOTH);
}


void stream(void)
{
    int i,j,im1,ip1,jm1,jp1,i0;
    for (j=0; j<nj; j++) {
	jm1=j-1;
	jp1=j+1;
	if (j==0) jm1=0;
	if (j==(nj-1)) jp1=nj-1;
	for (i=1; i<ni; i++) {
	    im1 = i-1;
	    ip1 = i+1;
	    if (i==0) im1=0;
	    if (i==(ni-1)) ip1=ni-1;
		i0 = I2D(ni, i, j);
	    tmpf1[i0] = f1[I2D(ni,im1,j)];
	    tmpf2[i0] = f2[I2D(ni,i,jp1)];
	    tmpf3[i0] = f3[I2D(ni,ip1,j)];
	    tmpf4[i0] = f4[I2D(ni,i,jm1)];
		tmpf5[i0] = f5[I2D(ni,im1,jp1)];
		tmpf6[i0] = f6[I2D(ni,ip1,jp1)];
		tmpf7[i0] = f7[I2D(ni,ip1,jm1)];
		tmpf8[i0] = f8[I2D(ni,im1,jm1)];
	}
    }

    for (j=0; j<nj; j++) {
	for (i=1; i<ni; i++) {
	    i0 = I2D(ni,i,j);
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
    int i,j,i0;
    float ro, rovx, rovy, vx(0.0f), vy(0.0f), v_sq_term;
    float f0eq, f1eq, f2eq, f3eq, f4eq, f5eq, f6eq, f7eq, f8eq;
    float rtau, rtau1;

    rtau = 1.f/tau;
    rtau1 = 1.f - rtau;

    for (j=0; j<nj; j++) {
	for (i=0; i<ni; i++) {

	    i0 = I2D(ni,i,j);

	    // Do the summations needed to evaluate the density and components of velocity
	    ro = f0[i0] + f1[i0] + f2[i0] + f3[i0] + f4[i0] + f5[i0] + f6[i0] + f7[i0] + f8[i0];
	    rovx = f1[i0] - f3[i0] + f5[i0] - f6[i0] - f7[i0] + f8[i0];
	    rovy = f2[i0] - f4[i0] + f5[i0] + f6[i0] - f7[i0] - f8[i0];
	    vx = rovx/ro;
	    vy = rovy/ro;

		vel_u[I2D(ni, i, j)] = vx;
		vel_v[I2D(ni, i, j)] = vy;

	    // Also load the velocity magnitude into plotvar - this is what we will
	    // display using OpenGL later
	    plotvar[i0] = sqrt(vx*vx + vy*vy);

	    v_sq_term = 1.5f*(vx*vx + vy*vy);

	    // Evaluate the local equilibrium f values in all directions
	    f0eq = ro * faceq1 * (1.f - v_sq_term);
	    f1eq = ro * faceq2 * (1.f + 3.f*vx + 4.5f*vx*vx - v_sq_term);
	    f2eq = ro * faceq2 * (1.f + 3.f*vy + 4.5f*vy*vy - v_sq_term);
	    f3eq = ro * faceq2 * (1.f - 3.f*vx + 4.5f*vx*vx - v_sq_term);
	    f4eq = ro * faceq2 * (1.f - 3.f*vy + 4.5f*vy*vy - v_sq_term);
	    f5eq = ro * faceq3 * (1.f + 3.f*(vx + vy) + 4.5f*(vx + vy)*(vx + vy) - v_sq_term);
	    f6eq = ro * faceq3 * (1.f + 3.f*(-vx + vy) + 4.5f*(-vx + vy)*(-vx + vy) - v_sq_term);
	    f7eq = ro * faceq3 * (1.f + 3.f*(-vx - vy) + 4.5f*(-vx - vy)*(-vx - vy) - v_sq_term);
	    f8eq = ro * faceq3 * (1.f + 3.f*(vx - vy) + 4.5f*(vx - vy)*(vx - vy) - v_sq_term);

	    // Simulate collisions by "relaxing" toward the local equilibrium
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
    int i,j,i0;
    float f1old,f2old,f3old,f4old,f5old,f6old,f7old,f8old;
    
    for (j=0;j<nj;j++){
	for (i=0;i<ni;i++){
	    i0=I2D(ni,i,j);
	    if (solid[i0]==0) {
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
    int i0,i1,i;

    for (i=0; i<ni; i++){
	i0 = I2D(ni,i,0);
	i1 = I2D(ni,i,nj-1);
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
    float f1new, f5new, f8new, vx_term;

    vx_term = 1.f + 3.f*vxin +3.f*vxin*vxin;
    f1new = roout * faceq2 * vx_term;
    f5new = roout * faceq3 * vx_term;
    f8new = f5new;

    for (j=0; j<nj; j++){
      i0 = I2D(ni,0,j);
      f1[i0] = f1new;
      f5[i0] = f5new;
      f8[i0] = f8new;
    }

}

void ex_BC_crude(void)
{
    int i0, i1, j;

    for (j=0; j<nj; j++){
	i0 = I2D(ni,ni-1,j);
	i1 = i0 - 1;
	f3[i0] = f3[i1];
	f6[i0] = f6[i1];
	f7[i0] = f7[i1];
    }
}

void apply_BCs(void)
{
    per_BC();

    solid_BC();
	 	
    in_BC();

    ex_BC_crude();
}

static void draw_vector_field(float * u, float * v, float lineWidth, float r, float g, float b)
{
	int i, j;
	float x, y, h;

	h = 1.0f / ni;

	glColor3f(r, g, b);
	glLineWidth(lineWidth);

	glBegin(GL_LINES);

	for (i = 1; i <= nj; i++) {
		x = (i - 0.5f)*h;
		for (j = 1; j <= ni; j++) {
			y = (j - 0.5f)*h;

			glVertex2f(x, y);
			glVertex2f(x + u[j * (ni + 2) + i] * streamline_length / ni,
				y + v[j * (ni + 2) + i] * streamline_length / ni);
		}
	}

	glEnd();
}

void display(void)
{
    int i,j,ip1,jp1,i0,icol,i1,i2,i3,i4,isol;
    float minvar,maxvar,frac;
	float x, y, h, d00, d01, d10, d11;

    minvar=0.0;
    maxvar=0.2;

	//stream();
	//apply_BCs();
    //collide();

	pre_display();

	glColor3f(0.0f, 1.0f, 0.0f);
	glPointSize(10.0f);
	glBegin(GL_POINTS);
	glVertex2f(5, 5);
	glEnd();
	glutSwapBuffers();
}

void resize(int w, int h)
{
	glutSetWindow(win_id);
	glutReshapeWindow(w, h);

	win_x = w;
	win_y = h;
}
    
void mouse(int button, int state, int x, int y)
{
    float xx,yy;

    if ((button == GLUT_LEFT_BUTTON) && (state == GLUT_DOWN)) {
        draw_solid_flag = 0;
        xx=x;
        yy=y;
        ipos_old=xx/width*ni;
        jpos_old=(height-yy)/height*nj;
    }

    if ((button == GLUT_RIGHT_BUTTON) && (state == GLUT_DOWN)) {
        draw_solid_flag = 1;
        xx=x;
        yy=y;
        ipos_old=xx/width*ni;
        jpos_old=(height-yy)/height*nj;
    }
}

void mouse_motion(int x, int y)
{
    float xx,yy,frac;
    int ipos,jpos,i,j,i1,i2,j1,j2, jlast, jnext;
    xx=x;
    yy=y;
    ipos=(int)(xx/width*(float)ni);
    jpos=(int)((height-yy)/height*(float)nj);

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
    
    jlast=j1;

    for (i=i1;i<=i2;i++){
        if (i1 != i2) {
            frac=(float)(i-i1)/(float)(i2-i1);
            jnext=(int)(frac*(j2-j1))+j1;
        }
        else {
            jnext=j2;
        }
        if (jnext >= jlast) {
            solid[I2D(ni,i,jlast)]=draw_solid_flag;
            for (j=jlast; j<=jnext; j++){
                solid[I2D(ni,i,j)]=draw_solid_flag;
            }
        }
        else {
            solid[I2D(ni,i,jlast)]=draw_solid_flag;
            for (j=jnext; j<=jlast; j++){
                solid[I2D(ni,i,j)]=draw_solid_flag;
            }
        }
        jlast = jnext;
    }

    
    ipos_old=ipos;
    jpos_old=jpos;
}

void shutdown()
{
	//float *f0, *f1, *f2, *f3, *f4, *f5, *f6, *f7, *f8;
	//float *tmpf0, *tmpf1, *tmpf2, *tmpf3, *tmpf4, *tmpf5, *tmpf6, *tmpf7, *tmpf8;
	//float *vel_u, *vel_v;
	//float *cmap, *plotvar;
	//int *solid;
	//unsigned int *cmap_rgba, *plot_rgba;  //rgba arrays for plotting

	if (!f0) free(f0);
	if (!f1) free(f1);
	if (!f2) free(f2);
	if (!f3) free(f3);
	if (!f4) free(f4);
	if (!f5) free(f5);
	if (!f6) free(f6);
	if (!f7) free(f7);
	if (!f8) free(f8);

	if (!tmpf0) free(tmpf0);
	if (!tmpf1) free(tmpf1);
	if (!tmpf2) free(tmpf2);
	if (!tmpf3) free(tmpf3);
	if (!tmpf4) free(tmpf4);
	if (!tmpf5) free(tmpf5);
	if (!tmpf6) free(tmpf6);
	if (!tmpf7) free(tmpf7);
	if (!tmpf8) free(tmpf8);

	if (!vel_u) free(vel_u);
	if (!vel_v) free(vel_v);

	if (!cmap) free(cmap);
	if (!plotvar) free(plotvar);
	if (!solid) free(solid);
	if (!cmap_rgba) free(cmap_rgba);
	if (!plot_rgba) free(plot_rgba);
}

static void idle_func(void){
	glutSetWindow(win_id);
	glutPostRedisplay();
}

int main(int argc, char **argv)
{
	int array_size_2d, totpoints, i;
	float rcol, gcol, bcol;

	// The following parameters are usually read from a file, but
	// hard code them for the demo:
	ni = 32;
	nj = 32;
	vxin = 0.04;
	roout = 1.0;
	tau = 0.51;
	streamline_length = 10.0f;
	// End of parameter list

	// Write parameters to screen
	printf("ni = %d\n", ni);
	printf("nj = %d\n", nj);
	printf("vxin = %f\n", vxin);
	printf("roout = %f\n", roout);
	printf("tau = %f\n", tau);


	totpoints = (ni + 2)*(nj + 2);
	array_size_2d = (ni + 2)*(nj + 2)*sizeof(float);

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

	plotvar = (float*)malloc(array_size_2d);

	plot_rgba = (unsigned int*)malloc(ni*nj*sizeof(unsigned int));

	solid = (int*)malloc(ni*nj*sizeof(int));

	vel_u = (float*)malloc(array_size_2d);
	vel_v = (float*)malloc(array_size_2d);

	faceq1 = 4.f / 9.f;
	faceq2 = 1.f / 9.f;
	faceq3 = 1.f / 36.f;
	for (i = 0; i<totpoints; i++) {
		f0[i] = faceq1 * roout * (1.f - 1.5f*vxin*vxin);
		f1[i] = faceq2 * roout * (1.f + 3.f*vxin + 4.5f*vxin*vxin - 1.5f*vxin*vxin);
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


	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);

	glutInitWindowPosition((glutGet(GLUT_SCREEN_WIDTH) - win_x) / 2,
		(glutGet(GLUT_SCREEN_HEIGHT) - win_y) / 2);
	glutInitWindowSize(win_x, win_y);
	win_id = glutCreateWindow("LBM Sim");

	glutMouseFunc(mouse);
	glutMotionFunc(mouse_motion);
	glutReshapeFunc(resize);
	glutIdleFunc(idle_func);
	glutDisplayFunc(display);
	glutMainLoop();
	shutdown();
	return 0;
}


