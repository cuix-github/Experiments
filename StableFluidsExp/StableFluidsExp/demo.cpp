#include <stdlib.h>
#include <stdio.h>
#include <glut.h>
#include "helpers.h"
#include "solver.h"

static int N;
static int nx;
static int ny;
static float dt, diff, visc;
static float force, source;
static float vort_conf_coef;
static float temp;
static int dvel;

static float * u, *v, *u_prev, *v_prev;
static float * fx, *fy;
static float * psi, *wn, *dw, *w_bar, *w_star;
static float * du, *dv;
static float * dens, *dens_prev;
static float * t, * t0;
static Particle * particles;

static int win_id;
static int win_x, win_y;
static int mouse_down[3];
static int omx, omy, mx, my;
static bool pause = false;
static float streamline_length = 10.0f;
static int stop_frame = 0;
static int frame_counter = 0;
static int numParticles = 0;
float world_scale = 0;

static void free_data(void)
{
	if (particles) free(particles);
	if (u) free(u);
	if (v) free(v);
	if (du) free(du);
	if (dv) free(dv);
	if (u_prev) free(u_prev);
	if (v_prev) free(v_prev);
	if (psi) free(psi);
	if (wn) free(wn);
	if (dw) free(dw);
	if (w_bar) free(w_bar);
	if (w_star) free(w_star);
	if (dens) free(dens);
	if (dens_prev) free(dens_prev);
	if (fx) free(fx);
	if (fy) free(fy);
	if (t) free(t);
	if (t0) free(t0);
}

static void clear_data(void)
{
	int i, size = (N + 2)*(N + 2);

	for (i = 0; i<size; i++) {
		fx[i] = fy[i] =
			u[i] = v[i] = u_prev[i] = v_prev[i] =
			dens[i] = dens_prev[i] =
			psi[i] =
			du[i] = dv[i] =
			wn[i] = dw[i] = w_bar[i] = w_star[i] = t[i] = t0[i] = 0.0f;
	}

	float r, g, b;
	for (int i = 0; i != numParticles; i++){
		particles[i].x = (N / 2 + VFXEpoch::RandomI(-40, 40)) * world_scale;
		particles[i].y = (VFXEpoch::RandomI(0, 30)) * world_scale;
		r = VFXEpoch::RandomF(0.5, 0.7);
		g = VFXEpoch::RandomF(0.6, 0.8);
		b = VFXEpoch::RandomF(0.8, 1.0);
		particles[i].color = vec3(r, g, b);
	}

	frame_counter = 0;
}

static int allocate_data(void)
{
	int size = (N + 2)*(N + 2);

	fx = (float *)malloc(size*sizeof(float));
	fy = (float *)malloc(size*sizeof(float));
	u = (float *)malloc(size*sizeof(float));
	v = (float *)malloc(size*sizeof(float));
	du = (float *)malloc(size*sizeof(float));
	dv = (float *)malloc(size*sizeof(float));
	u_prev = (float *)malloc(size*sizeof(float));
	v_prev = (float *)malloc(size*sizeof(float));
	psi = (float *)malloc(size*sizeof(float));
	wn = (float *)malloc(size*sizeof(float));
	dw = (float *)malloc(size*sizeof(float));
	w_bar = (float *)malloc(size*sizeof(float));
	w_star = (float *)malloc(size*sizeof(float));
	dens = (float *)malloc(size*sizeof(float));
	dens_prev = (float *)malloc(size*sizeof(float));
	t = (float *)malloc(size*sizeof(float));
	t0 = (float *)malloc(size*sizeof(float));
	particles = (Particle*)malloc(numParticles*sizeof(Particle));

	if (!fx || !fy ||
		!psi ||
		!wn || !dw || !w_bar || !w_star ||
		!u || !v || !du || !dv || !u_prev || !v_prev || !t || !t0 ||
		!dens || !dens_prev) {
		fprintf(stderr, "cannot allocate data\n");
		return (0);
	}

	// Initialize particle start position
	float r, g, b;
	for (int i = 0; i != numParticles; i++){
		particles[i].x = (N / 2 + VFXEpoch::RandomI(-40, 40)) * world_scale;
		particles[i].y = (VFXEpoch::RandomI(0, 30)) * world_scale;
		r = VFXEpoch::RandomF(0.5, 0.7);
		g = VFXEpoch::RandomF(0.6, 0.8);
		b = VFXEpoch::RandomF(0.8, 1.0);
		particles[i].color = vec3(r, g, b);
	}

	return (1);
}

static void pre_display(void)
{
	glViewport(0, 0, win_x, win_y);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, 1.0, 0.0, 1.0);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// Make the pixel looks round.
	glEnable(GL_POINT_SMOOTH);
}

static void post_display(void)
{
	glutSwapBuffers();
}

static void draw_vector_field(float * u, float * v, float lineWidth, float r, float g, float b)
{
	int i, j;
	float x, y, h;

	h = 1.0f / N;

	glColor3f(r, g, b);
	glLineWidth(lineWidth);

	glBegin(GL_LINES);

	for (i = 1; i <= N; i++) {
		x = (i - 0.5f)*h;
		for (j = 1; j <= N; j++) {
			y = (j - 0.5f)*h;

			glVertex2f(x, y);
			glVertex2f(x + u[j * (N + 2) + i] * streamline_length / N,
				y + v[j * (N + 2) + i] * streamline_length / N);
		}
	}

	glEnd();
}

static void draw_particles(float * u, float * v, float pointSize)
{
	glPointSize(pointSize);

	glBegin(GL_POINTS);
	for (int i = 0; i != numParticles; i++){
		if (particles[i].x < 0 || particles[i].x> 1 ||
			particles[i].y < 0 || particles[i].y> 1){
			particles[i].x = (N / 2 + VFXEpoch::RandomI(-40, 40)) * world_scale;
			particles[i].y = (VFXEpoch::RandomI(0, 30)) * world_scale;
		}
		glColor3f(particles[i].color.x, particles[i].color.y, particles[i].color.z);
		glVertex2f(particles[i].x, particles[i].y);
	}
	glEnd();
}

static void draw_scalar_field(float * field, int r, int g, int b)
{
	int i, j;
	float x, y, h, d00, d01, d10, d11;

	h = 1.0f / N;

	glBegin(GL_QUADS);

	for (i = 0; i <= N; i++) {
		x = (i - 0.5f)*h;
		for (j = 0; j <= N; j++) {
			y = (j - 0.5f)*h;

			d00 = dens[IX(i, j)];
			d01 = dens[IX(i, j + 1)];
			d10 = dens[IX(i + 1, j)];
			d11 = dens[IX(i + 1, j + 1)];

			glColor3f(d00, d00, d00); glVertex2f(x, y);
			glColor3f(d10, d10, d10); glVertex2f(x + h, y);
			glColor3f(d11, d11, d11); glVertex2f(x + h, y + h);
			glColor3f(d01, d01, d01); glVertex2f(x, y + h);
		}
	}

	glEnd();
}

static void get_from_UI(float * d, float * u, float * v)
{
	int i, j, size = (N + 2)*(N + 2);

	for (i = 0; i<size; i++) {
		u[i] = v[i] = d[i] = 0.0f;
	}

	if (!mouse_down[0] && !mouse_down[2]) return;

	i = (int)((mx / (float)win_x)*N + 1);
	j = (int)(((win_y - my) / (float)win_y)*N + 1);

	if (i<1 || i>N || j<1 || j>N) return;

	if (mouse_down[0]) {
		u[IX(i, j)] = force * (mx - omx);
		v[IX(i, j)] = force * (omy - my);
	}

	if (mouse_down[2]) {
		d[IX(i, j)] = source;
	}

	omx = mx;
	omy = my;

	return;
}

static void key_func(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'c':
	case 'C':
		clear_data();
		break;

	case 'q':
	case 'Q':
		free_data();
		exit(0);
		break;

	case 'v':
	case 'V':
		dvel = !dvel;
		break;

	case ' ':
		pause = !pause;
		break;
	}
}

static void mouse_func(int button, int state, int x, int y)
{
	omx = mx = x;
	omx = my = y;

	mouse_down[button] = state == GLUT_DOWN;
}

static void motion_func(int x, int y)
{
	mx = x;
	my = y;
}

static void reshape_func(int width, int height)
{
	glutSetWindow(win_id);
	glutReshapeWindow(width, height);

	win_x = width;
	win_y = height;
}

static void idle_func(void)
{
	get_from_UI(dens_prev, u_prev, v_prev);
	int idxX = N / 2;
	int idxY = 3;

	v_prev[IX(idxX, idxY)] = force;
	t0[IX(idxX, idxY)] = temp;
	t0[IX(idxX + 1, idxY)] = temp;
	t0[IX(idxX - 1, idxY)] = temp;
	t0[IX(idxX + 2, idxY)] = temp;
	t0[IX(idxX - 2, idxY)] = temp;
	dens_prev[IX(idxX, idxY)] = source;

	if (!pause){
		if (frame_counter != stop_frame)
		{
			IVOCKAdvance(N, particles, numParticles, fx, fy, psi, du, dv, wn, dw, w_bar, w_star, u, v, u_prev, v_prev, t, t0, visc, dt);
			computeBuoyancy(N, v, dens, t, 0.1f, 0.3f, dt);
			
			// TODO: Fix bugs
			computeVortConf(N, u, v, dt, vort_conf_coef);
			project(N, u, v, u_prev, v_prev);
			MoveScalarProperties(N, t, t0, u, v, 0.0f, dt);
			MoveScalarProperties(N, dens, dens_prev, u, v, diff, dt);
			frame_counter++;
		}
	}
	glutSetWindow(win_id);
	glutPostRedisplay();
}

static void display_func(void)
{
	if (!pause){
		pre_display();
		//draw_scalar_field(dens, 1.0f, 1.0f, 1.0f);
		//draw_scalar_field(t, 1.0f, 1.0f, 1.0f);
		draw_vector_field(u, v, 1.0, 0.0f, 1.0f, 0.0f);
		//draw_vector_field(du, dv, 1.0f, 0.0f, 1.0f, 0.0f);
		//draw_particles(u, v, 1.0f);
		post_display();
	}
}

static void open_glut_window(void)
{
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);

	glutInitWindowPosition((glutGet(GLUT_SCREEN_WIDTH) - win_x) / 2,
		(glutGet(GLUT_SCREEN_HEIGHT) - win_y) / 2);
	glutInitWindowSize(win_x, win_y);
	win_id = glutCreateWindow("Smoke");

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers();
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers();

	pre_display();

	glutKeyboardFunc(key_func);
	glutMouseFunc(mouse_func);
	glutMotionFunc(motion_func);
	glutReshapeFunc(reshape_func);
	glutIdleFunc(idle_func);
	glutDisplayFunc(display_func);
}

int main(int argc, char ** argv)
{
	N = 192;
	dt = 0.01f;
	diff = 0.0f;
	visc = 0.0f;
	force = 100.0f;
	source = 100.0f;
	temp = 50.0f;
	stop_frame = -1;
	numParticles = 10000;
	world_scale = 1.0 / N;
	vort_conf_coef = 0.15f;
	streamline_length = 5.0f;
	cout << "Default values of the simualtion: " << endl;
	cout << "Dim = " << N << " x " << N << endl;
	cout << "Time step = " << dt << endl;
	cout << "Diffuse = " << diff << endl;
	cout << "Viscosity = " << visc << endl;
	cout << "Force = " << force << endl;
	cout << "Source = " << source << endl;
	cout << "Vorticity Control coefficient = " << vort_conf_coef << endl;
	cout << "Number of Particles = " << numParticles << endl;

	if (!allocate_data()) exit(1);
	clear_data();

	win_x = 640;
	win_y = 720;
	open_glut_window();

	glutMainLoop();

	exit(0);
}