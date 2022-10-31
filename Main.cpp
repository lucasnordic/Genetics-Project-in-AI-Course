#include "Genetics.h"

//___________________________________________________
/*
* Variables
*/
int window_x = 1280;
int window_y = 720;
Genetics G{ window_x, window_y, 25};

//___________________________________________________
/*
* Function declarations
*/
void initGFX();					// Initialize gfx before starting loop
void loop();					// Function to use when displaying gfx in glut window
void ChangeSize(int w, int h);	// for glutReshapeFunc

//___________________________________________________
/*
* Main Method
*/
int main(int argc, char **argv){
	// srand(time(NULL)); // random number generation // used with RandomBetween()

	// Initialize glut
	glutInit(&argc, argv);
	glutInitWindowSize(window_x, window_y);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutCreateWindow("Traveling Salesman Genetic Algorithm");
	glutReshapeFunc(ChangeSize); // window reshape callback // triggered immediatelly on first display and on resize...

	// Initialize gfx
	initGFX();

	// Initialize Gene Pool(s)
	G.InitCitiesXY();
	G.InitFirstGen();

	// set draw function
	glutDisplayFunc(loop);
	glutIdleFunc(loop);

	// Main Loop
	glutMainLoop();
    return 0;
}

//___________________________________________________
/*
* Functions
*/
void initGFX() {
	glClearColor(0.0, 0.0, 0.0, 1.0); // Black background
	glPointSize(8.0); // pixel size of a point
	glLineWidth(2);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
}

void loop() {
	// Do new/next generation
	G.NewGeneration();

	// Draw
	glClear(GL_COLOR_BUFFER_BIT);
	G.DrawCitiesPaths();
	glutSwapBuffers(); // Swap buffers if GLUT_DOUBLE
}

void ChangeSize(int w, int h) {
	window_x = w; window_y = h;
	glViewport(0, 0, window_x, window_y);
	glMatrixMode(GL_PROJECTION); glLoadIdentity(); // GL_PROJECTION represents camera aperture, far-field, near-field...
	glOrtho(0., double(window_x), double(window_y), 0., -1., 1.);
	glMatrixMode(GL_MODELVIEW); glLoadIdentity(); // GL_MODELVIEW represents camera pos,dir,up/down...
}