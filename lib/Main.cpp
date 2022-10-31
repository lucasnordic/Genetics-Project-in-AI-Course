#include <sstream>  // for string streams
#include <string>	// for string

#include "Genetics.h"

//___________________________________________________
/*
* Variables
*/
int window_x = 1280;
int window_y = 920;
Genetics G{ window_x, window_y-200, 25};

//___________________________________________________
/*
* Function declarations
*/
void ChangeSize(int w, int h);	// for glutReshapeFunc
void initGFX();					// Initialize gfx before starting loop
void loop();					// Function to use when displaying gfx in glut window
void drawText(double x, double y, std::string text);// draw text
void drawOptionsInfoBar();		// on screen customizable options && info

//___________________________________________________
/*
* Main Method
*/
int main(int argc, char **argv){
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
void ChangeSize(int w, int h) {
	window_x = w; window_y = h;
	glViewport(0, 0, window_x, window_y);
	glMatrixMode(GL_PROJECTION); glLoadIdentity(); // GL_PROJECTION represents camera aperture, far-field, near-field...
	glOrtho(0., double(window_x), double(window_y), 0., -1., 1.);
	glMatrixMode(GL_MODELVIEW); glLoadIdentity(); // GL_MODELVIEW represents camera pos,dir,up/down...
}

void initGFX() {
	glClearColor(0.0, 0.0, 0.0, 1.0); // Black background
	glPointSize(10); // pixel size of a point
	glLineWidth(5);
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
	drawOptionsInfoBar();

	glutSwapBuffers(); // Swap buffers if GLUT_DOUBLE
}

void drawOptionsInfoBar() {
	// 1. Background
	glColor3ub(0, 0, 155);
	glLineWidth(100);
	//glBegin(GL_RECTAN);
	For(i, 2) {
		if (i == 0) { glVertex2dv(new double[] { 0, 820 }); }
		else { glVertex2dv(new double[] { (double)window_x, 820 }); }
	}
	glEnd();
	glLineWidth(5);

	glColor3ub(0, 127, 255);
	glBegin(GL_LINE_LOOP);
	For(i, 2) { 
		if (i == 0) { glVertex2dv(new double[] { 0, 720 }); }
		else		{ glVertex2dv(new double[] { (double)window_x, 720 }); }
	}
	glEnd();

	// 2. Info
	std::ostringstream gen1; gen1.precision(9); gen1 << "Best start fitness:    " << G.mBestStartChromLength;
	std::ostringstream gen0; gen0.precision(9); gen0 << "Generation:            " << G.mGenId;
	std::ostringstream gen3; gen3.precision(9); gen3 << "Best current fitness:  " << G.mChromDistances[0];
	std::ostringstream gen2; gen2.precision(4); gen2 << "% better than start:   " << 100 - (100 * G.mChromDistances[0] / G.mBestStartChromLength) << "%";
	
	std::ostringstream gen4;
	if (G.mGenId < G.maxGen)  gen4 << "Time Alive:            " << std::chrono::duration_cast<std::chrono::seconds>(G.mCurrTime);
	if (G.mGenId >= G.maxGen) gen4 << "Time until best chrom: " << std::chrono::duration_cast<std::chrono::seconds>(G.mTimeEnd);
	
	std::string running{};
	if (G.mRunning) { running = "True"; }
	else { running = "False"; }
	std::ostringstream gen5;					gen5 << "Is active:             " << running;

	drawText(100, 760, gen1.str());
	drawText(100, 800, gen0.str());
	drawText(100, 820, gen3.str());
	drawText(100, 840, gen2.str());
	drawText(100, 860, gen4.str());
	drawText(100, 880, gen5.str());

	// 3. Adjustable

}

void drawText(double x, double y, std::string text) {
	glColor3ub(215, 215, 0);
	//GFX_Text(mW - 250, mH - 45, "Generation: %d", mGenIdx);
	//double x = GoalFunc(0);
	//GFX_Text(mW - 250, mH - 25, "Best Path:  %3.1lf (%3.1lf%%)",
	//	x, 100. * x / mStartValueBestChrom);
	glRasterPos2d(x, y);
	For(i, text.length()) glutBitmapCharacter(GLUT_BITMAP_8_BY_13, text[i]);
}