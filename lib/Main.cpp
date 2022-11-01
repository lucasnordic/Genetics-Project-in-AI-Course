#include <sstream>  // for string streams
#include <string>	// for string
#include <map>		// for Genetics
#include <utility>
#include <future>	// for async, future

#include "Genetics.h"

//___________________________________________________
/*
* Variables, Structs...
*/
int						window_x = 1280;	
int						window_y = 920;					
std::map <int,Genetics>	gArray{};			// array of Genetics
bool					leftBtn, rightBtn;	// user input
std::future<bool>		a3;		// Async func for G.DrawCities()
std::future<bool>		a4;		// Async func for G.DrawPaths()

//___________________________________________________
/*
* Function declarations
*/
void changeSize(int w, int h);			// for glutReshapeFunc
void initGFX();							// Initialize gfx before starting loop
void loop();							// Function to use when displaying gfx in glut window
//__
void doGenetics();						// Spawn new or remake(new/copy) Genetics 
void specialKeys(int key, int x, int y);// for user input
void drawOptionsInfoBar();				// bottom info bar
void drawGeneticsInfo();				// on screen customizable options && info
void drawText(double x, double y, std::string text, int r, int g, int b);// draw text


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

	// Callback registration
	glutReshapeFunc(changeSize); // window reshape callback // triggered on first display and on resize...
	glutDisplayFunc(loop);		 // Set draw function
	glutIdleFunc(loop);			 // Set draw function
	glutSpecialFunc(specialKeys);

	// Initialize gfx
	initGFX();									

	// Main Loop
	glutMainLoop();
    return 0;
}


//___________________________________________________
/*
* Functions
*/
void changeSize(int w, int h) {
	window_x = w; window_y = h;
	glViewport(0, 0, window_x, window_y);
	glMatrixMode(GL_PROJECTION); glLoadIdentity();	// GL_PROJECTION represents camera aperture, far-field, near-field...
	glOrtho(0., double(window_x), double(window_y), 0., -1., 1.);
	glMatrixMode(GL_MODELVIEW); glLoadIdentity();	// GL_MODELVIEW represents camera pos,dir,up/down...
}

void initGFX() {
	glClearColor(0.0, 0.0, 0.0, 1.0);	// Black background
	glPointSize(10);					// pixel size of a point
	glLineWidth(5);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
}

void loop() {
	bool arrEmpty;
	Genetics* pGen;

	doGenetics();						// Spawn new or remake(new/copy) Genetics 

	// Draw -----------------------------------------------------------
	glClear(GL_COLOR_BUFFER_BIT);

	drawOptionsInfoBar();				// 1. Draw bottom info bar

	arrEmpty = gArray.empty();			// Arr mpty?
	if (a3._Is_ready() && !arrEmpty) { 
		pGen = &gArray.find(0)->second; 
		pGen->DrawCities();				// 2. Draw - Cities
	}
	if (a4._Is_ready() && !arrEmpty) {
		pGen = &gArray.find(0)->second;
		pGen->NewGeneration();			// ------ Do new Gen ------
		pGen->DrawPaths();				// 3. Draw - Paths
	}

	drawGeneticsInfo();					// 4. Draw - Current Genetics info

	glutSwapBuffers();					// Swap buffers if GLUT_DOUBLE
	// end Draw --------------------------------------------------------
}

void doGenetics() {
	Genetics O;	
	Genetics* pGen;
	bool arrEmpty = gArray.size() == 0;

	if (arrEmpty) {
		gArray.emplace(0, Genetics(window_x, window_y - 200, 25));					// Construct new Genetics in array
		pGen = &gArray.find(0)->second;												// Create a pointer to it
		a3 = std::async([](Genetics* pGen) { return pGen->InitCitiesXY(); }, pGen);	// Async lambda, wait for InitCi...
		a4 = std::async([](Genetics* pGen) { return pGen->InitFirstGen(); }, pGen);	// Async lambda, wait for InitCi... complete
	} 
	else if (!arrEmpty && leftBtn) {
		O = gArray.find(0)->second;													// old genetics
		gArray.erase(0);															// erase old genetics
		gArray.emplace(0, Genetics(window_x, window_y - 200, 25));					// Construct new Genetics in array
		pGen = &gArray.find(0)->second;												// Create a pointer to it
		For(i, O.mMCI) For(j, 2) pGen->mCityXY[i][j] = O.mCityXY[i][j];				// copy city x/y to new G
		a4 = std::async([](Genetics* pGen) { return pGen->InitFirstGen(); }, pGen);	// Async Process, wait for InitCi... complete
		leftBtn = false;
	}
	else if (!arrEmpty && rightBtn) {
		gArray.erase(0);	// erase old genetics
		rightBtn = false;
	}
}

void specialKeys(int key, int x, int y) {
	bool a4Rdy = a4._Is_ready();

	if ((!leftBtn || !rightBtn) && !a4Rdy) return; // if btn(s) pressed or draw not done, return;

	if		(key == GLUT_KEY_LEFT)	{ leftBtn = true; } 
	else if (key == GLUT_KEY_RIGHT) { rightBtn = true; }
}

void drawGeneticsInfo() {
	Genetics G;
	std::ostringstream gen1; gen1 << "Best start fitness:    ";
	std::ostringstream gen0; gen0 << "Generation:            ";
	std::ostringstream gen3; gen3 << "Best current fitness:  ";
	std::ostringstream gen2; gen2 << "% better than start:   ";
	std::ostringstream gen4; gen4 << "Time alive:            ";
	std::ostringstream gen8; gen8 << "Time until best chrom: ";
	std::ostringstream gen5; gen5 << "Is active:             ";
	bool arrEmpty = gArray.size() == 0;
	int tlo = 50;	// text x offset
	int tyo = -10;	// text y offset
	int r = 215, g = 215, b = 0;

	if (arrEmpty) {
		gen1 << "?"; gen0 << "?"; gen3 << "?"; gen2 << "?"; gen4 << "?"; gen8 << "?"; gen5 << "?"; // not working
	}
	else if (!arrEmpty) {
		G = gArray.find(0)->second;
		double percentWorse = 100 * G.mChromDistances[0] / G.mBestStartChromLength;
		double percentBetter = 100 - percentWorse;

		// Add Current Genetics info
		gen1.precision(9); gen1 << G.mBestStartChromLength;
		gen0.precision(9); gen0 << G.mGenId;
		gen3.precision(9); gen3 << G.mChromDistances[0];
		gen2.precision(4); gen2 << percentBetter << "%";
		gen2.precision(9); gen4 << G.mCurrTime;
		gen2.precision(9); gen8 << G.mTimeEnd;
		std::string running{}; if (G.mRunning) { running = "True"; } else { running = "False"; }
		gen5 << running;

		// Conditional color formatting
		r = ((int)percentWorse) * 2.55;
		g = ((int)percentBetter) * 2.55; // scale to 255;
	}

	// Draw text
	drawText(tlo, 790+tyo, gen1.str(), 215, 215, 0); // best start
	drawText(tlo, 810+tyo, gen0.str(), 215, 215, 0); // gen
	drawText(tlo, 830+tyo, gen3.str(), r, g, 0); // best curr
	drawText(tlo, 850+tyo, gen2.str(), r, g, 0); // % better
	drawText(tlo, 870+tyo, gen4.str(), 215, 215, 0); // alive
	drawText(tlo, 890+tyo, gen8.str(), 215, 215, 0); // best born
	drawText(tlo, 910+tyo, gen5.str(), 215, 215, 0); // active gen
}

void drawOptionsInfoBar() {
	int tlo = 50;  // text left offset on x
	int tyo = -10; // text y offset
	std::ostringstream gen6; gen6 << "left arrow  <-  ==  start a new population with the same city x/y";
	std::ostringstream gen7; gen7 << "right arrow ->  ==  start a new population with new city x/y";

	// 1. Background
	// Horizontal line
	glColor3ub(0, 127, 255);
	glBegin(GL_LINE_LOOP);
	For(i, 2) {
		if (i == 0) { glVertex2dv(new double[] { 0, 720 }); }
		else { glVertex2dv(new double[] { (double)window_x, 720 }); }
	}
	glEnd();

	// Vertical line
	glBegin(GL_LINE_LOOP);
	For(i, 2) {
		if (i == 0) { glVertex2dv(new double[] { (double)window_x / 3, 720 }); }
		else { glVertex2dv(new double[] { (double)window_x / 3, (double)window_y }); }
	}
	glEnd();

	// 2. Info
	drawText(tlo, 760 + tyo, "INFO", 215, 215, 0);
	drawText(((double)window_x / 3) + tlo, 760 + tyo, "KEYBOARD COMMANDS", 215, 215, 0);
	drawText(((double)window_x / 3) + tlo, 790 + tyo, gen6.str(), 215, 215, 0);
	drawText(((double)window_x / 3) + tlo, 810 + tyo, gen7.str(), 215, 215, 0);

	// 3. Other info
	if (gArray.empty() || !a4._Is_ready()) {
		drawText(((double)window_x / 2) - 100, ((double)window_y - 200) / 2, "GENERATING NEW POPULATION", 215, 215, 0);
	}
}

void drawText(double x, double y, std::string text, int r, int g, int b) {
	if (r < 0 || g < 0 || b < 0) { r = 255; g = 255; b = 0; } // stop negative number and color as default

	glColor3ub(r, g, b);
	glRasterPos2d(x, y);
	For(i, text.length()) glutBitmapCharacter(GLUT_BITMAP_8_BY_13, text[i]);
}