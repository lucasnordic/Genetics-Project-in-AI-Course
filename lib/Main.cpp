#include <sstream>  // for string streams
#include <string>	// for string
#include <map>		// for Genetics
#include <utility>
#include <future>	// for async, future
#include <chrono>	// for millisecond
#include <vector>	// for runInfo, list of runs in percentage
#include <numeric>	// for reduce()

#include "Genetics.h"

//___________________________________________________
/*
* Variables, Structs...
*/
bool					leftBtn, rightBtn, upBtn, downBtn; // user input
int						window_x = 1280;// 16:9
int						window_y = 920;	// 720 == 16:9 == no bottom bar
std::map <int,Genetics>	gArray{};		// array of a single Genetics
std::future<bool>		a3;				// Async func for G.DrawCities()
std::future<bool>		a4;				// Async func for G.DrawPaths()

struct {								
	int					amountRuns = 0;
	bool				testRun = true;
	double				avPercImprov = 0;
	std::vector<double>	listOfPerc{};
} runInfo;				// information of a test run of Genetics
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
void test30Cities();
float averageOfVect(std::vector<double> const& v);

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
* Main Functions
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

	doGenetics();						// 0. Spawn new or remake(new/copy) Genetics 

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

	if (arrEmpty) {										// 1. NEW GENETICS, NEW CITIES
		gArray.emplace(0, 
			Genetics(window_x, window_y - 200, 25));	// Construct new Genetics in array
		pGen = &gArray.find(0)->second;					// Create a pointer to it
		
		a3 = std::async([](Genetics* pGen) { 			// Async lambda, wait for InitCi...
			return pGen->InitCitiesXY(); }, pGen);
		a4 = std::async([](Genetics* pGen) {			// Async lambda, wait for InitCi... complete
			return pGen->InitPopulation(pGen->mch,true); }, pGen);

		runInfo.amountRuns++;
		runInfo.testRun = false;
	} 
	else if (!arrEmpty && leftBtn) {					// 2. NEW GENETICS, COPY CITY
		// 2.1 keep old, erase old
		O = gArray.find(0)->second;						// old genetics
		gArray.erase(0);								// erase old genetics

		// 2.2 Place new Genetics
		gArray.emplace(0, Genetics(window_x, window_y - 200, 25));	
		pGen = &gArray.find(0)->second;								
		For(i, O.mci) For(j, 2) {
			pGen->cityXY[i][j] = O.cityXY[i][j];		// copy city x/y to new G
		}

		a4 = std::async([](Genetics* pGen) {						
			return pGen->InitPopulation(pGen->mch, true);}, pGen);	

		// 2.4 calculate percentage improvement, add to list of percentages, calculate averages of runs
		if (runInfo.amountRuns > 0) {
			double percentBetter = 100 - (100 * O.chromDistances[0] / O.bestStartChromLength);
			runInfo.listOfPerc.push_back(percentBetter);
			if (runInfo.listOfPerc.size() > 0) runInfo.avPercImprov = averageOfVect(runInfo.listOfPerc);
		}

		leftBtn = false;
		runInfo.amountRuns++;
	}
	else if (!arrEmpty && rightBtn && !runInfo.testRun) {// 3. DELETE GENETICS AND CITIES
		gArray.erase(0);

		rightBtn = false;
		runInfo.amountRuns++;
	}
	else if (!arrEmpty && upBtn) {						// 4. MUTATE POPULATION
		pGen = &gArray.find(0)->second;
		int chroms = pGen->mch;

		for (int i = 0; i < chroms; i++) {				// for chroms
			for (int j = 0; j < pGen->mci/4; j++) {		// for half the amount of genes in chrom
				pGen->Mutate(i);				
			}
		}

		upBtn = false;
	}
	else if (!arrEmpty && downBtn) {					// 5. Do test run 
		// reset run info, set test run to true
		runInfo.amountRuns = 0;		
		runInfo.avPercImprov = 0.;	
		runInfo.listOfPerc.clear();	
		runInfo.testRun = true;

		leftBtn = true;
		downBtn = false;
	}

	if (a4._Is_ready() && runInfo.testRun == true) {	// 5.1 Start new test
		pGen = &gArray.find(0)->second;
		if (pGen->genId >= pGen->maxGen) {leftBtn = true;}

		test30Cities();	// Test preset city x/y;
	}
}

//___________________________________________________
/*
* Functions
*/
void specialKeys(int key, int x, int y) {
	bool a4Rdy = a4._Is_ready();

	if ((!leftBtn || !rightBtn) && !a4Rdy) return; // if btn(s) pressed or draw not done, return;

	if		(key == GLUT_KEY_LEFT)	{ leftBtn = true; } 
	else if (key == GLUT_KEY_RIGHT) { rightBtn = true; }
	else if (key == GLUT_KEY_UP) { upBtn = true; }
	else if (key == GLUT_KEY_DOWN) { 
		if (!runInfo.testRun)downBtn = true;
		runInfo.testRun = false; 
	}
}

void drawGeneticsInfo() {
	Genetics G;
	std::ostringstream gen1; gen1	<< "Best start fitness:    ";
	std::ostringstream gen0; gen0	<< "Generation:            ";
	std::ostringstream gen3; gen3	<< "Best current fitness:  ";
	std::ostringstream gen2; gen2	<< "% better than start:   ";
	std::ostringstream gen4; gen4	<< "Time alive:            ";
	std::ostringstream gen8; gen8	<< "Time until best chrom: ";
	std::ostringstream gen5; gen5	<< "Average% and TestRuns: ";
	bool arrEmpty = gArray.size() == 0;
	int tlo = 50;	// text x offset
	int tyo = -10;	// text y offset
	int r = 215, g = 215, b = 0;

	if (arrEmpty) {
		gen1 << "0"; gen0 << "0"; gen3 << "0"; gen2 << "0"; gen4 << "0"; gen8 << "0"; gen5 << "0"; // not working
	}
	else if (!arrEmpty) {
		G = gArray.find(0)->second;
		double percentWorse = 100 * G.chromDistances[0] / G.bestStartChromLength;
		double percentBetter = 100 - percentWorse;

		// Add Current Genetics info
		gen1.precision(9); gen1 << G.bestStartChromLength;
		gen0.precision(9); gen0 << G.genId;
		gen3.precision(9); gen3 << G.chromDistances[0];
		gen2.precision(4); gen2 << percentBetter << "%";
		gen2.precision(9); gen4 << G.timeCurr.count();
		gen2.precision(9); gen8 << G.timeEnd.count();
		gen5.precision(4); gen5 << runInfo.avPercImprov << "%   Run: " << runInfo.amountRuns;

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
	std::ostringstream gen6; gen6 << "left  arrow  <-  ==  start a new population with the same city x/y";
	std::ostringstream gen7; gen7 << "right arrow  ->  ==  start a new population with new city x/y";
	std::ostringstream gen8; gen8 << "up    arrow      ==  chernobyl mode || mutate the population x times";
	std::ostringstream gen9; 
	if (!gArray.empty()) {
		gen9 << "down  arrow      ==  load 30 preset cities && restart after " 
			<< gArray.find(0)->second.maxGen << " gen; !!M_CI == 30!!";
	}
	else {
		gen9 << "down  arrow      ==  load 30 preset cities && restart after x gen; !!M_CI == 30!!";
	}

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
	drawText(((double)window_x / 3) + tlo, 830 + tyo, gen8.str(), 215, 115, 0);
	drawText(((double)window_x / 3) + tlo, 850 + tyo, gen9.str(), 215, 115, 100);

	// 3. Other info
	if (gArray.empty() || !a4._Is_ready()) {
		drawText(((double)window_x / 2) - 100, ((double)window_y - 200) / 2, 
			"GENERATING NEW POPULATION", 215, 215, 0);
	}
	if (runInfo.testRun ) {
		drawText(((double)window_x / 2) - 60,
			((double)window_y - 240) / 2,
			"TEST RUN ACTIVE",
			215, 115, 100);
	}
}

void drawText(double x, double y, std::string text, int r, int g, int b) {
	if (r < 0 || g < 0 || b < 0) { r = 255; g = 255; b = 0; } // stop negative number and color as default

	glColor3ub(r, g, b);
	glRasterPos2d(x, y);
	For(i, text.length()) glutBitmapCharacter(GLUT_BITMAP_8_BY_13, text[i]);
}

float averageOfVect(std::vector<double> const& v) {
	if (v.empty()) {
		return 0;
	}

	auto const count = static_cast<double>(v.size());
	return std::reduce(v.begin(), v.end()) / count;
}

void test30Cities() {
	Genetics *pGen = &gArray.find(0)->second;
	pGen->cityXY [0][0] = {544.00000000000000};   pGen->cityXY[0][1] = {202.00000000000000};
	pGen->cityXY [1][0] = {69.000000000000000};    pGen->cityXY[1][1] = { 625.00000000000000 };
	pGen->cityXY [2][0] = { 666.00000000000000 };  pGen->cityXY[2][1] = { 36.000000000000000 };
	pGen->cityXY [3][0] = { 474.00000000000000 };  pGen->cityXY[3][1] = { 122.00000000000000 };
	pGen->cityXY [4][0] = { 356.00000000000000 }; pGen->cityXY [4][1] = { 510.00000000000000 };
	pGen->cityXY [5][0] = { 1226.0000000000000 }; pGen->cityXY [5][1] = { 336.00000000000000 };
	pGen->cityXY [6][0] = { 576.00000000000000 }; pGen->cityXY [6][1] = { 390.00000000000000 };
	pGen->cityXY [7][0] = { 325.00000000000000 }; pGen->cityXY [7][1] = { 590.00000000000000 };
	pGen->cityXY [8][0] = { 568.00000000000000 }; pGen->cityXY [8][1] = { 110.00000000000000 };
	pGen->cityXY [9][0] = { 1185.0000000000000 }; pGen->cityXY [9][1] = { 79.000000000000000 };
	pGen->cityXY[10][0] = { 68.000000000000000 }; pGen->cityXY[10][1] = { 407.00000000000000 };
	pGen->cityXY[11][0] = { 187.00000000000000 }; pGen->cityXY[11][1] = { 183.00000000000000 };
	pGen->cityXY[12][0] = { 1103.0000000000000 }; pGen->cityXY[12][1] = { 430.00000000000000 };
	pGen->cityXY[13][0] = { 663.00000000000000 }; pGen->cityXY[13][1] = { 428.00000000000000 };
	pGen->cityXY[14][0] = { 780.00000000000000 }; pGen->cityXY[14][1] = { 604.00000000000000 };
	pGen->cityXY[15][0] = { 1138.0000000000000 }; pGen->cityXY[15][1] = { 291.00000000000000 };
	pGen->cityXY[16][0] = { 336.00000000000000 }; pGen->cityXY[16][1] = { 684.00000000000000 };
	pGen->cityXY[17][0] = { 1224.0000000000000 }; pGen->cityXY[17][1] = { 311.00000000000000 };
	pGen->cityXY[18][0] = { 1235.0000000000000 }; pGen->cityXY[18][1] = { 399.00000000000000 };
	pGen->cityXY[29][0] = { 151.00000000000000 }; pGen->cityXY[29][1] = { 581.00000000000000 };
	pGen->cityXY[20][0] = { 772.00000000000000 }; pGen->cityXY[20][1] = { 219.00000000000000 };
	pGen->cityXY[21][0] = { 484.00000000000000 }; pGen->cityXY[21][1] = { 321.00000000000000 };
	pGen->cityXY[22][0] = { 618.00000000000000 }; pGen->cityXY[22][1] = { 361.00000000000000 };
	pGen->cityXY[23][0] = { 389.00000000000000 }; pGen->cityXY[23][1] = { 518.00000000000000 };
	pGen->cityXY[24][0] = { 289.00000000000000 }; pGen->cityXY[24][1] = { 145.00000000000000 };
	pGen->cityXY[25][0] = { 819.00000000000000 }; pGen->cityXY[25][1] = { 520.00000000000000 };
	pGen->cityXY[26][0] = { 287.00000000000000 }; pGen->cityXY[26][1] = { 120.00000000000000 };
	pGen->cityXY[27][0] = { 339.00000000000000 }; pGen->cityXY[27][1] = { 140.00000000000000 };
	pGen->cityXY[28][0] = { 75.000000000000000 }; pGen->cityXY[28][1] = { 602.00000000000000 };
	pGen->cityXY[29][0] = { 71.000000000000000 }; pGen->cityXY[29][1] = { 34.000000000000000 };
}