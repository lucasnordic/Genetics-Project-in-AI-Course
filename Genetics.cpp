#include "Genetics.h"

/*
* Constructors
*/
Genetics::Genetics() {
	/*
	* mGenePool[mCityXY.size()]
	* a = table[0][2]
	* b = table[2][0]
	*    0  1  2  3
	* 0 [ ][ ][a][ ]
	* 1 [ ][ ][ ][ ]
	* 2 [b][ ][ ][ ]
	*/
	mChromDistances = new double[M_CH];
	mGenePool = new int[M_CH][M_CI];
	mCityXY = new double[M_CI][2];
	maxWidth = 1280;
	maxHeight = 720;
	mRandSeed = rand(); //Random seed for each gene pool
	useMPC = true;
	mFromIdx = 0;
	mToIdx = 0;
	mWidth = 0;
}

Genetics::Genetics(int w, int h) {
	mChromDistances = new double[M_CH];
	mGenePool = new int[M_CH][M_CI];
	mCityXY = new double[M_CI][2];
	maxWidth = w;
	maxHeight = h;
	mRandSeed = rand(); //Random seed for each gene pool
	useMPC = true;
	mFromIdx = 0;
	mToIdx = 0;
	mWidth = 0;
}

/*
* Public methods
*/
void Genetics::InitCitiesXY() {
	For(i, M_CI) For(j, 2) {
		double rand = 0.0;
		if (j == 0) {			// if we want to change x value
			double distToSide = 25.0;
			rand = random(distToSide, maxWidth - distToSide);
		}
		else {					// else, y value
			double distToSide = 25.0;
			rand = random(distToSide, maxHeight - distToSide);
		}
		mCityXY[i][j] = rand;	// set x or y to random value
	}
}

void Genetics::InitFirstGen() {
	int mutAmo = 1337;						// amount of mutations
	mGenerationId = 0;
	For(i, M_CH) {							// For all chromosomes
		For(j, M_CI) mGenePool[i][j] = j;	// Initialize the gene pool with 0...N-1
		For(j, mutAmo) Mutate(i);			// mutate chromosomes
	}
	Genetics::EvalDistances();				// Evaluate all distances
}

void Genetics::NewGeneration() {
	For(i, GFX_GENS_PER_FRAME) EvalNextGen();
}

void Genetics::DrawCitiesPaths() {
	// Draw lines between cities from the chromosome with least length
	glColor3ub(55, 55, 215);
	glBegin(GL_LINE_LOOP);	// Function for displaying glVertex2dv();
	For(i, M_CI) {			// For all cities
		double xyArr[2]{};
		double* pArr{ mCityXY[mGenePool[0][i]] };
		xyArr[0] = *pArr; pArr++;
		xyArr[1] = *pArr; pArr--;
		glVertex2dv(xyArr);	// specify x,y for draw
	}
	glEnd();

	// Draw points for the generated cities
	glColor3ub(215, 215, 15);
	glBegin(GL_POINTS);		// Function for displaying glVertex2dv();
	For(i, M_CI) {			// For all cities
		double xyArr[2]{};
		double* pArr{ mCityXY[i] };
		xyArr[0] = *pArr; pArr++;
		xyArr[1] = *pArr; pArr--;
		glVertex2dv(xyArr); // Specify x,y for draw
	}
	glEnd();
}

/*
* Private methods
*/
void Genetics::Mutate(int chromID) {
	// Choose two random indexes
	int a = randomInt(0, (M_CI-1));
	int b = randomInt(0, (M_CI-1));

	// Don't allow same index
	while (a == b) { b = randomInt(0, (M_CI-1)); };

	// Switch values
	int temp = mGenePool[chromID][a];
	mGenePool[chromID][a] = mGenePool[chromID][b];
	mGenePool[chromID][b] = temp;
}

//double Genetics::RandomBetween(double min, double max) {
//	return min + (rand()) / ((RAND_MAX / (max - min))); // Not ideal, but simpler
//}

void Genetics::EvalDistances() {
	double l{0.0};	// for storing lowest distance
	For(i, M_CH) {  // For all Chromosomes
		mChromDistances[i] = CalcTotalDistance(i); 

		if (mChromDistances[i] < l || l == 0.0) l = mChromDistances[i];
	}

	mBestStartChromLength = l; 
}

double Genetics::CalcTotalDistance(int chromID) {
	double result{ 0.0 };
	double disX, disY;		// distance between two cities x,y coordinates

	For(i, M_CI) {
		int city0 = mGenePool[chromID][i];
		int city1;
		if ((i + 1) == M_CI) { city1 = mGenePool[chromID][0];}	// if we are at last index, point to first city
		else { city1 = mGenePool[chromID][i + 1]; }				// else point at next city index
		disX = mCityXY[city0][0] - mCityXY[city1][0];
		disY = mCityXY[city0][1] - mCityXY[city1][1];
		result = result + sqrt(disX * disX + disY * disY);		// final calculation for Euclidean distance
	}

	return result; // return summed distances
}

void Genetics::ReshapeGenePool() {

}

void Genetics::EvalNextGen() {

}

void Genetics::MPC_RepairChrom(int chromIdx) {
	// Step 1: find duplicates
	int duplicatesList[M_CI];
	int numberOfDuplicatesN = 0;
	// TODO; improve O(n)
	for (int i = 0; i < M_CI; i++)					// for all genes
	{
		int cCh1 = mGenePool[chromIdx][i];				// gene to compare
		int cCh2;									// second gene to compare

		if (i >= mFromIdx && i <= mToIdx) {
			for (int j = i + 1; j < M_CI - 1; j++)			// for all switched letters on left
			{
				cCh2 = mGenePool[chromIdx][j];
				if (cCh1 == cCh2 && cCh1 != -1)
				{
					mGenePool[chromIdx][j] = -1; duplicatesList[numberOfDuplicatesN] = j; numberOfDuplicatesN++;
				}
			}
		}
		else {
			for (int j = i + 1; j < M_CI - 1; j++)			// for all switched letters on left
			{
				cCh2 = mGenePool[chromIdx][j];
				if (cCh1 == cCh2 && cCh1 != -1)
				{
					mGenePool[chromIdx][i] = -1; duplicatesList[numberOfDuplicatesN] = i; numberOfDuplicatesN++;
				}
			}
		}
	}

	// Step 2: Sort the duplicates - TODO; replace sorting
	for (size_t i = 0; i < numberOfDuplicatesN; i++)
	{
		for (size_t j = i + 1; j < numberOfDuplicatesN; j++)
		{
			if (duplicatesList[j] < duplicatesList[i]) {
				int temp = duplicatesList[i];
				duplicatesList[i] = duplicatesList[j];
				duplicatesList[j] = temp;
			}
		}
	}

	// Step 3: Find out which numbers that are present in the chrom except for -1
	bool numbersIncluded[M_CI];
	For(i, M_CI) numbersIncluded[i] = false;
	for (size_t i = 0; i < M_CI; i++)
	{
		int c = mGenePool[chromIdx][i]; //  current
		if (c >= 0) {
			numbersIncluded[c] = true;
		};
	}

	// Step 4: Find out which numbers that are missing
	int missingNumbersList[M_CI];
	int numberOfMissingNumbersN = 0;
	for (size_t i = 0; i < M_CI; i++)
	{
		if (numbersIncluded[i] == false) {
			missingNumbersList[numberOfMissingNumbersN] = i;
			numberOfMissingNumbersN++;
		}
	}

	// Step 5: Fill holes with the missing numbers by ascending order
	for (size_t i = 0; i < numberOfMissingNumbersN; i++)
	{
		int curr = missingNumbersList[i];
		mGenePool[chromIdx][duplicatesList[i]] = curr;
	}
}