#include "Genetics.h"

/*
* Constructors
*/
Genetics::Genetics() {
	mChromDistances = new double[M_CH];
	mGenePool = new int[M_CH][M_CI];
	mIdxArr = new int[M_CH];
	mCityXY = new double[M_CI][2];

	mMaxWidth = 1280;
	mMaxHeight = 720;
	mDistToSide = 50.0;
	mBestStartChromLength = 0;
	mPrevGenBestChromLength = 0;
	mUseMPC = true;
	mGenId = 0;
}

Genetics::Genetics(int w, int h, int d) {
	mChromDistances = new double[M_CH];
	mGenePool = new int[M_CH][M_CI];
	mIdxArr = new int[M_CH];
	mCityXY = new double[M_CI][2];

	mMaxWidth = w;
	mMaxHeight = h;
	mDistToSide = d;
	mBestStartChromLength = 0;
	mPrevGenBestChromLength = 0;
	mUseMPC = true;
	mGenId = 0;
}

/*
* Public methods
*/
bool Genetics::InitCitiesXY() {
	For(i, M_CI) For(j, 2) {
		double rand = 0.0;
		if (j == 0) {				// if x value
			rand = randomInt(mDistToSide, mMaxWidth - mDistToSide);
		}
		else {						// else, y value
			rand = randomInt(mDistToSide, mMaxHeight - mDistToSide);
		}
		mCityXY[i][j] = rand;		// set x or y to random value
	}

	return true;
}

bool Genetics::InitFirstGen() {
	mTimeStart = std::chrono::high_resolution_clock::now();

	For(i, M_CH) {					// 1.1 For all chromosomes
		For(j, M_CI) {				// 1.2 For their genes
			mGenePool[i][j] = j;	// 2.1 Initialize the gene pool with {0,1...N-1}
		}

		For(j, INIT_MUTATIONS) Mutate(i); // 2.2 mutate chromosomes

		mIdxArr[i] = i;				// Initialize with {0,1...N-1}
	}

	Genetics::EvalDistances();		// Evaluate all distances

	return true;
}

void Genetics::NewGeneration() {
	For(i, GFX_GENS_PER_FRAME) DoNextGen();
}

void Genetics::DrawCities() {
	// Draw points for the generated cities
	glColor3ub(55, 215, 55);
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

void Genetics::DrawPaths() {
	// Draw lines between cities from the chromosome with best fitness value
	glColor3ub(15, 75, 175);
	glBegin(GL_LINE_LOOP);			// Function for displaying glVertex2dv();
	For(i, M_CI) {					// For all cities
		double xyArr[2]{};
		double* pArr{ mCityXY[mGenePool[0][i]] };
		xyArr[0] = *pArr; pArr++;
		xyArr[1] = *pArr; pArr--;
		glVertex2dv(xyArr);			// specify x,y for draw
	}
	glEnd();
}

/*
* Private methods
*/
void Genetics::SwapInt(int& from, int& to) {
	int temp = from;
	from = to;
	to = temp;
}

bool Genetics::CheckChrom(int id) {
	std::multiset <int> temp = {};
	For(i, M_CI) { temp.insert(mGenePool[id][i]); }
	For(i, M_CI) {
		int count = temp.count(i);
		if (count > 1) { // a similar gene is found, not good.
			std::cout << "Error, the chromosome has duplicate gene: " << i << std::endl;
			return false;
		}
	}

	return true;
}

void Genetics::Mutate(int chromID) {
	// Choose two random indexes
	int a = randomInt(0, (M_CI-1));
	int b = randomInt(0, (M_CI-1));

	// Don't allow same index
	while (a == b) { b = randomInt(0, (M_CI-1)); };

	// Switch values
	SwapInt(mGenePool[chromID][a], mGenePool[chromID][b]);
}

void Genetics::EvalDistances() {
	double l{0.0};									// lowest fitness value
	mPrevGenBestChromLength = mChromDistances[0];

	For(i, M_CH) {									// For all Chromosomes
		mChromDistances[i] = CalcTotalDistance(i);	// Calculate fitness value

		if (mChromDistances[i] < l || l == 0.0) l = mChromDistances[i];
	}

	if (mGenId == 0) { mBestStartChromLength = l; }		// store best initial chrom fitness
}

double Genetics::CalcTotalDistance(int chromID) {
	double result{ 0.0 };
	double disX, disY;		// distance between two cities x,y coordinates

	For(i, M_CI) {			// 1. For all genes in chromosome, calc total dist between points
		int city0 = mGenePool[chromID][i];
		int city1;
		if ((i + 1) == M_CI) { city1 = mGenePool[chromID][0];}	// if at last index, point first city
		else { city1 = mGenePool[chromID][i + 1]; }				// else point next city
		disX = mCityXY[city0][0] - mCityXY[city1][0];
		disY = mCityXY[city0][1] - mCityXY[city1][1];
		result = result + sqrt(disX * disX + disY * disY);		// Add euclidean distance of two points
	}

	return result;			// 2. Return total distance (fitness value)
}

bool Genetics::AreChromsEqual(int chromID0, int chromID1) {
	bool foundSimilarGene = false;
	For(i, M_CI) {
		int c0 = mGenePool[chromID0][i]; int c1 = mGenePool[chromID1][i];
		if (c0 != c1) return false;
	}

	return true;
}

// TODO: Improve O(n); Ideas: use a map with struct{id,genes[],fitness}, sorted by fitness
void Genetics::ReshapeGenePool() {
	// 1. Sort indexes of chroms based on fitness; low -> high
	std::sort(mIdxArr, mIdxArr + M_CH, sort_indices_d(mChromDistances));  
	if (mGenId == 0) mBestStartChromLength = mChromDistances[mIdxArr[0]]; // set best fitness at start
	
	// 2. Copy gene pool and distances arrays to keep original ordering
	double tempD[M_CH]{};
	int prevGP[M_CH][M_CI]{};
	for (int i = 0; i < M_CH; i++) {
		tempD[i] = mChromDistances[i];
		for (int j = 0; j < M_CI; j++)
		{
			prevGP[i][j] = mGenePool[i][j];
		}
	}

	// 3. Sort gene pool and distances by highest fitness ordering (stored in mIdxArr)
	//	   , in order to place parents(highest fitness) on top and children on bottom
	std::unordered_multiset<double> searched = {}; // found fitness values
	int idx = 0; int idx2 = (M_CH - 1); // idx = (0)+; idx2 = (N-1)-;
	For(i, M_CH) {													// 3.0 for all genes, distances
		int id = mIdxArr[i]; // id of fit value
		double d = tempD[id]; // original dist

		auto search = searched.find(d);								// conditional iterator
		if (search != searched.end()) {								// if fitness value already found
			For(j, M_CI) SwapInt(mGenePool[idx2][j], prevGP[id][j]);// 3.1 swap genes
			mChromDistances[idx2] = *search;						// 3.2 swap distances
			if (idx2 >= 0) idx2--;
		}
		else {
			For(j, M_CI) SwapInt(mGenePool[idx][j], prevGP[id][j]);	// 3.1 swap genes
			mChromDistances[idx] = d;								// 3.2 swap distances
			if (idx < M_CH) idx++;
		}

		searched.insert(d);	// insert to list of searched dist
	}
}

void Genetics::CopyParents() {
	int c = Genetics::PARENTS;					// child
	int p = 0;									// parent
	For_c(c, M_CH) {							// For all children chromosomes
		For(i, M_CI) {							// For all child genes
			mGenePool[c][i] = mGenePool[p][i];	// Copy parent genes to child
		}		

		p++;
	}
}

void Genetics::TransformChildren() {
	for (int c = PARENTS; c < M_CH; c++) {				// 1. for all children
		int rand = randomInt(0, GFX_GENS_PER_FRAME + 1);

		if (rand == 0) { Mutate(c); }					// 2.1 if (1 in GFX_GE...) Mutate()
		else if (c != M_CH - 1) { MPC(c, c + 1); }		// 2.2 else do MPC() on two children
		else { Mutate(c); }								// 2.3 if last index, simply Mutate()
	}
}

void Genetics::DoNextGen() {
	if (mGenId >= MAX_GEN) { mRunning = false; return; }	// when reaching max gen limit, stop
	if (TIME_TO_SLEEP > 0) std::this_thread::sleep_for(std::chrono::milliseconds(TIME_TO_SLEEP));

	// save current time
	auto now = std::chrono::high_resolution_clock::now();
	mCurrTime = std::chrono::duration_cast<std::chrono::milliseconds>(now - mTimeStart); // store time every gen
	
	// Main Genetic Functions
	ReshapeGenePool();				// Reshape the gene pool
	if (mChromDistances[0] < mPrevGenBestChromLength) mTimeEnd = mCurrTime; // store time until best chrom
	CopyParents();	
	TransformChildren();			// Transform clones/children
	EvalDistances();				// Re-evaluate distances/fitness values
	mGenId++;						// increase generation
}

void Genetics::MPC(int cId0, int cId1) {
	if (!mUseMPC) return;
	if (M_CH <= MPC_MIN_GENES) return;// prevent wrong mpc cutting if chrom too small

	// 1. Set range for gene swapping
	int width = M_CI / 2;				  // max allowable width = half a chrom
	int fromId = randomInt(0, width);	  // cut start index
	int toId = fromId + width - 1;	  // cut end index

	// 2. For the genes in range, swap them
	for(int i = fromId; i <= toId; i++){
		/*SwapInt(mGenePool[cId0][i], mGenePool[cId1][i]);*/
		int tmp = mGenePool[cId1][i];
		mGenePool[cId1][i] = mGenePool[cId0][i];
		mGenePool[cId0][i] = tmp;
	}

	// 3. Repair the damaged genes(duplicate genes) || Repair the damaged chromosome
	MPC_RepairChrom(cId0, fromId, toId); MPC_RepairChrom(cId1, fromId, toId);
}

void Genetics::MPC_RepairChrom(int chromId, int fromId, int toId) {
	// Step 1: find duplicates
	int duplicatesList[M_CI]{};
	For(i, M_CI) duplicatesList[i] = -2;		// initiate with -2 since we don't consider negative
	int numberOfDuplicatesN = 0;
	for (int i = 0; i < M_CI; i++)				// for all genes
	{
		int cCh1 = mGenePool[chromId][i];		// gene to compare
		int cCh2;								// second gene to compare

		if (i >= fromId && i <= toId) {			// [0,from,to,0] if gene is affected by MPC 
			for (int j = i + 1; j < M_CI; j++)	// for all genes
			{
				cCh2 = mGenePool[chromId][j];
				if (cCh1 == cCh2 && cCh1 != -1)	// if a similar gene is found
				{
					mGenePool[chromId][j] = -1;	// mark as duplicate 
					duplicatesList[numberOfDuplicatesN] = j; 
					numberOfDuplicatesN++;
				}
			}
		}
		else {									// [x,0,0,x] if gene is on left or right of MPC
			for (int j = i + 1; j < M_CI; j++)	// for all genes
			{
				cCh2 = mGenePool[chromId][j];
				if (cCh1 == cCh2 && cCh1 != -1)
				{
					mGenePool[chromId][i] = -1;
					duplicatesList[numberOfDuplicatesN] = i; 
					numberOfDuplicatesN++;
				}
			}
		}
	}

	// Step 2: Sort the duplicates || Selection sort | Time Complexity O(n^2)!
	for (size_t i = 0; i < numberOfDuplicatesN; i++)
	{
		for (size_t j = i + 1; j < numberOfDuplicatesN; j++)
		{
			if (duplicatesList[j] < duplicatesList[i]) {
				SwapInt(duplicatesList[i], duplicatesList[j]);
			}
		}
	}

	// Step 3: Find out which numbers that are present in the chrom except for -1
	bool numbersIncluded[M_CI]{false};
	for (size_t i = 0; i < M_CI; i++)
	{
		int c = mGenePool[chromId][i];
		if (c >= 0) {
			numbersIncluded[c] = true;
		};
	}

	// Step 4: Find out which numbers that are missing
	int missingNumbersList[M_CI]{};
	For(i, M_CI) missingNumbersList[i] = -1;
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
		mGenePool[chromId][duplicatesList[i]] = curr;
	}
}