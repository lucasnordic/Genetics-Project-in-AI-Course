#include "Genetics.h"

/*
* Constructors
*/
Genetics::Genetics() {
	mChromDistances = new double[M_CH];
	mGenePool = new int[M_CH][M_CI];
	mIdxArr = new int[M_CH];
	mCityXY = new double[M_CI][2];
	similarChroms = new bool[M_CH] { false };

	mMaxWidth = 1280;
	mMaxHeight = 720;
	mDistToSide = 50.0;
	mBestStartChromLength = 0;
	//mRandSeed = rand(); //Random seed for each gene pool
	//mRandSeed = 0; //Random seed for each gene pool
	mUseMPC = true;
	mFromId = 0;
	mToId = 0;
	mWidth = 0;
	mGenId = 0;
}

Genetics::Genetics(int w, int h, int d) {
	mChromDistances = new double[M_CH];
	mGenePool = new int[M_CH][M_CI];
	mIdxArr = new int[M_CH];
	mCityXY = new double[M_CI][2];
	similarChroms = new bool[M_CH]{ false };

	mMaxWidth = w;
	mMaxHeight = h;
	mDistToSide = d;
	mBestStartChromLength = 0;
	//mRandSeed = rand(); //Random seed for each gene pool
	//mRandSeed = 0; //Random seed for each gene pool
	mUseMPC = true;
	mFromId = 0;
	mToId = 0;
	mWidth = 0;
	mGenId = 0;
}

/*
* Public methods
*/
void Genetics::InitCitiesXY() {
	For(i, M_CI) For(j, 2) {
		double rand = 0.0;
		if (j == 0) {			// if we want to change x value
			rand = randomInt(mDistToSide, mMaxWidth - mDistToSide);
		}
		else {					// else, y value
			rand = randomInt(mDistToSide, mMaxHeight - mDistToSide);
		}
		mCityXY[i][j] = rand;	// set x or y to random value
	}
}

void Genetics::InitFirstGen() {
	int mutAmo = 1337;				// amount of mutations
	For(i, M_CH) {					// For all chromosomes
		//Chrom c = { i };
		For(j, M_CI) {
			//c.genes[j] = j;
			mGenePool[i][j] = j;	// Initialize the gene pool with {0,1...N-1}
			//mGenes.insert(multimap_d_2dArray::value_type(mChromDistances[i], mGenePool[i]));
			//mGenes.insert(multimap_d_2dArray::value_type(0.0, mGenePool[i]));
		}
		//For(i, M_CI) 
		For(j, mutAmo) Mutate(i);	// mutate chromosomes
	}
	For(i, M_CH) mIdxArr[i] = i;	// Initialize with {0,1...N-1}
	Genetics::EvalDistances();		// Evaluate all distances
}

void Genetics::NewGeneration() {
	For(i, GFX_GENS_PER_FRAME) DoNextGen();
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

		//delete pArr;
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

		//delete pArr;
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

void Genetics::SwapDouble(double& from, double& to) {
	double temp = from;
	from = to;
	to = temp;
}

void Genetics::Mutate(int chromID) {
	// Choose two random indexes
	int a = randomInt(0, (M_CI-1));
	int b = randomInt(0, (M_CI-1));

	// Don't allow same index
	while (a == b) { b = randomInt(0, (M_CI-1)); };

	// Switch values
	SwapInt(mGenePool[chromID][a], mGenePool[chromID][b]);
	//int temp = mGenePool[chromID][a];
	//mGenePool[chromID][a] = mGenePool[chromID][b];
	//mGenePool[chromID][b] = temp;

	if (mGenId != 0) CheckChrom(chromID);
}

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

// TODO: Rework this so it checks if [31204] == [40213] && [31204] == [43120] (?)
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
	std::fill_n(similarChroms, M_CH, false);							  // reset the similar array
	
	// 1.1 Sort indexes based on fitness; low -> high
	std::sort(mIdxArr, mIdxArr + M_CH, sort_indices_d(mChromDistances));
	if (mGenId == 0) mBestStartChromLength = mChromDistances[mIdxArr[0]]; // set best fitness at start
	
	// 1.2 Copy gene pool and distances arrays to keep original ordering
	double prevD[M_CH]{};
	int prevGP[M_CH][M_CI]{};
	for (int i = 0; i < M_CH; i++) {
		prevD[i] = mChromDistances[i];
		for (int j = 0; j < M_CI; j++)
		{
			prevGP[i][j] = mGenePool[i][j];
		}
	}
	//std::unordered_multiset<int> ufArr1 = {2,3,1,3,4};				// found fitness values
	//auto ser = ufArr1.find(1);
	// 1.3 Sort arrays by mIdxArr ordering
	std::unordered_multiset<double> searched = { }; // found fitness values
	//std::unordered_multiset<double> duplicates = { }; // found fitness values
	int ufArr3[M_CH][M_CI] = {}; // found fitness values
	//double testest = 2635.33179999999993;
	//ufArr.insert(testest);
	//std::unordered_multiset<double> ufArr = { };				
	//auto ser2 = ufArr.find(2635.33179999999993);
	//if (ser2 != ufArr.end()) {
	//	auto test = *ser2;
	//	std::cout.precision(18);
	//	std::cout << "Found     " << (*ser2) << '\n';
	//}
	//else {
	//	std::cout.precision(18);
	//	std::cout << "Not found " << "2" << '\n';
	//}

	int idx = 0;
	int idx2 = (M_CH - 1);
	double prevf = 0.0;
	double prev = 0.0;
	For(i, M_CH) {
		int currId = mIdxArr[i];
		double *cDist = &mChromDistances[i];
		bool found = false;
		double testy = mChromDistances[i];
		double prevDistance = prevD[currId];

		auto search = searched.find(prevD[currId]);				// conditional iterator
		if (search != searched.end()) {
			//found = true;
			//auto test = *search;
			//std::cout.precision(17);
			//if (*search == prevf) std::cout << "Found     " << "." << '\n';
			//if (*search != prevf) std::cout << "Found     " << (*search) << '\n';
			//prevf = *search;
			
			//duplicates.insert(*search);
			For(j, M_CI) SwapInt(mGenePool[idx2][j], prevGP[currId][j]);// 1.3.2 swap genes
			mChromDistances[idx2] = *search;
			if (idx2 >= 0) idx2--;
		}
		else {
			SwapDouble(mChromDistances[idx], prevD[currId]);					  // 1.3.1 swap distances
			For(j, M_CI) SwapInt(mGenePool[idx][j], prevGP[currId][j]);// 1.3.2 swap genes
			if (idx <= M_CH) idx++;

			//std::cout.precision(17);
			//if (prevD[currId] == prev) std::cout << "Not found " << "." << '\n';
			//if (prevD[currId] != prev) std::cout << "Not found " << prevD[currId] << '\n';
			//prev = prevD[currId];
		}
		
		//if (search != ufArr.end()) found = true;		// if fitness value found	
		searched.insert(prevDistance);					// insert to list

		//if (!found) {
		//	SwapDouble(mChromDistances[idx], prevD[currId]);					  // 1.3.1 swap distances
		//	For(j, M_CI) SwapInt(mGenePool[i][j], prevGP[currId][j]);// 1.3.2 swap genes
		//	idx++;
		//}
	}
	//int testtes = M_CH - idx2;
	//int i = 0;
	//for (auto iter = duplicates.begin(); iter != duplicates.end(); ++iter) {
	//	mChromDistances[idx + i] = *iter;
	//	i++;
	//}
	//for (int i = 0; i < M_CH; i++){ mChromDistances[idx+i] = ufArr2 }
	//For_c(idx, M_CH) { mChromDistances[idx] = ufArr2[]; }
	
	// 2.1 Mark the chroms that are similar
	//For(i, M_CH - 1) {				// For all chroms until next to last
	//	for (int j = i + 1; j < M_CH - 1; j++) {
	//		double dis0 = mChromDistances[i];
	//		double dis1 = mChromDistances[j];
	//		bool sameChrom = false;
	//		bool sameFit = false;
	//		
	//		// 2.2 Compare fitness value
	//		sameFit = (dis0 == dis1);
	//		if (!sameFit) break;	// we don't need to look at chroms if not same fitness value

	//		// 2.3 Compare their genes only if same fitness value
	//		sameChrom = AreChromsEqual(i, j);
	//		if (sameFit && sameChrom) {
	//			// 2.4 if they are 100% equal, save to similarChroms
	//			similarChroms[i] = true;
	//		}
	//	}
	//}
	For(i, M_CH) CheckChrom(i);
}

void Genetics::CopyParents() {
	int cChromId = Genetics::PARENTS;	// child
	int pChromId = 0;		// parent
	For_c(cChromId, M_CH) {	// For all children
		For(i, M_CI) {		// For all chromosomes
			//if (similarChroms[i]) return;
			
			//int* c = &mGenePool[cChromId][mIdxArr[i]];
			//int p = mGenePool[pChromId][mIdxArr[i]];
			//*c = p;			// copy (p)arent chroms to (c)hild
			mGenePool[cChromId][i] = mGenePool[pChromId][i];
		}		

		pChromId++;
	}
}

void Genetics::TransformChildren() {
	for (int c = PARENTS; c < M_CH; c++) {				// for all children
		int rand = randomInt(0, GFX_GENS_PER_FRAME + 1);// +1 to get correct num gen

		if (rand == 0) { Mutate(c); }					// 1 in GFX_GE... chance of mutation
		else if (c != M_CH - 1) { MPC(c, c + 1); }		// send in two children
		else { Mutate(c); }								// simply mutate last child
	}
}

void Genetics::DoNextGen() {
	if (mGenId >= MAX_GEN) return;	// stop when reaching max gen limit
	if (TIME_TO_SLEEP > 0) std::this_thread::sleep_for(std::chrono::milliseconds(TIME_TO_SLEEP));

	// Debugging
	std::cout.precision(18);
	std::cout << "Gen:" << mGenId << " Best:" << mChromDistances[0];
	//For(j, M_CI) std::cout << mGenePool[0][j];
	std::cout << std::endl;
	
	// Main generation generation methods
	ReshapeGenePool();				// Reshape the gene pool
	CopyParents();					// clone chroms from parents to children
	TransformChildren();			// Transform clones/children
	EvalDistances();				// Re-evaluate distances/fitness values

	mGenId++;						// increase generation
}

void Genetics::MPC(int cId0, int cId1) {
	int original[M_CI] = {};
	int original1[M_CI] = {};
	For(i, M_CI) { original[i] = mGenePool[cId0][i]; original1[i] = mGenePool[cId1][i]; }
	//int test[M_CI] = { 3,4,2,5,0,1 };
	//int test1[M_CI] = { 4,5,3,1,2,0 };
	//For(i, M_CI) { mGenePool[cId0][i] = test[i]; mGenePool[cId1][i] = test1[i]; }

	if (!mUseMPC) return;
	if (M_CH <= 4) return;	// prevent wrong mpc cutting if chrom 2 small

	// 1. Set range for gene swapping
	mWidth = M_CI / 2;		// max allowable width = half a chrom
	mFromId = randomInt(0, mWidth);
	mToId = mFromId + mWidth - 1;	

	mFromId = 1;
	mToId = 3;

	// 2. For the genes in range, swap them
	for(int i = mFromId; i <= mToId; i++){ 
		/*SwapInt(mGenePool[cId0][i], mGenePool[cId1][i]);*/
		int tmp = mGenePool[cId1][i];
		mGenePool[cId1][i] = mGenePool[cId0][i];
		mGenePool[cId0][i] = tmp;
	}

	// 3. Repair the damaged genes || Repair the damaged chromosome. (duplicates)
	MPC_RepairChrom(cId0); 
	MPC_RepairChrom(cId1);
}

void Genetics::MPC_RepairChrom(int chromId) {
	// Step 1: find duplicates
	int duplicatesList[M_CI]{};
	For(i, M_CI) duplicatesList[i] = -2; // initiate with -2 since we don't consider negative
	int numberOfDuplicatesN = 0;

	// TODO; improve O(n)
	for (int i = 0; i < M_CI; i++)					// for all genes
	{
		int cCh1 = mGenePool[chromId][i];			// gene to compare
		int cCh2;									// second gene to compare

		if (i >= mFromId && i <= mToId) {
			for (int j = i + 1; j < M_CI; j++)	// for all switched letters on left
			{
				cCh2 = mGenePool[chromId][j];
				if (cCh1 == cCh2 && cCh1 != -1)
				{
					mGenePool[chromId][j] = -1;
					duplicatesList[numberOfDuplicatesN] = j; 
					numberOfDuplicatesN++;
				}
			}
		}
		else {
			for (int j = i + 1; j < M_CI; j++)	// for all switched letters on left
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
	bool numbersIncluded[M_CI]{false};
	for (size_t i = 0; i < M_CI; i++)
	{
		int c = mGenePool[chromId][i]; // current
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

	CheckChrom(chromId);
}

bool Genetics::CheckChrom(int id) {
	std::multiset <int> temp = {};
	For(i, M_CI) { temp.insert(mGenePool[id][i]); }
	For(i, M_CI) { 
		int count = temp.count(i);
		if (count > 1) {
			std::cout << "ALERT!" << std::endl;
			return false;
		}
	}

	return true;
}