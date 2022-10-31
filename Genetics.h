#ifndef GENETICS_H_2022_10_20
#define GENETICS_H_2022_10_20

//#include <stdio.h>
//#include <stdarg.h>
#include <Math.h>
#include <time.h>
#include <thread>		// for SLEEP
#include <chrono>		// for SLEEP
#include <iostream>
#include <random>		// for uniform_x_distribution
#include <algorithm>	// for sorting
#include <unordered_set>
#include <set>

#include <GL/glut.h>

#define For(i,N) for (int (i) = 0; (i) < (N); (i)++)	// 
#define For_c(i,N) for ((i); (i) < (N); (i)++)			// user assigned i value

class Genetics {
	static const int M_CI = 20;				// default 20 ; > 3; MAX_CITIES
	static const int M_CH = 100;		    // default 100; MAX_CHROMOSOMES
	static const int MAX_GEN = 300;		// default 300;
	static const int MPCS_PER_MUT = 3;		// default 3  ;how many MPC's happen for each mutation
	static const int GFX_GENS_PER_FRAME = 1;// default 50 ;
	static const int TIME_TO_SLEEP = 0;		// default 0  ;milliseconds 
	static const int PARENTS = .7 * M_CH;	// default .7
public:
	Genetics();
	Genetics(int w, int h, int d);
	void	InitCitiesXY();						// 1. Initialize city x,y
	void	InitFirstGen();						// 2. Initialize first gen
	void	NewGeneration();					// 3.1 generate new generation
	void	DrawCitiesPaths();					// 3.2 draw cities and paths(for best Chromosome)
	//___________________________________________________
	//___________________________________________________
	double	(*mCityXY)[2];						// [MAX_CITIES][X/Y]
	//unsigned int mRandSeed;					// randomized seed
	bool	mUseMPC;							// whether to use Multi Point Crossover
	double	mBestStartChromLength;				// best chrom, at beginning; to eval algo success
	int		mGenId;						
    int		mMaxWidth, mMaxHeight, mDistToSide;	// Max width/height to generate cities;
	//int		mMaxChromGeneComp;					// TODO: ?
	//		TODO: Find a better way to sort arrays below
	int		(*mGenePool)[M_CI];			// [MAX_CHROMOSOMES][MAX_CITIES] // Gene Pool
	double* mChromDistances;			// store fitness values
	/*std::vector<std::vector<int>> mGenePool;
	std::vector<double> mChromDistances;*/
	bool*	similarChroms;				// for storing idx of similar chromosomes
	//struct Chrom { int id; int genes[M_CI]{};};
	//std::multimap<double[M_CH], int[M_CH][M_CI]> mGenes;
	//typedef std::multimap<double, int[M_CI]> multimap_d_2dArray;
	//multimap_d_2dArray mGenes;
	int*	mIdxArr;				// stores sorted indexes of mChromDistances, based on fitness {3,1...N-1}
private:
	//___________________________________________________
	// Modified from https://stackoverflow.com/questions/3909272/sorting-two-corresponding-arrays
	class sort_indices_d {			// sort indices (0,1,2...N-1) based on double arr
	private:
		double* mpArr;
	public:
		sort_indices_d(double* parr) : mpArr(parr) {}
		bool operator()(int i, int j) const { return mpArr[i] < mpArr[j]; }
	};
	//___________________________________________________
	template<typename T> T randomInt(T range_from, T range_to) { // For generating rand int numb
		std::random_device                  rand_dev;
		std::mt19937                        generator(rand_dev());
		std::uniform_int_distribution<T>    distr(range_from, range_to);
		return								distr(generator);
	}
	//___________________________________________________
	int		mFromId, mToId, mWidth;					// For MPC width, to, from
	//___________________________________________________
	//		Helper methods
	void	SwapInt(int& from, int& to);			// swap 2 int
	void	SwapDouble(double& from, double& to);	// swap 2 double
	//___________________________________________________
	//		Gene methods
	void	Mutate(int chromID);						// Mutate a Chromosome
	void	EvalDistances();							// Evaluate goal function/distances
	double	CalcTotalDistance(int chromID);				// FITNESS FUNCTION || GOAL FUNCTION || total distance
	void	DoNextGen();								// Evaluate Next generation
	void	ReshapeGenePool();							// Check for similar chroms; sort by fitness
	bool	AreChromsEqual(int chromID0, int chromID1);	// Check if chroms are equal;
	void	CopyParents();								// Create children from parents
	void	TransformChildren();						// Mutate and do MPC on children
	void	MPC(int chromID0, int chromID1);			// Multi Point Crossover
	void	MPC_RepairChrom(int chromId);				// Repair the chrome after MPC
	//-------------------------------------------------------
	bool	CheckChrom(int id);
};

#endif