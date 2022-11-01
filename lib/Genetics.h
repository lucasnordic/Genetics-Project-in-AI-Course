#ifndef GENETICS_H_2022_10_20
#define GENETICS_H_2022_10_20

#include <Math.h>
#include <time.h>
#include <thread>		// for TIME_TO_SLEEP
#include <chrono>		// for TIME_TO_SLEEP
#include <iostream>		// debugging
#include <random>		// for uniform_x_distribution
#include <algorithm>	// for sorting
#include <unordered_set>
#include <set>			

#include <GL/glut.h>

#define For(i,N) for (int (i) = 0; (i) < (N); (i)++)	
#define For_c(i,N) for ((i); (i) < (N); (i)++) // user assigned i value

class Genetics {
	static const int M_CI = 20;					// default 20  ;> 3; MAX_CITIES
	static const int M_CH = 100;				// default 100 ;MAX_CHROMOSOMES
	static const int MAX_GEN = 1337;			// default 300
	static const int MPCS_PER_MUT = 3;			// default 3   ;how many MPC's happen for each mutation
	static const int MPC_MIN_GENES = 4;			// default 4   ;minimum allowed genes for MPC
	static const int GFX_GENS_PER_FRAME = 2;	// default 1
	static const int PARENTS = .5 * M_CH;		// default .5 * M_CH ;Percentage of parents
	static const int TIME_TO_SLEEP = 0;			// default 0   ;milliseconds ;slow down generation
	static const int INIT_MUTATIONS = M_CI*20;  // default M_CI*20 ;N mutations on initial pop
public:
	Genetics();
	//~Genetics(); // TODO: Make memory leaks a thing of the past
	Genetics(int w, int h, int d);
	//___________________________________________________
	bool	InitCitiesXY();						// 1. Initialize city x,y
	bool	InitFirstGen();						// 2. Initialize first gen
	void	NewGeneration();					// 3  generate new generation
	void	DrawCities();						// 4.1 draw cities
	void	DrawPaths();						// 4.2 draw paths( for best Chromosome )
	//___________________________________________________
	int		maxGen = MAX_GEN;
	int		mMCI = M_CI;
	int		mMCH = M_CH;
	//___________________________________________________
	double	(*mCityXY)[2];						// [MAX_CITIES][X/Y]
	bool	mUseMPC;							// whether to use Multi Point Crossover
	double	mBestStartChromLength;				// best chrom, at beginning; to eval algo success
	double	mPrevGenBestChromLength;
	int		mGenId;						
    int		mMaxWidth, mMaxHeight, mDistToSide;	// Max width/height to generate cities;
	int		(*mGenePool)[M_CI];					// [MAX_CHROMOSOMES][MAX_CITIES] // Gene Pool
	double* mChromDistances;					// store fitness values
	bool	mRunning = true;					// if DoNextGen is happening();
	std::chrono::time_point<std::chrono::high_resolution_clock> mTimeStart{};	// gen start
	std::chrono::milliseconds mCurrTime{};										// time alive
	std::chrono::milliseconds mTimeEnd{};										// gen end
private:
	//___________________________________________________
	//Modified from https://stackoverflow.com/questions/3909272/sorting-two-corresponding-arrays
	class sort_indices_d {						// sort indices (0,1,2...N-1) based on double arr
	private:
		double* mpArr;
	public:
		sort_indices_d(double* parr) : mpArr(parr) {}
		bool operator()(int i, int j) const { 
			return mpArr[i] < mpArr[j]; 
		}
	};
	//___________________________________________________
	//For generating rand int numb
	template<typename T> T randomInt(T range_from, T range_to) {
		std::random_device                  rand_dev;
		std::mt19937                        generator(rand_dev());
		std::uniform_int_distribution<T>    distr(range_from, range_to);
		return								distr(generator);
	}
	//___________________________________________________
	int*	mIdxArr;					// stores INDEXES of mChromDistances, based on fitness score low -> high
	//___________________________________________________
	//		Helper methods
	void	SwapInt(int& from, int& to);				// swap 2 int
	bool	CheckChrom(int id);							// Helper function to check chrom health
	//___________________________________________________
	//		Gene methods
	void	Mutate(int chromID);						// Mutate a Chromosome
	void	EvalDistances();							// Evaluate goal function/distances
	double	CalcTotalDistance(int chromID);				// FITNESS FUNCTION || GOAL FUNCTION || total distance
	void	DoNextGen();								// Evaluate Next generation
	void	ReshapeGenePool();							// Sort by fitness; Check for similar chroms
	bool	AreChromsEqual(int chromID0, int chromID1);	// Check if chroms are equal;
	void	CopyParents();								// Create children from parents
	void	TransformChildren();						// Mutate and do MPC on children
	void	MPC(int chromID0, int chromID1);			// Multi Point Crossover
	void	MPC_RepairChrom(int chromId, int fromId, int toId);	// Repair the chrome after MPC
	//-------------------------------------------------------
};

#endif