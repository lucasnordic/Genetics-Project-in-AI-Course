#ifndef GENETICS_H_2022_10_20
#define GENETICS_H_2022_10_20

#include <Math.h>
#include <stdio.h>
#include <time.h>
//#include <stdarg.h>
//#include <thread>
#include <iostream>
#include <random> // for random num generation

#include <GL/glut.h>

#define For(i,N) for (int (i) = 0; (i) < (N); (i)++)

class Genetics {
	static const int M_CI = 10;				// default 100 | MAX_CITIES
	static const int M_CH = 20;		        // > 3; default 20 | MAX_CHROMOSOMES
	static const int MAX_GEN = 60000;		// default 300
	//static const int MPCS_PER_MUTATION = 3;	// default 3
	//static const int PARENTS_N = .6 * M_CH;	// default .7
	//___________________________________________________
	static const int GFX_GENS_PER_FRAME = 1;// default 50
	static const int TIME_TO_SLEEP = 0;		// default 0 - milliseconds 
	//___________________________________________________
public:
	Genetics();
	Genetics(int w, int h);
	void InitCitiesXY();					// 1. Initialize city x,y
	void InitFirstGen();					// 2. Initialize first gen
	void NewGeneration();                   // 3.1 generate new generation
	void DrawCitiesPaths();					// 3.2 draw cities and paths(for best Chromosome)
	//___________________________________________________
	double	(*mCityXY)[2];					//[MAX_CITIES][X/Y]
	unsigned int mRandSeed;                 // randomized seed
	bool	useMPC;							// whether to use Multi Point Crossover
	double	mBestStartChromLength;			//best chrom with shortest length, at beginning; to evaluate algo success
	int		mFromIdx, mToIdx, mWidth;
	//___________________________________________________
	int		(*mGenePool)[M_CI];				//[MAX_CHROMOSOMES][MAX_CITIES] // Gene Pool
	double* mChromDistances;
	int		mBestCurrentChromIdx;			//best current Chromosomes index
	int		mGenerationId;						
	//___________________________________________________
    int  maxWidth, maxHeight;				// max width/height to generate cities
private:
	template<typename T> T randomInt(T range_from, T range_to) { // For generating rand int numb
		std::random_device                  rand_dev;
		std::mt19937                        generator(rand_dev());
		std::uniform_int_distribution<T>    distr(range_from, range_to);
		return distr(generator);
	}
	template<typename T> T random(T range_from, T range_to) {	// For generating rand numb
		std::random_device                  rand_dev;
		std::mt19937                        generator(rand_dev());
		std::uniform_real_distribution<T>    distr(range_from, range_to);
		return distr(generator);
	}
	void	Mutate(int chromID);						// Mutate a Chromosome
	double	CalcTotalDistance(int chromID);			// GOAL FUNCTION || total distance for a chromosome
	void	EvalDistances();							// Evaluate goal function/distances
	void	ReshapeGenePool();
	void	EvalNextGen();								// Evaluate Next generation
	//double RandomBetween(double min, double max);
    //-------------------------------------------------------
	void	MPC_RepairChrom(int chromIdx);				// Repair the chrome after MPC
};

#endif