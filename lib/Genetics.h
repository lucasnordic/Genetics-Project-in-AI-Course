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
#include <map>		
#include <iomanip>		// for setprecision

#include <GL/glut.h>

#define For(i,N) for (int (i) = 0; (i) < (N); (i)++)	
#define For_c(i,N) for ((i); (i) < (N); (i)++) // user assigned i value

class Genetics {
	static const int MAX_GEN = 13370;			// default 300
	static const int M_CI = 30;					// default 20  ;> 3; MAX_CITIES
	static const int M_CH = 200;				// default 100 ;MAX_CHROMOSOMES
	static const int PARENTS = .7 * M_CH;		// default .7 * M_CH ;Percentage of parents
	static const int MPCS_PER_MUT = 3;			// default 3   ;how many MPC's happen for each mutation
	static const int INIT_MUTATIONS = M_CI*30;  // default M_CI*20 ;N mutations on initial pop
	static const int CHILD_MUT = 1000;			// default 1   ;gene mut for each mutating child chrom, each gen
	static const int MPC_MIN_GENES = 4;			// default 4   ;minimum allowed genes for MPC
	static const int GFX_GENS_PER_FRAME = 2;	// default 1
	static const int TIME_TO_SLEEP = 0;			// default 0   ;milliseconds ;slow down generation
public:
	Genetics();
	//~Genetics(); // TODO: Make memory leaks a thing of the past
	Genetics(int w, int h, int d);
	//___________________________________________________
	bool	InitCitiesXY();								 // 1. Initialize city x,y
	bool	InitPopulation(int amountChrom, bool newPop);// 2. Initialize first gen
	void	NewGeneration();							 // 3. Generate new generation
	void	DrawCities();								 // 4.1 draw cities
	void	DrawPaths();								 // 4.2 draw paths( for best Chromosome )
	//		Gene functions
	void	Mutate(int chromID);						// Mutate a Chromosome
	//___________________________________________________
	int			maxGen = MAX_GEN;
	int	mci = M_CI;
	int	mch = M_CH;
	//___________________________________________________
	double	(*cityXY)[2];						// [MAX_CITIES][X/Y]
	bool	useMPC;								// whether to use Multi Point Crossover
	//bool	useRoulette = false;				// use roulette wheel selection?
	double	bestStartChromLength;				// best chrom, at beginning; to eval algo success
	double	prevGenBestChromLength;				// best from previous gen(s)
	int		genId;								// loop count
    int		maxWidth, maxHeight, distToSide;	// Max width/height to generate cities;
	int		(*genePool)[M_CI];					// [MAX_CHROMOSOMES][MAX_CITIES] // Gene Pool
	double* chromDistances;						// store fitness values
	bool	running = true;						// if DoNextGen is happening();
	std::chrono::milliseconds timeCurr{};		// time alive
	std::chrono::milliseconds timeEnd{};		// gen end
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
	struct {

	};
	//___________________________________________________
	std::chrono::time_point<std::chrono::high_resolution_clock> mTimeStart{};	// gen start
	int*	mIdxArr;					// stores INDEXES of chromDistances, based on fitness score low -> high
	//___________________________________________________
	//		Helper functions
	void	SwapInt(int& from, int& to);				// swap 2 int
	bool	CheckChrom(int id);							// Helper function to check chrom health
	//___________________________________________________
	//		Gene functions
	void	EvalDistances(int amountChrom);				// Evaluate goal function/distances
	double	CalcTotalDistance(int chromID);				// FITNESS FUNCTION || GOAL FUNCTION || total distance
	bool	AreChromsEqual(int chromID0, int chromID1);	// Check if chroms are equal;
	//___________________________________________________
	void	DoNextGen();								// Evaluate Next generation
	void	ReshapeGenePool();							// Sort by fitness; Check for similar chroms
	//___________________________________________________
	void	CopyParents();								// Create children from parents
	//void	RouletteSelection();						// Select parent by roulette wheel
	//___________________________________________________
	void	TransformChildren();						// Mutate and do MPC on children
	void	MPC(int chromID0, int chromID1);			// Multi Point Crossover
	void	MPC_RepairChrom(int chromId, int fromId, int toId);	// Repair the chrome after MPC
	//-------------------------------------------------------
};

#endif