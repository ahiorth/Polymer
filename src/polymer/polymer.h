#ifndef POLYMER_H
#define POLYMER_H

#include <vector>
#include <random>
#include <string>
#include <iostream>
#include <fstream>
#include "omp.h"

class polymer_t
{
	//
	//solver class for calculating polymer degradation NC, ..., NNa is the number of Carbon,
	//hydrogen, oxygen, nitrogen, kalium, and natrium in the polumer monomer
	
	double Mw_C_ = 12.0107;     //g / mol
	double Mw_H_ = 1.0079;      //g / mol
	double Mw_O_ = 15.9994;     //g / mol
	double Mw_N_ = 28.014;      //g / mol
	double Mw_K_ = 39.098;      //g / mol
	double Mw_Na_ = 22.989769;  //g / mol

public:
	polymer_t(double Mw=6e6, double Nb = 1, int Nb_mono_deg = 1, int NC = 3, int NH = 5, int NO = 1, int NN = 1, int NK = 0, int NNa = 0);
	double Mw_poly_;
	double Mw_mono_;
	int    Nb_mono_;
	int    Nb_mono_deg_;
	int    Nb_poly_; // Number of bindings in original polymer
	double Mn_; // number average molecular weight
	double Mw_; // mass average molecular weight
	double PD_=1;
	int N_deg_monomers_ = 0;    // cummulative number of degraded monomers
	int N_deg_step_ = 0;        // number of degraded per step
	bool fully_degraded_ = false;
	int No_childs_ = 1;
	double M_child2_; // int overflows ....
	std::vector<int> idx_; // length of polymer chain with a change
	std::vector<int> dn_;  // change in chain length

	std::pair<int, int> pos_;

	std::vector<int> cum_num_bindings_;

	double degrade(int l);
	void write_polymer_bindings();
};

class polymer_solution_t
{
	//	takes as input an array of polymer objects, and degrades them
	//	polymer_types : list of polymer objects
	//	no_polymer_types : number of each polymer in the solution
public:
	polymer_solution_t(std::vector<polymer_t> polymer_types, int DT, int Tf);
	double Na_ = 6.0221409e+23; //Avogadros number
	int ACTIVE_ = 1;
	std::vector<polymer_t> polymer_types_;
	std::vector<int> active_;
	int time_end_;
	int DT_;
	int clock_ = 0;
	int DEBUG_ = 0;
	int MAKE_HIST_ = 0; // write histograms file - memory intensive
	int MAX_POLYMER_LENGTH_;
	std::vector<double> mass_degraded_t_;
	std::vector<double> r_;  // average number  of bounds broken
	std::vector<double> Mn_; // number  average molecular weight
	std::vector<double> Mw_; // mass    average molecular weight
	std::vector<int> Tn_;
	std::vector<int> cum_num_bindings_; 
	std::vector<std::vector<int>> hist_;
	double degrade_polymer();
	void degrade_polymer_solution();
	void write_polymer_hist(std::vector<int> &hist, int clock);
	std::ofstream hist_fnames_;
	std::vector<int> Ndist_;    // Number distribution
	std::vector<std::vector<int>> NdistT_; // number distrubtion at DT intervals
	std::pair<int, int> pos_;
	void get_polymer_fractions();
	void write_time_series();
	void write_time_denst(std::vector<int> &hist, int clock);

};

class MC_sampling_t
{
public:
	MC_sampling_t(std::vector<polymer_t> polymer_types, int DT, int Tf);
	const int NO_MC_RUNS_ = 100;
	int MAX_TIME_ = 0; 
	std::vector<polymer_solution_t> solMC_; 
	polymer_solution_t sol_FINAL_; // collect average values
	void simulate();
	void get_average();
	void get_max_time();

};

#endif