#ifndef POLYMER_H
#define POLYMER_H

#include <vector>
#include <random>
#include <string>
#include <iostream>
#include <fstream>

class polymer_t
{
	//
	//solver class for calculating polymer degradation NC, ..., NNa is the number of Carbon,
	//hydrogen, oxygen, nitrogen, kalium, and natrium in the polumer monomer
	polymer_t(double Mw, double Nb, int Nb_mono_deg, int NC, int NH, int NO, int NN, int NK, int NNa);

	double Mw_C_ = 12.0107;     //g / mol
	double Mw_H_ = 1.0079;      //g / mol
	double Mw_O_ = 15.9994;     //g / mol
	double Mw_N_ = 28.014;      //g / mol
	double Mw_K_ = 39.098;      //g / mol
	double Mw_Na_ = 22.989769;  //g / mol

public:
	double Mw_poly_;
	double Mw_mono_;
	int    Nb_mono_;
	int    Nb_mono_deg_;
	int    Nb_poly_;
	std::vector<int> N_polyT_; // vector of the individual length of polymer chain after time T
	int N_deg_monomers_ = 0;    // cummulative number of degraded monomers
	int N_deg_step_ = 0;        // number of degraded per step
	bool fully_degraded_ = false;

	double degrade();
};

class polymer_solution_t
{
	//	takes as input an array of polymer objects, and degrades them
	//	polymer_types : list of polymer objects
	//	no_polymer_types : number of each polymer in the solution
	double Na_ = 6.0221409e+23; //Avogadros number
	int ACTIVE_ = 1;
	int no_polymers_;
	std::vector<polymer_t> polymer_types_;
	std::vector<int> active_;
	int time_end_;
	int DT_;
	int clock_ = 0;
	std::vector<double> mass_degraded_t_;
	std::vector<double> Mn_;// number average molecular weight
	std::vector<double> Mw_;// mass average molecular weight
	std::vector<int> Tn_;
	double degrade_polymer();
	double degrade_polymer_solution();
	void write_polymer_hist();
	std::ofstream hist_fnames_;
	
	polymer_solution_t(std::vector<polymer_t> polymer_types, int DT, int Tf);
	std::vector<double> get_polymer_fractions();

};

#endif