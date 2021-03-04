//
#include "polymer.h"
#include "util.h"


polymer_t::polymer_t(double Mw, double Nb, int Nb_mono_deg, int NC, int NH, int NO, int NN, int NK, int NNa)
{
	Mn_=Mw_=Mw_poly_ = Mw;         //g / mol
	Mw_mono_ = NC * Mw_C_ + NH * Mw_H_ + NO * Mw_O_ + NK * Mw_K_ + NNa * Mw_Na_;
	// debug:
	//Mw_mono_ = 10.;
	Nb_mono_ = Nb;
	Nb_mono_deg_ = Nb_mono_deg;
	Nb_poly_ = (int)Mw_poly_ / Mw_mono_ * Nb_mono_ - 1;
	std::cout << "No of bonds " << Nb_poly_ << std::endl;
	M_child2_ = ((double) (Nb_poly_+1.)) * (Nb_poly_+1.); //Molecular weights squared

	idx_.resize(3, 0); // one can only lead to three changes
	dn_.push_back(-1); dn_.push_back(1); dn_.push_back(1);
	cum_num_bindings_.push_back(Nb_poly_);

}

void polymer_t::write_polymer_bindings()
{
	std::vector<int>::iterator it;
	int sum = 0;
	std::cout << "-----------------------------------------------" << std::endl;
	for (int i = 0; i < N_deg_monomers_; ++i)
		std::cout << 0 << std::endl;
	int prev = 0;
	for (it = cum_num_bindings_.begin(); it != cum_num_bindings_.end(); ++it)
	{
		std::cout << *it-prev << "\t" << std::endl;
		sum += (*it-prev);
		prev = *it;
	}

	std::cout << std::endl;
	std::cout << "total number of bindings " << sum << std::endl;
	std::cout << "-----------------------------------------------" << std::endl;
}

double polymer_t::degrade(int idx)
{
	// degrade bindings number l
	// all elements less or equal to a monomer are removed
	pos_ = findInRange(cum_num_bindings_, idx);
	int l    = pos_.first; // position in pos_
	int size = pos_.second; // chain_length
	int l2 = cum_num_bindings_[l] - idx;
	int l1 = size - l2;
	
	// polymer molecule of size size is split into l2 and l1-1
	idx_[0] = size; idx_[1] = l1 - 1; idx_[2] = l2;
	for (int i = l; i < cum_num_bindings_.size(); ++i)
		cum_num_bindings_[i]-=l2+1;

	N_deg_step_ = 0;
	if (l1 - 1 < Nb_mono_deg_)
	{
		N_deg_step_ += 1;
		cum_num_bindings_.erase(cum_num_bindings_.begin() + l); //remove polymer from list
	}
	else
	{
		int prev = (l == 0) ? (0):(cum_num_bindings_[l - 1]);
		cum_num_bindings_[l] = prev+l1-1; 
	}
	if (l2 < Nb_mono_deg_)
	{
		N_deg_step_ += 1;
	}
	else
	{
		int back = (cum_num_bindings_.size() > 0) ? (cum_num_bindings_.back()) : (0);
		cum_num_bindings_.push_back(l2 + back);
	}

	
	if (cum_num_bindings_.empty()) // list empty, no polymer left
	{
		fully_degraded_ = true;
	}
	
	N_deg_monomers_ += N_deg_step_;
	No_childs_++;
	M_child2_ -= 2. * l1 * (l2 + 1);

	// update number average:
	Mw_ = M_child2_ * Mw_mono_ * Mw_mono_ / Mw_poly_;
	Mn_ = Mw_poly_ / ((double)No_childs_);
	return N_deg_step_;
}



polymer_solution_t::polymer_solution_t(std::vector<polymer_t> polymer_types, int DT, int Tf)
{
	polymer_types_ = polymer_types;
	DT_ = DT;
	time_end_ = Tf;
	std::vector<polymer_t>::iterator it;
	
	MAX_POLYMER_LENGTH_ = polymer_types_[0].Nb_poly_;
	for (int i = 1; i < polymer_types_.size(); ++i)
		MAX_POLYMER_LENGTH_ = std::max(polymer_types_[i-1].Nb_poly_, polymer_types_[i].Nb_poly_);

	Ndist_.resize(MAX_POLYMER_LENGTH_+1, 0);
	int back = 0;
	for (it = polymer_types_.begin(); it != polymer_types_.end(); ++it)
	{
		active_.push_back(it - polymer_types_.begin());
		Ndist_[it->Nb_poly_]++;
		cum_num_bindings_.push_back(back+it->Nb_poly_);
		back = cum_num_bindings_.back();
	}
}
double polymer_solution_t::degrade_polymer()
{
	//
	// pick a chain in solution, degraded it and return the mass of the degraded part
	//
	//std::default_random_engine generator;
	//std::uniform_int_distribution<int> distribution(0, active_.size()-1);
	// pick a random degraded polymer,
	//int idx = distribution(generator);


	int chain_no = (cum_num_bindings_.back() <2)?(1):(rand() % (cum_num_bindings_.back()-1)+1);
	pos_ = findInRange(cum_num_bindings_, chain_no);
	int idx = pos_.first; // position of polymer
	polymer_t *pol;
	pol = &polymer_types_[idx];
	if (DEBUG_)
	{
		std::cout << "break binding " << chain_no << std::endl;
	}
	int prev = (idx == 0) ? (0) : (cum_num_bindings_[idx - 1]);
	double dmp = pol->degrade(chain_no- prev);
	
	Ndist_[pol->idx_[0]] += pol->dn_[0];
	Ndist_[pol->idx_[1]] += pol->dn_[1];
	Ndist_[pol->idx_[2]] += pol->dn_[2];


	if (pol->fully_degraded_)
	{
		if(DEBUG_) std::cout << "polymer " << idx << " fully degraded" << std::endl;
	}
	for (int i = idx; i < cum_num_bindings_.size(); ++i) // remove a chain 
		cum_num_bindings_[i] --;
	return dmp * pol->Mw_mono_ / Na_;
}

void polymer_solution_t::degrade_polymer_solution()
{
	hist_fnames_.open("polymer_fnames.out", std::ofstream::out);
	while ((cum_num_bindings_.back() > 0) && (clock_ < time_end_))
	{
		if (clock_ == 0 || clock_ % DT_ == 0)
		{

			
			get_polymer_fractions();
			if (DEBUG_)std::cout << "Molar weight avereage: "<<Mn_.back() << "\t" << Mw_.back() << "\t" << Mw_.back() / Mn_.back() << "\t" << r_.back()<< std::endl;
			if (DEBUG_)std::cout << "time step " << clock_ << std::endl;
			NdistT_.push_back(Ndist_);
		}
		
		clock_ += 1;
		mass_degraded_t_.push_back(degrade_polymer());
		
	}
	hist_fnames_.close();
}

void polymer_solution_t::write_polymer_hist(std::vector<int> &hist, int clock)
{
	std::ofstream hist_out;
	std::string fname = "hist" + std::to_string(clock) + ".out";
	hist_fnames_ << fname << std::endl;
	hist_out.open(fname, std::ofstream::out);
	std::vector<int>::iterator it;
	for (it = hist.begin(); it != hist.end(); ++it)
		hist_out << *it << "\n";
	hist_out.close();
}

void polymer_solution_t::write_time_denst(std::vector<int> &Ndist, int clock)
{
	std::ofstream hist_out;
	std::string fname = "denst" + std::to_string(clock) + ".out";
	hist_out.open(fname, std::ofstream::out);
	
	for (int i=0;i<MAX_POLYMER_LENGTH_+1;++i)
		hist_out << Ndist[i]<<"\n";
	hist_out.close();
}
void polymer_solution_t::get_polymer_fractions()
{
	// Have to calculate the average of each polymer molecules
	std::vector<int> Ncc_dist;
	std::vector<int>::iterator iti;
	std::vector<polymer_t>::iterator it;

	double Mw = 0.;
	double Mn = 0.;
	double Mw_norm = 0.;
	double Mn_norm = 0.;
	double No_cut = 0.;
	double No_cut_norm = 0.;

	if(MAKE_HIST_) Ncc_dist.reserve(polymer_types_.size());
	for (it = polymer_types_.begin(); it != polymer_types_.end(); it++)
	{
		if (it->N_deg_monomers_ > 0)
		{
			for(int j=0;j<it->N_deg_monomers_;++j)
				Ncc_dist.push_back(0);
		}
		if (MAKE_HIST_)
		{
			int prev = 0;
			for (iti = it->cum_num_bindings_.begin(); iti != it->cum_num_bindings_.end(); iti++)
			{
				Ncc_dist.push_back(*iti - prev);
				prev = *iti;
			}
		}
	
		int index = std::distance(polymer_types_.begin(), it);
		if (DEBUG_) it->write_polymer_bindings();
		if (DEBUG_)std::cout << "pol"<< index<<"\t"<<it->Mn_ << "\t" << it->Mw_ << "\t" << it->Mw_ / it->Mn_<<"\n";
		Mn += it->Mn_*(it->No_childs_);
		Mw += it->Mw_* (it->Mw_poly_);
		Mn_norm += (double) it->No_childs_;
		Mw_norm += (double) it->Mw_poly_;
		No_cut += it->No_childs_-1.;
		No_cut_norm += (double) it->Nb_poly_;
	}
	// Hack calculate from distributions
	/*double y = 0.;
	double y_norm = 0.;
	std::vector<double>w;
	double mx=0.,w_norm=0.;
	for (int i = 0; i < Ncc_dist.size(); ++i)
	{
		y_norm += (1. + Ncc_dist[i]);
		y += (1.+Ncc_dist[i]) * (1.+Ncc_dist[i]);
	}
	for (int i = 0; i < Ncc_dist.size(); ++i)
	{
		w.push_back((1. + Ncc_dist[i]) / y_norm);
	}

	for (int i = 0; i < Ncc_dist.size(); ++i)
	{
		mx += w[i] * (1. + Ncc_dist[i]);
	}
	mx *= polymer_types_[0].Mw_mono_;
	Mw_.push_back(mx);


	y = y / y_norm * polymer_types_[0].Mw_mono_;
	*/
	Mn_.push_back(Mn / Mn_norm);
	Mw_.push_back(Mw / Mw_norm);
	Tn_.push_back(clock_);
	r_.push_back(No_cut/ No_cut_norm);
	

	if (MAKE_HIST_)hist_.push_back(Ncc_dist);
}


void polymer_solution_t::write_time_series()
{
	std::ofstream off;
	off.open("polymer_t.out");
	off << "time\tMn\tMw\tPD\tFractionCut\n";
	for (int i = 0; i < Tn_.size(); ++i)
	{
		off << Tn_[i] << "\t" << Mn_[i] << "\t" << Mw_[i] << "\t" << Mw_[i] / Mn_[i] << "\t" << r_[i] << "\n";
		std::cout << " Time " << i << " finnished!" << std::endl;
		
	}
	off.close();
	if (MAKE_HIST_)
	{
#pragma omp parallel for
		for (int i = 0; i < Tn_.size(); ++i)
		{
			write_polymer_hist(hist_[i], (int)Tn_[i]);
			std::cout << " Hist " << i << " finnished!" << std::endl;
		}
	}
#pragma omp parallel for
	for (int i = 0; i < Tn_.size(); ++i)
	{
		write_time_denst(NdistT_[i], (int)Tn_[i]);
		std::cout << " Denst " << i << " finnished!" << std::endl;
	}
	std::cout << "Number of bonds " << polymer_types_[0].Nb_poly_ << std::endl;
}

MC_sampling_t::MC_sampling_t(std::vector<polymer_t> a, int DT, int Tf):sol_FINAL_(a,DT,Tf)
{
#pragma omp parallel for
	for (int i = 0; i < NO_MC_RUNS_; ++i)
	{
		solMC_.push_back(polymer_solution_t(a, DT, Tf));
	}
}

void MC_sampling_t::simulate()
{
#pragma omp parallel for
	for (int i = 0; i < NO_MC_RUNS_; ++i)
	{
		std::cout << " Mc no" << i << " started " << std::endl;
		solMC_[i].degrade_polymer_solution();
		std::cout << " MC no " << i << " finnished!" << std::endl;
	}
	get_max_time();
}

void MC_sampling_t::get_max_time()
{
	MAX_TIME_ = solMC_[0].Tn_.size();
	for (int i = 1; i < NO_MC_RUNS_; ++i)
		MAX_TIME_ = std::min(solMC_[i - 1].Tn_.size(), solMC_[i].Tn_.size());
}
void MC_sampling_t::get_average()
{
	std::vector<double> Mn2, Mw2, r2;
	double N = (double)NO_MC_RUNS_;
	Mn2.resize(MAX_TIME_, 0.);
	Mw2.resize(MAX_TIME_, 0.);
	r2.resize(MAX_TIME_, 0.);
	sol_FINAL_.Tn_.resize(MAX_TIME_, 0.);
	sol_FINAL_.Mn_.resize(MAX_TIME_, 0.);
	sol_FINAL_.Mw_.resize(MAX_TIME_, 0.);
	sol_FINAL_.r_.resize(MAX_TIME_, 0.);
	sol_FINAL_.NdistT_.resize(MAX_TIME_);
	if (sol_FINAL_.MAKE_HIST_)
		sol_FINAL_.hist_.resize(MAX_TIME_);

	for (int t = 0; t < MAX_TIME_; ++t)
	{
		Mn2[t] = Mw2[t] = r2[t] = sol_FINAL_.Mn_[t] = sol_FINAL_.Mw_[t] = sol_FINAL_.r_[t] = 0.;
		sol_FINAL_.NdistT_[t].resize(sol_FINAL_.MAX_POLYMER_LENGTH_ + 1, 0);
		sol_FINAL_.Tn_[t] = solMC_[0].Tn_[t];
		for (int i = 0; i < NO_MC_RUNS_; ++i)
		{
			sol_FINAL_.Mn_[t] += solMC_[i].Mn_[t];
			sol_FINAL_.Mw_[t] += solMC_[i].Mw_[t];
			sol_FINAL_.r_[t] += solMC_[i].r_[t];
			Mn2[t] += solMC_[i].Mn_[t] * solMC_[i].Mn_[t];
			Mw2[t] += solMC_[i].Mw_[t] * solMC_[i].Mw_[t];
			r2[t] += solMC_[i].r_[t] * solMC_[i].r_[t];
	
			if(sol_FINAL_.MAKE_HIST_)
				sol_FINAL_.hist_[t].insert(sol_FINAL_.hist_[t].end(), solMC_[i].hist_[t].begin(), solMC_[i].hist_[t].end());
			for (int no = 0; no < solMC_[i].NdistT_[t].size(); ++no)
			{
				sol_FINAL_.NdistT_[t][no] += solMC_[i].NdistT_[t][no];
			}
		}
		sol_FINAL_.Mn_[t] /= N;
		sol_FINAL_.Mw_[t] /= N;
		sol_FINAL_.r_[t]  /= N;
		Mn2[t] = Mn2[t] / N - sol_FINAL_.Mn_[t];
		Mw2[t] = Mw2[t] / N - sol_FINAL_.Mw_[t];
		r2[t]  = r2[t] / N - sol_FINAL_.r_[t];
		std::cout << " Time " << t << " finnished!" << std::endl;
	}
	sol_FINAL_.write_time_series();
}
	
