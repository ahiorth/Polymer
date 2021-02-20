//
#include "polymer.h"


polymer_t::polymer_t(double Mw = 6e6, double Nb = 1, int Nb_mono_deg = 1, int NC = 3, int NH = 5, int NO = 1, int NN = 1, int NK = 0, int NNa = 0)
{
	Mw_poly_ = Mw;         //g / mol
	Mw_mono_ = NC * Mw_C_ + NH * Mw_H_ + NO * Mw_O_ + NK * Mw_K_ + NNa * Mw_Na_;
	Nb_mono_ = Nb;
	Nb_mono_deg = Nb_mono_deg;
	Nb_poly_ = (int)Mw_poly_ / Mw_mono_ * Nb_mono_ - 1;
	N_polyT_.push_back(Nb_poly_);
}

double polymer_t::degrade()
{
	//
	// pick a random position in a polymer backbone to be degraded
	// all elements less or equal to a monomer are removed
	
	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution(0,N_polyT_.size());
	// pick a random degraded polymer,
	int l = distribution(generator);
	// ... and place in that polymer
	std::uniform_int_distribution<int> distribution2(0, N_polyT_[l]);
	int l1 = distribution2(generator);
	int l2 = N_polyT_[l] - l1;
	N_polyT_[l] = l1 - 1;
	N_deg_step_ = 0;
	if (N_polyT_[l] <= Nb_mono_deg_)
	{
		N_deg_step_ += l1;
		N_polyT_.erase(N_polyT_.begin() + l); //remove polymer from list
		if (l2 <= Nb_mono_deg_)
		{
			N_deg_step_ += l;
		}
		else
		{
			N_polyT_.push_back(l2);// # add to list
		}
		if (N_polyT_.empty()) // list empty, no polymer left
		{
			fully_degraded_ = true;
		}
		N_deg_monomers_ += N_deg_step_;
	}
	return N_deg_step_;
}

polymer_solution_t::polymer_solution_t(std::vector<polymer_t> polymer_types, int DT = 10000, int Tf = 1000000)
{
	polymer_types_ = polymer_types;
	DT_ = DT;
	time_end_ = Tf;
	hist_fnames_.open("polymer_fnames.out");
}
double polymer_solution_t::degrade_polymer()
{
	//
	// pick a polymer in solution, degraded it and return the mass of the degraded part
	//
	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution(0, active_.size());
	// pick a random degraded polymer,
	int idx = distribution(generator);
	polymer_t *pol;
	pol = &polymer_types_[active_[idx]];
	double dmp = pol->degrade();
	if (pol->fully_degraded_)
	{
		active_.erase(active_.begin() + idx);
		std::cout << "polymer " << idx << " fully degraded" << std::endl;
	}
	return dmp * pol->Mw_mono_ / Na_;
}

double polymer_solution_t::degrade_polymer_solution()
{
	while (active_.size() > 0 & clock_ < time_end_)
	{
		if (clock_ == 0 || clock_ % DT_ == 0)
		{
			write_polymer_hist();
		}
		mass_degraded_t_.push_back(degrade_polymer());
		clock_ += 1;
	}
}

void polymer_solution_t::write_polymer_hist()
{
	std::ofstream hist_out;
	std::string fname = "hist" + std::to_string(clock_) + ".out";
	hist_fnames_ << fname << std::endl;
	hist_out.open(fname);
}

std::vector<double> polymer_solution_t::get_polymer_fractions()
{
	std::vector<double> Ncc_dist;
	double Mw, Mn;
	std::vector<polymer_t>::iterator it;
	int no = 0;
	Mw = Mn = 0;
	Ncc_dist.reserve(polymer_types_.size());
	for (it = polymer_types_.begin(); it < polymer_types_.end(); it++)
	{
		if (it->N_deg_monomers_ > 0)
			Ncc_dist.push_back(it->N_deg_monomers_);
		Ncc_dist.insert(Ncc_dist.end(), it->N_polyT_.begin(), it->N_polyT_.end());
		double mn = 0;
		for (auto& n : it->N_polyT_)

			mn += n;
		Mn += mn * it->Mw_mono_ / it->N_polyT_.size();
		Mw += mn * it->Mw_mono_*it->Mw_mono_ / it->N_polyT_.size();
	}

	Mn

		Mn += sum(mm)*pol.Mw_mono / len(mm)
		self.Mn.append(Mn / len(self.polymer_types))
		self.Tn.append(self.clock)
		return Mw
}
//std::default_random_engine generator;
//std::uniform_int_distribution<int> distribution(1, 6);
//int dice_roll = distribution(generator);  // generates number in the range 1..6 
