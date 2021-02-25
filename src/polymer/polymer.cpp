//
#include "polymer.h"


polymer_t::polymer_t(double Mw, double Nb, int Nb_mono_deg, int NC, int NH, int NO, int NN, int NK, int NNa)
{
	Mn_=Mw_=Mw_poly_ = Mw;         //g / mol
	Mw_mono_ = NC * Mw_C_ + NH * Mw_H_ + NO * Mw_O_ + NK * Mw_K_ + NNa * Mw_Na_;
	// debug:
	//Mw_mono_ = 10.;
	Nb_mono_ = Nb;
	Nb_mono_deg_ = Nb_mono_deg;
	Nb_poly_ = (int)Mw_poly_ / Mw_mono_ * Nb_mono_ - 1;
	M_child2_ = ((double) (Nb_poly_+1)) * (Nb_poly_+1); //Molecular weights squared
	N_polyT_.push_back(Nb_poly_);
	M_polyT_.push_back(Mw_poly_);
}

void polymer_t::write_polymer_bindings()
{
	std::vector<int>::iterator it;
	int sum = 0;
	std::cout << "-----------------------------------------------" << std::endl;
	for (int i = 0; i < N_deg_monomers_; ++i)
		std::cout << 0 << std::endl;
	for (it = N_polyT_.begin(); it != N_polyT_.end(); ++it)
	{
		std::cout << *it << "\t" << std::endl;
		sum += *it;
	}

	std::cout << std::endl;
	std::cout << "total number of bindings " << sum << std::endl;
	std::cout << "-----------------------------------------------" << std::endl;
}

double polymer_t::degrade()
{
	//
	// pick a random position in a polymer backbone to be degraded
	// all elements less or equal to a monomer are removed
	
	//std::default_random_engine generator;
	//std::uniform_int_distribution<int> distribution(0,N_polyT_.size()-1);
	// pick a random degraded polymer,
	//int l = distribution(generator);
	int l = rand() % N_polyT_.size();
	// ... and place in that polymer
	//std::uniform_int_distribution<int> distribution2(1, N_polyT_[l]-1);
	//int l1 = distribution2(generator);
	int l1 = (N_polyT_[l]>1)?(rand() % (N_polyT_[l] - 1)+1):(1);
	int l2 = N_polyT_[l] - l1;
	N_polyT_[l] = l1 - 1;
	M_polyT_[l] = Mw_mono_ * l1;
	N_deg_step_ = 0;
	if (N_polyT_[l] < Nb_mono_deg_)
	{
		N_deg_step_ += 1;
		N_polyT_.erase(N_polyT_.begin() + l); //remove polymer from list
		M_polyT_.erase(M_polyT_.begin() + l); //remove polymer from list
	}
	if (l2 < Nb_mono_deg_)
	{
		N_deg_step_ += 1;
	}
	else
	{
		N_polyT_.push_back(l2);// # add to list
		M_polyT_.push_back(Mw_mono_*(l2+1.));
	}
	if (N_polyT_.empty()) // list empty, no polymer left
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



polymer_solution_t::polymer_solution_t(std::vector<polymer_t> polymer_types, int DT = 10000, int Tf = 1000000)
{
	polymer_types_ = polymer_types;
	DT_ = DT;
	time_end_ = Tf;


	std::vector<polymer_t>::iterator it;
	for (it = polymer_types_.begin(); it != polymer_types_.end(); ++it)
		active_.push_back(it-polymer_types_.begin());
}
double polymer_solution_t::degrade_polymer()
{
	//
	// pick a polymer in solution, degraded it and return the mass of the degraded part
	//
	//std::default_random_engine generator;
	//std::uniform_int_distribution<int> distribution(0, active_.size()-1);
	// pick a random degraded polymer,
	//int idx = distribution(generator);
	int idx=(active_.size()<2? 0:rand() % active_.size());
	polymer_t *pol;
	pol = &polymer_types_[active_[idx]];
	double dmp = pol->degrade();
	if (pol->fully_degraded_)
	{
		active_.erase(active_.begin() + idx);
		if(DEBUG_) std::cout << "polymer " << idx << " fully degraded" << std::endl;
	}
	return dmp * pol->Mw_mono_ / Na_;
}

void polymer_solution_t::degrade_polymer_solution()
{
	hist_fnames_.open("polymer_fnames.out", std::ofstream::out);
	while ((active_.size() > 0) && (clock_ < time_end_))
	{
		if (clock_ == 0 || clock_ % DT_ == 0)
		{

			
			get_polymer_fractions();
			write_polymer_hist();
			std::cout << "Molar weight avereage: "<<Mn_.back() << "\t" << Mw_.back() << "\t" << Mw_.back() / Mn_.back() << std::endl;
			std::cout << "time step " << clock_ << std::endl;
		}
		mass_degraded_t_.push_back(degrade_polymer());
		clock_ += 1;
	}
	hist_fnames_.close();
}

void polymer_solution_t::write_polymer_hist()
{
	std::ofstream hist_out;
	std::string fname = "hist" + std::to_string(clock_) + ".out";
	hist_fnames_ << fname << std::endl;
	hist_out.open(fname, std::ofstream::out);
	std::vector<int>::iterator it;
	for (it = hist_.begin(); it != hist_.end(); ++it)
		hist_out << *it << "\n";
	hist_out.close();
}
void polymer_solution_t::get_polymer_fractions()
{
	// Have to calculate the average of each polymer molecules
	std::vector<int> Ncc_dist;
	std::vector<polymer_t>::iterator it;

	double Mw = 0.;
	double Mn = 0.;
	double Mw_norm = 0.;
	double Mn_norm = 0.;
	double No_cut = 0.;
	double No_cut_norm = 0.;

	Ncc_dist.reserve(polymer_types_.size());
	for (it = polymer_types_.begin(); it < polymer_types_.end(); it++)
	{
		if (it->N_deg_monomers_ > 0)
		{
			Ncc_dist.push_back(it->N_deg_monomers_);
		}
		Ncc_dist.insert(Ncc_dist.end(), it->N_polyT_.begin(), it->N_polyT_.end());
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

	Mn_.push_back(Mn / Mn_norm);
	Mw_.push_back(Mw / Mw_norm);
	Tn_.push_back(clock_);
	r_.push_back(No_cut/ No_cut_norm);
	

	hist_= Ncc_dist;
}


void polymer_solution_t::write_time_series()
{
	std::ofstream off;
	off.open("polymer_t.out");
	off << "time\tMn\tMw\tPD\tFractionCut\n";
	for (int i = 0; i < Tn_.size(); ++i)
		off << Tn_[i] << "\t" << Mn_[i] << "\t" << Mw_[i] << "\t" << Mw_[i]/Mn_[i]<<"\t"<<r_[i]<<"\n";
	off.close();
}
//std::default_random_engine generator;
//std::uniform_int_distribution<int> distribution(1, 6);
//int dice_roll = distribution(generator);  // generates number in the range 1..6 
