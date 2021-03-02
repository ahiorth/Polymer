//
#include "../polymer/polymer.h"

int main(int argc, char* argv[])
{
	double Mw1 = 6e6;
	int	N1 = 1000;
	int	Tf = 1000000000;
	int	DT = 100000;

	std::vector<polymer_t> a;
	for (int i = 0; i < N1; ++i)
		a.push_back(polymer_t(Mw1));

	MC_sampling_t MC_run(a, DT, Tf);
	MC_run.simulate();
	MC_run.get_average();

//	polymer_solution_t sol(a, DT, Tf);
//	sol.degrade_polymer_solution();
//	sol.write_time_series();

};
