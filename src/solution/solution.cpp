//
#include "../polymer/polymer.h"

int main(int argc, char* argv[])
{
	double Mw1 = 6e3;
	int	N1 = 10000;
	int	Tf = 1000000;
	int	DT = 1000;

	std::vector<polymer_t> a;
	for (int i = 0; i < N1; ++i)
		a.push_back(polymer_t(Mw1));

	//double Mw2 = 1e3;
	//int N2 = (int)(Mw1 / Mw2 * N1);
	//std::vector<polymer_t> a;
	//for (int i = 0; i < N2; ++i)
	//	a.push_back(polymer_t(Mw2));

	polymer_solution_t sol(a, DT, Tf);

	sol.degrade_polymer_solution();
	sol.write_time_series();
	std::cout << "here is solution" << std::endl;
};