#include "../gspan.h"
#include <random>

std::random_device rd;
std::mt19937 mt(rd());
//std::mt19937 mt(12345);

bool Gspan::prob_prune(Ctree& node) {
	vector<double> probabilities{1 - prune_ratio, prune_ratio};
	std::discrete_distribution<int> rdm(probabilities.begin(), probabilities.end());
	return rdm(mt);
}
