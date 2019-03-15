#include "../gspan.h"
#include <random>

std::random_device rd;
std::mt19937 mt(rd());
//std::mt19937 mt(12345);

bool Gspan::prob_prune(Ctree& node) {
	double prune_support = gdata.size() * prune_ratio;
	if (prune_support <= support(node.g2tracers)) {
		return 0;
	}
	vector<double> probabilities{support(node.g2tracers) / prune_support, 1 - (support(node.g2tracers) / prune_support)};
	std::discrete_distribution<int> rdm(probabilities.begin(), probabilities.end());
	return rdm(mt);
}
