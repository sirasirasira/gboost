#include "../gspan.h"
#include <random>

std::random_device rd;
std::mt19937 mt(rd());
//std::mt19937 mt(12345);

bool Gspan::prob_prune(Ctree& node) {
	double opt_gain = fabs(opt_pat.gain);
	double node_gain = node.max_gain;
	double prune_bound = opt_gain * prune_ratio;
	if (prune_bound + opt_gain <= node_gain) {
		return 0;
	}
	vector<double> probabilities{(node_gain - opt_gain) / prune_bound, 1 - ((node_gain - opt_gain) / prune_bound)};
	std::discrete_distribution<int> rdm(probabilities.begin(), probabilities.end());
	return rdm(mt);
}
