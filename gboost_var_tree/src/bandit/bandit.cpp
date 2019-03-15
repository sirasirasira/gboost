#include "../gspan.h"
#include <random>

using std::map; using std::vector;
std::random_device rd;
//std::mt19937 mt(rd());
std::mt19937 mt(12345);

void Gspan::edge_grow(GraphToTracers& g2tracers, int n) {
	if (pattern.size() > maxpat) {
		return;
	}
	if (support(g2tracers) < minsup) {
		return;
	}
	if (!is_min()) {
		return;
	}

	PairSorter b_heap;
	map<int,PairSorter,greater<int> > f_heap;

	int maxtoc = scan_gspan(g2tracers,b_heap,f_heap);

	// make sum of edge
	int sum_of_edge = 0;
	sum_of_edge += b_heap.size();
	for (auto it = f_heap.begin(); it != f_heap.end(); ++it) {
		sum_of_edge += it->second.size();
	}
	if (sum_of_edge == 0) {
		return;
	}

	// make count map
	vector<int> count(sum_of_edge, 0);
	if (make_count_b(n, count, b_heap, f_heap, maxtoc)) {
		return;
	}

	// projecting...
	DFSCode  dcode;
	int c = 0;

	for (auto eitr = b_heap.begin(); eitr != b_heap.end(); ++eitr) {
		if (count[c] != 0) {
			dcode.labels = Triplet(-1,eitr->first.b,-1);
			dcode.time.set(maxtoc, eitr->first.a);
			pattern.push_back(dcode);
			edge_grow(eitr->second, count[c]);
			pattern.pop_back();
		}
		c++;
	}
	for (auto vitr = f_heap.begin(); vitr != f_heap.end(); ++vitr) {
		for (auto eitr = vitr->second.begin(); eitr != vitr->second.end(); ++eitr) {
			if (count[c] != 0) {
				dcode.labels = Triplet(-1,eitr->first.a,eitr->first.b);
				dcode.time.set(vitr->first,maxtoc+1);
				pattern.push_back(dcode);
				edge_grow(eitr->second, count[c]);
				pattern.pop_back();
			}
			c++;
		}
	}
	return;
} 

void Gspan::sira_search() {
	map<Triplet,GraphToTracers> heap;
	for (unsigned int gid = 0; gid < gdata.size(); ++gid) {
		EdgeTracer cursor;
		Triplet wild_edge;
		Graph& g = gdata[gid];

		for (unsigned int v=0; v<g.size(); ++v) {
			for (vector<Edge>::iterator e = g[v].begin(); e != g[v].end(); ++e) {
				if (e->labels.x > e->labels.z) {
					continue;
				}
				cursor.set(v,e->to,e->id,0);
				heap[e->labels][gid].push_back(cursor);
				if (wildcard_r>0) {
					wild_edge = e->labels;
					wild_edge.z =999;
					heap[wild_edge][gid].push_back(cursor);
					wild_edge = e->labels.reverse();
					wild_edge.z =999;
					cursor.set(e->to,v,e->id,0);
					heap[wild_edge][gid].push_back(cursor);
				}
			}
		}
	}
	pattern.resize(1);

	// random edge select
	vector<int> count(heap.size(), 0);
	if (make_count_b(flownum, count, heap)) {
		return;
	}
	
	int c = 0;
	for (auto it = heap.begin(); it != heap.end(); ++it) {
		//std::cout << count[c] << std::endl;
		if (count[c] == 0) {
			c++;
			continue;
		}
		if (support(it->second) < minsup) {
			c++;
			continue;
		}
		pattern[0].labels = it->first;
		pattern[0].time.set(0,1);
		edge_grow(it->second, count[c]);
		pattern.resize(1);
		c++;
	}
	
	// select in opt_list
	if (opt_list.size() > 1) {
		std::uniform_int_distribution<unsigned int> rdm(0, opt_list.size() - 1);
		opt_pat = opt_list[rdm(mt)];
	}
}

double Gspan::valuation() {
	double value;
	value = lambda * gain_ + (1 - lambda) * bound;
	return value;
}

bool Gspan::make_count_b(unsigned int n, vector<int>& count, map<Triplet, GraphToTracers>& heap) {
	vector<double> probabilities(count.size(), 0.0);
	vector<double> bounds(count.size(), 0.0);
	unsigned int i = 0;
	bool prune = true;
	for (auto it = heap.begin(); it != heap.end(); ++it) {
		pattern[0].labels = it->first;
		pattern[0].time.set(0,1);
		can_prune(it->second);
		probabilities[i] = valuation();
		bounds[i] = bound;
		pattern.resize(1);
		i++;
	}

	for (i = 0; i < count.size(); i++) {
		if (bounds[i] < fabs(opt_pat.gain)) {
			probabilities[i] = 0.0;
		} else {
			prune = false;
		}
	}

	std::discrete_distribution<int> rdm(probabilities.begin(), probabilities.end());
	for (i = 0; i < n; i++) {
		count[rdm(mt)]++;
	} 
	/*
	for (auto itr = count.begin(); itr != count.end(); itr++) {
		std::cout << *itr << std::endl;
	}
	*/
	return prune;
}

bool Gspan::make_count_b(unsigned int n, vector<int>& count, PairSorter& b_heap, map<int, PairSorter, greater<int>>& f_heap, int maxtoc) {
	vector<double> probabilities(count.size(), 0.0);
	vector<double> bounds(count.size(), 0.0);
	unsigned int i = 0;
	bool prune = true;
	for (auto eitr = b_heap.begin(); eitr != b_heap.end(); ++eitr) {
		DFSCode  dcode;
		dcode.labels = Triplet(-1,eitr->first.b,-1);
		dcode.time.set(maxtoc, eitr->first.a);
		pattern.push_back(dcode);
		can_prune(eitr->second);
		probabilities[i] = valuation();
		bounds[i] = bound;
		pattern.pop_back();
		i++;
	}
	for (auto vitr = f_heap.begin(); vitr != f_heap.end(); ++vitr) {
		for (auto eitr = vitr->second.begin(); eitr != vitr->second.end(); ++eitr) {
			DFSCode  dcode;
			dcode.labels = Triplet(-1,eitr->first.a,eitr->first.b);
			dcode.time.set(vitr->first,maxtoc+1);
			pattern.push_back(dcode);
			can_prune(eitr->second);
			probabilities[i] = valuation();
			bounds[i] = bound;
			pattern.pop_back();
			i++;
		}
	}

	for (i = 0; i < count.size(); i++) {
		if (bounds[i] < fabs(opt_pat.gain)) {
			probabilities[i] = 0.0;
		} else {
			prune = false;
		}
	}

	std::discrete_distribution<int> rdm(probabilities.begin(), probabilities.end());
	for (i = 0; i < n; i++) {
		count[rdm(mt)]++;
	}
	return prune;
}
