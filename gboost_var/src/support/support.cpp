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
	if (can_prune(g2tracers)) {
		return;
	}

	map<int,PairSorter,greater<int> > b_heap;
	map<int,PairSorter,greater<int> > f_heap;
	map<Pair_IT, int> b_count;
	map<Pair_IT, int> f_count;

	int maxtoc = scan_gspan(g2tracers,b_heap,f_heap);

	// make sum of edge
	int sum_of_edge = 0;
	for (auto it = b_heap.begin(); it != b_heap.end(); ++it) {
		sum_of_edge += it->second.size();
	}
	for (auto it = f_heap.begin(); it != f_heap.end(); ++it) {
		sum_of_edge += it->second.size();
	}
	if (sum_of_edge == 0) {
		return;
	}

	// make count map
	vector<int> count(sum_of_edge, 0);
	make_count_s(n, count, b_heap, f_heap);

	// projecting...
	DFSCode  dcode;
	int c = 0;

	for (auto vitr = b_heap.begin(); vitr != b_heap.end(); ++vitr) {
		for (auto eitr = vitr->second.begin(); eitr != vitr->second.end(); ++eitr) {
			if (count[c] != 0) {
				dcode.labels = Triplet(-1,eitr->first.a,-1);
				dcode.time.set(vitr->first, eitr->first.b);
				pattern.push_back(dcode);
				edge_grow(eitr->second, count[c]);
				pattern.pop_back();
			}
			c++;
		}
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
	make_count_s(flownum, count, heap);
	
	int c = 0;
	for (auto it = heap.begin(); it != heap.end(); ++it) {
		//std::cout << support(it->second) << ":" << count[c] << std::endl;
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
	
	//TODO
	// select in opt_list
	if (opt_list.size() > 1) {
		std::uniform_int_distribution<unsigned int> rdm(0, opt_list.size() - 1);
		opt_pat = opt_list[rdm(mt)];
	}
}

void Gspan::make_count_s(unsigned int n, vector<int>& count, map<Triplet, GraphToTracers>& heap) {
	//vector<int> probabilities;
	vector<double> probabilities(count.size(), 0.0);
	unsigned int i = 0;
	for (auto it = heap.begin(); it != heap.end(); ++it) {
		probabilities[i] = support(it->second);
		i++;
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
}

void Gspan::make_count_s(unsigned int n, vector<int>& count, map<int, PairSorter, greater<int>>& b_heap, map<int, PairSorter, greater<int>>& f_heap) {
	vector<double> probabilities(count.size(), 0.0);
	unsigned int i = 0;
	for (auto vitr = b_heap.begin(); vitr != b_heap.end(); ++vitr) {
		for (auto eitr = vitr->second.begin(); eitr != vitr->second.end(); ++eitr) {
			probabilities[i] = support(eitr->second);
			i++;
		}
	}
	for (auto vitr = f_heap.begin(); vitr != f_heap.end(); ++vitr) {
		for (auto eitr = vitr->second.begin(); eitr != vitr->second.end(); ++eitr) {
			probabilities[i] = support(eitr->second);
			i++;
		}
	}

	std::discrete_distribution<int> rdm(probabilities.begin(), probabilities.end());
	for (i = 0; i < n; i++) {
		count[rdm(mt)]++;
	}
}
