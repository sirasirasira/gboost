#include "gspan.h"
#include <sstream>
#include <set>
#include <strstream>

bool Gspan::can_prune(GraphToTracers& g2tracers){
	vnum++;

	double gain=0.0;
	double upos=0.0;
	double uneg=0.0;

	gain=-wbias;
	upos=-wbias;
	uneg=wbias;

	for (GraphToTracers::iterator it=g2tracers.begin();it!=g2tracers.end();++it) {
		int gid = it->first;
		gain += 2 * corlab[gid] * weight[gid];
		if (corlab[gid]>0) {
			upos += 2 * weight[gid];
		} else {
			uneg += 2 * weight[gid];
		}
	}

	//std::cout << gain << " " ;
	
	bound = std::max(upos,uneg);
	gain_ = fabs(gain);
	//std::cout << max_gain << std::endl;
	
	if (fabs(opt_pat.gain) > bound ) {
		return true;
	}

	if ( gain_ > fabs(opt_pat.gain)
			|| (fabs(gain_ - fabs(opt_pat.gain)) < 1e-10 && pattern.size() < opt_pat.size)) {
		opt_pat.gain = gain;
		opt_pat.optimalplace = g2tracers;
		opt_pat.size = pattern.size();
		opt_pat.dfs = pattern;

		opt_list.resize(0);
		opt_list.push_back(opt_pat);
	} else if (gain_ == fabs(opt_pat.gain) && pattern.size() == opt_pat.size) {
		DPat tmp;
		tmp.gain = gain;
		tmp.optimalplace = g2tracers;
		tmp.size = pattern.size();
		tmp.dfs = pattern;

		opt_list.push_back(tmp);
	}

	return false;
}

vector<int> Gspan::update(DPat& op){
	vector<int> locvec; 
	for (GraphToTracers::iterator it=op.optimalplace.begin();it!=op.optimalplace.end();++it) {
		locvec.push_back(it->first);
	}

    std::ostrstream ostrs;
    ostrs << op.dfs;
    ostrs << std::ends;
    op.dfscode = ostrs.str();

	return locvec;
}
