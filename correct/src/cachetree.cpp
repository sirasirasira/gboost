#include "gspan.h"
#include <sstream>
#include <set>
#include <strstream>

void Ctree::print(){
	std::cout << g2tracers.size() << ":" << pat.rebuild() <<std::endl;
	if(children.size()!=0){
		for(list<Ctree*>::iterator it=children.begin();it!=children.end();++it){
			std::cout << "-";
			(*it)->print();
		}
	}
}
void CRoot::print(){
	for(list<Ctree*>::iterator it=one_edge_graphs.begin();it!=one_edge_graphs.end();++it){
		(*it)->print();
	}
}
void Tdelete(Ctree* tree){
	if(tree->children.size()!=0){
		for(list<Ctree*>::iterator it=tree->children.begin();it!=tree->children.end();++it){
			Tdelete(*it);
		}
	}
	delete tree;
}
void Rdelete(CRoot* croot){
	if(croot->one_edge_graphs.size()!=0){
		for(list<Ctree*>::iterator it=croot->one_edge_graphs.begin();it!=croot->one_edge_graphs.end();++it){
			Tdelete(*it);
		}
	}
	delete croot;
	//std::cout << "delete complete ! "<< std::endl;
}

void Gspan::Crun(){
	if(first_flag==true){
		first_tree_make();
		first_flag=false;
		std::cout <<"CacheTree Node size : " << TNnum <<std::endl;
		return;
	}
	CacheTree_search();
	std::cout <<"CacheTree Node size : " << TNnum <<std::endl;
}

void Gspan::first_tree_make(){
	/***   init CRoot         ***/
	croot = new CRoot;
	TNnum = 1;
	/****  construct CRoot   ****/
	map<Triplet,GraphToTracers> heap;
	for(unsigned int gid = 0; gid < gdata.size(); ++gid){
		EdgeTracer cursor;
		Triplet wild_edge;
		Graph& g = gdata[gid];

		for(unsigned int v=0; v<g.size(); ++v){
			for(vector<Edge>::iterator e = g[v].begin(); e != g[v].end(); ++e){
				if (e->labels.x <= e->labels.z){
					cursor.set(v,e->to,e->id,0);
					heap[e->labels][gid].push_back(cursor);
				}
			}
		}
	}
	pattern.resize(1);
	for(map<Triplet,GraphToTracers>::iterator it = heap.begin(); it != heap.end(); ++it){
		if(support(it->second) < minsup){ continue;}
		pattern[0].labels = it->first;
		pattern[0].time.set(0,1);

		Ctree *node = new Ctree;
		TNnum++;
		nDFSCode d;
		d.set(pattern[0],NULL);
		node->pat = d;
		GraphToTracers& g2 = node->g2tracers;
		g2 = it->second;
		node->children.resize(0);
		croot->one_edge_graphs.push_back(node);
		edge_grow(*node);
		pattern.resize(1);
	}
}

void Gspan::edge_grow(Ctree& nnode){
	if(can_prune(nnode)) { return;}

	if(pattern.size() == maxpat){return;}

	PairSorter b_heap;
	map<int,PairSorter,greater<int> > f_heap;
	int maxtoc = scan_gspan(nnode.g2tracers,b_heap,f_heap);

	// projecting...
	DFSCode  dcode;
	for(PairSorter::iterator it = b_heap.begin(); it != b_heap.end(); ++it){	
		dcode.labels = Triplet(-1,it->first.b,-1);
		dcode.time.set(maxtoc, it->first.a);
		pattern.push_back(dcode);

		if(support(it->second) < minsup){pattern.pop_back(); continue;}
		if(!is_min()){pattern.pop_back(); continue;}

		//new Node projecting ...
		Ctree *node = new Ctree;
		TNnum++;
		nDFSCode d;
		d.set(dcode,&nnode.pat);
		node->pat = d;
		GraphToTracers& g2 = node->g2tracers;
		g2 = it->second;
		nnode.children.push_back(node);

		edge_grow(*node);
		pattern.pop_back();
	}

	for(map<int,PairSorter,greater<int> >::iterator it = f_heap.begin(); it != f_heap.end(); ++it){	
		for(PairSorter::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2){		
			dcode.labels = Triplet(-1,it2->first.a,it2->first.b);
			dcode.time.set(it->first,maxtoc+1);
			pattern.push_back(dcode);

			if(support(it2->second) < minsup){pattern.pop_back(); continue;}
			if(!is_min()){pattern.pop_back(); continue;}

			//new Node projecting ...
			Ctree *node = new Ctree;
			TNnum++;
			nDFSCode d;
			d.set(dcode,&nnode.pat);
			node->pat = d;

			GraphToTracers& g2 = node->g2tracers;
			g2 = it2->second;
			nnode.children.push_back(node);

			edge_grow(*node);
			pattern.pop_back();
		}
	}
}

bool Gspan::can_prune(Ctree& node){
	vnum++;

	double gain=0.0;
	double upos=0.0;
	double uneg=0.0;

	gain=-wbias;
	upos=-wbias;
	uneg=wbias;

	for(GraphToTracers::iterator it=node.g2tracers.begin();it!=node.g2tracers.end();++it){
		int gid = it->first;
		gain += 2 * corlab[gid] * weight[gid];
		if(corlab[gid]>0){
			upos += 2 * weight[gid];
		}else{
			uneg += 2 * weight[gid];
		}
	}
	node.max_gain = std::max(upos,uneg);
	if(fabs(opt_pat.gain) > node.max_gain ){ return true;}

	double gain_abs = fabs(gain);
	if(gain_abs > fabs(opt_pat.gain) || (fabs(gain_abs - fabs(opt_pat.gain))<1e-10 && pattern.size() < opt_pat.size)){
		opt_pat.gain = gain;
		opt_pat.optimalplace = &node;
		opt_pat.size = pattern.size();
	}
	return false;
}

void Gspan::CacheTree_search(){
	//std::cout << "cache tree search..." << std::endl;
	pattern.resize(0);
	for(list<Ctree*>::iterator it = croot->one_edge_graphs.begin();it != croot->one_edge_graphs.end();++it){
		pattern.push_back((*it)->pat.dcode);
		node_search(*(*it));
		pattern.pop_back();
	}
}

void Gspan::node_search(Ctree& node){
	//std::cout << "node search..." << std::endl;
	if(can_prune(node)) { return;}
	if(pattern.size() >= maxpat){ return;}
	if(node.children.size() == 0){
		can_grow_search(&node);
	}

	for(list<Ctree*>::iterator it = node.children.begin();it != node.children.end();++it){
		pattern.push_back((*it)->pat.dcode);
		node_search(**it);
		pattern.pop_back();
	}
}

void Gspan::can_grow_search(Ctree* n){
	pattern = n->pat.rebuild();
	//std::cout << pattern << std::endl;
	//edge_grow(**it);

	if(fabs(opt_pat.gain) > n->max_gain ){return;}

	PairSorter b_heap;
	map<int,PairSorter,greater<int> > f_heap;
	int maxtoc = scan_gspan(n->g2tracers,b_heap,f_heap);

	// projecting...
	DFSCode  dcode;
	for(PairSorter::iterator it = b_heap.begin(); it != b_heap.end(); ++it){	
		dcode.labels = Triplet(-1,it->first.b,-1);
		dcode.time.set(maxtoc, it->first.a);
		pattern.push_back(dcode);

		if(support(it->second) < minsup){pattern.pop_back(); return;}
		if(!is_min()){pattern.pop_back(); return;}

		//new Node projecting ...
		Ctree *node = new Ctree;
		TNnum++;
		nDFSCode d;
		d.set(dcode,&n->pat);
		node->pat = d;
		GraphToTracers& g2 = node->g2tracers;
		g2 = it->second;
		n->children.push_back(node);

		edge_grow(*node);
		pattern.pop_back();
	}

	for(map<int,PairSorter,greater<int> >::iterator it = f_heap.begin(); it != f_heap.end(); ++it){	
		for(PairSorter::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2){		
			dcode.labels = Triplet(-1,it2->first.a,it2->first.b);
			dcode.time.set(it->first,maxtoc+1);
			pattern.push_back(dcode);

			if(support(it2->second) < minsup){pattern.pop_back(); return;}
			if(!is_min()){pattern.pop_back(); return;}

			//new Node projecting ...
			Ctree *node = new Ctree;
			TNnum++;
			nDFSCode d;
			d.set(dcode,&n->pat);
			node->pat = d;

			GraphToTracers& g2 = node->g2tracers;
			g2 = it2->second;
			n->children.push_back(node);

			edge_grow(*node);
			pattern.pop_back();
		}
	}
}

void update(DPat& op){
	op.locsup.clear();
	for(GraphToTracers::iterator it=op.optimalplace->g2tracers.begin();it!=op.optimalplace->g2tracers.end();++it){
		op.locsup.push_back(it->first);
	}
	std::ostrstream ostrs;
	ostrs << op.optimalplace->pat.rebuild();
	ostrs << std::ends;
	op.dfscode = ostrs.str();

}
