#ifndef GSPAN_H_
#define GSPAN_H_

#include <map>
#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <ext/hash_map>
#include <math.h>

using std::map;
using std::list;
using std::string;
using std::vector;
using std::greater;
using __gnu_cxx::hash_map;

struct Triplet {
public:
  int x, y, z;
  Triplet reverse();
  explicit Triplet(){}
  explicit Triplet(int _x,int _y,int _z): x(_x),y(_y),z(_z){}
};

struct Pair { public: int a, b; void set(int, int); };
inline void Pair::set(int _a, int _b){ a = _a; b = _b; }

struct VertexPair: public Pair { public: int id; };
struct Edge   { public: int to, id; Triplet labels; };
struct DFSCode { public: Triplet labels; Pair time; };

typedef vector< vector<Edge> > AdjacentList;

class Graph: public AdjacentList {
 public:
  vector<int> label;
  int num_of_edges;
  int class_label;
  void resize(unsigned int s){ AdjacentList::resize(s); label.resize(s); }
};

class EdgeTracer {
 public:
  VertexPair  vpair;
  EdgeTracer* predec;
  explicit EdgeTracer(){};
  void set(int,int,int,EdgeTracer*);
};

inline void EdgeTracer::set(int a,int b,int id, EdgeTracer* pr){
  vpair.a = a;
  vpair.b = b;
  vpair.id = id;
  predec = pr;
}

typedef list<EdgeTracer> Tracers;
typedef map<int,Tracers> GraphToTracers;
inline void debugGraphToTracers(GraphToTracers);
typedef map<Pair,GraphToTracers> PairSorter;

inline Triplet Triplet::reverse(){ return Triplet(z,y,x); }

inline bool operator< (const Triplet& left, const Triplet& right){
  if (left.x!=-1 && right.x!=-1 && left.x != right.x) return (left.x < right.x);
  if (left.y!=-1 && right.y!=-1 && left.y != right.y) return (left.y < right.y);
  return (left.z < right.z);
}

inline bool operator<= (const Triplet& left, const Triplet& right){
  return !(right < left);
}

inline bool operator== (const Triplet& left, const Triplet& right){
  return (left.x==right.x && left.y==right.y && left.z==right.z);
}

inline bool operator< (const Pair& left, const Pair& right){
  if (left.a != right.a) return (left.a < right.a);
  return (left.b < right.b);
}

inline bool operator<= (const Pair& left, const Pair& right){
  return !(right < left);
}

inline bool operator== (const Pair& left, const Pair& right){
  return (left.a==right.a && left.b==right.b);
}

struct Pair_IT {
	public:
		int i;
		Pair p;
};

inline bool operator< (const Pair_IT& left, const Pair_IT& right){
  if (left.i != right.i) return (left.i < right.i);
  return (left.p < right.p);
}

inline bool operator<= (const Pair_IT& left, const Pair_IT& right){
  return !(right < left);
}

inline bool operator== (const Pair_IT& left, const Pair_IT& right){
  return (left.i==right.i && left.p==right.p);
}

inline bool operator== (const DFSCode& left, const DFSCode& right){
  if(left.time.a != right.time.a) return false;
  if(left.time.b != right.time.b) return false;	
  if(left.labels.x != right.labels.x) return false;
  if(left.labels.y != right.labels.y) return false;
  return (left.labels.z == right.labels.z);	
}

inline bool operator!= (const DFSCode& x, const DFSCode& y){
  return !(x==y);
}
inline bool operator< (const DFSCode& left, const DFSCode& right){
  if(left.time.a != right.time.a) return left.time.a > right.time.a;
  
  if(left.time.b != right.time.b) return left.time.b < right.time.b;
  
  if(left.labels.x != right.labels.x) return left.labels < right.labels;
  if(left.labels.y != right.labels.y) return left.labels.y < right.labels.y;
  return (left.labels.z < right.labels.z);	
}
inline std::ostream& operator<< (std::ostream& os, const vector<DFSCode> pattern){
  if(pattern.empty()) return os;
  os << "(" << pattern[0].labels.x << ") " << pattern[0].labels.y << " (0f" << pattern[0].labels.z << ")";
  /*
  for(unsigned int i=1; i<pattern.size(); ++i){
    if(pattern[i].time.a < pattern[i].time.b){
      os << " " << pattern[i].labels.y << " (" << pattern[i].time.a << "f" << pattern[i].labels.z << ")";
    }else{
      os << " " << pattern[i].labels.y << " (b" << pattern[i].time.b << ")";
    }
  }
  */
  for(unsigned int i=1; i<pattern.size(); ++i){
	  if (pattern[i].labels.x == -1 and pattern[i].labels.z == -1) {
		  os << " " << pattern[i].labels.y << " (" << pattern[i].time.a << "b" << pattern[i].time.b << ")";
	  } else {
		  os << " " << pattern[i].labels.y << " (" << pattern[i].time.a << "f" << pattern[i].labels.z << ")";
	  }
  }
  return os;
}

class nDFSCode {
 public:
  DFSCode  dcode;
  nDFSCode* pred;
  explicit nDFSCode(){};
  void set(DFSCode,nDFSCode*);
  vector<DFSCode> rebuild();
};

inline void nDFSCode::set(DFSCode _dcode, nDFSCode* pr){
  dcode = _dcode;
  pred = pr;
}

inline vector<DFSCode> nDFSCode::rebuild(){
  vector<DFSCode> p;
  if(pred != NULL){
    p  = pred->rebuild();
  }
  p.push_back(dcode);
  return p;
}

const vector<Graph> readGraphs(std::istream&);

struct DPat{//discrimination pattern
  std::string dfscode;
  GraphToTracers optimalplace;
  unsigned int size;
  vector<DFSCode> dfs;
  vector<int> locsup;
  double gain;
};


class Gspan {
 private:
  bool is_min();
  unsigned int support(GraphToTracers&);
  void scan_rm(vector<DFSCode>&, vector<int>&);
  bool min_checker(vector<DFSCode>&, Graph&, Tracers&);
  int  scan_gspan(GraphToTracers&, map<int,PairSorter,greater<int>>&, map<int,PairSorter,greater<int>>&);
 public:
  bool out_instances;
  unsigned int minsup;
  unsigned int maxpat;
  unsigned int flownum;
  unsigned int wildcard_r;
  double lambda;
  vector<Graph>   gdata;
  vector<DFSCode> pattern;
  void set_data(std::istream& is) {
    gdata = readGraphs(is);
  };
  void report(GraphToTracers&);

  //lpboost
  int vnum;
  int search_vnum;
  int scanum;
  int search_scanum;
  vector<double> weight;
  vector<double> corlab; 
  unsigned int max_itr;
  DPat opt_pat;
  vector<DPat> opt_list;
  double wbias;
  double bound;
  double gain_;
  double nu;
  double conv_epsilon;
  bool can_prune(GraphToTracers&);
  void lpboost();
  void edge_grow(GraphToTracers&, int);
  void make_count_r(unsigned int, vector<int>&);
  void make_count_s(unsigned int, vector<int>&, map<Triplet, GraphToTracers>&);
  void make_count_s(unsigned int, vector<int>&, map<int, PairSorter, greater<int>>&, map<int, PairSorter, greater<int>>&);
  bool make_count_b(unsigned int, vector<int>&, map<Triplet, GraphToTracers>&);
  bool make_count_b(unsigned int, vector<int>&, map<int, PairSorter, greater<int>>&, map<int, PairSorter, greater<int>>&, int);
  double valuation();

  //caching
  bool first_flag;
  unsigned int TNnum;
  void sira_search();
  vector<int> update(DPat&);
};

Graph toGraph(vector<DFSCode>&);

inline void init(EdgeTracer& et,int a,int b,int id){
  et.vpair.a = a;
  et.vpair.b = b;
  et.vpair.id = id;
  et.predec = 0;
}

#endif /*GSPAN_H_*/
