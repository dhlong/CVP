#ifndef __NETWORK_H__
#define __NETWORK_H__

#include "cvputility.h"
#include "dijkstra.h"

using namespace std;

// File format for reading
enum FileFormat { GENFLOT, TNTP, NETGEN };

struct NetworkArc {
  Vertex head, tail; 
  Real cap, cost;
  NetworkArc(Vertex o, Vertex d, Real C, Real c = 0);
};

// Graph contains a list of arcs
struct Graph {
  vector<NetworkArc> arcs;  
  Graph();
  Graph(const Graph &);
  int getNVertex() const;  
  void write_pajek(const char* filename);
};

typedef pair<Vertex, Real> Sink;
typedef pair<Vertex, Real> Source;

// Single commodity network with multiple sources and sinks
struct Network : public Graph {
  vector<Source> supply;
  vector<Sink> demand;

  Network();
  Network(const Network &n);
  Network(int n, Real cap); // generate nxn grid network

  // construct by reading file
  Network(const char* filename, FileFormat format);
};


// Commodity flow
struct CommoFlow{
  int origin, destination;
  Real demand;
  CommoFlow(Vertex o, Vertex d, Real D);
};


// Multi-commodity network
struct MultiCommoNetwork : public Graph {
  vector <CommoFlow> commoflows;
  MultiCommoNetwork();
  MultiCommoNetwork(const MultiCommoNetwork &net);

  // construct by reading file
  MultiCommoNetwork(const char* filename, FileFormat format);
};

class ShortestPathOracle{
 private:
  MultiCommoNetwork net;
  int V, A, K;

  AdjacentList adjl;
  vector< vector<arc_t> > indexarcl;
  vector< vector<arc_t> > indexadjl;
  vector< vertex_t * > trace;
  vector< char * > vb;
  vector< int > nv;

  vertex_t *heap;
  index_t  *pos;
  cost_t *d;

  bool has_solved;

  void solve();

 public:
  ShortestPathOracle(const MultiCommoNetwork &n);
  ~ShortestPathOracle();

  void set_cost(vertex_t u, vertex_t v, cost_t c){
    arc_t a = indexadjl[u][v];
    if(a>=0) adjl.costs[a] = c;
    has_solved = false;
  }

  void get_flows(Vector &sp);
};

#endif

