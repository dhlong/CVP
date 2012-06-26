#include "network.h"

//////////////////////////////////////////////////////////////
////////// Arc and Commodity Flow constructors
//////////
//////////////////////////////////////////////////////////////

NetworkArc::NetworkArc(Vertex o, Vertex d, Real C, Real c):
	head(o), tail(d), cap(C), cost(c) 
{}
  
CommoFlow::CommoFlow(Vertex o, Vertex d, Real D): 
	origin(o), destination(d), demand(D) 
{}


//////////////////////////////////////////////////////////////
////////// Graph constructors and functions
//////////
//////////////////////////////////////////////////////////////

int Graph::getNVertex() const{
	int V_=0;
	FOR(i, arcs.size()){
		if(V_ < arcs[i].head) V_ = arcs[i].head;
		if(V_ < arcs[i].tail) V_ = arcs[i].tail;
	}
	return V_+1;
}

void Graph::write_pajek(const char* filename){
	fstream f(filename, fstream::out);
	int V = getNVertex();
	f<<"*Vertices "<<V<<endl;
	FOR(i,V) f<<i+1<<" \""<<i+1<<"\""<<endl;
	f<<"*Arcs*"<<endl;
	FOR(i,arcs.size()) f<<arcs[i].head+1<<" "<<arcs[i].tail+1<<endl;
	f.close();
}

//////////////////////////////////////////////////////////////
////////// Single commodity Network constructors
//////////
//////////////////////////////////////////////////////////////


Network::Network(): Graph(), supply(), demand() {
}

// constructor by reading file
Network::Network(const char* filename, FileFormat format)
	: Graph(), supply(), demand() {

	if(format== NETGEN){
		fstream f(filename, fstream::in);
		char line[300];
		int status=0, V, E;
		Vertex i, j;
		Real r, c;

		if(!f) error_handle("File does not exist");
    
		while(!f.eof()){
			f.getline(line,300);
			if (strstr(line, "PROBLEM NUMBER") != NULL){
				sscanf(line, 
				       "%*[^0123456789]%d"
				       "%*[^0123456789]%d"
				       "%*[^0123456789]%d",
				       &i, &V, &E);
				break;
			}
		}

		while(!f.eof()){
			f.getline(line,300);
			if(strstr(line, "SUPPLY") != NULL) status = 1;
			else if(strstr(line, "ARCS") != NULL) status = 2;
			else if(strstr(line, "DEMAND") != NULL) status = 3;
			else if(strstr(line, "END") != NULL) break;
			else{
				if(status == 1){
					sscanf(line,"%d %lf",&i,&r);
					supply.push_back(make_pair(i-1,r));
				} 
				else if(status == 3){
					sscanf(line,"%d %lf",&i,&r);
					demand.push_back(make_pair(i-1,r));
				}
				else if(status == 2){
					sscanf(line,"%d%d%lf%lf",&i,&j,&r,&c);
					arcs.push_back(NetworkArc(i-1,j-1,c));
				}
			}
		}

		f.close();
	}
	else{
		error_handle("File format not supported for single commodity nework");
	}   
}

// constructor by grid generation
// generate a nxn grid network with 1 supply vertex at (0,0) 
// and 1 demand vertex at (n-1, n-1)
// and supply = demand = cap
Network::Network(int n, Real cap){
	int dx[2] = {0, 1}, dy[2] = {1, 0};
	int a, b, ii, jj;
	Real c;

	FOR(i,n) FOR(j,n) FOR(d,2){
		ii = i+dx[d], jj = j+dy[d];
		if(ii<0 || jj<0 || ii>=n || jj>=n) continue;

		a = i*n + j; b = ii*n + jj;
		c = (a==0 || b == n*n-1)? cap/2 : cap/4;
		arcs.push_back(NetworkArc(a,b,c));
	}

	supply.push_back(make_pair(0,cap));
	demand.push_back(make_pair(n*n-1,cap));
}

//////////////////////////////////////////////////////////////
////////// Multi-commodity Network constructors
//////////
//////////////////////////////////////////////////////////////

// constructor by reading file
MultiCommoNetwork::MultiCommoNetwork(const char * filename, FileFormat format)
	: Graph(), commoflows()
{
	cout<<"Reading "<<filename<<endl;
	if(format == GENFLOT){
		cout<<"GENFLOT format"<<endl;
		fstream f(filename, fstream::in);
		int V, A, K, o, d, r, c;
		f>>V>>A>>K;
		FOR(j,A) f>>o>>d>>c>>r, arcs.push_back(NetworkArc(o-1,d-1,r,c));
		FOR(i,K) f>>o>>d>>r, commoflows.push_back(CommoFlow(o-1,d-1,r));
		f.close();
	}

	else if(format == TNTP){
		cout<<"TNTP format"<<endl;
		string tripsfilename(filename);
		size_t pos = tripsfilename.find("_net");
		tripsfilename.replace(pos, 4, "_trips");
		cout<<tripsfilename<<endl;
    
		fstream fnet(filename, fstream::in);
		fstream fcommo(tripsfilename.c_str(), fstream::in);
    
		if(!fnet || !fcommo) error_handle("File does not exist.");
    
		char line[300], *str;
		int A=-1, V=-1;
    
		while(!fnet.eof()) {
			fnet.getline(line,300);
    
			str = strstr(line, "<NUMBER OF LINKS>");
			if(str!=NULL) {
				stringstream ss(str+17, stringstream::in);
				ss>>A;
			}

			str = strstr(line, "<NUMBER OF NODES>");
			if(str!=NULL){
				stringstream ss(str+17, stringstream::in);
				ss>>V;
			}

			if(strstr(line, "~")!=NULL) break;
		}
    
		if(A==-1 || V==-1) 
			error_handle("Input error: "
			             "no specification of number of links or nodes");

		int head, tail;
		Real cost, cap, tmp, demand;

		FOR(i, A){
			fnet.getline(line,300);
			stringstream ss(line, stringstream::in);
			ss>>head>>tail>>cap>>tmp>>cost;
			arcs.push_back(NetworkArc(head-1, tail-1, cap, cost));
		}

		head = -1;
		while(!fcommo.eof()){
			fcommo.getline(line,300);
			if(strstr(line,"<")!=NULL) continue;

			if((str = strstr(line,"Origin")) != NULL){
				stringstream ss(line+7,stringstream::in);
				ss>>head;
				continue;
			}

			if(head == -1) continue;
      
			stringstream ss(line, stringstream::in);
			while(!ss.eof()){
				ss>>tail; ss.ignore(300,':');
				if(ss.eof()) break;
				ss>>demand; ss.ignore(300,';');
				commoflows.push_back(CommoFlow(head-1,tail-1,demand));
			}
		}
  
		fnet.close(); fcommo.close();
	}

	else{
		error_handle("File format not supported for multi-commodinty network");
	}
}


ShortestPathOracle::ShortestPathOracle(const MultiCommoNetwork &n):
	net(n), V(n.getNVertex()), A(n.arcs.size()), K(n.commoflows.size()),
	indexarcl(V), indexadjl(V), trace(V), vb(V), nv(V, 0)
{
	adjl.V = V; adjl.A = A;
	malloc_adjl(&adjl);

	FOR(i, V) trace[i] = MALLOC(vertex_t, V);
	FOR(i, V) vb[i] = MALLOC(char, V);
  
	FOR(i,V)  
		indexarcl[i] = vector<arc_t>(V, -1), 
		indexadjl[i] = vector<arc_t>(V, -1);
  
	FOR(i, V+1) adjl.n_arcs[i] = 0;
	FOR(a, A)   adjl.n_arcs[net.arcs[a].head]++;
	FOR(i, V)   adjl.n_arcs[i+1] += adjl.n_arcs[i];
    
	FOR(a, A) {
		int aa =  --adjl.n_arcs[net.arcs[a].head];
		int u = net.arcs[a].head, v = net.arcs[a].tail;
		indexarcl [u] [v] = a;
		indexadjl [u] [v] = aa;
		adjl.adjacent_vertices[aa] = net.arcs[a].tail;
	}
  
	FOR(i, V) FOR(j, V) vb[i][j]= 0;
	FOR(k, net.commoflows.size()) {
		vb[net.commoflows[k].origin][net.commoflows[k].destination] = 1;
		nv[net.commoflows[k].origin]++;
	}  

	heap = MALLOC(vertex_t, V);
	pos  = MALLOC(index_t,  V);
	d = MALLOC(cost_t, V);
  

	has_solved = false;
}

ShortestPathOracle::~ShortestPathOracle(){
	FOR(i, V){
		FREE(vb[i]);
		FREE(trace[i]);
	}
	FREE(heap);
	FREE(pos);
	FREE(d);
	free_adjl(&adjl);
}

extern fstream iteration_report;

void ShortestPathOracle::solve(){
	FOR(i, V) if(nv[i]>0) 
		dijkstra(adjl, i, vb[i], nv[i], heap, pos, d, trace[i]);
	has_solved = true;
}

void ShortestPathOracle::get_flows(Vector &sp, bool use_tmp) {
	sp *= 0.0;
	vector< list<arc_t> > *paths = use_tmp ? new vector< list<arc_t> >(A) : NULL;
	if(!has_solved) solve();

	FOR(k, K){
		int u = net.commoflows[k].origin, v = net.commoflows[k].destination;
		Real demand = net.commoflows[k].demand;
		if(use_tmp)
			while(v>=0 && v!=u){
				(*paths)[indexarcl[trace[u][v]] [v]].push_back(k);
				v = trace[u][v];
			}
		else
			while(v>=0 && v!=u){
				sp.insert(indexarcl[trace[u][v]][v]*K + k) = demand;
				v = trace[u][v];
			}
	}
  
	if(use_tmp){
		FOR(a, A)
			for(list<arc_t>::iterator it = (*paths)[a].begin(); it != (*paths)[a].end(); ++it)
				sp.insert(a*K + *it) = net.commoflows[*it].demand;
		delete paths;
	}
}
