//#include <cmath>
//#include <vector>
//#include <algorithm>
#include <Rcpp.h>
#include "rcpplique_header.h"

using namespace Rcpp;

//' Rcpp Implementation of Basic Greedy Clique Partitioner
//'
//' This function assumes the matrix is symmetric and its diagonal is zeroed out
//' and performs no checks for correctness of the input.
//' Use \code{greedyCliquePartitioner()} for a robust version.
//'
//' @param W_r A numeric matrix, assumed to be symmetric with zeroes on its diagonal.
//' @return An integer vector encoding the clique membership of each node. Zeroes encode singletons.
//' @export 
// [[Rcpp::export]]
IntegerVector rcpplique_gw_basic(NumericMatrix W_r){
	int N_obj = W_r.ncol();
	std::vector<std::vector<float>> W (N_obj);
	for (auto &row : W) row.resize(N_obj);
	std::vector<std::vector<int>> A (N_obj);
	for (auto &row : A) row.resize(N_obj);

	for (int i=0; i< N_obj; ++i)
		for (int j=0; j< N_obj; ++j)
			{
			W[i][j]=W_r(i,j);
			A[i][j]= (W[i][j]>0);
			}	
	std::vector<int> memsh (N_obj,-1);	//-1 is label of` to be checked node
	std::vector<int> v_i (N_obj);         //vertex index (original)
	std::vector<float> degs (N_obj,0);         //vertex index (original)
	for (int i=0; i<N_obj; ++i){ v_i[i]=i;
	for (int j=0; j<N_obj; ++j)  degs[i]+=W[i][j];
		if (degs[i]==0) memsh[i]=0;
	}
	int counter=0;
	int label=0;
	reduce_to_unexplored_subgraph(A,W,v_i,degs,memsh);
	while (!v_i.empty()){ //main loop
			counter+=1;
			label+=1;               //TO DO: what if max is 0, handle this.
			std::vector<int> max_W_ij= max_edge( W, degs);
			for (auto i: max_W_ij) memsh[v_i[i]]=label;
			std::vector<float> hub_cons= W[max_W_ij[0]];
			std::vector<int>  cl_neig= A[max_W_ij[0]];
			update_neigh(cl_neig, A[max_W_ij[1]]); //diag of A is zero so here also max edge nodes are set to 0 such that they are not considered in next steps
			bool no_neigh=true;
			for (int i=0; i<hub_cons.size(); ++i) if (!cl_neig[i]) hub_cons[i]=0; else {if (no_neigh) no_neigh=false;};
			while (!no_neigh){   
				counter+=1;
				if (counter % 1000 == 0) {checkUserInterrupt(); Rcout<<v_i.size()<<"\n";};
				int newcomer= which_max(hub_cons);
				memsh[v_i[ newcomer]]=label; 
				update_neigh(cl_neig, A[newcomer]); 
				no_neigh=true;
				for (int i=0; i<hub_cons.size(); ++i) if (!cl_neig[i]) hub_cons[i]=0; else {if (no_neigh) no_neigh=false;};
			}; 
		reduce_to_unexplored_subgraph(A,W,v_i,degs,memsh);
		}
	for (int i=0; i<memsh.size(); ++i) if (memsh[i]==-1) memsh[i]=0;
	IntegerVector membership(memsh.size());
	for (int i=0; i<memsh.size(); ++i) membership[i]=memsh[i];
	return membership;
	}
