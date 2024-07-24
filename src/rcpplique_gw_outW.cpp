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
IntegerVector rcpplique_gw_outW(NumericMatrix W_r){
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
	std::vector<float> degs (N_obj,0);       
	for (int i=0; i<N_obj; ++i){ v_i[i]=i;
	for (int j=0; j<N_obj; ++j)   degs[i]+=W[i][j];
		if (degs[i]==0) memsh[i]=0;
				    }
	int counter=0;
	int label=0;
	reduce_to_unexplored_subgraph(A,W,v_i,degs,memsh);
	while (!v_i.empty()){ //main loop
			counter+=1;
			label+=1;
			std::vector<int> max_W_ij= max_edge( W, degs);
			for (auto i: max_W_ij) memsh[v_i[i]]=label;
			std::vector<int>  cl_neig= A[max_W_ij[0]];
			std::vector<float> out_deg(cl_neig.size()); 
			for (int i=0; i< v_i.size(); ++i){
				if ( i==max_W_ij[0] ) {out_deg[i] = -1;} else {
						      out_deg[i] = degs[i];
						      if (A[i][max_W_ij[0]]) out_deg[i]= out_deg[i] - W[i][max_W_ij[0]]; 
					}
				}
			int current_size=1;
			bool no_neigh=true;
			update_neigh_and_outStr( cl_neig,
				A[max_W_ij[1]], 
				out_deg,
				W[max_W_ij[1]],	
				no_neigh);
			current_size++;
			while (!no_neigh){
				counter+=1;
				if (counter % 1000 == 0) {checkUserInterrupt(); Rcout<<v_i.size()<<"\n";};
				int newcomer= which_max(out_deg);
				memsh[v_i[ newcomer]]=label; 
				no_neigh=true;
				update_neigh_and_outStr( cl_neig,
						A[newcomer], 
						out_deg,	
						W[newcomer],
						no_neigh);
			current_size++;
			}; 
		reduce_to_unexplored_subgraph(A,W,v_i,degs,memsh);
		}
	for (int i=0; i<memsh.size(); ++i) if (memsh[i]==-1) memsh[i]=0;
	IntegerVector membership(memsh.size());
	for (int i=0; i<memsh.size(); ++i) membership[i]=memsh[i];
	return membership;
	}
