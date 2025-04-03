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
IntegerVector rcpplique_gw_sumPC(NumericMatrix W_r){
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
	std::vector<float> strs (N_obj,0);        
	std::vector<int> degs (N_obj,0);        
	for (int i=0; i<N_obj; ++i){ v_i[i]=i;
	for (int j=0; j<N_obj; ++j)  {degs[i]+=A[i][j]; strs[i]+=W[i][j];}
		if (degs[i]==0) memsh[i]=0;
				    }
	int counter=0;
	int label=0;
	reduce_to_unexplored_subgraph2(A,W,v_i,strs,degs,memsh);
	while (!v_i.empty()){ //main loop
			counter+=1;
			//Rcout<<v_i.size()<<"\n";
			label+=1;
			int max_V= which_max(degs);
			memsh[v_i[max_V]]=label;
			//Rcout<<"building clique "<<label<<" \n";
			std::vector<int>  cl_neig= A[max_V];
			//score of the potential newcomer is 
			//sum of strengths of connections to the clique and its neighbors
			std::vector<float> neigh_score= W[max_V];
			for (int j=0; j<neigh_score.size();++j)
			       	if (cl_neig[j]==1)	
				for (int k=0; k<neigh_score.size();++k) 
				if (cl_neig[k]==1) neigh_score[j]+= W[k][j];
			//by construction of reduce_to_unexplored_subgraph 
			//there has to be at least one neighbor
			bool no_neigh=false; 
	//		Rcout<<"added vertex "<<max_V<<"\n";
			while (!no_neigh){
				//Rcout<<"adding newcomer \n";
				counter+=1;
				if (counter % 1000 == 0) {checkUserInterrupt(); Rcout<<v_i.size()<<"\n";};
				int newcomer= which_max(neigh_score);
	//			if (counter<10) Rcout<<"added vertex "<<newcomer<<"\n";
				memsh[v_i[ newcomer]]=label; 
				no_neigh=true;
				update_neigh_and_potential_score( cl_neig,
						A[newcomer], 
						neigh_score,	
						newcomer, 
						no_neigh, W);
	//			int current_size=0;
	//			for (int h=0; h<cl_neig.size();++h) if (cl_neig[h]) current_size++;
	//			Rcout<<" current clique neighborhood size is "<<current_size<<"\n";
		//		for (int i=0; i<neigh_score.size();++i) if( neigh_score[i]<0 ) Rcout<<neigh_score[i]<<" "<<cl_neig[i]<<"\n";
			}; 
		reduce_to_unexplored_subgraph2(A,W,v_i,strs,degs,memsh);
		}
	for (int i=0; i<memsh.size(); ++i) if (memsh[i]==-1) memsh[i]=0;
	IntegerVector membership(memsh.size());
	for (int i=0; i<memsh.size(); ++i) membership[i]=memsh[i];
	Rcout<<counter;
	return membership;
	}
