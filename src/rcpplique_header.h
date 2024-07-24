#ifndef RCPPLIQUE_H
#define RCPPLIQUE_H

#include <Rcpp.h>

int n_left(const std::vector<int> &mem);

bool is_neg1(int x);
bool a_is_empty(const std::vector<int> &x);	
bool w_is_empty(const std::vector<float> &x);

void reduce_to_unexplored_subgraph(
	std::vector<std::vector<int>> &A,
	std::vector<std::vector<float>> &W,
	std::vector<int> &V,
	std::vector<float> &DEG,
        const std::vector<int> &mem);

int which_max (const std::vector<float> &W_i);
int which_max (const std::vector<int> &A_i);

std::vector <int> max_edge (const std::vector<std::vector<float>> &W,
			    const std::vector<float> &degrees
			    );

void update_neigh( std::vector<int> &current_neigh,
				 std::vector<int> &new_neigh);

void update_neigh_and_avg( std::vector<int> &current_neigh,
				 std::vector<int> &new_neigh,
				std::vector<float> &avg_con,
				 std::vector<float> &new_con,
				bool &no_cl_neigh);
void update_neigh_and_max( std::vector<int> &current_neigh,
				 std::vector<int> &new_neigh,
				std::vector<float> &max_con,
				 std::vector<float> &new_con,
				bool &no_cl_neigh);
void update_neigh_and_outDeg( std::vector<int> &current_neigh,
				 std::vector<int> &new_neigh,
				std::vector<int> &out_deg,
				bool &no_cl_neigh);
#endif
