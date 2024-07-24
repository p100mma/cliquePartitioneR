#include <Rcpp.h>
using namespace Rcpp;

int n_left(const std::vector<int> &mem){
	int nl=0;
	for (auto x: mem) if (x==-1) nl+=1;
	return nl;
}

bool is_neg1(int x) { return x==-1; };
bool a_is_empty(const std::vector<int> &x) { return x.empty(); };	
bool w_is_empty(const std::vector<float> &x) { return x.empty(); };

void reduce_to_unexplored_subgraph(
	std::vector<std::vector<int>> &A,
	std::vector<std::vector<float>> &W,
	std::vector<int> &V,
	std::vector<float> &DEG,
        const std::vector<int> &mem){
	std::vector<int> indexes;	//inside V, DEG, W, A, potentially to stay
	std::vector<int> indexes2go;	
	for (int j=0; j< V.size(); j++) if (mem[V[j]]==-1) indexes.push_back(j); 
					else {
					V[j]=-1;
					DEG[j]=-1;
					W[j].clear();
					A[j].clear();
					for (int k=0; k< W.size(); k++)
						{
						if (!W[k].empty()){
						W[k][j]=-1;
						A[k][j]=-1;
							  	  }
						}
					      }
	for (int i: indexes){ bool i_stays=false;
		for(int j: indexes) if (A[i][j]>0) {i_stays=true; break;};
		if (!i_stays) { indexes2go.push_back(i);
				V[i]=-1;
				DEG[i]=-1;
				W[i].clear();
				A[i].clear();
				for (int k=0; k< W.size(); k++)
					{
					if (!W[k].empty()){
					W[k][i]=-1;
					A[k][i]=-1;
							  }
					}
				}
					}
	V.erase(std::remove_if(V.begin(),V.end(), is_neg1), V.end());
	DEG.erase(std::remove_if(DEG.begin(),DEG.end(), is_neg1), DEG.end());
	W.erase(std::remove_if(W.begin(),W.end(), w_is_empty), W.end());
	A.erase(std::remove_if(A.begin(),A.end(), a_is_empty), A.end());
	for (int i=0; i<W.size();++i)
		if (!A[i].empty())
		{
		W[i].erase(std::remove_if(W[i].begin(),W[i].end(), is_neg1), W[i].end());
		A[i].erase(std::remove_if(A[i].begin(),A[i].end(), is_neg1), A[i].end());
		}
}


int which_max (const std::vector<float> &W_i){
	return	std::distance(W_i.begin(),std::max_element(W_i.begin(),W_i.end()));
}

int which_max (const std::vector<int> &A_i){
	return	std::distance(A_i.begin(),std::max_element(A_i.begin(),A_i.end()));
}

std::vector <int> max_edge (const std::vector<std::vector<float>> &W,
			    const std::vector<float> &degrees
			    ){
	std::vector <float> colMaxes (W.size());
	for (int j=0; j< W.size(); ++j)
		colMaxes[j]=*(std::max_element(W[j].begin(), W[j].end()));
	int max_j =which_max(colMaxes);
	int max_i =which_max(W[max_j]);
	std::vector <int> edge (2);
	if (degrees[max_j] > degrees[max_i] ) edge= { max_j, max_i }; 
	else edge= {max_i, max_j};
	return edge;
}

void update_neigh( std::vector<int> &current_neigh,
				 std::vector<int> &new_neigh){
for (int i=0; i< current_neigh.size(); i++) 
	current_neigh[i] = current_neigh[i] && new_neigh[i];
}

void update_neigh_and_avg( std::vector<int> &current_neigh,
				 std::vector<int> &new_neigh,
				std::vector<float> &avg_con,
				 std::vector<float> &new_con,
				bool &no_cl_neigh){
update_neigh(current_neigh,new_neigh);
for (int i=0; i< avg_con.size(); i++) 
	if (current_neigh[i]){
	no_cl_neigh=false;
	avg_con[i] = avg_con[i] + new_con[i];
	} else
	avg_con[i] = 0;
}

void update_neigh_and_max( std::vector<int> &current_neigh,
				 std::vector<int> &new_neigh,
				std::vector<float> &max_con,
				 std::vector<float> &new_con,
				bool &no_cl_neigh){
update_neigh(current_neigh,new_neigh);
for (int i=0; i< max_con.size(); i++) 
	if (current_neigh[i]){
	no_cl_neigh=false;
	max_con[i] = std::max( max_con[i], new_con[i]);
	} else
	max_con[i] = 0;
}

void update_neigh_and_outDeg( std::vector<int> &current_neigh,
				 std::vector<int> &new_neigh,
				std::vector<int> &out_deg,
				bool &no_cl_neigh){
update_neigh(current_neigh,new_neigh);
for (int i=0; i< out_deg.size(); i++) 
	if (current_neigh[i]){
	no_cl_neigh=false;
	out_deg[i]= new_neigh[i] ? out_deg[i]-1 : out_deg[i];	
	} else
	out_deg[i]= -1;	
}

void update_neigh_and_outStr( std::vector<int> &current_neigh,
				 std::vector<int> &new_neigh,
				std::vector<float> &out_deg,
				std::vector<float> newcomer_con,
				bool &no_cl_neigh){
update_neigh(current_neigh,new_neigh);
for (int i=0; i< out_deg.size(); i++) 
	if (current_neigh[i]){
	no_cl_neigh=false;
	out_deg[i]= new_neigh[i] ? out_deg[i]- newcomer_con[i] : out_deg[i];	
	} else
	out_deg[i]= -1;	
}
