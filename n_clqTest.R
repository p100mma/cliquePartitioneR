args=commandArgs(trailingOnly=TRUE)
ndiv= as.integer(args[[1]])
score=as.integer(args[[2]])



library(HCCsim)
library(cliquePartitioneR)
source('cliquePartitioneR/thr_seeke.R')
if (score==1) csf= complexity_clqScore else csf= n_clqScore

data(brca)

S<- similarity_matrix(fastPearsonData(brca))
S_ltr<- S[lower.tri(S)]

time_start= Sys.time()
 
grid_search_exhaustive( S=S, n_divisions=ndiv, 
	clqScoreFun=csf, 
	expansion_mode="basic")-> init_search 

precision_search_piecewiseLinear(max_neighborhood=init_search$max_neighborhood, S=S, clqScoreFun=csf, tol=0.1)-> precis_search
print("table of cliques from precis search")
print(table(precis_search$partition))
print("n cliques before precis")
print(length(unique(init_search$partitions[[ init_search$max_idx]] )))
print("n cliques after precis")
print(length(unique(precis_search$partition)))
print("best score before precis")
print(init_search$Y_t[[init_search$max_idx]])
print("score & thr: precis")
print(precis_search$opt_res) 

print(Sys.time()-time_start)
