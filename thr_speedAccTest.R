args=commandArgs(trailingOnly=TRUE)
ndiv= as.integer(args[[1]])  #n steps in exhaustive thr search
score=as.integer(args[[2]])
prec= as.integer(args[[3]])  #if 1 then do precision search after exhaustive
which.dataset= as.integer(args[[4]])  # 1- BRCA; 2 - KIRC; 3- AGP



library(HCCsim)
library(cliquePartitioneR)
library(MDFS)
source('thr_seeke.R')
if (score==1) csf= complexity_clqScore else csf= n_clqScore

ds_names=c('gene_expr_data.rds', 'KIRC_gene_expr_data.rds', 'agp.rds')

X= readRDS(ds_names[[which.dataset]])

if (which.dataset==1) {
#	data(brca)
#	S<- similarity_matrix(fastPearsonData(brca))
	S<- similarity_matrix(fastPearsonData(X))
} else if (which.dataset==2) {
X<-readRDS("KIRC_gene_expr_data.rds")
V_X=apply(X,2,var)
X=as.matrix(X)
hV= V_X >= quantile(V_X, 0.25)
X=X[,hV]
gc()
S<- similarity_matrix(fastPearsonData(X))
} else if (which.dataset==3) {
X<-readRDS("agp.rds")
nrow1<-nrow(X)
td_rowsums<- rowSums(X)
max_rowsum<-3636 #obtained by Piotr from theory
small_rows<- td_rowsums < max_rowsum 
X<- X[!small_rows, ]
tdrownames<- rownames(X)
X<-X/rowSums(X)
rownames(X)<- tdrownames
#binarize
min_nonzero_counts<-20 #minimum nonzero counts to do any statistics with a taxon
Xbinary<-0*X
for(i in 1:ncol(X)) Xbinary[,i]<-(X[,i]>median(X[,i]))
Xbinary<-Xbinary[,colSums(Xbinary)>=min_nonzero_counts]
#remove Archaea and Eukaryotes
archaea_eucaryotes<-c(1:9, 1022:1028)
Xbinary[,-archaea_eucaryotes]->Xbinary
pairwiseMIsimilarityDiscrete(Xbinary)$S-> S
}


S_ltr<- S[lower.tri(S)]

exh_time_start= Sys.time()
 
grid_search_exhaustive( S=S, n_divisions=ndiv, 
	clqScoreFun=csf, 
	expansion_mode="basic")-> init_search 

exh_time= Sys.time() - exh_time_start

if (prec==0)
saveRDS( 
	list( 
	      time= exh_time,
	      partition=init_search$partitions[[ init_search$max_idx ]],
	      maximum= init_search$X_t[[ init_search$max_idx]],
	      objective=init_search$Y_t[[init_search$max_idx]],
	      X_t= init_search$X_t,
	      Y_t= init_search$Y_t
	    ),
	 sprintf("data=%d_ndiv=%d_score=%d_prec=%d.rds",
		  which.dataset, ndiv, score, prec) )

if (prec==1) {

prec_time_start=Sys.time()

precision_search_piecewiseLinear(max_neighborhood=init_search$max_neighborhood, S=S, clqScoreFun=csf, tol=0.1)-> precis_search
prec_time= Sys.time() - prec_time_start

saveRDS(
	list(
	      exh_time=exh_time,
	      exh_partition=init_search$partitions[[ init_search$max_idx ]],
	      exh_maximum= init_search$X_t[[ init_search$max_idx]],
	      exh_objective=init_search$Y_t[[init_search$max_idx]],
	      tot_time= exh_time + prec_time,
	      prec_partition= precis_search$partition,
	      prec_maximum=precis_search$opt_res$maximum,
	      prec_objective=precis_search$opt_res$objective),
	 sprintf("data=%d_ndiv=%d_score=%d_prec=%d.rds",
		  which.dataset, ndiv, score, prec) 
	    )
	      }
