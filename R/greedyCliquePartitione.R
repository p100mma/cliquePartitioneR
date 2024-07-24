#' Greedily partition input weighted graph into heavily connected cliques
#'
#' @details
#' Algorithm partitions nodes of the network by the following procedure. First, it picks the heaviest edge
#' available as a starting point for the 1st clique. Then the clique is expanded one node at a time,
#' by joining a non-clique breaking node according to the \code{expansion_mode} argument.
#' Clique is expanded until exhaustion, and a new clique is initialized on the graph excluding the cliques already
#' built.
#' 
#' If \code{expansion_mode} is set to \code{basic}, then the newcomers are neighbors of a seed node of the clique,
#' that is one of the nodes belonging to initializer edge (the one with the bigger weighted degree). Those
#' neighbors are chosen in the order of their connection strength with the seed node.
#'
#' If \code{expansion_mode} is set to \code{average} or \code{max}, then the newcomers are chosen by maximizing 
#' aggregated strength of the connections of newcomer to the members of the clique. Aggregation function used
#' for \code{average} is mean, while for \code{max}, the aggregated connection strength of newcomer is set to
#' the max over all of its connections to the clique.
#' 
#' Lastly, if \code{expansion_mode} is set to \code{outN} or \code{outW}, then newcomers are chosen by maximizing
#' their outer connectivity. For \code{outN}, procedure picks newcomer with maximal number of connections to nodes
#' outside of the clique. For \code{outW}, sum of weights of the connections of the newly joined node to the nodes
#' outside of the clique is maximized.
#'
#' @param W A symmetric matrix encoding weights between nodes, i.e. \code{W[i,j]} is the weight of the edge between nodes \code{i,j}. Weights are interpreted as similarities. \code{W[i,j]==0} implies no connection between \code{i,j}.
#' @param check_symmetry Logical telling whether to check symmetry of the input matrix. Not done by default.
#' @param expansion_mode Charecter string encoding one of possible ways the algorithm shall build a clique based on the strength of the connections of the neighbors.
#' @param unique_singletonLabels Logical telling whether isolated nodes should carry unique labels in the output partition. By default, each singleton gets a label of \code{0}.
#' @return A \code{list} containing following two elements:
#' \itemize{
#' \item \code{membership} -- a integer vector of length \code{ncol(W)} giving labels of cliques in which each node belongs, i.e. \code{membership[[i]]} is label of the clique of node \code{i}
#' \item \code{alg_params} --- a list containing parameters of used greedy clique partitioner algorithm
#' }
#' @examples
#' input_graph= matrix(runif(40*40)*sample(c(0,1),40*40,replace=TRUE, prob=c(0.7,0.3)),ncol=40)
#' input_graph= input_graph %*% t(input_graph) #generate symmetric matrix
#' greedyCliquePartitioner( input_graph)
#' @useDynLib cliquePartitioneR, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @export
greedyCliquePartitioner<- function(W,
				   check_symmetry=FALSE,
				   expansion_mode="basic",
				   unique_singletonLabels=FALSE) {


if (!(expansion_mode %in% c('basic', 'average', 'max', 'outN', 'outW')) )
	stop("unknown value of expansion_mode")

#check correctness of W matrix
if (length(W)==1) return(1)
if (is.null(dim(W))) stop("W should have a dim attribute")
stopifnot(length(dim(W))==2)
stopifnot(!is.null(ncol(W)))
stopifnot(!is.null(nrow(W)))
stopifnot(ncol(W)==nrow(W))
N<- ncol(W)
#check for trivial graphs
if (N==0){ message("W describes empty graph, returning NULL)")
		 return(NULL)
		}
if (N==1) {
	   return(1)
	  }

#symmetry check
if (check_symmetry) stopifnot(isSymmetric(W))

#0 out diagonal by convention
diag(W)<-0

if (expansion_mode=='basic')
	partition<-rcpplique_gw_basic(W)
if (expansion_mode=='average')
	partition<-rcpplique_gw_avg(W)
if (expansion_mode=='max')
	partition<-rcpplique_gw_max(W)
if (expansion_mode=='outN')
	partition<-rcpplique_gw_outN(W)
if (expansion_mode=='outW')
	partition<-rcpplique_gw_outW(W)

if (unique_singletonLabels)
	partition<- uniqueSingletonLabels(partition)


result<- list(membership=partition,
	      alg_params= list( expansion_mode=expansion_mode )
	      )

return(result)
}
