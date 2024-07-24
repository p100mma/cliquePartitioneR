

unlapply<- function(...) unlist(lapply(...))

#' Order values in the decreasing manner
#'
#' @param x A vector to be ordered
#' @return A permutation of indexes of \code{x} that sorts the vector in the decreasing manner.
#' @examples
#' y<- runif(10)*10
#' ord_y<- max2minOrder(y)
#' print(y)
#' print(ord_y) 
#' print(y[ord_y])
#' @export
max2minOrder<- function(x)
{
  stopifnot(is.numeric(x))
  order(-x)
}

#' A minimal, efficient table function for one dimensional vectors
#' 
#' Emulates \code{table} in an efficient manner for a narrow case of one dimensional atomic vectors. Note that no correctness checks are done.
#' 
#' @param values A 1 dimensional atomic vector
#' @return A \code{list} wirth the following elements: 
#' \itemize{
#' \item \code{value} a vector of unique values in \code{values}.
#' \item \code{count} a vector of counts of each element of \code{value} component, i.e. \code{fastTable(x)$count[i]} is frequency of \code{fastTable(x)$value[i]} in \code{x}.
#' 	   }
#' @examples
#' fastTable( sample(c(1,2,3,4),100,replace=TRUE,prob=c(0.2,0.1,0.5,0.2)))
#' @export
fastTable<- function(values)
{
    uqs<-unique(values)
    histo<-unlist(lapply(uqs, function(x) sum(values==x)))
    list(value=uqs,
         count=histo)
}

#' Check if an adjacency or weight matrix encodes a clique
#'
#' Assuming a input matrix encodes ajdacency relation or weights where \code{0} means no connection, function checks whether the nodes described by the matrix are forming a valid clique. Can also be used by specifying a subgraph of the input graph.
#'
#' @param WorA A numeric matrix containing adjacency relationship information between nodes. \code{WorA[i,j]>0} means that nodes \code{i,j} are connected.
#' @param subgr An indexing vector specifying a subraph of the nodes for which to check if they are forming a clique. If \code{NULL}, check is done using all nodes in the input graph.
#' @return \code{TRUE} if the graph in question is a clique, \code{FALSE} otherwise.
#' @examples
#' W<- matrix(runif(40*40), ncol=40) + 1
#' W= W %*% t(W)
#' isClique(W)
#' W[3,1]<-W[1,3]<-0
#' isClique(W)
#' isClique(W, 4:ncol(W))
#' @export
isClique<-function(WorA,subgr=NULL){
    if (!is.null(subgr))  WorA<- WorA[subgr,subgr,
                                        drop=FALSE]
    diag(WorA)<-0
    sum((WorA==0))==nrow(WorA)
}

#' Check if membership vector encodes proper cliques for a given adjacency relationship encoding matrix
#'
#' Assuming a input matrix encodes ajdacency relation or weights where \code{0} means no connection, function checks whether the membership vector provided properly partitions nodes into a set of valid cliques.
#' 
#' @param WorA A numeric matrix containing adjacency relationship information between nodes. \code{WorA[i,j]>0} means that nodes \code{i,j} are connected.
#' @param cl_mem An integer membership vector for wihch \code{cl_mem[i]} gives a label of the clique to which node \code{i} belongs. A value of \code{0} can possibly encode omitted nodes, for which the check is not performed.
#' @return \code{TRUE} if all groups of nodes described by \code{cl_mem} are valid cliques (including pairs and singletons), \code{FALSE} otherwise.
#' @examples
#' W<- matrix(runif(40*40), ncol=40) +1
#' W= W %*% t(W)
#' membership= c( rep(1,3), rep(2, ncol(W) - 3) )
#' areCliques(W, membership)
#' W[3,1]<-W[1,3]<-0
#' areCliques(W, membership)
#' membership[1:3]<-0
#' areCliques(W, membership)
#' @export
areCliques<-function(WorA,cl_mem){
non0<- unique(cl_mem[cl_mem!=0])
all(
unlapply(non0, function(lbl) isClique(WorA, cl_mem==lbl))
    )
}

#' Relabel singleton nodes from \code{0} to unique values
#'
#' Function is interpreted in the contex where input vector describes partition of nodes into groups. In general it performs a simple substitution of \code{0} values by unique integers.
#' 
#' @param cl_mem An integer membership vector for which \code{cl_mem[i]} gives a label of the group to which node \code{i} belongs. A value of \code{0} by assumption encodes singleton nodes. 
#' @return An altered integer membership vector in which each singleton (formerly encoded by \code{0}) gets a new, unique numerical label.
#' @examples
#' mem<- c(0,0,0,1,0,1,2,3,0,1,2,2,2,1,4,4,0,0)
#' uniqueSingletonLabels(mem)
#' @export
uniqueSingletonLabels<- function(cl_mem){
  
  singletonMask<- cl_mem==0
  n_clusters<- max(unique(cl_mem))
  cl_mem[singletonMask] = (n_clusters+1):(sum(singletonMask)+n_clusters)
  return(cl_mem)
}

#' Relabel small groups to \code{0}
#'
#' Based on input label vector, zero out labels corresponding to groups smaller than given minimal group size.
#'
#' @param labelV An integer membership vector for which \code{cl_mem[i]} gives a label of the group to which node \code{i} belongs. A value of \code{0} by assumption encodes nodes that will be not taken into the account by the relabeling operation.
#' @param mS Integer scalar giving minimum group size that is not relabeled to \code{0}.
#' @return An altered integer membership vector in which each label corresponding to the group of size \code{ < mS } is changed to \code{0}.
#' @examples
#' mem<- c(0,0,0,1,0,1,2,3,0,1,2,2,2,1,4,4,0,0)
#' zeroOutSmall(mem,3)
#' @export
zeroOutSmall<- function(labelV, mS) {
rment<- labelV
non0labs<- unique(labelV[labelV!=0])
non0cnts<- unlapply(non0labs, function(lab) sum(labelV==lab))
for (l in seq_along(non0labs)) 
    if (non0cnts[[l]] < mS )
       rment[ labelV==non0labs[[l]] ]=0
return(rment)
}


#' Relabel group labels to consecutive integers
#'
#' Based on input integer label vector, normalize labels such that if there are \code{M} unique labels in total, each of them gets mapped to one of the integers in \code{1:M}.
#'
#' @param labelV An integer membership vector for which \code{cl_mem[i]} gives a label of the group to which node \code{i} belongs. A value of \code{0} by assumption encodes nodes that will be not taken into the account by the relabeling operation.
#' @return An altered integer membership vector in which each label corresponding to the group of size \code{ < mS } is changed to \code{0}.
#' @examples
#' mem<- c(0,0,0,1,0,1,67,3,0,1,67,67,5,5,4,4,0,0)
#' tidyUpLabels(mem)
#' @export
tidyUpLabels<- function(labelV) {
uqn0<-unique(labelV[labelV!=0])
rmentL<- seq_along(uqn0)
rment<- rep(0, length(labelV))
for (l in seq_along(uqn0))
   rment[ labelV == uqn0[[l]] ]= rmentL[[l]]
return(rment)
}

