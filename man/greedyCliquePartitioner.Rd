% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/greedyCliquePartitione.R
\name{greedyCliquePartitioner}
\alias{greedyCliquePartitioner}
\title{Greedily partition input weighted graph into heavily connected cliques}
\usage{
greedyCliquePartitioner(
  W,
  check_symmetry = FALSE,
  expansion_mode = "basic",
  unique_singletonLabels = FALSE
)
}
\arguments{
\item{W}{A symmetric matrix encoding weights between nodes, i.e. \code{W[i,j]} is the weight of the edge between nodes \code{i,j}. Weights are interpreted as similarities. \code{W[i,j]==0} implies no connection between \code{i,j}.}

\item{check_symmetry}{Logical telling whether to check symmetry of the input matrix. Not done by default.}

\item{expansion_mode}{Charecter string encoding one of possible ways the algorithm shall build a clique based on the strength of the connections of the neighbors.}

\item{unique_singletonLabels}{Logical telling whether isolated nodes should carry unique labels in the output partition. By default, each singleton gets a label of \code{0}.}
}
\value{
A \code{list} containing following two elements:
\itemize{
\item \code{membership} -- a integer vector of length \code{ncol(W)} giving labels of cliques in which each node belongs, i.e. \code{membership[[i]]} is label of the clique of node \code{i}
\item \code{alg_params} --- a list containing parameters of used greedy clique partitioner algorithm
}
}
\description{
Greedily partition input weighted graph into heavily connected cliques
}
\details{
Algorithm partitions nodes of the network by the following procedure. First, it picks the heaviest edge
available as a starting point for the 1st clique. Then the clique is expanded one node at a time,
by joining a non-clique breaking node according to the \code{expansion_mode} argument.
Clique is expanded until exhaustion, and a new clique is initialized on the graph excluding the cliques already
built.

If \code{expansion_mode} is set to \code{basic}, then the newcomers are neighbors of a seed node of the clique,
that is one of the nodes belonging to initializer edge (the one with the bigger weighted degree). Those
neighbors are chosen in the order of their connection strength with the seed node.

If \code{expansion_mode} is set to \code{average} or \code{max}, then the newcomers are chosen by maximizing 
aggregated strength of the connections of newcomer to the members of the clique. Aggregation function used
for \code{average} is mean, while for \code{max}, the aggregated connection strength of newcomer is set to
the max over all of its connections to the clique.

Lastly, if \code{expansion_mode} is set to \code{outN} or \code{outW}, then newcomers are chosen by maximizing
their outer connectivity. For \code{outN}, procedure picks newcomer with maximal number of connections to nodes
outside of the clique. For \code{outW}, sum of weights of the connections of the newly joined node to the nodes
outside of the clique is maximized.
}
\examples{
input_graph= matrix(runif(40*40)*sample(c(0,1),40*40,replace=TRUE, prob=c(0.7,0.3)),ncol=40)
input_graph= input_graph \%*\% t(input_graph) #generate symmetric matrix
greedyCliquePartitioner( input_graph)
}
