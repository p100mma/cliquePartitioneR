\name{cliquePartitioneR-package}
\alias{cliquePartitioneR-package}
\alias{cliquePartitioneR}
\docType{package}
\title{
Partitioning Undirected Weighted Graphs Into Heavily Connected Cliques
}
\description{
Collection of functions for partitioning undirected weighted graphs into cliques,
         in a manner that takes into the account the weights of the edges. Designed to work 
	 mostly with dense graphs. Clique paritioning of such graphs
	 that also exhibit a correlation of weighted and plain degree is a meaningful way
	 of simplifying the analysis of such graphs, especially if clique parititioning
	 takes weights into the account. For the case of such graphs, simple greedy algorithms
	 can be used due to the weighted - unweighted degree correlations. Routines for such
	 algorithms were implemented in C++.
}
\details{
As an input one must provide a symmetric numeric matrix encoding similarities between nodes, where 0 means no similarity. 
Main function to use is \code{greedyCliquePartitioner}. 
By convention, clique assignments are encoded as membership vectors.
Singleton nodes not belonging to any cliques are given a label of \code{0}.
Package contains also auxilary functions for normalizing labels and for checking if given membership vector encodes a valid clique partition of a given input graph.
}
\author{
Piotr Stomma 

Maintainer: Piotr Stomma <p.stomma@uwb.edu.pl>
}
\keyword{ clique }
