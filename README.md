# cliqueParitioneR 1.0

Collection of functions for partitioning undirected weighted graphs into cliques,
in a manner that takes into the account the weights of the edges. Designed to work 
	 mostly with dense graphs, even fully connected. Clique paritioning of such graphs
	 that also exhibit a correlation of weighted and plain degree is a meaningful way
	 of simplifying the analysis of such graphs, especially if clique parititioning
	 takes weights into the account. For the case of such graphs, simple greedy algorithms
	 can be used due to the weighted - unweighted degree correlations. Routines for such
	 algorithms were implemented in C++.
