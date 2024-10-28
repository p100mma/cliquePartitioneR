
n_clqScore<- function(clique_labels){
	length(unique(clique_labels[clique_labels!=0]))
}

complexity_clqScore<- function(clique_labels){
	clique_labels<- uniqueSingletonLabels(clique_labels)
	fastTable(clique_labels)-> tabulation
	# - [log n_cliques! + \sum_cliques log(clique size!)]
	-(lfactorial(length(tabulation$value)) + sum(lfactorial(tabulation$count)))
}


thr_clq_f<- function(t, S, ordS_ltr, clqScoreFun,
		     expansion_mode= expansion_mode
		     ){
  stopifnot((t >= ordS_ltr[[1]]) &&(t <= ordS_ltr[[length(ordS_ltr)]]) )
  diag(S)<-0
  if (t %in% ordS_ltr) t_exact=TRUE else t_exact=FALSE
  if (t_exact) {
    S[S< t]=0
    greedyCliquePartitioner(S, expansion_mode = expansion_mode, unique_singletonLabels = FALSE)$membership -> clqs
    clqScoreFun(clqs) -> y_t
  #  prev_y=-9999
  } else {
       r.s_IDX<-min( which(ordS_ltr > t) )
       r.s<- ordS_ltr[[r.s_IDX]]
       l.s<- ordS_ltr[[r.s_IDX -1 ]]
       S[S<l.s]=0
       greedyCliquePartitioner(S, expansion_mode = expansion_mode, unique_singletonLabels = FALSE)$membership -> clqs
       clqScoreFun(clqs) -> y_l
       S[S<r.s]=0
       greedyCliquePartitioner(S, expansion_mode = expansion_mode, unique_singletonLabels = FALSE)$membership -> clqs
       clqScoreFun(clqs) -> y_r
       A= (y_r - y_l)/(r.s-l.s)
       b= y_l - A*l.s
       y_t= A*t + b
       print(c(t, y_t))
       #if (prev_y==y_t) n_eq_tries= n_eq_tries +1 
  }
  return(y_t)
    }

grid_search_exhaustive<- function( S, n_divisions, clqScoreFun, expansion_mode="basic") {

	  
	S_ltr<- S[lower.tri(S)]
	ord_S_ltr<- sort(unique(S_ltr))

	sr=range(ord_S_ltr)
	seq(from=sr[[1]],to=sr[[length(sr)]], length.out=n_divisions)-> seqq
	X_t=seqq
	Y_t=vector()
	partitions=list()

	Sx=S
	for (x_t in X_t){
	  Sx[Sx < x_t]=0
	  greedyCliquePartitioner(Sx, expansion_mode = expansion_mode, unique_singletonLabels = FALSE)$membership -> clqs
	  Y_t[[length(Y_t)+1]]<- clqScoreFun(clqs)
	  partitions[[ length(partitions)+1]]= clqs
	}

	max_idx<-which.max(Y_t)
	if (max_idx==1) {l_t= X_t[[1]]; r_t=X_t[[2]] 
	} else if (max_idx== length(Y_t)) {
	  l_t= X_t[[length(Y_t)-1]]
	  r_t= X_t[[length(Y_t)]]
	} else {
	  l_t= X_t[[ max_idx -1]]
	  r_t= X_t[[ max_idx +1]]
	}
	max_neighborhood= ord_S_ltr[ (ord_S_ltr >= l_t) & (ord_S_ltr <= r_t) ]
	list( X_t= X_t, Y_t= Y_t, max_idx=max_idx, max_neighborhood= max_neighborhood, partitions=partitions)
	}

precision_search_piecewiseLinear<- function(max_neighborhood, S, clqScoreFun, expansion_mode="basic", tol=0.1) {
	
	rs= range(max_neighborhood)
	l_t= rs[[1]]
	r_t= rs[[2]]
	opt_res=optimize(f = thr_clq_f, interval=c(l_t, r_t), maximum = TRUE,
         S=S, ordS_ltr= max_neighborhood, expansion_mode=expansion_mode, clqScoreFun=clqScoreFun, tol=tol)
	S[ S< opt_res$maximum ] =0
	list(opt_res=opt_res,
	     partition=greedyCliquePartitioner(S, expansion_mode = expansion_mode, unique_singletonLabels = FALSE)$membership
	     )
}


