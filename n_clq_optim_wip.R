


library(HCCsim)
library(cliquePartitioneR)
data(brca)

S<- similarity_matrix(fastPearsonData(brca))
S_ltr<- S[lower.tri(S)]

n_clq_f<- function(t, S, ordS_ltr){
  stopifnot((t >= ordS_ltr[[1]]) &&(t <= ordS_ltr[[length(ordS_ltr)]]) )
  diag(S)<-0
  if (t %in% ordS_ltr) t_exact=TRUE else t_exact=FALSE
  n_eq_tries=0
  if (t_exact) {
    S[S< t]=0
    greedyCliquePartitioner(S, expansion_mode = "basic", unique_singletonLabels = FALSE)$membership -> clqs
    length(unique(clqs)) -> y_t
  #  prev_y=-9999
  } else {
       r.s_IDX<-min( which(ordS_ltr > t) )
       r.s<- ordS_ltr[[r.s_IDX]]
       l.s<- ordS_ltr[[r.s_IDX -1 ]]
       S[S<l.s]=0
       greedyCliquePartitioner(S, expansion_mode = "basic", unique_singletonLabels = FALSE)$membership -> clqs
       length(unique(clqs)) -> y_l
       S[S<r.s]=0
       greedyCliquePartitioner(S, expansion_mode = "basic", unique_singletonLabels = FALSE)$membership -> clqs
       length(unique(clqs)) -> y_r
       A= (y_r - y_l)/(r.s-l.s)
       b= y_l - A*l.s
       y_t= A*t + b
       print(c(t, y_t))
       #if (prev_y==y_t) n_eq_tries= n_eq_tries +1 
  }
  return(y_t)
    }
  
ord_S_ltr<- sort(unique(S_ltr))



sr=range(ord_S_ltr)
seq(from=sr[[1]],to=sr[[length(sr)]], length.out=2)-> seqq
X_t=seqq
Y_t=vector()

Sx=S
for (x_t in X_t){
  Sx[Sx < x_t]=0
  greedyCliquePartitioner(Sx, expansion_mode = "basic", unique_singletonLabels = FALSE)$membership -> clqs
  Y_t[[length(Y_t)+1]]<-length(unique(clqs))
}

max_idx<-which.max(Y_t)
print(Y_t[[max_idx]])
if (max_idx==1) {l_t= X_t[[1]]; r_t=X_t[[2]] 
} else if (max_idx== length(Y_t)) {
  l_t= X_t[[length(Y_t)-1]]
  r_t= X_t[[length(Y_t)]]
} else {
  l_t= X_t[[ max_idx -1]]
  r_t= X_t[[ max_idx +1]]
}

ord_S_ltr = ord_S_ltr[ (ord_S_ltr >= l_t) & (ord_S_ltr <= r_t) ]
rs= range(ord_S_ltr)
l_t= rs[[1]]
r_t= rs[[2]]

optimize(f = n_clq_f, interval=c(l_t, r_t), maximum = TRUE,
         S=S, ordS_ltr= ord_S_ltr, tol= 0.1)->res_opt
res_opt$objective
