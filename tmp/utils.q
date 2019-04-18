\l ../kdtree.q
\d .clust

/takes in data and returns splitting dimensions + next splitting dimension
dim:{til[count first x],0}

/takes in kd-tree and num clust, returns true is num clust in t>desired num clust
nclust:{[cl;t]cl<count distinct exec clustIdx from t where valid}

/takes in a kd-tree, returns clust of 2 closest pts
closClust:{
 v:select from x where valid;
 a:first select from v where closDist=min closDist; /cl with min closDist
 b:first exec clust from v where idx=a`closIdx;     /cl closest to a
 select from x where clust in(b,a`clust)}           /cl of a and b

/takes in a kd-tree and cluster, returns idx of nearest neighbours
nnidx:{[t;cl]exec initi except cl`initi from t where valid,closIdx in cl`idx} 

/tree, closDist+closIdx to cl, nn, cl idxs, dist func, link func, initial idx
/upd closest dets for cl and dist for nn, returns upd tree
upd:{[t;cd;ci;nn;idxs;df;lf;ii]
 t:update closIdx:ci,closDist:cd from t where clust=first idxs;
 nni:exec idx from t where valid,initi in nn;
 t:kd.distC[df;lf]/[t;nni];
 {[c;t;j]update clustIdx:c from t where initi=j,valid}[enlist idxs]/[t;ii]}