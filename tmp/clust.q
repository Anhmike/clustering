\l kdtree.q
\l kd-util.q
\l utils.q
\d .clust

/hierarchical clustering
/* d  = data
/* cl = number of clusters
/* df = distance function/metric
/* lf = linkage function

hc:{[d;cl;df;lf]
 initcl[d;cl;df;lf;0b]}

/CURE algorithm
/* r = number of representative points
/* c = compression

cure:{[d;cl;r;c]
 initcl[d;cl;r;c;1b]}

/points in each cluster
initcl:{[d;cl;x1;x2;b]
 t:$[bb:(`$string x2)in`complete`average;createtab[d;x1];
    b;kd.createTree[d;sd:dim d;`e2dist;`single];
    kd.createTree[d;sd:dim d;x1;x2]];
 p:$[bb;{x<count distinct clt from y}[cl]hcca[x1;x2]/t;
  x2~`single;nclust[cl]hcsin[d;x1;x2;sd]/t;
  nclust[cl]cluster[d;x1;x2;sd;b]/t];
 d distinct exec clustIdx from p where valid}

/kd-tree with the two closest clusters merged and distances/indices updated
cluster:{[d;x1;x2;sd;b;t]
 v:val t;
 cl:closClust v;
 rep:$[b;curerep[d;cl;x1;x2];hcrep[d;cl;x2]];
 $[b;(df:`e2dist;lf:`single);(df:x1;lf:x2)];
 nn:nnidx[v;cl];
 t:kd.deleteN/[t;idxs:rep 2];
 dist:distc[df;val t;rp:rep 0];
 t:upd[t;dist;idxs;lf;ii:rep 1;rp;sd];
 nni:exec idx from t where initi in nn,valid;
 recalc[df;lf;t;nni;idxs;ii]}

/hc single link - kd-tree with the two closest clusters merged and distances/indices updated
hcsin:{[d;df;lf;sd;t]
 cl:closClust t;
 i0:first idxs:distinct raze cl`clustIdx;
 t:update clust:i0 from t where idx in cl`idx;ii:idxs;
 recalc[df;lf;t;cl`idx;idxs;ii]}