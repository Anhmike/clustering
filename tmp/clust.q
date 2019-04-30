\l kdtree.q
\l kd-util.q
\l utils.q
\d .clust

/* d  = data
/* cl = number of clusters
/* df = distance function/metric
/* lf = linkage function

hc:{[d;cl;df;lf]
 initcl[d;cl;df;lf;0b]}

/* r = number of representative points
/* c = compression

cure:{[d;cl;r;c]
 initcl[d;cl;r;c;1b]}

initcl:{[d;cl;x1;x2;b]
 t:$[b;kd.createTree[d;sd:dim d;edist2;`single];
  kd.createTree[d;sd:dim d;x1;x2]];
 p:$[x2~`single;nclust[cl]hcsin[d;x1;x2;sd]/t;
  nclust[cl]cluster[d;x1;x2;sd;b]/t];
 d distinct exec clustIdx from p where valid}

cluster:{[d;r;c;sd;b;t]
 v:val t;
 cl:closClust v;
 rep:$[b;curerep[d;cl;r;c];hcrep[d;cl;c]];
 $[b;(df:edist2;lf:`single);(df:r;lf:c)];
 nn:nnidx[v;cl];
 t:kd.deleteN/[t;idxs:rep 2];
 dist:distc[df;val t;rp:rep 0];
 t:upd[t;dist;idxs;df;lf;ii:rep 1;rep 0;sd];
 nni:exec idx from t where initi in nn,valid;
 t:kd.distC[df;lf]/[t;nni];
 {[c;t;j]update clustIdx:c from t where initi=j,valid}[enlist idxs]/[t;ii]}

hcsin:{[d;df;lf;sd;t]
 cl:closClust t;
 i0:first idxs:distinct raze cl`clustIdx;
 t:update clust:i0 from t where idx in cl`idx;
 nni:exec idx from t where closIdx in cl`idx;
 t:kd.distC[df;lf]/[t;nni];
 {[c;t;j]update clustIdx:c from t where initi=j,valid}[enlist idxs]/[t;idxs]}