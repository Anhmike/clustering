\l kdtree.q
\l kd-util.q
\l utils.q
\d .clust

/data, num clust, dist fnc, linkage fnc
hc:{[d;cl;df;lf]
 initcl[d;cl;df;lf;0b]
 }


/data, num clust, num rep pts, comp
cure:{[d;cl;r;c]
 initcl[d;cl;r;c;1b]}



initcl:{[d;cl;x1;x2;b]
 t:$[b;kd.createTree[d;sd:dim d;edist2;`single];
  kd.createTree[d;sd:dim d;x1;x2]];
 p:nclust[cl]cluster[d;x1;x2;sd;b]/t;
 d distinct exec clustIdx from p where valid}

cluster:{[d;r;c;sd;b;t]
 v:val t;
 cl:closClust v;
 nn:nnidx[v;cl];
 rep:$[b;curerep[d;cl;r;c];hcrep[d;cl;c]];
 t:kd.deleteN/[t;idxs:rep[2]];
 $[b;(df:edist2;lf:`single);(df:r;lf:c)];
 dist:distc[val t;rp:rep[0];df];
 t:upd[t;dist;idxs;df;lf;ii:rep[1];rep[0];sd];
 nni:exec idx from t where initi in nn,valid;
 t:kd.distC[df;lf]/[t;nni];
 {[c;t;j]update clustIdx:c from t where initi=j,valid}[enlist idxs]/[t;ii]}

