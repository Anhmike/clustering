\l ../kdtree.q
\l utils.q
\d .clust

/data, dist fnc, linkage fnc, splitting dimension, kd-tree
hcclust:{[d;df;lf;sd;t]
 cl:closClust t;
 nn:nnidx[t;cl];

 rp:ld[lf;0][d;idxs:distinct raze cl`clustIdx]; /rep pts
 $[lf~`centroid;[cent:1b;ii:first idxs];[cent:0b;ii:idxs]]; /if cent, new cl-> 1row

 t:kd.deleteN/[t;idxs];

 v:select from t where valid;
 dist:{x each y}[df]each flip rp-\:/:v`rep;         /dist between rp&other cl
 cd:dm im:imin dm:ld[lf;1]dist;                     /d of closest pt to cl
 ci:v[`idx]im;                                      /idx of closest pt to cl
 t:kd.insertKd[sd]/[t;rp;ii;i0:first idxs];         /insert new cl
 updp:v[`idx]n:raze where each{x>y}[v`closDist]each dist; /pts w/ closDist>d to rp
 if[0<count updp;t:{[n;p;dm;t;k]                    /if closer pts, upd tree
  update closDist:dm[n k],closIdx:enlist max t`idx from t where idx=p k
  }[n;updp;dm]/[t;til count updp]];

 upd[t;cd;ci;nn;idxs;df;lf;$[cent;i0;idxs]]}

/data, num clust, dist fnc, linkage fnc
hc:{[d;cl;df;lf]
 sd:dim d;
 t:kd.createTree[d;sd;df;lf];
 p:nclust[cl]hcclust[d;df;lf;sd]/t;
 d distinct exec clustIdx from p where valid
 }