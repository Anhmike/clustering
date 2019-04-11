\l ../kdtree.q
\d .clust

/data, dist fnc, linkage fnc, splitting dimension, kd-tree
hc.clust:{[d;df;lf;sd;t]
 a:first select from t where valid,closDist=min closDist;            /cl with min closDist
 b:first exec clust from t where idx=a`closIdx;                      /cl closest to a
 cl:select from t where clust in(b,a`clust);                         /cl of a and b
 nn:exec initi except cl`initi from t where valid,closIdx in cl`idx; /nearest neighbours to a&b
 rp:ld[lf;0][d;idxs:distinct raze cl`clustIdx];                      /rep pts
 t:kd.deleteN/[t;idxs];                                              /delete rp from current tree
 v:select from t where valid;                                        /get remaining pts in tree
 dist:{x each y}[df]each flip rp-\:/:v`rep;                          /dist between rp&remaining pts
 link:ld[lf;1]each dd:flip dist;                                     /idx of clos rp to each other pt
 cd:dm im:imin dm:{x y}'[dd;link];                                   /d of closest pt to cl
 ci:v[`idx]im;                                                       /idx of closest pt to cl
 $[lf~`centroid;[cent:1b;ii:first idxs];[cent:0b;ii:idxs]];          /if cent lf store new cl as 1row
 t:kd.insertKd[sd]/[t;rp;ii;i0:first idxs];                          /insert new cl
 updp:v[`idx]n:raze where each{x>y}[v`closDist]each dist;            /pts to upd w/ closDist>d to rp
 t:$[0=count updp;t;{[n;p;dm;t;k]                                    /if closer pts, upd tree
  update closDist:dm[n k],closIdx:enlist max t`idx from t where idx=p k
  }[n;updp;dm]/[t;til count updp]];
 t:update closIdx:ci,closDist:cd from t where clust=i0;              /upd clos idx&d of all pts in cl
 nni:exec idx from t where valid,initi in nn;                        /get idx of nearest neighbours
 t:kd.distC/[t;nni;df];                                              /calc new d to cl
 {[c;t;j]                                                            /upd clustIdx to init idx of
  update clustIdx:c from t where initi=j,valid                       /all pts now in cl
  }[enlist idxs]/[t;$[cent;i0;idxs]]
 }

/data, no clust, dist fnc, linkage fnc
hc.hc:{[d;c;df;lf]
 sd:til[count first d],0;  /get splitting dims + next dim, e.g. 2D:0 1 0, 3D:0 1 2 0, etc
 t:kd.createTree[d;sd;df]; /build initial kd-tree and then run HCA until no clust=c
 r:{[c;t]c<count distinct exec clustIdx from t where valid}[c]hc.clust[d;df;lf;sd]/t;
 d distinct exec clustIdx from r where valid  /return pts in each cl
 }