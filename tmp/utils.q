\l kdtree.q
\d .clust

/takes in data and returns splitting dimensions + next splitting dimension
dim:{count first x}

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

/insert new pt into the tree
newpt:{[t;idxs;rp;sd;ii;df;lf;nn]
 t:kd.deleteN/[t;idxs];
 v:select from t where valid;
 dist:{x each y}[df]each flip rp-\:/:v`rep;         /dist between rp&other cl
 cd:dm im:imin dm:ld[lf;1]dist;                     /d of closest pt to cl
 ci:v[`idx]im;                                      /idx of closest pt to cl
 t:kd.insertKd[sd]/[t;rp;ii;first idxs];         /insert new cl
 updp:v[`idx]n:raze where each{x>y}[v`closDist]each dist; /pts w/ closDist>d to rp;
 if[0<count updp;t:{[n;p;dm;t;k]                    /if closer pts, upd tree
  update closDist:dm[n k],closIdx:enlist max t`idx from t where idx=p k
  }[n;updp;dm]/[t;til count updp]];
 upd[t;cd;ci;nn;idxs;df;lf;ii]}

curerep:{[d;cl;r;c]
 mean:avg pts:d idxs:distinct raze cl`clustIdx; /mean of rep pts
 maxFromMean:idxs imax sum each{x*x}mean-/:d idxs; /max pts from mean
 rp:distinct d ii:r{[samp;idxs;nn]              /get most spread out rep pts
  nn,maxI imax{[samp;nn;maxI]
   min{[samp;nn;maxI]
    sum x*x:samp[maxI]-samp[nn]
   }[samp;maxI]each nn
  }[samp;nn]each maxI:idxs except nn
 }[d;idxs]/maxFromMean;
 rp:(rp*1-c)+\:c*mean; /apply comp to rp
 (rp;ii;idxs)
 }

hcrep:{[d;cl;lf]
 rp:ld[lf;0][d;idxs:distinct raze cl`clustIdx]; /rep pts
 ii:$[lf~`centroid;[first idxs];idxs]; /if cent, new cl-> 1row
 (rp;ii;idxs)
 } 
