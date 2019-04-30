\d .clust

/takes in data and returns splitting dimensions + next splitting dimension
dim:{count first x}

val:{select from x where valid};

/takes in kd-tree and num clust, returns true is num clust in t>desired num clust
nclust:{[cl;t]cl<count distinct exec clustIdx from t where valid}

/takes in a kd-tree, returns clust of 2 closest pts
closClust:{ 
 a:first select from x where closDist=min closDist; /cl with min closDist
 b:first exec clust from x where idx=a`closIdx;     /cl closest to a
 select from x where clust in(b,a`clust)}           /cl of a and b

/takes in a kd-tree and cluster, returns idx of nearest neighbours
nnidx:{[t;cl]exec initi except cl`initi from t where closIdx in cl`idx} 

/tree, closDist+closIdx to cl, nn, cl idxs, dist func, link func, initial idx
/upd closest dets for cl and dist for nn, returns upd tree
upd:{[t;dist;idxs;df;lf;ii;rp;sd]
 v:val t;
 cd:dm im:imin dm:ld[lf;1]dist;                     
 ci:v[`idx]im;                                      
 t:kd.insertKd[sd]/[t;rp;ii;first idxs];
 updp:v[`idx]n:raze where each{x>y}[v`closDist]each dist; 
 if[0<count updp;t:{[n;p;dm;t;k]   
  update closDist:dm[n k],closIdx:enlist max t`idx from t where idx=p k
  }[n;updp;dm]/[t;til count updp]]; 
 update closIdx:ci,closDist:cd from t where clust=first idxs
 }

distc:{[v;rp;df]{x each y}[df]each flip rp-\:/:v`rep}

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

