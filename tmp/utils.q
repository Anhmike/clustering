
\d .clust

/splitting dimensions of co-ordinates + next splitting dimension, e.g. (x,y) = 0 1 0
dim:{count first x}

/valid entries of a kd-tree
val:{select from x where valid}

/true if number of clusters in a kd-tree > desired number of clusters (cl)
nclust:{x<count distinct exec clustIdx from y where valid}

/2 closest clusters in a kd-tree
closClust:{ 
 a:first select from x where closDist=min closDist;
 b:first exec clust from x where idx=a`closIdx;
 select from x where clust in(b,a`clust)}

/index of nearest neighbours in kd-tree to cluster
nnidx:{[t;cl]exec initi except cl`initi from t where closIdx in cl`idx} 

/updated kd-tree with new closest distance/idx for cluster and nearest neighbours
/* t    = tree
/* dist = closest distances to cluster
/* idxs = closest index to cluster
/* nn   = nearest neighbours to cluster
/* lf   = linkage function
/* ii   = initial index
/* rp   = representative points
/* sd   = splitting dimensions

upd:{[t;dist;idxs;lf;ii;rp;sd]
 v:val t;
 cd:dm im:imin dm:ld[lf;1]dist;                     
 ci:v[`idx]im;                                      
 t:kd.insertKd[sd]/[t;rp;ii;first idxs];
 updp:v[`idx]n:raze where{x>y}'[v`closDist;dm]; 
 if[0<count updp;t:{[n;p;dm;t;k]   
  update closDist:dm[n k],closIdx:enlist max t`idx from t where idx=p k
  }[n;updp;dm]/[t;til count updp]]; 
 update closIdx:ci,closDist:cd from t where clust=first idxs
 }

/distance calulation (x) between valid points in tree (y) and points in cluster (z)
/distc:{{x each y}[x]each flip z-\:/:y`rep}
distc:{dd[x]@'/:z-\:/:y`rep}

/representative points for a cluster using CURE
/* d  = data points
/* cl = cluster
/* r  = number of representative points
/* c  = compression

curerep:{[d;cl;r;c]
 mean:avg pts:d idxs:distinct raze cl`clustIdx;
 maxFromMean:idxs imax sum each{x*x}mean-/:d idxs;
 rp:distinct d ii:r{[samp;idxs;nn]     /get most spread out pts as rp
  nn,maxI imax{[samp;nn;maxI]
   min{[samp;nn;maxI]
    sum x*x:samp[maxI]-samp[nn]
   }[samp;maxI]each nn
  }[samp;nn]each maxI:idxs except nn
 }[d;idxs]/maxFromMean;
 rp:(rp*1-c)+\:c*mean;                 /apply comp to rp
 (rp;ii;idxs)}

/representative points for cluster using hierarchical clustering
hcrep:{[d;cl;lf]
 rp:ld[lf;0][d;idxs:distinct raze cl`clustIdx];
 ii:$[lf~`centroid;first idxs;idxs]; /for centroid lf insert clust as 1 row (mean rp)
 (rp;ii;idxs)}

/updated tree with nearest neighbour distances recalculated
recalc:{[df;lf;t;nni;idxs;ii]
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

updw:{[lf;df;cd;t]
 p:sum exec count[i]*first rep by rep from t where clt=min cd;
 t:update rep:count[i]#enlist[p%count[i]]by clt from t where clt=min cd;
 ct:0!select n:count i,first rep,nn:any nni in cd by clt from t;
 nd:pointdist[df;lf;;ct]each uc:select from ct where nn;
 {[t;x;y]![t;enlist(=;`clt;x);0b;`nnd`nni!y[0],y 1]}/[t;uc`clt;nd]
 }

cldist:{[df;lf;x;y](d i;y[`clt]i:imin d:ldw[lf]each dd[df]@'/:raze each x-/:\:/:y`rep)}
pointdist:{[df;lf;x;y](d first i j;first ci j:where x[`clt]<>ci:y[`clt]i:2#iasc d:dd[df]each x[`rep]-/:y`rep)}

ldw:{z%(1%y)+1%x}

