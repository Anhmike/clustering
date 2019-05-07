\d .clust

/splitting dimensions of co-ordinates + next splitting dimension, e.g. (x,y) = 0 1 0
i.dim:{count first x}

/valid entries of a kd-tree
i.val:{select from x where valid}

/true if number of clusters in a kd-tree > desired number of clusters (cl)
i.cn1:{x<exec count distinct clt from y}
i.cn2:{x<count distinct exec cltidx from y where valid}

/same output
i.cl2tab:{`idx xasc flip`idx`clt!raze each(x;(count each x)#'min each x:exec distinct cltidx from x where valid)}
i.rtab:{update pts:x from @[i.cl2tab;y;{[x;y]`idx`clt`pts#y}[;y]]}


/2 closest clusters in a kd-tree
i.closclust:{ 
 a:first select from x where nnd=min nnd;
 b:first exec clt from x where idx=a`nni;
 select from x where clt in(b,a`clt)}

/index of nearest neighbours in kd-tree to cluster
i.nnidx:{[t;cl]exec initi except cl`initi from t where nni in cl`idx} 

/updated kd-tree with new closest distance/idx for cluster and nearest neighbours
/* t    = tree
/* dist = closest distances to cluster
/* idxs = closest index to cluster
/* nn   = nearest neighbours to cluster
/* lf   = linkage function
/* ii   = initial index
/* rp   = representative points
/* sd   = splitting dimensions

i.upd:{[t;dist;idxs;lf;ii;rp;sd]
 v:i.val t;
 cd:dm im:kd.i.imin dm:kd.i.ld[lf;1]dist;                     
 ci:v[`idx]im;                                      
 t:kd.insertcl[sd]/[t;rp;ii;first idxs];
 updp:v[`idx]n:raze where{x>y}'[v`nnd;dm]; 
 if[0<count updp;t:{[n;p;dm;t;k]   
  update nnd:dm[n k],nni:enlist max t`idx from t where idx=p k
  }[n;updp;dm]/[t;til count updp]]; 
 update nni:ci,nnd:cd from t where clt=first idxs
 }

/distance calulation (x) between valid points in tree (y) and points in cluster (z)
i.distc:{kd.i.dd[x]@'/:z-\:/:y`pts}

/representative points for a cluster using CURE - get most spread out and apply compression
/* d  = data points
/* cl = cluster
/* r  = number of representative points
/* c  = compression

i.curerep:{[d;cl;r;c]
 mean:avg d idxs:distinct raze cl`cltidx;
 maxp:idxs kd.i.imax sum each{x*x}mean-/:d idxs;
 rp:d ii:r{z,m kd.i.imax{min{sum k*k:x[z]-x[y]}[x;y]each z}[x;;z]each m:y except z}[d;idxs]/maxp;
 rp:(rp*1-c)+\:c*mean;
 (rp;ii;idxs)}

/representative points for cluster using hierarchical clustering
i.hcrep:{[d;cl;lf]
 rp:kd.i.ld[lf;0][d;idxs:distinct raze cl`cltidx];
 ii:$[lf~`centroid;first idxs;idxs];
 (rp;ii;idxs)}

/updated tree with nearest neighbour distances recalculated
i.recalc:{[df;lf;t;nni;idxs;ii]
 t:kd.distcalc[df;lf]/[t;nni];
 {[c;t;j]update cltidx:c from t where initi=j,valid}[enlist idxs]/[t;ii]}

/insert new pt into the tree
i.newpt:{[t;idxs;rp;sd;ii;df;lf;nn]
 t:kd.deleteN/[t;idxs];
 v:select from t where valid;
 dist:{x each y}[df]each flip rp-\:/:v`pts;
 cd:dm im:kd.i.imin dm:kd.i.ld[lf;1]dist;
 ci:v[`idx]im;
 t:kd.insertcl[sd]/[t;rp;ii;first idxs];
 updp:v[`idx]n:raze where each{x>y}[v`nnd]each dist;
 if[0<count updp;t:{[n;p;dm;t;k]
  update nnd:dm[n k],nni:enlist max t`idx from t where idx=p k
  }[n;updp;dm]/[t;til count updp]];
 upd[t;cd;ci;nn;idxs;df;lf;ii]}

/initial cluster table for complete/average linkage
i.buildtab:{
 d:{(d i;i:first 1_iasc d:kd.i.dd[z]each x-/:y)}[;x;y]each x;
 flip`idx`pts`clt`nni`nnd!(i;x;i:til count x;d[;1];d[;0])}

i.hcupd:{[df;lf;t;cl]
 t:select from t where clt<>cl[`clt];
 dm:$[lf=`ward;raze kd.i.ld[lf;1][cl`n]'[t`n;i.distc[df;t;enlist cl`pts]];kd.i.ld[lf;1]i.distc[df;t;cl`pts]];
 (cl`clt;`cd`ci!(dm;t`clt)@\:kd.i.imin dm)
 }
