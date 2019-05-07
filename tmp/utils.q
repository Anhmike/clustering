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
 rp:d r{z,m kd.i.imax{min{sum k*k:x[z]-x[y]}[x;y]each z}[x;;z]each m:y except z}[d;idxs]/maxp;
 rp:(rp*1-c)+\:c*mean;
 (rp;idxs;cl`initi)}

/representative points for cluster using hierarchical clustering
i.hcrep:{[d;cl;lf]
 rp:kd.i.ld[lf;0][d;idxs:distinct raze cl`cltidx];
 (rp;idxs;cl`initi)}

/initial cluster table for complete/average linkage
i.buildtab:{
 d:{(d i;i:first 1_iasc d:kd.i.dd[z]each x-/:y)}[;x;y]each x;
 flip`idx`pts`clt`nni`nnd!(i;x;i:til count x;d[;1];d[;0])}

i.hcupd:{[df;lf;t;cl]
 t:select from t where clt<>cl[`clt];
 dm:$[lf=`ward;raze kd.i.ld[lf;1][cl`n]'[t`n;i.distc[df;t;enlist cl`pts]];kd.i.ld[lf;1]i.distc[df;t;cl`pts]];
 (cl`clt;`cd`ci!(dm;t`clt)@\:kd.i.imin dm)
 }

