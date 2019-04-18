\l ../kdtree.q
\l utils.q
\d .clust

/data, num rep pts, comp, splitting dim, kd-tree
cureclust:{[d;r;c;sd;t]
 cl:closClust t;
 nn:nnidx[t;cl];

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
 
 t:kd.deleteN/[t;idxs];

 v:select from t where valid;
 dist:{min sum each k*k:x-\:y}[rp]each exec rep from v; /dist between rp&other cl
 ci:v[`idx]where dist=cd:min dist;                  /d&idx of closest pt to cl
 t:kd.insertKd[sd]/[t;rp;ii;first idxs];            /insert new cl
 updp:v[`idx]n:where dist<v`closDist;               /pts w/ closDist>d to rp
 if[0<count updp;t:{[rd;t;n;updp]                   /if closer pts, upd tree
  update closDist:rd[n],closIdx:enlist max[t`idx]from t where idx=updp
  }[dist]/[t;n;updp]];

 upd[t;cd;first ci;nn;idxs;edist2;`centroid;ii]}

/data, num clust, num rep pts, comp  
cure:{[d;cl;r;c]
 sd:dim d;                                     
 t:kd.createTree[d;sd;edist2;`centroid]; /df & lf hard coded
 p:nclust[cl]cureclust[d;r;c;sd]/t;
 d distinct exec clustIdx from p where valid
 }


















