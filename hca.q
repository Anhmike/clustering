\l ../kdtree.q
\d .kd

/min and max indices of a list
imax:{x?max x};
imin:{x?min x};

/run average hierarchical clustering alogrithm for data in a kd tree
hc.clust:{[sample;dim;kd]
 old:select from kd where valid,closDist=min closDist,idx in(idx[0];closIdx[0]);
 cpts:exec initi from kd where valid,closIdx in old`idx;
 mean:avg sample idxs:distinct raze old`clustIdx;
 updkd:tree.deleteN/[kd;idxs];
 v:select from updkd where valid;
 d:sum each x*x:mean-/:v`rep;
 p:v[`idx]n:where d<v`closDist;
 m:v[`idx]d?dmin:min d;
 updkd:tree.insertKd[dim]/[updkd;enlist mean;i0;i0:first idxs];
 updkd:$[0=count p;updkd;
  {[n;p;d;kd;k]update closDist:d[n k],closIdx:enlist max kd`idx from kd where idx=p k
   }[n;p;d]/[updkd;til count p]];
 updkd:update closIdx:first m,closDist:dmin from updkd where clust=i0;
 newc:exec idx from updkd where clust=i0;
 newp:exec idx except newc from updkd where valid,initi in cpts;
 updkd:tree.distC/[updkd;newp];
 update clustIdx:enlist idxs from updkd where valid,initi=i0
 }

/cluster a sample of data into n clusters using average heirarchical clustering
hc.hc:{[sample;nclust]
 dim:(til count first sample),0;
 t:tree.createTree[sample;dim];
 hctab:{[nclust;kd]nclust<count distinct exec clust from kd where valid}[nclust]hc.clust[sample;dim]/t;
 distinct exec clustIdx from hctab where valid
 }