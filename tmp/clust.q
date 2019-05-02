\l kdtree.q
\l kd-util.q
\l utils.q
\d .clust

/hierarchical clustering
/* d  = data
/* cl = number of clusters
/* df = distance function/metric
/* lf = linkage function

hc:{[d;cl;df;lf]
 t:$[b:lf in`complete`average`ward;createtab[d;df];kd.createTree[d;sd:dim d;df;lf]];
 $[lf~`ward;cn1[cl]nnc[df;lf]/t;b;cn1[cl]hcca[df;lf]/t;   /ward/complete/average
  lf~`single;cn2[cl]hcsin[d;df;lf;sd]/t;                  /single
  cn2[cl]cluster[d;df;lf;sd;0b]/t]}                       /centroid

/CURE algorithm
/* r = number of representative points
/* x2 = compression

cure:{[d;cl;r;c]
 t:kd.createTree[d;sd:dim d;`e2dist;`single];
 cn2[cl]cluster[d;r;c;sd;1b]/t}

/kd-tree with the two closest clusters merged and distances/indices updated
cluster:{[d;x1;x2;sd;b;t]
 v:val t;
 cl:closClust v;
 rep:$[b;curerep[d;cl;x1;x2];hcrep[d;cl;x2]];
 $[b;(df:`e2dist;lf:`single);(df:x1;lf:x2)];
 nn:nnidx[v;cl];
 t:kd.deleteN/[t;idxs:rep 2];
 dist:distc[df;val t;rp:rep 0];
 t:upd[t;dist;idxs;lf;ii:rep 1;rp;sd];
 nni:exec idx from t where initi in nn,valid;
 recalc[df;lf;t;nni;idxs;ii]}

/hc single link - kd-tree with the two closest clusters merged and distances/indices updated
hcsin:{[d;df;lf;sd;t]
 cl:closClust t;
 i0:first idxs:distinct raze cl`clustIdx;
 t:update clust:i0 from t where idx in cl`idx;ii:idxs;
 recalc[df;lf;t;cl`idx;idxs;ii]}

/initial cluster table for complete/average linkage
createtab:{
 d:{(d i;i:first 1_iasc d:dd[z]each x-/:y)}[;x;y]each x;
 flip`ind`rep`clust`nni`nnd!(i;x;i:til count x;d[;1];d[;0])}

/clustered data points for complete/average linkage
hcca:{[df;lf;t]
 cd:c,(t c:imin t`nnd)`nni;
 t:update clust:min cd from t where clust=max cd;
 nn:exec rep by clust from t where nni in cd;
 dd:distc[df]'[tc:{[x;y]select rep,clust from x where clust<>y}[t]each k:key nn;value nn];
 cd:dm@'im:imin each dm:{$[1=y;raze;ld[x;1]]z}[lf]'[value count each nn;dd];
 im:(tc@\:`clust)@'im;
 {[t;x;y;z]![t;enlist(=;`clust;x);0b;`nnd`nni!y,z]}/[t;k;cd;im]}

/clustered data points for ward linkage
nnc:{[df;lf;x]
 cd:c,d:(x c:imin x`nnd)`nni;
 x:update clust:min cd from x where clust=max cd;
 updw[lf;df;cd;x]}
