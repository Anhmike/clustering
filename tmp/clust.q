\l kdtree.q
\l kd-util.q
\l utils.q
\l ccode/cure.q
\d .clust

/hierarchical clustering
/* d  = data
/* cl = number of clusters
/* df = distance function/metric
/* lf = linkage function

hc:{[d;cl;df;lf]
 t:$[b:lf in`complete`average`ward;i.buildtab[d;df];kd.buildtree[d;til count d;sd:i.dim d;df]];
 i.rtab[d]$[lf~`ward;i.cn1[cl]algow[df;lf]/@[t;`nnd;%;2];b;i.cn1[cl]algoca[df;lf]/t;
  lf~`single;i.cn2[cl]algos[d;df;lf;sd]/t;i.cn2[cl]algocc[d;df;df;lf;sd;0b]/t]}

/CURE algorithm
/* r = number of representative points
/* c = compression

cure:{[d;cl;df;r;c;b;s]
 $[b;[cst:.cure.cure[r;c;cl;flip d];
 ([]idx:til count d;clt:{where y in'x}[cst]each til count d;pts:d)];
 [t:kd.buildtree[d;til count d;sd:i.dim d;df];
 $[s;{select from x where valid};i.rtab[d]]i.cn2[cl]algocc[d;df;r;c;sd;1b]/t]]}

/DBSCAN algorithm
/* d = data
/* p = minimum number of points per cluster
/* e = epsilon value

dbscan:{[d;p;e]
 dm:i.epdistmat[e;flip d]'[d;k:til count d];
 t:([]idx:k;dist:dm;clt:0N;valid:1b);
 i.rtabdb[d]{0N<>x 1}algodb[p]/(t;0;0)}

/CURE/centroid - merge two closest clusters and update distances/indices
/* x1 = r (CURE) or df (centroid)
/* x2 = c (CURE) or lf (centroid)
/* sd = splitting dimension
/* b  = 1b (CURE) or 0b (centroid)

algocc:{[d;df;x1;x2;sd;b;t]
 cl:i.closclust i.val t;
 rep:$[b;(i.curerep[d;idxs;x1;x2];idxs:distinct raze cl`cltidx;cl`initi);i.hcrep[d;cl;x2]];
 t:kd.insertcl[sd]/[t;rp;ii:first idxs;(count rp:rep 0)#enlist idxs:rep 1];
 t:kd.deletecl[df]/[t;rep 2];
 kd.distcalc[df]/[t;exec idx from t where clt in ii,valid]}

/Single - merge two closest clusters and update distances/indices
algos:{[d;df;lf;sd;t]
 cl:i.closclust t;
 i0:first idxs:distinct raze cl`cltidx;
 t:update clt:i0 from t where idx in cl`idx;
 t:kd.distcalc[df]/[t;cl`idx];
 {[c;t;j]update cltidx:c from t where initi=j,valid}[enlist idxs]/[t;idxs]}

/Complete/average - merge two closest clusters and update distances/indices
algoca:{[df;lf;t]
 cd:c,(t c:kd.i.imin t`nnd)`nni;
 t:update clt:min cd from t where clt=max cd;
 nn:0!select pts by clt from t where nni in cd;
 du:i.hcupd[df;lf;t]each nn;
 {[t;x]![t;enlist(=;`clt;x 0);0b;`nnd`nni!x 1]}/[t;du]}

/Ward - merge two closest clusters and update distances/indices
algow:{[df;lf;t]
 cd:c,d:(t c:kd.i.imin t`nnd)`nni;
 t:update clt:min cd from t where clt=max cd;
 p:sum exec count[i]*first pts by pts from t where clt=min cd;
 t:update pts:count[i]#enlist[p%count[i]]by clt from t where clt=min cd;
 ct:0!select n:count i,first pts,nn:any nni in cd by clt from t;
 du:i.hcupd[df;lf;ct]each select from ct where nn;
 {[t;x]![t;enlist(=;`clt;x 0);0b;`nnd`nni!x 1]}/[t;du]}

/DBSCAN
algodb:{[p;l]
 cl:{0<>sum type each x 1}i.dbclust[c:l 2;p]/(l 0;l 1);
 nc:first exec idx from t:cl 0 where valid;
 (t;nc;1+c)}

/cluster new points
clustnew:{
 cl:$[z;raze whichcl[x;exec idx from x where dir=2]each y;
  x[`clt]{.clust.kd.i.imin sum each k*k:y-/:x}[x`pts]each y];
 ([]pts:y;clt:cl)}
whichcl:{ind:@[;2]{0<count x 1}kd.bestdist[x;z;0n;`e2dist]/(0w;y;y;y);exec clt from x where idx=ind}