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

 t:$[b:lf in`complete`average`ward;i.buildtab[d;df];kd.buildtree[d;sd:i.dim d;df;lf]];
 $[lf~`ward;i.cn1[cl]algow[df;lf]/@[t;`nnd;%;2];b;i.cn1[cl]algoca[df;lf]/t; /`ward`complete`average

  lf~`single;i.cn2[cl]algos[d;df;lf;sd]/t;i.cn2[cl]algocc[d;df;lf;sd;0b]/t]} /`single`centroid

/CURE algorithm
/* r = number of representative points
/* c = compression

cure:{[d;cl;r;c]
 t:kd.buildtree[d;sd:i.dim d;`e2dist;`single];
 i.cn2[cl]algocc[d;r;c;sd;1b]/t}

/CURE/centroid - merge two closest clusters and update distances/indices
/* x1 = r (CURE) or df (centroid)
/* x2 = c (CURE) or lf (centroid)
/* sd = splitting dimension
/* b  = 1b (CURE) or 0b (centroid)

algocc:{[d;x1;x2;sd;b;t]
 cl:i.closclust v:i.val t;
 rep:$[b;i.curerep[d;cl;x1;x2];i.hcrep[d;cl;x2]];
 nn:i.nnidx[v;cl];
 t:kd.deletecl/[t;idxs:rep 2];
 if[b;x1:`e2dist;x2:`single];
 dist:i.distc[x1;i.val t;rp:rep 0];
 t:i.upd[t;dist;idxs;x2;ii:rep 1;rp;sd];
 nni:exec idx from t where initi in nn,valid;
 i.recalc[x1;x2;t;nni;idxs;ii]}

/Single - merge two closest clusters and update distances/indices
algos:{[d;df;lf;sd;t]
 cl:i.closclust t;
 i0:first idxs:distinct raze cl`cltidx;
 t:update clt:i0 from t where idx in cl`idx;
 i.recalc[df;lf;t;cl`idx;idxs;idxs]}

/Complete/average - merge two closest clusters and update distances/indices
algoca:{[df;lf;t]
 cd:c,(t c:kd.i.imin t`nnd)`nni;
 t:update clt:min cd from t where clt=max cd;
 nn:0!select pts by clt from t where nni in cd;
 du:i.hcupd[df;lf;t]each nn;
 {[t;x]![t;enlist(=;`clt;x 0);0b;`nnd`nni!value x 1]}/[t;du]}

/Ward - merge two closest clusters and update distances/indices
algow:{[df;lf;t]
 cd:c,d:(t c:kd.i.imin t`nnd)`nni;
 t:update clt:min cd from t where clt=max cd;
 p:sum exec count[i]*first pts by pts from t where clt=min cd;
 t:update pts:count[i]#enlist[p%count[i]]by clt from t where clt=min cd;
 ct:0!select n:count i,first pts,nn:any nni in cd by clt from t;
 du:i.hcupd[df;lf;ct]each select from ct where nn;
 {[t;x]![t;enlist(=;`clt;x 0);0b;`nnd`nni!value x 1]}/[t;du]}