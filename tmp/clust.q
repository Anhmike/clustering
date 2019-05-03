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
 t:$[b:lf in`complete`average`ward;i.buildtab[d;df];kd.buildtree[d;sd:i.dim d;df]];
 $[lf~`ward;i.cn1[cl]algow[df;lf]/@[t;`nnd;%;2];b;i.cn1[cl]algoca[df;lf]/t;  /`ward`complete`average
  lf~`single;i.cn2[cl]algos[d;df;lf;sd]/t;i.cn2[cl]algocc[d;df;lf;sd;0b;df]/t]} /`single`centroid

/CURE algorithm
/* r = number of representative points
/* c = compression

cure:{[d;cl;r;c;df]
 t:kd.buildtree[d;sd:i.dim d;df];
 i.cn2[cl]algocc[d;r;c;sd;1b;df]/t}

/CURE/centroid - merge two closest clusters and update distances/indices
/* x1 = r (CURE) or df (centroid)
/* x2 = c (CURE) or lf (centroid)
/* sd = splitting dimension
/* b  = 1b (CURE) or 0b (centroid)

algocc:{[d;x1;x2;sd;b;df;t]
 cl:i.closclust i.val t;
 rep:$[b;i.curerep[d;cl;x1;x2];i.hcrep[d;cl;x2]];
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
 dd:distc[df]'[tc:{[x;y]select rep,clust from x where clt<>y}[t]each k:key nn;value nn];
 cd:dm@'im:imin each dm:{$[1=y;raze;ld[x;1]]z}[lf]'[value count each nn;dd];
 im:(tc@\:`clust)@'im;
 {[t;x;y;z]![t;enlist(=;`clt;x);0b;`nnd`nni!y,z]}/[t;k;cd;im]}

/Ward - merge two closest clusters and update distances/indices
algow:{[df;lf;t]
 cd:c,d:(t c:kd.i.imin t`nnd)`nni;
 t:update clt:min cd from t where clt=max cd;
 p:sum exec count[i]*first pts by pts from t where clt=min cd;
 t:update pts:count[i]#enlist[p%count[i]]by clt from t where clt=min cd;
 ct:0!select n:count i,first pts,nn:any nni in cd by clt from t;
 du:i.hcupd[df;lf;ct]each select from ct where nn;
 {[t;x]![t;enlist(=;`clt;x 0);0b;`nnd`nni!value x 1]}/[t;du]}
