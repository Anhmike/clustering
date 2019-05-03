\d .clust

/build kd-tree with 1st data point as root and remaining points subsequently added
/* d  = data points
/* sd = splitting dimensions
/* df = distance function/metric
/* lf = linkage function

kd.buildtree:{[d;sd;df;lf]
 r:flip`idx`initi`pts`dir`dim`par`clt`cltidx`valid!(0;0;enlist @[d;0];2;0;0;0;enlist 0 0;1b);
 t:kd.insertcl[sd]/[r;1_d;cl;cl:1_til count d];
 t:update cltidx:enlist each til count d from t;
 kd.distcalc[df;lf]/[t;t`initi]}

/insert new cluster checking if L or R of root and looking at sd of each node
/* t  = kd-tree
/* ii = initial index
/* cl = cluster index

kd.insertcl:{[sd;t;d;ii;cl]
 nsd:{0<=first @[x;0]`idx}kd.i.nodedir[d;t]/enlist t 0;
 kd.i.insertn[t;nsd 1;nsd;d;ii;cl;sd]}

/distance calculation between clusters in kd-tree and single node (pt)
kd.distcalc:{[df;lf;t;pt]
 idpts:select par,clt,pts from t where idx=pt,valid;
 dist:{0<count x 1}kd.bestdist[t;first idpts`pts;first idpts`clt;df;lf]/
  (0W;raze[idpts[`par],raze exec idx from t where par=pt,valid]except pt;pt;pt);
 update nnd:dist 0,nni:dist 2 from t where idx=pt}

/returns list of best distance, points to search, closest index, searched indices
/* p  = index of cluster in the kd-tree
/* bd = current best distance from p to the closest cluster

kd.bestdist:{[t;p;cl;df;lf;bd]
 nn:bd 1;
 newn:select pts,idx from t where idx in nn,valid,clt<>cl;
 newD:imins,newn[`idx] ii?imins:min ii:{[df;p;x]kd.i.dd[df] x-p}[df;p]each newn`pts;
 if[(newD[0]<bd 0)&count[newn]<>0;bd[0]:newD 0;bd[2]:newD 1];
 axisD:raze[kd.i.splitdim[t;bd;p;df]each nn]except bd[3]:bd[3],nn;
 (bd 0;distinct axisD;bd 2;bd 3)}

/updated kd-tree with point X removed by moving points up the tree until X has no children
kd.deletecl:{[t;X]
 delCl:{0<>count select from x[0]where par=x[1],valid}
  kd.delnode/(t;first exec idx from t where initi=X,valid);
 update valid:0b from first delCl where idx=last delCl}

/updated kd-tree and next point to be deleted
kd.delnode:{
 t:x 0;X:x 1;
 nsd:$[ii:0=count ll:exec idx from t where dir=1,par=X,valid;
  first exec idx from t where par=X,dir=0,valid;first ll];
 if[ii;t:update dir:1 from t where idx=nsd];
 child:raze{[t;x]0<>count exec idx from t where par=first x,valid}[t]
  kd.i.branches[t]\nsd;
 newNode:kd.i.mindim[t;X;child];
 tree:kd.i.updatet[t;newNode;X];
 (update nni:X from tree where nni=first newNode`idx,valid;first newNode`idx)}