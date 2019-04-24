\d .clust

/insert 1st clust as root, then remaining clusts and indices
kd.createTree:{[d;sd;df;lf]
 r:flip`idx`initi`rep`dir`dim`parent`clust`clustIdx`valid!
  (0;0;enlist @[d;0];2;0;0;0;enlist 0 0;1b);
 t:kd.insertKd[sd]/[r;1_d;cl;cl:1_til count d];
 t:update clustIdx:enlist each til count d from t;
 kd.distC[df;lf]/[t;t`initi]
 }

/insert new clust into tree - check if L or R of init clust, insert looking at sd of each node, upd info of new node in tree
kd.insertKd:{[sd;t;d;l;cl]
 nsd:{0<=first @[x;0]`idx}nodedir[d;t]/enlist t 0;
 insertt[t;nsd 1;nsd;d;l;cl;sd]
 }

/calculate distances between clusters using tree
kd.distCalc:{[t;p;cl;df;lf;bd]
 nn:bd 1;
 newn:select rep,idx from t where idx in nn,valid,clust<>cl;
 newD:imins,newn[`idx] ii?imins:min ii:{[p;x]sum m*m:x-p}[p]each newn`rep;
 if[(newD[0]<bd[0])&count[newn]<>0;bd[0]:newD[0];bd[2]:newD[1]];
 axisD:(raze splitdim[t;bd;p;df]each nn)except bd[3]:bd[3],nn;
 (bd 0;distinct axisD;bd 2;bd 3)
 }
 
/upd distances and new closDist/idx in tree
kd.distC:{[df;lf;t;pt] 
 idpts:select parent,clust,rep from t where idx=pt,valid;
 dist:{0<count x[1]}
  kd.distCalc[t;first idpts`rep;first idpts`clust;df;lf]/
  (0W;(raze idpts[`parent],raze exec idx from t where parent=pt,valid)except pt;pt;pt);
 update closDist:dist[0],closIdx:dist[2] from t where idx=pt
 }

/delete node from tree
kd.delN:{
 t:x[0];X:x[1];
 nsd:$[ii:0=count ll:exec idx from t where dir=1,parent=X,valid;
  first exec idx from t where parent=X,dir=0,valid;
  first ll];
 if[ii;t:update dir:1 from t where idx=nsd];
 child:raze{[t;x]0<>count exec idx from t where parent=first x,valid}[t]
  branches[t]\nsd;
 newNode:mindim[t;X;child];
 tree:updatet[t;newNode;X];
 (update closIdx:X from tree where closIdx=first newNode`idx,valid;first newNode`idx)
 }

/ upd tree once pt deleted - repeat process moving pts up the tree until pt is reached with no children then delete node with no children
kd.deleteN:{[t;X]
 delCl:{0<>count select from x[0]where parent=x[1],valid}kd.delN/(t;first exec idx from t where initi=X,valid);
 update valid:0b from first delCl where idx=last delCl
 }
