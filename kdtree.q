\d .kd

/utils
imax:{x?max x}
imin:{x?min x}

/distance metrics
mdist:{sum abs x}
edist2:{sum x*x}
edist:{sqrt edist2 x}

/linkage
ld:`single`complete`centroid!(({x y};imin);({x y};imax);({enlist avg x y};imin))

/insert 1st clust as root, then remaining clusts and indices
/*  (d)ata points
/*  (s)plitting (d)imension
tree.createTree:{[d;sd;f]
 r:flip`idx`initi`rep`dir`dim`parent`clust`clustIdx`valid!
  (0;0;enlist @[d;0];2;0;0;0;enlist 0 0;1b);
 t:tree.insertKd[sd]/[r;1_d;cl;cl:1_til count d];
 t:update clustIdx:enlist each til count d from t;
 tree.distC/[t;t`initi;f]
 }

/insert new clust into tree - check if L or R of init clust, insert looking at sd of each node, upd info of new node in tree
/*  (t)ree
/*  l: current idx
/*  current (cl)uster
tree.insertKd:{[sd;t;d;l;cl]
 nsd:{0<=first @[x;0]`idx}{[d;t;nn]a:@[nn;0];
  sd:$[d[first a`dim]>raze[a`rep]first a`dim;1;0];
  i:select from t where dir=sd,parent=first a`idx,valid;
  (i;a;sd)}[d;t]/enlist t 0;      /a:prev pt, i:next pt to do split on
 p:nsd 1;
 $[not 0b in t`valid;
  t upsert([]idx:1+max t`idx;initi:l;clust:l;rep:enlist d;
   dim:sd 1+p`dim;valid:1b;parent:p`idx;dir:nsd 2);
  update idx:1+max t`idx,initi:l,clust:cl,rep:enlist d,valid:1b,dim:sd 1+p`dim,
   dir:nsd 2,parent:p`idx from t where idx=first exec idx from t where not valid]
 }

/calculate distances between clusters using tree
/*  (p)oint to query
/*  current (b)est (d)istance
tree.distCalc:{[t;p;cl;f;bd]
 newn:nn where{[cl;t;x]cl<>first exec clust from t where idx=x,valid}[cl;t]each nn:bd 1;
 newD:imins,newn ii?imins:min ii:{[t;p;f;x]f(first exec rep from t where idx=x,valid)-p}[t;p;f]each newn;
 $[(newD[0]<bd[0])&count[newn]<>0;(bd[0]:newD[0];bd[2]:newD[1]);];
 axisD:(raze{[t;bd;p;f;nn]
  ll:select from t where idx=nn,valid;
  nsd:$[(qdim:p[dd])<rdim:first[ll`rep]dd:first ll`dim;0;1];
  $[bd[0]>=f[rdim-qdim];exec idx from t where parent=nn,valid;
   exec idx from t where parent=nn,dir=nsd,valid],ll`parent
  }[t;bd;p;f]each nn)except bd[3]:bd[3],nn; 
 (bd[0];distinct axisD;bd[2];bd[3])
 }
 
/upd distances and new closDist/idx in tree
tree.distC:{[t;pt;f] 
  idpts:select parent,clust,rep from t where idx=pt,valid;
  dist:{0<>count @[x;1]}tree.distCalc[t;first idpts`rep;first idpts`clust;f]/
   (0W;(raze idpts[`parent],raze exec idx from t where parent=pt,valid)except pt;pt;pt);
 update closDist:dist[0],closIdx:dist[2] from t where idx=pt
 }

/delete node from tree
/*  X: pt to be deleted
tree.delN:{
 t:x[0];
 X:x[1];
 nsd:$[ii:0=count ll:exec idx from t where dir=1,parent=X,valid;
  first exec idx from t where parent=X,dir=0,valid;
  first ll];
 mindim:raze{[t;x]0<>count exec idx from t where parent=first x,valid}[t]
  {[t;x]raze exec idx from t where parent in x,valid}[t]\nsd;
 newP:mindim $[ii;imax;imin]raze({[t;x]first exec rep from t where idx=x,valid
  }[t]each mindim)[;first exec dim from t where idx=X,valid];
 newNode:select from t where idx=newP,valid;
 tree:update rep:newNode`rep,initi:newNode`initi,closDist:newNode`closDist,clust:newNode`clust,clustIdx:newNode`clustIdx, closIdx:newNode`closIdx from t where idx=X,valid;
 (update closIdx:X from tree where closIdx=first newNode`idx,valid;newP)
 }

/ upd tree once pt deleted - repeat process moving pts up the tree until pt is reached with no children then delete node with no children
tree.deleteN:{[t;X]
 delCl:{0<>count select from x[0]where parent=x[1],valid}tree.delN/(t;first exec idx from t where initi=X,valid);
 update valid:0b from first delCl where idx=last delCl
 }