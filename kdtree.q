\d .kd

/build tree - insert first clust as root then insert the remain clusters and their indices
tree.createTree:{[sample;diml]
 root:flip `idx`initi`rep`dir`dim`parent`clust`clustIdx`valid!(0;0;enlist @[sample;0];2;0;0;0;enlist (0 0);1b);
 kds:tree.insertKd[diml]/[root;1_sample;cl;cl:1_til count sample];
 kds:update clustIdx:enlist each til count[sample] from kds;
 tree.distC/[kds;kds`initi]
 }

/insert new clust into tree
tree.insertKd:{[diml;kd;samp;L;cl]
 /check if clust is to L or R of initial clust in tree
 dirn:{0<=first @[x;0]`idx}{[samp;kd;nn]
  /a is prev pt, i is next pt to do split on
  a:@[nn;0];
  /insert clust by looking at splitting dimension of each node
  dir1:$[samp[first a`dim]>(raze a`rep)[first a`dim];1;0];
  i:select from kd where dir=dir1,parent=first a`idx,valid;
  (i;a;dir1)}[samp;kd]/(i;i:kd[0]);
 par:dirn[1];
 /update the info of new node into tree
 $[not 0b in kd[`valid];
  kd upsert([]idx:max[kd`idx]+1;initi:L;clust:L;rep:enlist samp;dim:diml par[`dim]+1;valid:1b;parent:par`idx;dir:dirn[2]);
  update idx:(max kd`idx)+1,initi:L,clust:cl,rep:enlist samp,valid:1b,dim:diml par[`dim]+1,dir:dirn[2],parent:par`idx from kd where idx=first exec idx from kd where not valid]
 }

/calculate dist between clusters using tree
tree.distCalc:{[kd;query;cl;bestD]
 /nodes to search that arent in the same cluster;bd[1]=next pts to search
 newn:nn where{[cl;kd;x](first exec clust from kd where idx=x,valid)<>cl}[cl;kd]each nn:bestD[1];
 /get minimum dist of all searched nodes
 newD:imins,newn ii?imins:min ii:{[kd;query;x]sum m*m:(first exec rep from kd where idx=x,valid)-query}[kd;query]each newn;
 /if new dist is less than current best dist, then that becomes new best dist
 $[(newD[0]<bestD[0])&count[newn]<>0;(bestD[0]:newD[0];bestD[2]:newD[1]);];
 /calc dist between node and pts based on split dim
 /if <= bestD then search the leafs
 /go up the tree, to L/R depending on position
 axisD:(raze{[kd;bestD;query;nn]
  ll:select from kd where idx=nn,valid;
  dirn:$[(qdim:query[dd])<rdim:(first ll[`rep])[dd:first ll`dim];0;1];
  $[bestD[0]>=m*m:rdim-qdim;exec idx from kd where parent=nn,valid;exec idx from kd where parent=nn,dir=dirn,valid],ll`parent
  }[kd;bestD;query]each nn)except bestD[3]:bestD[3],nn; /bestD[3]=nodes already searched
 (bestD[0];distinct axisD;bestD[2];bestD[3])
 }
 
/upd distances
tree.distC:{[kd;pt] 
  idpts:select parent,clust,rep from kd where idx=pt,valid;
  dist:{0<>count @[x;1]}tree.distCalc[kd;first idpts`rep;first idpts`clust]/
   (0W;(raze idpts[`parent],raze exec idx from kd where parent=pt,valid)except pt;pt;pt);
 /update new closDist and idx in tree
 update closDist:dist[0],closIdx:dist[2] from kd where idx=pt
 }

/delete node from tree
tree.delN:{
 /tree and pt (X) to be deleted
 kd:x[0];
 X:x[1];
 dirn:$[ii:0=count ll:exec idx from kd where dir=1,parent=X,valid;
  first exec idx from kd where parent=X,dir=0,valid;
  first ll];
 /if pt has no R child, L child replaces
 mindim:raze{[kd;x]0<>count exec idx from kd where parent=first x,valid}[kd]
  {[kd;x]raze exec idx from kd where parent in x,valid}[kd]\dirn;
 /min/max value of rep pts from mindim based on splitting dim
 newP:mindim $[ii;imax;imin]raze({[kd;x]first exec rep from kd where idx=x,valid
  }[kd]each mindim)[;first exec dim from kd where idx=X,valid];
 /get info from tree of new position
 newNode:select from kd where idx=newP,valid;
 /newNode replaces the deleted pt, but dim and idx stay the same
 tree:update rep:newNode`rep,initi:newNode`initi,closDist:newNode`closDist,clust:newNode`clust,clustIdx:newNode`clustIdx, closIdx:newNode`closIdx from kd where idx=X,valid;
 /any pt that had newP as closest idx has to upd closIdx to new position
 (update closIdx:X from tree where closIdx=first newNode`idx,valid;newP)
 }

/upd tree once pt deleted
tree.deleteN:{[kd;X]
 /repeat process moving pts up the tree until pt is reached with no children
 delCl:{0<>count select from x[0] where parent=x[1],valid}tree.delN/(kd;first exec idx from kd where initi=X,valid);
 /delete the node with no children
 update valid:0b from first delCl where idx=last delCl
 }
