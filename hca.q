\l ../kdtree.q
\d .kd

/
2D example:
- .hc.hc will run until d has been classified into 4 clusters at which point it will return the points in each cluster

q)d
3.487966 2.617258 
3.052439 2.939565 
2.541804 2.855116 
2.993872 2.741651 
..
q).kd.hc.hc[d;4;.kd.edist2;`centroid]
(2.89733  2.502805   3.487966 2.617258   3.05243..
(-3.148617 3.378614  -3.140483 3.32115   -3.1571..
(3.468053 -3.205529  3.473971 -2.881551  3.43966..
(-2.819643 -2.775133 -2.650737 -2.697036 -2.6885..
\

/run hierarchical clustering alogrithm for data in a kd tree
/* (d)ata
/* (s)plitting (d)imension
/* (d)istance (f)unction: .kd.edist .kd.edist2 .kd.mdist
/* (l)inkage (f)unction: `single`complete`centroid
/* kd-(t)ree

hc.clust:{[d;sd;df;lf;t]
 a:first select from t where valid,closDist=min closDist;
 b:first exec clust from t where idx=a`closIdx;          
 c:select from t where clust in(b,a`clust);             
 cpts:exec initi except c`initi from t where valid,closIdx in c`idx;
 cent:$[lf~`centroid;1b;0b];
 rep:ld[lf;0][d;idxs:distinct raze c`clustIdx];
 updkd:tree.deleteN/[t;idxs];
 v:select from updkd where valid;
 aa:dm im:imin dm:{x y}'[dist;rm:ld[lf;1]each dist:flip dd:{x each y}[df]each flip rep-\:/:v`rep];
 ci:v[`idx]im;
 p:v[`idx]n:raze where each{x>y}[v`closDist]each dd;
 updid:$[cent;first idxs;idxs];
 updkd:tree.insertKd[sd]/[updkd;rep;updid;i0:first idxs];
 updkd:$[0=count p;updkd;
  {[n;p;dm;t;k]update closDist:dm[n k],closIdx:enlist max t`idx from t where idx=p k
   }[n;p;dm]/[updkd;til count p]];
 updkd:update closIdx:ci,closDist:aa from updkd where clust=i0;
 newp:exec idx from updkd where valid,initi in cpts;
 updkd:tree.distC/[updkd;newp;df];
 {[c;t;X]update clustIdx:c from t where initi=X,valid}[enlist idxs]/[updkd;$[cent;i0;idxs]]
 }

/cluster a sample of data into n clusters using average heirarchical clustering
/* (d)ata
/* no of (c)lusters
/* (d)istance (f)unction: .kd.edist .kd.edist2 .kd.mdist
/* (l)inkage (f)unction: `single`complete`centroid

hc.hc:{[d;c;df;lf]
 sd:til[count first d],0;
 t:tree.createTree[d;sd;df];
 r:{[c;t]c<count distinct exec clustIdx from t where valid}[c]hc.clust[d;sd;df;lf]/t;
 d distinct exec clustIdx from r where valid
 }

/NB: affinity and linkage functions contained in kdtree.q