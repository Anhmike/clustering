\l kdtree.q
\l kd-util.q
\l utils.q
\d .clust

/data, dist fnc, linkage fnc, splitting dimension, kd-tree
hcclust:{[d;df;lf;sd;t]
 cl:closClust t;
 nn:nnidx[t;cl];

 rep:hcrep[d;cl;lf];

 npt:newpt[t;rep[2];rep[0];sd;rep[1];df;lf;nn]
 }

/data, num clust, dist fnc, linkage fnc
hc:{[d;cl;df;lf]
 sd:dim d;
 t:kd.createTree[d;sd;df;lf];
 p:nclust[cl]hcclust[d;df;lf;sd]/t;
 d distinct exec clustIdx from p where valid
 }
