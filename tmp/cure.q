\l kdtree.q
\l kd-util.q
\l utils.q
\d .clust

/data, num rep pts, comp, splitting dim, kd-tree
cureclust:{[d;r;c;sd;t]
 cl:closClust t;
 nn:nnidx[t;cl];
 
 rep:curerep[d;cl;r;c];

 npt:newpt[t;rep[2];rep[0];sd;rep[1];edist2;`single;nn]}

/data, num clust, num rep pts, comp  
cure:{[d;cl;r;c]
 sd:dim d;                                     
 t:kd.createTree[d;sd;edist2;`centroid]; /df & lf hard coded
 p:nclust[cl]cureclust[d;r;c;sd]/t;
 d distinct exec clustIdx from p where valid
 }


















