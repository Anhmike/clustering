\d .clust

// The following is a second pass implementation of K-means clustering.
// issues with performance were hampering the use of the first pass when compared
// to the python equivalent (sklearn.cluster.KMeans)

/* (m)atrix
/* (n)umber of clusters
/*  k repetitions
/* (b)oolean indicating if kpp(1b) or random initialization(0b) is used

kmeans2:{[m;n;k;b]
 init:$[b;kpp2[m;n];randinit2[m;n]];
 centers:k{value avg each x group mindist2[x;y]}[m]/init;
 mindist2[m;centers]}
kpp2:{[m;n](n-1){y,x iwrand2[1]{x x?min x}each flip{sqrt sum x*x-:y}[flip x]'[y]}[m]/1?m}
 
// utils
mindist2:{{k:@[x;where x=0;:;0n];k?min k}each flip{sqrt sum x*x-:y}[flip x]'[y]}
randinit2:{neg[y]?x}
iwrand2:{[n;w]s binr n?last s:sums w}

// The following is an implementation where the assumption that the initial table is
// manipulated into the form ```value flip``` rather than ```flip value flip```
// the assumption being that for large tables this will be more performant

kmeans3:{[m;n;k;b;f]
 init:$[b;kpp3[m;n];randinitf3[m;n]];
 centers:k{{avg each x@\:y}[x]each value group mindist3[x;y;z]}[m;;f]/init;
 mindist3[m;centers;f]}

kmeanminibatch:{[m;n;k;b;f;p]
 init:$[b;kpp3[m;n];randinitf3[m;n]];
 centers:k{[x;y;z;k]{avg each x@\:y}[x]each value group mindist3[randinit3[x;floor k*count x 0];y;z]}[m;;f;p]/init;
 mindist3[m;centers;f]}

// utils
kpp3:{kpp2[flip x;y]}
mindist3:{{k:@[x;where x=0;:;0n];k?min k}each(,'/)z[x]'[y]}
randinitf3:{flip randinit3[x;y]}
randinit3:{x@\:neg[y]?til count x 0}

// Distance metrics
euclidean:{sqrt sum x*x-:y}
manhattan:{sum abs x-y}
chebyshev:{min abs x-y} 
