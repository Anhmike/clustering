\d .clust

// The following is an implementation where the assumption that the initial table is
// manipulated into the form ```value flip``` rather than ```flip value flip```
// the assumption being that for large tables this will be more performant

/* dt = data
/* n = number of clusters
/* k = number of iterations
/* i = initialisation type
/* d = distance metric

kmeans:{[dt;n;k;i;d]
 dt:typecast dt;
 init:$[i;kpp[dt;n];randinitf[dt;n]];
 centers:k{{avg each x@\:y}[x]each value group mindist[x;y;z]}[dt;;d]/init;
 mindist[dt;centers;d]}

/* p = percentage of data in each batch
kmeanminibatch:{[dt;n;k;i;d;p]
 dt:typecast dt;
 init:$[i;kpp[dt;n];randinitf[dt;n]];
 centers:k{[x;y;z;k]{avg each x@\:y}[x]each value group mindist[randinit[x;floor k*count x 0];y;z]}[dt;;d;p]/init;
 mindist[dt;centers;d]}

// utils
kpp:{kpp2[flip x;y]}
kpp2:{[m;n](n-1){y,x iwrand[1]{x x?min x}each flip{sqrt sum x*x-:y}[flip x]'[y]}[m]/1?m}
mindist:{{k:@[x;where x=0;:;0n];k?min k}each(,'/)z[x]'[y]}
randinitf:{flip randinit[x;y]}
randinit:{x@\:neg[y]?til count x 0}
typecast:{$[98=type x;value flip x;99=type x;value x;0=type x;x;'`type]}
iwrand:{[n;w]s binr n?last s:sums w}

// Distance metrics
euclidean:{sqrt sum x*x-:y}
manhattan:{sum abs x-y}
chebyshev:{min abs x-y} 


// Testing section

kmeansmerge:{[dt;n;k;i;d;t]
 dt:typecast dt;
 pt1:$[i;kpp[dt;n];randinitf[dt;n]];            / first seeded points
 cntrs:$[(::)~t;
           normcntr[dt;d;pt1;k];
           $[99=type t;
             $[(`mini=key t)0;minicntr[dt;d;pt1;k;value t];'`$"inappropriate dict"];
               '`$"must be dict"];
           '`$"inappropriate input"];
 mindist[dt;cntrs;d]}

normcntr:{[x;y;z;k]k{{avg each x@\:y}[x]each value group mindist[x;y;z]}[x;;y]/z}
minicntr:{[x;y;z;k;p]k{{avg each x@\:y}[x]each value group mindist[randinit[x;floor k*count x 0];y;z]}[x;;y;p]/z}
