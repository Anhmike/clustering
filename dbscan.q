\d .clust

dbscan.clust:{[c;p;l] /cl idx, minpts, (table;pts idx to search)
 t:l 0;
 s:l 1;
 ncl:{[t;p;s]
  raze{[x;y;z]if[y<=count cl:x[z]`dist;:exec idx from x where idx in cl,valid] /find clusters w/ npts>=p
   }[t;p]each exec idx from t where idx in s,valid
  }[t;p]each s;
 t:update clust:c,valid:0b from t where idx in distinct raze s;
 (t;ncl)
 }

dbscan.calc:{[p;l] /minpts, (table;pts idx to search;cl idx)
 t:l 0;
 s:l 1;
 c:l 2;
 cl:{0<>sum type each x 1}dbscan.clust[c;p]/(t;s);
 nc:first exec idx from t:cl 0 where valid;
 (t;nc;1+c)
 }

dbscan.dbscan:{[d;p;e] /data, minpts, epsilon
 dm:{[x;y;z;n]where y>=@[;n;:;0w]sum x*x-:z}[flip d;e]'[d;k:til count d]; /dist matrix - all dist within eps
 t:([]idx:k;dist:dm;clust:0N;valid:1b); /cluster table
 nn:{0N<>x 1}dbscan.calc[p]/(t;0;0);
 {[x;y]exec idx from x 0 where clust=y}[nn]each distinct nn[0]`clust
 }