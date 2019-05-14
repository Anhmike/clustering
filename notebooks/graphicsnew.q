plt:.p.import`matplotlib.pyplot

/plotting function
/* x = algo e.g.`hc`ward`cure`dbscan
/* y = data
/* z = inputs for the cluster functions

plot:{$[x in`ward`dbscan;plotwdb[x;y;z];plothkc[x;y;z]]}

/plot ward or dbscan
plotwdb:{
 s:.z.t;
 r:$[b:x~`ward;.ml.clust.hc;.ml.clust.dbscan][y;]. z;
 t:.z.t-s;
 $[2<count first y;[fig:plt[`:figure][];ax::fig[`:add_subplot][111;`projection pykw"3d"]];ax::plt];
 {ax[`:scatter]. flip x}each exec pts by clt from r;
 plt[`:title]"df/lf: e2dist/",string[x]," - ",string t;
 plt[`:show][];}

/plot hierarchical, kmeans or cure
plothkc:{
 $[b::2<count first y;fig::plt[`:figure][];
  [subplots::plt[`:subplots]. ud[x;0];fig::subplots[@;0];axarr::subplots[@;1]]];
 fig[`:set_size_inches;18.5;8.5];
 fig[`:subplots_adjust][`hspace pykw .5];
 {[a;d;c;f;i]
  if[b;ax:fig[`:add_subplot][;;i+1;`projection pykw"3d"]. ud[a]0];
  s:.z.t;
  r:ud[a;1][d;c]. f;
  t:.z.t-s;
  if[not b;ax:$[a~`hc;[j:cross[til 3;til 4]i;axarr[@;j 0][@;j 1]];axarr[@;i]]];
  {x[`:scatter]. flip y}[ax]each exec pts by clt from r;
  ax[`:set_title]ud[a;2;f]," - ",string t;
  }[x;y;z[;0]0]'[d;til count d:1_ 'z];
 plt[`:show][];}

/utils dictionary for plothkc
ud:`hc`cure`kmeans!(enlist each(3 4;1 3;1 3)),'
 (.ml.clust.hc;.ml.clust.ccure;.ml.clust.kmeans),'
 ({"df/lf: ",string[x 0],"/",string x 1};{"df/C: ",string[x 0],"/",string[x 3],"b"};{"df: ",string x 2})