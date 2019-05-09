\l clust.q
\d .clust
plt:.p.import`matplotlib.pyplot

plot:{
 subplots::plt[`:subplots]. $[b::2~count first z;3 4;1 3];
 fig::subplots[@;0];
 axarr::subplots[@;1];
 fig[`:set_size_inches;18.5;8.5];
 fig[`:subplots_adjust][`hspace pykw .5];
 {[d;c;f;i]
  s:.z.t;
  r:$[b;hc;cure][d;c;] . f;
  t:.z.t-s;
  j:cross[til 3;til 4]i;
  box:$[b;axarr[@;j 0][@;j 1];axarr[@;i]];
  {x[`:scatter][;]. flip y}[box]each exec pts by clt from r;
  box[`:set_title]"df/",$[b;["lf: ",string[f 0],"/",string f 1];
      ["C: ",string[f 0],"/",string[f 3],"b"]]," - ",string t;
  }[x;y]'[z;til count z];
 plt[`:show][];
 }

plotwdb:{
 s:.z.t;
 r:$[b:x~`ward;hc;dbscan][y;]. z;
 t:.z.t-s;
 {plt[`:scatter][;]. flip x}each exec pts by clt from r;
 plt[`:title]"df/lf: e2dist/",$[b;"ward";"dbscan"]," - ",string t;
 plt[`:show][];}

plotwdb3d:{
 s:.z.t;
 r:$[b:x~`ward;hc;dbscan][y;]. z;
 t:.z.t-s;
 fig:plt[`:figure][];
 ax::fig[`:add_subplot][111;`projection pykw"3d"];
 {ax[`:scatter][;;]. flip x}each exec pts by clt from r;
 plt[`:title]"df/lf: e2dist/,",string[x]," - ",string t;
 plt[`:show][];}