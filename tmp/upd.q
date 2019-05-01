\d .clust
 
dd:`edist`e2dist`mdist!({sqrt x wsum x};{x wsum x};{sum abs x}) /create complete disctionary of distance functions
ld:`complete`single`avg`ward`centroid!(max;min;avg;{z%(1%y)+1%x};(::)) /create complete dictionary of linkage functions  

imin:{x?min x}
imax:{x?max x}

cldist:{[df;lf;x;y](d i;y[`clt]i:imin d:ld[lf]each dd[df]@'/:raze each x-/:\:/:y`pts)}
/cldist:{[df;lf;x;y;z](d i;x[`clt]i:imin d:ld[lf]each dd[df]@'/:raze each z-/:\:/:exec pts from x:select clt,pts from x where clt<>y)}
/pointdist:{[df;lf;x;y](d i;y[`clt]i:first 1_iasc d:$[lf=`ward;ld[lf][x`n]'[y`n;];]dd[df]each x[`pts]-/:y`pts)}
pointdist:{[df;lf;x;y](d first i j;first ci j:where x[`clt]<>ci:y[`clt]i:2#iasc d:$[lf=`ward;ld[lf][x`n]'[y`n;];]dd[df]each x[`pts]-/:y`pts)}

/
createtab:{
 d:{(d i;i:first 1_iasc d:dd[z]each x-/:y)}[;x;y]each x;
 flip`ind`pts`clt`nni`nnd!(i;x;i:til count x;d[;1];d[;0])
 }
\

createtab:{
 t:flip`ind`pts`clt!(i;x;i:til count x);
 d:pointdist[y;`;;t]each t;
 update nni:d[;1],nnd:d[;0]from t
 }


hc.clust:{[d;c;df;lf]
 if[(lf=`ward) and not df in `edist`e2dist;0N!"Ward's method requires euclidian distances";:()];
 t:createtab[d;df];
 if[lf=`ward;t:@[t;`nnd;%;2]];
 {x<>count exec distinct clt from y}[c]nnc[df;lf]/t
 }


nnc:{[df;lf;x]
 /cd:c,d:(t c:$[count s;s;]0)`nni;
 /cd:c,d:(x c:exec rand clt from x where nnd=min nnd)`nni;
 cd:c,d:(x c:imin x`nnd)`nni;
 /if[d<>s 1;:(d,s;t)];
 x:update clt:min cd from x where clt=max cd;
 updf[lf;df;cd;x]
 }

updw:{[lf;df;cd;t]
 p:sum exec count[i]*first pts by pts from t where clt=min cd;
 t:update pts:count[i]#enlist[p%count[i]]by clt from t where clt=min cd;
 /cl:(%/)sum exec count[i]*first pts,count i by pts from t where clt=min cd;
 /t:update pts:count[i]#enlist cl by clt from t where clt=min cd;
 ct:0!select n:count i,first pts,nn:any nni in cd by clt from t;
 nd:pointdist[df;lf;;ct]each uc:select from ct where nn;
 {[t;x;y]![t;enlist(=;`clt;x);0b;`nnd`nni!y[0],y 1]}/[t;uc`clt;nd]
 }

updc:{[lf;df;cd;t]
 nn:exec pts by clt from t where nni in cd;
 /nc:value count each nn:exec pts by clt from t where nni in cd;
 nd:cldist[df;lf]'[value nn;{[x;y]0!select pts by clt from x where clt<>y}[t]each k:key nn];
 {[t;x;y]![t;enlist(=;`clt;x);0b;`nnd`nni!y[0],y 1]}/[t;k;nd] 
 }

upds:{[lf;df;cd;t]
 t:update nni:min cd from t where nni in cd;
 cl:exec pts by clt from t where clt=min cd;
 cld:cldist[df;lf;first value cl;0!select pts by clt from t where clt<>min cd];
 update nnd:cld[0],nni:cld[1]from t where clt=min cd 
 }


updf:`complete`avg`single`centroid`ward!updc[`complete],updc[`avg],upds[`single],updw[`centroid],updw[`ward]

