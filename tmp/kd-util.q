\d .clust

/utils
imax:{x?max x}
imin:{x?min x}
edist2:{x wsum x}
/distance and linkage dictionaries
dd:`e2dist`edist`mdist!({x wsum x};{sqrt edist2 x};{sum abs x})
ld:`single`complete`average`centroid!
 ({x y};{x y};{x y};{enlist avg x y}),'
 ({{x y}'[x;imin each x]};{{x y}'[x;imax each x]};{avg x};{raze x}),'
 (min;max;avg;min)

/updated kd-tree with new node inserted
/* t   = kd-tree
/* p   = parent node dimension and index
/* nsd = new splitting dimension
/* d   = data points
/* ii  = initial index
/* ci  = cluster index
/* sd  = splitting dimension

insertt:{[t;p;nsd;d;ii;ci;sd]
 $[not 0b in t`valid;
  t upsert([]idx:1+max t`idx;initi:ii;clust:ii;rep:enlist d;
   dim:(1+p`dim)mod sd;valid:1b;parent:p`idx;dir:nsd 2);
  update idx:1+max t`idx,initi:ii,clust:ci,rep:enlist d,valid:1b,dim:(1+p`dim)mod sd,
   dir:nsd 2,parent:p`idx from t where idx=first exec idx from t where not valid]}

/list of next node to split on, previous node and splitting dimensions
nodedir:{[d;t;nn]
 a:nn 0;
 sd:$[d[first a`dim]>raze[a`rep]first a`dim;1;0];
 i:select from t where dir=sd,parent=first a`idx,valid;
 (i;a;sd)}

/children of node X
branches:{[t;X]raze exec idx from t where parent in X,valid}

/tree node of child with minimum dimension
mindim:{[t;X;child]
 newP:child imin raze({[t;x]first exec rep from t where idx=x,valid
  }[t]each child)[;first exec dim from t where idx=X,valid];
 select from t where idx=newP,valid}


/updated kd-tree with new node n inserted
updatet:{[t;n;X]
 update rep:n`rep,initi:n`initi,closDist:n`closDist,clust:n`clust,
  clustIdx:n`clustIdx,closIdx:n`closIdx from t where idx=X,valid}

/new splitting dimension of node looking at parent
splitdim:{[t;bd;p;df;nn]
 a:select from t where idx=nn,valid;
 nsd:$[(qdim:p d)<rdim:first[a`rep]d:first a`dim;0;1];
 $[bd[0]>=dd[df]rdim-qdim;exec idx from t where parent=nn,valid;
  exec idx from t where parent=nn,dir=nsd,valid],a`parent}

