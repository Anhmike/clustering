\d .clust

/utils
imax:{x?max x}
imin:{x?min x}

/distance metrics
mdist:{sum abs x}
edist2:{x wsum x}
edist:{sqrt edist2 x}

/linkage dictionary
ld:`single`complete`average`centroid!
 ({x y};{x y};{x y};{enlist avg x y}),'
 ({{x y}'[dd;imin each dd:flip x]};{{x y}'[dd;imax each dd:flip x]};{avg x};{raze x}),'
 (min;max;avg;min)

insertt:{[t;p;nsd;d;l;cl;sd]
 $[not 0b in t`valid;
  t upsert([]idx:1+max t`idx;initi:l;clust:l;rep:enlist d;
   dim:(1+p`dim) mod sd;valid:1b;parent:p`idx;dir:nsd 2);
  update idx:1+max t`idx,initi:l,clust:cl,rep:enlist d,valid:1b,dim:(1+p`dim) mod sd,
   dir:nsd 2,parent:p`idx from t where idx=first exec idx from t where not valid]}

nodedir:{[d;t;nn]a:@[nn;0];
  sd:$[d[first a`dim]>raze[a`rep]first a`dim;1;0];
  i:select from t where dir=sd,parent=first a`idx,valid;
  (i;a;sd)}      /a=prev pt, i=next pt to do split on

branches:{[t;x]raze exec idx from t where parent in x,valid}

mindim:{[t;X;child]
 newP:child imin raze({[t;x]first exec rep from t where idx=x,valid
  }[t]each child)[;first exec dim from t where idx=X,valid];
 select from t where idx=newP,valid}

updatet:{[t;newNode;X]
 update rep:newNode`rep,initi:newNode`initi,closDist:newNode`closDist,
   clust:newNode`clust,clustIdx:newNode`clustIdx, closIdx:newNode`closIdx from t where idx=X,valid}

splitdim:{[t;bd;p;df;nn]
  ll:select from t where idx=nn,valid;
  nsd:$[(qdim:p[dd])<rdim:first[ll`rep]dd:first ll`dim;0;1];
  $[bd[0]>=df[rdim-qdim];exec idx from t where parent=nn,valid;
   exec idx from t where parent=nn,dir=nsd,valid],ll`parent
  }
