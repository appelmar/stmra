#
# Original version (c) 2018 Matthias Katzfuss, see https://github.com/katzfuss-group/MRA_JASA/tree/master/R
# Modifications (c) 2020 Marius Appel

##### smaller help functions

## calculate log determinant from cholesky factor of a matrix
log.det=function(mat.chol) 2*sum(log(diag(mat.chol)))


## return numbered index from tree index (i.e., inverse to indices)
num.ind=function(tree.ind,J){
  if(length(tree.ind)==0) num.index=1
  else {
    m = length(tree.ind)
    J = J[1:m]
    l = cumprod(c(1, rev(J)))
    J = c(1,J)
    nprv  = sum(cumprod(J)[1:m])
    num.index= nprv + 1 + sum(rev(tree.ind-1) * l[1:m])
  }
  return(num.index)
}



## calculate A-B*inv(C)*D' using cholesky of C
woodbury=function(A,B,C.chol,D) {
  A -   if(length(B)<=1 | length(C.chol)==0 | length(D)<=1)  0 else
    t(base::forwardsolve(C.chol,t(B)))%*%base::forwardsolve(C.chol,t(D))
}

## calculate B'*inv(C)*D using cholesky of C
quadform=function(B,C.chol,D) {
  if(length(B)<=1 | length(C.chol)==0 | length(D)<=1)  0 else #matrix(0,nrow=ncol(B),ncol=ncol(D)) else
    t(base::forwardsolve(C.chol,B))%*%forwardsolve(C.chol,D)
}


MRA.fast=function(theta,cov.fun,data,knots,indices,J, pred.locs=NULL) {

  pred=(!is.null(pred.locs)) # do prediction if locs are given (o/w return likelihood)

  ## extract dimensions and other constants
  n.ind=length(indices)
  M=length(indices[[n.ind]])


  indres=vector("list",M+1)
  for(m in 0:M) {
    if (m==0)  indres[[m+1]] = 1
    else {
      indres[[m+1]] = (indres[[m]][length(indres[[m]])]+1):(indres[[m]][length(indres[[m]])]+prod(J[1:m]))
    }
  }

  ## initialize
  Kc.B=vector("list",n.ind)
  if(pred) postmean=postvar=B.tilde=preds=vector("list",n.ind)
  R.prior.chol=vector("list",n.ind)
  A.tilde.cur=w.tilde.cur=A.tilde.prev=w.tilde.prev=vector("list",n.ind)
  loglik.j=numeric(length=n.ind)

  ## going from coarsest to finest resolution
  for(m in 0:M){
    for (ind in indres[[m+1]]) { # loop over all regions at current resolution
      inds=indices[[ind]] # full (j) index

      # create prior quantities
      V=Kc.B[[ind]]=vector("list",m+1)
      for(l in 0:m){
        V[[l+1]]=vector("list",m+1)
        ind.l=if(l==0) 1 else num.ind(inds[1:l],J)
        for(k in l:m){
          ind.k=if(k==0) 1 else num.ind(inds[1:k],J)
          V[[l+1]][[k+1]]= if(l==0) cov.fun(knots[[ind]],knots[[ind.k]],theta) else
            V[[l]][[k+1]]-tp(Kc.B[[ind]][[l]],Kc.B[[ind.k]][[l]])
        }
        if(l<m) {
          Kc.B[[ind]][[l+1]]=sol(R.prior.chol[[ind.l]],t(V[[l+1]][[l+1]]))
        } else {
          R.prior=V[[m+1]][[m+1]]
          R.prior.chol[[ind]]=cholesky(R.prior)
        }
      }

      # begin inference for regions at finest resolution M
      if(m==M) {

        # pre-compute solves
        Sic.B=vector("list",M+1)
        for(l in 0:M) Sic.B[[l+1]]=sol(R.prior.chol[[ind]],V[[l+1]][[l+1]])
        Sic.y=sol(R.prior.chol[[ind]],data[[ind]])

        # inference quantities
        w.tilde.prev[[ind]]=lapply(Sic.B,function(x) tp(x,Sic.y))
        A.tilde.prev[[ind]]=vector("list",m+1)
        for(l in 0:m) {
          A.tilde.prev[[ind]][[l+1]]=vector("list",m+1)
          for(k in l:m)  A.tilde.prev[[ind]][[l+1]][[k+1]] = tp(Sic.B[[l+1]],Sic.B[[k+1]])
        }

        # quantities for prediction or likelihood evaluation
        if(pred) {

          # calculate B.p and L
          V.p=Kc.Bp=V.pp=vector("list",m+1)
          for(l in 0:M){
            V.p[[l+1]]=vector("list",m+1)
            ind.l=if(l==0) 1 else num.ind(inds[1:l],J)
            for(k in l:M){
              ind.k=if(k==0) 1 else num.ind(inds[1:k],J)
              V.p[[l+1]][[k+1]]= if(l==0) cov.fun(pred.locs[[ind]],knots[[ind.k]],theta) else
                V.p[[l]][[k+1]]-tp(Kc.Bp[[l]],Kc.B[[ind.k]][[l]])
            }
            Kc.Bp[[l+1]]=sol(R.prior.chol[[ind.l]],t(V.p[[l+1]][[l+1]])) # Sic.L for l=M
            V.pp[[l+1]]=if(l==0) cov.fun(pred.locs[[ind]],pred.locs[[ind]],theta) else
              V.pp[[l]]-tp(Kc.Bp[[l]],Kc.Bp[[l]])
          }

          # initialize prediction inference
          if (is.matrix(pred.locs[[ind]])) {
            postmean[[ind]]=postvar[[ind]]=matrix(nrow=nrow(pred.locs[[ind]]),ncol=M+1)
          }
          else {
            postmean[[ind]]=postvar[[ind]]=matrix(nrow=length(pred.locs[[ind]]),ncol=M+1)
          }
          postmean[[ind]][,M+1]=tp(Kc.Bp[[M+1]],Sic.y)
          postvar[[ind]][,M+1]=diag(V.pp[[M+1]]-tp(Kc.Bp[[M+1]],Kc.Bp[[M+1]]))
          if(M>0) {
            B.tilde[[ind]]=vector("list",M+1); B.tilde[[ind]][[M+1]]=vector("list",M)
            for(k in 0:(M-1)) B.tilde[[ind]][[M+1]][[k+1]]=V.p[[k+1]][[k+1]]-tp(Kc.Bp[[M+1]],Sic.B[[k+1]])
          }

        } else {
          loglik.j[ind] = log.det(R.prior.chol[[ind]]) + tp(Sic.y,Sic.y)
        }

      }
    }
  }

  rm(V,Kc.B)
  if(pred) rm(V.p,V.pp,Kc.Bp)


  ## posterior inference (from finest to coarsest resolution)
  R.post.chol=Kc.A=Kc.w=w.mm=vector("list",n.ind)

  if(M>0) {
    for(m in seq(M-1,0,by=-1)){
      children=numeric(length=J[m+1])
      A.tilde.cur=w.tilde.cur=vector("list",n.ind)

      for(ind in indres[[m+1]]){

        inds=indices[[ind]] # full (j) index

        # sum up over children tildes
        for(j in 1:J[m+1]) children[j]=num.ind(as.numeric(c(inds,j)),J)
        w=vector("list",m+1)
        A=vector("list",m+1)
        for(l in 0:m){
          w[[l+1]]=Reduce('+',sapply(w.tilde.prev[children],`[`,l+1))
          A[[l+1]]=vector("list",m+1)
          for(k in l:m)  A[[l+1]][[k+1]]=Reduce('+',l.ex(A.tilde.prev[children],l+1,k+1))
        }

        # calculate cholesky of K.inv and save relevant w
        R.post = R.prior.chol[[ind]]%*%t(R.prior.chol[[ind]]) + A[[m+1]][[m+1]]
        R.post.chol[[ind]]=cholesky(R.post)
        w.mm[[ind]]=w[[m+1]]

        # pre-compute the solves required later
        Kc.w[[ind]]=sol(R.post.chol[[ind]],w.mm[[ind]])
        Kc.A[[ind]]=lapply(sapply(A,`[`,m+1),function(x) sol(R.post.chol[[ind]],t(x)))

        if(m>0) {
          # calculate w.tilde and A.tilde
          w.tilde.cur[[ind]]=mapply('-',w,lapply(Kc.A[[ind]],function(x) tp(x,Kc.w[[ind]])),SIMPLIFY=FALSE)
          A.tilde.cur[[ind]]=vector("list",m+1)
          for(l in 0:m) {
            A.tilde.cur[[ind]][[l+1]]=vector("list",m+1)
            for(k in l:m) A.tilde.cur[[ind]][[l+1]][[k+1]] = A[[l+1]][[k+1]] - tp(Kc.A[[ind]][[l+1]],Kc.A[[ind]][[k+1]])
          }
        }

        # likelihood evaluation
        if(!pred) loglik.j[ind] = log.det(R.post.chol[[ind]]) - log.det(R.prior.chol[[ind]]) - tp(Kc.w[[ind]],Kc.w[[ind]])

      }

      A.tilde.prev=A.tilde.cur
      w.tilde.prev=w.tilde.cur

    }
  }


  # spatial prediction
  if(pred) {
    for(ind in indres[[M+1]]){ # only for finest resolution
      if(M>0){
        if(M>1) { for(k in seq(M-1,1,by=-1)){
          ind.k=num.ind(indices[[ind]][1:k],J)
          Kc.Btilde=sol(R.post.chol[[ind.k]],t(B.tilde[[ind]][[k+2]][[k+1]]))
          B.tilde[[ind]][[k+1]]=vector("list",k)
          for(l in seq(k-1,0,by=-1)) B.tilde[[ind]][[k+1]][[l+1]]=B.tilde[[ind]][[k+2]][[l+1]]-tp(Kc.Btilde,Kc.A[[ind.k]][[l+1]])
        } }
        for(k in 0:(M-1)) {
          ind.k=if(k==0) 1 else num.ind(indices[[ind]][1:k],J)
          Kc.Btilde.cur=sol(R.post.chol[[ind.k]],t(B.tilde[[ind]][[k+2]][[k+1]]))
          postmean[[ind]][,k+1]=tp(Kc.Btilde.cur,Kc.w[[ind.k]])
          postvar[[ind]][,k+1]=diag(tp(Kc.Btilde.cur,Kc.Btilde.cur)) # colSums(Kc.Btilde.cur^2)
        }
      }
      preds[[ind]]=cbind(rowSums(postmean[[ind]]),rowSums(postvar[[ind]]))
    }
  }

  return( if(pred) preds else sum(loglik.j))

}
















