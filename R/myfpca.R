myfpca <- function(MatY, J, t.fit, t.pred, tbas){
  nt <- dim(MatY)[1]
  timebasis <- create.bspline.basis(c(0,1), nbasis=tbas)
  timefdPar <- fdPar(timebasis)
  daytempfd=smooth.basis(t.fit, MatY, timebasis)$fd
  pcafd <- pca.fd(daytempfd, nharm=J, timefdPar)
  harmfd <- pcafd[[1]]
  fdmat.all <- eval.fd(c(t.fit,t.pred), harmfd)
  fdmat=fdmat.all[1:nt,]
  
  etan=apply(fdmat^2,2,sum)
  values=pcafd$values[1:J]*etan/nt
  vectors=fdmat*t(matrix(sqrt(nt/etan),J,nt))
  fd.pre=fdmat.all[(nt+1):(nt+length(t.pred)),]*t(matrix(sqrt(nt/etan),J,length(t.pred)))
  return(list(values=values, vectors=vectors, pred.vec=fd.pre))
}
