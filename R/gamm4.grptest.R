#'Test the equality of nonlinear curves and surface estimations by semiparametric methods with correlated data
#'
#'This function tests the equality of nonlinear curves and surface estimations with correlated data based on L2 distance. The semiparametric estimation
#'uses 'gamm4' package with a compond symmetry correlation structure to adjust correlated observations.
#'The specific model considered here is
#'
#'y_ij= m_i(x_ij) + b_i + e_ij,
#'
#'where m_i(.), are semiparametric smooth functions; b_i are subject-specific random intercept; e_ij are subject-specific errors. The errors e_ij do not have to be independent N(0, sigma^2) errors. The errors can be heteroscedastic, i.e., e_ij = sigma_i(x_ij) * u_ij, where u_ij are independent identically distributed errors with mean 0 and variance 1.
#'
#'We are interested in the problem of testing the equality of the regression curves (when x is one-dimensional) or surfaces (when x is two-dimensional),
#'
#'H_0: m_1(.) = m_2(.) = ...v.s.H_1: otherwise
#'
#'The problem can also be viewed as the test of the equality in the one-sample problem for functional data.
#'@param formula A GAM formula.  This is like the formula for a glm except that smooth terms (s and t2 but not te) can be added to the right hand side of the formula. Note that ids for smooths and fixed smoothing parameters are not supported.
#'@param random An optional formula specifying the random effects structure in lmer style.
#'@param test An indicator of variable for testing nonlinear curves or surface estimations
#'@param data A data frame or list containing the model response variable and covariates required by the formula.
#'@param N.boot the number of bootstrap replicates. This should be a single positive integer.
#'@param m the number of the sampling points for the Monte-Carlo integration.
#'@param parallel Parallel computation of semiparametric estimations with bootstrap samples for getting test statistics under null hypothesis
#'@details A bootstrap algorithm is applied to test the equality of semiparametric curves or surfaces based on L2 distance.
#'@seealso \code{\link[gamm4]{gamm4}} \code{\link{gam.grptest}} \code{\link{plot.gamtest}}
#'@importFrom gamm4 gamm4
#'@importFrom parallel detectCores makeCluster stopCluster
#'@importFrom foreach foreach %dopar% %do%
#'@importFrom doSNOW registerDoSNOW
#'@import mgcv
#'@import stats
#'@importFrom Matrix bdiag
#'@examples
#'\dontrun{
#'#Test the equality of three nonlinear curves
#'m1 <- 120 #number of subjects in group 1
#'m2 <- 100 #number of subjects in group 2
#'m3 <- 110 #number of subjects in group 3
#'n1 <- 3 #number of repeated measurements for each subject in group 1
#'n2 <- 4 #number of repeated measurements for each subject in group 2
#'n3 <- 2 #number of repeated measurements for each subject in group 3
#'sigma1 <- 0.3
#'sigma2 <- 0.2
#'sigma3 <- 0.2
#'sigma.noise1 <- sigma.noise2 <- sigma.noise3 <- 0.1
#'f1 <- function(u) sin(2*pi*u)
#'f2 <- f3 <- function(u) sin(2*pi*u)+u/3
#'N1 <- m1*n1
#'N2 <- m2*n2
#'N3 <- m3*n3
#'
#'x11 <- runif(N1,0,1)
#'b1i <- rnorm(m1,0,sigma1)
#'b1 <- rep(b1i, each=n1)
#'id1 <- rep(1:m1, each=n1)
#'y1 <- f1(x11) + b1 + rnorm(N1, 0, sigma.noise1)
#'
#'x21 <- runif(N2,0,1)
#'b2i <- rnorm(m2,0,sigma2)
#'b2 <- rep(b2i,each=n2)
#'id2 <- rep((m1+1):(m1+m2),each=n2)
#'y2 <- f2(x21) + b2 + rnorm(N2,0,sigma.noise2)
#'
#'x31 <- runif(N3,0,1)
#'b3i <- rnorm(m3,0,sigma3)
#'b3 <- rep(b3i,each=n3)
#'id3 <- rep((m1+m2+1):(m1+m2+m3),each=n3)
#'y3 <- f3(x31) + b3 + rnorm(N3,0,sigma.noise2)
#'
#'dat <- data.frame(rbind(cbind(id1, x11,y1,1), cbind(id2, x21, y2,2), cbind(id3, x31, y3,3)))
#'colnames(dat)=c('id','x', 'y','grp')
#'testout <- gamm4.grptest(formula=y~s(x,k=6,bs="cr"), test=~grp,
#'                         random=~(1|id), data=dat, N.boot=200, m=225, parallel = TRUE)
#'testout
#'plot(testout)
#'
#'dat0 <- data.frame(rbind(cbind(id3, x31, y3, 3), cbind(id2, x21, y2, 2)))
#'colnames(dat0)=c('id', 'x', 'y', 'grp')
#'testout0 <- gamm4.grptest(formula=y~s(x,k=6,bs="cr"), test=~grp,
#'        random=~(1|id), data=dat0, N.boot=200, m=225, parallel= TRUE)
#'testout0$p.value
#'plot(testout0, test.statistic = TRUE)
#'
#'########
#'## Semiparametric test the equality for regression surfaces with longitudinal data
#'## Simulate data sets
#'f1 <- function(u,v) 2*u^2+3*v^2
#'f2 <- function(u,v) 2*u^2+3*v^2+sin(2*pi*u)
#'
#'m1 <- 100 #number of subjects in group 1
#'n1 <- 4 #number of repeated measurements for each subject in group 1
#'m2 <- 120 #number of subjects in group 2
#'n2 <- 3 #number of repeated measurements for each subject in group 2
#'N1 <- m1*n1
#'N2 <- m2*n2
#'sigma1 <- 0.2
#'sigma2 <- 0.15
#'sigma.noise1 <- 0.04
#'sigma.noise2 <- 0.05
#'x11 <- runif(N1,0,1)
#'x12 <- runif(N1,0,1)
#'b1i <- rnorm(m1,0,sigma1)
#'b1 <- rep(b1i,each=n1)
#'id1 <- rep(1:m1,each=n1)
#'y1 <- f1(x11,x12) + b1 + rnorm(N1,0, sigma.noise1)
#'
#'x21 <- runif(N2,0,1)
#'x22 <- runif(N2,0,1)
#'b2i <- rnorm(m2,0,sigma2)
#'b2 <- rep(b2i,each=n2)
#'id2 <- rep((m1+1):(m1+m2),each=n2)
#'y2 <- f2(x21,x22) + b2 + rnorm(N2,0,sigma.noise2)
#'
#'y3 <- f1(x21,x22) + b2 + rnorm(N2,0,sigma.noise2)
#'
#'dat <- data.frame(rbind(cbind(id1, x11, x12,y1,1), cbind(id2, x21, x22, y2,2)))
#'colnames(dat)=c('id','x1','x2', 'y','grp')
#'
#'test.spline1 <- gamm4.grptest(formula=y~t2(x1,x2), test=~grp,
#'                random=~(1|id), data=dat, N.boot=200, m=225, parallel=TRUE)
#'plot(test.spline1)
#'plot(test.spline1, type="plotly.persp")
#'plot(test.spline1, type="plotly.persp", data.pts=TRUE)
#'
#'dat0 <- data.frame(rbind(cbind(id1, x11, x12 , y1, 1), cbind(id2, x21, x22, y3, 2)))
#'colnames(dat0)=c('id','x1','x2', 'y','grp')
#'
#'test.spline0 <- gamm4.grptest(y~t2(x1,x2), test=~grp,
#'                random=~(1|id), data=dat0, N.boot=200, m=225, parallel=TRUE)
#'test.spline0
#'plot(test.spline0, test.statistic = FALSE)
#'plot(test.spline0)
#'plot(test.spline0, type="plotly.persp")
#'
#'########
#'## Data analyses with internal "outchild" dataset
#'
#'data("outchild")
#'outchild1016 <- outchild[(outchild$age<=16 & outchild$age>10),]
#'child.repw <- outchild1016[(outchild1016$RACE==1),]
#'child.reptest1 <- gamm4.grptest(HEIGHT~s(age), random=~(1|SID), 
#'                                test=~SEX, data=child.repw, parallel = TRUE)
#'child.reptest1
#'plot(child.reptest1)
#'plot(child.reptest1,test.statistic = FALSE)
#'
#'child.reptest2 <- gamm4.grptest(WEIGHT~t2(age,HEIGHT), random=~(1|SID), 
#'                               test=~SEX, data = child.repw, parallel = TRUE)
#'plot(child.reptest2,type="plotly.persp")
#'plot(child.reptest2,type="contour")}
#'@export

gamm4.grptest <- function(formula,random,test,data,N.boot=200,m=225,parallel=TRUE){
  gp <-  interpret.gam0(formula) # interpret the formula
  rand.term <- terms.formula(random)
  test.term <- terms.formula(test)
  data.bind <- data.frame(id=data[,as.character(attr(rand.term,"variables")[[2]])[3]],
                          x=data[,gp$pred.names], y=data[,gp$response],
                          group=as.factor(data[,as.character(attr(test.term,"variables")[[2]])]))
  data.bind <- na.omit(data.bind)
  data.bind <- data.bind[order(data.bind$group,data.bind$id),]
  mydataname <- c(as.character(attr(rand.term,"variables")[[2]])[3],gp$pred.names,gp$response,as.character(attr(test.term,"variables")[[2]]))
  
  if(length(gp$pred.names)==1) {
    x = as.matrix(data.bind$x)
    names(data.bind) <- c("id","x", "y", "group")
    knots.s <- gp$smooth.spec[[1]]$bs.dim
  } else{
    x = as.matrix(data.bind[,2:3])
    names(data.bind) <- c("id","x1", "x2", "y", "group")
    knots.s <- c(gp$smooth.spec[[1]]$margin[[1]]$bs.dim,gp$smooth.spec[[1]]$margin[[2]]$bs.dim)
  }
  id <- data.bind$id
  y <- data.bind$y
  group <- as.factor(data.bind$group)
  smooth.class <- class(gp$smooth.spec[[1]])
  
  if (ncol(x) == 1) {
    method <- "Test the equality of curves based on L2 distance"
  } else {
    if(ncol(x) == 2) {
      method <- "Test the equality of surfaces based on L2 distance"
    } else stop("The predictor 'x' should be one or two dimensional!!")
  }
  if (!is.numeric(x)) stop(paste0("argument",gp$pred.names, "must be numeric!"))
  if (!is.numeric(y)) stop(paste0("argument",gp$response, "must be numeric!"))
  
  # if (nrow(x) != length(id) | nrow(x) != length(y) | nrow(x) != length(group))
  #   stop("'id', 'x', 'y' and 'group' have different lengths!")
  # if (is.unsorted(group)) stop("sort by 'group' and 'id'!")
  
  g <- unique(data.bind$group)
  gn <- length(g)
  ny <- length(y)
  
  if (ncol(x)==1) {
    if (smooth.class=="tp.smooth.spec"){ #s()
      if (unique(knots.s)<0){
        fit0 <- gamm4(y ~ s(x), random=~(1|id), data=data.bind)
        fit.sub2 <- lapply(g,function(g)  fit <- gamm4(y ~ s(x),random=~(1|id), data=data.bind[group==g,]))
      }else{
        fit0 <- gamm4(y ~ s(x,k=get("knots.s",parent.frame(n=5))), random=~(1|id),data=data.bind)
        fit.sub2 <- lapply(g,function(g)  fit <- gamm4(y ~ s(x,k=get("knots.s",parent.frame(n=5))), random=~(1|id),data=data.bind[group==g,]))
      }
    }else if(smooth.class=="cr.smooth.spec"){ #s(,bs="cr")
      if (unique(knots.s)<0){
        fit0 <- gamm4(y ~ s(x,bs="cr"), random=~(1|id), data=data.bind)
        fit.sub2 <- lapply(g,function(g)  fit <- gamm4(y ~ s(x,bs="cr"),random=~(1|id),data=data.bind[group==g,]))
      }else{
        # print(parent.env(environment()))
        # print((environment()))
        # print(environment(knots.s))
        fit0 <- gamm4(y ~ s(x,bs="cr",k=get("knots.s",parent.frame(n=5))), random=~(1|id),data=data.bind)
        fit.sub2 <- lapply(g,function(g)  fit <- gamm4(y ~ s(x,bs="cr",k=get("knots.s",parent.frame(n=5))),random=~(1|id),data=data.bind[group==g,]))
      }
    }else if(smooth.class=="ps.smooth.spec"){ #s(,bs="ps")
      if (unique(knots.s)<0){
        fit0 <- gamm4(y ~ s(x,bs="ps"), random=~(1|id), data=data.bind)
        fit.sub2 <- lapply(g,function(g)  fit <- gamm4(y ~ s(x,bs="ps"),random=~(1|id),data=data.bind[group==g,]))
      }else{
        fit0 <- gamm4(y ~ s(x,bs="ps",k=get("knots.s",parent.frame(n=5))), random=~(1|id),data=data.bind)
        fit.sub2 <- lapply(g,function(g)  fit <- gamm4(y ~ s(x,bs="ps",k=get("knots.s",parent.frame(n=5))), random=~(1|id),data=data.bind[group==g,]))
      }
    }
  }else{
    if (smooth.class=="tp.smooth.spec"){ #s()
      if (is.null(unique(knots.s))){
        fit0 <- gamm4(y ~ s(x1,x2),random=~(1|id), data=data.bind)
        fit.sub2 <- lapply(g,function(g,x1,x2)  fit <- gamm4(y ~ s(x1,x2),random=~(1|id),data=data.bind[group==g,]))
      }else{
        fit0 <- gamm4(y ~ s(x1,x2,k=get("knots.s",parent.frame(n=5))),random=~(1|id), data=data.bind)
        fit.sub2 <- lapply(g,function(g,x1,x2)  fit <- gamm4(y ~ s(x1,x2,k=get("knots.s",parent.frame(n=5))),random=~(1|id), data=data.bind[group==g,]))
      }
    }else if(smooth.class=="tensor.smooth.spec"){ #ti()
      if (unique(knots.s)<0){
        fit0 <- gamm4(y ~ te(x1,x2),random=~(1|id), data=data.bind)
        fit.sub2 <- lapply(g,function(g,x1,x2)  fit <- gamm4(y ~ te(x1,x2),random=~(1|id),data=data.bind[group==g,]))
      }else{
        fit0 <- gamm4(y ~ te(x1,x2,k=get("knots.s",parent.frame(n=5))),random=~(1|id), data=data.bind)
        fit.sub2 <- lapply(g,function(g,x1,x2)  fit <- gamm4(y ~ te(x1,x2,k=get("knots.s",parent.frame(n=5))),random=~(1|id), data=data.bind[group==g,]))
      }
    }else if(smooth.class=="t2.smooth.spec"){ #ti()
      if (unique(knots.s)<0){
        fit0 <- gamm4(y ~ t2(x1,x2),random=~(1|id), data=data.bind)
        fit.sub2 <- lapply(g,function(g,x1,x2)  fit <- gamm4(y ~ t2(x1,x2),random=~(1|id),data=data.bind[group==g,]))
      }else{
        fit0 <- gamm4(y ~ t2(x1,x2,k=get("knots.s",parent.frame(n=5))),random=~(1|id), data=data.bind)
        fit.sub2 <- lapply(g,function(g,x1,x2)  fit <- gamm4(y ~ t2(x1,x2,k=get("knots.s",parent.frame(n=5))),random=~(1|id), data=data.bind[group==g,]))
      }
    }
  }
  
  ## pairwise difference
  pwdiff <- function(i, mat) {
    z <- mat[, i-1] - mat[, i:ncol(mat), drop = FALSE]
    colnames(z) <- paste(colnames(mat)[i-1], colnames(z), sep = "-")
    z
  }
  
  # Compute test statistics
  if (ncol(x)==1) {
    u.min <- max(unlist(lapply(g,function(x) min(data.bind$x[data.bind$group==x]))))
    u.max <- min(unlist(lapply(g,function(x) max(data.bind$x[data.bind$group==x]))))
    u <- runif(m, min=u.min, max=u.max)
    fit.sub.x <- matrix(unlist(lapply(fit.sub2,function(x) predict(x$gam,data.frame(x=u)))),nrow=m)
    if (length(g)==2){
      fit.sub.x.diff <- sapply(2:ncol(fit.sub.x), pwdiff, fit.sub.x)
    } else {
      fit.sub.x.diff <- do.call("cbind", sapply(2:ncol(fit.sub.x), pwdiff, fit.sub.x))
    }
    T.spline <- sum(apply(fit.sub.x.diff^2, 2, mean))
  } else{
    u1.min <- max(unlist(lapply(g,function(x) min(data.bind$x1[data.bind$group==x]))))
    u1.max <- min(unlist(lapply(g,function(x) max(data.bind$x1[data.bind$group==x]))))
    u1 <- runif(m, min=u1.min, max=u1.max)
    v1.min <- max(unlist(lapply(g,function(x) min(data.bind$x2[data.bind$group==x]))))
    v1.max <- min(unlist(lapply(g,function(x) max(data.bind$x2[data.bind$group==x]))))
    v1 <- runif(m, min=v1.min, max=v1.max)
    # fit.sub.x <- matrix(unlist(lapply(fit.sub2,function(x) predict(x,data.frame(u=u,v=v)))),nrow=length(dat$y)) #not consider resampling u,v
    fit.sub.x <- matrix(unlist(lapply(fit.sub2,function(x) predict(x$gam,data.frame(x1=u1,x2=v1)))),nrow=m)
    if (length(g)==2){
      fit.sub.x.diff <- sapply(2:ncol(fit.sub.x), pwdiff, fit.sub.x)
    } else {
      fit.sub.x.diff <- do.call("cbind", sapply(2:ncol(fit.sub.x), pwdiff, fit.sub.x))
    }
    T.spline <- sum(apply(fit.sub.x.diff^2, 2, mean))
  }
  
  #############################################################
  #bootstrap under null
  ##  Bootstrap
  if (ncol(x)==1) {
    hat.eta <- data.bind$y-matrix(unlist(lapply(g,function(i)  predict(fit.sub2[[i]]$gam,data.frame(x=data.bind$x[data.bind$group==i]))))) #all hat.eta
  }else{
    hat.eta <- data.bind$y-matrix(unlist(lapply(g,function(i)  predict(fit.sub2[[i]]$gam,data.frame(x1=data.bind$x1[data.bind$group==i],x2=data.bind$x2[data.bind$group==i]))))) #all hat.eta
  }
  # maxrepnum <- lapply(g, function(i) max(table(data.bind$id[data.bind$group==i]))) #maximum number of repeated measurements
  repnum <- lapply(g, function(i) as.data.frame(table(data.bind$id[data.bind$group==i]))$Freq) #number of repeated measurements
  
  bootfunc <- function(x,nboot=1){
    bootres0 <- matrix(rep(NA,nboot*length(x)),ncol=nboot)
    for (i in 1:nboot) bootres0[,i] <- x[sample.int(length(x),nrow(x),replace=TRUE),]
    return(bootres0)
  }
  
  temp_func<-  function(
    var.id, #estimation of random effects
    var.res,  #estimation of random effects
    repnumgrp,
    my.hat.eta)
  {
    # varcov <- matrix(sigmab, ncol = ni, nrow = ni) + diag(sigma2, ni)
    covmat_list<-lapply(repnumgrp,function(s){
      smallSigma <- var.res*diag(s)+var.id*rep(1,s)%*%t(rep(1,s))
    })
    hat.R <- as.matrix(Matrix::bdiag(covmat_list))  #bdiag gets sparse matrix
    hat.L <- t(chol(hat.R)) #cholesky decomposition hat.L%*%t(hat.L)=hat.R
    hat.e <- solve(hat.L) %*% my.hat.eta
    tilde.e <- hat.e-mean(hat.e)  #whitened residual tilde.e
    tilde.e.star <- bootfunc(tilde.e, nboot=N.boot)
    tilde.eta.star <- hat.L%*%tilde.e.star
    return(tilde.eta.star)
  }
  
  tilde.eta.star <- lapply(g,function(i) temp <- temp_func(
    var.id = as.data.frame(summary(fit.sub2[[i]]$mer)[[13]])$vcov[1],
    var.res =  as.data.frame(summary(fit.sub2[[i]]$mer)[[13]])$vcov[3],
    repnumgrp = repnum[[i]],
    my.hat.eta = cbind(data.bind, hat.eta)[group==i,'hat.eta'])) #get tilde.eta.star for each group
  
  if (ncol(x)==1) {
    y.boot <- matrix(rep(predict(fit0$gam,data.frame(x=data.bind$x)),N.boot),ncol=N.boot) + do.call("rbind", tilde.eta.star)#create sample with y*
  } else{
    y.boot <- matrix(rep(predict(fit0$gam,data.frame(x1=data.bind$x1,x2=data.bind$x2)),N.boot),ncol=N.boot) + do.call("rbind", tilde.eta.star)#create sample with y*
  }
  
  T.spline.boot1 <- function(y.boot,knots.s){#, nvar=4
    data1 <- data.frame(id,x,y.boot,group)
    # data1 <- data.frame(matrix(data, ncol=nvar))
    names(data1)=c('id','x', 'y','group')
    cols = c(1, 2, 3)
    data1[,cols] = apply(data1[,cols], 2, function(x) as.numeric(as.character(x)))
    if (smooth.class=="tp.smooth.spec"){ #s()
      if (unique(knots.s)<0){
        fit.sub3 <- lapply(g,function(g) fit <- gamm4(y ~ s(x),random=~(1|id),data=data1[data1$group==g,]))
      }else{
        fit.sub3 <- lapply(g,function(g) fit <- gamm4(y ~ s(x, k=get("knots.s",parent.frame(n=5))), random=~(1|id),data=data1[data1$group==g,]))
      }
    }else if (smooth.class=="cr.smooth.spec"){ #s(,bs="cr")
      if (unique(knots.s)<0){
        fit.sub3 <- lapply(g,function(g) fit <- gamm4(y ~ s(x,bs="cr"),random=~(1|id),data=data1[data1$group==g,]))
      }else{
        fit.sub3 <- lapply(g,function(g) fit <- gamm4(y ~ s(x,bs="cr",k=get("knots.s",parent.frame(n=5))), random=~(1|id),data=data1[data1$group==g,]))
      }
    }else if(smooth.class=="ps.smooth.spec"){ #s(,bs="ps")
      if (unique(knots.s)<0){
        fit.sub3 <- lapply(g,function(g) fit <- gamm4(y ~ s(x,bs="ps"),random=~(1|id),data=data1[data1$group==g,]))
      }else{
        fit.sub3 <- lapply(g,function(g) fit <- gamm4(y ~ s(x,bs="ps",k=get("knots.s",parent.frame(n=5))), random=~(1|id),data=data1[data1$group==g,]))
      }
    }
    fit.sub.x0 <- matrix(unlist(lapply(fit.sub3,function(x) predict(x$gam,data.frame(x=u)))),nrow=m)
    if (length(g)==2){
      fit.sub.x.diff0 <- sapply(2:ncol(fit.sub.x0), pwdiff, fit.sub.x0)
    } else {
      fit.sub.x.diff0 <- do.call("cbind", sapply(2:ncol(fit.sub.x0), pwdiff, fit.sub.x0))
    }
    T.splineout <- sum(apply(fit.sub.x.diff0^2, 2, mean))
    T.spline0 <- ifelse("1" %in% unlist(lapply(g,function(x) fit.sub3[[x]]$mer@optinfo$conv$opt)),NA,T.splineout) #if not converge,NA for T.spline0
    return(T.spline0)
  }
  
  T.spline.boot2 <- function(y.boot,knots.s){#, nvar=5
    data1 <- data.frame(id,x,y.boot,group)
    # data1 <- data.frame(matrix(data, ncol=nvar))
    names(data1)=c('id','x1', 'x2', 'y','group')
    cols = c(1, 2, 3, 4)
    data1[,cols] = apply(data1[,cols], 2, function(x) as.numeric(as.character(x)))
    
    if (smooth.class=="tp.smooth.spec"){ #s()
      if (is.null(unique(knots.s))){
        fit.sub3 <- lapply(g,function(g,x1,x2)  fit <- gamm4(y ~ s(x1,x2),random=~(1|id),data=data1[data1$group==g,]))
      }else{
        fit.sub3 <- lapply(g,function(g,x1,x2)  fit <- gamm4(y ~ s(x1,x2,k=get("knots.s",parent.frame(n=5))),random=~(1|id),data=data1[data1$group==g,]))
      }
    }else if(smooth.class=="t2.smooth.spec"){  #t2()
      if (unique(knots.s)<0){
        fit.sub3 <- lapply(g,function(g,x1,x2)  fit <- gamm4(y ~ t2(x1,x2),random=~(1|id),data=data1[data1$group==g,]))
      }else{
        fit.sub3 <- lapply(g,function(g,x1,x2)  fit <- gamm4(y ~ t2(x1,x2,k=get("knots.s",parent.frame(n=5))),random=~(1|id),data=data1[data1$group==g,]))
      }
    }else if(smooth.class=="tensor.smooth.spec"){  #t2()
      if (unique(knots.s)<0){
        fit.sub3 <- lapply(g,function(g,x1,x2)  fit <- gamm4(y ~ te(x1,x2),random=~(1|id),data=data1[data1$group==g,]))
      }else{
        fit.sub3 <- lapply(g,function(g,x1,x2)  fit <- gamm4(y ~ te(x1,x2,k=get("knots.s",parent.frame(n=5))),random=~(1|id),data=data1[data1$group==g,]))
      }
    }
    
    fit.sub.x0 <- matrix(unlist(lapply(fit.sub3,function(x) predict(x$gam,data.frame(x1=u1,x2=v1)))),nrow=m)
    
    if (length(g)==2){
      fit.sub.x.diff0 <- sapply(2:ncol(fit.sub.x0), pwdiff, fit.sub.x0)
    } else {
      fit.sub.x.diff0 <- do.call("cbind", sapply(2:ncol(fit.sub.x0), pwdiff, fit.sub.x0))
    }
    T.splineout <- sum(apply(fit.sub.x.diff0^2, 2, mean))
    T.spline0 <- ifelse("1" %in% unlist(lapply(g,function(x) fit.sub3[[x]]$mer@optinfo$conv$opt)),NA,T.splineout) #if not converge,NA for T.spline0
    return(T.spline0)
  }
  
  #################
  #start parallel
  if(parallel==FALSE){
    if (ncol(x)==1) {
      T.spline.boot <- apply(y.boot, 2, function(t) T.spline.boot1(t, knots.s))
    } else {
      T.spline.boot <- apply(y.boot, 2, function(t) T.spline.boot2(t, knots.s))
    }
  }else{
    apply.T.spline.boot1<-function(y.boot,T.spline.boot1,N.boot){
      y.boot.col <- i <- NULL
      foreach(y.boot.col=iblkcol(y.boot), .combine="c", .packages=c("foreach","gamm4"), .export=c("x","id","group","g","u","m","pwdiff","smooth.class","knots.s")) %dopar% {
        foreach(i = 1:ncol(y.boot.col)) %do% T.spline.boot1(array(y.boot.col[,i]),knots.s)
      }
    }
    apply.T.spline.boot2<-function(y.boot,T.spline.boot2,N.boot){
      y.boot.col <- i <- NULL
      foreach(y.boot.col=iblkcol(y.boot), .combine="c", .packages=c("foreach","gamm4"), .export=c("x","id","group","g","u1","v1","m","pwdiff","smooth.class","knots.s")) %dopar% {
        foreach(i = 1:ncol(y.boot.col)) %do% T.spline.boot2(array(y.boot.col[,i]),knots.s)
      }
    }
    no.cores <- detectCores()-1
    iblkcol <- function(a, chunks=no.cores) {
      n <- ncol(a)
      i <- 1
      
      nextEl <- function() {
        if (chunks <= 0 || n <= 0) stop('StopIteration')
        m <- ceiling(n / chunks)
        r <- seq(i, length=m)
        i <<- i + m
        n <<- n - m
        chunks <<- chunks - 1
        a[,r, drop=FALSE]
      }
      
      obj <- list(nextElem=nextEl)
      class(obj) <- c('abstractiter', 'iter')
      obj
    }
    myCl <- makeCluster(detectCores()-1)
    doSNOW::registerDoSNOW(myCl)
    
    if (ncol(x)==1) {
      T.spline.boot <- apply.T.spline.boot1(y.boot, T.spline.boot1,N.boot)
    } else {
      T.spline.boot <- apply.T.spline.boot2(y.boot, T.spline.boot2,N.boot)
    }
    stopCluster(myCl)
    #end parallel
    #############
  }
  
  if(length(T.spline.boot[!is.na(unlist(T.spline.boot))])<100) stop("Insufficient number of converged bootstrap samples. Increase number of bootstrap 'N.boot='!")
  
  pval <- (1+sum(T.spline.boot[!is.na(unlist(T.spline.boot))]>T.spline))/(1+length(T.spline.boot[!is.na(T.spline.boot)]))
  output <- list(statistic = T.spline, T.boot = unlist(T.spline.boot), p.value = pval, group = gn, fit = fit.sub2,
                 edf = fit.sub2[[1]]$edf, data = data.bind, mydataname = mydataname,
                 method = method, fcn = "gamm4.grptest")
  class(output) <- "gamtest"
  return(output)
}

