######################################
#cross-sectional
######################################
wild.boot <- function(x, nboot=1){
  if (!is.numeric(x)) stop("argument 'x1' must be numeric")
  x <- as.vector(x)
  nx <- length(x)
  if (nboot < 1) stop("'nboot' has to be greater than zero")
  if (nboot==1) {
    a <- rbinom(nx,1,prob=(5+sqrt(5))/10)
    w <- (1-sqrt(5))/2*a+(1+sqrt(5))/2*(1-a)
    x.wb <- w*x
    return(x.wb)
  }
  if (nboot >1 ){
    a0 <- as.matrix(rep(nx, times= nboot))
    a <- apply(a0, 1, function(x) {rbinom(x,1,prob=(5+sqrt(5))/10)})
    w <- (1-sqrt(5))/2*a+(1+sqrt(5))/2*(1-a)
    x.wb <- w*x
    return(x.wb)
  }
}

interpret.gam0 <- function(gf,textra=NULL)
  # interprets a gam formula of the generic form:
  #   y~x0+x1+x3*x4 + s(x5)+ s(x6,x7) ....
  # and returns:
  # 1. a model formula for the parametric part: pf (and pfok indicating whether it has terms)
  # 2. a list of descriptors for the smooths: smooth.spec
  # this is function does the work, and is called by in interpret.gam
{ p.env <- environment(gf) # environment of formula
tf <- terms.formula(gf,specials=c("s","te","ti","t2")) # specials attribute indicates which terms are smooth

terms <- attr(tf,"term.labels") # labels of the model terms 
nt <- length(terms) # how many terms?

if (attr(tf,"response") > 0) {  # start the replacement formulae
  response <- as.character(attr(tf,"variables")[2])
} else { 
  response <- NULL
}
sp <- attr(tf,"specials")$s     # array of indices of smooth terms 
tp <- attr(tf,"specials")$te    # indices of tensor product terms
tip <- attr(tf,"specials")$ti   # indices of tensor product pure interaction terms
t2p <- attr(tf,"specials")$t2   # indices of type 2 tensor product terms
off <- attr(tf,"offset") # location of offset in formula

## have to translate sp, tp, tip, t2p so that they relate to terms,
## rather than elements of the formula...
vtab <- attr(tf,"factors") # cross tabulation of vars to terms
if (length(sp)>0) for (i in 1:length(sp)) {
  ind <- (1:nt)[as.logical(vtab[sp[i],])]
  sp[i] <- ind # the term that smooth relates to
}
if (length(tp)>0) for (i in 1:length(tp)) {
  ind <- (1:nt)[as.logical(vtab[tp[i],])]
  tp[i] <- ind # the term that smooth relates to
} 
if (length(tip)>0) for (i in 1:length(tip)) {
  ind <- (1:nt)[as.logical(vtab[tip[i],])]
  tip[i] <- ind # the term that smooth relates to
} 
if (length(t2p)>0) for (i in 1:length(t2p)) {
  ind <- (1:nt)[as.logical(vtab[t2p[i],])]
  t2p[i] <- ind # the term that smooth relates to
} ## re-referencing is complete

k <- kt <- kti <- kt2 <- ks <- kp <- 1 # counters for terms in the 2 formulae
len.sp <- length(sp)
len.tp <- length(tp)
len.tip <- length(tip)
len.t2p <- length(t2p)
ns <- len.sp + len.tp + len.tip + len.t2p # number of smooths
pav <- av <- rep("",0)
smooth.spec <- list()
mgcvat <- "package:mgcv" %in% search() ## is mgcv in search path?
if (nt) for (i in 1:nt) { # work through all terms
  if (k <= ns&&((ks<=len.sp&&sp[ks]==i)||(kt<=len.tp&&tp[kt]==i)||
                (kti<=len.tip&&tip[kti]==i)||(kt2<=len.t2p&&t2p[kt2]==i))) { # it's a smooth
    ## have to evaluate in the environment of the formula or you can't find variables 
    ## supplied as smooth arguments, e.g. k <- 5;gam(y~s(x,k=k)), fails,
    ## but if you don't specify namespace of mgcv then stuff like 
    ## loadNamespace('mgcv'); k <- 10; mgcv::interpret.gam(y~s(x,k=k)) fails (can't find s)
    ## eval(parse(text=terms[i]),envir=p.env,enclos=loadNamespace('mgcv')) fails??
    ## following may supply namespace of mgcv explicitly if not on search path...
    if (mgcvat) st <- eval(parse(text=terms[i]),envir=p.env) else {
      st <- try(eval(parse(text=terms[i]),envir=p.env),silent=TRUE)
      if (inherits(st,"try-error")) st <- 
          eval(parse(text=terms[i]),enclos=p.env,envir=loadNamespace('mgcv'))
    }
    if (!is.null(textra)) { ## modify the labels on smooths with textra
      pos <- regexpr("(",st$lab,fixed=TRUE)[1]
      st$label <- paste(substr(st$label,start=1,stop=pos-1),textra,
                        substr(st$label,start=pos,stop=nchar(st$label)),sep="")
    }
    smooth.spec[[k]] <- st
    if (ks<=len.sp&&sp[ks]==i) ks <- ks + 1 else # counts s() terms
      if (kt<=len.tp&&tp[kt]==i) kt <- kt + 1 else # counts te() terms
        if (kti<=len.tip&&tip[kti]==i) kti <- kti + 1 else # counts ti() terms
          kt2 <- kt2 + 1                           # counts t2() terms
    k <- k + 1      # counts smooth terms 
  } else {          # parametric
    av[kp] <- terms[i] ## element kp on rhs of parametric
    kp <- kp+1    # counts parametric terms
  }
}    
if (!is.null(off)) { ## deal with offset 
  av[kp] <- as.character(attr(tf,"variables")[1+off])
  kp <- kp+1          
}

pf <- paste(response,"~",paste(av,collapse=" + "))
if (attr(tf,"intercept")==0) {
  pf <- paste(pf,"-1",sep="")
  if (kp>1) pfok <- 1 else pfok <- 0
} else { 
  pfok <- 1;if (kp==1) { 
    pf <- paste(pf,"1"); 
  }
}

fake.formula <- pf

if (length(smooth.spec)>0) 
  for (i in 1:length(smooth.spec)) {
    nt <- length(smooth.spec[[i]]$term)
    ff1 <- paste(smooth.spec[[i]]$term[1:nt],collapse="+")
    fake.formula <- paste(fake.formula,"+",ff1)
    if (smooth.spec[[i]]$by!="NA") {
      fake.formula <- paste(fake.formula,"+",smooth.spec[[i]]$by)
      av <- c(av,smooth.spec[[i]]$term,smooth.spec[[i]]$by)
    } else av <- c(av,smooth.spec[[i]]$term)
  }
fake.formula <- as.formula(fake.formula,p.env)
if (length(av)) {
  pred.formula <- as.formula(paste("~",paste(av,collapse="+")))
  pav <- all.vars(pred.formula) ## trick to strip out 'offset(x)' etc...
  pred.formula <- reformulate(pav) 
} else  pred.formula <- ~1
ret <- list(pf=as.formula(pf,p.env),pfok=pfok,smooth.spec=smooth.spec,
            fake.formula=fake.formula,response=response,fake.names=av,
            pred.names=pav,pred.formula=pred.formula)
class(ret) <- "split.gam.formula"
ret
} ## interpret.gam0

#'Test the equality of nonlinear curves and surface estimations by semiparametric method
#'
#'This function tests the equality of nonlinear curves and surface estimations based on L2 distance. The semiparametric estimation
#'uses 'mgcv' package.
#'The specific model considered here is
#'
#'y_ij= m_i(x_ij) + e_ij,
#'
#'where m_i(.), are semiparametric smooth functions; e_ij are subject-specific errors. The errors e_ij do not have to be independent N(0, sigma^2) errors. The errors can be heteroscedastic, i.e., e_ij = sigma_i(x_ij) * u_ij, where u_ij are independent identically distributed errors with mean 0 and variance 1.
#'
#'We are interested in the problem of testing the equality of the regression curves (when x is one-dimensional) or surfaces (when x is two-dimensional),
#'
#'H_0: m_1(.) = m_2(.) = ... v.s. H_1: otherwise
#'
#'The problem can also be viewed as the test of the equality in the one-sample problem for functional data.
#'
#'@param formula A GAM formula.  This is like the formula for a glm except that smooth terms (s and t2 but not te) can be added to the right hand side of the formula.
#'@param test An indicator of variable for testing nonlinear curves or surface estimations
#'@param data A data frame or list containing the model response variable and covariates required by the formula.
#'@param N.boot the number of bootstrap replicates. This should be a single positive integer.
#'@param m the number of the sampling points for the Monte-Carlo integration.
#'@param parallel Parallel computation of semiparametric estimations with bootstrap samples for getting test statistics under null hypothesis.
#'@details A bootstrap algorithm is applied to test the equality of semiparametric curves or surfaces based on L2 distance.
#'@seealso \code{\link[mgcv]{gam}} \code{\link{gamm4.grptest}} \code{\link{plot.gamtest}} \code{\link{T.L2c}}
#'@import mgcv
#'@importFrom parallel detectCores makeCluster stopCluster
#'@importFrom foreach foreach %dopar% %do%
#'@importFrom doSNOW registerDoSNOW
#'@import stats
#'@examples
#'n1 <- 200
#'x1 <- runif(n1,min=0, max=3)
#'sd1 <- 0.2
#'e1 <- rnorm(n1,sd=sd1)
#'y1 <- sin(2*x1) + cos(2*x1) + e1
#'
#'n2 <- 120
#'x2 <- runif(n2, min=0, max=3)
#'sd2 <- 0.25
#'e2 <- rnorm(n2, sd=sd2)
#'y2 <- sin(2*x2) + cos(2*x2) + x2 + e2
#'
#'data.bind <- rbind(cbind(x1,y1,1), cbind(x2,y2,2))
#'data.bind <- data.frame(data.bind)
#'colnames(data.bind)=c('x','y','group')
#'
#'t1 <- gam.grptest(y~s(x,bs="cr"), test=~group, data=data.bind, parallel=FALSE)
#'t1
#'plot(t1)
#'
#'########
#'\dontrun{
#'## Semiparametric test the equality for regression surfaces
#'## Simulate data sets
#'
#'n1 <- 500
#'x11 <- runif(n1,min=0, max=3)
#'x12 <- runif(n1,min=0, max=3)
#'sd1 <- 0.2
#'e1 <- rnorm(n1,sd=sd1)
#'y1 <- 2*x11^2 + 3*x12^2  + e1
#'
#'n2 <- 420
#'x21 <- runif(n2, min=0, max=3)
#'x22 <- runif(n2, min=0, max=3)
#'sd2 <- 0.25
#'e2 <- rnorm(n2, sd=sd2)
#'y2 <- 2*x21^2 + 3*x22^2 + 6*sin(2*pi*x21) + e2
#'
#'n3 <- 550
#'x31 <- runif(n3,min=0, max=3)
#'x32 <- runif(n3,min=0, max=3)
#'sd3 <- 0.2
#'e3 <- rnorm(n3,sd=sd1)
#'y3 <- 2*x31^2 + 3*x32^2  + e3
#'
#'data.bind <- rbind(cbind(x11, x12 ,y1,1), cbind(x21, x22, y2,2), cbind(x31, x32, y3,3))
#'data.bind <- data.frame(data.bind)
#'colnames(data.bind)=c('x1','x2', 'y','group')
#'
#'tspl <- gam.grptest(y~s(x1,x2), test=~group, data=data.bind, N.boot=200, m=225, parallel=FALSE)
#'tspl$p.value #p-value
#'plot(tspl, test.statistic = TRUE)
#'plot(tspl, type="contour")
#'plot(tspl, type="persp")
#'plot(tspl, type="plotly.persp")
#'plot(tspl, type="plotly.persp",data.pts=TRUE)}
#'
#'########
#'## Data analyses with internal "outchild" dataset
#'
#'data("outchild")
#'child<- outchild[order(outchild$SID,outchild$age),]
#'bs <- aggregate(.~SID, child, FUN=head, 1)
#'
#'childcur <- bs[,c("SEX","WEIGHT","age")]
#'test.grpsex1 <- gam.grptest(WEIGHT~s(age), test=~SEX, data=childcur)
#'test.grpsex1
#'plot(test.grpsex1)
#'plot(test.grpsex1,test.statistic=TRUE)
#'
#'childsurf <- bs[,c("SEX","HEIGHT","WEIGHT","age")]
#'test.grpsex2 <- gam.grptest(WEIGHT~s(HEIGHT,age), test=~SEX, data=childsurf)
#'test.grpsex2
#'plot(test.grpsex2)
#'plot(test.grpsex2, type="plotly.persp")
#'plot(test.grpsex2, type="plotly.persp",data.pts=TRUE)
#'@export

gam.grptest <- function(formula,test,data,N.boot=200,m=225,parallel=FALSE) {
  gp <-  interpret.gam0(formula) # interpret the formula
  test.term <- terms.formula(test)
  data.bind <- data.frame(x=data[,gp$pred.names], y=data[,gp$response],
                          group=as.factor(data[,as.character(attr(test.term,"variables")[[2]])]))
  data.bind <- na.omit(data.bind)
  data.bind <- data.bind[order(data.bind$group),]
  mydataname <- c(gp$pred.names,gp$response,as.character(attr(test.term,"variables")[[2]]))

  if (length(gp$pred.names)==1) {
    names(data.bind) <- c("x", "y", "group")
    data.bind <- data.bind[order(data.bind$group,data.bind$x),]
    knots.s <- gp$smooth.spec[[1]]$bs.dim
    x = as.matrix(data.bind$x)
  } else {
    names(data.bind) <- c("x1", "x2", "y", "group")
    data.bind <- data.bind[order(data.bind$group,data.bind$x1,data.bind$x2),]
    knots.s <- c(gp$smooth.spec[[1]]$margin[[1]]$bs.dim,gp$smooth.spec[[1]]$margin[[2]]$bs.dim)
    x = as.matrix(data.bind[,1:2])
  }

  y <- data.bind$y
  group <- as.factor(data.bind$group)
  smooth.class <- class(gp$smooth.spec[[1]])

  ## CheckValidity
  if (ncol(x) == 1) {
    method <- "Test the equality of curves based on L2 distance"
  } else {
    if(ncol(x) == 2) {
      method <- "Test the equality of surfaces based on L2 distance"
    } else stop("The predictor 'x' should be one or two dimensional!!")
  }

  if (!is.numeric(x)) stop(paste0("argument",gp$pred.names, "must be numeric!"))
  if (!is.numeric(y)) stop(paste0("argument",gp$response, "must be numeric!"))

  # if(nrow(x) != length(y) | nrow(x) != length(group))
  #   stop("'x', 'y' and 'group' have different lengths!")

  g <- unique(group)
  gn <- length(g)
  ny <- length(y)
  if(gn > ny/3) stop("check if there is error in the 'group' variable!")
  if(ny < 3*gn) stop("not enough observations!")


  if (ncol(x)==1) {
    if (smooth.class=="tp.smooth.spec"){ #s()
      if (unique(knots.s)<0){
        fit0 <- gam(y ~ s(x), data=data.bind,select = TRUE)
        fit.sub2 <- lapply(g,function(g)  fit <- gam(y ~ s(x),data=data.bind[group==g,],select = TRUE))
      }else{
        fit0 <- gam(y ~ s(x,k=get("knots.s",parent.frame(n=5))), data=data.bind,select = TRUE)
        fit.sub2 <- lapply(g,function(g)  fit <- gam(y ~ s(x,k=get("knots.s",parent.frame(n=5))), data=data.bind[group==g,],select = TRUE))
      }
    }else if(smooth.class=="cr.smooth.spec"){ #s(,bs="cr")
      if (unique(knots.s)<0){
        fit0 <- gam(y ~ s(x,bs="cr"), data=data.bind,select = TRUE)
        fit.sub2 <- lapply(g,function(g)  fit <- gam(y ~ s(x,bs="cr"), data=data.bind[group==g,],select = TRUE))
      }else{
        fit0 <- gam(y ~ s(x,bs="cr",k=get("knots.s",parent.frame(n=5))), data=data.bind,select = TRUE)
        fit.sub2 <- lapply(g,function(g)  fit <- gam(y ~ s(x,bs="cr",k=get("knots.s",parent.frame(n=5))), data=data.bind[group==g,],select = TRUE))
      }
    }else if(smooth.class=="ps.smooth.spec"){ #s(,bs="ps")
      if (unique(knots.s)<0){
        fit0 <- gam(y ~ s(x,bs="ps"), data=data.bind,select = TRUE)
        fit.sub2 <- lapply(g,function(g)  fit <- gam(y ~ s(x,bs="ps"), data=data.bind[group==g,],select = TRUE))
      }else{
        fit0 <- gam(y ~ s(x,bs="ps",k=get("knots.s",parent.frame(n=5))), data=data.bind,select = TRUE)
        fit.sub2 <- lapply(g,function(g)  fit <- gam(y ~ s(x,bs="ps",k=get("knots.s",parent.frame(n=5))), data=data.bind[group==g,],select = TRUE))
      }
    }
  }else{
    if (smooth.class=="tp.smooth.spec"){ #s()
      if (is.null(unique(knots.s))){
        fit0 <- gam(y ~ s(x1,x2), data=data.bind,select = TRUE)
        fit.sub2 <- lapply(g,function(g,x1,x2)  fit <- gam(y ~ s(x1,x2),data=data.bind[group==g,],select = TRUE))
      }else{
        fit0 <- gam(y ~ s(x1,x2,k=get("knots.s",parent.frame(n=5))), data=data.bind,select = TRUE)
        fit.sub2 <- lapply(g,function(g,x1,x2)  fit <- gam(y ~ s(x1,x2,k=get("knots.s",parent.frame(n=5))), data=data.bind[group==g,],select = TRUE))
      }
    }else if(smooth.class=="tensor.smooth.spec"){  #te()
      if (is.null(unique(knots.s))){
        fit0 <- gam(y ~ te(x1,x2), data=data.bind,select = TRUE)
        fit.sub2 <- lapply(g,function(g,x1,x2)  fit <- gam(y ~ te(x1,x2),data=data.bind[group==g,],select = TRUE))
      }else{
        fit0 <- gam(y ~ te(x1,x2,k=get("knots.s",parent.frame(n=5))), data=data.bind, select = TRUE)
        fit.sub2 <- lapply(g,function(g,x1,x2)  fit <- gam(y ~ te(x1,x2,k=get("knots.s",parent.frame(n=5))), data=data.bind[group==g,],select = TRUE))
      }
    }else if(smooth.class=="t2.smooth.spec"){ #ti()
      if (is.null(unique(knots.s))){
        fit0 <- gam(y ~ t2(x1,x2), data=data.bind,select = TRUE)
        fit.sub2 <- lapply(g,function(g,x1,x2)  fit <- gam(y ~ t2(x1,x2),data=data.bind[group==g,],select = TRUE))
      }else{
        fit0 <- gam(y ~ t2(x1,x2,k=get("knots.s",parent.frame(n=5))), data=data.bind, select = TRUE)
        fit.sub2 <- lapply(g,function(g,x1,x2)  fit <- gam(y ~ t2(x1,x2,k=get("knots.s",parent.frame(n=5))), data=data.bind[group==g,],select = TRUE))
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
    fit.sub.x <- matrix(unlist(lapply(fit.sub2,function(x) predict(x,data.frame(x=u)))),nrow=m)
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
    fit.sub.x <- matrix(unlist(lapply(fit.sub2,function(x) predict(x,data.frame(x1=u1,x2=v1)))),nrow=m)
    if (length(g)==2){
      fit.sub.x.diff <- sapply(2:ncol(fit.sub.x), pwdiff, fit.sub.x)
    } else {
      fit.sub.x.diff <- do.call("cbind", sapply(2:ncol(fit.sub.x), pwdiff, fit.sub.x))
    }
    T.spline <- sum(apply(fit.sub.x.diff^2, 2, mean))
  }

  #############################################################
  #bootstrap under null
  ## Wild Bootstrap
  y.boot <- matrix(rep(fit0$fitted,N.boot),length(fit0$fitted)) + wild.boot(fit0$res, nboot=N.boot) #create sample with y* from null model

  T.spline.boot1 <- function(y.boot,knots.s){
    data1 <- data.frame(x,y.boot,group)
    names(data1)=c('x', 'y','group')
    cols = c(1, 2)
    data1[,cols] = apply(data1[,cols], 2, function(x) as.numeric(as.character(x)))

    if (smooth.class=="tp.smooth.spec"){ #s()
      if (unique(knots.s)<0){
        fit.sub3 <- lapply(g,function(g) fit <- gam(y ~ s(x),data=data1[data1$group==g,],select = TRUE))
      }else{
        fit.sub3 <- lapply(g,function(g) fit <- gam(y ~ s(x,k=get("knots.s",parent.frame(n=5))), data=data1[data1$group==g,],select = TRUE))
      }
    }else if (smooth.class=="cr.smooth.spec"){ #s(,bs="cr")
      if (unique(knots.s)<0){
        fit.sub3 <- lapply(g,function(g) fit <- gam(y ~ s(x,bs="cr"),data=data1[data1$group==g,],select = TRUE))
      }else{
        fit.sub3 <- lapply(g,function(g) fit <- gam(y ~ s(x,bs="cr",k=get("knots.s",parent.frame(n=5))), data=data1[data1$group==g,],select = TRUE))
      }
    }else if (smooth.class=="ps.smooth.spec"){ #s(,bs="cr")
      if (unique(knots.s)<0){
        fit.sub3 <- lapply(g,function(g) fit <- gam(y ~ s(x,bs="ps"),data=data1[data1$group==g,],select = TRUE))
      }else{
        fit.sub3 <- lapply(g,function(g) fit <- gam(y ~ s(x,bs="ps",k=get("knots.s",parent.frame(n=5))), data=data1[data1$group==g,],select = TRUE))
      }
    }
    fit.sub.x0 <- matrix(unlist(lapply(fit.sub3,function(x) predict(x,data.frame(x=u)))),nrow=m)
    if (length(g)==2){
      fit.sub.x.diff0 <- sapply(2:ncol(fit.sub.x0), pwdiff, fit.sub.x0)
    } else {
      fit.sub.x.diff0 <- do.call("cbind", sapply(2:ncol(fit.sub.x0), pwdiff, fit.sub.x0))
    }
    T.spline0 <- sum(apply(fit.sub.x.diff0^2, 2, mean))
    return(T.spline0)
  }

  T.spline.boot2 <- function(y.boot,knots.s){#, nvar=5
    data1 <- data.frame(x,y.boot,group)
    names(data1)=c('x1', 'x2', 'y','group')
    cols = c(1, 2, 3)
    data1[,cols] = apply(data1[,cols], 2, function(x) as.numeric(as.character(x)))

    if (smooth.class=="tp.smooth.spec"){ #s()
      if (is.null(unique(knots.s))){
        fit.sub3 <- lapply(g,function(g,x1,x2)  fit <- gam(y ~ s(x1,x2),data=data1[data1$group==g,],select = TRUE))
      }else{
        fit.sub3 <- lapply(g,function(g,x1,x2)  fit <- gam(y ~ s(x1,x2,k=get("knots.s",parent.frame(n=5))),data=data1[data1$group==g,],select = TRUE))
      }
    }else if(smooth.class=="tensor.smooth.spec"){  #te()
      if (is.null(unique(knots.s))){
        fit.sub3 <- lapply(g,function(g,x1,x2)  fit <- gam(y ~ te(x1,x2),data=data1[data1$group==g,],select = TRUE))
      }else{
        fit.sub3 <- lapply(g,function(g,x1,x2)  fit <- gam(y ~ te(x1,x2,k=get("knots.s",parent.frame(n=5))),data=data1[data1$group==g,],select = TRUE))
      }
    }else if(smooth.class=="t2.smooth.spec"){  #t2()
      if (is.null(unique(knots.s))){
        fit.sub3 <- lapply(g,function(g,x1,x2)  fit <- gam(y ~ t2(x1,x2),data=data1[data1$group==g,],select = TRUE))
      }else{
        fit.sub3 <- lapply(g,function(g,x1,x2)  fit <- gam(y ~ t2(x1,x2,k=get("knots.s",parent.frame(n=5))),data=data1[data1$group==g,],select = TRUE))
      }
    }
    fit.sub.x0 <- matrix(unlist(lapply(fit.sub3,function(x) predict(x,data.frame(x1=u1,x2=v1)))),nrow=m)
    if (length(g)==2){
      fit.sub.x.diff0 <- sapply(2:ncol(fit.sub.x0), pwdiff, fit.sub.x0)
    } else {
      fit.sub.x.diff0 <- do.call("cbind", sapply(2:ncol(fit.sub.x0), pwdiff, fit.sub.x0))
    }
    T.spline0 <- sum(apply(fit.sub.x.diff0^2, 2, mean))
    return(T.spline0)
  }

  if(parallel==FALSE){
    if (ncol(x)==1) {
      T.spline.boot <- apply(y.boot, 2, function(t) T.spline.boot1(t,knots.s))
    } else {
      T.spline.boot <- apply(y.boot, 2, function(t) T.spline.boot2(t,knots.s))
    }
  }else{
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
    apply.T.spline.boot1 <- function(y.boot, T.spline.boot1, N.boot){
      y.boot.col <- i <- NULL
      foreach(y.boot.col = iblkcol(y.boot), .combine="c", .packages=c("foreach","gamm4"), 
              .export=c("x","group","g","u","m","pwdiff","smooth.class","knots.s")) %dopar% {
        foreach(i = 1:ncol(y.boot.col)) %do% T.spline.boot1(array(y.boot.col[,i]),knots.s)
      }
    }
    apply.T.spline.boot2 <- function(y.boot, T.spline.boot2, N.boot){
      y.boot.col <- i <- NULL
      foreach(y.boot.col = iblkcol(y.boot), .combine="c", .packages=c("foreach","gamm4"), 
              .export=c("x","group","g","u1","v1","m","pwdiff","smooth.class","knots.s")) %dopar% {
        foreach(i = 1:ncol(y.boot.col)) %do% T.spline.boot2(array(y.boot.col[,i]),knots.s)
      }
    }
    myCl <- makeCluster(detectCores()-1)
    registerDoSNOW(myCl)

    if (ncol(x)==1) {
      T.spline.boot <- apply.T.spline.boot1(y.boot, T.spline.boot1, N.boot)
    } else {
      T.spline.boot <- apply.T.spline.boot2(y.boot, T.spline.boot2, N.boot)
    }
    stopCluster(myCl)
    #end parallel
    #############
  }

  pval <- (1+sum(T.spline.boot>T.spline))/(1+N.boot)

  output <- list(statistic = T.spline, T.boot = unlist(T.spline.boot), p.value = pval, group = gn, fit = fit.sub2,
                 edf = fit.sub2[[1]]$edf, data = data.bind, mydataname = mydataname,
                 method = method, fcn="gam.grptest")
  class(output) <- "gamtest"
  return(output)
}


#'Test the equality of nonparametric curves or surfaces based on L2 distance
#'
#'This function tests the equality of nonparametric curves and surface estimations based on L2 distance.
#'The specific model considered here is
#'
#'y_ij= m_i(x_ij) + e_ij,
#'
#'where m_i(.), are semiparametric smooth functions; e_ij are subject-specific errors. The errors e_ij do not have to be independent N(0, sigma^2) errors. The errors can be heteroscedastic, i.e., e_ij = sigma_i(x_ij) * u_ij, where u_ij are independent identically distributed errors with mean 0 and variance 1.
#'
#'We are interested in the problem of testing the equality of the regression curves (when x is one-dimensional) or surfaces (when x is two-dimensional),
#'
#'H_0: m_1(.) = m_2(.) = ... v.s. H_1: otherwise
#'
#'The problem can also be viewed as the test of the equality in the one-sample problem for functional data.
#'
#'@param formula A regression formula.  This is like the formula for a lm.
#'@param test An indicator of variable for testing nonparametric curves or surface estimations
#'@param data A data frame or list containing the model response variable and covariates required by the formula.
#'@param N.boot the number of bootstrap replicates. This should be a single positive integer.
#'@param degree the degree of the local polynomials to be used. It can ben 0, 1 or 2.
#'@param criterion the criterion for automatic smoothing parameter selection: ``aicc'' denotes bias-corrected AIC criterion, ``gcv'' denotes generalized cross-validation.
#'@param family if ``gaussian'' fitting is by least-squares, and if ``symmetric'' a re-descending M estimator is used with Tukey's biweight function.
#'@param m the number of the sampling points for the Monte-Carlo integration.
#'@param user.span the user-defined parameter which controls the degree of smoothing.
#'@param ... other options from ``loess'' package.
#'@details A bootstrap algorithm is applied to test the equality of semiparametric curves or surfaces based on L2 distance.
#'@seealso \code{\link{gam.grptest}}
#'@import stats
#'@import mgcv
#'@examples 
#'n1 <- 200
#'x1 <- runif(n1,min=0, max=3)
#'sd1 <- 0.2
#'e1 <- rnorm(n1,sd=sd1)
#'y1 <- sin(2*x1) + cos(2*x1) + e1
#'
#'n2 <- 120
#'x2 <- runif(n2, min=0, max=3)
#'sd2 <- 0.25
#'e2 <- rnorm(n2, sd=sd2)
#'y2 <- sin(2*x2) + cos(2*x2) + x2 + e2
#'
#'dat <- data.frame(rbind(cbind(x1,y1,1), cbind(x2,y2,2)))
#'colnames(dat)=c('x','y','group')
#'
#'t1 <- T.L2c(formula=y~x,test=~group,data=dat)
#'t1$p.value
#'########
#'## Semiparametric test the equality for regression surfaces
#'## Simulate data sets
#'
#'n1 <- 200
#'x11 <- runif(n1,min=0, max=3)
#'x12 <- runif(n1,min=0, max=3)
#'sd1 <- 0.2
#'e1 <- rnorm(n1,sd=sd1)
#'y1 <- 2*x11^2 + 3*x12^2  + e1
#'
#'n2 <- 120
#'x21 <- runif(n2, min=0, max=3)
#'x22 <- runif(n2, min=0, max=3)
#'sd2 <- 0.25
#'e2 <- rnorm(n2, sd=sd2)
#'y2 <- 2*x21^2 + 3*x22^2 + sin(2*pi*x21) + e2
#'
#'n3 <- 150
#'x31 <- runif(n3,min=0, max=3)
#'x32 <- runif(n3,min=0, max=3)
#'sd3 <- 0.2
#'e3 <- rnorm(n3,sd=sd1)
#'y3 <- 2*x31^2 + 3*x32^2  + e3
#'
#'data.bind <- data.frame(rbind(cbind(x11, x12 ,y1,1), cbind(x21, x22, y2,2), cbind(x31, x32, y3,3)))
#'colnames(data.bind)=c('x1','x2', 'y','group')
#'
#'T.L2c(formula=y~x1+x2,test=~group,data=data.bind)
#'@export

T.L2c <- function(formula,test,data,N.boot=200,degree=1, criterion=c("aicc", "gcv"),
                  family = c("gaussian", "symmetric"), m=225, user.span=NULL, ...){
  data.bind <- data[complete.cases(data),]
  gp <-  interpret.gam0(formula)
  x <- data.bind[,gp$pred.names]
  y <- data.bind[,gp$response]
  test.term <- terms.formula(test)
  group <- data.bind[,as.character(attr(test.term,"variables")[[2]])]

  criterion <- match.arg(criterion)
  family <- match.arg(family)
  x <- as.matrix(x)
  if (ncol(x) == 1) {
    method <- "Test the equality of curves based on L2 distance"
  } else {
    if(ncol(x) == 2) {
      method <- "Test the equality of surfaces based on L2 distance"
    } else stop("The predictor 'x' should be one or two dimensional!!")
  }

  ## CheckValidity
  if (!is.numeric(x)) stop("argument 'x' must be numeric!")
  if (!is.numeric(y)) stop("argument 'y' must be numeric!")
  if (any(is.na(x))) stop("'x' contains missing values!")
  if (any(is.na(y))) stop("'y' contains missing values!")
  if (any(is.na(group))) stop("'group' contains missing values!")
  if (!is.null(user.span) && (length(user.span) != 1 || !is.numeric(user.span))) stop("argument 'user.span' must be a numerical number!")
  if(nrow(x) != length(y) | nrow(x) != length(group))
    stop("'x', 'y' and 'group' have different lengths!")

  g <- unique(group)
  gn <- length(g)
  ny <- length(y)
  if(gn > ny/3) stop("check if there is error in the 'group' variable!")
  if(ny < 3*gn) stop("not enough observations!")

  data.bind <- data.frame(x=x, y=y, group=group)
  if (ncol(x) == 1) {
    names(data.bind) <- c("x", "y", "group")
  } else { names(data.bind) <- c("x1", "x2", "y", "group") }

  opt.span <- function(model, criterion=c("aicc", "gcv"), span.range=c(.05, .95)){
    as.crit <- function (x) {
      span <- x$pars$span
      traceL <- x$trace.hat
      sigma2 <- sum(x$residuals^2 ) / (x$n-1)
      aicc <- log(sigma2) + 1 + 2* (2*(traceL+1)) / (x$n-traceL-2)
      gcv <- x$n*sigma2 / (x$n-traceL)^2
      result <- list(span=span, aicc=aicc, gcv=gcv)
      return(result)
    }
    criterion <- match.arg(criterion)
    fn <- function(span) {
      mod <- update(model, span=span)
      as.crit(mod)[[criterion]]
    }
    result <- optimize(fn, span.range)
    return(list(span=result$minimum, criterion=result$objective))
  }

  loc.fit.sub <- function(g, data, dim=c("one", "two"), degree=1, criterion=c("aicc", "gcv"), family = c("gaussian", "symmetric"), user.span=NULL, ...){
    dim <- match.arg(dim)
    opt.span.sub <- function(model, criterion=c("aicc", "gcv"), span.range=c(.05, .95)){
      as.crit <- function (x) {
        span <- x$pars$span
        traceL <- x$trace.hat
        sigma2 <- sum(x$residuals^2 ) / (x$n-1)
        aicc <- log(sigma2) + 1 + 2* (2*(traceL+1)) / (x$n-traceL-2)
        gcv <- x$n*sigma2 / (x$n-traceL)^2
        result <- list(span=span, aicc=aicc, gcv=gcv)
        return(result)
      }
      criterion <- match.arg(criterion)
      fn <- function(span) {
        mod <- update(model, span=span)
        as.crit(mod)[[criterion]]
      }
      result <- optimize(fn, span.range)
      return(list(span=result$minimum, criterion=result$objective))
    }

    subdata <- subset(data, group==g)

    if (dim=="one") {
      if (is.null(user.span)) {
        loc0 <- loess(y ~ x, degree=degree,family = family, data=subdata)
        span1 <- opt.span.sub(loc0, criterion=criterion)$span
      } else {
        span1 <- user.span
      }
      loc1 <- loess(y ~ x, degree=degree, span=span1, family = family, data=subdata,...)
    } else {
      if (is.null(user.span)) {
        loc0 <- loess(y ~ x1 + x2, degree=degree,family = family, data=subdata)
        span1 <- opt.span.sub(loc0, criterion=criterion)$span
      } else {
        span1 <- user.span
      }
      loc1 <- loess(y ~ x1 + x2, degree=degree, span=span1, family = family, data=subdata,...)
    }
    return(loc1)
  }

  ## Fit the curves or surfaces
  if (ncol(x)==1) {
    if (is.null(user.span)) {
      fit0 <- loess(y ~ x, degree=degree, family = family, data=data.bind, ...)
      span1 <- opt.span(fit0, criterion=criterion)$span
      fit <- loess(y ~ x, degree=degree, span=span1, family = family, data=data.bind, ...)
      fit.sub <- lapply(g, loc.fit.sub, data=data.bind, dim="one", degree=degree, criterion=criterion, family = family, ...)
    } else {
      span1 <- user.span
      fit <- loess(y ~ x, degree=degree, span=span1, family = family, data=data.bind, ...)
      fit.sub <- lapply(g, loc.fit.sub, data=data.bind, dim="one", degree=degree, criterion=criterion, family = family, user.span=span1, ...)
    }
  } else {
    if (is.null(user.span)) {
      fit0 <- loess(y ~ x1 + x2, degree=degree,family = family, data.bind, ...)
      span1 <- opt.span(fit0, criterion=criterion)$span
      fit <- loess(y ~ x1 + x2, degree=degree, span=span1, family = family, data=data.bind,...)
      fit.sub <- lapply(g, loc.fit.sub, data=data.bind, dim="two", degree=degree, criterion=criterion, family = family, ...)
    } else {
      span1 <- user.span
      fit <- loess(y ~ x1 + x2, degree=degree, span=span1, family = family, data=data.bind,...)
      fit.sub <- lapply(g, loc.fit.sub, data=data.bind, dim="two", degree=degree, criterion=criterion, family = family, user.span=span1, ...)
    }
  }


  ## Wild Bootstrap
  y.boot <- matrix(rep(fit$fitted,N.boot),fit$n) + wild.boot(fit$res, nboot=N.boot)
  if (ncol(x)==1) {
    x.boot <- matrix(rep(fit$x,N.boot),fit$n)
  } else {x.boot <- matrix(rep(fit$x,N.boot), 2*fit$n)}
  group.boot <- matrix(rep(data.bind$group,N.boot),fit$n)
  data.bind.boot <- rbind(x.boot, y.boot, group.boot)

  ## pairwise difference
  pwdiff <- function(i, mat) {
    z <- mat[, i-1] - mat[, i:ncol(mat), drop = FALSE]
    colnames(z) <- paste(colnames(mat)[i-1], colnames(z), sep = "-")
    z
  }

  ## Compute test statistics
  # find the range to calculate the integration
  if (ncol(x) == 1){
    u.min <- max(unlist(lapply(fit.sub, function(x) min(x$x))))
    u.max <- min(unlist(lapply(fit.sub, function(x) max(x$x))))
    u <- runif(m, min=u.min, max=u.max)
    fit.sub.u <- matrix(unlist(lapply(fit.sub, function(x) predict(x, data.frame(x = u)))),nrow=m)
    if (length(g)==2){
      fit.sub.u.diff <- sapply(2:ncol(fit.sub.u), pwdiff, fit.sub.u)
    } else {
      fit.sub.u.diff <- do.call("cbind", sapply(2:ncol(fit.sub.u), pwdiff, fit.sub.u))
    }
    T.L2 <- sum(apply(fit.sub.u.diff^2, 2, mean))
  } else {
    u1.min <- max(unlist(lapply(fit.sub, function(x) min(x$x[,1]))))
    u1.max <- min(unlist(lapply(fit.sub, function(x) max(x$x[,1]))))
    u1 <- runif(m, min=u1.min, max=u1.max)
    u2.min <- max(unlist(lapply(fit.sub, function(x) min(x$x[,2]))))
    u2.max <- min(unlist(lapply(fit.sub, function(x) max(x$x[,2]))))
    u2 <- runif(m, min=u2.min, max=u2.max)
    fit.sub.u <- matrix(unlist(lapply(fit.sub, function(x) predict(x, data.frame(x1 = u1, x2 = u2)))),nrow=m)
    if (length(g)==2){
      fit.sub.u.diff <- sapply(2:ncol(fit.sub.u), pwdiff, fit.sub.u)
    } else {
      fit.sub.u.diff <- do.call("cbind", sapply(2:ncol(fit.sub.u), pwdiff, fit.sub.u))
    }
    T.L2 <- sum(apply(fit.sub.u.diff^2, 2, mean))
  }

  span0 <- fit$pars$span
  span.sub <- unlist(lapply(fit.sub, function(x)  x$pars$span))
  g.span0 <- cbind(g,span.sub)

  T.L2.boot1 <- function(data, span, g.span, u, nvar=3, degree=1, family = c("gaussian", "symmetric")){
    data1 <- matrix(data, ncol=nvar)
    data1 <- data.frame(data1)
    colnames(data1)=c('x','y','group')

    loc.fit.sub0 <- function(g.span, data, degree=degree, family = family, ...){
      loc1 <- loess(y ~ x, degree=degree, span=g.span[2], subset=(group==g.span[1]), family = family, data=data, ...)
      return(loc1)
    }
    fit.sub <- apply(g.span, 1, loc.fit.sub0, data=data1, degree=degree, family=family, ...)
    fit.sub.u <- matrix(unlist(lapply(fit.sub, function(x) predict(x, data.frame(x = u)))),nrow=length(u))
    if (length(g)==2){
      fit.sub.u.diff <- sapply(2:ncol(fit.sub.u), pwdiff, fit.sub.u)
    } else {
      fit.sub.u.diff <- do.call("cbind", sapply(2:ncol(fit.sub.u), pwdiff, fit.sub.u))
    }
    T.L2 <- sum(apply(fit.sub.u.diff^2, 2, mean))
    return(T.L2)
  }


  T.L2.boot2 <- function(data, span, g.span, u1, u2, nvar=4, degree=1, family = c("gaussian", "symmetric")){
    data1 <- matrix(data, ncol=nvar)
    data1 <- data.frame(data1)
    colnames(data1)=c('x1', 'x2', 'y','group')

    loc.fit.sub0 <- function(g.span, data, degree=degree, family = family, ...){
      loc1 <- loess(y ~ x1 + x2, degree=degree, span=g.span[2], subset=(group==g.span[1]), family = family, data=data, ...)
      return(loc1)
    }
    fit.sub <- apply(g.span, 1, loc.fit.sub0, data=data1, degree=degree, family=family, ...)
    fit.sub.u <- matrix(unlist(lapply(fit.sub, function(x) predict(x, data.frame(x1 = u1, x2 = u2)))), nrow=length(u1))
    if (length(g)==2){
      fit.sub.u.diff <- sapply(2:ncol(fit.sub.u), pwdiff, fit.sub.u)
    } else {
      fit.sub.u.diff <- do.call("cbind", sapply(2:ncol(fit.sub.u), pwdiff, fit.sub.u))
    }

    T.L2 <- sum(apply(fit.sub.u.diff^2, 2, mean))
    return(T.L2)
  }

  if (ncol(x)==1) {
    T.L2.boot <- apply(data.bind.boot, 2, T.L2.boot1, span=span0, g.span=g.span0, u=u, degree=degree, family=family, ...)
  } else { T.L2.boot <- apply(data.bind.boot, 2, T.L2.boot2, span=span0, g.span=g.span0, u1=u1, u2=u2, degree=degree, family=family, ...)}

  pval <- (1+sum(T.L2.boot>T.L2))/(1+N.boot)

  output <- list(statistic=T.L2, T.boot=T.L2.boot, p.value = pval, group=gn, fit=fit.sub, spans=span.sub, degree=degree, criterion=criterion, family = family, data=data.bind, method=method, fcn="T.L2c")
  class(output) <- "gamtest"
  return(output)
}
