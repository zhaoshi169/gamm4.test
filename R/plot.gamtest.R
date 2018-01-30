#'Plot a gamtest Object
#'
#'This function plots the semiparametric estimation of nonlinear curves and surface.
#'
#'@param x A gamtest object.
#'@param test.statistic If TRUE, plot the density of the test statistic under null hypothesis; if FALSE, plot the estimated curves/surfaces.
#'@param test.stat.type must have "test.statistic=TRUE". Default is "test.stat.type=density". If "test.stat.type=hist", plot the histogram of the test statistic under null.
#'@param main The title of the plot.
#'@param n The number of points that are used to draw the curves or surfaces in the plot.
#'@param legend.position the position of legend in the plot: "topright", "topleft", "bottomright", "bottomleft", etc.
#'@param se.est If TRUE, plot the pointwise 95\% confidence intervals of curves; if FALSE, don't plot the pointwise confidence intervals.
#'@param type Only used for ploting surfaces. If "type=contour",then contour plot from vis.gam function in mgcv package; if "type=persp", then plot persp from vis.gam function;
#'if "type=plotly.persp", then plot persp from plotly package.
#'@param data.pts Only used for curve plotting. If TRUE, plot raw data points. If FALSE, only plot estimated curves/surfaces.
#'@param one.frame Only used when "type=plotly.persp". If "one.frame=TRUE", plot surface estimations in one frame; otherwise, in separate frames.
#'@param data.pts.3d If TRUE, plot 3D scatterplot by group using "plotly" package.
#'@param ... Other options from package ``mgcv'' ``vis.gam()'' function.
#'@details This function is to plot a gamtest object. If "test.statistic=TRUE", a density plot of the test statistic under null hypothesis will be generated;
#'if "test.statistic=FALSE", the estimated curves/surfaces for all groups are drawn.
#'@seealso \code{\link[mgcv]{gam}}
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
#'t1 <- gam.grptest(y~s(x,bs="cr"),test=~group,data=data.bind)
#'t1
#'plot(t1)
#'plot(t1,test.statistic=TRUE)
#'
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
#'y2 <- 2*x21^2 + 3*x22^2 + 4*sin(2*pi*x21) + e2
#'
#'n3 <- 150
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
#'tspl <- gam.grptest(y~te(x1,x2),test=~group,data=data.bind,N.boot=200,m=225,parallel=FALSE)
#'tspl$p.value #p-value
#'plot(tspl)
#'plot(tspl,test.statistic = TRUE)
#'plot(tspl,type="plotly.persp")
#'
#'@importFrom mgcv vis.gam
#'@importFrom plotly plot_ly add_surface add_trace %>% add_markers layout
#'@importFrom graphics hist legend matplot points par plot
#'@import stats
#'@export
plot.gamtest <- function(x, test.statistic=FALSE, test.stat.type="density", main="", n=256, legend.position="topright", se.est=FALSE, data.pts=FALSE, type="contour",one.frame=TRUE, data.pts.3d=FALSE,...){
  if (test.statistic) {
    par(mfrow=c(1,1))
    if (test.stat.type=="density"){
      plot(density(x$T.boot, from=0), type = "l", lwd=1.5, main=main, xlab="Test Statistic", ylab="Density")
      text <- paste(" T = ", formatC(x$statistic, digits = 4),"\n","p-value = ", formatC(x$p.value, digits = 4))
      legend(x = legend.position, legend = text)
    }else if(test.stat.type=="hist"){
      hist(x$T.boot, main=main, xlab="Test Statistic")
      text <- paste(" T = ", formatC(x$statistic, digits = 4),"\n","p-value = ", formatC(x$p.value, digits = 4))
      legend(x = legend.position, legend = text)
    }
  } else {
    fit.sub <- x$fit
    gn <- length(fit.sub)
    data.bind <- x$data

    if (x$fcn=="gamm4.grptest"){
      if (ncol(x$data)==4) {
        # if (type %in% c("contour","persp","plotly.persp")) stop(paste(type,"is only used for surface comparisons!"))
        par(mfrow=c(1,1))
        u.min <- max(unlist(lapply(levels(data.bind$group),function(x) min(data.bind$x[data.bind$group==x]))))
        u.max <- min(unlist(lapply(levels(data.bind$group),function(x) max(data.bind$x[data.bind$group==x]))))
        u <- seq(from=u.min, to=u.max, length.out=n)
        fit.sub.u <- matrix(unlist(lapply(fit.sub,function(x) predict(x$gam,data.frame(x=u),se.fit=TRUE))),nrow=n)
        if (se.est==TRUE){
          lower <- fit.sub.u[,seq(1,gn*2,by=2)]-1.96*fit.sub.u[,seq(2,gn*2,by=2)]
          upper <- fit.sub.u[,seq(1,gn*2,by=2)]+1.96*fit.sub.u[,seq(2,gn*2,by=2)]
          matplot(u, cbind(fit.sub.u[,seq(1,gn*2,by=2)],lower,upper), lty=1:x$group, col=1:x$group, type="l", lwd=1.5, xlab="x", ylab="m(x)", xlim = range(data.bind$x), ylim = range(data.bind$y))
          if (data.pts==TRUE) {points(data.bind$x, data.bind$y, col=data.bind$group)}
        }else{
          xl <-x$mydataname[2]
          yl <- x$mydataname[3]
          matplot(u, fit.sub.u[,seq(1,gn*2,by=2)], lty=1:x$group, col=1:x$group, type="l", lwd=1.5, xlab=xl, ylab=yl, xlim = range(data.bind$x), ylim = range(data.bind$y))
          if (data.pts==TRUE) {points(data.bind$x, data.bind$y, col=data.bind$group)}
        }
        text <- paste(x$mydataname[4], ": ", levels(data.bind$group), sep="")
        legend(x = legend.position, legend = text, lty=1:x$group, col=1:x$group, lwd=1.5)
      }else if (ncol(x$data)==5){
        if (type=="contour"){
          par(mfrow=c(2,ifelse(gn>2,2,1)))
          xl <-x$mydataname[2]
          yl <- x$mydataname[3]
          if (data.pts==TRUE) {invisible(lapply(as.numeric(levels(data.bind$group)),function(x) {
            mgcv::vis.gam(fit.sub[[x]]$gam,plot.type="contour",xlab=xl,ylab=yl,zlim=c(min(data.bind$y),max(data.bind$y)))
            points(data.bind$x1[data.bind$group==x],data.bind$x2[data.bind$group==x])
          }))} #suppress returning of list
          else{
            invisible(lapply(as.numeric(levels(data.bind$group)),function(x) mgcv::vis.gam(fit.sub[[x]]$gam,plot.type="contour",xlab=xl,ylab=yl,zlim=c(min(data.bind$y),max(data.bind$y)))))
          }
        }else if(type=="persp"){
          par(mfrow=c(2,ifelse(gn>2,2,1)))
          xl <-x$mydataname[2]
          yl <- x$mydataname[3]
          invisible(lapply(as.numeric(levels(data.bind$group)),function(x) mgcv::vis.gam(fit.sub[[x]]$gam,plot.type="persp",xlab=xl,ylab=yl,zlim=c(min(data.bind$y),max(data.bind$y)),...)))
        }
        else if(type=="plotly.persp"){
          u1.min <- max(unlist(lapply(levels(data.bind$group),function(x) min(data.bind$x1[data.bind$group==x]))))
          u1.max <- min(unlist(lapply(levels(data.bind$group),function(x) max(data.bind$x1[data.bind$group==x]))))
          u1 <- seq(from=u1.min, to=u1.max,length.out=n)
          v1.min <- max(unlist(lapply(levels(data.bind$group),function(x) min(data.bind$x2[data.bind$group==x]))))
          v1.max <- min(unlist(lapply(levels(data.bind$group),function(x) max(data.bind$x2[data.bind$group==x]))))
          v1 <- seq(from=v1.min, to=v1.max,length.out=n)
          u1v1 <- expand.grid(u1,v1)
          fit.sub.uv <- matrix(unlist(lapply(fit.sub,function(x) predict(x$gam,data.frame(x1=u1v1[,1],x2=u1v1[,2])))),nrow=n)
          if (data.pts.3d == TRUE) {plotly::plot_ly(data.bind, x = ~x1, y = ~x2, z = ~y, color = ~group) %>% add_markers() }
          if (one.frame==TRUE){
            if (x$group==2){
              p <- plotly::plot_ly(colors=c("red","blue"), x = u1, y = v1, showscale = FALSE)%>%
              plotly::layout(scene = list(xaxis = list(title = x$mydataname[2]), yaxis = list(title = x$mydataname[3]),
                                  zaxis = list(title = x$mydataname[4])))
              p <- add_surface(p, z = ~fit.sub.uv[,1:n], surfacecolor=matrix(rep(0,n*2),nrow=n,ncol=n),cauto=F,cmax=(x$group-1),cmin=0,showscale=F)
              p <- add_surface(p, z = ~fit.sub.uv[,(1+n):(n*2)], surfacecolor=matrix(rep(1,n*2),nrow=n,ncol=n),cauto=F,cmax=(x$group-1),cmin=0,showscale=F)
              p
            }
            else if (x$group==3){
              p <- plotly::plot_ly(colors=c("red","blue","green"), x = u1, y = v1, showscale = FALSE)%>%
                plotly::layout(scene = list(xaxis = list(title = x$mydataname[2]), yaxis = list(title = x$mydataname[3]),
                                            zaxis = list(title = x$mydataname[4])))
              p <- add_surface(p, z = ~fit.sub.uv[,1:n], surfacecolor=matrix(rep(0,n*2),nrow=n,ncol=n),cauto=F,cmax=(x$group-1),cmin=0,showscale=F)
              p <- add_surface(p, z = ~fit.sub.uv[,(1+n):(n*2)], surfacecolor=matrix(rep(1,n*2),nrow=n,ncol=n),cauto=F,cmax=(x$group-1),cmin=0,showscale=F)
              p <- add_surface(p, z = ~fit.sub.uv[,(1+n*2):(n*3)], surfacecolor=matrix(rep(2,n*2),nrow=n,ncol=n),cauto=F,cmax=(x$group-1),cmin=0,showscale=F)
              p
            }
            else if (x$group==4){
              p <- plotly::plot_ly(colors=c("red","blue","green","grey"), x = u1, y = v1, showscale = FALSE)%>%
                plotly::layout(scene = list(xaxis = list(title = x$mydataname[2]), yaxis = list(title = x$mydataname[3]),
                                            zaxis = list(title = x$mydataname[4])))
              p <- add_surface(p, z = ~fit.sub.uv[,1:n], surfacecolor=matrix(rep(0,n*2),nrow=n,ncol=n),cauto=F,cmax=(x$group-1),cmin=0,showscale=F)
              p <- add_surface(p, z = ~fit.sub.uv[,(1+n):(n*2)], surfacecolor=matrix(rep(1,n*2),nrow=n,ncol=n),cauto=F,cmax=(x$group-1),cmin=0,showscale=F)
              p <- add_surface(p, z = ~fit.sub.uv[,(1+n*2):(n*3)], surfacecolor=matrix(rep(2,n*2),nrow=n,ncol=n),cauto=F,cmax=(x$group-1),cmin=0,showscale=F)
              p <- add_surface(p, z = ~fit.sub.uv[,(1+n*3):(n*4)], surfacecolor=matrix(rep(3,n*2),nrow=n,ncol=n),cauto=F,cmax=(x$group-1),cmin=0,showscale=F)
              p
            }
          }
        }else if (one.frame==FALSE){
          if (x$group==2){
            p <- plotly::plot_ly(x = u1, y = v1, showscale = FALSE)%>%
              plotly::layout(scene = list(xaxis = list(title = x$mydataname[2]), yaxis = list(title = x$mydataname[3]),
                                          zaxis = list(title = x$mydataname[4])))
            p1 <- add_surface(p, z = ~fit.sub.uv[,1:n],  showscale = FALSE)
            p2 <- add_surface(p, z = ~fit.sub.uv[,(1+n):(n*2)], showscale = FALSE)
            p1
            p2
          }
          else if (x$group==3){
            p <- plotly::plot_ly(x = u1, y = v1, showscale = FALSE)%>%
              plotly::layout(scene = list(xaxis = list(title = x$mydataname[2]), yaxis = list(title = x$mydataname[3]),
                                          zaxis = list(title = x$mydataname[4])))
            p1 <- add_trace(p, z = ~fit.sub.uv[,1:n], type='surface', showscale = FALSE)
            p2 <- add_trace(p, z = ~fit.sub.uv[,(1+n):(n*2)], type='surface', showscale = FALSE)
            p3 <- add_trace(p, z = ~fit.sub.uv[,(1+n*2):(n*3)], type='surface', showscale = FALSE)
            p1
            p2
            p3
          }
          else if (x$group==4){
            p <- plotly::plot_ly(x = u1, y = v1, showscale = FALSE)%>%
              plotly::layout(scene = list(xaxis = list(title = x$mydataname[2]), yaxis = list(title = x$mydataname[3]),
                                          zaxis = list(title = x$mydataname[4])))
            p1 <- add_trace(p, z = ~fit.sub.uv[,1:n], type='surface', showscale = FALSE)
            p2 <- add_trace(p, z = ~fit.sub.uv[,(1+n):(n*2)], type='surface', showscale = FALSE)
            p3 <- add_trace(p, z = ~fit.sub.uv[,(1+n*2):(n*3)], type='surface', showscale = FALSE)
            p4 <- add_trace(p, z = ~fit.sub.uv[,(1+n*3):(n*4)], type='surface', showscale = FALSE)
            p1
            p2
            p3
            p4
          }
        }
      }
    } else if (x$fcn=="gam.grptest"){
      if (ncol(x$data)==3) {
   # if (type %in% c("contour","persp","plotly.persp")) stop(paste(type,"is only used for surface comparisons!"))
      par(mfrow=c(1,1))
      u.min <- max(unlist(lapply(levels(data.bind$group),function(x) min(data.bind$x[data.bind$group==x]))))
      u.max <- min(unlist(lapply(levels(data.bind$group),function(x) max(data.bind$x[data.bind$group==x]))))
      u <- seq(from=u.min, to=u.max, length.out=n)
      fit.sub.u <- matrix(unlist(lapply(fit.sub,function(x) predict(x,data.frame(x=u),se.fit=TRUE))),nrow=n)
      if (se.est==TRUE){
        lower <- fit.sub.u[,seq(1,gn*2,by=2)]-1.96*fit.sub.u[,seq(2,gn*2,by=2)]
        upper <- fit.sub.u[,seq(1,gn*2,by=2)]+1.96*fit.sub.u[,seq(2,gn*2,by=2)]
        matplot(u, cbind(fit.sub.u[,seq(1,gn*2,by=2)],lower,upper), lty=1:x$group, col=1:x$group, type="l", lwd=1.5, xlab="x", ylab="m(x)", xlim = range(data.bind$x), ylim = range(data.bind$y))
        if (data.pts==TRUE) {points(data.bind$x, data.bind$y, col=data.bind$group)}
      }else{
        matplot(u, fit.sub.u[,seq(1,gn*2,by=2)], lty=1:x$group, col=1:x$group, type="l", lwd=1.5, xlab=x$mydataname[1], ylab=x$mydataname[2], xlim = range(data.bind$x), ylim = range(data.bind$y))
        if (data.pts==TRUE) {points(data.bind$x, data.bind$y, col=data.bind$group)}
      }
      text <- paste(x$mydataname[3], ": ", levels(data.bind$group), sep="")
      legend(x = legend.position, legend = text, lty=1:x$group, col=1:x$group, lwd=1.5)
      } else if (ncol(x$data)==4){
        if (type=="contour"){
          par(mfrow=c(2,ifelse(gn>2,2,1)))
          xl <- x$mydataname[1]
          yl <- x$mydataname[2]
          if (data.pts==TRUE) {invisible(lapply(as.numeric(levels(data.bind$group)),function(x) {
            mgcv::vis.gam(fit.sub[[x]],plot.type="contour",xlab=xl,ylab=yl,zlim=c(min(data.bind$y),max(data.bind$y)))
            points(data.bind$x1[data.bind$group==x],data.bind$x2[data.bind$group==x])
            }))} #suppress returning of list
          else{
            invisible(lapply(as.numeric(levels(data.bind$group)),function(x) mgcv::vis.gam(fit.sub[[x]],plot.type="contour",xlab=xl,ylab=yl,zlim=c(min(data.bind$y),max(data.bind$y)))))
          }
        }else if(type=="persp"){
          par(mfrow=c(2,ifelse(gn>2,2,1)))
          xl <- x$mydataname[1]
          yl <- x$mydataname[2]
          invisible(lapply(as.numeric(levels(data.bind$group)),function(x) mgcv::vis.gam(fit.sub[[x]],plot.type="persp",xlab=xl,ylab=yl,zlim=c(min(data.bind$y),max(data.bind$y)),...)))
        }
        else if(type=="plotly.persp"){
          u1.min <- max(unlist(lapply(levels(data.bind$group),function(x) min(data.bind$x1[data.bind$group==x]))))
          u1.max <- min(unlist(lapply(levels(data.bind$group),function(x) max(data.bind$x1[data.bind$group==x]))))
          u1 <- seq(from=u1.min, to=u1.max,length.out=n)
          v1.min <- max(unlist(lapply(levels(data.bind$group),function(x) min(data.bind$x2[data.bind$group==x]))))
          v1.max <- min(unlist(lapply(levels(data.bind$group),function(x) max(data.bind$x2[data.bind$group==x]))))
          v1 <- seq(from=v1.min, to=v1.max,length.out=n)
          u1v1 <- expand.grid(u1,v1)
          fit.sub.uv <- matrix(unlist(lapply(fit.sub,function(x) predict(x,data.frame(x1=u1v1[,1],x2=u1v1[,2])))),nrow=n)
          if (one.frame==TRUE){
            if (x$group==2){
              p <- plot_ly(colors=c("red","blue"),x = u1, y = v1, showscale = FALSE)%>%
                layout(scene = list(xaxis = list(title = x$mydataname[1]), yaxis = list(title = x$mydataname[2]),
                                    zaxis = list(title = x$mydataname[3])))
              p <- add_surface(p, z = ~fit.sub.uv[,1:n],surfacecolor=matrix(rep(0,n*2),nrow=n,ncol=n),cauto=F,cmax=(x$group-1),cmin=0,showscale=F)
              p <- add_surface(p, z = ~fit.sub.uv[,(1+n):(n*2)], surfacecolor=matrix(rep(1,n*2),nrow=n,ncol=n),cauto=F,cmax=(x$group-1),cmin=0,showscale=F)
             p
            }
             else if (x$group==3){
               p <- plot_ly(colors=c("red","blue","green"),x = u1, y = v1, showscale = FALSE)%>%
                 layout(scene = list(xaxis = list(title = x$mydataname[1]), yaxis = list(title = x$mydataname[2]),
                                     zaxis = list(title = x$mydataname[3])))
               p <- add_surface(p, z = ~fit.sub.uv[,1:n], surfacecolor=matrix(rep(0,n*2),nrow=n,ncol=n),cauto=F,cmax=(x$group-1),cmin=0,showscale=F)
               p <- add_surface(p, z = ~fit.sub.uv[,(1+n):(n*2)], surfacecolor=matrix(rep(1,n*2),nrow=n,ncol=n),cauto=F,cmax=(x$group-1),cmin=0,showscale=F)
               p <- add_surface(p, z = ~fit.sub.uv[,(1+n*2):(n*3)], surfacecolor=matrix(rep(2,n*2),nrow=n,ncol=n),cauto=F,cmax=(x$group-1),cmin=0,showscale=F)
               p
             }
            else if (x$group==4){
              p <- plot_ly(colors=c("red","blue","green","grey"),x = u1, y = v1, showscale = FALSE)%>%
                layout(scene = list(xaxis = list(title = x$mydataname[1]), yaxis = list(title = x$mydataname[2]),
                                    zaxis = list(title = x$mydataname[3])))
              p <- add_surface(p, z = ~fit.sub.uv[,1:n], surfacecolor=matrix(rep(0,n*2),nrow=n,ncol=n),cauto=F,cmax=(x$group-1),cmin=0,showscale=F)
              p <- add_surface(p, z = ~fit.sub.uv[,(1+n):(n*2)], surfacecolor=matrix(rep(1,n*2),nrow=n,ncol=n),cauto=F,cmax=(x$group-1),cmin=0,showscale=F)
              p <- add_surface(p, z = ~fit.sub.uv[,(1+n*2):(n*3)], surfacecolor=matrix(rep(2,n*2),nrow=n,ncol=n),cauto=F,cmax=(x$group-1),cmin=0,showscale=F)
              p <- add_surface(p, z = ~fit.sub.uv[,(1+n*3):(n*4)], surfacecolor=matrix(rep(3,n*2),nrow=n,ncol=n),cauto=F,cmax=(x$group-1),cmin=0,showscale=F)
              p
            }
         }
          }else if (one.frame==FALSE){
            if (x$group==2){
              p <- plot_ly(x = u1, y = v1, showscale = FALSE)%>%
                layout(scene = list(xaxis = list(title = x$mydataname[1]), yaxis = list(title = x$mydataname[2]),
                                    zaxis = list(title = x$mydataname[3])))
              p1 <- add_surface(p, z = ~fit.sub.uv[,1:n],  showscale = FALSE)%>%
                layout(scene = list(zaxis = list(title = "z")))
              p2 <- add_surface(p, z = ~fit.sub.uv[,(1+n):(n*2)], showscale = FALSE)
              p1
              p2
            }
            else if (x$group==3){
              p <- plot_ly(x = u1, y = v1, showscale = FALSE)%>%
                layout(scene = list(xaxis = list(title = x$mydataname[1]), yaxis = list(title = x$mydataname[2]),
                                    zaxis = list(title = x$mydataname[3])))
              p1 <- add_surface(p, z = ~fit.sub.uv[,1:n], showscale = FALSE)
              p2 <- add_surface(p, z = ~fit.sub.uv[,(1+n):(n*2)], showscale = FALSE)
              p3 <- add_surface(p, z = ~fit.sub.uv[,(1+n*2):(n*3)], showscale = FALSE)
              p1
              p2
              p3
            }
            else if (x$group==4){
              p <- plot_ly(x = u1, y = v1, showscale = FALSE)%>%
                layout(scene = list(xaxis = list(title = x$mydataname[1]), yaxis = list(title = x$mydataname[2]),
                                    zaxis = list(title = x$mydataname[3])))
              p1 <- add_surface(p, z = ~fit.sub.uv[,1:n], showscale = FALSE)
              p2 <- add_surface(p, z = ~fit.sub.uv[,(1+n):(n*2)], showscale = FALSE)
              p3 <- add_surface(p, z = ~fit.sub.uv[,(1+n*2):(n*3)], showscale = FALSE)
              p4 <- add_surface(p, z = ~fit.sub.uv[,(1+n*3):(n*4)], showscale = FALSE)
              p1
              p2
              p3
              p4
            }
          }
        }
      }
    }
  }


