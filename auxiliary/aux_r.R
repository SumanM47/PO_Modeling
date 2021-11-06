rm(list=ls())
library(spatstat)

matrix_duplicate <- function(M){
  if(ncol(M)!=2) stop("Function only works for matrices with 2 columns")
  n <- nrow(M)
  initvec <- sapply(1:n,function(i){which(M[i,1]==M[,1] & M[i,2]==M[,2])})
  dupvec <- unique(setdiff(initvec,1:n))
  return(dupvec)
}

## Function to create multitype point process
## Must have a parent and at least one offspring
## Otherwise use rmpoispp() in spatstat
## This function is based on rmpoispp()
##win : window in question; must be of class owin
##lambdapar : parent intensity; positive scalar
##name_par : label for parent process; NULL or character string
##noff : number of offspring species; positive scalar
##mu0vec : average offspring per parent; positive scalar or vector of length noff
##hvec : bandwidth that determines tightness of the cluster; positive scalar or vector of length noff
##names_off : labels for offspring species; NULL or character vector of length noff
##nxtra : number of other species, defaults to 0; non-negative scalar
##lambdaxtra : intensities of the additional processes; positive scalar or vector of length nxtra, missing allowed if nxtra=0
##names_xtra : labels for the other species; NULL or character vector of length nxtra
multiprocppp <- function(win,
                         lambdapar, name_par=NULL,
                         noff,mu0vec,hvec,names_off=NULL,
                         nxtra=0,lambdaxtra,names_xtra=NULL){
  ## Initial Checks
  if(lambdapar <= 0 | length(lambdapar) != 1) stop("lambdapar must be a positive scalar")
  if(noff <= 0) stop("noff must be a positive integer")
  if(length(mu0vec)!=1 & length(mu0vec)!=noff) stop("mu0vec must be a positive scalar or vector of length noff")
  if(length(hvec)!=1 & length(hvec)!=noff) stop("hvec must be a positive scalar or vector of length noff")
  if(sum(mu0vec<=0)!=0) stop("mu0vec must be positive")
  if(sum(hvec<=0)!=0) stop("hvec must be positive")
  if(!is.owin(win)) stop("win must be object of class owin")
  if(nxtra < 0) stop("nxtra must be a non-negative integer")
  if(nxtra > 0 & missing(lambdaxtra)) stop("lambdaxtra is missing. Provide a positive scalar or vector of length nxtra")
  if(nxtra > 0 & (length(lambdaxtra)!=1 & length(lambdaxtra)!=nxtra))  stop("Provide a positive scalar or vector of length nxtra")
  if(nxtra > 0 & (sum(lambdaxtra <=0)!=0)) stop("lambdaxtra must be positive")
  if(!is.null(name_par) & (length(name_par)!=1 | !is.character(name_par))) stop("name_par must a character string of size 1")
  if(!is.null(names_off) & (length(names_off)!=noff & !is.character(names_off))) stop("names_off must be a character vector of size noff")
  if(nxtra > 0 & !is.null(names_xtra) & (length(names_xtra)!=nxtra | !is.character(names_xtra))) stop("names_xtra must be a character vector of size nxtra")
  
  ## Setting labels if not already given
  if(is.null(name_par)) name_par <- "parent"
  if(is.null(names_off)){
    names_off <- paste0("offspring",1:noff,sep="")
  }
  if(nxtra > 0 & is.null(names_xtra)){
    names_xtra <- paste0("xtra",1:nxtra,sep="")
  }
  
  ## Managing inputs
  if(length(mu0vec)==1) mu0vec <- rep(mu0vec,noff)
  if(length(hvec)==1) hvec <- rep(hvec,noff)
  if(nxtra > 0 & length(lambdaxtra)==1) lambdaxtra <- rep(lambdaxtra,nxtra)
  
  require('spatstat')
  
  ## Creating a factor of labels for marking
  types <- c(name_par,names_off,names_xtra)
  factortype <- factor(types, levels = types)
  
  ## First create the parent process
  X <- rpoispp(lambdapar*(1-exp(-max(mu0vec))),win=win)
  parentx <- X$x; parenty <- X$y; np <- X$n
  X <- X %mark% factortype[1]
  
  ## Create the offspring process
  ## Specify parent locations in input for the lambda function
  ## Superimpose on the already created marked point process
  for (i in 1:noff) {
    
    csize <- rpois(np,mu0vec[i])
    numoff <- sum(csize)
    x0 <- rep.int(parentx,csize)
    y0 <- rep.int(parenty,csize)
    
    dd <- matrix(rnorm(2*numoff,0,hvec[i]),ncol=2)
    xy <- xy.coords(dd)
    dx <- xy$x
    dy <- xy$y
    xoff <- x0 + dx
    yoff <- y0 + dy
    
    retain <- inside.owin(xoff,yoff,win)
    
    xoff <- xoff[retain]
    yoff <- yoff[retain]
    
    Y <- ppp(xoff,yoff,window=win)
    
    Y <- Y %mark% factortype[1+i]
    X <- superimpose(X, Y, W = X$window, check = FALSE)
  }
  ## Generate from other processes, if any
  if(nxtra > 0){
    for(i in 1:nxtra){
      Y <- rpoispp(lambdaxtra[i],lamx=NULL,
                   win=win,
                   nsim=1,drop=TRUE,ex=NULL,warnwin=TRUE)
      Y <- Y %mark% factortype[1+noff+i]
      X <- superimpose(X, Y, W = X$window, check = FALSE)
    }
  }
  
  ## Permuting to make it look randomize; not needed
  permu <- sample(X$n)
  X <- X[permu]
  
  return(X)
}


check_within_Window <- function(x,y,w){
  if(w$type == "rectangle"){
    xb <- w$xrange; xlb <- min(xb); xub <- max(xb)
    yb <- w$yrange; ylb <- min(yb); yub <- max(yb)
    out <- as.vector((x <= xub & x >= xlb)*(y <= yub & y >= ylb))
  }
  if(w$type=="polygonal"){
    require(sp)
    polyx <- w$bdry[[1]]$x
    polyy <- w$bdry[[1]]$y
    #dimx <- dim(x)
    ind <- sp::point.in.polygon(x,y,polyx,polyy)
    out <- as.numeric(ind>0)
  }
  return(out)
}

getbignum2 <- function(Y_x,Y_y,C_x,C_y,h){
  out <-  .Call("bignum2func",as.numeric(Y_x),as.numeric(Y_y),as.numeric(C_x),as.numeric(C_y),as.numeric(h))
  return(out)
}

gethpars <- function(Y_x,Y_y,C_x,C_y){
  hsamp <- .Call("gethsamp",as.numeric(Y_x),as.numeric(Y_y),as.numeric(C_x),as.numeric(C_y))
  hbar <- mean(hsamp)
  hsd <- sd(hsamp)
  
  lhsig2 = log(1+(hsd/hbar)^2)
  lhmean = log(hbar) - 0.5*lhsig2
  
  return(list("h_init" = hbar,"lhbar" = lhmean,"lhsig" = sqrt(lhsig2)))
}

subdata <- function(dat, genus){
  dd <- dat
  ind <- which(dat$marks==as.character(genus))
  dd$n <- length(ind)
  dd$x <- dat$x[ind]
  dd$y <- dat$y[ind]
  dd$markformat <- "none"
  dd$marks <- NULL
  
  return(dd)
}


NS_MCMC_new_fixp <- function(obj,parent_genus,offspring_genus,
                             bdtype=1,
                             jitter=FALSE,
                             al=0.01,bl=0.01,
                             am=0.01,bm=0.01,
                             ao=0.01,bo=0.01,
                             hclimp = 0.05,
                             B=100,
                             iters=20000,burn=10000,thin=1,step_int=30,
                             store_res=FALSE,outfile=NULL){
  if(store_res & is.null(outfile)){stop("Please specify intermediary output storage file")}
  require(fields)
  require(spatstat)
  
  ## Create window according to boundary type
  ## bdtype=1 (default) is to use as is; 2 and 3 are convex and concave boundaries
  
  if(bdtype==2){
    cvxhull <- rev(chull(obj$x,obj$y))
    wcvx <- owin(poly=list(x=obj$x[cvxhull],y=obj$y[cvxhull]))
    obj$window <- wcvx
  }
  
  
  
  ## Get Window Info
  W <- obj$window
  xmax <- max(W$xrange)
  xmin <- min(W$xrange)
  
  ymax <- max(W$yrange)
  ymin <- min(W$yrange)
  maxd <- sqrt((xmax-xmin)^2 + (ymax-ymin)^2)
  hclim = hclimp*maxd;
  hsd <- hclim/qnorm(0.995)
  W_area <- area(W)
  
  ##Sorting out input
  
  alltaxa <- unique(obj$marks)
  atn <- length(alltaxa)
  
  not <- length(offspring_genus)
  np <- length(parent_genus)
  if(np!=1 & np!=not){stop("Number of parent taxa must be one or same as the number of offspring taxa")}
  notp <- length(unique(c(parent_genus,offspring_genus)))
  nxt <- atn - notp
  
  up <- unique(as.character(parent_genus))
  upno <- setdiff(up,as.character(offspring_genus))
  nupno <- length(upno)  
  upnoid <- sapply(upno,function(v){min(which(parent_genus==v))})
  
  if(length(parent_genus)==1){parent_genus <- rep(parent_genus,not)}
  
  Cpart_x <- Cpart_y <- list(not)
  Ypart_x <- Ypart_y <- list(not)
  n_Ypart <- rep(0,not); n_Cpart <- rep(0,not)
  
  ## Subsetting the data to parent, offspring and other taxa info
  for(ind in 1:not){
    parent_dat <- subdata(obj,as.character(parent_genus[ind]))
    Cpart_x[[ind]] <- parent_dat$x; Cpart_y[[ind]] <- parent_dat$y; n_Cpart[ind] <- parent_dat$n
  }
  
  
  ## List of offspring subdata
  for(ind in 1:not){
    od <- subdata(obj,as.character(offspring_genus[ind]))
    # offspring_datall[[ind]] <- od
    Ypart_x[[ind]] <- od$x; Ypart_y[[ind]] <- od$y; n_Ypart[ind] <- od$n
  }
  
  
  ## List of Other subdata
  if(nxt > 0){
    xtra_genus <- setdiff(alltaxa,union(parent_genus,offspring_genus))
    Opart_x <- Opart_y <- list(nxt)
    n_Opart <- rep(0,nxt)
    
    for (ind in 1:nxt){
      xd <- subdata(obj,as.character(xtra_genus[ind]))
      # offspring_datall[[ind]] <- od
      Opart_x[[ind]] <- xd$x; Opart_y[[ind]] <- xd$y; n_Opart[ind] <- xd$n
    }
  }
  
  
  hvec <- rep(0,not)
  Bignum1vec <- Bignum2vec <- rep(0,not)
  for(ind in 1:not){
    Y_x <- Ypart_x[[ind]]; Y_y <- Ypart_y[[ind]]; n_Ys <- n_Ypart[ind]
    C_x <- Cpart_x[[ind]]; C_y <- Cpart_y[[ind]]; n_Cs <- n_Cpart[ind]
    hl <- gethpars(Y_x,Y_y,C_x,C_y) 
    hs <- hl$h_init
    if(jitter){
      hs <- jitter(hs)
    }
    hvec[ind] <- hs
    Bigmat_x <- C_x + hs*matrix(rnorm(n_Cs*B),n_Cs,B); Bigmat_y <- C_y + hs*matrix(rnorm(n_Cs*B),n_Cs,B)
    Bignum1vec[ind] <- sum(check_within_Window(Bigmat_x,Bigmat_y,W))/B
    
    Bignum2vec[ind] <- getbignum2(Y_x,Y_y,C_x,C_y,hs)
  }
  
  mu0vec <- n_Ypart[1:not]/Bignum1vec
  
  ll.lh <- -mu0vec*Bignum1vec + Bignum2vec
  
  ## Bookkeping
  v.h <- 0.0075
  v.beta <- 0.5
  
  ## Storage
  keep.lambdaC <- matrix(0,(iters-burn)/thin,nupno)
  keep.mu0 <- matrix(0,(iters-burn)/thin,not)
  keep.h <- matrix(0,(iters-burn)/thin,not)
  if(nxt > 0){
    keep.lambdaO <- matrix(0,(iters-burn)/thin,nxt)
  }
  if(store_res){
    temp.h <- temp.mu0 <- matrix(0,10000,not); temp.lambdaC <- matrix(0,10000,nupno)
    if(nxt > 0){
      temp.lambdaO <- matrix(0,10000,nxt)
    }
    store.h <- store.mu0 <- store.lambdaC <- store.lambdaO <- NULL
  }
  
  pratio <- xtra_ep <- 0
  
  ## GO!
  for(i in 1:iters){
    
    
    
    ## Update lambdaC
    lambdaC <- rgamma(nupno,al + n_Cpart[upnoid], bl + W_area)
    
    if(nxt > 0){
      ## Update lambdaO
      lambdaO <- rgamma(nxt,ao + n_Opart, bo + W_area)
    }
    
    
    ## Update alpha_j
    mu0vec <- rgamma(not,am + n_Ypart, bm + Bignum1vec)
    
    ## Update h_j
    for(j in 1:not){
      Y_x <- Ypart_x[[j]]; Y_y <- Ypart_y[[j]]; n_Ys <- n_Ypart[j]
      C_x <- Cpart_x[[j]]; C_y <- Cpart_y[[j]]; n_Cs <- n_Cpart[j]
      Bignum1s <- Bignum1vec[j]; Bignum2s <- Bignum2vec[j]
      hs <- hvec[j]
      
      ll.lhs <- -mu0vec[j]*Bignum1s + Bignum2s
      ll.lh[j] <- ll.lhs
      
      can_hs <- hs + v.h*rnorm(1)
      
      if(can_hs > 0){
        pratio <- -0.5*((can_hs/hsd)^2 - (hs/hsd)^2)
        Bigmat_x <- C_x + can_hs*matrix(rnorm(n_Cs*B),n_Cs,B); Bigmat_y <- C_y + can_hs*matrix(rnorm(n_Cs*B),n_Cs,B)
        can_Bignum1 <- sum(check_within_Window(Bigmat_x,Bigmat_y,W))/B

        can_Bignum2 <- getbignum2(Y_x,Y_y,C_x,C_y,can_hs)
        
        can_ll.lhs <- -mu0vec[j]*can_Bignum1 + can_Bignum2
        
        a_lhs <- can_ll.lhs - ll.lhs + pratio
        if(log(runif(1)) < a_lhs){
          Bignum1vec[j] <- can_Bignum1
          Bignum2vec[j] <- can_Bignum2
          hvec[j] <- can_hs
          ll.lh[j] <- can_ll.lhs
        }
      }
    }
    
    
    if(i > burn & (i-burn)%%thin == 0){
      index <- (i-burn)/thin
      keep.lambdaC[index,] <- lambdaC
      keep.mu0[index,] <- mu0vec
      keep.h[index,] <- hvec
      if(nxt > 0){
        keep.lambdaO[index,] <- lambdaO
      }
    }
    if(store_res){
      ind2 <- 1 + ((i-1)%%10000)
      temp.lambdaC[ind2,] <- lambdaC
      temp.mu0[ind2,] <- mu0vec
      temp.h[ind2,] <- hvec
      if(nxt > 0){temp.lambdaO[ind2,] <- lambdaO}
      if(i%%10000==0){
        store.lambdaC <- rbind(store.lambdaC,temp.lambdaC)
        store.mu0 <- rbind(store.mu0,temp.mu0)
        store.h <- rbind(store.h,temp.h)
        if(nxt > 0){
          store.lambdaO <- rbind(store.lambdaO,temp.lambdaO)
        }
        temp.h <- temp.mu0 <- matrix(0,10000,not); temp.lambdaC <- matrix(0,10000,nupno)
        if(nxt > 0){
          temp.lambdaO <- matrix(0,10000,nxt)
        }
        lst <- list("alllambdaC" = store.lambdaC,"allmu0"=store.mu0,"allh"=store.h,"ha"=ha,"hb"=hb,"alllambdaO"=store.lambdaO)
        save(lst,file=outfile)
      }
    }
    
  }
  if(nxt==0){keep.lambdaO <- NULL}
  lst <- list("lambdaC"=keep.lambdaC,"mu0"=keep.mu0,"h"=keep.h,"ha"=ha,"hb"=hb,"lambdaO"=keep.lambdaO)
  return(lst)
}