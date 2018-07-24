sbde <- function(y, nsamp = 1e3, thin = 10, cens = rep(0,length(y)), wt = rep(1,length(y)), incr = 0.01, par = "prior", nknots = 6, hyper = list(sig = c(.1,.1), lam = c(6,4), kap = c(0.1,0.1,1)), prox.range = c(.2,.95), acpt.target = 0.15, ref.size = 3, blocking = "single", temp = 1, expo = 2, blocks.mu, blocks.S, fix.nu = FALSE, fbase = c("t", "t+", "gpd", "unif"), verbose = TRUE){
	
    fbase.choice <- match(fbase[1], c("t", "t+", "gpd", "unif"))
    if(is.na(fbase.choice)) stop("Only 't', 't+', 'gpd' or 'unif' is allowed for the choice of fbase")
    if(fbase.choice == 1){
        s0 <- function(nu = Inf) return(qt(0.9, df = nu))
        log_f0 <- function(x, nu = Inf) return(dt(x*s0(nu), df = nu, log = TRUE) + log(s0(nu)))
        f0 <- function(x, nu = Inf) return(dt(x * s0(nu), df = nu) * s0(nu))
        F0 <- function(x, nu = Inf) return(pt(x * s0(nu), df = nu))
    } else if (fbase.choice == 2) {
        if(any(y < 0)) stop("'t+' base option may be used only for positive valued data")
        s0 <- function(nu = Inf) return(qt(0.95, df = nu))
        log_f0 <- function(x, nu = Inf) return(log(x > 0) * (log(2.0) + dt(x * s0(nu), df = nu, log = TRUE) + log(s0(nu))))
        f0 <- function(x, nu = Inf) return((x > 0) * 2.0 * dt(x*s0(nu), df = nu)*s0(nu))
        F0 <- function(x, nu = Inf) return((x > 0) * 2.0 * (pt(x*s0(nu), df = nu) - 0.5))
    } else if (fbase.choice == 3) {
        if(any(y < 0)) stop("'gpd' base option may be used only for positive valued data")
        s0 <- function(nu = 1) return(1) ##return(qgpd(0.5, shape = 1/nu))
        log_f0 <- function(x, nu = 1) return(dgpd(x * s0(nu), shape = 1/nu, log = TRUE) + log(s0(nu)))
        f0 <- function(x, nu = 1) return(dgpd(x * s0(nu), shape = 1/nu) * s0(nu))
        F0 <- function(x, nu = 1) return(pgpd(x * s0(nu), shape = 1/nu))
    } else {
        fix.nu <- 1
        s0 <- function(nu = Inf) return(1)
        log_f0 <- function(x, nu = Inf) return(dunif(x, -1,1, log = TRUE))
        f0 <- function(x, nu = Inf) return(dunif(x, -1,1))
        F0 <- function(x, nu = Inf) return(punif(x, -1,1))
    }
    base.bundle <- list(log_f0 = log_f0, f0 = f0, F0 = F0, s0 = s0)
    
    n <- length(y)
    tau.g <- seq(0, 1, incr)
    L <- length(tau.g)
    mid <- which.min(abs(tau.g - 0.5))
	
	tau.kb <- seq(0,1,len = nknots) 
	tau.k <- tau.kb
    delta <- tau.g[2]
    tau.0 <- tau.g[mid]
	
	a.sig <- hyper$sig; if(is.null(a.sig)) a.sig <- c(.1, .1)
	a.lam <- hyper$lam; if(is.null(a.lam)) a.lam <- c(6, 4)
	a.kap <- hyper$kap; if(is.null(a.kap)) a.kap <- c(0.1, 0.1, 1); a.kap <- matrix(a.kap, nrow = 3); nkap <- ncol(a.kap); a.kap[3,] <- log(a.kap[3,])
	hyper.reduced <- c(a.sig, c(a.kap))
	
	prox.grid <- proxFn(max(prox.range), min(prox.range), 0.5)
	ngrid <- length(prox.grid)
	lamsq.grid <- lamFn(prox.grid)^2
	prior.grid <- -diff(pbeta(c(1, (prox.grid[-1] + prox.grid[-ngrid])/2, 0), a.lam[1], a.lam[2]))
	lp.grid <- log(prior.grid)
	
	d.kg <- abs(outer(tau.k, tau.g, "-"))^expo	
	d.kk <- abs(outer(tau.k, tau.k, "-"))^expo
	gridmats <- matrix(NA, nknots*(L + nknots)+2, ngrid)
	K0 <- 0
    t1 <- Sys.time()
	for(i in 1:ngrid){
		K.grid <- exp(-lamsq.grid[i] * d.kg); K.knot <- exp(-lamsq.grid[i] * d.kk);	diag(K.knot) <- 1 + 1e-10	
		R.knot <- chol(K.knot); A.knot <- solve(K.knot, K.grid)
		gridmats[,i] <- c(c(A.knot), c(R.knot), sum(log(diag(R.knot))), lp.grid[i])		
		K0 <- K0 + prior.grid[i] * K.knot
	}
    t2 <- Sys.time()
	
	niter <- nsamp * thin
	dimpars <- c(n, L, mid - 1, nknots, ngrid, ncol(a.kap), niter, thin, nsamp)
	
    if(par == "prior") par <- rep(0, nknots+3)
    if(fix.nu) par[nknots+3] <- nuFn.inv(fix.nu)
    
    if(fbase.choice == 2 | fbase.choice == 3){
        if(substr(blocking, 1, 3) == "std") {
            blocking <- "std1"
        }
    }
	npar <- nknots+3
	if(blocking == "single"){
		blocks <- list(rep(TRUE, npar))
	} else if(blocking == "single2"){
		blocks <- list(rep(TRUE, npar), rep(FALSE, npar))
		blocks[[2]][nknots + 1:3] <- TRUE
	} else if(blocking == "single3"){
		blocks <- list(rep(TRUE, npar), rep(FALSE, npar), rep(FALSE, npar))
        blocks[[2]][1:nknots] <- TRUE
		blocks[[3]][nknots + 1:3] <- TRUE
    } else if(blocking == "std0"){ # same as single2
        blocks <- list(rep(TRUE, npar), rep(FALSE, npar))
        blocks[[2]][nknots + 1:3] <- TRUE
    } else if(blocking == "std1"){ # same as single2
		blocks <- replicate(2, rep(FALSE, npar), simplify = FALSE)
		blocks[[1]] <- rep(TRUE, npar)
		blocks[[2]][nknots + 1:3] <- TRUE
    } else if(blocking == "std2"){ # same as single3
		blocks <- replicate(3, rep(FALSE, npar), simplify = FALSE)
		blocks[[1]] <- rep(TRUE, npar)
		blocks[[2]][nknots + 1] <- TRUE
		blocks[[3]][nknots+1 + 1:2] <- TRUE
	} else if(blocking == "std3"){
		blocks <- replicate(3, rep(FALSE, npar), simplify = FALSE)
		blocks[[1]][1:nknots] <- TRUE
		blocks[[2]][nknots + 1] <- TRUE
		blocks[[3]][nknots + 1 + 1:2] <- TRUE
	} else if(blocking == "std4"){
		blocks <- replicate(3, rep(FALSE, npar), simplify = FALSE)
		blocks[[1]][c(1:nknots, nknots + 1)] <- TRUE
		blocks[[2]][nknots + 1] <- TRUE
		blocks[[3]][nknots + 1 + 1:2] <- TRUE
	} else if(blocking == "std5"){
		blocks <- replicate(4, rep(FALSE, npar), simplify = FALSE)
		blocks[[1]][c(1:nknots, nknots + 1)] <- TRUE
		blocks[[2]][nknots + 1] <- TRUE
		blocks[[3]][nknots+1 + 1:2] <- TRUE
		blocks[[4]][1:npar] <- TRUE
	} else {
		blocks <- replicate(npar, rep(FALSE, npar), simplify = FALSE)
		for(i in 1:npar) blocks[[i]][i] <- TRUE
	}
	
	
	nblocks <- length(blocks)
	if(fix.nu) for(j in 1:nblocks) blocks[[j]][nknots+3] <- FALSE
    if(fbase.choice == 2 | fbase.choice == 3) for(j in 1:nblocks) blocks[[j]][nknots+1] <- FALSE

	
	blocks.ix <- c(unlist(lapply(blocks, which))) - 1
	blocks.size <- sapply(blocks, sum)
	if(missing(blocks.mu)) blocks.mu <- rep(0, sum(blocks.size))
	if(missing(blocks.S)){
		blocks.S <- lapply(blocks.size, function(q) diag(1, q))
		if(substr(blocking, 1, 3) == "std"){
			blocks.S[[1]][1:nknots, 1:nknots] <- K0
			if(as.numeric(substr(blocking, 4,5)) > 1){
                m0 <- quantile(y, prob = tau.0)
                dd <- density(y, n = 1, from = m0, to = m0)
                blocks.S[[2]] <- tau.0 * (1 - tau.0)/(n * dd$y^2)
				blocks.S[[3]] <- matrix(c(1, 0, 0, .1), 2, 2)
			}
			if(as.numeric(substr(blocking, 4,5)) == 5){
				slist <- list(); length(slist) <- 3
				slist[[1]] <- K0
                m0 <- quantile(y, prob = tau.0)
                dd <- density(y, n = 1, from = m0, to = m0)
                slist[[2]] <- tau.0 * (1 - tau.0)/(n * dd$y^2)
				slist[[3]] <- matrix(c(1, 0, 0, .1), 2, 2)
				blocks.S[[4]] <- as.matrix(bdiag(slist))
			}
		}
		blocks.S <- unlist(blocks.S)
	}
		
	imcmc.par <- c(nblocks, ref.size, verbose, max(10, niter/1e4), rep(0, nblocks))
	dmcmc.par <- c(temp, 0.999, rep(acpt.target, nblocks), 2.38 / sqrt(blocks.size))
	
	
	tm.c <- system.time(oo <- .C("SBDE", par = as.double(par), y = as.double(y), cens = as.integer(cens), wt = as.double(wt),
								 hyper = as.double(hyper.reduced), dim = as.integer(dimpars), gridmats = as.double(gridmats),
								 tau.g = as.double(tau.g), muV = as.double(blocks.mu), SV = as.double(blocks.S), blocks = as.integer(blocks.ix), 
								 blocks.size = as.integer(blocks.size), dmcmcpar = as.double(dmcmc.par), 
								 imcmcpar = as.integer(imcmc.par), parsamp = double(nsamp * length(par)), 
								 acptsamp = double(nsamp * nblocks), lpsamp = double(nsamp), fbase.choice = as.integer(fbase.choice)))
	if(verbose) cat("elapsed time:", round(tm.c[3]), "seconds\n")
	
	oo$y <- y; oo$gridmats <- gridmats; oo$prox <- prox.grid; oo$runtime <- tm.c[3]
    oo$base.bundle <- base.bundle
	class(oo) <- "sbde"
	return(oo)
}


update.sbde <- function(object, nadd, append = TRUE, ...){
    niter <- object$dim[7]; thin <- object$dim[8]; nsamp <- object$dim[9]
    if(missing(nadd)) nadd <- nsamp
    par <- object$par; npar <- length(par)
    dimpars <- object$dim
    dimpars[7] <- nadd * thin
    dimpars[9] <- nadd
    nblocks <- object$imcmcpar[1]
    object$imcmcpar[4] <- max(10, nadd * thin/1e4)
    
    tm.c <- system.time(oo <- .C("SBDE", par = as.double(par), y = as.double(object$y), cens = as.integer(object$cens), wt = as.double(object$wt),
        hyper = as.double(object$hyper), dim = as.integer(dimpars), gridmats = as.double(object$gridmats),
        tau.g = as.double(object$tau.g), muV = as.double(object$muV), SV = as.double(object$SV), blocks = as.integer(object$blocks),
        blocks.size = as.integer(object$blocks.size), dmcmcpar = as.double(object$dmcmcpar),
        imcmcpar = as.integer(object$imcmcpar), parsamp = double(nadd * npar),
        acptsamp = double(nadd * nblocks), lpsamp = double(nadd), fbase.choice = as.integer(object$fbase.choice)))
    cat("elapsed time:", round(tm.c[3]), "seconds\n")

    oo$y <- object$y; oo$gridmats <- object$gridmats; oo$prox <- object$prox; oo$runtime <- object$runtime+tm.c[3]
    oo$base.bundle <- object$base.bundle
    if(append){
        oo$dim[7] <- niter + nadd * thin
        oo$dim[9] <- nsamp + nadd
        oo$parsamp <- c(object$parsamp, oo$parsamp)
        oo$acptsamp <- c(object$acptsamp, oo$acptsamp)
        oo$lpsamp <- c(object$lpsamp, oo$lpsamp)
    }
    class(oo) <- "sbde"
    return(oo)
}

coef.sbde <- function(object, burn.perc = 0.5, nmc = 200, reduce = TRUE, ...){
    niter <- object$dim[7]
    nsamp <- object$dim[9]
    pars <- matrix(object$parsamp, ncol = nsamp)
    ss <- unique(round(nsamp * seq(burn.perc, 1, len = nmc + 1)[-1]))
    
    n <- object$dim[1]; p <- 0; L <- object$dim[2]; mid <- object$dim[3] + 1; nknots <- object$dim[4]; ngrid <- object$dim[5]

    parametric.list <- list(gam0 = pars[nknots + 1,ss], sigma = sigFn(pars[nknots + 2,ss], a.sig), nu = nuFn(pars[nknots + 3,ss]))
    gamsignu <- t(sapply(parametric.list, quantile, pr = c(0.5, 0.025, 0.975)))
    dimnames(gamsignu)[[2]] <- c("Estimate", "Lo95%", "Up95%")
    invisible(list(parametric = gamsignu))
}

summary.sbde <- function(object, ntrace = 1000, burn.perc = 0.5, plot.dev = TRUE, more.details = FALSE, ...){
    thin <- object$dim[8]
    nsamp <- object$dim[9]
    pars <- matrix(object$parsamp, ncol = nsamp)
    ss <- unique(pmax(1, round(nsamp * (1:ntrace/ntrace))))
    post.burn <- (ss > nsamp * burn.perc)
    dimpars <- object$dim
    dimpars[7] <- length(ss)
    
    n <- object$dim[1]; p <- 0; ngrid <- object$dim[5]
    sm <- .C("DEV", pars = as.double(pars[,ss]), y = as.double(object$y), cens = as.integer(object$cens), wt = as.double(object$wt), hyper = as.double(object$hyper), dim = as.integer(dimpars), gridmats = as.double(object$gridmats), tau.g = as.double(object$tau.g), devsamp = double(length(ss)), llsamp = double(length(ss)*n), pgsamp = double(length(ss)*ngrid), dist = as.integer(object$fbase.choice))
    deviance <- sm$devsamp
    ll <- matrix(sm$llsamp, ncol = length(ss))
    fit.waic <- waic(ll[,post.burn])
    pg <- matrix(sm$pgsamp, ncol = length(ss))
    prox.samp <- object$prox[apply(pg[1:ngrid,], 2, function(pr) sample(length(pr), 1, prob = pr))]
    
    if(more.details) par(mfrow = c(2,2), mar = c(5,4,3,2)+.1)
    if(plot.dev){
        plot(thin * ss, deviance, ty = "l", xlab = "Markov chain iteration", ylab = "Deviance", bty = "n", main = "Fit trace plot", ...)
        grid(col = "gray")
    }
    
    if(more.details){
        ngrid <- length(object$prox)
        prior.grid <- exp(object$gridmats[nrow(object$gridmats),])
        lam.priorsamp <- lamFn(sample(object$prox, ntrace, replace = TRUE, prob = prior.grid))
        lam.prior.q <- quantile(lam.priorsamp, pr = c(.025, .5, .975))
        lam.samp <- lamFn(prox.samp)
        a <- min(lamFn(object$prox))
        b <- diff(range(lamFn(object$prox))) * 1.2
        plot(thin * ss, lam.samp, ty = "n", ylim = a + c(0, b), bty = "n", ann = FALSE, axes = FALSE)
        axis(1)
        for(i in 1:1){
            abline(h = b * (i-1) + lamFn(object$prox), col = "gray")
            abline(h = b * (i - 1) + lam.prior.q, col = "red", lty = c(2,1,2))
            lines(thin * ss, b * (i-1) + lam.samp, lwd = 1, col = 4)
            if(i %% 2) axis(2, at = b * (i-1) + lamFn(object$prox[c(1,ngrid)]), labels = round(object$prox[c(1,ngrid)],2), las = 1, cex.axis = 0.6)
            mtext(substitute(beta[index], list(index = i - 1)), side = 4, line = 0.5, at = a + b * (i - 1) + 0.4*b, las = 1)
        }
        title(xlab = "Markov chain iteration", ylab = "Proxmity posterior", main = "Mixing over GP scaling")
        
        theta <- as.mcmc(t(matrix(object$parsamp, ncol = nsamp)[,ss[post.burn]]))
        gg <- geweke.diag(theta, .1, .5)
        zvals <- gg$z
        
        pp <- 2 * (1 - pnorm(abs(zvals)))
        plot(sort(pp), ylab = "Geweke p-values", xlab = "Parameter index (reordered)", main = "Convergence diagnosis", ty = "h", col = 4, ylim = c(0, 0.3), lwd = 2)
        abline(h = 0.05, col = 2, lty = 2)
        abline(a = 0, b = 0.1 / length(pp), col = 2, lty = 2)
        mtext(c("BH-10%", "5%"), side = 4, at = c(0.1, 0.05), line = 0.1, las = 0, cex = 0.6)
        
        npar <- length(object$par)
        image(1:npar, 1:npar, cor(theta), xlab = "Parameter index", ylab = "Parameter index", main = "Parameter correlation")
        
    }
    invisible(list(deviance = deviance, pg = pg, prox = prox.samp, ll = ll, waic = fit.waic))
}

predict.sbde <- function(object, burn.perc = 0.5, nmc = 200, yRange = range(object$y), yLength = 401, ...){
    thin <- object$dim[8]
    nsamp <- object$dim[9]
    pars <- matrix(object$parsamp, ncol = nsamp)
    ss <- unique(round(nsamp * seq(burn.perc, 1, len = nmc + 1)[-1]))
    dimpars <- object$dim
    dimpars[7] <- length(ss)

    yGrid <- seq(yRange[1], yRange[2], len = yLength)
    n <- object$dim[1]; p <- 0; ngrid <- object$dim[5]
    dimpars[1] <- yLength
    
    pred <- .C("PRED", pars = as.double(pars[,ss]), yGrid = as.double(yGrid), hyper = as.double(object$hyper), dim = as.integer(dimpars), gridmats = as.double(object$gridmats), tau.g = as.double(object$tau.g), ldenssamp = double(length(ss)*yLength), dist = as.integer(object$fbase.choice))
    dens <- matrix(exp(pred$ldenssamp), ncol = length(ss))
    return(list(y = yGrid, fsamp = dens, fest = t(apply(dens, 1, quantile, pr = c(.025, .5, .975)))))
}


waic <- function(logliks, print = TRUE){
	lppd <- sum(apply(logliks, 1, logmean))
	p.waic.1 <- 2 * lppd - 2 * sum(apply(logliks, 1, mean))
	p.waic.2 <- sum(apply(logliks, 1, var))
	waic.1 <- -2 * lppd + 2 * p.waic.1
	waic.2 <- -2 * lppd + 2 * p.waic.2
	if(print) cat("WAIC.1 =", round(waic.1, 2), ", WAIC.2 =", round(waic.2, 2), "\n")
	invisible(c(WAIC1 = waic.1, WAIC2 = waic.2))
}




lamFn <- function(prox) return(sqrt(-100*log(prox)))
nuFn <- function(z) return(0.5 + 5.5*exp(z/2)) 
nuFn.inv <- function(nu) return(2*log((nu - 0.5)/5.5))
sigFn <- function(z, a.sig) return(exp(z/2)) 
sigFn.inv <- function(s, a.sig) return(2 * log(s))
unitFn <- function(u) return(pmin(1 - 1e-10, pmax(1e-10, u)))

sum.sq <- function(x) return(sum(x^2))
extract <- function(lo, vn) return(lo[[vn]])
logmean <- function(lx) return(max(lx) + log(mean(exp(lx - max(lx)))))
logsum <- function(lx) return(logmean(lx) + log(length(lx)))
shrinkFn <- function(x) return(1) ##(1/(1 + log(x)))
trape <- function(x, h, len = length(x)) return(c(0, cumsum(.5 * (x[-1] + x[-len]) * (h[-1] - h[-len]))))

getBands <- function(b, col = 2, lwd = 1, plot = TRUE, add = FALSE, x = seq(0,1,len=nrow(b)), remove.edges = TRUE, ...){
	colRGB <- col2rgb(col)/255
	colTrans <- rgb(colRGB[1], colRGB[2], colRGB[3], alpha = 0.2)
	b.med <- apply(b, 1, quantile, pr = .5)
	b.lo <- apply(b, 1, quantile, pr = .025)
	b.hi <- apply(b, 1, quantile, pr = 1 - .025)
	L <- nrow(b)
	ss <- 1:L; ss.rev <- L:1
	if(remove.edges){
		ss <- 2:(L-1); ss.rev <- (L-1):2
	}
	if(plot){
		if(!add) plot(x[ss], b.med[ss], ty = "n", ylim = range(c(b.lo[ss], b.hi[ss])), ...)  
		polygon(x[c(ss, ss.rev)], c(b.lo[ss], b.hi[ss.rev]), col = colTrans, border = colTrans)
		lines(x[ss], b.med[ss], col = col, lwd = lwd)  
	}
	invisible(cbind(b.lo, b.med, b.hi))
}


klGP <- function(lam1, lam2, nknots = 11){
	tau <- seq(0, 1, len = nknots)
	dd <- outer(tau, tau, "-")^2
	K1 <- exp(-lam1^2 * dd); diag(K1) <- 1 + 1e-10; R1 <- chol(K1); log.detR1 <- sum(log(diag(R1)))
	K2 <- exp(-lam2^2 * dd); diag(K2) <- 1 + 1e-10; R2 <- chol(K2); log.detR2 <- sum(log(diag(R2)))
	return(log.detR2-log.detR1 - 0.5 * (nknots - sum(diag(solve(K2, K1)))))
}



proxFn <- function(prox.Max, prox.Min, kl.step = 1){
	prox.grid <- prox.Max
	j <- 1
	while(prox.grid[j] > prox.Min){
		prox1 <- prox.grid[j]
		prox2 <- prox.Min
		kk <- klGP(lamFn(prox1), lamFn(prox2))
		while(kk > kl.step){
			prox2 <- (prox1 + prox2)/2
			kk <- klGP(lamFn(prox1), lamFn(prox2))
		}
		j <- j + 1
		prox.grid <- c(prox.grid, prox2)
	}
	return(prox.grid)
}

transform.grid <- function(w, ticks, dists){
    return((1-dists) * w[ticks] + dists * w[ticks+1])
}
