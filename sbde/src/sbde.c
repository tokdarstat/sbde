#include "R.h"
#include "Rmath.h"
#include "R_ext/Applic.h"

double log2(double x);
double *vect(int n);
int *ivect(int n);
double **mymatrix(int nr, int nc);
void Rprintvec(char *a, char *format, double *x, int n);
void Rprintmat(char *a, char *format, double **x, int m, int n, int flip);
void Rprintveci(char *a, char *format, int *x, int n);
double sumsquares(double *x, int n);
double inprod(double *x, double *y, int n);
double rnormtrunc(double mu, double sigma, double lo, double hi);
double vmax(double *x, int n);
double vmin(double *x, int n);
double logsum(double *lx, int n);
double logmean(double *lx, int n);
double sum(double *x, int n);
int rdraw(int n, double *lprob, int inlog);
void locator_string(int *ix, int n, char *a);
void locator_string_inverse(char *a, int *ix);
void mmprod(double **a, double **b, double **c, int m, int k, int n, int atrans, int btrans, int ctrans);
void mvprod(double **a, double *b, double *c, int m, int k, int atrans);
void set_lower_tri_zero(double **A, int n, int m);
void spchol(double **R, int N, double tol, int *pivot, int *rank, int max_rank,
			double gpcov(int, int, double*, int*, double **, int), double *covpar, int *include,
			double **K0, int addK, int max_scan, double *d, int dopivoting, 
			int dostopping, int padzero, double *lpen, double *d2);
void chol(double **R, int N, double tol, int *pivot, int *rank, int max_rank,
		  double *d, double **A, int dopivoting, int padzero, double eps);
void trisolve(double **R, int m, double *b, double *x, int transpose);
void triprod(double **R, int m, int n, double *x, double *b, int transpose);
void adMCMC(int niter, int thin, int npar, double *par, double **mu, double ***S, double *acpt_target, 
			double decay, int refresh, double *lm, double temp, int nblocks, int **blocks, int *blocks_size, 
			int *refresh_counter, int verbose, int ticker, double lpFn(double *, double), 
			double *parsamp, double *acptsamp, double *lpsamp);
void transform_grid(double *w, double *v, int *ticks, double *dists);
double part_trape(double target, double baseline, double a, double b, double Delta);

//global integers
int n, p, L, mid, m, nkap, ngrid, shrink, dist;

//value assigned global constants
double *taugrid, *akap, *bkap, *lpkap, asig, bsig, shrinkFactor, ***Agrid, ***Rgrid, *ldRgrid, *lpgrid, **x, *y, *wt;
int *cens;

//memory assigned global variables
double *lb;
double **wgrid, *llgrid, *zknot, *lw;
double *w0, *np_density, *np_cumu_density, *gam;

double sigFn(double z){
	//return sqrt(1.0 / qgamma(pnorm(z, 0.0, 1.0, 1, 1), asig, 1.0/bsig, 1, 1));	
	return exp(z/2.0);
} 	
double sigFn_inv(double s) {
	//return qnorm(pgamma(1.0 / s*s, asig, 1.0/bsig, 1, 1), 0.0, 1.0, 1, 1);
	return 2.0*log(s);
}
double nuFn(double z)  {
	return 0.5 + 5.5*exp(z/2.0);
}
double nuFn_inv(double nu) {
	return 2.0*log((nu - 0.5)/5.5);
}

double dgpd(double x, double nu, int in_log){
    double val = -(nu + 1.0) * log1p(x/nu);
    if(!in_log) val = exp(val);
    return val;
}
double pgpd(double x, double nu){
    return 1.0 - exp(-nu * log1p(x/nu));
}
double qgpd(double p, double nu){
    return nu * expm1(-log1p(-p)/nu);
}


double s0(double nu){
    double val;
    switch(dist){
        case 2:
            val = qt(0.95, nu, 1, 0);
            break;
        case 3:
            val = qgpd(0.5, nu);
            break;
        case 4:
            val = 1.0;
            break;
        default:
            val = qt(0.9, nu, 1, 0);
            break;
    }
    return val;
}

double f0(double x, double nu) {
    double val;
    switch (dist) {
        case 2:
            val = 2.0 * dt(x * s0(nu), nu, 0) * s0(nu);
            break;
        case 3:
            val = dgpd(x * s0(nu), nu, 0) * s0(nu);
            break;
        case 4:
            val = dunif(x, -1.0, 1.0, 0);
            break;
        default:
            val = dt(x * s0(nu), nu, 0) * s0(nu);
            break;
    }
    return val;
}

double log_f0(double x, double nu) {
    double val;
    switch (dist) {
        case 2:
            val = log(2.0) + dt(x * s0(nu), nu, 1) + log(s0(nu));
            break;
        case 3:
            val = dgpd(x * s0(nu), nu, 1) + log(s0(nu));
            break;
        case 4:
            val = dunif(x, -1.0, 1.0, 1);
            break;
        default:
            val = dt(x * s0(nu), nu, 1) + log(s0(nu));
            break;
    }
    return val;
}

double F0(double x, double nu) {
    double val;
    switch (dist) {
        case 2:
            val = 2.0 * (pt(x * s0(nu), nu, 1, 0) - 0.5);
            break;
        case 3:
            val = pgpd(x * s0(nu), nu);
            break;
        case 4:
            val = punif(x, -1.0, 1.0, 1, 0);
            break;
        default:
            val = pt(x * s0(nu), nu, 1, 0);
            break;
    }
    return val;
}

void trape(double *x, double *h, int n, double *c, int reverse){
	int i, j = 0;
	c[0] = 0.0;
	if(reverse){
		for(i = 1; i < n; i++){
			c[i] = c[i-1] + 0.5 * (h[j] - h[j-1]) * (x[j] + x[j-1]); 
			j--;
		}
	} else {
		for(i = 1; i < n; i++){
			c[i] = c[i-1] + 0.5 * (h[j+1] - h[j]) * (x[j] + x[j+1]); 
			j++;
		}
	}	
}

int locate(double u){
    int pos = floor((L - 1) * u);
    if(pos == L - 1) pos = L - 1;
    return pos;
}

double ppFn0(double *wknot, double *w, double *postgrid){
	int i, l;
	double akapm, zss;
	for(i = 0; i < ngrid; i++){
		mvprod(Agrid[i], wknot, wgrid[i], L, m, 0);
		trisolve(Rgrid[i], m, wknot, zknot, 1);
		zss = sumsquares(zknot, m);
		akapm = akap[0] + 0.5 * (double)m;
		llgrid[i] = -akapm * log1p(0.5 * zss / bkap[0]);
        //akapm = 0.1 + 0.5 * (double)m;
        //llgrid[i] = -akapm * log1p(0.5 * zss / 0.1);
		postgrid[i] = llgrid[i] - ldRgrid[i] + lpgrid[i];
	}
	
	double lps = logsum(postgrid, ngrid);
	for(i = 0; i < ngrid; i++) postgrid[i] = exp(postgrid[i] - lps);
	for(l = 0; l < L; l++) for(w[l] = 0.0, i = 0; i < ngrid; i++) w[l] += wgrid[i][l] * postgrid[i];
	return lps;
}	

double logpostFn(double *par, double temp, int llonly, double *ll, double *pg){
    
    int i, l, reach = 0;
    
    double lps0 = ppFn0(par, w0, pg);
    reach += m;
    
    double w0max = vmax(w0, L);
    for(l = 0; l < L; l++) np_density[l] = exp(w0[l] - w0max);
    trape(np_density, taugrid, L, np_cumu_density, 0);
    double totmass = np_cumu_density[L-1];
    for(l = 0; l < L; l++) np_density[l] /= totmass;
    
    for(i = 0; i < n; i++) ll[i] = 0.0;

    if(temp > 0.0){
        double gam0 = par[reach++];
        double sigma = sigFn(par[reach++]);
        double log_sigma = log(sigma);
        double nu = nuFn(par[reach++]);
        
        int loc_grid = 0;
        double y_std = 0.0, y_trans = 0.5;
        for(i = 0; i < n; i++){
            y_std = (y[i] - gam0) / sigma;
            y_trans = F0(y_std, nu);
            ll[i] = log_f0(y_std, nu) - log_sigma;
            loc_grid = locate(y_trans);
            ll[i] += log(((taugrid[loc_grid + 1] - y_trans) * np_density[loc_grid] + (y_trans - taugrid[loc_grid]) * np_density[loc_grid + 1])/(taugrid[loc_grid+1] - taugrid[loc_grid]));
        }
    }
    double lp = temp * inprod(ll, wt, n);
    if(!llonly) lp += lps0 + dlogis(par[m+2], 0.0, 1.0, 1);
    
    return lp;
}



double *llvec, *pgvec;
double lpFn(double *par, double temp){
    return logpostFn(par, temp, 0, llvec, pgvec);
}
double lpFn1(double *par, double *pgvec){
	return logpostFn(par, 1.0, 1, llvec, pgvec);
}

void SBDE(double *par, double *yVar, int *status, double *weights, double *hyper, int *dim, double *gridpars, double *tauG, double *muVar, double *SVar, int *blocksVar, int *blocks_size, double *dmcmcpar, int *imcmcpar, double *parsamp, double *acptsamp, double *lpsamp, int *distribution){
	
	int i, k, l;
	
	int reach = 0;
    n = dim[reach++]; p = 0; L = dim[reach++]; mid = dim[reach++];
	m = dim[reach++]; ngrid = dim[reach++]; nkap = dim[reach++];
    dist = distribution[0];
	int niter = dim[reach++], thin = dim[reach++], npar = m+3;
	taugrid = tauG;
	
	asig = hyper[0]; bsig = hyper[1];
	akap = vect(nkap); bkap = vect(nkap); lpkap = vect(nkap);
	for(reach = 2, i = 0; i < nkap; i++){
		akap[i] = hyper[reach++];
		bkap[i] = hyper[reach++];
		lpkap[i] = hyper[reach++];
	}
	
	reach = 0; 
	Agrid = (double ***)R_alloc(ngrid, sizeof(double **));
	Rgrid = (double ***)R_alloc(ngrid, sizeof(double **));
	ldRgrid = vect(ngrid);
	lpgrid = vect(ngrid);
	
	for(i = 0; i < ngrid; i++){
		Agrid[i] = mymatrix(L, m);
		for(l = 0; l < L; l++) for(k = 0; k < m; k++) Agrid[i][l][k] = gridpars[reach++];
		
		Rgrid[i] = mymatrix(m, m);
		for(k = 0; k < m; k++) for(l = 0; l < m; l++) Rgrid[i][l][k] = gridpars[reach++];
		
		ldRgrid[i] = gridpars[reach++];
		lpgrid[i] = gridpars[reach++];
	}
	
	y = yVar;
	cens = status;
    wt = weights;
	
	lb = vect(10);
	wgrid = mymatrix(ngrid, L);
	lw = vect(nkap);
	llgrid = vect(ngrid);
	zknot = vect(m);
	w0 = vect(L);
	np_density = vect(L);
	np_cumu_density = vect(L);
	llvec = vect(n);
	pgvec = vect(ngrid);
	
	int b;
	int nblocks = imcmcpar[0], refresh = imcmcpar[1], verbose = imcmcpar[2], ticker = imcmcpar[3];
	int *refresh_counter = imcmcpar + 4;
	
	double temp = dmcmcpar[0], decay = dmcmcpar[1];
	double *acpt_target = dmcmcpar + 2, *lm = dmcmcpar + 2 + nblocks;
	
	double **mu = (double **)R_alloc(nblocks, sizeof(double *));
	double ***S = (double ***)R_alloc(nblocks, sizeof(double **));
	int **blocks = (int **)R_alloc(nblocks, sizeof(int *));
	
	int mu_point = 0, S_point = 0, blocks_point = 0;
	for(b = 0; b < nblocks; b++){
		mu[b] = muVar + mu_point; mu_point += blocks_size[b];
		S[b] = (double **)R_alloc(blocks_size[b], sizeof(double *));
		for(i = 0; i < blocks_size[b]; i++){
			S[b][i] = SVar + S_point;
			S_point += blocks_size[b];
		}
		blocks[b] = blocksVar + blocks_point; blocks_point += blocks_size[b];
	}
	adMCMC(niter, thin, npar, par, mu, S, acpt_target, decay, refresh, lm, temp, 
		   nblocks, blocks, blocks_size, refresh_counter, verbose, ticker, 
		   lpFn, parsamp, acptsamp, lpsamp);
	
	
}


void DEV(double *par, double *yVar, int *status, double *weights, double *hyper, int *dim, double *gridpars, double *tauG, double *devsamp, double *llsamp, double *pgsamp, int *distribution){
    
    int i, k, l;
    
    int reach = 0;
    n = dim[reach++]; p = 0; L = dim[reach++]; mid = dim[reach++];
    m = dim[reach++]; ngrid = dim[reach++]; nkap = dim[reach++];
    int niter = dim[reach++], npar = (m+1) + 2;
    
    dist = distribution[0];
    
    reach = 0;
    taugrid = tauG;
    asig = hyper[0]; bsig = hyper[1];
    akap = vect(nkap); bkap = vect(nkap); lpkap = vect(nkap);
    for(reach = 2, i = 0; i < nkap; i++){
        akap[i] = hyper[reach++];
        bkap[i] = hyper[reach++];
        lpkap[i] = hyper[reach++];
    }
    
    reach = 0;
    Agrid = (double ***)R_alloc(ngrid, sizeof(double **));
    Rgrid = (double ***)R_alloc(ngrid, sizeof(double **));
    ldRgrid = vect(ngrid);
    lpgrid = vect(ngrid);
    
    for(i = 0; i < ngrid; i++){
        Agrid[i] = mymatrix(L, m);
        for(l = 0; l < L; l++) for(k = 0; k < m; k++) Agrid[i][l][k] = gridpars[reach++];
        
        Rgrid[i] = mymatrix(m, m);
        for(k = 0; k < m; k++) for(l = 0; l < m; l++) Rgrid[i][l][k] = gridpars[reach++];
        
        ldRgrid[i] = gridpars[reach++];
        lpgrid[i] = gridpars[reach++];
    }
    
    y = yVar;
    cens = status;
    wt = weights;
    
    lb = vect(10);
    wgrid = mymatrix(ngrid, L);
    lw = vect(nkap);
    llgrid = vect(ngrid);
    zknot = vect(m);
    w0 = vect(L);
    np_density = vect(L);
    np_cumu_density = vect(L);
    
    reach = 0;
    int iter, reach2 = 0, reach3 = 0;
    for(iter = 0; iter < niter; iter++){
        devsamp[iter] = -2.0 * logpostFn(par + reach, 1.0, 1, llsamp + reach2, pgsamp + reach3);
        reach += npar; reach2 += n; reach3 += ngrid;
    }
}




void PRED(double *par, double *yGrid, double *hyper, int *dim, double *gridpars, double *tauG, double *logdenssamp, int *distribution){
    
    int i, k, l;
    
    dist = distribution[0];
    int reach = 0;
    n = dim[reach++]; p = 0; L = dim[reach++]; mid = dim[reach++];
    m = dim[reach++]; ngrid = dim[reach++]; nkap = dim[reach++];
    int niter = dim[reach++], npar = (m+1) + 2;
    
    reach = 0;
    taugrid = tauG;
    asig = hyper[0]; bsig = hyper[1];
    akap = vect(nkap); bkap = vect(nkap); lpkap = vect(nkap);
    for(reach = 2, i = 0; i < nkap; i++){
        akap[i] = hyper[reach++];
        bkap[i] = hyper[reach++];
        lpkap[i] = hyper[reach++];
    }
    
    reach = 0;
    Agrid = (double ***)R_alloc(ngrid, sizeof(double **));
    Rgrid = (double ***)R_alloc(ngrid, sizeof(double **));
    ldRgrid = vect(ngrid);
    lpgrid = vect(ngrid);
    
    for(i = 0; i < ngrid; i++){
        Agrid[i] = mymatrix(L, m);
        for(l = 0; l < L; l++) for(k = 0; k < m; k++) Agrid[i][l][k] = gridpars[reach++];
        
        Rgrid[i] = mymatrix(m, m);
        for(k = 0; k < m; k++) for(l = 0; l < m; l++) Rgrid[i][l][k] = gridpars[reach++];
        
        ldRgrid[i] = gridpars[reach++];
        lpgrid[i] = gridpars[reach++];
    }
    
    y = yGrid;
    cens = ivect(n); for(i = 0; i < n; i++) cens[i] = 0;
    wt = vect(n); for(i = 0; i < n; i++) wt[i] = 1.0;
    
    lb = vect(10);
    wgrid = mymatrix(ngrid, L);
    lw = vect(nkap);
    llgrid = vect(ngrid);
    zknot = vect(m);
    w0 = vect(L);
    np_density = vect(L);
    np_cumu_density = vect(L);
    
    reach = 0;
    double lldummy, *pgdummy = vect(ngrid);
    int iter, reach2 = 0, reach3 = 0;
    for(iter = 0; iter < niter; iter++){
        lldummy = logpostFn(par + reach, 1.0, 1, logdenssamp + reach2, pgdummy);
        reach += npar; reach2 += n; reach3 += ngrid;
    }
}


// ------ mydefs ------ //

double log2(double x){
	return log(x)/log(2.0);
}

double * vect(int n){
	return (double *)R_alloc(n, sizeof(double));
}

int * ivect(int n){
	return (int *)R_alloc(n, sizeof(int));
}

double ** mymatrix(int nr, int nc){
	int   i;
	double **m;	
	m = (double **) R_alloc(nr, sizeof(double *));
	for (i = 0; i < nr; i++)
		m[i] = (double *) R_alloc(nc, sizeof(double));
	return m;
}


void Rprintvec(char *a, char *format, double *x, int n){
	int i;
	Rprintf("%s", a);
	for(i = 0; i < n; i++)
		Rprintf(format, x[i]);
	Rprintf("\n");
}

void Rprintmat(char *a, char *format, double **x, int m, int n, int flip){
	int i, j;
	Rprintf("%s\n", a);
	for(i = 0; i < m; i++){
		for(j = 0; j < n; j++)
			Rprintf(format, x[i][j]);
		Rprintf("\n");
	}
}


void Rprintveci(char *a, char *format, int *x, int n){
	int i;
	Rprintf("%s", a);
	for(i = 0; i < n; i++)
		Rprintf(format, x[i]);
	Rprintf("\n");
}


double sumsquares(double *x, int n){
	double ss = 0.0;
	int i;
	for(i = 0; i < n; i++)
		ss += x[i] * x[i];
	return ss;
}

double inprod(double *x, double *y, int n){
	double ip = 0.0;
	int i;
	for(i = 0; i < n; i++)
		ip += x[i] * y[i];
	return ip;
}


double rnormtrunc(double mu, double sigma, double lo, double hi){
	double u = runif(0.0, 1.0);
	double p = u * pnorm(hi, mu, sigma, 1, 0) + (1.0 - u) * pnorm(lo, mu, sigma, 1, 0);
	if(p <= 0.0) p = 1.0e-10;
	if(p >= 1.0) p = 1.0 - 1.0e-10;
	return qnorm(p, mu, sigma, 1, 0);
}



double matern(double x, double phi, int kappa){ 
	/*Returns the Matern function for x, phi and kappa.*/ 
	
	/* Variables */ 
	double ans, cte; 
	double uphi=x/phi; 
	
	/* Matern */ 
	
	if (uphi==0) return 1; 
	else{ 
		if (kappa==0.5) 
			ans = exp(-uphi); 
		else { 
			cte = R_pow(2, (-(kappa-1)))/gammafn(kappa); 
			ans = cte * R_pow(uphi, kappa) * bessel_k(uphi,kappa,1); 
		} 
	} 
	
	return ans; 
} 

double vmax(double *x, int n){
	int i;
	double xmax = x[0];
	for(i = 1; i < n; i++) if(x[i] > xmax) xmax = x[i];
	return xmax;
}

double logsum(double *lx, int n){
	double lxmax = vmax(lx, n), a = 0.0;
	int i;
	for(i = 0; i < n; i++) a += exp(lx[i] - lxmax);
	return lxmax + log(a);
}

double logmean(double *lx, int n){
	return logsum(lx, n) - log((double)n);
} 

double sum(double *x, int n){
	double a = 0.0;
	int i;
	for(i = 0; i < n; i++) a += x[i];
	return a;
}

int rdraw(int n, double *prob, int inlog){
	double psum, u = runif(0.0, 1.0), cprob;
	int j = 0;
	
	if(inlog){
		psum = logsum(prob, n);	
		cprob = exp(prob[0] - psum);
		while(u > cprob && j < n - 1){
			j++;
			cprob += exp(prob[j] - psum);
		}
	} else {
		psum = sum(prob, n);
		cprob = prob[0] / psum;
		while(u > cprob && j < n - 1){
			j++;
			if(prob[j] > 0.0) cprob += prob[j] / psum;
		}
	}
	return j;
}


void locator_string(int *ix, int n, char *a){
	const char *fmt[2];	fmt[0] = "%d";	fmt[1] = ".%d";	
	int i, skip = 0;
	for(i = 0; i < n; i++){
		if(ix[i]){
			sprintf(a + skip, fmt[skip > 0], i + 1);
			skip = strlen(a);
		}
	}
}


void locator_string_inverse(char *a, int *ix){
	const char s[2] = ".";
	char *token;
	token = strtok(a, s);
	while( token != NULL ) {
		ix[atoi(token) - 1] = 1;
		token = strtok(a, s);
	}
}


void mmprod(double **a, double **b, double **c, int m, int k, int n, int atrans, int btrans, int ctrans){
	int i, j, l;
	if(!ctrans){
		if(atrans && btrans){
			for(i = 0; i < m; i++)
				for(j = 0; j < n; j++)
					for(c[i][j] = 0.0, l = 0; l < k; l++) c[i][j] += a[l][i] * b[j][l];
		} else if (!atrans && btrans){
			for(i = 0; i < m; i++)
				for(j = 0; j < n; j++)
					for(c[i][j] = 0.0, l = 0; l < k; l++) c[i][j] += a[i][l] * b[j][l];		
		} else if (atrans && !btrans){
			for(i = 0; i < m; i++)
				for(j = 0; j < n; j++)
					for(c[i][j] = 0.0, l = 0; l < k; l++) c[i][j] += a[l][i] * b[l][j];		
		} else {
			for(i = 0; i < m; i++)
				for(j = 0; j < n; j++)
					for(c[i][j] = 0.0, l = 0; l < k; l++) c[i][j] += a[i][l] * b[l][j];
		}
	} else {
		if(atrans && btrans){
			for(i = 0; i < m; i++)
				for(j = 0; j < n; j++)
					for(c[i][j] = 0.0, l = 0; l < k; l++) c[j][i] += a[l][i] * b[j][l];
		} else if (!atrans && btrans){
			for(i = 0; i < m; i++)
				for(j = 0; j < n; j++)
					for(c[i][j] = 0.0, l = 0; l < k; l++) c[j][i] += a[i][l] * b[j][l];		
		} else if (atrans && !btrans){
			for(i = 0; i < m; i++)
				for(j = 0; j < n; j++)
					for(c[i][j] = 0.0, l = 0; l < k; l++) c[j][i] += a[l][i] * b[l][j];		
		} else {
			for(i = 0; i < m; i++)
				for(j = 0; j < n; j++)
					for(c[i][j] = 0.0, l = 0; l < k; l++) c[j][i] += a[i][l] * b[l][j];
		}		
	}
}

void mvprod(double **a, double *b, double *c, int m, int k, int atrans){
	int i, l;
	if(atrans){
		for(i = 0; i < m; i++)
			for(c[i] = 0.0, l = 0; l < k; l++) c[i] += a[l][i] * b[l];
	} else {
		for(i = 0; i < m; i++)
			for(c[i] = 0.0, l = 0; l < k; l++) c[i] += a[i][l] * b[l];
	}
}

double vmin(double *x, int n){
	int i;
	double xmin = x[0];
	for(i = 1; i < n; i++) if(x[i] < xmin) xmin = x[i];
	return xmin;
}


//----- spchol ----//

void set_lower_tri_zero(double **A, int n, int m ){
	int i, j;
	for(i = 0; i < n; i++)
		for(j = i + 1; j < m; j++)
			A[j][i] = 0.0;
}


void spchol(double **R, int N, double tol, int *pivot, int *rank, int max_rank,
			double gpcov(int, int, double*, int*, double**, int), double *covpar, int *include,
			double **K0, int addK, int max_scan, double *d, int dopivoting, 
			int dostopping, int padzero, double *lpen, double *d2){
	
	
	// sparse cholesky factorization with pivoting and diagonal augmentation
	// accepts an empty matrix R and a function gpcov(i, j, ...) that is called 
	// to compute the (i,j)-th element of the original covariance matrix
	
	
	set_lower_tri_zero(R, N, max_rank);
	
	int i, a, l;
	double u, b;
	
	if(dopivoting){
		for(i = 0; i < N; i++)
			pivot[i] = i;
	}
	for(i = 0; i < N; i++) d[i] = gpcov(pivot[i], pivot[i], covpar, include, K0, addK);
	for(i = 0; i < max_scan; i++) d2[i] = lpen[pivot[i]];
	
	int k = 0, max_diag;
	for(max_diag = k, i = k + 1; i < max_scan; i++)
		if(d[i] > d[max_diag] * exp(d2[max_diag] - d2[i]))
			max_diag = i;
	tol *= d[max_diag];
	int flag = (k < max_rank);
	if(dostopping)
		flag = (d[max_diag] > tol);	
	
	while(flag){
		if(dopivoting){
			if(max_diag > k){
				a = pivot[k];
				pivot[k] = pivot[max_diag];
				pivot[max_diag] = a;
				
				b = d[k];
				d[k] = d[max_diag];
				d[max_diag] = b;
				
				b = d2[k];
				d2[k] = d2[max_diag];
				d2[max_diag] = b;
				
				for(i = 0; i < k; i++){
					b = R[i][k];
					R[i][k] = R[i][max_diag];
					R[i][max_diag] = b;
				}
			}
		}
		
		R[k][k] = sqrt(d[k]);
		
		for(i = k + 1; i < N; i++){
			u = gpcov(pivot[i], pivot[k], covpar, include, K0, addK);
			for(R[k][i] = u, l = 0; l < k; l++)
				R[k][i] -= R[l][i] * R[l][k];
			R[k][i] /= R[k][k];
			d[i] -= R[k][i] * R[k][i];
		}
		
		k++;
		flag = (k < max_rank);
		if(flag && dostopping){
			for(max_diag = k, i = k + 1; i < max_scan; i++)
				if(d[i] > d[max_diag] * exp(d2[max_diag] - d2[i]))
					max_diag = i;
			flag = (d[max_diag] > tol);	
		}
	}
	
	rank[0] = k;
	if(padzero){
		for(l = k; l < N; l++)
			d[l] = 0.0;
	}
}

void chol(double **R, int N, double tol, int *pivot, int *rank, int max_rank,
		  double *d, double **A, int dopivoting, int padzero, double eps){
	
	set_lower_tri_zero(R, N, max_rank);
	
	int i, a, l;
	double u, b;
	
	for(i = 0; i < N; i++){
		pivot[i] = i;
		d[i] = A[i][i] + eps * (1.0 + A[i][i]);
	}
	
	int k = 0, max_diag;
	for(max_diag = k, i = k + 1; i < N; i++)
		if(d[i] > d[max_diag])
			max_diag = i;
	int flag = (d[max_diag] > tol);	
	
	
	while(flag){
		if(dopivoting){
			if(max_diag > k){
				a = pivot[k];
				pivot[k] = pivot[max_diag];
				pivot[max_diag] = a;
				
				b = d[k];
				d[k] = d[max_diag];
				d[max_diag] = b;
				
				for(i = 0; i < k; i++){
					b = R[i][k];
					R[i][k] = R[i][max_diag];
					R[i][max_diag] = b;
				}
			}
		}
		
		R[k][k] = sqrt(d[k]);
		
		for(i = k + 1; i < N; i++){
			u = A[pivot[i]][pivot[k]];
			for(R[k][i] = u, l = 0; l < k; l++)
				R[k][i] -= R[l][i] * R[l][k];
			R[k][i] /= R[k][k];
			d[i] -= R[k][i] * R[k][i];
		}
		
		k++;
		flag = (k < max_rank);
		if(flag){
			for(max_diag = k, i = k + 1; i < N; i++)
				if(d[i] > d[max_diag])
					max_diag = i;
			flag = (d[max_diag] > tol);	
		}
	}
	
	rank[0] = k;
	if(padzero){
		for(l = k; l < N; l++)
			d[l] = 0.0;
	}
}



void trisolve(double **R, int m, double *b, double *x, int transpose){
	
	int i, j;
	if(transpose){
		for(j = 0; j < m; j++){
			for(x[j] = b[j], i = 0; i < j; i++)
				x[j] -= x[i] * R[i][j];
			x[j] /= R[j][j];
		}
	} else {	
		for(j = m - 1; j >= 0; j--){
			for(x[j] = b[j], i = j + 1; i < m; i++) 
				x[j] -= R[j][i] * x[i];
			x[j] /= R[j][j];
		}
	}	  
}



void triprod(double **R, int m, int n, double *x, double *b, int transpose){
	
	int i, j;
	if(transpose){
		for(i = 0; i < m; i++)
			for(b[i] = 0.0, j = 0; j <= i; j++)
				b[i] += R[j][i] * x[j];
		for(; i < n; i++)
			for(b[i] = 0.0, j = 0; j < m; j++)
				b[i] += R[j][i] * x[j];
	} else{
		for(i = 0; i < m; i++)
			for(b[i] = 0.0, j = i; j < n; j++)
				b[i] += R[i][j] * x[j];
	}
}

//------ mcmc ------//

void adMCMC(int niter, int thin, int npar, double *par, double **mu, double ***S, double *acpt_target, 
			double decay, int refresh, double *lm, double temp, int nblocks, int **blocks, int *blocks_size, 
			int *refresh_counter, int verbose, int ticker, double lpFn(double *, double), 
			double *parsamp, double *acptsamp, double *lpsamp){
	
	int b, i, j;
	double ***R = (double ***)R_alloc(nblocks, sizeof(double **));
	int **blocks_pivot = (int **)R_alloc(nblocks, sizeof(double *)), *blocks_rank = ivect(nblocks);
	double *blocks_d = vect(npar);
	for(b = 0; b < nblocks; b++){
		R[b] = mymatrix(blocks_size[b], blocks_size[b]);
		blocks_pivot[b] = ivect(blocks_size[b]);
		chol(R[b], blocks_size[b], 1.0e-8, blocks_pivot[b], blocks_rank + b, blocks_size[b], blocks_d, S[b], 0, 0, 1.0e-10);
	}
	
	int *chunk_size = ivect(nblocks);
	double *acpt_chunk = vect(nblocks);
	double **parbar_chunk = (double **)R_alloc(nblocks, sizeof(double *));
	double *alpha = vect(nblocks);
	double *frac = vect(nblocks);
	//double *ppick = vect(nblocks);
	
	for(b = 0; b < nblocks; b++){
		chunk_size[b] = 0;
		acpt_chunk[b] = 0.0;
		parbar_chunk[b] = vect(blocks_size[b]);		
		for(i = 0; i < blocks_size[b]; i++) parbar_chunk[b][i] = 0.0;		
		frac[b] = sqrt(1.0 / ((double)refresh_counter[b] + 1.0));
	}
	
	double lpval, lpvalnew, lp_diff;
	double **parstore = mymatrix(2, npar), *par_incr = vect(npar), *zsamp = vect(npar);
	int ipar = 0, iparnew = 1;
	for(i = 0; i < npar; i++) parstore[0][i] = par[i];
	
	lpval = lpFn(parstore[ipar], temp);
	if(verbose) Rprintf("Initial lp = %g\n", lpval);
	int iter, store_lp = 0, store_par = 0, store_acpt = 0;
	double *alpha_run = vect(nblocks), chs, lambda;
	for(b = 0; b < nblocks; b++) alpha_run[b] = 0.0;
	//	int run_counter = 0;
	int *run_counter = ivect(nblocks);
	for(b = 0; b < nblocks; b++) run_counter[b] = 0;
	
	GetRNGstate();
	for(iter = 0; iter < niter; iter++){
		//		for(b = 0; b < nblocks; b++) ppick[b] = fabs(log(1.0e-6 + acpt_chunk[b]) - log(acpt_target[0]));
		for(b = 0; b < nblocks; b++){
			//b = floor(nblocks * runif(0.0, 1.0));
			//		b = rdraw(nblocks, ppick, 0);
			chunk_size[b]++;
			for(i = 0; i < blocks_size[b]; i++) zsamp[i] = rnorm(0.0, 1.0);
			triprod(R[b], blocks_size[b], blocks_size[b], zsamp, par_incr, 1);
			for(i = 0; i < npar; i++) parstore[iparnew][i] = parstore[ipar][i];
			lambda = lm[b] * sqrt(3.0 / rgamma(3.0, 1.0));
			for(i = 0; i < blocks_size[b]; i++) parstore[iparnew][blocks[b][i]] += lambda * par_incr[i];
			lpvalnew = lpFn(parstore[iparnew], temp);
			lp_diff = lpvalnew - lpval;
			alpha[b] = exp(lp_diff); if(alpha[b] > 1.0) alpha[b] = 1.0;
			if(log(runif(0.0, 1.0)) < lp_diff){      
				ipar = iparnew;
				iparnew = !ipar;
				lpval = lpvalnew;
			}			
			alpha_run[b] = ((double)run_counter[b] * alpha_run[b] + alpha[b]) / ((double)(run_counter[b] + 1.0));
			run_counter[b]++;
		}
		if((iter + 1) % thin == 0){
			lpsamp[store_lp++] = lpval;
			for(i = 0; i < npar; i++) parsamp[store_par++] = parstore[ipar][i];
			for(b = 0; b < nblocks; b++) acptsamp[store_acpt++] = alpha[b];
		}
		
		if(verbose){
			if(niter < ticker || (iter + 1) % (niter / ticker) == 0){
				Rprintf("iter = %d, lp = %g ", iter + 1, lpval);
				Rprintvec("acpt = ", "%0.2f ", alpha_run, nblocks);
				for(b = 0; b < nblocks; b++){
					alpha_run[b] = 0.0;
					run_counter[b] = 0;
				}
			}
		}
		
		for(b = 0; b < nblocks; b++){
			chs = (double) chunk_size[b]; if(chs < 1.0) chs = 1.0;
			acpt_chunk[b] = acpt_chunk[b] + (alpha[b] - acpt_chunk[b]) / chs;
			for(i = 0; i < blocks_size[b]; i++)
				parbar_chunk[b][i] += (parstore[ipar][blocks[b][i]] - parbar_chunk[b][i]) / chs;  
			
			if(chunk_size[b] == refresh * blocks_size[b]){
				refresh_counter[b]++;
				frac[b] = sqrt(1.0 / ((double)refresh_counter[b] + 1.0));
				for(i = 0; i < blocks_size[b]; i++){
					for(j = 0; j < i; j++){
						S[b][i][j] = (1.0 - frac[b]) * S[b][i][j] + frac[b] * (parbar_chunk[b][i] - mu[b][i]) * (parbar_chunk[b][j] - mu[b][j]);
						S[b][j][i] = S[b][i][j];
					}
					S[b][i][i] = (1.0 - frac[b]) * S[b][i][i] + frac[b] * (parbar_chunk[b][i] - mu[b][i]) * (parbar_chunk[b][i] - mu[b][i]);
				}
				chol(R[b], blocks_size[b], 1.0e-8, blocks_pivot[b], blocks_rank + b, blocks_size[b], blocks_d, S[b], 0, 0, 1.0e-10);
				for(i = 0; i < blocks_size[b]; i++) mu[b][i] += frac[b] * (parbar_chunk[b][i] - mu[b][i]);
				lm[b] *= exp(frac[b] * (acpt_chunk[b] - acpt_target[b]));
				acpt_chunk[b] = 0.0;
				for(i = 0; i < blocks_size[b]; i++) parbar_chunk[b][i] = 0.0;
				chunk_size[b] = 0;
			}
		}
	}
	PutRNGstate();
	for(i = 0; i < npar; i++) par[i] = parstore[ipar][i];
}


void transform_grid(double *w, double *v, int *ticks, double *dists){
    int l;
    for(l = 0; l < L; l++) v[l] = (1.0 - dists[l]) * w[ticks[l]] + dists[l] * w[ticks[l]+1];
}

double part_trape(double target, double baseline, double a, double b, double Delta){
    double h = (-a*Delta + sqrt(a*a * Delta*Delta + 2.0*Delta*(b-a)*(target - baseline)))/((b-a)*Delta);
    return ((1.0 - h)*a + h*b);
}
