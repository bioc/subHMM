/* History: Aug 13 2018 Initial coding
*/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Memory.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define DEBUG 0

#define CHECK_MEM(obj) if (obj == NULL) {error("Memory");}
#define TWOPI 6.2831853071795862
#define LOG2 0.69314718055994529
#define LOG0ARG 1.0e-300
#define LOG0 -690.77552789821368
#define MAX_EP5 1.0e300
#define DMISSTEST_LT -999999998.0
#define NPARM 6
#define LARGE_NEG -1.0e100

/* For optimizer */
#define stepredn	0.2
#define acctol	0.0001
#define reltest	10.0
#define expTo0       -9.0e100

#define IARG_mainJ 0
#define IARG_tauJ 1
#define IARG_J 2
#define IARG_n 3
#define IARG_maxiter 4
#define IARG_print 5
#define IARG_DEBUG 6

#define DARG_pr0 0
#define DARG_lpr0 1
#define DARG_sigma20_lr 2
#define DARG_lsigma20_lr 3
#define DARG_lpsi0 4
#define DARG_sigma20_lor 5
#define DARG_lsigma20_lor 6
#define DARG_lsz0 7
#define DARG_P5_LB 8
#define DARG_V0_LB 9
#define DARG_V0_UB 10
#define DARG_v0 11
#define DARG_loglike_eps 12
#define DARG_parm_eps 13
#define DARG_reltol 14

void akh_ca(int*, int*, int*, int*, double*, double*, double*, double*,
            double*, double*, double*, double*, double*, double*, int*,
            double*, double*, double*, double*, int*, double*, double*);
void bkh_cb(int*, int*, double*, double*, double*, double*, double*,
            double*, double*, double*, double*, int*, double*, double*,
            int*, double*, double*);
void locFunc_hatp0(int*, int*, int*, int*, double*, double*, double*, 
                   double*, double*, double*, double*, double*, double*, 
                   double*, int*, double*, double*, double*, double*, 
                   double*, int*, double*, double*, double*, double*, double*);
void C_hes_logL1(int*, int*, double*, double*, double*, double*, double*, 
                 double*, int*, double*, double*, int*, double*);
void C_myoptimC(double*, int*, double*, double*, double*, int*, int*, int*, 
                double*, double*, int*, double*, double*, double*, 
                double*, double*, double*, double*, double*, int*);
void C_EM_alg(int*, double*, double*, double*, double*, double*, double*, 
              double*, double*, double*, double*, double*, double*, double*, 
              double*, double*, double*, int*, int*, double*, double*, 
              double*, double*, double*);

static const R_CMethodDef callMethods[] = {
  {"akh_ca", (DL_FUNC)&akh_ca, 22},
  {"bkh_cb", (DL_FUNC)&bkh_cb, 17},
  {"locFunc_hatp0", (DL_FUNC)&locFunc_hatp0, 26},
  {"C_hes_logL1", (DL_FUNC)&C_hes_logL1, 13},
  {"C_myoptimC", (DL_FUNC)&C_myoptimC, 20},
  {"C_EM_alg", (DL_FUNC)&C_EM_alg, 24},
  {NULL, NULL, 0}
};

/* Used for debugging
void writeVec(fid, vec, len, ncol)
FILE *fid;
double *vec;
int len, ncol;
{
  int i, m=0;

  for (i=0; i<len-1; i++) {
    fprintf(fid, "%20.11g\t", vec[i]);
    m++;
    if (ncol && (m == ncol)) {
      m = 0;
      fprintf(fid, "\n");
    }
  }
  fprintf(fid, "%20.11g\n", vec[len-1]);

} 
void setVecToVal(vec, n, val)
double *vec, val;
int n;
{
  int i;
  for (i=0; i<n; i++) vec[i] = val;

}
void setMatToVal(mat, nr, nc, val)
double **mat, val;
int nr, nc;
{
  int i, j;
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) mat[i][j] = val;
  }
}
void print_dVec(vec, n, n2, name)
double *vec;
int n, n2;
char name[10];
{
  int i, m=0;
  Rprintf("%s \n", name);
  for (i=0; i<n; i++) {
    Rprintf(" %g ", vec[i]);
    m++;
    if (n2 && (m == n2)) {
      Rprintf("\n");
      m = 0;
    }
  }
  Rprintf("\n \n");
}
void print_dMat(mat, nr, nc, name)
double **mat;
int nr, nc;
char name[10];
{
  int i, j;
  Rprintf("%s \n", name);
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) Rprintf(" %g ", mat[i][j]);
    Rprintf("\n");
  }
  Rprintf("\n");
}
void print_sumAbsVec(vec, n) 
double *vec;
int n;
{
  int i;
  double sum=0.0;
  for (i=0; i<n; i++) sum += fabs(vec[i]);
  Rprintf("absSum=%g %40.8f\n", sum, sum);

}
void print_sumAbsMat(mat, nr, nc) 
double **mat;
int nr, nc;
{
  int i, j;
  double sum=0.0;
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) sum += fabs(mat[i][j]);
  }
  Rprintf("absSum=%g %40.8f\n", sum, sum);

}
*/


static char * cVec_alloc(n, initFlag, initVal)
int n, initFlag;
char initVal;
{
  int i;
  char *ret, *p;

  ret = (char *) malloc(n*sizeof(char));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }

  return(ret);

} /* END: cVec_alloc */

/* Function to allocate memory for a double vector */
static double * dVec_alloc(n, initFlag, initVal)
int n, initFlag;
double initVal;
{
  int i;
  double *ret, *p;

  ret = (double *) malloc(n*sizeof(double));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }

  return(ret);

} /* END: dVec_alloc */

static double ** ptrdVec_alloc(n, initFlag)
int n, initFlag;
{
  int i;
  double **ret, **p;

  ret = (double **) malloc(n*sizeof(double *));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = NULL;
  }

  return(ret);

} /* END: dVec_alloc */


/* Function to allocate a double matrix */
static double ** dMat_alloc(nrow, ncol, initFlag, initVal)
int nrow, ncol, initFlag;
double initVal;
{
  double **mat, **ptr;
  int i;

  mat = (double **) malloc(nrow*sizeof(double *));
  CHECK_MEM(mat);
  for (i=0, ptr=mat; i<nrow; i++, ptr++) *ptr = dVec_alloc(ncol, initFlag, initVal);

  return(mat);

} /* END: dMat_alloc */

/* Function to allocate a double array (3-dim) */
/*
static double *** dArray_alloc(n1, n2, n3, initFlag, initVal)
int n1, n2, n3, initFlag;
double initVal;
{
  double ***mat;
  int i, j;

  mat = (double ***) malloc(n1*sizeof(double **));
  CHECK_MEM(mat);
  for (i=0; i<n1; i++) {
    mat[i] = (double **) malloc(n2*sizeof(double *));
    CHECK_MEM(mat[i]);
    for (j=0; j<n2; j++) mat[i][j] = dVec_alloc(n3, initFlag, initVal);
  }
     
  return(mat);

}  */

/* Function to free a matrix */
static void matrix_free(x, n)
void **x;
int n;
{
  int i;
  for (i=0; i<n; i++) {
    if (x[i]) free(x[i]);
  }
  free(x);

} /* END: matrix_free */

/* Function to free an aray (3d) */
/*
static void array_free(x, n1, n2)
void ***x;
int n1, n2;
{
  int i, j;
  for (i=0; i<n1; i++) {
    for (j=0; j<n2; j++) free(x[i][j]);
    free(x[i]);
  }
  free(x);

}  */


/* Function to fill in a matrix from a vector (by column) */
static void fillMat(vec, nr, nc, out)
double *vec, **out;
int nr, nc;
{
  int i, j, col=0, ii;

  ii = 0;
  for (j=0; j<nc; j++) {
    for (i=0; i<nr; i++) {
      out[i][col] = vec[ii];
      ii++;
    }
    col++;
  }

} /* END: fillMat */

static double maxAbsDiff(v1, v2, n)
double *v1, *v2;
int n;
{
  int i;
  double val, maxval=-1.0, *p1, *p2;

  for (i=0, p1=v1, p2=v2; i<n; i++, p1++, p2++) {
    val = fabs(*p1 - *p2);
    if (val > maxval) maxval = val;
  }

  return(maxval);

} /* END: maxAbsDiff */

static void copy_dVec(from, into, n)
double *from, *into;
int n;
{
  int i;
  double *p1, *p2;

  for (i=0, p1=from, p2=into; i<n; i++, p1++, p2++) *p2 = *p1;

} /* END: copy_dVec */

static double sumvec(vec, n)
double *vec;
int n;
{
  int i;
  double ret=0.0, *p;

  for (i=0, p=vec; i<n; i++, p++) ret += *p;

  return(ret);

}

static void logvec(vec, n, ret)
double *vec, *ret;
int n;
{
  int i;
  double *p1, *p2;

  for (i=0, p1=vec, p2=ret; i<n; i++, p1++, p2++) *p2 = log(*p1);

}

static void reverseCumSum(vec, n, ret)
double *vec, *ret; 
int n;
{
  int i;
  double sum;

  if (ret == vec) error("INTERNAL CODING ERROR in reverseCumSum");
  sum    = sumvec(vec, n);

  ret[0] = sum;
  for (i=1; i<n; i++) ret[i] = ret[i-1] - vec[i-1]; 

}

/* The pzzks array always get summed over the 3rd dimension, so we really only
   need it to be 2 dimensions. */
static void get_pzzks(n, J, hatp0, akh_all, bkh_all, tmp0, tmp1, tmp2, tmp3, 
             t1, t2, lr, logor2, mu0_lrs, notMissing, ncpVecJ, vec, tvec, ret)
int n, J, *notMissing;
double t1, t2, tmp0, tmp1, tmp2, tmp3, **hatp0, **akh_all, **bkh_all;
double *logor2, *lr, *mu0_lrs, *ncpVecJ, *vec, *tvec, **ret;
{
  int k, i, j;
  double *vd1, xchisq, lrk, d1, d2, *bvec, *MAT, *avec, *pd;

  vd1 = dVec_alloc(J, 0, 0.0); 
  MAT = dVec_alloc(J, 0, 0.0); 

  for (k=1; k<n-1; k++) {
    /* Get tmat[k, ] */
    if (notMissing[k]) {
      xchisq = t2*logor2[k];
      for (i=0; i<J; i++) {
        vd1[i] = t2*dnchisq(xchisq, 1, ncpVecJ[i], 0);
      }
    } else {
      for (i=0; i<J; i++) vd1[i] = 1.0;
    }

    lrk  = lr[k];
    bvec = bkh_all[k];
    for (i=0; i<J; i++) {
      d1     = lrk-mu0_lrs[i];
      d2     = d1*d1;
      MAT[i] = bvec[i]*(tmp1/(tmp2*pow(tmp0+d2/t1, tmp3))*vd1[i]);
    }

    /* pzzks[,,k-1] <- tvec[k-1]*(hatp0*(akh_all[k-1,]%*%mat)) */

    d1  = tvec[k-1];
    pd = akh_all[k-1]; 
    for (i=0; i<J; i++) {
      d2 = pd[i];
      for (j=0; j<J; j++) {
        ret[i][j] += d1*hatp0[i][j]*d2*MAT[j];
      }
    }
  }


  avec = akh_all[n-2];
  for (i=0; i<J; i++) {
    d1 = avec[i];
    for (j=0; j<J; j++) {
      ret[i][j] += hatp0[i][j]*d1*vec[j];
    }
  }


  free(vd1);
  free(MAT);

} /* END: get_pzzks */

static double sum_pzzks(pzzks, n11, n12, n21, n22)
double **pzzks;
int n11, n12, n21, n22;
{
  int i, j;
  double sum=0.0;

  for (i=n11; i<=n12; i++) {
    for (j=n21; j<=n22; j++) sum += pzzks[i][j];
    
  }

  return(sum);

} /* END: sum_pzzks */

/* To sum subsets from matrix: pzzks[(n11:n12)[-(m11:m12)], n21:n22] */
static double sum_pzzks_subset1(pzzks, n11, n12, m11, m12, n21, n22)
double **pzzks;
int n11, n12, m11, m12, n21, n22;
{
  int i, j, ind=0;
  double sum=0.0;

  for (i=n11; i<=n12; i++) {
    if ((ind < m11) || (ind > m12)) {
      for (j=n21; j<=n22; j++) sum += pzzks[i][j];
    }
    ind++;
  }

  return(sum);

} /* END: sum_pzzks_subset1 */


/* Function to compute 
  hatp0_m  <- c(1-hatp0_m01, hatp0_m01, hatp0_m10, 1-hatp0_m10)
*/
static void get_hatp0_m(pzzks, n, J, mainJ, ret)
double **pzzks, *ret;
int n, J, mainJ;
{
  /* 
    hatp0_m01 <- sum(pzzks[1:mainJ,(mainJ+1):J,])/sum(pzzks[1:mainJ,,])
    hatp0_m10 <- sum(pzzks[(mainJ+1):J, 1:mainJ,])/sum(pzzks[(mainJ+1):J,,])
    hatp0_m   <- c(1-hatp0_m01, hatp0_m01, hatp0_m10, 1-hatp0_m10)
  */
 
  double d1, d2, d3, d4, hatp0_m01, hatp0_m10;

  d1 = sum_pzzks(pzzks, 0,     mainJ-1, mainJ, J-1);
  d2 = sum_pzzks(pzzks, 0,     mainJ-1, 0,     J-1);
  d3 = sum_pzzks(pzzks, mainJ, J-1,     0,     mainJ-1);
  d4 = sum_pzzks(pzzks, mainJ, J-1,     0,     J-1);

  hatp0_m01 = d1/d2;
  hatp0_m10 = d3/d4;

  ret[0] = 1-hatp0_m01;
  ret[1] = hatp0_m01;
  ret[2] = hatp0_m10;
  ret[3] = 1-hatp0_m10;

  return;

} /* END: get_hatp0_m */

static void get_hatp0_M(pzzks, n, J, mainJ, tauJ, ret)
double **pzzks, **ret;
int n, J, mainJ, tauJ;
{
  int i1, g, l, i2, i3, i41, i42, i51, i52, gm1, lm1, Jm1;
  double num, den, d1, d2, d3, d4, d5, d6;

  /* ret is mainJxmainJ, computed by column */

  Jm1  = J - 1;
  i1   = 1 + mainJ;
  for (g=1; g<=mainJ; g++) {
    gm1 = g - 1;
    i2  = gm1*tauJ;
    i3  = g*tauJ;
    i41 = i1 + i2;
    i42 = i1 + i3 - 1;
    for (l=1; l<=mainJ; l++) {
      lm1 = l - 1;
      i51 = i1 + lm1*tauJ;
      i52 = i1 + l*tauJ - 1;

      /* Offset array indices by 1 */
      d1  = sum_pzzks(pzzks, lm1,   lm1,   gm1,   gm1);
      d2  = sum_pzzks(pzzks, lm1,   lm1,   i41-1, i42-1);
      d3  = sum_pzzks(pzzks, i51-1, i52-1, gm1,   gm1);
      d4  = sum_pzzks(pzzks, i51-1, i52-1, i41-1, i42-1);
      d5  = sum_pzzks(pzzks, lm1,   lm1,   0,     Jm1);
      d6  = sum_pzzks(pzzks, i51-1, i52-1, 0,     Jm1);
      num = d1 + d2 + d3 + d4;
      den = d5 + d6; 
      ret[lm1][gm1] = num/den;
    }
  }

  return;

} /* END: get_hatp0_M */

static double sum_pzks1(pzks1, n11, n12, n21, n22)
double **pzks1;
int n11, n12, n21, n22;
{
  int i, j;
  double sum=0.0;

  for (i=n11; i<=n12; i++) {
    for (j=n21; j<=n22; j++) sum += pzks1[i][j];
  }

  return(sum);

} /* sum_pzks1 */

static void get_hatp0_t0(pzzks, n, J, mainJ, tauJ, pzks1, ret)
double **pzzks, **pzks1, *ret;
int n, J, mainJ, tauJ;
{
  int g, t, tmpi, i1, i2, i3, i4, mainJm1, i4m1, Jm1, i2m1, i3m1, i1ptm1;
  double *pret, dd1, dd2, dd3, denom, num, d2, d3;

  /* ret here is a vector of length mainJ*tauJ */

  pret    = ret;
  mainJm1 = mainJ - 1;
  Jm1     = J - 1;
  for (g=1; g<=mainJ; g++) {
    tmpi  = (g - 1)*tauJ;
    i1    = mainJ + tmpi;
    i2    = tmpi + 1;
    i3    = g*tauJ;
    i4    = mainJ + i3;
    i4m1  = i4 - 1; 
    i2m1  = i2 - 1;
    i3m1  = i3 - 1;
    dd1   = sum_pzks1(pzks1, 0, 0, i1, i4m1);
    dd2   = sum_pzzks(pzzks, 0,  mainJm1,   i1,   i4m1);
    dd3   = sum_pzzks_subset1(pzzks, mainJ, Jm1, i2m1, i3m1, mainJ+i2m1, i4m1);
    denom = dd1 + dd2 + dd3;

    for (t=1; t<=tauJ; t++) {
      i1ptm1  = i1 + t - 1;
      d2      = sum_pzzks(pzzks, 0, mainJm1, i1ptm1, i1ptm1);
      d3      = sum_pzzks_subset1(pzzks, mainJ, Jm1, i2m1, i3m1, i1ptm1, i1ptm1);
      num     = pzks1[0][i1ptm1] + d2 + d3;
      *pret++ = num/denom;
    }
  }

  return;

} /* END: get_hatp0_t0 */

static void get_r0s(pzks1, n, J, mainJ, tauJ, hatp0_t0, ret)
double **pzks1, *hatp0_t0, *ret;
int n, J, mainJ, tauJ;
{
  double d0, d1, d2, *r0, u01, u02;
  int i, j, i21, i22, ii, k;

  r0  = dVec_alloc(mainJ, 0, 0.0);
  d0  = sum_pzks1(pzks1, 0, 0, 0, J-1);

  /* Compute r0 
    r0 <- sapply(1:mainJ, function(j) sum(pzks1[1,c(j,(mainJ+1+tauJ*(j-1)):(mainJ+tauJ*j))])/sum(pzks1[1,]))
  */
  for (j=1; j<=mainJ; j++) {
    d1      = pzks1[0][j-1];
    i21     = mainJ+1+tauJ*(j-1);
    i22     = mainJ+tauJ*j; 
    d2      = sum_pzks1(pzks1, 0, 0, i21-1, i22-1);
    r0[j-1] = (d1 + d2)/d0; 
  }

  /* u0 <- c(sum(pzks1[1,1:mainJ]), sum(pzks1[1,(mainJ+1):J]))/sum(pzks1[1,]) */
  u01 = sum_pzks1(pzks1, 0, 0, 0, mainJ-1)/d0;
  u02 = sum_pzks1(pzks1, 0, 0, mainJ, J-1)/d0;

  /* 
   r0s <- c(r0*u0[1], c(sapply(1:mainJ, function(j) sapply(1:tauJ, function(t) r0[j]*u0[2]*hatp0_t0[(j-1)*tauJ+t]))))
  */

  for (j=0; j<mainJ; j++) ret[j] = r0[j]*u01;
  ii = mainJ;
  for (i=0; i<mainJ; i++) {
    d1 = r0[i]*u02;
    k  = i*tauJ;
    for (j=0; j<tauJ; j++) {
      ret[ii] = d1*hatp0_t0[k+j];
      ii++;
    }
  }


  free(r0);

  return;

} /* END: get_r0s */

/*
static void get_hatp0_mgtype(pzzks, n, J, mainJ, tauJ, r0s, ret)
double **pzzks, *r0s, *ret;
int n, J, mainJ, tauJ;
{
  double d2, d3, den, *pret;
  int i, j, i21, i22, ii;

  pret = ret;
  for (j=1; j<=mainJ; j++) {
    i21 = 1+mainJ+tauJ*(j-1);
    i22 = mainJ+tauJ*j;
    d3  = sum_pzzks(pzzks, 0, J-1, i21-1, i22-1);
    d2  = 0.0;
    for (i=i21-1; i<i22; i++) d2 += r0s[i];
    den = d2 + d3;
    ii  = i21-1;
    for (i=0; i<tauJ; i++) {
      d3      = sum_pzzks(pzzks, 0, J-1, ii, ii);
      *pret++ = (r0s[ii] + d3)/den;
      ii++;
    }
  }

  return;

} 
*/

static void init_hatp0(hatp0, mainJ, hatp0_M, hatp0_m, J)
double **hatp0, **hatp0_M, *hatp0_m;
int mainJ, J;
{
  /* copy first mainJ elements
  for(i in 1:mainJ){
    hatp0[i,] <- c(hatp0_M[i,]*hatp0_m[1], hatp01[i,] )
  }
  */

  int i, j;
  double d;

  /* init to 0 */
  for (i=0; i<J; i++) {
    for (j=0; j<J; j++) hatp0[i][j] = 0.0;
  }

  d = hatp0_m[0];
  for (i=0; i<mainJ; i++) {
    for (j=0; j<mainJ; j++) hatp0[i][j] = hatp0_M[i][j]*d;
  }

} /* END: init_hatp0 */

static void update_hatp0_01(hatp0, mainJ, tauJ, hatp0_M, hatp0_m, hatp0_t0)
double **hatp0, **hatp0_M, *hatp0_m, *hatp0_t0;
int mainJ, tauJ;
{
  int h, hh, i1, i2, i, hm1, hhm1;
  double d2, tmp, d3;

  d2 = hatp0_m[1];
  for (h=1; h<=mainJ; h++) {
    hm1 = h - 1;
    for (hh=1; hh<=mainJ; hh++) {
      hhm1 = hh - 1;
      d3   = hatp0_M[hm1][hhm1];
      i1   = hhm1*tauJ + 1 - 1;
      i2   = hh*tauJ - 1;
      for (i=i1; i<=i2; i++) {
        tmp = d3*d2*hatp0_t0[i];
        hatp0[hm1][mainJ + i] = tmp; 
      }
    }
  }

} /* END: update_hatp0_01 */

static void update_hatp0_11(hatp0, mainJ, tauJ, hatp0_M, hatp0_m, hatp0_t0)
double **hatp0, **hatp0_M, *hatp0_m, *hatp0_t0;
int mainJ, tauJ;
{
  int h, hh, i1, i2, i, j, hm1, hhm1, mtJ, row, col;
  double d4, d3, d2, **hatp_11s, tmp;

  /* 
  pt <- diag(tauJ)
  for(h in 1:mainJ){
    hatp_11s <- matrix(0,tauJ,mainJ*tauJ)
    hatp_11s[,((h-1)*tauJ+1):(h*tauJ)] <- hatp0_M[h,h]*hatp0_m[4]*pt
    for(hh in c(1:mainJ)[-h]){
      hatp_11s[,((hh-1)*tauJ+1):(hh*tauJ)] <- hatp0_M[h,hh]*hatp0_m[4]*matrix(c(hatp0_t0[((hh-1)*tauJ+1):(hh*tauJ)]), byrow=T, nrow=tauJ, ncol=tauJ)
    }
    hatp11 <- rbind(hatp11, hatp_11s)
  }
  hatp0[(mainJ+1):J,(mainJ+1):J] <- hatp11
  */

  mtJ      = mainJ*tauJ;
  hatp_11s = dMat_alloc(tauJ, mtJ, 0, 0.0);
  d4       = hatp0_m[3];
  row      = mainJ;

  for (h=1; h<=mainJ; h++) {
    /* init hatp_11s to 0 */
    for (j=0; j<tauJ; j++) {
      for (i=0; i<mtJ; i++) hatp_11s[j][i] = 0.0;
    }

    hm1 = h - 1;
    d2  = hatp0_M[hm1][hm1]*d4;
    i1   = hm1*tauJ + 1 - 1;
    i2   = h*tauJ - 1;   
    j    = 0;
    for (i=i1; i<=i2; i++) {
      hatp_11s[j][i] = d2;
      j++;
    }

    for (hh=1; hh<=mainJ; hh++) {  
      if (hh != h) {
        hhm1 = hh - 1;
        i1   = hhm1*tauJ+1 - 1;
        i2   = hh*tauJ - 1;   
        d3   = hatp0_M[hm1][hhm1];
        d2   = d3*d4;
        for (j=0; j<tauJ; j++) {
          for (i=i1; i<=i2; i++) {
            tmp = d2*hatp0_t0[i];
            hatp_11s[j][i] = tmp; 
          }
        }
      }
    }

    /* copy hatp_11s to hatp0 */
    for (i=0; i<tauJ; i++) {
      col = mainJ;
      for (j=0; j<mtJ; j++) {
        hatp0[row][col] = hatp_11s[i][j];
        col++;  
      }
      row++; 
    }
  }

  matrix_free((void **) hatp_11s, tauJ);

} /* END: update_hatp0_11 */

static void update_hatp0_10(hatp0, mainJ, tauJ, hatp0_M, hatp0_m)
double **hatp0, **hatp0_M, *hatp0_m;
int mainJ, tauJ;
{
  /* 
   for(i in 1:mainJ){
    temp <- matrix(c(hatp0_M[i,]*hatp0_m[3]), nrow=1)
    temp1 <- matrix(rep(temp,tauJ), byrow=T, ncol=mainJ)
    hatp10 <- rbind(hatp10, temp1)
  }
  hatp0[(mainJ+1):J,1:mainJ] <- hatp10
  */

  int i, j, k, row;
  double d3, *vec, *prow;

  vec = dVec_alloc(mainJ, 0, 0.0);

  d3  = hatp0_m[2];
  row = mainJ;

  for (i=0; i<mainJ; i++) {
    prow = hatp0_M[i];
    for (j=0; j<mainJ; j++) vec[j] = prow[j]*d3;
    /* Copy same vector for each set of tauJ rows */
    for (j=0; j<tauJ; j++) {
      for (k=0; k<mainJ; k++) hatp0[row][k] = vec[k];
      row++; 
    }
  }

  free(vec);

} /* END: update_hatp0_10 */

static void update_hatp0(hatp0, mainJ, tauJ, hatp0_M, hatp0_m, hatp0_t0, J)
double **hatp0, **hatp0_M, *hatp0_m, *hatp0_t0;
int mainJ, tauJ, J;
{
  init_hatp0(hatp0, mainJ, hatp0_M, hatp0_m, J);
  update_hatp0_01(hatp0, mainJ, tauJ, hatp0_M, hatp0_m, hatp0_t0);
  update_hatp0_11(hatp0, mainJ, tauJ, hatp0_M, hatp0_m, hatp0_t0);
  update_hatp0_10(hatp0, mainJ, tauJ, hatp0_M, hatp0_m);

} /* END: update_hatp0 */

static void copyMatIntoVecByCol(mat, nr, nc, ret)
double **mat, *ret;
int nr, nc;
{
  int i, j;
  double *p;

  p = ret;
  for (i=0; i<nc; i++) {
    for (j=0; j<nr; j++) {
      *p = mat[j][i];
      p++;
    }
  }

} /* END: copyMatIntoVecByCol */

void locFunc_hatp0(pn, pJ, pmainJ, ptauJ, pt1, pt2, tvec, ncpVecJ,
                   ptmp0, ptmp1, ptmp2, ptmp3, logor2, lr, notMissing,
                   akh_allVec, bkh_allVec, mu0_lrs, vec, pzks1Vec,
                   ret_code, ret_r0s, hatp0Vec, hatp0_Mvec, hatp0_m, hatp0_t0)
int *pn, *pJ, *pmainJ, *ptauJ, *notMissing, *ret_code;
double *pt1, *pt2, *tvec, *ncpVecJ, *ptmp0, *ptmp1, *ptmp2, *ptmp3;
double *logor2, *lr, *akh_allVec, *bkh_allVec, *mu0_lrs, *hatp0Vec, *vec;
double *pzks1Vec, *hatp0_Mvec, *hatp0_m, *hatp0_t0;
double *ret_r0s;
{
  int n, J, mainJ, tauJ;
  double t1, t2, tmp0, tmp1, tmp2, tmp3;
  double **akh_all, **bkh_all, **hatp0, **pzks1;
  double **pzzks, **hatp0_M;

  *ret_code = -1;
  n         = *pn;
  J         = *pJ;
  mainJ     = *pmainJ;
  tauJ      = *ptauJ;
  t1        = *pt1;
  t2        = *pt2;
  tmp0      = *ptmp0;
  tmp1      = *ptmp1;
  tmp2      = *ptmp2;
  tmp3      = *ptmp3;

  akh_all      = dMat_alloc(n, J, 0, 0.0);
  bkh_all      = dMat_alloc(n, J, 0, 0.0);
  hatp0        = dMat_alloc(J, J, 0, 0.0);
  pzks1        = dMat_alloc(n, J, 0, 0.0);
  pzzks        = dMat_alloc(J, J, 1, 0.0);
  hatp0_M      = dMat_alloc(mainJ, mainJ, 0, 0.0);

  fillMat(akh_allVec, n, J, akh_all);
  fillMat(bkh_allVec, n, J, bkh_all);
  fillMat(hatp0Vec, J, J, hatp0);
  fillMat(pzks1Vec, n, J, pzks1);

  get_pzzks(n, J, hatp0, akh_all, bkh_all, tmp0, tmp1, tmp2, tmp3, 
    t1, t2, lr, logor2, mu0_lrs, notMissing, ncpVecJ, vec, tvec, pzzks);

  get_hatp0_m(pzzks, n, J, mainJ, hatp0_m);
  get_hatp0_M(pzzks, n, J, mainJ, tauJ, hatp0_M);
  get_hatp0_t0(pzzks, n, J, mainJ, tauJ, pzks1, hatp0_t0);
  get_r0s(pzks1, n, J, mainJ, tauJ, hatp0_t0, ret_r0s); 

  /* Must be called last */
  update_hatp0(hatp0, mainJ, tauJ, hatp0_M, hatp0_m, hatp0_t0, J);

  /* Copy hatp0 into return vector by column */
  copyMatIntoVecByCol(hatp0, J, J, hatp0Vec);

  /* Copy hatp0_M into return vector by column */
  copyMatIntoVecByCol(hatp0_M, mainJ, mainJ, hatp0_Mvec);

  matrix_free((void **) akh_all, n); 
  matrix_free((void **) bkh_all, n);
  matrix_free((void **) hatp0, J);
  matrix_free((void **) pzks1, n);
  matrix_free((void **) pzzks, J);
  matrix_free((void **) hatp0_M, mainJ);

  *ret_code = 0;
  return;

} /* END: locFunc_hatp0 */

static void rowVecMat(vec, mat, nr, nc, ret)
double *vec, **mat, *ret;
int nr, nc;
{
  int i, j;
  double sum;

  for (i=0; i<nc; i++) {
    sum = 0.0;
    for (j=0; j<nr; j++) {
      sum += vec[j]*mat[j][i];
    }
    ret[i] = sum;
  }

  return;

} /* END: rowVecMat */

static void matColVec(vec, mat, nr, nc, ret)
double *vec, **mat, *ret;
int nr, nc;
{
  int i, j;
  double sum;

  for (i=0; i<nr; i++) {
    sum = 0.0;
    for (j=0; j<nc; j++) {
      sum += vec[j]*mat[i][j];
    }
    ret[i] = sum;
  }

  return;

} /* END: matColVec */


void akh_ca(pn, pJ, pmainJ, ptauJ, 
                   ppt5v0, ppt5v0POWpt5v0, pv0Plus1Over2, pgamma_v0Over2,
                   ptmp1, ptmp2, ptmp3, tmp4, logor2, lr, notMissing,
                   mu0_lrs, a1h, hatp0Vec, pc1, 
                   ret_code, ret_akh_all, ret_ca_all)
int *pn, *pJ, *pmainJ, *ptauJ, *notMissing, *ret_code;
double *ppt5v0, *ppt5v0POWpt5v0, *pv0Plus1Over2, *pgamma_v0Over2;
double *ptmp1, *ptmp2, *ptmp3, *tmp4;
double *logor2, *lr, *a1h, *mu0_lrs, *hatp0Vec, *pc1;
double *ret_akh_all, *ret_ca_all;
{
  int n, J, k, ii, j, notMiss;
  double tmp1, tmp2, tmp3, c1, tmp5, cka, lrk, logor2k, sum_akd;
  double **hatp0, pt5v0, v0Plus1Over2, gamma_v0Over2;
  double *a1h_hatp0, *akd, d1, d2;

  *ret_code     = -1;
  n             = *pn;
  J             = *pJ;
  tmp1          = *ptmp1;
  tmp2          = *ptmp2;
  tmp3          = *ptmp3;
  c1            = *pc1;
  pt5v0         = *ppt5v0;
  v0Plus1Over2  = *pv0Plus1Over2;
  gamma_v0Over2 = *pgamma_v0Over2;

  /* akh_all.new[1,] <- a1h, ca_all.new[1]   <- c1 */
  ret_ca_all[0] = c1;
  for (k=0; k<J; k++) ret_akh_all[k] = a1h[k];

  hatp0     = dMat_alloc(J, J, 0, 0.0);
  fillMat(hatp0Vec, J, J, hatp0);
  a1h_hatp0 = dVec_alloc(J, 0, 0.0);
  akd       = dVec_alloc(J, 0, 0.0);

  /* Initialize akd with a1h */
  for (j=0; j<J; j++) akd[j] = a1h[j];

  ii = J;
  for (k=1; k<n; k++) {
    notMiss = notMissing[k];
    lrk     = lr[k];
    logor2k = logor2[k];
    sum_akd = 0.0;

    /* Compute a1h*hatp0 */
    rowVecMat(akd, hatp0, J, J, a1h_hatp0);

    for (j=0; j<J; j++) {
      d1   = lrk-mu0_lrs[j];
      d1   = d1*d1;
      d2   = pow(pt5v0+(d1/tmp2), v0Plus1Over2);
      tmp5 = tmp1/(gamma_v0Over2*d2);
      
      if (!notMiss) {
        d1 = a1h_hatp0[j]*tmp5;
      } else { 
        d1 = a1h_hatp0[j]*(tmp5*(tmp3*dnchisq(logor2k*tmp3, 1, tmp4[j], 0)));
      }
      sum_akd += d1;
      akd[j]  = d1;
    }
    cka             = 1.0/sum_akd;
    ret_ca_all[k]   = cka;
    for  (j=0; j<J; j++) {
      d1              = cka*akd[j];
      akd[j]          = d1;
      ret_akh_all[ii] = d1;
      ii++;
    }
  }

  matrix_free((void **) hatp0, J);
  free(a1h_hatp0);
  free(akd);

  *ret_code = 0;
  return;

} /* END: akh_ca */

void C_hes_logL1(pn, pJ, parm, mu0_lrs, mu0_lors,
               ppi, logor2, lr, notMissing, r0s, hatp0Vec,  
               ret_code, ret_loglike)
int *pn, *pJ, *notMissing, *ret_code;
double *parm, *mu0_lrs, *mu0_lors, *ppi, *logor2, *lr, *r0s, *hatp0Vec; 
double *ret_loglike;
{
  int n, J, k, j, notMiss;
  double lrk, logor2k, sum_akd, dpi, lsigma20_lr, lsigma20_lor, v0, sigma20_lr;
  double **hatp0, p5v0, p5v0p5v0, v01, gv01, gv0, tmp1, tmp2, tmp3, cka, *a1d;
  double *a1h_hatp0, sum, akd, d1, d2, d3, *ncpVec, logL_Y, c1, sum_log_ca_all, dtmp;

#if DEBUG
Rprintf("BEGIN: C_hes_logL1\n");
#endif

  *ret_code     = -1;
  n             = *pn;
  J             = *pJ;
  dpi           = *ppi;
  lsigma20_lr   = parm[3-1];
  lsigma20_lor  = parm[4-1];
  v0            = parm[5-1];
  sigma20_lr    = exp(lsigma20_lr);
  p5v0          = 0.5*v0;
  p5v0p5v0      = pow(p5v0, p5v0);
  v01           = (v0 + 1.0)/2.0;
  gv01          = gammafn(v01);
  gv0           = gammafn(v0/2.0);
  tmp1          = p5v0p5v0*gv01*pow(2.0*dpi*sigma20_lr, -0.5);
  tmp2          = 2.0*sigma20_lr;
  tmp3          = exp(-lsigma20_lor);
  sum_akd       = 0.0;

  ncpVec    = dVec_alloc(J, 0, 0.0);
  a1d       = dVec_alloc(J, 0, 0.0);
  a1h_hatp0 = dVec_alloc(J, 0, 0.0);
  hatp0     = dMat_alloc(J, J, 0, 0.0);
  fillMat(hatp0Vec, J, J, hatp0);

#if DEBUG
Rprintf("Compute ncp\n");
#endif

  for (k=0; k<J; k++) ncpVec[k] = mu0_lors[k]*mu0_lors[k]*tmp3*r0s[k];

#if DEBUG
Rprintf("Compute a1d\n");
#endif

  /* Compute a1d */
  notMiss = notMissing[0];
  lrk     = lr[0];
  logor2k = logor2[0];
  sum     = 0.0;
  for (k=0; k<J; k++) {
    d1 = lrk - mu0_lrs[k];
    d1 = d1*d1;
    d2 = gv0*pow(p5v0 + d1/tmp2, v01);
    d3 = (tmp1/d2)*r0s[k];  
    if (notMiss) d3 = d3*tmp3*dnchisq(logor2k*tmp3, 1.0, ncpVec[k], 0);

    a1d[k] = d3;
    sum += d3;
  }


#if DEBUG
Rprintf("Re-define ncp\n");
#endif

  c1             = 1.0/sum;
  sum_log_ca_all = log(c1);
  for (k=0; k<J; k++) {
    a1d[k]    = c1*a1d[k]; 
    ncpVec[k] = mu0_lors[k]*mu0_lors[k]*tmp3;
  }

#if DEBUG
dtmp = ncpVec[0];
for (k=1; k<J; k++) {
  d1 = ncpVec[k];
  if (d1 > dtmp) dtmp = d1;
}
Rprintf("Max ncp = %g\n", dtmp);
Rprintf("Main loop\n");
#endif


  for (k=1; k<n; k++) {
    notMiss = notMissing[k];
    logor2k = logor2[k];
    lrk     = lr[k];
    dtmp    = logor2k*tmp3;
    
    rowVecMat(a1d, hatp0, J, J, a1h_hatp0);
    sum_akd = 0.0;
    for (j=0; j<J; j++) {
      d1  = lrk - mu0_lrs[j];
      d1  = d1*d1;
      d2  = gv0*pow(p5v0 + d1/tmp2, v01);
      akd = a1h_hatp0[j]*(tmp1/d2);  
      if (notMiss) akd = akd*tmp3*dnchisq(dtmp, 1.0, ncpVec[j], 0);
      sum_akd += akd;
      a1d[j] = akd;
    }
    cka = 1.0/sum_akd;

    sum_log_ca_all += log(cka);
    sum_akd = 0.0;
    for (j=0; j<J; j++) {
      d1       = a1d[j]*cka;
      a1d[j]   = d1;
      sum_akd += d1;
    }
  }

#if DEBUG
Rprintf("end main loop\n");
#endif


  logL_Y = -sum_log_ca_all + log(sum_akd);
  *ret_loglike = logL_Y;
  
  free(ncpVec);
  free(a1d);
  free(a1h_hatp0);
  matrix_free((void **) hatp0, J);

  *ret_code = 0;

#if DEBUG
Rprintf("END: C_hes_logL1\n");
#endif


  return;

} /* END: C_hes_logL1 */

void bkh_cb(pn, pJ, ptmp0, ptmp1, ptmp2, ptmp3, ptmp4, ptmp5, tmp6,
    logor2, lr, notMissing, mu0_lrs, hatp0Vec,  
                   ret_code, ret_bkh_all, ret_cb_all)
int *pn, *pJ, *notMissing, *ret_code;
double *ptmp0, *ptmp1, *ptmp2, *ptmp3, *ptmp4, *ptmp5, *tmp6;
double *logor2, *lr, *mu0_lrs, *hatp0Vec;
double *ret_bkh_all, *ret_cb_all;
{
  int n, J, k, ii, kp1, j;
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, lrk1, logor2k1;
  double **hatp0, d1, d2, *bnh, *vec, *bkd, ckb;
  
  *ret_code = -1;
  n         = *pn;
  J         = *pJ;
  tmp0      = *ptmp0;
  tmp1      = *ptmp1;
  tmp2      = *ptmp2;
  tmp3      = *ptmp3;
  tmp4      = *ptmp4;
  tmp5      = *ptmp5;

  hatp0     = dMat_alloc(J, J, 0, 0.0);
  bnh       = dVec_alloc(J, 0, 0.0);
  vec       = dVec_alloc(J, 0, 0.0);
  bkd       = dVec_alloc(J, 0, 0.0);
  fillMat(hatp0Vec, J, J, hatp0);

  /* We are working backwards, so get results for final index, row */
  d1 = 1.0/((double) J);
  ret_cb_all[n-1] = d1;
  ii = (n-1)*J;
  for (k=0; k<J; k++) {
    ret_bkh_all[ii] = d1;
    ii++;
    bnh[k] = d1;
  }

  for (k=n-2; k>-1; k--) {
    kp1      = k + 1;
    logor2k1 = logor2[kp1];
    lrk1     = lr[kp1];

    if (!notMissing[kp1]) {
      for (j=0; j<J; j++) {
        d1     = lrk1 - mu0_lrs[j];
        d2     = pow(tmp0 + d1*d1/tmp3, tmp4);
        vec[j] = tmp1/(tmp2*d2)*bnh[j];
      }
    } else {
      for (j=0; j<J; j++) {
        d1     = lrk1 - mu0_lrs[j];
        d2     = pow(tmp0 + d1*d1/tmp3, tmp4);
        vec[j] = tmp1/(tmp2*d2)*(tmp5*dnchisq(logor2k1*tmp5, 1.0, tmp6[j], 0))*bnh[j];
      }
    }
    matColVec(vec, hatp0, J, J, bkd);
    d1  = 0.0;
    for (j=0; j<J; j++) d1 += bkd[j];
    ckb = 1.0/d1;

    ret_cb_all[k] = ckb;
    ii            = k*J;   
    for (j=0; j<J; j++) {
      d1 = ckb*bkd[j];
      ret_bkh_all[ii] = d1;
      ii++;
      bnh[j] = d1;
    }
  } 

  free(bkd);
  free(bnh);
  free(vec);
  matrix_free((void **) hatp0, J);

  *ret_code = 0;
  return;
  
} /* END: bkh_cb */



static double ** Lmatrix(int n)
{
    int   i;
    double **m;

    m = (double **) malloc(n*sizeof(double *));
    for (i = 0; i < n; i++)
	m[i] = (double *) malloc((i + 1)*sizeof(double));
    return m;
}





static double logBase2(d)
double d;
{
  double ret;
  ret = log(d)/LOG2;
  return(ret);

} /* END: logBase2 */

static void mu0_lr_func(nj, mainJ, lpsi0, lpr0, lsz0, ctzs, ctms, ret)
int nj, mainJ;
double lpsi0, lpr0, lsz0, *ctzs, *ctms, *ret;
{
  int i;
  double pr0, sz0, tmp1;
  
  pr0 = 1.0/(1.0 + exp(lpr0));
  sz0 = 1.0/(1.0 + exp(lsz0));
  for (i=0; i<mainJ; i++) ret[i] = logBase2(2.0+pr0*(ctzs[i]-2.0)) - lpsi0;
  for (i=mainJ; i<nj; i++) {
    tmp1   = ctms[i - mainJ];
    ret[i] = logBase2(2.0+(tmp1-2.0)*pr0+(ctzs[i]-tmp1)*sz0*pr0) - lpsi0;
  } 
  
  return;

} /* END: mu0_lr_func */    

static void mu0_lor_func(nj, mainJ, lpr0, lsz0, mz, pz, mz_sub, pz_sub, ret)
int nj, mainJ;
double lpr0, lsz0, *mz, *pz, *mz_sub, *pz_sub, *ret;
{
  int i, j;
  double pr0, sz0, tmp1, tmp2, arg1, arg2, oneMinusPr0, szpr;
  
  pr0         = 1.0/(1.0 + exp(lpr0));
  sz0         = 1.0/(1.0 + exp(lsz0));
  oneMinusPr0 = 1.0 - pr0;
  szpr        = sz0*pr0;

  for (i=0; i<mainJ; i++) {
    arg1 = mz[i]*pr0 + oneMinusPr0;
    if (arg1 < LOG0ARG) {
      tmp1 = LOG0; 
    } else {
      tmp1 = log(arg1);
    }
    arg2 = pz[i]*pr0 + oneMinusPr0;
    if (arg2 < LOG0ARG) {
      tmp2 = LOG0; 
    } else {
      tmp2 = log(arg2);
    }
    ret[i] = tmp1 - tmp2;
  }

  for (i=mainJ; i<nj; i++) {
    j    = i - mainJ;
    tmp1 = mz_sub[j];
    arg1 = 1.0 + (mz[i] - tmp1)*szpr + (tmp1 - 1.0)*pr0;   
    if (arg1 < LOG0ARG) {
      tmp1 = LOG0; 
    } else {
      tmp1 = log(arg1);
    }
    tmp2 = pz_sub[j];
    arg2 = 1.0 + (pz[i] - tmp2)*szpr + (tmp2 - 1.0)*pr0;
    if (arg2 < LOG0ARG) {
      tmp2 = LOG0; 
    } else {
      tmp2 = log(arg2);
    }
    ret[i] = tmp1 - tmp2;
  } 
  
  return;

} /* END: mu0_lor_func */   

static double negLoglike(parm, nparm, p5L, p6L, p6U, n, J, mainJ, lr, logor2,
    notMissing, pzks1, ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors)
int nparm, n, J, *notMissing, mainJ;
double *parm, p5L, p6L, p6U, *lr, *logor2, **pzks1, *ctzs, *ctms, *mz, *pz, *mz_sub, *pz_sub;
double *mu0_lrs, *mu0_lors;
{
  int i, nomiss, j;
  double p1, p2, p3, p4, p5, p6, p62, sigma20_lr, ep5, d1, d2, d3, mat, ret=0.0, p6p1, lri;
  double *vec, pzks1ij, logor2i, darg;

  p1         = parm[0];
  p2         = parm[1];
  p3         = parm[2];
  p4         = parm[3];
  p5         = parm[4];
  p6         = parm[5];
  sigma20_lr = exp(p4);
  
  /* Reparameterize */
  p62  = p6*p6;
  p6   = p6L + (p6U - p6L)*p62/(1.0 + p62);

  p5   = p5L + p5*p5;
  ep5  = exp(-p5);
  if (ep5 > MAX_EP5) error("ERROR: set bounds on parm[5]");
  p6p1 = p6 + 1.0;

  mu0_lr_func(J, mainJ, p1, p2, p3, ctzs, ctms, mu0_lrs);
  mu0_lor_func(J, mainJ, p2, p3, mz, pz, mz_sub, pz_sub, mu0_lors);


  /* Define the ncp and store in mu0_lors */
  for (i=0; i<J; i++) {
    d1          = mu0_lors[i];
    mu0_lors[i] = d1*d1*ep5;
  }

  d1 = -0.5*log(TWOPI*sigma20_lr) + log(gammafn(p6p1/2.0)) - log(gammafn(p6/2.0)) + 0.5*p6*log(p6/2.0);
  d2 = 0.5*p6p1;
  for (i=0; i<n; i++) {
    lri     = lr[i];
    logor2i = logor2[i];
    vec     = pzks1[i];
    nomiss  = notMissing[i];
    darg    = logor2i*ep5;
    for (j=0; j<J; j++) {
      d3      = lri - mu0_lrs[j];
      pzks1ij = vec[j];
      mat     = (d1 - d2*log(0.5*(p6 + d3*d3/sigma20_lr)))*pzks1ij;
      if (nomiss) mat = mat - pzks1ij*(p5 - dnchisq(darg, 1.0, mu0_lors[j], 1));
      ret += mat;
    }
  }

  ret = -ret;

  return(ret);

} /* END: negLoglike */

static void gradient(eta, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors, ret)
double *eta, p5L, p6L, p6U, *lr, *logor2, **pzks1, *ctzs, *ctms, *mz, *pz, *mz_sub, *pz_sub, *mu0_lrs, *mu0_lors;
int nparm, nobs, J, mainJ, *notMissing;
double *ret;
{
  int i;
  double fxplush, fxminush, h, save, h2, *ptr, *ptrret;

  /* default step size */
  h = 1e-3;
  h2 = 2.0*h;

  for (i=0, ptr=eta, ptrret=ret; i<nparm; i++, ptr++, ptrret++) {
    save = *ptr;  

    *ptr = save + h;  
    /*fxplush = negloglike(eta, op);*/
    fxplush = negLoglike(eta, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors);

    *ptr = save - h;  
    /*fxminush = negloglike(eta, op);*/
    fxminush =  negLoglike(eta, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors);
    *ptr = save;

    *ptrret = (fxplush - fxminush)/h2;
  }

} /* END: gradient */


static void myvmmin(int n0, double *b, double *Fmin, 
      int maxit, int trace, double abstol, double reltol, int nREPORT, 
      int *fncount, int *grcount, int *fail,
      double p5L, double p6L, double p6U, int nobs, int J, int mainJ, double *lr, double *logor2, int *notMissing, 
      double **pzks1, double *ctzs, double *ctms, double *mz, double *pz, 
      double *mz_sub, double *pz_sub, double *mu0_lrs, double *mu0_lors)
{
    int accpoint, enough;
    double *g, *t, *X, *c, **B;
    int   count, funcount, gradcount;
    double f, gradproj;
    int   i, j, ilast, iter = 0;
    double s, steplength;
    double D1, D2;
    int   n, *l, nparm=n0;

    if (maxit <= 0) {
	*fail = 0;
	/* *Fmin = negloglike(b, op);*/
       *Fmin =  negLoglike(b, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors);
	*fncount = *grcount = 0;
	return;
    }

    l = (int *) R_alloc(n0, sizeof(int));
    n = 0;
    /*for (i = 0; i < n0; i++) if (mask[i]) l[n++] = i;*/
    for (i = 0; i < n0; i++) l[n++] = i;
    g = (double *) R_alloc(n0, sizeof(double));
    t = (double *) R_alloc(n, sizeof(double));
    X = (double *) R_alloc(n, sizeof(double));
    c = (double *) R_alloc(n, sizeof(double));
    B = Lmatrix(n);
    /*f = negloglike(b, op);*/
    f = negLoglike(b, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors);
 
    if (!isfinite(f))
	Rprintf("initial value in 'vmmin' is not finite");
    if (trace) Rprintf("initial  value %f \n", f);
    *Fmin = f;
    funcount = gradcount = 1;
    /*gradient(b, op, g);*/
    gradient(b, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors, g);
    iter++;
    ilast = gradcount;

    do {
	if (ilast == gradcount) {
	    for (i = 0; i < n; i++) {
		for (j = 0; j < i; j++) B[i][j] = 0.0;
		B[i][i] = 1.0;
	    }
	}
	for (i = 0; i < n; i++) {
	    X[i] = b[l[i]];
	    c[i] = g[l[i]];
	}
	gradproj = 0.0;
	for (i = 0; i < n; i++) {
	    s = 0.0;
	    for (j = 0; j <= i; j++) s -= B[i][j] * g[l[j]];
	    for (j = i + 1; j < n; j++) s -= B[j][i] * g[l[j]];
	    t[i] = s;
	    gradproj += s * g[l[i]];
	}

	if (gradproj < 0.0) {	/* search direction is downhill */
	    steplength = 1.0;
	    accpoint = 0;
	    do {

		count = 0;
		for (i = 0; i < n; i++) {
		    b[l[i]] = X[i] + steplength * t[i];
		    if (reltest + X[i] == reltest + b[l[i]]) /* no change */
			count++;
		}
		if (count < n) {
		    /*f = negloglike(b, op);*/
                  f = negLoglike(b, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1, 
                              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors);

		    funcount++;
		    accpoint = isfinite(f) &&
			(f <= *Fmin + gradproj * steplength * acctol);
		    if (!accpoint) {
			steplength *= stepredn;
		    }
		}

	    } while (!(count == n || accpoint));
	    enough = (f > abstol) && 
		fabs(f - *Fmin) > reltol * (fabs(*Fmin) + reltol);

	    /* stop if value if small or if relative change is low */
	    if (!enough) {
		count = n;
		*Fmin = f;
	    }

	    if (count < n) {/* making progress */
		*Fmin = f;
		/*gradient(b, op, g);*/
              gradient(b, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors, g);

		gradcount++;
		iter++;
		D1 = 0.0;
		for (i = 0; i < n; i++) {
		    t[i] = steplength * t[i];
		    c[i] = g[l[i]] - c[i];
		    D1 += t[i] * c[i];
		}
		if (D1 > 0) {
		    D2 = 0.0;
		    for (i = 0; i < n; i++) {
			s = 0.0;
			for (j = 0; j <= i; j++)
			    s += B[i][j] * c[j];
			for (j = i + 1; j < n; j++)
			    s += B[j][i] * c[j];
			X[i] = s;
			D2 += s * c[i];
		    }
		    D2 = 1.0 + D2 / D1;
		    for (i = 0; i < n; i++) {
			for (j = 0; j <= i; j++)
			    B[i][j] += (D2 * t[i] * t[j]
					- X[i] * t[j] - t[i] * X[j]) / D1;
		    }
		} else {	/* D1 < 0 */
		    ilast = gradcount;
		}
	    } else {	/* no progress */
		if (ilast < gradcount) {
		    count = 0;
		    ilast = gradcount;
		}
	    }
	} else {		/* uphill search */
	    count = 0;
	    if (ilast == gradcount) count = n;
	    else ilast = gradcount;
	    /* Resets unless has just been reset */
	}
	if (trace && (iter % nREPORT == 0))
	    Rprintf("iter%4d value %f\n", iter, f);
	if (iter >= maxit) break;
	if (gradcount - ilast > 2 * n)
	    ilast = gradcount;	/* periodic restart */

    } while (count != n || ilast != gradcount);
    if (trace) {
	Rprintf("final  value %f \n", *Fmin);
	if (iter < maxit) Rprintf("converged\n");
	else Rprintf("stopped after %i iterations\n", iter);
    }
    *fail = (iter < maxit) ? 0 : 1;
    *fncount = funcount;
    *grcount = gradcount;
    /*
    if (op->debug) {
      Rprintf("fncount = %d   grcount = %d \n", funcount, gradcount);
    }
    */
}


/* Main function for optimizing loglike */
void C_myoptimC(parm, pnparm, pp5L, pp6L, pp6U, pnobs, pJ, pmainJ, lr, logor2, notMissing,
             pzks1Vec, ctzs, ctms, mz, pz, mz_sub, pz_sub, preltol,
             ret_code)
double *parm, *pp5L, *pp6L, *pp6U, *lr, *logor2, *pzks1Vec, *ctzs, *ctms;
double *mz, *pz, *mz_sub, *pz_sub, *preltol;
int *pnparm, *pJ,*pnobs, *pmainJ, *notMissing, *ret_code; 
{
  int nparm, nobs, J, mainJ;
  int nREPORT, fail, maxit, trace, fncount, grcount;
  double p5L, p6L, p6U, abstol, reltol, Fmin, **pzks1, *mu0_lrs, *mu0_lors;

  *ret_code = -1;
  nparm     = *pnparm;
  nobs      = *pnobs;
  p5L       = *pp5L;
  p6L       = *pp6L;
  p6U       = *pp6U;
  J         = *pJ;
  mainJ     = *pmainJ;
  reltol    = *preltol;

  nREPORT   = 1;
  abstol    = -1.0e100;
  fail      = 1;
  fncount   = 0;
  grcount   = 0;
  maxit     = 100;
  trace     = 0;

  mu0_lrs  = dVec_alloc(J, 0, 0.0);
  mu0_lors = dVec_alloc(J, 0, 0.0);
  pzks1    = dMat_alloc(nobs, J, 0, 0.0);
  fillMat(pzks1Vec, nobs, J, pzks1);

  myvmmin(nparm, parm, &Fmin, maxit, trace, abstol, reltol, nREPORT, 
      &fncount, &grcount, &fail,
      p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, 
      pzks1, ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors);

  free(mu0_lrs);
  free(mu0_lors);
  matrix_free((void **) pzks1, nobs);

  *ret_code = fail;

  return;

} /* END: C_myoptimC */

/*
mu0_lr_func(nj, mainJ, lpsi0, lpr0, lsz0, ctzs, ctms, ret)
mu0_lor_func(nj, mainJ, lpr0, lsz0, mz, pz, mz_sub, pz_sub, ret)

hes_logL1(pn, pJ, parm, mu0_lrs, mu0_lors,
               ppi, logor2, lr, notMissing, r0s, hatp0Vec,  
               ret_code, ret_loglike)


static void optim_hessian(eta, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors, ret)
double *eta, p5L, p6L, p6U, *lr, *logor2, **pzks1, *ctzs, *ctms, *mz, *pz, *mz_sub, *pz_sub, *mu0_lrs, *mu0_lors;
int nparm, nobs, J, mainJ, *notMissing;
double *ret;
{
  int i, nparms, j, ni;
  double *gxplush, *gxminush, h, save, h2, *ptr;

  gxplush  = dVec_alloc(nparms, 0, 0.0);
  gxminush = dVec_alloc(nparms, 0, 0.0);

  h = 1e-3;
  h2 = 2.0*h;
  for (i=0, ptr=eta; i<nparms; i++, ptr++) {
    save = *ptr;  

    *ptr = save + h;  
    gradient(eta, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors, gxplush);

    *ptr = save - h;  
    gradient(eta, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors, gxminush);

    *ptr = save;

    ni = i*nparms;
    for (j=i; j<nparms; j++) ret[ni + j] = (gxplush[j] - gxminush[j])/h2;
    for (j=i+1; j<nparms; j++) ret[j*nparms + i] = ret[ni + j];
  }

  free(gxplush);
  free(gxminush);

}  END: optim_hessian */

/**************************************************************************
***************************************************************************
***************************************************************************
***************************************************************************
***************************************************************************/

struct mystr {
  int n;
  int mainJ;
  int tauJ;
  int J;
  int maxiter;
  int print;
  char *notMissing;
  int nparms;
  int converged;
  int niter;
  int debug;
  int *maxcol1, *maxcol2, *maxcol3;

  double *ctzs;
  double *ctms;
  double *mz;
  double *pz;
  double *mz_sub;
  double *pz_sub;
  double pr0;
  double lpr0;
  double sigma20_lr;
  double lsigma20_lr;
  double lpsi0;
  double sigma20_lor;
  double lsigma20_lor;
  double lsz0;
  double P5_LB;
  double V0_LB;
  double V0_UB;
  double v0;
  double loglike_eps;
  double parm_eps;
  double reltol;

  double pt5v0;
  double pt5v0POWpt5v0;
  double v0Plus1Over2;
  double gamma_pt5v0; 
  double gamma_v0Plus1Over2; /* gamma((v0+1)/2) */
  double twoPiSig20lrPOWmpt5; /* (2*pi*sigma20_lr)^-0.5 */
  double TMP1; /* pt5v0POWpt5v0*gamma((v0+1)/2)*(2*pi*sigma20_lr)^-0.5 */
  double sigma20_lor_inv; /* exp(-lsigma20_lor) */
  
  double *mu0_lrs;
  double *mu0_lors;
  double *parms;
  double *parms0;
  double *a1d;
  double sum_a1d;
  double *lr;
  double *logor2;
  double *r0s;
  double *ncpVec;
  double *akh_all;  /* nxJ matrix by rows */
  double *ca_all;
  double *log_ca_all;
  double *a1h;
  double **hatp0;
  double **hatp0_M;
  double *hatp0_m;
  double *hatp0_t0;
  double sum_akh_all_n;
  double *r0;
  double *u0;
  double ll0, ll1;
  double *cb_all;
  double cb_all_n; /* last element before log-transforming */
  double *bkh_pks;
  double **pzzks;

  double *tmpvec1_J, *tmpvec2_J, *tmpvec3_J;
  double *tmpvec1_n, *tmpvec2_n, *tmpvec3_n;
};
typedef struct mystr MYSTR;


static void getRowPtrs(rowMatAsVec, nr, nc, ret)
double *rowMatAsVec, **ret;
int nr, nc;
{
  int i;
  
  for (i=0; i<nr; i++) ret[i] = &rowMatAsVec[i*nc];

} 

static int getNotMissingVec(logor2, n, ret)
double *logor2;
int n;
char *ret;
{
  int i, nmiss=0;
  
  for (i=0; i<n; i++) {
    if (logor2[i] < DMISSTEST_LT) {
      ret[i] = 0;
      nmiss++;
    } else {
      /* not missing */
      ret[i] = 1;
    }
  }

  return(nmiss);

} 

static void updateParms(parms, mystr)
double *parms;
MYSTR *mystr;
{
  double v0;

  mystr->lpsi0               = parms[0];
  mystr->lpr0                = parms[1];
  mystr->lsz0                = parms[2];
  mystr->lsigma20_lr         = parms[3];
  mystr->lsigma20_lor        = parms[4];
  v0                         = parms[5];
  mystr->v0                  = v0;
  mystr->sigma20_lr          = exp(mystr->lsigma20_lr);
  mystr->sigma20_lor         = exp(mystr->lsigma20_lor);
  mystr->pr0                 = 1.0/(1.0 + exp(mystr->lpr0));
  mystr->pt5v0               = 0.5*v0;
  mystr->pt5v0POWpt5v0       = pow(mystr->pt5v0, mystr->pt5v0);
  mystr->v0Plus1Over2        = 0.5*(v0 + 1.0);
  mystr->gamma_pt5v0         = gammafn(mystr->pt5v0);
  mystr->gamma_v0Plus1Over2  = gammafn(mystr->v0Plus1Over2);
  mystr->twoPiSig20lrPOWmpt5 = 1.0/sqrt(TWOPI*mystr->sigma20_lr);
  mystr->sigma20_lor_inv     = 1.0/mystr->sigma20_lor;
  mystr->TMP1                = (mystr->pt5v0POWpt5v0)*(mystr->gamma_v0Plus1Over2)*(mystr->twoPiSig20lrPOWmpt5);

  mu0_lr_func(mystr->J, mystr->mainJ, mystr->lpsi0, mystr->lpr0, mystr->lsz0,
              mystr->ctzs, mystr->ctms, mystr->mu0_lrs);
  mu0_lor_func(mystr->J, mystr->mainJ, mystr->lpr0, mystr->lsz0, mystr->mz, 
               mystr->pz, mystr->mz_sub, mystr->pz_sub, mystr->mu0_lors);

} /* END: updateParms */


static void mystr_init(mystr, iargs, dargs, lr, logor2, mz, mz_sub, pz, pz_sub,
     ctzs, ctms, hatp0_M, hatp0_m, hatp0_t0, r0s, r0, u0, ret_parms, ret_mu0_lrs,
     ret_mu0_lors, ret_pzks1, ret_hatp0)
MYSTR *mystr;
int *iargs;
double *dargs, *lr, *logor2, *hatp0_M, *hatp0_m, *hatp0_t0, *r0s, *r0, *u0;
double *ret_parms, *ret_mu0_lrs, *ret_mu0_lors;
double *mz, *mz_sub, *pz, *pz_sub, *ctzs, *ctms, *ret_pzks1, *ret_hatp0;
{
  int n, J, tauJ, mainJ, debug;

  debug               = iargs[IARG_DEBUG];
  if (debug) Rprintf("Begin: mystr_init\n");

  mystr->mainJ        = iargs[IARG_mainJ];
  mystr->tauJ         = iargs[IARG_tauJ];
  mystr->J            = iargs[IARG_J];
  mystr->n            = iargs[IARG_n];
  mystr->maxiter      = iargs[IARG_maxiter];
  mystr->print        = iargs[IARG_print];
  mystr->debug        = debug;

  mystr->pr0          = dargs[DARG_pr0];
  mystr->lpr0         = dargs[DARG_lpr0];
  mystr->sigma20_lr   = dargs[DARG_sigma20_lr];
  mystr->lsigma20_lr  = dargs[DARG_lsigma20_lr];
  mystr->sigma20_lor  = dargs[DARG_sigma20_lor];
  mystr->lsigma20_lor = dargs[DARG_lsigma20_lor];
  mystr->lpsi0        = dargs[DARG_lpsi0];
  mystr->lsz0         = dargs[DARG_lsz0];
  mystr->P5_LB        = dargs[DARG_P5_LB];
  mystr->V0_LB        = dargs[DARG_V0_LB];
  mystr->V0_UB        = dargs[DARG_V0_UB];
  mystr->v0           = dargs[DARG_v0];
  mystr->loglike_eps  = dargs[DARG_loglike_eps];
  mystr->parm_eps     = dargs[DARG_parm_eps];
  mystr->reltol       = dargs[DARG_reltol];

  mystr->lr           = lr;
  mystr->logor2       = logor2;
  mystr->mz           = mz;
  mystr->mz_sub       = mz_sub;
  mystr->pz           = pz;
  mystr->pz_sub       = pz_sub;
  mystr->ctzs         = ctzs;
  mystr->ctms         = ctms;
  mystr->hatp0_m      = hatp0_m;
  mystr->hatp0_t0     = hatp0_t0;
  mystr->r0s          = r0s;
  mystr->r0           = r0;
  mystr->u0           = u0;
  mystr->parms        = ret_parms;
  mystr->mu0_lrs      = ret_mu0_lrs;
  mystr->mu0_lors     = ret_mu0_lors;

  n     = mystr->n;
  J     = mystr->J;
  tauJ  = mystr->tauJ;
  mainJ = mystr->mainJ;

  mystr->hatp0_M = ptrdVec_alloc(mainJ, 0);
  getRowPtrs(hatp0_M, mainJ, mainJ, mystr->hatp0_M);
  
  mystr->hatp0 = ptrdVec_alloc(J, 0);
  getRowPtrs(ret_hatp0, J, J, mystr->hatp0);
  update_hatp0(mystr->hatp0, mainJ, tauJ, mystr->hatp0_M, hatp0_m, hatp0_t0, J);

  mystr->a1d        = dVec_alloc(J, 0, 0.0);
  mystr->a1h        = dVec_alloc(J, 0, 0.0);
  mystr->ncpVec     = dVec_alloc(J, 0, 0.0);
  mystr->akh_all    = dVec_alloc(n*J, 0, 0.0);
  mystr->ca_all     = dVec_alloc(n, 0, 0.0);
  mystr->log_ca_all = dVec_alloc(n, 0, 0.0);
  mystr->cb_all     = dVec_alloc(n, 0, 0.0);
  mystr->bkh_pks    = ret_pzks1;
  mystr->pzzks      = dMat_alloc(J, J, 1, 0.0);

  mystr->tmpvec1_J  = dVec_alloc(J, 0, 0.0);
  mystr->tmpvec2_J  = dVec_alloc(J, 0, 0.0);
  mystr->tmpvec3_J  = dVec_alloc(J, 0, 0.0);
  mystr->tmpvec1_n  = dVec_alloc(n, 0, 0.0);
  mystr->tmpvec2_n  = dVec_alloc(n, 0, 0.0);
  mystr->tmpvec3_n  = dVec_alloc(n, 0, 0.0);

  mystr->notMissing = cVec_alloc(n, 0, 0);
  getNotMissingVec(logor2, n, mystr->notMissing);

  mystr->parms[0]   = mystr->lpsi0;
  mystr->parms[1]   = mystr->lpr0;
  mystr->parms[2]   = mystr->lsz0;
  mystr->parms[3]   = mystr->lsigma20_lr;
  mystr->parms[4]   = mystr->lsigma20_lor; 
  mystr->parms[5]   = mystr->v0; 
  updateParms(mystr->parms, mystr);
  mystr->parms0     = dVec_alloc(NPARM, 0, 0.0);
  copy_dVec(mystr->parms, mystr->parms0, NPARM);

  mystr->sum_a1d       = 0.0;
  mystr->sum_akh_all_n = 0.0;
  mystr->converged     = 0;
  mystr->niter         = 0;
  mystr->ll0           = 0.0;
  mystr->ll1           = 0.0;
  mystr->cb_all_n      = 0.0;

  if (debug) Rprintf("End: mystr_init\n");

}

static void mystr_free(mystr)
MYSTR *mystr;
{
  if (mystr->debug) Rprintf("Begin: mystr_free\n");

  free(mystr->hatp0_M);
  free(mystr->hatp0);
  free(mystr->a1d);
  free(mystr->ncpVec);
  free(mystr->akh_all);
  free(mystr->ca_all);
  free(mystr->log_ca_all);
  free(mystr->a1h);
  free(mystr->cb_all);
  free(mystr->parms0);
  free(mystr->notMissing);
  matrix_free((void **) mystr->pzzks, mystr->J);
  free(mystr->tmpvec1_J);
  free(mystr->tmpvec2_J);
  free(mystr->tmpvec3_J);
  free(mystr->tmpvec1_n);
  free(mystr->tmpvec2_n);
  free(mystr->tmpvec3_n);

  if (mystr->debug) Rprintf("End: mystr_free\n");

} 
 
static double get_a1d(mystr, a1d)
MYSTR *mystr;
double *a1d;
{
  if (mystr->debug) Rprintf("Begin: get_a1d\n");

  int k, J=mystr->J;
  double *mu0_lors=mystr->mu0_lors, *mu0_lrs=mystr->mu0_lrs, *r0s=mystr->r0s, lrk, logor2k;
  double *ncpVec=mystr->ncpVec, tmp3=mystr->sigma20_lor_inv, *lr=mystr->lr;
  double *logor2=mystr->logor2, gv0=mystr->gamma_pt5v0, v01=mystr->v0Plus1Over2;
  double tmp2, p5v0=mystr->pt5v0, sum, d1, d2, d3, tmp1;
  char *notMissing=mystr->notMissing, notMiss;

  for (k=0; k<J; k++) ncpVec[k] = mu0_lors[k]*mu0_lors[k]*tmp3;

  notMiss = notMissing[0];
  lrk     = lr[0];
  logor2k = logor2[0];
  sum     = 0.0;
  tmp2    = 2.0*mystr->sigma20_lr;
  tmp1    = mystr->TMP1;

  for (k=0; k<J; k++) {
    d1 = lrk - mu0_lrs[k];
    d1 = d1*d1;
    d2 = gv0*pow(p5v0 + d1/tmp2, v01);
    d3 = (tmp1/d2)*r0s[k];  
/*Rprintf("%d %g %g %g\n", k+1, tmp1, d2, r0s[k]);*/
    if (notMiss) d3 = d3*tmp3*dnchisq(logor2k*tmp3, 1.0, ncpVec[k], 0);

    a1d[k] = d3;
    sum += d3;
  }
  if (mystr->debug) Rprintf("End: get_a1d\n");
 

  return(sum);

} 

static void get_akh_ca(mystr)
MYSTR *mystr;
{
  if (mystr->debug) Rprintf("Begin: get_akh_ca\n");

  int n=mystr->n, J=mystr->J, k, ii, j;
  double tmp1=mystr->TMP1, tmp2, tmp3, c1, tmp5, cka, lrk, logor2k, sum_akd;
  double pt5v0, v0Plus1Over2, gamma_v0Over2, *lr=mystr->lr, *logor2=mystr->logor2;
  double *a1h_hatp0, *akd, d1, d2, *ca_all=mystr->ca_all, *a1h=mystr->a1h, *a1d=mystr->a1d;
  double **hatp0=mystr->hatp0, *mu0_lrs=mystr->mu0_lrs, *mu0_lors=mystr->mu0_lors, *tmp4;
  double *akh_all=mystr->akh_all, *pd;
  char notMiss, *notMissing=mystr->notMissing;

  tmp2          = 2.0*mystr->sigma20_lr;
  tmp3          = mystr->sigma20_lor_inv;
  c1            = 1.0/(mystr->sum_a1d);
  pt5v0         = mystr->pt5v0;
  v0Plus1Over2  = mystr->v0Plus1Over2;
  gamma_v0Over2 = mystr->gamma_pt5v0;
  akd           = mystr->tmpvec1_J;
  a1h_hatp0     = mystr->tmpvec2_J;
  tmp4          = mystr->tmpvec3_J;

  ca_all[0]     = c1;
  for (k=0; k<J; k++) {
    d1         = c1*a1d[k];
    a1h[k]     = d1;
    akh_all[k] = d1;
    akd[k]     = d1;
    d2         = mu0_lors[k];
    tmp4[k]    = d2*d2*tmp3;
  }

  ii = J;
  for (k=1; k<n; k++) {
    notMiss = notMissing[k];
    lrk     = lr[k];
    logor2k = logor2[k];
    sum_akd = 0.0;
    rowVecMat(akd, hatp0, J, J, a1h_hatp0);

    for (j=0; j<J; j++) {
      d1   = lrk-mu0_lrs[j];
      d1   = d1*d1;
      d2   = pow(pt5v0+(d1/tmp2), v0Plus1Over2);
      tmp5 = tmp1/(gamma_v0Over2*d2);
      
      if (!notMiss) {
        d1 = a1h_hatp0[j]*tmp5;
      } else { 
        d1 = a1h_hatp0[j]*(tmp5*(tmp3*dnchisq(logor2k*tmp3, 1, tmp4[j], 0)));
      }
      sum_akd += d1;
      akd[j]  = d1;
    }
    cka         = 1.0/sum_akd;
    ca_all[k]   = cka;
    for  (j=0; j<J; j++) {
      d1          = cka*akd[j];
      akd[j]      = d1;
      akh_all[ii] = d1;
      ii++;
    }
  }

  /* Compute sum(akh_all[n, ]) */
  pd = akh_all;
  mystr->sum_akh_all_n = sumvec(&pd[(n-1)*J], J);
  if (mystr->debug) Rprintf("End: get_akh_ca\n");

} 


static void get_bkh_cb(mystr)
MYSTR *mystr;
{
  if (mystr->debug) Rprintf("Begin: get_bkh_cb\n");

  int n, J, k, ii, kp1, j;
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, *tmp6, lrk1, logor2k1;
  double **hatp0, d1, d2, *bnh, *vec, *bkd, ckb;
  double *mu0_lors, *mu0_lrs, *cb_all, *bkh_all, *lr, *logor2;
  char *notMissing;
  
  n          = mystr->n;
  J          = mystr->J;
  tmp0       = mystr->pt5v0;
  tmp1       = mystr->TMP1;
  tmp2       = mystr->gamma_pt5v0;
  tmp3       = 2.0*mystr->sigma20_lr;
  tmp4       = mystr->v0Plus1Over2;
  tmp5       = mystr->sigma20_lor_inv;
  tmp6       = mystr->ncpVec;
  mu0_lors   = mystr->mu0_lors;
  mu0_lrs    = mystr->mu0_lrs;
  cb_all     = mystr->cb_all;
  bkh_all    = mystr->bkh_pks;
  bnh        = mystr->tmpvec1_J;
  bkd        = mystr->tmpvec2_J;
  vec        = mystr->tmpvec3_J;
  lr         = mystr->lr;
  logor2     = mystr->logor2;
  notMissing = mystr->notMissing;
  hatp0      = mystr->hatp0;
  
  for (k=0; k<J; k++) {
    d1      = mu0_lors[k];
    tmp6[k] = d1*d1*tmp5;
  }

  /* We are working backwards, so get results for final index, row */
  d1 = 1.0/((double) J);
  cb_all[n-1] = d1;
  ii = (n-1)*J;
  for (k=0; k<J; k++) {
    bkh_all[ii] = d1;
    ii++;
    bnh[k] = d1;
  }

  for (k=n-2; k>-1; k--) {
    kp1      = k + 1;
    logor2k1 = logor2[kp1];
    lrk1     = lr[kp1];
    if (!notMissing[kp1]) {
      for (j=0; j<J; j++) {
        d1     = lrk1 - mu0_lrs[j];
        d2     = pow(tmp0 + d1*d1/tmp3, tmp4);
        vec[j] = tmp1/(tmp2*d2)*bnh[j];
      }
    } else {
      for (j=0; j<J; j++) {
        d1     = lrk1 - mu0_lrs[j];
        d2     = pow(tmp0 + d1*d1/tmp3, tmp4);
        vec[j] = tmp1/(tmp2*d2)*(tmp5*dnchisq(logor2k1*tmp5, 1.0, tmp6[j], 0))*bnh[j];
      }
    }

    matColVec(vec, hatp0, J, J, bkd);
    d1  = 0.0;
    for (j=0; j<J; j++) d1 += bkd[j];
    ckb = 1.0/d1;

    cb_all[k] = ckb;
    ii        = k*J;   
    for (j=0; j<J; j++) {
      d1 = ckb*bkd[j];
      bkh_all[ii] = d1;
      ii++;
      bnh[j] = d1;
    }
  } 

  /* Save last element, vector will be log transformed */
  mystr->cb_all_n = cb_all[n-1];

  if (mystr->debug) Rprintf("End: get_bkh_cb\n");

} 

static void get_pzks1(mystr)
MYSTR *mystr;
{
  if (mystr->debug) Rprintf("Begin: get_pzks1\n");

   /*
    pzks             <- akh_all*bkh_all/sum(akh_all[n,])
    lca_all          <- log(ca_all)
    lcb_all          <- log(cb_all)
    pzks1            <- matrix(data=NA, nrow=n, ncol=J)
    pzks1[n ,]       <- pzks[n,]/cb_all[n]
    tmp1             <- reverseCumSum(lca_all)
    tmp2             <- reverseCumSum(lcb_all)
    pzks1[1:(n-1), ] <- exp(tmp1[2:n] - tmp2[1:(n-1)])*pzks[1:(n-1), ]
   */

  int i, j, n, J, nm1, ii;
  double *akh_all, *bkh_all, sum_akh_all_n, cb_all_n;
  double *vec1, *vec2, d1, d2, d3, *tmpvec1_n, *tmpvec2_n;

  /* cb_all has been log transformed */

  n             = mystr->n;
  nm1           = n - 1;
  J             = mystr->J;
  akh_all       = mystr->akh_all;
  bkh_all       = mystr->bkh_pks; /* Currently contains bkh_all */
  sum_akh_all_n = mystr->sum_akh_all_n;
  cb_all_n      = mystr->cb_all_n; 
  tmpvec1_n     = mystr->tmpvec1_n; /* reverseCumSum(log_ca_all, n, tmpvec1_n) */
  tmpvec2_n     = mystr->tmpvec2_n; /* reverseCumSum(cb_all, n, tmpvec2_n) */

  /* pzks1[n ,]       <- pzks[n,]/cb_all[n] 
     store in bkh_pks array */
  i    = nm1*J;
  vec1 = &akh_all[i];
  vec2 = &bkh_all[i];
  d1   = sum_akh_all_n*cb_all_n;

  for (j=0; j<J; j++) {
    vec2[j] = vec1[j]*vec2[j]/d1;
  }

  ii = 0;
  for (i=0; i<nm1; i++) {
    d1 = exp(tmpvec1_n[i+1] - tmpvec2_n[i]);
    d2 = d1/sum_akh_all_n;

    for (j=0; j<J; j++) {
      d3          = akh_all[ii]*bkh_all[ii];
      bkh_all[ii] = d2*d3;
      ii++;
    } 
  }

  if (mystr->debug) Rprintf("End: get_pzks1\n");

}

static void get_pzzks2(mystr)
MYSTR *mystr;
{
  if (mystr->debug) Rprintf("Begin: get_pzzks2\n");

  int k, i, j, n, J, nm1;
  double *vd1, xchisq, lrk, d1, d2, d3, d4, d5, d6, *bvec, *MAT, *avec, *pd, t2, tmp0, tmp1, tmp2, tmp3;
  double *logor2, *ncpVec, *mu0_lors, *lr, *bkh_all, *mu0_lrs, t3, *akh_all;
  double *tvec, *tveca, *tvecb, *ca_all, *vec, lrn, t1, **hatp0, **ret;
  char *notMissing, notMiss;

  vd1        = mystr->tmpvec1_J;
  MAT        = mystr->tmpvec2_J;
  vec        = mystr->tmpvec3_J;
  notMissing = mystr->notMissing;
  logor2     = mystr->logor2;
  lr         = mystr->lr;
  t1         = 2.0*mystr->sigma20_lr;
  t2         = mystr->sigma20_lor_inv;
  t3         = mystr->sum_akh_all_n;
  n          = mystr->n;
  J          = mystr->J;
  mu0_lors   = mystr->mu0_lors;
  mu0_lrs    = mystr->mu0_lrs;
  ncpVec     = mystr->ncpVec;
  bkh_all    = mystr->bkh_pks;
  tmp0       = mystr->pt5v0;
  tmp1       = mystr->TMP1;
  tmp2       = mystr->gamma_pt5v0;
  tmp3       = mystr->v0Plus1Over2;
  tveca      = mystr->tmpvec1_n;  /* reverseCumSum(lca_all) */
  tvecb      = mystr->tmpvec2_n;  /* reverseCumSum(lcb_all) */
  tvec       = mystr->tmpvec3_n;
  ca_all     = mystr->ca_all;
  nm1        = n - 1;
  hatp0      = mystr->hatp0;
  akh_all    = mystr->akh_all;
  ret        = mystr->pzzks;
 
  /* Initialize */
  for (i=0; i<J; i++) {
    for (j=0; j<J; j++) ret[i][j] = 0.0;
  }

  for (j=0; j<J; j++) {
    d1        = mu0_lors[j];
    ncpVec[j] = d1*d1*t2;
  }

  for (i=0; i<n; i++) {
    j       = 1 + i;
    d1      = exp(tveca[j+1] - tvecb[j]);
    d2      = ca_all[j]/t3;
    tvec[i] = d1*d2;
  }


  d1      = ca_all[nm1]*tmp1;
  lrn     = lr[nm1]; 
  d6      = logor2[nm1]*t2;
  notMiss = notMissing[nm1];
  for (j=0; j<J; j++) {
    d2 = lrn - mu0_lrs[j];
    d2 = d2*d2/t1;
    d3 = tmp0 + d2;
    d4 = tmp2*pow(d3, tmp3);
    d5 = d1/d4;
    if (notMiss) d5 = d5*t2*dnchisq(d6, 1, ncpVec[j], 0);
    vec[j] = d5/t3;
  }

  for (k=1; k<n-1; k++) {
    /* Get tmat[k, ] */
    if (notMissing[k]) {
      xchisq = t2*logor2[k];
      for (i=0; i<J; i++) {
        vd1[i] = t2*dnchisq(xchisq, 1, ncpVec[i], 0);
      }
    } else {
      for (i=0; i<J; i++) vd1[i] = 1.0;
    }

    lrk  = lr[k];
    bvec = &bkh_all[k*J];
    for (i=0; i<J; i++) {
      d1     = lrk-mu0_lrs[i];
      d2     = d1*d1;
      MAT[i] = bvec[i]*(tmp1/(tmp2*pow(tmp0+d2/t1, tmp3))*vd1[i]);
    }

    /* pzzks[,,k-1] <- tvec[k-1]*(hatp0*(akh_all[k-1,]%*%mat)) */

    d1  = tvec[k-1];
    pd = &akh_all[(k-1)*J]; 
    for (i=0; i<J; i++) {
      d2 = pd[i];
      for (j=0; j<J; j++) {
        ret[i][j] += d1*hatp0[i][j]*d2*MAT[j];
      }
    }
  }

  avec = &akh_all[(n-2)*J];
  for (i=0; i<J; i++) {
    d1 = avec[i];
    for (j=0; j<J; j++) {
      ret[i][j] += hatp0[i][j]*d1*vec[j];
    }
  }

  if (mystr->debug) Rprintf("End: get_pzzks2\n");


} /* END: get_pzzks2 */

static double sum_pzks1Vec(pzks1Vec, n11, n12, n21, n22, J)
double *pzks1Vec;
int n11, n12, n21, n22, J;
{
  int i, j;
  double sum=0.0, *vec;

  j = J - 1;
  if ((n21 > j) || (n22 > j)) error("INTERNAL CODING ERROR in sum_pzks1Vec");

  for (i=n11; i<=n12; i++) {
    vec = &pzks1Vec[i*J]; 
    for (j=n21; j<=n22; j++) sum += vec[j];
  }

  return(sum);

} /* sum_pzks1Vec */


static void get_r0s2(pzks1, n, J, mainJ, tauJ, hatp0_t0, r0, ret)
double *pzks1, *hatp0_t0, *r0, *ret;
int n, J, mainJ, tauJ;
{
  double d0, d1, d2, u01, u02;
  int i, j, i21, i22, ii, k;

  d0  = sum_pzks1Vec(pzks1, 0, 0, 0, J-1, J);

  /* Compute r0 
    r0 <- sapply(1:mainJ, function(j) sum(pzks1[1,c(j,(mainJ+1+tauJ*(j-1)):(mainJ+tauJ*j))])/sum(pzks1[1,]))
  */
  for (j=1; j<=mainJ; j++) {
    /*d1      = pzks1[0][j-1];*/
    d1      = pzks1[j-1]; 
    i21     = mainJ+1+tauJ*(j-1);
    i22     = mainJ+tauJ*j; 
    d2      = sum_pzks1Vec(pzks1, 0, 0, i21-1, i22-1, J);
    r0[j-1] = (d1 + d2)/d0; 
  }

  /* u0 <- c(sum(pzks1[1,1:mainJ]), sum(pzks1[1,(mainJ+1):J]))/sum(pzks1[1,]) */
  u01 = sum_pzks1Vec(pzks1, 0, 0, 0, mainJ-1, J)/d0;
  u02 = sum_pzks1Vec(pzks1, 0, 0, mainJ, J-1, J)/d0;

  /* 
   r0s <- c(r0*u0[1], c(sapply(1:mainJ, function(j) sapply(1:tauJ, function(t) r0[j]*u0[2]*hatp0_t0[(j-1)*tauJ+t]))))
  */


  for (j=0; j<mainJ; j++) ret[j] = r0[j]*u01;

  ii = mainJ;
  for (i=0; i<mainJ; i++) {
    d1 = r0[i]*u02;
    k  = i*tauJ;
    for (j=0; j<tauJ; j++) {
      ret[ii] = d1*hatp0_t0[k+j];
      ii++;
    }
  }


} /* END: get_r0s2 */

static void get_hatp0_t0_2(pzzks, n, J, mainJ, tauJ, pzks1Vec, ret)
double **pzzks, *pzks1Vec, *ret;
int n, J, mainJ, tauJ;
{
  int g, t, tmpi, i1, i2, i3, i4, mainJm1, i4m1, Jm1, i2m1, i3m1, i1ptm1;
  double *pret, dd1, dd2, dd3, denom, num, d2, d3;

  /* ret here is a vector of length mainJ*tauJ */

  pret    = ret;
  mainJm1 = mainJ - 1;
  Jm1     = J - 1;
  for (g=1; g<=mainJ; g++) {
    tmpi  = (g - 1)*tauJ;
    i1    = mainJ + tmpi;
    i2    = tmpi + 1;
    i3    = g*tauJ;
    i4    = mainJ + i3;
    i4m1  = i4 - 1; 
    i2m1  = i2 - 1;
    i3m1  = i3 - 1;
    dd1   = sum_pzks1Vec(pzks1Vec, 0, 0, i1, i4m1, J);
    dd2   = sum_pzzks(pzzks, 0,  mainJm1,   i1,   i4m1);
    dd3   = sum_pzzks_subset1(pzzks, mainJ, Jm1, i2m1, i3m1, mainJ+i2m1, i4m1);
    denom = dd1 + dd2 + dd3;

    for (t=1; t<=tauJ; t++) {
      i1ptm1  = i1 + t - 1;
      d2      = sum_pzzks(pzzks, 0, mainJm1, i1ptm1, i1ptm1);
      d3      = sum_pzzks_subset1(pzzks, mainJ, Jm1, i2m1, i3m1, i1ptm1, i1ptm1);
      num     = pzks1Vec[i1ptm1] + d2 + d3;
      *pret++ = num/denom;
    }
  }

  return;

} /* END: get_hatp0_t0_2 */


void get_hatp0(mystr)
MYSTR *mystr;
{
  if (mystr->debug) Rprintf("Begin: get_hatp0\n");

  int n, J, mainJ, tauJ;
  double *pzks1, **pzzks;

  n         = mystr->n;
  J         = mystr->J;
  mainJ     = mystr->mainJ;
  tauJ      = mystr->tauJ;
  pzzks     = mystr->pzzks;
  pzks1     = mystr->bkh_pks;

  /* pzks1 is now stored in bkh_pks ??? */

  get_hatp0_m(pzzks, n, J, mainJ, mystr->hatp0_m);
  get_hatp0_M(pzzks, n, J, mainJ, tauJ, mystr->hatp0_M);
  get_hatp0_t0_2(pzzks, n, J, mainJ, tauJ, pzks1, mystr->hatp0_t0);
  get_r0s2(pzks1, n, J, mainJ, tauJ, mystr->hatp0_t0, mystr->tmpvec1_J, mystr->r0s); 

  /* Must be called last */
  update_hatp0(mystr->hatp0, mainJ, tauJ, mystr->hatp0_M, mystr->hatp0_m, mystr->hatp0_t0, J);

  if (mystr->debug) Rprintf("End: get_hatp0\n");

} /* END: get_hatp0 */


static double get_loglike(log_ca_all, n, sum_akh_all_n) 
double *log_ca_all, sum_akh_all_n;
int n;
{
  double ret, sum;
  
  sum = sumvec(log_ca_all, n);
  ret = -sum + log(sum_akh_all_n);

  return(ret);

} 

static double negLoglike2(parm, nparm, p5L, p6L, p6U, n, J, mainJ, lr, logor2,
    notMissing, pzks1Vec, ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors)
int nparm, n, J, mainJ;
double *parm, p5L, p6L, p6U, *lr, *logor2, *pzks1Vec;
double *mu0_lrs, *mu0_lors, *ctzs, *ctms, *mz, *pz, *mz_sub, *pz_sub;
char *notMissing;
{
  int i, j;
  double p1, p2, p3, p4, p5, p6, p62, sigma20_lr, ep5, d1, d2, d3, mat, ret=0.0, p6p1, lri;
  double *vec, pzks1ij, logor2i, darg;
  char nomiss;

  p1         = parm[0];
  p2         = parm[1];
  p3         = parm[2];
  p4         = parm[3];
  p5         = parm[4];
  p6         = parm[5];
  sigma20_lr = exp(p4);
  
  /* Reparameterize */
  p62  = p6*p6;
  p6   = p6L + (p6U - p6L)*p62/(1.0 + p62);

  p5   = p5L + p5*p5;
  ep5  = exp(-p5);
  if (ep5 > MAX_EP5) error("ERROR: set bounds on parm[5]");
  p6p1 = p6 + 1.0;

  mu0_lr_func(J, mainJ, p1, p2, p3, ctzs, ctms, mu0_lrs);
  mu0_lor_func(J, mainJ, p2, p3, mz, pz, mz_sub, pz_sub, mu0_lors);


  /* Define the ncp and store in mu0_lors */
  for (i=0; i<J; i++) {
    d1          = mu0_lors[i];
    mu0_lors[i] = d1*d1*ep5;
  }

  d1 = -0.5*log(TWOPI*sigma20_lr) + log(gammafn(p6p1/2.0)) - log(gammafn(p6/2.0)) + 0.5*p6*log(p6/2.0);
  d2 = 0.5*p6p1;
  for (i=0; i<n; i++) {
    lri     = lr[i];
    logor2i = logor2[i];
    vec     = &pzks1Vec[i*J];
    nomiss  = notMissing[i];
    darg    = logor2i*ep5;
    for (j=0; j<J; j++) {
      d3      = lri - mu0_lrs[j];
      pzks1ij = vec[j];
      mat     = (d1 - d2*log(0.5*(p6 + d3*d3/sigma20_lr)))*pzks1ij;
      if (nomiss) mat = mat - pzks1ij*(p5 - dnchisq(darg, 1.0, mu0_lors[j], 1));
      ret += mat;
    }
  }
  ret = -ret;
/*Rprintf("%20.11f\n", ret);*/
  return(ret);

} /* END: negLoglike2 */

static void gradient2(eta, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1Vec, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors, ret)
double *eta, p5L, p6L, p6U, *lr, *logor2, *pzks1Vec, *mu0_lrs, *mu0_lors;
int nparm, nobs, J, mainJ;
double *ret, *ctzs, *ctms, *mz, *pz, *mz_sub, *pz_sub;
char *notMissing;
{
  int i;
  double fxplush, fxminush, h, save, h2, *ptr, *ptrret;

  /* default step size */
  h = 1e-3;
  h2 = 2.0*h;

  for (i=0, ptr=eta, ptrret=ret; i<nparm; i++, ptr++, ptrret++) {
    save = *ptr;  

    *ptr = save + h;  
    fxplush = negLoglike2(eta, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1Vec, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors);

    *ptr = save - h;  
    fxminush =  negLoglike2(eta, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1Vec, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors);
    *ptr = save;

    *ptrret = (fxplush - fxminush)/h2;
  }

} /* END: gradient2 */


static void myvmmin2(int n0, double *b, double *Fmin, 
      int maxit, int trace, double abstol, double reltol, int nREPORT, 
      int *fncount, int *grcount, int *fail,
      double p5L, double p6L, double p6U, int nobs, int J, int mainJ, double *lr, double *logor2, char *notMissing, 
      double *pzks1Vec, double *ctzs, double *ctms, double *mz, double *pz, 
      double *mz_sub, double *pz_sub, double *mu0_lrs, double *mu0_lors)
{
    int accpoint, enough;
    double *g, *t, *X, *c, **B;
    int   count, funcount, gradcount;
    double f, gradproj;
    int   i, j, ilast, iter = 0;
    double s, steplength;
    double D1, D2;
    int   n, *l, nparm=n0;

    if (maxit <= 0) {
	*fail = 0;
       *Fmin =  negLoglike2(b, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1Vec, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors);
	*fncount = *grcount = 0;
	return;
    }

    l = (int *) R_alloc(n0, sizeof(int));
    n = 0;
    /*for (i = 0; i < n0; i++) if (mask[i]) l[n++] = i;*/
    for (i = 0; i < n0; i++) l[n++] = i;
    g = (double *) R_alloc(n0, sizeof(double));
    t = (double *) R_alloc(n, sizeof(double));
    X = (double *) R_alloc(n, sizeof(double));
    c = (double *) R_alloc(n, sizeof(double));
    B = Lmatrix(n);
    f = negLoglike2(b, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1Vec, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors);
 
    if (!isfinite(f))
	Rprintf("initial value in 'vmmin' is not finite");
    if (trace) Rprintf("initial  value %f \n", f);
    *Fmin = f;
    funcount = gradcount = 1;
    gradient2(b, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1Vec, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors, g);
    iter++;
    ilast = gradcount;

    do {
	if (ilast == gradcount) {
	    for (i = 0; i < n; i++) {
		for (j = 0; j < i; j++) B[i][j] = 0.0;
		B[i][i] = 1.0;
	    }
	}
	for (i = 0; i < n; i++) {
	    X[i] = b[l[i]];
	    c[i] = g[l[i]];
	}
	gradproj = 0.0;
	for (i = 0; i < n; i++) {
	    s = 0.0;
	    for (j = 0; j <= i; j++) s -= B[i][j] * g[l[j]];
	    for (j = i + 1; j < n; j++) s -= B[j][i] * g[l[j]];
	    t[i] = s;
	    gradproj += s * g[l[i]];
	}

	if (gradproj < 0.0) {	/* search direction is downhill */
	    steplength = 1.0;
	    accpoint = 0;
	    do {

		count = 0;
		for (i = 0; i < n; i++) {
		    b[l[i]] = X[i] + steplength * t[i];
		    if (reltest + X[i] == reltest + b[l[i]]) /* no change */
			count++;
		}
		if (count < n) {
                  f = negLoglike2(b, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1Vec, 
                              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors);

		    funcount++;
		    accpoint = isfinite(f) &&
			(f <= *Fmin + gradproj * steplength * acctol);
		    if (!accpoint) {
			steplength *= stepredn;
		    }
		}

	    } while (!(count == n || accpoint));
	    enough = (f > abstol) && 
		fabs(f - *Fmin) > reltol * (fabs(*Fmin) + reltol);

	    /* stop if value if small or if relative change is low */
	    if (!enough) {
		count = n;
		*Fmin = f;
	    }

	    if (count < n) {/* making progress */
		*Fmin = f;
              gradient2(b, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1Vec, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors, g);

		gradcount++;
		iter++;
		D1 = 0.0;
		for (i = 0; i < n; i++) {
		    t[i] = steplength * t[i];
		    c[i] = g[l[i]] - c[i];
		    D1 += t[i] * c[i];
		}
		if (D1 > 0) {
		    D2 = 0.0;
		    for (i = 0; i < n; i++) {
			s = 0.0;
			for (j = 0; j <= i; j++)
			    s += B[i][j] * c[j];
			for (j = i + 1; j < n; j++)
			    s += B[j][i] * c[j];
			X[i] = s;
			D2 += s * c[i];
		    }
		    D2 = 1.0 + D2 / D1;
		    for (i = 0; i < n; i++) {
			for (j = 0; j <= i; j++)
			    B[i][j] += (D2 * t[i] * t[j]
					- X[i] * t[j] - t[i] * X[j]) / D1;
		    }
		} else {	/* D1 < 0 */
		    ilast = gradcount;
		}
	    } else {	/* no progress */
		if (ilast < gradcount) {
		    count = 0;
		    ilast = gradcount;
		}
	    }
	} else {		/* uphill search */
	    count = 0;
	    if (ilast == gradcount) count = n;
	    else ilast = gradcount;
	    /* Resets unless has just been reset */
	}
	if (trace && (iter % nREPORT == 0))
	    Rprintf("iter%4d value %f\n", iter, f);
	if (iter >= maxit) break;
	if (gradcount - ilast > 2 * n)
	    ilast = gradcount;	/* periodic restart */

    } while (count != n || ilast != gradcount);
    if (trace) {
	Rprintf("final  value %f \n", *Fmin);
	if (iter < maxit) Rprintf("converged\n");
	else Rprintf("stopped after %i iterations\n", iter);
    }
    *fail = (iter < maxit) ? 0 : 1;
    *fncount = funcount;
    *grcount = gradcount;
    /*
    if (op->debug) {
      Rprintf("fncount = %d   grcount = %d \n", funcount, gradcount);
    }
    */

} /* END: myvmmin2 */

static int call_optim(mystr)
MYSTR *mystr;
{
  if (mystr->debug) Rprintf("Begin: call_optim\n");

  int nparm, nobs, J, mainJ;
  int nREPORT, fail, maxit, trace, fncount, grcount;
  double p5L, p6L, p6U, abstol, reltol, Fmin, *pzks1, *mu0_lrs, *mu0_lors, *parm, p5, p6;

  fail      = -1;
  nparm     = NPARM;
  nobs      = mystr->n;
  p5L       = mystr->P5_LB;
  p6L       = mystr->V0_LB;
  p6U       = mystr->V0_UB;
  J         = mystr->J;
  mainJ     = mystr->mainJ;
  reltol    = mystr->reltol;
  parm      = mystr->parms;

  nREPORT   = 1;
  abstol    = -1.0e100;
  fail      = 1;
  fncount   = 0;
  grcount   = 0;
  maxit     = 100;
  trace     = 0;

  mu0_lrs  = mystr->tmpvec1_J;
  mu0_lors = mystr->tmpvec2_J;
  pzks1    = mystr->bkh_pks;

  /* Reparameterize */
  p6      = parm[5];
  parm[5] = sqrt((p6 - p6L)/(p6U - p6));
  parm[4] = sqrt(parm[4] - p5L);


  myvmmin2(nparm, parm, &Fmin, maxit, trace, abstol, reltol, nREPORT, 
      &fncount, &grcount, &fail,
      p5L, p6L, p6U, nobs, J, mainJ, mystr->lr, mystr->logor2, mystr->notMissing, 
      pzks1, mystr->ctzs, mystr->ctms, mystr->mz, mystr->pz, mystr->mz_sub, 
      mystr->pz_sub, mu0_lrs, mu0_lors);

  p6      = parm[5];
  p6      = p6*p6;
  parm[5] = p6L + (p6U - p6L)*(p6/(1 + p6));
  p5      = parm[4];
  parm[4] = p5L + p5*p5; 

  if (mystr->debug) Rprintf("End: call_optim\n");

  return(fail);

} /* END: call_optim */

static void print_parms(mystr)
MYSTR *mystr;
{
  if (mystr->debug) Rprintf("Begin: print_parms\n");

  double d1, d2, d3, d4, d5, d6;
  Rprintf("ploidy  purity  clonal.prop  logR.var  logOR.var  df\n");

  d1 = pow(2.0, mystr->lpsi0);
  d2 = mystr->pr0;
  d3 = 1.0/(1.0 + exp(mystr->lsz0));
  d4 = mystr->sigma20_lr;
  d5 = mystr->sigma20_lor;
  d6 = mystr->v0;
  Rprintf("%g %g %g %g %g %g\n", d1, d2, d3, d4, d5, d6);

  if (mystr->debug) Rprintf("End: print_parms\n");

}

static void print_init_loglike(mystr)
MYSTR *mystr;
{
  if (mystr->debug) Rprintf("Begin: print_init_loglike\n");

  Rprintf("Initial loglikelihood = %g\n", mystr->ll0);
  if (mystr->print > 1) print_parms(mystr);

  if (mystr->debug) Rprintf("End: print_init_loglike\n");

}


static int checkStop(ll0, ll1, iter, mystr)
double ll0, ll1;
int iter;
MYSTR *mystr;
{
  if (mystr->debug) Rprintf("Begin: checkStop\n");

  double v1, v2;
  int ret=0;

  v1  = fabs(ll0 - ll1);
  v2  = maxAbsDiff(mystr->parms0, mystr->parms, NPARM);
  if ((v1 < mystr->loglike_eps) || (v2 < mystr->parm_eps)) ret=1;
  
  if (mystr->print) {
    Rprintf("Iter %d loglike = %g loglike diff = %g max parm diff = %g\n", iter, ll1, v1, v2);
    if (mystr->print > 1) print_parms(mystr);
  }

  if (mystr->debug) Rprintf("End: checkStop\n");

  return(ret);

} /* END: checkStop */

static int EM_alg_iter(mystr)
MYSTR *mystr;
{
  if (mystr->debug) Rprintf("Begin: EM_alg_iter\n");

  int iter, rc=0, ret=-1, n;
  double ll0, ll1, *parms0, *parms1;

  ll0    = mystr->ll0;
  parms0 = mystr->parms0;
  parms1 = mystr->parms;
  n      = mystr->n;

  for (iter=1; iter<=mystr->maxiter; iter++) {

    get_bkh_cb(mystr);
    reverseCumSum(mystr->log_ca_all, n, mystr->tmpvec1_n);
    logvec(mystr->cb_all, n, mystr->cb_all);
    reverseCumSum(mystr->cb_all, n, mystr->tmpvec2_n);
    get_pzzks2(mystr); 
    get_pzks1(mystr);  
    get_hatp0(mystr);

    rc = call_optim(mystr); 
    if (rc) error("optimizer did not converge");
    updateParms(mystr->parms, mystr);

    mystr->sum_a1d = get_a1d(mystr, mystr->a1d);
    get_akh_ca(mystr);
    logvec(mystr->ca_all, n, mystr->log_ca_all);
    ll1 = get_loglike(mystr->log_ca_all, n, mystr->sum_akh_all_n);
    if (!R_FINITE(ll1)) error("Non-finite log-likelihood");
    if (checkStop(ll0, ll1, iter, mystr)) {
      mystr->converged = 1;
      mystr->niter     = iter;
      ret              = 0;
      break;
    }
    copy_dVec(parms1, parms0, NPARM);
    ll0 = ll1;
  }
  if (mystr->debug) Rprintf("End: EM_alg_iter\n");

  return(ret);

}

static void EM_alg(mystr)
MYSTR *mystr;
{
  if (mystr->debug) Rprintf("Begin: EM_alg\n");
  int n=mystr->n;

  mystr->sum_a1d = get_a1d(mystr, mystr->a1d); 
  get_akh_ca(mystr);  /* Also computes sum(akh_all[n]) */
  logvec(mystr->ca_all, mystr->n, mystr->log_ca_all);
  mystr->ll0 = get_loglike(mystr->log_ca_all, n, mystr->sum_akh_all_n);
  if (!R_FINITE(mystr->ll0)) error("Non-finite initial log-likelihood");
  if (mystr->print) print_init_loglike(mystr);
  EM_alg_iter(mystr);
  if (mystr->debug) Rprintf("End: EM_alg\n");


} /* END: EM_alg */

/*
static void get_maxcol(mystr)
MYSTR *mystr;
{
  int i, k, mainJ, J, n, ii, *maxcol1, *maxcol2, *maxcol3, m2, m3;
  double *pzks1, val, maxval2, maxval3;

  pzks1   = mystr->bkh_pks;
  maxcol1 = mystr->maxcol1;
  maxcol2 = mystr->maxcol2;
  maxcol3 = mystr->maxcol3;
  mainJ   = mystr->mainJ;
  J       = mystr->J;
  n       = mystr->n;
  ii      = 0;
  
  for (i=0; i<n; i++) {
    maxval2 = LARGE_NEG;
    maxval3 = LARGE_NEG;
    m2      = 0;
    m3      = mainJ;
    for (k=0; k<mainJ; k++) {
      val = pzks1[ii];
      ii++;
      if (val > maxval2) {
        maxval2 = val;
        m2      = k;
      }
    }
    for (k=mainJ; k<J; k++) {
      val = pzks1[ii];
      ii++;
      if (val > maxval3) {
        maxval3 = val;
        m3      = k;
      }
    }
    if (maxval2 > maxval3) {
      maxcol1[i] = m2;
    } else {
      maxcol1[i] = m3;
    }
    maxcol2[i] = m2;
    maxcol3[i] = m3;
  }
}  
*/

/* hatp0_M, hatp0_m, hatp0_t0, r0s are input and output */
void C_EM_alg(iargs, dargs, lr, logor2, mz, mz_sub, pz, pz_sub,
     ctzs, ctms, r0, u0, hatp0_M, hatp0_m, hatp0_t0, r0s, 
     ret_parms, ret_conv, ret_niter, ret_mu0_lrs, ret_mu0_lors, ret_loglike,
     ret_pzks1, ret_hatp0)
int *iargs, *ret_conv, *ret_niter;
double *mz, *mz_sub, *pz, *pz_sub, *ctzs, *ctms;
double *dargs, *lr, *logor2, *hatp0_M, *hatp0_m, *hatp0_t0, *r0s, *ret_hatp0;
double *ret_parms, *ret_mu0_lrs, *ret_mu0_lors, *r0, *u0, *ret_loglike, *ret_pzks1;
{
  MYSTR mystr;

  mystr_init(&mystr, iargs, dargs, lr, logor2, mz, mz_sub, pz, pz_sub,
             ctzs, ctms, hatp0_M, hatp0_m, hatp0_t0, r0s, r0, u0, 
             ret_parms, ret_mu0_lrs, ret_mu0_lors, ret_pzks1, ret_hatp0);

  EM_alg(&mystr);

  /* Return: parms, conv, ll1, hatp0_M, hatp0_m, hatp0_t0, mu0_lrs, mu0_lors */
  *ret_conv    = mystr.converged;
  *ret_niter   = mystr.niter;
  *ret_loglike = mystr.ll1;
  copyMatIntoVecByCol(mystr.hatp0_M, mystr.mainJ, mystr.mainJ, hatp0_M);  

  mystr_free(&mystr);

  return;

}  /* END: C_EM_alg */



void R_init_subHMM(DllInfo *dll)
{
    R_registerRoutines(dll, callMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
}

