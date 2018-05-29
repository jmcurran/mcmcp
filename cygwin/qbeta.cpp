#include <math.h>
#include <float.h>
#include <errno.h>
#include <iostream>
#include <fstream>
using namespace std;

double lgamma(double arg);

double lbeta(double a, double b)
{
	return lgamma(a)+lgamma(b)-lgamma(a+b);
}

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998--2001  The R Development Core Team
 *  based on code (C) 1979 and later Royal Statistical Society
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
 */

double fmax2(double x, double y)
{
	return (x < y) ? y : x;
}

double fmin2(double x, double y)
{
	return (x < y) ? x : y;
}

double log1p(double x){
	if (x == 0.) 
		return 0.;/* speed */
	return log(1+x);
}




/*  SYNOPSIS
 *
 *	#include <Rmath.h>
 *	double expm1(double x);
 *
 *  DESCRIPTION
 *
 *	Compute the Exponential minus 1
 *
 *			exp(x) - 1
 *
 *      accurately also when x is close to zero, i.e. |x| << 1
 *
 *  NOTES
 *
 *	As log1p(), this is a standard function in some C libraries,
 *	particularly GNU and BSD (but is neither ISO/ANSI C nor POSIX).
 *
 *  We supply a substitute for the case when there is no system one.
 */

double expm1(double x)
{
    double y, a = fabs(x);

    if (a < DBL_EPSILON) return x;
    if (a > 0.697) return exp(x) - 1;  /* negligible cancellation */

    if (a > 1e-8)
	y = exp(x) - 1;
    else /* Taylor expansion, more accurate in this range */
	y = (x / 2 + 1) * x;

    /* Newton step for solving   log(1 + y) = x   for y : */
    /* WARNING: does not work for y ~ -1: bug in 1.5.0 */
    y -= (1 + y) * (log1p (y) - x);
    return y;
}

 /* double pbeta_raw(double x, double pin, double qin, int lower_tail)
 *
 *  DESCRIPTION
 *
 *	Returns distribution function of the beta distribution.
 *	( = The incomplete beta ratio I_x(p,q) ).
 *
 *  NOTES
 *
 *	This routine is a translation into C of a Fortran subroutine
 *	by W. Fullerton of Los Alamos Scientific Laboratory.
 *
 *  REFERENCE
 *
 *	Bosten and Battiste (1974).
 *	Remark on Algorithm 179, CACM 17, p153, (1974).
 */

/* This is called from	qbeta(.) in a root-finding loop --- be FAST! */

double pbeta_raw(double x, double pin, double qin, int lower_tail)
{
    double ans, c, finsum, p, ps, p1, q, term, xb, xi, y;
    int n, i, ib, swap_tail;

    const double eps = .5*DBL_EPSILON;
    const double sml = DBL_MIN;
    const double lneps = log(eps);
    const double lnsml = log(sml);

    /* swap tails if x is greater than the mean */
    if (pin / (pin + qin) < x) {
	swap_tail = 1;
	y = 1 - x;
	p = qin;
	q = pin;
    }
    else {
	swap_tail = 0;
	y = x;
	p = pin;
	q = qin;
    }

    if ((p + q) * y / (p + 1) < eps) {

	/* tail approximation */

	xb = p * log(fmax2(y, sml)) - log(p) - lbeta(p, q);
	if (xb > lnsml && y != 0) {
	    ans = (swap_tail == lower_tail) ? -expm1(xb) : exp(xb);
	} else {
	    ans = (swap_tail == lower_tail) ? 1. : 0;
	}
    }
    else {
	/* MM: __ FIXME __ : This takes forever (or ends wrongly)
	  when (one or) both p & q  are huge

	  ./pf.c  now has a cheap fix -- can use that here, but better
	  "get it right"  (PD to R-core on 20 Feb 2000)a
	*/

	/* evaluate the infinite sum first.  term will equal */
	/* y^p / beta(ps, p) * (1 - ps)-sub-i * y^i / fac(i) */

	/* Ly := log(y) */
	double Ly = swap_tail ? log1p(-x) : log(y);

	ps = q - floor(q);
	xb = p * Ly;
	if (ps == 0)
	    ps = 1; /*==> lbeta(ps,p)= log Beta(1,p) = log(1/p) = -log(p) */
	else
	    xb -= (lbeta(ps, p) + log(p));
	ans = 0;
	if (xb >= lnsml) {
	    ans = exp(xb);
	    term = ans * p;
	    if (ps != 1) {
		n = fmax2(lneps/Ly, 4.0);
		for(i=1 ; i <= n ; i++) {
		    xi = i;
		    term *= (xi - ps) * y / xi;
		    ans += term / (p + xi);
		}
	    }
	}

	/* now evaluate the finite sum, maybe. */

	if (q > 1) {

	    double liy;/* == log(1-y) */
	    if(swap_tail) {
		c = 1./x;/* == 1/(1 - y) */
		liy = log(x);
	    }
	    else {
		c = 1./(1. - y);
		liy = log1p(-y);
	    }
	    xb = p * Ly + q * liy - lbeta(p, q) - log(q);
	    ib = fmax2(xb / lnsml, 0.0);
	    term = exp(xb - ib * lnsml);
	    p1 = q * c / (p + q - 1);

	    finsum = 0;
	    n = q;
	    if (q == n)
		n--;
	    for(i= 1; i <= n; i++) {
		if (p1 <= 1 && term / eps <= finsum)
		    break;
		xi = i;
		term = (q - xi + 1) * c * term / (p + q - xi);
		if (term > 1) {
		    ib--;
		    term *= sml;
		}
		if (ib == 0)
		    finsum += term;
	    }
	    ans += finsum;
	}
	if (swap_tail == lower_tail)
	    ans = 1 - ans;
	ans = fmax2(fmin2(ans, 1.), 0.);
    }
    return ans;
}

/* Reference:
 * Cran, G. W., K. J. Martin and G. E. Thomas (1977).
 *	Remark AS R19 and Algorithm AS 109,
 *	Applied Statistics, 26(1), 111-114.
 * Remark AS R83 (v.39, 309-310) and the correction (v.40(1) p.236)
 *	have been incorporated in this version.
 */



/* set the exponent of accu to -2r-2 for r digits of accuracy */
/*---- NEW ---- -- still fails for p = 1e11, q=.5*/

#define fpu 3e-308
/* acu_min:  Minimal value for accuracy 'acu' which will depend on (a,p);
	     acu_min >= fpu ! */
#define acu_min 1e-300
#define lower fpu
#define upper 1-2.22e-16

#define const1 2.30753
#define const2 0.27061
#define const3 0.99229
#define const4 0.04481

static volatile double xtrunc;/* not a real global .. delicate though! */

double qbeta(double alpha, double p, double q, int lower_tail, int log_p)
{
    int swap_tail, i_pb, i_inn;
    double a, adj, logbeta, g, h, pp, p_, prev, qq, r, s, t, tx, w, y, yprev;
    double acu;
    volatile double xinbta;

    /* define accuracy and initialize */

    xinbta = alpha;

    /* test for admissibility of parameters */

	if((log_p && alpha>0) || (!log_p && (alpha<0 || alpha>1)))
		return -DBL_MAX*(1-1e-15);
    
    if(p < 0. || q < 0.)
		return -DBL_MAX*(1-1e-15);

	p_ = log_p ? (lower_tail ? exp(alpha) : -expm1(alpha)) : (lower_tail ? (alpha) : (1 - (alpha))); /* lower_tail prob (in any case) */

    if (p_ == 0. || p_ == 1.)
	return p_;

    logbeta = lbeta(p, q);

    /* change tail if necessary;  afterwards   0 < a <= 1/2	 */
    if (p_ <= 0.5) {
	a = p_;	pp = p; qq = q; swap_tail = 0;
    } else { /* change tail, swap  p <-> q :*/
	a = (!lower_tail && !log_p)? alpha : 1 - p_;
	pp = q; qq = p; swap_tail = 1;
    }

    /* calculate the initial approximation */

    r = sqrt(-2 * log(a));
    y = r - (const1 + const2 * r) / (1. + (const3 + const4 * r) * r);
    if (pp > 1 && qq > 1) {
	r = (y * y - 3.) / 6.;
	s = 1. / (pp + pp - 1.);
	t = 1. / (qq + qq - 1.);
	h = 2. / (s + t);
	w = y * sqrt(h + r) / h - (t - s) * (r + 5. / 6. - 2. / (3. * h));
	xinbta = pp / (pp + qq * exp(w + w));
    } else {
	r = qq + qq;
	t = 1. / (9. * qq);
	t = r * pow(1. - t + y * sqrt(t), 3.0);
	if (t <= 0.)
	    xinbta = 1. - exp((log1p(-a)+ log(qq) + logbeta) / qq);
	else {
	    t = (4. * pp + r - 2.) / t;
	    if (t <= 1.)
		xinbta = exp((log(a * pp) + logbeta) / pp);
	    else
		xinbta = 1. - 2. / (t + 1.);
	}
    }

    /* solve for x by a modified newton-raphson method, */
    /* using the function pbeta_raw */

    r = 1 - pp;
    t = 1 - qq;
    yprev = 0.;
    adj = 1;
    /* Sometimes the approximation is negative! */
    if (xinbta < lower)
	xinbta = 0.5;
    else if (xinbta > upper)
	xinbta = 0.5;

    /* Desired accuracy should depend on  (a,p)
     * This is from Remark .. on AS 109, adapted.
     * However, it's not clear if this is "optimal" for IEEE double prec.

     * acu = fmax2(acu_min, pow(10., -25. - 5./(pp * pp) - 1./(a * a)));

     * NEW: 'acu' accuracy NOT for squared adjustment, but simple;
     * ---- i.e.,  "new acu" = sqrt(old acu)

     */
    acu = fmax2(acu_min, pow(10., -13 - 2.5/(pp * pp) - 0.5/(a * a)));
    tx = prev = 0.;	/* keep -Wall happy */

    for (i_pb=0; i_pb < 1000; i_pb++) {
	y = pbeta_raw(xinbta, pp, qq, /*lower_tail = */ true);
	/* y = pbeta_raw2(xinbta, pp, qq, logbeta) -- to SAVE CPU; */
#ifdef IEEE_754
	if(!R_FINITE(y))
#else
	    if (errno)
#endif
		return -DBL_MAX*(1-1e-15);

	y = (y - a) *
	    exp(logbeta + r * log(xinbta) + t * log1p(-xinbta));
	if (y * yprev <= 0.)
	    prev = fmax2(fabs(adj),fpu);
	g = 1;
	for (i_inn=0; i_inn < 1000;i_inn++) {
	    adj = g * y;
	    if (fabs(adj) < prev) {
		tx = xinbta - adj; /* trial new x */
		if (tx >= 0. && tx <= 1) {
		    if (prev <= acu)	goto L_converged;
		    if (fabs(y) <= acu) goto L_converged;
		    if (tx != 0. && tx != 1)
			break;
		}
	    }
	    g /= 3;
	}
	xtrunc = tx;	/* this prevents trouble with excess FPU */
				/* precision on some machines. */
	if (xtrunc == xinbta)
	    goto L_converged;
	xinbta = tx;
	yprev = y;
    }
    /*-- NOT converged: Iteration count --*/
 //   ML_ERROR(ME_PRECISION);

 L_converged:
    if (swap_tail)
	return 1 - xinbta;
    return xinbta;
}

#ifndef _goldenratio
#define _goldenratio 0.618033988749895
#define _oneminusgoldenratio 0.3819660112501052
#endif

void fitbeta(double m, double dUpperQuantile, double *a, double *b, double dP, double dLowerLim, double dUpperLim)
{
	// assume coefficients are in [1,100] - should be safe enough
	// using Golden section search

	double dA = dLowerLim;
	double dB = dUpperLim;

	double alpha_lower = dA + (dB-dA)*_oneminusgoldenratio;
	double alpha_upper = dA + (dB-dA)*_goldenratio;

	double beta_lower=(alpha_lower*(1-m)+2*m-1)/m;
	double beta_upper=(alpha_upper*(1-m)+2*m-1)/m;

	double q_lower=fabs(qbeta(dP,alpha_lower,beta_lower,true,false)-dUpperQuantile);
	double q_upper=fabs(qbeta(dP,alpha_upper,beta_upper,true,false)-dUpperQuantile);

	while(q_lower>1e-3&&q_upper>1e-3){
		if(q_lower>q_upper){
			dA = alpha_lower;
			alpha_lower = alpha_upper;
			q_lower = q_upper;

			alpha_upper = dA + (dB-dA)*_goldenratio;
			beta_upper=(alpha_upper*(1-m)+2*m-1)/m;
			q_upper=fabs(qbeta(dP,alpha_upper,beta_upper,true,false)-dUpperQuantile);
		}else{
			dB = alpha_upper;
			alpha_upper = alpha_lower;
			q_upper = q_lower;

			alpha_lower = dA + (dB-dA)*_oneminusgoldenratio;
			beta_lower=(alpha_lower*(1-m)+2*m-1)/m;
			q_lower=fabs(qbeta(dP,alpha_lower,beta_lower,true,false)-dUpperQuantile);
		}
	}

	if(q_lower<1e-3)
	{
		*a = alpha_lower;
		*b = beta_lower;
	}
	/* else */
	*a = alpha_upper;
	*b = beta_upper;

}
