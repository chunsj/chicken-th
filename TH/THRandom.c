#include "THGeneral.h"
#include "THRandom.h"

#ifndef _WIN32
#include <fcntl.h>
#include <unistd.h>
#endif

/* Code for the Mersenne Twister random generator.... */
#define n _MERSENNE_STATE_N
#define m _MERSENNE_STATE_M

/* Creates (unseeded) new generator*/
static THGenerator* THGenerator_newUnseeded()
{
  THGenerator *self = THAlloc(sizeof(THGenerator));
  memset(self, 0, sizeof(THGenerator));
  self->left = 1;
  self->seeded = 0;
  self->normal_is_valid = 0;
  return self;
}

/* Creates new generator and makes sure it is seeded*/
THGenerator* THGenerator_new()
{
  THGenerator *self = THGenerator_newUnseeded();
  THRandom_seed(self);
  return self;
}

THGenerator* THGenerator_copy(THGenerator *self, THGenerator *from)
{
    memcpy(self, from, sizeof(THGenerator));
    return self;
}

void THGenerator_free(THGenerator *self)
{
  THFree(self);
}

int THGenerator_isValid(THGenerator *_generator)
{
  if ((_generator->seeded == 1) &&
    (_generator->left > 0 && _generator->left <= n) && (_generator->next <= n))
    return 1;

  return 0;
}

#ifndef _WIN32
static unsigned long readURandomLong()
{
  int randDev = open("/dev/urandom", O_RDONLY);
  unsigned long randValue;
  if (randDev < 0) {
    THError("Unable to open /dev/urandom");
  }
  ssize_t readBytes = read(randDev, &randValue, sizeof(randValue));
  if (readBytes < sizeof(randValue)) {
    THError("Unable to read from /dev/urandom");
  }
  close(randDev);
  return randValue;
}
#endif // _WIN32

unsigned long THRandom_seed(THGenerator *_generator)
{
#ifdef _WIN32
  unsigned long s = (unsigned long)time(0);
#else
  unsigned long s = readURandomLong();
#endif
  THRandom_manualSeed(_generator, s);
  return s;
}

/* The next 4 methods are taken from http:www.math.keio.ac.jpmatumotoemt.html
   Here is the copyright:
   Some minor modifications have been made to adapt to "my" C... */

/*
   A C-program for MT19937, with initialization improved 2002/2/10.
   Coded by Takuji Nishimura and Makoto Matsumoto.
   This is a faster version by taking Shawn Cokus's optimization,
   Matthe Bellew's simplification, Isaku Wada's double version.

   Before using, initialize the state by using init_genrand(seed)
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.keio.ac.jp/matumoto/emt.html
   email: matumoto@math.keio.ac.jp
*/

/* Macros for the Mersenne Twister random generator... */
/* Period parameters */
/* #define n 624 */
/* #define m 397 */
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UMASK 0x80000000UL /* most significant w-r bits */
#define LMASK 0x7fffffffUL /* least significant r bits */
#define MIXBITS(u,v) ( ((u) & UMASK) | ((v) & LMASK) )
#define TWIST(u,v) ((MIXBITS(u,v) >> 1) ^ ((v)&1UL ? MATRIX_A : 0UL))
/*********************************************************** That's it. */

void THRandom_manualSeed(THGenerator *_generator, unsigned long the_seed_)
{
  int j;

  /* This ensures reseeding resets all of the state (i.e. state for Gaussian numbers) */
  THGenerator *blank = THGenerator_newUnseeded();
  THGenerator_copy(_generator, blank);
  THGenerator_free(blank);

  _generator->the_initial_seed = the_seed_;
  _generator->state[0] = _generator->the_initial_seed & 0xffffffffUL;
  for(j = 1; j < n; j++)
  {
    _generator->state[j] = (1812433253UL * (_generator->state[j-1] ^ (_generator->state[j-1] >> 30)) + j);
    /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
    /* In the previous versions, mSBs of the seed affect   */
    /* only mSBs of the array state[].                        */
    /* 2002/01/09 modified by makoto matsumoto             */
    _generator->state[j] &= 0xffffffffUL;  /* for >32 bit machines */
  }
  _generator->left = 1;
  _generator->seeded = 1;
}

unsigned long THRandom_initialSeed(THGenerator *_generator)
{
  return _generator->the_initial_seed;
}

void THRandom_nextState(THGenerator *_generator)
{
  unsigned long *p = _generator->state;
  int j;

  _generator->left = n;
  _generator->next = 0;

  for(j = n-m+1; --j; p++)
    *p = p[m] ^ TWIST(p[0], p[1]);

  for(j = m; --j; p++)
    *p = p[m-n] ^ TWIST(p[0], p[1]);

  *p = p[m-n] ^ TWIST(p[0], _generator->state[0]);
}

unsigned long THRandom_random(THGenerator *_generator)
{
  unsigned long y;

  if (--(_generator->left) == 0)
    THRandom_nextState(_generator);
  y = *(_generator->state + (_generator->next)++);

  /* Tempering */
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);

  return y;
}

/* generates a random number on [0,1)-double-interval */
static double __uniform__(THGenerator *_generator)
{
  /* divided by 2^32 */
  return (double)THRandom_random(_generator) * (1.0/4294967296.0);
}

/*********************************************************

 Thanks *a lot* Takuji Nishimura and Makoto Matsumoto!

 Now my own code...

*********************************************************/

double THRandom_uniform(THGenerator *_generator, double a, double b)
{
  return(__uniform__(_generator) * (b - a) + a);
}

double THRandom_normal(THGenerator *_generator, double mean, double stdv)
{
  THArgCheck(stdv > 0, 2, "standard deviation must be strictly positive");

  /* This is known as the Box-Muller method */
  if(!_generator->normal_is_valid)
  {
    _generator->normal_x = __uniform__(_generator);
    _generator->normal_y = __uniform__(_generator);
    _generator->normal_rho = sqrt(-2. * log(1.0-_generator->normal_y));
    _generator->normal_is_valid = 1;
  }
  else
    _generator->normal_is_valid = 0;

  if(_generator->normal_is_valid)
    return _generator->normal_rho*cos(2.*M_PI*_generator->normal_x)*stdv+mean;
  else
    return _generator->normal_rho*sin(2.*M_PI*_generator->normal_x)*stdv+mean;
}

double THRandom_exponential(THGenerator *_generator, double lambda)
{
  return(-1. / lambda * log(1-__uniform__(_generator)));
}

double THRandom_cauchy(THGenerator *_generator, double median, double sigma)
{
  return(median + sigma * tan(M_PI*(__uniform__(_generator)-0.5)));
}

/* Faut etre malade pour utiliser ca.
   M'enfin. */
double THRandom_logNormal(THGenerator *_generator, double mean, double stdv)
{
  THArgCheck(stdv > 0, 2, "standard deviation must be strictly positive");
  return(exp(THRandom_normal(_generator, mean, stdv)));
}

double THRandom_gamma(THGenerator *state, double shape)
{
  if (shape < 1.) {
    double u, v, x, y;
    while (1) {
      u = __uniform__(state);
      v = THRandom_exponential(state, 1);
      if (u <= 1.0 - shape) {
        x = pow(u, 1./shape);
        if (x <= v) {
          return x;
        }
      }
      else {
        y = -log((1 - u)/shape);
        x = pow(1.0 - shape + shape*y, 1./shape);
        if (x <= (v + y)) {
          return x;
        }
      }
    }
  }
  else if (shape > 1.) {
    double d = shape - (1./3.);
    double c = 1./sqrt(9. * d);
    double u, v, x = 0;
    do {
      x = THRandom_normal(state, 0, 1);
      v = (1 + c * x) * (1 + c * x) * (1 + c * x);
      u = __uniform__(state);
    } while (v <= 0. ||
             (((log(u) >= 0.5 * x * x + d * (1 - v + log(v)))) &&
              (u < 1.0 - 0.0331*(x*x)*(x*x))));
    return d * v;
  }
  else {
    return THRandom_exponential(state, 1);
  }
}

double THRandom_gamma2 (THGenerator *state, double a, double scale)
{
  const static double sqrt32 = 5.656854;
  const static double exp_m1 = 0.36787944117144232159;/* exp(-1) = 1/e */
  /* Coefficients q[k] - for q0 = sum(q[k]*a^(-k))
   * Coefficients a[k] - for q = q0+(t*t/2)*sum(a[k]*v^k)
   * Coefficients e[k] - for exp(q)-1 = sum(e[k]*q^k)
   */
  const static double q1 = 0.04166669;
  const static double q2 = 0.02083148;
  const static double q3 = 0.00801191;
  const static double q4 = 0.00144121;
  const static double q5 = -7.388e-5;
  const static double q6 = 2.4511e-4;
  const static double q7 = 2.424e-4;

  const static double a1 = 0.3333333;
  const static double a2 = -0.250003;
  const static double a3 = 0.2000062;
  const static double a4 = -0.1662921;
  const static double a5 = 0.1423657;
  const static double a6 = -0.1367177;
  const static double a7 = 0.1233795;

  /* State variables [FIXME for threading!] :*/
  static double aa = 0.;
  static double aaa = 0.;
  static double s, s2, d;    /* no. 1 (step 1) */
  static double q0, b, si, c;/* no. 2 (step 4) */

  double e, p, q, r, t, u, v, w, x, ret_val;

  if (isnan(a) || isnan(scale)) {
    THError("gamma2: shape and scale should be valid positive numbers.");
    return nan("");
  }
  if (a <= 0.0 || scale <= 0.0) {
    if(scale == 0. || a == 0.) return 0.;
    THError("gamma2: shape and scale should be positive values.");
    return nan("");
  }

  if (a < 1.) { /* GS algorithm for parameters a < 1 */
    e = 1.0 + exp_m1 * a;
    for(;;) {
      p = e * __uniform__(state);
      if (p >= 1.0) {
        x = -log((e - p) / a);
        if (THRandom_exponential(state, 1) >= (1.0 - a) * log(x))
          break;
      } else {
        x = exp(log(p) / a);
        if (THRandom_exponential(state, 1) >= x)
          break;
      }
    }
    return scale * x;
  }

  /* --- a >= 1 : GD algorithm --- */

  /* Step 1: Recalculations of s2, s, d if a has changed */
  if (a != aa) {
    aa = a;
    s2 = a - 0.5;
    s = sqrt(s2);
    d = sqrt32 - s * 12.0;
  }
  /* Step 2: t = standard normal deviate,
     x = (s,1/2) -normal deviate. */

  /* immediate acceptance (i) */
  t = THRandom_normal(state, 0, 1);
  x = s + 0.5 * t;
  ret_val = x * x;
  if (t >= 0.0)
    return scale * ret_val;

  /* Step 3: u = 0,1 - uniform sample. squeeze acceptance (s) */
  u = __uniform__(state);
  if (d * u <= t * t * t)
    return scale * ret_val;

  /* Step 4: recalculations of q0, b, si, c if necessary */

  if (a != aaa) {
    aaa = a;
    r = 1.0 / a;
    q0 = ((((((q7 * r + q6) * r + q5) * r + q4) * r + q3) * r
           + q2) * r + q1) * r;

      /* Approximation depending on size of parameter a */
      /* The constants in the expressions for b, si and c */
      /* were established by numerical experiments */

      if (a <= 3.686) {
        b = 0.463 + s + 0.178 * s2;
        si = 1.235;
        c = 0.195 / s - 0.079 + 0.16 * s;
      } else if (a <= 13.022) {
        b = 1.654 + 0.0076 * s2;
        si = 1.68 / s + 0.275;
        c = 0.062 / s + 0.024;
      } else {
        b = 1.77;
        si = 0.75;
        c = 0.1515 / s;
      }
  }
  /* Step 5: no quotient test if x not positive */

  if (x > 0.0) {
    /* Step 6: calculation of v and quotient q */
    v = t / (s + s);
    if (fabs(v) <= 0.25)
      q = q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v
                                + a3) * v + a2) * v + a1) * v;
    else
      q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);


      /* Step 7: quotient acceptance (q) */
      if (log(1.0 - u) <= q)
        return scale * ret_val;
  }

  for(;;) {
    /* Step 8: e = standard exponential deviate
     *	u =  0,1 -uniform deviate
     *	t = (b,si)-double exponential (laplace) sample */
    e = THRandom_exponential(state, 1);
    u = __uniform__(state);
    u = u + u - 1.0;
    if (u < 0.0)
      t = b - si * e;
    else
      t = b + si * e;
    /* Step	 9:  rejection if t < tau(1) = -0.71874483771719 */
    if (t >= -0.71874483771719) {
      /* Step 10:	 calculation of v and quotient q */
      v = t / (s + s);
      if (fabs(v) <= 0.25)
        q = q0 + 0.5 * t * t *
          ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v
            + a2) * v + a1) * v;
      else
        q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
        /* Step 11:	 hat acceptance (h) */
        /* (if q not positive go to step 8) */
        if (q > 0.0) {
          w = expm1(q);
          /*  ^^^^^ original code had approximation with rel.err < 2e-7 */
          /* if t is rejected sample again at step 8 */
          if (c * fabs(u) <= w * exp(e - 0.5 * t * t))
            break;
        }
    }
  } /* repeat .. until  `t' is accepted */
  x = s + 0.5 * t;
  return scale * x * x;
}

double THRandom_beta(THGenerator *_generator, double a, double b)
{
  double output;

  if (a <= 1 && b <= 1) {
    double U, V, X, Y;

    while (1) {
      U = __uniform__(_generator);
      V = __uniform__(_generator);
      X = pow(U, 1/a);
      Y = pow(V, 1/a);

      if ((X + Y) <= 1) {
        if ((X + Y) > 0) {
          output = X / (X + Y);
          break;
        }
        else {
          double logX = log(U)/a;
          double logY = log(V)/b;
          double logM = logX > logY ? logX : logY;
          logX -= logM;
          logY -= logM;

          output = exp(logX - log(exp(logX) + exp(logY)));
          break;
        }
      }
    }
  }
  else {
    double Ga = THRandom_gamma(_generator, a);
    double Gb = THRandom_gamma(_generator, b);
    output = Ga / (Ga + Gb);
  }

  return output;
}

int THRandom_geometric(THGenerator *_generator, double p)
{
  THArgCheck(p > 0 && p < 1, 1, "must be > 0 and < 1");
  return((int)(log(1-__uniform__(_generator)) / log(p)) + 1);
}

int THRandom_bernoulli(THGenerator *_generator, double p)
{
  THArgCheck(p >= 0 && p <= 1, 1, "must be >= 0 and <= 1");
  return(__uniform__(_generator) <= p);
}

int THRandom_binomial(THGenerator *_generator, int nin, double pp)
{
  /* FIXME: These should become THREAD_specific globals : */

  static double c, fm, npq, p1, p2, p3, p4, qn;
  static double xl, xll, xlr, xm, xr;

  static double psave = -1.0;
  static int nsave = -1;
  static int mv;

  double f, f1, f2, u, v, w, w2, x, x1, x2, z, z2;
  double p, q, np, g, r, al, alv, amaxp, ffm, ynorm;
  int i, ix, k, nn;

  THArgCheck(pp >= 0 && pp <= 1, 1, "must be >= 0 and <= 1");

  r = nin;

  if (r == 0 || pp == 0.) return 0;
  if (pp == 1.) return r;

  if (r >= INT_MAX)/* evade integer overflow,
                      and r == INT_MAX gave only even values */ {
    THError("binomial: too big n");
    return (int)nan("");
  }
  /* else */
  nn = (int) r;

  p = fmin(pp, 1. - pp);
  q = 1. - p;
  np = nn * p;
  r = p / q;
  g = r * (nn + 1);

  /* Setup, perform only when parameters change [using static (globals): */

  /* FIXING: Want this thread safe
     -- use as little (thread globals) as possible
  */
  if (pp != psave || nn != nsave) {
    psave = pp;
    nsave = nn;
    if (np < 30.0) {
      /* inverse cdf logic for mean less than 30 */
      qn = pow(q, nn);
      goto L_np_small;
    } else {
      ffm = np + p;
      mv = (int) ffm;
      fm = mv;
      npq = np * q;
      p1 = (int)(2.195 * sqrt(npq) - 4.6 * q) + 0.5;
      xm = fm + 0.5;
      xl = xm - p1;
      xr = xm + p1;
      c = 0.134 + 20.5 / (15.3 + fm);
      al = (ffm - xl) / (ffm - xl * p);
      xll = al * (1.0 + 0.5 * al);
      al = (xr - ffm) / (xr * q);
      xlr = al * (1.0 + 0.5 * al);
      p2 = p1 * (1.0 + c + c);
      p3 = p2 + c / xll;
      p4 = p3 + c / xlr;
    }
  } else if (nn == nsave) {
    if (np < 30.0)
      goto L_np_small;
  }

  /*-------------------------- np = n*p >= 30 : ------------------- */
  for(;;) {
    u = __uniform__(_generator) * p4;
    v = __uniform__(_generator);
    /* triangular region */
    if (u <= p1) {
      ix = (int)(xm - p1 * v + u);
      goto finis;
    }
    /* parallelogram region */
    if (u <= p2) {
      x = xl + (u - p1) / c;
      v = v * c + 1.0 - fabs(xm - x) / p1;
      if (v > 1.0 || v <= 0.)
        continue;
      ix = (int) x;
    } else {
      if (u > p3) {	/* right tail */
        ix = (int)(xr - log(v) / xlr);
        if (ix > nn)
          continue;
        v = v * (u - p3) * xlr;
      } else {/* left tail */
        ix = (int)(xl + log(v) / xll);
        if (ix < 0)
          continue;
        v = v * (u - p2) * xll;
      }
    }
    /* determine appropriate way to perform accept/reject test */
    k = abs(ix - mv);
    if (k <= 20 || k >= npq / 2 - 1) {
      /* explicit evaluation */
      f = 1.0;
      if (mv < ix) {
        for (i = mv + 1; i <= ix; i++)
          f *= (g / i - r);
      } else if (mv != ix) {
        for (i = ix + 1; i <= mv; i++)
          f /= (g / i - r);
      }
      if (v <= f)
        goto finis;
    } else {
      /* squeezing using upper and lower bounds on log(f(x)) */
      amaxp = (k / npq) * ((k * (k / 3. + 0.625) + 0.1666666666666) / npq + 0.5);
      ynorm = -k * k / (2.0 * npq);
      alv = log(v);
      if (alv < ynorm - amaxp)
        goto finis;
      if (alv <= ynorm + amaxp) {
        /* stirling's formula to machine accuracy */
        /* for the final acceptance/rejection test */
        x1 = ix + 1;
        f1 = fm + 1.0;
        z = nn + 1 - fm;
        w = nn - ix + 1.0;
        z2 = z * z;
        x2 = x1 * x1;
        f2 = f1 * f1;
        w2 = w * w;
        if (alv <= xm * log(f1 / x1) + (n - mv + 0.5) * log(z / w) + (ix - mv) * log(w * p / (x1 * q)) + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / f2) / f2) / f2) / f2) / f1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / z2) / z2) / z2) / z2) / z / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / x2) / x2) / x2) / x2) / x1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / w2) / w2) / w2) / w2) / w / 166320.)
          goto finis;
      }
    }
  }

 L_np_small:
  /*---------------------- np = n*p < 30 : ------------------------- */

  for(;;) {
    ix = 0;
    f = qn;
    u = __uniform__(_generator);
    for(;;) {
      if (u < f)
        goto finis;
      if (ix > 110)
        break;
      u -= f;
      ix++;
      f *= (g / ix - r);
    }
  }
 finis:
  if (psave > 0.5)
    ix = n - ix;
  return (int)ix;
}

double afc(int i)
{
  // If (i > 7), use Stirling's approximation, otherwise use table lookup.
  const static double al[8] =
    {
      0.0,/*ln(0!)=ln(1)*/
      0.0,/*ln(1!)=ln(1)*/
      0.69314718055994530941723212145817,/*ln(2) */
      1.79175946922805500081247735838070,/*ln(6) */
      3.17805383034794561964694160129705,/*ln(24)*/
      4.78749174278204599424770093452324,
      6.57925121201010099506017829290394,
      8.52516136106541430016553103634712
      /* 10.60460290274525022841722740072165, approx. value below =
         10.6046028788027; rel.error = 2.26 10^{-9}

         FIXME: Use constants and if(n > ..) decisions from ./stirlerr.c
         -----  will be even *faster* for n > 500 (or so)
      */
    };

  if (i < 0) {
    THError("afc: i should be positive number.");
    return -1; // unreached
  }
  if (i <= 7)
    return al[i];
  // else i >= 8 :
  double di = i, i2 = di*di;
  return (di + 0.5) * log(di) - di + 0.918938533204672741780329736406 +
    (0.0833333333333333 - 0.00277777777777778 / i2) / di;
}


int THRandom_hypergeometric(THGenerator *_generator, int nn1in, int nn2in, int kkin)
{
  int nn1, nn2, kk;
  int ix; // return value (coerced to double at the very end)
  char setup1, setup2;

  /* These should become 'thread_local globals' : */
  static int ks = -1, n1s = -1, n2s = -1;
  static int mv, minjx, maxjx;
  static int k, n1, n2; // <- not allowing larger integer par
  static double tn;

  // II :
  static double w;
  // III:
  static double a, d, s, xl, xr, kl, kr, lamdl, lamdr, p1, p2, p3;

  /* check parameter validity */

  if (nn1in < 0 || nn2in < 0 || kkin < 0 || kkin > nn1in + nn2in){
    THError("hypergeometric: nr, nb, k should be positive; k <= nr + nb");
    return (int)nan("");
  }

  nn1 = (int)nn1in;
  nn2 = (int)nn2in;
  kk  = (int)kkin;

  /* if new parameter values, initialize */
  if (nn1 != n1s || nn2 != n2s) {
    setup1 = 1;	setup2 = 1;
  } else if (kk != ks) {
    setup1 = 0;	setup2 = 1;
  } else {
    setup1 = 0;	setup2 = 0;
  }
  if (setup1) {
    n1s = nn1;
    n2s = nn2;
    tn = nn1 + nn2;
    if (nn1 <= nn2) {
      n1 = nn1;
      n2 = nn2;
    } else {
      n1 = nn2;
      n2 = nn1;
    }
  }
  if (setup2) {
    ks = kk;
    if (kk + kk >= tn) {
      k = (int)(tn - kk);
    } else {
      k = kk;
    }
  }
  if (setup1 || setup2) {
    mv = (int) ((k + 1.) * (n1 + 1.) / (tn + 2.));
    minjx = 0 < k - n2 ? k - n2 : 0;
    maxjx = n1 < k ? n1 : k;
  }
  /* generate random variate --- Three basic cases */

  if (minjx == maxjx) { /* I: degenerate distribution ---------------- */
    ix = maxjx;
    goto L_finis; // return appropriate variate

  } else if (mv - minjx < 10) { // II: (Scaled) algorithm HIN (inverse transformation) ----
    const static double scale = 1e25; // scaling factor against (early) underflow
    const static double con = 57.5646273248511421;
    // 25*log(10) = log(scale) { <==> exp(con) == scale }
    if (setup1 || setup2) {
      double lw; // log(w);  w = exp(lw) * scale = exp(lw + log(scale)) = exp(lw + con)
      if (k < n2) {
        lw = afc(n2) + afc(n1 + n2 - k) - afc(n2 - k) - afc(n1 + n2);
      } else {
        lw = afc(n1) + afc(     k     ) - afc(k - n2) - afc(n1 + n2);
      }
      w = exp(lw + con);
    }
    double p, u;
  L10:
    p = w;
    ix = minjx;
    u = __uniform__(_generator) * scale;
    while (u > p) {
      u -= p;
      p *= ((double) n1 - ix) * (k - ix);
      ix++;
      p = p / ix / (n2 - k + ix);
      if (ix > maxjx)
        goto L10;
      // FIXME  if(p == 0.)  we also "have lost"  => goto L10
    }
  } else { /* III : H2PE Algorithm --------------------------------------- */

    double u,v;

    if (setup1 || setup2) {
      s = sqrt((tn - k) * k * n1 * n2 / (tn - 1) / tn / tn);

      /* remark: d is defined in reference without int. */
      /* the truncation centers the cell boundaries at 0.5 */

      d = (int) (1.5 * s) + .5;
      xl = mv - d + .5;
      xr = mv + d + .5;
      a = afc(mv) + afc(n1 - mv) + afc(k - mv) + afc(n2 - k + mv);
      kl = exp(a - afc((int) (xl)) - afc((int) (n1 - xl))
               - afc((int) (k - xl))
               - afc((int) (n2 - k + xl)));
      kr = exp(a - afc((int) (xr - 1))
               - afc((int) (n1 - xr + 1))
               - afc((int) (k - xr + 1))
               - afc((int) (n2 - k + xr - 1)));
      lamdl = -log(xl * (n2 - k + xl) / (n1 - xl + 1) / (k - xl + 1));
      lamdr = -log((n1 - xr + 1) * (k - xr + 1) / xr / (n2 - k + xr));
      p1 = d + d;
      p2 = p1 + kl / lamdl;
      p3 = p2 + kr / lamdr;
    }
    int n_uv = 0;
  L30:
    u = __uniform__(_generator) * p3;
    v = __uniform__(_generator);
    n_uv++;
    if(n_uv >= 10000) {
      THError("hypergeometric: too much rejections.");
      return (int)nan("");
    }

    if (u < p1) {		/* rectangular region */
      ix = (int) (xl + u);
    } else if (u <= p2) {	/* left tail */
      ix = (int) (xl + log(v) / lamdl);
      if (ix < minjx)
        goto L30;
      v = v * (u - p1) * lamdl;
    } else {		/* right tail */
      ix = (int) (xr - log(v) / lamdr);
      if (ix > maxjx)
        goto L30;
      v = v * (u - p2) * lamdr;
    }

    /* acceptance/rejection test */
    char reject = 1;

    if (mv < 100 || ix <= 50) {
      /* explicit evaluation */
      /* The original algorithm (and TOMS 668) have
         f = f * i * (n2 - k + i) / (n1 - i) / (k - i);
         in the (m > ix) case, but the definition of the
         recurrence relation on p134 shows that the +1 is
         needed. */
      int i;
      double f = 1.0;
      if (mv < ix) {
        for (i = mv + 1; i <= ix; i++)
          f = f * (n1 - i + 1) * (k - i + 1) / (n2 - k + i) / i;
      } else if (mv > ix) {
        for (i = ix + 1; i <= mv; i++)
          f = f * i * (n2 - k + i) / (n1 - i + 1) / (k - i + 1);
      }
      if (v <= f) {
        reject = 0;
      }
    } else {

      const static double deltal = 0.0078;
      const static double deltau = 0.0034;

      double e, g, r, t, y;
      double de, dg, dr, ds, dt, gl, gu, nk, nm, ub;
      double xk, xm, xn, y1, ym, yn, yk, alv;

      /* squeeze using upper and lower bounds */
      y = ix;
      y1 = y + 1.0;
      ym = y - mv;
      yn = n1 - y + 1.0;
      yk = k - y + 1.0;
      nk = n2 - k + y1;
      r = -ym / y1;
      s = ym / yn;
      t = ym / yk;
      e = -ym / nk;
      g = yn * yk / (y1 * nk) - 1.0;
      dg = 1.0;
      if (g < 0.0)
        dg = 1.0 + g;
      gu = g * (1.0 + g * (-0.5 + g / 3.0));
      gl = gu - .25 * (g * g * g * g) / dg;
      xm = mv + 0.5;
      xn = n1 - mv + 0.5;
      xk = k - mv + 0.5;
      nm = n2 - k + xm;
      ub = y * gu - mv * gl + deltau
        + xm * r * (1. + r * (-0.5 + r / 3.0))
        + xn * s * (1. + s * (-0.5 + s / 3.0))
        + xk * t * (1. + t * (-0.5 + t / 3.0))
        + nm * e * (1. + e * (-0.5 + e / 3.0));
        /* test against upper bound */
        alv = log(v);
        if (alv > ub) {
          reject = 1;
        } else {
          /* test against lower bound */
          dr = xm * (r * r * r * r);
          if (r < 0.0)
            dr /= (1.0 + r);
          ds = xn * (s * s * s * s);
          if (s < 0.0)
            ds /= (1.0 + s);
          dt = xk * (t * t * t * t);
          if (t < 0.0)
            dt /= (1.0 + t);
          de = nm * (e * e * e * e);
          if (e < 0.0)
            de /= (1.0 + e);
          if (alv < ub - 0.25 * (dr + ds + dt + de)
              + (y + mv) * (gl - gu) - deltal) {
            reject = 0;
          }
          else {
            /* * Stirling's formula to machine accuracy
             */
            if (alv <= (a - afc(ix) - afc(n1 - ix)
                        - afc(k - ix) - afc(n2 - k + ix))) {
              reject = 0;
            } else {
              reject = 1;
            }
          }
        }
    } // else
    if (reject)
      goto L30;
  }


 L_finis:
  /* return appropriate variate */

  if (kk + kk >= tn) {
    if (nn1 > nn2) {
      ix = kk - nn2 + ix;
    } else {
      ix = nn1 - ix;
    }
  } else {
    if (nn1 > nn2)
      ix = kk - ix;
  }
  return ix;

}

int THRandom_poisson(THGenerator *_generator, double mu)
{
#define a0	-0.5
#define a1	 0.3333333
#define a2	-0.2500068
#define a3	 0.2000118
#define a4	-0.1661269
#define a5	 0.1421878
#define a6	-0.1384794
#define a7	 0.1250060

#define one_7	0.1428571428571428571
#define one_12	0.0833333333333333333
#define one_24	0.0416666666666666667

  /* Factorial Table (0:9)! */
  const static double fact[10] =
    {
      1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880.
    };

  /* These are static --- persistent between calls for same mu : */
  static int l, mv;

  static double b1, b2, c, c0, c1, c2, c3;
  static double pp[36], p0, p, q, s, d, omega;
  static double big_l;/* integer "w/o overflow" */
  static double muprev = 0., muprev2 = 0.;/*, muold	 = 0.*/

  /* Local Vars  [initialize some for -Wall]: */
  double del, difmuk= 0., E= 0., fk= 0., fx, fy, g, px, py, t, u= 0., v, x;
  double pois = -1.;
  int k, kflag, big_mu, new_big_mu = 0;

  if (mu <= 0.)
    return 0;

  big_mu = mu >= 10.;
  if(big_mu)
    new_big_mu = 0;

  if (!(big_mu && mu == muprev)) {/* maybe compute new persistent par.s */

    if (big_mu) {
      new_big_mu = 1;
      /* Case A. (recalculation of s,d,l	because mu has changed):
       * The poisson probabilities pk exceed the discrete normal
       * probabilities fk whenever k >= m(mu).
       */
      muprev = mu;
      s = sqrt(mu);
      d = 6. * mu * mu;
      big_l = floor(mu - 1.1484);
      /* = an upper bound to m(mu) for all mu >= 10.*/
    }
    else { /* Small mu ( < 10) -- not using normal approx. */

      /* Case B. (start new table and calculate p0 if necessary) */

      /*muprev = 0.;-* such that next time, mu != muprev ..*/
      if (mu != muprev) {
        muprev = mu;
        mv = (1 > mu ? 1 : (int)mu);
        l = 0; /* pp[] is already ok up to pp[l] */
        q = p0 = p = exp(-mu);
      }

      for(;;) {
        /* Step U. uniform sample for inversion method */
        u = __uniform__(_generator);
        if (u <= p0)
          return 0.;

        /* Step T. table comparison until the end pp[l] of the
           pp-table of cumulative poisson probabilities
           (0.458 > ~= pp[9](= 0.45792971447) for mu=10 ) */
        if (l != 0) {
          for (k = (u <= 0.458) ? 1 : (l < mv ? l : mv);  k <= l; k++)
            if (u <= pp[k])
              return k;
          if (l == 35) /* u > pp[35] */
            continue;
        }
        /* Step C. creation of new poisson
           probabilities p[l..] and their cumulatives q =: pp[k] */
        l++;
        for (k = l; k <= 35; k++) {
          p *= mu / k;
          q += p;
          pp[k] = q;
          if (u <= q) {
            l = k;
            return k;
          }
        }
        l = 35;
      } /* end(repeat) */
    }/* mu < 10 */

  } /* end {initialize persistent vars} */

  /* Only if mu >= 10 : ----------------------- */

  /* Step N. normal sample */
  g = mu + s * THRandom_normal(_generator, 0, 1);/* norm_rand() ~ N(0,1), standard normal */

  if (g >= 0.) {
    pois = floor(g);
    /* Step I. immediate acceptance if pois is large enough */
    if (pois >= big_l)
      return pois;
    /* Step S. squeeze acceptance */
    fk = pois;
    difmuk = mu - fk;
    u = __uniform__(_generator); /* ~ U(0,1) - sample */
    if (d * u >= difmuk * difmuk * difmuk)
      return pois;
  }

  /* Step P. preparations for steps Q and H.
     (recalculations of parameters if necessary) */

  if (new_big_mu || mu != muprev2) {
    /* Careful! muprev2 is not always == muprev
       because one might have exited in step I or S
    */
    muprev2 = mu;
    omega = 0.398942280401432677939946059934 / s;
    /* The quantities b1, b2, c3, c2, c1, c0 are for the Hermite
     * approximations to the discrete normal probabilities fk. */

    b1 = one_24 / mu;
    b2 = 0.3 * b1 * b1;
    c3 = one_7 * b1 * b2;
    c2 = b2 - 15. * c3;
    c1 = b1 - 6. * b2 + 45. * c3;
    c0 = 1. - b1 + 3. * b2 - 15. * c3;
    c = 0.1069 / mu; /* guarantees majorization by the 'hat'-function. */
  }

  if (g >= 0.) {
    /* 'Subroutine' F is called (kflag=0 for correct return) */
    kflag = 0;
    goto Step_F;
  }


  for(;;) {
    /* Step E. Exponential Sample */

    E = THRandom_exponential(_generator, 1);	/* ~ Exp(1) (standard exponential) */

    /*  sample t from the laplace 'hat'
        (if t <= -0.6744 then pk < fk for all mu >= 10.) */
    u = 2 * __uniform__(_generator) - 1.;
    t = 1.8 + ((u >= 0) ? fabs(E) : -fabs(E));
    if (t > -0.6744) {
      pois = floor(mu + s * t);
      fk = pois;
      difmuk = mu - fk;

      /* 'subroutine' F is called (kflag=1 for correct return) */
      kflag = 1;

    Step_F: /* 'subroutine' F : calculation of px,py,fx,fy. */

      if (pois < 10) { /* use factorials from table fact[] */
        px = -mu;
        py = pow(mu, pois) / fact[(int)pois];
      }
      else {
        /* Case pois >= 10 uses polynomial approximation
           a0-a7 for accuracy when advisable */
        del = one_12 / fk;
        del = del * (1. - 4.8 * del * del);
        v = difmuk / fk;
        if (fabs(v) <= 0.25)
          px = fk * v * v * (((((((a7 * v + a6) * v + a5) * v + a4) *
                                v + a3) * v + a2) * v + a1) * v + a0)
            - del;
          else /* |v| > 1/4 */
            px = fk * log(1. + v) - difmuk - del;
            py = 0.398942280401432677939946059934 / sqrt(fk);
      }
      x = (0.5 - difmuk) / s;
      x *= x;/* x^2 */
      fx = -0.5 * x;
      fy = omega * (((c3 * x + c2) * x + c1) * x + c0);
      if (kflag > 0) {
        /* Step H. Hat acceptance (E is repeated on rejection) */
        if (c * fabs(u) <= py * exp(px + E) - fy * exp(fx + E))
          break;
      } else
        /* Step Q. Quotient acceptance (rare case) */
        if (fy - u * fy <= py * exp(px - fx))
          break;
    }/* t > -.67.. */
  }
  return pois;
}
