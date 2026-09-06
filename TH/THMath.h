#ifndef _THMATH_H
#define _THMATH_H

static inline double TH_sigmoid(double value) {
  if (value < 0) {
    if (value < -744.44)
      return 0.0f;
    else {
      double k = exp(value);
      double r = k / (1 + k);
      return r;
    }
  }
  else {
    if (value > 36.7368)
      return 1.0f;
    else
        return 1.0 / (1.0 + exp(-value));
  }
}

static inline double TH_frac(double x) {
  return x - trunc(x);
}

static inline double TH_rsqrt(double x) {
  return 1.0 / sqrt(x);
}

static inline double TH_lerp(double a, double b, double weight) {
  return a + weight * (b-a);
}

static inline double TH_lgamma(double a) {
  if ((a >= 1) && (a <= 100)) {
    extern double _lgamma_cache[100];
    int na = round(a);
    if (na == a) {
      if (_lgamma_cache[na-1] <= 0) {
        _lgamma_cache[na-1] = lgamma(a);
      }
      return _lgamma_cache[na-1];
    }  else {
      return lgamma(a);
    }
  } else {
    return lgamma(a);
  }
}

static inline float TH_sigmoidf(float value) {
  if (value < 0) {
    if (value < -103.27)
      return 0.0f;
    else {
      float k = expf(value);
      float r = k / (1.0f + k);
      return r;
    }
  }
  else {
    if (value > 16.6355)
      return 1.0f;
    else
      return 1.0f / (1.0f + expf(-value));
  }
}

static inline float TH_fracf(float x) {
  return x - truncf(x);
}

static inline float TH_rsqrtf(float x) {
  return 1.0f / sqrtf(x);
}

static inline float TH_lerpf(float a, float b, float weight) {
  return a + weight * (b-a);
}

static inline float TH_lgammaf(float a) {
  if ((a >= 1) && (a <= 100)) {
    extern double _lgamma_cache[100];
    int na = round(a);
    if (na == a) {
      if (_lgamma_cache[na-1] <= 0) {
        _lgamma_cache[na-1] = lgamma(a);
      }
      return (float)_lgamma_cache[na-1];
    }  else {
      return lgammaf(a);
    }
  } else {
    return lgammaf(a);
  }
}

#endif // _THMATH_H
