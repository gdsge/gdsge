// normal pdf, erfc and normal cdf
template<class T>
T normpdf(T x, double mu, double sigma) 
{
    T r = 0.3989422804/sigma*exp(-.5*(x-mu)*(x-mu)/(sigma*sigma));
    return r;
};


adouble kf_erfc(adouble x)
{
  const double p0 = 220.2068679123761;
  const double p1 = 221.2135961699311;
  const double p2 = 112.0792914978709;
  const double p3 = 33.912866078383;
  const double p4 = 6.37396220353165;
  const double p5 = .7003830644436881;
  const double p6 = .03526249659989109;
  const double q0 = 440.4137358247522;
  const double q1 = 793.8265125199484;
  const double q2 = 637.3336333788311;
  const double q3 = 296.5642487796737;
  const double q4 = 86.78073220294608;
  const double q5 = 16.06417757920695;
  const double q6 = 1.755667163182642;
  const double q7 = .08838834764831844;
  adouble expntl, z, p;
  z = fabs(x) * 1.41421356237;
  // if (z > 37.) return x > 0.? 0. : 2.;
  if (z>37)
  {
    if (x>0)
    {
      adouble r=0.0;
      return r;
    }
    else
    {
      adouble r=2.0;
      return r;
    }
  }
  expntl = exp(z * z * - .5);
  if (z < 10. / 1.41421356237) // for small z
      p = expntl * ((((((p6 * z + p5) * z + p4) * z + p3) * z + p2) * z + p1) * z + p0)
      / (((((((q7 * z + q6) * z + q5) * z + q4) * z + q3) * z + q2) * z + q1) * z + q0);
  else p = expntl / 2.506628274631001 / (z + 1. / (z + 2. / (z + 3. / (z + 4. / (z + .65)))));
  if (x>0)
  {
    adouble r = p*2.;
    return r;
  }
  else
  {
    adouble r = (1-p)*2.;
    return r;
  }
};

double kf_erfc(double x)
{
  const double p0 = 220.2068679123761;
  const double p1 = 221.2135961699311;
  const double p2 = 112.0792914978709;
  const double p3 = 33.912866078383;
  const double p4 = 6.37396220353165;
  const double p5 = .7003830644436881;
  const double p6 = .03526249659989109;
  const double q0 = 440.4137358247522;
  const double q1 = 793.8265125199484;
  const double q2 = 637.3336333788311;
  const double q3 = 296.5642487796737;
  const double q4 = 86.78073220294608;
  const double q5 = 16.06417757920695;
  const double q6 = 1.755667163182642;
  const double q7 = .08838834764831844;
  double expntl, z, p;
  z = fabs(x) * 1.41421356237;
  // if (z > 37.) return x > 0.? 0. : 2.;
  if (z>37)
  {
    if (x>0)
      return 0.;
    else
      return 2.;
  }
  expntl = exp(z * z * - .5);
  if (z < 10. / 1.41421356237) // for small z
      p = expntl * ((((((p6 * z + p5) * z + p4) * z + p3) * z + p2) * z + p1) * z + p0)
      / (((((((q7 * z + q6) * z + q5) * z + q4) * z + q3) * z + q2) * z + q1) * z + q0);
  else p = expntl / 2.506628274631001 / (z + 1. / (z + 2. / (z + 3. / (z + 4. / (z + .65)))));
  if (x>0)
  {
    double r = p*2.;
    return r;
  }
  else
  {
    double r = (1-p)*2.;
    return r;
  }
};

template<class T>
T normcdf(T x, double mu, double sigma)
{
  T r = .5*kf_erfc(-(x-mu)/(1.41421356237*sigma));
  return r;
};

// Log gamma function from https://github.com/lh3/samtools/blob/master/bcftools/kfunc.c
inline double kf_lgamma(double z)
{
    double x = 0;
    x += 0.1659470187408462e-06 / (z+7);
    x += 0.9934937113930748e-05 / (z+6);
    x -= 0.1385710331296526     / (z+5);
    x += 12.50734324009056      / (z+4);
    x -= 176.6150291498386      / (z+3);
    x += 771.3234287757674      / (z+2);
    x -= 1259.139216722289      / (z+1);
    x += 676.5203681218835      / z;
    x += 0.9999999999995183;
    double r = log(x) - 5.58106146679532777 - z + (z-0.5) * log(z+6.5);
    return r;
};


#define KF_GAMMA_EPS 1e-14
// Incomplete lower gamma function from https://github.com/lh3/samtools/blob/master/bcftools/kfunc.c
double _kf_gammap(double s, double z)
{
    double sum, x;
    int k;
    for (k = 1, sum = x = 1.; k < 100; ++k) {
        sum += (x *= z / (s + k));
        if (x / sum < KF_GAMMA_EPS) break;
    }
    double r = exp(s * log(z) - z - kf_lgamma(s + 1.) + log(sum));
    return r;
};

adouble _kf_gammap(double s, adouble z)
{
    adouble sum, x;
    int k;
    for (k = 1, sum = x = 1.; k < 100; ++k) {
        sum += (x *= z / (s + k));
        if (x / sum < KF_GAMMA_EPS) break;
    }
    adouble r = exp(s * log(z) - z - kf_lgamma(s + 1.) + log(sum));
    return r;
};

// gamma pdf from wiki https://en.wikipedia.org/wiki/Gamma_distribution
double gampdf(double xx, double alpha, double scale)
{
    double beta = 1/scale;
    double r = pow(beta,alpha)/exp(kf_lgamma(alpha)) * pow(xx,alpha-1)*exp(-beta*xx);
    return r;
};

adouble gampdf(adouble xx, double alpha, double scale)
{
    double beta = 1/scale;
    adouble r = pow(beta,alpha)/exp(kf_lgamma(alpha)) * pow(xx,alpha-1)*exp(-beta*xx);
    return r;
};

// gamma cdf from wiki https://en.wikipedia.org/wiki/Gamma_distribution
double gamcdf(double xx, double alpha, double scale)
{
    double beta = 1/scale;
    double r = _kf_gammap(alpha,beta*xx);
    return r;
};

adouble gamcdf(adouble xx, double alpha, double scale)
{
    double beta = 1/scale;
    adouble r = _kf_gammap(alpha,beta*xx);
    return r;
};