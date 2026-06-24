
// Min
inline double MIN(double a, double b)
{
  if (a<b)
    return a;
  else
    return b;
}

inline adouble MIN(adouble aa, adouble bb)
{
  double a = aa.value();
  double b = bb.value();
  
  if (a<b)
    return aa;
  else
    return bb;
}

inline adouble MIN(adouble aa, double b)
{
  double a = aa.value();
  
  if (a<b)
  {
    return aa;
  }
  else
  {
    adouble bb = b;
    return bb;
  }
}

inline adouble MIN(double a, adouble bb)
{
  double b = bb.value();
  
  if (a<b)
  {
    adouble aa = a;
    return aa;
  }
  else
  {
    return bb;
  }
}

// Max
inline double MAX(double a, double b)
{
  if (a>b)
    return a;
  else
    return b;
}

inline adouble MAX(adouble aa, adouble bb)
{
  double a = aa.value();
  double b = bb.value();
  
  if (a>b)
    return aa;
  else
    return bb;
}

inline adouble MAX(adouble aa, double b)
{
  double a = aa.value();
  
  if (a>b)
  {
    return aa;
  }
  else
  {
    adouble bb = b;
    return bb;
  }
}

inline adouble MAX(double a, adouble bb)
{
  double b = bb.value();
  
  if (a>b)
  {
    adouble aa = a;
    return aa;
  }
  else
  {
    return bb;
  }
}


// Min
inline double min(double a, double b)
{
  if (a<b)
    return a;
  else
    return b;
}

inline adouble min(adouble aa, adouble bb)
{
  double a = aa.value();
  double b = bb.value();
  
  if (a<b)
    return aa;
  else
    return bb;
}

inline adouble min(adouble aa, double b)
{
  double a = aa.value();
  
  if (a<b)
  {
    return aa;
  }
  else
  {
    adouble bb = b;
    return bb;
  }
}

inline adouble min(double a, adouble bb)
{
  double b = bb.value();
  
  if (a<b)
  {
    adouble aa = a;
    return aa;
  }
  else
  {
    return bb;
  }
}

// Max
inline double max(double a, double b)
{
  if (a>b)
    return a;
  else
    return b;
}

inline adouble max(adouble aa, adouble bb)
{
  double a = aa.value();
  double b = bb.value();
  
  if (a>b)
    return aa;
  else
    return bb;
}

inline adouble max(adouble aa, double b)
{
  double a = aa.value();
  
  if (a>b)
  {
    return aa;
  }
  else
  {
    adouble bb = b;
    return bb;
  }
}

inline adouble max(double a, adouble bb)
{
  double b = bb.value();
  
  if (a>b)
  {
    adouble aa = a;
    return aa;
  }
  else
  {
    return bb;
  }
}
