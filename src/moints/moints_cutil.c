
#if defined(CRAY_T3E) || defined(CRAY_T3D)
int ONBITMASK( int *len )
#else
int onbitmask_( int *len )
#endif
{
  unsigned int		mask;

  mask = ~((~0) << *len);
  return ((int)mask);
}

  

