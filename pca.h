/* Some parameters and functions needed by LFA, from GROMACS.
 * By Zhiyong Zhang, 02/02/2004.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "xdrfile.h"

#define XX 0
#define YY 1
#define ZZ 2
#define STRLEN 200
#define snew(ptr,nelem) (ptr)=save_calloc(#ptr,(nelem),sizeof(*(ptr)))
#define sfree(ptr) save_free(#ptr,(ptr))
#define sqr(x) ((x)*(x))

void *save_calloc(char *name,unsigned nelem,unsigned elsize);
void save_free(char *name,void *ptr);

extern char *fgets2(char *s, int n, FILE *stream);

static inline void cprod(const rvec a,const rvec b,rvec c)
{
  c[XX]=a[YY]*b[ZZ]-a[ZZ]*b[YY];
  c[YY]=a[ZZ]*b[XX]-a[XX]*b[ZZ];
  c[ZZ]=a[XX]*b[YY]-a[YY]*b[XX];
}

static inline void oprod(rvec a,rvec b,rvec c)
{
  c[XX]=a[YY]*b[ZZ]-a[ZZ]*b[YY];
  c[YY]=a[ZZ]*b[XX]-a[XX]*b[ZZ];
  c[ZZ]=a[XX]*b[YY]-a[YY]*b[XX];
}

//static inline float distance2(const rvec a,const rvec b)
static inline float distance2(rvec a,rvec b)
{
  return sqr(b[XX]-a[XX]) + sqr(b[YY]-a[YY]) + sqr(b[ZZ]-a[ZZ]);
}

static inline float iprod(const rvec a,const rvec b)
{
  return (a[XX]*b[XX]+a[YY]*b[YY]+a[ZZ]*b[ZZ]);
}

static inline float norm(const rvec a)
{
  return (float)sqrt(a[XX]*a[XX]+a[YY]*a[YY]+a[ZZ]*a[ZZ]);
}

static inline float gmx_angle(const rvec a, const rvec b)
{
    rvec w;
    float wlen,s;

    cprod(a,b,w);

    wlen  = norm(w);
    s     = iprod(a,b);

    return atan2(wlen,s);
}

static inline void clear_mat(matrix a)
{
  a[XX][XX]=a[XX][YY]=a[XX][ZZ]=0.0;
  a[YY][XX]=a[YY][YY]=a[YY][ZZ]=0.0;
  a[ZZ][XX]=a[ZZ][YY]=a[ZZ][ZZ]=0.0;
}
