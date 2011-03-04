#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "rnruni.h"
#include "correlation.h"
#define _FFTW3_TINY 1.0e-5



double* lamda_sqrt_t(const double *alpha_r,
                     const int nt, const double dt)
{
  fftw_complex *lamda=
    (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*2*nt);

  fftw_complex *corr=
    (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*2*nt);

  double *rlt = (double*)malloc(sizeof(double)*2*nt);

  assert( lamda && corr && rlt);

  bzero(lamda,sizeof(fftw_complex)*2*nt);
  bzero(corr,sizeof(fftw_complex)*2*nt);
  bzero(rlt,sizeof(double)*2*nt);

  fftw_plan fftpln = fftw_plan_dft_1d (2*nt, corr, lamda,
                                       FFTW_FORWARD, FFTW_MEASURE);

  int i;
  for(i=0;i<nt;i++)
    corr[i]=alpha_r[i];
  i=nt;
  corr[i]=alpha_r[i-1]/2.0;
  for(i=nt+1;i<2*nt;i++)
    corr[i]=alpha_r[2*nt-i];

  fftw_execute(fftpln);
//   FILE *fp=fopen("lamda.dat","w");
  for(i=0;i<nt*2;i++)
  {
//     fprintf(fp,"%d\t%f+%fi\n",i,creal(lamda[i]),cimag(lamda[i]));
    if(creal(lamda[i])< -_FFTW3_TINY ||
        fabs(cimag(lamda[i]))> _FFTW3_TINY)
    {
       fprintf(stderr,"correlation not positive: %d\t%e\t%e\n",
               i,creal(lamda[i]),cimag(lamda[i]));
      lamda[i]=0.0+0.0*I;
    }
    if(creal(lamda[i])<0) lamda[i]=0.0+0.0*I;
    rlt[i] = sqrt(creal(lamda[i]));
  }
//   fclose(fp);
  fftw_free(lamda);
  fftw_free(corr);
  fftw_destroy_plan(fftpln);
  return rlt;
}

fftw_plan  fftw3_plan_for_corr_rand;
fftw_complex *fftw3_in_for_corr_rand;
fftw_complex *fftw3_out_for_corr_rand;

int corr_rand_gener_prep(const int nt)

{
#define plan fftw3_plan_for_corr_rand
#define in fftw3_in_for_corr_rand
#define out fftw3_out_for_corr_rand
  in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*2*nt);
  out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*2*nt);
  plan = fftw_plan_dft_1d (2*nt, in, out,
                           FFTW_BACKWARD, FFTW_MEASURE);
  assert(in&&out);
  return 0;
#undef plan
#undef in
#undef out
}

int corr_rand_gener_free()
{
#define plan fftw3_plan_for_corr_rand
#define in fftw3_in_for_corr_rand
#define out fftw3_out_for_corr_rand
  fftw_free(in);
  fftw_free(out);
  fftw_destroy_plan(plan);
  return 0;
#undef plan
#undef in
#undef out
}


int corr_rand_gener(const int nt, double** rand, double *lamda, int flag)
{
#define plan fftw3_plan_for_corr_rand
#define in fftw3_in_for_corr_rand
#define out fftw3_out_for_corr_rand
  double pt, ar, ai, u1, u2;
  double ivfct = sqrt(0.5/nt);
  double *ptr = (double*)in;
  double *ptt = (double*)out;

if(flag==0){
  int i;
  for(i=0; i<2*nt; i++)
  {
    if(lamda[i]<_FFTW3_TINY && lamda[i]>-_FFTW3_TINY)
    {
      ptr[2*i]=0;
      ptr[2*i+1]=0;
      continue;
    }
#define _NOISE_TYPE 1
//#dfine _NOISE_TYPE 2
//1 Gaussian noise
//2 binary noise
#if _NOISE_TYPE==1
    do
    {
      u1 = rnruni();
    }
    while (u1<=1.e-8);
    do
    {
      u2 = rnruni();
    }
    while (u2<=1.e-8);
    u2 *= 2*M_PI;
    pt = sqrt(-2.*log(u1))*ivfct*lamda[i];
    ar = pt*cos(u2);
    ai = pt*sin(u2);
    //		in[i]=ar+ai*I;
    ptr[2*i] = ar;
    ptr[2*i+1] = ai;
#elif _NOISE_TYPE==2
    double tip=ivfct*lamda[i];
    u1 = rnruni();
    u2 = rnruni();
    ptr[2*i]=(u1>=0.5)? tip:-tip;
    ptr[2*i+1]=(u2<=0.5)? -tip:tip;
#endif
  }

  fftw_execute(plan);
  for(i=0;i<2*nt;i++) ptr[i]=ptt[2*i];
  ptr=ptr+2*nt;
  for(i=0;i<2*nt;i++) ptr[i]=ptt[2*i+1];
  ptr=(double*) in;
  *rand=ptr;
  return 0;
}
else if(flag==1)  { *rand=ptr+nt; return 0;}
else if(flag==2)  { *rand=ptr+2*nt;return 0;}
else if(flag==3)  { *rand=ptr+3*nt; return 0;}

#undef plan
#undef in
#undef out
#undef _NOISE_TYPE
  return 0;
}






int corr_rand_gener_bare(const int nt, double** rand, double *lamda)
{
#define plan fftw3_plan_for_corr_rand
#define in fftw3_in_for_corr_rand
#define out fftw3_out_for_corr_rand
  double /*pt, ar, ai,*/ u1, u2;
  double ivfct = sqrt(0.5/nt);
  double *ptr = (double*)in;
//  double *ptt = (double*)out;

  int i;
  for(i=0; i<2*nt; i++)
  {
#define _TRUNCATE_MEMORY 0
//0:not truncate memory, 1: truncate memory
#if _TRUNCATE_MEMORY==1
    if(lamda[i]<_FFTW3_TINY && lamda[i]>-_FFTW3_TINY)
    {
      ptr[2*i]=0;
      ptr[2*i+1]=0;
      continue;
    }
#endif
#undef _TRUNCATE_MEMORY
#define _NOISE_TYPE 2
    //#define _NOISE_TYPE 2
    //1 Gaussian noise
    //2 binary noise
#if _NOISE_TYPE==1
    do
    {
      u1 = rnruni();
    }
    while (u1<=1.e-8);
    do
    {
      u2 = rnruni();
    }
    while (u2<=1.e-8);
    u2 *= 2*M_PI;
    pt = sqrt(-2.*log(u1))*ivfct*lamda[i];
    ar = pt*cos(u2);
    ai = pt*sin(u2);
    //		in[i]=ar+ai*I;
    ptr[2*i] = ar;
    ptr[2*i+1] = ai;
#elif _NOISE_TYPE==2
    double tip=ivfct*lamda[i];
    u1 = rnruni();
    u2 = rnruni();
    ptr[2*i]=(u1>=0.5)? tip:-tip;
    ptr[2*i+1]=(u2<=0.5)? -tip:tip;
#endif

  }

  fftw_execute(plan);
  *rand=(double*)out;
  return 0;

#undef plan
#undef in
#undef out
#undef _NOISE_TYPE
}

int whitenoise(double sig,double *r1, double *r2)
{
  double pt, u1, u2;
  do
  {
    u1 = rnruni();
  }
  while (u1<=1.e-8);
  do
  {
    u2 = rnruni();
  }
  while (u2<=1.e-8);
  u2 *= 2*M_PI;
  pt = sqrt(-2.*log(u1))*sig;
  *r1 = pt*cos(u2);
  *r2 = pt*sin(u2);
  return 0;
}

