#ifndef _CORRELATION_H
#define _CORRELATION_H
#include <fftw3.h>

//generate noise with positive with correlation function
double* lamda_sqrt_t(const double *alpha_r, 
		const int nt, const double dt); 
int corr_rand_gener_bare(const int nt, double** rand, double *lamda);
int corr_rand_gener_prep(const int nt);
int corr_rand_gener(const int nt,double **rand, double *lamda, int flag);
int corr_rand_gener_free();
int whitenoise(double sig,double *u1, double *u2);
extern fftw_plan fftw3_plan_for_corr_rand;
extern fftw_complex *fftw3_in_for_corr_rand;
extern fftw_complex *fftw3_out_for_corr_rand;

/*
How to use:
    for positive definite spectrum (will be explained later)
    1. Calculate the memory function and put the results into an array pointed by alpha_r, i.e., alpha_r[i]=memory(t_i),t_i=\delta_t*i.
    2. Put alpha_r into lamda_sqrt_t(const double *alpha,const int nt, const double dt), where the spectrum (the Fourier transform of the correlation matrix is diagonal, and positve definite, the square root of diagonal elements collected in a array are called spectrum here) of alpha_r[i] is calculated and put into a array pointed by the return value of the function.
    3. Call corr_rand_gener_pre(const int nt), which alloc the memory for fft algorithm. 
    4. Call corr_rand_gener(const int nt,double **rand, double *lamda, int flag);the colored noise will be generated and put into an array allocted in step 3 and pointed by fftw3_out_for_corr_rand. lamda is the spectrum array returned by the last step.
    5. Call corr_rand_gener_free to free the memory.
    
    the subroutine corr_rand_gener() generates 4 random vector per time (because of the algorithm). So there is a flag, if flag==0, the function genrate new 4 random vectors. Else, the function
    simply return the vector generated before. 
    
example:
int main(int argc, char *argv[])
{
  double alpha[1024];
  double *rand;
  for(int i=0;i<1024;i++)
  alpha[i]=exp(-dt*i);
  double *lamda=lamda_sqrt_t(alpha_r, 1024, 0.01);
  corr_rand_gener_prep(1024);
  corr_rand_gener(1024,&rand, lamda, 0);
}
*/

#endif
