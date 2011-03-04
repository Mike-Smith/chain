#ifndef _INTEGRATION_H
#define _INTEGRATION_H

#include <complex.h>

typedef _Complex double complex99;
//typedef double complex99;

//typedef void (data_ana)(double,int,complex99*,double*);
/*the parametre are, respectively
 t, the auto variable
 n, the number of dimentions
 y(t), the n dimentional variable at time t
 rlt, the interested value calculated from t,y[t]
 */
//typedef void (derive)(double,complex99*,int,complex99[],complex99[]);
/*
 t, the auto variable
 rdm, the noise at time t
 n, the number of dimentions
 y(t), the n dimentional variable at time t
 dydt(t), the derive of y(t);
 */

void euler_one_step(complex99 y[], complex99 dydx[], int n,
	double x, double h, complex99* rdm, complex99 yout[],
	void (*drv)(double,complex99*,int,complex99[],complex99[]));

void euler_propagator(complex99 vstart[], int nvar, double x1,
		double x2, int nstep, complex99 *noise, double *rlt,
		void (*drv)(double,complex99*,int,complex99[],complex99[]),
		void (*dt_an)(double,int,complex99*,double*)) ;

extern complex99 *rk_dym_d,*rk_dyt_d,*rk_yt_d;
extern int rk_kmax, rk_kount;
//temprory data for rkdumb
void rk4_one_step(complex99 y[], complex99 dydx[], int n,
		double x, double h, complex99* rdm, complex99 yout[],
		void (*derivs)(double,complex99*,int,complex99[],complex99[]));
int rk2_white ( double y[], int ny, double x1,
                      double x2, int nx, double *rlt, double *rdm, int nobv,
                      void ( *derivs ) ( int, double, double*,double [],double [] ),
                      void ( *data_ana ) ( int, double, double [], double *, int ) );
void rkdumb(complex99 vstart[], int nvar, double x1, double x2,
		int nstep, complex99 *noise, double *rlt,
		void (*derivs)(double,complex99*,int,complex99[],complex99[]),
		void (*data_ana)(double,int,complex99*,double*));

extern complex99 *dydx, *dym, *dyt, *yt;
int alloc_mem_for_propagator(int ny);
int free_mem_for_propagator();
int propagator(complex99 y[], int ny, double x1,
		double x2, int nx, double *rlt, complex99 *rdm,
		void (*derivs)(int, double,complex99*,complex99 [],complex99 []),
		void (*data_ana)(int, double, complex99 [], double *)) ;

int propagator_purtub(complex99 y[], int ny, double x1,
		double x2,int nx, complex99 *rdm,
		void (*derivs)(int, double,complex99*,complex99 [],complex99 []),
		complex99 *rho0,complex99 *rho1);

void propagator_derive(int ny, double x, complex99 *rdm,
		complex99* y, complex99* dydx);
void propagator_data_ana(int ny, double x, complex99 *y,double *rlt);

int accumulate_rlt(double *rlt_tmp,
		double *stat,int *good_orbit, int nt);
int accumulate_rlt_mult(double *rlt_tmp,
			 double *stat,int *good_orbit, int nt, int obvs);
//propagator for real
int propagator_real(double y[], int ny, double x1,
		    double x2, int nx, double *rlt, double *rdm, int nobv,
		    void (*derivs)(int, double, double*,double [],double []),
		    void (*data_ana)(int, double, double [], double *, int)) ;
int propagator_real_with_integration_kernal(double y[], int ny, double x1,
					    double x2, int nx, double *rlt, double *rdm, int nobv,
					    void (*derivs)(int, double, double*,double [],double []),
					    void (*data_ana)(int, double, double [], double *, int),
					    void (*dump_traj)(int,int, double*,double *));
int euler_real(double y[], int ny, double x1,
	       double x2, int nx, double *rlt, double *rdm, int nobv,
	       void (*derivs)(int, double, double*,double [],double []),
	       void (*data_ana)(int, double, double [], double *, int));
int alloc_mem_for_propagator_real(int ny);
int free_mem_for_propagator_real();
void data_analyse_example_real(int nvar,double x, double *y, double* rlt,int nobvs);
void propagator_derive_example_real(int ny, double x, double *rdm, double *y, double *dydx);
extern double *dydx_r, *dym_r, *dyt_r, *yt_r;
int propagator_real_save_intermidiate_data(double *y, FILE *fp,int ny);

#endif
