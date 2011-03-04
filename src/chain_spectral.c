#include <sys/time.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <fftw3.h>
#include "rnruni.h"
#include "integration.h"
#include "correlation.h"
#include "mytime.h"
#include "chain.h"

//#include <acml_mv.h>
// extern double T_l=0.010;//2.0;//Kelvin;//The temperature of left bath.
// extern double T_r=0.0;//1.0;//Kelvin;//The temperature of right bath.
// extern double tau_l=40.0;//fs//The
// extern double tau_r=40.0;//fs
// extern double epselon_l=5.0e-2;//fs^-1
// extern double epselon_r=5.0e-2;//fs^-1
// extern double mass=12.0;//12.0;//amu
// extern double k_b=8.314472477042e-3;
//
// extern int num_of_atoms=6;
// extern double x_eq=1.54e+2;//pm
// extern double D=3.84e+02/(6.0*6.0);//0.5*0.2700;//3.678e+02/(0.01*0.01);//unit of energy
// extern double a=1.875e-2*(6.0);//1.0;//1.875e-2*(0.01);//pm^-1
// //double k_force_constant=2*D*a^2=0.2700;//only for harmonic case
// extern double *memory_of_p1;
// extern double *memory_of_pN;
// extern int N_memory_step_left;
// extern int N_memory_step_right;
extern double *integration_kernal_l;
extern double *integration_kernal_r;
//extern double *integration_kernal;
extern double integration_kernal_time_step_length;
extern double *momentum_r;
extern double *momentum_l;
//extern double *momentum;
extern int momentum_for_current_time_r;
extern int momentum_for_current_time_l;
extern int truncation_time_step_number_l;
extern int truncation_time_step_number_r;
inline double propagator_auxiliary_integration(double * integration_kernal,
    int truncation_time_step_number, double *momentum, int momentum_for_current_time);

int chain_spectral(double *dt_of_rlt, int *nt_of_rlt,
                   double** rlt_mpi, int *good_orb, int *n_obvs, char* output_file_name,
                   int argc, char* argv[])
{

  int nt  = 1024*128;
  double dt = 0.01;
  double tf = nt*dt;


  {//output parameters to datafile
    myprint(output_file_name,"num_of_atoms","%d",&num_of_atoms);
    myprint(output_file_name,"D","%f",&D);
    myprint(output_file_name,"a","%f",&a);
    myprint(output_file_name,"Temperature left","%f",&T_l);
    myprint(output_file_name,"Temperature right","%f",&T_r);
    myprint(output_file_name,"memory left","%f",&tau_l);
    myprint(output_file_name,"memory right","%f",&tau_r);
    myprint(output_file_name,"epselon_l","%f",&epselon_l);
    myprint(output_file_name,"epselon_r","%f",&epselon_r);
    myprint(output_file_name,"mass","%f",&mass);
    myprint(output_file_name,"nt","%d",&nt);
    myprint(output_file_name,"dt","%f",&dt);
    myprint(output_file_name,"tf","%f",&tf);
    myprint(output_file_name,"function","%s","chain_spectral");
    myprint(output_file_name,"spectral","%s","\\eta \\omega 1/(1+\\tau^2*\\omega^2)");
    myprint(output_file_name,"comment","%s","2 baths with different spectral, anharmonic chain should be better because of frequency mix feature of anharmonic system");
  }
  //calculate \gamma(t) from \hat{gamma}(\omega);
  double para_gamma[2];
  para_gamma[0]=epselon_l;para_gamma[1]=tau_l;
  draw_spectral("spectral_left.dat",para_gamma,spectral_UO);
  double *alpha_l=gamma_t_from_spectral_truncate(2*nt, 0.5*dt, tau_l, 2, para_gamma, spectral_UO);//\gamma(t) for left bath
  para_gamma[0]=epselon_r;para_gamma[1]=tau_r;
  draw_spectral("spectral_right.dat",para_gamma,spectral_UO);
  double *alpha_r = gamma_t_from_spectral_truncate(2*nt, 0.5*dt, tau_r, 2, para_gamma, spectral_UO);//\gamma(t) for right bath
  {//output gamma_l(r)
    FILE *fp=fopen("gamma_l.dat","w");
    int i;
    for(i=0;i<2*nt;i++)
    {
      fprintf(fp,"%f  %e\n",i*dt*0.5,alpha_l[i]);
    }
    fclose(fp);
    fp=fopen("gamma_r.dat","w");
    for(i=0;i<2*nt;i++)
    {
      fprintf(fp,"%f  %e\n",i*dt*0.5,alpha_r[i]);
    }
    fclose(fp);
  }
  //open memory for integration kernal in Langevin equation
  truncation_time_step_number_l=truncation_time_step_number_estimate(0.5*dt, 2*nt, alpha_l,tau_l);
  truncation_time_step_number_r=truncation_time_step_number_estimate(0.5*dt, 2*nt, alpha_r,tau_r);
//   truncation_time_step_number_l=(int)(truncation_time_step_number_l/2);//alpha_l(r) use time step 0.5*dt. but the momentum_l(r) will use time step dt, so only half memory is needed.
//   truncation_time_step_number_r=(int)(truncation_time_step_number_r/2);//alpha_l(r) use time step 0.5*dt. but the momentum_l(r) will use time step dt, so only half memory is needed.
  {//output number of truncation time steps for integration kernal in Langevin equation
    myprint(output_file_name,"truncation_time_step_number_l","%d",&truncation_time_step_number_l);
    myprint(output_file_name,"truncation_time_step_number_r","%d",&truncation_time_step_number_r);
  }
  integration_kernal_l=(double*)malloc(sizeof(double)*(truncation_time_step_number_l+truncation_time_step_number_r));
  assert(integration_kernal_l);
  integration_kernal_r=integration_kernal_l+truncation_time_step_number_l;
  {//set values for integration_kernal_l(r)
    int i;
    for(i=0;i<truncation_time_step_number_l;i++)
//       integration_kernal_l[i]=alpha_l[2*i];
      integration_kernal_l[i]=alpha_l[i];

    for(i=0;i<truncation_time_step_number_r;i++)
//       integration_kernal_r[i]=alpha_r[2*i];
      integration_kernal_r[i]=alpha_r[i];

  }
  {//output integration_kernal_l(r)
    int i;
    FILE *fp=fopen("integration_kernal_l.dat","w");
    for(i=0;i<truncation_time_step_number_l;i++)
//       fprintf(fp,"%f  %e\n",i*dt,integration_kernal_l[i]);
      fprintf(fp,"%f  %e\n",i*dt/2.0,integration_kernal_l[i]);

    fclose(fp);
    fp=fopen("integration_kernal_r.dat","w");
    for(i=0;i<truncation_time_step_number_r;i++)
//       fprintf(fp,"%f  %e\n",i*dt,integration_kernal_r[i]);
      fprintf(fp,"%f  %e\n",i*dt/2.0,integration_kernal_r[i]);
    fclose(fp);
  }
//   integration_kernal_time_step_length=dt;
  integration_kernal_time_step_length=dt/2.0;
  momentum_l=(double*)malloc(sizeof(double)*(truncation_time_step_number_l+truncation_time_step_number_r));
  assert(momentum_l);
  bzero(momentum_l,sizeof(double)*(truncation_time_step_number_l+truncation_time_step_number_r));
  momentum_r=momentum_l+truncation_time_step_number_l;
  momentum_for_current_time_l=0;
  momentum_for_current_time_r=0;


  double para_t_m_k[3];
  para_t_m_k[0]=T_l;para_t_m_k[1]=mass;para_t_m_k[2]=k_b;
  alpha_t_from_gamma_t(2*nt, 0.5*dt, 3, para_t_m_k, alpha_l);
  para_t_m_k[0]=T_r;para_t_m_k[1]=mass;para_t_m_k[2]=k_b;
  alpha_t_from_gamma_t(2*nt, 0.5*dt, 3, para_t_m_k, alpha_r);
  {//output \alpha_l(r)=k_b T/m \gamma(t)
    FILE *fp=fopen("alpha_l.dat","w");
    int i;
    for(i=0;i<2*nt;i++)
    {
      fprintf(fp,"%f  %e\n",i*dt*0.5,alpha_l[i]);
    }
    fclose(fp);
    fp=fopen("alpha_r.dat","w");
    for(i=0;i<2*nt;i++)
    {
      fprintf(fp,"%f  %e\n",i*dt*0.5,alpha_r[i]);
    }
    fclose(fp);
  }
  double *lamda_l = lamda_sqrt_t(alpha_l,2*nt,0.5*dt);
  double *lamda_r = lamda_sqrt_t(alpha_r,2*nt,0.5*dt);
  assert(lamda_l&&lamda_r);
  {//output lamda_l(r)=Fourier[alhpa_l(r)]
    FILE *fp=fopen("lambda_l.dat","w");
    int i;
    double dw=2.0*3.1415926/(nt*dt);
    for(i=0;i<2*nt;i++)
    {
      fprintf(fp,"%f  %e\n",i*dw,lamda_l[i]);
    }
    fclose(fp);
    fp=fopen("lambda_r.dat","w");
    for(i=0;i<2*nt;i++)
    {
      fprintf(fp,"%f  %e\n",i*dw,lamda_r[i]);
    }
    fclose(fp);
  }

  free(alpha_l);
  free(alpha_r);
  corr_rand_gener_prep(2*nt);//open cache for FFTW3 operation680
  corr_rand_gener_2bath_malloc(2*nt);//open memory to store 4 random vector with correlation alpha_l and 4 rand vectors with correlation alpha_r

  //handle the rlt data
  int nobvs=1;//how many observables will be calculated in each step. used in (void*)data_ana(...)
  double *rlt = (double*)malloc(sizeof(double)*2*nt*nobvs);
  double *rlt_tmp = (double*)malloc(sizeof(double)*nt*nobvs);
  assert(rlt&&rlt_tmp);
  bzero(rlt,sizeof(double)*2*nt*nobvs);
  bzero(rlt_tmp,sizeof(double)*nt*nobvs);
  //handle propagator
  int num_of_terms=(num_of_atoms+2)*2;
  double *rho=(double*)malloc(sizeof(double)*num_of_terms);
  alloc_mem_for_propagator_real(num_of_terms);
  //handle the loop
  int i, N=4*2600,good_orbit=0;
  run_time time_loop_begin,time_loop_end;
  mytime(&time_loop_begin);
  {
    myprint(output_file_name,"num of orbs","%d",&N);
    myprint(output_file_name,"loop begin time","time",&(time_loop_begin));
  }
  for(i=0;i<N;i++)
  {
    {
      int one_percent=(int)(N/100);
      if (one_percent==0) one_percent=1;
      if(i%one_percent==0)
        printf("\n%d percent finished\n", (int)(100*i/N));
    }
    double *rand=0;
    corr_rand_gener_2bath(2*nt, &rand, lamda_r, lamda_l, i%4);
    {//output random noise
//             int i;
//             for(i=0;i<nt;i++)
//             {
// 	      rlt_tmp[i]=rand[4*i];
//             }
    }
    { //    test corr
      //       int j;
      //       int Nn=1024;//1024*56;
      //       int TT=nt-2*1024;
      //       for(j=0;j<TT;j++)
      // 	rlt_tmp[j]=rand[1]*rand[1+4*j];//4*j comes from 2 (left+right) multiplies by 2(dt/2 instead of dt)
    }

    initial_chain(rho,num_of_atoms);
    propagator_real_with_integration_kernal(rho, num_of_terms, 0.0, tf,
                                            nt, rlt_tmp,rand,nobvs,
//                                             propagator_derive_chain_white_harm,
//                                             propagator_derive_chain_UO_harm,
//                                             propagator_derive_chain_UO_mors,
                                            propagator_derive_chain_general_spectral,
//                                             data_analyse_chain_heat_flux_fast,
//                                             data_analyse_chain_temperature_profile
//                                             data_analyse_chain_heat_flux_harm_mult
//                                             data_analyse_chain_heat_flux_harm
                                            data_analyse_chain_heat_flux_mors,
//                                             data_analyse_chain_position_mult
//                                             data_analyse_chain_heat_flux_fast_mult
//                                             data_analyse_chain_total_energy_mors
                                            dump_momentum_of_1st_and_Nth_atoms
// 					    dump_random_force_of_1st_and_Nth_atoms
                                           );

    accumulate_rlt_mult(rlt_tmp,rlt,&good_orbit, nt,nobvs);

  }

  mytime(&time_loop_end);
  {
    myprint(output_file_name,"loop end time","time",&(time_loop_end));
  }
  //	free_mem_for_propagator_real();

  corr_rand_gener_free();
  free(lamda_r);
  free(lamda_l);
  free(momentum_l);
  free(integration_kernal_l);

  corr_rand_gener_2bath_free();
  free(rlt_tmp);
  // 	free(memory_of_p1);
  // 	free(memory_of_pN);
  //	free(equi);
  // put result back to main()
  *good_orb=good_orbit;
  *rlt_mpi = rlt;
  *nt_of_rlt = nt;
  *dt_of_rlt = dt;
  *n_obvs=nobvs;

  //output
  return 0;
}

double *integration_kernal_l;
double *integration_kernal_r;
// double *integration_kernal;
double integration_kernal_time_step_length;
double *momentum_r;
double *momentum_l;
// double *momentum;
int momentum_for_current_time_r;
int momentum_for_current_time_l;
int truncation_time_step_number_l;
int truncation_time_step_number_r;



/**
 * This function decide how many times needed to be stored to carry on the integration part of Langevin equation.
 * @param dt 
 * @param (* spect)( int , double *, double ) 
 * @param time_scale 
 * @return 
 */
int truncation_time_step_number_estimate(double dt, int nt, double *alpha,double time_scale)
{
  int time_scale_step=(int)(time_scale/dt);
  double alpha_at_time_scale;
  {
    if(fabs(alpha[0])/2>fabs(alpha[time_scale_step]))
      alpha_at_time_scale=fabs(alpha[0])/2;
    else
      alpha_at_time_scale=fabs(alpha[time_scale_step]);
  }
  double rate=0.0005;
  int i;
  for(i=nt-1; i>0; i--)
  {
    if(fabs(alpha[i]/alpha_at_time_scale)>rate)
      break;
  }
  return i;
}

/**
 * The function calculate the generalized friction term in langevin equation (the integration). The trajectory before current time t is stored in 
 * a ring-like data-structure, data was write into the ring sequentially. The underlying data structure is linear, so I need to manage the data with
 * several pointers.
 * for the first atom
 * y_1=\int_0^{t}{dt'\gamma(t-t')p(t')}=\int_0^{t}{d\tau \gamma(\tau) p(t-\tau)}; or in the descrete form
 * y_1=\delta\sum_{n=0}^{N-1}{\gamma_n*p_{N-1-n}}, where \delta*(N-1)=t is the truncation time. Notice that this is a naive way the carry on the integration (rectangle summation).
 * \gamma_n=\gamma(n*\delta); p_{N-n}=p(\delta*(N-n))
 * the data-structure of integration_kernal and momentum:
 * integration_kernal=[\gamma_0, \gamma_1, \gamma_2,...\gamma_{N-1}], where N is the truncation number
 * momentum=[p_{N-1-k},p_{N-1-k+1},p_{N-1-k+2},...,p_{N-1-k+k=N-1},p_{0},p_{1},p_{2},...p_{N-2-k}], where k is momentum_for_current_time.
 * see picture ring_data_structure_for_friction_kernal.jpg for a schematic discription
 * @param integration_kernal : \gamma(t) in friction term in the Langevin equation; see propagator_derive_chain_UO.jpg
 * @param truncation_time_step_number : because \gamma(t) decays quickly. there's no need to integrate from 0 to t, instead, integration inteval is 
 * chosen to be [t0,t]. t-t0 is truncation length, which is descreted into N steps, and N is truncation_time_step_number.
 * @param momentum :the pointer pointed to the head of the underlying data-structure for p_1(t) and p_N(t)
 * @param momentum_for_current_time : the actual position for p_1(t) and p_N(t)
 * @return the value of the integration, i.e. y_1 and y_N in propagator_derive_chain_UO.jpg
 */
inline double propagator_auxiliary_integration(double * integration_kernal,
    int truncation_time_step_number, double *momentum, int momentum_for_current_time)
{
  double rlt=0.0;

#define k momentum_for_current_time
#define N truncation_time_step_number
#define gamma integration_kernal
#define p momentum
  //all the macro define above is for consistance with ring_data_structure_for_friction_kernal.jpg
  double *p_curr=p+k; //p_curr=p': see pic:ring_data_structure_for_friction_kernal.jpg
  int i;
  for (i=0;i<k+1;i++)//\gamma_0 p_k+\gamma_1 p_{k-1}+... +\gamma_k p_0
  {
    rlt+=(gamma[i]*p[k-i]);
  }
  for (i=k+1;i<N;i++)
  {
    rlt+=(gamma[i]*p_curr[N-i]);
  }
#undef k
#undef N
#undef gamma
#undef p
  return rlt*integration_kernal_time_step_length;
}

/**
 * all the same as inline double propagator_auxiliary_integration(double * integration_kernal, int truncation_time_step_number, double *momentum, int momentum_for_current_time),
 * except Simpson rule is used instead of naive rectangle summation.
 * @param integration_kernal 
 * @param truncation_time_step_number 
 * @param momentum 
 * @param momentum_for_current_time 
 * @return 
 */
inline double propagator_auxiliary_integration_simpson(double * integration_kernal,
					       int truncation_time_step_number, double *momentum, int momentum_for_current_time)
{
  double rlt=0.0;
  if(truncation_time_step_number<9)
  {
    printf("propagator_auxiliary_integration_simpson:truncation_time_step_number should be larger than 9\n");
    exit(0);
  }
#define k momentum_for_current_time
#define N truncation_time_step_number
#define gamma integration_kernal
#define p momentum
  //all the macro define above is for consistance with ring_data_structure_for_friction_kernal.jpg
/*  double *p_curr=p+k; //p_curr=p': see pic:ring_data_structure_for_friction_kernal.jpg*/
  int i;
  double coef[5]={17.0,59.0,43.0,49.0,48.0};//the coefs are required by Simpson's rule, see http://en.wikipedia.org/wiki/Simpson's_rule
  for (i=0;i<N;i++)
  {
    int index;
    double c;
    if(i<=k) index=k-i;
    else index= k+N-i;
    if(i<4) c=coef[i];
    else if (i>N-5) c=coef[N-1-i];
    else c=coef[4];
    rlt+=c*(gamma[i]*p[index]);
  }
  rlt=(rlt*integration_kernal_time_step_length/48.0);

#undef k
#undef N
#undef gamma
#undef p
  return rlt;

}

/**
 * For a general form of langevin equation, there's a integration to be calculated for every time step. So I need to record the trajectory before current time step.
 * @param ny 
 * @param x 
 * @param rdm 
 * @param y 
 * @param dydx 
 */
void propagator_derive_chain_general_spectral(int ny, double x, double *rdm,
    double* y, double* dydx)
{
  /*
    y[0]---------------x_0 //not used in this function
    y[1]---------------y_1 //the integration term in Langevin equation. check propagator_derive_chain_UO.jpg
    y[2]---------------x_1 
    y[3]---------------p_1
    y[2*i]-------------x_{i}
    y[2*i+1]-----------p_{i}
    y[2*N]-------------x_N
    y[2*N+1]-----------p_N
    y[2*(N+1)]---------x_{N+1} // not used in this function
  y[2*(N+1)+1]-------y_N //the integration term in Langevin equation. check propagator_derive_chain_UO.jpg
  */
  int N=ny/2-2;
  double r1;
  bzero(dydx,sizeof(double)*ny);
  int i;
#define fastexp exp
  for(i=1;i<N;i++)
  {
    dydx[2*i]=y[2*i+1]/mass;//    (d/dt)x_i=p_i/mq
    r1=y[2*i+2]-y[2*i]-x_eq;
    r1=fastexp(-a*r1);
    dydx[2*i+1]+=(-r1*(r1-1)*2.0*a*D);
    dydx[2*(i+1)+1]-=(-r1*(r1-1)*2.0*a*D);
  }
  dydx[2*N]=y[2*N+1]/mass;

  {//this is for check the dumped trajectory
//       char name[64]={};
//       FILE *fp=fopen("dum_traj.dat","w");
//       int i;
//       for(i=0;i<truncation_time_step_number_l+truncation_time_step_number_r;i++)
//       {
// 	fprintf(fp,"%e\n",momentum_l[i]);
//       }
//       fclose(fp);
//       fp=fopen("dum_integration_kernal.dat","w");
//       for(i=0;i<truncation_time_step_number_l+truncation_time_step_number_r;i++)
//       {
// 	fprintf(fp,"%e\n",integration_kernal_l[i]);
//       }
//       fclose(fp);
  }
//   dydx[1]=propagator_auxiliary_integration(integration_kernal_l, truncation_time_step_number_l, momentum_l, momentum_for_current_time_l);
//   dydx[2*N+3]=propagator_auxiliary_integration(integration_kernal_r, truncation_time_step_number_r, momentum_r, momentum_for_current_time_r);//see propagator_derive_chain_UO.jpg
  dydx[1]=propagator_auxiliary_integration_simpson(integration_kernal_l, truncation_time_step_number_l, momentum_l, momentum_for_current_time_l);
  dydx[2*N+3]=propagator_auxiliary_integration_simpson(integration_kernal_r, truncation_time_step_number_r, momentum_r, momentum_for_current_time_r);//see propagator_derive_chain_UO.jpg

  r1=y[2]-y[0]-x_eq;
  r1=fastexp(-a*r1);
  dydx[3]-=(-r1*(r1-1)*2.0*a*D);
  r1=y[2*(N+1)]-y[2*N]-x_eq;
  r1=fastexp(-a*r1);
  dydx[2*N+1]+=(-r1*(r1-1)*2.0*a*D);
  dydx[3]=dydx[3]-y[1]+mass*rdm[1];
  dydx[2*N+1]=dydx[2*N+1]-y[2*N+3]+mass*rdm[0];

#undef fastexp
  //test x_0 periodic and x_{N+1} damping:
  //  dydx[0]=0.1*sin(x);
  //  dydx[2*N+1]-=y[2*N+1];
}

/**
 * Basically the same as 4th order Runge-Kutta method, except it stores the trajectories of the 1st and Nth atom for calculating the friction kernal in Langevin equation. To inhance the speed, a
 * ring-like data structure is used to store the trajectories. In this function, a pointer need to be managed (k in pic:ring_data_structure_for_friction_kernal.jpg) to keep the data structure.
 * @param y[] 
 * @param ny 
 * @param x1 
 * @param x2 
 * @param nx 
 * @param rlt 
 * @param rdm 
 * @param nobv 
 * @param (* derivs)( int , double , double *, double [], double [] ) 
 * @param (* data_ana)( int , double , double [], double *, int ) 
 * @param (* dump_traj)(int ny, int k, double *y, double* data) k is the current time step, dum_traj uses this number to determine the dumping site of the data, 
 * data[] points to the data array. The
 * realization is basically the same as data_ana
 * @return 
 */
int propagator_real_with_integration_kernal(double y[], int ny, double x1,
    double x2, int nx, double *rlt, double *rdm, int nobv,
    void (*derivs)(int, double, double*,double [],double []),
    void (*data_ana)(int, double, double [], double *, int),
    void (*dump_traj)(int,int, double*,double *))
{


  double h=(x2-x1)/nx;
  double hh=h*0.5, h6=h/6.0, x, xh;
  int xn, i;
#define SAVE_INTERMIDIATE_DATA 0
#if SAVE_INTERMIDIATE_DATA==1
  FILE *fp = fopen("tmp_data_for_propagator_real.dat","a");
#endif

  for(xn=0, x=x1; xn<nx; xn++)
  {
    xh=x+hh;
    (*derivs)(ny, x, rdm, y, dydx_r);
    //    (*data_ana)(ny,x,y,rlt,nobv);//to calculate the observables here is to avoid to calculate the force twice(one in propagation, one in data_ana).
    for(i=0;i<ny;i++) yt_r[i]=y[i]+hh*dydx_r[i];
    (*derivs)(ny, xh, rdm+2, yt_r, dyt_r);
    for(i=0;i<ny;i++) yt_r[i]=y[i]+hh*dyt_r[i];
    (*derivs)(ny, xh, rdm+2, yt_r, dym_r);
    for(i=0;i<ny;i++)
    {
      yt_r[i]=y[i]+h*dym_r[i];
      dym_r[i]+=dyt_r[i];
    }
    (*derivs)(ny, x+h, rdm+4, yt_r, dyt_r);

    //move the the new point
//     for(i=0;i<ny;i++)
//       y[i]+=h6*(dydx_r[i]+dyt_r[i]+2.0*dym_r[i]);
    //to dump the data in the middle of the step, the equation above is divided into 2 parts
    for(i=0;i<ny;i++)
      y[i]+=h6*(dydx_r[i]+dym_r[i]);
    (*dump_traj)(ny,2*xn,y,momentum_l);//dump the traj point at the middle of the step.
    {//this is for check the dumped trajectory
//       FILE *fp=fopen("dum_traj.dat","w");
//       int i;
//       for(i=0;i<truncation_time_step_number_l+truncation_time_step_number_r;i++)
//       {
// 	fprintf(fp,"%e\n",momentum_l[i]);
//       }
//       fclose(fp);
    }
    for(i=0;i<ny;i++)
      y[i]+=h6*(dyt_r[i]+dym_r[i]);
    (*dump_traj)(ny,2*xn+1,y,momentum_l);//momentum_r=momentum_l+truncation_time_step_number_l, so there's no need to transfer this parameter.
    {//this is for check the dumped trajectory
//       FILE *fp=fopen("dum_traj.dat","w");
//       int i;
//       for(i=0;i<truncation_time_step_number_l+truncation_time_step_number_r;i++)
//       {
// 	fprintf(fp,"%e\n",momentum_l[i]);
//       }
//       fclose(fp);
    }
    x+=h;
    rdm+=4;
    //rk4 beyond
    (*data_ana)(ny,x,y,rlt,nobv);
    rlt+=nobv;
//    (*dump_traj)(ny,xn,y,momentum_l);//momentum_r=momentum_l+truncation_time_step_number_l, so there's no need to transfer this parameter.
//     (*dump_traj)(ny,xn,rdm,momentum_l);//dump random forces


#if SAVE_INTERMIDIATE_DATA==1
    propagator_real_save_intermidiate_data(y,fp,ny);
#endif

  }
#if SAVE_INTERMIDIATE_DATA==1
  fclose(fp);
#endif
#undef SAVE_INTERMIDIDATE_DATA
  return 0;

}

void dump_position_of_1st_and_Nth_atoms(int ny, int xn, double *y, double *data)
{
  double *atom_l=data;
  double *atom_r=data+truncation_time_step_number_l;
  int k_l=xn%truncation_time_step_number_l;
  int k_r=xn%truncation_time_step_number_r;
  momentum_for_current_time_r=k_r;
  momentum_for_current_time_l=k_l;
  //position of atom_1 is y[2]
  //position of atom_N is y[2*N]
  int N=ny/2-2;
  atom_l[k_l]=y[2];
  atom_r[k_r]=y[2*N];


}
void dump_momentum_of_1st_and_Nth_atoms(int ny, int xn, double *y, double *data)
{
  double *atom_l=data;
  double *atom_r=data+truncation_time_step_number_l;
  int k_l=xn%truncation_time_step_number_l;
  int k_r=xn%truncation_time_step_number_r;
  momentum_for_current_time_r=k_r;
  momentum_for_current_time_l=k_l;
  //momentum of atom_1 is y[3]
  //momentum of atom_N is y[2*N+1]
  int N=ny/2-2;
  atom_l[k_l]=y[3];
  atom_r[k_r]=y[2*N+1];


}
void dump_random_force_of_1st_and_Nth_atoms(int ny, int xn, double *y, double *data)
{
  double *atom_l=data;
  double *atom_r=data+truncation_time_step_number_l;
  int k_l=xn%truncation_time_step_number_l;
  int k_r=xn%truncation_time_step_number_r;
  momentum_for_current_time_r=k_r;
  momentum_for_current_time_l=k_l;
  atom_l[k_l]=y[0];
  atom_r[k_r]=y[1];


}
/**
 *This function calculate response function \gamma(t) from \hat{\gamma}(\omega), using N discrete points in \omega axe. The problem is that we need to increase N to check convergence, which will be done in function alpha_r_t_spectral(...). See fft_spectral_response_1(2).jpg
 * @param N Number of points used to take FFT tranformation
 * @param dt time steps
 * @param n_para number of parameters used to calculate \hat{\gamma}(\omega)
 * @param para parameters used to calculate \hat{\gamma}(\omega)
 * @param (* spect)( int , double *, double ) 
 * @return a N-array containing \gamma(ti), where ti=i*dt
 */
double *gamma_t_from_spectral_test(int N, double dt, int n_para, double *para, double (*spect)(int, double*, double))
{
  fftw_complex *rf=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);//inverse fourier transform (according to fftw3 document) \gamma(t)
  //http://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html#The-1d-Discrete-Fourier-Transform-_0028DFT_0029

  fftw_complex *sd=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);//spectral density \hat{\gamma}(\omega0

  double dw=2.0*3.1415926/(N*dt);

  int i;
  //  FILE *fp=fopen("spectral.dat","w");
  for(i=0;i<N;i++)
  {
    double wi=i*dw;
    sd[i]=(*spect)(n_para,para,wi);
    //    fprintf(fp,"%f   %e\n",wi,sd[i]);
  }
  fftw_plan fftpln = fftw_plan_dft_1d (N, sd, rf, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(fftpln);

  double *rlt=(double*)malloc(sizeof(double)*N);
  double fct=2.0/(N*dt);
  for(i=0;i<N;i++)
  {
    rlt[i]=fct*creal(rf[i]);
  }
  fftw_free(rf);
  fftw_free(sd);
  fftw_destroy_plan(fftpln);
  //  fclose(fp);
  return rlt;

}
/**
 *this function calculate the response function \gamma(t) according to the spectral function \hat｛\gamma｝(\omega)=\eta\omega cutoff(\omega), where cutoff(\omega) is a cutoff function guarantee J(\omega) goes to zero when \omega goes to infinity. The relation between \gamma(t) and \hat{\gamma}(\omega) (Eq.4 of THE JOURNAL OF CHEMICAL PHYSICS 128, 224710 2008) \hat{\gamma}(\omega)=\int{\exp(-i\omega t)\gamma(t)dt}, so we can take FFT directly to \hat{\gamma}(\omega) to calculate \gamma(t). gamma_t_from_spectral_test directly carries on FFT. This function uses gamma_t_from_spectral_test multiple times with increasing N to check and garantee convergence.
 * @param nt 
 * @param dt 
 * @param para 
 * @return 
 */
double *gamma_t_from_spectral(int nt, double dt, double tc, int n_para, double *para, double (*spect)(int, double*, double))
{
  int N=nt;
  double time_scale=tc;
  int convergence=0;
  double *rlt=gamma_t_from_spectral_test(N, dt, n_para, para, spect);
  double error=0.0;
  while(convergence==0 && N<1024*1024*16)
  {
    N=N*2;
    double *rlt1=gamma_t_from_spectral_test(N, dt, n_para, para, spect);
    int i_time_scale=(int)(time_scale/dt);
    error=(rlt1[i_time_scale]-rlt[i_time_scale])/rlt[i_time_scale];
    //    printf("error at tc:%f\n",fabs(error));
    if (fabs(error)<0.001)
    {
      convergence=1;
    }
    //    printf("convergence:%d\n",convergence);
    free(rlt);
    rlt=rlt1;
  }
  if (convergence==0)
    printf("gamma_t_from_spectral: precision convergence not reached");
  //  printf("%d points in spectral function is used, truncation frequency is %e, precision at tc is %e\n",N,2.0*3.14159/dt,error);
  double *response_function=(double*)malloc(sizeof(double)*nt);
  int i;
  for(i=0;i<nt;i++)
  {
    response_function[i]=rlt[i];
  }
  free(rlt);
  return response_function;
}
/**
 * There are considerable numeric errors with FFT scheme when t is large enough. For example, when t=6000, the result of FFT tranform scheme is 10 times large as the exact result, although, both virtually zero. So I decide to truncate the result, that is, to set alpha(t)=0 when t>n_max*t_c, where tc is characteristic time scale and n_max is a integer.  This function will use  gamma_t_from_spectral(...) to calculate \gamma(t) up to t_max=n_nax*tc. \gamma(t) above t_max is set to zero.
 * @param nt number of time steps
 * @param dt length of time steps
 * @param tc tc is the charactoristic time scale of the response function. Above tmax=n_max*tc, \gamma(t) is set to zero. n_max is determined
 * recursively in the function.
 * @param n_para number of parameters 
 * @param para  parameters
 * @param (* spect)( int , double *, double ) 
 * @return 
 */
double *gamma_t_from_spectral_truncate(int nt, double dt, double tc, int n_para, double *para, double (*spect)(int, double*, double))
{
//  int n_max=2;
  int convergence=0;
  int nc=(int)(tc/dt);
  int N=nc;
  double *rlt=NULL;
  while(convergence==0 && N<1024*1024*16)
  {
    N=N*2;
    rlt=gamma_t_from_spectral(N, dt, tc, n_para, para, spect);
    double ratio=rlt[nc-1]/rlt[N-1];
    //    printf("length:%f  ratio:%f\n", N*dt, ratio);
    if(fabs(ratio)>1000)// if ratio is large enough, the truncation time is long enough.
    {
      convergence=1;
    }
    else
    {
      free(rlt);
    }
  }
  if (convergence==0)
    printf("gamma_t_from_spectral_truncate: time length convergence not reached\n");
  //  printf("truncation time is %f\n",N*dt);

  //fill zero to the rest of the array.
  double *rlt1=(double*)malloc(sizeof(double)*nt);
  int i;
  for(i=0;i<N&&i<nt;i++)
  {
    rlt1[i]=rlt[i];
  }
  free(rlt);
  for(i=N;i<nt;i++)
  {
    rlt1[i]=0.0;
  }
  return rlt1;
}

double *alpha_t_from_gamma_t(int nt, double dt, int n_para, double *para, double *gamma_t)
{
  //para[0]=temperature
  //para[1]=mass
  //para[2]=k_b
  double T=para[0];
  double mass=para[1];
  double k_b=para[2];
  double fct=k_b*T/mass;
  int i;
  for(i=0;i<nt;i++)
  {
    gamma_t[i]=gamma_t[i]*fct;
  }
  return gamma_t;

}
int  draw_spectral(char *name, double *para,double (*spect)(int, double*, double))
{
  FILE *fp=fopen(name,"w");
  int i;
  double dw=0.1;
  int N=1024;
  for(i=0;i<N;i++)
  {
    double w=i*dw;
    double now=spect(2,para,w);
    fprintf(fp,"%f  %e\n",w,now);
  }
  fclose(fp);
  return 0;
}

/**
 * This function defines the Omega case spectral density function \hat{\gamma}(\omega)=\eta*\omega*\exp*(-\omega/\omega_c)
 * @param n_para Number of parameters 
 * @param para pointer to a paramters
 * @param o \omega
 * @return the spectral density  \hat{\gamma}(\omega)
 */
double spectral_OM(int n_para, double *para, double o)
{
  //omega case spectral density function with exponential cutoff function
  double eta=para[0];
  double tc=para[1];//tc=1/\omega_c
  return  eta*o/(1+o*o*tc*tc);
  //  return  eta*o*exp(-o*tc);
}

double spectral_UO(int n_para, double *para, double o)
{
  /*
    para[0]=characteristic spectral density: \epsilon_n in THE JOURNAL OF CHEMICAL PHYSICS 128, 224710 (2008);
    para[1]=characteristic time scale:\tau^{n}_{c};
  */
  double en=para[0];
  double tn=para[1];
  return 2.0*en/(1+o*o*tn*tn);
}

/**
 * This function defines spectral density function for test purples
 * @param n_para Number of parameters
 * @param para pointer to parameters
 * @param o \omega
 * @return the spectral density \hat{\gamma}(\omega)
 */
double spectral_test(int n_para, double *para, double o)
{
  //exponential spectral density
  //\hat{\gamma}(\omega)=\eta*\exp(-\omega/\omega_c)
  double eta=para[0];//eta
  double W=para[1];//W is omega_c
  return eta*exp(-o/W);
}
