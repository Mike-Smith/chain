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
double gamma_para_b_l=0.0;
double gamma_para_b_r=0.0;


/*
The parameters:
To make the results for chain_expcos roughly comparable to those of chain (with UO noise), I follow some corresponses.
In expcos noise case, \gamma(\omega)= 2*c*a*(w^2+a^2+b^2)/((a^2+w^2-2*w*b+b^2)*(a^2+w^2+2*w*b+b^2))
\gamma(t)=c*exp(-a*t)*cos(b*t)
In UO noise, \gamma(\omega)=2*e/(1+w^2*tau^2)
\gamma(t)=e/tau*exp(-t/tau)
1. integrations of \gamma(\omega) for UO case and expcos case are roughly the same. which means c=e/tau
2. b and a are determined by the position and hieght of the peak.
3. The width of \gamma(\omega) for UO case is roughly 1/tau=0.025-0.1fs^-1, So for expcos case, I choose the height of the peak so that the width of \gamma(\omega) to be roughly the same. The results are also displayed in find_a_b.m. It's lucky characteristic frequency of the bond is 0.15fs^-1, so the width of \gamma(\omega) is small enough to distinguish omega and 2omega.
*/

int chain_expcos(double *dt_of_rlt, int *nt_of_rlt,
                 double** rlt_mpi, int *good_orb, int *n_obvs, char* output_file_name,
                 int argc, char* argv[])
{

  if(argc!=8)
  {
    printf("wrong input parameters\n");
    exit(0);
  }
  else
  {
    tau_l=atof(argv[2]);//1/a
    tau_r=atof(argv[3]);
    epselon_l=atof(argv[4]);//c
    epselon_r=atof(argv[5]);
    gamma_para_b_l=atof(argv[6]);//b
    gamma_para_b_r=atof(argv[7]);

  }

  int nt  = 1024*1024;

  double dt = 0.1;//123
  double tf = nt*dt;
  /*  gamma_para_b_l=1.0;
    gamma_para_b_r=1.0;
    tau_l=20.0;//fs//The
    tau_r=2.0;//fs
    epselon_l=5e-2;//5.0e-2*2.0/3.1415926;//fs^-1
    epselon_r=5e-1;//5.0e-2*2.0/3.1415926;//fs^-1*/
  T_l=300.0;
  T_r=0.0;
  num_of_atoms=10;

  D=3.84e+02/(6.0*6.0);//0.5*0.2700;//3.678e+02/(0.01*0.01);//unit of energy
  a=1.875e-2*(6.0);//1.0;//1.875e-2*(0.01);//pm^-1

  {

      myprint(output_file_name,"num_of_atoms","%d",&num_of_atoms);
      myprint(output_file_name,"D","%f",&D);
      myprint(output_file_name,"a","%f",&a);
      myprint(output_file_name,"Temperature left","%f",&T_l);
      myprint(output_file_name,"Temperature right","%f",&T_r);
      myprint(output_file_name,"memory left","%f",&tau_l);
      myprint(output_file_name,"memory right","%f",&tau_r);
      myprint(output_file_name,"epselon_l","%f",&epselon_l);
      myprint(output_file_name,"epselon_r","%f",&epselon_r);
      myprint(output_file_name,"gamma_para_b_l","%f",&gamma_para_b_l);
      myprint(output_file_name,"gamma_para_b_r","%f",&gamma_para_b_r);
      myprint(output_file_name,"mass","%f",&mass);
      myprint(output_file_name,"nt","%d",&nt);
      myprint(output_file_name,"dt","%f",&dt);
      myprint(output_file_name,"tf","%f",&tf);
      myprint(output_file_name,"function","%s","chain_expcos.c");
      myprint(output_file_name,"explain","%s","the same as chain.c except it uses \\gamma(t)=c*exp(-at)cos(bt) as response function");

  }

  //return the memory function of the form \gamma(t)=c*exp(-a*t)*cos(b*t)
  //where a=1/\tau=1/para[4], b=gamma_para_b=para[5],c=\epselon=para[3],
  double para[6];
  para[0]=T_l;para[1]=mass;para[2]=k_b;para[3]=epselon_l;para[4]=tau_l;para[5]=gamma_para_b_l;
  double *alpha_l = alpha_r_t_expcos(2*nt, 0.5*dt, para);//generate memory function
  para[0]=T_r;para[1]=mass;para[2]=k_b;para[3]=epselon_r;para[4]=tau_r;para[5]=gamma_para_b_r;
  double *alpha_r = alpha_r_t_expcos(2*nt, 0.5*dt, para);



  double *lamda_l = lamda_sqrt_t(alpha_l,2*nt,0.5*dt);
  double *lamda_r = lamda_sqrt_t(alpha_r,2*nt,0.5*dt);
  {//output lamda_l(r)=Fourier[alhpa_l(r)]
    // #if  __WRITE_DATA_TO_FILE__ == 2
    //0 nowhere, 1 stdout, 2 datafile
    if(__WRITE_DATA_TO_FILE__ == 2 )
    {
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
    // #endif
  }

  free(alpha_l);
  free(alpha_r);
  corr_rand_gener_prep(2*nt);//open cache for FFTW3 operation

  double *rand=0;
  corr_rand_gener_2bath_malloc(2*nt);//open memory to store 4 random vector with correlation alpha_l and 4 rand vectors with correlation alpha_r



  //handle the rlt data
  int nobvs=1;//num_of_atoms;//how many observables will be calculated in each step.
  double *rlt = (double*)malloc(sizeof(double)*2*nt*nobvs);
  double *rlt_tmp = (double*)malloc(sizeof(double)*nt*nobvs);
  assert(rlt&&rlt_tmp);
  bzero(rlt,sizeof(double)*2*nt*nobvs);
  bzero(rlt_tmp,sizeof(double)*nt*nobvs);
  //handle propagator
  int num_of_terms=(num_of_atoms+3)*2;
  double *rho=(double*)malloc(sizeof(double)*num_of_terms);
  alloc_mem_for_propagator_real(num_of_terms);
  //handle the loop

  int i,good_orbit=0, N=4*10;
  run_time time_loop_begin,time_loop_end;
  mytime(&time_loop_begin);
  {
    myprint(output_file_name,"num of orbs","%d",&N);
    myprint(output_file_name,"loop begin time","time",&(time_loop_begin));
  }
  for(i=0;i<N;i++)
  {
    {
      //       int one_percent=(int)(N/100);
      //       if (one_percent==0) one_percent=1;
      //       if(i%one_percent==0)
      //         printf("\n%d percent finished\n", (int)(100*i/N));
    }

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
    initial_chain_expcos(rho,num_of_atoms);

    propagator_real(rho, num_of_terms, 0.0, tf,
                    nt, rlt_tmp,rand,nobvs,
//                      propagator_derive_chain_expcos_mors,
                       propagator_derive_chain_expcos_harm,
//                    data_analyse_chain_heat_flux_mors_expcos
                    data_analyse_chain_heat_flux_harm_expcos
                   );

    accumulate_rlt_mult(rlt_tmp,rlt,&good_orbit, nt,nobvs);

  }

  mytime(&time_loop_end);
  {
    myprint(output_file_name,"loop end time","time",&(time_loop_end));
  }
  //	free_mem_for_propagator_real();
#if _NOISE==1
  corr_rand_gener_free();
  free(lamda_r);
  free(lamda_l);
#endif
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

double *alpha_r_t_expcos(int nt, double dt, double *para)
{//  para[0]=T; para[1]=mass; para[2]=k_b; para[3]=epselon; para[4]=tau; para[5]=gamma_para_b;see Eq. 2 in JCP 128, 224710 (2008)
  //return the memory function of the form \gamma(t)=c*exp(-a*t)*cos(b*t)
  //where a=1/\tau=1/para[4], b=gamma_para_b=para[5],c=\epselon=para[3],
  int i;
  double a=1.0/para[4],b=para[5],c=para[3],T=para[0],m=para[1],k=para[2];
  //  printf("%f\t%f\t%f\t%f\t%f\n",para[0],para[1],para[2],para[3],para[4]);
  FILE *fp=0;
  char fn[1024];
  if(__WRITE_DATA_TO_FILE__ == 2)
  {
    // #if  __WRITE_DATA_TO_FILE__ == 2
    //0 nowhere, 1 stdout, 2 datafile
    sprintf(fn,"alpha_r_t_expcos%1.1f%1.1f%1.1f%1.1f%1.1f%1.1f.dat",para[0],para[1],para[2],para[3],para[4],para[5]);
    fp=fopen(fn,"w");
  }
  // #endif
  double *rlt = malloc(sizeof(double)*nt);
  for(i=0;i<nt;i++)
  {//rlt=(k_b*T/m)c*exp(-a*t)*cos(b*t)
    double t=dt*i;
    rlt[i]=(k*T/m)*c*exp(-a*t)*cos(b*t);
    if(__WRITE_DATA_TO_FILE__ == 2)
    {
      // #if  __WRITE_DATA_TO_FILE__ == 2
      //0 nowhere, 1 stdout, 2 datafile
      fprintf(fp,"%f\t%f\n",t,rlt[i]);
      // #endif
    }
  }
  // #if  __WRITE_DATA_TO_FILE__ == 2
  if(__WRITE_DATA_TO_FILE__ == 2)
  {

    //0 nowhere, 1 stdout, 2 datafile
    fclose(fp);
  }
  // #endif
  // #if  __WRITE_DATA_TO_FILE__ == 2
  if(__WRITE_DATA_TO_FILE__ == 2)
  {
    //0 nowhere, 1 stdout, 2 datafile
    //output the spectral function correspond to this gamma(t)
    sprintf(fn,"gamma_o_expcos%1.1f%1.1f%1.1f%1.1f%1.1f%1.1f.dat",para[0],para[1],para[2],para[3],para[4],para[5]);
    fp=fopen(fn,"w");
    double dw=2*3.1415926/(nt*dt);
    for(i=0;i<nt;i++)
    {//
      double w=i*dw;
      double gw=2.0*c*a*(w*w+a*a+b*b)/((a*a+w*w+b*b-2.0*w*b)*(a*a+w*w+b*b+2.0*w*b));
      fprintf(fp,"%f\t%e\n",w,gw);
    }
    fclose(fp);
  }
  // #endif
  return rlt;
}
int initial_chain_expcos(double *y,int num_of_atoms)
{
  /*
  y[0]---------------x_0
  y[1]---------------y_1// y1(0)=0
  y[2]---------------x_1
  y[3]---------------p_1
  y[2*i]-------------x_{i}
  y[2*i+1]-----------p_{i}
  y[2*N]-------------x_N
  y[2*N+1]-----------p_N
  y[2*(N+1)]---------x_{N+1}
  y[2*(N+1)+1]-------y_N// yN(0)=0;
  y[2*(N+1)+2]-------y_b// yb(0)=ya(0)=0: see expcos.jpg for definition of y_a and y_b
  y[2*(N+1)+3]-------y_a
  N is num_of_atoms
  */


  y[0]=0;
  y[1]=0;
  y[2*num_of_atoms+2]=(num_of_atoms+1)*x_eq;
  y[2*num_of_atoms+3]=0;
  y[2*num_of_atoms+4]=0;
  y[2*num_of_atoms+5]=0;
  //all at equil position still
  int i;
  for(i=1;i<=num_of_atoms;i++)
  {
    y[2*i]=i*x_eq;
    y[2*i+1]=0;
  }
  return 0;

}
void propagator_derive_chain_expcos_harm(int ny, double x, double *rdm,
    double* y, double* dydx)
{//halmiltonian chain
  /*
    y[0]---------------x_0
    y[1]---------------y_1
    y[2]---------------x_1
    y[3]---------------p_1
    y[2*i]-------------x_{i}
    y[2*i+1]-----------p_{i}
    y[2*N]-------------x_N
    y[2*N+1]-----------p_N
    y[2*(N+1)]---------x_{N+1}
    y[2*(N+1)+1]-------y_N
    y[2*(N+1)+2]-------y_b// yb(0)=ya(0)=0: see expcos.jpg for definition of y_a and y_b
    y[2*(N+1)+3]-------y_a
  */
  int N=ny/2-3;
  int i;
  double r1;
  bzero(dydx,sizeof(double)*ny);
  double k_cons=D*2.0*a*a;
  for(i=1;i<N;i++)
  {
    dydx[2*i]=y[2*i+1]/mass;//    (d/dt)x_i=p_i/m
    r1=y[2*i+2]-y[2*i]-x_eq;//r1=x_{i+1}-x_i-x_eq
    dydx[2*i+1]+=r1*k_cons;//(d/dt)p_i+=k*(x_{i+1}-x_i), see eq.1 and eq.2@ hamiltonian.jpg
    dydx[2*(i+1)+1]-=r1*k_cons;
  }
  dydx[2*N]=y[2*N+1]/mass;
  r1=y[2]-y[0]-x_eq;
  dydx[3]-=r1*k_cons;
  r1=y[2*(N+1)]-y[2*N]-x_eq;
  dydx[2*N+1]+=r1*k_cons;
  dydx[3]=dydx[3]-y[1]+mass*rdm[1];
  dydx[2*N+1]=dydx[2*N+1]-y[2*N+3]+mass*rdm[0];
#define _A_l ((1.0)/(tau_l))
#define _A_r ((1.0)/(tau_r))
#define _B_l (gamma_para_b_l)
#define _B_r (gamma_para_b_r)
#define _C_l (epselon_l)
#define _C_r (epselon_r)
  dydx[1]=_C_l*y[3]-_A_l*y[1]-_B_l*y[2*N+5];
  dydx[2*N+3]=_C_r*y[2*N+1]-_A_r*y[2*N+3]-_B_r*y[2*N+4];//see propagator_derive_chain_UO.jpg
  dydx[2*N+5]=-_A_l*y[2*N+5]+_B_l*y[1];
  dydx[2*N+4]=-_A_r*y[2*N+4]+_B_r*y[2*N+3];
#undef _A_l
#undef _A_r
#undef _B_l
#undef _B_r
#undef _C_l
#undef _C_r
  //test x_0 periodic and x_{N+1} damping:
  //  dydx[0]=0.1*sin(x);
  //  dydx[2*N+1]-=y[2*N+1];
}
void propagator_derive_chain_expcos_mors(int ny, double x, double *rdm,
    double* y, double* dydx)
{
  /*
    y[0]---------------x_0
    y[1]---------------y_1
    y[2]---------------x_1
    y[3]---------------p_1
    y[2*i]-------------x_{i}
    y[2*i+1]-----------p_{i}
    y[2*N]-------------x_N
    y[2*N+1]-----------p_N
    y[2*(N+1)]---------x_{N+1}
    y[2*(N+1)+1]-------y_N
    y[2*(N+1)+2]-------y_b// yb(0)=ya(0)=0: see expcos.jpg for definition of y_a and y_b
    y[2*(N+1)+3]-------y_a
  */
  int N=ny/2-3;
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
  r1=y[2]-y[0]-x_eq;
  r1=fastexp(-a*r1);
  dydx[3]-=(-r1*(r1-1)*2.0*a*D);
  r1=y[2*(N+1)]-y[2*N]-x_eq;
  r1=fastexp(-a*r1);
  dydx[2*N+1]+=(-r1*(r1-1)*2.0*a*D);
  dydx[3]=dydx[3]-y[1]+mass*rdm[1];
  dydx[2*N+1]=dydx[2*N+1]-y[2*N+3]+mass*rdm[0];
#define _A_l ((1.0)/(tau_l))
#define _A_r ((1.0)/(tau_r))
#define _B_l (gamma_para_b_l)
#define _B_r (gamma_para_b_r)
#define _C_l (epselon_l)
#define _C_r (epselon_r)
  dydx[1]=_C_l*y[3]-_A_l*y[1]-_B_l*y[2*N+5];
  dydx[2*N+3]=_C_r*y[2*N+1]-_A_r*y[2*N+3]-_B_r*y[2*N+4];//see propagator_derive_chain_UO.jpg
  dydx[2*N+5]=-_A_l*y[2*N+5]+_B_l*y[1];
  dydx[2*N+4]=-_A_r*y[2*N+4]+_B_r*y[2*N+3];
#undef _A_l
#undef _A_r
#undef _B_l
#undef _B_r
#undef _C_l
#undef _C_r
#undef fastexp
  //test x_0 periodic and x_{N+1} damping:
  //  dydx[0]=0.1*sin(x);
  //  dydx[2*N+1]-=y[2*N+1];
}
void data_analyse_chain_heat_flux_mors_expcos(int nvar,double x,
    double *y, double* rlt, int nobv)
{//potential energy
  /*
  y[0]---------------x_0
  y[1]---------------y_1
  y[2]---------------x_1
  y[3]---------------p_1
  y[2*i]-------------x_{i}
  y[2*i+1]-----------p_{i}
  y[2*N]-------------x_N
  y[2*N+1]-----------p_N
  y[2*(N+1)]---------x_{N+1}
  y[2*(N+1)+1]-------y_N
  y[2*(N+1)+2]-------y_b// yb(0)=ya(0)=0: see expcos.jpg for definition of y_a and y_b
  y[2*(N+1)+3]-------y_a
  N is num_of_atoms
  */
  int N=nvar/2-3;
  int i;
  double J=0.0;
  double F=0.0;
  //  printf("%f ",x);
  for(i=1;i<N;i++)
  {
    double r1=y[2*i+2]-y[2*i]-x_eq;//r1=x_{i+1}-x_i-x_eq
    //     double r2=y[2*i]-y[2*i-2]-x_eq;//r2=x_{i}-x_{i-1}-x_eqD
    r1=exp(-a*r1);
    //     r2=exp(-a*r2);
    F=-2*a*D*r1*(r1-1);//see eq.1 and eq.2@ hamiltonian.jpg
    double Jk=F*(y[2*i+1]+y[2*(i+1)+1])/mass;//v=p/m
    J+=Jk;
    //    printf("%08f ",Jk);
  }
  //  printf("\n");
  *rlt=-J/(2*(N-1));
}
void data_analyse_chain_heat_flux_harm_expcos(int nvar,double x,
    double *y, double* rlt, int nobv)
{//potential energy
  /*
  y[0]---------------x_0
  y[1]---------------y_1
  y[2]---------------x_1
  y[3]---------------p_1
  y[2*i]-------------x_{i}
  y[2*i+1]-----------p_{i}
  y[2*N]-------------x_N
  y[2*N+1]-----------p_N
  y[2*(N+1)]---------x_{N+1}
  y[2*(N+1)+1]-------y_N
  y[2*(N+1)+2]-------y_b// yb(0)=ya(0)=0: see expcos.jpg for definition of y_a and y_b
  y[2*(N+1)+3]-------y_a
  N is num_of_atoms
  */
  int N=nvar/2-3;
  int i;
  double J=0.0;
  double F=0.0;
  //  printf("%f ",x);
  double k_cons=2*a*a*D;
  for(i=1;i<N;i++)
  {
    double r1=y[2*i+2]-y[2*i]-x_eq;//r1=x_{i+1}-x_i-x_eq
    F=-r1*k_cons;//see eq.1 and eq.2@ hamiltonian.jpg
    double Jk=F*(y[2*i+1]+y[2*(i+1)+1])/mass;//v=p/m
    J+=Jk;
    //    printf("%08f ",Jk);
  }
  //  printf("\n");
  *rlt=J/(2*(N-1));
}

