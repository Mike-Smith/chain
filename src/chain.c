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



double T_l=300.0;//2.0;//Kelvin;//The temperature of left bath.
double T_r=0.0;//1.0;//Kelvin;//The temperature of right bath.
double tau_l=10.0;//fs//The
double tau_r=10.0;//fs
double epselon_l=5.0e-2;//fs^-1
double epselon_r=5.0e-2;//fs^-1
double mass=12.0;//12.0;//amu
double k_b=8.314472477042e-3;//Boltzman constant
int num_of_atoms=10;
double x_eq=1.54e+2;//pm
double D=3.84e+02/(6.0*6.0);//0.5*0.2700;//3.678e+02/(0.01*0.01);//unit of energy
double a=1.875e-2*(6.0);//1.0;//1.875e-2*(0.01);//pm^-1
//double k_force_constant=2*D*a^2=0.2700;//only for harmonic case
double *memory_of_p1;
double *memory_of_pN;
int N_memory_step_left;
int N_memory_step_right;

//???
int chain(double *dt_of_rlt, int *nt_of_rlt,
		double** rlt_mpi, int *good_orb, int *n_obvs, char* output_file_name, int argc, char *argv[])
{

  if(argc!=(2))
//2 for main -stdout, 4 for 4 parameters for \gamma(t)
  {
    printf("wrong input parameters\n");
    exit(0);
  }
  else
  {
 //   num_of_atoms=atof(argv[2]);
//     tau_l=atof(argv[2]);
//     tau_r=atof(argv[3]);
//     epselon_l=atof(argv[4]);
//     epselon_r=atof(argv[5]);
  }

  int nt  = 1024*128;
  double dt = 0.1;
  double tf = nt*dt;


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
		myprint(output_file_name,"mass","%f",&mass);
		myprint(output_file_name,"nt","%d",&nt);
		myprint(output_file_name,"dt","%f",&dt);
		myprint(output_file_name,"tf","%f",&tf);
		myprint(output_file_name,"function","%s","chain.c");

	}

	double para[5];
	para[0]=T_l;para[1]=mass;para[2]=epselon_l;para[3]=tau_l;para[4]=k_b;
	double *alpha_l = alpha_r_t_UO(2*nt, 0.5*dt, para);//generate memory function
	para[0]=T_r;para[1]=mass;para[2]=epselon_r;para[3]=tau_r;para[4]=k_b;
	double *alpha_r = alpha_r_t_UO(2*nt, 0.5*dt, para);
//	{
//	  FILE *fp=fopen("rlt_ch.dat","w");
//	  int i;
//	  for(i=0;i<2*nt;i++)
//	  {
//	    fprintf(fp,"%f  %e\n",i*dt/2.0,alpha_l[i]);
//	  }
//	  fclose(fp);
//	}

	//alloc memory to handle the integral part of the equation of motion.
// 	N_memory_step_left=memory_truncate(alpha_l,nt);
// 	printf("%d\n",N_memory_step_left);//    print l.split()


// 	exit(0);
// 	N_memory_step_right=memory_truncate(alpha_r,nt);
// 	memory_of_p1=(double*)malloc(sizeof(double)*N_memory_step_left);//alloc memory to store p1 before current time to calculate the integration
// 	memory_of_pN=(double*)malloc(sizeof(double)*N_memory_step_right);//alloc memory to store pN before current time to calculate the integration
// 	assert(memory_of_p1&&memory_of_pN);
// 	bzero(memory_of_p1,sizeof(double)*N_memory_step_left);
// 	bzero(memory_of_pN,sizeof(double)*N_memory_step_right);


	double *lamda_l = lamda_sqrt_t(alpha_l,2*nt,0.5*dt);
	double *lamda_r = lamda_sqrt_t(alpha_r,2*nt,0.5*dt);
	free(alpha_l);
	free(alpha_r);
	corr_rand_gener_prep(2*nt);//open cache for FFTW3 operation680

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
	int num_of_terms=(num_of_atoms+2)*2;
	double *rho=(double*)malloc(sizeof(double)*num_of_terms);
	alloc_mem_for_propagator_real(num_of_terms);
//handle the loop

	int i, N=4*4*100,good_orbit=0;
// 	//initial randomly set momentum
//   	initial_chain(rho,num_of_atoms);
//  	double *equi=(double*)malloc(sizeof(double)*num_of_terms);
//  	for(i=0;i<10;i++)
//  	{
//  	  corr_rand_gener_2bath(2*nt, &rand, lamda_r, lamda_l, i%4);
//  	  propagator_real(rho, num_of_terms, 0.0, tf,
// 			  nt, rlt_tmp,rand,nobvs,
// //			  propagator_derive_chain_UO_harm,
// 			  propagator_derive_chain_UO_mors,
// // 			  data_analyse_chain
// 			  data_analyse_chain_heat_flux_mors
// 			 );
//  // 	  int kk;
//  //  	  for(kk=0;kk<nt;kk++)
//  //  	    printf("%f\n",rlt_tmp[kk]);
//  	}
//  	memcpy(equi,rho,sizeof(double)*num_of_terms);
// 	//\initial randomly set momentum
	run_time time_loop_begin,time_loop_end;
	mytime(&time_loop_begin);
	{
	  myprint(output_file_name,"num of orbs","%d",&N);
	  myprint(output_file_name,"loop begin time","time",&(time_loop_begin));
	}
	for(i=0;i<N;i++)
	{
// 		if(i%100==0)
// 		  printf("\n%d orb left\n", (int)((N-i)));

//		corr_rand_gener(2*nt,&rand,lamda_l,i%4);
		corr_rand_gener_2bath(2*nt, &rand, lamda_r, lamda_l, i%4);

// 		white_rand_gener_2bath(2*nt, dt/2.0, T_r, T_l, &rand, i%4);

//		bzero(rand,sizeof(double)*2*nt);
		//output noise
// 		int j;
// 		for( j=0;j<2*nt;j++)
// 		  printf("%f\t%f\t%f\n",j*dt/2.0,rand[2*j],rand[2*j+1]);
// 		exit(0);
//test corr
//  		int j;
// 		int Nn=1024;//1024*56;
// 		int TT=nt-2*1024;
// 		for(j=0;j<TT;j++)zotero://attachment/13087/
// 		  rlt_tmp[j]=rand[2*Nn+1]*rand[2*Nn+1+4*j];//*dt*dt/4.0;//*rand[1];

//		rho[0]=1;rho[1]=1;
 		initial_chain(rho,num_of_atoms);
// 		memcpy(rho,equi,sizeof(double)*num_of_terms);

		propagator_real(rho, num_of_terms, 0.0, tf,
				nt, rlt_tmp,rand,nobvs,
//				propagator_derive_chain_white_harm,
//  				propagator_derive_chain_UO_harm,
 				propagator_derive_chain_UO_mors,
// 				data_analyse_chain_heat_flux_fast
//  				data_analyse_chain_temperature_profile
//				data_analyse_chain_heat_flux_harm_mult
// 				data_analyse_chain_heat_flux_harm
				data_analyse_chain_heat_flux_mors
// 				data_analyse_chain_position_mult
// 				data_analyse_chain_heat_flux_fast_mult
// 				data_analyse_chain_total_energy_mors
                              );

//  		euler_real(rho, num_of_terms, 0.0, tf,
//  				nt, rlt_tmp,rand,nobvs,
//  				propagator_derive_chain_white_harm,
// 				propagator_derive_chain_UO_harm,
// 				propagator_derive_chain_UO_mors,
// 				data_analyse_chain_heat_flux_fast
//   				data_analyse_chain_temperature_profile
// 				data_analyse_chain_position_mult
// 				data_analyse_chain_heat_flux_harm_mult
// 				data_analyse_chain_heat_flux_harm
//  				data_analyse_chain_heat_flux_fast_mult
// 				data_analyse_chain_total_energy_harm
//  			  );

// 		int j;
// 		rlt_tmp[0]=1;
// 		for(j=0;j<nt-1;j++){
// 		  double mid=rlt_tmp[j]-(1+rand[4*j+1])*(1+rand[4*j+1])*rlt_tmp[j]*dt/2;
// 		  rlt_tmp[j+1]=mid-(1+rand[4*j+3])*(1+rand[4*j+3])*mid*dt/2;
// 		}
		// 		accumulate_rlt(rlt_tmp,rlt,&good_orbit, nt);zotero://attachment/13087/
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

double *alpha_r_t_UO(int nt, double dt, double *para)
{//return the memory function
  //para[0]=T, para[1]=m, para[2]=\epselon, para[3]=\tau_c,para[4]=k_b. see Eq. 2 in JCP 128, 224710 (2008)
  int i;
//  printf("%f\t%f\t%f\t%f\t%f\n",para[0],para[1],para[2],para[3],para[4]);
//  char fn[1024];
//  sprintf(fn,"alpha_r_t_UO%1e%1e%1e%1e%1e.dat",para[0],para[1],para[2],para[3],para[4]);
//  FILE *fp=fopen(fn,"w");
  double *rlt = malloc(sizeof(double)*nt);
  for(i=0;i<nt;i++)
  {//rlt=(k_b T_n /m) (\epselon/\tau) exp(-t/\tau)
    rlt[i]=(para[4]*para[0]/para[1])*(para[2]/para[3])*exp(-dt*i/para[3]);
//    fprintf(fp,"%f\t%f\n",dt*i,rlt[i]);
  }
//   fclose(fp);
  return rlt;
}



double *corr_rand_gener_2bath_mem;
int corr_rand_gener_2bath_malloc(int nt)
{
  corr_rand_gener_2bath_mem=(double*)malloc(sizeof(double)*8*nt);
  return 0;
}
int corr_rand_gener_2bath_free()
{
  free(corr_rand_gener_2bath_mem);
  return 0;
}

int white_rand_gener_2bath(const int nt, double dt, double T_right, double T_left, double** rand, int flag)
{
    //8 white noises will be generated, they are
  //right[0][0],right[0][1],right[0][2]...--the 1st random for right bath
  //right[1][0],right[1][1],right[1][2]...--the 2nd random for right bath
  //right[2][0],right[2][1],right[2][2]...--the 3rd random for right bath
  //right[3][0],right[3][1],right[3][2]...--the 4th random for right bath
  //left[0][0], left[0][1], left[0][2]...--the 1st random vector for left bath
  //left[1][0], left[1][1], left[1][2]...--the 1st random vector for left bath
  //left[2][0], left[2][1], left[2][2]...--the 1st random vector for left bath
  //left[3][0], left[3][1], left[3][2]...--the 1st random vector for left bath
  //The data above will be stored in pnew array with the sequence
  //right[0][0] left[0][0] right[0][1] left[0][1] ......right[0][nt-1] left[0][nt-1]
  //                                                                         ^
  //                                                                         |
  //                                                                    Number 2*nt element:pnew[2*nt-1]
 //right[1][0] left[1][0] right[1][1] left[1][1] ......right[1][nt-1] left[1][nt-1]
  //                                                                         ^
  //                                                                         |
  //                                                                    Number 4*nt element:pnew[4*nt-1]
   //right[2][0] left[2][0] right[2][1] left[2][1] ......right[2][nt-1] left[2][nt-1]
  //                                                                         ^
  //                                                                         |
  //                                                                    Number 6*nt element:pnew[6*nt-1]
  //right[3][0] left[3][0] right[3][1] left[3][1] ......right[3][nt-1] left[3][nt-1]
  //                                                                         ^
  //                                                                         |
  //                                                                    Number 8*nt element : pnew[8*nt-1]

  //the structure of pfft is
  //color[0][0] color[1][0] color[0][1] color[1][1] .........color[0][nt-1] color[1][nt-1]
  //                                                                                 ^
  //---------------------------------------------------------------------------------|
  //                                                                               pfft[2*nt-1]
  //color[2][0] color[3][0] color[2][1] color[3][1] .........color[2][nt-1] color[3][nt-1]
  //                                                                                 ^
  //---------------------------------------------------------------------------------|
  //                                                                               pfft[4*nt-1]
  //to adapt pfft to pnew, it's neccessary to rearrange the sequence
  double *pnew=corr_rand_gener_2bath_mem;
  assert(pnew);
  double sig_right=sqrt(2.0*k_b*T_right*epselon_r/(mass*dt));
  double sig_left=sqrt(2.0*k_b*T_left*epselon_l/(mass*dt));
  if(flag==0)
  {
    int i;
    for(i=0;i<4*nt;i++)
    {
      double pt, u1, u2;
      do
      {
	u1 = rnruni();
      }
      while (u1<=1.e-15);
      do
      {
	u2 = rnruni();
      }
      while (u2<=1.e-15);
      u2 *= 2*M_PI;
      pt = sqrt(-2.*log(u1));
      pnew[2*i] = pt*cos(u2)*sig_right;
      pnew[2*i+1] = pt*sin(u2)*sig_left;
    }
    *rand=corr_rand_gener_2bath_mem;
  }
  else if(flag==1)
  {
    *rand=corr_rand_gener_2bath_mem+nt;
  }
  else if(flag==2)
  {
    *rand=corr_rand_gener_2bath_mem+2*nt;
  }
  else if(flag==3)
  {
    *rand=corr_rand_gener_2bath_mem+3*nt;
  }
  return 0;

}
int corr_rand_gener_2bath(const int nt, double** rand, double *lamda_r, double *lamda_l, int flag)
{
  //8 color noises will be generated, they are
  //right[0][0],right[0][1],right[0][2]...--the 1st random for right bath
  //right[1][0],right[1][1],right[1][2]...--the 2nd random for right bath
  //right[2][0],right[2][1],right[2][2]...--the 3rd random for right bath
  //right[3][0],right[3][1],right[3][2]...--the 4th random for right bath
  //left[0][0], left[0][1], left[0][2]...--the 1st random vector for left bath
  //left[1][0], left[1][1], left[1][2]...--the 1st random vector for left bath
  //left[2][0], left[2][1], left[2][2]...--the 1st random vector for left bath
  //left[3][0], left[3][1], left[3][2]...--the 1st random vector for left bath
  //The data above will be stored in pnew array with the sequence
  //right[0][0] left[0][0] right[0][1] left[0][1] ......right[0][nt-1] left[0][nt-1]
  //                                                                         ^
  //                                                                         |
  //                                                                    Number 2*nt element:pnew[2*nt-1]
 //right[1][0] left[1][0] right[1][1] left[1][1] ......right[1][nt-1] left[1][nt-1]
  //                                                                         ^
  //                                                                         |
  //                                                                    Number 4*nt element:pnew[4*nt-1]
   //right[2][0] left[2][0] right[2][1] left[2][1] ......right[2][nt-1] left[2][nt-1]
  //                                                                         ^
  //                                                                         |
  //                                                                    Number 6*nt element:pnew[6*nt-1]
  //right[3][0] left[3][0] right[3][1] left[3][1] ......right[3][nt-1] left[3][nt-1]
  //                                                                         ^
  //                                                                         |
  //                                                                    Number 8*nt element : pnew[8*nt-1]

  //the structure of pfft is
  //color[0][0] color[1][0] color[0][1] color[1][1] .........color[0][nt-1] color[1][nt-1]
  //                                                                                 ^
  //---------------------------------------------------------------------------------|
  //                                                                               pfft[2*nt-1]
  //color[2][0] color[3][0] color[2][1] color[3][1] .........color[2][nt-1] color[3][nt-1]
  //                                                                                 ^
  //---------------------------------------------------------------------------------|
  //                                                                               pfft[4*nt-1]
  //to adapt pfft to pnew, it's neccessary to rearrange the sequence

  double *pfft=(double*)fftw3_out_for_corr_rand;
  double *pnew=corr_rand_gener_2bath_mem;
  assert(pfft);assert(pnew);
  if(flag==0)
  {
    corr_rand_gener_bare(nt, &pfft, lamda_r);
    int i;
    for(i=0;i<2*nt;i++)
      pnew[2*i]=pfft[2*i];
    pnew=corr_rand_gener_2bath_mem+4*nt;
    for(i=0;i<2*nt;i++)
      pnew[2*i]=pfft[2*i+1];
    corr_rand_gener_bare(nt, &pfft, lamda_l);
    pnew=corr_rand_gener_2bath_mem;
    for(i=0;i<2*nt;i++)
      pnew[2*i+1]=pfft[2*i];
    pnew=corr_rand_gener_2bath_mem+4*nt;
    for(i=0;i<2*nt;i++)
      pnew[2*i+1]=pfft[2*i+1];
    *rand=corr_rand_gener_2bath_mem;
  }
  else if(flag==1)
  {
    *rand=corr_rand_gener_2bath_mem+nt;
  }
  else if(flag==2)
  {
    *rand=corr_rand_gener_2bath_mem+2*nt;
  }
  else if(flag==3)
  {
    *rand=corr_rand_gener_2bath_mem+3*nt;
  }
  return 0;
}



void propagator_derive_chain_UO_mors(int ny, double x, double *rdm,
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
*/
  int N=ny/2-2;
  //calculate the pair twice
//   int i;
//   for(i=1;i<=N;i++)
//   {
//     dydx[2*i]=y[2*i+1]/mass;//    (d/dt)x_i=p_i/m
//     double r1=y[2*i+2]-y[2*i]-x_eq;//r1=x_{i+1}-x_i-x_eq
//     double r2=y[2*i]-y[2*i-2]-x_eq;//r2=x_{i}-x_{i-1}-x_eqD
//     r1=fastexp(-a*r1);
//     r2=fastexp(-a*r2);
//     dydx[2*i+1]=(r2*(r2-1)-r1*(r1-1))*2*a*D;//see eq.1 and eq.2@ hamiltonian.jpg
//   }
//   dydx[3]=dydx[3]-y[1]+mass*rdm[1];
//   dydx[2*N+1]=dydx[2*N+1]-y[2*N+3]+mass*rdm[0];
//   dydx[1]=(epselon_l/tau_l)*y[3]-y[1]/tau_l;
//   dydx[2*N+3]=(epselon_r/tau_r)*y[2*N+1]-y[2*N+3]/tau_r;//see propagator_derive_chain_UO.jpg
  //calculate force pair once
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
  dydx[1]=(epselon_l/tau_l)*y[3]-y[1]/tau_l;
  dydx[2*N+3]=(epselon_r/tau_r)*y[2*N+1]-y[2*N+3]/tau_r;//see propagator_derive_chain_UO.jpg
#undef fastexp
  //test x_0 periodic and x_{N+1} damping:
//  dydx[0]=0.1*sin(x);
//  dydx[2*N+1]-=y[2*N+1];
}
void propagator_derive_chain_UO_harm(int ny, double x, double *rdm,
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
*/
  int N=ny/2-2;
  int i;
  //calculate force pair twice
//   for(i=1;i<=N;i++)
//   {
//     dydx[2*i]=y[2*i+1]/mass;//    (d/dt)x_i=p_i/m
//     double r1=y[2*i+2]-y[2*i]-x_eq;//r1=x_{i+1}-x_i-x_eq
//     double r2=y[2*i]-y[2*i-2]-x_eq;//r2=x_{i}-x_{i-1}-x_eqD
//     dydx[2*i+1]=(r1-r2)*D*2.0*a*a;
//   }
//   dydx[3]=dydx[3]-y[1]+mass*rdm[1];
//   dydx[2*N+1]=dydx[2*N+1]-y[2*N+3]+mass*rdm[0];
//   dydx[1]=(epselon_l/tau_l)*y[3]-y[1]/tau_l;
//   dydx[2*N+3]=(epselon_r/tau_r)*y[2*N+1]-y[2*N+3]/tau_r;//see propagator_derive_chain_UO.jpg

  //calculate force pair only once.
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
  dydx[1]=(epselon_l/tau_l)*y[3]-y[1]/tau_l;
  dydx[2*N+3]=(epselon_r/tau_r)*y[2*N+1]-y[2*N+3]/tau_r;//see propagator_derive_chain_UO.jpg

  //test x_0 periodic and x_{N+1} damping:
//  dydx[0]=0.1*sin(x);
//  dydx[2*N+1]-=y[2*N+1];
}
void propagator_derive_chain_white_harm(int ny, double x, double *rdm,
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
*/
  int N=ny/2-2;
  int i;
  //calculate force pair twice
//   for(i=1;i<=N;i++)
//   {
//     dydx[2*i]=y[2*i+1]/mass;//    (d/dt)x_i=p_i/m
//     double r1=y[2*i+2]-y[2*i]-x_eq;//r1=x_{i+1}-x_i-x_eq
//     double r2=y[2*i]-y[2*i-2]-x_eq;//r2=x_{i}-x_{i-1}-x_eqD
//     dydx[2*i+1]=(r1-r2)*D*2.0*a*a;
//   }
//   dydx[3]=dydx[3]-y[1]+mass*rdm[1];
//   dydx[2*N+1]=dydx[2*N+1]-y[2*N+3]+mass*rdm[0];
//   dydx[1]=(epselon_l/tau_l)*y[3]-y[1]/tau_l;
//   dydx[2*N+3]=(epselon_r/tau_r)*y[2*N+1]-y[2*N+3]/tau_r;//see propagator_derive_chain_UO.jpg

  //calculate force pair only once.
  double r1;
  bzero(dydx,sizeof(double)*ny);
  double k_cons=D*2.0*a*a;
  for(i=1;i<N;i++)
  {
    dydx[2*i]=y[2*i+1]/mass;//    (d/dt)x_i=p_i/m
    r1=y[2*i+2]-y[2*i]-x_eq;//r1=x_{i+1}-x_i-x_eq
    dydx[2*i+1]+=r1*k_cons;//see eq.1 and eq.2@ hamiltonian.jpg
    dydx[2*(i+1)+1]-=r1*k_cons;
  }
  dydx[2*N]=y[2*N+1]/mass;
  r1=y[2]-y[0]-x_eq;
  dydx[3]-=r1*k_cons;
  r1=y[2*(N+1)]-y[2*N]-x_eq;
  dydx[2*N+1]+=r1*k_cons;
  dydx[3]=dydx[3]-epselon_l*y[3]+mass*rdm[1];
  dydx[2*N+1]=dydx[2*N+1]-epselon_r*y[2*N+1]+mass*rdm[0];


  //test x_0 periodic and x_{N+1} damping:
//  dydx[0]=0.1*sin(x);
//  dydx[2*N+1]-=y[2*N+1];
}

int initial_chain(double *y,int num_of_atoms)
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
  N is num_of_atoms
  */

  y[0]=0;
  y[1]=0;
  y[2*num_of_atoms+2]=(num_of_atoms+1)*x_eq;
  y[2*num_of_atoms+3]=0;
//all at equil position still
  int i;
  for(i=1;i<=num_of_atoms;i++)
  {
    y[2*i]=i*x_eq;
    y[2*i+1]=0;
  }

////   y[2]+=x_eq/30.0;//a perturbation at the initial state

//  //assign by temperature
//   y[0]=0;
//   y[1]=0;
//   y[2*num_of_atoms+2]=(num_of_atoms+1)*x_eq;
//   y[2*num_of_atoms+3]=0;
//   double dT=(T_l-T_r)/num_of_atoms;
//   int i;
//   for(i=1;i<=num_of_atoms;i++)
//   {
//     double T=T_l-i*dT;
//     double sig1=sqrt(k_b*T*mass);
// //    double sig2=sqrt(k_b*T/D*0.5);
//     double a1,a2;
//     whitenoise(sig1,&a1,&a2);
//     y[2*i+1]=a1;
// //    whitenoise(sig2,&a1,&a2);
//     y[2*i]=i*x_eq;//+a1;
//   }
  return 0;
}

void data_analyse_chain(int nvar,double x,
			       double *y, double* rlt, int nobv)
{
//  *rlt = y[21]*y[21]/(2*mass);//average kinetic energy on atom 10
  int i=2;//the flux on i;
  double F=-dyt_r[2*i+1];
  *rlt=F*(y[2*i+1]+y[2*(i+1)+1])/mass;

}

void data_analyse_chain_temperature_profile(int nvar,double x,
			double *y, double* rlt, int nobv)
{
//  int N=nvar/2-2;
  int i;
  for(i=1;i<=nobv;i++)
  {
    rlt[i-1]=(y[2*i+1]*y[2*i+1])/(mass*k_b);
  }

}
void data_analyse_chain_kenergy(int nvar,double x,
			double *y, double* rlt, int nobv)
{//kinetic energy
  double ke=0;
  int N=nvar/2-2;
  int i;
  for(i=1;i<=N;i++)
  {
    ke+=(y[2*i+1]*y[2*i+1]);
  }
  ke=ke/(2*mass);
  *rlt=ke;
}
void data_analyse_chain_penergy_mors(int nvar,double x,
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
  N is num_of_atoms
    */
  double v=0;
  int N=nvar/2-2;
  int i;
  for(i=0;i<=N;i++)
  {
    double r=exp(-a*(y[2*(i+1)]-y[2*i]-x_eq))-1;
    v+=D*r*r;
  }
  *rlt=v;
}
void data_analyse_chain_penergy_harm(int nvar,double x,
				     double *y, double* rlt, int nobv)
{//potential energy
    /*
  y[0]---------------x_0
  y[1]---------------y_1
  y[2]---------------x_1
  y[3]---------------p_1
  y[2*i]-------------x_{i}
  y[2*i+ 1]-----------p_{i}
  y[2*N]-------------x_N
  y[2*N+1]-----------p_N
  y[2*(N+1)]---------x_{N+1}
  y[2*(N+1)+1]-------y_N
  N is num_of_atoms
    */
  double v=0;
  int N=nvar/2-2;
  int i;
  for(i=0;i<=N;i++)
  {
    double r=y[2*(i+1)]-y[2*i]-x_eq;
    v+=0.5*D*r*r;
  }
  *rlt=v;
}
void data_analyse_chain_total_energy_mors(int nvar,double x,
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
  N is num_of_atoms
    */
  double ve=0;
  double ke=0;
  data_analyse_chain_kenergy(nvar,x,y,&ke,nobv);
  data_analyse_chain_penergy_mors(nvar,x,y, &ve,nobv);
  *rlt=ve+ke;
}
void data_analyse_chain_total_energy_harm(int nvar,double x,
					  double *y, double* rlt,  int nobv)
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
  N is num_of_atoms
    */
  double ve=0;
  double ke=0;
  data_analyse_chain_kenergy(nvar,x,y,&ke,nobv);
  data_analyse_chain_penergy_harm(nvar,x,y, &ve,nobv);
  *rlt=ve+ke;
}
void data_analyse_chain_heat_flux_mors(int nvar,double x,
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
  N is num_of_atoms
    */
  int N=nvar/2-2;
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

void data_analyse_chain_heat_flux_harm(int nvar,double x,
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
  N is num_of_atoms
    */
  int N=nvar/2-2;
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

void data_analyse_chain_heat_flux_harm_mult(int nvar,double x,
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
  N is num_of_atoms
    */
  int N=nvar/2-2;
  int i;
//  double J=0.0;
  double F=0.0;
//  printf("%f ",x);
  double k_cons=2*a*a*D;
  for(i=1;i<N;i++)
  {
    double r1=y[2*i+2]-y[2*i]-x_eq;//r1=x_{i+1}-x_i-x_eq
    F=-r1*k_cons;//see eq.1 and eq.2@ hamiltonian.jpg
    rlt[i-1]=F*(y[2*i+1]+y[2*(i+1)+1])/(2.0*mass);//v=p/m

//    printf("%08f ",Jk);
  }
//  printf("\n");
}

// void data_analyse_chain_heat_flux_fast(int nvar,double x,
// 				  double *y, double* rlt, int nobv)
// {//to use dp/dt to calculate force
//   //double *dydx_r, *dym_r, *dyt_r, *yt_r;
//   //To use this, also change the     (*data_ana)(ny,x,y,rlt) in
//  //propagator_real() to (*data_ana)(ny,x,dyt_r,rlt);
//   //potential energy
//     /*
//   y[0]---------------x_0
//   y[1]---------------y_1
//   y[2]---------------x_1
//   y[3]---------------p_1
//   y[2*i]-------------x_{i}
//   y[2*i+1]-----------p_{i}
//   y[2*N]-------------x_N
//   y[2*N+1]-----------p_N
//   y[2*(N+1)]---------x_{N+1}
//   y[2*(N+1)+1]-------y_N
//   N is num_of_atoms
//     */
//   int N=nvar/2-2;
//   int i;
//   double J=0.0;
//   double F=0.0;
// //   printf("%f ",x);
//   for(i=1;i<N;i++)
//   {
//     F=-dydx_r[2*i+1];
//     double Jk=F*(y[2*i+1]+y[2*(i+1)+1])/mass;
// //     printf("%08f ",Jk);
//     J+=Jk;
//   }
// //   printf("\n");
//   *rlt=J/(2*(N-1));
// }

// void data_analyse_chain_heat_flux_fast_mult(int nvar,double x,
// 				       double *y, double* rlt, int nobv)
// {//to use dp/dt to calculate force
//   //double *dydx_r, *dym_r, *dyt_r, *yt_r;
//   //To use this, also change the     (*data_ana)(ny,x,y,rlt) in
//  //propagator_real() to (*data_ana)(ny,x,dyt_r,rlt);
//   //potential energy
//     /*
//   y[0]---------------x_0
//   y[1]---------------y_1
//   y[2]---------------x_1
//   y[3]---------------p_1
//   y[2*i]-------------x_{i}
//   y[2*i+1]-----------p_{i}
//   y[2*N]-------------x_N
//   y[2*N+1]-----------p_N
//   y[2*(N+1)]---------x_{N+1}
//   y[2*(N+1)+1]-------y_N
//   N is num_of_atoms
//     */
//   int i;
//   for(i=1;i<nobv;i++)
//   {
//     double F=-dydx_r[2*i+1];
//     *(rlt+i-1)=F*(y[2*i+1]+y[2*(i+1)+1])/mass;
//   }
//   *(rlt+nobv-1)=0;
// }
void data_analyse_chain_position_mult(int nvar,double x,
					    double *y, double* rlt, int nobv)
{//to use dp/dt to calculate force
  //double *dydx_r, *dym_r, *dyt_r, *yt_r;
  //To use this, also change the     (*data_ana)(ny,x,y,rlt) in
 //propagator_real() to (*data_ana)(ny,x,dyt_r,rlt);
  //potential energy
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
  N is num_of_atoms
    */
  int i;
  for(i=1;i<=nobv;i++)
  {
    *(rlt+i-1)=y[2*i];
  }
}

int memory_truncate(double *alpha, int nt)
{
  int i;
  int rlt=nt;
  double max=alpha[0];
  for(i=0;i<nt;i++)
  {
    if(max<alpha[i]) max=alpha[i];
    if(fabs(alpha[0]/alpha[2*i])>1.0e2)
    {
      rlt=i;
    }
  }
  return rlt;
}
