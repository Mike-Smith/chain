#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "math.h"
#include "integration.h"


complex99 OMEGA;
complex99 GAMMA;
complex99 *rk_dym_d,*rk_dyt_d,*rk_yt_d;


void euler_one_step ( complex99 y[], complex99 dydx[], int n,
                      double x, double h, complex99* rdm, complex99 yout[],
                      void ( drv ) ( double,complex99*,int,complex99[],complex99[] ) )
{

  //-----add to debug
  /*
  	int i;
  	for(i=0;i<n;i++)
  		y[i]=(complex99)i;
  	(*drv)(x,rdm,n,y,dydx);
  	FILE *fp = fopen("check_derive.dat","w");
  	for(i=0;i<n;i++)
  		fprintf(fp,"%d\t%f\t%f\t%f\t%f\n",i/4,creal(y[i]),cimag(y[i]),
  				creal(dydx[i]),cimag(dydx[i]));
  	fclose(fp);
  	exit(0);
  */
  //-----add to debug
  ( *drv ) ( x,rdm,n,y,dydx );
  int i;
  for ( i=0;i<n;i++ )
    yout[i]=y[i]+h*dydx[i];
  return;
}



void euler_propagator ( complex99 vstart[], int nvar, double x1,
                        double x2, int nstep, complex99 *noise, double *rlt,
                        void ( *drv ) ( double,complex99*,int,complex99[],complex99[] ),
                        void ( *dt_an ) ( double,int,complex99*,double* ) )
{
  complex99 *y  = ( complex99 * ) malloc ( sizeof ( complex99 ) *nvar );
  complex99 *dydx = ( complex99 * ) malloc ( sizeof ( complex99 ) *nvar );
  if ( !y || !dydx )
  {
    perror ( "failed to malloc arrays.\n" );
    exit ( 1 );
  }

  ( *dt_an ) ( x1, nvar, vstart, rlt );

  memcpy ( y, vstart, sizeof ( complex99 ) *nvar );

  double x=x1;
  double h= ( x2-x1 ) /nstep;


  double *prlt = rlt;
  complex99 *rdm = noise;

  int k;
  for ( k=0;k<nstep;k++ )
  {
    euler_one_step ( y,dydx,nvar,x,h,rdm,vstart,drv );
    x+=h;
    rdm++;
    ( *dt_an ) ( x, nvar, vstart, prlt );
    prlt++;
    memcpy ( y, vstart, sizeof ( complex99 ) *nvar );
  }
  free ( y );
  free ( dydx );
}



void rk4_one_step ( complex99 y[], complex99 dydx[], int n,
                    double x, double h, complex99* rdm, complex99 yout[],
                    void ( *derivs ) ( double,complex99*,int,complex99[],complex99[] ) )
{
  int i;
  double xh,hh,h6;

  hh=h*0.5;
  h6=h/6.0;
  xh=x+hh;
  ( *derivs ) ( x,rdm,n,y,dydx );
  for ( i=0;i<n;i++ ) rk_yt_d[i]=y[i]+hh*dydx[i];
  ( *derivs ) ( xh,rdm+1,n,rk_yt_d,rk_dyt_d );
  for ( i=0;i<n;i++ ) rk_yt_d[i]=y[i]+hh*rk_dyt_d[i];
  ( *derivs ) ( xh,rdm+1,n,rk_yt_d,rk_dym_d );
  for ( i=0;i<n;i++ )
  {
    rk_yt_d[i]=y[i]+h*rk_dym_d[i];
    rk_dym_d[i] += rk_dyt_d[i];
  }
  ( *derivs ) ( x+h,rdm+2,n,rk_yt_d,rk_dyt_d );
  for ( i=0;i<n;i++ )
    yout[i]=y[i]+h6* ( dydx[i]+rk_dyt_d[i]+2.0*rk_dym_d[i] );
}


int rk_kmax, rk_kount;

void rkdumb ( complex99 vstart[], int nvar, double x1, double x2,
              int nstep, complex99 *noise, double *rlt,
              void ( *derivs ) ( double,complex99*,int,complex99[],complex99[] ),
              void ( *data_ana ) ( double,int,complex99*,double* ) )
{
  //in rkdumb, noise should be twice as many as in euler.
  int k;
  double x,h;
  complex99 *v,*dv;

  rk_dym_d = ( complex99 * ) malloc ( sizeof ( complex99 ) *nvar );
  rk_dyt_d = ( complex99 * ) malloc ( sizeof ( complex99 ) *nvar );
  rk_yt_d  = ( complex99 * ) malloc ( sizeof ( complex99 ) *nvar );
  if ( !rk_dym_d || !rk_dyt_d || !rk_yt_d )
  {
    perror ( "failed to malloc arrays.\n" );
    exit ( 1 );
  }

  if ( rk_kmax )
  {
    rk_kount = 0;
    ( *data_ana ) ( x1,nvar,vstart,rlt );
    rk_kount++;
  }
  if ( x1 == x2 ) return;

  v  = ( complex99 * ) malloc ( sizeof ( complex99 ) *nvar );
  dv = ( complex99 * ) malloc ( sizeof ( complex99 ) *nvar );
  if ( !v || !dv )
  {
    perror ( "failed to malloc arrays.\n" );
    exit ( 1 );
  }
  memcpy ( v, vstart, sizeof ( complex99 ) *nvar );
  x=x1;
  h= ( x2-x1 ) /nstep;
  complex99 *rdm=noise;
  for ( k=0;k<nstep;k++ )
  {
    rk4_one_step ( v,dv,nvar,x,h,rdm,vstart,derivs );
    ( *data_ana ) ( x,nvar,vstart,rlt );
    x += h;
    rdm+=2;
    rlt++;
    rk_kount++;
    memcpy ( v,vstart,sizeof ( complex99 ) *nvar );
  }

  free ( dv );
  free ( v );
  free ( rk_dym_d );
  free ( rk_dyt_d );
  free ( rk_yt_d );
}
void derive_example ( double x, complex99* rdm ,int nvar, complex99 y[],
                      complex99 dydx[] )
{
  *dydx= ( * ( y+1 ) ) *I;
  * ( dydx+1 ) = ( *y ) *I;
}
complex99 *dydx, *dym, *dyt, *yt;
int alloc_mem_for_propagator ( int ny )
{
  dydx= ( complex99* ) malloc ( sizeof ( complex99 ) *ny );
  dym= ( complex99* ) malloc ( sizeof ( complex99 ) *ny );
  dyt= ( complex99* ) malloc ( sizeof ( complex99 ) *ny );
  yt = ( complex99* ) malloc ( sizeof ( complex99 ) *ny );
  assert ( dydx&&dym&&dyt&&yt );
  return 0;
}

int free_mem_for_propagator()
{
  free ( dydx );
  free ( dym );
  free ( dyt );
  free ( yt );
  return 0;
}

int propagator ( complex99 y[], int ny, double x1,
                 double x2, int nx, double *rlt, complex99 *rdm,
                 void ( *derivs ) ( int, double,complex99*,complex99 [],complex99 [] ),
                 void ( *data_ana ) ( int, double, complex99 [], double * ) )
{


  double h= ( x2-x1 ) /nx;
  double hh=h*0.5, h6=h/6, x, xh;
  int xn, i;
#define SAVE_INTERMIDIATE_DATA 0
#if SAVE_INTERMIDIATE_DATA==1
  FILE *fp = fopen ( "tmp_data_for_rkdumb.dat","w" );
#endif

  for ( xn=0, x=x1+h; xn<nx; xn++, x+=h, rlt++ )
  {
    xh=x+hh;
    ( *derivs ) ( ny, x, rdm, y, dydx );
    for ( i=0;i<ny;i++ ) yt[i]=y[i]+hh*dydx[i];
    ( *derivs ) ( ny, xh, rdm+1, yt, dyt );
    for ( i=0;i<ny;i++ ) yt[i]=y[i]+hh*dyt[i];
    ( *derivs ) ( ny, xh, rdm+1, yt, dym );
    for ( i=0;i<ny;i++ )
    {
      yt[i]=y[i]+h*dym[i];
      dym[i]+=dyt[i];
    }
    ( *derivs ) ( ny, x+h, rdm+2, yt, dyt );
    rdm+=2;
    for ( i=0;i<ny;i++ )
      y[i]+=h6* ( dydx[i]+dyt[i]+2.0*dym[i] );
    //rk4 beyond
    ( *data_ana ) ( ny,x,y,rlt );
#if SAVE_INTERMIDIATE_DATA==1
    fprintf ( fp, "%1.10e\t%1.10e\n", x, *rlt );
#endif

  }
#if SAVE_INTERMIDIATE_DATA==1
  fclose ( fp );
#endif
#undef SAVE_INTERMIDIDATE_DATA
  return 0;
#undef SAVE_INTERMIDIATE_DATA
}

int propagator_purtub ( complex99 y[], int ny, double x1,
                        double x2,int nx, complex99 *rdm,
                        void ( *derivs ) ( int, double,complex99*,complex99 [],complex99 [] ),
                        complex99 *rho0,complex99 *rho1 )
{
  //not begin writing yet
  return 0;

}
void propagator_derive ( int ny, double x, complex99 *rdm,
                         complex99* y, complex99* dydx )
{
  dydx[0]=- ( 1+I ) *y[0];
}
void propagator_data_ana ( int ny, double x, complex99 *y,double *rlt )
{
  *rlt = creal ( *y );
}

void data_analyse_example ( double x, int nvar,
                            complex99 *y, double* rlt )
{
  *rlt = conj ( *y ) * ( *y )-conj ( * ( y+1 ) ) * ( * ( y+1 ) );
}




int accumulate_rlt ( double *rlt_tmp,
                     double *stat,int *good_orbit, int nt )
{
  //accumulate data of every simulation, the result will be stored in stat to calculate average value and statistical variance later.
  int i;
  // //
  // // 	for(i=0;i<nt;i++)
  // // 	{
  // // 		if( ! (rlt_tmp[i]<100000 && rlt_tmp[i]>-100000))
  // // 			break;// result explodes, drop it.
  // // 	}
  // 	if(i!=nt) return 1;//result explodes, drop it.

  for ( i=0;i<nt;i++ )
  {
    double tmp=rlt_tmp[i];
    stat[2*i]+=tmp;
    stat[2*i+1]+=tmp*tmp;
  }
  ( *good_orbit ) ++;
  return 0;
}
int accumulate_rlt_mult ( double *rlt_tmp,
                          double *stat,int *good_orbit, int nt, int obvs )
{
  //the multiple-observable version of  accumulate_rlt
  int i;
  int NN=obvs*nt;
  // //
  // // 	for(i=0;i<nt;i++)
  // // 	{
  // // 		if( ! (rlt_tmp[i]<100000 && rlt_tmp[i]>-100000))
  // // 			break;// result explodes, drop it.
  // // 	}
  // 	if(i!=nt) return 1;//result explodes, drop it.

  for ( i=0;i<NN;i++ )
  {
    double tmp=rlt_tmp[i];
    stat[2*i]+=tmp;
    stat[2*i+1]+=tmp*tmp;
  }
  ( *good_orbit ) ++;
  return 0;
}
double *dydx_r, *dym_r, *dyt_r, *yt_r;
int propagator_real ( double y[], int ny, double x1,
                      double x2, int nx, double *rlt, double *rdm, int nobv,
                      void ( *derivs ) ( int, double, double*,double [],double [] ),
                      void ( *data_ana ) ( int, double, double [], double *, int ) )
{


  double h= ( x2-x1 ) /nx;
  double hh=h*0.5, h6=h/6.0, x, xh;
  int xn, i;
#define SAVE_INTERMIDIATE_DATA 0
#if SAVE_INTERMIDIATE_DATA==1
  FILE *fp = fopen ( "tmp_data_for_propagator_real.dat","a" );
#endif

  for ( xn=0, x=x1; xn<nx; xn++ )
  {
    xh=x+hh;
    ( *derivs ) ( ny, x, rdm, y, dydx_r );
    //    (*data_ana)(ny,x,y,rlt,nobv);//to calculate the observables here is to avoid to calculate the force twice.
    for ( i=0;i<ny;i++ ) yt_r[i]=y[i]+hh*dydx_r[i];
    ( *derivs ) ( ny, xh, rdm+2, yt_r, dyt_r );
    for ( i=0;i<ny;i++ ) yt_r[i]=y[i]+hh*dyt_r[i];
    ( *derivs ) ( ny, xh, rdm+2, yt_r, dym_r );
    for ( i=0;i<ny;i++ )
    {
      yt_r[i]=y[i]+h*dym_r[i];
      dym_r[i]+=dyt_r[i];
    }
    ( *derivs ) ( ny, x+h, rdm+4, yt_r, dyt_r );

    //move the the new point
    for ( i=0;i<ny;i++ )
      y[i]+=h6* ( dydx_r[i]+dyt_r[i]+2.0*dym_r[i] );
    x+=h;
    rdm+=4;
    //rk4 beyond
    ( *data_ana ) ( ny,x,y,rlt,nobv );
    rlt+=nobv;
#if SAVE_INTERMIDIATE_DATA==1
    propagator_real_save_intermidiate_data ( y,fp,ny );
#endif

  }
#if SAVE_INTERMIDIATE_DATA==1
  fclose ( fp );
#endif
#undef SAVE_INTERMIDIDATE_DATA
  return 0;

}



int rk2_white ( double y[], int ny, double x1,
                      double x2, int nx, double *rlt, double *rdm, int nobv,
                      void ( *derivs ) ( int, double, double*,double [],double [] ),
                      void ( *data_ana ) ( int, double, double [], double *, int ) )
{

//this is from the paper Physical Review A, Volume 45, Number 2
//  double h= ( x2-x1 ) /nx;
//  double hh=h*0.5, h6=h/6.0, x, xh;
//  int xn, i;
//#define SAVE_INTERMIDIATE_DATA 0
//#if SAVE_INTERMIDIATE_DATA==1
//  FILE *fp = fopen ( "tmp_data_for_propagator_real.dat","a" );
//#endif
//
//  for ( xn=0, x=x1; xn<nx; xn++ )
//  {
//    xh=x+hh;
//    ( *derivs ) ( ny, x, rdm, y, dydx_r );
//    //    (*data_ana)(ny,x,y,rlt,nobv);//to calculate the observables here is to avoid to calculate the force twice.
//    for ( i=0;i<ny;i++ ) yt_r[i]=y[i]+hh*dydx_r[i];
//    ( *derivs ) ( ny, xh, rdm+2, yt_r, dyt_r );
//    for ( i=0;i<ny;i++ ) yt_r[i]=y[i]+hh*dyt_r[i];
//    ( *derivs ) ( ny, xh, rdm+2, yt_r, dym_r );
//    for ( i=0;i<ny;i++ )
//    {
//      yt_r[i]=y[i]+h*dym_r[i];
//      dym_r[i]+=dyt_r[i];
//    }
//    ( *derivs ) ( ny, x+h, rdm+4, yt_r, dyt_r );
//
//    //move the the new point
//    for ( i=0;i<ny;i++ )
//      y[i]+=h6* ( dydx_r[i]+dyt_r[i]+2.0*dym_r[i] );
//    x+=h;
//    rdm+=4;
//    //rk4 beyond
//    ( *data_ana ) ( ny,x,y,rlt,nobv );
//    rlt+=nobv;
//#if SAVE_INTERMIDIATE_DATA==1
//    propagator_real_save_intermidiate_data ( y,fp,ny );
//#endif
//
//  }
//#if SAVE_INTERMIDIATE_DATA==1
//  fclose ( fp );
//#endif
//#undef SAVE_INTERMIDIDATE_DATA
  return 0;

}
int euler_real ( double y[], int ny, double x1,
                 double x2, int nx, double *rlt, double *rdm, int nobv,
                 void ( *derivs ) ( int, double, double*,double [],double [] ),
                 void ( *data_ana ) ( int, double, double [], double *, int ) )
{


  double h= ( x2-x1 ) /nx;
  double hh=h*0.5, x, xh;
  int xn, i;
#define SAVE_INTERMIDIATE_DATA 0
#if SAVE_INTERMIDIATE_DATA==1
  FILE *fp = fopen ( "tmp_data_for_propagator_real.dat","a" );
#endif

  for ( xn=0, x=x1; xn<nx; xn++ )
  {
    xh=x+hh;
    ( *derivs ) ( ny, x, rdm, y, dydx_r );
    ( *data_ana ) ( ny,x,y,rlt,nobv );//to calculate the observables here is to avoid to calculate the force twice.
    for ( i=0;i<ny;i++ ) y[i]=y[i]+hh*dydx_r[i];
    ( *derivs ) ( ny, xh, rdm+2, y, dydx_r );

    //move to the new point
    x+=h;
    for ( i=0;i<ny;i++ ) y[i]=y[i]+hh*dydx_r[i];
    rdm=rdm+4;
    rlt+=nobv;

#if SAVE_INTERMIDIATE_DATA==1
    propagator_real_save_intermidiate_data ( y,fp,ny );
#endif

  }
#if SAVE_INTERMIDIATE_DATA==1
  fclose ( fp );
#endif
#undef SAVE_INTERMIDIDATE_DATA
  return 0;

}

int propagator_real_save_intermidiate_data ( double *y, FILE *fp,int ny )
{
  int ll;
  for ( ll=0;ll<ny;ll++ )
    fprintf ( fp, "%1.5f ", y[ll] );
  fprintf ( fp,"\n" );
  return 0;
}
int alloc_mem_for_propagator_real ( int ny )
{
  dydx_r= ( double* ) malloc ( sizeof ( double ) *ny );
  dym_r= ( double* ) malloc ( sizeof ( double ) *ny );
  dyt_r= ( double* ) malloc ( sizeof ( double ) *ny );
  yt_r = ( double* ) malloc ( sizeof ( double ) *ny );
  assert ( dydx_r&&dym_r&&dyt_r&&yt_r );
  return 0;
}

int free_mem_for_propagator_real()
{
  free ( dydx_r );
  free ( dym_r );
  free ( dyt_r );
  free ( yt_r );
  return 0;
}

int alloc_mem_for_rk2_white ( int ny )
{
	  dydx_r= ( double* ) malloc ( sizeof ( double ) *ny );
	  dym_r= ( double* ) malloc ( sizeof ( double ) *ny );
	  dyt_r= ( double* ) malloc ( sizeof ( double ) *ny );
	  yt_r = ( double* ) malloc ( sizeof ( double ) *ny );
	  assert ( dydx_r&&dym_r&&dyt_r&&yt_r );
	  return 0;
}
int free_mem_for_rk2_white()
{
  free ( dydx_r );
  free ( dym_r );
  free ( dyt_r );
  free ( yt_r );
  return 0;
}
void data_analyse_example_real ( int nvar,double x,
                                 double *y, double* rlt, int nobv )
{
  *rlt = y[1];
}
void propagator_derive_example_real ( int ny, double x, double *rdm,
                                      double* y, double* dydx )
{
  dydx[0]=- ( 1+2* ( *rdm ) ) * ( 1+2* ( *rdm ) ) *y[0];
  dydx[1]=0;
}

//  int propagator_verlet ( double y[], int ny, double x1,
//                          double x2, int nx, double *rlt, double *rdm, int nobv,
//                          void ( *update_xp ) ( int,int, double, double*,double [],double [] ),
//                          void ( *data_ana ) ( int, double, double [], double *, int ) )
//  {
//    //see wikipedia:Verlet integration
//    // there is substancial difference between Verlet and RK-4. because x and p are treated differently in
//    // Verlet while the same in RK-4. So I need to rearrange the code a little bit.
//    //To follow the frame of the program, Verlet scheme is represent in the following form:
//    //I use x and p instead of r and v
//    // x(t+dt)=x(t)+x'(t) dt + 1/2*p'(t)/m dt^2    --------1
//    // p(t+dt/2)=p(t) + p'(t)*dt/2                  -------2
//    // p'(t+dt)=f(x(t+dt))                          -------3
//    // p(t+dt)=p(t+dt/2)+p'(t+dt)*dt/2               ------4
//    // As x and p are treated differently, there must be a protocol about how to store and process the data,
//    // I don't want to rewrite other part of the programm. So I'll follow the protocol in chain.c. The difference is that now propagator_verlet
//    // also has to know the protocol, while propagator_real doesn't know.
//    /*
//    y[0]---------------x_0 //not used in this function
//    y[1]---------------y_1 //the integration term in Langevin equation. check propagator_derive_chain_UO.jpg
//    y[2]---------------x_1
//    y[3]---------------p_1
//    y[2*i]-------------x_{i}
//    y[2*i+1]-----------p_{i}
//    y[2*N]-------------x_N
//    y[2*N+1]-----------p_N
//    y[2*(N+1)]---------x_{N+1} // not used in this function
//    y[2*(N+1)+1]-------y_N //the integration term in Langevin equation. check propagator_derive_chain_UO.jpg
//    */
//    //pseudo code
//    // calculate x'(0), p'(0) and store the result in dydx_r
//    //
//    // equation 1 and 2               --------loop start
//    // update p'(t+dt) by x(t+dt)
//    // equation 4
//    // update x'(t+dt) by p(t+dt)
//    // goto loop start                   ------------loop end
//    double h= ( x2-x1 ) /nx;
//    double hh=h*0.5, x=x1, xh;
//    int xn, i;
//    int N=ny/2-2;
//    ( *update_xp ) ( 1,ny,x,rdm,y,dydx_r ); //calculate p'(0)
//    ( *update_xp ) ( 0,ny,x,rdm,y,dydx_r ); //calculate x'(0)
//    for ( xn=0, x=x1; xn<nx; xn++ )
//    {
//      for(i=1;i<=N;i++)
//      {
//        y[2*i]=y[2*i]+dydx_r[2*i]*h+0.5*dydx_r[2*i+1]*h*h;
//        y[2*i+1]=y[2*i+1]+0.5*dydx_r[2*i+1]*h;
//      }
//      xh=x+hh;
//    }
//    return 0;
//  }

/*
#include <stdio.h>
int main(void)
{
	derive *drv = derive_example;
	data_ana *dt_an = data_analyse_example;

	complex99 y0[2]={1.0,0.0};
	double *rlt = (double*)malloc(sizeof(double)*10000);

	euler_propagator(y0,2,0,10,10000,rlt,drv,dt_an);

	int i; double t=0.0,dt=10.0/10000.0;
	for(i=0; i<10000; i++,t+=dt)
		printf("%f\t%f\n",t,rlt[i]);
	return 0;
}
*/
//test the propagator
/*	complex99 y[1]={1+1*I};
	double *rlt = (double*)malloc(sizeof(double)*1024);
	alloc_mem_for_propagator(1);
	propagator(y, 1, 0, 10.0, 1024, rlt,
			propagator_derive,propagator_data_ana) ;
	free_mem_for_propagator();
	FILE *fp = fopen("rlt.dat", "w");
	assert(fp);
	int i;
	double t=0;
	dt = 10.0/1024;
	for(i=0,t=dt;i<1024;i++,t+=dt){
			fprintf(fp, "%f\t%f\n",
					t,rlt[i]);
	}
	free(rlt);
	fclose(fp);
*/
