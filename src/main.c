
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rnruni.h"
#include "mytime.h"
#include "chain.h"

//0 nowhere, 1 stdout, 2 datafile



int main(int argc, char * argv[])
{

  if (argc==1)
  {
  	 printf("need input parameter(s)");
  	 exit(0);
  }
  if(0==strcmp(argv[1],"-no_output")) __WRITE_DATA_TO_FILE__=0;
  else if (0==strcmp(argv[1],"-stdout")) __WRITE_DATA_TO_FILE__=1;
  else if(0==strcmp(argv[1],"-data_file")) __WRITE_DATA_TO_FILE__=2;
  else if(0==strcmp(argv[1],"-digested")) __WRITE_DATA_TO_FILE__=3;
  else
  {
    printf("wrong input parameters\n");
    exit(0);
  }
  //--------------------------------------------------------------------
  run_time pro_begin;
  mytime(&pro_begin);
  izy+=pro_begin.year*pro_begin.month*pro_begin.day*pro_begin.hour
       *pro_begin.minute*pro_begin.sec;
  char output_file_name[256]={};
  // #if __WRITE_DATA_TO_FILE__ == 0
  if (__WRITE_DATA_TO_FILE__ == 0)
  {
    //    output_file_name[]={};
  }
  // #elif __WRITE_DATA_TO_FILE__ == 1
  else if(__WRITE_DATA_TO_FILE__ == 1)
  {
    strcpy(output_file_name,"stdout");
  }
  // #elif __WRITE_DATA_TO_FILE__ == 2
  else if(__WRITE_DATA_TO_FILE__ == 2 || __WRITE_DATA_TO_FILE__ == 3)
  {
    //     char output_file_name[256]={};
    {
      //     char directory[256]={};
      //     getcwd (directory, 256);
      //specific output data file name by the starting time of the programme.
      char name_tag[256]={};
      if(argv[1]!=0) strcpy(name_tag,argv[1]);
      sprintf(output_file_name,"z%d_%d_%d_%d_%d_%d_%s",
              pro_begin.year,pro_begin.month,
              pro_begin.day,pro_begin.hour,
              pro_begin.minute,pro_begin.sec,name_tag);
      printf("%s\n",output_file_name);
    }
  }
  // #else
  else
  {
    printf("you use wrong main.c");
    exit(0);
  }
  // #endif


  {
    int nt=0;
    double dt=0;
    int good_orb;
    int nobvs=0;
    double *rlt=0;//(double*)malloc(sizeof(double)*nt*2);


//    chain(&dt, &nt, &rlt,&good_orb, &nobvs, output_file_name,argc,argv);
//    chain_spectral(&dt, &nt, &rlt,&good_orb, &nobvs, output_file_name,argc, argv);
    chain_expcos(&dt, &nt, &rlt,&good_orb, &nobvs, output_file_name,argc,argv);
//      chain_ana(&dt, &nt, &rlt,&good_orb, &nobvs, output_file_name,argc,argv);

    int i;
    for(i=0;i<nt*nobvs;i++)
    {
      rlt[2*i]/=good_orb;
      //the average of all good orbs
      rlt[2*i+1]=
        sqrt(fabs(rlt[2*i+1]/good_orb-rlt[2*i]*rlt[2*i])
             /good_orb);
      //the stand errors of the results
    }
    char fn[128]={};
    FILE *fp=0;
    // #if __WRITE_DATA_TO_FILE__ == 1
    if (__WRITE_DATA_TO_FILE__ == 1)
    {
      fp=stdout;
      // #endif
    }
    // #if __WRITE_DATA_TO_FILE__ == 2
    if(__WRITE_DATA_TO_FILE__ == 2 || __WRITE_DATA_TO_FILE__==3)
    {
      strcpy(fn,output_file_name);
      fp = fopen(strcat(fn,".dat"),"a");
    }
    // #endif
    myprint(output_file_name,"#number of good orbit","%d",&good_orb);

    //time average
    double *time_av=(double*)malloc(sizeof(double)*nobvs*2);
    bzero(time_av,sizeof(double)*nobvs*2);
    int nskip=nt/4;//skip some time steps so that the steady state is achieved.
    myprint(output_file_name,"#number of skipped time steps","%d",&nskip);
    if(nskip>=nt)
    {
      printf("nskip should be less than nt\n");
      exit(0);
    }
    for(i=nskip;i<nt;i++)
    {
      int j;
      for(j=0;j<nobvs*2;j++)
      {
        time_av[j]+=rlt[i*2*nobvs+j];
      }
    }
    // #if __WRITE_DATA_TO_FILE__ != 0
    if(__WRITE_DATA_TO_FILE__ != 0)
    {
      fprintf(fp,"\n#observable(s)");
      for(i=0;i<nobvs*2;i++)
      {
        fprintf(fp,"\t%e",time_av[i]/(nt-nskip));
        //       printf("%e\t",time_av[i]/(nt-nskip));
      }
      fprintf(fp,"\n");
    }
    // #elif __WRITE_DATA_TO_FILE__ == 0
    else if(__WRITE_DATA_TO_FILE__ == 0)
    {
      for(i=0;i<nobvs*2;i++)
      {
        printf("%e\t",time_av[i]/(nt-nskip));
      }
    }
    // #else
    else
    {
      printf("__WRITE_DATA_TO_FILE__ == %d\n",__WRITE_DATA_TO_FILE__);
      // #endif
    }
    free(time_av);

    if(__WRITE_DATA_TO_FILE__== 1 || __WRITE_DATA_TO_FILE__==2)
    {
      // #if __WRITE_DATA_TO_FILE__ != 0
      double t=dt;
      //the data is stored in the array in this sequence:
      //t1                                        t2                                          t3
      //o1  v1  o2   v2   o3   v3   o4     v4     o1   v1   o2   v2    o3   v3    o4    v4   o1     v1    o2     v2     o3    v3.....
      //o1 is observable1 at time t1, v1 is its variance. After all observables and variance for t1 is listed, comes observables and their variance of t2, and so on.
      int n_pick=nt/(1024*8*32);
      if(n_pick==0) n_pick=1;
      for(i=0,t=0;i<nt;i+=n_pick,t+=dt*n_pick)
      {
        fprintf(fp,"\n%f",t);
        int j;
        for(j=0;j<nobvs;j++)
        {
          fprintf(fp,"\t%e",rlt[i*2*nobvs+2*j]);
          fprintf(fp,"\t%e",rlt[i*2*nobvs+2*j+1]);
        }
      }
      fclose(fp);
      free(rlt);
      // #endif
    }
  }


  return 0;
}

