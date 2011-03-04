
#include "mytime.h"
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <complex.h>
int __WRITE_DATA_TO_FILE__ =2;
//0 nowhere, 1 stdout, 2 datafile
int  mytime(run_time *a)
{
  const char* (wday[])={"Sun","Mon","Tue","Wed","Thu","Fri","Sat"};
  time_t timep;
  struct tm *p;
  struct timeval tv;
  time(&timep);
  p=localtime(&timep);

  a->year = -100+p->tm_year;
  a->month = 1+p->tm_mon;
  a->day = p->tm_mday;
  strcpy(a->weekday, wday[p->tm_wday]);
  a->hour = p->tm_hour;
  a->minute = p->tm_min;
  a->sec = p->tm_sec;
  if(gettimeofday(&tv, NULL)==0)
  {
    a->sys_sec = tv.tv_sec;
    a->sys_usec = tv.tv_usec;
  }
  else
  {
    a->sys_sec = 0;
    a->sys_usec = 0;
  }
  return 0;

}

typedef _Complex double complex99;
int myprint(char* output_file_name, char* comment, char *format, void * para)
{
  char empty[]={};
  if (0==strcmp(output_file_name,empty))
  {
    return 0;
  }
  FILE *fp;
  if (0==strcmp(output_file_name,"stdout"))
  {
    fp=stdout;
  }
  else{
    char fn[128]={};
    strcpy(fn,output_file_name);
    fp = fopen(strcat(fn,".dat"),"a");
  }
  fprintf(fp,"\n#%s\t",comment);
  if(strcmp(format,"%f")==0)
  {
    double *p_para = (double*)para;
    fprintf(fp,"%f",*p_para);
  }
  else if(strcmp(format,"%d")==0)
  {
    int *p_para = (int*)para;
    fprintf(fp,"%d",*p_para);
  }
  else if(strcmp(format,"%s")==0)
  {
    char *p_para = (char*)para;
    fprintf(fp,"%s",p_para);
  }
  else if(strcmp(format,"complex99")==0)
  {
    complex99 *p_para = (complex99*)para;
    fprintf(fp,"%f+%fi",creal(*p_para),cimag(*p_para));
  }
  else if(strcmp(format,"time")==0)
  {
    run_time *p_para = (run_time*)para;
    fprintf(fp,"y.m.d\th:m:s\t200%d.%d.%d.\t%d:%d:%d",
            p_para->year,p_para->month,p_para->day,
            p_para->hour,p_para->minute,p_para->sec);
  }
  else
  {
    fprintf(fp,"myprint:format not recog");
  }
  
  if (0!=strcmp(output_file_name,"stdout"))
  {
    fclose(fp);
  }
  return 0;

}
