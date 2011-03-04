#include "rnruni.h"
#define azy 16807
#define mzy 2147483647
#define qzy 127773
#define rzy 2836
#define convzy (1.0/(mzy-1))

long int izy=1983475888;

inline double rnruni(void){
/************************************************************************
*   Random numbers uniformly distributed between 0 and 1. izy is seed.  *
*   (coded by Y. Zhang, ITP, Acdemia Sincica, 2000.10. zhangy@itp.ac.cn *
*   can be replaced by any other suitable generator for random numbers) *
*   Put the following line before main():                               *
*   double ranzy(void); long izy=123456789;                             *
*************************************************************************/

  int lzy;
  double rlt;

  lzy=izy/qzy;
  izy=azy*(izy-qzy*lzy)-rzy*lzy;
  if (izy<0) izy+=mzy;
  rlt = convzy*(izy-1);

  return rlt;
}

#undef azy 
#undef mzy 
#undef qzy 
#undef rzy
#undef convzy 
//-nuy
/*
int main(int argc, char *argv[])
{
	int N=1000;
	double x=0;
	double y=0;
	for(int i=0; i<N; i++)
	{
		y=rnruni();
		x+=y;
		printf("%1.10f\n",y); 
	}
	x/=N;
	printf("average is %1.10f\n", x);
	return 0;

}
*/




