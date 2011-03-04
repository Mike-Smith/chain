#ifndef _MY_TIME_H
#define _MY_TIME_H
extern int __WRITE_DATA_TO_FILE__;
//0 nowhere, 1 stdout, 2 datafile

typedef struct RUN_TIME{
	int year;
	int month;
	int day;
	int hour;
	int minute;
	int sec;
	char weekday[4];
	long int sys_sec;
	long int sys_usec;
	char comment[128];
} run_time;


int myprint(char* output_file_name, char* comment, char *format, void * para);
int mytime(run_time *a);
#endif


