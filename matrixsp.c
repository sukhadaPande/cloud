#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <pthread.h>
#include <unistd.h>
#include <sys/time.h>


#include <immintrin.h>
#include <stdio.h>


const int no_loops = 1000000;
void perform_FLOPS(int);
void * fpoperations(void *);


#define MAX 4 

#define MAX_THREAD 4

int matA[MAX][MAX]; 
int matB[MAX][MAX]; 
int matC[MAX][MAX]; 

int step_i = 0; 

double randomNumberInRange(int n) {
	return ((n/9.0) * (2.0)) - 1.0;
}



int main(int argc, char** argv)
{
        int no_threads[4] = {1,2,4,8},i;//no of threads array
		for (int i = 0; i < MAX; i++) { 
			for (int j = 0; j < MAX; j++) 
			{
			matA[i][j] = randomNumberInRange(j); 
			matB[i][j] = randomNumberInRange(j); 
			} 
		} 
        size_t threads_size=sizeof(no_threads)/sizeof(int);//size of thread
        printf("no of threads\t\t\tNo of Operations\t\tOperation\t\t\tIOPS/FLOPS\t\tTime(Seconds)\n\n");
        for(i=0;i<threads_size;i++)
        {
                perform_FLOPS(no_threads[i]);//performing double instructions
        }
        return 0;

}

void perform_FLOPS(int no_threads)
{
        struct timeval startTime, endTime;
        double total_time,no_flops;
        // no of iterations
        long double no_operations= 1000000000000;       // no of iterations
        pthread_t *threads = (pthread_t*) malloc(no_threads * sizeof (pthread_t));
        gettimeofday(&startTime, NULL);//getting current time
        int i;
        for (int i = 0; i < no_threads; i++)
        {
		                pthread_create(&threads[i], NULL, fpoperations, (void*)(intptr_t)no_threads);//creating mutliple thread$

        }

        for (int i = 0; i < no_threads; i++)
        {
                pthread_join(threads[i], NULL);//wait till all threads complete execution
        }
        gettimeofday(&endTime, NULL);//getting current time
        total_time =(double) (endTime.tv_sec - startTime.tv_sec)+(double) (endTime.tv_usec - startTime.tv_usec) / 1000000 ;       
		no_flops = ((double)no_operations / (total_time )) / 1e9;//number of flops in giga flops
        printf("%i\t\t\t%Lf\t\t%s\t\t\t%f\t\t%f\n\n", no_threads,no_operations,"Double Operations", no_flops, total_time);        
		free(threads);
}

void * fpoperations(void * arg)
{
		int no_threads = (int)(intptr_t)arg;
		int core = step_i++; 
 
	for (int i = core * MAX / no_threads; i < (core + 1) * MAX / no_threads; i++) 
		for (int j = 0; j < MAX; j++) 
			for (int k = 0; k < MAX; k++) 
				matC[i][j] += matA[i][k] * matB[k][j]; 
		
}