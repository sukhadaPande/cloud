#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/time.h>
#include <immintrin.h>


void * spoperations(void *);

int main(int argc, char** argv)
{
	printf("Start of program\n");
	int threads_no[4] = {1,2,4,8};
	size_t size_of_thread=sizeof(threads_no)/sizeof(int);
	int i;
	
	for(i=0;i<size_of_thread;i++)
	{
		//calculate_GLOPS(threads_no[i]);
	double flops , gflops,effiency;
	long double no_operations= 1000000000000;
	double total_time,final_time,theoreticalvalue = 79.3;
	int t = threads_no[i];
	printf(" Thread %i" ,t); 
	pthread_t *threads = (pthread_t*) malloc(threads_no[i] * sizeof (pthread_t));
	
	int i;
	struct timeval start_time, end_time;
	gettimeofday(&start_time, NULL);
	for (i = 0; i < threads_no[i]; i++) 
		{
			pthread_create(&threads[i], NULL, spoperations, (void*)(intptr_t)threads_no[i]);
		}
	
	for (i = 0; i < threads_no[i]; i++) 
		{
			pthread_join(threads[i], NULL);
		}
	gettimeofday(&end_time, NULL);
	total_time =(double) (end_time.tv_sec - start_time.tv_sec)+(double) (end_time.tv_usec - start_time.tv_usec)/ 1000000; 
	
	flops = (no_operations / (total_time ));
	gflops = flops/ 1000000000;
	effiency = (gflops/theoreticalvalue)* 100;
	free(threads);
	
	printf("No of threads %i\n" ,t);
	printf("GFLOPS %f\n" , gflops);
	printf("Time taken %f\n",total_time);
	printf("Efficiency %f\n", effiency);

	}
	return 0;
	
}


void * spoperations(void * arg)
{
	int no_threads = (int)(intptr_t)arg;
		__m256 evens = _mm256_set_ps(2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0);
		__m256 odds = _mm256_set_ps(1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0);
		
		long double no_of_operations=1000000000000/8;
        int outerloop = (no_of_operations)/1000000;
        int result1= (outerloop)/no_threads;
		

		for (int z = 0; z < (1000000) ; z++)
		 {
                for(int y =0; y <(result1);y++)
                {
					
				__m256 result = _mm256_sub_ps(evens, odds);
				
                }
        }
}