#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <pthread.h>
#include <unistd.h>
#include <sys/time.h>


#include <immintrin.h>
#include <stdio.h>
#define SIZE 1000




float matA[SIZE][SIZE]; 
float matB[SIZE][SIZE]; 
float transB[SIZE][SIZE];
float matC[SIZE][SIZE];
int step = 0;  

void GFLOPS(int);
void * matrixmsp(void *);



float randomNo(int n) {
	return ((n/9.0) * (2.0)) - 1.0;
}


int main(int argc, char** argv)
{
		int count =0;
		int value;
        int no_threads[4] = {1,2,4,8},i;//no of threads array
		for (int i = 0; i < SIZE; i++) { 
			for (int j = 0; j < SIZE; j++) 
			{
			value = count++;
			matA[i][j] = randomNo(rand() % 10); 
			matB[i][j] = randomNo(rand() % 10); 
			} 
		} 
        size_t threads_size=sizeof(no_threads)/sizeof(int);
		 
		for (int i = 0; i <SIZE; i++) 
			for (int j = 0; j < SIZE; j++) 
				transB[i][j] = matB[j][i];

        for(i=0;i<threads_size;i++)
        {
                GFLOPS(no_threads[i]);
        }
        return 0;

}

void GFLOPS(int no_threads)
{
        struct timeval startTime, endTime;
        double total_time,gflops,efficiency,flops,theoriticalValue = 79.3;
        // no of iterations
        long double no_operations= 2000000000000;       	
        
        gettimeofday(&startTime, NULL);
		pthread_t *threads = (pthread_t*) malloc(no_threads * sizeof (pthread_t));
        int i;
        for (int i = 0; i < no_threads; i++)
        {
		    pthread_create(&threads[i], NULL, matrixmsp, (void*)(intptr_t)no_threads);

        }

        for (int i = 0; i < no_threads; i++)
        {
            pthread_join(threads[i], NULL);
        }
        gettimeofday(&endTime, NULL);
        total_time =(double) (endTime.tv_sec - startTime.tv_sec)+(double) (endTime.tv_usec - startTime.tv_usec) / 1000000 ;
		flops =  ((double)no_operations / (total_time ));
		gflops = (flops) / 1000000000;
		efficiency = gflops/theoriticalValue;
       

		printf("No of threads %i\n" ,no_threads);
		printf("GFLOPS %f\n" , gflops);
		printf("Time taken %f\n",total_time);
		printf("Efficiency %f\n", efficiency);
		free(threads);
}

void * matrixmsp(void * arg)
{
		int no_threads = (int)(intptr_t)arg;
		int incremental_step = step++; 
 
		for (int i = incremental_step * SIZE / no_threads; i < (incremental_step + 1) * SIZE / no_threads; i++) 
			for (int j = 0; j < SIZE; j++) 
				for (int k = 0; k < SIZE; k++) 
					matC[i][j] += matA[i][k] * transB[j][k]; 
		
}