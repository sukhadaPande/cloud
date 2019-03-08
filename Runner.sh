gcc CPUbenchmark.c -mavx -pthread -o CPUbenchmark
gcc MatrixMultiplication.c -pthread -o2 matrixmultiplication


./CPUbenchmark
./matrixmultiplication