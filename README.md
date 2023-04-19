# MulticoreProj
## How to Build

On CIMS, run 

    module load gcc-12.2
    
    gcc -O3 -std=c99 -Wall -fopenmp -o forecasting forecasting.c -lm

## How to Run

    ./forecasting thread_number

## Result
The output file will be in the same execute path, 
 "output" has the predicted data with all combination of alpha and beta, optimal_output has the optimal predicted data
