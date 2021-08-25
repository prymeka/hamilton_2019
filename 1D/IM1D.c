// the file IM1D_functions needs to be the same directory as this file
// 1-Dimensional Ising model with update paths: random, every 2nd, every 3rd, order
// running the model and outputting the data; the constants and functions are in the header file

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// constants and functions 
#include "IM1D_Functions.h" 

int main(){
    
    srand(time(NULL));
    int bins_number = ceil( (double)MCS/(double)BINS_SIZE ); // number of bins used to average over 

    double beta = 0.6; // temperature 

    // openning file 
    FILE *fptr;
    char name[FILENAME_MAX];
    sprintf(name, "Data_1D_%.2f.csv", beta);
    fptr = fopen(name, "w");
    fprintf(fptr, "beta=%.2f\n", beta);
    // columns from left: separation, avg random, sd random, avg 2nd, sd 2nd, 
    // avg 3rd, sd 3rd, avg order, sd order
    fprintf(fptr, "separation,avg_random,sd_random,avg_2nd,sd_2nd,avg_3rd,sd_3rd,avg_order,sd_order\n");

    // looping 10 times and saving all to one file
    for( int i=0; i<10; i++ ){

        if ( i == 0 ){ printf( "The process has been started...\n" ); }
        else { printf( "\n" ); }

        // arrays that will hold all the data 
        double avg_r[SEPARATION], standard_deviation_r[SEPARATION];
        double avg_2[SEPARATION], standard_deviation_2[SEPARATION];
        double avg_3[SEPARATION], standard_deviation_3[SEPARATION];
        double avg_o[SEPARATION], standard_deviation_o[SEPARATION]; 
        
        // running the program to collect data for all update paths 
        Run_Random( beta, bins_number, avg_r, standard_deviation_r );
        printf( "Random Completed - %d/10...\n", i+1 );
        Run_2ND( beta, bins_number, avg_2, standard_deviation_2 );
        printf( "Every 2nd Completed - %d/10...\n", i+1 );
        Run_3RD( beta, bins_number, avg_3, standard_deviation_3 );
        printf( "Every 3rd Completed - %d/10...\n", i+1 );
        Run_Order( beta, bins_number, avg_o, standard_deviation_o );
        printf( "Order Completed - %d/10...\n", i+1 );
        
        // outputing data into a csv file
        for ( int d=0; d<SEPARATION; d++ ){
            fprintf(fptr,"%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf  \n", d, avg_r[d], standard_deviation_r[d], avg_2[d], standard_deviation_2[d], avg_3[d], standard_deviation_3[d], avg_o[d], standard_deviation_o[d]);
        } 

    }  
    // closing file
    fclose(fptr);

    return 0; 
}
