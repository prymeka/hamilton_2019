// the file IM3D_functions needs to be the same directory as this file
// 3-Dimensional Ising model with update paths: random, order, Hilbert curve, Lebesque curve
// running the model and outputting the data; the constants and functions are in the header file

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// constants and functions 
#include "IM3D_Functions.h" 

int main(){
    
    srand(time(NULL));
    int bins_number = ceil( (double)MCS/(double)BINS_SIZE ); // number of bins used to average over 

    double beta = 0.6; // temperature 

    // openning file 
    FILE *fptr;
    char name[FILENAME_MAX];
    sprintf(name, "Data_3D_%.2f.csv", beta);
    fptr = fopen(name, "w");
    fprintf(fptr, "beta=%.2f\n", beta);
    // columns from left: separation, avg random, sd random, avg order, sd order, avg hilbert, sd hilbert, 
    // avg lebesgue, sd lebesgue
    fprintf(fptr, "separation,avg_random,sd_random,avg_order,sd_order,avg_hilbert,sd_hilbert,avg_lebesgue,sd_lebesgue\n");

    // looping 10 times and saving all to one file
    for( int i=0; i<10; i++ ){

        if ( i == 0 ){ printf( "The process has been started...\n" ); }
        else { printf( "\n" ); }

        // arrays that will hold all the data 
        double avg_r[SEPARATION], standard_deviation_r[SEPARATION];
        double avg_o[SEPARATION], standard_deviation_o[SEPARATION];
        double avg_h[SEPARATION], standard_deviation_h[SEPARATION];
        double avg_l[SEPARATION], standard_deviation_l[SEPARATION]; 
        
        // running the program to collect data for all update paths 
        Run_Random( beta, bins_number, avg_r, standard_deviation_r );
        printf( "Random Completed - %d/10...\n", i+1 );
        Run_Order( beta, bins_number, avg_o, standard_deviation_o );
        printf( "Order Completed - %d/10...\n", i+1 );
        Run_Hilbert( beta, bins_number, avg_h, standard_deviation_h );
        printf( "Hilbert Completed - %d/10...\n", i+1 );
        Run_Lebesgue( beta, bins_number, avg_l, standard_deviation_l );
        printf( "Lebesgue Completed - %d/10...\n", i+1 );
        
        // outputing data into a csv file
        for ( int d=0; d<SEPARATION; d++ ){
            fprintf( fptr,"%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", d, avg_r[d], standard_deviation_r[d], avg_o[d], standard_deviation_o[d], avg_h[d], standard_deviation_h[d], avg_l[d], standard_deviation_l[d] );
        } 

    }  
    // closing file
    fclose(fptr);

    return 0; 
}
