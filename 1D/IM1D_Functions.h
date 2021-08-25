// this file needs to be in the same directory as the main file
// constants and functions to run the 2-dimensional model for all update orders

const int N = 1000; //total number of sites
const int MCS = 10000; // total number of states to be generated
const int BINS_SIZE = 100; // size of bins to average over in order to smooth out fluctuations
const int SEPARATION = 11; // correlation will be calc. for separation of 0 to SEPARATION-1


void InitialiseSigma( int sigma[] );
int ChoosePosition_Random();
int ChoosePosition_2ND( int c );
int ChoosePosition_3RD( int c );
int ChoosePosition_Order( int c );
int DeltaU( int sigma[], int x );
int TestFlip( int e, double beta );
void InitializeCorrelation( int bins_number, double correl_data[][SEPARATION] );
void Correlation( int a, int sigma[], double correl_data[][SEPARATION] );
void Average( int bins_number, double correl_data[][SEPARATION], double avg[] );
void StandardDeviation( int bins_number, double correl_data[][SEPARATION], double avg[], double standard_deviation[] );
int Run_Random( double beta, int bins_number, double avg[SEPARATION], double standard_deviation[SEPARATION] );
int Run_2ND( double beta, int bins_number, double avg[SEPARATION], double standard_deviation[SEPARATION] );
int Run_3RD( double beta, int bins_number, double avg[SEPARATION], double standard_deviation[SEPARATION] );
int Run_Order( double beta, int bins_number, double avg[SEPARATION], double standard_deviation[SEPARATION] );


void InitialiseSigma( int sigma[] ){
    // initialise array holding all spins by randomly assigning +/- 1
    for ( int i=0; i<N; i++ ){
        if (( double)rand()/(double)RAND_MAX <= 0.5 ){
            sigma[i] = 1;
        }
        else{
            sigma[i] = -1;
        }
    }
}

int ChoosePosition_Random(){
    // return random coordinate of a point
    return rand()%N;
}

int ChoosePosition_2ND( int c ){
    // return coordinate of the c-th point from 0 going every second point
    if ( 2*c < N ){
        return 2*c;
    }
    else if ( N%2 == 0 ){
        return 2*c-N+1;
    }
    else{
        return 2*c-N;
    }
}

int ChoosePosition_3RD( int c ){
    // return coordinate of the c-th point from 0 going every third point
    if ( N%3 == 0 ){
        if ( 3*c < N ){
            return 3*c;
        }
        else if ( 3*c-N+1 < N ){
            return 3*c-N+1;
        }
        else{
            return 3*c-2*N+2;
        }
    }
    else if ( N%3 == 1 ){
        if ( 3*c < N ){
            return 3*c;
        }
        else if ( 3*c-N-1 < N ){
            return 3*c-N-1;
        }
        else{
            return 3*c-2*N+1;
        }
    }
    else{
        if ( 3*c < N ){
            return 3*c;
        }
        else if ( 3*c-N < N ){
            return 3*c-N;
        }
        else{
            return 3*c-2*N;
        }
    }
}

int ChoosePosition_Order( int c ){
    // return coordinate of c-th point from 0
    return c;
}

int DeltaU( int sigma[], int x ){
    // return energy difference after flip
    int j;
    double s;

    if ( x != 0 && x != N-1 ){
        j = sigma[x-1] + sigma[x+1];
    }
    else if ( x == 0 ){
        j = sigma[N-1] + sigma[1];
    }
    else if ( x == N-1 ){
        j = sigma[N-2] + sigma[0];
    }

    s = 2 * sigma[x] * j;
    return s;
}

int TestFlip( int e, double beta ){
    // test whether the site should be flipped: yes = return 0; no = return 1
    if ( e <= 0 ){
        return 0;
    }
    else if ( exp(-(double)e*beta) >= 1 ){
        return 0;
    }
    else if ( (double)rand()/(double)RAND_MAX < exp(-(double)e*beta) ){
        return 0;
    }
    else{
        return 1;
    }
}

void InitializeCorrelation( int bins_number, double correl_data[][SEPARATION] ){
    // initialise the correletaion array to 0s
    for ( int i=0; i<bins_number; i++ ){
        for( int d=0; d<SEPARATION; d++ ){
            correl_data[i][d] = 0;
        }
    }
}

void Correlation( int a, int sigma[], double correl_data[][SEPARATION] ){
    // calculate the correletion between a site and a site 'd' away for all x, y, d<11 and save to array
    for ( int d=0; d<SEPARATION; d++ ){
        for ( int i=0; i<N; i++ ){
            if ( i+d < N ){
                correl_data[a][d] += sigma[i]*sigma[i+d]/(double)(N*BINS_SIZE);
            }
            else{
                correl_data[a][d] += sigma[i]*sigma[i+d-N]/(double)(N*BINS_SIZE);
            }
        }
    }
}

void Average( int bins_number, double correl_data[][SEPARATION], double avg[] ){
    // calculate the average correlation over all bins for each separation value
    for ( int d=0; d<SEPARATION; d++ ){
        avg[d] = 0;
        for ( int i=0; i<bins_number; i++ ){
            avg[d] += correl_data[i][d]/bins_number;
        }
    }
}

void StandardDeviation( int bins_number, double correl_data[][SEPARATION], double avg[], double standard_deviation[] ){
    // find standard deviation over all bins for each separation value
    double sd_sum[SEPARATION];

    for ( int d=0; d<SEPARATION; d++ ){
        sd_sum[d] = 0;
        for ( int i=0; i<bins_number; i++ ){
            sd_sum[d] += (correl_data[i][d] - avg[d])*(correl_data[i][d] - avg[d]);
        }
        standard_deviation[d] = sqrt(sd_sum[d]/(bins_number-1));
    }
}

int Run_Random( double beta, int bins_number, double avg[SEPARATION], double standard_deviation[SEPARATION] ){
    // run the metropolis algorithm, find avg. and s.d. for random
    int sigma[N];
    InitialiseSigma( sigma );

    double correl_data[bins_number][SEPARATION];
    InitializeCorrelation( bins_number, correl_data );

    int x, e;
    for ( int a=0; a<bins_number; a++ ){
        for ( int b=0; b<BINS_SIZE; b++ ){
            for ( int c=0; c<N; c++ ){
                x = ChoosePosition_Random();
                e = DeltaU( sigma, x );
                if ( TestFlip( e, beta ) == 0 ){
                    sigma[x] = -sigma[x];
                }
            }
            Correlation( a, sigma, correl_data );
        }
    }
    Average( bins_number, correl_data, avg );
    StandardDeviation( bins_number, correl_data, avg, standard_deviation );
}

int Run_2ND( double beta, int bins_number, double avg[SEPARATION], double standard_deviation[SEPARATION] ){
    // run the metropolis algorithm, find avg. and s.d. for every 2nd
    int sigma[N];
    InitialiseSigma( sigma );

    double correl_data[bins_number][SEPARATION];
    InitializeCorrelation( bins_number, correl_data );

    int x, e;
    for ( int a=0; a<bins_number; a++ ){
        for ( int b=0; b<BINS_SIZE; b++ ){
            for ( int c=0; c<N; c++ ){
                x = ChoosePosition_2ND( c );
                e = DeltaU( sigma, x );
                if ( TestFlip( e, beta ) == 0 ){
                    sigma[x] = -sigma[x];
                }
            }
            Correlation( a, sigma, correl_data );
        }
    }
    Average( bins_number, correl_data, avg );
    StandardDeviation( bins_number, correl_data, avg, standard_deviation );
}

int Run_3RD( double beta,  int bins_number, double avg[SEPARATION], double standard_deviation[SEPARATION] ){
    // run the metropolis algorithm, find avg. and s.d. for every 3rd
    int sigma[N];
    InitialiseSigma( sigma );

    double correl_data[bins_number][SEPARATION];
    InitializeCorrelation( bins_number, correl_data );

    int x, e;
    for ( int a=0; a<bins_number; a++ ){
        for ( int b=0; b<BINS_SIZE; b++ ){
            for ( int c=0; c<N; c++ ){
                x = ChoosePosition_3RD( c );
                e = DeltaU( sigma, x );
                if ( TestFlip( e, beta ) == 0 ){
                    sigma[x] = -sigma[x];
                }
            }
            Correlation( a, sigma, correl_data );
        }
    }
    Average( bins_number, correl_data, avg );
    StandardDeviation( bins_number, correl_data, avg, standard_deviation );
}

int Run_Order( double beta, int bins_number, double avg[SEPARATION], double standard_deviation[SEPARATION] ){
    // run the metropolis algorithm, find avg. and s.d. for order
    int sigma[N];
    InitialiseSigma( sigma );

    double correl_data[bins_number][SEPARATION];
    InitializeCorrelation( bins_number, correl_data );

    int x, e;
    for ( int a=0; a<bins_number; a++ ){
        for ( int b=0; b<BINS_SIZE; b++ ){
            for ( int c=0; c<N; c++ ){
                x = ChoosePosition_Order( c );
                e = DeltaU( sigma, x );
                if ( TestFlip( e, beta ) == 0 ){
                    sigma[x] = -sigma[x];
                }
            }
            Correlation( a, sigma, correl_data );
        }
    }
    Average( bins_number, correl_data, avg );
    StandardDeviation( bins_number, correl_data, avg, standard_deviation );
}
