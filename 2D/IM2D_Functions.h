// this file needs to be in the same directory as the main file
// constants and functions to run the 2-dimensional model for all update paths

const int SIZE=64; // size of the lattice (total number of sites is SIZE^2)
const int MCS = 10000; // total number of states to be generated
const int BINS_SIZE = 100; // size of bins to average over in order to smooth out fluctuations
const int SEPARATION = 11; // correlation will be calc. for separation of 0 to SEPARATION-1

// struct used in ChoosePosition_Random and Order to return position
typedef struct{
    int x;
    int y;
} position;


void InitialiseSigma( int sigma[][SIZE] );
position ChoosePosition_Random();
position ChoosePosition_Order( int c );
void Hilbert( int x, int y, int width, int initial1, int initial2, int Position[][2] );
void Lebesgue( int x, int y, int width, int Position[][2] );
void Gcurve( int x, int y, int width, int Position[][2] );
double DeltaU( int sigma[][SIZE], int x, int y );
int TestFlip( int e, double beta );
void InitializeCorrelation( int bins_number, double correl_data[][SEPARATION] );
void Correlation( int a, int N, int sigma[][SIZE], double correl_data[][SEPARATION] );
void Average( int bins_number, double correl_data[][SEPARATION], double avg[] );
void StandardDeviation( int bins_number, double correl_data[][SEPARATION], double avg[], double standard_deviation[] );
void Run_Random( double beta, int bins_number, double avg[SEPARATION], double standard_deviation[SEPARATION] );
void Run_Order( double beta, int bins_number, double avg[SEPARATION], double standard_deviation[SEPARATION] );
void Run_Hilbert( double beta, int bins_number, double avg[SEPARATION], double standard_deviation[SEPARATION] );
void Run_Lebesgue( double beta, int bins_number, double avg[SEPARATION], double standard_deviation[SEPARATION] );
void Run_Gcurve( double beta, int bins_number, double avg[SEPARATION], double standard_deviation[SEPARATION] );


void InitialiseSigma( int sigma[][SIZE] ){
    // initialise array holding all spins by randomly assigning +/- 1
    for ( int a=0; a<SIZE; a++ ){
        for ( int b=0; b<SIZE; b++ ){
            if ( (double)rand()/(double)RAND_MAX >= 0.5 ){
                sigma[b][a] = 1;
            }
            else{
                sigma[b][a] = -1;
            }
        }
    }
}

position ChoosePosition_Random(){
    // return random values of x and y as struct
    position p = {rand()%SIZE, rand()%SIZE};
    return p;
}

position ChoosePosition_Order( int c ){
    // return c-th point (value of x and y as struct) going row-by-row
    int x_temp = c;
    while ( x_temp >= SIZE ){
        x_temp -= SIZE;
    }

    int y_temp = c/SIZE;

    position p = {x_temp, y_temp};
    return p;
}

void Hilbert( int x, int y, int width, int initial1, int initial2, int Position[][2] ){
    // update Position array with points according to Hilbert curve
    if ( width == 1 ){
        static int count_h = -1;
        count_h++;
        if ( count_h == SIZE*SIZE ){ count_h = 0; }

        Position[count_h][0] = x;
        Position[count_h][1] = y;
        return;
    }

    width /= 2;
    Hilbert( x+initial1*width,     y+initial1*width,     width,   initial1,   1-initial2, Position );
    Hilbert( x+initial2*width,     y+(1-initial2)*width, width,   initial1,   initial2, Position );
    Hilbert( x+(1-initial1)*width, y+(1-initial1)*width, width,   initial1,   initial2, Position );
    Hilbert( x+(1-initial2)*width, y+initial2*width,     width,   1-initial1, initial2, Position );
}

void Lebesgue( int x, int y, int width, int Position[][2] ){
    // update Position array with points according to Lebesgue curve
    if ( width == 1 ){
        static int count_l = -1;
        count_l++;
        if ( count_l == SIZE*SIZE ){ count_l = 0; }

        Position[count_l][0] = x;
        Position[count_l][1] = y;
        return;
    }

    width /= 2;
    Lebesgue( x, y, width, Position );
    Lebesgue( x+1*width, y, width, Position );
    Lebesgue( x, y+1*width, width, Position );
    Lebesgue( x+1*width, y+1*width, width, Position );
}

void Gcurve( int x, int y, int width, int Position[][2] ){
    // update Position array with points according to Gcurve curve
    if ( width == 1 ){
        static int count_g = -1;
        count_g++;
        if ( count_g == SIZE*SIZE ){ count_g = 0; }

        Position[count_g][0] = x;
        Position[count_g][1] = y;
        return;
    }

    width /= 4;
    Gcurve( x, y, width, Position );
    Gcurve( x,y+1*width, width, Position );
    Gcurve( x,y+2*width, width, Position );
    Gcurve( x,y+3*width, width, Position );
    Gcurve( x+1*width, y+3*width, width, Position );
    Gcurve( x+2*width, y+3*width, width, Position );
    Gcurve( x+3*width, y+3*width, width, Position );
    Gcurve( x+3*width, y+2*width, width, Position );
    Gcurve( x+2*width, y+2*width, width, Position );
    Gcurve( x+1*width, y+2*width, width, Position );
    Gcurve( x+1*width, y+1*width, width, Position );
    Gcurve( x+1*width, y, width, Position );
    Gcurve( x+2*width, y, width, Position );
    Gcurve( x+2*width, y+1*width, width, Position );
    Gcurve( x+3*width, y+1*width, width, Position );
    Gcurve( x+3*width, y, width, Position );
}

double DeltaU( int sigma[][SIZE], int x, int y ){
    // return energy difference after flip
    int up, down, right, left;
    double s;

    if ( y==SIZE-1 ){
        up = 0;
    }
    else{
        up = y+1;
    }
    if ( y==0 ){
        down = SIZE-1;
    }
    else{
        down = y-1;
    }
    if ( x==SIZE-1 ){
        right = 0;
    }
    else{
        right = x+1;
    }
    if ( x==0 ){
        left = SIZE-1;
    }
    else{
        left = x-1;
    }

    s = 2 * sigma[x][y] * (sigma[x][up] + sigma[x][down] + sigma[right][y] + sigma[left][y]);
    return s;
}

int TestFlip( int e, double beta ){
    // test whether the site should be flipped: yes = return 0; no = return 1
    if ( e <= 0 ){
        return 0;
    }
    else if ( (double)rand()/(double)RAND_MAX < exp( -e*beta ) ){
        return 0;
    }
    else{
        return 1;
    }
}

void InitializeCorrelation( int bins_number, double correl_data[][SEPARATION] ){
    // initialise the correletaion array to 0s
    for ( int a=0; a<bins_number; a++ ){
        for ( int b=0; b<SEPARATION; b++ ){
            correl_data[a][b] = 0;
        }
    }
}

void Correlation( int a, int N, int sigma[][SIZE], double correl_data[][SEPARATION] ){
    // calculate the correletion between a site and a site 'd' away for all x, y, d<11 and save to array
    double norm = N*BINS_SIZE;

    for ( int d=0; d<SEPARATION; d++ ){
        for ( int y=0; y<SIZE; y++ ){
            for ( int x=0; x<SIZE; x++ ){
                if ( x+d < SIZE ){
                    correl_data[a][d] += (double)(sigma[x][y]*sigma[x+d][y])/norm;
                }
                else{
                    correl_data[a][d] += (double)(sigma[x][y]*sigma[x+d-SIZE][y])/norm;
                }
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
    }
    for ( int d=0; d<SEPARATION; d++ ){
        standard_deviation[d] = sqrt(sd_sum[d]/(bins_number-1));
    }
}

void Run_Random( double beta, int bins_number, double avg[SEPARATION], double standard_deviation[SEPARATION] ){
    // run the metropolis algorithm, find avg. and s.d. for random
    int N = SIZE*SIZE;

    int sigma[SIZE][SIZE];
    InitialiseSigma( sigma );

    double correl_data[bins_number][SEPARATION];
    InitializeCorrelation( bins_number, correl_data );

    int e, x, y;
    for ( int a=0; a<bins_number; a++ ){
        for ( int b=0; b<BINS_SIZE; b++ ){
            for ( int c=0; c<N; c++ ){
                position p = ChoosePosition_Random();
                x = p.x;
                y = p.y;
                e = DeltaU( sigma, x, y );
                if ( TestFlip( e, beta ) == 0 ){
                    sigma[x][y] = -sigma[x][y];
                }
            }
            Correlation( a, N, sigma, correl_data );
        }
    }

    Average(bins_number, correl_data, avg );
    StandardDeviation( bins_number, correl_data, avg, standard_deviation );
}

void Run_Order( double beta, int bins_number, double avg[SEPARATION], double standard_deviation[SEPARATION] ){
    // run the metropolis algorithm, find avg. and s.d. for order
    int N = SIZE*SIZE;

    int sigma[SIZE][SIZE];
    InitialiseSigma( sigma );

    double correl_data[bins_number][SEPARATION];
    InitializeCorrelation( bins_number, correl_data );

    int e, x, y;
    for ( int a=0; a<bins_number; a++ ){
        for ( int b=0; b<BINS_SIZE; b++ ){
            for ( int c=0; c<N; c++ ){
                position p = ChoosePosition_Order( c );
                x = p.x;
                y = p.y;
                e = DeltaU( sigma, x, y );
                if ( TestFlip( e, beta ) == 0 ){
                    sigma[x][y] = -sigma[x][y];
                }
            }
            Correlation( a, N, sigma, correl_data );
        }
    }

    Average(bins_number, correl_data, avg );
    StandardDeviation( bins_number, correl_data, avg, standard_deviation );
}

void Run_Hilbert( double beta, int bins_number, double avg[SEPARATION], double standard_deviation[SEPARATION] ){
    // run the metropolis algorithm, find avg. and s.d. for Hilbert curve
    int N = SIZE*SIZE;

    int sigma[SIZE][SIZE];
    InitialiseSigma( sigma );

    double correl_data[bins_number][SEPARATION];
    InitializeCorrelation( bins_number, correl_data );

    int Position[N][2];
    Hilbert( 0, 0, SIZE, 0, 0, Position );

    int e, x, y;
    for ( int a=0; a<bins_number; a++ ){
        for ( int b=0; b<BINS_SIZE; b++ ){
            for ( int c=0; c<N; c++ ){
                x = Position[c][0];
                y = Position[c][1];
                e = DeltaU( sigma, x, y );
                if ( TestFlip( e, beta ) == 0 ){
                    sigma[x][y] = -sigma[x][y];
                }
            }
            Correlation( a, N, sigma, correl_data );
        }
    }

    Average(bins_number, correl_data, avg );
    StandardDeviation( bins_number, correl_data, avg, standard_deviation );
}

void Run_Lebesgue( double beta, int bins_number, double avg[SEPARATION], double standard_deviation[SEPARATION] ){
    // run the metropolis algorithm, find avg. and s.d. for Lebesgue curve
    int N = SIZE*SIZE;

    int sigma[SIZE][SIZE];
    InitialiseSigma( sigma );

    double correl_data[bins_number][SEPARATION];
    InitializeCorrelation( bins_number, correl_data );

    int Position[N][2];
    Lebesgue( 0, 0, SIZE, Position );

    int e, x, y;
    for ( int a=0; a<bins_number; a++ ){
        for ( int b=0; b<BINS_SIZE; b++ ){
            for ( int c=0; c<N; c++ ){
                x = Position[c][0];
                y = Position[c][1];
                e = DeltaU( sigma, x, y );
                if ( TestFlip( e, beta ) == 0 ){
                    sigma[x][y] = -sigma[x][y];
                }
            }
            Correlation( a, N, sigma, correl_data );
        }
    }

    Average(bins_number, correl_data, avg );
    StandardDeviation( bins_number, correl_data, avg, standard_deviation );
}

void Run_Gcurve( double beta, int bins_number, double avg[SEPARATION], double standard_deviation[SEPARATION] ){
    // run the metropolis algorithm, find avg. and s.d. for Gcurve curve
    int N = SIZE*SIZE;

    int sigma[SIZE][SIZE];
    InitialiseSigma( sigma );

    double correl_data[bins_number][SEPARATION];
    InitializeCorrelation( bins_number, correl_data );

    int Position[N][2];
    Gcurve( 0, 0, SIZE, Position );

    int e, x, y;
    for ( int a=0; a<bins_number; a++ ){
        for ( int b=0; b<BINS_SIZE; b++ ){
            for ( int c=0; c<N; c++ ){
                x = Position[c][0];
                y = Position[c][1];
                e = DeltaU( sigma, x, y );
                if ( TestFlip( e, beta ) == 0 ){
                    sigma[x][y] = -sigma[x][y];
                }
            }
            Correlation( a, N, sigma, correl_data );
        }
    }

    Average(bins_number, correl_data, avg );
    StandardDeviation( bins_number, correl_data, avg, standard_deviation );
}
