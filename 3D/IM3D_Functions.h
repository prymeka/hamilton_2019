// this file needs to be in the same directory as the main file
// constants and functions to run the 3-dimensional model for all update orders

const int SIZE=16; // size of the lattice (total number of sites is SIZE^3)
const int MCS = 10000; // total number of states to be generated
const int BINS_SIZE = 100; // size of bins to average over in order to smooth out fluctuations
const int SEPARATION = 11; // correlation will be calc. for separation of 0 to SEPARATION-1

// struct used in ChoosePosition_Random and Order to return position
typedef struct{
    int x;
    int y;
    int z;
} position;


void InitializeSigma( int sigma[][SIZE][SIZE] );
position ChoosePosition_Random();
position ChoosePosition_Order( int c );
void Hilbert( int s, int x, int y, int z, int dx, int dy, int dz, int dx2, int dy2, int dz2, int dx3, int dy3, int dz3, int Position[][3] );
void Lebesgue( int x, int y, int z, int width, int Position[][3] );
double DeltaU( int sigma[][SIZE][SIZE], int x, int y, int z );
int TestFlip( int e, double beta );
void InitializeCorrelation( int bins_number, double correl_data[][SEPARATION] );
void Correlation( int a, int N, int sigma[][SIZE][SIZE], double correl_data[][SEPARATION] );
void Average( int bins_number, double correl_data[][SEPARATION], double avg[] );
void StandardDeviation( int bins_number, double correl_data[][SEPARATION], double avg[], double standard_deviation[] );
void Run_Random( double beta, int bins_number, double avg[SEPARATION], double standard_deviation[SEPARATION] );
void Run_Order( double beta, int bins_number, double avg[SEPARATION], double standard_deviation[SEPARATION] );
void Run_Hilbert( double beta, int bins_number, double avg[SEPARATION], double standard_deviation[SEPARATION] );
void Run_Lebesgue( double beta, int bins_number, double avg[SEPARATION], double standard_deviation[SEPARATION] );


void InitializeSigma( int sigma[][SIZE][SIZE] ){
    // initialise array holding all spins by randomly assigning +/- 1
    for ( int a=0; a<SIZE; a++ ){
        for ( int b=0; b<SIZE; b++ ){
            for ( int c=0; c<SIZE; c++ ){
                if ( (double)rand()/(double)RAND_MAX >= 0.5 ){
                    sigma[a][b][c] = 1;
                }
                else{
                    sigma[a][b][c] = -1;
                }
            }
        }
    }
}

position ChoosePosition_Random(){
    // return random values of x, y and z as struct
    position p = {rand()%SIZE, rand()%SIZE, rand()%SIZE};
    return p;
}

position ChoosePosition_Order( int c ){
    // return c-th point (value of,  y, and z as struct) going row-by-row
    int x = c%SIZE;
    int y = c;
    while ( y >= SIZE*SIZE )
    {
        y -= SIZE*SIZE;
    }
    y = y/SIZE;
    int z = c/(SIZE*SIZE);

    position p = {x, y, z};
    return p;
}

void Hilbert( int s, int x, int y, int z, int dx, int dy, int dz, int dx2, int dy2, int dz2, int dx3, int dy3, int dz3, int Position[][3] ){
    // update Position array with points according to Hilbert curve
    if( s == 1 ){
        static int count_h = -1;
        count_h++;
        if ( count_h == SIZE*SIZE*SIZE ){ count_h = 0; }

        Position[count_h][0] = x;
        Position[count_h][1] = y;
        Position[count_h][2] = z;
        return;
    }
    else{
        s/=2;
        if( dx<0 ){ x-=s*dx; }
        if( dy<0 ){ y-=s*dy; }
        if( dz<0 ){ z-=s*dz; }
        if( dx2<0 ){ x-=s*dx2; }
        if( dy2<0 ){ y-=s*dy2; }
        if( dz2<0 ){ z-=s*dz2; }
        if( dx3<0 ){ x-=s*dx3; }
        if( dy3<0 ){ y-=s*dy3; }
        if( dz3<0 ){ z-=s*dz3; }
        Hilbert( s, x, y, z, dx2, dy2, dz2, dx3, dy3, dz3, dx, dy, dz, Position );
        Hilbert( s, x+s*dx, y+s*dy, z+s*dz, dx3, dy3, dz3, dx, dy, dz, dx2, dy2, dz2, Position );
        Hilbert( s, x+s*dx+s*dx2, y+s*dy+s*dy2, z+s*dz+s*dz2, dx3, dy3, dz3, dx, dy, dz, dx2, dy2, dz2, Position );
        Hilbert( s, x+s*dx2, y+s*dy2, z+s*dz2, -dx, -dy, -dz, -dx2, -dy2, -dz2, dx3, dy3, dz3, Position );
        Hilbert( s, x+s*dx2+s*dx3, y+s*dy2+s*dy3, z+s*dz2+s*dz3, -dx, -dy, -dz, -dx2, -dy2, -dz2, dx3, dy3, dz3, Position );
        Hilbert( s, x+s*dx+s*dx2+s*dx3, y+s*dy+s*dy2+s*dy3, z+s*dz+s*dz2+s*dz3, -dx3, -dy3, -dz3, dx, dy, dz, -dx2, -dy2, -dz2, Position );
        Hilbert( s, x+s*dx+s*dx3, y+s*dy+s*dy3, z+s*dz+s*dz3, -dx3, -dy3, -dz3, dx, dy, dz, -dx2, -dy2, -dz2, Position );
        Hilbert( s, x+s*dx3, y+s*dy3, z+s*dz3, dx2, dy2, dz2, -dx3, -dy3, -dz3, -dx, -dy, -dz, Position );
    }
}

void Lebesgue( int x, int y, int z, int width, int Position[][3] ){
    // update Position array with points according to Lebesgue curve
    if( width == 1 ){
        static int count_l = -1;
        count_l++;
        if ( count_l == SIZE*SIZE*SIZE ){ count_l = 0; }

        Position[count_l][0] = x;
        Position[count_l][1] = y;
        Position[count_l][2] = z;
        return;
    }

    width /= 2;
    Lebesgue( x, y, z, width, Position );
    Lebesgue( x+width, y, z, width, Position );
    Lebesgue( x, y+width, z, width, Position );
    Lebesgue( x+width, y+width, z, width, Position );
    Lebesgue( x, y, z+width, width, Position );
    Lebesgue( x+width, y, z+width, width, Position );
    Lebesgue( x, y+width, z+width, width, Position );
    Lebesgue( x+width, y+width, z+width, width, Position );
}

double DeltaU( int sigma[][SIZE][SIZE], int x, int y, int z ){
    // return energy difference after flip
    int up, down, right, left, front, back;
    double s;

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
    if ( y==SIZE-1 ){
        front = 0;
    }
    else{
        front = y+1;
    }
    if ( y==0 ){
        back = SIZE-1;
    }
    else{
        back = y-1;
    }
    if ( z==SIZE-1 ){
        up = 0;
    }
    else{
        up = z+1;
    }
    if ( z==0 ){
        down = SIZE-1;
    }
    else{
        down = z-1;
    }

    s = 2 * sigma[x][y][z] * (sigma[right][y][z] + sigma[left][y][z] + sigma[x][front][z] + sigma[x][back][z] + sigma[x][y][up] + sigma[x][y][down]);
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

void Correlation( int a, int N, int sigma[][SIZE][SIZE], double correl_data[][SEPARATION] ){
    // calculate the correletion between a site and a site 'd' away for all x, y, d<11 and save to array
    double norm = N*BINS_SIZE;

    for ( int d=0; d<SEPARATION; d++ ){
        for ( int z=0; z<SIZE; z++ ){
            for ( int y=0; y<SIZE; y++ ){
                for ( int x=0; x<SIZE; x++ ){
                    if ( x+d < SIZE ){
                        correl_data[a][d] += (double)(sigma[x][y][z]*sigma[x+d][y][z])/norm;
                    }
                    else{
                        correl_data[a][d] += (double)(sigma[x][y][z]*sigma[x+d-SIZE][y][z])/norm;
                    }
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
    int N = SIZE*SIZE*SIZE;

    int sigma[SIZE][SIZE][SIZE];
    InitializeSigma( sigma );

    double correl_data[bins_number][SEPARATION];
    InitializeCorrelation( bins_number, correl_data );

    int e, x, y, z;
    for ( int a=0; a<bins_number; a++ ){
        for ( int b=0; b<BINS_SIZE; b++ ){
            for ( int c=0; c<N; c++ ){
                position p = ChoosePosition_Random();
                x = p.x;
                y = p.y;
                z = p.z;
                e = DeltaU( sigma, x, y, z );
                if ( TestFlip( e, beta ) == 0 ){
                    sigma[x][y][z] = -sigma[x][y][z];
                }
            }
            Correlation( a, N, sigma, correl_data );
        }
    }
    Average( bins_number, correl_data, avg );
    StandardDeviation( bins_number, correl_data, avg, standard_deviation );
}

void Run_Order( double beta, int bins_number, double avg[SEPARATION], double standard_deviation[SEPARATION] ){
    // run the metropolis algorithm, find avg. and s.d. for order
    int N = SIZE*SIZE*SIZE;

    int sigma[SIZE][SIZE][SIZE];
    InitializeSigma( sigma );

    double correl_data[bins_number][SEPARATION];
    InitializeCorrelation( bins_number, correl_data );

    int e, x, y, z;
    for ( int a=0; a<bins_number; a++ ){
        for ( int b=0; b<BINS_SIZE; b++ ){
            for ( int c=0; c<N; c++ ){
                position p = ChoosePosition_Order( c );
                x = p.x;
                y = p.y;
                z = p.z;
                e = DeltaU( sigma, x, y, z );
                if ( TestFlip( e, beta ) == 0 ){
                    sigma[x][y][z] = -sigma[x][y][z];
                }
            }
            Correlation( a, N, sigma, correl_data );
        }
    }
    Average( bins_number, correl_data, avg );
    StandardDeviation( bins_number, correl_data, avg, standard_deviation );
}

void Run_Hilbert( double beta, int bins_number, double avg[SEPARATION], double standard_deviation[SEPARATION] ){
    // run the metropolis algorithm, find avg. and s.d. for Hilbert curve
    int N = SIZE*SIZE*SIZE;

    int sigma[SIZE][SIZE][SIZE];
    InitializeSigma( sigma );

    double correl_data[bins_number][SEPARATION];
    InitializeCorrelation( bins_number, correl_data );

    int Position[N][3];
    Hilbert( SIZE, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, Position );

    int e, x, y, z;
    for ( int a=0; a<bins_number; a++ ){
        for ( int b=0; b<BINS_SIZE; b++ ){
            for ( int c=0; c<N; c++ ){
                x = Position[c][0];
                y = Position[c][1];
                z = Position[c][2];
                e = DeltaU( sigma, x, y, z );
                if ( TestFlip( e, beta ) == 0 ){
                    sigma[x][y][z] = -sigma[x][y][z];
                }
            }
            Correlation( a, N, sigma, correl_data );
        }
    }
    Average( bins_number, correl_data, avg );
    StandardDeviation( bins_number, correl_data, avg, standard_deviation );
}

void Run_Lebesgue( double beta, int bins_number, double avg[SEPARATION], double standard_deviation[SEPARATION] ){
    // run the metropolis algorithm, find avg. and s.d. for Lebesgue curve
    int N = SIZE*SIZE*SIZE;

    int sigma[SIZE][SIZE][SIZE];
    InitializeSigma( sigma );

    double correl_data[bins_number][SEPARATION];
    InitializeCorrelation( bins_number, correl_data );

    int Position[N][3];
    Lebesgue( 0, 0, 0, SIZE, Position );

    int e, x, y, z;
    for ( int a=0; a<bins_number; a++ ){
        for ( int b=0; b<BINS_SIZE; b++ ){
            for ( int c=0; c<N; c++ ){
                x = Position[c][0];
                y = Position[c][1];
                z = Position[c][2];
                e = DeltaU( sigma, x, y, z );
                if ( TestFlip( e, beta ) == 0 ){
                    sigma[x][y][z] = -sigma[x][y][z];
                }
            }
            Correlation( a, N, sigma, correl_data );
        }
    }
    Average( bins_number, correl_data, avg );
    StandardDeviation( bins_number, correl_data, avg, standard_deviation );
}
