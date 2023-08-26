#include "General.h"

//_______________________________________________________________________________
// Print an error message and exit the program
//_______________________________________________________________________________
void ErrorMessage ( const char* message )
{
    std::cout << message << std::endl;
    
    exit(1);
}

//_______________________________________________________________________________
// Calculates error function
//_______________________________________________________________________________
double Erf ( double x )
{
    double t = 1.0/(1.0 + 0.5*Absolute(x));
    double tau = t*exp ( - x*x 
                         - 1.26551223 
                         + 1.00002368*t 
                         + 0.37409196*t*t 
                         + 0.09678418*t*t*t 
                         - 0.18628806*t*t*t*t 
                         + 0.27886807*t*t*t*t*t 
                         - 1.13520398*t*t*t*t*t*t 
                         + 1.48851587*t*t*t*t*t*t*t 
                         - 0.82215223*t*t*t*t*t*t*t*t 
                         + 0.17087277*t*t*t*t*t*t*t*t*t );
    
    return (x >= 0.0 ? (1.0-tau) : (tau-1.0));
}

//_______________________________________________________________________________
// This program generates a random integer number between x0 and x1
//_______________________________________________________________________________
int RandomNumber ( Int x0,
                   Int x1 )
{
    return x0 + rand()%((x1-x0)+1);
}

//_______________________________________________________________________________
// This program generates a random double between 0.0 and 1.0
//_______________________________________________________________________________
double RandomNumber ( )
{
    return (double)rand()/((double)(RAND_MAX));
}

//_______________________________________________________________________________
// This program generates a normally distributed random double with a specified 
// mean and a standard deviation.
//_______________________________________________________________________________
double NormallyDistributedRandomNumber ( Double Mean,
                                         Double SD )
{
    double U1, U2, W, multiplier;
    static double X1, X2;
    static int call = 0;
    
    if (call == 1)
    {
        call = !call;
        return (Mean + SD * (double) X2);
    }
    
    do
    {
        U1 = -1.0 + ((double) rand() / RAND_MAX) * 2.0;
        U2 = -1.0 + ((double) rand() / RAND_MAX) * 2.0;
        W = pow (U1,2) + pow (U2,2);
    } while (W >= 1 || W == 0);
    
    multiplier = sqrt((-2.0 * log (W)) / W);
    X1 = U1 * multiplier;
    X2 = U2 * multiplier;
    
    call = !call;
    
    return (Mean + SD * (double) X1);
}

//_______________________________________________________________________________
// Calculation of square root of (a*a + b*b)
//_______________________________________________________________________________
double SquareRootSquaredSum ( double a,
                              double b )
{
    return (Absolute(a) > Absolute(b) ? Absolute(a)*sqrt(1.0+b*b/(a*a)) : (Absolute(b) == 0.0 ? 0.0 : Absolute(b)*sqrt(1.0+a*a/(b*b))));
}

//_______________________________________________________________________________
// Calculation of square root of (a + i*b) where i = sqrt(-1.0)
// Real and complex part will be returned in a and b respectively.
//_______________________________________________________________________________
void SquareRoot ( double &a,
                  double &b )
{
    a = sqrt(2.0*(SquareRootSquaredSum(a,b)-a))*0.5;
    b /= (2.0*a);
    Swap(a,b);
}

//_______________________________________________________________________________
// This function returns the roots of a quadratic equation of the form 
// ax^2 + bx + c = 0. Here 'a', 'b' and 'c' are assumed to be real. 
// Roots are returned in 'x0' and 'x1'. If 'real' is true then both 
// of those roots are real and those are 'x0' and 'x1'. Otherwise roots 
// are (x0 + i x1) and (x0 - i x1) where i = sqrt(-1).
//_______________________________________________________________________________
void roots ( double &x0,
             double &x1,
             double a,
             double b,
             double c,
             bool &real )
{
    double d = (b*b - 4.0*a*c), denominator = 1.0/(2.0*a);
    
    if (d >= 0)
    {
        real = true;
        d = sqrt(d);
        x0 = (-b + d)*denominator;
        x1 = (-b - d)*denominator;
    }
    else
    {
        real = false;
        d = sqrt(-d);
        x0 = -b*denominator;
        x1 = d*denominator;
    }
}

//_______________________________________________________________________________
// This function returns the roots of a cubic equation of the form 
// ax^3 + bx^2 + cx + d = 0. 
// Here 'a', 'b', 'c' and 'd' are assumed to be real.
//_______________________________________________________________________________
void roots ( double *x,
             double a,
             double b,
             double c,
             double d,
             bool &real )
{
    
}

//_______________________________________________________________________________
// Periodic index
//_______________________________________________________________________________
int PeriodicIndex ( int i,
                    int N )
{
    return (i > 0 ? (i % N) : PeriodicIndex(N - (Absolute(i) % N),N));
}

//_______________________________________________________________________________
// Periodic modulus
//_______________________________________________________________________________
double PeriodicModulus ( double x,
                         double Lx )
{
    while ((x >= Lx) || (x < 0))
    {
        if (x >= Lx)    
            x -= Lx;
        else        
            x += Lx;
    }
    
    return x;
}

//_______________________________________________________________________________
// Periodic distance between two points
//_______________________________________________________________________________
double PeriodicDistance ( double xi,
                          double xo,
                          double Lx )
{
    return Minimum(Absolute(xo-xi),(Lx-Absolute(xo-xi)));
}

double PeriodicDistance ( double xi,
                          double yi,
                          double xo,
                          double yo,
                          double Lx,
                          double Ly )
{
    double lx = PeriodicDistance(xi,xo,Lx), ly = PeriodicDistance(yi,yo,Ly);
    
    return sqrt(lx*lx + ly*ly);
}

double PeriodicDistance ( double xi,
                          double yi,
                          double zi,
                          double xo,
                          double yo,
                          double zo,
                          double Lx,
                          double Ly,
                          double Lz )
{
    double lx = PeriodicDistance(xi,xo,Lx);
    double ly = PeriodicDistance(yi,yo,Ly);
    double lz = PeriodicDistance(zi,zo,Lz);
    
    return sqrt(lx*lx + ly*ly + lz*lz);
}

//_______________________________________________________________________________
// Periodic displacement
//_______________________________________________________________________________
double PeriodicDisplacement ( double xf,
                              double xi,
                              double Lx )
{
    return ( (xf-xi) >  Lx*0.5 ? (xf-xi)-Lx :
           ( (xf-xi) < -Lx*0.5 ? (xf-xi)+Lx : (xf-xi) ) );
}

//_______________________________________________________________________________
// Decimal to Fraction
//_______________________________________________________________________________
double ContinuedFraction ( LLint *A,
                           int i )
{
        long int Numerator = 1, Denominator = A[i];
        
        for (int j = i; j > 0; j--)
        {
            Numerator += Denominator*A[j-1];
            Swap(Numerator,Denominator);
        }
        
        return (double)(Numerator)/((double)(Denominator));
}

void DecimalToFraction ( double Number,
                         LLint &Numerator,
                         LLint &Denominator,
                         double Tolerance )
{
    if (Absolute(Number) < Tolerance)
    {
        Numerator = 0;
        Denominator = 1;
        std::cout << Number << " = 0/1";
    }
    else
    {
        int i = 0;
        double sign = Sign(Number);
        
        Number *= sign;
        
        // Create an array of variable length or a linked list.
        LLint *A;
        double Number0 = Number;
        
        int MemoryLength = 100;
        
        Number = 1.0/Number;
        
        Allocate(A,MemoryLength);

        bool flag = true;
        
        while (flag)
        {
            A[i] = floor(Number);
            Number -= A[i];
            
            if (Absolute(Number0 - ContinuedFraction(A,i)) < Tolerance)
                flag = false;
            else
                Number = 1.0/Number;
            i++;
            
            if (i == MemoryLength)
                ErrorMessage("Allocated memory is not enough to store the continued fraction!");
        }
        
        Numerator = 1;
        Denominator = A[i-1];
        
        for (int j = (i-1); j > 0; j--)
        {
            Numerator += Denominator*A[j-1];
            Swap(Numerator,Denominator);
        }
        
        Numerator *= (LLint)(sign);
        Deallocate(A,MemoryLength);
    }
}

void DecimalToFraction ( double Number,
                         LLint &FullNumber,
                         LLint &Numerator,
                         LLint &Denominator,
                         double Tolerance )
{
    if (Absolute(Number) > 1.0)
    {
        FullNumber = (LLint)(int(Number));
        
        Number -= (double)(FullNumber);
        Number *= Sign(Number);
    }
    else
        FullNumber = 0;
    
    if (Absolute(Number) < Tolerance)
    {
        Numerator = 0;
        Denominator = 1;
        std::cout << Number << " = 0/1";
    }
    else
    {
        int i = 0;
        double sign = Sign(Number);
        
        Number *= sign;
        
        // Create an array of variable length or a linked list.
        LLint *A;
        double Number0 = Number;
        
        int MemoryLength = 100;
        
        Number = 1.0/Number;
        
        Allocate(A,MemoryLength);
        
        bool flag = true;
        
        while (flag)
        {
            A[i] = floor(Number);
            Number -= A[i];
            
            if (Absolute(Number0 - ContinuedFraction(A,i)) < Tolerance)
                flag = false;
            else
                Number = 1.0/Number;
            i++;
            
            if (i == MemoryLength)
                ErrorMessage("Allocated memory is not enough to store the continued fraction!");
        }
        
        Numerator = 1;
        Denominator = A[i-1];
        
        for (int j = (i-1); j > 0; j--)
        {
            Numerator += Denominator*A[j-1];
            Swap(Numerator,Denominator);
        }
        
        Numerator *= (LLint)(sign);
        Deallocate(A,MemoryLength);
    }
}

//_______________________________________________________________________________
// Prime factorization
//_______________________________________________________________________________
void primeFactors ( Int N )
{
    int n = N;
    
    // Print the number of 2s that divide n 
    while ( n % 2 == 0) 
    { 
        std::cout << 2;
        
        if (n/2 != 1)
            std::cout << " x ";
        
        n /= 2;
    }
    
    // n must be odd at this point. So we can skip one element (Note i += 2)
    for (int i = 3; i <= sqrt(n); i += 2) 
    {
        // While i divides n, print i and divide n 
        while (n % i == 0) 
        { 
            std::cout << i;
            
            if (n/i != 1)
                std::cout << " x ";
            
            n /= i;
        }
    } 
    
    // This condition is to handle the case when n is a prime number greater than 2 
    if (n > 2)
        std::cout << n << std::endl;
    
    std::cout << std::endl;
}
