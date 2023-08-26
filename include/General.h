#pragma once

#ifndef General_H
#define General_H

#include <bits/stdc++.h>

using Int  = const int;
using Double  = const double;
using LLint  = long long int;

const double SQRT2   = 1.4142135623730950;
const double SQRT3   = 1.7320508075688772;
const double PI      = 3.1415926535897932;
const double EPSILON = 1.000000000000E-12;

template <class Type> static Type Absolute(Type m);
template <class Type> static void Swap(Type &a, Type &b);

template <class Type> static void CheckNaN(Type);
template <class Type> static void CheckNaN(Type*, Int, Int = 0);
template <class Type> static void CheckNaN(Type**, Int, Int, Int = 0, Int = 0);
template <class Type> static void CheckNaN(Type***, Int, Int, Int, Int = 0, Int = 0, Int = 0);
template <class Type> static void CheckNaN(Type****, Int, Int, Int, Int, Int = 0, Int = 0, Int = 0, Int = 0);
template <class Type> static void CheckNaN(Type*****, Int, Int, Int, Int, Int, Int = 0, Int = 0, Int = 0, Int = 0, Int = 0);

template <class Type> static void CheckINF(Type);
template <class Type> static void CheckINF(Type*, Int, Int = 0);
template <class Type> static void CheckINF(Type**, Int, Int, Int = 0, Int = 0);
template <class Type> static void CheckINF(Type***, Int, Int, Int, Int = 0, Int = 0, Int = 0);
template <class Type> static void CheckINF(Type****, Int, Int, Int, Int, Int = 0, Int = 0, Int = 0, Int = 0);
template <class Type> static void CheckINF(Type*****, Int, Int, Int, Int, Int, Int = 0, Int = 0, Int = 0, Int = 0, Int = 0);

template <class Type> static void Allocate(Type*&, Int, Int = 0);
template <class Type> static void Allocate(Type**&, Int, Int, Int = 0, Int = 0);
template <class Type> static void Allocate(Type***&, Int, Int, Int, Int = 0, Int = 0, Int = 0);
template <class Type> static void Allocate(Type****&, Int, Int, Int, Int, Int = 0, Int = 0, Int = 0, Int = 0);
template <class Type> static void Allocate(Type*****&, Int, Int, Int, Int, Int, Int = 0, Int = 0, Int = 0, Int = 0, Int = 0);

template <class Type> static void Deallocate(Type*, Int, Int = 0);
template <class Type> static void Deallocate(Type**, Int, Int, Int = 0, Int = 0);
template <class Type> static void Deallocate(Type***, Int, Int, Int, Int = 0, Int = 0, Int = 0);
template <class Type> static void Deallocate(Type****, Int, Int, Int, Int, Int = 0, Int = 0, Int = 0, Int = 0);
template <class Type> static void Deallocate(Type*****, Int, Int, Int, Int, Int, Int = 0, Int = 0, Int = 0, Int = 0, Int = 0);

template <class Type> static Type Median(Type a0, Type a1, Type a2);

template <class Type> static Type Maximum(Type a0, Type a1);
template <class Type> static Type Maximum(Type a0, Type a1, Type a2);
template <class Type> static Type Maximum(Type a0, Type a1, Type a2, Type a3);
template <class Type> static Type Maximum(Type a0, Type a1, Type a2, Type a3, Type a4);
template <class Type> static Type Maximum(Type a0, Type a1, Type a2, Type a3, Type a4, Type a5);
template <class Type> static Type Maximum(Type a0, Type a1, Type a2, Type a3, Type a4, Type a5, Type a6);
template <class Type> static Type Maximum(Type a0, Type a1, Type a2, Type a3, Type a4, Type a5, Type a6, Type a7);
template <class Type> static Type Maximum(Type a0, Type a1, Type a2, Type a3, Type a4, Type a5, Type a6, Type a7, Type a8);
template <class Type> static Type Maximum(Type a0, Type a1, Type a2, Type a3, Type a4, Type a5, Type a6, Type a7, Type a8, Type a9);

template <class Type> static Type Minimum(Type a0, Type a1);
template <class Type> static Type Minimum(Type a0, Type a1, Type a2);
template <class Type> static Type Minimum(Type a0, Type a1, Type a2, Type a3);
template <class Type> static Type Minimum(Type a0, Type a1, Type a2, Type a3, Type a4);
template <class Type> static Type Minimum(Type a0, Type a1, Type a2, Type a3, Type a4, Type a5);
template <class Type> static Type Minimum(Type a0, Type a1, Type a2, Type a3, Type a4, Type a5, Type a6);
template <class Type> static Type Minimum(Type a0, Type a1, Type a2, Type a3, Type a4, Type a5, Type a6, Type a7);
template <class Type> static Type Minimum(Type a0, Type a1, Type a2, Type a3, Type a4, Type a5, Type a6, Type a7, Type a8);
template <class Type> static Type Minimum(Type a0, Type a1, Type a2, Type a3, Type a4, Type a5, Type a6, Type a7, Type a8, Type a9);

template <class Type> static Type MinMod(Type a0, Type a1, Type a2);
template <class Type> static Type MinMod(Type a0, Type a1, Type a2, Type a3);
template <class Type> static Type MinMod(Type a0, Type a1, Type a2, Type a3, Type a4);
template <class Type> static Type MinMod(Type a0, Type a1, Type a2, Type a3, Type a4, Type a5);
template <class Type> static Type MinMod(Type a0, Type a1, Type a2, Type a3, Type a4, Type a5, Type a6);
template <class Type> static Type MinMod(Type a0, Type a1, Type a2, Type a3, Type a4, Type a5, Type a6, Type a7);
template <class Type> static Type MinMod(Type a0, Type a1, Type a2, Type a3, Type a4, Type a5, Type a6, Type a7, Type a8);
template <class Type> static Type MinMod(Type a0, Type a1, Type a2, Type a3, Type a4, Type a5, Type a6, Type a7, Type a8, Type a9);

template <class Type> static Type Sign(Type a);
template <class Type> static Type Sign(Type a, Type b);

template <class Type> static void PrintMatrix(Type **A, Int d0, Int d1, int = 5);

void ErrorMessage(const char* message);

double Erf(double x);

int RandomNumber(Int x0, Int x1);
double RandomNumber();

double NormallyDistributedRandomNumber(Double = 0.0, Double = 1.0);

double SquareRootSquaredSum(double a, double b);

void SquareRoot(double &a, double &b);

void roots(double &x0, double &x1, double a, double b, double c, bool &real);
void roots(double *x, double a, double b, double c, double d, bool &real);

int PeriodicIndex(int i, int N);

double PeriodicModulus(double x, double Lx);

double PeriodicDistance(double xi, double xo, double Lx);
double PeriodicDistance(double xi, double yi, double xo, double yo, double Lx, double Ly);
double PeriodicDistance(double xi, double yi, double zi, double xo, double yo, double zo, double Lx, double Ly, double Lz);

double PeriodicDisplacement(double xf, double xi, double Lx);

double ContinuedFraction(LLint *A, int i);

void DecimalToFraction(double Number, LLint &Numerator, LLint &Denominator, double = 1.0E-10);
void DecimalToFraction ( double Number, LLint &FullNumber, LLint &Numerator, LLint &Denominator, double = 1.0E-10);

void primeFactors(Int N);

//________________________________________________________________________________________________
// This function returns absolute value
//________________________________________________________________________________________________
template <class Type> static Type Absolute ( Type m )
{
    return (m >= 0.0 ? m : -m);
}

//________________________________________________________________________________________________
// This function swaps two numbers
//________________________________________________________________________________________________
template <class Type> static void Swap ( Type &a,
                                         Type &b )
{
    a = a + b;
    b = a - b;
    a = a - b;
}

//_______________________________________________________________________________
// Checks for NaN
//_______________________________________________________________________________
template <class Type> static void CheckNaN ( Type A )
{
    if (A != A)
        ErrorMessage("NaN encountered!");
}

template <class Type> static void CheckNaN ( Type *A,
                                             Int d1,
                                             Int offset1 )
{
    for (int i = -offset1; i < (d1-offset1); ++i)
        if (A[i] != A[i])
            ErrorMessage("NaN encountered!");
}

template <class Type> static void CheckNaN ( Type **A,
                                             Int d1,
                                             Int d2,
                                             Int offset1,
                                             Int offset2 )
{
    for (int i = -offset1; i < (d1-offset1); ++i)
        for (int j = -offset2; j < (d2-offset2); ++j)
            if (A[i][j] != A[i][j])
                ErrorMessage("NaN encountered!");
}

template <class Type> static void CheckNaN  ( Type ***A,
                                              Int d1,
                                              Int d2,
                                              Int d3,
                                              Int offset1,
                                              Int offset2,
                                              Int offset3 )
{
    for (int i = -offset1; i < (d1-offset1); ++i)
        for (int j = -offset2; j < (d2-offset2); ++j)
            for (int k = -offset3; k < (d3-offset3); ++k)
                if (A[i][j][k] != A[i][j][k])
                    ErrorMessage("NaN encountered!");
}

template <class Type> static void CheckNaN ( Type ****A,
                                             Int d1,
                                             Int d2,
                                             Int d3,
                                             Int d4,
                                             Int offset1,
                                             Int offset2,
                                             Int offset3,
                                             Int offset4 )
{
    for (int i = -offset1; i < (d1-offset1); ++i)
        for (int j = -offset2; j < (d2-offset2); ++j)
            for (int k = -offset3; k < (d3-offset3); ++k)
                for (int l = -offset4; l < (d4-offset4); ++l)
                    if (A[i][j][k][l] != A[i][j][k][l])
                        ErrorMessage("NaN encountered!");
}

template <class Type> static void CheckNaN ( Type *****A,
                                             Int d1,
                                             Int d2,
                                             Int d3,
                                             Int d4,
                                             Int d5,
                                             Int offset1,
                                             Int offset2,
                                             Int offset3,
                                             Int offset4,
                                             Int offset5 )
{
    for (int i = -offset1; i < (d1-offset1); ++i)
        for (int j = -offset2; j < (d2-offset2); ++j)
            for (int k = -offset3; k < (d3-offset3); ++k)
                for (int l = -offset4; l < (d4-offset4); ++l)
                    for (int p = -offset5; p < (d5-offset5); ++p)
                        if (A[i][j][k][l][p] != A[i][j][k][l][p])
                            ErrorMessage("NaN encountered!");
}

//_______________________________________________________________________________
// Checks for INF
//_______________________________________________________________________________
template <class Type> static void CheckINF ( Type A )
{
    if ( isinf(A) )
        ErrorMessage("INF encountered!");
}

template <class Type> static void CheckINF ( Type *A,
                                             Int d1,
                                             Int offset1 )
{
    for (int i = -offset1; i < (d1-offset1); ++i)
        if ( isinf(A[i]) )
            ErrorMessage("INF encountered!");
}

template <class Type> static void CheckINF ( Type **A,
                                             Int d1,
                                             Int d2,
                                             Int offset1,
                                             Int offset2 )
{
    for (int i = -offset1; i < (d1-offset1); ++i)
        for (int j = -offset2; j < (d2-offset2); ++j)
            if ( isinf(A[i][j]) )
                ErrorMessage("INF encountered!");
}

template <class Type> static void CheckINF ( Type ***A,
                                             Int d1,
                                             Int d2,
                                             Int d3,
                                             Int offset1,
                                             Int offset2,
                                             Int offset3 )
{
    for (int i = -offset1; i < (d1-offset1); ++i)
        for (int j = -offset2; j < (d2-offset2); ++j)
            for (int k = -offset3; k < (d3-offset3); ++k)
                if ( isinf(A[i][j][k]) )
                    ErrorMessage("INF encountered!");
}

template <class Type> static void CheckINF ( Type ****A,
                                             Int d1,
                                             Int d2,
                                             Int d3,
                                             Int d4,
                                             Int offset1,
                                             Int offset2,
                                             Int offset3,
                                             Int offset4 )
{
    for (int i = -offset1; i < (d1-offset1); ++i)
        for (int j = -offset2; j < (d2-offset2); ++j)
            for (int k = -offset3; k < (d3-offset3); ++k)
                for (int l = -offset4; l < (d4-offset4); ++l)
                    if ( isinf(A[i][j][k][l]) )
                        ErrorMessage("INF encountered!");
}

template <class Type> static void CheckINF ( Type *****A,
                                             Int d1,
                                             Int d2,
                                             Int d3,
                                             Int d4,
                                             Int d5,
                                             Int offset1,
                                             Int offset2,
                                             Int offset3,
                                             Int offset4,
                                             Int offset5 )
{
    for (int i = -offset1; i < (d1-offset1); ++i)
        for (int j = -offset2; j < (d2-offset2); ++j)
            for (int k = -offset3; k < (d3-offset3); ++k)
                for (int l = -offset4; l < (d4-offset4); ++l)
                    for (int p = -offset5; p < (d5-offset5); ++p)
                        if ( isinf(A[i][j][k][l][p]) )
                            ErrorMessage("INF encountered!");
}

//_______________________________________________________________________________
// Allocate in contiguous memory locations.
//_______________________________________________________________________________
template <class Type> static void Allocate ( Type *&m,
                                             Int d1,
                                             Int offset1 )
{
    m = new (std::nothrow) Type [d1];
    
    if (m == 0)
        ErrorMessage("Error: Memory can not be allocated!");
    
    for (int i{}; i < d1; ++i)
        m[i] = (Type)(0.0);
    
    m += offset1;
}

template <class Type> static void Allocate ( Type **&m,
                                             Int d1,
                                             Int d2,
                                             Int offset1,
                                             Int offset2 )
{
    m = new (std::nothrow) Type* [d1];
    
    if (m == 0)
        ErrorMessage("Error: Memory can not be allocated!");
    
    m[0] = new (std::nothrow) Type [d1*d2];
    
    if (m[0] == 0)
        ErrorMessage("Error: Memory can not be allocated!");
    
    for (int i{}; i < d1; ++i)
        m[i] = &(m[0][d2*i]);
    
    for (int i{}; i < d1; ++i)
        for (int j{}; j < d2; ++j)
            m[i][j] = (Type)(0.0);
    
    m += offset1;
    
    for (int i{-offset1}; i < (d1-offset1); ++i)
        m[i] += offset2;
}

template <class Type> static void Allocate ( Type ***&m,
                                             Int d1,
                                             Int d2,
                                             Int d3,
                                             Int offset1,
                                             Int offset2,
                                             Int offset3 )
{
    m = new (std::nothrow) Type** [d1];
    
    if (m == 0)
        ErrorMessage("Error: Memory can not be allocated!");
    
    for (int i{}; i < d1; ++i)
    {
        m[i] = new (std::nothrow) Type* [d2];
        
        if (m[i] == 0)
            ErrorMessage("Error: Memory can not be allocated!");
    }
    
    m[0][0] = new (std::nothrow) Type [d1*d2*d3];
    
    if (m[0][0] == 0)
        ErrorMessage("Error: Memory can not be allocated!");
    
    for (int i{}; i < d1; ++i)
        for (int j{}; j < d2; ++j)
            m[i][j] = &(m[0][0][i*d2*d3+j*d3]);
    
    for (int i{}; i < d1; ++i)
        for (int j{}; j < d2; ++j)
            for (int k{}; k < d3; ++k)
                m[i][j][k] = (Type)(0.0);
    
    m += offset1;
    
    for (int i{-offset1}; i < (d1-offset1); ++i)
        m[i] += offset2;
    
    for (int i{-offset1}; i < (d1-offset1); ++i)
        for (int j{-offset2}; j < (d2-offset2); ++j)
            m[i][j] += offset3;
}

template <class Type> static void Allocate ( Type ****&m,
                                             Int d1,
                                             Int d2,
                                             Int d3,
                                             Int d4,
                                             Int offset1,
                                             Int offset2,
                                             Int offset3,
                                             Int offset4 )
{
    m = new (std::nothrow) Type*** [d1];
    
    if (m == 0)
        ErrorMessage("Error: Memory can not be allocated!");
    
    for (int i{}; i < d1; ++i)
    {
        m[i] = new (std::nothrow) Type** [d2];
        
        if (m[i] == 0)
            ErrorMessage("Error: Memory can not be allocated!");

        for (int j{}; j < d2; ++j)
        {
            m[i][j] = new (std::nothrow) Type* [d3];
            
            if (m[i][j] == 0)
                ErrorMessage("Error: Memory can not be allocated!");
        }
    }
    
    m[0][0][0] = new (std::nothrow) Type [d1*d2*d3*d4];
    
    if (m[0][0][0] == 0)
        ErrorMessage("Error: Memory can not be allocated!");
    
    for (int i{}; i < d1; ++i)
        for (int j{}; j < d2; ++j)
            for (int k{}; k < d3; ++k)
                m[i][j][k] = &(m[0][0][0][i*d2*d3*d4 + j*d3*d4 + k*d4]);
    
    for (int i{}; i < d1; ++i)
        for (int j{}; j < d2; ++j)
            for (int k{}; k < d3; ++k)
                for (int l{}; l < d4; ++l)
                    m[i][j][k][l] = (Type)(0.0);
    
    m += offset1;
    
    for (int i{-offset1}; i < (d1-offset1); ++i)
        m[i] += offset2;
    
    for (int i{-offset1}; i < (d1-offset1); ++i)
        for (int j = -offset2; j < (d2-offset2); ++j)
            m[i][j] += offset3;
    
    for (int i{-offset1}; i < (d1-offset1); ++i)
        for (int j{-offset2}; j < (d2-offset2); ++j)
            for (int k{-offset3}; k < (d3-offset3); ++k)
                m[i][j][k] += offset4;
}

template <class Type> static void Allocate ( Type *****&m,
                                             Int d1,
                                             Int d2,
                                             Int d3,
                                             Int d4,
                                             Int d5,
                                             Int offset1,
                                             Int offset2,
                                             Int offset3,
                                             Int offset4,
                                             Int offset5 )
{
    m = new (std::nothrow) Type**** [d1];
    
    if (m == 0)
        ErrorMessage("Error: Memory can not be allocated!");
    
    for (int i{}; i < d1; ++i)
    {
        m[i] = new (std::nothrow) Type*** [d2];
        
        if (m[i]==0)
            ErrorMessage("Error: Memory can not be allocated!");

        for (int j{}; j < d2; ++j)
        {
            m[i][j] = new (std::nothrow) Type** [d3];
            
            if (m[i][j] == 0)
                ErrorMessage("Error: Memory can not be allocated!");
            
            for (int k{}; k < d3; ++k)
            {
                m[i][j][k] = new (std::nothrow) Type* [d4];
                
                if (m[i][j][k] == 0)
                    ErrorMessage("Error: Memory can not be allocated!");
            }
        }
    }
    
    m[0][0][0][0] = new (std::nothrow) Type [d1*d2*d3*d4*d5];
    
    if (m[0][0][0][0] == 0)
        ErrorMessage("Error: Memory can not be allocated!");
    
    for (int i{}; i < d1; ++i)
        for (int j{}; j < d2; ++j)
            for (int k{}; k < d3; ++k)
                for (int l{}; l < d4; ++l)
                    m[i][j][k][l] = &(m[0][0][0][0][i*d2*d3*d4*d5 + j*d3*d4*d5 + k*d4*d5 + l*d5]);
    
    for (int i{}; i < d1; ++i)
        for (int j{}; j < d2; ++j)
            for (int k{}; k < d3; ++k)
                for (int l{}; l < d4; ++l)
                    for (int p{}; p < d5; ++p)
                        m[i][j][k][l][p] = (Type)(0.0);
    
    m += offset1;
    
    for (int i{-offset1}; i < (d1-offset1); ++i)
        m[i] += offset2;
    
    for (int i{-offset1}; i < (d1-offset1); ++i)
        for (int j{-offset2}; j < (d2-offset2); ++j)
            m[i][j] += offset3;
    
    for (int i{-offset1}; i < (d1-offset1); ++i)
        for (int j{-offset2}; j < (d2-offset2); ++j)
            for (int k{-offset3}; k < (d3-offset3); ++k)
                m[i][j][k] += offset4;
    
    for (int i{-offset1}; i < (d1-offset1); ++i)
        for (int j{-offset2}; j < (d2-offset2); ++j)
            for (int k{-offset3}; k < (d3-offset3); ++k)
                for (int l{-offset4}; l < (d4-offset4); ++l)
                    m[i][j][k][l] += offset5;
}

//_______________________________________________________________________________
// Deallocate
//_______________________________________________________________________________
template <class Type> static void Deallocate ( Type *m,
                                               Int d1,
                                               Int offset1 )
{
    m -= offset1;
    
    delete [] m;
    
    m = nullptr;
}

template <class Type> static void Deallocate ( Type **m,
                                               Int d1,
                                               Int d2,
                                               Int offset1,
                                               Int offset2 )
{
    for (int i{-offset1}; i < (d1-offset1); ++i)
        m[i] -= offset2;
    
    m -= offset1;
    
    delete [] m[0];
    
    delete [] m;
    
    m = nullptr;
}

template <class Type> static void Deallocate ( Type ***m,
                                               Int d1,
                                               Int d2,
                                               Int d3,
                                               Int offset1,
                                               Int offset2,
                                               Int offset3 )
{
    for (int i{-offset1}; i < (d1-offset1); ++i)
        for (int j{-offset2}; j < (d2-offset2); ++j)
            m[i][j] -= offset3;
    
    for (int i{-offset1}; i < (d1-offset1); ++i)
        m[i] -= offset2;
    
    m -= offset1;
    
    delete [] m[0][0];
    
    for (int i{}; i < d1; ++i)
        delete [] m[i];
    
    delete [] m;
    
    m = nullptr;
}

template <class Type> static void Deallocate ( Type ****m,
                                               Int d1,
                                               Int d2,
                                               Int d3,
                                               Int d4,
                                               Int offset1,
                                               Int offset2,
                                               Int offset3,
                                               Int offset4 )
{
    for (int i{-offset1}; i < (d1-offset1); ++i)
        for (int j{-offset2}; j < (d2-offset2); ++j)
            for (int k{-offset3}; k < (d3-offset3); ++k)
                m[i][j][k] -= offset4;
    
    for (int i{-offset1}; i < (d1-offset1); ++i)
        for (int j{-offset2}; j < (d2-offset2); ++j)
            m[i][j] -= offset3;
    
    for (int i{-offset1}; i < (d1-offset1); ++i)
        m[i] -= offset2;
    
    m -= offset1;
    
    delete [] m[0][0][0];
    
    for (int i{}; i < d1; ++i)
    {
        for (int j{}; j < d2; ++j)
            delete [] m[i][j];
        
        delete [] m[i];
    }
    
    delete [] m;
    
    m = nullptr;
}

template <class Type> static void Deallocate ( Type *****m,
                                               Int d1,
                                               Int d2,
                                               Int d3,
                                               Int d4,
                                               Int d5,
                                               Int offset1,
                                               Int offset2,
                                               Int offset3,
                                               Int offset4,
                                               Int offset5 )
{
    for (int i{-offset1}; i < (d1-offset1); ++i)
        for (int j{-offset2}; j < (d2-offset2); ++j)
            for (int k{-offset3}; k < (d3-offset3); ++k)
                for (int l{-offset4}; l < (d4-offset4); ++l)
                    m[i][j][k][l] -= offset5;
    
    for (int i{-offset1}; i < (d1-offset1); ++i)
        for (int j{-offset2}; j < (d2-offset2); ++j)
            for (int k{-offset3}; k < (d3-offset3); ++k)
                m[i][j][k] -= offset4;
    
    for (int i{-offset1}; i < (d1-offset1); ++i)
        for (int j{-offset2}; j < (d2-offset2); ++j)
            m[i][j] -= offset3;
    
    for (int i{-offset1}; i < (d1-offset1); ++i)
        m[i] -= offset2;
    
    m -= offset1;
    
    delete [] m[0][0][0][0];
    
    for (int i{}; i < d1; ++i)
    {
        for (int j{}; j < d2; ++j)
        {
            for (int k{}; k < d3; ++k)
                delete [] m[i][j][k];
            
            delete [] m[i][j];
        }
        
        delete [] m[i];
    }
    
    delete [] m;
    
    m = nullptr;
}

//_______________________________________________________________________________
// Median
//_______________________________________________________________________________
template <class Type> static Type Median ( Type a0,
                                           Type a1,
                                           Type a2 )
{
    return  (Minimum(a0,a1,a2) == a0 ? Minimum(a1,a2) 
          : (Maximum(a0,a1,a2) == a0 ? Maximum(a1,a2) : a0));
}

//_______________________________________________________________________________
// Maximum
//_______________________________________________________________________________
template <class Type> static Type Maximum ( Type a0,
                                            Type a1 )
{
        return (a0 > a1 ? a0 : a1);
}

template <class Type> static Type Maximum ( Type a0,
                                            Type a1,
                                            Type a2 )
{
    return (a0 > a1 ? Maximum(a0,a2) 
                    : Maximum(a1,a2));
}

template <class Type> static Type Maximum ( Type a0,
                                            Type a1,
                                            Type a2,
                                            Type a3 )
{
    return (a0 > a1 ? Maximum(a0,a2,a3) 
                    : Maximum(a1,a2,a3));
}

template <class Type> static Type Maximum ( Type a0,
                                            Type a1,
                                            Type a2,
                                            Type a3,
                                            Type a4 )
{
    return (a0 > a1 ? Maximum(a0,a2,a3,a4) 
                    : Maximum(a1,a2,a3,a4));
}

template <class Type> static Type Maximum ( Type a0,
                                            Type a1,
                                            Type a2,
                                            Type a3,
                                            Type a4,
                                            Type a5 )
{
    return (a0 > a1 ? Maximum(a0,a2,a3,a4,a5) 
                    : Maximum(a1,a2,a3,a4,a5));
}

template <class Type> static Type Maximum ( Type a0,
                                            Type a1,
                                            Type a2,
                                            Type a3,
                                            Type a4,
                                            Type a5,
                                            Type a6 )
{
    return (a0 > a1 ? Maximum(a0,a2,a3,a4,a5,a6) 
                    : Maximum(a1,a2,a3,a4,a5,a6)); 
}

template <class Type> static Type Maximum ( Type a0,
                                            Type a1,
                                            Type a2,
                                            Type a3,
                                            Type a4,
                                            Type a5,
                                            Type a6,
                                            Type a7 )
{
    return (a0 > a1 ? Maximum(a0,a2,a3,a4,a5,a6,a7) 
                    : Maximum(a1,a2,a3,a4,a5,a6,a7));
}

template <class Type> static Type Maximum ( Type a0,
                                            Type a1,
                                            Type a2,
                                            Type a3,
                                            Type a4,
                                            Type a5,
                                            Type a6,
                                            Type a7,
                                            Type a8 )
{
    return (a0 > a1 ? Maximum(a0,a2,a3,a4,a5,a6,a7,a8) 
                    : Maximum(a1,a2,a3,a4,a5,a6,a7,a8)); 
}

template <class Type> static Type Maximum ( Type a0,
                                            Type a1,
                                            Type a2,
                                            Type a3,
                                            Type a4,
                                            Type a5,
                                            Type a6,
                                            Type a7,
                                            Type a8,
                                            Type a9 )
{
    return (a0 > a1 ? Maximum(a0,a2,a3,a4,a5,a6,a7,a8,a9) 
                    : Maximum(a1,a2,a3,a4,a5,a6,a7,a8,a9));
}

//_______________________________________________________________________________
// Minimum
//_______________________________________________________________________________
template <class Type> static Type Minimum ( Type a0,
                                            Type a1 )
{
    return (a0 < a1 ? a0 : a1);
}

template <class Type> static Type Minimum ( Type a0,
                                            Type a1,
                                            Type a2 )
{
    return (a0 < a1 ? Minimum(a0,a2) 
                    : Minimum(a1,a2));
}

template <class Type> static Type Minimum ( Type a0,
                                            Type a1,
                                            Type a2,
                                            Type a3 )
{
    return (a0 < a1 ? Minimum(a0,a2,a3) 
                    : Minimum(a1,a2,a3));
}

template <class Type> static Type Minimum ( Type a0,
                                            Type a1,
                                            Type a2,
                                            Type a3,
                                            Type a4 )
{
    return (a0 < a1 ? Minimum(a0,a2,a3,a4) 
                    : Minimum(a1,a2,a3,a4));
}

template <class Type> static Type Minimum ( Type a0,
                                            Type a1,
                                            Type a2,
                                            Type a3,
                                            Type a4,
                                            Type a5 )
{
    return (a0 < a1 ? Minimum(a0,a2,a3,a4,a5) 
                    : Minimum(a1,a2,a3,a4,a5));
}

template <class Type> static Type Minimum ( Type a0,
                                            Type a1,
                                            Type a2,
                                            Type a3,
                                            Type a4,
                                            Type a5,
                                            Type a6 )
{
    return (a0 < a1 ? Minimum(a0,a2,a3,a4,a5,a6) 
                    : Minimum(a1,a2,a3,a4,a5,a6));
}

template <class Type> static Type Minimum ( Type a0,
                                            Type a1,
                                            Type a2,
                                            Type a3,
                                            Type a4,
                                            Type a5,
                                            Type a6,
                                            Type a7 )
{
    return (a0 < a1 ? Minimum(a0,a2,a3,a4,a5,a6,a7) 
                    : Minimum(a1,a2,a3,a4,a5,a6,a7));
}

template <class Type> static Type Minimum ( Type a0,
                                            Type a1,
                                            Type a2,
                                            Type a3,
                                            Type a4,
                                            Type a5,
                                            Type a6,
                                            Type a7,
                                            Type a8 )
{
    return (a0 < a1 ? Minimum(a0,a2,a3,a4,a5,a6,a7,a8) 
                    : Minimum(a1,a2,a3,a4,a5,a6,a7,a8));
}

template <class Type> static Type Minimum ( Type a0,
                                            Type a1,
                                            Type a2,
                                            Type a3,
                                            Type a4,
                                            Type a5,
                                            Type a6,
                                            Type a7,
                                            Type a8,
                                            Type a9 )
{
    return (a0 < a1 ? Minimum(a0,a2,a3,a4,a5,a6,a7,a8,a9) 
                    : Minimum(a1,a2,a3,a4,a5,a6,a7,a8,a9));
}

//_______________________________________________________________________________
// MinMod
//_______________________________________________________________________________
template <class Type> static Type MinMod ( Type a0,
                                           Type a1 )
{
        return Median(a0,a1,(Type)(0.0));
}

template <class Type> static Type MinMod ( Type a0,
                                           Type a1,
                                           Type a2 )
{
    return (( (a0 > 0.0) 
             && (a1 > 0.0) 
             && (a2 > 0.0) ) ? Minimum(a0,a1,a2) 
             : (( (a0 < 0.0) 
                && (a1 < 0.0) 
                && (a2 < 0.0)) ? Maximum(a0,a1,a2) : 0.0));
}

template <class Type> static Type MinMod ( Type a0,
                                           Type a1,
                                           Type a2,
                                           Type a3 )
{
    return (( (a0 > 0.0) 
             && (a1 > 0.0) 
             && (a2 > 0.0) 
             && (a3 > 0.0)) ? Minimum(a0,a1,a2,a3) 
             : (( (a0 < 0.0) 
                && (a1 < 0.0) 
                && (a2 < 0.0) 
                && (a3 < 0.0)) ? Maximum(a0,a1,a2,a3) : 0.0));
}

template <class Type> static Type MinMod ( Type a0,
                                           Type a1,
                                           Type a2,
                                           Type a3,
                                           Type a4 )
{
    return (( (a0 > 0.0) 
             && (a1 > 0.0) 
             && (a2 > 0.0) 
             && (a3 > 0.0) 
             && (a4 > 0.0)) ? Minimum(a0,a1,a2,a3,a4) 
             : (( (a0 < 0.0) 
                && (a1 < 0.0) 
                && (a2 < 0.0) 
                && (a3 < 0.0) 
                && (a4 < 0.0)) ? Maximum(a0,a1,a2,a3,a4) : 0.0));
}

template <class Type> static Type MinMod ( Type a0,
                                           Type a1,
                                           Type a2,
                                           Type a3,
                                           Type a4,
                                           Type a5 )
{
    return (( (a0 > 0.0) 
             && (a1 > 0.0) 
             && (a2 > 0.0) 
             && (a3 > 0.0) 
             && (a4 > 0.0) 
             && (a5 > 0.0)) ? Minimum(a0,a1,a2,a3,a4,a5) 
             : (( (a0 < 0.0) 
                && (a1 < 0.0) 
                && (a2 < 0.0) 
                && (a3 < 0.0) 
                && (a4 < 0.0) 
                && (a5 < 0.0)) ? Maximum(a0,a1,a2,a3,a4,a5) : 0.0));
}

template <class Type> static Type MinMod ( Type a0,
                                           Type a1,
                                           Type a2,
                                           Type a3,
                                           Type a4,
                                           Type a5,
                                           Type a6 )
{
    return (( (a0 > 0.0) 
             && (a1 > 0.0) 
             && (a2 > 0.0) 
             && (a3 > 0.0) 
             && (a4 > 0.0) 
             && (a5 > 0.0) 
             && (a6 > 0.0)) ? Minimum(a0,a1,a2,a3,a4,a5,a6) 
             : (( (a0 < 0.0) 
                && (a1 < 0.0) 
                && (a2 < 0.0) 
                && (a3 < 0.0) 
                && (a4 < 0.0) 
                && (a5 < 0.0) 
                && (a6 < 0.0)) ? Maximum(a0,a1,a2,a3,a4,a5,a6) : 0.0));
}

template <class Type> static Type MinMod ( Type a0,
                                           Type a1,
                                           Type a2,
                                           Type a3,
                                           Type a4,
                                           Type a5,
                                           Type a6,
                                           Type a7 )
{
    return (( (a0 > 0.0) 
             && (a1 > 0.0) 
             && (a2 > 0.0) 
             && (a3 > 0.0) 
             && (a4 > 0.0) 
             && (a5 > 0.0) 
             && (a6 > 0.0) 
             && (a7 > 0.0)) ? Minimum(a0,a1,a2,a3,a4,a5,a6,a7) 
             : (( (a0 < 0.0) 
                && (a1 < 0.0) 
                && (a2 < 0.0) 
                && (a3 < 0.0) 
                && (a4 < 0.0) 
                && (a5 < 0.0) 
                && (a6 < 0.0) 
                && (a7 < 0.0)) ? Maximum(a0,a1,a2,a3,a4,a5,a6,a7) : 0.0));
}

template <class Type> static Type MinMod ( Type a0,
                                           Type a1,
                                           Type a2,
                                           Type a3,
                                           Type a4,
                                           Type a5,
                                           Type a6,
                                           Type a7,
                                           Type a8 )
{
    return (( (a0 > 0.0) 
             && (a1 > 0.0) 
             && (a2 > 0.0) 
             && (a3 > 0.0) 
             && (a4 > 0.0) 
             && (a5 > 0.0) 
             && (a6 > 0.0) 
             && (a7 > 0.0) 
             && (a8 > 0.0)) ? Minimum(a0,a1,a2,a3,a4,a5,a6,a7,a8) 
             : (( (a0 < 0.0) 
                && (a1 < 0.0) 
                && (a2 < 0.0) 
                && (a3 < 0.0) 
                && (a4 < 0.0) 
                && (a5 < 0.0) 
                && (a6 < 0.0) 
                && (a7 < 0.0) 
                && (a8 < 0.0)) ? Maximum(a0,a1,a2,a3,a4,a5,a6,a7,a8) : 0.0));
}

template <class Type> static Type MinMod ( Type a0,
                                           Type a1,
                                           Type a2,
                                           Type a3,
                                           Type a4,
                                           Type a5,
                                           Type a6,
                                           Type a7,
                                           Type a8,
                                           Type a9 )
{
    return (( (a0 > 0.0) 
             && (a1 > 0.0) 
             && (a2 > 0.0) 
             && (a3 > 0.0) 
             && (a4 > 0.0) 
             && (a5 > 0.0) 
             && (a6 > 0.0) 
             && (a7 > 0.0) 
             && (a8 > 0.0) 
             && (a9 > 0.0)) ? Minimum(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9) 
             : (( (a0 < 0.0) 
                && (a1 < 0.0) 
                && (a2 < 0.0) 
                && (a3 < 0.0) 
                && (a4 < 0.0) 
                && (a5 < 0.0) 
                && (a6 < 0.0) 
                && (a7 < 0.0) 
                && (a8 < 0.0) 
                && (a9 < 0.0)) ? Maximum(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9) : 0.0));
}

//_______________________________________________________________________________
// Sign
//_______________________________________________________________________________
template <class Type> static Type Sign ( Type a )
{
    return (a > EPSILON ? 1.0 : (a < -EPSILON  ? -1.0 : 0.0));
}

//_______________________________________________________________________________
// Returns |a| or -|a| if sign of b is + or -, respectively.
//_______________________________________________________________________________
template <class Type> static Type Sign ( Type a,
                                         Type b )
{
    return (b >= 0.0 ? Absolute(a) : -Absolute(a));
}

//_______________________________________________________________________________
// Print a 2D array of dimension d0 x d1 on terminal
//_______________________________________________________________________________
template <class Type> static void PrintMatrix ( Type **A,
                                                Int d0,
                                                Int d1,
                                                int StripLength )
{
    if ( ( d0 <= 0 ) || ( d1 <= 0 ) )
        std::cout << std::endl << "  (None)" << std::endl;
    else
    {
        int ColumnMax;
        
        //  Print the columns of the matrix, in strips of 5.
        for (int j = 0; j < d1; j += StripLength)
        {
            std::cout << std::endl << "  Col:    ";
            
            ColumnMax = Minimum(j+StripLength,d1);
            
            for (int i = j; i < ColumnMax; i++ )
                std::cout << std::setw(7) << i << "       ";
            
            std::cout << std::endl << "  Row" << std::endl << std::endl;
            
            for (int i = 0; i < d0; i++)
            {
                std::cout << std::setw(5) << i << ": ";
                
                for (int k = j; k < ColumnMax; k++ )
                    std::cout << std::setw(12) << ( Absolute(A[i][k]) < EPSILON ? (Type)(0.0) : A[i][k] ) << "  ";
                
                std::cout << std::endl;
            }
        }
    }
}
#endif
