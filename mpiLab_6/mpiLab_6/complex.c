#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct COMPLEX complex;

struct COMPLEX {
    double real;
    double imaginary;
};

complex complexGen(int i) {
    srand(i + (int)time(NULL));
    
    complex num;
    num.real = rand() % 11 - 5;
    num.imaginary = rand() % 11 - 5;
    
    return num;
}

complex complexSum(complex a, complex b) {
    complex c;
    c.real = a.real + b.real;
    c.imaginary = a.imaginary + b.imaginary;
    return c;
}

complex complexSub(complex a, complex b) {
    complex c;
    c.real = a.real - b.real;
    c.imaginary = a.imaginary - b.imaginary;
    return c;
}

complex complexMult(complex a, complex b) {
    complex c;
    c.real = (a.real * b.real - a.imaginary * b.imaginary);
    c.imaginary = (a.real * b.imaginary + a.imaginary * b.real);
    return c;
}

complex complexDivision(complex a, complex b) {
    complex c;
    c.real = (a.real * b.real + a.imaginary * b.imaginary) / (pow(b.real, 2) + pow(b.imaginary, 2));
    c.imaginary = (b.real * a.imaginary + a.real * b.imaginary) / (pow(b.real, 2) + pow(b.imaginary, 2));
    return c;
}

void complexPrint(complex num) {
    if (num.imaginary >= 0)
        printf("%.1f + %.1fi\t", num.real, num.imaginary);
    else
        printf("%.1f - %.1fi\t", num.real, fabs(num.imaginary));
}

complex toComplex(int intNum) {
    complex complexNum;
    complexNum.real = intNum;
    complexNum.imaginary = 0;
    return complexNum;
}

int toInt(complex complexNum) {
    return complexNum.real;
}
