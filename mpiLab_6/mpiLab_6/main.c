#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define BASE 10 // notation
#define N 4

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

typedef struct LONG_NUMBER longNum;
typedef unsigned long int lengthNum;

struct LONG_NUMBER {
    int *values; // array with number digits written in reverse order
    lengthNum length;
};

longNum longNumGen(int length) {
    longNum num;
    num.length = length;
    num.values = malloc(num.length * sizeof(int));

    for (lengthNum i = 0; i < num.length; i++)
        num.values[i] = rand() % 5;
    
    return num;
}

longNum schönhageStrasseFunction(longNum, longNum);

void longNumPrint(longNum num) {
    // will not output zeros before the most significant bit
    int a = 0;
    for (lengthNum i = num.length - 1; i > 0; i--) {
        if (num.values[i] && a == 0)
            a = 1;
        if (a)
            printf("%d", num.values[i]);
    }
    printf("%d", num.values[0]);
    puts("");
}

longNum longNumSum(longNum a, longNum b) {
    longNum c;
    
    if (a.length >= b.length) {
        c.length = a.length + 1;
        c.values = malloc(c.length * sizeof(int));
        
        c.values[a.length - 1] = a.values[a.length - 1];
        c.values[a.length] = 0;
        
        for (lengthNum i = 0; i < b.length; i++)
            c.values[i] = a.values[i] + b.values[i];
    } else {
        c.length = b.length + 1;
        c.values = malloc(c.length * sizeof(int));
        
        c.values[b.length - 1] = b.values[b.length - 1];
        c.values[b.length] = 0;
        
        for (lengthNum i = 0; i < a.length; i++)
            c.values[i] = a.values[i] + b.values[i];
    }
    
    return c;
}

longNum longNumSub(longNum a, longNum b) {
    if (a.length >= b.length) {
        for (lengthNum i = 0; i < b.length; i++)
            a.values[i] -= b.values[i];
        return a;
    } else {
        for (lengthNum i = 0; i < a.length; i++)
            b.values[i] -= a.values[i];
        return b;
    }
}
 
void normalization(longNum num) {
    for (lengthNum i = 0; i < num.length - 1; i++) {
        // if the number is greater than the maximum, then the transfer is performed,
        // otherwise we make a loan
        if (num.values[i] >= BASE) {
            int carryover = num.values[i] / BASE;
            num.values[i + 1] += carryover;
            num.values[i] -= carryover * BASE;
        } else if (num.values[i] < 0) {
            int carryover = (num.values[i] + 1) / BASE - 1;
            num.values[i + 1] += carryover;
            num.values[i] -= carryover * BASE;
        }
    }
}

longNum bringing(longNum num, lengthNum len) {
    longNum newNum;
    newNum.length = len;
    newNum.values = malloc(newNum.length * sizeof(int));
    
    for (int i = 0; i < num.length; i++)
        newNum.values[i] = num.values[i];
    
    for (int i = (int)num.length; i < newNum.length; i++)
        newNum.values[i] = 0;
    
    return newNum;
}

complex *fftM(complex *as, lengthNum size) {
    lengthNum n = size;
    int k = 0;
    while ((1 << k) < n) k++;
    
    int *rev = malloc(n * sizeof(int));
    rev[0] = 0;
    int high1 = -1;
    for (int i = 1; i < n; i++) {
        if ((i & (i - 1)) == 0)
            high1++;
        rev[i] = rev[i ^ (1 << high1)];
        rev[i] |= (1 << (k - high1 - 1));
    }
    
    complex *roots = malloc(n * sizeof(complex));
    for (int i = 0; i < n; i++) {
        double alpha = 2 * M_PI * i / n;
        roots[i].real = cos(alpha);
        roots[i].imaginary = sin(alpha);
    }
    
    complex *cur = malloc(n * sizeof(complex));
    for (int i = 0; i < n; i++)
        cur[i] = as[rev[i]];
    
    for (int len = 1; len < n; len <<= 1) {
        complex *ncur = malloc(n * sizeof(complex));
        lengthNum rstep = n / (len * 2);
        for (int pdest = 0; pdest < n;) {
            int p1 = pdest;
            for (int i = 0; i < len; i++) {
                complex val = complexMult(roots[i * rstep], cur[p1 + len]);
                ncur[pdest] = complexSum(cur[p1], val);
                ncur[pdest + len] = complexSub(cur[p1], val);
                pdest++;
                p1++;
            }
            pdest += len;
        }
        for (int i = 0; i < n; i++) {
            complex tmp;
            tmp = cur[i];
            cur[i] = ncur[i];
            ncur[i] = tmp;
        }
    }
    return cur;
}

complex *fftM_rev(complex *as, lengthNum n) {
    complex *res = fftM(as, n);
    
    complex tmp;
    tmp.real = n;
    tmp.imaginary = 0;
    
    for (int i = 0; i < n; i++)
        res[i] = complexDivision(res[i], tmp);
    
    return res;
}

longNum schönhageStrassenFunction(longNum a, longNum b) {
    // reduction to one rank
    if (a.length > b.length)
        b = bringing(b, a.length);
    else if (a.length < b.length)
        a = bringing(a, b.length);
    
    lengthNum size = a.length;
    
    longNum c;
    c.length = size + size;
    c.values = malloc(c.length * sizeof(int));
    for (int i = 0; i < c.length; i++)
        c.values[i] = 0;
    
    // convert a to complex
    complex *aComplex = malloc(size * sizeof(complex));
    for (int i = 0; i < size; i++)
        aComplex[i] = toComplex(a.values[i]);
    
    // convert b to complex
    complex *bComplex = malloc(size * sizeof(complex));
    for (int i = 0; i < size; i++)
        bComplex[i] = toComplex(b.values[i]);
    
    // apply ftt
    complex *aFtt = fftM(aComplex, size);
    complex *bFtt = fftM(bComplex, size);
    
    complex *cComplex = malloc(c.length * sizeof(complex));
    for (int i = 0; i < c.length; i++)
        cComplex[i] = toComplex(c.values[i]);
    
    for (int i = 0; i < size; i++)
        cComplex[i] = complexMult(aFtt[i], bFtt[i]);
    
    // apply iftt
    complex *cIftt = fftM_rev(cComplex, c.length);
    
    // convert cIftt to int
    for (int i = 0; i < c.length; i++)
        c.values[i] = toInt(cIftt[i]);
    
    normalization(c);
     
    return c;
}

int main(int argc, char *argv[]) {
    srand((int)time(NULL));

    // length of numbers
    int longNumLen = atoi(argv[1]);

    int size, rank;
    MPI_Status status;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // create group
    MPI_Group MPI_GROUP_ORIGINAL, MPI_GROUP_SCHÖNHAGE_STRASSEN;
    MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_ORIGINAL);
    int inGroup[size];
    for (int i = 0; i < size; i++)
        inGroup[i] = i;
    MPI_Group_incl(MPI_GROUP_ORIGINAL, size, inGroup, &MPI_GROUP_SCHÖNHAGE_STRASSEN);

    // create communicator
    MPI_Comm MPI_COMM_SCHÖNHAGE_STRASSEN;
    
    int index[size];
    for (int i = 0; i < size; i++)
        index[i] += i + 1;
    
    int edges[size];
    for (int i = 0; i < size; i++)
        edges[i] = i + 1;
    
    MPI_Graph_create(MPI_COMM_WORLD, size, index, edges, 1, &MPI_COMM_SCHÖNHAGE_STRASSEN);
    
    for (int i = 0; i < size; i++) {
        if (rank == 0 && i == 0) {
            // create special datatype for send
            MPI_Datatype TYPE0;
            MPI_Type_contiguous(longNumLen, MPI_INT, &TYPE0);
            MPI_Type_commit(&TYPE0);

            // generate long number
            longNum m_0 = longNumGen(longNumLen);
            printf("Process %d:\n m_0: ", rank);
            longNumPrint(m_0);

            MPI_Send(&m_0.length, 1, MPI_UNSIGNED_LONG, 1, 1, MPI_COMM_SCHÖNHAGE_STRASSEN);
            MPI_Send(m_0.values, 1, TYPE0, 1, 0, MPI_COMM_SCHÖNHAGE_STRASSEN);

            // special datatype invalidation
            MPI_Type_free(&TYPE0);
        } else if (rank > 0 && rank == i) {
            // first we take the number of rank
            lengthNum num = 0;
            MPI_Recv(&num, 1, MPI_UNSIGNED_LONG, i - 1, 1, MPI_COMM_SCHÖNHAGE_STRASSEN, &status);
            longNum recv_m_i = longNumGen((int)num);
            recv_m_i.length = num;

            // create special datatype for send
            MPI_Datatype TYPE_RECV_I;
            MPI_Type_contiguous((int)recv_m_i.length, MPI_INT, &TYPE_RECV_I);
            MPI_Type_commit(&TYPE_RECV_I);

            MPI_Recv(recv_m_i.values, 1, TYPE_RECV_I, i - 1, 0, MPI_COMM_SCHÖNHAGE_STRASSEN, &status);

            // special datatype invalidation
            MPI_Type_free(&TYPE_RECV_I);

            printf("\nProcess %d:\n recv_m_%d: ", rank, status.MPI_SOURCE);
            longNumPrint(recv_m_i);

            // generate long number
            longNum m_i = longNumGen(longNumLen);
            printf(" m_%d: ", rank);
            longNumPrint(m_i);

            longNum res_m = schönhageStrasseFunction(m_i, recv_m_i);
            printf(" m_%d x m_%d: ", status.MPI_SOURCE, rank);
            longNumPrint(res_m);

            if (rank != (size - 1)) {
                // create special datatype for send
                MPI_Datatype TYPE_I;
                MPI_Type_contiguous((int)res_m.length, MPI_INT, &TYPE_I);
                MPI_Type_commit(&TYPE_I);

                MPI_Send(&res_m.length, 1, MPI_UNSIGNED_LONG, i + 1, 1, MPI_COMM_SCHÖNHAGE_STRASSEN);
                MPI_Send(res_m.values, 1, TYPE_I, i + 1, 0, MPI_COMM_SCHÖNHAGE_STRASSEN);

                // special datatype invalidation
                MPI_Type_free(&TYPE_I);
            }
        }
    }

    // group and communicator invalidation
    MPI_Group_free(&MPI_GROUP_SCHÖNHAGE_STRASSEN);
    MPI_Comm_free(&MPI_COMM_SCHÖNHAGE_STRASSEN);

    MPI_Finalize();
    
    return 0;
}








longNum schönhageStrasseFunction(longNum a, longNum b) {
    // reduction to one rank
    if (a.length > b.length)
        b = bringing(b, a.length);
    else if (a.length < b.length)
        a = bringing(a, b.length);
    
    longNum c;
    c.length = a.length + b.length;
    c.values = malloc(c.length * sizeof(int));
    
    if (a.length < N) {
        memset(c.values, 0, c.length * sizeof(int));
        
        for (lengthNum i = 0; i < a.length; i++)
            for (lengthNum j = 0; j < b.length; j++)
                c.values[i + j] += a.values[i] * b.values[j];
    } else {
        // least significant part of a
        longNum aPart1;
        aPart1.values = a.values;
        aPart1.length = (a.length + 1) / 2;
        
        // the leading part of the number a
        longNum aPart2;
        aPart2.values = a.values + aPart1.length;
        aPart2.length = a.length / 2;
        
        // least significant part of b
        longNum bPart1;
        bPart1.values = b.values;
        bPart1.length = (b.length + 1) / 2;
        
        // the leading part of the number b
        longNum bPart2;
        bPart2.values = b.values + bPart1.length;
        bPart2.length = b.length / 2;
        
        // the sum of the parts of a
        longNum aPartsSum = longNumSum(aPart1, aPart2);
        normalization(aPartsSum);
         
        // the sum of the parts of b
        longNum bPartsSum = longNumSum(bPart1, bPart2);
        normalization(bPartsSum);
         
        // product of sums of parts
        longNum productSumsParts = schönhageStrasseFunction(aPartsSum, bPartsSum);
         
        // significant part
        longNum productSignPart = schönhageStrasseFunction(aPart1, bPart1);
         
        // the leading part
        longNum productLeadPart = schönhageStrasseFunction(aPart2, bPart2);
         
        // finding the sum of the middle terms
        longNum sumMidTerms = longNumSub(longNumSub(productSumsParts, productSignPart), productLeadPart);
         
        // summation of a polynomial
        memcpy(c.values, productSignPart.values, productSignPart.length * sizeof(int));
        memcpy(c.values + productSignPart.length, productLeadPart.values, productLeadPart.length * sizeof(int));
         
        for (lengthNum i = 0; i < sumMidTerms.length; i++)
            c.values[aPart1.length + i] += sumMidTerms.values[i];
    }
    normalization(c);
     
    return c;
}
