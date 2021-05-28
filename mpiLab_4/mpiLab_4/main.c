#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#define M 3

typedef struct COMPLEX complex;
typedef struct COMPLEX_MATRIX complex_matrix;

struct COMPLEX {
    double real;
    double imaginary;
};

struct COMPLEX_MATRIX {
    complex matrix[M][M];
};

complex complexGen(int i) {
    srand(i + (int)time(NULL));
    
    complex num;
    num.real = rand() % 11 - 5;
    num.imaginary = rand() % 11 - 5;
    
    return num;
}

complex_matrix matrixGen(int num) {
    complex_matrix arr;
    for (int i = 0; i < M; i++)
        for (int j = 0; j < M; j++) {
            // provides random
            // sleep(i + j);
            arr.matrix[i][j] = complexGen(num + i + j);
        }
    return arr;
}

complex complexSum(complex a, complex b) {
    complex c;
    c.real = a.real + b.real;
    c.imaginary = a.imaginary + b.imaginary;
    return c;
}

complex complexMult(complex a, complex b) {
    complex c;
    c.real = (a.real * b.real - a.imaginary * b.imaginary);
    c.imaginary = (a.real * b.imaginary + a.imaginary * b.real);
    return c;
}

complex_matrix matrixMult(complex_matrix a, complex_matrix b) {
    complex_matrix c;
    for (int i = 0; i < M; i++)
        for (int j = 0; j < M; j++) {
            c.matrix[i][j].real = 0;
            c.matrix[i][j].real = 0;
            for (int k = 0; k < M; k++)
                c.matrix[i][j] = complexSum(c.matrix[i][j], complexMult(a.matrix[i][k], b.matrix[k][j]));
                
        }
    return c;
}

void complexPrint(complex num) {
    if (num.imaginary >= 0)
        printf("%.1f + %.1fi\t", num.real, num.imaginary);
    else
        printf("%.1f - %.1fi\t", num.real, fabs(num.imaginary));
}

void matrixPrint(complex_matrix arr) {
    puts("");
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++)
            complexPrint(arr.matrix[i][j]);
        puts("");
    }
    puts("");
}

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // create type for complex number
    MPI_Datatype _MPI_COMPLEX;
    MPI_Type_contiguous(2, MPI_DOUBLE, &_MPI_COMPLEX);
    MPI_Type_commit(&_MPI_COMPLEX);
    
    // create type for complex matrix
    MPI_Datatype MPI_COMPLEX_MATRIX;
    MPI_Type_contiguous(M * M, _MPI_COMPLEX, &MPI_COMPLEX_MATRIX);
    MPI_Type_commit(&MPI_COMPLEX_MATRIX);
    
    for (int i = 0; i < size; i++) {
        if (rank == 0 && i == 0) {
            complex_matrix m_0 = matrixGen(i);
            printf("Process %d:\nm_0: (%d)", rank, i);
            matrixPrint(m_0);
            
            MPI_Send(&m_0, 1, MPI_COMPLEX_MATRIX, 1, 0, MPI_COMM_WORLD);
        } else if (rank > 0 && rank == i) {
            complex_matrix recv_m_i;
            MPI_Recv(&recv_m_i, 1, MPI_COMPLEX_MATRIX, i - 1, 0, MPI_COMM_WORLD, &status);
            
            printf("Process %d:\nrecv_m_%d: (%d)", rank, status.MPI_SOURCE, i);
            matrixPrint(recv_m_i);
            
            complex_matrix m_i = matrixGen(i);
            printf("m_%d: (%d)", rank, i);
            matrixPrint(m_i);
            
            complex_matrix res_m = matrixMult(recv_m_i, m_i);
            printf("m_%d x m_%d: (%d)", status.MPI_SOURCE, rank, i);
            matrixPrint(res_m);
            
            if (rank != (size - 1))
                MPI_Send(&res_m, 1, MPI_COMPLEX_MATRIX, i + 1, 0, MPI_COMM_WORLD);
        }
    }
    
    // type invalidation
    MPI_Type_free(&_MPI_COMPLEX);
    MPI_Type_free(&MPI_COMPLEX_MATRIX);
    
    MPI_Finalize();
    
    return 0;
}

//    void printSize(MPI_Datatype type) {
//        int size = 0;
//        MPI_Type_size(type, &size);
//        printf("size = %d bytes\n", size);
//    }
//
//    void printExtent(MPI_Datatype type) {
//        MPI_Aint lb = 0, extent = 0;
//        MPI_Type_get_extent(type, &lb, &extent);
//        printf("lb = %ld\nextent = %ld bytes\n", lb, extent);
//        puts("");
//    }
//
//    if (rank == 0) {
//        int arr[2][2] = {{1, 2}, {1, 1}};
//        MPI_Send(&arr, 4, MPI_INT, 1, 0, MPI_COMM_WORLD);
//    } else if (rank == 1) {
//        int brr[2][2];
//        MPI_Recv(&brr, 4, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
//        for (int i = 0; i < 2; i++) {
//            for (int j = 0; j < 2; j++)
//                printf("%d ", brr[i][j]);
//            puts("");
//        }
//    }
//
//    MPI_Datatype TWO_INT;
//    MPI_Type_contiguous(2, MPI_INT, &TWO_INT);
//    MPI_Type_commit(&TWO_INT);
//    printSize(TWO_INT);
//    puts("");
//
//    MPI_Datatype SIX_INT;
//    MPI_Type_vector(2, 3, 4, MPI_INT, &SIX_INT);
//    MPI_Type_commit(&SIX_INT);
//    printSize(SIX_INT);
//    printExtent(SIX_INT);
//
//    MPI_Datatype FIVE_INT;
//    int blocklengths[] = {1, 4};
//    int displacements[] = {0, 2};
//    MPI_Type_indexed(2, blocklengths, displacements, MPI_INT, &FIVE_INT);
//    MPI_Type_commit(&FIVE_INT);
//    printSize(FIVE_INT);
//    printExtent(FIVE_INT);
//
//    MPI_Datatype STRUCT;
//    MPI_Datatype types[] = {MPI_INT, MPI_SHORT, MPI_CHAR};
//    int blocklen[] = {1, 6, 4};
//    MPI_Aint displace[] = {0, 12, 26};
//    MPI_Type_create_struct(3, blocklen, displace, types, &STRUCT);
//    MPI_Type_commit(&STRUCT);
//    printSize(STRUCT);
//    printExtent(STRUCT);
//
//
//    MPI_Type_free(&TWO_INT);
//    MPI_Type_free(&SIX_INT);
//    MPI_Type_free(&FIVE_INT);
//    MPI_Type_free(&STRUCT);
