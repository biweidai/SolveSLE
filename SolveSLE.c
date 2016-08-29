#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mkl_lapacke.h"

#define N 1999                                                                                 //the order of the Matrix 
#define address "/home/biwei/eigenvector"                                                      //the address to print the eigenvector
#define l 0                                                                                     
#define rmin 0.0
#define rmax 20.0 
#define step (double)(rmax-rmin)/(N+1.0)
#define b 1                                                                                     //rmax=b*r[N-1] boundary condition

extern void print_eigenvalues( char* desc, MKL_INT n, double* wr, double* wi );
extern void print_eigenvectors(MKL_INT n, double* wr, double* wi, double* v, MKL_INT ldv );

double r[N];
double Matrix[N*N]={0};
double vr[N*N];

void ConstructMatrix();

int main()
{
	
	int i;
	for(i=0;i<N;i++)                                                                       //initialize r[]
	{
		r[i]=rmin+step*(i+1);
	}

	ConstructMatrix();
/*
	printf("Matrix:\n%f %f %f %f\n",Matrix[0],Matrix[1],Matrix[2],Matrix[3]);              //print the nonzero term of the matrix 
	for(i=1;i<(N-1);i++)
	{
		printf("%f %f %f\n",Matrix[i*N+i-1],Matrix[i*N+i],Matrix[i*N+i+1]);
	}
	printf("%f %f\n",Matrix[N*N-2],Matrix[N*N-1]);	
*/        
	/* Locals */
        MKL_INT n = N, lda = N, ldvl = N, ldvr = N, info;
        /* Local arrays */
        double wr[N], wi[N], *vl=NULL; 
       
        /* Solve eigenproblem */
        info = LAPACKE_dgeev( LAPACK_ROW_MAJOR, 'N', 'V', n, Matrix, lda, wr, wi,
                        vl, ldvl, vr, ldvr );
        /* Check for convergence */
        if( info > 0 ) {
                printf( "The algorithm failed to compute eigenvalues.\n" );
                exit( 1 );
        }
        /* Print eigenvalues */
        print_eigenvalues( "Eigenvalues", n, wr, wi );
        print_eigenvectors( n, wr, wi, vr, ldvr );
        exit( 0 );
}

void SecondDerivative()                                                                    //construct the matrix of second derivative term
{
	Matrix[0]-=2*(1.0+r[0])*(1.0+r[0])/(step*step)*log(1.0+r[0]);
	Matrix[1]+=5*(1.0+r[0])*(1.0+r[0])/(step*step)*log(1.0+r[0]);
	Matrix[2]-=4*(1.0+r[0])*(1.0+r[0])/(step*step)*log(1.0+r[0]);
	Matrix[3]+=(1.0+r[0])*(1.0+r[0])/(step*step)*log(1.0+r[0]);
	int i;
	for(i=1;i<N-1;i++)
	{
		Matrix[i-1+i*N]-=(1.0+r[i])*(1.0+r[i])/(step*step)*log(1.0+r[i]);	
		Matrix[i+i*N]+=(1.0+r[i])*(1.0+r[i])/(step*step)*log(1.0+r[i])*2.0;	
		Matrix[i+1+i*N]-=(1.0+r[i])*(1.0+r[i])/(step*step)*log(1.0+r[i]);	
	}
	Matrix[N*N-2]-=(1.0+r[N-1])*(1.0+r[N-1])/(step*step)*log(1.0+r[N-1]);
	Matrix[N*N-1]-=(1.0+r[N-1])*(1.0+r[N-1])/(step*step)*log(1.0+r[N-1])*(b-2.0);
}

void FirstDerivative()                                                                    // construct the matrix of first derivative term
{
	Matrix[0]+=2*(1.0+r[0])/(step);
	Matrix[1]-=(1.0+r[0])/(step);
	Matrix[2]-=2*(1.0+r[0])/(step);
	Matrix[3]+=(1.0+r[0])/(step);
	int i;
	for(i=1;i<N-1;i++)
	{
		Matrix[i-1+i*N]+=(1.0+r[i])/(step);
		Matrix[i+1+i*N]-=(1.0+r[i])/(step);
	}
	Matrix[N*N-2]+=(1.0+r[N-1])/(step);
	Matrix[N*N-1]-=b*(1.0+r[N-1])/(step);
}

void ConstructMatrix()
{
	SecondDerivative();

	FirstDerivative();

	int i;
	for(i=0;i<N;i++) 
	{
		Matrix[i*N+i]+=1.0+l*(l+1.0)*(1.0+r[i])*(1.0+r[i])*log(1.0+r[i])/(r[i]*r[i]);
	}
}

/* Auxiliary routine: printing eigenvalues */
void print_eigenvalues( char* desc, MKL_INT n, double* wr, double* wi ) {
        MKL_INT j;
        printf( "\n %s\n", desc );
   for( j = 0; j < n; j++ ) {
      if( wi[j] == (float)0.0 ) {
	if(wr[j]>50||wr[j]<0) continue;
         printf( " %6.15f", wr[j] );
      } else {
         printf( " (%6.2f,%6.2f)", wr[j], wi[j] );
      }
   }
   printf( "\n" );
}

/* Auxiliary routine: printing eigenvectors */
void print_eigenvectors( MKL_INT n, double* wr, double* wi, double* v, MKL_INT ldv )   //print the first several eigenvectors to /home/biwei/eigenvector 
{
        MKL_INT i, j;
	FILE *fp;
	double t;
	fp=fopen(address,"w");
	for(i=0;i<N;i++)
		fprintf(fp,"%f ",r[i]);
	fprintf(fp,"\n");
	for( i = 0; i < n; i++ ) 
	{
		if(wi[i] == (float)0.0 && wr[i]>0 && wr[i]<50)
		{
			t=1.0/v[i];	
			fprintf(fp,"%3.3f\n",wr[i]);
			for(j=0;j<n;j++)
				fprintf(fp,"%3.5f ",v[ldv*j+i]*t);
                 
			fprintf(fp,"\n");
		}	

	}
}
