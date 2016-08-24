#include <stdlib.h>
#include "f2c.h"
#include <stdio.h>
#include <math.h>
#include "clapack.h"
#define RowNumber 100
#define l 0
#define rmin 0.01
#define rmax 101.01
#define step (rmax-rmin)/(RowNumber+1)
#define a 0
#define b 0

float r[RowNumber],Matrix[RowNumber][RowNumber]={0};

void ConstructMatrix();

void ComputeEvalue();

int main()
{
	float Eigenvalue[RowNumber];
	for(int i=0;i<RowNumber;i++)
	{
		r[i]=rmin+step*(i+1);
	}

	ConstructMatrix();

	//ComputeEvalue(Eigenvalue);

	for(int i=0;i<RowNumber;i++)
		printf("%f ",Eigenvalue[i]);
	printf("\n");

	return 0;
}

void SecondDerivative()
{
	Matrix[0][0]-=(1+r[0])*(1+r[0])/(step*step)*log(1+r[0])*(a-2);
	Matrix[0][1]-=(1+r[0])*(1+r[0])/(step*step)*log(1+r[0]);
	for(int i=1;i<RowNumber-1;i++)
	{
		Matrix[i][i-1]-=(1+r[i])*(1+r[i])/(step*step)*log(1+r[i]);	
		Matrix[i][i]+=(1+r[i])*(1+r[i])/(step*step)*log(1+r[i])*2;	
		Matrix[i][i+1]-=(1+r[i])*(1+r[i])/(step*step)*log(1+r[i]);	
	}
	Matrix[RowNumber-1][RowNumber-2]-=(1+r[RowNumber-1])*(1+r[RowNumber])/(step*step)*log(1+r[RowNumber-1]);
	Matrix[RowNumber-1][RowNumber-1]-=(1+r[RowNumber-1])*(1+r[RowNumber])/(step*step)*log(1+r[RowNumber-1])*(b-2);
}

void FirstDerivative()
{
	Matrix[0][0]+=(1+r[0])/step*a;
	Matrix[0][1]-=(1+r[0])/step;
	for(int i=1;i<RowNumber-1;i++)
	{
		Matrix[i][i-1]+=(1+r[i])/step;
		Matrix[i][i+1]-=(1+r[i])/step;
	}
	Matrix[RowNumber-1][RowNumber-2]+=(1+r[RowNumber-1])/step;
	Matrix[RowNumber-1][RowNumber-1]-=(1+r[RowNumber-1])/step*b;
}

void ConstructMatrix()
{
	SecondDerivative();
	FirstDerivative();

	for(int i=0;i<RowNumber;i++)
	{
		Matrix[i][i]+=1+l*(l+1)*(1+r[i])*(1+r[i])/(r[i]*r[i])*log(1+r[i]);
	}
	
}

