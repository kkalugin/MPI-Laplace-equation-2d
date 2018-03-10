#include <stdio.h>
#include <mpi.h>
#include <math.h>



int main(int argc, char *argv[]){
    double* A = NULL;
    double* B = NULL;
	double* initA = NULL;
	int* sendcount = NULL;
	int* displs = NULL;
		int* pSend = NULL;
		int* pDispls = NULL;
		double* ptrA = NULL;
		double* ptrB = NULL;

    int rank, size;
    int Asize = 12;
	int rowCount;
	int ifBound;
	int maxiter;
	double Norm, MaxNorm;
	double eps;
	MPI_Status stat;

	
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	
	if(!(Asize%size)){
		ifBound = (rank%(size-1))==0;
		rowCount = Asize / size + 2 - ifBound;
		eps = 0.001;
		maxiter = 100;
		}
	else{
		printf("\n\nERROR\n\n");
		return 0;
	}
		
  
	
    if(rank==0){
        initA = new double [Asize * Asize];
		sendcount = new int [size];
        displs = new int [size];
		ptrA = initA;
		
		for(int i = 0; i < Asize; i++)
			*ptrA++ = 1;
        for(int i = 1; i < Asize - 1; i++){
			*ptrA++ = 1;
            for(int j = 1; j < Asize - 1; j++)
                *ptrA++ = 0;
			*ptrA++ = 1;
		}
		for(int i = 0; i < Asize; i++)
			*ptrA++ = 1;
		
		pSend = sendcount;
		pDispls = displs;
		
		*pSend++ = rowCount * Asize;
		*pDispls++ = 0;
		*pDispls = (rowCount - 2) * Asize;
		for(int i = 1; i <  size - 1; i++){
			*pSend++ = (rowCount + 1) * Asize;
			*++pDispls = displs[1] + i * (rowCount - 1) * Asize;
		}
		*pSend = rowCount * Asize;
		
		printf("\n%d iterations to go\n", maxiter);
	}
	
	A = new double [rowCount  * Asize];
	B = new double [ (Asize-2) * (rowCount-2)];
	
    MPI_Scatterv(initA, sendcount, displs, MPI_DOUBLE, A, rowCount * Asize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	do{
		MaxNorm = -1.0;
		ptrB = B;
		ptrA = A + Asize + 1;
		for(int i = 0; i < rowCount - 2; i++){
			for(int j = 0; j < Asize - 2; j++){
				*ptrB = 1.0/4 * (*(ptrA-1) + *(ptrA + 1) + *(ptrA-Asize) + *(ptrA+Asize));
				ptrA;
				Norm = fabs(*ptrB - *ptrA);
				if(Norm > MaxNorm)
					MaxNorm = Norm;
				ptrB++;
				ptrA++;
			}
			ptrA += 2;
		}
						
		ptrB = B;
		ptrA = A + Asize + 1;
		for(int i = 0; i < rowCount - 2; i++){
			for(int j = 0; j < Asize - 2; j++)
				*ptrA++ = *ptrB++;
			ptrA += 2;
		}
		
		if(ifBound){
		 /* FIXME */
		 if(rank==0){
			 MPI_Send(B + (Asize-2) * (rowCount-3), Asize-2, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
			 MPI_Recv(A + Asize * (rowCount-1) + 1, Asize-2, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &stat);
		 }
		 else{
			MPI_Send(B, Asize-2, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
			MPI_Recv(A + 1, Asize-2, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &stat);
		 }	
		}		
		else{
			MPI_Send(B, Asize-2, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
			MPI_Send(B + (Asize-2) * (rowCount-3), Asize-2, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
			MPI_Recv(A + 1, Asize-2, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &stat);
			MPI_Recv(A + Asize * (rowCount-1) + 1, Asize-2, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &stat);
		}
		
		MPI_Reduce(&Norm, &MaxNorm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Bcast(&MaxNorm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		maxiter--;
	
	}while(MaxNorm > eps && maxiter);
	

	MPI_Gather(A + (rank!=0) * Asize, Asize * Asize / size, MPI_DOUBLE, initA, Asize * Asize / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
	if(rank==0){
		ptrA = initA;
		printf("\n\nSolution found\n%d iteration left. Current error is %.4f\n", maxiter, MaxNorm);
        printf("\nMatrix A is:\n");
        for(int i = 0; i < Asize; i++){
            for(int j = 0; j < Asize; j++)
                printf("%.3f\t", *ptrA++);				
            printf("\n");
        }
        delete []initA;		
    }
    delete []A;
    delete []B;


    MPI_Finalize();
}
