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
    int Asize = 1000;
	int rowCount;
	int ifBound;
	int curriter, stopiter;
	double Norm, MaxNorm;
	double eps;
	MPI_Status stat;


    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	// N / size must be integer!
	if(!(Asize%size)){
		ifBound = (rank%(size-1))==0;
		rowCount = Asize / size + 2 - ifBound;
		eps = 0.0001;
		curriter = 0;
		stopiter = 2000;
		}
	else{
		printf("\n\nERROR! N / size must be integer!\n\n");
		return 0;
	}



	// Initialization
    if(rank==0){
		printf("\nInitialization of solving Laplase's equasion. Domain size is %d*%d, number of processors is %d\n\n", Asize, Asize, size);

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

		printf("\nCalculation started. %d iterations to go\n\n\n", stopiter);
		/*
		for(int i = 0; i <  size; i++)
			printf("\ndispls is %d, sendcound is%d\n", displs[i], sendcount[i]);
*/
	}

	A = new double [rowCount  * Asize];
	B = new double [ (Asize-2) * (rowCount-2)];
	

    MPI_Scatterv(initA, sendcount, displs, MPI_DOUBLE, A, rowCount * Asize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	do{
		MaxNorm = -1.0;
		ptrB = B;
		ptrA = A + Asize + 1;

		// Caculate current ineration
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

		// Copy current answer to next iteration
		ptrB = B;
		ptrA = A + Asize + 1;
		for(int i = 0; i < rowCount - 2; i++){
			for(int j = 0; j < Asize - 2; j++)
				*ptrA++ = *ptrB++;
			ptrA += 2;
		}

		// Fill remaining values				
		if(rank < size - 1)
			MPI_Send(B + (Asize-2) * (rowCount-3), Asize-2, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
		if(rank > 0)
			MPI_Recv(A + 1, Asize-2, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &stat);
		if(rank > 0)
			MPI_Send(B, Asize-2, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
		if(rank < size - 1)
			MPI_Recv(A + Asize * (rowCount-1) + 1, Asize-2, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &stat);


		// Find MaxNorm among all nodes.
		MPI_Reduce(&MaxNorm, &Norm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Bcast(&Norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		curriter++;
		if(!(curriter % 50) && (rank == 0)){
           printf("Completed %d iteration, value of T(20,20) is %.3f\n", curriter, *(A + 20*Asize + 20));
        }
	}while((Norm > eps) && (curriter < stopiter));

	

	// Collect the answer
	MPI_Gather(A + (rank!=0) * Asize, Asize * Asize / size, MPI_DOUBLE, initA, Asize * Asize / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// Ptint out answer
	if(rank == 0){
		ptrA = initA;
		if(curriter==stopiter)		
			printf("\n\nSolution NOT found. 0 iteration out of %d left.\nCurrent error is %.5f. Target error is %.5f\n", stopiter, Norm, eps);
		else			
			printf("\n\nSolution found and took %d iteration.\nCurrent error is %.5f. Target error is %.5f\n", curriter, Norm, eps);
/*
        printf("\nMatrix A is:\n");
        for(int i = 0; i < Asize; i++){
            for(int j = 0; j < Asize; j++)
                printf("%.3f\t", *ptrA++);
            printf("\n");
        }
*/
        delete []initA;
		delete []sendcount;
        delete []displs;
    }
    delete []A;
    delete []B;


    MPI_Finalize();
}
