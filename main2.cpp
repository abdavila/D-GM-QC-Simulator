#include <iostream>
#include <vector>
#include "dgm.h"
#include "complex.h"
#include "gates.h"
#include <unistd.h>

using namespace std;

#define NUM_AMOSTRAS 30
#define NUM_THREADS 4

#define EXC_RANGE 5


bool print_statistic(vector <float> amostra){
	float med, desv;
	med = desv = 0;
	int number = amostra.size();

	sort(amostra.begin(), amostra.end());

	for (int i = EXC_RANGE; i < number - EXC_RANGE; i++){
//		cout << amostra[i] << endl;
		med += amostra[i];
	}
	med /= (number - EXC_RANGE*2);

	for (int i = EXC_RANGE; i < number - EXC_RANGE; i++){
		desv += pow((amostra[i] - med),2.0);
	}

	desv = sqrt(desv/(number-1)) / med * 100.0;

	//if (desv > 1) return false;

	cout << "\nMEDIA: " << med << "\t\tDESV: " << desv << endl;
	return true;
}

int main(int argc, char** argv){
	struct timeval timev, tvBegin, tvEnd;
	float t;

	int qubits = atoi(argv[1]);
	int m = atoi(argv[2]);
	int threads = atoi(argv[3]);
	//int had_pos = atoi(argv[2]);
	//int factor = atoi(argv[2]);
	//int stress = 1;
	//if (argc == 4)
	//	stress = atoi(argv[3]);

	//string function = "ID,ID,ID,ID,ID";
	vector <string> function;


//	for (int i = 0; i < stress; i++){
//		function.push_back(Hadamard(qubits, had_pos, 1));
//		function.push_back(Pauli_Z(qubits, had_pos, 1));
		//function.push_back(Pauli_Z(qubits, had_pos, 1));
//	}
	
	//float complex *state = (float complex*) calloc(pow(2, qubits), sizeof(float complex));
	//state[0] = 1;

	//struct timeval timev, tvBegin, tvEnd;
	//gettimeofday(&tvBegin, NULL);

	//long *spectro;

	//spectro = SpectroExecute(qubits, function);

	DGM dgm;

	dgm.qubits = qubits;
	dgm.exec_type = t_CPU;
	dgm.n_threads = threads;
	dgm.multi_gpu = 1;
	//dgm.qubits = qubits;

	vector <string> had = Controlled1(5, 0, 1, "H");

	dgm.allocateMemory();
	dgm.setMemoryValue(0);

	dgm.executeFunction(Hadamard(qubits, 0, qubits));

	gettimeofday(&tvBegin, NULL);
	dgm.measure(m);
	gettimeofday(&tvEnd, NULL);
	timeval_subtract(&timev, &tvEnd, &tvBegin);
	t = timev.tv_sec + (timev.tv_usec / 1000000.0);

	cout << t << endl;

	gettimeofday(&tvBegin, NULL);
	dgm.measure_coalesced(m);
	gettimeofday(&tvEnd, NULL);
	timeval_subtract(&timev, &tvEnd, &tvBegin);
	t = timev.tv_sec + (timev.tv_usec / 1000000.0);

	cout << t << endl;
	

/*
	long mem_size = pow(2,qubits);

	cout << (pow(2,qubits)) << endl;	

	for (long i = 0; i < mem_size; i++){
		for (long j = 0; j < mem_size; j++){
			//fwrite(&spectro[j*mem_size + i], sizeof(unsigned char), 1, image);
			cout << spectro[j*mem_size + i] << " ";
		}
		cout << endl;
	}
*/


	
	//gettimeofday(&tvEnd, NULL);
	//timeval_subtract(&timev, &tvEnd, &tvBegin);
	//float t = timev.tv_sec + (timev.tv_usec / 1000000.0);
	//cout << t << endl;

	//free(state);
	//printMem(state, 4);
	return 0;
}
