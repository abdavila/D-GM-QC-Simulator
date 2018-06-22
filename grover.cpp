#include "dgm.h"
#include "common.h"
#include "gates.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

string ControledZ(int qubits);
string Oracle1(long qubits, long int value);

#define NUM_AMOSTRAS 1
#define NUM_THREADS 4

#define EXC_RANGE 0

float print_statistic(vector <float> amostra){
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

	cout << med << "\t" << desv;

	return med;

}


float HadamardNQB(long qubits, long num_of_it, int qb_pos, int multi_gpu, int qbs_region, int coalesc, int tam_block, int rept){
	DGM dgm;

	dgm.exec_type = t_GPU;
	dgm.n_threads = 1;
	dgm.multi_gpu = multi_gpu;
	dgm.qubits = qubits;

	dgm.qbs_region = qbs_region;
	dgm.coalesc = coalesc;
	dgm.tam_block = tam_block;
	dgm.rept = rept;
	//dgm.multi_gpu = 2;

	dgm.allocateMemory();
	dgm.setMemoryValue(0);

	dgm.setFunction(Hadamard(qubits, qb_pos, 1), num_of_it);


	struct timeval timev, tvBegin, tvEnd;
	float t;


	gettimeofday(&tvBegin, NULL);
	dgm.execute(1);
	gettimeofday(&tvEnd, NULL);
	timeval_subtract(&timev, &tvEnd, &tvBegin);
	t = timev.tv_sec + (timev.tv_usec / 1000000.0);

	dgm.freeMemory();

	return t;
}

long* GroverSpectro(long qubits, long value, int spec_block){
	DGM dgm;

	dgm.qubits = qubits;
	dgm.exec_type = t_SPEC;
	dgm.spec_block = spec_block;
	
	dgm.allocateSpectro();

	
	string H = Hadamard(qubits, 0, qubits);
	string orcl = Oracle1(qubits, value);
	string CZ = ControledZ(qubits);

	vector <string> grover_step;

	grover_step.push_back(orcl);
	
	for (int i = 1; i < qubits; i++){
		grover_step.push_back(Hadamard(qubits, i, 1));
		grover_step.push_back(Pauli_X(qubits, i, 1));
	}
	
	grover_step.push_back(CZ);

	for (int i = qubits-1; i >= 1; i--){
		grover_step.push_back(Pauli_X(qubits, i, 1));
		grover_step.push_back(Hadamard(qubits, i, 1));
	}


	int num_of_it = 1; //min(50, (int) (M_PI/4.0*sqrt(1<<(qubits-1))));

	dgm.setFunction(H);
	dgm.setFunction(grover_step, num_of_it, false);

	dgm.execute(1);

	return dgm.spectro;

}

float Grover(long qubits, long value, int type, int multi_gpu, int qbs_region, int coalesc, int tam_block, int rept){
	DGM dgm;

	dgm.exec_type = type;
	dgm.n_threads = multi_gpu;
	dgm.multi_gpu = multi_gpu;
	dgm.qubits = qubits;

	dgm.qbs_region = qbs_region;
	dgm.coalesc = coalesc;
	dgm.tam_block = tam_block;
	dgm.rept = rept;

	dgm.allocateMemory();
	dgm.setMemoryValue(1<<(qubits-1));

	
	string H = Hadamard(qubits, 0, qubits);
	string orcl = Oracle1(qubits, value);
	string CZ = ControledZ(qubits);

	vector <string> grover_step;

	grover_step.push_back(orcl);
	
	for (int i = 1; i < qubits; i++){
		grover_step.push_back(Hadamard(qubits, i, 1));
		grover_step.push_back(Pauli_X(qubits, i, 1));
	}
	
	grover_step.push_back(CZ);

	for (int i = qubits-1; i >= 1; i--){
		grover_step.push_back(Pauli_X(qubits, i, 1));
		grover_step.push_back(Hadamard(qubits, i, 1));
	}


	int num_of_it = (int) (M_PI/4.0*sqrt(1<<(qubits-1)));
	long result = 0;

	dgm.setFunction(H);
	dgm.setFunction(grover_step, num_of_it, false);

	dgm.execute(1);


	for (long i = 1; i < qubits; i++)
		result = (result << 1) + dgm.measure(i);

	dgm.freeMemory();

    //cout << result << endl;

	//cout << num_of_it << endl;
	//if (result == value) return true;
	return 0;

}

string Oracle1(long qubits, long int value){
	vector <string> t(qubits);
	int ctrl_v;

	for (int i = qubits - 1; i >= 1; i--){
		ctrl_v = value&1;
		value = value >> 1;
		if (ctrl_v) t[i] = "Control1(1)";
		else t[i] = "Control1(0)";
	}
	t[0] = "Target1(X)";

	return concatena(t, qubits);

}

string ControledZ(int qubits){
	vector <string> cz;
	cz.push_back("ID");
	for (int i = 0; i < qubits-2; i++) cz.push_back("Control1(1)");
	cz.push_back("Target1(Z)");

	return concatena(cz, qubits);
}

int main(int argc, char **argv){
	srand(time(NULL));
	if (argc < 3) return 0;

	struct timeval timev, tvBegin, tvEnd;
	float t;

	int num_qubits = atoi(argv[1]);
	int exec_type = atoi(argv[2]);
	int num_threads = atoi(argv[3]);
	int value = 100;
	//int spec_block = 1024;

	//if (argc == 3) spec_block = atoi(argv[2]);
//	int it = atoi(argv[2]);
//	int qb_range = 1;
//	if (argc == 4) qb_range = atoi(argv[3]) + 1;

	gettimeofday(&tvBegin, NULL);
	Grover(num_qubits, value, exec_type, num_threads, 9, 4, 256, 2);
	//Grover(qb, value, t_GPU, 1, qbs_region, coalesc, tam_block, rept);
	gettimeofday(&tvEnd, NULL);
	timeval_subtract(&timev, &tvEnd, &tvBegin);
	t = timev.tv_sec + (timev.tv_usec / 1000000.0);

	cout << t << endl;



	/*
	vector <float> amostras;


	long *spectro;
	spectro = GroverSpectro(num_qubits, value, spec_block);

	long size = spec_block;
	cout << size << endl;	

	for (long i = 0; i < size; i++){
		for (long j = 0; j < size; j++){
			
			if (spectro[j*size + i]){
				cout << long(spectro[j*size + i]) << " ";
			}
			else{
				cout << "0 ";	
			}
			
		}
		cout << endl;
	}
	*/


	/*
	cout << "##### CPU #####" << endl;
	for (int qb = num_qubits; qb < num_qubits + qb_range; qb++){
		cout << "QUBITS: " << qb << endl;
		for (int a = 0; a < NUM_AMOSTRAS; a++){
			gettimeofday(&tvBegin, NULL);
			Grover(qb, value, t_CPU, 1);
			gettimeofday(&tvEnd, NULL);
			timeval_subtract(&timev, &tvEnd, &tvBegin);
			t = timev.tv_sec + (timev.tv_usec / 1000000.0);
			amostras.push_back(t);
		}
		print_statistic(amostras);
		amostras.clear();
	}


	cout << "##### PAR_CPU #####" << endl;
	for (int qb = num_qubits; qb < num_qubits + qb_range; qb++){
		cout << "QUBITS: " << qb << endl;
		for (int th = 1; th <= NUM_THREADS; th*=2){
			cout << "Threads: " << th << endl;
			for (int a = 0; a < NUM_AMOSTRAS; a++){
				gettimeofday(&tvBegin, NULL);
				Grover(qb, value, t_PAR_CPU, th);
				gettimeofday(&tvEnd, NULL);
				timeval_subtract(&timev, &tvEnd, &tvBegin);
				t = timev.tv_sec + (timev.tv_usec / 1000000.0);
				amostras.push_back(t);
			}
			print_statistic(amostras);
			amostras.clear();
		}
	}
	*/

	//float matrix_result[7][4];
/*
	cout << "##### GPU #####" << endl;
	for (int qb = num_qubits; qb < num_qubits + qb_range; qb++){
		cout << "QUBITS: " << qb << endl;

		for (int coalesc = 4; coalesc <= 4; coalesc++){
//			cout << "COALESC: " << coalesc << endl;
			int qbs_region = 9;

			//int qbs_region = coalesc+1; qbs_region <= 9; qbs_region++){
			int tam_block = 256;
			int rept = 4;
				
			amostras.clear();
			for (int a = 0; a < NUM_AMOSTRAS; a++){
				gettimeofday(&tvBegin, NULL);
				Grover(qb, value, t_GPU, 1, qbs_region, coalesc, tam_block, rept);
				gettimeofday(&tvEnd, NULL);
				timeval_subtract(&timev, &tvEnd, &tvBegin);
				t = timev.tv_sec + (timev.tv_usec / 1000000.0);
				amostras.push_back(t);
			}
			print_statistic(amostras);
			cout << "\t";

			amostras.clear();
			for (int a = 0; a < NUM_AMOSTRAS; a++){
				gettimeofday(&tvBegin, NULL);
				Grover(qb, value, t_GPU, 2, qbs_region, coalesc, tam_block, rept);
				gettimeofday(&tvEnd, NULL);
				timeval_subtract(&timev, &tvEnd, &tvBegin);
				t = timev.tv_sec + (timev.tv_usec / 1000000.0);
				amostras.push_back(t);
			}
			print_statistic(amostras);
			cout << endl;
			//}
//			cout << "####################" << endl;
		}
		cout << "####################" << endl; 
	}
*/
	/*

	cout << "##### GPU #####" << endl;
	for (int qb = num_qubits; qb < num_qubits + qb_range; qb++){
		cout << "QUBITS: " << qb << endl;

		for (int coalesc = 0; coalesc <= 6; coalesc++){
			//cout << "COALESC: " << coalesc << endl;
			int qbs_region = 7;

			//for (int qbs_region = coalesc+1; qbs_region <= 9; qbs_region++){
			//	cout << "QBS_REGION: " << qbs_region << endl;

			int tam_block = 256;


			//for (int rept = 1; rept <= 8 && ((tam_block * rept) <= 2048); rept*=2){
			int rept = 4;
			//cout << "REPT: " << rept << endl;
			for (int a = 0; a < NUM_AMOSTRAS; a++){
				t = HadamardNQB(qb, it, 0, 1, qbs_region, coalesc, tam_block, rept);
				amostras.push_back(t);
			}
			cout << print_statistic(amostras) << "\t";
			//matrix_result[coalesc][0] = print_statistic(amostras);
			amostras.clear();
			for (int a = 0; a < NUM_AMOSTRAS; a++){
				t = HadamardNQB(qb, it, 1, 1, qbs_region, coalesc, tam_block, rept);
				amostras.push_back(t);
			}
			cout << print_statistic(amostras) << "\t";
			//matrix_result[coalesc][1] = print_statistic(amostras);
			amostras.clear();
			for (int a = 0; a < NUM_AMOSTRAS; a++){
				t = HadamardNQB(qb, it, 0, 2, qbs_region, coalesc, tam_block, rept);
				amostras.push_back(t);
			}
			cout << print_statistic(amostras) << "\t";
			//matrix_result[coalesc][2] = print_statistic(amostras);
			amostras.clear();
			for (int a = 0; a < NUM_AMOSTRAS; a++){
				t = HadamardNQB(qb, it, 1, 2, qbs_region, coalesc, tam_block, rept);
				amostras.push_back(t);
			}
			cout << print_statistic(amostras) << endl;
			//matrix_result[coalesc][3] = print_statistic(amostras);
			amostras.clear();
			//}
			//cout << "####################" << endl;
			//}

			//}
		}
		cout << "####################" << endl; 

	}
	*/
}
