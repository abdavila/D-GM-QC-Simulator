#include "dgm.h"
#include "common.h"
#include "gates.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

void Fuzzy(){
	
	vector <qubit> entrada;

	float x1, x2, y1, y2;

	x1 = 1.0/3;
	x2 = 1.0/2;
	y1 = 1.0/3;
	y2 = 1.0/2;

	entrada.push_back(qubit(sqrt(x1)));
	entrada.push_back(qubit(sqrt(x2)));
	entrada.push_back(qubit(sqrt(y1)));
	entrada.push_back(qubit(sqrt(y2)));
	
	entrada.push_back(ZERO);
	entrada.push_back(ZERO);
	entrada.push_back(ZERO);

	entrada[0].print();
	entrada[1].print();
	entrada[2].print();
	entrada[3].print();

	int qubits = entrada.size();


	vector <int> ctrls;
	ctrls.push_back(0);
	ctrls.push_back(1);
	ctrls.push_back(2);
	ctrls.push_back(3);
	ctrls.push_back(4);

	vector <string> alg;
//	alg.push_back(Toffoli(qubits, 0, 1, 2));


	alg.push_back(Toffoli(qubits, 0, 2, 4));
	alg.push_back(Pauli_X(qubits, 1));
	alg.push_back(Pauli_X(qubits, 3));
	alg.push_back(Toffoli(qubits, 1, 3, 5));
	alg.push_back(Pauli_X(qubits, 1));
	alg.push_back(Pauli_X(qubits, 3));
	alg.push_back(Pauli_X(qubits, 5));

	alg.push_back(Pauli_X(qubits, 5));
	alg.push_back(Pauli_X(qubits, 4));
	alg.push_back(Toffoli(qubits, 4, 5, 6));
	alg.push_back(Pauli_X(qubits, 5));
	alg.push_back(Pauli_X(qubits, 4));

	//alg.push_back(Hadamard(qubits, 0));
	//alg.push_back(Controlled1(qubits, 0, 1, "R1"));
	//alg.push_back(Hadamard(qubits, 1));	
	//alg.push_back(Controlled1(qubits, 2, 0, "R2"));

	//alg.push_back(Controlled1(qubits, 2, 1, "R1"));
	//alg.push_back(Hadamard(qubits, 2));

	for (int i = 0; i < alg.size(); i++) cout << alg[i] << endl;
	cout << endl;

	DGM dgm;
	dgm.qubits = qubits;
	//dgm.exec_type = t_PAR_CPU;
	//dgm.num_threads = threads;
	dgm.setMemory(genMem(entrada));

	for (int i = 0; i < alg.size(); i++) dgm.executeFunction(alg[i]);

	vector <int> mq;

	printMem(dgm.r_mem, qubits);

	cout << endl;

	//mq.push_back(4);
	//mq.push_back(5);
	//mq.push_back(6);
	//mq.push_back(2);
	//mq.push_back(3);
	//dgm.measure(mq);

	//printMem(dgm.r_mem, qubits);
	cout << endl;

	mq.clear();
	mq.push_back(2);
	dgm.measure(mq);

return;
	cout << endl;

	mq.clear();
	mq.push_back(5);
	dgm.measure(mq);

	cout << endl;

	mq.clear();
	mq.push_back(6);
	dgm.measure(mq);

}


int main(int argc, char **argv){
	Fuzzy();
}

