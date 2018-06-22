#ifndef _GATES_H_
#define _GATES_H_

#include "common.h"
#include <vector>
#include <map>
#include <string>
#include <complex.h>
#include <math.h>

using namespace std;

string CNot(int qubits, int ctrl, int target, int cv = 1);
string Toffoli(int qubits, int ctrl1, int ctrl2, int target, int cv = 3);
string Controlled1(int qubits, int ctrl, int target, string op, int cv = 1);
string Controlled2(int qubits, int ctrl1, int ctrl2, int target, string op, int cv = 3);
string ControlledN(int qubits, vector<int> ctrls, int target, string op, int cv = -1);
string Pauli_X(int qubits, int reg, int width = 1);
string Pauli_Z(int qubits, int reg, int width = 1);
string Hadamard(int qubits, int reg, int width = 1);

string concatena(vector <string> vec, int size, bool rev = false);


class Gates{
public:
	static map <string, float complex*> list;
	Gates();
	~Gates();
	void init();
	float complex* getMatrix(string gateName);
	bool addGate(string name, float complex* matrix);
	bool addGate(string name, float complex a0, float complex a1, float complex a2, float complex a3);
};

#endif
