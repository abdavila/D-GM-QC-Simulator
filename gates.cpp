#include "gates.h"
#include <iostream>
#include <cmath>

using namespace std;

map <string, float complex*> Gates::list;

Gates::Gates(){
	init();
}

Gates::~Gates(){}

void Gates::init(){
	if (Gates::list.size() == 0){
		addGate("H", 1.0/sqrt(2), 1.0/sqrt(2), 1.0/sqrt(2), -1.0/sqrt(2));
		addGate("X", 0.0, 1.0, 1.0, 0.0);
		addGate("Z", 1.0, 0.0, 0.0, -1.0);
		addGate("Y", 0.0, 1.0*I, -1.0*I, 0.0);
		addGate("R1", 1.0, 0, 0.0, cpowf(M_E, M_PI*I));
		addGate("R2", 1.0, 0, 0.0, cpowf(M_E, M_PI*I/2.0));
		addGate("R3", 1.0, 0, 0.0, cpowf(M_E, M_PI*I/4.0));
	}
}

float complex* Gates::getMatrix(string gateName){
	if (Gates::list.find(gateName) == Gates::list.end()) return 0;

	return Gates::list[gateName];
}

bool Gates::addGate(string name, float complex* matrix){
	if (Gates::list.find(name) != Gates::list.end()) return false;

	Gates::list[name] = matrix;
	return true;
}

bool Gates::addGate(string name, float complex a0, float complex a1, float complex a2, float complex a3){
//        cout << "GATES" << endl;

	if (Gates::list.find(name) != Gates::list.end()) return false;

//        cout << "GATES" << endl;

	float complex* matrix = new float complex[4];
	matrix[0] = a0;
	matrix[1] = a1;
	matrix[2] = a2;
	matrix[3] = a3;

	/*
	cout << "Gate Adicionado " << name << endl;
	cout << "Matrix" << endl;
	cout << crealf(matrix[0]) << " , " << cimagf(matrix[0]) << "\t" << crealf(matrix[1]) << " , " << cimagf(matrix[1]) << endl;
	cout << crealf(matrix[2]) << " , " << cimagf(matrix[2]) << "\t" << crealf(matrix[3]) << " , " << cimagf(matrix[3]) << endl << endl;
	*/

	Gates::list[name] = matrix;
	return true;
}

/////////////////////////////////////////////////////////////////////

string CNot(int qubits, int ctrl, int target, int cv){
	vector <string> cn (qubits, "ID");
	cn[ctrl] = "Control1(0)";
	if (cv) cn[ctrl] = "Control1(1)";
	cn[target] = "Target1(X)";

	return concatena(cn, qubits);
}

string Toffoli(int qubits, int ctrl1, int ctrl2, int target, int cv){
	vector <string> tf (qubits, "ID");
	tf[ctrl1] = "Control1(0)";
	if (cv>>1) tf[ctrl1] = "Control1(1)";
	tf[ctrl2] = "Control1(0)";
	if (cv&1) tf[ctrl2] = "Control1(1)";
	tf[target] = "Target1(X)";

	return concatena(tf, qubits);
}

string Controlled1(int qubits, int ctrl, int target, string op, int cv){
	vector <string> c1 (qubits, "ID");
	c1[ctrl] = "Control1(0)";
	if (cv) c1[ctrl] = "Control1(1)";
	c1[target] = "Target1(" + op + ")";

	return concatena(c1, qubits);
}


string Controlled2(int qubits, int ctrl1, int ctrl2, int target, string op, int cv){
	vector <string> c2 (qubits, "ID");
	c2[ctrl1] = "Control1(0)";
	if (cv&2) c2[ctrl1] = "Control1(1)";
	c2[ctrl2] = "Control1(0)";
	if (cv&1) c2[ctrl2] = "Control1(1)";
	c2[target] = "Target1(" + op + ")";

	return concatena(c2, qubits);
}

string ControlledN(int qubits, vector <int> ctrls, int target, string op, int cv){
	if (cv == -1) cv = pow(2,ctrls.size()) - 1;
	vector <string> c (qubits, "ID");

	for (int i = ctrls.size() - 1; i >=0; i--){
		c[ctrls[i]] = "Control1(0)";
		if (cv & 1) c[ctrls[i]] = "Control1(1)";
		cv = cv >> 1;
	}

	c[target] = "Target1(" + op + ")";

	return concatena(c, qubits);
}


string Pauli_X(int qubits, int reg, int width){
	vector <string> px (qubits, "ID");
	for (int i = 0; i < width; i++) px[i+reg] = "X";

	return concatena(px, qubits);
}

string Pauli_Z(int qubits, int reg, int width){
	vector <string> px (qubits, "ID");
	for (int i = 0; i < width; i++) px[i+reg] = "Z";

	return concatena(px, qubits);
}


string Hadamard(int qubits, int reg, int width){
	vector <string> hn (qubits, "ID");
	for (int i = 0; i < width; i++) hn[i+reg] = "H";

	return concatena(hn, qubits);
}


string concatena(vector <string> vec, int size, bool rev){
	string s;
	if (!rev){
		s = vec[0];
		for (int i = 1; i < size; i++)
			s += "," + vec[i];
	}
	else{
		s = vec[size-1];
		for (int i = size - 2; i >= 0; i--)
			s += "," + vec[i];
	}

	return s;
}
