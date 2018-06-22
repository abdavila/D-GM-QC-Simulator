#include "genMem.h"
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <fstream>

using namespace std;

void tensorProduct(vector <qubit> &list, float complex *state, int pos, int num, float complex value){
	if (num == list.size()){
		state[pos] = value;
		return;
	}

	if (list[num].zero != 0) tensorProduct(list, state, pos << 1, num+1, (value * list[num].zero));
	if (list[num].one != 0) tensorProduct(list, state, (pos << 1) + 1, num+1, (value * list[num].one));
}

float complex* genMem(vector <qubit> entrada){
	int mem_size = pow(2,entrada.size());
	
	float complex *state = (float complex *) calloc(mem_size, sizeof(float)*2);
	
	tensorProduct(entrada, state, 0, 0, 1.0);

	return state;
}

void saveAsXML(int qubits, float complex *state, string name){
	string arq;
	int mem_size = pow(2,qubits);

	ofstream out(name.c_str());

	out << "<memoria>" << endl;
	out << "<posicao dimensao=\"1\" tamanho=\""<< qubits << "\"></posicao>" << endl;	  
	
	out << "<valores>Complex";
	for (int i = 1; i < qubits; i++){out << ",Complex";};
	out << "</valores>" << endl;

	out << "<dados>";

	out << "(" << crealf(state[0]);
	if (cimagf(state[0]) >= 0) out << "+";
	out << cimagf(state[0]) << "j)";

	for (int i = 1; i < mem_size; i++){
		out << ", (" << crealf(state[i]);
		if (cimagf(state[i]) >= 0) out << "+";
		out << cimagf(state[i]) << "j)";
	}

	out << "</dados>\n</memoria>" << endl;
	out.close();
}



