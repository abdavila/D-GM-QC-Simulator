#ifndef _GEN_MEM_H_
#define _GEN_MEM_H_

#include <iostream>
#include <vector>
#include <complex.h>
#include <string>
#include <math.h>
#include "common.h"

using namespace std;

struct qubit{
	float complex zero;
	float complex one;
	qubit(float complex z, float complex o){zero = z; one = o;};
	qubit(float complex o){one = o; zero = cpowf(1-cpowf(o,2.0), 0.5);};
	void print(){cout << crealf(zero) << " " << crealf(one) << endl;}
};

static qubit ZERO(1,0);
static qubit ONE(0,1);

void tensorProduct(vector <qubit> &list, float complex *state, int pos, int num, float complex value);
float complex* genMem(vector <qubit> entrada);

void saveAsXML(int qubits, float complex *state, string name);


#endif

