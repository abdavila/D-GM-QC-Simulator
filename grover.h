#ifndef _GEN_MEM_H_
#define _GEN_MEM_H_

#include "common.h"
#include "gates.h"

using namespace std;


float HadamardNQB(long qubits, long num_of_it, int qb_pos, int multi_gpu, int qbs_region, int coalesc, int tam_block, int rept);



string ControledZ(int qubits);
string Oracle1(long qubits, long int value);

vector<string> GroverOrcl1(long qubits, long value){


	
}
