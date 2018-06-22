#ifndef _DGM_H_
#define _DGM_H_

#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
//#include <complex.h>
#include <sys/time.h>
//#include <pthread.h>
#include "common.h"
#include "gates.h"
//#include "genMem.h"

using namespace std;

long* SpectroExecute(int qubits, vector<string> function);

float complex* GenericExecute(float complex *state, string function, int qubits, int type, int threads, int factor);
float complex*  GenericExecute(float complex *state, vector<string> function, int qubits, int type, int threads, int factor);

extern "C" bool setDevice(int num = 0);

//extern "C" float complex* GpuExecution0(float complex* r_memory, float complex* w_memory, PT **pts, int qubits, long MAX_PT, long MAX_QB, int it);
extern "C" float complex* GpuExecutionWrapper(float complex* r_memory, PT **pts, int qubits, int multi_gpu, int qbs_region, int tam_block, int rept, int coalesc, int num_it);
extern "C" float complex* GpuExecution(float complex* r_memory, float complex* w_memory, PT **pts, int qubits, float *total_time, long MAX_PT, long MAX_QB, int it);
extern "C" float complex* GpuExecution2(float complex* r_memory, PT **pts, int pts_size, int qubits, long MAX_PT, int it);
extern "C" float complex* GpuExecution3(float complex* r_memory, float complex* w_memory, int sub_size, int shift_write, PT *pt, int qubits, long MAX_PT, long MAX_QB, int it);
//extern "C" float complex* PCpuExecution1(float complex* r_memory, float complex* w_memory, PT **pts, int qubits, long n_threads);
//extern "C" float complex* PCpuExecution2(float complex* r_memory, float complex* w_memory, PT **pts, int qubits, long n_threads);
float complex* PCpuExecution1(float complex* r_memory, float complex* w_memory, PT **pts, int qubits, long n_threads, int it);
float complex* PCpuExecution3(float complex* r_memory, float complex* w_memory, PT **pts, int qubits, long n_threads);
void PCpuExecution1_0(float complex* r_memory, float complex* w_memory, PT *pt, long mem_size, long n_threads, long desloc = 0);
void PCpuExecution1_1(float complex* r_memory, float complex* w_memory, PT *pt, long mem_size, long n_threads, long desloc = 0);
void PCpuExecution1_2(float complex* r_memory, float complex* w_memory, PT *pt, long mem_size, long n_threads, long desloc = 0);
void PCpuExecution1_3(float complex* r_memory, float complex* w_memory, PT *pt, long mem_size, long n_threads, long desloc = 0);

void PCpuExecution2_1(float complex* r_memory, float complex* w_memory, PT *pt, long mem_size, long n_threads, long desloc = 0);
void PCpuExecution2_11(float complex* r_memory, float complex* w_memory, PT *pt, long mem_size, long n_threads, long desloc = 0);
void PCpuExecution3_1(float complex* r_memory, float complex* w_memory, PT *pt, long mem_size, long n_threads, long desloc = 0);

void PCpuExecution4_0(float complex *r_memory, PT *pt, long mem_size, long n_threads, long desloc = 0);
void PCpuExecution4_1(float complex* r_memory, PT *pt, long mem_size, long n_threads, long desloc = 0);
void PCpuExecution4_2(float complex* r_memory, PT *pt, long mem_size, long n_threads, long desloc = 0);
void PCpuExecution4_3(float complex* r_memory, PT *pt, long mem_size, long n_threads, long desloc = 0);

inline long LINE (long pos, long shift){
	return ((pos >> shift) & 1) * 2;
}
inline long BASE (long pos, long shift){
	return pos & (~(1 << shift));
}

//void *cpufunction(void *param);

enum {
	t_CPU,
	t_PAR_CPU,
	t_GPU,
	t_HYBRID,
	t_SPEC
};

/*
struct param{
	int id;
	long start;
	long end;
	PT **pts;

	float complex *r_mem;
	float complex *w_mem;

	pthread_barrier_t *barrier;
};
*/

class Group{
public:
	vector <string> ops;
	vector <long> pos_ops;
	vector <bool> ctrl;
	vector <long> pos_ctrl;

	Group(){};
	bool isAfected(int pos, int afect);
};

class DGM{
public:
	vector <string> diag;
	long MAX_QB, MAX_PT, qb_afected;
	long n_threads;
	long factor;

	int exec_type;

	int multi_gpu;
	int qbs_region;
	int tam_block;
	int rept;
	int coalesc;

	long *spectro;
	int spec_block;
	long spec_region;


	//string function;
	vector <PT*> vec_pts;
	PT** pts;
	long qubits;

	float measure_value;

	float elapsed_time;
	struct timeval timev;

	float complex *r_mem;
	float complex *w_mem;

	DGM();
	~DGM();

	bool en_print;

	void printPTs();
	void erase();
	void setExecType(int type);

	void allocateMemory();
	void setMemory(float complex *mem);
	void freeMemory();
	void setMemoryValue(int pos);
	void allocateSpectro();
	int measure(int q_pos);
	int measure_coalesced(int k);
	void measure(vector<int> q_pos);
	void colapse(int q_pos, int value);

	void setFunction(string function, int it = 1, bool er = true);
	void setFunction(vector<string> steps, int it = 1, bool er = true);
	map <long, Group> genGroups(string step);
	void genPTs(map<long, Group> &gps, vector <PT*> &step_pts);
	//bool isAfected(Group &g, long pos);
	void genMatrix(float complex* matrix, vector<float complex*> &matrices, long tam, long current, long line, long column, float complex cmplx);

	//float complex getOpValue(string op, long l, long c);
	/*
	setMemory();
	setQbLimit();
	setSubSize();
	*/
	void executeFunction(string function, int it = 1);
	void executeFunction(vector<string> steps, int it = 1);
	float complex* execute(int it);
	float complex* HybridExecution(float complex *r_memory, float complex *w_memory, int it);
	float complex* CpuExecution1(float complex *r_memory, float complex *w_memory, int it);
	void CpuExecution1_1(float complex *r_memory, float complex *w_memory, PT *pt, long mem_size);
	void CpuExecution1_2(float complex *r_memory, float complex *w_memory, PT *pt, long mem_size);
	void CpuExecution1_3(float complex *r_memory, float complex *w_memory, PT *pt, long mem_size);

	void CpuExecution2_1(float complex *r_memory, float complex *w_memory, PT *pt, long mem_size);

	void CpuExecution_spectro(int it);
	void CpuExecution1_spectro_1(PT *pt, long mem_size);
	void CpuExecution1_spectro_2(PT *pt, long mem_size);
	void CpuExecution1_spectro_3(PT *pt, long mem_size);
	//float complex* CpuExecution1_2(float complex *r_memory, float complex *w_memory);
	//float complex* CpuExecution2(float complex *r_memory, float complex *w_memory);
	//float complex* CpuExecution2_1(float complex *r_memory, float complex *w_memory);
};

#endif
