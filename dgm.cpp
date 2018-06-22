#include <iostream>
#include "dgm.h"
#include <omp.h>
#include <unistd.h>
#include <cstdio>
#include <iterator>

#define CHUNK 2048

void Tokenize(const string& str, vector<string>& tokens, const string& delimiters = ",")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////

long* SpectroExecute(int qubits, vector<string> function){
	DGM dgm;
	dgm.exec_type = t_SPEC;
	dgm.qubits = qubits;

	dgm.allocateSpectro();

	dgm.executeFunction(function);

	return dgm.spectro;
}


///////////////////////////////////////////////////////////////////////////////////////////////

float complex* GenericExecute(float complex *state, string function, int qubits, int type, int threads, int factor = 0){
	DGM dgm;
	dgm.exec_type = type;
	dgm.n_threads = threads;
	dgm.qubits = qubits;
	dgm.factor = factor;

	dgm.setMemory(state);

	dgm.executeFunction(function);

	state = dgm.r_mem;

	dgm.r_mem = NULL;

	free(dgm.w_mem);

	return state;
}

float complex* GenericExecute(float complex *state, vector<string> function, int qubits, int type, int threads, int factor = 0){
	DGM dgm;
	dgm.exec_type = type;
	dgm.n_threads = threads;
	dgm.qubits = qubits;
	dgm.factor = factor;
	dgm.setMemory(state);

	//struct timeval t, tBegin, tEnd;
	//float te;

	//gettimeofday(&tBegin, NULL);

	dgm.executeFunction(function);

	printMem(dgm.r_mem, 4);

	//gettimeofday(&tEnd, NULL);
	//timeval_subtract(&t, &tEnd, &tBegin);
	//te = t.tv_sec + (t.tv_usec / 1000000.0);

	//cout << "Time: " << te << endl;

	//state = dgm.w_mem;
	//dgm.r_mem = NULL;

	//printMem(state, 4);
	//printMem(dgm.w_mem, 4);

	free(dgm.w_mem);
	free(dgm.r_mem);


	return state;
}

///////////////////////////////////////////////////////////////////////////////////////////////

DGM::DGM(){
	MAX_QB = QB_LIMIT;
	MAX_PT = PT_TAM;

	pts = NULL;
	r_mem = w_mem = NULL;
	en_print = false;
	exec_type = t_CPU;
	factor = 1;
	multi_gpu = 1;

	spec_block = 1;
}

DGM::~DGM(){erase();}

void DGM::setExecType(int type){
	exec_type = type;
}

void DGM::printPTs(){
	for (int i = 0; i < vec_pts.size() -1; i++){
		vec_pts[i]->print();
	}
}

void DGM::erase(){
    if (!pts) return;

    long i = 0;
    while (pts[i] != NULL){
    	pts[i]->destructor();
    	free(pts[i]);
    	i++;
    }

    vec_pts.clear();
    pts = NULL;
}

void DGM::allocateMemory(){
	r_mem = (float complex*) calloc(pow(2, qubits), sizeof(float complex));
	w_mem = (float complex*) calloc(pow(2, qubits), sizeof(float complex));
}

void DGM::setMemory(float complex* mem){
	freeMemory();
	r_mem = mem;
	w_mem = (float complex*) calloc(pow(2, qubits), sizeof(float complex));
}

void DGM::freeMemory(){
	if (r_mem) free(r_mem);
	if (w_mem) free(w_mem);
	r_mem = w_mem = NULL;
}

void DGM::setMemoryValue(int pos){
	r_mem[pos] = 1;
}

void DGM::allocateSpectro(){
	spectro = (long *) calloc(spec_block*spec_block, sizeof(long));
	spec_region = (pow(2,qubits))/spec_block;
}

int DGM::measure(int q_pos){
	long size = pow(2.0, qubits);

	long shift = (qubits - 1 - q_pos);

	int count_one, count_zero, num_pb;
	float zero, one, norm, r;
	one = zero = 0;

	omp_set_num_threads(n_threads);	
	#pragma omp parallel for reduction(+:zero,one)
	for (long i = 0; i < size; i++){
        if ((i >> shift)&1) one += pow(crealf(r_mem[i]), 2.0) + pow(cimagf(r_mem[i]), 2.0);
		else zero += pow(crealf(r_mem[i]), 2.0) + pow(cimagf(r_mem[i]), 2.0);
	}

	printf("%f %f\n", zero, one);
	

	long m;
	/*
	srand (time(NULL));
	count_one = 0;
	count_zero = 0;	
	num_pb = 1;

	for (int i = 0; i < num_pb; i++){
		r = (double) rand() / RAND_MAX;
		if (zero > r) count_zero++;
		else count_one++;
	}

	if (count_one > count_zero){
		measure_value = one;
		norm = sqrt(one);
		m = 1;
	}
	else{
		measure_value = zero;
		norm = sqrt(zero);
		m = 0;
	}

	long mask;
	mask = pow(2,shift) - 1;
	#pragma omp for
	for (long i = 0; i < size/2; i++){
		long pos = (i << 1) - (i&mask);
		r_mem[pos] = r_mem[pos | (m << shift)]/norm;
		r_mem[pos | (1<<shift)] = 0.0;
	}
	*/

	return m;
}

int DGM::measure_coalesced(int k){
	long size = pow(2.0, qubits);

	long gap = pow(2, (k+1));
	long stride = pow(2, k);

	float zero, one;
	zero = one = 0;

	omp_set_num_threads(n_threads);
	#pragma omp parallel for reduction(+:zero,one)
	for (long i = 0; i < size; i+=gap){
		for (long j = i; j < (i+stride); j++){
			zero += pow(crealf(r_mem[j]), 2.0) + pow(cimagf(r_mem[j]), 2.0);
			one += pow(crealf(r_mem[j+stride]), 2.0) + pow(cimagf(r_mem[j+stride]), 2.0);
		}
	}

	return 0;
}




void DGM::colapse(int q_pos, int value){
	long size = pow(2.0, qubits);
	long mask = (qubits - 1 - q_pos);

	float m, norm, r;
	m = 0;

	for (long i = 0; i < size; i++)
        if (((i >> mask)&1) == value) m += pow(crealf(r_mem[i]), 2.0) + pow(cimagf(r_mem[i]), 2.0);

	cout << m << endl;

	m = sqrt(m);
	for (long i = 0; i < size; i++){
		if (((i >> mask)&1) == value) r_mem[i] = r_mem[i]/m;
		else r_mem[i] = 0.0;
	}
}

void DGM::measure(vector<int> q_pos){
	long mask = 0;

	for (int i =0; i < q_pos.size(); i++) mask = mask | (1<<(qubits - 1 - q_pos[i]));

	map <long, float> m;

	long size = pow(2.0, qubits);

	for (long i =0; i < size; i++) m[i&mask] += pow(crealf(r_mem[i]), 2.0) + pow(cimagf(r_mem[i]), 2.0);

	for (map<long,float>::iterator it=m.begin(); it!=m.end(); ++it)
		cout << it->first << " =>  " << it->second << '\n';

	return;

	float norm = sqrt(m[0]);
	long pos;
	long shift = (qubits - 1 - q_pos[0]);
	mask = pow(2,shift) - 1;
	for (long i = 0; i < size/2; i++){
		pos = (i << 1) - (i&mask);
		r_mem[pos] = r_mem[pos]/norm;
		r_mem[pos | (1<<shift)] = 0.0;
	}

}

void DGM::setFunction(string function, int it, bool er){
	//if (er) erase();
	vector <string> steps;

	Tokenize(function, steps, ";");

	setFunction(steps, it, er);
}

void DGM::setFunction(vector <string> steps, int it, bool er){
	if (er) erase();
	else vec_pts.pop_back();


	vector <PT*> step_pts, vec_tmp;
	map<long, Group> gps;

	for (long j = 0; j< it; j++)
	for (long i = 0; i < steps.size(); i++){
		gps = genGroups(steps[i]);
		genPTs(gps, step_pts);

		//cout << "STEP: " << i << endl;

		if (i%2)
			sort(step_pts.begin(), step_pts.end(), increasing);
		else
			sort(step_pts.begin(), step_pts.end(), decreasing);

		//for (long j = 0; j < step_pts.size(); j++)
		//	cout << step_pts[j]->start << endl;

		vec_pts.insert(vec_pts.end(), step_pts.begin(), step_pts.end());
	}

	
		//vec_pts.insert(vec_pts.end(), vec_tmp.begin(), vec_tmp.end());


	//cout << "Vec Size 2: " << vec_pts.size() << endl;


	vec_pts.push_back(NULL);

	pts = &vec_pts[0];
}

map <long, Group> DGM::genGroups(string step){
	vector <string> ops;
	Tokenize(step, ops); //separa os operadores usando "," como delimitador
	qubits = ops.size();

	size_t found_c, found_t, p;
	string str;
	long pos, ctrl_value, ctrl_num;
	
	map<long, Group> gps;

	char * pEnd;
	pos = 0;
	vector<string>::iterator it;
	for (it = ops.begin() ; it != ops.end(); ++it){ //percorre os operadores
		str = *it;
		//cout << str << endl;
		found_c = str.find("Control"); //tamanho 7
		found_t = str.find("Target");  //tamanho 6
		p = str.find("(") + 1;

		if (found_c != string::npos){ //Controle
			ctrl_num = strtol(str.c_str()+7, &pEnd, 10);
			ctrl_value = strtol(str.c_str()+p, &pEnd, 10);

			gps[ctrl_num].ctrl.push_back(ctrl_value); //adicona o valor do controle
			gps[ctrl_num].pos_ctrl.push_back(pos);  //e a sua posição ao map relacionado ao controle
		}
		else if(found_t != string::npos){ //Target
			ctrl_num = strtol(str.c_str()+6, &pEnd, 10);
			str = str.substr(p, str.size()-p-1);

			gps[ctrl_num].ops.push_back(str);     //adicona o operador
			gps[ctrl_num].pos_ops.push_back(pos); //e a sua posição ao map relacionado ao target
		}
		else{ //operador normal
			if (str != "ID"){ //se for ID ignora
				gps[0].ops.push_back(str);     //adiciona o operador
				gps[0].pos_ops.push_back(pos); //e a sua posição ao map '0'
			}
		}
		pos++;
	}
	
	return gps;
}

void DGM::genPTs(map<long, Group> &gps, vector <PT*> &step_pts){
	step_pts.clear();
	Gates gates;

	map<long,Group>::iterator it;	
	Group gp;
	PT* pt;
	long ctrl_mask, ctrl_value, ctrl_count;
	long start, tam, size;
	
	for (it = gps.begin(); it != gps.end(); ++it){ //percorre os grupos
		gp = it->second;
		size = gp.ops.size();
		
		ctrl_count = gp.ctrl.size();
		ctrl_value = ctrl_mask = 0;

		for (long i = 0; i < ctrl_count; i++){ //gera a mascara e o valor do controle (em binário)
			gp.pos_ctrl[i] =  qubits - gp.pos_ctrl[i] - 1;
			ctrl_mask += (1 << gp.pos_ctrl[i]);
			if (gp.ctrl[i]) ctrl_value += (1 << gp.pos_ctrl[i]);
		}

		for (int p = 0; p < size; p++){
			
			pt = (PT*) malloc(sizeof(PT));
			pt->affected = false;

			pt->qubits = 1;
			pt->start = qubits - gp.pos_ops[p];
			pt->end = pt->start - 1;
			pt->mat_size = 2;
			
			pt->matrix = gates.getMatrix(gp.ops[p]);

			pt->ctrl_value = ctrl_value;
			pt->ctrl_mask = ctrl_mask;
			pt->ctrl_count = ctrl_count;

			if (ctrl_count){
				pt->ctrl_pos = (long*)malloc(sizeof(long) * ctrl_count);
				copy(gp.pos_ctrl.begin(), gp.pos_ctrl.end(), pt->ctrl_pos);
			}

			step_pts.push_back(pt);
		}
	}
}

/*
bool DGM::isAfected(Group &gp, long pos){
	string str = gp.ops[pos];
	size_t p = str.find("(");		
	str = str.substr(0, p);

	return (gp.pos_ops[pos] < qb_afected) && (find(diag.begin(), diag.end(), str) == diag.end());
}
*/

void DGM::genMatrix(float complex* matrix, vector<float complex*> &matrices, long tam, long current, long line, long column, float complex cmplx){
	if (cmplx == 0.0) return;

	if (current == tam){ //percorreu até a ultima matriz
		matrix[line*(1<<tam) + column] = cmplx;
		return;
	}

	for (long l = 0; l < 2; l++)
		for (long c = 0; c < 2; c++)
			genMatrix(matrix, matrices, tam, current+1, (line<<1)|l, (column<<1)|c, cmplx * matrices[current][l*2+c]);
}


void DGM::executeFunction(vector <string> function, int it){
	setFunction(function);
	execute(it);
}

void DGM::executeFunction(string function, int it){
	if (function == "") return;

	setFunction(function);
	execute(it);
}

float complex* DGM::execute(int it){
    //float *s_time = new float[1];
    float complex* result;

	switch (exec_type){
		case t_HYBRID:
			result = HybridExecution(r_mem, w_mem, it);
			break;
		case t_CPU:
			result = CpuExecution1(r_mem, w_mem, it);
			break;
		case t_SPEC:
			CpuExecution_spectro(it);
			break;
		case t_PAR_CPU:
			result = PCpuExecution1(r_mem, w_mem, pts, qubits, n_threads, it);
			break;
		case t_GPU:
	    	result = GpuExecutionWrapper(r_mem, pts, qubits, qbs_region, multi_gpu, tam_block, rept, coalesc, it);
			break;
		default:
			exit(1);
	}
	    
    if (result == w_mem) swap(r_mem, w_mem);
	
    //delete s_time;
    return result;

}

float complex* DGM::HybridExecution(float complex *r_memory, float complex *w_memory, int it){
	long mem_size, sub_size, sub_write, sub_qubits, num_parts, i, count;
	bool aff;
	float complex *r_cpu, *w_cpu, *aux;

	mem_size = pow(2.0, qubits);
	num_parts = pow(2.0, ceil(log2(n_threads + 1)))*factor;
	sub_size = mem_size/num_parts;
	sub_qubits = log2(sub_size);

	omp_set_nested(1);
	//setDevice();

	cout << "AQUI " << num_parts << " " << sub_qubits << " " << sub_size << endl;

	//omp_lock_t writelock;
	//omp_init_lock(&writelock);

	//omp_set_num_threads(4);

	int pause;

	struct timeval tcpu, tgpu, tcBegin, tcEnd, tgBegin, tgEnd;
	float tc, tg;
	gettimeofday(&tcBegin, NULL);
	
	
	i = 0;

	while (pts[i] != NULL){
		sub_write = -1;
		count = 1;

		if (pts[i]->start <= sub_qubits)
			while ((pts[i+count] != NULL) && (pts[i+count]->start <= sub_qubits)){
				count++;
			}
		
		//cout << "\nCount: " << count << "  ";

		#pragma omp parallel sections num_threads(2)
		{
			#pragma omp section 					//CPU EXECUTION
			{
				long cpu_write;
				
				#pragma omp critical
				{
					sub_write++;
					cpu_write = sub_write;
				}
				
				while (cpu_write < num_parts){
					long desloc = sub_size*cpu_write;
					r_cpu = r_memory;
					w_cpu = w_memory;

					if (count == 1){
						PCpuExecution1_0(r_memory, w_memory+desloc, pts[i], sub_size, n_threads, desloc);
					}
					else{
						for (int p = 0; p < count; p++){
							PCpuExecution4_0(r_memory+desloc, pts[i+p], sub_size, n_threads, desloc);
						}
					}
					//printf("CPU: %ld\n", cpu_write);

					#pragma omp critical
					{
						sub_write++;
						cpu_write = sub_write;
					}
				}
			}
			#pragma omp section 					//GPU EXECUTION
			{
				/*
				long gpu_write;

				#pragma omp critical
				{
					sub_write++;
					gpu_write = sub_write;
				}
				
				while (gpu_write < num_parts){
					long desloc = sub_size*gpu_write;

					//cout << "\nPTS: " << i << " \nDesloc: " << desloc << endl;

					//gettimeofday(&tgBegin, NULL);
					if (pts[i]->start <= sub_qubits){
						GpuExecution2(r_memory+desloc, pts+i, count, sub_qubits, MAX_PT, 1);
					}
					else{
						GpuExecution3(r_memory, w_memory+desloc, sub_size, desloc, pts[i], qubits, MAX_PT, MAX_QB, 1);
					}
					//gettimeofday(&tgEnd, NULL);
					//timeval_subtract(&tgpu, &tgEnd, &tgBegin);
					//tg = tgpu.tv_sec + (tgpu.tv_usec / 1000000.0);
					//printf("GPU: %ld\n", gpu_write);

					//printf("GPU: %ld\n", gpu_write);
					
			    	#pragma omp critical
					{
						sub_write++;
						gpu_write = sub_write;
					}
				}
				*/
			}
		}
		if (count == 1) swap_ptr(&r_memory, &w_memory);
		i+=count;
	}

	//printMem(r_memory, 4);

	return r_memory;
}


float complex* DGM::CpuExecution1(float complex *r_memory, float complex *w_memory, int it){
	long mem_size = pow(2.0, qubits);

    for (int x = 0; x < it; x++){
	    long i = 0;
	    while (pts[i] != NULL){
	    	//printMem(r_memory, 4);
	    	//char c;
	    	//cin >> c;
	    	//getchar();
	    	long mt = pts[i]->matrixType();

	    	switch (mt){
	    		case DENSE:
	    			//ALTERA AQUI
	    			//CpuExecution1_1(r_memory, w_memory, pts[i], mem_size);
	    			CpuExecution2_1(r_memory, w_memory, pts[i], mem_size);
	    			break;
	    		case DIAG_PRI:
	    			CpuExecution1_2(r_memory, w_memory, pts[i], mem_size);
	    			break;
	    		case DIAG_SEC:
	    			CpuExecution1_3(r_memory, w_memory, pts[i], mem_size);
	    			break;
	    		default:
	    			exit(1);
	    	}
		    i++;
		    swap_ptr(&r_memory, &w_memory);
	    }
    }
    //printMem(r_memory, 4);
	
	swap_ptr(&r_memory, &w_memory);
	
	return w_memory;
}


void DGM::CpuExecution1_1(float complex *r_memory, float complex *w_memory, PT *pt, long mem_size){ //Denso
	long base, shift, l;
	
	shift = pt->end;

	//cout << "###########################################" << endl;
		
	if (!pt->ctrl_count) 	//operador não controlado
		for (long pos = 0; pos < mem_size; pos++){
			l = ((pos >> shift) & 1) * 2;
			base = pos & (~(1 << shift));

			w_memory[pos] = pt->matrix[l] * r_memory[base] + 
							pt->matrix[l+1] * r_memory[base | (1 << shift)];
							
		}
	else					//operador controlado
		for (long pos = 0; pos < mem_size; pos++){
			if ((pos & pt->ctrl_mask) == pt->ctrl_value){
				l = ((pos >> shift) & 1) * 2;
				base = pos & (~(1 << shift));
	
				w_memory[pos] = pt->matrix[l] * r_memory[base] + 
								pt->matrix[l+1] * r_memory[base | (1 << shift)];
			}
			else
				w_memory[pos] = r_memory[pos];
		}

}

void DGM::CpuExecution1_2(float complex *r_memory, float complex *w_memory, PT *pt, long mem_size){ //Diagonal Principal
	long shift = pt->end;
		
	if (!pt->ctrl_count)	//operador não controlado
		for (long pos = 0; pos < mem_size; pos++)
			w_memory[pos] = pt->matrix[((pos >> shift) & 1) * 3] * r_memory[pos];
	else					//operador controlado
		for (long pos = 0; pos < mem_size; pos++)
			if ((pos & pt->ctrl_mask) == pt->ctrl_value)
				w_memory[pos] = pt->matrix[((pos >> shift) & 1) * 3] * r_memory[pos];
			else
				w_memory[pos] = r_memory[pos];

}

void DGM::CpuExecution1_3(float complex *r_memory, float complex *w_memory, PT *pt, long mem_size){ //Diagonal Secundária
	long l, shift = pt->end;
		
	if (!pt->ctrl_count)	//operador não controlado
		for (long pos = 0; pos < mem_size; pos++)
			w_memory[pos] = pt->matrix[((pos >> shift) & 1) + 1] * r_memory[pos^(1<<shift)];
	else					//operador controlado
		for (long pos = 0; pos < mem_size; pos++)	
			if ((pos & pt->ctrl_mask) == pt->ctrl_value)
				w_memory[pos] = pt->matrix[((pos >> shift) & 1) + 1] * r_memory[pos^(1<<shift)];
			else
				w_memory[pos] = r_memory[pos];
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void DGM::CpuExecution_spectro(int it){
	long mem_size = pow(2.0, qubits);

    for (int x = 0; x < it; x++){
	    long i = 0;
	    while (pts[i] != NULL){
	    	long mt = pts[i]->matrixType();

	    	switch (mt){
	    		case DENSE:
	    			//ALTERA AQUI
	    			//CpuExecution1_1(r_memory, w_memory, pts[i], mem_size);
	    			CpuExecution1_spectro_1(pts[i], mem_size);
	    			break;
	    		case DIAG_PRI:
	    			CpuExecution1_spectro_2(pts[i], mem_size);
	    			break;
	    		case DIAG_SEC:
	    			CpuExecution1_spectro_3(pts[i], mem_size);
	    			break;
	    		default:
	    			exit(1);
	    	}
		    i++;
	    }
    }
    //printMem(r_memory, 4);
}

void DGM::CpuExecution1_spectro_1(PT *pt, long mem_size){ //Diagonal Principal
	long pos0, pos1, base, shift, l;
	
	shift = 1 << pt->end;
	long range = mem_size / 2;

	for (long pos = 0; pos < range; pos++){
		pos0 = (pos * 2) - (pos & (shift-1));
		pos1 = pos0 | shift;

		pos0/=spec_region;
		pos1/=spec_region;

		if ((!pt->ctrl_count) || ((pos0 & pt->ctrl_mask) == pt->ctrl_value)){
			spectro[pos0*spec_block + pos0]++;
			spectro[pos0*spec_block + pos1]++;
			spectro[pos1*spec_block + pos0]++;
			spectro[pos1*spec_block + pos1]++;
		}
	}
}


void DGM::CpuExecution1_spectro_2(PT *pt, long mem_size){ //Diagonal Principal
	long shift = pt->end;

	for(long pos = 0; pos < mem_size; pos++)
		if (!pt->ctrl_count || ((pos & pt->ctrl_mask) == pt->ctrl_value)){
			long p = pos/spec_region; 
			spectro[p*spec_block + p]++;
		}
}

void DGM::CpuExecution1_spectro_3(PT *pt, long mem_size){ //Diagonal Secundária
	long l, shift = pt->end;

	for(long pos = 0; pos < mem_size; pos++)
		if (!pt->ctrl_count || ((pos & pt->ctrl_mask) == pt->ctrl_value)){
			long p0 = pos/spec_region;
			long p1 = (pos^(1<<shift))/spec_region;

			spectro[p0*spec_block + p1]++;
		}
}


/*
float complex* DGM::CpuExecution2(float complex *r_memory, float complex *w_memory, int it){
	long mem_size = pow(2.0, qubits);

    for (int x = 0; x < it; x++){
	    long i = 0;
	    while (pts[i] != NULL){
	    	long mt = pts[i]->matrixType();

	    	switch (mt){
	    		case DENSE:
	    			CpuExecution2_1(r_memory, w_memory, pts[i], mem_size);
	    			break;
	    		case DIAG_PRI:
	    			CpuExecution2_2(r_memory, w_memory, pts[i], mem_size);
	    			break;
	    		case DIAG_SEC:
	    			CpuExecution2_3(r_memory, w_memory, pts[i], mem_size);
	    			break;
	    		default:
	    			exit(1);
	    	}
		    i++;
		    swap_ptr(&r_memory, &w_memory);
	    }
    }
	
	swap_ptr(&r_memory, &w_memory);
	
	return w_memory;
}
*/

void DGM::CpuExecution2_1(float complex *r_memory, float complex *w_memory, PT *pt, long mem_size){ //Diagonal Principal
	long pos0, pos1, base, shift, l;
	
	shift = 1 << pt->end;
	mem_size /= 2;
		
	if (!pt->ctrl_count) 	//operador não controlado
		for (long pos = 0; pos < mem_size; pos++){
			pos0 = (pos * 2) - (pos & (shift-1));
			pos1 = pos0 | shift;

			w_memory[pos0] = pt->matrix[0] * r_memory[pos0] + pt->matrix[1] * r_memory[pos1];
			w_memory[pos1] = pt->matrix[2] * r_memory[pos0] + pt->matrix[3] * r_memory[pos1];
		}
	else					//operador controlado
		for (long pos = 0; pos < mem_size; pos++){
			pos0 = (pos * 2) - (pos & (shift-1));
			pos1 = pos0 | shift;
			if ((pos0 & pt->ctrl_mask) == pt->ctrl_value){
				w_memory[pos0] = pt->matrix[0] * r_memory[pos0] + pt->matrix[1] * r_memory[pos1];
				w_memory[pos1] = pt->matrix[2] * r_memory[pos0] + pt->matrix[3] * r_memory[pos1];			
			}
			else{
				w_memory[pos0] = r_memory[pos0];
				w_memory[pos1] = r_memory[pos1];
			}
		}

}


float complex* PCpuExecution1(float complex *r_memory, float complex *w_memory, PT **pts, int qubits, long n_threads, int it){
	long mem_size, i;
	
	mem_size = pow(2.0, qubits);

	//cout << pts[0]->end << endl;

    for (int x = 0; x < it; x++){
	    for (i = 0; pts[i] != NULL; i++){
	    	PCpuExecution1_0(r_memory, w_memory, pts[i], mem_size, n_threads);
		    
		    //swap_ptr(&r_memory, &w_memory);
	    }
    }
	swap_ptr(&r_memory, &w_memory);
	
	return w_memory;
}

void PCpuExecution1_0(float complex *r_memory, float complex *w_memory, PT *pt, long mem_size, long n_threads, long desloc){
    long mt = pt->matrixType();

	switch (mt){
		case DENSE:
			//ALTERA AQUI
			PCpuExecution2_11(r_memory, w_memory, pt, mem_size, n_threads, desloc);
			//PCpuExecution2_1(r_memory, w_memory, pt, mem_size, n_threads, desloc);
			break;
		case DIAG_PRI:
			PCpuExecution1_2(r_memory, w_memory, pt, mem_size, n_threads, desloc);
			break;
	   	case DIAG_SEC:
	   		PCpuExecution1_3(r_memory, w_memory, pt, mem_size, n_threads, desloc);
	   		break;
	   	default:
	   		exit(1);
	}
}

void PCpuExecution1_1(float complex *r_memory, float complex *w_memory, PT *pt, long mem_size, long n_threads, long desloc){ //Operador denso
	long base, shift, i, l, p0, p1;
	
	omp_set_num_threads(n_threads);
	
	i = 0;
	shift = pt->end;
		
	if (!pt->ctrl_count)
		#pragma omp parallel for private(l,base) //schedule(static, CHUNK)
		for (long pos = 0; pos < mem_size; pos++){
			l = (((pos+desloc) >> shift) & 1) * 2; //LINE((pos+desloc), shift); 
			base = (pos+desloc) & (~(1 << shift)); //BASE((pos+desloc), shift);
	
			w_memory[pos] = pt->matrix[l] * r_memory[base] + 
							pt->matrix[l+1] * r_memory[base | (1 << shift)];
		}
	else
		#pragma omp parallel for private(l,base) //schedule(static, CHUNK)
		for (long pos = 0; pos < mem_size; pos++){
			if (((pos+desloc) & pt->ctrl_mask) == pt->ctrl_value){
				l = LINE(pos+desloc, shift);
				base = BASE(pos+desloc, shift);
	
				w_memory[pos] = pt->matrix[l] * r_memory[base] + 
								pt->matrix[l+1] * r_memory[base | (1 << shift)];
			}
			else
				w_memory[pos] = r_memory[pos+desloc];
		}
}

void PCpuExecution1_2(float complex *r_memory, float complex *w_memory, PT *pt, long mem_size, long n_threads, long desloc){ //Diagonal Principal
	long shift, i;
	
	omp_set_num_threads(n_threads);
	
	shift = pt->end;
		
	if (!pt->ctrl_count)
		#pragma omp parallel for //schedule(static, CHUNK)
		for (long pos = 0; pos < mem_size; pos++)
			w_memory[pos] = pt->matrix[(((pos+desloc) >> shift) & 1) * 3] * r_memory[pos+desloc];
	else
		#pragma omp parallel for //schedule(static, CHUNK)
		for (long pos = 0; pos < mem_size; pos++)
			if (((pos+desloc) & pt->ctrl_mask) == pt->ctrl_value)
				w_memory[pos] = pt->matrix[(((pos+desloc) >> shift) & 1) * 3] * r_memory[pos+desloc];
			else
				w_memory[pos] = r_memory[pos+desloc];

}

void PCpuExecution1_3(float complex *r_memory, float complex *w_memory, PT *pt, long mem_size, long n_threads, long desloc){ //Diagonal Secundária
	omp_set_num_threads(n_threads);

	long shift = pt->end;
		
	if (!pt->ctrl_count)
		#pragma omp parallel for //schedule(static, CHUNK)
		for (long pos = 0; pos < mem_size; pos++)
			w_memory[pos] = pt->matrix[LINE(pos+desloc, shift)] * r_memory[(pos+desloc)^(1<<shift)];
	else
		#pragma omp parallel for //schedule(static, CHUNK)
		for (long pos = 0; pos < mem_size; pos++)
			if (((pos+desloc) & pt->ctrl_mask) == pt->ctrl_value)
				w_memory[pos] = pt->matrix[LINE(pos+desloc, shift)] * r_memory[(pos+desloc)^(1<<shift)];
			else
				w_memory[pos] = r_memory[pos+desloc];
}


void PCpuExecution2_1(float complex* r_memory, float complex* w_memory, PT *pt, long mem_size, long n_threads, long desloc){
	long pos, pos0, pos1, base, shift, i;
	float complex tmp;
	
	mem_size /= 2;

	omp_set_num_threads(n_threads);

	shift = 1 << pt->end;
			
	if (!pt->ctrl_count)
		#pragma omp parallel for private(pos0,pos1,tmp)
		for (long pos = 0; pos < mem_size; pos++){
			pos0 = (pos * 2) - (pos & (shift-1));
			pos1 = pos0 | shift;

			tmp = pt->matrix[0] * r_memory[pos0] + pt->matrix[1] * r_memory[pos1];
			r_memory[pos1] = pt->matrix[2] * r_memory[pos0] + pt->matrix[3] * r_memory[pos1];
			r_memory[pos0] = tmp;
		}
	else
		#pragma omp parallel for private(pos0,pos1,tmp)
		for (long pos = 0; pos < mem_size; pos++){
			pos0 = (pos * 2) - (pos & (shift-1));
			pos1 = pos0 | shift;
			if ((pos0 & pt->ctrl_mask) == pt->ctrl_value){
				tmp = pt->matrix[0] * r_memory[pos0] + pt->matrix[1] * r_memory[pos1];
				r_memory[pos1] = pt->matrix[2] * r_memory[pos0] + pt->matrix[3] * r_memory[pos1];
				r_memory[pos0] = tmp;
			}
			
		}
}

void PCpuExecution2_11(float complex* r_memory, float complex* w_memory, PT *pt, long mem_size, long n_threads, long desloc){
	long pos, pos0, pos1, base, shift, i;
	float complex tmp;
	
	//mem_size /= 2;

	omp_set_num_threads(n_threads);

	shift = 1 << pt->end;

	int gap = pow(2, pt->end + 1);
	int stride = pow(2, pt->end);


	//cout << "Mem size: " << mem_size << endl;
	//cout << "Iterations: " << ((mem_size/gap)*stride) << "   ---   " << (mem_size/gap) << " " << stride << endl;
	
	if (!pt->ctrl_count)
		#pragma omp parallel for private(tmp)
		for (int i = 0; i < mem_size; i+=gap){
			//#pragma omp parallel for private (tmp)
			for (int p = 0; p < stride; p++){
				tmp = pt->matrix[2] * r_memory[i+p] + pt->matrix[3] * r_memory[i+p+stride];
				r_memory[i+p] = pt->matrix[0] * r_memory[i+p] + pt->matrix[1] * r_memory[i+p+stride];
				r_memory[i+p+stride] = tmp;
			}
		}
	else
		//#pragma omp parallel for private(tmp)
		for (long i = 0; i < mem_size; i+=gap){
			for (long p = i; p < i+stride; p++){
				if ((p & pt->ctrl_mask) == pt->ctrl_value){
					tmp = pt->matrix[2] * r_memory[p] + pt->matrix[3] * r_memory[p+stride];
					r_memory[p] = pt->matrix[0] * r_memory[p] + pt->matrix[1] * r_memory[p+stride];
					r_memory[p+stride] = tmp;
				}
			}
		}
}

void PCpuExecution3_1(float complex *r_memory, float complex *w_memory, PT *pt, long mem_size, long n_threads, long desloc){ //Operador denso
	long base, shift, i, l;
	
	omp_set_num_threads(n_threads);
	
	i = 0;
	shift = pt->end;
		
	if (!pt->ctrl_count)
		#pragma omp parallel for private(l,base) //schedule(static, CHUNK)
		for (long pos = 0; pos < mem_size; pos+=n_threads){
			for (long p = 0; p< n_threads; p++){
				l = (((pos+p+desloc) >> shift) & 1) * 2; //LINE((pos+desloc), shift); 
				base = (pos+p+desloc) & (~(1 << shift)); //BASE((pos+desloc), shift);
		
				w_memory[pos+p] = pt->matrix[l] * r_memory[base] + 
								pt->matrix[l+1] * r_memory[base | (1 << shift)];
			}
		}
	else
		#pragma omp parallel for private(l,base) //schedule(static, CHUNK)
		for (long pos = 0; pos < mem_size; pos+=n_threads){
			for (long p = 0; p< n_threads; p++){
				if (((pos+p+desloc) & pt->ctrl_mask) == pt->ctrl_value){
					l = LINE(pos+p+desloc, shift);
					base = BASE(pos+p+desloc, shift);
		
					w_memory[pos+p] = pt->matrix[l] * r_memory[base] + 
									pt->matrix[l+1] * r_memory[base | (1 << shift)];
				}
				else
					w_memory[pos] = r_memory[pos+desloc];
			}
		}
}

void PCpuExecution4_0(float complex *r_memory, PT *pt, long mem_size, long n_threads, long desloc){
    long mt = pt->matrixType();

	switch (mt){
		case DENSE:
			PCpuExecution4_1(r_memory, pt, mem_size, n_threads, desloc);
			break;
		case DIAG_PRI:
			PCpuExecution4_2(r_memory, pt, mem_size, n_threads, desloc);
			break;
	   	case DIAG_SEC:
	   		PCpuExecution4_3(r_memory, pt, mem_size, n_threads, desloc);
	   		break;
	   	default:
	   		exit(1);
	}
}

void PCpuExecution4_1(float complex* r_memory, PT *pt, long mem_size, long n_threads, long desloc){
	long pos, pos0, pos1, base, shift, i;
	float complex temp0;
	
	mem_size /= 2;

	omp_set_num_threads(n_threads);

	shift = 1 << pt->end;
			
	if (!pt->ctrl_count)
		#pragma omp parallel for private(pos0,pos1,temp0)
		for (long pos = 0; pos < mem_size; pos++){
			pos0 = (pos * 2) - (pos & (shift-1));
			pos1 = pos0 | shift;

			temp0 = r_memory[pos0];
				
			r_memory[pos0] = pt->matrix[0] * r_memory[pos0] + pt->matrix[1] * r_memory[pos1];
			r_memory[pos1] = pt->matrix[2] * temp0 + pt->matrix[3] * r_memory[pos1];
		}
	else
		#pragma omp parallel for private(pos0,pos1,temp0)
		for (long pos = 0; pos < mem_size; pos++){
			pos0 = (pos * 2) - (pos & (shift-1));
			pos1 = pos0 | shift;

			if (((pos0 + desloc) & pt->ctrl_mask) == pt->ctrl_value){
				temp0 = r_memory[pos0];
				r_memory[pos0] = pt->matrix[0] * r_memory[pos0] + pt->matrix[1] * r_memory[pos1];
				r_memory[pos1] = pt->matrix[2] * temp0 + pt->matrix[3] * r_memory[pos1];			
			}
		}
}

void PCpuExecution4_2(float complex *r_memory, PT *pt, long mem_size, long n_threads, long desloc){ //Diagonal Principal
	omp_set_num_threads(n_threads);
	
	long shift = pt->end;
		
	if (!pt->ctrl_count)
		#pragma omp parallel for //schedule(static, CHUNK)
		for (long pos = 0; pos < mem_size; pos++)
			r_memory[pos] = pt->matrix[((pos >> shift) & 1) * 3] * r_memory[pos];
	else
		#pragma omp parallel for //schedule(static, CHUNK)
		for (long pos = 0; pos < mem_size; pos++)
			if (((pos+desloc) & pt->ctrl_mask) == pt->ctrl_value)
				r_memory[pos] = pt->matrix[((pos >> shift) & 1) * 3] * r_memory[pos];
			else
				r_memory[pos] = r_memory[pos];

}

void PCpuExecution4_3(float complex *r_memory, PT *pt, long mem_size, long n_threads, long desloc){ //Diagonal Secundária
	omp_set_num_threads(n_threads);

	long shift = pt->end;
		
	if (!pt->ctrl_count)
		#pragma omp parallel for //schedule(static, CHUNK)
		for (long pos = 0; pos < mem_size; pos++)
			r_memory[pos] = pt->matrix[LINE(pos, shift)] * r_memory[pos^(1<<shift)];
	else
		#pragma omp parallel for //schedule(static, CHUNK)
		for (long pos = 0; pos < mem_size; pos++)
			if (((pos+desloc) & pt->ctrl_mask) == pt->ctrl_value)
				r_memory[pos] = pt->matrix[LINE(pos, shift)] * r_memory[pos^(1<<shift)];
}


float complex* PCpuExecution3(float complex *r_memory, float complex *w_memory, PT **pts, int qubits, long n_threads){
	long mem_size, i;
	
	mem_size = pow(2.0, qubits)/n_threads;

	omp_set_num_threads(n_threads);
	
	#pragma omp parallel
	{
		long base, shift, i, l;
		long start = omp_get_thread_num() * mem_size;
		long end = start + mem_size;
		PT *pt;

		i = 0;
		while (pts[i] != NULL){
			pt = pts[i];
			shift = pt->end;
			for (long pos = start; pos < end; pos+=2){
				l = ((pos >> shift) & 1) * 2;
				base = pos & (~(1 << shift));
		
				w_memory[pos] = pt->matrix[l] * r_memory[base] + pt->matrix[l+1] * r_memory[base | (1 << shift)];

				l = (((pos+1) >> shift) & 1) * 2;
				base = (pos+1) & (~(1 << shift));
		
				w_memory[pos+1] = pt->matrix[l] * r_memory[base] + pt->matrix[l+1] * r_memory[base | (1 << shift)];
			}
			i++;
			#pragma omp barrier
			#pragma omp single
			{
				swap_ptr(&r_memory, &w_memory);
			}
		}
	}
	
	swap_ptr(&r_memory, &w_memory);
	
	return w_memory;
}

