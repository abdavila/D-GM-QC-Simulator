#include "common.h"
#include <math.h>

float round_precision = 0.000000001;

PT::PT(){
	matrix = NULL;
	ctrl_pos = ctrl_rest = NULL;
}

void PT::destructor(){
	if ((mat_size != 1) && !matrix) free(matrix);
	if (!ctrl_pos) free(ctrl_pos);
	if (!ctrl_rest) free(ctrl_rest);
}
	
long PT::ctrlAffect(long qubit){
	long i;
	for (i = 0; i < ctrl_count; i++) if (ctrl_pos[i] < qubit) return i;

	return ctrl_count;
}

long PT::matrixType(){
	if ((matrix[1] == 0.0) && (matrix[2] == 0.0)) return DIAG_PRI;
	else if ((matrix[0] == 0.0) && (matrix[3] == 0.0)) return DIAG_SEC;
	
	return DENSE;
}

/*
void PT::ctrlRest(long afect){
	long p, v, aux;
	long n = pow(2,ctrl_count-afect);
	ctrl_rest_count = n - 1;
	//if (ctrl_rest != NULL) free(ctrl_rest);
	if (n) ctrl_rest = (long*)malloc(sizeof(long)*n - 1);
	
	p = 0;
	for (long i = 0; i < n; i++){
		v = 0;
		aux = i;
		for (long c = ctrl_count - 1; c >= afect; c--){
			v = v | ((aux & 1) << ctrl_pos[c]);
			aux = aux >> 1;
		}
		if (v != ctrl_value){
			ctrl_rest[p] = v;
			p++;
		}
	}
}
*/
void PT::setArgs(long *arg, long affect){
	long ctrl_affect = ctrlAffect(affect);
	//ctrlRest(ctrl_afect);
				
	//arg[MAT_SIZE] = mat_size;
	arg[SHIFT] = end;
	arg[CTRL_MASK] = ctrl_mask;
	arg[CTRL_VALUE] = ctrl_value;
	//arg[CTRL_COUNT] = ctrl_count - ctrl_affect;
	//arg[CTRL_REST] = ctrl_rest_count;
		
	//arg[MAT_START] = 0;
	//arg[MAT_END] = arg[MAT_SIZE];
	//arg[ACUMM] = 0;	
}

void PT::setArgs_soft(long *arg, long affect){
	//long ctrl_affect = ctrlAffect(affect);
	//ctrlRest(ctrl_afect);
				
	//arg[MAT_SIZE] = mat_size;
	arg[SHIFT] = end;
	arg[CTRL_MASK] = ctrl_mask;
	arg[CTRL_VALUE] = ctrl_value;
	//arg[CTRL_COUNT] = ctrl_count - ctrl_affect;
	//arg[CTRL_REST] = ctrl_rest_count;
		
	//arg[MAT_START] = 0;
	//arg[MAT_END] = arg[MAT_SIZE];
	//arg[ACUMM] = 0;	
}

void PT::setArgsGPU(long *arg, int region_start, int region_size, int coalesc){
	long mask_coalesc = (1 << coalesc) - 1;
	long mask_region = ((1<<(region_size-coalesc)) - 1) << region_start;
	long mask_rest = ~(mask_coalesc | mask_region);



	arg[SHIFT] = end;
	if (end >= coalesc)	arg[SHIFT] = end - region_start + coalesc;

	if (ctrl_mask){
		arg[CTRL_MASK] = ctrl_mask & mask_rest;
		arg[CTRL_VALUE] = ctrl_value & mask_rest;

		arg[CTRL_REG_MASK] = (ctrl_mask & mask_coalesc) | ((ctrl_mask & mask_region) >> (region_start - coalesc));


		arg[CTRL_REG_VALUE] = (ctrl_value & mask_coalesc) | ((ctrl_value & mask_region) >> (region_start - coalesc));
	}
	else
	{
		arg[CTRL_MASK] = arg[CTRL_VALUE] = arg[CTRL_REG_MASK] = arg[CTRL_REG_VALUE] = 0;
	}
	/*
	if (ctrl_mask){
		printf("\n\nEnd: %d\nregion_start: %ld\nregion_size: %ld\ncoalesc: %ld\nCtrl-(value: %ld, mask: %ld, count: %ld)\n**mask**:%ld\n", end, region_start, region_size, coalesc, ctrl_value, ctrl_mask, ctrl_count, mask);
		printf("Global Ctrl-(value: %ld, mask: %ld)\n", arg[CTRL_VALUE], arg[CTRL_MASK]);
		printf("Region Ctrl-(value: %ld, mask: %ld)\n----------------------------------------\n\n", arg[CTRL_REG_VALUE], arg[CTRL_REG_MASK]);
	}
	*/

}
	
void PT::print(){
	printf("qubits: %d\nmat_size: %d\nstart: %d\nend: %d\nAffect: %d\nCtrl-(value: %ld, mask: %ld, count: %ld)\n", qubits, mat_size, start, end, affected, ctrl_value, ctrl_mask, ctrl_count);
	
	for (int i = 0; i < ctrl_count; i++)
		printf("%d: %ld\n", i, ctrl_pos[i]);

	/*
	printf("Rest:\n");
	for (int i = 0; i < ctrl_rest_count; i++){
		printf("%ld\n", ctrl_rest[i]);
	}
	*/
}
	
void PT::printMatrix(){
	for (int i = 0; i < mat_size; i++){
		for (int j = 0; j < mat_size; j++)
			printf("%d: %.4f, %.4f  \t", i*mat_size+j,crealf(matrix[i*mat_size+j]), cimagf(matrix[i*mat_size+j]));
		printf("\n");
	}

}

void printMem(float complex* mem, int qubits){
	long size = pow(2, qubits);
	float range = 1.0/pow(2,21);

	float real, imag, f;
	for (long i = 0; i < size; i++){
		real = imag = 0;
		f = fabs(crealf(mem[i]));		
		//if (f > round_precision)
			real = crealf(mem[i]);

		f = fabs(cimagf(mem[i]));
		//if (f > round_precision)
			imag = cimagf(mem[i]);
		//if (real != 0 || imag != 0)
			printf("%ld:\t%.6f %.6f\n", i, real, imag);
	}
}

void printMemExp(float complex* mem, int qubits, int reg1, int reg2, long n){
	long size = pow(2, qubits);
	float range = 1.0/pow(2,21);

	long mask = pow(2, n) - 1;

	float real, imag, f;
	long last_X = 0;
	for (long i = 0; i < size; i++){
		real = imag = 0;
		f = fabs(crealf(mem[i]));
	   	if (f > round_precision)
			real = crealf(mem[i]);

		f = fabs(cimagf(mem[i]));
		if (f > round_precision)
			imag = cimagf(mem[i]);

		if ((imag != 0) || (real != 0)){
			long X = (i>>(qubits-reg1-n))&mask;
			long Exp = (i>>(qubits-reg2-n))&mask;
			printf("%ld\t>>  X: %ld\tExp: %ld\tDif: %ld\t\t\tV: %f %f\n", i, X, Exp, X-last_X, crealf(mem[i]), cimagf(mem[i]));
			last_X = X;
			//printf("%ld >>> X: %ld\tB: %ld\tOver: b%ld %ld\tV: %f %f\n", i, (i>>(n+2))&mask, (i>>1)&(mask), (i&1), ((i>>(n+1))&1), crealf(mem[i]), cimagf(mem[i]));
		}
	}
}

void printMemCheckExp(float complex* mem, int qubits, long width, long a, long N){
	long size = pow(2, qubits);
	float range = 1.0/pow(2,21);

	long mask = pow(2, width) - 1;

	float real, imag, f;
	long last_X = 0;
	for (long i = 0; i < size; i++){
		real = imag = 0;
		f = fabs(crealf(mem[i]));
	   	if (f > round_precision)
			real = crealf(mem[i]);

		f = fabs(cimagf(mem[i]));
		if (f > round_precision)
			imag = cimagf(mem[i]);

		if ((imag != 0) || (real != 0)){
			long X = (i>>(2*width+2));
			long Exp = (i>>(width+2))&mask;
			bool ok = true;

			if (modular_pow(a,X,N) != Exp) ok = "Erro";

			printf("%ld\t>>  X: %ld\tExp: %ld\tDif: %ld\t\t", i, X, Exp, X-last_X);
			last_X = X;
			if (modular_pow(a,X,N) != Exp) printf("Errado\n");
			else printf("\n");

			//printf("%ld >>> X: %ld\tB: %ld\tOver: b%ld %ld\tV: %f %f\n", i, (i>>(n+2))&mask, (i>>1)&(mask), (i&1), ((i>>(n+1))&1), crealf(mem[i]), cimagf(mem[i]));
		}
	}
}

long modular_pow(long base, long exponent, long modulus){
    long result = 1;
    base = base % modulus;
	exponent = exponent;
    while (exponent){
        if (exponent&1) result = (result * base) % modulus;
        exponent = exponent >> 1;
        base = (base * base) % modulus;
	}
    return result;
}


bool increasing(const PT *pt1, const PT *pt2){
	if (pt1->affected == pt2->affected)
		return pt1->start > pt2->start;
	else
		return pt2->affected;
}

bool decreasing(const PT *pt1, const PT *pt2){
	if (pt1->affected == pt2->affected)
		return pt1->start < pt2->start;
	else
		return pt1->affected;
}

////////////////////////////////////////////////////////////
void swap_value(int *v1, int *v2){
	int aux = *v1;
	*v1 = *v2;
	*v2 = aux;
}

void swap_ptr(float **ptr1, float **ptr2){
	float *aux = *ptr1;
	*ptr1 = *ptr2;
	*ptr2 = aux;
}

void swap_ptr(float complex **ptr1, float complex **ptr2){
	float complex *aux = *ptr1;
	*ptr1 = *ptr2;
	*ptr2 = aux;
}
//////////////////////////////////////////////////////////

int timeval_subtract(struct timeval *result, struct timeval *t2, struct timeval *t1)
{
    long int diff = (t2->tv_usec + 1000000 * t2->tv_sec) - (t1->tv_usec + 1000000 * t1->tv_sec);
    result->tv_sec = diff / 1000000;
    result->tv_usec = diff % 1000000;

    return (diff<0);
}


