#include <cuComplex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"


bool error();
static int inst = 0;
static int call_count = 0;
static int call_peer_count = 0;


struct host_partial_transform{
	PT *pt;
	PT *pts;
	int num_pts;
	cuFloatComplex *matrix, *read_memory, *write_memory;
	long read_size, write_size;

	void malloc(long size, long pt_size){
		cudaMalloc((void**)&write_memory, size);
		cudaMalloc((void**)&read_memory, size);
		cudaMalloc((void**)&matrix, pow(4, pt_size)*sizeof(float complex));
	};
	void malloc(long r_size, long w_size, long pt_size){
		read_size = r_size;
		write_size = w_size;
		
		cudaMalloc((void**)&read_memory, read_size*sizeof(float complex));
		cudaMalloc((void**)&write_memory, write_size*sizeof(float complex));
		cudaMalloc((void**)&matrix, pow(4, pt_size)*sizeof(float complex));
	};

	void malloc_read(long size, long pt_size){
		write_memory = NULL;
		//cudaMalloc((void**)&write_memory, size);
		cudaMalloc((void**)&read_memory, size);
		cudaMalloc((void**)&matrix, pow(4, pt_size)*sizeof(float complex));	
	}

	void free(){
		if (write_memory) cudaFree(write_memory);
		cudaFree(read_memory);
		cudaFree(matrix); 
	};
	
	void swap(){
		cuFloatComplex **ptr1, **ptr2;
		ptr1 = &read_memory;
		ptr2 = &write_memory;
		
		
		cuFloatComplex *aux = *ptr1;
		*ptr1 = *ptr2;
		*ptr2 = aux;
	}
};
typedef host_partial_transform HPT;

struct DEV_OP{
	long arg[TAM_ARG];
	cuFloatComplex matrix[4];
};


extern "C" bool setDevice(int num = 0){
	return cudaFree(0);
}

extern "C" bool enablePeerAccess(){
	cudaSetDevice(0);
	cudaDeviceEnablePeerAccess(1, 0);
    
    cudaSetDevice(1);
    cudaDeviceEnablePeerAccess(0, 0);

    cudaGetLastError();

    return true;
}

__constant__ long c_arg[1][1];
__constant__ cuFloatComplex cmatrix[1][1];

__constant__ DEV_OP op[OPS_BLOCK];

__constant__ cuFloatComplex *gpu_pointer[4];


inline int GET_BLOCK_ID(PT *pt, int coalesc, int qbs_region){
	return (pt->end - coalesc)/(qbs_region-coalesc);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//extern "C"
template <int t_TAM_BLOCK, int t_REPT, int t_COALESC>
__global__ void ApplyValuesC01(int const b_pos, int const n_bits, int const count, int const shift, int const block_shift){
	long pos, block = (blockIdx.x + block_shift) * t_REPT;
	
	int i, c, thId = threadIdx.x;

	__shared__ cuFloatComplex s[t_REPT][t_TAM_BLOCK*2];

	long block_base[t_REPT];

	for (i = 0; i < t_REPT; i++){
		block_base[i] = (block + i) << t_COALESC;
		block_base[i] = (block_base[i] >> b_pos) << (b_pos + n_bits) | (block_base[i] & ((1 << b_pos) - 1)); //OPEN SPACE
	}
	
	
	for (i = 0; i < t_REPT; i++){
		pos = block_base[i] | ((thId >> t_COALESC) << b_pos) | (thId & ((1 << t_COALESC)-1));
		s[i][thId] = gpu_pointer[pos/shift][pos%shift];

		pos = pos | (1 << (b_pos+n_bits-1));
		s[i][thId+t_TAM_BLOCK] = gpu_pointer[pos/shift][pos%shift];
	}



	int pos0, pos1, op_bit;	
	cuFloatComplex tmp;
	
	for (c = 0; c < count; c++){
		__syncthreads();

		//if (print) printf("Ctrl G: M -> %ld    V -> %ld\n", op[c].arg[CTRL_MASK], op[c].arg[CTRL_VALUE]);
		//if (print) printf("Ctrl R: M -> %ld    V -> %ld\n", op[c].arg[CTRL_REG_MASK], op[c].arg[CTRL_REG_VALUE]);
		op_bit = 1 << op[c].arg[SHIFT];

		pos0 = (thId * 2) - (thId & (op_bit - 1));
		pos1 = pos0 | op_bit;
		
		for (i = 0; i < t_REPT; i++){
			if (((block_base[i] & op[c].arg[CTRL_MASK]) == op[c].arg[CTRL_VALUE]) && ((pos0 & op[c].arg[CTRL_REG_MASK]) == op[c].arg[CTRL_REG_VALUE])){

				//long ctrl = block_base | ((pos0 >> COALESC) << b_pos) | (pos0 & ((1 << COALESC)-1));

					//if (print) printf("threadIdx: %d  - opbit: %d  ---    pos: %d    e   %d\n", threadIdx.x, op_bit, pos0, pos1);

					tmp = cuCaddf(cuCmulf(s[i][pos0], op[c].matrix[0]), cuCmulf(s[i][pos1], op[c].matrix[1]));
					s[i][pos1] = cuCaddf(cuCmulf(s[i][pos0], op[c].matrix[2]), cuCmulf(s[i][pos1], op[c].matrix[3]));			
					s[i][pos0] = tmp;
			}
		}
	}
	__syncthreads();


	for (i = 0; i < t_REPT; i++){
		pos = block_base[i] | ((thId >> t_COALESC) << b_pos) | (thId & ((1 << t_COALESC)-1));
		gpu_pointer[pos/shift][pos%shift] = s[i][thId];

		pos = pos | (1 << (b_pos+n_bits-1));
		gpu_pointer[pos/shift][pos%shift] = s[i][thId+t_TAM_BLOCK];
	}
	
}



template <int t_TAM_BLOCK, int t_REPT, int t_COALESC>	
void GpuExecution01(float* r_memory, PT **pts, int qubits, int qbs_region, int multi_gpu, int num_it){
	//printf("%d  --  %d  --  %d  --  %d\n", t_TAM_BLOCK, t_REPT, t_COALESC, qbs_region);
	DEV_OP operators[OPS_BLOCK];

	inst = 0;

	dim3 block, dim;

	long mem_size = pow(2.0, qubits);
	long mem_desloc = mem_size/multi_gpu;
	
	//long rept = REPT;		//número de substate cada bloco fica responsável
	long nth = mem_size/multi_gpu/t_REPT/2;//2;	// /2 porque cada thread fica responsável por duas posições & /2 pelas 2 GPUS

	long malloc_size = (mem_size * (sizeof(float complex)))/multi_gpu;

	
	block.x = t_TAM_BLOCK;
	(nth > block.x)? dim.x = nth/block.x : block.x = nth;


	int block_region_size = log(block.x)/log(2) + 1;

	if (block_region_size < qbs_region){
		printf("ERRO: Região do bloco menor que a região de qubits\n");
		exit(1);
	}

	cuFloatComplex *gpu_mem[4];


	if (multi_gpu > 1){
		for (int d = 0; d < multi_gpu; d++){
			cudaSetDevice(d);
			for (int j = 0; j < multi_gpu; j++)
				if (d!=j) cudaDeviceEnablePeerAccess(j, 0);
		}
		cudaGetLastError();
	}
	

	for (int d = 0; d < multi_gpu; d++){
		cudaSetDevice(d);
		cudaMalloc(&gpu_mem[d], malloc_size); error();
		cudaMemcpy(gpu_mem[d], r_memory + (mem_desloc*2)*d, malloc_size, cudaMemcpyHostToDevice); error();
	}

	for (int d = 0; d < multi_gpu; d++){
		cudaSetDevice(d);
		cudaMemcpyToSymbol(gpu_pointer, gpu_mem, multi_gpu*sizeof(cuFloatComplex*)); error();
	}

	int i;
	for (int it = 0; it < num_it; it++){
		i = 0;

		while (pts[i]!= NULL){
			int qbs_block_id, region_start, is_peer, max_end, c = 0; //CONTADOR
			is_peer = 0;

			while (pts[i+c] != NULL &&
				   pts[i+c]->end < t_COALESC &&
				   c < OPS_BLOCK)
			{
				c++;
			}

			max_end = t_COALESC;

			if (pts[i+c] != NULL &&
				c < OPS_BLOCK)
			{
				qbs_block_id = GET_BLOCK_ID(pts[i+c], t_COALESC, qbs_region);

				do
				{
					max_end = max(max_end, pts[i+c]->end);
					c++;
				}
				while (pts[i+c] != NULL &&
					   qbs_block_id == GET_BLOCK_ID(pts[i+c], t_COALESC, qbs_region) &&
					   c < OPS_BLOCK);
			}

			region_start = max(t_COALESC, (max_end - (block_region_size - t_COALESC) + 1));

			is_peer = ((region_start + (block_region_size - t_COALESC)) > (qubits-multi_gpu+1));

			//printf("COUNT %d\nREGION_START %d\nREGION_BLOCK %d\nSHIFT %ld\nDIMX %d\n", c, region_start, (block_region_size - COALESC), mem_desloc, dim.x);
					

			for (int j = 0; j < c; j++){
				memcpy(operators[j].matrix, pts[i+j]->matrix, 4*sizeof(float complex));
				pts[i+j]->setArgsGPU(operators[j].arg, region_start, block_region_size, t_COALESC);
			}

			if (is_peer){
				for (int d = 0; d < multi_gpu; d++){
					cudaSetDevice(d);
					cudaDeviceSynchronize();
				}
			}

			for (int d = 0; d < multi_gpu; d++){
				cudaSetDevice(d);
				cudaMemcpyToSymbol(op, operators, c*sizeof(DEV_OP));
			}
			
			for (int d = 0; d < multi_gpu; d++){
				cudaSetDevice(d);
				ApplyValuesC01<t_TAM_BLOCK, t_REPT, t_COALESC><<<dim,block>>>(region_start, (block_region_size - t_COALESC), c, mem_desloc, dim.x * d);
			}

			
			if (is_peer){
				for (int d = 0; d < multi_gpu; d++){
					cudaSetDevice(d);
					cudaDeviceSynchronize();
				}
			}

			call_count++;
			if (is_peer) call_peer_count++;

			i += c;
		}
	}

	for (int d = 0; d < multi_gpu; d++){
		cudaMemcpy(r_memory + (mem_desloc*2)*d, gpu_mem[d], malloc_size, cudaMemcpyDeviceToHost); error();
		cudaFree(gpu_mem[d]); error();
	}

	//printf("Kernel Calls %d\nPeer Calls %d\n", call_count, call_peer_count);
}



template <int t_COALESC>
void GEWrapper2(float* r_memory, PT **pts, int qubits, int qbs_region, int multi_gpu, int tam_block, int rept, int num_it){
	switch(tam_block){
		case 32:
			switch(rept){
				case 1:
					GpuExecution01<32, 1, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				case 2:
					GpuExecution01<32, 2, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				case 4:
					GpuExecution01<32, 4, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				case 8:
					GpuExecution01<32, 8, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				case 16:
					GpuExecution01<32, 16, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				case 32:
					GpuExecution01<32, 32, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				default:
					printf("Invalid REPT");
			}
			break;
		case 64:
			switch(rept){
				case 1:
					GpuExecution01<64, 1, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				case 2:
					GpuExecution01<64, 2, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				case 4:
					GpuExecution01<64, 4, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				case 8:
					GpuExecution01<64, 8, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				case 16:
					GpuExecution01<64, 16, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				case 32:
					GpuExecution01<64, 32, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				default:
					printf("Invalid REPT");
			}
			break;
		case 128:
			switch(rept){
				case 1:
					GpuExecution01<128, 1, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				case 2:
					GpuExecution01<128, 2, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				case 4:
					GpuExecution01<128, 4, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				case 8:
					GpuExecution01<128, 8, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				case 16:
					GpuExecution01<128, 16, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				default:
					printf("Invalid REPT");
			}
			break;
		case 256:
			switch(rept){
				case 1:
					GpuExecution01<256, 1, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				case 2:
					GpuExecution01<256, 2, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				case 4:
					GpuExecution01<256, 4, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				case 8:
					GpuExecution01<256, 8, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				default:
					printf("Invalid REPT");
			}
			break;
		case 512:
			switch(rept){
				case 1:
					GpuExecution01<512, 1, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				case 2:
					GpuExecution01<512, 2, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				case 4:
					GpuExecution01<512, 4, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				default:
					printf("Invalid REPT");
			}
			break;
		case 1024:
			switch(rept){
				case 1:
					GpuExecution01<1024, 1, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				case 2:
					GpuExecution01<1024, 2, t_COALESC>(r_memory, pts, qubits, qbs_region, multi_gpu, num_it);
					break;
				default:
					printf("Invalid REPT");
			}
			break;
		default:
			printf("Invalid TAM_BLOCK");
	}
}


extern "C" float* GpuExecutionWrapper(float* r_memory, PT **pts, int qubits, int qbs_region, int multi_gpu, int tam_block, int rept, int coalesc, int num_it){
	switch(coalesc){
		case 0:
			GEWrapper2<0>(r_memory, pts, qubits, qbs_region, multi_gpu, tam_block, rept, num_it);
			break;
		case 1:
			GEWrapper2<1>(r_memory, pts, qubits, qbs_region, multi_gpu, tam_block, rept, num_it);
			break;
		case 2:
			GEWrapper2<2>(r_memory, pts, qubits, qbs_region, multi_gpu, tam_block, rept, num_it);
			break;
		case 3:
			GEWrapper2<3>(r_memory, pts, qubits, qbs_region, multi_gpu, tam_block, rept, num_it);
			break;
		case 4:
			GEWrapper2<4>(r_memory, pts, qubits, qbs_region, multi_gpu, tam_block, rept, num_it);
			break;
		case 5:
			GEWrapper2<5>(r_memory, pts, qubits, qbs_region, multi_gpu, tam_block, rept, num_it);
			break;
		case 6:
			GEWrapper2<6>(r_memory, pts, qubits, qbs_region, multi_gpu, tam_block, rept, num_it);
			break;
		case 7:
			GEWrapper2<7>(r_memory, pts, qubits, qbs_region, multi_gpu, tam_block, rept, num_it);
			break;
		case 8:
			GEWrapper2<8>(r_memory, pts, qubits, qbs_region, multi_gpu, tam_block, rept, num_it);
			break;
		default:
			printf("Invalid COALESC");
	}

	return r_memory;
}

bool error(){
	inst++;
	cudaError_t e;
	e = cudaGetLastError();
	if (e == cudaSuccess) return false;
	printf("inst: %d\nerror: %d - %s\n", inst, e, cudaGetErrorString (e));
	exit(1);
	return true;
}
