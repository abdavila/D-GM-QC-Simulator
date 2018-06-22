# Build tools
NVCC = nvcc $(ARCH)
CXX = g++-5
GCC = gcc-5
ARCH = -arch=sm_52

#QBS_REGION = 4
#D = -D QBS_REGION=$(QBS_REGION)
OPS_BLOCK=200

# here are all the objects
GPUOBJS = kernel.o
OBJS = dgm.o common.o gates.o genMem.o

# make and compile

#gpu: main.o $(OBJS) $(GPUOBJS)
#	$(NVCC) -o gpu.out main.o $(OBJS) $(GPUOBJS) -Xcompiler "-fopenmp"

gpu: main2.o final.o
	$(NVCC) -o gpu.out main2.o final.o -Xcompiler "-fopenmp"

cpu: main2.o $(OBJS)
	$(CXX) -o cpu.out main2.o $(OBJS) -fopenmp

cpu_qft: shor.o $(OBJS)
	$(CXX) -o cpu_qft.out shor.o $(OBJS) -fopenmp


final.o: $(OBJS) $(GPUOBJS)
	ld -r $(OBJS) $(GPUOBJS) -o final.o

shor: shor.o $(OBJS) $(GPUOBJS)
	$(NVCC) -o shor.out shor.o $(OBJS) $(GPUOBJS) -Xcompiler "-fopenmp"

grover: grover.o $(OBJS) $(GPUOBJS)
	$(NVCC) -o grover.out grover.o $(OBJS) $(GPUOBJS) -Xcompiler "-fopenmp"

fuzzy: fuzzy.o $(OBJS) $(GPUOBJS)
	$(NVCC) -o fuzzy.out fuzzy.o $(OBJS) $(GPUOBJS) -Xcompiler "-fopenmp"

grover.o: grover.cpp
	$(CXX) -c grover.cpp -fopenmp 

shor.o: shor.cpp
	$(CXX) -c shor.cpp -fopenmp 

fuzzy.o: fuzzy.cpp
	$(CXX) -c fuzzy.cpp -fopenmp 

kernel.o: kernel.cu
	$(NVCC) -c -D OPS_BLOCK=$(OPS_BLOCK) kernel.cu
	
main.o: main.cpp
	$(CXX) -c main.cpp

main2.o: main2.cpp
	$(CXX) -c main2.cpp

dgm.o: dgm.cpp
	$(CXX) -c dgm.cpp -fopenmp -O3 -fcx-limited-range

gates.o: gates.cpp
	$(CXX) -c gates.cpp

genMem.o: genMem.cpp
	$(CXX) -c genMem.cpp
	
common.o: common.c
	$(CXX) -c common.c 
	
clean:
	rm *.o *.out
