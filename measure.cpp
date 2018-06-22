#include <iostream>
#include <vector>
#include <algorithm>
#include <complex>
#include <unistd.h>
#include <omp.h>
#include <sys/time.h>


using namespace std;

#define NUM_AMOSTRAS 12
#define NUM_THREADS 4

#define EXC_RANGE 1

int timeval_subtract(struct timeval *result, struct timeval *t2, struct timeval *t1)
{
    long int diff = (t2->tv_usec + 1000000 * t2->tv_sec) - (t1->tv_usec + 1000000 * t1->tv_sec);
    result->tv_sec = diff / 1000000;
    result->tv_usec = diff % 1000000;

    return (diff<0);
}


bool print_statistic(vector <float> amostra){
	float med, desv;
	med = desv = 0;
	int number = amostra.size();

	sort(amostra.begin(), amostra.end());

	for (int i = EXC_RANGE; i < number - EXC_RANGE; i++){
//		cout << amostra[i] << endl;
		med += amostra[i];
	}
	med /= (number - EXC_RANGE*2);

	for (int i = EXC_RANGE; i < number - EXC_RANGE; i++){
		desv += pow((amostra[i] - med),2.0);
	}

	desv = sqrt(desv/(number-1)) / med * 100.0;

	//if (desv > 1) return false;

	cout << med << "\t" << desv << endl;
	return true;
}



int measure(complex<float> *state, int qubits, int k, int width, int n_threads){
	long size = pow(2.0, qubits);
	int n = 1 << width;
	float m[n] = {0};

	int mask = (n - 1);
	cout << size << "  " << n <<  "   " << mask << endl;

	//float one, zero;

	//omp_set_num_threads(n_threads);
	//#pragma omp parallel for reduction(+:m[0:n]) //if (cond)
	for (long i = 0; i < size; i++){
		cout << pow(state[i],2.0) << endl;
		m[(i>>k)&mask] += pow(state[i].real(),2.0) + pow(state[i].imag(),2.0);
	}

	cout << m[0] << "  " << m[1] << endl;


	return 0;
}



int main(int argc, char** argv){
	struct timeval timev, tvBegin, tvEnd;
	float t;


	int qubits = atoi(argv[1]);

	long mem_size = pow(2.0, qubits);

	complex<float> *state = new complex<float>[mem_size];
	cout << 1.0/mem_size << endl;
	complex<float> value = sqrt(1.0/mem_size);
	cout << value << endl;
	for (int i = 0; i < mem_size; i++){
		state[i] = value;
//		cout << state[i] << endl;
	}
	cout << pow(state[0],2.0) << endl;

	gettimeofday(&tvBegin, NULL);
	measure(state, qubits, atoi(argv[2]), atoi(argv[3]), 1);
	gettimeofday(&tvEnd, NULL);
	timeval_subtract(&timev, &tvEnd, &tvBegin);
	t = timev.tv_sec + (timev.tv_usec / 1000000.0);

	delete state;


	return 0;

	vector <float> amostras;

	for (int qubits = 18; qubits < 29; qubits++){
		cout << "Qubits " << qubits << endl;
		long mem_size = pow(2.0, qubits);
		complex<float> *state = new complex<float>[mem_size];


		cout << 1.0/mem_size << endl;
		complex<float> value = sqrt(1.0/mem_size);
		cout << value << endl;
		for (int i = 0; i < mem_size; i++){
			state[i] = value;
	//		cout << state[i] << endl;
		}
		cout << pow(state[0], 2.0);
		for (int n_threads = 16; n_threads <= 16; n_threads*=2){
			for (int width = 1; width <= 1; width++){
				amostras.clear();
				for (int a = 0; a < NUM_AMOSTRAS; a++){
					gettimeofday(&tvBegin, NULL);
					measure(state, qubits, 0, width, n_threads);
					gettimeofday(&tvEnd, NULL);
					timeval_subtract(&timev, &tvEnd, &tvBegin);
					t = timev.tv_sec + (timev.tv_usec / 1000000.0);

					amostras.push_back(t);
				}
				print_statistic(amostras);
			}
			cout << endl;
		}
		cout << endl;
		delete state;
	}

	return 0;
}
