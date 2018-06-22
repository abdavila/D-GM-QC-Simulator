#include <iostream>
#include <vector>
#include "dgm.h"
#include "complex.h"
#include <unistd.h>
#include <papi.h>

using namespace std;

#define NUM_AMOSTRAS 30
#define NUM_THREADS 4

#define EXC_RANGE 5

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

	cout << "\nMEDIA: " << med << "\t\tDESV: " << desv << endl;
	return true;

}

#define NUM_EVENTS_L1 5
#define NUM_EVENTS_L2 7
#define NUM_EVENTS_L2_2 7
#define NUM_EVENTS_L3 9
//#define NUM_EVENTS_TLB 2
#define NUM_EVENTS_OTHERS 7


//int Events_L1[NUM_EVENTS_L1] = {PAPI_L1_DCM,PAPI_L1_ICM,PAPI_L1_TCM,PAPI_L1_LDM,PAPI_L1_STM};
//int Events_L2[NUM_EVENTS_L2] = {PAPI_L2_DCM,PAPI_L2_ICM,PAPI_L2_TCM,PAPI_L2_STM,PAPI_L2_DCH,PAPI_L2_DCA,PAPI_L2_DCR};
//int Events_L2_2[NUM_EVENTS_L2_2] = {PAPI_L2_DCW,PAPI_L2_ICH,PAPI_L2_ICA,PAPI_L2_ICR,PAPI_L2_TCA,PAPI_L2_TCR,PAPI_L2_TCW};
//int Events_L3[NUM_EVENTS_L3] = {PAPI_L3_TCM,PAPI_L3_DCA,PAPI_L3_DCR,PAPI_L3_DCW,PAPI_L3_ICA,PAPI_L3_ICR,PAPI_L3_TCA,PAPI_L3_TCR,PAPI_L3_TCW};
//int Events_TLB[NUM_EVENTS_TLB] = {PAPI_TLB_DM,PAPI_TLB_IM};
//int Events_Others[NUM_EVENTS_OTHERS] = {PAPI_TLB_DM,PAPI_TLB_IM,PAPI_FP_INS,PAPI_FP_OPS,PAPI_STL_ICY,PAPI_TOT_CYC,PAPI_TOT_INS};

//const char* const Events_Names_L1[] = {"PAPI_L1_DCM","PAPI_L1_ICM","PAPI_L1_TCM","PAPI_L1_LDM","PAPI_L1_STM"};
//const char* const Events_Names_L2[] = {"PAPI_L2_DCM","PAPI_L2_ICM","PAPI_L2_TCM","PAPI_L2_STM","PAPI_L2_DCH","PAPI_L2_DCA","PAPI_L2_DCR","PAPI_L2_DCW","PAPI_L2_ICH","PAPI_L2_ICA","PAPI_L2_ICR","PAPI_L2_TCA","PAPI_L2_TCR","PAPI_L2_TCW"};
//const char* const Events_Names_L3[] = {"PAPI_L3_TCM","PAPI_L3_DCA","PAPI_L3_DCR","PAPI_L3_DCW","PAPI_L3_ICA","PAPI_L3_ICR","PAPI_L3_TCA","PAPI_L3_TCR","PAPI_L3_TCW"};
//const char* const Events_Names_TLB[] = {"PAPI_TLB_DM","PAPI_TLB_IM"};
//const char* const Events_Names_OTHERS[] = {"PAPI_FP_INS","PAPI_FP_OPS","PAPI_STL_ICY","PAPI_TOT_CYC","PAPI_TOT_INS"};

//long_long values_L1[NUM_EVENTS_L1];																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																															
//long_long values_L2[NUM_EVENTS_L2];
//long_long values_L3[NUM_EVENTS_L3];
//long_long values_TLB[NUM_EVENTS_TLB];
//long_long values_OTHERS[NUM_EVENTS_OTHERS];
/*
void print_events(vector <int> events, vector <long_long> values){
	char name[PAPI_MAX_STR_LEN];

	for (int i = 0; i < events.size(); i++){
		PAPI_event_code_to_name(events[i], name);
		cout << name << "\t" << values[i] << endl;
	}
}
*/

int main(int argc, char** argv){
	int qubits = atoi(argv[1]);
	//int had_pos = atoi(argv[2]);
	//int factor = atoi(argv[2]);
	//int stress = atoi(argv[3]);

	/*
	int retval;

	retval = PAPI_library_init(PAPI_VER_CURRENT);

	if (!PAPI_is_initialized()){
		cout << "SHIT " << retval << endl;
		exit(1);
	}

	vector <int*> list;
	list.push_back(Events_L1);
	list.push_back(Events_L2);
	list.push_back(Events_L2_2);
	list.push_back(Events_L3);
//	list.push_back(Events_TLB);
	list.push_back(Events_Others);

	vector <int> size_list;
	size_list.push_back(NUM_EVENTS_L1);
	size_list.push_back(NUM_EVENTS_L2);
	size_list.push_back(NUM_EVENTS_L2_2);
	size_list.push_back(NUM_EVENTS_L3);
//	size_list.push_back(NUM_EVENTS_TLB);
	size_list.push_back(NUM_EVENTS_OTHERS);


	vector <int> events;
	vector <long_long> values;

	int i = 0 | PAPI_PRESET_MASK;
	PAPI_event_info_t info;

	printf("Name\t\t\t\t\t\t   Description\n");
	do {
		retval = PAPI_get_event_info(i, &info);
		if (retval == PAPI_OK && info.count > 0) {
			printf("%-30s %s\n", info.symbol, info.long_descr);
		}
	} while (PAPI_enum_event(&i, PAPI_ENUM_ALL) == PAPI_OK);
	*/

	//string function = "ID,ID,ID,ID,ID";
	//string function = "H,H,H,H,H";

	string had = "ID";
	if (had_pos == (qubits-1)) had = "H";

	for (int i = 1; i < qubits; i++){
		if ((qubits-i-1) == had_pos) had += ",H";
		else had += ",ID";
	}

	cout << had << endl;

	vector<string> f;

	for (int i = 0; i < stress; i++) f.push_back(had);

	float complex *state = (float complex*) calloc(pow(2, qubits), sizeof(float complex));
	state[0] = 1;

	state = GenericExecute(state, f, qubits, t_CPU, 4, factor);

	return 0;

	/*

	cout << "\n\nLets Get Started\n\n" << endl;
	string had;// = "H,H";//,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H,H";
	vector <float> amostra;

	//for (had_pos = qubits-1; had_pos >= 0; had_pos--){
		had = "ID";
		if (had_pos == (qubits-1)) had = "H";

		for (int i = 1; i < qubits; i++){
			if ((qubits-i-1) == had_pos) had += ",H";
			else had += ",ID";
		}

		cout << "HAD_POS: " << had_pos << endl;
		cout << had << endl;

		vector <string> function;

		for (int i = 0; i < stress; i++) function.push_back(had);

		struct timeval t, tBegin, tEnd;
		float te;

		amostra.clear();

		for (int a = 0; a < list.size(); a++){
			float complex *state = (float complex*) calloc(pow(2, qubits), sizeof(float complex));
			state[0] = 1;

			gettimeofday(&tBegin, NULL);

			events.clear();
			events.assign(list[a], list[a]+size_list[a]);

			values.resize(events.size());


//			int PA = PAPI_start_counters(&events[0], events.size());

//			if (PA != PAPI_OK){
//				cout << PA << endl;
//				exit(1);
//			}

			state = GenericExecute(state, function, qubits, t_CPU, 4, 1);

//			if (PAPI_stop_counters(&values[0], values.size()) != PAPI_OK) exit(1);

//			print_events(events, values);

			gettimeofday(&tEnd, NULL);
			timeval_subtract(&t, &tEnd, &tBegin);
			te = t.tv_sec + (t.tv_usec / 1000000.0);

			//cout << te << endl;

			amostra.push_back(te);
			//free(state);
		}

		cout << "\n" << endl;

		//if (!print_statistic(amostra)){
		//	cout << "RETRY " << had_pos << endl;
			//had_pos++;
		//}
		//printf("Time: %f\n", te);
	//}

	//printMem(state, 4);
	*/
	return 0;
}
