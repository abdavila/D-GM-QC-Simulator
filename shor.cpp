#include "dgm.h"
#include "common.h"
#include "gates.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

#define NUM_AMOSTRAS 14

#define EXC_RANGE 2

//////////////////////////////////////////////////////
//string concatena(vector <string> vec, int size, bool rev = false);
string int2str(int number);
long mul_inv(long a, long b);
long modular_pow(long base, long exponent, long modulus);
//////////////////////////////////////////////////////

void ApplyQFT(int qubits, int type, int threads);

vector <string> QFT(int qubits, int reg, int over, int width);
vector <string> QFT2(int qubits, int reg, int width);
vector <string> RQFT(int qubits, int reg, int over, int width);
vector <string> CSwapR(int qubits, int ctrl, int reg1, int reg2, int width);
vector <string> SwapOver(int qubits, int reg, int width);
string genRot(int qubits, int reg, long value);

vector <string> CU(int qubits, int ctrl, int reg1, int reg2, int width, long a, long N);

vector <string> CMultMod(int qubits, int ctrl, int reg1, int reg2, int over, int over_bool, int width, long a, long N);
vector <string> CRMultMod(int qubits, int ctrl, int reg1, int reg2, int over, int over_bool, int width, long a, long N);
vector <string> C2AddMod(int qubits, int ctrl1, int ctrl2, int reg, int over, int over_bool, int width, long a, long N);
vector <string> C2SubMod(int qubits, int ctrl1, int ctrl2, int reg, int over, int over_bool, int width, long a, long N);

vector <string> CAdjust(int qubits, int ctrl, int reg, int width, long N);

string C2AddF(int qubits, int ctrl1, int ctrl2, int reg, int over, long num, int width);
string CAddF(int qubits, int ctrl1, int reg, int over, long num, int width);
string AddF(int qubits, int reg, int over, long num, int width);
vector <string> AddF(int qubits, int reg, int over, long num, int width, bool controlled);

string C2SubF(int qubits, int ctrl1, int ctrl2, int reg, int over, long num, int width);
string CSubF(int qubits, int ctrl1, int reg, int over, long num, int width);
string SubF(int qubits, int reg, int over, long num, int width);
vector <string> SubF(int qubits, int reg, int over, long num, int width, bool controlled);

//string CNot(int qubits, int ctrl, int target, int cv = 1);
//string Toffoli(int qubits, int ctrl1, int ctrl2, int target, int cv = 3);
string Pauly_X(int qubits, int reg, int width);
string HadN(int qubits, int reg, int width);

float print_statistic(vector <float> amostra){
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

	cout << med << "\t" << desv;

	return med;

}


int
revert_bits(int res, int n){
	int c = 0;

	for (int i =0; i < n; i++){
		c = (c<<1) | (res&1);
		res = res >> 1;
	}
	return c;
}

int
quantum_ipow(int a, int b)
{
  int i;
  int r=1;

  for(i=0; i<b ;i++)
    r*=a;

  return r;
}

/* Calculate the greatest common divisor with Euclid's algorithm */

int
quantum_gcd(int u, int v)
{
  int r;

  while(v)
    {
	  r = v;
	  v = u % v;
	  u = r;
      //r = u % v;
      //u = v;
      //v = r;
    }
  return u;
}


void
quantum_frac_approx(int *a, int *b, int width)
{
  float f = (float) *a / *b;
  float g=f;
  int i, num2=0, den2=1, num1=1, den1=0, num=0, den=0;

  do
    {
      i = (int) (g+0.000005);

      g -= i-0.000005;
      g = 1.0/g;

      if (i * den1 + den2 > 1<<width)
	break;

      num = i * num1 + num2;
      den = i * den1 + den2;

      num2 = num1;
      den2 = den1;
      num1 = num;
      den1 = den;

    } while(fabs(((double) num / den) - f) > 1.0 / (2 * (1 << width)));

  *a = num;
  *b = den;

  return;
}


void ApplyQFT(int qubits, int type, int multi_gpu, int qbs_region, int coalesc, int tam_block, int rept){
	DGM dgm;
	dgm.exec_type = type;
	dgm.multi_gpu = multi_gpu;

	dgm.qbs_region = qbs_region;
	dgm.coalesc = coalesc;
	dgm.tam_block = tam_block;
	dgm.rept = rept;

	dgm.qubits = qubits;
	dgm.allocateMemory();
	dgm.setMemoryValue(0);

	vector<string> qft = QFT2(qubits,0,qubits);

	dgm.executeFunction(qft);
}

//N - Number to ne factored
//type - Execution Type
//threads - Number of threads to be used in case of a parallel execution on CPU
void Shor(long N, int type,  int multi_gpu, int qbs_region, int coalesc, int tam_block, int rept){
	long a, n, mod_a, mod_inv_a, aux, m, res;

	int qubits, qft_qb, reg1, reg2, over, over_bool;
	int f1, f2, factor;

	aux = N;
	a = n = 0;
	while (aux){
		n++;
		aux = aux >> 1;
	}
	qubits = 2*n+3;

	qft_qb = 0;
	reg1 = 1;
	reg2 = n+2;
	over = n+1;
	over_bool = qubits - 1;


	while((quantum_gcd(N, a) > 1) || (a < 2)){
		a = rand() % N;
	}

	//cout << "Seed:\t" << a << endl;

	DGM dgm;
	dgm.exec_type = type;
	dgm.multi_gpu = multi_gpu;

	dgm.qbs_region = qbs_region;
	dgm.coalesc = coalesc;
	dgm.tam_block = tam_block;
	dgm.rept = rept;

	dgm.qubits = qubits;
	dgm.allocateMemory();
	dgm.setMemoryValue((1<<(n+2)));

	string X0 = Pauly_X(qubits, 0, 1);
	string H0 = HadN(qubits, qft_qb, 1);

//	cout << "AQUI" << endl;

	res = 0;
	int L = 2*n-1;
	long inv_a = mul_inv(a,N);

	vector <string> func, f;

	for (int i = L; i >= 0; i--){
		mod_a = modular_pow(a, pow(2,i), N);
		mod_inv_a = modular_pow(inv_a, pow(2,i), N);

//	        cout << "AQUI " << i << endl;

		func.clear();

		//cout << mod_a << " " << mod_inv_a << endl;

//	        cout << "AQUI W" << endl;

		//if (mod_a != 1){
		func.push_back(H0);
//	        cout << "AQUI" << endl;

		f = CMultMod(qubits, qft_qb, reg1, reg2, over, over_bool, n, mod_a, N);
//        	cout << "AQUI" << endl;


		func.insert(func.end(), f.begin(), f.end());
		f = CSwapR(qubits, qft_qb, reg1, reg2, n);
		func.insert(func.end(), f.begin(), f.end());
//	        cout << "AQUI" << endl;

		f = CRMultMod(qubits, qft_qb, reg1, reg2, over, over_bool, n, mod_a, N);
		func.insert(func.end(), f.begin(), f.end());
		func.push_back(H0);
		//}
		//else cout << 1 << endl;
//	        cout << "AQUI X" << endl;

		//printMemExp(dgm.r_mem, qubits, reg1, reg2, n);

		if (res) func.push_back(genRot(qubits, qft_qb, res));

		//cout << "Passo " << i << " " << func.size() << endl;
		//for (int l = 0; l < func.size(); l++){
		//	cout << func[l] << endl;
		//}

//	        cout << "AQUI Y" << endl;

		dgm.executeFunction(func);

//	        cout << "AQUI Z" << endl;

		m = dgm.measure(qft_qb);
//		cout << m << endl;
		res = (res << 1) | m;
	}

	//return;

//       cout << "AQUI" << endl;


  int c = revert_bits(res, 2*n);

  //cout << c << "   " << res << endl;

  if(c==0)
  {
      //printf("Fail - Measured Zero.\n");
      return;
  }

  int q = 1<<(2*n);

  //printf("Measured %i (%f), ", c, (float)c/q);

  quantum_frac_approx(&c, &q, n);

  //printf("fractional approximation is %i/%i.\n", c, q);


	int r = q;
	int i = 1;
	while ((r*i) < (1<<n)){
		if (modular_pow(a, r*i, N) == 1){
			q = r * i;
			break;
		}
		i++;
	}
	if (q >= N) q = r;


/*
  if((q % 2 == 1) && (2*q<(1<<n)))
    {
      printf("Odd denominator, trying to expand by 2.\n");
      q *= 2;
    }

  if(q % 2 == 1)
    {
      printf("Odd period, try again.\n");
      return;
    }
*/

  //printf("Possible period is %i.\n", q);

//long modular_pow(long base, long exponent, long modulus)
  //f1 = quantum_ipow(a, q/2) + 1 % N;
 // f1 = modular_pow(a, q/2, N) + 1;
  //f2 = quantum_ipow(a, q/2) - 1 % N;
//  f2 = modular_pow(a, q/2, N) - 1;

	i = modular_pow(a, q/2, N);
	f1 = quantum_gcd(N, i+1);
	f2 = quantum_gcd(N, i-1);

	if(f1>f2)
		factor=f1;
	else
		factor=f2;

    if((factor < N) && (factor > 1))
    {
    	//printf("Sucess\n");
    	//printf("%ld = %i * %i\n", N, factor, (int)N/factor);
		return;
    }

	if (r!=q){
	  i = modular_pow(a, r/2, N);
		f1 = quantum_gcd(N, i+1);
		f2 = quantum_gcd(N, i-1);

		if(f1>f2)
			factor=f1;
		else
			factor=f2;

  		if((factor < N) && (factor > 1)){
  			//printf("Sucess\n");
			//printf("R: %ld = %i * %i\n", N, factor, (int)N/factor);
			return;
		}
	}
	//printf("FAIL\n");

}

string genRot(int qubits, int reg, long value){
	vector <string> func(qubits, "ID");
	string name;

	int k = 2;
	float complex rot, eps;
	eps = M_E;

	rot = 1;
	while (value){
		if (value&1) rot *= cpowf(eps, -2*M_PI*I/pow(2.0, k));
		value = value >> 1;
		k++;
	}

	if (rot != 1){
		Gates g;
		name = "Rot_" + int2str(value);
		g.addGate(name, 1.0, 0.0, 0.0, rot);
		func[reg] = name;

		return concatena(func, qubits);
	}

	return "";
}

vector <string> CU(int qubits, int ctrl, int reg1, int reg2, int width, long a, long N){
	vector <string> m, rm, sw, u;
/*
	m = CMultMod(qubits, ctrl, reg1, reg2, width, a, N);
	rm = CRMultMod(qubits, ctrl, reg1, reg2, width, mul_inv(a,N), N);
	sw = CSwapR(qubits, ctrl, reg1, reg2+1, width);


	u = m;
	u.insert(u.end(), sw.begin(), sw.end());
	u.insert(u.end(), rm.begin(), rm.end());
*/
	return u;
}

vector<string> CMultMod(int qubits, int ctrl, int reg1, int reg2, int over, int over_bool, int width, long a, long N){
//        cout << "MULT" << endl;

	int ctrl2;
	vector <string> qft = QFT(qubits, reg2, over, width);
//        cout << "MULT" << endl;


	string HN = HadN(qubits, reg2, width);

//        cout << "MULT" << endl;

	vector <string> rqft = RQFT(qubits, reg2, over, width);

//        cout << "MULT" << endl;

	//////////////////////////////////////////////////////////////

//        cout << "MULT" << endl;


	vector <string> mult_mod;
	vector <string> am;
	mult_mod.push_back(HadN(qubits, over, 1));
	mult_mod.push_back(HN);

	ctrl2 = reg1 + width - 1;
	for (int i = 0; i < width; i++){
		am = C2AddMod(qubits, ctrl, ctrl2-i, reg2, over, over_bool, width, a, N);
		//for (int j = 0; j < am.size(); j++) cout << am[j] << endl;
		//exit(1);
		mult_mod.insert(mult_mod.end(), am.begin(), am.end());

		a = (a*2)%N;
	}

//        cout << "MULT" << endl;


	mult_mod.insert(mult_mod.end(), rqft.begin(), rqft.end());

	return mult_mod;

}

vector<string> CRMultMod(int qubits, int ctrl, int reg1, int reg2, int over, int over_bool, int width, long a, long N){
	int ctrl2;
	vector <string> qft = QFT(qubits, reg2, over, width);
	vector <string> rqft = RQFT(qubits, reg2, over, width);

	//////////////////////////////////////////////////////////////

	vector <string> mult_mod;
	vector <string> am;

	ctrl2 = reg1 + width - 1;
	for (int i = 0; i < width; i++){
		am = C2SubMod(qubits, ctrl, ctrl2-i, reg2, over, over_bool, width, a, N);
		mult_mod.insert(mult_mod.begin(), am.begin(), am.end());

		a = (a*2)%N;
	}

	mult_mod.insert(mult_mod.begin(), qft.begin(), qft.end());
	mult_mod.insert(mult_mod.end(), rqft.begin(), rqft.end());

	return mult_mod;
}

vector <string> C2AddMod(int qubits, int ctrl1, int ctrl2, int reg, int over, int over_bool, int width, long a, long N){
	vector<string> qft = QFT(qubits, reg, over, width);
	vector<string> rqft = RQFT(qubits, reg, over, width);

	string c2_add_a = C2AddF(qubits, ctrl1, ctrl2, reg, over, a, width);
	string c2_sub_a = C2SubF(qubits, ctrl1, ctrl2, reg, over, a, width);

	string sub_N = SubF(qubits, reg, over, N, width);
	string c_add_N = CAddF(qubits, over_bool, reg, over, N, width);

	string n_over = Pauly_X(qubits, over, 1);
	string c_over = CNot(qubits, over, over_bool);

	vector <string> func;

	func.push_back(c2_add_a);
	func.push_back(sub_N);
	func.insert(func.end(), rqft.begin(), rqft.end());
	func.push_back(c_over);
	func.insert(func.end(), qft.begin(), qft.end());
	func.push_back(c_add_N);
	func.push_back(c2_sub_a);
	func.insert(func.end(), rqft.begin(), rqft.end());
	func.push_back(n_over);
	func.push_back(c_over);
	func.push_back(n_over);
	func.insert(func.end(), qft.begin(), qft.end());
	func.push_back(c2_add_a);

	return func;
}

vector <string> C2SubMod(int qubits, int ctrl1, int ctrl2, int reg, int over, int over_bool, int width, long a, long N){
	vector <string> qft = QFT(qubits, reg, over, width);
	vector <string> rqft = RQFT(qubits, reg, over, width);

	string c2_add_a = C2AddF(qubits, ctrl1, ctrl2, reg, over, a, width);
	string c2_sub_a = C2SubF(qubits, ctrl1, ctrl2, reg, over, a, width);

	string add_N = AddF(qubits, reg, over, N, width);
	string c_add_N = CAddF(qubits, over_bool, reg, over, N, width);
	string c_sub_N = CSubF(qubits, over_bool, reg, over, N, width);

	string n_over = Pauly_X(qubits, over, 1);
	string c_over = CNot(qubits, over, over_bool);

	vector <string> func;

	func.push_back(c2_sub_a);
	func.insert(func.end(), rqft.begin(), rqft.end());
	func.push_back(n_over);
	func.push_back(c_over);
	func.push_back(n_over);
	func.insert(func.end(), qft.begin(), qft.end());
	func.push_back(c2_add_a);
	func.push_back(c_sub_N);
	func.insert(func.end(), rqft.begin(), rqft.end());
	func.push_back(c_over);
	func.insert(func.end(), qft.begin(), qft.end());
	func.push_back(add_N);
	func.push_back(c2_sub_a);

	return func;
}

//////////////////////////////////////////////////////////////////////////
/*
string CNot(int qubits, int ctrl, int target, int cv){
	vector <string> cn (qubits, "ID");
	cn[ctrl] = "Control1(0)";
	if (cv) cn[ctrl] = "Control1(1)";
	cn[target] = "Target1(X)";

	return concatena(cn, qubits);
}

string Toffoli(int qubits, int ctrl1, int ctrl2, int target, int cv){
	vector <string> tf (qubits, "ID");
	tf[ctrl1] = "Control1(0)";
	if (cv>>1) tf[ctrl1] = "Control1(1)";
	tf[ctrl2] = "Control1(0)";
	if (cv&1) tf[ctrl2] = "Control1(1)";
	tf[target] = "Target1(X)";

	return concatena(tf, qubits);
}
*/
string Pauly_X(int qubits, int reg, int width){
	vector <string> px (qubits, "ID");
	for (int i = 0; i < width; i++) px[i+reg] = "X";

	return concatena(px, qubits);
}


string HadN(int qubits, int reg, int width){
	vector <string> hn (qubits, "ID");
	for (int i = 0; i < width; i++) hn[i+reg] = "H";

	return concatena(hn, qubits);
}

//////////////////////////////////////////////////////////////////////////

string CAddF(int qubits, int ctrl1, int reg, int over, long num, int width){
	vector <string> caf = AddF(qubits, reg, over, num, width, true);

	caf[ctrl1] = "Control1(1)";

	return concatena(caf, qubits);
}


string C2AddF(int qubits, int ctrl1, int ctrl2, int reg, int over, long num, int width){
	vector <string> caf = AddF(qubits, reg, over, num, width, true);

	caf[ctrl1] = "Control1(1)";
	caf[ctrl2] = "Control1(1)";

	return concatena(caf, qubits);
}

string AddF(int qubits, int reg, int over, long num, int width){
	return concatena(AddF(qubits, reg, over, num, width, false), qubits);
}

vector <string> AddF(int qubits, int reg, int over, long num, int width, bool controlled){
	int size = width+1;
	vector <float complex> rot (size, 1);
	float complex  c;

	Gates g;

	long aux = num;

	float complex eps = M_E;

	for (int i = 0; i < size; i++){
		if (aux&1)
			for (int j = i; j < size; j++)
				rot[j] *= cpowf(eps, 2*M_PI*I/pow(2.0, j-i+1));
		aux = aux >> 1;
	}

	vector<string> add(qubits, "ID");
	string name;

	aux = reg+width-1;
	c = 1;
	for (int i = 0; i < size; i++){
		if (rot[i] != c){
			name = "ADD_" + int2str(num) + "_" + int2str(i);
			g.addGate(name, 1.0, 0.0, 0.0, rot[i]);
			if (controlled) name = "Target1(" + name + ")";
			add[aux-i] = name;
		}
	}

	name = add[reg-1];
	add[reg-1] = "ID";
	add[over] = name;

	return add;
}

string CSubF(int qubits, int ctrl1, int reg, int over, long num, int width){
	vector <string> csf = SubF(qubits, reg, over, num, width, true);

	csf[ctrl1] = "Control1(1)";

	return concatena(csf, qubits);
}

string C2SubF(int qubits, int ctrl1, int ctrl2, int reg, int over, long num, int width){
	vector <string> csf = SubF(qubits, reg, over, num, width, true);

	csf[ctrl1] = "Control1(1)";
	csf[ctrl2] = "Control1(1)";

	return concatena(csf, qubits);
}

string SubF(int qubits, int reg, int over, long num, int width){
	return concatena(SubF(qubits, reg, over, num, width, false), qubits);
}

vector <string> SubF(int qubits, int reg, int over, long num, int width, bool controlled){
	long size = width+1;
	vector <float complex> rot (size, 1);
	float complex  c;

	Gates g;

	long aux = num;

	float complex eps = M_E;
	for (int i = 0; i < size; i++){
		if (aux&1)
			for (int j = i; j < size; j++)
				rot[j] *= cpowf(eps, -2*M_PI*I/pow(2.0, j-i+1));
		aux = aux >> 1;
	}

	vector<string> sub(qubits, "ID");
	string name;

	aux = reg+width-1;
	c = 1;
	for (int i = 0; i < size; i++){
		if (rot[i] != c){
			name = "SUB_" + int2str(num) + "_" + int2str(i);
			g.addGate(name, 1.0, 0.0, 0.0, rot[i]);
			if (controlled) name = "Target1(" + name + ")";
			sub[aux-i] = name;
		}
	}

	name = sub[reg-1];
	sub[reg-1] = "ID";
	sub[over] = name;

	return sub;
}

vector <string> QFT(int qubits, int reg, int over, int width){

//        cout << "QFT" << endl;

	string s, name;
	vector <string> qft;

	Gates g;
	float complex c;
	for (int i = 1; i <= width+1; i++){
                name = "R" + int2str(i);

//		cout << "QFT " << i << endl;
		c = M_E;
		c = cpowf(c, 2*M_PI*I/pow(2.0, i));

//                cout << "QFT " << i << endl;
		g.addGate(name, 1.0, 0.0, 0.0, c);
//                cout << "QFT " << i << endl;
	}

//        cout << "QFT" << endl;


	vector <string> base (qubits, "ID");

	qft.push_back(HadN(qubits, over, 1));
	for (int j = 0; j < width; j++){
		base[j+reg] = "Control1(1)";
		base[over] = "Target1(R" + int2str(j+2) + ")";

		s = concatena(base, qubits);
		qft.push_back(s);
		base[j+reg] = "ID";
	}
	base[over] = "ID";

	for (int i = 0; i < width; i++){
		qft.push_back(HadN(qubits, i+reg, 1));

		for (int j = i+1; j < width; j++){
			base[j+reg] = "Control1(1)";
			base[i+reg] = "Target1(R" + int2str(j-i+1) + ")";

			s = concatena(base, qubits);
			qft.push_back(s);

			base[j+reg] = "ID";
		}
		base[i+reg] = "ID";
	}

//        cout << "QFT" << endl;


	return qft;
}

vector <string> QFT2(int qubits, int reg, int width){
	string s;
	vector <string> qft;

	Gates g;
	float complex c;
	for (int i = 1; i <= width+1; i++){
		c = M_E;
		c = cpowf(c, 2*M_PI*I/pow(2.0, i));
		g.addGate("R-" + int2str(i), 1.0, 0.0, 0.0, c);
	}

	vector <string> base (qubits, "ID");

	for (int i = 0; i < width; i++){
		base[i+reg] = "H";
		s = concatena(base, qubits);
		qft.push_back(s);

		for (int j = i+1; j < width; j++){
			base[j+reg] = "Control1(1)";
			base[i+reg] = "Target1(R-" + int2str(j-i+1) + ")";

			s = concatena(base, qubits);
			qft.push_back(s);

			base[j+reg] = "ID";
		}
		base[i+reg] = "ID";
	}

	return qft;
}

vector <string> RQFT(int qubits, int reg, int over, int width){
	string s;
	vector <string> rqft;

	Gates g;
	float complex c;
	for (int i = 1; i <= width+1; i++){
		c = M_E;
		c = cpowf(c, -2*M_PI*I/pow(2.0, i));
		g.addGate("R'" + int2str(i), 1.0, 0.0, 0.0, c);
	}

	vector <string> base (qubits, "ID");


	rqft.push_back(HadN(qubits, over, 1));
	for (int j = 0; j < width; j++){
		base[j+reg] = "Control1(1)";
		base[over] = "Target1(R'" + int2str(j+2) + ")";

		s = concatena(base, qubits);
		rqft.push_back(s);
		base[j+reg] = "ID";
	}
	base[over] = "ID";

	for (int i = 0; i < width; i++){
		base[i+reg] = "H";
		s = concatena(base, qubits);
		rqft.push_back(s);

		for (int j = i+1; j < width; j++){
			base[j+reg] = "Control1(1)";
			base[i+reg] = "Target1(R'" + int2str(j-i+1) + ")";

			s = concatena(base, qubits);
			rqft.push_back(s);

			base[j+reg] = "ID";
		}
		base[i+reg] = "ID";
	}
	reverse(rqft.begin(), rqft.end());

	return rqft;
}

vector <string> CSwapR(int qubits, int ctrl, int reg1, int reg2, int width){
	vector <string> sw;
	vector <string>	base (qubits, "ID");
	string s1, s2;

	for (int i = 0; i < width; i++){
		base[ctrl] = "Control1(1)";
		base[i+reg1] = "Target1(X)";
		base[i+reg2] = "Control1(1)";
		s1 = concatena(base, qubits);

		base[ctrl] = "ID";
		base[i+reg1] = "Control1(1)";
		base[i+reg2] = "Target1(X)";
		s2 = concatena(base, qubits);

		base[i+reg1] = base[i+reg2] = "ID";

		sw.push_back(s2);
		sw.push_back(s1);
		sw.push_back(s2);
	}

	return sw;
}

vector <string> SwapOver(int qubits, int reg, int width){
	vector <string> so;

	for(int i=0; i<width/2; i++){
		so.push_back(CNot(qubits, reg+width-i-1, reg+i));
		so.push_back(CNot(qubits, reg+i, reg+width-i-1));
		so.push_back(CNot(qubits, reg+width-i-1, reg+i));
    }

	return so;
}

/////////////////////////////////////////////////////
/*
string concatena(vector <string> vec, int size, bool rev){
	string s;
	if (!rev){
		s = vec[0];
		for (int i = 1; i < size; i++)
			s += "," + vec[i];
	}
	else{
		s = vec[size-1];
		for (int i = size - 2; i >= 0; i--)
			s += "," + vec[i];
	}

	return s;
}
*/

long mul_inv(long a, long b){
	long b0 = b, t, q;
	long x0 = 0, x1 = 1;
	if (b == 1) return 1;
	while (a > 1) {
		q = a / b;
		t = b, b = a % b, a = t;
		t = x0, x0 = x1 - q * x0, x1 = t;
	}
	if (x1 < 0) x1 += b0;

	return x1;
}

string int2str(int number){
	stringstream ss;
	ss << number;

	string str = ss.str();

	return str;
}

//////////////////////////////////////////////////////

int main(int argc, char** argv){
	srand (time(NULL));

	struct timeval timev, tvBegin, tvEnd;
	float t;
	vector <float> amostras;

	//QFT

	int qbs_region = 9;
	int tam_block = 256;
	int rept = 4;
	int coalesc = 4;

	for (int qubits = atoi(argv[1]); qubits <= atoi(argv[2]); qubits++){
		for (int a = 0; a < NUM_AMOSTRAS; a++){
			gettimeofday(&tvBegin, NULL);
			ApplyQFT(qubits, t_GPU, 1, qbs_region, coalesc, tam_block, rept);
			gettimeofday(&tvEnd, NULL);
			timeval_subtract(&timev, &tvEnd, &tvBegin);
			t = timev.tv_sec + (timev.tv_usec / 1000000.0);

			amostras.push_back(t);
		}
		print_statistic(amostras);
		cout << "\t";
		amostras.clear();

		for (int a = 0; a < NUM_AMOSTRAS; a++){
			gettimeofday(&tvBegin, NULL);
			ApplyQFT(qubits, t_GPU, 2, qbs_region, coalesc, tam_block, rept);
			gettimeofday(&tvEnd, NULL);
			timeval_subtract(&timev, &tvEnd, &tvBegin);
			t = timev.tv_sec + (timev.tv_usec / 1000000.0);

			amostras.push_back(t);
		}
		print_statistic(amostras);
		cout << "\n";
		amostras.clear();

		cout << "####################################" << endl << endl;
	}
	
	/* // SHOR
	vector <int> vet;

	//vet.push_back(57);
	//vet.push_back(119);
	//vet.push_back(253);
	//vet.push_back(485);
	vet.push_back(1017);
	vet.push_back(2045);

	int qbs_region = 9;
	int tam_block = 256;
	int rept = 4;

	for (int pos = 0; pos < vet.size(); pos++){
		cout << "VALUE: " << vet[pos] << endl;

		for (int coalesc = 4; coalesc <=4; coalesc++){

			for (int a = 0; a < NUM_AMOSTRAS; a++){
				gettimeofday(&tvBegin, NULL);
				Shor(vet[pos], t_GPU, 1, qbs_region, coalesc, tam_block, rept);
				gettimeofday(&tvEnd, NULL);
				timeval_subtract(&timev, &tvEnd, &tvBegin);
				t = timev.tv_sec + (timev.tv_usec / 1000000.0);

				amostras.push_back(t);
			}
			print_statistic(amostras);
			cout << "\t";
			amostras.clear();

			for (int a = 0; a < NUM_AMOSTRAS; a++){
				gettimeofday(&tvBegin, NULL);
				Shor(vet[pos], t_GPU, 2, qbs_region, coalesc, tam_block, rept);
				gettimeofday(&tvEnd, NULL);
				timeval_subtract(&timev, &tvEnd, &tvBegin);
				t = timev.tv_sec + (timev.tv_usec / 1000000.0);

				amostras.push_back(t);
			}
			print_statistic(amostras);
			cout << "\n";
			amostras.clear();
		}
		cout << "####################################" << endl << endl;
	}
	*/

	return 0;
}
