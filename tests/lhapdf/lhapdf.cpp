#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/Utils.h"
#include <iostream>
#include <vector>
#include <string>
using namespace LHAPDF;
using namespace std;
 
 
int main() {

	const PDF* pdf = mkPDF("CT18NLO", 0);
	vector<int> pids = pdf->flavors();
	vector<string> spids;
	for (int pid : pids)
		spids.emplace_back(lexical_cast<string>(pid));

	for (const string& spid : spids)
		cout << spid << ", ";
	cout << endl;
	
	delete pdf;
	return 0;
}
