#include <iostream>
#include "../kmer.h"
#include "../compact_sequence.h"
#include "../unitig_graph.h"

using namespace std;

int main() {
	cout << sizeof(CompactSequence) << endl;
	cout << sizeof(UnitigGraphVertex) << endl;
	cout << sizeof(Kmer<1>) << endl;
	cout << sizeof(Kmer<2>) << endl;
	cout << sizeof(Kmer<3>) << endl;
	cout << sizeof(Kmer<4>) << endl;
	return 0;
}