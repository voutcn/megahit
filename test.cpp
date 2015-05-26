#include <cstdio>
#include "sequence_package.h"

int main() {
	SequencePackage<> sp;
	sp.clear();
	sp.AppendSeq("acggt", 5);
	sp.AppendSeq("aaaaaaaaaggggaaaaaaaaaaaaaataaaaaaaaaaa", 39);

	for (int i = 0; i < sp.size(); ++i) {
		printf("Length: %d; ", sp.length(i));
		for (int j = 0; j < sp.length(i); ++j) {
			printf("%c", "ACGT"[sp.get_base(i, j)]);
		}
		puts("");
	}

	return 0;
}