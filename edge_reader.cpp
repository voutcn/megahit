#include <cstdio>
#include <cstdlib>
#include <cstdint>

using namespace std;

int main(int argc, char **argv) {
	if (argc < 2) {
		fprintf(stderr, "Usage %s <k_size> <edge_file>\n", argv[0]);
		exit(1);
	}

	int k = atoi(argv[1]);
	FILE *edge_file = fopen(argv[2], "rb");

	unsigned words_per_edge = ((k + 1) * 2 + 2 + 31) / 32; // DivCeiling
	uint32_t* edge = (uint32_t*) malloc(sizeof(uint32_t) * words_per_edge);

	while (fread(edge, sizeof(uint32_t), words_per_edge, edge_file) == words_per_edge) {
		uint32_t dollar_flag = edge[(k+1)/16] >> (k+1)%16*2 & 3;
		if (dollar_flag & 1) {
			// first char of the edge is $
			putchar('$');
		} else {
			putchar("ACGT"[edge[0] & 3]);
		}

		for (int i = 1; i < k; ++i) {
			putchar("ACGT"[edge[i/16] >> i%16*2 & 3]);
		}

		if (dollar_flag & 2) {
			// last char of the edge is $
			putchar('$');
		} else {
			putchar("ACGT"[edge[k/16] >> k%16*2 & 3]);
		}
		puts("");
	}

	return 0;
}