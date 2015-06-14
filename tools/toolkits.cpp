#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main_contig2fastg(int argc, char** argv);
int main_read_stat(int argc, char **argv);
int main_trim_lowq_tail(int argc, char **argv);
int main_filter_by_len(int argc, char **argv);

void show_help(const char *program_name) {
	fprintf(stderr, "Usage: %s <sub_program> [sub options]\n"
					"    sub-programs:\n"
	                "       contig2fastg          convert MEGAHIT's k*.contigs.fa to fastg format that can be viewed by Bandage\n"
	                "       readstat              calculate read stats (# of reads, bases, longest, shortest, average)\n"
	                "       trim                  trim low quality tail of fastq reads\n"
	                "       filterbylen           filter contigs by length\n",
	                program_name);
}

int main(int argc, char **argv) {
	if (argc < 2) {
		show_help(argv[0]);
		exit(1);
	}

	if (strcmp(argv[1], "contig2fastg") == 0) {
		return main_contig2fastg(argc - 1, argv + 1);
	} else if (strcmp(argv[1], "readstat") == 0) {
		return main_read_stat(argc - 1 , argv + 1);
	} else if (strcmp(argv[1], "trim") == 0) {
		return main_trim_lowq_tail(argc - 1, argv + 1);
	} else if (strcmp(argv[1], "filterbylen") == 0) {
		return main_filter_by_len(argc - 1, argv + 1);
	} else {
		show_help(argv[0]);
		return 1;
	}
}