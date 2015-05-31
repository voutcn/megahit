#include <string>

#include "omp.h"

#include "options_description.h"
#include "local_assembler.h"
#include "utils.h"

struct local_asm_opt_t {
	std::string contig_file;
	std::string read_file;

    int kmin;
    int kmax;
    int step;
    int seed_kmer;

    int min_contig_len;
    int sparsity;
    double similarity;
    double min_mapping_len;

    int num_threads;

    local_asm_opt_t() {
        kmin = 11;
        kmax = 41;
        step = 6;
        seed_kmer = 31;
        min_contig_len = 200;
        sparsity = 8;
        similarity = 0.95;
        min_mapping_len = 75;
        num_threads = 0;
    }
} opt;

void ParseOption(int argc, char *argv[]) {
    OptionsDescription desc;

    desc.AddOption("contig_file", "c", opt.contig_file, "contig file");
    desc.AddOption("read_file", "r", opt.read_file, "read file");
    desc.AddOption("kmin", "", opt.kmin, "");
    desc.AddOption("kmax", "", opt.kmax, "");
    desc.AddOption("step", "", opt.step, "");
    desc.AddOption("seed_kmer", "", opt.seed_kmer, "kmer size for seeding alignments");
    desc.AddOption("min_contig_len", "", opt.min_contig_len, "");
    desc.AddOption("min_mapping_len", "", opt.min_mapping_len, "");
    desc.AddOption("sparsity", "", opt.sparsity, "sparsity of hash mapper");
    desc.AddOption("similarity", "", opt.similarity, "alignment similarity threshold");
    desc.AddOption("num_threads", "t", opt.num_threads, "");

    try {
        desc.Parse(argc, argv);
        if (opt.contig_file == "") {
            throw std::logic_error("no contig file!");
        }
        if (opt.read_file == "") {
        	throw std::logic_error("no read file!");
        }
        if (opt.num_threads == 0) {
            opt.num_threads = omp_get_max_threads();
        }
    } catch (std::exception &e) {
        std::cerr << e.what() << std::endl;
        std::cerr << "Usage: " << argv[0] << " -c contigs.fa -r reads.fq > out.local_contig.fa" << std::endl;
        std::cerr << "options:" << std::endl;
        std::cerr << desc << std::endl;
        exit(1);
    }
}

static AutoMaxRssRecorder recorder;

int main(int argc, char **argv) {
	ParseOption(argc, argv);

    LocalAssembler la(opt.min_contig_len, opt.seed_kmer, opt.sparsity);
    la.set_kmer(opt.kmin, opt.kmax, opt.step);
    la.set_mapping_threshold(opt.similarity, opt.min_mapping_len);
    la.set_num_threads(opt.num_threads);

	la.ReadContigs(opt.contig_file.c_str());
    la.BuildHashMapper();
	la.AddReadLib(opt.read_file.c_str(), LocalAssembler::kFastx, true);
	la.EstimateInsertSize();
    la.MapToContigs();
    la.LocalAssemble();

	return 0;
}