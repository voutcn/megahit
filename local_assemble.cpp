#include <string>

#include "options_description.h"
#include "local_assembler.h"

struct local_asm_opt_t {
	std::string contig_file;
	std::string read_file;
} opt;

void ParseOption(int argc, char *argv[]) {
    OptionsDescription desc;

    desc.AddOption("contig_file", "c", opt.contig_file, "contig file");
    desc.AddOption("read_file", "r", opt.read_file, "read file");

    try {
        desc.Parse(argc, argv);
        if (opt.contig_file == "") {
            throw std::logic_error("no contig file!");
        }
        if (opt.read_file == "") {
        	throw std::logic_error("no read file!");
        }
    } catch (std::exception &e) {
        std::cerr << e.what() << std::endl;
        std::cerr << "Usage: " << argv[0] << " -c contigs.fa -r reads.fq" << std::endl;
        std::cerr << "options:" << std::endl;
        std::cerr << desc << std::endl;
        exit(1);
    }
}

int main(int argc, char **argv) {
	LocalAssembler la(200, 31);
	ParseOption(argc, argv);
	la.ReadContigs(opt.contig_file.c_str());
    la.BuildHashMapper();
	la.AddReadLib(opt.read_file.c_str(), LocalAssembler::kFastx, true);
	la.EstimateInsertSize();
    la.MapToContigs();

	return 0;
}