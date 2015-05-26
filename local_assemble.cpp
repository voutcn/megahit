#include "local_assembler.h"

struct local_asm_opt_t {

};

int main(int argc, char **argv) {
	LocalAssembler la(200, 31);
	la.ReadContigs("/nas1/dhli/dbg/data/7.iowa.prairie/iowa_prairie.contigs.k27_10_87.300.reupload.fna.gz");
	la.AddReadLib("-", LocalAssembler::kFastx, false);

	return 0;
}