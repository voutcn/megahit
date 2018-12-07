#include "read_lib_functions-inl.h"
#include "utils.h"

void DisplayHelp(const char *program) {
    fprintf(stderr, "Usage %s <read_lib_file> <out_prefix>\n", program);
}

int main_build_lib(int argc, char **argv) {
    AutoMaxRssRecorder recorder;

    if (argc < 3) {
        DisplayHelp(argv[0]);
        exit(1);
    }

    bool is_reverse = false;
    bool verbose = true;
    ReadAndWriteMultipleLibs(argv[1], is_reverse, argv[2], verbose);

    return 0;
}