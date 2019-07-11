#include "sequence/io/sequence_lib.h"
#include "utils/utils.h"

void DisplayHelp(const char *program) {
  pfprintf(stderr, "Usage {s} <read_lib_file> <out_prefix>\n", program);
}

int main_build_lib(int argc, char **argv) {
  AutoMaxRssRecorder recorder;

  if (argc < 3) {
    DisplayHelp(argv[0]);
    exit(1);
  }
  SequenceLibCollection::Build(argv[1], argv[2]);

  return 0;
}