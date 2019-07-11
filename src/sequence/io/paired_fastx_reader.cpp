//
// Created by vout on 5/11/19.
//

#include "paired_fastx_reader.h"

int64_t PairedFastxReader::Read(SeqPackage *pkg, int64_t max_num,
                                int64_t max_num_bases, bool reverse) {
  int64_t num_bases = 0;
  for (int64_t i = 0; i < max_num; i += 2) {
    auto r0 = readers_[0]->ReadNext();
    auto r1 = readers_[1]->ReadNext();
    if (r0 && r1) {
      int b0 = 0, e0 = r0->seq.l;
      int b1 = 0, e1 = r1->seq.l;

      if (trim_n_) {
        FastxReader::TrimN(r0->seq.s, r0->seq.l, &b0, &e0);
        FastxReader::TrimN(r1->seq.s, r1->seq.l, &b1, &e1);
      }

      if (reverse) {
        pkg->AppendReversedStringSequence(r0->seq.s + b0, e0 - b0);
        pkg->AppendReversedStringSequence(r1->seq.s + b1, e1 - b1);
      } else {
        pkg->AppendStringSequence(r0->seq.s + b0, e0 - b0);
        pkg->AppendStringSequence(r1->seq.s + b1, e1 - b1);
      }

      num_bases += e0 - b0 + e1 - b1;

      if (num_bases >= max_num_bases) {
        return i + 2;
      }
    } else {
      return i;
    }
  }
  return max_num;
}
