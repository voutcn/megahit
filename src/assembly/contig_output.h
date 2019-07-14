//
// Created by vout on 11/21/18.
//

#ifndef MEGAHIT_CONTIG_OUTPUT_H
#define MEGAHIT_CONTIG_OUTPUT_H

#include <cassert>
#include <cstdint>
#include <string>

#include "sequence/io/contig/contig_writer.h"

class UnitigGraph;

void OutputContigs(UnitigGraph &graph, ContigWriter *contig_writer,
                   ContigWriter *final_contig_writer, bool change_only,
                   uint32_t min_standalone);

#endif  // MEGAHIT_CONTIG_OUTPUT_H
