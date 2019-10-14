//
// Created by vout on 11/21/18.
//

#ifndef MEGAHIT_LOW_DEPTH_REMOVER_H
#define MEGAHIT_LOW_DEPTH_REMOVER_H

#include <cstdint>

class UnitigGraph;

uint32_t RemoveLowDepth(UnitigGraph &graph, double min_depth);
bool RemoveLocalLowDepth(UnitigGraph &graph, double min_depth, uint32_t max_len,
                         uint32_t local_width, double local_ratio,
                         bool permanent_rm, uint32_t *num_removed);
uint32_t IterateLocalLowDepth(UnitigGraph &graph, double min_depth,
                              uint32_t min_len, uint32_t local_width,
                              double local_ratio, bool permanent_rm = false);

#endif  // MEGAHIT_LOW_DEPTH_REMOVER_H
