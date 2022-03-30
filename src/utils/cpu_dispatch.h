//
// Created by vout on 4/14/19.
//

#ifndef MEGAHIT_CPU_DISPATCH_H
#define MEGAHIT_CPU_DISPATCH_H

inline bool HasPopcnt() { return false; }
inline bool HasBmi2() { return false; }

#endif  // MEGAHIT_CPU_DISPATCH_H
