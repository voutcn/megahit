//
// Created by vout on 4/14/19.
//

#ifndef MEGAHIT_CPU_DISPATCH_H
#define MEGAHIT_CPU_DISPATCH_H

inline bool HasPopcnt() {
  unsigned eax, ebx, ecx, edx;
#ifdef _MSC_VER
  int cpuid[4];
  __cpuid(cpuid, 1);
  eax = cpuid[0], ebx = cpuid[1], ecx = cpuid[2], edx = cpuid[3];
#else
  asm volatile("cpuid\n\t"
               : "=a"(eax), "=b"(ebx), "=c"(ecx), "=d"(edx)
               : "0"(1));
#endif
  return ecx >> 23U & 1U;
}

inline bool HasBmi2() {
  unsigned eax, ebx, ecx, edx;
#ifdef _MSC_VER
  int cpuid[4];
  __cpuidex(cpuid, 7, 0);
  eax = cpuid[0], ebx = cpuid[1], ecx = cpuid[2], edx = cpuid[3];
#else
  asm volatile("cpuid\n\t"
               : "=a"(eax), "=b"(ebx), "=c"(ecx), "=d"(edx)
               : "0"(7), "2"(0));
#endif
  return ebx >> 8U & 1U;
}

#endif  // MEGAHIT_CPU_DISPATCH_H
