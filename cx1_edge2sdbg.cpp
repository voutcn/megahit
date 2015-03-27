#include "cx1_edge2sdbg.h"

namespace cx1_edge2sdbg {

// helpers


// cx1 core functions
int64_t encode_lv1_diff_base(int64_t read_id, edge2sdbg_global_t &g);
void    read_edge_prepare(edge2sdbg_global_t &g); // num_items_, num_cpu_threads_ and num_output_threads_ must be set here
void*   lv0_calc_bucket_size(void*); // pthread working function
void    init_global_and_set_cx1(edge2sdbg_global_t &g);
void*   lv1_fill_offset(void*); // pthread working function
void*   lv2_extract_substr(void*); // pthread working function
void    lv2_sort(edge2sdbg_global_t &g);
void    lv2_pre_output_partition(edge2sdbg_global_t &g);
void*   lv2_output(void*); // pthread working function
void    lv2_post_output(edge2sdbg_global_t &g);
void    post_proc(edge2sdbg_global_t &g);

}