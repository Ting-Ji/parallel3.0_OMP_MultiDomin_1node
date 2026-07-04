#ifndef _WRITE_TEC_H
#define _WRITE_TEC_H
#include "Geo.h"
#include <fstream>
int write_tecplot(int num_total_node, int num_ele, int num_node_per_ele, int num_result_per_node, Point * node_array, int ** node_list, double ** node_result, int plot_result, int write_disp, char * out1, double amplitude);
#endif
