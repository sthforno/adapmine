#include "../include/common.h"
#include "../include/motif_generator.h"
#include "../include/pattern.h"
#include "../include/schedule.h"
#include "omp.h"
#include <../include/dataloader.h>
#include <../include/graph.h>

#include <algorithm>
#include <assert.h>
#include <iostream>
#include <string>

double test_pattern(Graph *g, Pattern &pattern) {

    bool is_pattern_valid;
    int performance_modeling_type;
    bool use_in_exclusion_optimize;

    performance_modeling_type = 1;
    use_in_exclusion_optimize = true;
    Schedule_IEP schedule_our(pattern, is_pattern_valid, performance_modeling_type, 1, use_in_exclusion_optimize, g->v_cnt, g->e_cnt);
    assert(is_pattern_valid == true);

    double t1, t2;
    double total_time = 0;

    int times = 1;

    printf("thread num: %d\n", omp_get_max_threads());

    for (int i = 0; i < times; ++i) {
        t1 = get_wall_time();
        long long ans_our = g->pattern_matching(schedule_our, true); // 进行匹配操作
        t2 = get_wall_time();

        printf("ans: %lld time: %.6lf\n", ans_our, t2 - t1);
        total_time += (t2 - t1);
        if (i == times - 1) {
            schedule_our.print_schedule();
        }
        fflush(stdout);
    }
    total_time /= times;
    printf("Counting time cost: %.6lf s\n", total_time);
    return total_time;
}

int main(int argc, char *argv[]) {
    Graph *g;
    DataLoader D;

    if (argc != 3) {
        printf("usage: %s graph_file clique_size\n", argv[0]);
        return 0;
    }

    DataType type = DataType::Patents;

    // load_data(Graph *&g, DataType type, const char *path, bool binary_input, int oriented_type)
    // bool ok = D.fast_load(g, argv[1]);
    bool ok = D.load_data(g, type, argv[1], 0, 0);
    if (!ok) {
        printf("Load data failed\n");
        return 0;
    }

    printf("Load data success!\n");
    fflush(stdout);

    // 构建团
    int size = atoi(argv[2]);
    Pattern pattern(size);
    for (int i = 0; i < size; i++) {
        for (int j = i + 1; j < size; j++) {
            pattern.add_edge(i, j);
        }
    }

    // 从数据图中移除一些无关的边
    reduce_edges_for_clique(*g);
    // 团匹配
    test_pattern(g, pattern);
    delete g;
    return 0;
}