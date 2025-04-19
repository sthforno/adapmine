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

int main(int argc, char *argv[]) {
    Graph *g;
    DataLoader D;

    if (argc != 3) {
        printf("usage: %s graph_file pattern_size\n", argv[0]);
        return 0;
    }
    DataType type = DataType::Patents;
    // bool ok = D.fast_load(g, argv[1]);
    bool ok = D.load_data(g, type, argv[1], 0, 0);
    if (!ok) {
        printf("Load data failed\n");
        return 0;
    }

    printf("Load data success!\n");
    fflush(stdout);
    int size = atoi(argv[2]);

    printf("thread num: %d\n", omp_get_max_threads());

    // g->triangle_statistics();
    // g->four_clique_statistics();

    if (size == 3)
        g->motif_counting_3();
    else
        g->motif_counting(size);
    delete g;
    return 0;
}
