#undef NDEBUG
#include <common.h>
#include <dataloader.h>
#include <graph.h>
#include <schedule.h>
#include <vertex_set.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <string>

#include <omp.h>
#include <sys/time.h>

// 频繁子图挖掘，处理带标签图，关键在于挖掘超出某个阈值的子图。
int main(int argc, char *argv[]) {
    LabeledGraph *g;
    DataLoader D;

    // assert(argc == 4);
    // const std::string type = argv[1];
    const std::string path = argv[1];
    const int max_edge = atoi(argv[2]);
    const int min_support = atoi(argv[3]);

    DataType my_type;

    GetDataType(my_type, "Patents");

    if (my_type == DataType::Invalid) {
        printf("Dataset not found!\n");
        return 0;
    }

    g = new LabeledGraph();
    assert(D.load_labeled_data(g, my_type, path.c_str()) == true);

    printf("thread num: %d\n", omp_get_max_threads());

    double total_time = 0;
    int times = 1;
    for (int i = 0; i < times; i++) {
        double this_time = 0.0;
        g->fsm(max_edge, min_support, &this_time);
        total_time += this_time;
    }
    total_time /= times;
    printf("Counting time cost: %.6lf s\n", total_time);
    return 0;
}
