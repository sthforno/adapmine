#include "../include/approx.h"
#include "../include/common.h"
#include "../include/graph.h"
#include "../include/motif_generator.h"
#include "../include/schedule_IEP.h"
#include "../include/vertex_set.h"
#include "timeinterval.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <queue>
#include <random>
#include <sstream>
#include <sys/time.h>
#include <unistd.h>

// 计算大部分顶点的度数阈值
int threshhold_computing(Graph *g) {
    int *record_degree = (int *)malloc(sizeof(int) * g->v_cnt);
    memset(record_degree, 0, sizeof(int) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        record_degree[g->vertex[i + 1] - g->vertex[i]]++;
    }
    int threshhold = 1;
    double ratio = 0;
    while (ratio < 0.9) {
        double temp_ratio = double(record_degree[threshhold]) / double(g->v_cnt);
        threshhold++;
        ratio += temp_ratio;
    }

    free(record_degree);
    return threshhold;
}

struct filter {
    uint32_t *edge_project;
    uint32_t *record;
    uint32_t *degree;
    uint32_t edge_size;
};

// 简单无放回采样
void four_clique_egonet_sample(e_index_t process_edge, Graph *g, uint64_t &count) {
    v_index_t left_vtx = g->edge[process_edge];
    v_index_t right_vtx = g->edge_from[process_edge];

    VertexSet third_candidate;
    third_candidate.init();
    // 构建第三个顶点的候选集合
    VertexSet tmp_vset_first;
    tmp_vset_first.init(g->vertex[left_vtx + 1] - g->vertex[left_vtx], g->edge + g->vertex[left_vtx]);
    VertexSet tmp_vset_second;
    tmp_vset_second.init(g->vertex[right_vtx + 1] - g->vertex[right_vtx], g->edge + g->vertex[right_vtx]);
    third_candidate.intersection(tmp_vset_first, tmp_vset_second);

    // 寻找第四个顶点的候选解大小

    int *loop_data_ptr = third_candidate.get_data_ptr();
    for (int i = 0; i < third_candidate.get_size(); i++) {
        int vertex_third = loop_data_ptr[i];
        if (vertex_third < right_vtx || g->vertex[vertex_third + 1] - g->vertex[vertex_third] < 3) {
            continue;
        }
        VertexSet tmp_vtx_forth;
        tmp_vtx_forth.init(g->vertex[vertex_third + 1] - g->vertex[vertex_third], g->edge + g->vertex[vertex_third]);
        VertexSet forth_candidate;
        forth_candidate.init();
        forth_candidate.intersection(tmp_vtx_forth, third_candidate);
        count += forth_candidate.get_size();
    }
}

void four_diamond_egonet_sample_sub_neighbo(e_index_t process_edge, Graph *g, uint64_t &count) {
    v_index_t left_vtx = g->edge[process_edge];
    v_index_t right_vtx = g->edge_from[process_edge];

    VertexSet third_candidate;
    third_candidate.init();
    // 构建第三个顶点的候选集合
    VertexSet tmp_vset_first;
    tmp_vset_first.init(g->vertex[left_vtx + 1] - g->vertex[left_vtx], g->edge + g->vertex[left_vtx]);
    VertexSet tmp_vset_second;
    tmp_vset_second.init(g->vertex[right_vtx + 1] - g->vertex[right_vtx], g->edge + g->vertex[right_vtx]);
    third_candidate.intersection(tmp_vset_first, tmp_vset_second);

    // 寻找第四个顶点的候选解大小

    int *loop_data_ptr = third_candidate.get_data_ptr();
    for (int i = 0; i < third_candidate.get_size(); i++) {
        // 处理每一个第三个顶点
        int vertex_third = loop_data_ptr[i];
        VertexSet tmp_vtx_forth;
        tmp_vtx_forth.init(g->vertex[vertex_third + 1] - g->vertex[vertex_third], g->edge + g->vertex[vertex_third]);

        VertexSet forth_candidate;
        forth_candidate.init();
        // 第四个顶点为第三个顶点候选集的大小，减去第三个顶点的邻接边表
        int forth_candidate_size = forth_candidate.unordered_subtraction_size(third_candidate, tmp_vtx_forth);
        count += forth_candidate_size;
    }
}

void four_diamond_egonet_sample_sub_cur(e_index_t process_edge, Graph *g, uint64_t &count) {
    v_index_t left_vtx = g->edge[process_edge];
    v_index_t right_vtx = g->edge_from[process_edge];

    VertexSet third_candidate;
    third_candidate.init();
    // 构建第三个顶点的候选集合
    VertexSet tmp_vset_first;
    tmp_vset_first.init(g->vertex[left_vtx + 1] - g->vertex[left_vtx], g->edge + g->vertex[left_vtx]);
    VertexSet tmp_vset_second;
    tmp_vset_second.init(g->vertex[right_vtx + 1] - g->vertex[right_vtx], g->edge + g->vertex[right_vtx]);
    third_candidate.intersection(tmp_vset_first, tmp_vset_second);

    // 寻找第四个顶点的候选解大小
    int num = third_candidate.get_size();

    count += (num * (num - 1) / 2);
}

// 混合采样
void four_clique_egonet_sample_hybird(e_index_t process_edge, Graph *g, uint64_t &count, uint64_t &sample_count, uint32_t &sample_times,
                                      double random_num) {
    v_index_t left_vtx = g->edge[process_edge];
    v_index_t right_vtx = g->edge_from[process_edge];

    VertexSet third_candidate;
    third_candidate.init();
    // 构建第三个顶点的候选集合
    VertexSet tmp_vset_first;
    tmp_vset_first.init(g->vertex[left_vtx + 1] - g->vertex[left_vtx], g->edge + g->vertex[left_vtx]);
    VertexSet tmp_vset_second;
    tmp_vset_second.init(g->vertex[right_vtx + 1] - g->vertex[right_vtx], g->edge + g->vertex[right_vtx]);
    third_candidate.intersection(tmp_vset_first, tmp_vset_second);

    if (third_candidate.get_size() < 10) {
        // 寻找第四个顶点的候选解大小
        int *loop_data_ptr = third_candidate.get_data_ptr();
        for (int i = 0; i < third_candidate.get_size(); i++) {
            int vertex_third = loop_data_ptr[i];
            if (vertex_third < right_vtx) {
                continue;
            }
            VertexSet tmp_vtx_forth;
            tmp_vtx_forth.init(g->vertex[vertex_third + 1] - g->vertex[vertex_third], g->edge + g->vertex[vertex_third]);
            VertexSet forth_candidate;
            forth_candidate.init();
            forth_candidate.intersection(tmp_vtx_forth, third_candidate);
            count += forth_candidate.get_size();
        }
    } else {
        int *loop_data_ptr = third_candidate.get_data_ptr();
        int vertex_third = loop_data_ptr[int(random_num * third_candidate.get_size())];
        VertexSet tmp_vtx_forth;
        tmp_vtx_forth.init(g->vertex[vertex_third + 1] - g->vertex[vertex_third], g->edge + g->vertex[vertex_third]);
        VertexSet forth_candidate;
        forth_candidate.init();
        forth_candidate.intersection(tmp_vtx_forth, third_candidate);
        sample_count += forth_candidate.get_size() * third_candidate.get_size();
        sample_times++;
    }
}

void four_clique_egonet_record(e_index_t process_edge, Graph *g, uint64_t &count, filter *filter_data) {
    v_index_t left_vtx = g->edge[process_edge];
    v_index_t right_vtx = g->edge_from[process_edge];

    VertexSet third_candidate;
    third_candidate.init();
    // 构建第三个顶点的候选集合
    VertexSet tmp_vset_first;
    tmp_vset_first.init(g->vertex[left_vtx + 1] - g->vertex[left_vtx], g->edge + g->vertex[left_vtx]);
    VertexSet tmp_vset_second;
    tmp_vset_second.init(g->vertex[right_vtx + 1] - g->vertex[right_vtx], g->edge + g->vertex[right_vtx]);
    third_candidate.intersection(tmp_vset_first, tmp_vset_second);

    // 寻找第四个顶点的候选解大小

    int *loop_data_ptr = third_candidate.get_data_ptr();
    for (int i = 0; i < third_candidate.get_size(); i++) {
        int vertex_third = loop_data_ptr[i];
        if (vertex_third < right_vtx || filter_data->degree[vertex_third] < 3) {
            continue;
        }
        VertexSet tmp_vtx_forth;
        tmp_vtx_forth.init(g->vertex[vertex_third + 1] - g->vertex[vertex_third], g->edge + g->vertex[vertex_third]);
        VertexSet forth_candidate;
        forth_candidate.init();
        forth_candidate.intersection(tmp_vtx_forth, third_candidate);
        count += forth_candidate.get_size();
    }
    filter_data->record[process_edge] = 1;
}

void four_clique_egonet_record_hybird(e_index_t process_edge, Graph *g, uint64_t &count, filter *filter_data) {
    v_index_t left_vtx = g->edge[process_edge];
    v_index_t right_vtx = g->edge_from[process_edge];

    VertexSet third_candidate;
    third_candidate.init();
    // 构建第三个顶点的候选集合
    VertexSet tmp_vset_first;
    tmp_vset_first.init(g->vertex[left_vtx + 1] - g->vertex[left_vtx], g->edge + g->vertex[left_vtx]);
    VertexSet tmp_vset_second;
    tmp_vset_second.init(g->vertex[right_vtx + 1] - g->vertex[right_vtx], g->edge + g->vertex[right_vtx]);
    third_candidate.intersection(tmp_vset_first, tmp_vset_second);

    // 寻找第四个顶点的候选解大小

    int *loop_data_ptr = third_candidate.get_data_ptr();
    for (int i = 0; i < third_candidate.get_size(); i++) {
        int vertex_third = loop_data_ptr[i];
        if (vertex_third < right_vtx) {
            continue;
        }

        VertexSet tmp_vtx_forth;
        tmp_vtx_forth.init(g->vertex[vertex_third + 1] - g->vertex[vertex_third], g->edge + g->vertex[vertex_third]);
        VertexSet forth_candidate;
        forth_candidate.init();
        forth_candidate.intersection(tmp_vtx_forth, third_candidate);
        count += forth_candidate.get_size();
    }
    filter_data->record[process_edge] = 1;
}

void four_clique_egonet_filter(e_index_t process_edge, Graph *g, uint64_t &count, filter *filter_data) {

    e_index_t process_edge_locate = process_edge;
    while (filter_data->record[process_edge] != -1) { // 直到找到一条没有被处理过的边
        process_edge = filter_data->edge_project[process_edge];
    }

    v_index_t left_vtx = g->edge[process_edge];
    v_index_t right_vtx = g->edge_from[process_edge];

    VertexSet third_candidate;
    third_candidate.init();
    // 构建第三个顶点的候选集合
    VertexSet tmp_vset_first;
    tmp_vset_first.init(g->vertex[left_vtx + 1] - g->vertex[left_vtx], g->edge + g->vertex[left_vtx]);
    VertexSet tmp_vset_second;
    tmp_vset_second.init(g->vertex[right_vtx + 1] - g->vertex[right_vtx], g->edge + g->vertex[right_vtx]);
    third_candidate.intersection(tmp_vset_first, tmp_vset_second);

    // 寻找第四个顶点的候选解大小

    int *loop_data_ptr = third_candidate.get_data_ptr();
    for (int i = 0; i < third_candidate.get_size(); i++) {
        int vertex_third = loop_data_ptr[i];
        VertexSet tmp_vtx_forth;
        tmp_vtx_forth.init(g->vertex[vertex_third + 1] - g->vertex[vertex_third], g->edge + g->vertex[vertex_third]);
        VertexSet forth_candidate;
        forth_candidate.init();
        forth_candidate.intersection(tmp_vtx_forth, third_candidate);
        count += forth_candidate.get_size();
    }
    // 更新投射
    filter_data->record[process_edge] = -1; // 更新实际边状态
    if (process_edge_locate != filter_data->edge_size)
        filter_data->edge_project[process_edge_locate] = filter_data->edge_project[filter_data->edge_size]; // 将当前位置记录最后一条边所在的位置
    filter_data->edge_size--;                                                                               // 更新剩余边数量
}

//==================================================================================================================================================//
//============精确算法=============================================================================================================================//
//==================================================================================================================================================//
void demo_pattern_test_exact(Graph *g) {
    uint64_t count = 0;
    g->edge_from = (v_index_t *)malloc(sizeof(v_index_t) * g->e_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        for (int j = g->vertex[i]; j < g->vertex[i + 1]; j++) {
            g->edge_from[j] = i;
        }
    }

    for (e_index_t i = 0; i < g->e_cnt; ++i) {
        v_index_t vtx_left = g->edge_from[i];
        v_index_t vtx_right = g->edge[i];

        if (vtx_left > vtx_right) {
            continue;
        }
        // four_clique_egonet_sample(i, g, count);
        // four_diamond_egonet_sample_sub_neighbo(i, g, count); //
        four_diamond_egonet_sample_sub_cur(i, g, count); // 251755062 和sc答案一样
    }

    std::cout << "钻石计数为 " << count << std::endl;
}

//==================================================================================================================================================//
//==========有放回采样算法============================================================================================================================//
//==================================================================================================================================================//
void demo_pattern_test_sample(Graph *g) {

    double estimate[10] = {0};
    double error[10] = {0};
    uint64_t count = 0;
    g->edge_from = (v_index_t *)malloc(sizeof(v_index_t) * g->e_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        for (int j = g->vertex[i]; j < g->vertex[i + 1]; j++) {
            g->edge_from[j] = i;
        }
    }

    std::uniform_int_distribution<e_index_t> edge_dist(0, g->e_cnt - 1);
    std::random_device rd;
    std::default_random_engine gen(rd());
    int record = 1;
    int sample_time = 0;
    while (record < 7) {
        auto eid = edge_dist(gen);
        if (g->edge_from[eid] > g->edge[eid]) {
            continue;
        }
        four_clique_egonet_sample(eid, g, count);
        sample_time++;
        if (std::pow(10, record) == sample_time) {
            record++;
            estimate[record] = count * g->e_cnt / (2 * sample_time);
        }
    }

    for (int i = 1; i < 7; i++) {
        std::cout << "四团计数为 " << estimate[i] << std::endl;
        error[i] = (estimate[i] - (148834439.0 * 8.0)) / (148834439.0 * 8.0);
        std::cout << "误差为 " << error[i] << std::endl;
    }
}
// remark: 采样数一千误差百分白，一万误差百分之十八，十万误差百分之五，百万误差百分之一。  一共五百万边
// remark: sk 去掉一次对称性。  十：1 百：3 千：0.2 万：0.2 十万：0.06 百万：0.03     运行时间90秒
// remark: sk 去掉两次对称性。  十：1 百：0.91 千：0.44 万：0.64 十万：0.02 百万：0.0003     运行时间41秒
// remark: sk 去掉两次对称性。  十：1 百：0.70 千：0.26 万：0.16 十万：0.21 百万：0.07     运行时间41秒
// remark: sk 去掉两次对称性。  十：1 百：0.98 千：0.82 万：0.16 十万：0.08 百万：0.01     运行时间47秒
//==================================================================================================================================================//
//==========过滤后有放回采样=============================================================================================================================//
//==================================================================================================================================================//
void demo_pattern_test_filter_sample(Graph *g) {
    // 第一步，在稀疏部分过滤一遍
    // 可以考虑两种稀疏化方法，处理边的两个端点，以及处理边两个端点的端点

    // 稠密稀疏的threshhould的计算和选择
    int dense_threshold = g->e_cnt / g->v_cnt; // 稠密度数的阈值      frendster 28 twitter 35 livejouranl 17 youtube 5 patent 5 pokec 30 tp 28 sk 11
    int sparse_threshold = 0;
    sparse_threshold = threshhold_computing(g); // 稀疏度数的阈值
    // int threshold = std::sqrt(g->v_cnt);
    // int threshold = dense_threshold;
    int threshold = g->v_cnt / 10000;

    uint64_t count = 0;
    g->edge_from = (v_index_t *)malloc(sizeof(v_index_t) * g->e_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        for (int j = g->vertex[i]; j < g->vertex[i + 1]; j++) {
            g->edge_from[j] = i;
        }
    }

    // 过滤并记录剩余边
    filter filter_data;
    filter_data.edge_project = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
    filter_data.record = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
    filter_data.edge_size = 0;
    memset(filter_data.edge_project, -1, sizeof(uint32_t) * g->e_cnt);
    memset(filter_data.record, -1, sizeof(uint32_t) * g->e_cnt);

    filter_data.degree = (uint32_t *)malloc(sizeof(uint32_t) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        filter_data.degree[i] = g->vertex[i + 1] - g->vertex[i];
    }

    for (e_index_t i = 0; i < g->e_cnt; ++i) {
        v_index_t vtx_left = g->edge_from[i];
        v_index_t vtx_right = g->edge[i];
        if (filter_data.degree[vtx_left] < threshold && filter_data.degree[vtx_right] < threshold && vtx_left < vtx_right) {
            four_clique_egonet_record(i, g, count, &filter_data);
        }
    }
    for (e_index_t i = 0; i < g->e_cnt; ++i) {
        if (filter_data.record[i] == -1 && g->edge_from[i] < g->edge[i]) {
            filter_data.edge_project[filter_data.edge_size++] = i;
        }
    }

    // 随机采样
    uint64_t count_sample = 0;
    double estimate[10] = {0};
    double error[10] = {0};
    std::uniform_int_distribution<e_index_t> edge_dist(0, filter_data.edge_size - 1);
    std::random_device rd;
    std::default_random_engine gen(rd());
    int record = 1;
    int processedge;
    for (int i = 0; i <= 1000000; i++) {
        auto eid = edge_dist(gen);
        processedge = filter_data.edge_project[eid];
        four_clique_egonet_sample(processedge, g, count_sample);
        if (std::pow(10, record) == i) {
            estimate[record] = count_sample * filter_data.edge_size / i;
            record++;
        }
    }

    for (int i = 1; i < 7; i++) {
        std::cout << "四团计数为 " << estimate[i] + count << std::endl;
        error[i] = (estimate[i] + double(count) - (4986965.0 * 8.0)) / (4986965.0 * 8.0); // 4986965 148834439
        std::cout << "误差为 " << error[i] << std::endl;
    }

    free(filter_data.edge_project);
    free(filter_data.record);
    free(filter_data.degree);
}

// remark 百：100， 千：25，万：4，十万：0.9，百万误差百分之0.06   使用平方剩余边45万
// remark: sk 去掉一次对称性 十0.45 百0.16 千：0.07，万：0.03，十万 0.006，百万：0.004  运行时间442秒
// remark: sk 去掉两次对称性 十0.73 百0.75 千：0.43，万：0.12，十万 0.008，百万：0.002  运行时间50秒  15为threshhold 边数减少很少
// remark: sk 去掉两次对称性 十0.94 百0.62 千：0.98，万：0.35，十万 0.03，百万：0.009  运行时间50秒  15为threshhold 边数减少很少
// remark: sk 去掉两次对称性 十0.33 百1.44 千：0.12，万：0.13，十万 0.01，百万：0.006  运行时间67秒  100为threshhold
// remark: sk 去掉两次对称性 十0.81 百1.12 千：0.04，万：0.05，十万 0.07，百万：0.01  运行时间69秒  100为threshhold
// remark: sk 去掉两次对称性 十0.75 百0.41 千：0.07，万：0.08，十万 0.02，百万：0.007  运行时间87秒  100为threshhold
// remark: sk 去掉两次对称性 十2.67 百0.63 千：0.02，万：0.04，十万 0.02，百万：0.009  运行时间217秒  1400为threshhold

//==================================================================================================================================================//
//============过滤后无放回采样算法====================================================================================================================//
//==================================================================================================================================================//
void demo_pattern_test_filter_reduce_one(Graph *g) {
    // 第一步，在稀疏部分过滤一遍
    // 可以考虑两种稀疏化方法，处理边的两个端点，以及处理边两个端点的端点

    // 稠密稀疏的threshhould的计算和选择
    int dense_threshold = g->e_cnt / g->v_cnt; // 稠密度数的阈值      frendster 28 twitter 35 livejouranl 17 youtube 2 patent 5 pokec 30 tp 28 sk 11
    int sparse_threshold = 0;
    sparse_threshold = threshhold_computing(g); // 稀疏度数的阈值
    int threshold = std::sqrt(g->v_cnt);

    uint64_t count = 0;
    g->edge_from = (v_index_t *)malloc(sizeof(v_index_t) * g->e_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        for (int j = g->vertex[i]; j < g->vertex[i + 1]; j++) {
            g->edge_from[j] = i;
        }
    }

    // 过滤并记录剩余边
    filter filter_data;
    filter_data.edge_project = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
    filter_data.record = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
    filter_data.edge_size = 0;
    memset(filter_data.edge_project, -1, sizeof(uint32_t) * g->e_cnt);
    memset(filter_data.record, -1, sizeof(uint32_t) * g->e_cnt);

    int32_t *degree = (int32_t *)malloc(sizeof(int32_t) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        degree[i] = g->vertex[i + 1] - g->vertex[i];
    }

    for (e_index_t i = 0; i < g->e_cnt; ++i) {
        v_index_t vtx_left = g->edge_from[i];
        v_index_t vtx_right = g->edge[i];
        if (degree[vtx_left] < threshold && degree[vtx_right] < threshold && vtx_left < vtx_right) {
            four_clique_egonet_record(i, g, count, &filter_data);
        }
    }
    for (e_index_t i = 0; i < g->e_cnt; ++i) {
        if (filter_data.record[i] == -1 && g->edge_from[i] < g->edge[i]) {
            filter_data.edge_project[filter_data.edge_size++] = i;
        }
    }

    // 随机采样
    std::random_device rd;
    std::default_random_engine gen(rd());
    uint32_t *order = (uint32_t *)malloc(sizeof(uint32_t) * filter_data.edge_size);
    for (uint32_t i = 0; i < filter_data.edge_size; i++) {
        order[i] = i;
    }
    std::shuffle(order, order + filter_data.edge_size, gen);

    int record = 1;
    uint64_t count_sample = 0;
    double estimate[10] = {0};
    double error[10] = {0};
    for (int i = 0; i < filter_data.edge_size; i++) {
        auto eid = filter_data.edge_project[order[i]];
        four_clique_egonet_sample(eid, g, count_sample);
        if (std::pow(10, record) == i) {
            estimate[record] = count_sample * filter_data.edge_size / i;
            record++;
        }
    }
    estimate[record] = count_sample;
    for (int i = 1; i <= record; i++) {
        std::cout << "四团计数为 " << estimate[i] + count << std::endl;
        error[i] = (estimate[i] + double(count) - (148834439.0 * 8.0)) / (148834439.0 * 8.0);
        std::cout << "误差为 " << error[i] << std::endl;
    }

    free(filter_data.edge_project);
    free(filter_data.record);
    free(degree);
}

// remark sk 十：0.46 百：0.46 千：0.19 万：0.01 十万：0.001 百万：0  470秒
// remark sk 十：0.32 百：0.24 千：0.03 万：0.01 十万：0.001 百万： 0  470秒
// 采用度数剪枝
void demo_pattern_test_filter_reduce_two(Graph *g) {
    // 第一步，在稀疏部分过滤一遍
    // 可以考虑两种稀疏化方法，处理边的两个端点，以及处理边两个端点的端点

    // 稠密稀疏的threshhould的计算和选择
    int dense_threshold = g->e_cnt / g->v_cnt; // 稠密度数的阈值      frendster 28 twitter 35 livejouranl 17 youtube 2 patent 5 pokec 30 tp 28 sk 11
    int sparse_threshold = 0;
    sparse_threshold = threshhold_computing(g); // 稀疏度数的阈值
    int threshold = std::sqrt(g->v_cnt);

    uint64_t count = 0;
    g->edge_from = (v_index_t *)malloc(sizeof(v_index_t) * g->e_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        for (int j = g->vertex[i]; j < g->vertex[i + 1]; j++) {
            g->edge_from[j] = i;
        }
    }

    // 过滤并记录剩余边
    filter filter_data;
    filter_data.edge_project = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
    filter_data.record = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
    filter_data.edge_size = 0;
    memset(filter_data.edge_project, -1, sizeof(uint32_t) * g->e_cnt);
    memset(filter_data.record, -1, sizeof(uint32_t) * g->e_cnt);

    filter_data.degree = (uint32_t *)malloc(sizeof(uint32_t) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        filter_data.degree[i] = g->vertex[i + 1] - g->vertex[i];
    }

    for (e_index_t i = 0; i < g->e_cnt; ++i) {
        v_index_t vtx_left = g->edge_from[i];
        v_index_t vtx_right = g->edge[i];
        if (filter_data.degree[vtx_left] < threshold && filter_data.degree[vtx_right] < threshold && vtx_left < vtx_right) {
            if (filter_data.degree[vtx_left] < 3 || filter_data.degree[vtx_right] < 3) {
                filter_data.record[i] = 1;
            }
            four_clique_egonet_record(i, g, count, &filter_data);
        }
    }
    for (e_index_t i = 0; i < g->e_cnt; ++i) {
        if (filter_data.record[i] == -1 && g->edge_from[i] < g->edge[i]) {
            filter_data.edge_project[filter_data.edge_size++] = i;
        }
    }

    // 随机采样
    std::random_device rd;
    std::default_random_engine gen(rd());
    uint32_t *order = (uint32_t *)malloc(sizeof(uint32_t) * filter_data.edge_size);
    for (uint32_t i = 0; i < filter_data.edge_size; i++) {
        order[i] = i;
    }
    std::shuffle(order, order + filter_data.edge_size, gen);

    int record = 1;
    uint64_t count_sample = 0;
    double estimate[10] = {0};
    double error[10] = {0};
    for (int i = 0; i < filter_data.edge_size; i++) {
        auto eid = filter_data.edge_project[order[i]];
        four_clique_egonet_sample(eid, g, count_sample);
        if (std::pow(10, record) == i) {
            estimate[record] = count_sample * filter_data.edge_size / i;
            record++;
        }
    }
    estimate[record] = count_sample;
    for (int i = 1; i <= record; i++) {
        std::cout << "四团计数为 " << estimate[i] + count << std::endl;
        error[i] = (estimate[i] + double(count) - (148834439.0 * 8.0)) / (148834439.0 * 8.0);
        std::cout << "误差为 " << error[i] << std::endl;
    }

    free(filter_data.edge_project);
    free(filter_data.record);
    free(filter_data.degree);
}

// remark sk 十：0.46 百：0.46 千：0.19 万：0.01 十万：0.001 百万：0  470秒  无剪枝
// remark sk 十：0.32 百：0.24 千：0.03 万：0.01 十万：0.001 百万： 0  470秒  无剪枝
// remark 加上剪枝560秒

// 采用了不同的预测方式
void demo_pattern_test_filter_reduce_three(Graph *g) {
    // 稠密稀疏的threshhould的计算和选择
    int dense_threshold = g->e_cnt / g->v_cnt; // 稠密度数的阈值      frendster 28 twitter 35 livejouranl 17 youtube 2 patent 5 pokec 30 tp 28 sk 11
    int sparse_threshold = 0;
    sparse_threshold = threshhold_computing(g); // 稀疏度数的阈值
    int threshold = std::sqrt(g->v_cnt);

    uint64_t count = 0;
    g->edge_from = (v_index_t *)malloc(sizeof(v_index_t) * g->e_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        for (int j = g->vertex[i]; j < g->vertex[i + 1]; j++) {
            g->edge_from[j] = i;
        }
    }

    // 过滤并记录剩余边
    filter filter_data;
    filter_data.edge_project = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
    filter_data.record = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
    filter_data.edge_size = 0;
    memset(filter_data.edge_project, -1, sizeof(uint32_t) * g->e_cnt);
    memset(filter_data.record, -1, sizeof(uint32_t) * g->e_cnt);

    int32_t *degree = (int32_t *)malloc(sizeof(int32_t) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        degree[i] = g->vertex[i + 1] - g->vertex[i];
    }

    for (e_index_t i = 0; i < g->e_cnt; ++i) {
        v_index_t vtx_left = g->edge_from[i];
        v_index_t vtx_right = g->edge[i];
        if (degree[vtx_left] < threshold && degree[vtx_right] < threshold && vtx_left < vtx_right) {
            four_clique_egonet_record(i, g, count, &filter_data);
        }
    }
    for (e_index_t i = 0; i < g->e_cnt; ++i) {
        if (filter_data.record[i] == -1 && g->edge_from[i] < g->edge[i]) {
            filter_data.edge_project[filter_data.edge_size++] = i;
        }
    }

    // 随机采样
    std::random_device rd;
    std::default_random_engine gen(rd());
    uint32_t *order = (uint32_t *)malloc(sizeof(uint32_t) * filter_data.edge_size);
    for (uint32_t i = 0; i < filter_data.edge_size; i++) {
        order[i] = i;
    }
    std::shuffle(order, order + filter_data.edge_size, gen);

    int record = 1;
    uint64_t count_sample = 0;
    uint64_t count_sample_total = 0;
    uint64_t count_sample_temp = 0;
    double estimate[10] = {0};
    double error[10] = {0};
    for (int i = 0; i < filter_data.edge_size; i++) {
        auto eid = filter_data.edge_project[order[i]];
        count_sample_temp = count_sample;
        four_clique_egonet_sample(eid, g, count_sample);
        count_sample_total += ((count_sample - count_sample_temp) * (filter_data.edge_size - i) + count_sample); // 当前采样预测整体
        if (std::pow(10, record) == i) {
            estimate[record] = count_sample_total / i;
            record++;
        }
    }
    estimate[record] = count_sample_total / filter_data.edge_size;
    for (int i = 1; i <= record; i++) {
        std::cout << "四团计数为 " << estimate[i] + count << std::endl;
        error[i] = (estimate[i] + double(count) - (148834439.0 * 8.0)) / (148834439.0 * 8.0);
        std::cout << "误差为 " << error[i] << std::endl;
    }

    free(filter_data.edge_project);
    free(filter_data.record);
    free(degree);
}
// remark： 预测误差最后可以收敛为0 百：25 千：3  万：0.2  十万：0.006  百万：0
// remark：sk 十：0.52 百：0.34 千：0.22 万：0.0009 十万：0.008 百万：0.002  482秒
// remark：sk 十：0.51 百：0.52 千：0.19 万：0.09 十万：0.06 百万：0.008  482秒
void demo_pattern_test_filter_reduce_hybird(Graph *g) {
    // 第一步，在稀疏部分过滤一遍
    // 可以考虑两种稀疏化方法，处理边的两个端点，以及处理边两个端点的端点

    // 稠密稀疏的threshhould的计算和选择
    int dense_threshold = g->e_cnt / g->v_cnt; // 稠密度数的阈值      frendster 28 twitter 35 livejouranl 17 youtube 2 patent 5 pokec 30 tp 28 sk 11
    int sparse_threshold = 0;
    sparse_threshold = threshhold_computing(g); // 稀疏度数的阈值
    int threshold = g->v_cnt / 10000;

    uint64_t count = 0;
    g->edge_from = (v_index_t *)malloc(sizeof(v_index_t) * g->e_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        for (int j = g->vertex[i]; j < g->vertex[i + 1]; j++) {
            g->edge_from[j] = i;
        }
    }

    // 过滤并记录剩余边
    filter filter_data;
    filter_data.edge_project = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
    filter_data.record = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
    filter_data.edge_size = 0;
    memset(filter_data.edge_project, -1, sizeof(uint32_t) * g->e_cnt);
    memset(filter_data.record, -1, sizeof(uint32_t) * g->e_cnt);

    int32_t *degree = (int32_t *)malloc(sizeof(int32_t) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        degree[i] = g->vertex[i + 1] - g->vertex[i];
    }

    for (e_index_t i = 0; i < g->e_cnt; ++i) {
        v_index_t vtx_left = g->edge_from[i];
        v_index_t vtx_right = g->edge[i];
        if (degree[vtx_left] < threshold && degree[vtx_right] < threshold && vtx_left < vtx_right) {
            four_clique_egonet_record_hybird(i, g, count, &filter_data);
        }
    }
    for (e_index_t i = 0; i < g->e_cnt; ++i) {
        if (filter_data.record[i] == -1 && g->edge_from[i] < g->edge[i]) {
            filter_data.edge_project[filter_data.edge_size++] = i;
        }
    }

    // 随机采样
    std::random_device rd;
    std::default_random_engine gen(rd());
    std::uniform_int_distribution<uint32_t> edge_dist(0, filter_data.edge_size - 1);
    uint32_t *order = (uint32_t *)malloc(sizeof(uint32_t) * filter_data.edge_size);
    for (uint32_t i = 0; i < filter_data.edge_size; i++) {
        order[i] = i;
    }
    std::shuffle(order, order + filter_data.edge_size, gen);

    int record = 1;
    uint64_t sample_count = 0;
    uint32_t sample_times = 0;
    double estimate[10] = {0};
    double error[10] = {0};
    for (int i = 0; i < filter_data.edge_size; i++) {
        auto eid = filter_data.edge_project[order[i]];
        double random_num = double(edge_dist(gen)) / double(filter_data.edge_size);
        four_clique_egonet_sample_hybird(eid, g, count, sample_count, sample_times, random_num);
        if (std::pow(10, record) == sample_times) {
            estimate[record] = sample_count / sample_times;
            record++;
        }
    }
    estimate[record] = sample_count / sample_times;
    for (int i = 1; i <= record; i++) {
        std::cout << "四团计数为 " << estimate[i] + count << std::endl;
        error[i] = ((estimate[i] * sample_times * 2.0 / 3.0) + double(count) - (148834439.0 * 8.0)) / (148834439.0 * 8.0);
        std::cout << "误差为 " << error[i] << std::endl;
    }

    free(filter_data.edge_project);
    free(filter_data.record);
    free(degree);
}
// remark: sk, threshhold 169, 一轮过滤剩余五百万边（缩减约一半），500秒
// 十，0.7 百，0.16 千，0.42 万，0.05 十万，0.006 百万，0.003 五百万 0.001

// 混合采样简单实现280秒 对称性打破（都只使用一种）
// th 候选集大小为10 270秒 十，0.79 百，0.62 千，0.18 万，0.09 十万，0.04 百万，0.02 全部 0.006。

// 混合采样混合预测197秒 对称性打破（精确算法和采样使用不同的预测方式）
//  十，0.51 百，0.08 千，0.17 万，0.06 十万，0.06 百万，0.07 全部 0.003。
//  十，0.89 百，0.23 千，0.11 万，0.01 十万，0.02 百万，0.02 全部 0.02。

//==================================================================================================================================================//
//========过滤+无放回+es和ns的混合====================================================================================================================//
//==================================================================================================================================================//

//==================================================================================================================================================//
//===========调度方案构建============================================================================================================================//
//==================================================================================================================================================//

void scheduel_motif_generate_approx(Pattern P, uint32_t pattern_size) {
    MotifGenerator mg(pattern_size);
    std::vector<Pattern> motifs = mg.generate();
    for (int i = 0; i < motifs.size(); ++i) {
        Pattern p = motifs[i];
    }
}

// 输入：模式的邻接矩阵和大小
// 输出：调度方案

// 要求：1，匹配顺序+依赖关系 2，剪枝（越稠密越好）3，使用数值计算代替循环计算
// 4，限制（打破对称性，越早越好）5，交集提前计算（暗示哪些交集可以提前计算）
// 6，随机选择友好  7降低树的高度（等价于增加不想关的顶点的数量IEP？），

// 邻接矩阵对称。
struct schedule_approx {
    // 基本功能  1，匹配顺序 2，依赖关系 3，简单剪枝（高度顶点优先潜在的实现了顶点覆盖集优先和降低匹配树高度的功能） 4，数值计算代替循环计算
    // //现有实现依据，稠密的在前面以提前剪枝，与以匹配顶点更加相关的在前面，交集操作提前(与现有匹配顶点依赖关系最大的顶点优先匹配)以缩小范围,
    // 高度顶点放在后面的话，可能要执行多次 多顶点交集操作
    // (现有实现考虑的是局部最优而不是全局最优，使用性能模型可以？)
    uint32_t *order;            // 匹配顺序
    uint32_t *order_degree;     // 匹配顺序中顶点的度
    uint32_t *dependcy;         // 依赖的前面的顶点
    uint32_t *dependcy_ptr;     // 依赖的指针
    uint32_t compute_depth = 0; // 记录可以数值计算的深度
    uint32_t pattern_size;      // 记录模式的大小

    bool symmetry_label; // 为ns构建约束打破对称性

    // 为es构建约束打破对称性

    // 待实现功能   1，为es施加约束打破对称性，同时满足随机选择要求  2，a，复用计算结果 b，将计算提到循环之前（只有es需要）。
    uint32_t *cache_insetion_recorder; // 记录当前轮次可以预先计算哪些后续需要使用的结果
};

void restricts_generate_approx(const int *cur_adj_mat, std::vector<std::vector<std::pair<int, int>>> &restricts, int pattern_size) {
    Schedule_IEP schedule(cur_adj_mat, pattern_size);
    schedule.aggressive_optimize_get_all_pairs(restricts);
    int size = schedule.get_size();
    Graph *complete;
    DataLoader *D = new DataLoader();
    assert(D->load_data(complete, size + 1));
    long long ans = complete->pattern_matching(schedule, 1) / schedule.get_multiplicity();
    int thread_num = 1;
    for (int i = 0; i < restricts.size();) {
        Schedule_IEP cur_schedule(schedule.get_adj_mat_ptr(), schedule.get_size());
        cur_schedule.add_restrict(restricts[i]);
        long long cur_ans = complete->pattern_matching(cur_schedule, thread_num);
        if (cur_ans != ans) {
            restricts.erase(restricts.begin() + i);
        } else {
            ++i;
        }
    }

    delete complete;
    delete D;
}
void scheduel_generate_approx_for_ns(Pattern *P, schedule_approx *schedule) {
    // 先对矩阵进行排序
    int size = P->get_size();
    int *adj_mat = new int[size * size];
    memcpy(adj_mat, P->get_adj_mat_ptr(), size * size * sizeof(int));

    // 检测对称性打破约束
    std::vector<std::vector<std::pair<int, int>>> restricts_vector;
    restricts_vector.clear();
    restricts_generate_approx(adj_mat, restricts_vector, size);
    for (int i = 0; i < restricts_vector.size(); ++i) {
        for (int j = 0; j < restricts_vector[i].size(); ++j) {
            printf("%d %d ", restricts_vector[i][j].first, restricts_vector[i][j].second);
        }
        printf("\n");
    }

    int *vtx_ptr = new int[size + 2];
    int *vtx_end = new int[size * size];
    int *degree = new int[size];
    vtx_ptr[0] = 0;
    int iter_vtx = 1;
    vtx_ptr[iter_vtx] = vtx_ptr[0];
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (adj_mat[INDEX(i, j, size)]) {
                vtx_end[vtx_ptr[iter_vtx]++] = j;
            }
        }
        vtx_ptr[iter_vtx + 1] = vtx_ptr[iter_vtx];
        iter_vtx++;
    }
    for (int i = 0; i < size; i++) {
        degree[i] = vtx_ptr[i + 1] - vtx_ptr[i];
    }

    // 初始化调度方案
    schedule->order = new uint32_t[size];
    memset(schedule->order, 0, sizeof(uint32_t) * size);
    schedule->order_degree = new uint32_t[size];
    schedule->dependcy = new uint32_t[size * size];
    schedule->dependcy_ptr = new uint32_t[size];
    schedule->compute_depth = size;
    schedule->pattern_size = size;
    schedule->cache_insetion_recorder = new uint32_t[size];
    if (restricts_vector.size() == 0) {
        schedule->symmetry_label = false;
    } else {
        schedule->symmetry_label = true;
    }

    // 优先处理约束连接边
    int *already_assigned = new int[size];
    memset(already_assigned, 0, sizeof(int) * size);
    int max_weight = 0;
    if (schedule->symmetry_label) {
        if (restricts_vector.size() == 1) {
            if (restricts_vector[0].size() == 1 &&
                adj_mat[INDEX(restricts_vector[0][0].first, restricts_vector[0][0].second, size)]) { // 只有一对约束并且相连
                schedule->order[0] = restricts_vector[0][0].first;
                schedule->order[1] = restricts_vector[0][0].second;
            } else {
                for (int i = 0; i < restricts_vector[0].size(); i++) { // 查找相邻的顶点对中边数最大的
                    if (degree[restricts_vector[0][i].first] + degree[restricts_vector[0][i].second] > max_weight &&
                        adj_mat[INDEX(restricts_vector[0][i].first, restricts_vector[0][i].second, size)]) {
                        max_weight = degree[restricts_vector[0][i].first] + degree[restricts_vector[0][i].second];
                        schedule->order[0] = restricts_vector[0][i].first;
                        schedule->order[1] = restricts_vector[0][i].second;
                    }
                }
            }

        } else {
            for (int i = 0; i < restricts_vector.size(); i++) { // 有多组约束的情况
                for (int j = 0; j < restricts_vector[i].size(); j++) {
                    if (degree[restricts_vector[i][j].first] + degree[restricts_vector[i][j].second] > max_weight) {
                        if (adj_mat[INDEX(restricts_vector[i][j].first, restricts_vector[i][j].second, size)]) {
                            max_weight = degree[restricts_vector[i][j].first] + degree[restricts_vector[i][j].second];
                            schedule->order[0] = restricts_vector[i][j].first;
                            schedule->order[1] = restricts_vector[i][j].second;
                        }
                    }
                }
            }
        }
    }

    // 没找到则找最重的边
    max_weight = 0;
    if (schedule->order[0] == schedule->order[1]) {
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (adj_mat[INDEX(i, j, size)]) {
                    if (vtx_ptr[i + 1] - vtx_ptr[i] + vtx_ptr[j + 1] - vtx_ptr[j] > max_weight) {
                        max_weight = vtx_ptr[i + 1] - vtx_ptr[i] + vtx_ptr[j + 1] - vtx_ptr[j];
                        schedule->order[0] = i;
                        schedule->order[1] = j;
                    }
                }
            }
        }
    }

    if (vtx_ptr[schedule->order[0] + 1] - vtx_ptr[schedule->order[0]] < vtx_ptr[schedule->order[1] + 1] - vtx_ptr[schedule->order[1]]) {
        int tempvtx = schedule->order[0];
        schedule->order[0] = schedule->order[1];
        schedule->order[1] = tempvtx;
    }
    already_assigned[schedule->order[0]] = 1;
    already_assigned[schedule->order[1]] = 1;
    schedule->order_degree[0] = vtx_ptr[schedule->order[0] + 1] - vtx_ptr[schedule->order[0]];
    schedule->order_degree[1] = vtx_ptr[schedule->order[1] + 1] - vtx_ptr[schedule->order[1]];

    int already_assigned_cnt = 2;

    //  相邻的顶点并且度数最大       先标记相邻，然后标记与几个顶点相邻， 然后选择度数  交集优先操作
    //  如果发生与一个相邻的度数比较大？与两个相邻的度数反而比较小？

    int *touched_times = new int[size];
    int max_touched_times;

    int *candidate_vtx = new int[size];
    int candidate_vtx_cnt = 0;

    int *source_label = new int[size];

    int max_degree;

    /*确定顶点处理顺序 以及记录 可以进行数值计算的层次 */
    while (already_assigned_cnt != size) {

        // 根据已有分配决定后续分配
        memset(touched_times, 0, sizeof(int) * size);
        for (int i = 0; i < already_assigned_cnt; i++) {
            // 根据已有分配构建候选集
            for (int ptr = vtx_ptr[schedule->order[i]]; ptr < vtx_ptr[schedule->order[i] + 1]; ptr++) {
                int touched_vtx = vtx_end[ptr];
                if (already_assigned[touched_vtx] == 0)
                    touched_times[touched_vtx]++;
            }
        }

        // 寻找被访问次数最多的顶点集合
        max_touched_times = 0;
        for (int i = 0; i < size; i++) {
            if (touched_times[i] > max_touched_times) {
                max_touched_times = touched_times[i];
                candidate_vtx_cnt = 0;
                candidate_vtx[candidate_vtx_cnt++] = i;
            } else if (touched_times[i] == max_touched_times && touched_times[i] != 0) {
                candidate_vtx[candidate_vtx_cnt++] = i;
            }
        }

        // 判断是否可以使用数值计算代替
        if (candidate_vtx_cnt == size - already_assigned_cnt && size - already_assigned_cnt != 1) {
            // 判断是否共源,只要有一个不共源就不行
            memset(source_label, 0, sizeof(int) * size);
            for (int ptr = vtx_ptr[candidate_vtx[0]]; ptr < vtx_ptr[candidate_vtx[0] + 1]; ptr++) {
                source_label[vtx_end[ptr]] = 1;
            }
            bool flag = true;
            for (int i = 1; i < candidate_vtx_cnt; i++) {
                for (int ptr = vtx_ptr[candidate_vtx[i]]; ptr < vtx_ptr[candidate_vtx[i] + 1]; ptr++) {
                    if (source_label[vtx_end[ptr]] == 0) {
                        flag = false;
                    }
                }
                if (flag == false) {
                    break;
                }
            }

            // 如果共源，记录这些顶点，并且记录可以进行数值计算的位置
            if (flag) {
                schedule->compute_depth = already_assigned_cnt;
                for (int i = 0; i < candidate_vtx_cnt; i++) {
                    schedule->order[already_assigned_cnt] = candidate_vtx[i];
                    schedule->order_degree[already_assigned_cnt] = vtx_ptr[candidate_vtx[i] + 1] - vtx_ptr[candidate_vtx[i]];
                    already_assigned_cnt++;
                }
            }
        }
        if (already_assigned_cnt == size) {
            break;
        }

        // 从候选集中选择度数最大的进行分配
        max_degree = 0;
        for (int i = 0; i < candidate_vtx_cnt; i++) {
            if (vtx_ptr[candidate_vtx[i] + 1] - vtx_ptr[candidate_vtx[i]] > max_degree) {
                max_degree = vtx_ptr[candidate_vtx[i] + 1] - vtx_ptr[candidate_vtx[i]];
                schedule->order[already_assigned_cnt] = candidate_vtx[i];
            }
        }
        schedule->order_degree[already_assigned_cnt] =
            vtx_ptr[schedule->order[already_assigned_cnt] + 1] - vtx_ptr[schedule->order[already_assigned_cnt]];
        already_assigned[schedule->order[already_assigned_cnt]] = 1;
        already_assigned_cnt++;
    }

    /* 确定好顶点处理顺序后，构建依赖关系 */
    schedule->dependcy_ptr[0] = 0;
    for (int process_vtx = 2; process_vtx < size; process_vtx++) {
        for (int i = 0; i < process_vtx; i++) {
            if (adj_mat[INDEX(process_vtx, i, size)]) {
                schedule->dependcy[schedule->dependcy_ptr[process_vtx - 2]++] = i;
            }
        }
        schedule->dependcy_ptr[process_vtx - 1] = schedule->dependcy_ptr[process_vtx - 2];
    }

    free(touched_times);
    free(candidate_vtx);
    free(source_label);
    free(already_assigned);
    free(adj_mat);
    free(vtx_ptr);
    free(vtx_end);
    free(degree);
}
// 测试共享顶点是否有用
void locality_verify(Graph *g, schedule_approx *s) {

    g->edge_from = (v_index_t *)malloc(sizeof(v_index_t) * g->e_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        for (int j = g->vertex[i]; j < g->vertex[i + 1]; j++) {
            g->edge_from[j] = i;
        }
    }

    uint64_t batch[5] = {10, 20, 40, 80, 160};

    uint32_t *vertex = new uint32_t[g->v_cnt];
    uint32_t *edge = new uint32_t[g->e_cnt];
    uint64_t vertex_size;
    uint64_t edge_size;
    std::random_device rd;
    std::default_random_engine gen(rd());
    std::uniform_int_distribution<uint64_t> dist(0, g->v_cnt - 1);
    double sample_count;
    for (int exp = 0; exp < 5; exp++) {
        vertex_size = 0;
        edge_size = 0;
        for (int i = 0; i < batch[exp]; i++) {
            uint64_t random_vtx = dist(gen);
            vertex[vertex_size++] = random_vtx;
            for (int ptr = g->edge[random_vtx]; ptr < g->edge[random_vtx + 1]; ptr++) {
                edge[edge_size++] = ptr;
            }
        }
        std::shuffle(edge, edge + edge_size, gen);

        // 开始实验测试

        // 边的测试
        TimeInterval Time;
        Time.check();
        for (int i = 0; i < edge_size; i++) {
            auto eid = edge[i];
            pattern_sample_es(eid, g, s, sample_count);
        }
        Time.print("edge");

        // 顶点的测试
        Time.check();
        for (int i = 0; i < vertex_size; i++) {
            auto vtx = vertex[i];
            for (uint64_t ptr = g->edge[vtx]; ptr < g->edge[vtx + 1]; ptr++)
                pattern_sample_es(ptr, g, s, sample_count);
        }
        Time.print("vertex");
    }
}
void scheduel_generate_approx_for_es(Pattern *P, schedule_approx *schedule) {
    // 先对矩阵进行排序
    int size = P->get_size();
    int *adj_mat = new int[size * size];
    memcpy(adj_mat, P->get_adj_mat_ptr(), size * size * sizeof(int));

    std::vector<std::vector<std::pair<int, int>>> restricts_vector;
    restricts_vector.clear();
    restricts_generate_approx(adj_mat, restricts_vector, size);

    for (int i = 0; i < restricts_vector.size(); ++i) {
        for (int j = 0; j < restricts_vector[i].size(); ++j) {
            printf("%d %d\n", restricts_vector[i][j].first, restricts_vector[i][j].second);
        }
    }
    // 根据限制类型生成不同的限制组

    int *vtx_ptr = new int[size + 2];
    int *vtx_end = new int[size * size];
    vtx_ptr[0] = 0;
    int iter_vtx = 1;
    vtx_ptr[iter_vtx] = vtx_ptr[0];
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (adj_mat[INDEX(i, j, size)]) {
                vtx_end[vtx_ptr[iter_vtx]++] = j;
            }
        }
        vtx_ptr[iter_vtx + 1] = vtx_ptr[iter_vtx];
        iter_vtx++;
    }

    // 初始化调度方案
    schedule->order = new uint32_t[size];
    schedule->order_degree = new uint32_t[size];
    schedule->dependcy = new uint32_t[size * size];
    schedule->dependcy_ptr = new uint32_t[size + 1];
    schedule->compute_depth = size;
    schedule->pattern_size = size;
    schedule->cache_insetion_recorder = new uint32_t[size];

    // 寻找最重的边作为起始边
    int *already_assigned = new int[size];
    memset(already_assigned, 0, sizeof(int) * size);
    int max_weight = 0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (adj_mat[INDEX(i, j, size)]) {
                if (vtx_ptr[i + 1] - vtx_ptr[i] + vtx_ptr[j + 1] - vtx_ptr[j] > max_weight) {
                    max_weight = vtx_ptr[i + 1] - vtx_ptr[i] + vtx_ptr[j + 1] - vtx_ptr[j];
                    schedule->order[0] = i;
                    schedule->order[1] = j;
                }
            }
        }
    }

    if (vtx_ptr[schedule->order[0] + 1] - vtx_ptr[schedule->order[0]] < vtx_ptr[schedule->order[1] + 1] - vtx_ptr[schedule->order[1]]) {
        int tempvtx = schedule->order[0];
        schedule->order[0] = schedule->order[1];
        schedule->order[1] = tempvtx;
    }
    already_assigned[schedule->order[0]] = 1;
    already_assigned[schedule->order[1]] = 1;
    schedule->order_degree[0] = vtx_ptr[schedule->order[0] + 1] - vtx_ptr[schedule->order[0]];
    schedule->order_degree[1] = vtx_ptr[schedule->order[1] + 1] - vtx_ptr[schedule->order[1]];

    int already_assigned_cnt = 2;

    //  相邻的顶点并且度数最大       先标记相邻，然后标记与几个顶点相邻， 然后选择度数  交集优先操作
    //  如果发生与一个相邻的度数比较大？与两个相邻的度数反而比较小？

    int *touched_times = new int[size];
    int max_touched_times;

    int *candidate_vtx = new int[size];
    int candidate_vtx_cnt = 0;

    int *source_label = new int[size];

    int max_degree;

    /*确定顶点处理顺序 以及记录 可以进行数值计算的层次 */
    while (already_assigned_cnt != size) {

        // 根据已有分配决定后续分配
        memset(touched_times, 0, sizeof(int) * size);
        for (int i = 0; i < already_assigned_cnt; i++) {
            // 根据已有分配构建候选集
            for (int ptr = vtx_ptr[schedule->order[i]]; ptr < vtx_ptr[schedule->order[i] + 1]; ptr++) {
                int touched_vtx = vtx_end[ptr];
                if (already_assigned[touched_vtx] == 0)
                    touched_times[touched_vtx]++;
            }
        }

        // 寻找被访问次数最多的顶点集合
        max_touched_times = 0;
        for (int i = 0; i < size; i++) {
            if (touched_times[i] > max_touched_times) {
                max_touched_times = touched_times[i];
                candidate_vtx_cnt = 0;
                candidate_vtx[candidate_vtx_cnt++] = i;
            } else if (touched_times[i] == max_touched_times) {
                candidate_vtx[candidate_vtx_cnt++] = i;
            }
        }

        // 判断是否可以使用数值计算代替
        if (candidate_vtx_cnt == size - already_assigned_cnt && size - already_assigned_cnt != 1) {
            // 判断是否共源,只要有一个不共源就不行
            memset(source_label, 0, sizeof(int) * size);
            for (int ptr = vtx_ptr[candidate_vtx[0]]; ptr < vtx_ptr[candidate_vtx[0] + 1]; ptr++) {
                source_label[vtx_end[ptr]] = 1;
            }
            bool flag = true;
            for (int i = 1; i < candidate_vtx_cnt; i++) {
                for (int ptr = vtx_ptr[candidate_vtx[i]]; ptr < vtx_ptr[candidate_vtx[i] + 1]; ptr++) {
                    if (source_label[vtx_end[ptr]] == 0) {
                        flag = false;
                    }
                }
                if (flag == false) {
                    break;
                }
            }

            // 如果共源，记录这些顶点，并且记录可以进行数值计算的位置
            if (flag) {
                schedule->compute_depth = already_assigned_cnt;
                for (int i = 0; i < candidate_vtx_cnt; i++) {
                    schedule->order[already_assigned_cnt] = candidate_vtx[i];
                    schedule->order_degree[already_assigned_cnt] = vtx_ptr[candidate_vtx[i] + 1] - vtx_ptr[candidate_vtx[i]];
                    already_assigned_cnt++;
                }
            }
        }
        if (already_assigned_cnt == size) {
            break;
        }

        // 从候选集中选择度数最大的进行分配
        max_degree = 0;
        for (int i = 0; i < candidate_vtx_cnt; i++) {
            if (vtx_ptr[candidate_vtx[i] + 1] - vtx_ptr[candidate_vtx[i]] > max_degree) {
                max_degree = vtx_ptr[candidate_vtx[i] + 1] - vtx_ptr[candidate_vtx[i]];
                schedule->order[already_assigned_cnt] = candidate_vtx[i];
            }
        }
        schedule->order_degree[already_assigned_cnt] =
            vtx_ptr[schedule->order[already_assigned_cnt] + 1] - vtx_ptr[schedule->order[already_assigned_cnt]];
        already_assigned[schedule->order[already_assigned_cnt]] = 1;
        already_assigned_cnt++;
    }

    /* 确定好顶点处理顺序后，构建实际处理顶点顺序的依赖关系 */
    schedule->dependcy_ptr[0] = 0;
    schedule->dependcy_ptr[1] = 0;
    schedule->dependcy_ptr[2] = 0;
    for (int process_id = 2; process_id < size; process_id++) {
        schedule->dependcy_ptr[process_id + 1] = schedule->dependcy_ptr[process_id];
        int process_vtx = schedule->order[process_id];
        for (int i = 0; i < process_id; i++) {
            int pre_vtx = schedule->order[i];
            if (adj_mat[INDEX(process_vtx, pre_vtx, size)]) {
                schedule->dependcy[schedule->dependcy_ptr[process_id + 1]++] = i;
            }
        }
    }

    free(touched_times);
    free(candidate_vtx);
    free(source_label);
    free(already_assigned);
    free(adj_mat);
    free(vtx_ptr);
    free(vtx_end);
}

Pattern pattern_load(char *filename) {
    std::fstream file(filename);
    if (!file.is_open()) {
        printf("open file %s failed\n", filename);
    }

    std::string line;
    std::getline(file, line);

    int vtx_num, edge_num;
    std::istringstream iss(line);
    iss >> vtx_num >> edge_num;

    char *buffer = new char[vtx_num * vtx_num];
    for (int i = 0; i < vtx_num * vtx_num; i++) {
        buffer[i] = '0';
    }

    int vtxleft, vtxright;
    for (int i = 0; i < edge_num; i++) {
        std::getline(file, line);
        std::istringstream iss(line);
        iss >> vtxleft >> vtxright;
        buffer[INDEX(vtxleft, vtxright, vtx_num)] = '1';
    }

    Pattern pattern(vtx_num, buffer);
    return pattern;
}

int pattern_count_ns(Graph *g, schedule_approx *s) {

    // 接收到的两个顶点
    uint32_t left_vtx;
    uint32_t right_vtx;

    // 已经匹配的顶点集合
    VertexSet matched_vetex_set;
    uint32_t *matched_vetex = new uint32_t[s->pattern_size];
    matched_vetex_set.init();
    matched_vetex_set.push_back(left_vtx);
    matched_vetex_set.push_back(right_vtx);
    matched_vetex[0] = left_vtx;
    matched_vetex[1] = right_vtx;

    // 候选集初始化
    VertexSet candidate_vtx[s->pattern_size];
    for (int i = 0; i < s->pattern_size; i++) {
        candidate_vtx[i].init();
    }
    candidate_vtx[0].push_back(left_vtx);
    candidate_vtx[1].push_back(right_vtx);

    // 临时顶点集合用于计算，scalegpm不构建候选集，因此需要使用两个反复计算
    VertexSet tmp_vtx_set;
    uint32_t tmp_vtx_set_size;

    // 随机采样一条边之后继续进行处理
    for (int process_id = 2; process_id < s->pattern_size; process_id++) {

        // 构架候选集 直接计算交集或者邻居集合
        // dependency记录的应该已经匹配id的顺序

        // 获取依赖的第一个顶点集合
        int tmp_pre_vtx = matched_vetex[s->dependcy[s->dependcy_ptr[process_id]]];
        tmp_vtx_set_size = g->vertex[tmp_pre_vtx + 1] - g->vertex[tmp_pre_vtx];
        candidate_vtx[process_id].init(tmp_vtx_set_size, g->edge + g->vertex[tmp_pre_vtx]);

        // 与后面的做交集
        for (int j = s->dependcy_ptr[process_id] + 1; j < s->dependcy_ptr[process_id + 1]; j++) {
            tmp_pre_vtx = matched_vetex[s->dependcy[j]];
            tmp_vtx_set_size = g->vertex[tmp_pre_vtx + 1] - g->vertex[tmp_pre_vtx];
            tmp_vtx_set.init(tmp_vtx_set_size, g->edge + g->vertex[tmp_pre_vtx]);
            candidate_vtx[process_id].intersection_with(tmp_vtx_set);
        }

        // graphset没有显示的做差集，而是在需要的时候判断已匹配顶点集合中有没有当前选择顶点
        // 对于显示的做差集还是有需要
        // TODO ::实现显示的差集操作
        int scale = candidate_vtx[process_id].unordered_subtraction_size(candidate_vtx[process_id], matched_vetex_set);

        // 如果可以用数值计算代替集合计算
        if (s->compute_depth = process_id) {
            int size = candidate_vtx[process_id].get_size();
            int scale = 1;
            for (int i = 0; i < s->pattern_size - s->compute_depth; i++) {
                scale = scale * size / (i + 1);
                size--;
            }
            return scale;
        }

        // m个顶点中有n个可行顶点，等概率从m中选择顶点，是否这n个顶点有相同的概率被选择？
        // TODO 随机选择顶点
        int random_vtx;
        matched_vetex_set.push_back(random_vtx);
    }
}

int main(int argc, char *argv[]) {
    Graph *g;
    DataLoader D;

    if (argc != 3) {
        printf("usage: %s graph_file pattern_file\n", argv[0]);
        return 0;
    }
    DataType type = DataType::Patents;
    // bool ok = D.fast_load(g, argv[1]);
    // bool ok = D.load_data(g, type, argv[1], 0, 0);
    // if (!ok) {
    //     printf("Load data failed\n");
    //     return 0;
    // }

    printf("Load data success!\n");
    fflush(stdout);
    printf("thread num: %d\n", omp_get_max_threads());

    TimeInterval allTime;
    allTime.check();

    // 模式生成测试
    Pattern pattern = pattern_load(argv[2]);
    schedule_approx schedule;
    scheduel_generate_approx_for_ns(&pattern, &schedule);
    for (int i = 0; i < pattern.get_size(); i++) {
        printf("order: %d, degree: %d\n", schedule.order[i], schedule.order_degree[i]);
    }

    // demo_pattern_test_exact(g);
    // demo_pattern_test_sample(g);
    // demo_pattern_test_filter_sample(g);
    // demo_pattern_test_filter_reduce_one(g);
    // demo_pattern_test_filter_reduce_two(g);
    // demo_pattern_test_filter_reduce_hybird(g);
    allTime.print("Total time cost");

    // g->triangle_statistics();
    // g->four_clique_statistics();

    // delete g;
    return 0;
}