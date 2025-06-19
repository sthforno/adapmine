#include "../include/approx.h"
#include "../include/common.h"
#include "../include/graph.h"
#include "../include/motif_generator.h"
#include "../include/schedule_IEP.h"
#include "../include/vertex_set.h"
#include "timeinterval.h"
#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <mpi.h>
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
    uint32_t *order;        // 匹配顺序
    uint32_t *order_degree; // 匹配顺序中顶点的度
    int *dependcy;          // 依赖的前面的顶点                //复用计算结果，将依赖顶点转换为计算结果
    int *dependcy_ptr;      // 依赖的指针
    int *dependcy_label;    // 记录是否有重复利用的计算结果
    uint32_t compute_depth; // 记录可以数值计算的深度
    uint32_t pattern_size;  // 记录模式的大小

    bool symmetry_label; // 为ns构建约束打破对称性

    // 差集优化？
    int *sub_id;
    int *sub_ptr;

    // 为es构建约束打破对称性
    // 待实现功能   1，为es施加约束打破对称性，同时满足随机选择要求  2，b,复用部分计算结果 c，将计算提到循环之前（只有es需要）。
};

struct filter {

    uint32_t *record;
    uint32_t *degree;
    uint32_t edge_size;

    // 边上的采样估计和采样次数，以及是否被移除
    uint32_t *sample_times_per_edge;
    double *sample_count_per_edge;
    double *sample_square_per_edge;
    bool *delete_edge;

    // 区域内的采样估计和采样次数
    int num_of_thread;
    uint32_t **sample_times_per_region;
    double **sample_count_per_region;
    double **sample_square_per_region;

    // 记录顶点和边的颜色
    int *vertex_record;
    int *edge_record;

    // 原始的顶点和边
    int64_t *edge_project;
    int *vertex_project;

    // 不同颜色的边数
    uint32_t num_of_color;
    uint64_t *num_of_different_weight_of_edges;
    uint64_t *prefix_sum_of_edges;
};

// 初始化过滤的数据结构

int find_range(uint32_t threshhold[], const uint32_t &degree, int size) {
    for (int i = 0; i < size; i++) {
        if (degree <= threshhold[i]) {
            return i - 1;
        }
    }
}

void filter_init(filter *filter_data, Graph *g, schedule_approx *s, uint32_t threshhold[]) {

    // 记录哪些边可以采样
    filter_data->edge_project = (int64_t *)malloc(sizeof(int64_t) * g->e_cnt);
    filter_data->edge_record = (int *)malloc(sizeof(int) * g->e_cnt);
    memset(filter_data->edge_project, -1, sizeof(int64_t) * g->e_cnt);
    memset(filter_data->edge_record, -1, sizeof(int) * g->e_cnt);
    filter_data->degree = (uint32_t *)malloc(sizeof(uint32_t) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        filter_data->degree[i] = g->vertex[i + 1] - g->vertex[i];
    }

    // 标记边
    // 记录每种范围边数据的大小
    filter_data->num_of_different_weight_of_edges = (uint64_t *)malloc(sizeof(uint64_t) * filter_data->num_of_color);
    for (uint64_t i = 0; i < filter_data->num_of_color; ++i) {
        filter_data->num_of_different_weight_of_edges[i] = 0;
    }
    for (int i = 0; i < g->v_cnt; ++i) {
        if (filter_data->degree[i] >= s->order_degree[0]) {
            for (int64_t j = g->vertex[i]; j < g->vertex[i + 1]; ++j) {
                uint32_t temp_degree = filter_data->degree[g->edge[j]];
                if (temp_degree >= s->order_degree[1]) {
                    int record_edge_color = find_range(threshhold, temp_degree, filter_data->num_of_color + 1);
                    filter_data->edge_record[j] = record_edge_color;
                    filter_data->num_of_different_weight_of_edges[record_edge_color]++;
                }
            }
        }
    }

    // 投射边
    // 前缀和形式
    filter_data->prefix_sum_of_edges = new uint64_t[filter_data->num_of_color + 1];
    filter_data->prefix_sum_of_edges[0] = 0;
    for (uint32_t i = 1; i < filter_data->num_of_color + 1; i++) {
        filter_data->prefix_sum_of_edges[i] = filter_data->prefix_sum_of_edges[i - 1] + filter_data->num_of_different_weight_of_edges[i - 1];
    }
    uint64_t weight_edge_ptr[filter_data->num_of_color] = {0}; // 只在此处使用，将边按照其出端点的度数分类投射。
    for (uint i = 0; i < filter_data->num_of_color; i++) {
        weight_edge_ptr[i] = filter_data->prefix_sum_of_edges[i];
    }
    for (uint32_t i = 0; i < g->e_cnt; ++i) {
        if (filter_data->edge_record[i] != -1)
            filter_data->edge_project[weight_edge_ptr[filter_data->edge_record[i]]++] = i;
    }

    // 记录采样结果
    filter_data->sample_times_per_edge = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
    filter_data->sample_count_per_edge = (double *)malloc(sizeof(double) * g->e_cnt);
    filter_data->sample_square_per_edge = (double *)malloc(sizeof(double) * g->e_cnt);
    memset(filter_data->sample_times_per_edge, 0, sizeof(uint32_t) * g->e_cnt);
    memset(filter_data->sample_count_per_edge, 0, sizeof(double) * g->e_cnt);
    memset(filter_data->sample_square_per_edge, 0, sizeof(double) * g->e_cnt);

    filter_data->num_of_thread = omp_get_max_threads();

    filter_data->sample_times_per_region = (uint32_t **)malloc(sizeof(uint32_t *) * (filter_data->num_of_thread + 1));
    filter_data->sample_count_per_region = (double **)malloc(sizeof(double *) * (filter_data->num_of_thread + 1));
    filter_data->sample_square_per_region = (double **)malloc(sizeof(double *) * (filter_data->num_of_thread + 1));

    for (int i = 0; i <= filter_data->num_of_thread; i++) {
        filter_data->sample_times_per_region[i] = (uint32_t *)malloc(sizeof(uint32_t) * filter_data->num_of_color);
        filter_data->sample_count_per_region[i] = (double *)malloc(sizeof(double) * filter_data->num_of_color);
        filter_data->sample_square_per_region[i] = (double *)malloc(sizeof(double) * filter_data->num_of_color);
        memset(filter_data->sample_times_per_region[i], 0, sizeof(uint32_t) * filter_data->num_of_color);
        memset(filter_data->sample_count_per_region[i], 0, sizeof(double) * filter_data->num_of_color);
        memset(filter_data->sample_square_per_region[i], 0, sizeof(double) * filter_data->num_of_color);
    }
}

void free_filter_data(filter &filter_data) {
    free(filter_data.edge_project);
    free(filter_data.edge_record);
    free(filter_data.degree);
    free(filter_data.sample_count_per_edge);
    free(filter_data.sample_times_per_edge);
    free(filter_data.sample_square_per_edge);

    for (int i = 0; i <= filter_data.num_of_thread; i++) {
        free(filter_data.sample_count_per_region[i]);
        free(filter_data.sample_square_per_region[i]);
        free(filter_data.sample_times_per_region[i]);
    }
    free(filter_data.sample_count_per_region);
    free(filter_data.sample_square_per_region);
    free(filter_data.sample_times_per_region);
}

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

int intersection_size(const VertexSet &set0, const VertexSet &set1) {
    int i = 0;
    int j = 0;
    int size0 = set0.get_size();
    int size1 = set1.get_size();

    int data0 = set0.get_data(0);
    int data1 = set1.get_data(0);
    int size = 0;

    // 遍历两个数组
    while (i < size0 && j < size1) {
        data0 = set0.get_data(i);
        data1 = set1.get_data(j);
        if (data0 < data1)
            ++i;
        else if (data0 > data1)
            ++j;
        else {
            ++i;
            ++j;
            size++;
        }
    }

    return size;
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
    schedule->dependcy = new int[size * size];
    schedule->dependcy_ptr = new int[size + 1];
    schedule->dependcy_label = new int[size];
    memset(schedule->dependcy_label, -1, sizeof(int) * size);
    schedule->compute_depth = size;
    schedule->pattern_size = size;
    if (restricts_vector.size() == 0) {
        schedule->symmetry_label = false;
    } else {
        schedule->symmetry_label = true;
    }

    schedule->sub_id = new int[size * size];
    schedule->sub_ptr = new int[size + 1];

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

    /* 确定好顶点处理顺序后，构建依赖关系和差集信息 */
    schedule->dependcy_ptr[0] = 0;
    schedule->dependcy_ptr[1] = 0;
    schedule->dependcy_ptr[2] = 0;

    schedule->sub_ptr[0] = 0;
    schedule->sub_ptr[1] = 0;
    schedule->sub_ptr[2] = 0;
    for (int process_vtx = 2; process_vtx < size; process_vtx++) {
        schedule->dependcy_ptr[process_vtx + 1] = schedule->dependcy_ptr[process_vtx];
        schedule->sub_ptr[process_vtx + 1] = schedule->sub_ptr[process_vtx];
        for (int i = 0; i < process_vtx; i++) {
            if (adj_mat[INDEX(process_vtx, i, size)]) {
                schedule->dependcy[schedule->dependcy_ptr[process_vtx + 1]++] = i;
            } else {
                schedule->sub_id[schedule->sub_ptr[process_vtx + 1]++] = i;
            }
        }
    }

    /*构建好依赖关系之后，转化为记录已有计算结果，此处实现的是完全复用 */
    VertexSet tmpset[2];
    int similar_vtx, max_similar_size;
    for (int process_vtx = size - 1; process_vtx >= 2; process_vtx--) {

        tmpset[0].init(schedule->dependcy_ptr[process_vtx + 1] - schedule->dependcy_ptr[process_vtx],
                       schedule->dependcy + schedule->dependcy_ptr[process_vtx]);
        max_similar_size = 0;
        for (int pre_vtx = process_vtx - 1; pre_vtx >= 2; pre_vtx--) {
            tmpset[1].init(schedule->dependcy_ptr[pre_vtx + 1] - schedule->dependcy_ptr[pre_vtx],
                           schedule->dependcy + schedule->dependcy_ptr[pre_vtx]);
            int temp_similar_size = intersection_size(tmpset[0], tmpset[1]);
            if (temp_similar_size > max_similar_size && tmpset[1].get_size() == temp_similar_size) {
                max_similar_size = temp_similar_size;
                similar_vtx = pre_vtx;
            }
        }

        // 如果有可以复用的结果，对其进行复用
        if (max_similar_size >= 2) {
            tmpset[1].init(schedule->dependcy_ptr[similar_vtx + 1] - schedule->dependcy_ptr[similar_vtx],
                           schedule->dependcy + schedule->dependcy_ptr[similar_vtx]);

            int record_ptr = schedule->dependcy_ptr[process_vtx + 1];
            for (int ptr = tmpset[0].get_size() - 1; ptr >= 0; ptr--) {
                int temp_vtx = tmpset[0].get_data(ptr);
                if (!tmpset[1].has_data(temp_vtx)) {
                    schedule->dependcy[--record_ptr] = temp_vtx;
                }
            }
            while (record_ptr != schedule->dependcy_ptr[process_vtx]) {
                schedule->dependcy[--record_ptr] = -1;
            }
            schedule->dependcy[schedule->dependcy_ptr[process_vtx]] = similar_vtx;
            schedule->dependcy_label[process_vtx] = 0;
        }
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

template <typename T>
int binary_search(const T data[], int n, const T &target) {
    int mid, l = 0, r = n - 1;
    while (l <= r) {
        mid = (l + r) >> 1;
        if (data[mid] < target) {
            l = mid + 1;
        } else if (data[mid] > target) {
            r = mid - 1;
        } else {
            return mid;
        }
    }
    return -1;
}
void subtraction(Graph *g, VertexSet &set1, VertexSet &set2) {
    int *data = set1.get_data_ptr();
    int *label = new int[set1.get_size()];
    memset(label, 0, sizeof(int) * set1.get_size());

    for (int i = 0; i < set2.get_size(); i++) {
        int temp_vtx = set2.get_data(i);
        int loc = binary_search(data, set1.get_size(), temp_vtx);
        if (loc != -1) {
            label[loc] = -1;
        }
    }

    int size = 0;
    for (int i = 0; i < set1.get_size(); i++) {
        if (label[i] != -1) {
            data[size++] = data[i];
        }
    }

    int delet_vtx = set1.get_size() - size;
    for (int i = 0; i < delet_vtx; i++) {
        set1.pop_back();
    }

    free(label);
}

void subtraction_optimized(Graph *g, VertexSet &set1, int process_id, schedule_approx *s, uint32_t *matched_vtx) {
    int *data = set1.get_data_ptr();
    int *label = new int[set1.get_size()];
    memset(label, 0, sizeof(int) * set1.get_size());

    for (int ptr = s->sub_ptr[process_id]; ptr < s->sub_ptr[process_id + 1]; ptr++) {
        int pre_id = s->sub_id[ptr];
        int temp_vtx = matched_vtx[pre_id];
        int loc = binary_search(data, set1.get_size(), temp_vtx);
        if (loc != -1) {
            label[loc] = -1;
        }
    }

    int size = 0;
    for (int i = 0; i < set1.get_size(); i++) {
        if (label[i] != -1) {
            data[size++] = data[i];
        }
    }

    int delet_vtx = set1.get_size() - size;
    for (int i = 0; i < delet_vtx; i++) {
        set1.pop_back();
    }

    free(label);
}
// 传过来的edge是已经满足剪枝要求和对称性要求的了
void pattern_sample_ns(uint32_t edge, Graph *g, schedule_approx *s, double &sample_count, double &count_sample_square) {

    std::random_device rd;
    std::default_random_engine gen(rd());

    // 接收到的两个顶点
    uint32_t left_vtx = g->edge_from[edge];
    uint32_t right_vtx = g->edge[edge];

    // 已经匹配的顶点集合
    VertexSet matched_vetex_set;
    uint32_t *matched_vetex = new uint32_t[s->pattern_size];
    matched_vetex_set.init();
    matched_vetex_set.push_back(left_vtx);
    matched_vetex_set.push_back(right_vtx);
    matched_vetex[0] = left_vtx;
    matched_vetex[1] = right_vtx;

    double temp_count = 1;
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
        if (s->dependcy_label[process_id] == -1) {
            int tmp_pre_vtx = matched_vetex[s->dependcy[s->dependcy_ptr[process_id]]];
            int tmp_vtx_set_size = g->vertex[tmp_pre_vtx + 1] - g->vertex[tmp_pre_vtx];
            candidate_vtx[process_id].copy(tmp_vtx_set_size, g->edge + g->vertex[tmp_pre_vtx]);

            // 与后面的做交集
            VertexSet tmp_vtx_set;
            for (int j = s->dependcy_ptr[process_id] + 1; j < s->dependcy_ptr[process_id + 1]; j++) {
                tmp_pre_vtx = matched_vetex[s->dependcy[j]];
                tmp_vtx_set_size = g->vertex[tmp_pre_vtx + 1] - g->vertex[tmp_pre_vtx];
                tmp_vtx_set.init(tmp_vtx_set_size, g->edge + g->vertex[tmp_pre_vtx]);
                candidate_vtx[process_id].intersection_with(tmp_vtx_set);
            }
        } else {
            // 直接复用结果
            int reuse_id = s->dependcy[s->dependcy_ptr[process_id]];
            candidate_vtx[process_id].copy(candidate_vtx[reuse_id].get_size(), candidate_vtx[reuse_id].get_data_ptr());

            // 与后面的做交集
            VertexSet tmp_vtx_set;
            int tmp_pre_vtx;
            int tmp_vtx_set_size;
            for (int j = s->dependcy_ptr[process_id] + 1; j < s->dependcy_ptr[process_id + 1]; j++) {
                if (s->dependcy[j] != -1) {
                    tmp_pre_vtx = matched_vetex[s->dependcy[j]];
                    tmp_vtx_set_size = g->vertex[tmp_pre_vtx + 1] - g->vertex[tmp_pre_vtx];
                    tmp_vtx_set.init(tmp_vtx_set_size, g->edge + g->vertex[tmp_pre_vtx]);
                    candidate_vtx[process_id].intersection_with(tmp_vtx_set);
                }
            }
        }

        // graphset没有显示的做差集，而是在需要的时候判断已匹配顶点集合中有没有当前选择顶点。
        // 有一个类似的实现，anti顶点移除，对于candidat中的顶点到差集中去做二分查找。D*logd次搜索，D次移动。后者到前者中查找呢，d*logD次查找，D*d次移动。
        // 考虑后者，从被减去的顶点集合中的顶点往候选集中查找，然后标记位置，最后一轮移除，这样就是d*logD次查找,d次移动。
        // TODO ::实现显示的差集操作  我们的实现产生的大小更小？

        // 两种差集操作符合预期，关键在于候选集中会有-1顶点和重复顶点  关键在预测初始化时，就有-1边和重复边, 关键在于边表就有-1和重复边
        // VertexSet tmp_set;
        // tmp_set.copy(candidate_vtx[process_id].get_size(), candidate_vtx[process_id].get_data_ptr());
        // int tmp = tmp_set.unordered_subtraction_size(candidate_vtx[process_id], matched_vetex_set);

        // subtraction(g, candidate_vtx[process_id], matched_vetex_set);

        if (candidate_vtx[process_id].get_size() == 0) {
            return;
        }
        subtraction_optimized(g, candidate_vtx[process_id], process_id, s, matched_vetex);
        if (candidate_vtx[process_id].get_size() == 0) {
            return;
        }

        // if (candidate_vtx[process_id].get_size() != tmp) {
        //     printf("wrong");
        //     subtraction(g, tmp_set, matched_vetex_set);
        // }

        // 如果可以用数值计算代替集合计算
        if (s->compute_depth == process_id && s->compute_depth != s->pattern_size) {
            uint64_t size = candidate_vtx[process_id].get_size();
            uint32_t num = s->pattern_size - s->compute_depth;

            if (size < num) {
                return;
            }

            for (int i = 1; i <= num; i++) {
                temp_count = temp_count * size;
                size--;
            }

            // for (int i = 1; i <= num; i++) {
            //     temp_count = temp_count * size / i;
            //     size--;
            // }

            for (int i = 2; i < process_id; i++) {
                temp_count *= candidate_vtx[i].get_size();
            }

            sample_count += (temp_count * g->e_cnt);
            count_sample_square += (temp_count * g->e_cnt * temp_count * g->e_cnt);
            return;
        }

        // m个顶点中有n个可行顶点，等概率从m中选择顶点，是否这n个顶点有相同的概率被选择？
        // TODO 随机选择是否存在优化的地方？

        if (candidate_vtx[process_id].get_size() >= 1) {
            std::uniform_int_distribution<uint64_t> dist(0, candidate_vtx[process_id].get_size() - 1);
            int random_id = dist(gen);
            int random_vtx = candidate_vtx[process_id].get_data(random_id);
            matched_vetex_set.push_back(random_vtx);
            matched_vetex[process_id] = random_vtx;

            if (g->degree[random_vtx] < s->order_degree[process_id]) {
                return;
            }

            // for (int i = 0; i <= process_id; i++) {
            //     if (matched_vetex_set.has_data(matched_vetex[i])) {
            //         continue;
            //     } else {
            //         printf("wrong");
            //     }
            // }

        } else {
            return;
        }
    }

    for (int i = 0; i < s->pattern_size; i++) {
        temp_count = temp_count * candidate_vtx[i].get_size();
    }

    // std::ofstream of;
    // of.open("livej_triangle.txt", std::ios::app);
    // of << temp_count << std::endl;
    // of.close();
    sample_count += (temp_count * g->e_cnt);
    count_sample_square += (temp_count * temp_count * g->e_cnt * g->e_cnt);
}

void pattern_sample_ns_record(uint32_t edge, Graph *g, schedule_approx *s, filter *filter_data) {

    std::random_device rd;
    std::default_random_engine gen(rd());

    __sync_fetch_and_add(&filter_data->sample_times_per_edge[edge], 1);

    // 接收到的两个顶点
    uint32_t left_vtx = g->edge_from[edge];
    uint32_t right_vtx = g->edge[edge];

    // 已经匹配的顶点集合
    VertexSet matched_vetex_set;
    uint32_t *matched_vetex = new uint32_t[s->pattern_size];
    matched_vetex_set.init();
    matched_vetex_set.push_back(left_vtx);
    matched_vetex_set.push_back(right_vtx);
    matched_vetex[0] = left_vtx;
    matched_vetex[1] = right_vtx;

    double temp_count = 1;
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
        if (s->dependcy_label[process_id] == -1) {
            int tmp_pre_vtx = matched_vetex[s->dependcy[s->dependcy_ptr[process_id]]];
            int tmp_vtx_set_size = g->vertex[tmp_pre_vtx + 1] - g->vertex[tmp_pre_vtx];
            candidate_vtx[process_id].copy(tmp_vtx_set_size, g->edge + g->vertex[tmp_pre_vtx]);
            // 与后面的做交集
            VertexSet tmp_vtx_set;
            for (int j = s->dependcy_ptr[process_id] + 1; j < s->dependcy_ptr[process_id + 1]; j++) {
                tmp_pre_vtx = matched_vetex[s->dependcy[j]];
                tmp_vtx_set_size = g->vertex[tmp_pre_vtx + 1] - g->vertex[tmp_pre_vtx];
                tmp_vtx_set.init(tmp_vtx_set_size, g->edge + g->vertex[tmp_pre_vtx]);
                candidate_vtx[process_id].intersection_with(tmp_vtx_set);
            }
        } else {
            // 直接复用结果
            int reuse_id = s->dependcy[s->dependcy_ptr[process_id]];
            candidate_vtx[process_id].copy(candidate_vtx[reuse_id].get_size(), candidate_vtx[reuse_id].get_data_ptr());

            // 与后面的做交集
            VertexSet tmp_vtx_set;
            int tmp_pre_vtx;
            int tmp_vtx_set_size;
            for (int j = s->dependcy_ptr[process_id] + 1; j < s->dependcy_ptr[process_id + 1]; j++) {
                if (s->dependcy[j] != -1) {
                    tmp_pre_vtx = matched_vetex[s->dependcy[j]];
                    tmp_vtx_set_size = g->vertex[tmp_pre_vtx + 1] - g->vertex[tmp_pre_vtx];
                    tmp_vtx_set.init(tmp_vtx_set_size, g->edge + g->vertex[tmp_pre_vtx]);
                    candidate_vtx[process_id].intersection_with(tmp_vtx_set);
                }
            }
        }

        // graphset没有显示的做差集，而是在需要的时候判断已匹配顶点集合中有没有当前选择顶点。
        // 有一个类似的实现，anti顶点移除，对于candidat中的顶点到差集中去做二分查找。D*logd次搜索，D次移动。后者到前者中查找呢，d*logD次查找，D*d次移动。
        // 考虑后者，从被减去的顶点集合中的顶点往候选集中查找，然后标记位置，最后一轮移除，这样就是d*logD次查找,d次移动。
        // TODO ::实现显示的差集操作  我们的实现产生的大小更小？

        // 两种差集操作符合预期，关键在于候选集中会有-1顶点和重复顶点  关键在预测初始化时，就有-1边和重复边, 关键在于边表就有-1和重复边
        // VertexSet tmp_set;
        // tmp_set.copy(candidate_vtx[process_id].get_size(), candidate_vtx[process_id].get_data_ptr());
        // int tmp = tmp_set.unordered_subtraction_size(candidate_vtx[process_id], matched_vetex_set);

        // subtraction(g, candidate_vtx[process_id], matched_vetex_set);
        if (candidate_vtx[process_id].get_size() == 0) {
            return;
        }

        subtraction_optimized(g, candidate_vtx[process_id], process_id, s, matched_vetex);
        if (candidate_vtx[process_id].get_size() == 0) {
            return;
        }

        // if (candidate_vtx[process_id].get_size() != tmp) {
        //     printf("wrong");
        //     subtraction(g, tmp_set, matched_vetex_set);
        // }

        // 如果可以用数值计算代替集合计算
        if (s->compute_depth == process_id && s->compute_depth != s->pattern_size) {
            uint64_t size = candidate_vtx[process_id].get_size();
            uint32_t num = s->pattern_size - s->compute_depth;

            if (size < num) {
                return;
            }

            for (int i = 1; i <= num; i++) {
                temp_count = temp_count * size;
                size--;
            }

            for (int i = 2; i < process_id; i++) {
                temp_count *= candidate_vtx[i].get_size();
            }

            double temp_squra = temp_count * temp_count;

#pragma omp atomic
            filter_data->sample_count_per_edge[edge] += temp_count;

#pragma omp atomic
            filter_data->sample_square_per_edge[edge] += temp_squra;

            return;
        }

        // m个顶点中有n个可行顶点，等概率从m中选择顶点，是否这n个顶点有相同的概率被选择？
        // TODO 随机选择是否存在优化的地方？

        if (candidate_vtx[process_id].get_size() >= 1) {
            std::uniform_int_distribution<uint64_t> dist(0, candidate_vtx[process_id].get_size() - 1);
            int random_id = dist(gen);
            int random_vtx = candidate_vtx[process_id].get_data(random_id);
            matched_vetex_set.push_back(random_vtx);
            matched_vetex[process_id] = random_vtx;

            if (g->degree[random_vtx] < s->order_degree[process_id]) {
                return;
            }

        } else {
            return;
        }
    }

    for (int i = 0; i < s->pattern_size; i++) {
        temp_count = temp_count * candidate_vtx[i].get_size();
    }
    double temp_squra = temp_count * temp_count;
#pragma omp atomic
    filter_data->sample_count_per_edge[edge] += temp_count;

#pragma omp atomic
    filter_data->sample_square_per_edge[edge] += temp_squra;

    // __sync_fetch_and_add(&filter_data->sample_times_per_edge[edge], 1);
}

// 会移除无法产生匹配的边
void pattern_sample_ns_record_filter(uint32_t edge, Graph *g, schedule_approx *s, filter *filter_data) {

    std::random_device rd;
    std::default_random_engine gen(rd());

    __sync_fetch_and_add(&filter_data->sample_times_per_edge[edge], 1);

    // 接收到的两个顶点
    uint32_t left_vtx = g->edge_from[edge];
    uint32_t right_vtx = g->edge[edge];

    // 已经匹配的顶点集合
    VertexSet matched_vetex_set;
    uint32_t *matched_vetex = new uint32_t[s->pattern_size];
    matched_vetex_set.init();
    matched_vetex_set.push_back(left_vtx);
    matched_vetex_set.push_back(right_vtx);
    matched_vetex[0] = left_vtx;
    matched_vetex[1] = right_vtx;

    double temp_count = 1;
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
        if (s->dependcy_label[process_id] == -1) {
            int tmp_pre_vtx = matched_vetex[s->dependcy[s->dependcy_ptr[process_id]]];
            int tmp_vtx_set_size = g->vertex[tmp_pre_vtx + 1] - g->vertex[tmp_pre_vtx];
            candidate_vtx[process_id].copy(tmp_vtx_set_size, g->edge + g->vertex[tmp_pre_vtx]);
            // 与后面的做交集
            VertexSet tmp_vtx_set;
            for (int j = s->dependcy_ptr[process_id] + 1; j < s->dependcy_ptr[process_id + 1]; j++) {
                tmp_pre_vtx = matched_vetex[s->dependcy[j]];
                tmp_vtx_set_size = g->vertex[tmp_pre_vtx + 1] - g->vertex[tmp_pre_vtx];
                tmp_vtx_set.init(tmp_vtx_set_size, g->edge + g->vertex[tmp_pre_vtx]);
                candidate_vtx[process_id].intersection_with(tmp_vtx_set);
            }
        } else {
            // 直接复用结果
            int reuse_id = s->dependcy[s->dependcy_ptr[process_id]];
            candidate_vtx[process_id].copy(candidate_vtx[reuse_id].get_size(), candidate_vtx[reuse_id].get_data_ptr());

            // 与后面的做交集
            VertexSet tmp_vtx_set;
            int tmp_pre_vtx;
            int tmp_vtx_set_size;
            for (int j = s->dependcy_ptr[process_id] + 1; j < s->dependcy_ptr[process_id + 1]; j++) {
                if (s->dependcy[j] != -1) {
                    tmp_pre_vtx = matched_vetex[s->dependcy[j]];
                    tmp_vtx_set_size = g->vertex[tmp_pre_vtx + 1] - g->vertex[tmp_pre_vtx];
                    tmp_vtx_set.init(tmp_vtx_set_size, g->edge + g->vertex[tmp_pre_vtx]);
                    candidate_vtx[process_id].intersection_with(tmp_vtx_set);
                }
            }
        }

        // graphset没有显示的做差集，而是在需要的时候判断已匹配顶点集合中有没有当前选择顶点。
        // 有一个类似的实现，anti顶点移除，对于candidat中的顶点到差集中去做二分查找。D*logd次搜索，D次移动。后者到前者中查找呢，d*logD次查找，D*d次移动。
        // 考虑后者，从被减去的顶点集合中的顶点往候选集中查找，然后标记位置，最后一轮移除，这样就是d*logD次查找,d次移动。
        // TODO ::实现显示的差集操作  我们的实现产生的大小更小？

        // 两种差集操作符合预期，关键在于候选集中会有-1顶点和重复顶点  关键在预测初始化时，就有-1边和重复边, 关键在于边表就有-1和重复边
        // VertexSet tmp_set;
        // tmp_set.copy(candidate_vtx[process_id].get_size(), candidate_vtx[process_id].get_data_ptr());
        // int tmp = tmp_set.unordered_subtraction_size(candidate_vtx[process_id], matched_vetex_set);

        // subtraction(g, candidate_vtx[process_id], matched_vetex_set);

        if (candidate_vtx[process_id].get_size() == 0) {
            if (process_id == 2)
                filter_data->delete_edge[edge] = true;
            return;
        }
        subtraction_optimized(g, candidate_vtx[process_id], process_id, s, matched_vetex);

        if (candidate_vtx[process_id].get_size() == 0) {
            if (process_id == 2)
                filter_data->delete_edge[edge] = true;
            return;
        }

        // if (candidate_vtx[process_id].get_size() != tmp) {
        //     printf("wrong");
        //     subtraction(g, tmp_set, matched_vetex_set);
        // }

        // 如果可以用数值计算代替集合计算
        if (s->compute_depth == process_id && s->compute_depth != s->pattern_size) {
            uint64_t size = candidate_vtx[process_id].get_size();
            uint32_t num = s->pattern_size - s->compute_depth;

            if (size < num) {
                return;
            }

            for (int i = 1; i <= num; i++) {
                temp_count = temp_count * size;
                size--;
            }

            for (int i = 2; i < process_id; i++) {
                temp_count *= candidate_vtx[i].get_size();
            }

            double temp_squra = temp_count * temp_count;

#pragma omp atomic
            filter_data->sample_count_per_edge[edge] += temp_count;

#pragma omp atomic
            filter_data->sample_square_per_edge[edge] += temp_squra;

            return;
        }

        // m个顶点中有n个可行顶点，等概率从m中选择顶点，是否这n个顶点有相同的概率被选择？
        // TODO 随机选择是否存在优化的地方？

        if (candidate_vtx[process_id].get_size() >= 1) {
            std::uniform_int_distribution<uint64_t> dist(0, candidate_vtx[process_id].get_size() - 1);
            int random_id = dist(gen);
            int random_vtx = candidate_vtx[process_id].get_data(random_id);
            matched_vetex_set.push_back(random_vtx);
            matched_vetex[process_id] = random_vtx;

            if (g->degree[random_vtx] < s->order_degree[process_id]) {
                return;
            }

        } else {
            return;
        }
    }

    for (int i = 0; i < s->pattern_size; i++) {
        temp_count = temp_count * candidate_vtx[i].get_size();
    }
    double temp_squra = temp_count * temp_count;
#pragma omp atomic
    filter_data->sample_count_per_edge[edge] += temp_count;

#pragma omp atomic
    filter_data->sample_square_per_edge[edge] += temp_squra;

    // __sync_fetch_and_add(&filter_data->sample_times_per_edge[edge], 1);
}

// 同时更新域内的记录信息
void pattern_sample_ns_record_filter_region(uint32_t edge, Graph *g, schedule_approx *s, filter *filter_data, int thread_id, uint32_t threshhold[]) {

    std::random_device rd;
    std::default_random_engine gen(rd());

    // 寻找边的颜色领域
    int region_id = find_range(threshhold, filter_data->degree[g->edge[edge]], filter_data->num_of_color + 1);

    __sync_fetch_and_add(&filter_data->sample_times_per_edge[edge], 1); // 一定要在开始做，否则会产生错误
    __sync_fetch_and_add(&filter_data->sample_times_per_region[thread_id][region_id], 1);

    // 接收到的两个顶点
    uint32_t left_vtx = g->edge_from[edge];
    uint32_t right_vtx = g->edge[edge];

    // 已经匹配的顶点集合
    VertexSet matched_vetex_set;
    uint32_t *matched_vetex = new uint32_t[s->pattern_size];
    matched_vetex_set.init();
    matched_vetex_set.push_back(left_vtx);
    matched_vetex_set.push_back(right_vtx);
    matched_vetex[0] = left_vtx;
    matched_vetex[1] = right_vtx;

    double temp_count = 1;
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
        if (s->dependcy_label[process_id] == -1) {
            int tmp_pre_vtx = matched_vetex[s->dependcy[s->dependcy_ptr[process_id]]];
            int tmp_vtx_set_size = g->vertex[tmp_pre_vtx + 1] - g->vertex[tmp_pre_vtx];
            candidate_vtx[process_id].copy(tmp_vtx_set_size, g->edge + g->vertex[tmp_pre_vtx]);
            // 与后面的做交集
            VertexSet tmp_vtx_set;
            for (int j = s->dependcy_ptr[process_id] + 1; j < s->dependcy_ptr[process_id + 1]; j++) {
                tmp_pre_vtx = matched_vetex[s->dependcy[j]];
                tmp_vtx_set_size = g->vertex[tmp_pre_vtx + 1] - g->vertex[tmp_pre_vtx];
                tmp_vtx_set.init(tmp_vtx_set_size, g->edge + g->vertex[tmp_pre_vtx]);
                candidate_vtx[process_id].intersection_with(tmp_vtx_set);
            }
        } else {
            // 直接复用结果
            int reuse_id = s->dependcy[s->dependcy_ptr[process_id]];
            candidate_vtx[process_id].copy(candidate_vtx[reuse_id].get_size(), candidate_vtx[reuse_id].get_data_ptr());

            // 与后面的做交集
            VertexSet tmp_vtx_set;
            int tmp_pre_vtx;
            int tmp_vtx_set_size;
            for (int j = s->dependcy_ptr[process_id] + 1; j < s->dependcy_ptr[process_id + 1]; j++) {
                if (s->dependcy[j] != -1) {
                    tmp_pre_vtx = matched_vetex[s->dependcy[j]];
                    tmp_vtx_set_size = g->vertex[tmp_pre_vtx + 1] - g->vertex[tmp_pre_vtx];
                    tmp_vtx_set.init(tmp_vtx_set_size, g->edge + g->vertex[tmp_pre_vtx]);
                    candidate_vtx[process_id].intersection_with(tmp_vtx_set);
                }
            }
        }

        // graphset没有显示的做差集，而是在需要的时候判断已匹配顶点集合中有没有当前选择顶点。
        // 有一个类似的实现，anti顶点移除，对于candidat中的顶点到差集中去做二分查找。D*logd次搜索，D次移动。后者到前者中查找呢，d*logD次查找，D*d次移动。
        // 考虑后者，从被减去的顶点集合中的顶点往候选集中查找，然后标记位置，最后一轮移除，这样就是d*logD次查找,d次移动。
        // TODO ::实现显示的差集操作  我们的实现产生的大小更小？

        // 两种差集操作符合预期，关键在于候选集中会有-1顶点和重复顶点  关键在预测初始化时，就有-1边和重复边, 关键在于边表就有-1和重复边
        // VertexSet tmp_set;
        // tmp_set.copy(candidate_vtx[process_id].get_size(), candidate_vtx[process_id].get_data_ptr());
        // int tmp = tmp_set.unordered_subtraction_size(candidate_vtx[process_id], matched_vetex_set);

        // subtraction(g, candidate_vtx[process_id], matched_vetex_set);

        if (candidate_vtx[process_id].get_size() == 0) {
            if (process_id == 2)
                filter_data->delete_edge[edge] = true;
            return;
        }
        subtraction_optimized(g, candidate_vtx[process_id], process_id, s, matched_vetex);

        if (candidate_vtx[process_id].get_size() == 0) {
            if (process_id == 2)
                filter_data->delete_edge[edge] = true;
            return;
        }

        // if (candidate_vtx[process_id].get_size() != tmp) {
        //     printf("wrong");
        //     subtraction(g, tmp_set, matched_vetex_set);
        // }

        // 如果可以用数值计算代替集合计算
        if (s->compute_depth == process_id && s->compute_depth != s->pattern_size) {
            uint64_t size = candidate_vtx[process_id].get_size();
            uint32_t num = s->pattern_size - s->compute_depth;

            if (size < num) {
                return;
            }

            for (int i = 1; i <= num; i++) {
                temp_count = temp_count * size;
                size--;
            }

            for (int i = 2; i < process_id; i++) {
                temp_count *= candidate_vtx[i].get_size();
            }

            double temp_squra = temp_count * temp_count;

#pragma omp atomic
            filter_data->sample_count_per_edge[edge] += temp_count;

#pragma omp atomic
            filter_data->sample_square_per_edge[edge] += temp_squra;

#pragma omp atomic
            filter_data->sample_count_per_region[thread_id][region_id] += temp_count;

#pragma omp atomic
            filter_data->sample_square_per_region[thread_id][region_id] += temp_squra;

            return;
        }

        // m个顶点中有n个可行顶点，等概率从m中选择顶点，是否这n个顶点有相同的概率被选择？
        // TODO 随机选择是否存在优化的地方？

        if (candidate_vtx[process_id].get_size() >= 1) {
            std::uniform_int_distribution<uint64_t> dist(0, candidate_vtx[process_id].get_size() - 1);
            int random_id = dist(gen);
            int random_vtx = candidate_vtx[process_id].get_data(random_id);
            matched_vetex_set.push_back(random_vtx);
            matched_vetex[process_id] = random_vtx;

            if (g->degree[random_vtx] < s->order_degree[process_id]) {
                return;
            }

        } else {
            return;
        }
    }

    for (int i = 0; i < s->pattern_size; i++) {
        temp_count = temp_count * candidate_vtx[i].get_size();
    }
    double temp_squra = temp_count * temp_count;
#pragma omp atomic
    filter_data->sample_count_per_edge[edge] += temp_count;

#pragma omp atomic
    filter_data->sample_square_per_edge[edge] += temp_squra;

#pragma omp atomic
    filter_data->sample_count_per_region[thread_id][region_id] += temp_count;

#pragma omp atomic
    filter_data->sample_square_per_region[thread_id][region_id] += temp_squra;
}

uint64_t recur_match(Graph *g, schedule_approx *s, VertexSet *candidate_vtx, VertexSet &matched_vetex_set, uint32_t *matched_vetex, int process_id) {

    uint64_t temp_count = 0;

    // 获取依赖的第一个顶点集合
    if (s->dependcy_label[process_id] == -1) {
        int tmp_pre_vtx = matched_vetex[s->dependcy[s->dependcy_ptr[process_id]]];
        int tmp_vtx_set_size = g->vertex[tmp_pre_vtx + 1] - g->vertex[tmp_pre_vtx];
        candidate_vtx[process_id].copy(tmp_vtx_set_size, g->edge + g->vertex[tmp_pre_vtx]);
        // 与后面的做交集
        VertexSet tmp_vtx_set;
        for (int j = s->dependcy_ptr[process_id] + 1; j < s->dependcy_ptr[process_id + 1]; j++) {
            tmp_pre_vtx = matched_vetex[s->dependcy[j]];
            tmp_vtx_set_size = g->vertex[tmp_pre_vtx + 1] - g->vertex[tmp_pre_vtx];
            tmp_vtx_set.init(tmp_vtx_set_size, g->edge + g->vertex[tmp_pre_vtx]);
            candidate_vtx[process_id].intersection_with(tmp_vtx_set);
        }
    } else {
        // 直接复用结果
        int reuse_id = s->dependcy[s->dependcy_ptr[process_id]];
        candidate_vtx[process_id].copy(candidate_vtx[reuse_id].get_size(), candidate_vtx[reuse_id].get_data_ptr());

        // 与后面的做交接
        VertexSet tmp_vtx_set;

        int tmp_pre_vtx;
        int tmp_vtx_set_size;
        for (int j = s->dependcy_ptr[process_id] + 1; j < s->dependcy_ptr[process_id + 1]; j++) {
            if (s->dependcy[j] != -1) {
                tmp_pre_vtx = matched_vetex[s->dependcy[j]];
                tmp_vtx_set_size = g->vertex[tmp_pre_vtx + 1] - g->vertex[tmp_pre_vtx];
                tmp_vtx_set.init(tmp_vtx_set_size, g->edge + g->vertex[tmp_pre_vtx]);
                candidate_vtx[process_id].intersection_with(tmp_vtx_set);
            }
        }
    }

    // 如果当前层是最后一层
    if (process_id == s->pattern_size - 1) {
        temp_count = candidate_vtx[process_id].unordered_subtraction_size(candidate_vtx[process_id], matched_vetex_set);
        return temp_count;
    }

    // subtraction(g, candidate_vtx[process_id], matched_vetex_set);
    subtraction_optimized(g, candidate_vtx[process_id], process_id, s, matched_vetex);

    // 如果可以用数值计算代替集合计算
    // if (s->compute_depth == process_id) {
    //     int size = candidate_vtx[process_id].get_size();

    //     for (int i = 0; i < s->pattern_size - s->compute_depth; i++) {
    //         temp_count = temp_count * size / (i + 1);
    //         size--;
    //     }
    //     sample_count = sample_count + temp_count;
    //     return;
    // }

    // 对于每一个匹配顶点进行处理

    int *iter_data = candidate_vtx[process_id].get_data_ptr();
    for (int i = 0; i < candidate_vtx[process_id].get_size(); i++) {
        if (g->degree[iter_data[i]] < s->order_degree[process_id]) {
            continue;
        }
        matched_vetex_set.push_back(iter_data[i]);
        matched_vetex[process_id] = iter_data[i];
        temp_count = temp_count + recur_match(g, s, candidate_vtx, matched_vetex_set, matched_vetex, process_id + 1);
        matched_vetex_set.pop_back();
    }
    return temp_count;
}

void pattern_sample_es(uint32_t edge, Graph *g, schedule_approx *s, double &sample_count) {

    // 接收到的两个顶点
    uint32_t left_vtx = g->edge_from[edge];
    uint32_t right_vtx = g->edge[edge];

    // 已经匹配的顶点集合
    VertexSet matched_vetex_set;
    uint32_t *matched_vetex = new uint32_t[s->pattern_size];
    matched_vetex_set.init();
    matched_vetex_set.push_back(left_vtx);
    matched_vetex_set.push_back(right_vtx);
    matched_vetex[0] = left_vtx;
    matched_vetex[1] = right_vtx;

    int temp_count = 0;
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

    int process_id = 2;

    // 构架第三个顶点的候选集
    int tmp_pre_vtx = matched_vetex[s->dependcy[s->dependcy_ptr[process_id]]];
    tmp_vtx_set_size = g->vertex[tmp_pre_vtx + 1] - g->vertex[tmp_pre_vtx];
    candidate_vtx[process_id].copy(tmp_vtx_set_size, g->edge + g->vertex[tmp_pre_vtx]);

    // 与后面的做交集
    for (int j = s->dependcy_ptr[process_id] + 1; j < s->dependcy_ptr[process_id + 1]; j++) {
        tmp_pre_vtx = matched_vetex[s->dependcy[j]];
        tmp_vtx_set_size = g->vertex[tmp_pre_vtx + 1] - g->vertex[tmp_pre_vtx];
        tmp_vtx_set.init(tmp_vtx_set_size, g->edge + g->vertex[tmp_pre_vtx]);
        candidate_vtx[process_id].intersection_with(tmp_vtx_set);
    }

    // subtraction(g, candidate_vtx[process_id], matched_vetex_set);
    subtraction_optimized(g, candidate_vtx[process_id], process_id, s, matched_vetex);

    // 对于每一个匹配顶点进行处理

    int *iter_data = candidate_vtx[process_id].get_data_ptr();
    for (int i = 0; i < candidate_vtx[process_id].get_size(); i++) {
        if (g->degree[iter_data[i]] < s->order_degree[process_id]) {
            continue;
        }
        matched_vetex_set.push_back(iter_data[i]);
        matched_vetex[process_id] = iter_data[i];
        temp_count = temp_count + recur_match(g, s, candidate_vtx, matched_vetex_set, matched_vetex, process_id + 1);
        matched_vetex_set.pop_back();
    }

    sample_count += temp_count;
}

void edge_from_init(Graph *g) {
    g->edge_from = (v_index_t *)malloc(sizeof(v_index_t) * g->e_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        for (uint64_t j = g->vertex[i]; j < g->vertex[i + 1]; j++) {
            g->edge_from[j] = i;
        }
    }
}

void pattern_count_ns(Graph *g, schedule_approx *s) {

    // 初始化阈值和记录边的端点
    int dense_threshold = g->e_cnt / g->v_cnt; // 平均度数  frendster 28 twitter 35 livejouranl 17 youtube 2 patent 5 pokec 30 tp 28 sk 11
    int sparse_threshold = 0;
    sparse_threshold = threshhold_computing(g);
    int threshold = sparse_threshold;
    double count = 0;

    edge_from_init(g);

    // 过滤并记录剩余边
    filter filter_data;
    filter_data.edge_project = (int64_t *)malloc(sizeof(int64_t) * g->e_cnt);
    filter_data.record = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
    filter_data.edge_size = 0;
    memset(filter_data.edge_project, -1, sizeof(int64_t) * g->e_cnt);
    memset(filter_data.record, -1, sizeof(uint32_t) * g->e_cnt);
    filter_data.degree = (uint32_t *)malloc(sizeof(uint32_t) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        filter_data.degree[i] = g->vertex[i + 1] - g->vertex[i];
    }

    if (s->symmetry_label) {
        // for (e_index_t i = 0; i < g->e_cnt; ++i) {
        //     v_index_t vtx_left = g->edge_from[i];
        //     v_index_t vtx_right = g->edge[i];
        //     if (filter_data.degree[vtx_left] < threshold && filter_data.degree[vtx_right] < threshold && vtx_left < vtx_right) {
        //         if (filter_data.degree[vtx_left] < s->order_degree[0] || filter_data.degree[vtx_right] < s->order_degree[1]) {
        //             filter_data.record[i] = 1;
        //         }
        //         pattern_sample_es(i, g, s, count);
        //         filter_data.record[i] = 1;
        //     }
        // }
        for (e_index_t i = 0; i < g->e_cnt; ++i) {
            if (filter_data.record[i] == -1 && g->edge_from[i] < g->edge[i]) {
                filter_data.edge_project[filter_data.edge_size++] = i;
            }
        }
    } else {
        for (e_index_t i = 0; i < g->e_cnt; ++i) {
            v_index_t vtx_left = g->edge_from[i];
            v_index_t vtx_right = g->edge[i];
            if (filter_data.degree[vtx_left] < threshold && filter_data.degree[vtx_right] < threshold) {
                if (filter_data.degree[vtx_left] < s->order_degree[0] || filter_data.degree[vtx_right] < s->order_degree[1]) {
                    filter_data.record[i] = 1;
                }
                pattern_sample_es(i, g, s, count);
                filter_data.record[i] = 1;
            }
        }
        for (e_index_t i = 0; i < g->e_cnt; ++i) {
            if (filter_data.record[i] == -1) {
                filter_data.edge_project[filter_data.edge_size++] = i;
            }
        }
    }

    /*随机采样*/

    // 构建采样顺序
    std::random_device rd;
    std::default_random_engine gen(rd());
    uint32_t *order = (uint32_t *)malloc(sizeof(uint32_t) * filter_data.edge_size);
    for (uint32_t i = 0; i < filter_data.edge_size; i++) {
        order[i] = i;
    }
    std::shuffle(order, order + filter_data.edge_size, gen);

    // 记录采样结果
    int record = 0;
    double estimate[100] = {0};
    double error[100] = {0};

    // 混合采样
    // for (int i = 0; i < filter_data.edge_size; i++) {
    //     auto eid = filter_data.edge_project[order[i]];
    //     pattern_sample_ns(eid, g, s, count_sample);
    //     if (std::pow(10, record) == i) {
    //         estimate[record] = count_sample * filter_data.edge_size / i;
    //         record++;
    //     }
    // }
    // estimate[record] = count_sample;

    // std::uniform_int_distribution<int> dist(0, filter_data.edge_size - 1);

    std::uniform_int_distribution<uint64_t> dist(0, g->e_cnt - 1);

    // 随机采样100万次。
    int iter_num[10] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000};

    while (record < 7) {
        double global_ans = 0;
        // #pragma omp parallel reduction(+ : global_ans)
        {
            double count_sample = 0;
            double count_sample_square = 0;
            // #pragma omp for schedule(dynamic) nowait
            for (int i = 0; i < iter_num[record]; i++) {
                // auto eid = filter_data.edge_project[order[dist(gen)]];
                auto eid = dist(gen);
                pattern_sample_ns(eid, g, s, count_sample, count_sample_square);
                // pattern_sample_es(eid, g, s, count_sample);
            }
            global_ans += count_sample;
        }

        // estimate[record] = global_ans * filter_data.edge_size / (std::pow(10, record));
        estimate[record] = global_ans * g->e_cnt / std::pow(10, record);
        record++;
    } // bug除错了，应该每次都重新更新count_sample

    // for (int i = 0; i < 10; i++) {
    //     for (int j = 0; j < 1000000; j++) {
    //         auto eid = filter_data.edge_project[dist(gen)];
    //         pattern_sample_ns(eid, g, s, count_sample);
    //     }
    //     estimate[i] = count_sample * filter_data.edge_size / ((i + 1) * 1000000);
    // }

    // 输出结果
    for (int i = 0; i < 7; i++) {
        std::cout << "triangletriangle计数为 " << estimate[i] + count << std::endl;
        error[i] = (estimate[i] + double(count) - (53073844144.0 * 3.0)) / (53073844144.0 * 3.0); // house 53073844144  4986965  对三角形43378336949
        std::cout << "误差为 " << error[i] << std::endl;
    }

    free(filter_data.edge_project);
    free(filter_data.record);
    free(filter_data.degree);
}

void pattern_count_ns_with_symmetry_break(Graph *g, schedule_approx *s) {

    // 初始化阈值和记录边的端点
    int dense_threshold = g->e_cnt / g->v_cnt; // 平均度数  frendster 28 twitter 35 livejouranl 17 youtube 2 patent 5 pokec 30 tp 28 sk 11
    int sparse_threshold = 0;
    sparse_threshold = threshhold_computing(g);
    int threshold = sparse_threshold;
    double count = 0;

    edge_from_init(g);

    // 过滤并记录剩余边
    filter filter_data;
    filter_data.edge_project = (int64_t *)malloc(sizeof(int64_t) * g->e_cnt);
    filter_data.record = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
    filter_data.edge_size = 0;
    memset(filter_data.edge_project, -1, sizeof(int64_t) * g->e_cnt);
    memset(filter_data.record, -1, sizeof(uint32_t) * g->e_cnt);
    filter_data.degree = (uint32_t *)malloc(sizeof(uint32_t) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        filter_data.degree[i] = g->vertex[i + 1] - g->vertex[i];
    }

    if (s->symmetry_label) {
        // for (e_index_t i = 0; i < g->e_cnt; ++i) {
        //     v_index_t vtx_left = g->edge_from[i];
        //     v_index_t vtx_right = g->edge[i];
        //     if (filter_data.degree[vtx_left] < threshold && filter_data.degree[vtx_right] < threshold && vtx_left < vtx_right) {
        //         if (filter_data.degree[vtx_left] < s->order_degree[0] || filter_data.degree[vtx_right] < s->order_degree[1]) {
        //             filter_data.record[i] = 1;
        //         }
        //         pattern_sample_es(i, g, s, count);
        //         filter_data.record[i] = 1;
        //     }
        // }
        for (e_index_t i = 0; i < g->e_cnt; ++i) {
            if (filter_data.record[i] == -1 && g->edge_from[i] < g->edge[i]) {
                filter_data.edge_project[filter_data.edge_size++] = i;
            }
        }
    } else {
        for (e_index_t i = 0; i < g->e_cnt; ++i) {
            v_index_t vtx_left = g->edge_from[i];
            v_index_t vtx_right = g->edge[i];
            if (filter_data.degree[vtx_left] < threshold && filter_data.degree[vtx_right] < threshold) {
                if (filter_data.degree[vtx_left] < s->order_degree[0] || filter_data.degree[vtx_right] < s->order_degree[1]) {
                    filter_data.record[i] = 1;
                }
                pattern_sample_es(i, g, s, count);
                filter_data.record[i] = 1;
            }
        }
        for (e_index_t i = 0; i < g->e_cnt; ++i) {
            if (filter_data.record[i] == -1) {
                filter_data.edge_project[filter_data.edge_size++] = i;
            }
        }
    }

    /*随机采样*/

    // 构建采样顺序
    std::random_device rd;
    std::default_random_engine gen(rd());
    uint32_t *order = (uint32_t *)malloc(sizeof(uint32_t) * filter_data.edge_size);
    for (uint32_t i = 0; i < filter_data.edge_size; i++) {
        order[i] = i;
    }
    std::shuffle(order, order + filter_data.edge_size, gen);

    // 记录采样结果
    int record = 0;
    double estimate[100] = {0};
    double error[100] = {0};

    // 混合采样
    // for (int i = 0; i < filter_data.edge_size; i++) {
    //     auto eid = filter_data.edge_project[order[i]];
    //     pattern_sample_ns(eid, g, s, count_sample);
    //     if (std::pow(10, record) == i) {
    //         estimate[record] = count_sample * filter_data.edge_size / i;
    //         record++;
    //     }
    // }
    // estimate[record] = count_sample;

    std::uniform_int_distribution<uint64_t> dist(0, filter_data.edge_size - 1);

    // 随机采样100万次。
    int iter_num[10] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000};

    while (record < 7) {
        double global_ans = 0;
        // #pragma omp parallel reduction(+ : global_ans)
        {
            double count_sample = 0;
            double count_sample_square = 0;
            // #pragma omp for schedule(dynamic) nowait
            for (int i = 0; i < iter_num[record]; i++) {
                auto eid = filter_data.edge_project[order[dist(gen)]];
                pattern_sample_ns(eid, g, s, count_sample, count_sample_square);
            }
            global_ans += count_sample;
        }

        estimate[record] = global_ans * filter_data.edge_size / (std::pow(10, record));
        record++;
    } // bug除错了，应该每次都重新更新count_sample

    // for (int i = 0; i < 10; i++) {
    //     for (int j = 0; j < 1000000; j++) {
    //         auto eid = filter_data.edge_project[dist(gen)];
    //         pattern_sample_ns(eid, g, s, count_sample);
    //     }
    //     estimate[i] = count_sample * filter_data.edge_size / ((i + 1) * 1000000);
    // }

    // 输出结果
    for (int i = 0; i < 7; i++) {
        std::cout << "triangletriangle计数为 " << estimate[i] + count << std::endl;
        error[i] = (estimate[i] + double(count) - (53073844144.0)) / (53073844144.0); // house 53073844144  4986965  对三角形43378336949
        std::cout << "误差为 " << error[i] << std::endl;
    }

    free(filter_data.edge_project);
    free(filter_data.record);
    free(filter_data.degree);
}

void pattern_count_ns_partition(Graph *g, schedule_approx *s) {

    // 初始化阈值和记录边的端点
    int dense_threshold = g->e_cnt / g->v_cnt; // 平均度数  frendster 28 twitter 35 livejouranl 17 youtube 2 patent 5 pokec 30 tp 28 sk 11
    int sparse_threshold = 0;
    sparse_threshold = threshhold_computing(g);
    int threshold = sparse_threshold;
    double count = 0;

    edge_from_init(g);

    // 过滤并记录剩余边
    filter filter_data;
    filter_data.edge_project = (int64_t *)malloc(sizeof(int64_t) * g->e_cnt);
    filter_data.record = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
    filter_data.edge_size = 0;
    memset(filter_data.edge_project, -1, sizeof(int64_t) * g->e_cnt);
    memset(filter_data.record, -1, sizeof(uint32_t) * g->e_cnt);
    filter_data.degree = (uint32_t *)malloc(sizeof(uint32_t) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        filter_data.degree[i] = g->vertex[i + 1] - g->vertex[i];
    }

    if (s->symmetry_label) {
        // for (e_index_t i = 0; i < g->e_cnt; ++i) {
        //     v_index_t vtx_left = g->edge_from[i];
        //     v_index_t vtx_right = g->edge[i];
        //     if (filter_data.degree[vtx_left] < threshold && filter_data.degree[vtx_right] < threshold && vtx_left < vtx_right) {
        //         if (filter_data.degree[vtx_left] < s->order_degree[0] || filter_data.degree[vtx_right] < s->order_degree[1]) {
        //             filter_data.record[i] = 1;
        //         }
        //         pattern_sample_es(i, g, s, count);
        //         filter_data.record[i] = 1;
        //     }
        // }
        for (e_index_t i = 0; i < g->e_cnt; ++i) {
            if (filter_data.record[i] == -1 && g->edge_from[i] < g->edge[i]) {
                filter_data.edge_project[filter_data.edge_size++] = i;
            }
        }
    } else {
        for (e_index_t i = 0; i < g->e_cnt; ++i) {
            v_index_t vtx_left = g->edge_from[i];
            v_index_t vtx_right = g->edge[i];
            if (filter_data.degree[vtx_left] < threshold && filter_data.degree[vtx_right] < threshold) {
                if (filter_data.degree[vtx_left] < s->order_degree[0] || filter_data.degree[vtx_right] < s->order_degree[1]) {
                    filter_data.record[i] = 1;
                }
                pattern_sample_es(i, g, s, count);
                filter_data.record[i] = 1;
            }
        }
        for (e_index_t i = 0; i < g->e_cnt; ++i) {
            if (filter_data.record[i] == -1) {
                filter_data.edge_project[filter_data.edge_size++] = i;
            }
        }
    }

    /*随机采样*/

    // 构建采样顺序
    std::random_device rd;
    std::default_random_engine gen(rd());
    uint32_t *order = (uint32_t *)malloc(sizeof(uint32_t) * filter_data.edge_size);
    for (uint32_t i = 0; i < filter_data.edge_size; i++) {
        order[i] = i;
    }
    std::shuffle(order, order + filter_data.edge_size, gen);

    // 记录采样结果
    int record = 0;
    double estimate[100] = {0};
    double error[100] = {0};
    // 随机采样100万次。
    int iter_num[10] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000};
    std::uniform_int_distribution<uint64_t> dist1(0, (filter_data.edge_size - 1) / 5);
    std::uniform_int_distribution<uint64_t> dist2(0, (filter_data.edge_size - 1) / 5 * 1.5);
    std::uniform_int_distribution<uint64_t> dist3(0, (filter_data.edge_size - 1) / 5 * 0.5);
    std::uniform_int_distribution<uint64_t> dist4(0, (filter_data.edge_size - 1) / 5);
    std::uniform_int_distribution<uint64_t> dist5(0, (filter_data.edge_size - 1) / 5);

    while (record < 7) {
        uint64_t global_ans = 0;
        // #pragma omp parallel reduction(+ : global_ans)
        {
            double count_sample1 = 0;
            double count_sample2 = 0;
            double count_sample3 = 0;
            double count_sample4 = 0;
            double count_sample5 = 0;
            double count_sample_square = 0;
            // #pragma omp for schedule(dynamic) nowait
            for (int i = 0; i < iter_num[record] / 5; i++) {
                auto eid1 = filter_data.edge_project[order[dist1(gen)]];
                pattern_sample_ns(eid1, g, s, count_sample1, count_sample_square);
                auto eid2 = filter_data.edge_project[order[dist2(gen) + (filter_data.edge_size / 5)]];
                pattern_sample_ns(eid2, g, s, count_sample2, count_sample_square);
                auto eid3 = filter_data.edge_project[order[dist3(gen) + (filter_data.edge_size / 2)]];
                pattern_sample_ns(eid3, g, s, count_sample3, count_sample_square);
                auto eid4 = filter_data.edge_project[order[dist4(gen) + (filter_data.edge_size / 5 * 3)]];
                pattern_sample_ns(eid4, g, s, count_sample4, count_sample_square);
                auto eid5 = filter_data.edge_project[order[dist5(gen) + (filter_data.edge_size / 5 * 4)]];
                pattern_sample_ns(eid5, g, s, count_sample5, count_sample_square);
            }

            count_sample1 = count_sample1 * (filter_data.edge_size) / ((std::pow(10, record)));
            count_sample2 = count_sample2 * (filter_data.edge_size * 1.5) / ((std::pow(10, record)));
            count_sample3 = count_sample3 * (filter_data.edge_size * 0.5) / ((std::pow(10, record)));
            count_sample4 = count_sample4 * (filter_data.edge_size) / ((std::pow(10, record)));
            count_sample5 = count_sample5 * (filter_data.edge_size) / ((std::pow(10, record)));
            global_ans += count_sample1 + count_sample2 + count_sample3 + count_sample4 + count_sample5;
        }
        // estimate[record] = global_ans * filter_data.edge_size / (std::pow(10, record));
        estimate[record] = global_ans;
        record++;
    } // bug除错了，应该每次都重新更新count_sample

    // 输出结果
    for (int i = 0; i < 7; i++) {
        std::cout << "triangletriangle计数为 " << estimate[i] + count << std::endl;
        error[i] = (estimate[i] + double(count) - (53073844144.0)) / (53073844144.0); // house 53073844144  4986965  对三角形43378336949
        std::cout << "误差为 " << error[i] << std::endl;
    }

    free(filter_data.edge_project);
    free(filter_data.record);
    free(filter_data.degree);
}

void pattern_count_es(Graph *g, schedule_approx *s) {

    // 初始化阈值和记录边的端点
    int dense_threshold = g->e_cnt / g->v_cnt; // 平均度数  frendster 28 twitter 35 livejouranl 17 youtube 2 patent 5 pokec 30 tp 28 sk 11
    int sparse_threshold = 0;
    sparse_threshold = threshhold_computing(g);
    int threshold = std::sqrt(g->v_cnt);
    double count = 0;

    edge_from_init(g);

    // 过滤并记录剩余边
    filter filter_data;
    filter_data.edge_project = (int64_t *)malloc(sizeof(int64_t) * g->e_cnt);
    filter_data.record = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
    filter_data.edge_size = 0;
    memset(filter_data.edge_project, -1, sizeof(int64_t) * g->e_cnt);
    memset(filter_data.record, -1, sizeof(uint32_t) * g->e_cnt);
    filter_data.degree = (uint32_t *)malloc(sizeof(uint32_t) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        filter_data.degree[i] = g->vertex[i + 1] - g->vertex[i];
    }

    if (s->symmetry_label) {
        for (e_index_t i = 0; i < g->e_cnt; ++i) {
            v_index_t vtx_left = g->edge_from[i];
            v_index_t vtx_right = g->edge[i];
            if (filter_data.degree[vtx_left] < threshold && filter_data.degree[vtx_right] < threshold && vtx_left < vtx_right) {
                if (filter_data.degree[vtx_left] < s->order_degree[0] || filter_data.degree[vtx_right] < s->order_degree[1]) {
                    filter_data.record[i] = 1;
                    continue;
                }
                pattern_sample_es(i, g, s, count);
                filter_data.record[i] = 1;
            }
        }
        for (e_index_t i = 0; i < g->e_cnt; ++i) {
            if (filter_data.record[i] == -1 && g->edge_from[i] < g->edge[i]) {
                filter_data.edge_project[filter_data.edge_size++] = i;
            }
        }
    } else {
        for (e_index_t i = 0; i < g->e_cnt; ++i) {
            v_index_t vtx_left = g->edge_from[i];
            v_index_t vtx_right = g->edge[i];
            if (filter_data.degree[vtx_left] < threshold && filter_data.degree[vtx_right] < threshold) {
                if (filter_data.degree[vtx_left] < s->order_degree[0] || filter_data.degree[vtx_right] < s->order_degree[1]) {
                    filter_data.record[i] = 1;
                    continue;
                }
                pattern_sample_es(i, g, s, count);
                filter_data.record[i] = 1;
            }
        }
        for (e_index_t i = 0; i < g->e_cnt; ++i) {
            if (filter_data.record[i] == -1) {
                filter_data.edge_project[filter_data.edge_size++] = i;
            }
        }
    }

    /*随机采样*/

    // 构建采样顺序
    std::random_device rd;
    std::default_random_engine gen(rd());
    uint32_t *order = (uint32_t *)malloc(sizeof(uint32_t) * filter_data.edge_size);
    for (uint32_t i = 0; i < filter_data.edge_size; i++) {
        order[i] = i;
    }
    std::shuffle(order, order + filter_data.edge_size, gen);

    // 记录采样结果
    int record = 1;
    double count_sample = 0;
    double estimate[20] = {0};
    double error[20] = {0};

    // 混合采样
    for (int i = 0; i < filter_data.edge_size; i++) {
        auto eid = filter_data.edge_project[order[i]];
        pattern_sample_es(eid, g, s, count_sample);
        if (std::pow(10, record) == i) {
            estimate[record] = count_sample * filter_data.edge_size / i;
            record++;
        }
    }
    estimate[record] = count_sample;

    // 随机采样。
    // std::uniform_int_distribution<int> dist(0, filter_data.edge_size - 1);
    // for (int i = 0; i < 10; i++) {
    //     for (int j = 0; j < 1000; j++) {
    //         auto eid = filter_data.edge_project[order[dist(gen)]];
    //         pattern_sample_es(eid, g, s, count_sample);
    //     }
    //     estimate[i] = count_sample * filter_data.edge_size / ((i + 1) * 1000);
    // }

    // 输出结果
    for (int i = 1; i < 20; i++) {
        std::cout << "clique 计数为 " << estimate[i] + count << std::endl;
        error[i] = (estimate[i] + double(count) - (43378336949.0 * 3.0)) /
                   (43378336949.0 * 3.0); // house 53073844144 4clique 4986965 5clique 7211947   对三角形 43378336949
        std::cout << "误差为 " << error[i] << std::endl;
    }

    free(filter_data.edge_project);
    free(filter_data.record);
    free(filter_data.degree);
}

#define BUFFER_SIZE 2000
void pattern_count_ns_mpi_fix_sample_times(Graph *g, schedule_approx *s) {

    // 获取边的起点
    edge_from_init(g);
    // 分布式执行

    // 采样批量
    uint64_t sampling_step = 10000;

    // 全局变量
    uint64_t global_sampled_times = 0;
    uint64_t global_prob_inverse = 0;

    int comm_sz, my_rank, num_processes;

    MPI_Request recvrqst[BUFFER_SIZE];
    MPI_Status status;
    uint64_t recv_prob_inverse[BUFFER_SIZE] = {0};
    uint32_t headptr = 0, tailptr = 0;

    uint32_t sampling_times = 1000000;

    double count_sample = 0;
    double count_sample_square = 0;

#pragma omp parallel num_threads(2)
    {
#pragma omp barrier
#pragma omp master // 首先每个进程的主线程执行 初始化操作
        {
            int provided;
            MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &provided);
            MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
            MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
            MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
            MPI_Barrier(MPI_COMM_WORLD);
        }
#pragma omp barrier

#pragma omp master
        { // 当预测结果缓存中还有空间时，接收一个预测结果。预测结果缓存中有没处理的结果时，接收结果，更新到全局，然后更新采样次数。
          // 主节点主线程只管自己接收到了多少预测结果，以及相关的预测次数，不需要管实际的。    因此工作线程中也可以写为死循环。
            while (global_sampled_times < sampling_times) {
                if (!my_rank) {
                    if (tailptr - headptr < BUFFER_SIZE) {
                        MPI_Irecv(&recv_prob_inverse[tailptr % BUFFER_SIZE], 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,
                                  &recvrqst[tailptr % BUFFER_SIZE]);
                        tailptr++;
                    }
                    if (headptr < tailptr) // 如果缓存中还有数据 ，将缓存中的数据接收并且更新整体
                    {
                        if (global_sampled_times >= sampling_times)
                            break;
                        int testflag = 0;
                        MPI_Test(&recvrqst[headptr % BUFFER_SIZE], &testflag, MPI_STATUS_IGNORE);
                        if (testflag) {
#pragma omp atomic
                            global_prob_inverse += recv_prob_inverse[tailptr % BUFFER_SIZE];
#pragma omp atomic
                            global_sampled_times += sampling_step;
                            headptr++;
                        }
                    }
                }
                if (my_rank) // 非主进程的主线程的工作：做一些局部的统计工作，算ok了传递给主进程。
                    MPI_Recv(&global_sampled_times, 1, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD, &status);
                // break;
                // printf("caculate local result");
            }
            if (!my_rank) {
                for (uint32_t i = 1; i < num_processes; i++)
                    MPI_Send(&global_sampled_times, 1, MPI_UNSIGNED_LONG_LONG, i, 0, MPI_COMM_WORLD);
            }
        }

        // 每个进程的非主线程执行的任务，计算完更新主节点的主线程
        if (omp_get_thread_num()) {
            MPI_Request sendrqst[BUFFER_SIZE]; // 用于标记某个发送请求是否已经发送完毕

            // 每个线程 保存自己的预测数量和采样次数
            uint64_t local_prob_inverse[BUFFER_SIZE] = {0};
            uint32_t headptr = 0, tailptr = 0;

            // 随机数生成器
            std::random_device rd;
            std::default_random_engine gen(rd());
            std::uniform_int_distribution<uint64_t> dist(0, g->e_cnt - 1);

            // 如果采样次数不够，持续采样并且发送本地预测结果 ， 构建了一个缓存，缓存满了才进行发送，异步的进行发送
            // 每个线程都会去更新主节点？

            // 死循环，直到运算误差小于某个值
            // 与主节点的采样次数一致，会同步更新，因此会退出
            while (global_sampled_times < sampling_times) {
                if (tailptr - headptr < BUFFER_SIZE) // 发送缓存存在空间，继续采样并发送
                {
                    // 按照固定的batch进行测量
                    count_sample = 0;
                    for (int i = 0; i < sampling_step; i++) {
                        auto eid = dist(gen);
                        pattern_sample_ns(eid, g, s, count_sample, count_sample_square);
                    }
                    local_prob_inverse[tailptr % BUFFER_SIZE] = count_sample;

                    MPI_Isend(&local_prob_inverse[tailptr % BUFFER_SIZE], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
                              &sendrqst[tailptr % BUFFER_SIZE]); // 最后一个位置时标记，标记该操作是否已经完成
                    tailptr++;
                }
                if (headptr < tailptr) // 如果发送缓存不够用了，查看是否已经发送完成了。 一个异步的更新过程。
                {
                    int testflag = 0;
                    MPI_Test(&sendrqst[headptr % BUFFER_SIZE], &testflag, MPI_STATUS_IGNORE); // 异步的操作需要一个test来执行
                    if (testflag)
                        headptr++;
                    if (global_sampled_times >= sampling_times)
                        break;
                }

                // if (my_rank && headptr > 30) {
                //     break;
                // }
            }
            std::cout << "work thread " << "headptr " << headptr << " tailptr " << tailptr << std::endl;
            // for (int i = 0; i < 100; i++) {
            //     if (local_prob_inverse[i] != recv_prob_inverse[i]) {
            //         printf("wrong");
            //     }
            // }
        }
    }

    // 线程并行执行区域结束，这里只有主进程主线程会执行
    if (!my_rank) {
        // for (int i = 0; i < 100; i++) {
        //     if (recv_prob_inverse[i] == 0) {
        //         printf("wrong");
        //     }
        // }

        std::cout << "headptr " << headptr << " tailptr " << tailptr << std::endl;
        std::cout << "global summary global_sample_times: " << global_prob_inverse << " " << global_sampled_times << std::endl;
        uint64_t estimated_pattern_count = global_prob_inverse * g->e_cnt / (global_sampled_times);
        std::cout << "pattern_nums " << estimated_pattern_count << std::endl;
        double estimated_pattern_count_error = (estimated_pattern_count - (53073844144.0 * 3.0)) / (53073844144.0 * 3.0);
        std::cout << "estimated_error " << estimated_pattern_count_error << std::endl;
    }

    MPI_Finalize();
}

void pattern_count_ns_mpi_fix_error(Graph *g, schedule_approx *s) {

    // 获取边的起点
    edge_from_init(g);

    // 分布式执行

    // 采样批量
    uint64_t sampling_step = 10000;

    // 全局变量
    uint64_t global_sampled_times = 0;
    uint64_t global_sample_count = 0;
    uint64_t global_sample_squra = 0;

    bool error_bound_satisfied = false;
    double var = 0.0;
    double std_dev = 0.0;
    double real_error = 0.0;
    double error_bound = 0.05;

    int comm_sz, my_rank, num_processes;

    MPI_Status status;
    uint64_t recv_sample_count[BUFFER_SIZE] = {0};
    uint64_t recv_sample_squra[BUFFER_SIZE] = {0};
    MPI_Request recv_sample_rqst[BUFFER_SIZE];
    MPI_Request recv_square_rqst[BUFFER_SIZE];
    uint32_t headptr = 0, tailptr = 0;

    uint32_t sampling_times = 1000000;

    uint64_t count_sample = 0;
    uint64_t count_sample_square = 0;

#pragma omp parallel num_threads(2)
    {
#pragma omp barrier
#pragma omp master // 首先每个进程的主线程执行 初始化操作
        {
            int provided;
            MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &provided);
            MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
            MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
            MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
            MPI_Barrier(MPI_COMM_WORLD);
        }
#pragma omp barrier

#pragma omp master
        { // 当预测结果缓存中还有空间时，接收一个预测结果。预测结果缓存中有没处理的结果时，接收结果，更新到全局，然后更新采样次数。
          // 主节点主线程只管自己接收到了多少预测结果，以及相关的预测次数，不需要管实际的。    因此工作线程中也可以写为死循环。
            while (!error_bound_satisfied) {
                if (!my_rank) {
                    if (tailptr - headptr < BUFFER_SIZE) {
                        MPI_Irecv(&recv_sample_count[tailptr % BUFFER_SIZE], 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,
                                  &recv_sample_rqst[tailptr % BUFFER_SIZE]);
                        MPI_Irecv(&recv_sample_squra[tailptr % BUFFER_SIZE], 1, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD,
                                  &recv_square_rqst[tailptr % BUFFER_SIZE]);
                        tailptr++;
                    }
                    if (headptr < tailptr) // 如果缓存中还有数据 ，将缓存中的数据接收并且更新整体
                    {
                        int testflag1 = 0;
                        int testflag2 = 0;
                        MPI_Test(&recv_sample_rqst[headptr % BUFFER_SIZE], &testflag1, MPI_STATUS_IGNORE);
                        MPI_Test(&recv_square_rqst[headptr % BUFFER_SIZE], &testflag2, MPI_STATUS_IGNORE);
                        if (testflag1 && testflag2) {
                            // #pragma omp atomic
                            global_sample_count += recv_sample_count[tailptr % BUFFER_SIZE];
                            global_sample_squra += recv_sample_squra[tailptr % BUFFER_SIZE];
                            // #pragma omp atomic
                            global_sampled_times += sampling_step;
                            headptr++;

                            // 计算全体方差和每个值得标准差
                            var = (double)(global_sample_squra / global_sampled_times - std::pow((global_sample_count / global_sampled_times), 2));
                            std_dev = std::sqrt(var / global_sampled_times);

                            real_error = std_dev * 2.576 / (global_sample_count / global_sampled_times); // 99%的信心
                            if (real_error < error_bound) {
                                error_bound_satisfied = true;
                            }
                        }

                        if (error_bound_satisfied) {
                            for (uint32_t i = 1; i < num_processes; i++)
                                MPI_Send(&error_bound_satisfied, 1, MPI_C_BOOL, i, 0, MPI_COMM_WORLD);
                            break;
                        }
                    }
                }
                if (my_rank) // 非主进程的主线程的工作：做一些局部的统计工作，算ok了传递给主进程。
                    MPI_Recv(&error_bound_satisfied, 1, MPI_C_BOOL, 0, 0, MPI_COMM_WORLD, &status);
                // break;
                // printf("caculate local result");
            }
        }

        // 每个进程的非主线程执行的任务，计算完更新主节点的主线程
        if (omp_get_thread_num()) {
            MPI_Request send_count_rqst[BUFFER_SIZE]; // 用于标记某个发送请求是否已经发送完毕
            MPI_Request send_squre_rqst[BUFFER_SIZE];

            // 每个线程 保存自己的预测数量和采样次数
            uint64_t local_count_sample[BUFFER_SIZE] = {0};
            uint64_t local_count_square[BUFFER_SIZE] = {0};
            uint32_t headptr = 0, tailptr = 0;

            // 随机数生成器
            std::random_device rd;
            std::default_random_engine gen(rd());
            std::uniform_int_distribution<uint64_t> dist(0, g->e_cnt - 1);

            // 如果采样次数不够，持续采样并且发送本地预测结果 ， 构建了一个缓存，缓存满了才进行发送，异步的进行发送
            // 每个线程都会去更新主节点？

            // 死循环，直到运算误差小于某个值
            // 与主节点的采样次数一致，会同步更新，因此会退出
            while (!error_bound_satisfied) {
                if (tailptr - headptr < BUFFER_SIZE) // 发送缓存存在空间，继续采样并发送
                {
                    // 按照固定的batch进行测量
                    count_sample = 0;
                    count_sample_square = 0;
                    uint64_t sampling_times = 0;
                    for (int i = 0; i < sampling_step; i++) {
                        auto eid = dist(gen);
                        // pattern_sample_ns(eid, g, s, count_sample, count_sample_square);
                    }

                    local_count_sample[tailptr % BUFFER_SIZE] = count_sample;
                    local_count_square[tailptr % BUFFER_SIZE] = count_sample_square;

                    MPI_Isend(&local_count_sample[tailptr % BUFFER_SIZE], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
                              &send_count_rqst[tailptr % BUFFER_SIZE]); // 最后一个位置时标记，标记该操作是否已经完成
                    MPI_Isend(&local_count_square[tailptr % BUFFER_SIZE], 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD,
                              &send_squre_rqst[tailptr % BUFFER_SIZE]);

                    tailptr++;
                }
                if (headptr < tailptr) // 如果发送缓存不够用了，查看是否已经发送完成了。 一个异步的更新过程。
                {
                    int testflag1 = 0;
                    int testflag2 = 0;
                    MPI_Test(&send_count_rqst[headptr % BUFFER_SIZE], &testflag1, MPI_STATUS_IGNORE); // 异步的操作需要一个test来执行
                    MPI_Test(&send_squre_rqst[headptr % BUFFER_SIZE], &testflag2, MPI_STATUS_IGNORE);
                    if (testflag1 & testflag2)
                        headptr++;
                    if (error_bound_satisfied)
                        break;
                }
                // if (my_rank && headptr > 30) {
                //     break;
                // }
            }
            std::cout << "work thread " << "headptr " << headptr << " tailptr " << tailptr << std::endl;
            // for (int i = 0; i < 100; i++) {
            //     if (local_prob_inverse[i] != recv_prob_inverse[i]) {
            //         printf("wrong");
            //     }
            // }
        }
    }

    // 线程并行执行区域结束，这里只有主进程主线程会执行
    if (!my_rank) {
        // for (int i = 0; i < 100; i++) {
        //     if (recv_prob_inverse[i] == 0) {
        //         printf("wrong");
        //     }
        // }

        std::cout << "headptr " << headptr << " tailptr " << tailptr << std::endl;
        std::cout << "global summary global_sample_times: " << global_sample_count << " " << global_sampled_times << std::endl;
        uint64_t estimated_pattern_count = global_sample_count * g->e_cnt / (global_sampled_times);
        std::cout << "pattern_nums " << estimated_pattern_count << std::endl;
        double estimated_pattern_count_error = (estimated_pattern_count - (53073844144.0 * 3.0)) / (53073844144.0 * 3.0);
        std::cout << "estimated_error " << estimated_pattern_count_error << std::endl;
    }

    MPI_Finalize();
}

void pattern_count_ns_mpi_fix_error_filter(Graph *g, schedule_approx *s) {

    // 获取边的起点
    edge_from_init(g);
    /*过滤操作*/
    // 初始化阈值
    int sparse_threshold = g->e_cnt / g->v_cnt;    // 平均度数  frendster 28 twitter 35 livejouranl 17 youtube 2 patent 5 pokec 30 tp 28 sk 11
    int dense_threshold = threshhold_computing(g); // 百分之九十的顶点的度数阈值
    if (dense_threshold < sparse_threshold) {
        int temp_threshold = dense_threshold;
        dense_threshold = sparse_threshold;
        sparse_threshold = temp_threshold;
    }
    sparse_threshold = 0;

    // 过滤的数据结构。
    filter filter_data;
    filter_data.edge_project = (int64_t *)malloc(sizeof(int64_t) * g->e_cnt);
    filter_data.record = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
    memset(filter_data.edge_project, -1, sizeof(int64_t) * g->e_cnt);
    memset(filter_data.record, -1, sizeof(uint32_t) * g->e_cnt);
    filter_data.degree = (uint32_t *)malloc(sizeof(uint32_t) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        filter_data.degree[i] = g->vertex[i + 1] - g->vertex[i];
    }

    uint64_t num_of_light_edges = 0;
    uint64_t num_of_medium_edges = 0;
    uint64_t num_of_heavy_edges = 0;

    uint64_t light_edge_ptr = 0;
    uint64_t medium_edge_ptr = 0;
    uint64_t heavy_edge_ptr = 0;
    if (s->symmetry_label) {
        // 标记边的种类
        for (e_index_t i = 0; i < g->e_cnt; ++i) {
            v_index_t vtx_left = g->edge_from[i];
            v_index_t vtx_right = g->edge[i];
            if (vtx_left < vtx_right) {
                if (filter_data.degree[vtx_left] >= s->order_degree[0] && filter_data.degree[vtx_right] >= s->order_degree[1]) {
                    if (filter_data.degree[vtx_left] < sparse_threshold && filter_data.degree[vtx_right] < sparse_threshold) {
                        filter_data.record[i] = 0;
                        num_of_light_edges++;
                    } else if (filter_data.degree[vtx_left] < dense_threshold && filter_data.degree[vtx_right] < dense_threshold) {
                        filter_data.record[i] = 1;
                        num_of_medium_edges++;
                    } else {
                        filter_data.record[i] = 2;
                        num_of_heavy_edges++;
                    }
                }
            }
        }
        medium_edge_ptr = num_of_light_edges;
        heavy_edge_ptr = num_of_light_edges + num_of_medium_edges;
        // 投射记录
        for (e_index_t i = 0; i < g->e_cnt; ++i) {
            if (filter_data.record[i] == 0) {
                filter_data.edge_project[light_edge_ptr++] = i;
            } else if (filter_data.record[i] == 1) {
                filter_data.edge_project[medium_edge_ptr++] = i;
            } else if (filter_data.record[i] == 2) {
                filter_data.edge_project[heavy_edge_ptr++] = i;
            }
        }
    } else {
        for (e_index_t i = 0; i < g->e_cnt; ++i) {
            v_index_t vtx_left = g->edge_from[i];
            v_index_t vtx_right = g->edge[i];
            if (filter_data.degree[vtx_left] < sparse_threshold && filter_data.degree[vtx_right] < sparse_threshold) {
                filter_data.record[i] = 0;
                num_of_light_edges++;
            } else if (filter_data.degree[vtx_left] < dense_threshold && filter_data.degree[vtx_right] < dense_threshold) {
                filter_data.record[i] = 1;
                num_of_medium_edges++;
            } else {
                filter_data.record[i] = 2;
                num_of_heavy_edges++;
            }
        }
        medium_edge_ptr = num_of_light_edges;
        heavy_edge_ptr = num_of_light_edges + num_of_medium_edges;
        // 投射记录
        for (e_index_t i = 0; i < g->e_cnt; ++i) {
            if (filter_data.record[i] == 0) {
                filter_data.edge_project[light_edge_ptr++] = i;
            } else if (filter_data.record[i] == 1) {
                filter_data.edge_project[medium_edge_ptr++] = i;
            } else if (filter_data.record[i] == 2) {
                filter_data.edge_project[heavy_edge_ptr++] = i;
            }
        }
    }

    // 用于模式数量估计的数据结构，每个节点存储一份。被过滤的数据（精确）需要唯一计算，每个节点计算一部分。被保留的数据（近似）可以多次计算，每个节点保留自己的版本。
    // 与原始图对应
    filter_data.sample_count_per_edge = (double *)malloc(sizeof(double) * g->e_cnt);
    filter_data.sample_times_per_edge = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
    filter_data.sample_square_per_edge = (double *)malloc(sizeof(double) * g->e_cnt);
    memset(filter_data.sample_count_per_edge, 0, sizeof(double) * g->e_cnt);
    memset(filter_data.sample_times_per_edge, 0, sizeof(uint32_t) * g->e_cnt);
    memset(filter_data.sample_square_per_edge, 0, sizeof(double) * g->e_cnt);

    // 分布式执行
    // 工作流程
    // 1. 每个节点根据rank处理一个范围内的light edge。
    // 2. 处理完lightedge后，将本地信息发送到主节点。然后按照batch处理所有的medium和heavy edge，异步的发送统计信息。
    // 3. 如果处理的batch数量超过某个值，还没有受到收敛的信息，处理meidium edge，然后统计记录的heavy数据，凑个整数发送heavy edge的统计信息。
    // *主节点需要构建多个存储结构分别记录来自不同节点的信息，以便于后续的更新。
    // 4. 可以进一步优化，将heavy edge继续细分，最后实现无偏估计

    // 采样批量
    uint64_t batch_size = 100;

    // 每个节点判断是否满足要求的标志，每个线程都会使用到
    bool error_bound_satisfied = false;

    // 主线程判断是否满足要求使用到的变量，只有主线程会使用到
    double var = 0.0;
    double std_dev = 0.0;
    double real_error = 0.0;
    double error_bound = 0.01;

    // 计算节点相关信息
    int comm_sz, my_rank, num_processes;

    // 计算节点处理的边表范围
    uint32_t edge_range_start;
    uint32_t edge_range_end;
    uint32_t edge_range_size_per_node;
    uint32_t edge_range_size_per_thread;
    double local_node_count = 0;                 // 本地精确计数
    uint32_t local_node_count_compute_times = 0; // 计算次数

    // 缓存相关
    MPI_Status status;

    // 主节点接收 近似估计数据    数据类型无差别，可以异步执行。
    double recv_sample_count_and_squra[BUFFER_SIZE] = {0};
    MPI_Request recv_sample_rqst[BUFFER_SIZE];
    uint32_t approx_headptr = 0, approx_tailptr = 0;
    double approx_count = 0;
    double approx_squra = 0;
    uint64_t approx_sample_times = 0;

    // 主节点接收 精确计算数据
    double recv_exact_count[BUFFER_SIZE] = {0};
    MPI_Request recv_exact_rqst[BUFFER_SIZE];
    uint32_t exact_headptr = 0, exact_tailptr = 0;
    double exact_count = 0;
    uint64_t exact_count_received_times = 0;

    // 其他节点接收计算结果
    MPI_Request recv_errror_bound_rqst;

#pragma omp parallel num_threads(2)
    {
#pragma omp barrier
#pragma omp master
        {
            // 每个进程的主线程执行 初始化操作
            int provided;
            MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &provided);
            MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
            MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
            MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
            MPI_Barrier(MPI_COMM_WORLD);

            // 根据节点id划分边表范围
            edge_range_size_per_node = num_of_light_edges / num_processes;
            edge_range_size_per_thread = edge_range_size_per_node / (omp_get_max_threads() - 1);
            edge_range_start = my_rank * edge_range_size_per_node;
            if (my_rank == num_processes - 1)
                edge_range_end = num_of_light_edges;
            else
                edge_range_end = (my_rank + 1) * edge_range_size_per_node;
        }
#pragma omp barrier

#pragma omp master
        {
            while (!error_bound_satisfied) {
                if (!my_rank) {
                    /*接收精确计算的结果*/
                    // TODO 改为满足要求后就不执行了
                    // 接收来自其他节点的消息
                    if (exact_tailptr - exact_headptr < num_processes) {
                        MPI_Irecv(&recv_exact_count[exact_tailptr % BUFFER_SIZE], 1, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD,
                                  &recv_exact_rqst[exact_tailptr % BUFFER_SIZE]);
                        exact_tailptr = exact_tailptr + 1;
                    }
                    if (exact_headptr < exact_tailptr) {
                        int testflag = 0;
                        MPI_Test(&recv_exact_rqst[exact_headptr % BUFFER_SIZE], &testflag, MPI_STATUS_IGNORE);
                        if (testflag) {
                            exact_count += recv_exact_count[exact_headptr % BUFFER_SIZE];
                            exact_headptr = exact_headptr + 1;
                        }
                    }
                    // 接收来自本地的消息
                    if (local_node_count_compute_times == omp_get_max_threads() - 1) {
                        exact_count += local_node_count;
                        local_node_count_compute_times = 0;
                    }

                    /*接收近似计算的结果*/
                    // 当缓存中有数据的时候接收数据
                    if (approx_tailptr - approx_headptr < BUFFER_SIZE) {
                        MPI_Irecv(&recv_sample_count_and_squra[approx_tailptr % BUFFER_SIZE], 2, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,
                                  &recv_sample_rqst[approx_tailptr % BUFFER_SIZE]);
                        approx_tailptr = approx_tailptr + 2;
                    }
                    if (approx_headptr < approx_tailptr) // 如果缓存中还有数据 ，将缓存中的数据接收并且更新整体
                    {
                        int testflag = 0;
                        MPI_Test(&recv_sample_rqst[approx_headptr % BUFFER_SIZE], &testflag, MPI_STATUS_IGNORE);
                        if (testflag) {
                            // #pragma omp atomic
                            approx_count += recv_sample_count_and_squra[approx_tailptr % BUFFER_SIZE];
                            approx_squra += recv_sample_count_and_squra[approx_tailptr % BUFFER_SIZE + 1];
                            // #pragma omp atomic
                            approx_sample_times += batch_size;
                            approx_headptr = approx_headptr + 2;

                            // 计算全体方差和每个值得标准差
                            var = (double)(approx_squra / approx_sample_times - std::pow((approx_count / approx_sample_times), 2));
                            std_dev = std::sqrt(var / approx_sample_times);

                            real_error = std_dev * 2.576 / (approx_count / approx_sample_times); // 99%的信心
                            if (real_error < error_bound) {
                                error_bound_satisfied = true;
                            }
                        }
                        // 满足误差要求之后，更新各个节点的信息
                        if (error_bound_satisfied) {
                            for (uint32_t i = 1; i < num_processes; i++)
                                MPI_Send(&error_bound_satisfied, 1, MPI_C_BOOL, i, 0, MPI_COMM_WORLD);
                            break;
                        }
                    }
                }
                if (my_rank) // 非主进程的主线程的工作：做一些局部的统计工作，算ok了传递给主进程。
                {
                    // 如果局部精确已经计算完了，发送给主进程
                    if (local_node_count_compute_times == omp_get_max_threads() - 1) {
                        MPI_Send(&local_node_count, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
                        local_node_count_compute_times = 0;
                    }

                    // 如果不满足退出条件，尝试接收一下消息
                    if (!error_bound_satisfied) {
                        // 没有发送消息过来，会一直阻塞。
                        int error_bound_flag = 0;
                        MPI_Irecv(&error_bound_satisfied, 1, MPI_C_BOOL, 0, 0, MPI_COMM_WORLD, &recv_errror_bound_rqst);
                        MPI_Test(&recv_errror_bound_rqst, &error_bound_flag, MPI_STATUS_IGNORE);
                        if (error_bound_flag == 0) {
                            MPI_Cancel(&recv_errror_bound_rqst);
                            MPI_Wait(&recv_errror_bound_rqst, &status);
                        }
                    }
                }
            }
        }

        // 每个进程的非主线程执行的任务，计算完更新主节点的主线程
        if (omp_get_thread_num()) {
            // 变量初始化
            double count = 0;
            bool exact_computing = false;
            int computing_type = 0;
            uint32_t num_of_batch = 0;

            // 发送缓存相关
            MPI_Request send_sample_rqst[BUFFER_SIZE]; // 用于标记某个发送请求是否已经发送完毕
            double local_sample_count_and_squra[BUFFER_SIZE] = {0};
            uint32_t headptr = 0, tailptr = 0;

            // 随机数生成器
            std::random_device rd;
            std::default_random_engine gen(rd());
            // TODO:可能需要修改范围
            std::uniform_int_distribution<uint64_t> dist(0, num_of_heavy_edges + num_of_medium_edges - 1);

            // 执行，直到满足计算要求
            while (!error_bound_satisfied) {
                /*精确计算*/
                // 数据更新到本节点的主线程
                if (exact_computing == false) {
                    if (computing_type == 0) {
                        // 计算本地处理的light edge的范围，并进行计算计算。
                        uint32_t local_range_start = edge_range_start + (omp_get_thread_num() - 1) * edge_range_size_per_thread;
                        uint32_t local_range_end;
                        if (omp_get_thread_num() != omp_get_max_threads())
                            local_range_end = edge_range_start + omp_get_thread_num() * edge_range_size_per_thread;
                        else
                            local_range_end = edge_range_end;

                        for (int i = local_range_start; i < local_range_end; i++) {
                            uint32_t edge_id = filter_data.edge_project[i];
                            pattern_sample_es(edge_id, g, s, count);
                        }
#pragma omp atomic
                        local_node_count += count;

#pragma omp atomic
                        local_node_count_compute_times = local_node_count_compute_times + 1;

                        exact_computing = true;
                    } else if (computing_type == 1) {
                        exact_computing = true;
                    }
                }

                /*近似计算*/
                // 异步计算。在发送缓存中还有空间的时候进行计算并且发送数据。
                // tailptr - headptr是可能超出BUFFER_SIZE的
                if (tailptr - headptr < BUFFER_SIZE) {

                    // 按照固定的batch进行测量
                    double count_sample = 0;
                    double count_sample_square = 0;
                    for (uint64_t i = 0; i < batch_size; i++) {
                        auto eid = filter_data.edge_project[dist(gen) + num_of_light_edges];
                        pattern_sample_ns(eid, g, s, count_sample, count_sample_square);
                    }
                    num_of_batch++;
                    if (num_of_batch > 10000) {
                        exact_computing = false;
                        computing_type = 1;
                    }

                    local_sample_count_and_squra[tailptr % BUFFER_SIZE] = count_sample;
                    local_sample_count_and_squra[tailptr % BUFFER_SIZE + 1] = count_sample_square;

                    MPI_Isend(&local_sample_count_and_squra[tailptr % BUFFER_SIZE], 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
                              &send_sample_rqst[tailptr % BUFFER_SIZE]);
                    tailptr = tailptr + 2;
                }

                // 检测是否发送完成，更新发送缓存
                if (headptr < tailptr) {
                    int testflag = 0;
                    MPI_Test(&send_sample_rqst[headptr % BUFFER_SIZE], &testflag, MPI_STATUS_IGNORE); // 异步的操作需要一个test来执行
                    if (testflag)
                        headptr = headptr + 2;
                }
            }

            // 输出一下满足计算要求后，本地计算消耗的资源
            std::cout << "work thread " << "headptr " << headptr << " tailptr " << tailptr << std::endl;
            // for (int i = 0; i < 100; i++) {
            //     if (local_prob_inverse[i] != recv_prob_inverse[i]) {
            //         printf("wrong");
            //     }
            // }
        }
    }

    // 线程并行执行区域结束，这里只有主进程主线程会执行
    if (!my_rank) {
        // for (int i = 0; i < 100; i++) {
        //     if (recv_prob_inverse[i] == 0) {
        //         printf("wrong");
        //     }
        // }
        std::cout << "global summary global_sample_times: " << approx_count << " " << approx_sample_times << std::endl;
        double estimated_pattern_count = approx_count * (num_of_heavy_edges + num_of_medium_edges) / (approx_sample_times);
        std::cout << "pattern_nums " << estimated_pattern_count + exact_count << std::endl;
        double estimated_pattern_count_error = (estimated_pattern_count + exact_count - (26184567090560.0)) / (26184567090560.0);
        std::cout << "estimated_error " << estimated_pattern_count_error << std::endl;
    }

    // lj图 house 26184567090560 单线程3000秒，64线程几十秒 scalegpm 1.1秒
    MPI_Finalize();

    free(filter_data.edge_project);
    free(filter_data.record);
    free(filter_data.degree);
    free(filter_data.sample_count_per_edge);
    free(filter_data.sample_times_per_edge);
    free(filter_data.sample_square_per_edge);
}
void remove_anti_edge_vertices(VertexSet &out_buf, const VertexSet &in_buf, const Schedule_IEP &sched, const VertexSet &partial_embedding, int vp,
                               Graph *g) {

    out_buf.init();
    assert(&out_buf != &in_buf);
    auto d_out = out_buf.get_data_ptr();
    assert(d_out != nullptr);
    auto d_in = in_buf.get_data_ptr();
    // printf("d_out:%x d_in:%x\n", d_out, d_in);
    assert(d_out != d_in);
    int n_in = in_buf.get_size();
    int out_size = 0;

    for (int i = 0; i < n_in; i++) {
        bool produce_output = true;
        for (int u = 0; u < vp; ++u) {
            if (sched.get_adj_mat_ptr()[u * sched.get_size() + vp] != 0)
                continue;

            auto v = partial_embedding.get_data(u); // 获取位置为u的顶点
            e_index_t l, r;
            l = g->vertex[v];
            r = g->vertex[v + 1];
            int m = r - l; // m = |N(v)|

            if (binary_search(&g->edge[l], m, d_in[i])) {
                produce_output = false;
                break;
            }
        }
        if (produce_output) {
            d_out[out_size++] = d_in[i];
        }
    }
    out_buf.set_size(out_size);
}

void pattern_matching_aggressive_func(const Schedule_IEP &schedule, VertexSet *vertex_set, VertexSet &subtraction_set, VertexSet &tmp_set,
                                      long long &local_ans, int depth, int *ans_buffer, Graph *g) {

    // 获取需要处理的模式id
    int loop_set_prefix_id = schedule.get_loop_set_prefix_id(depth); // 从第一个模式顶点开始处理

    // 获取当前处理的模式顶点的候选集
    auto vset = &vertex_set[loop_set_prefix_id];
    int loop_size = (*vset).get_size();
    if (loop_size <= 0)
        return;

    // 处理当前顶点集合中的所有顶点
    int *loop_data_ptr = (*vset).get_data_ptr();
    // 如果是顶点诱导模式，可以进行优化
    if (schedule.is_vertex_induced) {
        VertexSet &diff_buf = vertex_set[schedule.get_total_prefix_num() + depth]; // 获取接近末尾的模式定对的候选集？？？？？？？？
        // printf("depth:%d loop_set_prefix_id = %d diff_buf: %d\n",depth,
        // loop_set_prefix_id, schedule.get_total_prefix_num() + depth);
        remove_anti_edge_vertices(diff_buf, vertex_set[loop_set_prefix_id], schedule, subtraction_set, depth,
                                  g); // 根据该候选集对于当前顶点进行移除？？？？？？？？？？？

        loop_data_ptr = diff_buf.get_data_ptr(); // 更新当前处理顶点的候选集
        loop_size = diff_buf.get_size();
        vset = &diff_buf;
    }
    // Case: in_exclusion_optimize_num > 1

    // 如果iep优化数大于1，则进行iep计算并且返回

    if (depth == schedule.get_size() - schedule.get_in_exclusion_optimize_num()) { // 存在iep优化潜力的情况下，处理到了最后一个顶点
        int last_pos = -1;
        long long val;

        for (int i = 0; i < schedule.in_exclusion_optimize_vertex_id.size(); ++i) {
            if (schedule.in_exclusion_optimize_vertex_flag[i]) {
                ans_buffer[i] = vertex_set[schedule.in_exclusion_optimize_vertex_id[i]].get_size() - schedule.in_exclusion_optimize_vertex_coef[i];
            } else {
                ans_buffer[i] = VertexSet::unordered_subtraction_size(vertex_set[schedule.in_exclusion_optimize_vertex_id[i]], subtraction_set);
            }
        }
        for (int pos = 0; pos < schedule.in_exclusion_optimize_coef.size(); ++pos) {
            if (pos == last_pos + 1)
                val = ans_buffer[schedule.in_exclusion_optimize_ans_pos[pos]];
            else {
                if (val != 0)
                    val = val * ans_buffer[schedule.in_exclusion_optimize_ans_pos[pos]];
            }
            if (schedule.in_exclusion_optimize_flag[pos]) {
                last_pos = pos;
                local_ans += val * schedule.in_exclusion_optimize_coef[pos];
            }
        }
        return;
    }
    // 如果iep优化数小于等于1
    if (depth == schedule.get_size() - 1) { // 处理到了最后一个顶点，不存在iep优化潜力的情况下
        // TODO : try more kinds of calculation.
        // For example, we can maintain an ordered set, but it will cost more to
        // maintain itself when entering or exiting recursion.
        if (schedule.get_total_restrict_num() > 0) {
            int min_vertex = g->v_cnt;
            for (int i = schedule.get_restrict_last(depth); i != -1; i = schedule.get_restrict_next(i))
                if (min_vertex > subtraction_set.get_data(schedule.get_restrict_index(i)))
                    min_vertex = subtraction_set.get_data(schedule.get_restrict_index(i));
            int size_after_restrict =
                std::lower_bound((*vset).get_data_ptr(), (*vset).get_data_ptr() + (*vset).get_size(), min_vertex) - (*vset).get_data_ptr();
            if (size_after_restrict > 0)
                local_ans += VertexSet::unordered_subtraction_size((*vset), subtraction_set, size_after_restrict);
        } else
            local_ans += VertexSet::unordered_subtraction_size((*vset), subtraction_set);
        return;
    }

    /*    #pragma omp critical
        {
            if( depth == 1) ++dep1_cnt;
            if( depth == 2) ++dep2_cnt;
            if( depth == 3) ++dep3_cnt;
        }*/
    // TODO : min_vertex is also a loop invariant

    // 前面考虑了处理到最后一个顶点的情况（有iep或者没有iep的情况下）。
    // ***这里考虑处理靠前的顶点的情况***//

    // 处理当前深度的所有限制，查看限制当前模式顶点的模式顶点匹配顶点最小的是哪一个
    int min_vertex = g->v_cnt;
    for (int i = schedule.get_restrict_last(depth); i != -1; i = schedule.get_restrict_next(i)) // 处理当前深度的所有限制
        if (min_vertex > subtraction_set.get_data(schedule.get_restrict_index(i)))
            min_vertex = subtraction_set.get_data(schedule.get_restrict_index(i));

    // 处理
    for (int i = 0; i < loop_size; ++i) {
        if (min_vertex <= loop_data_ptr[i]) // 如果限制当前模式顶点的模式顶点的匹配顶点小于当前顶点，退出循环， 打破对称性的一种方法。
            break;

        // 对于当前候选集中的顶点依次进行匹配，如果已经存在于差集中，取下一个
        int vertex = loop_data_ptr[i];
        if (subtraction_set.has_data(vertex))
            continue;
        e_index_t l, r;
        l = g->vertex[vertex];
        r = g->vertex[vertex + 1];

        // 根据当前匹配顶点，为后续每一个模式顶点重新构建候选集？  在构建边表的时候已经为后续每一个匹配施加限制了。
        bool is_zero = false;
        for (int prefix_id = schedule.get_last(depth); prefix_id != -1; prefix_id = schedule.get_next(prefix_id)) {
            vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &g->edge[l], (int)(r - l), prefix_id, vertex);
            // if( vertex_set[prefix_id].get_size() == 0 && prefix_id <
            // schedule.get_basic_prefix_num()) {
            if (vertex_set[prefix_id].get_size() == schedule.break_size[prefix_id]) {
                is_zero = true;
                break;
            }
        }
        if (is_zero)
            continue;
        // subtraction_set.insert_ans_sort(vertex);
        subtraction_set.push_back(vertex); // 将当前匹配的顶点插入差集
        pattern_matching_aggressive_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, depth + 1, ans_buffer, g); // 继续进行深度搜索
        subtraction_set.pop_back(); // 深度搜索处理完之后，弹出当前匹配顶点
    }
}

// void vertex_vs_edge_test(Graph *g, schedule_approx *s, Pattern p) {

//     // 记录需要处理的顶点
//     uint32_t *vtx_record = (uint32_t *)malloc(sizeof(uint32_t) * g->v_cnt);
//     memset(vtx_record, 0, sizeof(uint32_t) * g->v_cnt);

// edge_from_init(g);
//     /*过滤操作*/
//     // 初始化阈值
//     int sparse_threshold = g->e_cnt / g->v_cnt;    // 平均度数  frendster 28 twitter 35 livejouranl 17 youtube 2 patent 5 pokec 30 tp 28 sk 11
//     int dense_threshold = threshhold_computing(g); // 百分之九十的顶点的度数阈值
//     if (dense_threshold < sparse_threshold) {
//         int temp_threshold = dense_threshold;
//         dense_threshold = sparse_threshold;
//         sparse_threshold = temp_threshold;
//     }
//     sparse_threshold = 8;

//     // 过滤的数据结构。
//     filter filter_data;
//     filter_data.edge_project = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
//     filter_data.record = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
//     memset(filter_data.edge_project, -1, sizeof(uint32_t) * g->e_cnt);
//     memset(filter_data.record, -1, sizeof(uint32_t) * g->e_cnt);
//     filter_data.degree = (uint32_t *)malloc(sizeof(uint32_t) * g->v_cnt);
//     for (int i = 0; i < g->v_cnt; i++) {
//         filter_data.degree[i] = g->vertex[i + 1] - g->vertex[i];
//     }

//     uint64_t num_of_light_edges = 0;
//     uint64_t num_of_medium_edges = 0;
//     uint64_t num_of_heavy_edges = 0;
//     filter_data.light_edge_ptr = 0;
//     filter_data.medium_edge_ptr = 0;
//     filter_data.heavy_edge_ptr = 0;
//     for (e_index_t i = 0; i < g->e_cnt; ++i) {
//         v_index_t vtx_left = g->edge_from[i];
//         v_index_t vtx_right = g->edge[i];
//         if (filter_data.degree[vtx_left] < sparse_threshold || filter_data.degree[vtx_right] < sparse_threshold) {
//             filter_data.record[i] = 0;
//             num_of_light_edges++;
//             if (filter_data.degree[vtx_left] < sparse_threshold && filter_data.degree[vtx_right] < sparse_threshold) {
//                 vtx_record[vtx_left] = 1;
//                 vtx_record[vtx_right] = 1;
//             } else if (filter_data.degree[vtx_left] < sparse_threshold) {
//                 vtx_record[vtx_left] = 1;
//             } else if (filter_data.degree[vtx_right] < sparse_threshold) {
//                 vtx_record[vtx_right] = 1;
//             }
//         } else if (filter_data.degree[vtx_left] < dense_threshold || filter_data.degree[vtx_right] < dense_threshold) {
//             filter_data.record[i] = 1;
//             num_of_medium_edges++;
//             if (filter_data.degree[vtx_left] < dense_threshold && filter_data.degree[vtx_right] < dense_threshold) {
//                 vtx_record[vtx_left] = 2;
//                 vtx_record[vtx_right] = 2;
//             } else if (filter_data.degree[vtx_left] < dense_threshold) {
//                 vtx_record[vtx_left] = 2;
//             } else if (filter_data.degree[vtx_right] < dense_threshold) {
//                 vtx_record[vtx_right] = 2;
//             }
//         } else {
//             filter_data.record[i] = 2;
//             num_of_heavy_edges++;
//         }
//     }
//     filter_data.medium_edge_ptr = num_of_light_edges;
//     filter_data.heavy_edge_ptr = num_of_light_edges + num_of_medium_edges;
//     // 投射记录
//     for (e_index_t i = 0; i < g->e_cnt; ++i) {
//         if (filter_data.record[i] == 0) {
//             filter_data.edge_project[filter_data.light_edge_ptr++] = i;
//         } else if (filter_data.record[i] == 1) {
//             filter_data.edge_project[filter_data.medium_edge_ptr++] = i;
//         } else if (filter_data.record[i] == 2) {
//             filter_data.edge_project[filter_data.heavy_edge_ptr++] = i;
//         }
//     }

//     // 过滤边对于es的消耗
//     TimeInterval allTime;
//     allTime.check();
//     double count = 0;
//     for (int i = 0; i < num_of_light_edges; i++) {
//         uint32_t edge_id = filter_data.edge_project[i];
//         pattern_sample_es(edge_id, g, s, count);
//     }
//     allTime.print("es time cost");

//     bool is_pattern_valid;
//     bool use_in_exclusion_optimize = true;
//     Schedule_IEP schedule_iep(p, is_pattern_valid, 1, 1, use_in_exclusion_optimize, g->v_cnt, g->e_cnt, g->tri_cnt);
//     int *ans_buffer = new int[schedule_iep.in_exclusion_optimize_vertex_id.size()];  // 底层无交集顶点的大小
//     VertexSet *vertex_set = new VertexSet[schedule_iep.get_total_prefix_num() + 10]; // 构建了多个顶点集合作为候选集
//     VertexSet subtraction_set;                                                       // 构建一个差集
//     VertexSet tmp_set;                                                               // 构建一个临时集合
//     subtraction_set.init();
//     long long local_ans = 0;
//     // TODO : try different chunksize

//     // 对于每一个顶点进行处理

//     TimeInterval vtxTime;
//     vtxTime.check();

//     for (int vertex = 0; vertex < g->v_cnt; ++vertex) {
//         if (vtx_record[vertex] == 1) {
//             e_index_t l, r;
//             l = g->vertex[vertex];
//             r = g->vertex[vertex + 1];

//             // 为调度中的每一个顶点初始化候选顶点集合。
//             for (int prefix_id = schedule_iep.get_last(0); prefix_id != -1; prefix_id = schedule_iep.get_next(prefix_id)) {
//                 vertex_set[prefix_id].build_vertex_set(schedule_iep, vertex_set, &g->edge[l], (int)(r - l), prefix_id);
//             }
//             // subtraction_set.insert_ans_sort(vertex);
//             subtraction_set.push_back(vertex); // 将当前顶点放入差集中，也就是后续的所有顶点都不能够是该顶点
//                                                // if (schedule.get_total_restrict_num() > 0 && clique == false)

//             // 进行匹配操作

//             pattern_matching_aggressive_func(schedule_iep, vertex_set, subtraction_set, tmp_set, local_ans, 1, ans_buffer, g);

//             subtraction_set.pop_back();
//         }
//         // double end_time = get_wall_time();
//         // printf("my thread time %d %.6lf\n", omp_get_thread_num(), end_time -
//         // start_time);

//         // TODO : Computing multiplicty for a pattern
//     }
//     delete[] vertex_set;
//     vtxTime.print("vt time cost");
// }
// 通过rank来确定进程，通过tag来确认消息

// 工作流程
// 先采样，后过滤。
// 每个节点先记录 分层的边集和顶点集
// 第一步过滤不需要的边集

// 异步的进行采样流程，如果长时间不收敛，进行过滤

// 过滤操作，主线程重新扫一遍所有的数据然后，然后发送给主节点。后续的更新不扫了。

// 消息的tag

/*关键，要使用与从顶点出发顺序一样的调度顺序*/
void pattern_count_ns_mpi_fix_error_filter_vertex_after_sample(Graph *g, schedule_approx *s, const Schedule_IEP &schedule_iep, Pattern p) {
    // 获取边的起点

    edge_from_init(g);

    /*过滤操作*/
    // 初始化阈值
    int sparse_threshold = g->e_cnt / g->v_cnt;    // 平均度数  frendster 28 twitter 35 livejouranl 17 youtube 2 patent 5 pokec 30 tp 28 sk 11
    int dense_threshold = threshhold_computing(g); // 百分之九十的顶点的度数阈值
    if (dense_threshold < sparse_threshold) {
        int temp_threshold = dense_threshold;
        dense_threshold = sparse_threshold;
        sparse_threshold = temp_threshold;
    }

    // 过滤的数据结构。
    filter filter_data;
    filter_data.edge_project = (int64_t *)malloc(sizeof(int64_t) * g->e_cnt);
    filter_data.edge_record = (int *)malloc(sizeof(int) * g->e_cnt);
    filter_data.vertex_project = (int *)malloc(sizeof(int) * g->v_cnt);
    filter_data.vertex_record = (int *)malloc(sizeof(int) * g->v_cnt);
    memset(filter_data.edge_project, -1, sizeof(int64_t) * g->e_cnt);
    memset(filter_data.edge_record, -1, sizeof(int) * g->e_cnt);
    memset(filter_data.vertex_record, -1, sizeof(int) * g->v_cnt);
    memset(filter_data.vertex_project, -1, sizeof(int) * g->v_cnt);

    filter_data.degree = (uint32_t *)malloc(sizeof(uint32_t) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        filter_data.degree[i] = g->vertex[i + 1] - g->vertex[i];
    }

    uint64_t num_of_light_edges = 0;
    uint64_t num_of_medium_edges = 0;
    uint64_t num_of_heavy_edges = 0;
    uint64_t num_of_light_vertex = 0;
    uint64_t num_of_medium_vertex = 0;
    uint64_t num_of_heavy_vertex = 0;

    uint64_t light_edge_ptr = 0;
    uint64_t medium_edge_ptr = 0;
    uint64_t heavy_edge_ptr = 0;
    uint64_t light_vertex_ptr = 0;
    uint64_t medium_vertex_ptr = 0;
    uint64_t heavy_vertex_ptr = 0;

    for (uint64_t i = 0; i < g->v_cnt; ++i) {
        if (filter_data.degree[i] > s->order_degree[0]) {
            if (filter_data.degree[i] < sparse_threshold) {
                filter_data.vertex_record[i] = 0;
                num_of_light_vertex++;
                for (uint64_t j = g->vertex[i]; j < g->vertex[i + 1]; ++j) {
                    filter_data.edge_record[j] = 0;
                    num_of_light_edges++;
                }
            } else if (filter_data.degree[i] < dense_threshold) {
                filter_data.vertex_record[i] = 1;
                num_of_medium_vertex++;
                for (uint64_t j = g->vertex[i]; j < g->vertex[i + 1]; ++j) {
                    filter_data.edge_record[j] = 1;
                    num_of_medium_edges++;
                }
            } else {
                filter_data.vertex_record[i] = 2;
                num_of_heavy_vertex++;
                for (uint64_t j = g->vertex[i]; j < g->vertex[i + 1]; ++j) {
                    filter_data.edge_record[j] = 2;
                    num_of_heavy_edges++;
                }
            }
        }
    }
    // 投射记录
    medium_vertex_ptr = num_of_light_vertex;
    heavy_vertex_ptr = num_of_light_vertex + num_of_medium_vertex;
    for (uint32_t i = 0; i < g->v_cnt; ++i) {
        if (filter_data.vertex_record[i] == 0) {
            filter_data.vertex_project[light_vertex_ptr++] = i;
        } else if (filter_data.vertex_record[i] == 1) {
            filter_data.vertex_project[medium_vertex_ptr++] = i;
        } else if (filter_data.vertex_record[i] == 2) {
            filter_data.vertex_project[heavy_vertex_ptr++] = i;
        }
    }
    medium_edge_ptr = num_of_light_edges;
    heavy_edge_ptr = num_of_light_edges + num_of_medium_edges;
    for (e_index_t i = 0; i < g->e_cnt; ++i) {
        if (filter_data.edge_record[i] == 0) {
            filter_data.edge_project[light_edge_ptr++] = i;
        } else if (filter_data.edge_record[i] == 1) {
            filter_data.edge_project[medium_edge_ptr++] = i;
        } else if (filter_data.edge_record[i] == 2) {
            filter_data.edge_project[heavy_edge_ptr++] = i;
        }
    }

    // 用于模式数量估计的数据结构，每个节点存储一份。被过滤的数据（精确）需要唯一计算，每个节点计算一部分。被保留的数据（近似）可以多次计算，每个节点保留自己的版本。
    // 与原始图对应
    filter_data.sample_times_per_edge = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
    filter_data.sample_count_per_edge = (double *)malloc(sizeof(double) * g->e_cnt);
    filter_data.sample_square_per_edge = (double *)malloc(sizeof(double) * g->e_cnt);
    memset(filter_data.sample_times_per_edge, 0, sizeof(uint32_t) * g->e_cnt);
    memset(filter_data.sample_count_per_edge, 0, sizeof(double) * g->e_cnt);
    memset(filter_data.sample_square_per_edge, 0, sizeof(double) * g->e_cnt);

    // 分布式执行
    // 工作流程
    // 1. 每个节点根据rank处理一个范围内的light edge。
    // 2. 处理完lightedge后，将本地信息发送到主节点。然后按照batch处理所有的medium和heavy edge，异步的发送统计信息。
    // 3. 如果处理的batch数量超过某个值，还没有受到收敛的信息，处理meidium edge，然后统计记录的heavy数据，凑个整数发送heavy edge的统计信息。
    // *主节点需要构建多个存储结构分别记录来自不同节点的信息，以便于后续的更新。
    // 4. 可以进一步优化，将heavy edge继续细分，最后实现无偏估计

    // 采样批量
    uint64_t batch_size = 10000;

    // 每个节点判断是否满足要求的标志，每个线程都会使用到
    bool error_bound_satisfied = false;

    // 主线程判断是否满足要求使用到的变量，只有主线程会使用到
    double var = 0.0;
    double std_dev = 0.0;
    double real_error = 0.0;
    double error_bound = 0.01;

    // 计算节点相关信息
    int comm_sz, my_rank, num_processes;

    // 缓存相关
    MPI_Status status;

    // 主节点接收 近似估计数据    数据类型无差别，可以异步执行。
    double recv_sample_count_and_squra[BUFFER_SIZE] = {0};
    MPI_Request recv_sample_rqst[BUFFER_SIZE];
    uint32_t approx_headptr = 0, approx_tailptr = 0;
    double approx_count = 0;
    double approx_squra = 0;
    uint64_t approx_sample_times = 0;

    // 主节点接收 精确计算数据
    double recv_exact_count[BUFFER_SIZE] = {0};
    MPI_Request recv_exact_rqst[BUFFER_SIZE];
    uint32_t exact_headptr = 0, exact_tailptr = 0;
    double exact_count = 0;
    uint64_t exact_count_received_times = 0;

    // 精确计算相关  //分为三段
    bool exact_computing = false;
    bool local_exact_computing[64] = {false};
    int computing_type = 0;

    uint32_t vertex_range_size_per_node[3];
    uint32_t vertex_range_start_per_node[3];
    uint32_t vertex_range_end_per_node[3];
    uint32_t vertex_chunk_per_task = 1000;
    uint32_t vertex_range_ptr[3];

    double local_node_exact_count = 0;                 // 本地精确计数
    uint32_t local_node_exact_count_compute_times = 0; // 记录有多少个线程完成了精确计数

    uint64_t local_batch_sample_times = 0;

    // 其他节点 接收 误差要求，精确计算要求
    MPI_Request recv_errror_bound_rqst;
    MPI_Request recv_exact_comput_rqst;

#pragma omp parallel num_threads(8)
    {
#pragma omp barrier
#pragma omp master
        {
            // 每个进程的主线程执行 初始化操作
            int provided;
            MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &provided);
            MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
            MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
            MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
            MPI_Barrier(MPI_COMM_WORLD);

            // 根据节点id划分顶点处理范围
            vertex_range_size_per_node[0] = num_of_light_vertex / num_processes;
            vertex_range_start_per_node[0] = my_rank * vertex_range_size_per_node[0];
            vertex_range_ptr[0] = vertex_range_start_per_node[0];
            if (my_rank == num_processes - 1)
                vertex_range_end_per_node[0] = num_of_light_vertex;
            else
                vertex_range_end_per_node[0] = (my_rank + 1) * vertex_range_size_per_node[0];

            vertex_range_size_per_node[1] = num_of_medium_vertex / num_processes;
            vertex_range_start_per_node[1] = my_rank * vertex_range_size_per_node[1] + num_of_light_vertex;
            vertex_range_ptr[1] = vertex_range_start_per_node[1];
            if (my_rank == num_processes - 1)
                vertex_range_end_per_node[1] = num_of_medium_vertex + num_of_light_vertex;
            else
                vertex_range_end_per_node[1] = (my_rank + 1) * vertex_range_size_per_node[1] + num_of_light_vertex;

            vertex_range_size_per_node[2] = num_of_heavy_vertex / num_processes;
            vertex_range_start_per_node[2] = my_rank * vertex_range_size_per_node[2] + num_of_medium_vertex + num_of_light_vertex;
            vertex_range_ptr[2] = vertex_range_start_per_node[2];
            if (my_rank == num_processes - 1)
                vertex_range_end_per_node[2] = num_of_heavy_vertex + num_of_medium_vertex + num_of_light_vertex;
            else
                vertex_range_end_per_node[2] = (my_rank + 1) * vertex_range_size_per_node[2] + num_of_medium_vertex + num_of_light_vertex;
        }
#pragma omp barrier

#pragma omp master
        {
            while (!error_bound_satisfied) {
                // 两个发送，两个接收任务   发送精确计算要求 3号tag ，误差满足要求 0号tag，接收精确计算结果 2，近似计算结果 1
                if (!my_rank) {

                    /*接收 精确计算结果*/
                    if (exact_tailptr - exact_headptr < num_processes) {
                        MPI_Irecv(&recv_exact_count[exact_tailptr % BUFFER_SIZE], 1, MPI_DOUBLE, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD,
                                  &recv_exact_rqst[exact_tailptr % BUFFER_SIZE]);
                        exact_tailptr = exact_tailptr + 1;
                    }
                    if (exact_headptr < exact_tailptr) {
                        int testflag = 0;
                        MPI_Test(&recv_exact_rqst[exact_headptr % BUFFER_SIZE], &testflag, MPI_STATUS_IGNORE);
                        if (testflag) {
                            exact_count += recv_exact_count[exact_headptr % BUFFER_SIZE];
                            exact_headptr = exact_headptr + 1;
                        }
                    }
                    // 接收来自本地的消息
                    if (local_node_exact_count_compute_times == omp_get_max_threads() - 1) {
                        exact_count += local_node_exact_count;
                        local_node_exact_count_compute_times = 0;

                        // TODO::修改为全局的精确计算完成了
                        if (computing_type == 3) {
                            error_bound_satisfied = true;
                        }
                    }

                    /*接收 近似计算结果*/
                    // 当缓存中有数据的时候接收数据
                    if (approx_tailptr - approx_headptr < BUFFER_SIZE) {
                        MPI_Irecv(&recv_sample_count_and_squra[approx_tailptr % BUFFER_SIZE], 2, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD,
                                  &recv_sample_rqst[approx_tailptr % BUFFER_SIZE]);
                        approx_tailptr = approx_tailptr + 2;
                    }
                    if (approx_headptr < approx_tailptr) // 如果缓存中还有数据 ，将缓存中的数据接收并且更新整体
                    {
                        int testflag = 0;
                        MPI_Test(&recv_sample_rqst[approx_headptr % BUFFER_SIZE], &testflag, MPI_STATUS_IGNORE);
                        if (testflag) {
                            // #pragma omp atomic
                            approx_count += recv_sample_count_and_squra[approx_tailptr % BUFFER_SIZE];
                            approx_squra += recv_sample_count_and_squra[approx_tailptr % BUFFER_SIZE + 1];
                            // #pragma omp atomic
                            approx_sample_times += batch_size;
                            approx_headptr = approx_headptr + 2;

                            // 计算全体方差和每个值得标准差
                            var = (double)(approx_squra / approx_sample_times - std::pow((approx_count / approx_sample_times), 2));
                            std_dev = std::sqrt(var / approx_sample_times);

                            real_error = std_dev * 2.576 / (approx_count / approx_sample_times); // 99%的信心
                            if (real_error < error_bound) {
                                error_bound_satisfied = true;
                            }
                        }
                        // 满足误差要求，更新各个节点的信息
                        if (error_bound_satisfied) {
                            for (uint32_t i = 1; i < num_processes; i++)
                                MPI_Send(&error_bound_satisfied, 1, MPI_C_BOOL, i, 0, MPI_COMM_WORLD);
                            break;
                        }
                        // 长时间不满足，发送精确计算要求
                        if (approx_sample_times % 1000000 == 0 && approx_sample_times != 0) {
                            exact_computing = true;
                            for (int i = 0; i < omp_get_max_threads(); i++)
                                local_exact_computing[i] = true;

                            for (uint32_t i = 1; i < num_processes; i++)
                                MPI_Send(&exact_computing, 1, MPI_C_BOOL, i, 3, MPI_COMM_WORLD);
                        }
                    }
                }
                // 三个任务：接收精确计算要求，和误差满足退出要求，发送精确计算结果
                if (my_rank) {
                    // 如果有精确计算任务，查看是否完成并且发送，否则尝试接收精确计算任务要求
                    while (local_batch_sample_times % 10000 == 0 && local_batch_sample_times != 0) {
                        if (exact_computing) {
                            if (local_node_exact_count_compute_times == omp_get_max_threads() - 1) {
                                MPI_Send(&local_node_exact_count, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
                                local_node_exact_count_compute_times = 0;
                                exact_computing = false;
                            }
                        } else {
                            int exact_comput_flag = 0;
                            MPI_Irecv(&exact_computing, 1, MPI_C_BOOL, 0, 3, MPI_COMM_WORLD, &recv_exact_comput_rqst);
                            MPI_Test(&recv_exact_comput_rqst, &exact_comput_flag, MPI_STATUS_IGNORE);
                            if (exact_comput_flag) {
                                for (int i = 0; i < omp_get_max_threads(); i++) {
                                    local_exact_computing[i] = true;
                                }
                            } else {
                                MPI_Cancel(&recv_exact_comput_rqst);
                                MPI_Wait(&recv_exact_comput_rqst, MPI_STATUS_IGNORE);
                            }
                        }
                    }
                    // 如果不满足退出条件，尝试接收一下退出条件
                    if (!error_bound_satisfied) {
                        // 没有发送消息过来，会一直阻塞。
                        int error_bound_flag = 0;
                        MPI_Irecv(&error_bound_satisfied, 1, MPI_C_BOOL, 0, 0, MPI_COMM_WORLD, &recv_errror_bound_rqst);
                        MPI_Test(&recv_errror_bound_rqst, &error_bound_flag, MPI_STATUS_IGNORE);
                        if (error_bound_flag == 0) {
                            MPI_Cancel(&recv_errror_bound_rqst);
                            MPI_Wait(&recv_errror_bound_rqst, MPI_STATUS_IGNORE);
                        }
                    }
                }
            }
        }

        //==================================================================================================================================================//
        //==========计算任务===============================================================================================================================//
        //==================================================================================================================================================//

        // 两个任务：精确计算结果更新到主线程，近似计算结果异步发送
        if (omp_get_thread_num()) {

            std::cout << "current thread  " << omp_get_thread_num() << std::endl;

            // 发送缓存相关
            MPI_Request send_sample_rqst[BUFFER_SIZE]; // 标记发送请求是否已经发送完毕
            double local_sample_count_and_squra[BUFFER_SIZE] = {0};
            uint32_t headptr = 0, tailptr = 0;

            // 随机数生成器
            std::random_device rd;
            std::default_random_engine gen(rd());
            // TODO:可能需要修改范围
            uint64_t range = num_of_light_edges + num_of_medium_edges + num_of_heavy_edges;
            uint64_t start_edge_id = 0;
            // uint64_t range = num_of_medium_edges + num_of_heavy_edges;
            // uint64_t start_edge_id = num_of_light_edges;
            // uint64_t range = num_of_heavy_edges;
            // uint64_t start_edge_id = num_of_light_edges + num_of_medium_edges;
            std::uniform_int_distribution<uint64_t> dist(0, range - 1);

            while (!error_bound_satisfied) {
                /*精确计算 被动被主线程出发 */
                if (exact_computing && local_exact_computing[omp_get_thread_num()]) {
                    if (computing_type == 0) {
                        while (vertex_range_ptr[0] < vertex_range_end_per_node[0]) {

                            // 通过对于加法操作的原子操作，保证每个线程的vertex_range_ptr[0]的更新是线程安全的。
                            uint32_t local_range_start = __sync_fetch_and_add(&vertex_range_ptr[0], vertex_chunk_per_task);
                            // 还有顶点没有处理完
                            if (local_range_start < vertex_range_end_per_node[0]) {

                                uint32_t local_range_end;
                                if (local_range_start + vertex_chunk_per_task < vertex_range_end_per_node[0])
                                    local_range_end = local_range_start + vertex_chunk_per_task;
                                else
                                    local_range_end = vertex_range_end_per_node[0];

                                int *ans_buffer = new int[schedule_iep.in_exclusion_optimize_vertex_id.size()];  // 底层无交集顶点的大小
                                VertexSet *vertex_set = new VertexSet[schedule_iep.get_total_prefix_num() + 10]; // 构建了多个顶点集合作为候选集
                                VertexSet subtraction_set;                                                       // 构建一个差集
                                VertexSet tmp_set;                                                               // 构建一个临时集合
                                subtraction_set.init();
                                long long local_ans = 0;
                                for (int i = local_range_start; i < local_range_end; ++i) {
                                    uint32_t vertex = filter_data.vertex_project[i];
                                    e_index_t l, r;
                                    l = g->vertex[vertex];
                                    r = g->vertex[vertex + 1];

                                    // 为调度中的每一个顶点初始化候选顶点集合。
                                    for (int prefix_id = schedule_iep.get_last(0); prefix_id != -1; prefix_id = schedule_iep.get_next(prefix_id)) {
                                        vertex_set[prefix_id].build_vertex_set(schedule_iep, vertex_set, &g->edge[l], (int)(r - l), prefix_id);
                                    }
                                    subtraction_set.push_back(vertex);
                                    pattern_matching_aggressive_func(schedule_iep, vertex_set, subtraction_set, tmp_set, local_ans, 1, ans_buffer, g);
                                    subtraction_set.pop_back();
                                }
                                delete[] vertex_set;

#pragma omp atomic
                                local_node_exact_count += local_ans;
                            }
                        }
#pragma omp atomic
                        local_node_exact_count_compute_times++;
                        // 更新状态
                        computing_type = 1;
                        local_exact_computing[omp_get_thread_num()] = false; // 更新计算状态
                        range = num_of_medium_edges + num_of_heavy_edges;    // 更新采样范围
                        start_edge_id = num_of_light_edges;

                    } else if (computing_type == 1) {
                        while (vertex_range_ptr[1] < vertex_range_end_per_node[1]) {

                            // 通过对于加法操作的原子操作，保证每个线程的vertex_range_ptr[0]的更新是线程安全的。
                            uint32_t local_range_start = __sync_fetch_and_add(&vertex_range_ptr[1], vertex_chunk_per_task);
                            // 还有顶点没有处理完
                            if (local_range_start < vertex_range_end_per_node[1]) {

                                uint32_t local_range_end;
                                if (local_range_start + vertex_chunk_per_task < vertex_range_end_per_node[1])
                                    local_range_end = local_range_start + vertex_chunk_per_task;
                                else
                                    local_range_end = vertex_range_end_per_node[1];

                                int *ans_buffer = new int[schedule_iep.in_exclusion_optimize_vertex_id.size()];  // 底层无交集顶点的大小
                                VertexSet *vertex_set = new VertexSet[schedule_iep.get_total_prefix_num() + 10]; // 构建了多个顶点集合作为候选集
                                VertexSet subtraction_set;                                                       // 构建一个差集
                                VertexSet tmp_set;                                                               // 构建一个临时集合
                                subtraction_set.init();
                                long long local_ans = 0;
                                for (int i = local_range_start; i < local_range_end; ++i) {
                                    uint32_t vertex = filter_data.vertex_project[i];
                                    e_index_t l, r;
                                    l = g->vertex[vertex];
                                    r = g->vertex[vertex + 1];

                                    // 为调度中的每一个顶点初始化候选顶点集合。
                                    for (int prefix_id = schedule_iep.get_last(0); prefix_id != -1; prefix_id = schedule_iep.get_next(prefix_id)) {
                                        vertex_set[prefix_id].build_vertex_set(schedule_iep, vertex_set, &g->edge[l], (int)(r - l), prefix_id);
                                    }
                                    subtraction_set.push_back(vertex);
                                    pattern_matching_aggressive_func(schedule_iep, vertex_set, subtraction_set, tmp_set, local_ans, 1, ans_buffer, g);
                                    subtraction_set.pop_back();
                                }
                                delete[] vertex_set;

#pragma omp atomic
                                local_node_exact_count += local_ans;
                            }
                        }
#pragma omp atomic
                        local_node_exact_count_compute_times++;
                        // 更新状态
                        computing_type = 2;
                        local_exact_computing[omp_get_thread_num()] = false; // 更新计算状态
                        range = num_of_heavy_edges;                          // 更新采样范围
                        start_edge_id = num_of_light_edges + num_of_medium_edges;
                    } else if (computing_type == 2) {
                        while (vertex_range_ptr[2] < vertex_range_end_per_node[2]) {

                            // 通过对于加法操作的原子操作，保证每个线程的vertex_range_ptr[0]的更新是线程安全的。
                            uint32_t local_range_start = __sync_fetch_and_add(&vertex_range_ptr[2], vertex_chunk_per_task);
                            // 还有顶点没有处理完
                            if (local_range_start < vertex_range_end_per_node[2]) {

                                uint32_t local_range_end;
                                if (local_range_start + vertex_chunk_per_task < vertex_range_end_per_node[2])
                                    local_range_end = local_range_start + vertex_chunk_per_task;
                                else
                                    local_range_end = vertex_range_end_per_node[2];

                                int *ans_buffer = new int[schedule_iep.in_exclusion_optimize_vertex_id.size()];  // 底层无交集顶点的大小
                                VertexSet *vertex_set = new VertexSet[schedule_iep.get_total_prefix_num() + 10]; // 构建了多个顶点集合作为候选集
                                VertexSet subtraction_set;                                                       // 构建一个差集
                                VertexSet tmp_set;                                                               // 构建一个临时集合
                                subtraction_set.init();
                                long long local_ans = 0;
                                for (int i = local_range_start; i < local_range_end; ++i) {
                                    uint32_t vertex = filter_data.vertex_project[i];
                                    e_index_t l, r;
                                    l = g->vertex[vertex];
                                    r = g->vertex[vertex + 1];

                                    // 为调度中的每一个顶点初始化候选顶点集合。
                                    for (int prefix_id = schedule_iep.get_last(0); prefix_id != -1; prefix_id = schedule_iep.get_next(prefix_id)) {
                                        vertex_set[prefix_id].build_vertex_set(schedule_iep, vertex_set, &g->edge[l], (int)(r - l), prefix_id);
                                    }
                                    subtraction_set.push_back(vertex);
                                    pattern_matching_aggressive_func(schedule_iep, vertex_set, subtraction_set, tmp_set, local_ans, 1, ans_buffer, g);
                                    subtraction_set.pop_back();
                                }
                                delete[] vertex_set;

#pragma omp atomic
                                local_node_exact_count += local_ans;
                            }
                        }
#pragma omp atomic
                        local_node_exact_count_compute_times++;
                        // 更新状态
                        computing_type = 3;
                        local_exact_computing[omp_get_thread_num()] = false; // 更新计算状态
                        range = 1;                                           // 更新采样范围
                        start_edge_id = num_of_light_edges + num_of_medium_edges + num_of_heavy_edges - 1;
                    }
                }

                /*近似计算*/
                // 按照固定的batch计算
                std::uniform_int_distribution<uint64_t> dist(0, range - 1);
                if (tailptr - headptr < BUFFER_SIZE) {
                    double count_sample = 0;
                    double count_sample_square = 0;
                    for (uint64_t i = 0; i < batch_size; i++) {
                        auto random_id = dist(gen);
                        auto eid = filter_data.edge_project[random_id + start_edge_id];
                        pattern_sample_ns(eid, g, s, count_sample, count_sample_square);
                    }

#pragma omp atomic
                    local_batch_sample_times++;

                    local_sample_count_and_squra[tailptr % BUFFER_SIZE] = count_sample;
                    local_sample_count_and_squra[tailptr % BUFFER_SIZE + 1] = count_sample_square;

                    MPI_Isend(&local_sample_count_and_squra[tailptr % BUFFER_SIZE], 2, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD,
                              &send_sample_rqst[tailptr % BUFFER_SIZE]);
                    tailptr = tailptr + 2;
                }

                // 检测是否发送完成，更新发送缓存
                if (headptr < tailptr) {
                    int testflag = 0;
                    MPI_Test(&send_sample_rqst[headptr % BUFFER_SIZE], &testflag, MPI_STATUS_IGNORE); // 异步的操作需要一个test来执行
                    if (testflag)
                        headptr = headptr + 2;
                }
            }
        }
    }

    // 主进程主线程执行
    if (!my_rank) {
        std::cout << "global summary global_sample_times: " << approx_count << " " << approx_sample_times << std::endl;
        double estimated_pattern_count = approx_count * (num_of_heavy_edges + num_of_medium_edges + num_of_light_edges) / (approx_sample_times);
        std::cout << "pattern_nums " << estimated_pattern_count + exact_count << std::endl;
        double estimated_pattern_count_error = (estimated_pattern_count + exact_count - (26184567090560.0 * 2.0)) / (26184567090560.0 * 2.0);
        std::cout << "estimated_error " << estimated_pattern_count_error << std::endl;
    }

    // lj图 house 26184567090560 单线程3000秒，64线程几十秒 scalegpm 1.1秒
    MPI_Finalize();

    free(filter_data.edge_project);
    free(filter_data.edge_record);
    free(filter_data.vertex_project);
    free(filter_data.vertex_record);
    free(filter_data.degree);
    free(filter_data.sample_count_per_edge);
    free(filter_data.sample_times_per_edge);
    free(filter_data.sample_square_per_edge);
}
// 工作方式缺点：每个线程都会去更新主节点的主线程，主线的压力很大。
//  主线程会更新全局采样次数，而工作线程不会更新全局采样次数，因此工作线程会死循环执行，直到主线程跳出循环结束任务 实现错误
//  // 输入图，模式，采样次数和线程数
//  void MPI_pattern_counting(Graph &large_graph, Pattern &graph_pattern, uint64_t sampling_times, uint32_t num_threads) {
//      // 整体状态
//      uint32_t lost_edge_cases[100];

//     // 采样步骤
//     uint64_t sampling_step = 10000;

//     // 全局变量
//     uint64_t global_sampled_times = 0;
//     double global_prob_inverse = 0;

//     int comm_sz, my_rank, num_processes;

//     MPI_Request recvrqst[BUFFER_SIZE];
//     MPI_Status status;
//     double recv_prob_inverse[BUFFER_SIZE] = {0};
//     uint32_t headptr = 0, tailptr = 0;

//     auto start = chrono::high_resolution_clock::now();

// #pragma omp parallel num_threads(num_threads)
//     {
// #pragma omp barrier
// #pragma omp master // 首先每个进程的主线程执行 初始化操作
//         {
//             int provided;
//             MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &provided);
//             MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
//             MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
//             MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
//             MPI_Barrier(MPI_COMM_WORLD);
//         }
// #pragma omp barrier

//         // 每个进程的主线程执行任务，更新采样次数。对于主进程还更新预测数量
// #pragma omp master
//         {
//             while (global_sampled_times < sampling_times) {
//                 if (!my_rank) // 主进程的主线程
//                 {
//                     if (tailptr - headptr < BUFFER_SIZE) // 如果缓存数据较少，继续接收数据
//                     {
//                         MPI_Irecv(&recv_prob_inverse[tailptr % BUFFER_SIZE], 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,
//                                   &recvrqst[tailptr % BUFFER_SIZE]);
//                         tailptr++;
//                     }
//                     if (headptr < tailptr) // 如果缓存数据足够多了，执行更新操作。采样次数够了，退出任务，否择更新预测数量和采样次数
//                     {
//                         if (global_sampled_times >= sampling_times)
//                             break;
//                         int testflag = 0;
//                         MPI_Test(&recvrqst[headptr % BUFFER_SIZE], &testflag, MPI_STATUS_IGNORE);
//                         if (testflag) {
// #pragma omp atomic
//                             global_prob_inverse += recv_prob_inverse[tailptr % BUFFER_SIZE];
// #pragma omp atomic
//                             global_sampled_times += sampling_step;
//                             headptr++;
//                         }
//                     }
//                 }
//                 if (my_rank) // 非主进程的主线程执行，更新全局的采样次数
//                     MPI_Recv(&global_sampled_times, 1, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD, &status);
//             }
//         }

//         // 每个进程的其他线程执行的任务
//         if (omp_get_thread_num()) {
//             MPI_Request sendrqst[BUFFER_SIZE]; // 用于标记某个发送请求是否已经发送完毕
//             double local_prob_inverse[BUFFER_SIZE] = {0};
//             uint64_t local_sampled_times = 0;
//             uint32_t headptr = 0, tailptr = 0;

//             random_device sd;
//             default_random_engine this_rand_generator(sd());

//             while (global_sampled_times < sampling_times) // 采样数量不够时
//             {
//                 if (tailptr - headptr < BUFFER_SIZE) // 发送内存存在空间，继续采样并发送
//                 {
//                     tie(local_prob_inverse[tailptr % BUFFER_SIZE], local_sampled_times) =
//                         estimate_pattern(large_graph, graph_pattern, sampling_step, lost_edge_cases, this_rand_generator);
//                     MPI_Isend(&local_prob_inverse[tailptr % BUFFER_SIZE], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
//                               &sendrqst[tailptr % BUFFER_SIZE]); // 最后一个位置时标记，标记该操作是否已经完成
//                     tailptr++;
//                 }
//                 if (headptr < tailptr) // 如果缓存用完了，查看是否已经发送完成了。 一个异步的更新过程。
//                 {
//                     int testflag = 0;
//                     MPI_Test(&sendrqst[headptr % BUFFER_SIZE], &testflag, MPI_STATUS_IGNORE); // 异步的操作需要一个test来执行
//                     if (testflag)
//                         headptr++;
//                     if (global_sampled_times >= sampling_times)
//                         break;
//                 }
//             }
//         }
//     }

//     auto stop = chrono::high_resolution_clock::now();
//     auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);

//     // 从主进程上发送全局采样次数给所有进程，发送完就结束了？   为什么这一段程序在这里？ 存在依赖操作的应该按照顺序写
//     if (!my_rank) {
//         for (uint32_t i = 1; i < num_processes; i++)
//             MPI_Send(&global_sampled_times, 1, MPI_UNSIGNED_LONG_LONG, i, 0, MPI_COMM_WORLD);
//         cout << setprecision(16) << "raw results after merge nodes (total_prob_inverse, total_sampled_times): " << global_prob_inverse << "
//         "
//              << global_sampled_times << endl;
//         double estimated_pattern_count = global_prob_inverse / global_sampled_times;
//         cout << "pattern_nums time_consumed(us) sampling_times sampling/second\n"
//              << estimated_pattern_count << " " << duration.count() << " " << global_sampled_times << " "
//              << global_sampled_times * 1.0 / duration.count() * 1000000 << endl;
//     }

//     MPI_Finalize();
// }

void print_stata(Graph *g, Pattern &pattern, schedule_approx &schedule) {
    for (int i = 0; i < pattern.get_size(); i++) {
        printf("order: %d, degree: %d\n", schedule.order[i], schedule.order_degree[i]);
    }

    for (int i = 0; i < pattern.get_size(); i++) {
        printf("vertex: %d: ", i);
        for (int ptr = schedule.dependcy_ptr[i]; ptr < schedule.dependcy_ptr[i + 1]; ptr++) {
            printf(" dependcy: %d", schedule.dependcy[ptr]);
        }
        printf("\n");
    }
    for (int i = 0; i < pattern.get_size(); i++) {
        printf("vertex: %d: ", i);
        for (int ptr = schedule.sub_ptr[i]; ptr < schedule.sub_ptr[i + 1]; ptr++) {
            printf(" sub: %d", schedule.sub_id[ptr]);
        }
        printf("\n");
    }
}

void scheduel_generate_for_adap_ns(const Schedule_IEP &sched, schedule_approx *schedule) {
    // 先对矩阵进行排序
    int size = sched.get_size();
    int *adj_mat = new int[size * size];
    memcpy(adj_mat, sched.get_adj_mat_ptr(), size * size * sizeof(int));

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

    // 获取模式中连接的状态
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
    schedule->dependcy = new int[size * size];
    schedule->dependcy_ptr = new int[size + 1];
    schedule->dependcy_label = new int[size];
    memset(schedule->dependcy_label, -1, sizeof(int) * size);
    schedule->compute_depth = size;
    schedule->pattern_size = size;
    if (restricts_vector.size() == 0) {
        schedule->symmetry_label = false;
    } else {
        schedule->symmetry_label = true;
    }

    schedule->sub_id = new int[size * size];
    schedule->sub_ptr = new int[size + 1];

    // 适配新的模式
    int *already_assigned = new int[size];
    memset(already_assigned, 0, sizeof(int) * size);
    schedule->order[0] = 0;
    schedule->order[1] = 1;
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

    /* 确定好顶点处理顺序后，构建依赖关系和差集信息 */
    schedule->dependcy_ptr[0] = 0;
    schedule->dependcy_ptr[1] = 0;
    schedule->dependcy_ptr[2] = 0;

    schedule->sub_ptr[0] = 0;
    schedule->sub_ptr[1] = 0;
    schedule->sub_ptr[2] = 0;
    for (int process_vtx = 2; process_vtx < size; process_vtx++) {
        schedule->dependcy_ptr[process_vtx + 1] = schedule->dependcy_ptr[process_vtx];
        schedule->sub_ptr[process_vtx + 1] = schedule->sub_ptr[process_vtx];
        for (int i = 0; i < process_vtx; i++) {
            if (adj_mat[INDEX(process_vtx, i, size)]) {
                schedule->dependcy[schedule->dependcy_ptr[process_vtx + 1]++] = i;
            } else {
                schedule->sub_id[schedule->sub_ptr[process_vtx + 1]++] = i;
            }
        }
    }

    /*构建好依赖关系之后，转化为记录已有计算结果，此处实现的是完全复用 */
    VertexSet tmpset[2];
    int similar_vtx, max_similar_size;
    for (int process_vtx = size - 1; process_vtx >= 2; process_vtx--) {

        tmpset[0].init(schedule->dependcy_ptr[process_vtx + 1] - schedule->dependcy_ptr[process_vtx],
                       schedule->dependcy + schedule->dependcy_ptr[process_vtx]);
        max_similar_size = 0;
        for (int pre_vtx = process_vtx - 1; pre_vtx >= 2; pre_vtx--) {
            tmpset[1].init(schedule->dependcy_ptr[pre_vtx + 1] - schedule->dependcy_ptr[pre_vtx],
                           schedule->dependcy + schedule->dependcy_ptr[pre_vtx]);
            int temp_similar_size = intersection_size(tmpset[0], tmpset[1]);
            if (temp_similar_size > max_similar_size && tmpset[1].get_size() == temp_similar_size) {
                max_similar_size = temp_similar_size;
                similar_vtx = pre_vtx;
            }
        }

        // 如果有可以复用的结果，对其进行复用
        if (max_similar_size >= 2) {
            tmpset[1].init(schedule->dependcy_ptr[similar_vtx + 1] - schedule->dependcy_ptr[similar_vtx],
                           schedule->dependcy + schedule->dependcy_ptr[similar_vtx]);

            int record_ptr = schedule->dependcy_ptr[process_vtx + 1];
            for (int ptr = tmpset[0].get_size() - 1; ptr >= 0; ptr--) {
                int temp_vtx = tmpset[0].get_data(ptr);
                if (!tmpset[1].has_data(temp_vtx)) {
                    schedule->dependcy[--record_ptr] = temp_vtx;
                }
            }
            while (record_ptr != schedule->dependcy_ptr[process_vtx]) {
                schedule->dependcy[--record_ptr] = -1;
            }
            schedule->dependcy[schedule->dependcy_ptr[process_vtx]] = similar_vtx;
            schedule->dependcy_label[process_vtx] = 0;
        }
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

static double phi(double x) {
    // constants
    double a1 = 0.254829592;
    double a2 = -0.284496736;
    double a3 = 1.421413741;
    double a4 = -1.453152027;
    double a5 = 1.061405429;
    double p = 0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x) / sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0 / (1.0 + p * x);
    double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);

    return 0.5 * (1.0 + sign * y);
}

// binary search for near-match
static double inv_phi(double p) {
    double Z_EPSILON = 0.000001;
    double min = 0, max = 6;
    do {
        double guess = min + (max - min) / 2;
        double result = phi(guess);
        if (result > p)
            max = guess;
        if (result < p)
            min = guess;
    } while (min + Z_EPSILON < max);
    return min + (max - min) / 2;
}

static double confidence_to_z(double c) { return inv_phi(1 - ((1 - c) / 2)); }

void exact_matching(Graph *g, const Schedule_IEP &schedule_iep) {
#pragma omp parallel num_threads(7)
    {
        int *ans_buffer = new int[schedule_iep.in_exclusion_optimize_vertex_id.size()];  // 底层无交集顶点的大小
        VertexSet *vertex_set = new VertexSet[schedule_iep.get_total_prefix_num() + 10]; // 构建了多个顶点集合作为候选集
        VertexSet subtraction_set;                                                       // 构建一个差集
        VertexSet tmp_set;                                                               // 构建一个临时集合
        subtraction_set.init();
        long long local_ans = 0;

#pragma omp for
        for (uint32_t vertex = 0; vertex < g->v_cnt; ++vertex) {
            e_index_t l, r;
            l = g->vertex[vertex];
            r = g->vertex[vertex + 1];

            // 为调度中的每一个顶点初始化候选顶点集合。
            for (int prefix_id = schedule_iep.get_last(0); prefix_id != -1; prefix_id = schedule_iep.get_next(prefix_id)) {
                vertex_set[prefix_id].build_vertex_set(schedule_iep, vertex_set, &g->edge[l], (int)(r - l), prefix_id);
            }
            subtraction_set.push_back(vertex);
            pattern_matching_aggressive_func(schedule_iep, vertex_set, subtraction_set, tmp_set, local_ans, 1, ans_buffer, g);
            subtraction_set.pop_back();
        }
        delete[] vertex_set;
    }
}

// 阈值设置，百分之五十原则

void print_error(double real_error, uint64_t approx_sample_times, bool *flag, int recur_depth, double real_error_local[]) {

    std::cout << "error " << real_error << std::endl;
    for (int i = 0; i < recur_depth; i++)
        printf("real_error_local[%d] = %f\n", i, real_error_local[i]);

    if (real_error < 0.1 && flag[0]) {
        std::cout << "error bound 10% sample times " << approx_sample_times << std::endl;
        flag[0] = false;
    }
    if (real_error < 0.05 && flag[1]) {
        std::cout << "error bound 5% sample times " << approx_sample_times << std::endl;
        flag[1] = false;
    }
    if (real_error < 0.01 && flag[2]) {
        std::cout << "error bound 1% sample times " << approx_sample_times << std::endl;
        flag[2] = false;
    }
    if (real_error < 0.001 && flag[3]) {
        std::cout << "error bound 0.1% sample times " << approx_sample_times << std::endl;
        flag[3] = false;
    }
}

void adaptive_approximate_exact_fusion(Graph *g, schedule_approx *s, const Schedule_IEP &schedule_iep, Pattern p) {

    edge_from_init(g);

    // 平均度数  frendster 28 twitter 35 livejouranl 17 youtube 2 patent 5 pokec 30 tp 28 sk 11
    // 阈值划分10次 只保留千分之一的顶点 结果波动还很大，
    int recur_depth = 50;
    int vertex_depth = 0;

    uint32_t threshhold[recur_depth + 1] = {0};
    int *record_degree = (int *)malloc(sizeof(int) * g->v_cnt);
    memset(record_degree, 0, sizeof(int) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        record_degree[g->vertex[i + 1] - g->vertex[i]]++;
    }
    double ratio = 0;
    double ratio_threshhold = 0.5;
    double ratio_acum = 0.5;
    // threshhold[0] = 1;
    int iter_threshhoud = 1;
    uint64_t iter_record = 0;
    double change_ratio = 0.0;

    while (iter_record < g->v_cnt) {
        double temp_ratio = double(record_degree[iter_record]) / double(g->v_cnt);
        ratio += temp_ratio;
        iter_record++;
        if (ratio > ratio_threshhold) {
            threshhold[iter_threshhoud++] = iter_record;
            // change_ratio = double(threshhold[iter_threshhoud - 1]) / double(threshhold[iter_threshhoud - 2]);
            ratio_acum = ratio_acum * 0.5;
            ratio_threshhold = ratio_threshhold + ratio_acum;
            if (iter_threshhoud == recur_depth || iter_record > 20000 || iter_threshhoud > 15) { //|| change_ratio < 1.1
                break;
            }
        }
    }
    for (int i = 1; i < recur_depth; i++) {
        if (threshhold[i] == 0) {
            recur_depth = i;
            break;
        }
    }
    threshhold[recur_depth] = g->v_cnt;
    free(record_degree);

    // 一共有recur_depth中颜色， 但是查找阈值的时候需要用recur_depth +1.因为阈值第一位为0

    // // 查看阈值
    // for (int i = 0; i < recur_depth; i++) {
    //     printf("threshhold[%d] = %d\n", i, threshhold[i]);
    // }

    // 过滤的数据结构。
    filter filter_data;
    filter_data.edge_project = (int64_t *)malloc(sizeof(int64_t) * g->e_cnt);
    filter_data.edge_record = (int *)malloc(sizeof(int) * g->e_cnt);
    filter_data.vertex_project = (int *)malloc(sizeof(int) * g->v_cnt);
    filter_data.vertex_record = (int *)malloc(sizeof(int) * g->v_cnt);
    memset(filter_data.edge_project, -1, sizeof(int64_t) * g->e_cnt);
    memset(filter_data.edge_record, -1, sizeof(int) * g->e_cnt);
    memset(filter_data.vertex_record, -1, sizeof(int) * g->v_cnt);
    memset(filter_data.vertex_project, -1, sizeof(int) * g->v_cnt);
    filter_data.delete_edge = (bool *)malloc(sizeof(bool) * g->e_cnt);
    memset(filter_data.delete_edge, false, sizeof(bool) * g->e_cnt);

    filter_data.degree = (uint32_t *)malloc(sizeof(uint32_t) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        filter_data.degree[i] = g->vertex[i + 1] - g->vertex[i];
    }

    uint64_t discard_vertex = 0;

    // 标记顶点和边，以及记录数据
    uint64_t num_of_different_weight_of_edges[recur_depth] = {0};
    uint64_t num_of_different_weight_of_vertex[recur_depth] = {0};
    for (uint64_t i = 0; i < g->v_cnt; ++i) {
        if (filter_data.degree[i] >= s->order_degree[0]) {
            int record_vertex_color = find_range(threshhold, filter_data.degree[i], recur_depth + 1);
            filter_data.vertex_record[i] = record_vertex_color;
            num_of_different_weight_of_vertex[record_vertex_color]++;
            for (uint64_t j = g->vertex[i]; j < g->vertex[i + 1]; ++j) {
                uint32_t temp_degree = filter_data.degree[g->edge[j]];
                if (temp_degree >= s->order_degree[1]) {
                    int record_edge_color = find_range(threshhold, temp_degree, recur_depth + 1);
                    filter_data.edge_record[j] = record_edge_color;
                    num_of_different_weight_of_edges[record_edge_color]++;
                }
            }
        } else {
            discard_vertex++;
        }
    }

    // 投射顶点，记录所有的顶点
    uint64_t prefix_sum_of_vertex[recur_depth + 1] = {0};
    for (int i = 1; i < recur_depth + 1; i++) {
        prefix_sum_of_vertex[i] = prefix_sum_of_vertex[i - 1] + num_of_different_weight_of_vertex[i - 1];
    }
    uint64_t weight_vtx_ptr[recur_depth] = {0};
    for (int i = 1; i < recur_depth; i++) {
        weight_vtx_ptr[i] = prefix_sum_of_vertex[i];
    }
    for (uint32_t i = 0; i < g->v_cnt; ++i) {
        if (filter_data.vertex_record[i] != -1)
            filter_data.vertex_project[weight_vtx_ptr[filter_data.vertex_record[i]]++] = i;
    }

    // 投射边
    // 前缀和形式
    uint64_t prefix_sum_of_edges[recur_depth + 1] = {0};
    for (int i = 1; i < recur_depth + 1; i++) {
        prefix_sum_of_edges[i] = prefix_sum_of_edges[i - 1] + num_of_different_weight_of_edges[i - 1];
    }
    // 利用前缀和更新指针
    uint64_t weight_edge_ptr[recur_depth] = {0};
    for (int i = 0; i < recur_depth; i++) {
        weight_edge_ptr[i] = prefix_sum_of_edges[i];
    }
    for (uint32_t i = 0; i < g->e_cnt; ++i) {
        if (filter_data.edge_record[i] != -1)
            filter_data.edge_project[weight_edge_ptr[filter_data.edge_record[i]]++] = i;
    }

    // 记录采样结果
    filter_data.sample_times_per_edge = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
    filter_data.sample_count_per_edge = (double *)malloc(sizeof(double) * g->e_cnt);
    filter_data.sample_square_per_edge = (double *)malloc(sizeof(double) * g->e_cnt);
    memset(filter_data.sample_times_per_edge, 0, sizeof(uint32_t) * g->e_cnt);
    memset(filter_data.sample_count_per_edge, 0, sizeof(double) * g->e_cnt);
    memset(filter_data.sample_square_per_edge, 0, sizeof(double) * g->e_cnt);

    // 采样批量
    uint64_t batch_size = 1000;

    // 随机数生成器
    std::random_device rd;
    std::default_random_engine gen(rd());

    // 第一次采样的范围 ，*****后续要实时更新
    uint64_t range = prefix_sum_of_edges[recur_depth];
    // 采样记录
    uint64_t approx_sample_times = 0;
    uint64_t real_sample_times = 0;
    double approx_count = 0;
    double approx_squra = 0;
    uint64_t mask = 1ULL;
    mask <<= 63;
    bool flag[4] = {true, true, true, true};

    // 主线程判断是否满足要求使用到的变量，只有主线程会使用到

    double var_local[recur_depth] = {0};
    double std_dev_local[recur_depth] = {0};
    double estimate_error_local[recur_depth] = {0};
    uint64_t approx_sample_times_local[recur_depth] = {0};
    double approx_count_local[recur_depth] = {0};
    double approx_squra_local[recur_depth] = {0};
    double real_error_gloabl = 0.0;
    double real_error = 1.0;

    // 线程工作范围
    int num_of_thread = omp_get_max_threads();
    uint64_t start_edge_id[num_of_thread] = {0};
    uint64_t end_edge_id[num_of_thread];
    for (int i = 0; i < num_of_thread; i++) {
        end_edge_id[i] = prefix_sum_of_edges[recur_depth];
    }

    // 精确计算相关变量
    uint32_t processing_round = 0;
    uint32_t print_times = 0;
    bool exact_computing = false;
    bool get_ready[num_of_thread] = {false};
    uint64_t vertex_range_ptr = 0;
    uint64_t vertex_range_end = 0;
    uint64_t exact_count = 0;
    uint64_t vertex_chunk_per_task = 256;

    // 收敛速率检验相关
    double pre_real_error = 1.0;        // 记录前一个误差
    double pre_approx_sample_times = 0; // 记录前一个采样次数
    double coverge_range = 10000.0;     // 检查的范围差距根据处理速率不断调整。
    double coverge_range_ratio = 1.1;
    double coverge_range_ratio_motif5 = 1.5;

    // 分层的总结
    double weight_average = 0.0;
    double weight_squra = 0.0;
    double weight_std_error = 0.0;
    double weight_real_error = 0.0;

    double error_threshhold[5] = {0.0};
    double var_threshhold[recur_depth] = {0.0};
    double var_threshhold_pre[recur_depth] = {0.0};

    double multiplicity = schedule_iep.get_multiplicity();
    double iep_redun = schedule_iep.get_in_exclusion_optimize_redundancy();

    double var = 0.0;
    double real_average = 0.0;

    // 根据iep数和稀疏程度判断是否难以计算
    // double iep_num = schedule_iep.get_in_exclusion_optimize_num();
    // bool sparse = false;
    int start_part_id = 0;
    for (int i = 0; i < recur_depth; i++) {
        if (num_of_different_weight_of_edges[i] > 0) {
            start_part_id = i;
            break;
        }
    }

    // 误差 边界
#pragma omp parallel
    {
// 可以使用master来实现
#pragma omp master
        {
            while (real_error > 0.05 && approx_sample_times <= 100000000) {

                // 1. 统计全局误差和局部误差。
                // 2. 统计误差收敛速率
                // 3. 如果误差收敛速率够大，自适应采样调整范围
                // 4. 如果误差收敛速率过小，精确计算，重新调整采样范围。
                if (approx_sample_times > 1000 * print_times) {

                    //==================================================================================================================================================//
                    //=================误差统计=========================================================================================================================//
                    //==================================================================================================================================================//
                    /* 局部误差 */ // 记录平方和，计数和，局部采样次数
                    for (int i = start_part_id; i < recur_depth; i++) {
                        approx_sample_times_local[i] = 0.0;
                        approx_count_local[i] = 0.0;
                        approx_squra_local[i] = 0.0;
                    }
                    for (int i = start_part_id; i < recur_depth; i++) {
                        for (int ptr = prefix_sum_of_edges[i]; ptr < prefix_sum_of_edges[i + 1]; ptr++) {
                            auto eid = filter_data.edge_project[ptr];
                            approx_sample_times_local[i] += filter_data.sample_times_per_edge[eid];
                            approx_count_local[i] += filter_data.sample_count_per_edge[eid];
                            approx_squra_local[i] += filter_data.sample_square_per_edge[eid];
                        }
                    }
                    // 局部均值的偏差
                    for (int i = start_part_id; i < recur_depth; i++) {
                        var_local[i] = (double)(approx_squra_local[i] / approx_sample_times_local[i] -
                                                std::pow((approx_count_local[i] / approx_sample_times_local[i]), 2));
                    }

                    // 局部总和的标准偏差
                    // for (int i = 0; i < recur_depth; i++) {
                    //     double cur_errror = std::sqrt(var_local[i] / approx_sample_times_local[i]) * num_of_different_weight_of_edges[i];
                    //     std::cout << " 线程" << i << "偏差 " << cur_errror << std::endl;
                    // }

                    /* 统计误差 */
                    // 总体方差等于局部方差除以采样次数之和， 总体均值为局部均值之和
                    var = 0.0;
                    real_average = 0.0;
                    for (int i = start_part_id; i < recur_depth; i++) {
                        var +=
                            (var_local[i] * num_of_different_weight_of_edges[i] * num_of_different_weight_of_edges[i]) / approx_sample_times_local[i];
                        real_average += (approx_count_local[i] / approx_sample_times_local[i] * num_of_different_weight_of_edges[i]);
                    }
                    real_average += (double)exact_count * multiplicity / iep_redun;
                    double std_dev = std::sqrt(var);
                    real_error = std_dev * 2.576 / real_average; // 99%的信心

                    if (std::isnan(real_error)) {
                        real_error = 1.0;
                        continue;
                    }

                    // double realreal_error = (real_average - (7375094981 * 2)) / (7375094981 * 2); // 7375094981 5007311619323
                    // double realreal_error = (real_average - (13283084785.0 * multiplicity)) / (13283084785.0 * multiplicity);

                    // printf("real_average: %lf, approx_sample_times: %lu, estimate_error: %lf\n", real_average, approx_sample_times, real_error);

                    if (flag[0] && real_error < 0.1) {
                        printf("approx_count: %lf, approx_sample_times: %lu, real_error: %lf\n", real_average, approx_sample_times, real_error);
                        flag[0] = false;
                    }
                    if (flag[1] && real_error < 0.05) {
                        printf("approx_count: %lf, approx_sample_times: %lu, real_error: %lf\n", real_average, approx_sample_times, real_error);
                        flag[1] = false;
                    }

                    print_times++;

                    //==================================================================================================================================================//
                    //===========检查误差收敛速率===============================================================================================================================//
                    //==================================================================================================================================================//

                    /* 检验误差收敛速率 */
                    // 检查的间隔逐渐增加。当误差上升的时候，进行精确计算，更新边采样的范围
                    if (approx_sample_times - pre_approx_sample_times > coverge_range) {
                        pre_approx_sample_times = approx_sample_times;

                        if (coverge_range < 1000000) {
                            if (p.get_size() <= 4)
                                coverge_range *= coverge_range_ratio;
                            else if (p.get_size() <= 5)
                                coverge_range *= coverge_range_ratio_motif5;
                        }

                        if ((processing_round >= 5 && p.get_size() >= 5) || p.get_size() >= 6) {
                            coverge_range = 10000000;
                        }

                        if (real_error > pre_real_error || pre_real_error - real_error < 0.001) {
                            exact_computing = true;
                            pre_real_error = real_error;
                        } else {
                            pre_real_error = real_error;
                        }
                    }

                    // 1，更新精确计算范围。 2，打开精确计算开关。3，可能存在前一个精确计算还没有结束就进入精确计算的情况。
                    if (exact_computing) {
                        for (int i = 1; i < num_of_thread; i++) {
                            if (get_ready[i] == true) {
                                i--;
                            }
                        }

                        // 分配精确计算的范围
                        vertex_range_ptr = prefix_sum_of_vertex[processing_round];     // 顶点处理指针
                        vertex_range_end = prefix_sum_of_vertex[processing_round + 1]; // 顶点处理范围

                        // 更新精确计算的粒度
                        vertex_chunk_per_task = vertex_chunk_per_task / 2;
                        if (vertex_chunk_per_task <= 1) {
                            vertex_chunk_per_task = 1;
                        }

                        // 打开所有线程的精确计算开关
                        for (int i = 1; i < num_of_thread; i++) {
                            get_ready[i] = true;
                        }

                        // 最后一部分顶点，查看所有精确计算是否完成任务
                        if (processing_round + 1 == recur_depth) {
                            for (int i = 1; i < num_of_thread; i++) {
                                if (get_ready[i] == true) {
                                    i--;
                                }
                            }

                            real_error = 0.0;
                            break;
                        }

                        // 记录现有采样范围，以及更新边范围相关数据
                        uint64_t size = prefix_sum_of_edges[recur_depth];
                        for (int i = 0; i < recur_depth; i++) {
                            num_of_different_weight_of_edges[i] = 0;
                            prefix_sum_of_edges[i] = 0;
                        }

                        uint64_t temp_ptr = 0;
                        int pre_edge_color = 0;
                        int cur_edge_color = 0;
                        int cur_vtx_color = 0;
                        for (int i = 0; i < size; i++) {
                            uint64_t temp_edge = filter_data.edge_project[i];
                            if (filter_data.delete_edge[temp_edge]) {
                                continue;
                            }
                            uint64_t temp_degree = filter_data.degree[g->edge_from[temp_edge]];
                            cur_vtx_color = find_range(threshhold, temp_degree, recur_depth + 1);
                            // 起点度数够大，保留
                            if (cur_vtx_color > processing_round) {
                                // 记录边并且获取边的颜色
                                filter_data.edge_project[temp_ptr++] = temp_edge;
                                cur_edge_color = find_range(threshhold, filter_data.degree[g->edge[temp_edge]], recur_depth + 1);

                                // 更新相应颜色的范围
                                if (cur_edge_color != pre_edge_color) {
                                    pre_edge_color = cur_edge_color;
                                    num_of_different_weight_of_edges[cur_edge_color]++;
                                    prefix_sum_of_edges[cur_edge_color] = temp_ptr - 1;

                                } else {
                                    num_of_different_weight_of_edges[pre_edge_color]++;
                                }
                            }
                        }
                        prefix_sum_of_edges[recur_depth] = temp_ptr;
                        processing_round++;
                    }

                    //==================================================================================================================================================//
                    //=========最后调整工作线程的工作范围==================================================================================================================================//
                    //==================================================================================================================================================//

                    /* 自适应的修改采样范围 */
                    // 记录 偏差阈值 迭代上升。
                    for (int i = 0; i < recur_depth; i++) {
                        var_threshhold[i] = std::sqrt(var_local[i] / approx_sample_times_local[i]) * num_of_different_weight_of_edges[i] /
                                            (approx_count_local[i] / approx_sample_times_local[i]);
                        var_threshhold_pre[i] = var_threshhold[i];
                    }
                    std::sort(var_threshhold, var_threshhold + recur_depth);
                    // for (int i = 0; i < 5; i++) {
                    //     int var_threshhold_index = recur_depth - std::pow(2, (i + 1)) + 1;
                    //     if (var_threshhold_index >= 0)
                    //         error_threshhold[i] = var_threshhold[var_threshhold_index] - 1;
                    // }

                    int remain_thread = num_of_thread;
                    int allocate_thread_ptr = 0;
                    int allocate_thread = 0;
                    int allocate_threshhold = 0;
                    int adaptive_start = 0;
                    while (remain_thread > num_of_thread / 2) {
                        allocate_thread = remain_thread / 4;
                        // 查找当前应该分配的阈值的下标。
                        for (int i = 0; i < recur_depth; i++) {
                            if (var_threshhold_pre[i] > error_threshhold[allocate_threshhold]) {
                                adaptive_start = i;
                            }
                        }
                        allocate_threshhold++;
                        for (int i = allocate_thread_ptr; i < allocate_thread_ptr + allocate_thread; i++) {
                            start_edge_id[i] = prefix_sum_of_edges[adaptive_start];
                            end_edge_id[i] = prefix_sum_of_edges[recur_depth];
                        }
                        remain_thread -= allocate_thread;
                        allocate_thread_ptr += allocate_thread;
                    }
                    for (int i = allocate_thread_ptr; i < num_of_thread; i++) {
                        end_edge_id[i] = prefix_sum_of_edges[recur_depth];
                    }
                }
            }
            if (flag[2] && real_error <= 0.01) {
                printf("approx_count: %lf, approx_sample_times: %lu, real_error: %lf\n", real_average, approx_sample_times, real_error);
                flag[2] = false;
            }
        }
        if (omp_get_thread_num()) {
            while (real_error > 0.05 && approx_sample_times <= 100000000) {
                int cur_thread = omp_get_thread_num();
                // 如果需要精确计算
                if (get_ready[cur_thread]) {
                    while (vertex_range_ptr < vertex_range_end) {
                        // 获取当前线程的处理范围
                        uint32_t local_range_start = __sync_fetch_and_add(&vertex_range_ptr, vertex_chunk_per_task);
                        if (local_range_start < vertex_range_end) {
                            uint32_t local_range_end;
                            if (local_range_start + vertex_chunk_per_task < vertex_range_end)
                                local_range_end = local_range_start + vertex_chunk_per_task;
                            else
                                local_range_end = vertex_range_end;

                            /* 精确计算阶段 */
                            int *ans_buffer = new int[schedule_iep.in_exclusion_optimize_vertex_id.size()];  // 底层无交集顶点的大小
                            VertexSet *vertex_set = new VertexSet[schedule_iep.get_total_prefix_num() + 10]; // 构建了多个顶点集合作为候选集
                            VertexSet subtraction_set;                                                       // 构建一个差集
                            VertexSet tmp_set;                                                               // 构建一个临时集合
                            subtraction_set.init();
                            long long local_ans = 0;
                            for (int i = local_range_start; i < local_range_end; ++i) {
                                uint32_t vertex = filter_data.vertex_project[i];
                                e_index_t l, r;
                                l = g->vertex[vertex];
                                r = g->vertex[vertex + 1];

                                // 为调度中的每一个顶点初始化候选顶点集合。
                                for (int prefix_id = schedule_iep.get_last(0); prefix_id != -1; prefix_id = schedule_iep.get_next(prefix_id)) {
                                    vertex_set[prefix_id].build_vertex_set(schedule_iep, vertex_set, &g->edge[l], (int)(r - l), prefix_id);
                                }
                                subtraction_set.push_back(vertex);
                                pattern_matching_aggressive_func(schedule_iep, vertex_set, subtraction_set, tmp_set, local_ans, 1, ans_buffer, g);
                                subtraction_set.pop_back();
                            }
                            delete[] vertex_set;
                            delete[] ans_buffer;

#pragma omp atomic
                            exact_count += local_ans;
                        }
                    }

                    get_ready[cur_thread] = false;
                    exact_computing = false;
                }

                std::uniform_int_distribution<uint64_t> dist(start_edge_id[cur_thread],
                                                             end_edge_id[cur_thread] - 1); //***采样范围的更新可以不需要原子操作
                for (uint64_t i = 0; i < batch_size; i++) {
                    auto random_id = dist(gen);
                    auto eid = filter_data.edge_project[random_id]; // 一开始就处理最后一个区间，但是会不断过滤
                    pattern_sample_ns_record_filter(eid, g, s, &filter_data);
                    // pattern_sample_ns_record(eid, g, s, &filter_data);
                }

#pragma omp atomic
                approx_sample_times += batch_size;
            }
        }
    }

    free(filter_data.edge_project);
    free(filter_data.edge_record);
    free(filter_data.vertex_project);
    free(filter_data.vertex_record);
    free(filter_data.degree);

    free(filter_data.sample_count_per_edge);
    free(filter_data.sample_times_per_edge);
    free(filter_data.sample_square_per_edge);
    free(filter_data.delete_edge);
}

void approximate_exact_fusion_fifty(Graph *g, schedule_approx *s, const Schedule_IEP &schedule_iep, Pattern p) {

    edge_from_init(g);

    // 平均度数  frendster 28 twitter 35 livejouranl 17 youtube 2 patent 5 pokec 30 tp 28 sk 11
    // 阈值划分10次 只保留千分之一的顶点 结果波动还很大，
    int recur_depth = 50;
    int vertex_depth = 0;

    uint32_t threshhold[recur_depth + 1] = {0};
    int *record_degree = (int *)malloc(sizeof(int) * g->v_cnt);
    memset(record_degree, 0, sizeof(int) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        record_degree[g->vertex[i + 1] - g->vertex[i]]++;
    }
    double ratio = 0;
    double ratio_threshhold = 0.5;
    double ratio_acum = 0.5;
    // threshhold[0] = 1;
    int iter_threshhoud = 1;
    uint64_t iter_record = 0;
    double change_ratio = 0.0;

    while (iter_record < g->v_cnt) {
        double temp_ratio = double(record_degree[iter_record]) / double(g->v_cnt);
        ratio += temp_ratio;
        iter_record++;
        if (ratio > ratio_threshhold) {
            threshhold[iter_threshhoud++] = iter_record;
            // change_ratio = double(threshhold[iter_threshhoud - 1]) / double(threshhold[iter_threshhoud - 2]);
            ratio_acum = ratio_acum * 0.5;
            ratio_threshhold = ratio_threshhold + ratio_acum;
            if (iter_threshhoud == recur_depth || iter_record > 20000 || iter_threshhoud > 15) { //|| change_ratio < 1.1
                break;
            }
        }
    }
    for (int i = 1; i < recur_depth; i++) {
        if (threshhold[i] == 0) {
            recur_depth = i;
            break;
        }
    }
    threshhold[recur_depth] = g->v_cnt;
    free(record_degree);

    // 一共有recur_depth中颜色， 但是查找阈值的时候需要用recur_depth +1.因为阈值第一位为0

    // 查看阈值
    // for (int i = 0; i < recur_depth; i++) {
    //     printf("threshhold[%d] = %d\n", i, threshhold[i]);
    // }

    // 过滤的数据结构。
    filter filter_data;
    filter_data.edge_project = (int64_t *)malloc(sizeof(int64_t) * g->e_cnt);
    filter_data.edge_record = (int *)malloc(sizeof(int) * g->e_cnt);
    filter_data.vertex_project = (int *)malloc(sizeof(int) * g->v_cnt);
    filter_data.vertex_record = (int *)malloc(sizeof(int) * g->v_cnt);
    memset(filter_data.edge_project, -1, sizeof(int64_t) * g->e_cnt);
    memset(filter_data.edge_record, -1, sizeof(int) * g->e_cnt);
    memset(filter_data.vertex_record, -1, sizeof(int) * g->v_cnt);
    memset(filter_data.vertex_project, -1, sizeof(int) * g->v_cnt);
    filter_data.delete_edge = (bool *)malloc(sizeof(bool) * g->e_cnt);
    memset(filter_data.delete_edge, false, sizeof(bool) * g->e_cnt);

    filter_data.degree = (uint32_t *)malloc(sizeof(uint32_t) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        filter_data.degree[i] = g->vertex[i + 1] - g->vertex[i];
    }

    uint64_t discard_vertex = 0;

    // 标记顶点和边，以及记录数据
    uint64_t num_of_different_weight_of_edges[recur_depth] = {0};
    uint64_t num_of_different_weight_of_vertex[recur_depth] = {0};
    for (uint64_t i = 0; i < g->v_cnt; ++i) {
        if (filter_data.degree[i] >= s->order_degree[0]) {
            int record_vertex_color = find_range(threshhold, filter_data.degree[i], recur_depth + 1);
            filter_data.vertex_record[i] = record_vertex_color;
            num_of_different_weight_of_vertex[record_vertex_color]++;
            for (uint64_t j = g->vertex[i]; j < g->vertex[i + 1]; ++j) {
                uint32_t temp_degree = filter_data.degree[g->edge[j]];
                if (temp_degree >= s->order_degree[1]) {
                    int record_edge_color = find_range(threshhold, temp_degree, recur_depth + 1);
                    filter_data.edge_record[j] = record_edge_color;
                    num_of_different_weight_of_edges[record_edge_color]++;
                }
            }
        } else {
            discard_vertex++;
        }
    }

    // 投射顶点，记录所有的顶点
    uint64_t prefix_sum_of_vertex[recur_depth + 1] = {0};
    for (int i = 1; i < recur_depth + 1; i++) {
        prefix_sum_of_vertex[i] = prefix_sum_of_vertex[i - 1] + num_of_different_weight_of_vertex[i - 1];
    }
    uint64_t weight_vtx_ptr[recur_depth] = {0};
    for (int i = 1; i < recur_depth; i++) {
        weight_vtx_ptr[i] = prefix_sum_of_vertex[i];
    }
    for (uint32_t i = 0; i < g->v_cnt; ++i) {
        if (filter_data.vertex_record[i] != -1)
            filter_data.vertex_project[weight_vtx_ptr[filter_data.vertex_record[i]]++] = i;
    }

    // 投射边
    // 前缀和形式
    uint64_t prefix_sum_of_edges[recur_depth + 1] = {0};
    for (int i = 1; i < recur_depth + 1; i++) {
        prefix_sum_of_edges[i] = prefix_sum_of_edges[i - 1] + num_of_different_weight_of_edges[i - 1];
    }
    // 利用前缀和更新指针
    uint64_t weight_edge_ptr[recur_depth] = {0};
    for (int i = 0; i < recur_depth; i++) {
        weight_edge_ptr[i] = prefix_sum_of_edges[i];
    }
    for (uint32_t i = 0; i < g->e_cnt; ++i) {
        if (filter_data.edge_record[i] != -1)
            filter_data.edge_project[weight_edge_ptr[filter_data.edge_record[i]]++] = i;
    }

    // 记录采样结果
    filter_data.sample_times_per_edge = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
    filter_data.sample_count_per_edge = (double *)malloc(sizeof(double) * g->e_cnt);
    filter_data.sample_square_per_edge = (double *)malloc(sizeof(double) * g->e_cnt);
    memset(filter_data.sample_times_per_edge, 0, sizeof(uint32_t) * g->e_cnt);
    memset(filter_data.sample_count_per_edge, 0, sizeof(double) * g->e_cnt);
    memset(filter_data.sample_square_per_edge, 0, sizeof(double) * g->e_cnt);

    // 采样批量
    uint64_t batch_size = 1000;

    // 随机数生成器
    std::random_device rd;
    std::default_random_engine gen(rd());

    // 第一次采样的范围 ，*****后续要实时更新
    uint64_t range = prefix_sum_of_edges[recur_depth];
    // 采样记录
    uint64_t approx_sample_times = 0;
    uint64_t real_sample_times = 0;
    double approx_count = 0;
    double approx_squra = 0;
    uint64_t mask = 1ULL;
    mask <<= 63;
    bool flag[4] = {true, true, true, true};

    // 主线程判断是否满足要求使用到的变量，只有主线程会使用到

    double var_local[recur_depth] = {0};
    double std_dev_local[recur_depth] = {0};
    double estimate_error_local[recur_depth] = {0};
    uint64_t approx_sample_times_local[recur_depth] = {0};
    double approx_count_local[recur_depth] = {0};
    double approx_squra_local[recur_depth] = {0};
    double real_error_gloabl = 0.0;
    double real_error = 1.0;

    // 线程工作范围
    int num_of_thread = omp_get_max_threads();
    uint64_t start_edge_id[num_of_thread] = {0};
    uint64_t end_edge_id[num_of_thread];
    for (int i = 0; i < num_of_thread; i++) {
        end_edge_id[i] = prefix_sum_of_edges[recur_depth];
    }

    // 精确计算相关变量
    uint32_t processing_round = 0;
    uint32_t print_times = 0;
    bool exact_computing = true;
    bool get_ready[num_of_thread] = {false};
    uint64_t vertex_range_ptr = 0;
    uint64_t vertex_range_end = 0;
    uint64_t exact_count = 0;
    uint64_t vertex_chunk_per_task = 100;

    // 收敛速率检验相关
    double pre_real_error = 1.0;        // 记录前一个误差
    double pre_approx_sample_times = 0; // 记录前一个采样次数
    double coverge_range = 5000.0;      // 检查的范围差距根据处理速率不断调整。
    double coverge_range_ratio = 1.0;

    // 分层的总结
    double weight_average = 0.0;
    double weight_squra = 0.0;
    double weight_std_error = 0.0;
    double weight_real_error = 0.0;

    double error_threshhold[5] = {0.0};
    double var_threshhold[recur_depth] = {0.0};
    double var_threshhold_pre[recur_depth] = {0.0};

    double multiplicity = schedule_iep.get_multiplicity();
    double iep_redun = schedule_iep.get_in_exclusion_optimize_redundancy();

    double var = 0.0;
    double real_average = 0.0;

    // 误差 边界
#pragma omp parallel
    {
// 可以使用master来实现
#pragma omp master
        {
            while (real_error > 0.1 && approx_sample_times <= 100000000) {

                // 1，精确计算百分之五十后不再继续精确计算。
                if (approx_sample_times > 1000 * print_times) {

                    // 1，更新精确计算范围。 2，打开精确计算开关。3，可能存在前一个精确计算还没有结束就进入精确计算的情况。
                    if (exact_computing) {
                        for (int i = 1; i < num_of_thread; i++) {
                            if (get_ready[i] == true) {
                                i--;
                            }
                        }

                        // 分配精确计算的范围
                        vertex_range_ptr = prefix_sum_of_vertex[processing_round];     // 顶点处理指针
                        vertex_range_end = prefix_sum_of_vertex[processing_round + 1]; // 顶点处理范围
                        // 打开所有线程的精确计算开关
                        for (int i = 1; i < num_of_thread; i++) {
                            get_ready[i] = true;
                        }

                        // 查看所有精确计算是否完成任务
                        if (processing_round + 1 == 1) {
                            for (int i = 1; i < num_of_thread; i++) {
                                if (get_ready[i] == true) {
                                    i--;
                                }
                            }
                        }

                        // 记录现有采样范围，以及更新边范围相关数据
                        uint64_t size = prefix_sum_of_edges[recur_depth];
                        for (int i = 0; i < recur_depth; i++) {
                            num_of_different_weight_of_edges[i] = 0;
                            prefix_sum_of_edges[i] = 0;
                        }

                        uint64_t temp_ptr = 0;
                        int pre_edge_color = 0;
                        int cur_edge_color = 0;
                        int cur_vtx_color = 0;
                        for (int i = 0; i < size; i++) {
                            uint64_t temp_edge = filter_data.edge_project[i];
                            if (filter_data.delete_edge[temp_edge]) {
                                continue;
                            }
                            uint64_t temp_degree = filter_data.degree[g->edge_from[temp_edge]];
                            cur_vtx_color = find_range(threshhold, temp_degree, recur_depth + 1);
                            // 起点度数够大，保留
                            if (cur_vtx_color > processing_round) {
                                // 记录边并且获取边的颜色
                                filter_data.edge_project[temp_ptr++] = temp_edge;
                                cur_edge_color = find_range(threshhold, filter_data.degree[g->edge[temp_edge]], recur_depth + 1);

                                // 更新相应颜色的范围
                                if (cur_edge_color != pre_edge_color) {
                                    pre_edge_color = cur_edge_color;
                                    num_of_different_weight_of_edges[cur_edge_color]++;
                                    prefix_sum_of_edges[cur_edge_color] = temp_ptr - 1;

                                } else {
                                    num_of_different_weight_of_edges[pre_edge_color]++;
                                }
                            }
                        }
                        prefix_sum_of_edges[recur_depth] = temp_ptr;
                        processing_round++;
                    }

                    //==================================================================================================================================================//
                    //=================误差统计=========================================================================================================================//
                    //==================================================================================================================================================//
                    /* 局部误差 */ // 记录平方和，计数和，局部采样次数
                    for (int i = 0; i < recur_depth; i++) {
                        approx_sample_times_local[i] = 0.0;
                        approx_count_local[i] = 0.0;
                        approx_squra_local[i] = 0.0;
                    }
                    for (int i = 0; i < recur_depth; i++) {
                        for (int ptr = prefix_sum_of_edges[i]; ptr < prefix_sum_of_edges[i + 1]; ptr++) {
                            auto eid = filter_data.edge_project[ptr];
                            approx_sample_times_local[i] += filter_data.sample_times_per_edge[eid];
                            approx_count_local[i] += filter_data.sample_count_per_edge[eid];
                            approx_squra_local[i] += filter_data.sample_square_per_edge[eid];
                        }
                    }
                    // 局部均值的偏差
                    for (int i = 0; i < recur_depth; i++) {
                        var_local[i] = (double)(approx_squra_local[i] / approx_sample_times_local[i] -
                                                std::pow((approx_count_local[i] / approx_sample_times_local[i]), 2));
                    }

                    // 局部总和的标准偏差
                    // for (int i = 0; i < recur_depth; i++) {
                    //     double cur_errror = std::sqrt(var_local[i] / approx_sample_times_local[i]) * num_of_different_weight_of_edges[i];
                    //     std::cout << " 线程" << i << "偏差 " << cur_errror << std::endl;
                    // }

                    /* 统计误差 */
                    // 总体方差等于局部方差除以采样次数之和， 总体均值为局部均值之和
                    var = 0.0;
                    real_average = 0.0;
                    for (int i = 0; i < recur_depth; i++) {
                        var +=
                            (var_local[i] * num_of_different_weight_of_edges[i] * num_of_different_weight_of_edges[i]) / approx_sample_times_local[i];
                        real_average += (approx_count_local[i] / approx_sample_times_local[i] * num_of_different_weight_of_edges[i]);
                    }
                    real_average += (double)exact_count * multiplicity / iep_redun;
                    double std_dev = std::sqrt(var);
                    real_error = std_dev * 2.576 / real_average; // 99%的信心

                    // double realreal_error = (real_average - (7375094981 * 2)) / (7375094981 * 2); // 7375094981 5007311619323
                    // double realreal_error = (real_average - (13283084785.0 * multiplicity)) / (13283084785.0 * multiplicity);

                    // printf("approx_count: %lf, approx_sample_times: %lu, estimate_error: %lf\n", real_average, approx_sample_times, real_error);

                    if (flag[0] && real_error < 0.1) {
                        printf("approx_count: %lf, approx_sample_times: %lu, real_error: %lf\n", real_average, approx_sample_times, real_error);
                        flag[0] = false;
                    }
                    if (flag[1] && real_error < 0.05) {
                        printf("approx_count: %lf, approx_sample_times: %lu, real_error: %lf\n", real_average, approx_sample_times, real_error);
                        flag[1] = false;
                    }

                    print_times++;

                    //==================================================================================================================================================//
                    //=========最后调整工作线程的工作范围==================================================================================================================================//
                    //==================================================================================================================================================//

                    /* 自适应的修改采样范围 */
                    // 记录 偏差阈值 迭代上升。
                    for (int i = 0; i < recur_depth; i++) {
                        var_threshhold[i] = std::sqrt(var_local[i] / approx_sample_times_local[i]) * num_of_different_weight_of_edges[i] /
                                            (approx_count_local[i] / approx_sample_times_local[i]);
                        var_threshhold_pre[i] = var_threshhold[i];
                    }
                    std::sort(var_threshhold, var_threshhold + recur_depth);
                    // for (int i = 0; i < 5; i++) {
                    //     int var_threshhold_index = recur_depth - std::pow(2, (i + 1)) + 1;
                    //     if (var_threshhold_index >= 0)
                    //         error_threshhold[i] = var_threshhold[var_threshhold_index] - 1;
                    // }

                    int remain_thread = num_of_thread;
                    int allocate_thread_ptr = 0;
                    int allocate_thread = 0;
                    int allocate_threshhold = 0;
                    int adaptive_start = 0;
                    while (remain_thread > num_of_thread / 2) {
                        allocate_thread = remain_thread / 4;
                        // 查找当前应该分配的阈值的下标。
                        for (int i = 0; i < recur_depth; i++) {
                            if (var_threshhold_pre[i] > error_threshhold[allocate_threshhold]) {
                                adaptive_start = i;
                            }
                        }
                        allocate_threshhold++;
                        for (int i = allocate_thread_ptr; i < allocate_thread_ptr + allocate_thread; i++) {
                            start_edge_id[i] = prefix_sum_of_edges[adaptive_start];
                            end_edge_id[i] = prefix_sum_of_edges[recur_depth];
                        }
                        remain_thread -= allocate_thread;
                        allocate_thread_ptr += allocate_thread;
                    }
                    for (int i = allocate_thread_ptr; i < num_of_thread; i++) {
                        end_edge_id[i] = prefix_sum_of_edges[recur_depth];
                    }
                }
            }
            if (flag[2] && real_error <= 0.01) {
                printf("approx_count: %lf, approx_sample_times: %lu, real_error: %lf\n", real_average, approx_sample_times, real_error);
                flag[2] = false;
            }
        }
        if (omp_get_thread_num()) {
            while (real_error > 0.1 && approx_sample_times <= 100000000) {
                int cur_thread = omp_get_thread_num();
                // 如果需要精确计算
                if (get_ready[cur_thread]) {
                    while (vertex_range_ptr < vertex_range_end) {
                        // 获取当前线程的处理范围
                        uint32_t local_range_start = __sync_fetch_and_add(&vertex_range_ptr, vertex_chunk_per_task);
                        if (local_range_start < vertex_range_end) {
                            uint32_t local_range_end;
                            if (local_range_start + vertex_chunk_per_task < vertex_range_end)
                                local_range_end = local_range_start + vertex_chunk_per_task;
                            else
                                local_range_end = vertex_range_end;

                            /* 精确计算阶段 */
                            int *ans_buffer = new int[schedule_iep.in_exclusion_optimize_vertex_id.size()];  // 底层无交集顶点的大小
                            VertexSet *vertex_set = new VertexSet[schedule_iep.get_total_prefix_num() + 10]; // 构建了多个顶点集合作为候选集
                            VertexSet subtraction_set;                                                       // 构建一个差集
                            VertexSet tmp_set;                                                               // 构建一个临时集合
                            subtraction_set.init();
                            long long local_ans = 0;
                            for (int i = local_range_start; i < local_range_end; ++i) {
                                uint32_t vertex = filter_data.vertex_project[i];
                                e_index_t l, r;
                                l = g->vertex[vertex];
                                r = g->vertex[vertex + 1];

                                // 为调度中的每一个顶点初始化候选顶点集合。
                                for (int prefix_id = schedule_iep.get_last(0); prefix_id != -1; prefix_id = schedule_iep.get_next(prefix_id)) {
                                    vertex_set[prefix_id].build_vertex_set(schedule_iep, vertex_set, &g->edge[l], (int)(r - l), prefix_id);
                                }
                                subtraction_set.push_back(vertex);
                                pattern_matching_aggressive_func(schedule_iep, vertex_set, subtraction_set, tmp_set, local_ans, 1, ans_buffer, g);
                                subtraction_set.pop_back();
                            }
                            delete[] vertex_set;
                            delete[] ans_buffer;

#pragma omp atomic
                            exact_count += local_ans;
                        }
                    }

                    get_ready[cur_thread] = false;
                    exact_computing = false;
                }

                std::uniform_int_distribution<uint64_t> dist(start_edge_id[cur_thread],
                                                             end_edge_id[cur_thread] - 1); //***采样范围的更新可以不需要原子操作
                for (uint64_t i = 0; i < batch_size; i++) {
                    auto random_id = dist(gen);
                    auto eid = filter_data.edge_project[random_id]; // 一开始就处理最后一个区间，但是会不断过滤
                    pattern_sample_ns_record_filter(eid, g, s, &filter_data);
                    // pattern_sample_ns_record(eid, g, s, &filter_data);
                }

#pragma omp atomic
                approx_sample_times += batch_size;
            }
        }
    }

    free(filter_data.edge_project);
    free(filter_data.edge_record);
    free(filter_data.vertex_project);
    free(filter_data.vertex_record);
    free(filter_data.degree);

    free(filter_data.sample_count_per_edge);
    free(filter_data.sample_times_per_edge);
    free(filter_data.sample_square_per_edge);
    free(filter_data.delete_edge);
}

void approximate_exact_fusion_seventyfive(Graph *g, schedule_approx *s, const Schedule_IEP &schedule_iep, Pattern p) {

    edge_from_init(g);

    // 平均度数  frendster 28 twitter 35 livejouranl 17 youtube 2 patent 5 pokec 30 tp 28 sk 11
    // 阈值划分10次 只保留千分之一的顶点 结果波动还很大，
    int recur_depth = 50;
    int vertex_depth = 0;

    uint32_t threshhold[recur_depth + 1] = {0};
    int *record_degree = (int *)malloc(sizeof(int) * g->v_cnt);
    memset(record_degree, 0, sizeof(int) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        record_degree[g->vertex[i + 1] - g->vertex[i]]++;
    }
    double ratio = 0;
    double ratio_threshhold = 0.5;
    double ratio_acum = 0.5;
    // threshhold[0] = 1;
    int iter_threshhoud = 1;
    uint64_t iter_record = 0;
    double change_ratio = 0.0;

    while (iter_record < g->v_cnt) {
        double temp_ratio = double(record_degree[iter_record]) / double(g->v_cnt);
        ratio += temp_ratio;
        iter_record++;
        if (ratio > ratio_threshhold) {
            threshhold[iter_threshhoud++] = iter_record;
            // change_ratio = double(threshhold[iter_threshhoud - 1]) / double(threshhold[iter_threshhoud - 2]);
            ratio_acum = ratio_acum * 0.5;
            ratio_threshhold = ratio_threshhold + ratio_acum;
            if (iter_threshhoud == recur_depth || iter_record > 20000 || iter_threshhoud > 15) { //|| change_ratio < 1.1
                break;
            }
        }
    }
    for (int i = 1; i < recur_depth; i++) {
        if (threshhold[i] == 0) {
            recur_depth = i;
            break;
        }
    }
    threshhold[recur_depth] = g->v_cnt;
    free(record_degree);

    // 一共有recur_depth中颜色， 但是查找阈值的时候需要用recur_depth +1.因为阈值第一位为0

    // 查看阈值
    // for (int i = 0; i < recur_depth; i++) {
    //     printf("threshhold[%d] = %d\n", i, threshhold[i]);
    // }

    // 过滤的数据结构。
    filter filter_data;
    filter_data.edge_project = (int64_t *)malloc(sizeof(int64_t) * g->e_cnt);
    filter_data.edge_record = (int *)malloc(sizeof(int) * g->e_cnt);
    filter_data.vertex_project = (int *)malloc(sizeof(int) * g->v_cnt);
    filter_data.vertex_record = (int *)malloc(sizeof(int) * g->v_cnt);
    memset(filter_data.edge_project, -1, sizeof(int64_t) * g->e_cnt);
    memset(filter_data.edge_record, -1, sizeof(int) * g->e_cnt);
    memset(filter_data.vertex_record, -1, sizeof(int) * g->v_cnt);
    memset(filter_data.vertex_project, -1, sizeof(int) * g->v_cnt);
    filter_data.delete_edge = (bool *)malloc(sizeof(bool) * g->e_cnt);
    memset(filter_data.delete_edge, false, sizeof(bool) * g->e_cnt);

    filter_data.degree = (uint32_t *)malloc(sizeof(uint32_t) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        filter_data.degree[i] = g->vertex[i + 1] - g->vertex[i];
    }

    uint64_t discard_vertex = 0;

    // 标记顶点和边，以及记录数据
    uint64_t num_of_different_weight_of_edges[recur_depth] = {0};
    uint64_t num_of_different_weight_of_vertex[recur_depth] = {0};
    for (uint64_t i = 0; i < g->v_cnt; ++i) {
        if (filter_data.degree[i] >= s->order_degree[0]) {
            int record_vertex_color = find_range(threshhold, filter_data.degree[i], recur_depth + 1);
            filter_data.vertex_record[i] = record_vertex_color;
            num_of_different_weight_of_vertex[record_vertex_color]++;
            for (uint64_t j = g->vertex[i]; j < g->vertex[i + 1]; ++j) {
                uint32_t temp_degree = filter_data.degree[g->edge[j]];
                if (temp_degree >= s->order_degree[1]) {
                    int record_edge_color = find_range(threshhold, temp_degree, recur_depth + 1);
                    filter_data.edge_record[j] = record_edge_color;
                    num_of_different_weight_of_edges[record_edge_color]++;
                }
            }
        } else {
            discard_vertex++;
        }
    }

    // 投射顶点，记录所有的顶点
    uint64_t prefix_sum_of_vertex[recur_depth + 1] = {0};
    for (int i = 1; i < recur_depth + 1; i++) {
        prefix_sum_of_vertex[i] = prefix_sum_of_vertex[i - 1] + num_of_different_weight_of_vertex[i - 1];
    }
    uint64_t weight_vtx_ptr[recur_depth] = {0};
    for (int i = 1; i < recur_depth; i++) {
        weight_vtx_ptr[i] = prefix_sum_of_vertex[i];
    }
    for (uint32_t i = 0; i < g->v_cnt; ++i) {
        if (filter_data.vertex_record[i] != -1)
            filter_data.vertex_project[weight_vtx_ptr[filter_data.vertex_record[i]]++] = i;
    }

    // 投射边
    // 前缀和形式
    uint64_t prefix_sum_of_edges[recur_depth + 1] = {0};
    for (int i = 1; i < recur_depth + 1; i++) {
        prefix_sum_of_edges[i] = prefix_sum_of_edges[i - 1] + num_of_different_weight_of_edges[i - 1];
    }
    // 利用前缀和更新指针
    uint64_t weight_edge_ptr[recur_depth] = {0};
    for (int i = 0; i < recur_depth; i++) {
        weight_edge_ptr[i] = prefix_sum_of_edges[i];
    }
    for (uint32_t i = 0; i < g->e_cnt; ++i) {
        if (filter_data.edge_record[i] != -1)
            filter_data.edge_project[weight_edge_ptr[filter_data.edge_record[i]]++] = i;
    }

    // 记录采样结果
    filter_data.sample_times_per_edge = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
    filter_data.sample_count_per_edge = (double *)malloc(sizeof(double) * g->e_cnt);
    filter_data.sample_square_per_edge = (double *)malloc(sizeof(double) * g->e_cnt);
    memset(filter_data.sample_times_per_edge, 0, sizeof(uint32_t) * g->e_cnt);
    memset(filter_data.sample_count_per_edge, 0, sizeof(double) * g->e_cnt);
    memset(filter_data.sample_square_per_edge, 0, sizeof(double) * g->e_cnt);

    // 采样批量
    uint64_t batch_size = 1000;

    // 随机数生成器
    std::random_device rd;
    std::default_random_engine gen(rd());

    // 第一次采样的范围 ，*****后续要实时更新
    uint64_t range = prefix_sum_of_edges[recur_depth];
    // 采样记录
    uint64_t approx_sample_times = 0;
    uint64_t real_sample_times = 0;
    double approx_count = 0;
    double approx_squra = 0;
    uint64_t mask = 1ULL;
    mask <<= 63;
    bool flag[4] = {true, true, true, true};

    // 主线程判断是否满足要求使用到的变量，只有主线程会使用到

    double var_local[recur_depth] = {0};
    double std_dev_local[recur_depth] = {0};
    double estimate_error_local[recur_depth] = {0};
    uint64_t approx_sample_times_local[recur_depth] = {0};
    double approx_count_local[recur_depth] = {0};
    double approx_squra_local[recur_depth] = {0};
    double real_error_gloabl = 0.0;
    double real_error = 1.0;

    // 线程工作范围
    int num_of_thread = omp_get_max_threads();
    uint64_t start_edge_id[num_of_thread] = {0};
    uint64_t end_edge_id[num_of_thread];
    for (int i = 0; i < num_of_thread; i++) {
        end_edge_id[i] = prefix_sum_of_edges[recur_depth];
    }

    // 精确计算相关变量
    uint32_t processing_round = 0;
    uint32_t print_times = 0;
    bool exact_computing = true;
    bool get_ready[num_of_thread] = {false};
    uint64_t vertex_range_ptr = 0;
    uint64_t vertex_range_end = 0;
    uint64_t exact_count = 0;
    uint64_t vertex_chunk_per_task = 100;

    // 收敛速率检验相关
    double pre_real_error = 1.0;        // 记录前一个误差
    double pre_approx_sample_times = 0; // 记录前一个采样次数
    double coverge_range = 5000.0;      // 检查的范围差距根据处理速率不断调整。
    double coverge_range_ratio = 1.0;

    // 分层的总结
    double weight_average = 0.0;
    double weight_squra = 0.0;
    double weight_std_error = 0.0;
    double weight_real_error = 0.0;

    double error_threshhold[5] = {0.0};
    double var_threshhold[recur_depth] = {0.0};
    double var_threshhold_pre[recur_depth] = {0.0};

    double multiplicity = schedule_iep.get_multiplicity();
    double iep_redun = schedule_iep.get_in_exclusion_optimize_redundancy();

    double var = 0.0;
    double real_average = 0.0;

    // 误差 边界
#pragma omp parallel
    {
// 可以使用master来实现
#pragma omp master
        {
            while (real_error > 0.1 && approx_sample_times <= 100000000) {

                // 1，精确计算百分之五十后不再继续精确计算。
                if (approx_sample_times > 1000 * print_times) {

                    // 1，更新精确计算范围。 2，打开精确计算开关。3，可能存在前一个精确计算还没有结束就进入精确计算的情况。
                    if (exact_computing) {
                        for (int i = 1; i < num_of_thread; i++) {
                            if (get_ready[i] == true) {
                                i--;
                            }
                        }

                        // 分配精确计算的范围
                        vertex_range_ptr = prefix_sum_of_vertex[0]; // 顶点处理指针
                        vertex_range_end = prefix_sum_of_vertex[2]; // 顶点处理范围
                        processing_round++;
                        // 打开所有线程的精确计算开关
                        for (int i = 1; i < num_of_thread; i++) {
                            get_ready[i] = true;
                        }

                        // 查看所有精确计算是否完成任务
                        if (processing_round + 1 == 2) {
                            for (int i = 1; i < num_of_thread; i++) {
                                if (get_ready[i] == true) {
                                    i--;
                                }
                            }
                        }

                        // 记录现有采样范围，以及更新边范围相关数据
                        uint64_t size = prefix_sum_of_edges[recur_depth];
                        for (int i = 0; i < recur_depth; i++) {
                            num_of_different_weight_of_edges[i] = 0;
                            prefix_sum_of_edges[i] = 0;
                        }

                        uint64_t temp_ptr = 0;
                        int pre_edge_color = 0;
                        int cur_edge_color = 0;
                        int cur_vtx_color = 0;
                        for (int i = 0; i < size; i++) {
                            uint64_t temp_edge = filter_data.edge_project[i];
                            if (filter_data.delete_edge[temp_edge]) {
                                continue;
                            }
                            uint64_t temp_degree = filter_data.degree[g->edge_from[temp_edge]];
                            cur_vtx_color = find_range(threshhold, temp_degree, recur_depth + 1);
                            // 起点度数够大，保留
                            if (cur_vtx_color > processing_round) {
                                // 记录边并且获取边的颜色
                                filter_data.edge_project[temp_ptr++] = temp_edge;
                                cur_edge_color = find_range(threshhold, filter_data.degree[g->edge[temp_edge]], recur_depth + 1);

                                // 更新相应颜色的范围
                                if (cur_edge_color != pre_edge_color) {
                                    pre_edge_color = cur_edge_color;
                                    num_of_different_weight_of_edges[cur_edge_color]++;
                                    prefix_sum_of_edges[cur_edge_color] = temp_ptr - 1;

                                } else {
                                    num_of_different_weight_of_edges[pre_edge_color]++;
                                }
                            }
                        }
                        prefix_sum_of_edges[recur_depth] = temp_ptr;
                        processing_round++;
                    }

                    //==================================================================================================================================================//
                    //=================误差统计=========================================================================================================================//
                    //==================================================================================================================================================//
                    /* 局部误差 */ // 记录平方和，计数和，局部采样次数
                    for (int i = 0; i < recur_depth; i++) {
                        approx_sample_times_local[i] = 0.0;
                        approx_count_local[i] = 0.0;
                        approx_squra_local[i] = 0.0;
                    }
                    for (int i = 0; i < recur_depth; i++) {
                        for (int ptr = prefix_sum_of_edges[i]; ptr < prefix_sum_of_edges[i + 1]; ptr++) {
                            auto eid = filter_data.edge_project[ptr];
                            approx_sample_times_local[i] += filter_data.sample_times_per_edge[eid];
                            approx_count_local[i] += filter_data.sample_count_per_edge[eid];
                            approx_squra_local[i] += filter_data.sample_square_per_edge[eid];
                        }
                    }
                    // 局部均值的偏差
                    for (int i = 0; i < recur_depth; i++) {
                        var_local[i] = (double)(approx_squra_local[i] / approx_sample_times_local[i] -
                                                std::pow((approx_count_local[i] / approx_sample_times_local[i]), 2));
                    }

                    // 局部总和的标准偏差
                    // for (int i = 0; i < recur_depth; i++) {
                    //     double cur_errror = std::sqrt(var_local[i] / approx_sample_times_local[i]) * num_of_different_weight_of_edges[i];
                    //     std::cout << " 线程" << i << "偏差 " << cur_errror << std::endl;
                    // }

                    /* 统计误差 */
                    // 总体方差等于局部方差除以采样次数之和， 总体均值为局部均值之和
                    var = 0.0;
                    real_average = 0.0;
                    for (int i = 0; i < recur_depth; i++) {
                        var +=
                            (var_local[i] * num_of_different_weight_of_edges[i] * num_of_different_weight_of_edges[i]) / approx_sample_times_local[i];
                        real_average += (approx_count_local[i] / approx_sample_times_local[i] * num_of_different_weight_of_edges[i]);
                    }
                    real_average += (double)exact_count * multiplicity / iep_redun;
                    double std_dev = std::sqrt(var);
                    real_error = std_dev * 2.576 / real_average; // 99%的信心

                    // double realreal_error = (real_average - (7375094981 * 2)) / (7375094981 * 2); // 7375094981 5007311619323
                    // double realreal_error = (real_average - (13283084785.0 * multiplicity)) / (13283084785.0 * multiplicity);

                    // printf("approx_count: %lu, real_error: %lf, estimate_error: %lf\n", approx_sample_times, realreal_error, real_error);

                    if (flag[0] && real_error < 0.1) {
                        printf("approx_count: %lf, approx_sample_times: %lu, real_error: %lf\n", real_average, approx_sample_times, real_error);
                        flag[0] = false;
                    }
                    if (flag[1] && real_error < 0.05) {
                        printf("approx_count: %lf, approx_sample_times: %lu, real_error: %lf\n", real_average, approx_sample_times, real_error);
                        flag[1] = false;
                    }

                    print_times++;

                    //==================================================================================================================================================//
                    //=========最后调整工作线程的工作范围==================================================================================================================================//
                    //==================================================================================================================================================//

                    /* 自适应的修改采样范围 */
                    // 记录 偏差阈值 迭代上升。
                    for (int i = 0; i < recur_depth; i++) {
                        var_threshhold[i] = std::sqrt(var_local[i] / approx_sample_times_local[i]) * num_of_different_weight_of_edges[i] /
                                            (approx_count_local[i] / approx_sample_times_local[i]);
                        var_threshhold_pre[i] = var_threshhold[i];
                    }
                    std::sort(var_threshhold, var_threshhold + recur_depth);
                    // for (int i = 0; i < 5; i++) {
                    //     int var_threshhold_index = recur_depth - std::pow(2, (i + 1)) + 1;
                    //     if (var_threshhold_index >= 0)
                    //         error_threshhold[i] = var_threshhold[var_threshhold_index] - 1;
                    // }

                    int remain_thread = num_of_thread;
                    int allocate_thread_ptr = 0;
                    int allocate_thread = 0;
                    int allocate_threshhold = 0;
                    int adaptive_start = 0;
                    while (remain_thread > num_of_thread / 2) {
                        allocate_thread = remain_thread / 4;
                        // 查找当前应该分配的阈值的下标。
                        for (int i = 0; i < recur_depth; i++) {
                            if (var_threshhold_pre[i] > error_threshhold[allocate_threshhold]) {
                                adaptive_start = i;
                            }
                        }
                        allocate_threshhold++;
                        for (int i = allocate_thread_ptr; i < allocate_thread_ptr + allocate_thread; i++) {
                            start_edge_id[i] = prefix_sum_of_edges[adaptive_start];
                            end_edge_id[i] = prefix_sum_of_edges[recur_depth];
                        }
                        remain_thread -= allocate_thread;
                        allocate_thread_ptr += allocate_thread;
                    }
                    for (int i = allocate_thread_ptr; i < num_of_thread; i++) {
                        end_edge_id[i] = prefix_sum_of_edges[recur_depth];
                    }
                }
            }
            if (flag[2] && real_error <= 0.01) {
                printf("approx_count: %lf, approx_sample_times: %lu, real_error: %lf\n", real_average, approx_sample_times, real_error);
                flag[2] = false;
            }
        }
        if (omp_get_thread_num()) {
            while (real_error > 0.1 && approx_sample_times <= 100000000) {
                int cur_thread = omp_get_thread_num();
                // 如果需要精确计算
                if (get_ready[cur_thread]) {
                    while (vertex_range_ptr < vertex_range_end) {
                        // 获取当前线程的处理范围
                        uint32_t local_range_start = __sync_fetch_and_add(&vertex_range_ptr, vertex_chunk_per_task);
                        if (local_range_start < vertex_range_end) {
                            uint32_t local_range_end;
                            if (local_range_start + vertex_chunk_per_task < vertex_range_end)
                                local_range_end = local_range_start + vertex_chunk_per_task;
                            else
                                local_range_end = vertex_range_end;

                            /* 精确计算阶段 */
                            int *ans_buffer = new int[schedule_iep.in_exclusion_optimize_vertex_id.size()];  // 底层无交集顶点的大小
                            VertexSet *vertex_set = new VertexSet[schedule_iep.get_total_prefix_num() + 10]; // 构建了多个顶点集合作为候选集
                            VertexSet subtraction_set;                                                       // 构建一个差集
                            VertexSet tmp_set;                                                               // 构建一个临时集合
                            subtraction_set.init();
                            long long local_ans = 0;
                            for (int i = local_range_start; i < local_range_end; ++i) {
                                uint32_t vertex = filter_data.vertex_project[i];
                                e_index_t l, r;
                                l = g->vertex[vertex];
                                r = g->vertex[vertex + 1];

                                // 为调度中的每一个顶点初始化候选顶点集合。
                                for (int prefix_id = schedule_iep.get_last(0); prefix_id != -1; prefix_id = schedule_iep.get_next(prefix_id)) {
                                    vertex_set[prefix_id].build_vertex_set(schedule_iep, vertex_set, &g->edge[l], (int)(r - l), prefix_id);
                                }
                                subtraction_set.push_back(vertex);
                                pattern_matching_aggressive_func(schedule_iep, vertex_set, subtraction_set, tmp_set, local_ans, 1, ans_buffer, g);
                                subtraction_set.pop_back();
                            }
                            delete[] vertex_set;
                            delete[] ans_buffer;

#pragma omp atomic
                            exact_count += local_ans;
                        }
                    }

                    get_ready[cur_thread] = false;
                    exact_computing = false;
                }

                std::uniform_int_distribution<uint64_t> dist(start_edge_id[cur_thread],
                                                             end_edge_id[cur_thread] - 1); //***采样范围的更新可以不需要原子操作
                for (uint64_t i = 0; i < batch_size; i++) {
                    auto random_id = dist(gen);
                    auto eid = filter_data.edge_project[random_id]; // 一开始就处理最后一个区间，但是会不断过滤
                    pattern_sample_ns_record_filter(eid, g, s, &filter_data);
                    // pattern_sample_ns_record(eid, g, s, &filter_data);
                }

#pragma omp atomic
                approx_sample_times += batch_size;
            }
        }
    }

    free(filter_data.edge_project);
    free(filter_data.edge_record);
    free(filter_data.vertex_project);
    free(filter_data.vertex_record);
    free(filter_data.degree);

    free(filter_data.sample_count_per_edge);
    free(filter_data.sample_times_per_edge);
    free(filter_data.sample_square_per_edge);
    free(filter_data.delete_edge);
}

void approximate_exact_fusion_fifty_allrange(Graph *g, schedule_approx *s, const Schedule_IEP &schedule_iep, Pattern p) {

    edge_from_init(g);

    // 平均度数  frendster 28 twitter 35 livejouranl 17 youtube 2 patent 5 pokec 30 tp 28 sk 11
    // 阈值划分10次 只保留千分之一的顶点 结果波动还很大，
    int recur_depth = 50;
    int vertex_depth = 0;

    uint32_t threshhold[recur_depth + 1] = {0};
    int *record_degree = (int *)malloc(sizeof(int) * g->v_cnt);
    memset(record_degree, 0, sizeof(int) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        record_degree[g->vertex[i + 1] - g->vertex[i]]++;
    }
    double ratio = 0;
    double ratio_threshhold = 0.5;
    double ratio_acum = 0.5;
    // threshhold[0] = 1;
    int iter_threshhoud = 1;
    uint64_t iter_record = 0;
    double change_ratio = 0.0;

    while (iter_record < g->v_cnt) {
        double temp_ratio = double(record_degree[iter_record]) / double(g->v_cnt);
        ratio += temp_ratio;
        iter_record++;
        if (ratio > ratio_threshhold) {
            threshhold[iter_threshhoud++] = iter_record;
            // change_ratio = double(threshhold[iter_threshhoud - 1]) / double(threshhold[iter_threshhoud - 2]);
            ratio_acum = ratio_acum * 0.5;
            ratio_threshhold = ratio_threshhold + ratio_acum;
            if (iter_threshhoud == recur_depth || iter_record > 20000 || iter_threshhoud > 18) { //|| change_ratio < 1.1
                break;
            }
        }
    }
    for (int i = 1; i < recur_depth; i++) {
        if (threshhold[i] == 0) {
            recur_depth = i;
            break;
        }
    }
    threshhold[recur_depth] = g->v_cnt;
    free(record_degree);

    // 一共有recur_depth中颜色， 但是查找阈值的时候需要用recur_depth +1.因为阈值第一位为0

    // 过滤的数据结构。
    filter filter_data;
    filter_data.edge_project = (int64_t *)malloc(sizeof(int64_t) * g->e_cnt);
    filter_data.edge_record = (int *)malloc(sizeof(int) * g->e_cnt);
    filter_data.vertex_project = (int *)malloc(sizeof(int) * g->v_cnt);
    filter_data.vertex_record = (int *)malloc(sizeof(int) * g->v_cnt);
    memset(filter_data.edge_project, -1, sizeof(int64_t) * g->e_cnt);
    memset(filter_data.edge_record, -1, sizeof(int) * g->e_cnt);
    memset(filter_data.vertex_record, -1, sizeof(int) * g->v_cnt);
    memset(filter_data.vertex_project, -1, sizeof(int) * g->v_cnt);
    filter_data.delete_edge = (bool *)malloc(sizeof(bool) * g->e_cnt);
    memset(filter_data.delete_edge, false, sizeof(bool) * g->e_cnt);

    filter_data.degree = (uint32_t *)malloc(sizeof(uint32_t) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        filter_data.degree[i] = g->vertex[i + 1] - g->vertex[i];
    }

    uint64_t discard_vertex = 0;

    // 标记顶点和边，以及记录数据
    uint64_t num_of_different_weight_of_edges[recur_depth] = {0};
    uint64_t num_of_different_weight_of_vertex[recur_depth] = {0};
    for (uint64_t i = 0; i < g->v_cnt; ++i) {
        if (filter_data.degree[i] >= s->order_degree[0]) {
            int record_vertex_color = find_range(threshhold, filter_data.degree[i], recur_depth + 1);
            filter_data.vertex_record[i] = record_vertex_color;
            num_of_different_weight_of_vertex[record_vertex_color]++;
            for (uint64_t j = g->vertex[i]; j < g->vertex[i + 1]; ++j) {
                uint32_t temp_degree = filter_data.degree[g->edge[j]];
                if (temp_degree >= s->order_degree[1]) {
                    int record_edge_color = find_range(threshhold, temp_degree, recur_depth + 1);
                    filter_data.edge_record[j] = record_edge_color;
                    num_of_different_weight_of_edges[record_edge_color]++;
                }
            }
        } else {
            discard_vertex++;
        }
    }

    // 投射顶点，记录所有的顶点
    uint64_t prefix_sum_of_vertex[recur_depth + 1] = {0};
    for (int i = 1; i < recur_depth + 1; i++) {
        prefix_sum_of_vertex[i] = prefix_sum_of_vertex[i - 1] + num_of_different_weight_of_vertex[i - 1];
    }
    uint64_t weight_vtx_ptr[recur_depth] = {0};
    for (int i = 1; i < recur_depth; i++) {
        weight_vtx_ptr[i] = prefix_sum_of_vertex[i];
    }
    for (uint32_t i = 0; i < g->v_cnt; ++i) {
        if (filter_data.vertex_record[i] != -1)
            filter_data.vertex_project[weight_vtx_ptr[filter_data.vertex_record[i]]++] = i;
    }

    // 投射边
    // 前缀和形式
    uint64_t prefix_sum_of_edges[recur_depth + 1] = {0};
    for (int i = 1; i < recur_depth + 1; i++) {
        prefix_sum_of_edges[i] = prefix_sum_of_edges[i - 1] + num_of_different_weight_of_edges[i - 1];
    }
    // 利用前缀和更新指针
    uint64_t weight_edge_ptr[recur_depth] = {0};
    for (int i = 0; i < recur_depth; i++) {
        weight_edge_ptr[i] = prefix_sum_of_edges[i];
    }
    for (uint32_t i = 0; i < g->e_cnt; ++i) {
        if (filter_data.edge_record[i] != -1)
            filter_data.edge_project[weight_edge_ptr[filter_data.edge_record[i]]++] = i;
    }

    // 记录采样结果
    filter_data.sample_times_per_edge = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
    filter_data.sample_count_per_edge = (double *)malloc(sizeof(double) * g->e_cnt);
    filter_data.sample_square_per_edge = (double *)malloc(sizeof(double) * g->e_cnt);
    memset(filter_data.sample_times_per_edge, 0, sizeof(uint32_t) * g->e_cnt);
    memset(filter_data.sample_count_per_edge, 0, sizeof(double) * g->e_cnt);
    memset(filter_data.sample_square_per_edge, 0, sizeof(double) * g->e_cnt);

    // 采样批量
    uint64_t batch_size = 1000;

    // 随机数生成器
    std::random_device rd;
    std::default_random_engine gen(rd());

    // 第一次采样的范围 ，*****后续要实时更新
    uint64_t range = prefix_sum_of_edges[recur_depth];
    // 采样记录
    uint64_t approx_sample_times = 0;
    uint64_t real_sample_times = 0;
    double approx_count = 0;
    double approx_squra = 0;
    uint64_t mask = 1ULL;
    mask <<= 63;
    bool flag[4] = {true, true, true, true};

    // 主线程判断是否满足要求使用到的变量，只有主线程会使用到

    double var_local[recur_depth] = {0};
    double std_dev_local[recur_depth] = {0};
    double estimate_error_local[recur_depth] = {0};
    uint64_t approx_sample_times_local[recur_depth] = {0};
    double approx_count_local[recur_depth] = {0};
    double approx_squra_local[recur_depth] = {0};
    double real_error_gloabl = 0.0;
    double real_error = 1.0;

    // 线程工作范围
    int num_of_thread = omp_get_max_threads();
    uint64_t start_edge_id[num_of_thread] = {0};
    uint64_t end_edge_id[num_of_thread];
    for (int i = 0; i < num_of_thread; i++) {
        end_edge_id[i] = prefix_sum_of_edges[recur_depth];
    }

    // 精确计算相关变量
    uint32_t processing_round = 0;
    uint32_t print_times = 0;
    bool exact_computing = true;
    bool get_ready[num_of_thread] = {false};
    uint64_t vertex_range_ptr = 0;
    uint64_t vertex_range_end = 0;
    uint64_t exact_count = 0;
    uint64_t vertex_chunk_per_task = 100;

    // 收敛速率检验相关
    double pre_real_error = 1.0;        // 记录前一个误差
    double pre_approx_sample_times = 0; // 记录前一个采样次数
    double coverge_range = 5000.0;      // 检查的范围差距根据处理速率不断调整。
    double coverge_range_ratio = 1.0;

    // 分层的总结
    double weight_average = 0.0;
    double weight_squra = 0.0;
    double weight_std_error = 0.0;
    double weight_real_error = 0.0;

    double error_threshhold[5] = {0.0};
    double var_threshhold[recur_depth] = {0.0};
    double var_threshhold_pre[recur_depth] = {0.0};

    double multiplicity = schedule_iep.get_multiplicity();
    double iep_redun = schedule_iep.get_in_exclusion_optimize_redundancy();

    double var = 0.0;
    double real_average = 0.0;
    double sample_times = 0.0;

    // 误差 边界
#pragma omp parallel
    {
// 可以使用master来实现
#pragma omp master
        {
            while (real_error > 0.01 && approx_sample_times <= 100000000) {

                // 1，精确计算百分之五十后不再继续精确计算。
                if (approx_sample_times > 1000 * print_times) {

                    // 1，更新精确计算范围。 2，打开精确计算开关。3，可能存在前一个精确计算还没有结束就进入精确计算的情况。
                    if (exact_computing) {
                        for (int i = 1; i < num_of_thread; i++) {
                            if (get_ready[i] == true) {
                                i--;
                            }
                        }

                        // 分配精确计算的范围
                        vertex_range_ptr = prefix_sum_of_vertex[processing_round];     // 顶点处理指针
                        vertex_range_end = prefix_sum_of_vertex[processing_round + 1]; // 顶点处理范围
                        // 打开所有线程的精确计算开关
                        for (int i = 1; i < num_of_thread; i++) {
                            get_ready[i] = true;
                        }

                        // 查看所有精确计算是否完成任务
                        if (processing_round + 1 == 1) {
                            for (int i = 1; i < num_of_thread; i++) {
                                if (get_ready[i] == true) {
                                    i--;
                                }
                            }
                        }

                        // 记录现有采样范围，以及更新边范围相关数据
                        uint64_t size = prefix_sum_of_edges[recur_depth];
                        for (int i = 0; i < recur_depth; i++) {
                            num_of_different_weight_of_edges[i] = 0;
                            prefix_sum_of_edges[i] = 0;
                        }

                        uint64_t temp_ptr = 0;
                        int pre_edge_color = 0;
                        int cur_edge_color = 0;
                        int cur_vtx_color = 0;
                        for (int i = 0; i < size; i++) {
                            uint64_t temp_edge = filter_data.edge_project[i];
                            if (filter_data.delete_edge[temp_edge]) {
                                continue;
                            }
                            uint64_t temp_degree = filter_data.degree[g->edge_from[temp_edge]];
                            cur_vtx_color = find_range(threshhold, temp_degree, recur_depth + 1);
                            // 起点度数够大，保留
                            if (cur_vtx_color > processing_round) {
                                // 记录边并且获取边的颜色
                                filter_data.edge_project[temp_ptr++] = temp_edge;
                                cur_edge_color = find_range(threshhold, filter_data.degree[g->edge[temp_edge]], recur_depth + 1);

                                // 更新相应颜色的范围
                                if (cur_edge_color != pre_edge_color) {
                                    pre_edge_color = cur_edge_color;
                                    num_of_different_weight_of_edges[cur_edge_color]++;
                                    prefix_sum_of_edges[cur_edge_color] = temp_ptr - 1;

                                } else {
                                    num_of_different_weight_of_edges[pre_edge_color]++;
                                }
                            }
                        }
                        prefix_sum_of_edges[recur_depth] = temp_ptr;

                        for (int i = 0; i < num_of_thread; i++) {
                            end_edge_id[i] = prefix_sum_of_edges[recur_depth];
                        }

                        processing_round++;
                    }

                    //==================================================================================================================================================//
                    //=================误差统计=========================================================================================================================//
                    //==================================================================================================================================================//

                    /* 统计误差 */
                    // 总体方差等于局部方差除以采样次数之和， 总体均值为局部均值之和
                    var = 0.0;
                    real_average = 0.0;
                    sample_times = 0.0;
                    for (int i = 0; i < prefix_sum_of_edges[recur_depth]; i++) {
                        auto eid = filter_data.edge_project[i];
                        sample_times += filter_data.sample_times_per_edge[eid];
                        real_average += filter_data.sample_count_per_edge[eid];
                        var += filter_data.sample_square_per_edge[eid];
                    }
                    real_average /= sample_times;
                    real_average *= prefix_sum_of_edges[recur_depth]; // 总体估计

                    var = var / sample_times;
                    var = var * prefix_sum_of_edges[recur_depth] * prefix_sum_of_edges[recur_depth]; // 总体方差
                    var = var - real_average * real_average;
                    var = var / sample_times;
                    double std_dev = std::sqrt(var);

                    real_error = std_dev * 2.576 / real_average;
                    real_average += (double)exact_count * multiplicity / iep_redun;

                    // double realreal_error = (real_average - (7375094981 * 2)) / (7375094981 * 2); // 7375094981 5007311619323
                    // double realreal_error = (real_average - (13283084785.0 * multiplicity)) / (13283084785.0 * multiplicity);

                    // printf("approx_sample_times: %lu, real_error: %lf, estimate_error: %lf\n", approx_sample_times, realreal_error, real_error);

                    if (flag[0] && real_error < 0.1) {
                        printf("approx_count: %lf, approx_sample_times: %lu, real_error: %lf\n", real_average, approx_sample_times, real_error);
                        flag[0] = false;
                    }
                    if (flag[1] && real_error < 0.05) {
                        printf("approx_count: %lf, approx_sample_times: %lu, real_error: %lf\n", real_average, approx_sample_times, real_error);
                        flag[1] = false;
                    }

                    print_times++;
                }
            }
            if (flag[2] && real_error <= 0.01) {
                printf("approx_count: %lf, approx_sample_times: %lu, real_error: %lf\n", real_average, approx_sample_times, real_error);
                flag[2] = false;
            }
        }
        if (omp_get_thread_num()) {
            while (real_error > 0.01 && approx_sample_times <= 100000000) {
                int cur_thread = omp_get_thread_num();
                // 如果需要精确计算
                if (get_ready[cur_thread]) {
                    while (vertex_range_ptr < vertex_range_end) {
                        // 获取当前线程的处理范围
                        uint32_t local_range_start = __sync_fetch_and_add(&vertex_range_ptr, vertex_chunk_per_task);
                        if (local_range_start < vertex_range_end) {
                            uint32_t local_range_end;
                            if (local_range_start + vertex_chunk_per_task < vertex_range_end)
                                local_range_end = local_range_start + vertex_chunk_per_task;
                            else
                                local_range_end = vertex_range_end;

                            /* 精确计算阶段 */
                            int *ans_buffer = new int[schedule_iep.in_exclusion_optimize_vertex_id.size()];  // 底层无交集顶点的大小
                            VertexSet *vertex_set = new VertexSet[schedule_iep.get_total_prefix_num() + 10]; // 构建了多个顶点集合作为候选集
                            VertexSet subtraction_set;                                                       // 构建一个差集
                            VertexSet tmp_set;                                                               // 构建一个临时集合
                            subtraction_set.init();
                            long long local_ans = 0;
                            for (int i = local_range_start; i < local_range_end; ++i) {
                                uint32_t vertex = filter_data.vertex_project[i];
                                e_index_t l, r;
                                l = g->vertex[vertex];
                                r = g->vertex[vertex + 1];

                                // 为调度中的每一个顶点初始化候选顶点集合。
                                for (int prefix_id = schedule_iep.get_last(0); prefix_id != -1; prefix_id = schedule_iep.get_next(prefix_id)) {
                                    vertex_set[prefix_id].build_vertex_set(schedule_iep, vertex_set, &g->edge[l], (int)(r - l), prefix_id);
                                }
                                subtraction_set.push_back(vertex);
                                pattern_matching_aggressive_func(schedule_iep, vertex_set, subtraction_set, tmp_set, local_ans, 1, ans_buffer, g);
                                subtraction_set.pop_back();
                            }
                            delete[] vertex_set;
                            delete[] ans_buffer;

#pragma omp atomic
                            exact_count += local_ans;
                        }
                    }

                    get_ready[cur_thread] = false;
                    exact_computing = false;
                }

                std::uniform_int_distribution<uint64_t> dist(start_edge_id[cur_thread],
                                                             end_edge_id[cur_thread] - 1); //***采样范围的更新可以不需要原子操作
                for (uint64_t i = 0; i < batch_size; i++) {
                    auto random_id = dist(gen);
                    auto eid = filter_data.edge_project[random_id]; // 一开始就处理最后一个区间，但是会不断过滤
                    pattern_sample_ns_record_filter(eid, g, s, &filter_data);
                    // pattern_sample_ns_record(eid, g, s, &filter_data);
                }

#pragma omp atomic
                approx_sample_times += batch_size;
            }
        }
    }

    free(filter_data.edge_project);
    free(filter_data.edge_record);
    free(filter_data.vertex_project);
    free(filter_data.vertex_record);
    free(filter_data.degree);

    free(filter_data.sample_count_per_edge);
    free(filter_data.sample_times_per_edge);
    free(filter_data.sample_square_per_edge);
    free(filter_data.delete_edge);
}

void exact_count_version_two(Graph *g, schedule_approx *s, const Schedule_IEP &schedule_iep, Pattern p) {

    edge_from_init(g);

    // 平均度数  frendster 28 twitter 35 livejouranl 17 youtube 2 patent 5 pokec 30 tp 28 sk 11
    // 阈值划分10次 只保留千分之一的顶点 结果波动还很大，
    int recur_depth = 50;
    int vertex_depth = 0;

    uint32_t threshhold[recur_depth + 1] = {0};
    int *record_degree = (int *)malloc(sizeof(int) * g->v_cnt);
    memset(record_degree, 0, sizeof(int) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        record_degree[g->vertex[i + 1] - g->vertex[i]]++;
    }
    double ratio = 0;
    double ratio_threshhold = 0.5;
    double ratio_acum = 0.5;
    // threshhold[0] = 1;
    int iter_threshhoud = 1;
    uint64_t iter_record = 0;
    double change_ratio = 0.0;

    while (iter_record < g->v_cnt) {
        double temp_ratio = double(record_degree[iter_record]) / double(g->v_cnt);
        ratio += temp_ratio;
        iter_record++;
        if (ratio > ratio_threshhold) {
            threshhold[iter_threshhoud++] = iter_record;
            // change_ratio = double(threshhold[iter_threshhoud - 1]) / double(threshhold[iter_threshhoud - 2]);
            ratio_acum = ratio_acum * 0.5;
            ratio_threshhold = ratio_threshhold + ratio_acum;
            if (iter_threshhoud == recur_depth || iter_record > 20000 || iter_threshhoud > 18) { //|| change_ratio < 1.1
                break;
            }
        }
    }
    for (int i = 1; i < recur_depth; i++) {
        if (threshhold[i] == 0) {
            recur_depth = i;
            break;
        }
    }
    threshhold[recur_depth] = g->v_cnt;
    free(record_degree);

    // 一共有recur_depth中颜色， 但是查找阈值的时候需要用recur_depth +1.因为阈值第一位为0

    // 查看阈值
    // for (int i = 0; i < recur_depth; i++) {
    //     printf("threshhold[%d] = %d\n", i, threshhold[i]);
    // }

    // 过滤的数据结构。
    filter filter_data;
    filter_data.edge_project = (int64_t *)malloc(sizeof(int64_t) * g->e_cnt);
    filter_data.edge_record = (int *)malloc(sizeof(int) * g->e_cnt);
    filter_data.vertex_project = (int *)malloc(sizeof(int) * g->v_cnt);
    filter_data.vertex_record = (int *)malloc(sizeof(int) * g->v_cnt);
    memset(filter_data.edge_project, -1, sizeof(int64_t) * g->e_cnt);
    memset(filter_data.edge_record, -1, sizeof(int) * g->e_cnt);
    memset(filter_data.vertex_record, -1, sizeof(int) * g->v_cnt);
    memset(filter_data.vertex_project, -1, sizeof(int) * g->v_cnt);
    filter_data.delete_edge = (bool *)malloc(sizeof(bool) * g->e_cnt);
    memset(filter_data.delete_edge, false, sizeof(bool) * g->e_cnt);

    filter_data.degree = (uint32_t *)malloc(sizeof(uint32_t) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        filter_data.degree[i] = g->vertex[i + 1] - g->vertex[i];
    }

    uint64_t discard_vertex = 0;

    // 标记顶点和边，以及记录数据
    uint64_t num_of_different_weight_of_edges[recur_depth] = {0};
    uint64_t num_of_different_weight_of_vertex[recur_depth] = {0};
    for (uint64_t i = 0; i < g->v_cnt; ++i) {
        if (filter_data.degree[i] >= s->order_degree[0]) {
            int record_vertex_color = find_range(threshhold, filter_data.degree[i], recur_depth + 1);
            filter_data.vertex_record[i] = record_vertex_color;
            num_of_different_weight_of_vertex[record_vertex_color]++;
            for (uint64_t j = g->vertex[i]; j < g->vertex[i + 1]; ++j) {
                uint32_t temp_degree = filter_data.degree[g->edge[j]];
                if (temp_degree >= s->order_degree[1]) {
                    int record_edge_color = find_range(threshhold, temp_degree, recur_depth + 1);
                    filter_data.edge_record[j] = record_edge_color;
                    num_of_different_weight_of_edges[record_edge_color]++;
                }
            }
        } else {
            discard_vertex++;
        }
    }

    // 投射顶点，记录所有的顶点
    uint64_t prefix_sum_of_vertex[recur_depth + 1] = {0};
    for (int i = 1; i < recur_depth + 1; i++) {
        prefix_sum_of_vertex[i] = prefix_sum_of_vertex[i - 1] + num_of_different_weight_of_vertex[i - 1];
    }
    uint64_t weight_vtx_ptr[recur_depth] = {0};
    for (int i = 1; i < recur_depth; i++) {
        weight_vtx_ptr[i] = prefix_sum_of_vertex[i];
    }
    for (uint32_t i = 0; i < g->v_cnt; ++i) {
        if (filter_data.vertex_record[i] != -1)
            filter_data.vertex_project[weight_vtx_ptr[filter_data.vertex_record[i]]++] = i;
    }

    // 投射边
    // 前缀和形式
    uint64_t prefix_sum_of_edges[recur_depth + 1] = {0};
    for (int i = 1; i < recur_depth + 1; i++) {
        prefix_sum_of_edges[i] = prefix_sum_of_edges[i - 1] + num_of_different_weight_of_edges[i - 1];
    }
    // 利用前缀和更新指针
    uint64_t weight_edge_ptr[recur_depth] = {0};
    for (int i = 0; i < recur_depth; i++) {
        weight_edge_ptr[i] = prefix_sum_of_edges[i];
    }
    for (uint32_t i = 0; i < g->e_cnt; ++i) {
        if (filter_data.edge_record[i] != -1)
            filter_data.edge_project[weight_edge_ptr[filter_data.edge_record[i]]++] = i;
    }

    // 记录采样结果
    filter_data.sample_times_per_edge = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
    filter_data.sample_count_per_edge = (double *)malloc(sizeof(double) * g->e_cnt);
    filter_data.sample_square_per_edge = (double *)malloc(sizeof(double) * g->e_cnt);
    memset(filter_data.sample_times_per_edge, 0, sizeof(uint32_t) * g->e_cnt);
    memset(filter_data.sample_count_per_edge, 0, sizeof(double) * g->e_cnt);
    memset(filter_data.sample_square_per_edge, 0, sizeof(double) * g->e_cnt);

    // 采样批量
    uint64_t batch_size = 1000;

    // 随机数生成器
    std::random_device rd;
    std::default_random_engine gen(rd());

    // 第一次采样的范围 ，*****后续要实时更新
    uint64_t range = prefix_sum_of_edges[recur_depth];
    // 采样记录
    uint64_t approx_sample_times = 0;
    uint64_t real_sample_times = 0;
    double approx_count = 0;
    double approx_squra = 0;
    uint64_t mask = 1ULL;
    mask <<= 63;
    bool flag[4] = {true, true, true, true};

    // 主线程判断是否满足要求使用到的变量，只有主线程会使用到

    double var_local[recur_depth] = {0};
    double std_dev_local[recur_depth] = {0};
    double estimate_error_local[recur_depth] = {0};
    uint64_t approx_sample_times_local[recur_depth] = {0};
    double approx_count_local[recur_depth] = {0};
    double approx_squra_local[recur_depth] = {0};
    double real_error_gloabl = 0.0;
    double real_error = 1.0;

    // 线程工作范围
    int num_of_thread = omp_get_max_threads();
    uint64_t start_edge_id[num_of_thread] = {0};
    uint64_t end_edge_id[num_of_thread];
    for (int i = 0; i < num_of_thread; i++) {
        end_edge_id[i] = prefix_sum_of_edges[recur_depth];
    }

    // 精确计算相关变量
    uint32_t processing_round = 0;
    uint32_t print_times = 0;
    bool exact_computing = false;
    bool get_ready[num_of_thread] = {false};
    uint64_t vertex_range_ptr = 0;
    uint64_t vertex_range_end = 0;
    uint64_t exact_count = 0;
    uint64_t vertex_chunk_per_task = 10;

    // 收敛速率检验相关
    double pre_real_error = 1.0;        // 记录前一个误差
    double pre_approx_sample_times = 0; // 记录前一个采样次数
    double coverge_range = 5000.0;      // 检查的范围差距根据处理速率不断调整。
    double coverge_range_ratio = 1.0;

    // 分层的总结
    double weight_average = 0.0;
    double weight_squra = 0.0;
    double weight_std_error = 0.0;
    double weight_real_error = 0.0;

    double error_threshhold[5] = {0.0};
    double var_threshhold[recur_depth] = {0.0};
    double var_threshhold_pre[recur_depth] = {0.0};

    double multiplicity = schedule_iep.get_multiplicity();
    double iep_redun = schedule_iep.get_in_exclusion_optimize_redundancy();

    // 误差 边界
#pragma omp parallel
    {
// 可以使用master来实现
#pragma omp master
        {
            while (real_error > 0.01) {

                //==================================================================================================================================================//
                //=================误差统计=========================================================================================================================//
                //==================================================================================================================================================//

                /* 统计误差 */
                // 总体方差等于局部方差除以采样次数之和， 总体均值为局部均值之和

                // double realreal_error = (real_average - (13283084785.0 * multiplicity)) / (13283084785.0 * multiplicity);

                for (int i = 1; i < num_of_thread; i++) {
                    if (get_ready[i] == true) {
                        i--;
                    }
                }
                // 分配精确计算的范围
                vertex_range_ptr = prefix_sum_of_vertex[0];           // 顶点处理指针
                vertex_range_end = prefix_sum_of_vertex[recur_depth]; // 顶点处理范围
                // 打开所有线程的精确计算开关
                for (int i = 1; i < num_of_thread; i++) {
                    get_ready[i] = true;
                }

                // 最后一部分顶点，查看所有精确计算是否完成任务

                for (int i = 1; i < num_of_thread; i++) {
                    if (get_ready[i] == true) {
                        i--;
                    }
                }
                double real_average = 0.0;
                real_average += (double)exact_count * multiplicity / iep_redun;
                real_error = 0.0;
            }
        }
        if (omp_get_thread_num()) {
            while (real_error > 0.01) {
                int cur_thread = omp_get_thread_num();
                // 如果需要精确计算
                if (get_ready[cur_thread]) {
                    while (vertex_range_ptr < vertex_range_end) {
                        // 获取当前线程的处理范围
                        uint32_t local_range_start = __sync_fetch_and_add(&vertex_range_ptr, vertex_chunk_per_task);
                        if (local_range_start < vertex_range_end) {
                            uint32_t local_range_end;
                            if (local_range_start + vertex_chunk_per_task < vertex_range_end)
                                local_range_end = local_range_start + vertex_chunk_per_task;
                            else
                                local_range_end = vertex_range_end;

                            /* 精确计算阶段 */
                            int *ans_buffer = new int[schedule_iep.in_exclusion_optimize_vertex_id.size()];  // 底层无交集顶点的大小
                            VertexSet *vertex_set = new VertexSet[schedule_iep.get_total_prefix_num() + 10]; // 构建了多个顶点集合作为候选集
                            VertexSet subtraction_set;                                                       // 构建一个差集
                            VertexSet tmp_set;                                                               // 构建一个临时集合
                            subtraction_set.init();
                            long long local_ans = 0;
                            for (int i = local_range_start; i < local_range_end; ++i) {
                                uint32_t vertex = filter_data.vertex_project[i];
                                e_index_t l, r;
                                l = g->vertex[vertex];
                                r = g->vertex[vertex + 1];

                                // 为调度中的每一个顶点初始化候选顶点集合。
                                for (int prefix_id = schedule_iep.get_last(0); prefix_id != -1; prefix_id = schedule_iep.get_next(prefix_id)) {
                                    vertex_set[prefix_id].build_vertex_set(schedule_iep, vertex_set, &g->edge[l], (int)(r - l), prefix_id);
                                }
                                subtraction_set.push_back(vertex);
                                pattern_matching_aggressive_func(schedule_iep, vertex_set, subtraction_set, tmp_set, local_ans, 1, ans_buffer, g);
                                subtraction_set.pop_back();
                            }
                            delete[] vertex_set;
                            delete[] ans_buffer;

#pragma omp atomic
                            exact_count += local_ans;
                        }
                    }

                    get_ready[cur_thread] = false;
                    exact_computing = false;
                }
            }
        }
    }

    free(filter_data.edge_project);
    free(filter_data.edge_record);
    free(filter_data.vertex_project);
    free(filter_data.vertex_record);
    free(filter_data.degree);

    free(filter_data.sample_count_per_edge);
    free(filter_data.sample_times_per_edge);
    free(filter_data.sample_square_per_edge);
    free(filter_data.delete_edge);
}

void exact_count_test(Graph *g, schedule_approx *s, const Schedule_IEP &schedule_iep, Pattern p) {
    long long global_ans = 0;
    // printf("pattern_matching extra memory: %.3lf MB\n", thread_count * (schedule.get_total_prefix_num() + 10) * sizeof(int) *
    // (VertexSet::max_intersection_size * 2) / 1024.0 / 1024.0); fflush(stdout);
#pragma omp parallel reduction(+ : global_ans)
    {
        //   double start_time = get_wall_time();
        //   double current_time;
        int *ans_buffer = new int[schedule_iep.in_exclusion_optimize_vertex_id.size()];  // 底层无交集顶点的大小
        VertexSet *vertex_set = new VertexSet[schedule_iep.get_total_prefix_num() + 10]; // 构建了多个顶点集合作为候选集
        VertexSet subtraction_set;                                                       // 构建一个差集
        VertexSet tmp_set;                                                               // 构建一个临时集合
        subtraction_set.init();
        long long local_ans = 0;

        // 对于每一个顶点进行处理
#pragma omp for schedule(dynamic)
        for (int vertex = 0; vertex < g->v_cnt; ++vertex) {

            e_index_t l, r;
            l = g->vertex[vertex];
            r = g->vertex[vertex + 1];

            // 为调度中的每一个顶点初始化候选顶点集合。
            for (int prefix_id = schedule_iep.get_last(0); prefix_id != -1; prefix_id = schedule_iep.get_next(prefix_id)) {
                vertex_set[prefix_id].build_vertex_set(schedule_iep, vertex_set, &g->edge[l], (int)(r - l), prefix_id);
            }

            subtraction_set.push_back(vertex); // 将当前顶点放入差集中，也就是后续的所有顶点都不能够是该顶点
            pattern_matching_aggressive_func(schedule_iep, vertex_set, subtraction_set, tmp_set, local_ans, 1, ans_buffer, g);
            subtraction_set.pop_back();
            // printf("for %d %d\n", omp_get_thread_num(), vertex);
        }
        delete[] vertex_set;
        // TODO : Computing multiplicty for a pattern
        global_ans += local_ans;
        // printf("local_ans %d %lld\n", omp_get_thread_num(), local_ans);
    }
    global_ans = global_ans / schedule_iep.get_in_exclusion_optimize_redundancy();
    printf("exact_count %lld\n", global_ans);
}

void approximate_exact_fusion_eightyseven(Graph *g, schedule_approx *s, const Schedule_IEP &schedule_iep, Pattern p) {

    edge_from_init(g);

    // 平均度数  frendster 28 twitter 35 livejouranl 17 youtube 2 patent 5 pokec 30 tp 28 sk 11
    // 阈值划分10次 只保留千分之一的顶点 结果波动还很大，
    int recur_depth = 50;
    int vertex_depth = 0;

    uint32_t threshhold[recur_depth + 1] = {0};
    int *record_degree = (int *)malloc(sizeof(int) * g->v_cnt);
    memset(record_degree, 0, sizeof(int) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        record_degree[g->vertex[i + 1] - g->vertex[i]]++;
    }
    double ratio = 0;
    double ratio_threshhold = 0.5;
    double ratio_acum = 0.5;
    // threshhold[0] = 1;
    int iter_threshhoud = 1;
    uint64_t iter_record = 0;
    double change_ratio = 0.0;

    while (iter_record < g->v_cnt) {
        double temp_ratio = double(record_degree[iter_record]) / double(g->v_cnt);
        ratio += temp_ratio;
        iter_record++;
        if (ratio > ratio_threshhold) {
            threshhold[iter_threshhoud++] = iter_record;
            // change_ratio = double(threshhold[iter_threshhoud - 1]) / double(threshhold[iter_threshhoud - 2]);
            ratio_acum = ratio_acum * 0.5;
            ratio_threshhold = ratio_threshhold + ratio_acum;
            if (iter_threshhoud == recur_depth || iter_record > 20000 || iter_threshhoud > 15) { //|| change_ratio < 1.1
                break;
            }
        }
    }
    for (int i = 1; i < recur_depth; i++) {
        if (threshhold[i] == 0) {
            recur_depth = i;
            break;
        }
    }
    threshhold[recur_depth] = g->v_cnt;
    free(record_degree);

    // 一共有recur_depth中颜色， 但是查找阈值的时候需要用recur_depth +1.因为阈值第一位为0

    // 查看阈值
    // for (int i = 0; i < recur_depth; i++) {
    //     printf("threshhold[%d] = %d\n", i, threshhold[i]);
    // }

    // 过滤的数据结构。
    filter filter_data;
    filter_data.edge_project = (int64_t *)malloc(sizeof(int64_t) * g->e_cnt);
    filter_data.edge_record = (int *)malloc(sizeof(int) * g->e_cnt);
    filter_data.vertex_project = (int *)malloc(sizeof(int) * g->v_cnt);
    filter_data.vertex_record = (int *)malloc(sizeof(int) * g->v_cnt);
    memset(filter_data.edge_project, -1, sizeof(int64_t) * g->e_cnt);
    memset(filter_data.edge_record, -1, sizeof(int) * g->e_cnt);
    memset(filter_data.vertex_record, -1, sizeof(int) * g->v_cnt);
    memset(filter_data.vertex_project, -1, sizeof(int) * g->v_cnt);
    filter_data.delete_edge = (bool *)malloc(sizeof(bool) * g->e_cnt);
    memset(filter_data.delete_edge, false, sizeof(bool) * g->e_cnt);

    filter_data.degree = (uint32_t *)malloc(sizeof(uint32_t) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        filter_data.degree[i] = g->vertex[i + 1] - g->vertex[i];
    }

    uint64_t discard_vertex = 0;

    // 标记顶点和边，以及记录数据
    uint64_t num_of_different_weight_of_edges[recur_depth] = {0};
    uint64_t num_of_different_weight_of_vertex[recur_depth] = {0};
    for (uint64_t i = 0; i < g->v_cnt; ++i) {
        if (filter_data.degree[i] >= s->order_degree[0]) {
            int record_vertex_color = find_range(threshhold, filter_data.degree[i], recur_depth + 1);
            filter_data.vertex_record[i] = record_vertex_color;
            num_of_different_weight_of_vertex[record_vertex_color]++;
            for (uint64_t j = g->vertex[i]; j < g->vertex[i + 1]; ++j) {
                uint32_t temp_degree = filter_data.degree[g->edge[j]];
                if (temp_degree >= s->order_degree[1]) {
                    int record_edge_color = find_range(threshhold, temp_degree, recur_depth + 1);
                    filter_data.edge_record[j] = record_edge_color;
                    num_of_different_weight_of_edges[record_edge_color]++;
                }
            }
        } else {
            discard_vertex++;
        }
    }

    // 投射顶点，记录所有的顶点
    uint64_t prefix_sum_of_vertex[recur_depth + 1] = {0};
    for (int i = 1; i < recur_depth + 1; i++) {
        prefix_sum_of_vertex[i] = prefix_sum_of_vertex[i - 1] + num_of_different_weight_of_vertex[i - 1];
    }
    uint64_t weight_vtx_ptr[recur_depth] = {0};
    for (int i = 1; i < recur_depth; i++) {
        weight_vtx_ptr[i] = prefix_sum_of_vertex[i];
    }
    for (uint32_t i = 0; i < g->v_cnt; ++i) {
        if (filter_data.vertex_record[i] != -1)
            filter_data.vertex_project[weight_vtx_ptr[filter_data.vertex_record[i]]++] = i;
    }

    // 投射边
    // 前缀和形式
    uint64_t prefix_sum_of_edges[recur_depth + 1] = {0};
    for (int i = 1; i < recur_depth + 1; i++) {
        prefix_sum_of_edges[i] = prefix_sum_of_edges[i - 1] + num_of_different_weight_of_edges[i - 1];
    }
    // 利用前缀和更新指针
    uint64_t weight_edge_ptr[recur_depth] = {0};
    for (int i = 0; i < recur_depth; i++) {
        weight_edge_ptr[i] = prefix_sum_of_edges[i];
    }
    for (uint32_t i = 0; i < g->e_cnt; ++i) {
        if (filter_data.edge_record[i] != -1)
            filter_data.edge_project[weight_edge_ptr[filter_data.edge_record[i]]++] = i;
    }

    // 记录采样结果
    filter_data.sample_times_per_edge = (uint32_t *)malloc(sizeof(uint32_t) * g->e_cnt);
    filter_data.sample_count_per_edge = (double *)malloc(sizeof(double) * g->e_cnt);
    filter_data.sample_square_per_edge = (double *)malloc(sizeof(double) * g->e_cnt);
    memset(filter_data.sample_times_per_edge, 0, sizeof(uint32_t) * g->e_cnt);
    memset(filter_data.sample_count_per_edge, 0, sizeof(double) * g->e_cnt);
    memset(filter_data.sample_square_per_edge, 0, sizeof(double) * g->e_cnt);

    // 采样批量
    uint64_t batch_size = 1000;

    // 随机数生成器
    std::random_device rd;
    std::default_random_engine gen(rd());

    // 第一次采样的范围 ，*****后续要实时更新
    uint64_t range = prefix_sum_of_edges[recur_depth];
    // 采样记录
    uint64_t approx_sample_times = 0;
    uint64_t real_sample_times = 0;
    double approx_count = 0;
    double approx_squra = 0;
    uint64_t mask = 1ULL;
    mask <<= 63;
    bool flag[4] = {true, true, true, true};

    // 主线程判断是否满足要求使用到的变量，只有主线程会使用到

    double var_local[recur_depth] = {0};
    double std_dev_local[recur_depth] = {0};
    double estimate_error_local[recur_depth] = {0};
    uint64_t approx_sample_times_local[recur_depth] = {0};
    double approx_count_local[recur_depth] = {0};
    double approx_squra_local[recur_depth] = {0};
    double real_error_gloabl = 0.0;
    double real_error = 1.0;

    // 线程工作范围
    int num_of_thread = omp_get_max_threads();
    uint64_t start_edge_id[num_of_thread] = {0};
    uint64_t end_edge_id[num_of_thread];
    for (int i = 0; i < num_of_thread; i++) {
        end_edge_id[i] = prefix_sum_of_edges[recur_depth];
    }

    // 精确计算相关变量
    uint32_t processing_round = 0;
    uint32_t print_times = 0;
    bool exact_computing = true;
    bool get_ready[num_of_thread] = {false};
    uint64_t vertex_range_ptr = 0;
    uint64_t vertex_range_end = 0;
    uint64_t exact_count = 0;
    uint64_t vertex_chunk_per_task = 80;

    // 收敛速率检验相关
    double pre_real_error = 1.0;        // 记录前一个误差
    double pre_approx_sample_times = 0; // 记录前一个采样次数
    double coverge_range = 5000.0;      // 检查的范围差距根据处理速率不断调整。
    double coverge_range_ratio = 1.0;

    // 分层的总结
    double weight_average = 0.0;
    double weight_squra = 0.0;
    double weight_std_error = 0.0;
    double weight_real_error = 0.0;

    double error_threshhold[5] = {0.0};
    double var_threshhold[recur_depth] = {0.0};
    double var_threshhold_pre[recur_depth] = {0.0};

    double multiplicity = schedule_iep.get_multiplicity();
    double iep_redun = schedule_iep.get_in_exclusion_optimize_redundancy();

    double var = 0.0;
    double real_average = 0.0;

    // 误差 边界
#pragma omp parallel
    {
// 可以使用master来实现
#pragma omp master
        {
            while (real_error > 0.1 && approx_sample_times <= 100000000) {

                // 1，精确计算百分之五十后不再继续精确计算。
                if (approx_sample_times > 1000 * print_times) {

                    // 1，更新精确计算范围。 2，打开精确计算开关。3，可能存在前一个精确计算还没有结束就进入精确计算的情况。
                    if (exact_computing) {
                        for (int i = 1; i < num_of_thread; i++) {
                            if (get_ready[i] == true) {
                                i--;
                            }
                        }

                        // 分配精确计算的范围
                        vertex_range_ptr = prefix_sum_of_vertex[0]; // 顶点处理指针
                        vertex_range_end = prefix_sum_of_vertex[3]; // 顶点处理范围
                        processing_round = processing_round + 2;
                        // 打开所有线程的精确计算开关
                        for (int i = 1; i < num_of_thread; i++) {
                            get_ready[i] = true;
                        }

                        // 查看所有精确计算是否完成任务
                        if (processing_round + 1 == 3) {
                            for (int i = 1; i < num_of_thread; i++) {
                                if (get_ready[i] == true) {
                                    i--;
                                }
                            }
                        }

                        // 记录现有采样范围，以及更新边范围相关数据
                        uint64_t size = prefix_sum_of_edges[recur_depth];
                        for (int i = 0; i < recur_depth; i++) {
                            num_of_different_weight_of_edges[i] = 0;
                            prefix_sum_of_edges[i] = 0;
                        }

                        uint64_t temp_ptr = 0;
                        int pre_edge_color = 0;
                        int cur_edge_color = 0;
                        int cur_vtx_color = 0;
                        for (int i = 0; i < size; i++) {
                            uint64_t temp_edge = filter_data.edge_project[i];
                            if (filter_data.delete_edge[temp_edge]) {
                                continue;
                            }
                            uint64_t temp_degree = filter_data.degree[g->edge_from[temp_edge]];
                            cur_vtx_color = find_range(threshhold, temp_degree, recur_depth + 1);
                            // 起点度数够大，保留
                            if (cur_vtx_color > processing_round) {
                                // 记录边并且获取边的颜色
                                filter_data.edge_project[temp_ptr++] = temp_edge;
                                cur_edge_color = find_range(threshhold, filter_data.degree[g->edge[temp_edge]], recur_depth + 1);

                                // 更新相应颜色的范围
                                if (cur_edge_color != pre_edge_color) {
                                    pre_edge_color = cur_edge_color;
                                    num_of_different_weight_of_edges[cur_edge_color]++;
                                    prefix_sum_of_edges[cur_edge_color] = temp_ptr - 1;

                                } else {
                                    num_of_different_weight_of_edges[pre_edge_color]++;
                                }
                            }
                        }
                        prefix_sum_of_edges[recur_depth] = temp_ptr;
                        processing_round++;
                    }

                    //==================================================================================================================================================//
                    //=================误差统计=========================================================================================================================//
                    //==================================================================================================================================================//
                    /* 局部误差 */ // 记录平方和，计数和，局部采样次数
                    for (int i = 0; i < recur_depth; i++) {
                        approx_sample_times_local[i] = 0.0;
                        approx_count_local[i] = 0.0;
                        approx_squra_local[i] = 0.0;
                    }
                    for (int i = 0; i < recur_depth; i++) {
                        for (int ptr = prefix_sum_of_edges[i]; ptr < prefix_sum_of_edges[i + 1]; ptr++) {
                            auto eid = filter_data.edge_project[ptr];
                            approx_sample_times_local[i] += filter_data.sample_times_per_edge[eid];
                            approx_count_local[i] += filter_data.sample_count_per_edge[eid];
                            approx_squra_local[i] += filter_data.sample_square_per_edge[eid];
                        }
                    }
                    // 局部均值的偏差
                    for (int i = 0; i < recur_depth; i++) {
                        var_local[i] = (double)(approx_squra_local[i] / approx_sample_times_local[i] -
                                                std::pow((approx_count_local[i] / approx_sample_times_local[i]), 2));
                    }

                    // 局部总和的标准偏差
                    // for (int i = 0; i < recur_depth; i++) {
                    //     double cur_errror = std::sqrt(var_local[i] / approx_sample_times_local[i]) * num_of_different_weight_of_edges[i];
                    //     std::cout << " 线程" << i << "偏差 " << cur_errror << std::endl;
                    // }

                    /* 统计误差 */
                    // 总体方差等于局部方差除以采样次数之和， 总体均值为局部均值之和
                    var = 0.0;
                    real_average = 0.0;
                    for (int i = 0; i < recur_depth; i++) {
                        var +=
                            (var_local[i] * num_of_different_weight_of_edges[i] * num_of_different_weight_of_edges[i]) / approx_sample_times_local[i];
                        real_average += (approx_count_local[i] / approx_sample_times_local[i] * num_of_different_weight_of_edges[i]);
                    }
                    real_average += (double)exact_count * multiplicity / iep_redun;
                    double std_dev = std::sqrt(var);
                    real_error = std_dev * 2.576 / real_average; // 99%的信心

                    // double realreal_error = (real_average - (7375094981 * 2)) / (7375094981 * 2); // 7375094981 5007311619323
                    // double realreal_error = (real_average - (13283084785.0 * multiplicity)) / (13283084785.0 * multiplicity);

                    // printf("approx_count: %lu, real_error: %lf, estimate_error: %lf\n", approx_sample_times, realreal_error, real_error);

                    if (flag[0] && real_error < 0.1) {
                        printf("approx_count: %lf, approx_sample_times: %lu, real_error: %lf\n", real_average, approx_sample_times, real_error);
                        flag[0] = false;
                    }
                    if (flag[1] && real_error < 0.05) {
                        printf("approx_count: %lf, approx_sample_times: %lu, real_error: %lf\n", real_average, approx_sample_times, real_error);
                        flag[1] = false;
                    }

                    print_times++;

                    //==================================================================================================================================================//
                    //=========最后调整工作线程的工作范围==================================================================================================================================//
                    //==================================================================================================================================================//

                    /* 自适应的修改采样范围 */
                    // 记录 偏差阈值 迭代上升。
                    for (int i = 0; i < recur_depth; i++) {
                        var_threshhold[i] = std::sqrt(var_local[i] / approx_sample_times_local[i]) * num_of_different_weight_of_edges[i] /
                                            (approx_count_local[i] / approx_sample_times_local[i]);
                        var_threshhold_pre[i] = var_threshhold[i];
                    }
                    std::sort(var_threshhold, var_threshhold + recur_depth);
                    // for (int i = 0; i < 5; i++) {
                    //     int var_threshhold_index = recur_depth - std::pow(2, (i + 1)) + 1;
                    //     if (var_threshhold_index >= 0)
                    //         error_threshhold[i] = var_threshhold[var_threshhold_index] - 1;
                    // }

                    int remain_thread = num_of_thread;
                    int allocate_thread_ptr = 0;
                    int allocate_thread = 0;
                    int allocate_threshhold = 0;
                    int adaptive_start = 0;
                    while (remain_thread > num_of_thread / 2) {
                        allocate_thread = remain_thread / 4;
                        // 查找当前应该分配的阈值的下标。
                        for (int i = 0; i < recur_depth; i++) {
                            if (var_threshhold_pre[i] > error_threshhold[allocate_threshhold]) {
                                adaptive_start = i;
                            }
                        }
                        allocate_threshhold++;
                        for (int i = allocate_thread_ptr; i < allocate_thread_ptr + allocate_thread; i++) {
                            start_edge_id[i] = prefix_sum_of_edges[adaptive_start];
                            end_edge_id[i] = prefix_sum_of_edges[recur_depth];
                        }
                        remain_thread -= allocate_thread;
                        allocate_thread_ptr += allocate_thread;
                    }
                    for (int i = allocate_thread_ptr; i < num_of_thread; i++) {
                        end_edge_id[i] = prefix_sum_of_edges[recur_depth];
                    }
                }
            }
            if (flag[2] && real_error <= 0.01) {
                printf("approx_count: %lf, approx_sample_times: %lu, real_error: %lf\n", real_average, approx_sample_times, real_error);
                flag[2] = false;
            }
        }
        if (omp_get_thread_num()) {
            while (real_error > 0.1 && approx_sample_times <= 100000000) {
                int cur_thread = omp_get_thread_num();
                // 如果需要精确计算
                if (get_ready[cur_thread]) {
                    while (vertex_range_ptr < vertex_range_end) {
                        // 获取当前线程的处理范围
                        uint32_t local_range_start = __sync_fetch_and_add(&vertex_range_ptr, vertex_chunk_per_task);
                        if (local_range_start < vertex_range_end) {
                            uint32_t local_range_end;
                            if (local_range_start + vertex_chunk_per_task < vertex_range_end)
                                local_range_end = local_range_start + vertex_chunk_per_task;
                            else
                                local_range_end = vertex_range_end;

                            /* 精确计算阶段 */
                            int *ans_buffer = new int[schedule_iep.in_exclusion_optimize_vertex_id.size()];  // 底层无交集顶点的大小
                            VertexSet *vertex_set = new VertexSet[schedule_iep.get_total_prefix_num() + 10]; // 构建了多个顶点集合作为候选集
                            VertexSet subtraction_set;                                                       // 构建一个差集
                            VertexSet tmp_set;                                                               // 构建一个临时集合
                            subtraction_set.init();
                            long long local_ans = 0;
                            for (int i = local_range_start; i < local_range_end; ++i) {
                                uint32_t vertex = filter_data.vertex_project[i];
                                e_index_t l, r;
                                l = g->vertex[vertex];
                                r = g->vertex[vertex + 1];

                                // 为调度中的每一个顶点初始化候选顶点集合。
                                for (int prefix_id = schedule_iep.get_last(0); prefix_id != -1; prefix_id = schedule_iep.get_next(prefix_id)) {
                                    vertex_set[prefix_id].build_vertex_set(schedule_iep, vertex_set, &g->edge[l], (int)(r - l), prefix_id);
                                }
                                subtraction_set.push_back(vertex);
                                pattern_matching_aggressive_func(schedule_iep, vertex_set, subtraction_set, tmp_set, local_ans, 1, ans_buffer, g);
                                subtraction_set.pop_back();
                            }
                            delete[] vertex_set;
                            delete[] ans_buffer;

#pragma omp atomic
                            exact_count += local_ans;
                        }
                    }

                    get_ready[cur_thread] = false;
                    exact_computing = false;
                }

                std::uniform_int_distribution<uint64_t> dist(start_edge_id[cur_thread],
                                                             end_edge_id[cur_thread] - 1); //***采样范围的更新可以不需要原子操作
                for (uint64_t i = 0; i < batch_size; i++) {
                    auto random_id = dist(gen);
                    auto eid = filter_data.edge_project[random_id]; // 一开始就处理最后一个区间，但是会不断过滤
                    pattern_sample_ns_record_filter(eid, g, s, &filter_data);
                    // pattern_sample_ns_record(eid, g, s, &filter_data);
                }

#pragma omp atomic
                approx_sample_times += batch_size;
            }
        }
    }

    free(filter_data.edge_project);
    free(filter_data.edge_record);
    free(filter_data.vertex_project);
    free(filter_data.vertex_record);
    free(filter_data.degree);

    free(filter_data.sample_count_per_edge);
    free(filter_data.sample_times_per_edge);
    free(filter_data.sample_square_per_edge);
    free(filter_data.delete_edge);
}

// 阈值划分
void threshhold_compute(Graph *g, uint32_t threshhold[], int &recur_depth) {

    // 记录每种度数的顶点数量
    int *record_degree = (int *)malloc(sizeof(int) * g->v_cnt);
    memset(record_degree, 0, sizeof(int) * g->v_cnt);
    for (int i = 0; i < g->v_cnt; i++) {
        record_degree[g->vertex[i + 1] - g->vertex[i]]++;
    }

    double ratio = 0;
    double ratio_threshhold = 0.5;
    double ratio_acum = 0.5;
    // threshhold[0] = 1;
    int iter_threshhoud = 1;
    uint64_t iter_record = 0;
    double change_ratio = 0.0;

    // 处理每一种度数的顶点，而不是每一个顶点
    while (iter_record < g->v_cnt) {

        // 当前种类的顶点战整体的比值
        double temp_ratio = double(record_degree[iter_record]) / double(g->v_cnt);
        ratio += temp_ratio;
        iter_record++;
        if (ratio > ratio_threshhold) {
            threshhold[iter_threshhoud++] = iter_record;
            // change_ratio = double(threshhold[iter_threshhoud - 1]) / double(threshhold[iter_threshhoud - 2]);
            ratio_acum = ratio_acum * 0.5;
            ratio_threshhold = ratio_threshhold + ratio_acum;
            if (iter_threshhoud == recur_depth || iter_record > 20000 || iter_threshhoud > 16) { //|| change_ratio < 1.1
                break;
            }
        }
    }
    for (int i = 1; i < recur_depth; i++)
        if (threshhold[i] == 0 || (threshhold[i] - threshhold[i - 1] < 20 && i > 10)) {
            recur_depth = i;
            break;
        }

    // 阈值里面记录的是顶点的度数
    threshhold[recur_depth] = g->v_cnt;
    free(record_degree);
}

void adaptive_ns(Graph *g, schedule_approx *s) {

    edge_from_init(g);

    // 平均度数  frendster 28 twitter 35 livejouranl 17 youtube 2 patent 5 pokec 30 tp 28 sk 11
    // 初始化阈值
    // 阈值划分10次 只保留千分之一的顶点 结果波动还很大，
    int recur_depth = 50;
    // 划分的阈值记录
    uint32_t threshhold[recur_depth + 1] = {0};
    threshhold_compute(g, threshhold, recur_depth);

    // 过滤数据的初始化
    filter filter_data;
    filter_data.num_of_color = recur_depth;
    filter_init(&filter_data, g, s, threshhold);

    //==================================================================================================================================================//
    //==================================================================================================================================================//
    //==================================================================================================================================================//

    // 采样批量
    uint64_t batch_size = 1000;

    // 随机数生成器
    std::random_device rd;
    std::default_random_engine gen(rd());

    // 采样记录
    uint64_t approx_sample_times = 0;
    uint64_t real_sample_times = 0;
    double approx_count = 0;
    double approx_squra = 0;
    bool flag[3] = {true, true, true};

    // 主误差和局部误差
    double var = 0.0;
    double std_dev = 0.0;
    double real_error = 1.0;

    double var_local[recur_depth] = {0};
    double std_dev_local[recur_depth] = {0};
    double estimate_error_local[recur_depth] = {0};
    uint64_t approx_sample_times_local[recur_depth] = {0};
    double approx_count_local[recur_depth] = {0};
    double approx_squra_local[recur_depth] = {0};
    double real_error_gloabl = 0.0;

    double real_average = 0.0;

    uint32_t print_times = 0; // 输出次数

    // 自适应分配范围
    uint64_t start_edge_id[filter_data.num_of_thread] = {0};
    uint64_t end_edge_id[filter_data.num_of_thread];
    for (int i = 0; i < filter_data.num_of_thread; i++)
        end_edge_id[i] = filter_data.prefix_sum_of_edges[recur_depth];

    double error_threshhold[10] = {0.0};
    double var_threshhold[recur_depth] = {0.0};
    double var_threshhold_pre[recur_depth] = {0.0};

    // // 保证每条边只被处理依次
    // int *record = (int *)malloc(sizeof(int) * g->e_cnt);
    // memset(record, 0, sizeof(int) * g->e_cnt);

    uint64_t global_sample_times;

    bool record_flag[2] = {true, true};

    // 误差 边界
#pragma omp parallel
    {
        // 可以使用master来实现
        while (real_error > 0.01 && approx_sample_times <= g->e_cnt * 10) {
#pragma omp master
            {
                // 更新采样区域
                if (approx_sample_times > 1000 * print_times) {

                    /* 局部误差 */
                    // 统计每个分区的估计数据，从每个线程采样数据中收集。
                    for (uint32_t i = 0; i < filter_data.num_of_color; i++) {
                        filter_data.sample_times_per_region[0][i] = 0;
                        filter_data.sample_count_per_region[0][i] = 0;
                        filter_data.sample_square_per_region[0][i] = 0;
                    }
                    for (int i = 1; i < filter_data.num_of_thread; i++) {
                        for (uint32_t j = 0; j < filter_data.num_of_color; j++) {
                            filter_data.sample_times_per_region[0][j] += filter_data.sample_times_per_region[i][j];
                            filter_data.sample_count_per_region[0][j] += filter_data.sample_count_per_region[i][j];
                            filter_data.sample_square_per_region[0][j] += filter_data.sample_square_per_region[i][j];
                        }
                    }

                    // 局部方差对于全局影响的计算
                    for (uint32_t i = 0; i < filter_data.num_of_color; i++) {
                        var_local[i] =
                            (double)(filter_data.sample_square_per_region[0][i] / filter_data.sample_times_per_region[0][i] -
                                     std::pow((filter_data.sample_count_per_region[0][i] / filter_data.sample_times_per_region[0][i]), 2)) *
                            std::pow(filter_data.num_of_different_weight_of_edges[i], 2) / filter_data.sample_times_per_region[0][i];
                    }

                    /* 统计误差 */
                    // 总体方差等于局部方差除以采样次数之和， 总体均值为局部均值之和
                    var = 0.0;
                    real_average = 0.0;
                    for (int i = 0; i < filter_data.num_of_color; i++) {
                        var += var_local[i];
                        real_average += (filter_data.sample_count_per_region[0][i] / filter_data.sample_times_per_region[0][i] *
                                         filter_data.num_of_different_weight_of_edges[i]);
                    }
                    std_dev = std::sqrt(var);
                    real_error = std_dev * 2.576 / real_average; // 99%的信心

                    if (std::isnan(real_error)) {
                        real_error = 1.0;
                    }

                    if (flag[0] && real_error < 0.1) {
                        printf("approx_count: %lf, approx_sample_times: %lu, real_error: %lf\n", real_average, approx_sample_times, real_error);
                        flag[0] = false;
                    }
                    if (flag[1] && real_error < 0.05) {
                        printf("approx_count: %lf, approx_sample_times: %lu, real_error: %lf\n", real_average, approx_sample_times, real_error);
                        flag[1] = false;
                    }
                    print_times++;

                    /* 自适应的修改采样范围 */
                    // 记录 偏差阈值 迭代上升。
                    for (int i = 0; i < recur_depth; i++) {
                        var_threshhold[i] = std::sqrt(var_local[i] / approx_sample_times_local[i]) * filter_data.num_of_different_weight_of_edges[i] /
                                            (approx_count_local[i] / approx_sample_times_local[i]);
                        var_threshhold_pre[i] = var_threshhold[i];
                    }
                    std::sort(var_threshhold, var_threshhold + recur_depth);

                    // 取大小为前十的边界
                    for (int i = 0; i < 10; i++) {
                        int var_threshhold_index = recur_depth - 1 - i;
                        if (var_threshhold_index >= 0)
                            error_threshhold[i] = var_threshhold[var_threshhold_index] - 1;
                    }

                    int remain_thread = filter_data.num_of_thread;
                    int allocate_thread_ptr = 0;
                    int allocate_thread = 0;
                    int allocate_threshhold = 0;
                    int adaptive_start = 0;
                    while (remain_thread > filter_data.num_of_thread / 2) {
                        allocate_thread = remain_thread / 4;
                        // 查找当前应该分配的阈值的下标。
                        for (int i = 0; i < recur_depth; i++) {
                            if (var_threshhold_pre[i] > error_threshhold[allocate_threshhold]) {
                                adaptive_start = i;
                            }
                        }
                        allocate_threshhold++;
                        for (int i = allocate_thread_ptr; i < allocate_thread_ptr + allocate_thread; i++) {
                            start_edge_id[i] = filter_data.prefix_sum_of_edges[adaptive_start];
                            // end_edge_id[i] = prefix_sum_of_edges[adaptive_start + 1];
                        }
                        remain_thread -= allocate_thread;
                        allocate_thread_ptr += allocate_thread;
                    }
                    for (int i = allocate_thread_ptr; i < filter_data.num_of_thread; i++) {
                        end_edge_id[i] = filter_data.prefix_sum_of_edges[recur_depth];
                    }
                }
            }

            if (omp_get_thread_num()) {
                int cur_thread = omp_get_thread_num();
                // 如果需要精确计算
                std::uniform_int_distribution<uint64_t> dist(start_edge_id[cur_thread],
                                                             end_edge_id[cur_thread] - 1); //***采样范围的更新可以不需要原子操作
                for (uint64_t i = 0; i < batch_size; i++) {
                    auto random_id = dist(gen);
                    auto eid = filter_data.edge_project[random_id]; // 一开始就处理最后一个区间，但是会不断过滤
                                                                    // pattern_sample_ns_record(eid, g, s, &filter_data);
                                                                    // pattern_sample_ns_record_filter(eid, g, s, &filter_data);
                    pattern_sample_ns_record_filter_region(eid, g, s, &filter_data, cur_thread, threshhold);
                }

#pragma omp atomic
                approx_sample_times += batch_size;
            }
        }
    }
    if (flag[2] && real_error <= 0.01) {
        printf("approx_count: %lf, approx_sample_times: %lu, real_error: %lf\n", real_average, approx_sample_times, real_error);
        flag[2] = false;
    }

    free_filter_data(filter_data);
}

void partition_ns(Graph *g, schedule_approx *s) {

    edge_from_init(g);
    // 平均度数  frendster 28 twitter 35 livejouranl 17 youtube 2 patent 5 pokec 30 tp 28 sk 11
    // 初始化阈值
    // 阈值划分10次 只保留千分之一的顶点 结果波动还很大，

    // 阈值计算
    // 递归划分的最大轮次
    int recur_depth = 50;
    // 划分的阈值记录
    uint32_t threshhold[recur_depth + 1] = {0};
    threshhold_compute(g, threshhold, recur_depth);

    // 过滤数据的初始化
    filter filter_data;
    filter_data.num_of_color = recur_depth;
    filter_init(&filter_data, g, s, threshhold);

    // 采样批量
    uint64_t batch_size = 1000;

    // 随机数生成器
    std::random_device rd;
    std::default_random_engine gen(rd());

    // 采样记录
    uint64_t approx_sample_times = 0;
    uint64_t real_sample_times = 0;
    double approx_count = 0;
    double approx_squra = 0;
    bool flag[3] = {true, true, true};
    uint32_t print_times = 0; // 输出次数

    // 主误差和局部误差
    double var = 0.0;
    double std_dev = 0.0;
    double real_error = 1.0;

    double var_local[recur_depth] = {0};
    double std_dev_local[recur_depth] = {0};
    double estimate_error_local[recur_depth] = {0};
    uint64_t approx_sample_times_local[recur_depth] = {0};
    double approx_count_local[recur_depth] = {0};
    double approx_squra_local[recur_depth] = {0};
    double real_error_gloabl = 0.0;
    double real_average = 0.0;

    // 自适应分配范围
    uint64_t start_edge_id[filter_data.num_of_thread] = {0};
    uint64_t end_edge_id[filter_data.num_of_thread];
    for (int i = 0; i < filter_data.num_of_thread; i++)
        end_edge_id[i] = filter_data.prefix_sum_of_edges[recur_depth];

#pragma omp parallel
    {

        while (real_error > 0.01 && approx_sample_times <= g->e_cnt * 10) {
// 可以使用master来实现
#pragma omp master
            {
                // 更新采样区域
                if (approx_sample_times > 1000 * print_times) {

                    /* 局部误差 */
                    // 统计每个分区的估计数据，从每个线程采样数据中收集。
                    for (uint32_t i = 0; i < filter_data.num_of_color; i++) {
                        filter_data.sample_times_per_region[0][i] = 0;
                        filter_data.sample_count_per_region[0][i] = 0;
                        filter_data.sample_square_per_region[0][i] = 0;
                    }
                    for (int i = 1; i < filter_data.num_of_thread; i++) {
                        for (uint32_t j = 0; j < filter_data.num_of_color; j++) {
                            filter_data.sample_times_per_region[0][j] += filter_data.sample_times_per_region[i][j];
                            filter_data.sample_count_per_region[0][j] += filter_data.sample_count_per_region[i][j];
                            filter_data.sample_square_per_region[0][j] += filter_data.sample_square_per_region[i][j];
                        }
                    }

                    // 局部方差对于全局影响的计算
                    for (uint32_t i = 0; i < filter_data.num_of_color; i++) {
                        var_local[i] =
                            (double)(filter_data.sample_square_per_region[0][i] / filter_data.sample_times_per_region[0][i] -
                                     std::pow((filter_data.sample_count_per_region[0][i] / filter_data.sample_times_per_region[0][i]), 2)) *
                            std::pow(filter_data.num_of_different_weight_of_edges[i], 2) / filter_data.sample_times_per_region[0][i];
                    }

                    /* 统计误差 */
                    // 总体方差等于局部方差除以采样次数之和， 总体均值为局部均值之和
                    var = 0.0;
                    real_average = 0.0;
                    for (int i = 0; i < filter_data.num_of_color; i++) {
                        var += var_local[i];
                        real_average += (filter_data.sample_count_per_region[0][i] / filter_data.sample_times_per_region[0][i] *
                                         filter_data.num_of_different_weight_of_edges[i]);
                    }
                    std_dev = std::sqrt(var);
                    real_error = std_dev * 2.576 / real_average; // 99%的信心

                    if (std::isnan(real_error)) {
                        real_error = 1.0;
                    }

                    if (flag[0] && real_error < 0.1) {
                        printf("approx_count: %lf, approx_sample_times: %lu, real_error: %lf\n", real_average, approx_sample_times, real_error);
                        flag[0] = false;
                    }
                    if (flag[1] && real_error < 0.05) {
                        printf("approx_count: %lf, approx_sample_times: %lu, real_error: %lf\n", real_average, approx_sample_times, real_error);
                        flag[1] = false;
                    }
                    print_times++;
                }
            }

            if (omp_get_thread_num()) {
                int cur_thread = omp_get_thread_num();
                // 如果需要精确计算
                std::uniform_int_distribution<uint64_t> dist(start_edge_id[cur_thread],
                                                             end_edge_id[cur_thread] - 1); //***采样范围的更新可以不需要原子操作
                for (uint64_t i = 0; i < batch_size; i++) {
                    auto random_id = dist(gen);
                    auto eid = filter_data.edge_project[random_id]; // 一开始就处理最后一个区间，但是会不断过滤
                                                                    // pattern_sample_ns_record(eid, g, s, &filter_data);
                                                                    // pattern_sample_ns_record_filter(eid, g, s, &filter_data);
                    pattern_sample_ns_record_filter_region(eid, g, s, &filter_data, cur_thread, threshhold);
                }
#pragma omp atomic
                approx_sample_times += batch_size;
            }
        }
    }
    if (flag[2] && real_error <= 0.01) {
        printf("approx_count: %lf, global_sample_times: %lu, real_error: %lf\n", real_average, approx_sample_times, real_error);
        flag[2] = false;
    }

    free_filter_data(filter_data);
}

void scale_gpm_ns(Graph *g, schedule_approx *s) {

    edge_from_init(g);

    // 采样批量
    uint64_t batch_size = 1000;

    // 随机数生成器
    std::random_device rd;
    std::default_random_engine gen(rd());

    // 采样范围
    uint64_t range = g->e_cnt;

    // 采样记录
    uint64_t approx_sample_times = 0;
    double approx_count = 0;
    double approx_squra = 0;

    // 主线程判断是否满足要求使用到的变量，只有主线程会使用到
    double var = 0.0;
    double std_dev = 0.0;
    double real_error = 1.0;
    int thread_count = omp_get_max_threads();

    uint64_t print_times = 0;

    bool flag[3] = {true, true, true};

    // 误差 边界
#pragma omp parallel
    {
// 可以使用master来实现
#pragma omp master
        {
            while (real_error > 0.05 && approx_sample_times <= 100000000) {
                /* 收集一轮预估值，计算误差 */
                if (approx_sample_times > 1000 * print_times) {
                    var = (double)(approx_squra / approx_sample_times - std::pow((approx_count / approx_sample_times), 2));
                    std_dev = std::sqrt(var / approx_sample_times);
                    real_error = std_dev * 2.576 / (approx_count / approx_sample_times); // 99%的信心

                    // double real_count = 26184567090560.0 * 2.0;
                    // double ave = approx_count / approx_sample_times;
                    // double std_err = ave - real_count;
                    // double realreal_error = std_err / real_count;
                    // double realreal_error = ((approx_count / approx_sample_times) - (341906226 * 8)) / (341906226 * 8);
                    // (real_average - (5007311619323.0 * 2.0)) / (5007311619323.0 * 2.0);
                    // printf("approx_count: %lf, realreal_error: %lf, real_error: %lf\n", approx_count / approx_sample_times, realreal_error,
                    //        real_error);
                    // printf("real error %lf\n", real_error);

                    if (flag[0] && real_error < 0.1) {
                        printf("approx_count: %lf, approx_sample_times: %lu, real_error: %lf\n", approx_count / approx_sample_times,
                               approx_sample_times, real_error);
                        fflush(stdout);
                        flag[0] = false;
                    }
                    if (flag[1] && real_error < 0.05) {
                        printf("approx_count: %lf, approx_sample_times: %lu, real_error: %lf\n", approx_count / approx_sample_times,
                               approx_sample_times, real_error);
                        flag[1] = false;
                        fflush(stdout);
                    }
                    print_times++;
                }
                if (flag[2] && real_error <= 0.01) {
                    printf("approx_count: %lf, std_approx_sample_timesdev: %lu, real_error: %lf\n", approx_count / approx_sample_times,
                           approx_sample_times, real_error);
                    flag[2] = false;
                }
            }
        }
        if (omp_get_thread_num()) {
            std::uniform_int_distribution<uint64_t> dist(0, range - 1);
            double sample_count = 0;
            double sample_squra = 0;
            while (real_error > 0.05 && approx_sample_times <= 100000000) {
                sample_count = 0;
                sample_squra = 0;

                for (uint64_t i = 0; i < batch_size; i++) {
                    auto random_id = dist(gen);
                    // 一开始就处理最后一个区间，但是会不断过滤
                    pattern_sample_ns(random_id, g, s, sample_count, sample_squra);
                }

#pragma omp critical
                {
                    approx_squra += sample_squra;
                    approx_count += sample_count;
                    approx_sample_times += batch_size;
                }
            }
        }
    }
}

int main(int argc, char *argv[]) {
    Graph *g;
    DataLoader D;

    if (argc != 3) {
        printf("usage: %s graph_file pattern_file\n", argv[0]);
        return 0;
    }
    // DataType type = DataType::skitter;
    bool ok = D.fast_load(g, argv[1]);
    // bool ok = D.load_data(g, type, argv[1], 0, 0);
    if (!ok) {
        printf("Load data failed\n");
        return 0;
    }

    printf("Load data success!\n");

    printf("Graph %s\n", argv[1]);
    fflush(stdout);
    printf("thread num: %d\n", omp_get_max_threads());

    TimeInterval allTime;

    // 模式载入，生成近似调度模式
    Pattern pattern = pattern_load(argv[2]);
    bool is_pattern_valid;
    bool use_in_exclusion_optimize = true;
    Schedule_IEP schedule_iep(pattern, is_pattern_valid, 1, 1, use_in_exclusion_optimize, g->v_cnt, g->e_cnt, g->tri_cnt);
    schedule_approx schedule;
    // scheduel_generate_approx_for_ns(&pattern, &schedule);
    scheduel_generate_for_adap_ns(schedule_iep, &schedule);
    g->degree = (uint32_t *)malloc(sizeof(uint32_t) * g->v_cnt);
    uint64_t max_degree = 0;
    for (int i = 0; i < g->v_cnt; i++) {
        g->degree[i] = g->vertex[i + 1] - g->vertex[i];
        if (g->degree[i] > max_degree) {
            max_degree = g->degree[i];
        }
    }

    double errorz = confidence_to_z(0.99); // 百分之99置信度的边界

    VertexSet::max_intersection_size = max_degree + 1;
    int multiplicity = schedule_iep.get_multiplicity();

    allTime.check();
    // pattern_count_ns_mpi_fix_sample_times(g, &schedule);
    // pattern_count_ns_mpi_fix_error(g, &schedule);
    // pattern_count_ns_mpi_fix_error_filter(g, &schedule);
    // pattern_count_ns_mpi_fix_error_filter_vertex_after_sample(g, &schedule, schedule_iep, pattern);
    // exact_count_test(g, &schedule, schedule_iep, pattern);
    // exact_count_version_two(g, &schedule, schedule_iep, pattern);

    // approximate_exact_fusion_fifty(g, &schedule, schedule_iep, pattern);
    // approximate_exact_fusion_seventyfive(g, &schedule, schedule_iep, pattern);
    // approximate_exact_fusion_eightyseven(g, &schedule, schedule_iep, pattern);
    // approximate_exact_fusion_fifty_allrange(g, &schedule, schedule_iep, pattern);
    // adaptive_approximate_exact_fusion(g, &schedule, schedule_iep, pattern);
    adaptive_ns(g, &schedule);
    // partition_ns(g, &schedule);
    // scale_gpm_ns(g, &schedule); // scalegpm的复现

    // locality_verify(g, &schedule);
    // exact_matching(g, schedule_iep);
    // vertex_vs_edge_test(g, &schedule, pattern);
    // pattern_count_ns(g, &schedule);
    // pattern_count_ns_with_symmetry_break(g, &schedule);
    // pattern_count_ns_partition(g, &schedule);
    // pattern_count_es(g, &schedule);
    // demo_pattern_test_exact(g);
    // demo_pattern_test_sample(g);
    // demo_pattern_test_filter_sample(g);
    // demo_pattern_test_filter_reduce_one(g);
    // demo_pattern_test_filter_reduce_two(g);
    // demo_pattern_test_filter_reduce_hybird(g);
    allTime.print("Total time cost");

    // g->triangle_statistics();
    // g->four_clique_statistics();

    delete g;
    return 0;
}

// 实验结果讨论：com-youtube 4clique  采用顺序交集 320秒，采用复用部分计算结果 72 秒 刚好差四倍的原因？
