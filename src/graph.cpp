#include "../include/graph.h"
#include "../include/common.h"
#include "../include/graphmpi.h"
#include "../include/motif_generator.h"
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
#include <mpi.h>
#include <omp.h>
#include <queue>
#include <random>
#include <sys/time.h>
#include <unistd.h>

void Graph::build_reverse_edges() {
    edge_from = new int[e_cnt];
    for (int u = 0; u < v_cnt; ++u)
        for (int v = vertex[u]; v < vertex[u + 1]; ++v)
            edge_from[v] = u;
}

int Graph::intersection_size(v_index_t v1, v_index_t v2) {
    e_index_t l1, r1;
    get_edge_index(v1, l1, r1);
    e_index_t l2, r2;
    get_edge_index(v2, l2, r2);
    int ans = 0;
    while (l1 < r1 && l2 < r2) {
        if (edge[l1] < edge[l2]) {
            ++l1;
        } else {
            if (edge[l2] < edge[l1]) {
                ++l2;
            } else {
                ++l1;
                ++l2;
                ++ans;
            }
        }
    }
    return ans;
}

int Graph::intersection_size_mpi(v_index_t v1, v_index_t v2) {
    Graphmpi &gm = Graphmpi::getinstance();
    int ans = 0;
    if (gm.include(v2))
        return intersection_size(v1, v2);
    e_index_t l1, r1;
    get_edge_index(v1, l1, r1);
    int *data = gm.getneighbor(v2);
    for (e_index_t l2 = 0; l1 < r1 && ~data[l2];) {
        if (edge[l1] < data[l2]) {
            ++l1;
        } else if (edge[l1] > data[l2]) {
            ++l2;
        } else {
            ++l1;
            ++l2;
            ++ans;
        }
    }
    return ans;
}

int Graph::intersection_size_clique(v_index_t v1, v_index_t v2) {
    e_index_t l1, r1;
    get_edge_index(v1, l1, r1);
    e_index_t l2, r2;
    get_edge_index(v2, l2, r2);
    v_index_t min_vertex = v2;
    int ans = 0;
    if (edge[l1] >= min_vertex || edge[l2] >= min_vertex)
        return 0;
    while (l1 < r1 && l2 < r2) {
        if (edge[l1] < edge[l2]) {
            if (edge[++l1] >= min_vertex)
                break;
        } else {
            if (edge[l2] < edge[l1]) {
                if (edge[++l2] >= min_vertex)
                    break;
            } else {
                ++ans;
                if (edge[++l1] >= min_vertex)
                    break;
                if (edge[++l2] >= min_vertex)
                    break;
            }
        }
    }
    return ans;
}

long long Graph::triangle_counting() {
    long long ans = 0;
    for (int v = 0; v < v_cnt; ++v) {
        // for v in G
        e_index_t l, r;
        get_edge_index(v, l, r);
        for (e_index_t v1 = l; v1 < r; ++v1) {
            // for v1 in N(v)
            ans += intersection_size(v, edge[v1]);
        }
    }
    ans /= 6;
    return ans;
}

long long Graph::triangle_counting_mt() {
    long long ans = 0;
#pragma omp parallel
    {
        tc_mt(&ans);
    }
    return ans;
}

void Graph::tc_mt(long long *global_ans) {
    long long my_ans = 0;
#pragma omp for schedule(dynamic)
    for (int v = 0; v < v_cnt; ++v) {
        // for v in G
        e_index_t l, r;
        get_edge_index(v, l, r);
        for (e_index_t v1 = l; v1 < r; ++v1) {
            if (v <= edge[v1])
                break;
            // for v1 in N(v)
            my_ans += intersection_size_clique(v, edge[v1]);
        }
    }
#pragma omp critical
    {
        *global_ans += my_ans;
    }
}

long long Graph::triangle_counting_mpi() {
    /*int mynodel, mynoder;
    long long tot_ans;
    Graphmpi &gm = Graphmpi::getinstance();
#pragma omp parallel
    {
#pragma omp master
        {
            auto k = gm.init(thread_count, this);
            mynodel = k.first;
            mynoder = k.second;
        }
#pragma omp barrier //mynodel have to be calculated before running other threads
#pragma omp master
        {
            tot_ans = gm.runmajor();
        }
        long long thread_ans = 0;
#pragma omp for schedule(dynamic) nowait
        for(int v = mynodel; v < mynoder; v++) {
            // for v in G
            int l, r;
            get_edge_index(v, l, r);
            for(int v1 = l; v1 < r; ++v1) {
                //for v1 in N(v)
                thread_ans += intersection_size_mpi(v, edge[v1]);
            }
        }
        //gm.idle(thread_ans);
    }
    return tot_ans / 6ll;*/
    return -1;
}

void Graph::get_edge_index(v_index_t v, e_index_t &l, e_index_t &r) const {
    l = vertex[v];
    r = vertex[v + 1];
}

template <typename T>
bool binary_search(const T data[], int n, const T &target) {
    int mid, l = 0, r = n - 1;
    while (l <= r) {
        mid = (l + r) >> 1;
        if (data[mid] < target) {
            l = mid + 1;
        } else if (data[mid] > target) {
            r = mid - 1;
        } else {
            return true;
        }
    }
    return false;
}

void Graph::remove_anti_edge_vertices(VertexSet &out_buf, const VertexSet &in_buf, const Schedule_IEP &sched, const VertexSet &partial_embedding,
                                      int vp) {

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
            get_edge_index(v, l, r);
            int m = r - l; // m = |N(v)|

            if (binary_search(&edge[l], m, d_in[i])) {
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

void Graph::pattern_matching_func(const Schedule_IEP &schedule, VertexSet *vertex_set, VertexSet &subtraction_set, long long &local_ans, int depth,
                                  bool clique) {
    int loop_set_prefix_id = schedule.get_loop_set_prefix_id(depth);
    auto &vset = vertex_set[loop_set_prefix_id];
    int loop_size = vset.get_size();
    if (loop_size <= 0)
        return;
    int *loop_data_ptr = vset.get_data_ptr();
    /*if (clique == true)
      {
      int last_vertex = subtraction_set.get_last();
    // The number of this vertex must be greater than the number of last vertex.
    loop_start = std::upper_bound(loop_data_ptr, loop_data_ptr + loop_size,
    last_vertex) - loop_data_ptr;
    }*/
    if (schedule.is_vertex_induced) {
        VertexSet &diff_buf = vertex_set[schedule.get_total_prefix_num() + depth];
        diff_buf.init();
        remove_anti_edge_vertices(diff_buf, vertex_set[loop_set_prefix_id], schedule, subtraction_set, depth);
        loop_data_ptr = diff_buf.get_data_ptr();
        loop_size = diff_buf.get_size();
        vset = diff_buf;
    }
    if (depth == schedule.get_size() - 1) {
        // TODO : try more kinds of calculation.
        // For example, we can maintain an ordered set, but it will cost more to
        // maintain itself when entering or exiting recursion.
        if (clique == true)
            local_ans += loop_size;
        else if (loop_size > 0)
            local_ans += VertexSet::unordered_subtraction_size(vset, subtraction_set);
        return;
    }

    int last_vertex = subtraction_set.get_last();
    for (int i = 0; i < loop_size; ++i) {
        if (last_vertex <= loop_data_ptr[i] && clique == true)
            break;
        int vertex = loop_data_ptr[i];
        if (!clique)
            if (subtraction_set.has_data(vertex))
                continue;
        e_index_t l, r;
        get_edge_index(vertex, l, r);
        bool is_zero = false;
        for (int prefix_id = schedule.get_last(depth); prefix_id != -1; prefix_id = schedule.get_next(prefix_id)) {
            vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &edge[l], (int)(r - l), prefix_id, vertex,
                                                   clique); // 为每个顶点重新构建候选解
            if (vertex_set[prefix_id].get_size() == 0) {
                is_zero = true;
                break;
            }
        }
        if (is_zero)
            continue;
        // subtraction_set.insert_ans_sort(vertex);
        subtraction_set.push_back(vertex);
        pattern_matching_func(schedule, vertex_set, subtraction_set, local_ans, depth + 1, clique);
        subtraction_set.pop_back();
    }
}

// 模式匹配函数
long long Graph::pattern_matching(const Schedule_IEP &schedule, bool clique) {
    //    intersection_times_low = intersection_times_high = 0;
    //    dep1_cnt = dep2_cnt = dep3_cnt = 0;
    long long global_ans = 0;
    // printf("pattern_matching extra memory: %.3lf MB\n", thread_count * (schedule.get_total_prefix_num() + 10) * sizeof(int) *
    // (VertexSet::max_intersection_size * 2) / 1024.0 / 1024.0); fflush(stdout);
#pragma omp parallel reduction(+ : global_ans)
    {
        //   double start_time = get_wall_time();
        //   double current_time;
        int *ans_buffer = new int[schedule.in_exclusion_optimize_vertex_id.size()];  // 底层无交集顶点的大小
        VertexSet *vertex_set = new VertexSet[schedule.get_total_prefix_num() + 10]; // 构建了多个顶点集合作为候选集
        VertexSet subtraction_set;                                                   // 构建一个差集
        VertexSet tmp_set;                                                           // 构建一个临时集合
        subtraction_set.init();
        long long local_ans = 0;
        // TODO : try different chunksize

        // 对于每一个顶点进行处理
#pragma omp for schedule(dynamic) nowait
        for (int vertex = 0; vertex < v_cnt; ++vertex) {

            e_index_t l, r;
            get_edge_index(vertex, l, r); // 获取顶点的边集的索引

            // 为调度中的每一个顶点初始化候选顶点集合。
            for (int prefix_id = schedule.get_last(0); prefix_id != -1; prefix_id = schedule.get_next(prefix_id)) {
                vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &edge[l], (int)(r - l), prefix_id);
            }
            // subtraction_set.insert_ans_sort(vertex);
            subtraction_set.push_back(vertex); // 将当前顶点放入差集中，也就是后续的所有顶点都不能够是该顶点
            // if (schedule.get_total_restrict_num() > 0 && clique == false)

            // 进行匹配操作
            if (true) {
                pattern_matching_aggressive_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, 1, ans_buffer);
            } else
                pattern_matching_func(schedule, vertex_set, subtraction_set, local_ans, 1, clique);
            subtraction_set.pop_back();
            // printf("for %d %d\n", omp_get_thread_num(), vertex);
        }
        // double end_time = get_wall_time();
        // printf("my thread time %d %.6lf\n", omp_get_thread_num(), end_time -
        // start_time);
        delete[] vertex_set;
        // TODO : Computing multiplicty for a pattern
        global_ans += local_ans;
        // printf("local_ans %d %lld\n", omp_get_thread_num(), local_ans);
    }
    return global_ans / schedule.get_in_exclusion_optimize_redundancy();
}

void Graph::clique_matching_func(const Schedule_IEP &schedule, VertexSet *vertex_set, Bitmap *bs, long long &local_ans, int depth) {
    int loop_set_prefix_id = schedule.get_loop_set_prefix_id(depth);
    int loop_size = vertex_set[loop_set_prefix_id].get_size();
    int *loop_data_ptr = vertex_set[loop_set_prefix_id].get_data_ptr();

    for (int i = 0; i < loop_size; ++i) {
        int vertex = loop_data_ptr[i];
        e_index_t l, r;
        get_edge_index(vertex, l, r);
        int prefix_id = schedule.get_last(depth); // only one prefix
        if (depth == schedule.get_size() - 2)
            vertex_set[prefix_id].build_vertex_set_bs_only_size(schedule, vertex_set, bs, &edge[l], (int)(r - l), prefix_id, depth);
        else
            vertex_set[prefix_id].build_vertex_set_bs(schedule, vertex_set, bs, &edge[l], (int)(r - l), prefix_id, depth);

        if (depth == schedule.get_size() - 2) {
            local_ans += vertex_set[prefix_id].get_size();
            continue;
        }
        if (vertex_set[prefix_id].get_size() == 0)
            continue;
        clique_matching_func(schedule, vertex_set, bs, local_ans, depth + 1);
        int *_data = vertex_set[prefix_id].get_data_ptr(), _size = vertex_set[prefix_id].get_size();
        for (int j = 0; j < _size; j++)
            bs->dec(_data[j]);
    }
}

long long Graph::pattern_matching_edge_task(const Schedule_IEP &schedule, int edge_id, VertexSet *vertex_sets, VertexSet &partial_embedding,
                                            VertexSet &tmp_set, int *ans_buffer) {
    int v0, v1, l, r;
    v0 = edge_from[edge_id];
    v1 = edge[edge_id];
    if (schedule.get_restrict_last(1) != -1 && v0 <= v1)
        return 0;

    l = vertex[v0], r = vertex[v0 + 1];
    for (int prefix_id = schedule.get_last(0); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
        vertex_sets[prefix_id].build_vertex_set(schedule, vertex_sets, &edge[l], r - l, prefix_id);

    l = vertex[v1], r = vertex[v1 + 1];
    for (int prefix_id = schedule.get_last(1); prefix_id != -1; prefix_id = schedule.get_next(prefix_id)) {
        vertex_sets[prefix_id].build_vertex_set(schedule, vertex_sets, &edge[l], r - l, prefix_id);
        if (vertex_sets[prefix_id].get_size() == 0 && prefix_id < schedule.get_basic_prefix_num())
            return 0;
    }

    partial_embedding.push_back(v0);
    partial_embedding.push_back(v1);
    long long ans = 0;
    pattern_matching_aggressive_func(schedule, vertex_sets, partial_embedding, tmp_set, ans, 2, ans_buffer);
    partial_embedding.pop_back();
    partial_embedding.pop_back();
    return ans;
}

// 进行实际模式匹配
// 输入调度方案，所有候选集，差集（不能被匹配的顶点），临时集合，
void Graph::pattern_matching_aggressive_func(const Schedule_IEP &schedule, VertexSet *vertex_set, VertexSet &subtraction_set, VertexSet &tmp_set,
                                             long long &local_ans, int depth, int *ans_buffer) {

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
        remove_anti_edge_vertices(diff_buf, vertex_set[loop_set_prefix_id], schedule, subtraction_set,
                                  depth); // 根据该候选集对于当前顶点进行移除？？？？？？？？？？？

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
            int min_vertex = v_cnt;
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
    int min_vertex = v_cnt;
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
        get_edge_index(vertex, l, r);

        // 根据当前匹配顶点，为后续每一个模式顶点重新构建候选集？  在构建边表的时候已经为后续每一个匹配施加限制了。
        bool is_zero = false;
        for (int prefix_id = schedule.get_last(depth); prefix_id != -1; prefix_id = schedule.get_next(prefix_id)) {
            vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &edge[l], (int)(r - l), prefix_id, vertex);
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
        pattern_matching_aggressive_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, depth + 1, ans_buffer); // 继续进行深度搜索
        subtraction_set.pop_back(); // 深度搜索处理完之后，弹出当前匹配顶点
    }
}

long long Graph::pattern_matching_mpi(const Schedule_IEP &schedule, int thread_count, bool clique) {
    Graphmpi &gm = Graphmpi::getinstance();
    long long global_ans = 0;
#pragma omp parallel num_threads(thread_count)
    {
#pragma omp master
        {
            gm.init(thread_count, this);
        }
#pragma omp barrier // mynodel have to be calculated before running other
                    // threads
#pragma omp master
        {
            global_ans = gm.runmajor();
        }
        if (omp_get_thread_num()) {
            int *ans_buffer = new int[schedule.in_exclusion_optimize_vertex_id.size()];
            VertexSet *vertex_set = new VertexSet[schedule.get_total_prefix_num()];
            long long local_ans = 0;
            VertexSet subtraction_set;
            VertexSet tmp_set;
            subtraction_set.init();
            auto match_start_vertex = [&](int vertex, int *data, int size) {
                for (int prefix_id = schedule.get_last(0); prefix_id != -1; prefix_id = schedule.get_next(prefix_id)) {
                    vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, data, size, prefix_id);
                }
                // subtraction_set.insert_ans_sort(vertex);
                subtraction_set.push_back(vertex);
                pattern_matching_aggressive_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, 1, ans_buffer);
                subtraction_set.pop_back();
            };
            for (std::pair<int, int> range;;) {
                if ((range = gm.get_vertex_range()).first == -1)
                    break;
                // for (int vertex = v_cnt - range.second; vertex < v_cnt -
                // range.first; vertex++) {//backwards slower than forwards
                for (int vertex = range.first; vertex < range.second; vertex++) {
                    e_index_t l, r;
                    get_edge_index(vertex, l, r);
                    match_start_vertex(vertex, edge + l, r - l);
                }
            }
            delete[] vertex_set;
            gm.report(local_ans);
            gm.end();
        }
    }
    return global_ans;
}

void Graph::pattern_matching_aggressive_func_mpi(const Schedule_IEP &schedule, VertexSet *vertex_set, VertexSet &subtraction_set, VertexSet &tmp_set,
                                                 long long &local_ans, int depth) {
    int loop_set_prefix_id = schedule.get_loop_set_prefix_id(depth);
    int loop_size = vertex_set[loop_set_prefix_id].get_size();
    if (loop_size <= 0)
        return;
    int *loop_data_ptr = vertex_set[loop_set_prefix_id].get_data_ptr();
    if (depth == schedule.get_size() - 1) {
        // TODO : try more kinds of calculation.
        // For example, we can maintain an ordered set, but it will cost more to
        // maintain itself when entering or exiting recursion.
        if (schedule.get_total_restrict_num() > 0) {
            int min_vertex = v_cnt;
            for (int i = schedule.get_restrict_last(depth); i != -1; i = schedule.get_restrict_next(i))
                if (min_vertex > subtraction_set.get_data(schedule.get_restrict_index(i)))
                    min_vertex = subtraction_set.get_data(schedule.get_restrict_index(i));
            const VertexSet &vset = vertex_set[loop_set_prefix_id];
            int size_after_restrict = std::lower_bound(vset.get_data_ptr(), vset.get_data_ptr() + vset.get_size(), min_vertex) - vset.get_data_ptr();
            if (size_after_restrict > 0)
                local_ans += VertexSet::unordered_subtraction_size(vertex_set[loop_set_prefix_id], subtraction_set, size_after_restrict);
        } else
            local_ans += VertexSet::unordered_subtraction_size(vertex_set[loop_set_prefix_id], subtraction_set);
        return;
    }

    // TODO : min_vertex is also a loop invariant
    int min_vertex = v_cnt;
    for (int i = schedule.get_restrict_last(depth); i != -1; i = schedule.get_restrict_next(i))
        if (min_vertex > subtraction_set.get_data(schedule.get_restrict_index(i)))
            min_vertex = subtraction_set.get_data(schedule.get_restrict_index(i));
    for (int i = 0; i < loop_size; ++i) {
        if (min_vertex <= loop_data_ptr[i])
            break;
        int vertex = loop_data_ptr[i];
        if (subtraction_set.has_data(vertex))
            continue;
        int *data, size;
        Graphmpi &gm = Graphmpi::getinstance();
        // if (gm.include(vertex)) {
        if (true) {
            e_index_t l, r;
            get_edge_index(vertex, l, r);
            data = edge + l;
            size = r - l;
        } else {
            data = gm.getneighbor(vertex);
            size = gm.getdegree();
        }
        bool is_zero = false;
        for (int prefix_id = schedule.get_last(depth); prefix_id != -1; prefix_id = schedule.get_next(prefix_id)) {
            vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, data, size, prefix_id, vertex);
            if (vertex_set[prefix_id].get_size() == 0) {
                is_zero = true;
                break;
            }
        }
        if (is_zero)
            continue;
        // subtraction_set.insert_ans_sort(vertex);
        subtraction_set.push_back(vertex);
        int *ans_buffer = new int[schedule.in_exclusion_optimize_vertex_id.size()];
        pattern_matching_aggressive_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, depth + 1, ans_buffer);
        subtraction_set.pop_back();
    }
}

// 对于每个顶点的邻接表，只保留比自己大的顶点
void erase_edge(Graph &g) {
    int newe = 0;
    printf("start erasing edge\n");
    fflush(stdout);

    for (int v = 0; v < g.v_cnt; v++) {
        int l = newe;
        for (long long e = g.vertex[v]; e < g.vertex[v + 1] && e < g.e_cnt; e++) {
            if (g.edge[e] >= v)
                continue;
            g.edge[newe] = g.edge[e];
            newe++;
        }
        g.vertex[v] = l;
    }
    g.vertex[g.v_cnt] = newe;
    g.e_cnt = newe;
    printf("newe: %ld\n", g.e_cnt);
}

// reduce two directed edge to one
void reduce_edges_for_clique(Graph &g) {
    printf("Trying to reduce edge. the pattern is a clique.\n");
    erase_edge(g);
    printf("Finish reduce.\n");
}

void Graph::get_third_layer_size(const Schedule_IEP &schedule, int *count) const {
    //    intersection_times_low = intersection_times_high = 0;
    //    dep1_cnt = dep2_cnt = dep3_cnt = 0;
    uint32_t *edge_from = new uint32_t[e_cnt];
    for (uint32_t i = 0; i < v_cnt; ++i)
        for (e_index_t j = vertex[i]; j < vertex[i + 1]; ++j)
            edge_from[j] = i;
#pragma omp parallel num_threads(64)
    {
        VertexSet *vertex_set = new VertexSet[schedule.get_total_prefix_num()];
#pragma omp for schedule(dynamic) nowait
        for (int e = 0; e < e_cnt; ++e) {
            count[e] = 0;

            int v0 = edge[e];
            int v1 = edge_from[e];
            e_index_t l, r;
            get_edge_index(v0, l, r);

            for (int prefix_id = schedule.get_last(0); prefix_id != -1; prefix_id = schedule.get_next(prefix_id)) {
                vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &edge[l], (int)(r - l), prefix_id);
            }
            get_edge_index(v1, l, r);
            for (int prefix_id = schedule.get_last(1); prefix_id != -1; prefix_id = schedule.get_next(prefix_id)) {
                vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &edge[l], r - l, prefix_id);
                count[e] += vertex_set[prefix_id].get_size();
            }
        }

        delete[] vertex_set;
    }
}

// direct reorder, which doesn't take the cache into account
void Graph::reorder_edge(const Schedule_IEP &schedule, e_index_t *new_order, e_index_t *task_start, int total_devices) const {
    int total_chunks = (e_cnt - 1) / chunk_size + 1;

    int *third_layer_count = new int[e_cnt];
    // get_third_layer_size(schedule, third_layer_count);

    int *chunk_size_sum = new int[total_chunks];
    for (int i = 0; i < total_chunks; i++)
        chunk_size_sum[i] = 0;
    for (int i = 0; i < e_cnt; i++)
        chunk_size_sum[i / chunk_size] += third_layer_count[i];

    int *chunk_order = new int[total_chunks];
    for (int i = 0; i < total_chunks; i++)
        chunk_order[i] = i;
    // std::sort(chunk_order, chunk_order + total_chunks, [chunk_size_sum](int i, int j){
    //     return chunk_size_sum[i] > chunk_size_sum[j] || (chunk_size_sum[i] == chunk_size_sum[j] && i < j);
    // });

    int *new_chunk_order = new int[total_chunks];
    for (int i = 0; i < total_chunks; i++) {
        new_chunk_order[chunk_order[i]] = i;
    }

    for (int i = 0; i <= total_devices; i++) {
        task_start[i] = 0;
    }
    for (int i = 0; i < total_chunks; i++) {
        task_start[new_chunk_order[i] % total_devices + 1] += (i == total_chunks - 1 ? ((e_cnt - 1) % chunk_size + 1) : chunk_size);
    }
    for (int i = 1; i <= total_devices; i++) {
        task_start[i] += task_start[i - 1];
    }

    for (int e = 0; e < e_cnt; e++)
        new_order[e] = e;
    // std::sort(new_order, new_order + e_cnt, [total_devices, new_chunk_order](int i, int j){
    //     // with third layer consideration
    //     // int i_chunk = new_chunk_order[i / chunk_size], j_chunk = new_chunk_order[j / chunk_size];
    //     // without third layer consideration
    //     int i_chunk = i / chunk_size, j_chunk = j / chunk_size;
    //     int i_devices = i_chunk % total_devices, j_devices = j_chunk % total_devices;
    //     return i_devices < j_devices || (i_devices == j_devices && i < j);
    // });
    delete[] third_layer_count;
    delete[] chunk_size_sum;
    delete[] chunk_order;
    delete[] new_chunk_order;
}

void degree_orientation_init(Graph *original_g, Graph *&g) {
    g = new Graph();

    int *degree;
    int n = original_g->v_cnt;
    long long m = original_g->e_cnt;
    g->v_cnt = n;
    g->e_cnt = m >> 1;
    g->edge = new int[g->e_cnt];
    g->vertex = new int64_t[g->v_cnt + 1];

    degree = new int[n];
    for (int u = 0; u < n; ++u) {
        degree[u] = original_g->vertex[u + 1] - original_g->vertex[u];
    }

    int64_t ptr = 0;
    int max_degree = 0;
    for (int u = 0; u < n; ++u) {
        g->vertex[u] = ptr;
        for (int64_t i = original_g->vertex[u]; i < original_g->vertex[u + 1]; ++i) {
            int v = original_g->edge[i];
            if ((degree[u] < degree[v]) || (degree[u] == degree[v] && u < v)) {
                g->edge[ptr++] = v;
            }
        }
        if (max_degree < ptr - g->vertex[u]) {
            max_degree = ptr - g->vertex[u];
        }
    }
    g->vertex[n] = ptr;
    assert(ptr == g->e_cnt);
    // partition_num = (max_degree - 1) / 32 + 1;
    // assert(max_degree < THREADS_PER_BLOCK);
    printf("Maximum degree after orientation: %d\n", max_degree);

    delete[] degree;
}

void degeneracy_orientation_init(Graph *original_g, Graph *&g) {
    g = new Graph();

    // 这里用 vector<int> node_with_degree[x] 来装当前度数为 x
    // 的所有节点，做到了线性求 k-core 序
    int *degree;
    int *order;
    int *vector_ptr;
    bool *in_vector;
    int n = original_g->v_cnt;
    long long m = original_g->e_cnt;
    std::vector<std::vector<int>> node_with_degree(n, std::vector<int>());
    g->v_cnt = n;
    g->e_cnt = m >> 1;
    g->edge = new int[g->e_cnt];
    g->vertex = new int64_t[g->v_cnt + 1];

    degree = new int[n];
    order = new int[n];
    in_vector = new bool[n];
    vector_ptr = new int[n];

    for (int u = 0; u < n; ++u) {
        in_vector[u] = true;
        degree[u] = original_g->vertex[u + 1] - original_g->vertex[u];
        node_with_degree[degree[u]].push_back(u);
    }

    for (int i = 0; i < n; ++i) {
        vector_ptr[i] = 0;
    }

    // 精巧的实现
    int order_ptr = 0;
    for (int current_min_degree = 0; current_min_degree < n; ++current_min_degree) {
        bool back = false;
        for (int &i = vector_ptr[current_min_degree]; i < node_with_degree[current_min_degree].size(); ++i) {
            int u = node_with_degree[current_min_degree][i];
            if (in_vector[u]) {
                order[u] = order_ptr++;
                in_vector[u] = false;
                for (int64_t j = original_g->vertex[u]; j < original_g->vertex[u + 1]; ++j) {
                    int v = original_g->edge[j];
                    if (in_vector[v]) {
                        node_with_degree[--degree[v]].push_back(v);
                        if (degree[v] == current_min_degree - 1) {
                            back = true;
                        }
                    }
                }
            }
        }
        if (back) {
            // 由于 for 循环里还要 +1，这里实际上只减了 1
            current_min_degree -= 2;
        }
    }

    int64_t ptr = 0;
    int max_degree = 0;
    for (int u = 0; u < n; ++u) {
        g->vertex[u] = ptr;
        for (int64_t i = original_g->vertex[u]; i < original_g->vertex[u + 1]; ++i) {
            int v = original_g->edge[i];
            if (order[u] < order[v]) {
                g->edge[ptr++] = v;
            }
        }
        if (max_degree < ptr - g->vertex[u]) {
            max_degree = ptr - g->vertex[u];
        }
    }
    g->vertex[n] = ptr;
    assert(ptr == g->e_cnt);
    // partition_num = (max_degree - 1) / 32 + 1;
    printf("Maximum vertex degree after orientation: %d\n", max_degree);

    delete[] degree;
    delete[] order;
    delete[] in_vector;
    delete[] vector_ptr;
}

// 对于motif的测试，生成motif后，生成IEP，进行匹配
void Graph::motif_counting(int pattern_size) {

    double total_counting_time = 0;
    TimeInterval allTime, tmpTime;
    allTime.check();

    MotifGenerator mg(pattern_size);
    std::vector<Pattern> motifs = mg.generate();
    for (int i = 0; i < motifs.size(); ++i) {
        Pattern p = motifs[i];

        printf("pattern = \n");
        p.print();

        printf("max intersection size %d\n", VertexSet::max_intersection_size);
        bool is_pattern_valid;
        bool use_in_exclusion_optimize = true;

        // 生成调度方案 iep，限制，性能模型
        // 输入模式，模式可行标签，  顶点数，边数，三角形数。
        Schedule_IEP schedule_iep(p, is_pattern_valid, 1, 1, use_in_exclusion_optimize, v_cnt, e_cnt, tri_cnt);
        // Schedule_IEP schedule_iep(p, is_pattern_valid, 1, 1, use_in_exclusion_optimize, 41652230, 2936730364, 34824916864);
        if (!is_pattern_valid) {
            printf("pattern is invalid!\n");
            continue;
        }

        tmpTime.check();
        // 需要了解，怎么进行匹配
        long long ans = this->pattern_matching(schedule_iep);
        total_counting_time += tmpTime.print("Pattern counting time: ");

        printf("ans: %lld\n", ans);
    }

    printf("Counting time cost: %.6lf s\n", total_counting_time);
    allTime.print("Total time cost");
}

// 以边为中心，一套代码同时测试楔形和三角形

int32_t get_intersection_size(const v_index_t *a, int32_t na, const v_index_t *b, int32_t nb) {
    int32_t ans = 0;
    int32_t i = 0, j = 0;
    while (i < na) {
        if (i >= 1) {
            assert(a[i] > a[i - 1]);
        }
        while (j < nb && b[j] < a[i])
            j++;
        if (j == nb)
            break;
        if (a[i] == b[j])
            ans++;
        i++;
    }
    return ans;
}

// 对于大小为3的motif进行计数。
void Graph::motif_counting_3() {
    int64_t tri_cnt = 0, wedge_cnt = 0;

    // 记录每条边的起点
    uint32_t *edge_from = new uint32_t[e_cnt];
    for (uint32_t i = 0; i < v_cnt; ++i)
        for (e_index_t j = vertex[i]; j < vertex[i + 1]; ++j)
            edge_from[j] = i;

    auto t1 = std::chrono::system_clock::now();

    // 对于每一条边进行处理
#pragma omp parallel for schedule(dynamic) reduction(+ : tri_cnt) reduction(+ : wedge_cnt)
    for (e_index_t i = 0; i < e_cnt; i++) {
        v_index_t v0 = edge_from[i], v1 = edge[i];
        e_index_t l0 = vertex[v0], r0 = vertex[v0 + 1]; // 左侧端点的边表的两个指针

        // 直接计算与左侧顶点相关的楔形数量
        if (i == l0) {
            e_index_t d = r0 - l0;
            wedge_cnt += d * (d - 1) / 2;
        }

        // 避免重复计数
        if (v0 <= v1)
            continue;

        // 右侧顶点边表指针
        e_index_t l1 = vertex[v1], r1 = vertex[v1 + 1];

        // 通过两个顶点的邻居边表求交集来记录三角形的数量
        tri_cnt += get_intersection_size(&edge[l0], r0 - l0, &edge[l1], r1 - l1);
    }

    wedge_cnt -= tri_cnt;

    tri_cnt /= 3;
    auto t2 = std::chrono::system_clock::now();
    double time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() * 1e-6;

    printf("triangle: %ld wedge: %ld time:%.6lf\n", tri_cnt, wedge_cnt, time);
}

void Graph::triangle_statistics() {
    int64_t tri_cnt = 0, wedge_cnt = 0;

    // 记录每条边的起点
    uint32_t *edge_from = new uint32_t[e_cnt];
    for (uint32_t i = 0; i < v_cnt; ++i)
        for (e_index_t j = vertex[i]; j < vertex[i + 1]; ++j)
            edge_from[j] = i;

    auto t1 = std::chrono::system_clock::now();

    std::ofstream of;
    of.open("triangle_statistics.txt");

    uint32_t *count = new uint32_t[5000];
    memset(count, 0, sizeof(uint32_t) * 5000);

    int64_t max_count = 0;
    // 对于每一条边进行处理
    // #pragma omp parallel for schedule(dynamic) reduction(+ : tri_cnt) reduction(+ : wedge_cnt)
    for (e_index_t i = 0; i < e_cnt; i++) {
        v_index_t v0 = edge_from[i], v1 = edge[i];
        e_index_t l0 = vertex[v0], r0 = vertex[v0 + 1]; // 左侧端点的边表的两个指针

        // 直接计算与左侧顶点相关的楔形数量
        if (i == l0) {
            e_index_t d = r0 - l0;
            wedge_cnt += d * (d - 1) / 2;
        }

        // 避免重复计数
        if (v0 <= v1)
            continue;

        // 右侧顶点边表指针
        e_index_t l1 = vertex[v1], r1 = vertex[v1 + 1];

        // 通过两个顶点的邻居边表求交集来记录三角形的数量
        int64_t num = get_intersection_size(&edge[l0], r0 - l0, &edge[l1], r1 - l1);
        // of << num << std::endl;
        // if (max_count < num) {
        //     max_count = num;
        // }
        count[num]++;
        tri_cnt += num;
    }

    for (int i = 0; i < 5000; i++) {
        of << count[i] << std::endl;
    }

    of.close();

    wedge_cnt -= tri_cnt;

    tri_cnt /= 3;
    auto t2 = std::chrono::system_clock::now();
    double time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() * 1e-6;

    printf("triangle: %ld wedge: %ld time:%.6lf\n", tri_cnt, wedge_cnt, time);
}

void static_get_intersection_set(std::vector<u_int32_t> &candidate_set, const v_index_t *a, int32_t na, const v_index_t *b, int32_t nb) {
    int32_t i = 0, j = 0;
    while (i < na) { // 边表是有序的，一次处理
        if (i >= 1) {
            assert(a[i] > a[i - 1]);
        }
        while (j < nb && b[j] < a[i])
            j++;
        if (j == nb)
            break;
        if (a[i] == b[j])
            candidate_set.push_back(a[i]);
        i++;
    }
}

u_int64_t static_get_intersection_size(std::vector<u_int32_t> &candidate_set, const v_index_t *a, int32_t na) {
    int32_t i = 0, j = 0;
    u_int64_t ans = 0;
    while (i < na) {
        if (i >= 1) {
            assert(a[i] > a[i - 1]);
        }
        while (j < candidate_set.size() && candidate_set[j] < a[i])
            j++;
        if (j == candidate_set.size())
            break;
        if (a[i] == candidate_set[j])
            ans++;
        i++;
    }
    return ans;
}

void Graph::four_clique_statistics() {

    int64_t tri_cnt = 0, wedge_cnt = 0;

    // 记录每条边的起点
    uint32_t *edge_from = new uint32_t[e_cnt];
    for (uint32_t i = 0; i < v_cnt; ++i)
        for (e_index_t j = vertex[i]; j < vertex[i + 1]; ++j)
            edge_from[j] = i;

    auto t1 = std::chrono::system_clock::now();

    std::ofstream of;
    of.open("four_clique_statistics.txt");

    std::vector<uint32_t> third_candidate_set;

    // 随机生成器
    std::random_device rd;
    std::mt19937 gen = std::mt19937(rd());

    uint32_t *count = new uint32_t[5000];
    memset(count, 0, sizeof(uint32_t) * 5000);

    int64_t max_count = 0;
    // 对于每一条边进行处理
    // #pragma omp parallel for schedule(dynamic) reduction(+ : tri_cnt) reduction(+ : wedge_cnt)
    for (e_index_t i = 0; i < e_cnt; i++) {
        v_index_t v0 = edge_from[i], v1 = edge[i];
        e_index_t l0 = vertex[v0], r0 = vertex[v0 + 1];
        e_index_t l1 = vertex[v1], r1 = vertex[v1 + 1];

        third_candidate_set.clear();

        // 避免选取到重复的边
        if (v0 <= v1)
            continue;

        static_get_intersection_set(third_candidate_set, &edge[l0], r0 - l0, &edge[l1], r1 - l1);

        if (third_candidate_set.size() <= 1) {
            count[0]++;
            continue;
        }

        std::uniform_int_distribution<> dis(1, third_candidate_set.size());
        v_index_t v2 = third_candidate_set[dis(gen) - 1];
        e_index_t l2 = vertex[v2], r2 = vertex[v2 + 1];

        u_int64_t tempsize = static_get_intersection_size(third_candidate_set, &edge[l2], r2 - l2);

        u_int64_t num = tempsize * third_candidate_set.size();

        // if (num == 2) {
        //     std::cout << std::endl;
        // }
        // of << num << std::endl;
        // if (max_count < num) {
        //     max_count = num;
        // }
        if (num < 5000)
            count[num]++;
    }

    for (int i = 0; i < 5000; i++) {
        of << count[i] << std::endl;
    }

    of.close();

    wedge_cnt -= tri_cnt;

    tri_cnt /= 3;
    auto t2 = std::chrono::system_clock::now();
    double time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() * 1e-6;

    printf("triangle: %ld wedge: %ld time:%.6lf\n", tri_cnt, wedge_cnt, time);
}