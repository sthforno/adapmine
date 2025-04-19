#include "../include/schedule_IEP.h"
#include "../include/dataloader.h"
#include "../include/graph.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <stdexcept>

// 输入模式，顶点数，边数，三角形数 最后一个默认为false
Schedule_IEP::Schedule_IEP(const Pattern &pattern, bool &is_pattern_valid, int performance_modeling_type, int restricts_type,
                           bool use_in_exclusion_optimize, int v_cnt, e_index_t e_cnt, long long tri_cnt, bool vertex_induced) {
    // if (!use_in_exclusion_optimize) {
    //     throw std::logic_error("Schedule_IEP: must set use_in_exclusion_optimize to true.\n");
    // }

    // 如果是顶点诱导的，那么不能使用iep，记录具有iep冗余 这里默认不是顶点诱导的
    is_vertex_induced = vertex_induced;
    if (vertex_induced) {
        use_in_exclusion_optimize = false;
        in_exclusion_optimize_redundancy = 1;
    }

    // 性能模型需要使用到三角形数量
    if (performance_modeling_type == 1 && tri_cnt == -1) {
        throw std::logic_error("Fatal: Can not use performance modeling if not have triangle number of this dataset.\n");
    }

    // 构建adj_mat矩阵存储模式中存在的边
    is_pattern_valid = true;
    size = pattern.get_size();
    adj_mat = new int[size * size];
    memcpy(adj_mat, pattern.get_adj_mat_ptr(), size * size * sizeof(int));

    // The I-th loop consists of at most the intersection of i-1 VertexSet.
    // So the max number of prefix = 0 + 1 + ... + size-1 = size * (size-1) / 2

    // 最多可能存在的边数。
    int max_prefix_num = size * (size - 1) / 2 + 1;
    father_prefix_id = new int[max_prefix_num]; // 模式边相关，推测用于做限制            一共九个数组

    // 用于处理某个顶点及其所有相邻的边。
    last = new int[size];           // 模式顶点相关
    next = new int[max_prefix_num]; // 模式边相关。

    // 用于剪枝，如果某个候选集的大小 小于需要的某个大小的话
    break_size = new int[max_prefix_num]; // 模式边相关，

    // 应该是记录顶点的处理顺序。索引是深度，记录的是应该处理的顶点
    loop_set_prefix_id = new int[size]; // 模式顶点相关

    prefix = new Prefix[max_prefix_num]; // 模式边相关

    // 与约束相关的有三个，用于打破对称性。
    // 一个记录对于某个顶点的第一个约束的位置，一个记录当前约束下一个约束的位置，一个记录实际的约束是什么（约束当前顶点的顶点都大于当前顶点）。
    // 约束就是有向的边，因此约束的数量与边相关
    restrict_last = new int[size];            // 模式顶点相关
    restrict_next = new int[max_prefix_num];  // 模式边相关
    restrict_index = new int[max_prefix_num]; // 模式边相关

    memset(father_prefix_id, -1, max_prefix_num * sizeof(int));
    memset(loop_set_prefix_id, -1, sizeof(int) * size);
    memset(last, -1, size * sizeof(int));
    memset(next, -1, max_prefix_num * sizeof(int));
    memset(restrict_next, -1, max_prefix_num * sizeof(int));
    memset(break_size, -1, max_prefix_num * sizeof(int));
    memset(restrict_last, -1, size * sizeof(int));
    memset(restrict_next, -1, max_prefix_num * sizeof(int));
    memset(restrict_next, -1, max_prefix_num * sizeof(int));

    std::vector<std::pair<int, int>> best_pairs;
    best_pairs.clear();
    // Initialize adj_mat
    // If we use performance_modeling, we may change the order of vertex,
    // the best order produced by performance_modeling(...) is saved in best_order[]
    // Finally, we use best_order[] to relocate adj_mat
    //
    // 如果使用性能模型，可能会需要该边模式中顶点的顺序
    if (performance_modeling_type != 0) {

        // 生成所有可能的顶点处理顺序，与模式无关。然后根据模式移除不合法的处理顺序（后续处理的顶点，必然与至少前面一个顶点在模式中相连）
        unsigned int pow = 1;
        for (int i = 2; i <= size; ++i)
            pow *= i;

        std::vector<std::vector<int>> candidate_permutations;
        candidate_permutations.clear();

        bool use[size];
        for (int i = 0; i < size; ++i)
            use[i] = false;
        std::vector<int> tmp_vec;
        get_full_permutation(candidate_permutations, use, tmp_vec, 0);
        assert(candidate_permutations.size() == pow);

        remove_invalid_permutation(candidate_permutations);

        // 如果性能模型等于1，进一步筛选调度方案（保留底层不相交的的顶点最多的一系列模式）。
        if (performance_modeling_type == 1) {
            int max_val = 0;
            for (const auto &vec : candidate_permutations) {
                max_val = std::max(max_val, get_vec_optimize_num(vec));
            }
            k_val = max_val;
            std::vector<std::vector<int>> tmp;
            tmp.clear();
            for (const auto &vec : candidate_permutations)
                if (get_vec_optimize_num(vec) == max_val) {
                    tmp.push_back(vec);
                }
            candidate_permutations = tmp;
        }

        std::vector<int> best_order(size);
        double min_val = 1e18;
        bool have_best = false;
        const int *pattern_adj_mat = pattern.get_adj_mat_ptr();

        // printf("number of candidate permutations: %ld\n", candidate_permutations.size());

        // 对于所有合法的，并且具有最大IEP优化潜力的所有调度进行处理，判断增加了约束之后哪些最好。
        for (const auto &vec : candidate_permutations) {

            // 根据调度顺序重新排列矩阵
            copy_adj_mat_from(vec, pattern_adj_mat);

            setup_optimization_info(use_in_exclusion_optimize);
            std::vector<std::vector<std::pair<int, int>>> restricts_vector;
            restricts_vector.clear();

            // 根据模式矩阵，生成所有可能的约束
            if (restricts_type == 1) {
                restricts_generate(adj_mat, restricts_vector);
            } else {
                Schedule_IEP schedule(adj_mat, size);

                std::vector<std::pair<int, int>> pairs;
                schedule.GraphZero_aggressive_optimize(pairs);

                restricts_vector.clear();
                restricts_vector.push_back(pairs);
            }

            // 如果约束大小为0，说明这个模式不具有自同构性质，直接预测性能
            if (restricts_vector.size() == 0) {
                std::vector<std::pair<int, int>> Empty;
                Empty.clear();

                double val = 1e18;
                if (performance_modeling_type == 1) {
                    val = new_estimate_schedule_restrict(Empty, v_cnt, e_cnt, tri_cnt);
                }
                // if(performance_modeling_type == 1) {
                //     val = our_estimate_schedule_restrict(vec, Empty, v_cnt, e_cnt, tri_cnt);
                // }
                // else {
                //     if(performance_modeling_type == 2) {
                //         val = GraphZero_estimate_schedule_restrict(vec, Empty, v_cnt, e_cnt);
                //     }
                //     else {
                //         val = Naive_estimate_schedule_restrict(vec, Empty, v_cnt, e_cnt);
                //     }
                // }

                if (have_best == false || val < min_val) {
                    have_best = true;
                    min_val = val;
                    for (int i = 0; i < size; ++i)
                        best_order[i] = vec[i];
                    best_pairs = Empty;
                }
            }

            // 处理所有的约束，预测性能
            for (const auto &pairs : restricts_vector) {
                double val = 1e18;
                if (performance_modeling_type == 1) {
                    val = new_estimate_schedule_restrict(pairs, v_cnt, e_cnt, tri_cnt);
                }
                // if(performance_modeling_type == 1) {
                //     val = our_estimate_schedule_restrict(vec, pairs, v_cnt, e_cnt, tri_cnt);
                // }
                // else {
                //     if(performance_modeling_type == 2) {
                //         val = GraphZero_estimate_schedule_restrict(vec, pairs, v_cnt, e_cnt);
                //     }
                //     else {
                //         val = Naive_estimate_schedule_restrict(vec, pairs, v_cnt, e_cnt);
                //     }
                // }

                if (have_best == false || val < min_val) {
                    have_best = true;
                    min_val = val;
                    for (int i = 0; i < size; ++i)
                        best_order[i] = vec[i];
                    best_pairs = pairs;
                }
            }
        }

        // 最后保存最好的处理顺序
        copy_adj_mat_from(best_order, pattern_adj_mat);

        // printf("best order:\n");
        // for(auto i: best_order) {
        //     printf("%d ",i);
        // }
        // printf("\n");
    } else {
        // std::vector< int > I;
        // I.clear();
        // for(int i = 0; i < size; ++i) I.push_back(i);

        // 构建一个限制向量
        std::vector<std::vector<std::pair<int, int>>> restricts_vector;
        restricts_vector.clear();

        // 根据限制类型生成不同的限制组
        if (restricts_type != 0) {
            if (restricts_type == 1) {
                restricts_generate(adj_mat, restricts_vector);
            } else {
                std::vector<std::pair<int, int>> pairs;
                GraphZero_aggressive_optimize(pairs);

                restricts_vector.clear();
                restricts_vector.push_back(pairs);
            }
        }

        bool have_best = false;
        double min_val;

        // 处理每一组约束，寻找最优的约束方法
        for (const auto &pairs : restricts_vector) {
            double val;
            // if(restricts_type == 1) {
            //     val = our_estimate_schedule_restrict(I, pairs, v_cnt, e_cnt, tri_cnt);
            // }
            // else {
            //     val = GraphZero_estimate_schedule_restrict(I, pairs, v_cnt, e_cnt);
            // }
            if (restricts_type) {
                val = new_estimate_schedule_restrict(pairs, v_cnt, e_cnt, tri_cnt);
            }

            if (have_best == false || val < min_val) {
                have_best = true;
                min_val = val;
                best_pairs = pairs;
            }
        }
    }

    // 检查模式是否是连接的
    if (!check_connectivity()) {
        is_pattern_valid = false;
        throw std::runtime_error("pattern is not connected");
    }

    // 记录是否使用了iep优化
    setup_optimization_info(use_in_exclusion_optimize);

    // 记录更好的限制
    if (restricts_type != 0)
        add_restrict(best_pairs);

    // 记录是否存在iep冗余？？
    set_in_exclusion_optimize_redundancy();

    // print_schedule();
}

Schedule_IEP::Schedule_IEP(const int *_adj_mat, int _size) {
    size = _size;
    adj_mat = new int[size * size];

    memcpy(adj_mat, _adj_mat, size * size * sizeof(int));

    // The I-th loop consists of at most the intersection of i-1 VertexSet.
    // So the max number of prefix = 0 + 1 + ... + size-1 = size * (size-1) / 2
    int max_prefix_num = size * (size - 1) / 2 + 1;
    father_prefix_id = new int[max_prefix_num];
    last = new int[size];
    next = new int[max_prefix_num];
    break_size = new int[max_prefix_num];
    loop_set_prefix_id = new int[size];
    prefix = new Prefix[max_prefix_num];
    restrict_last = new int[size];
    restrict_next = new int[max_prefix_num];
    restrict_index = new int[max_prefix_num];
    memset(father_prefix_id, -1, max_prefix_num * sizeof(int));
    memset(last, -1, size * sizeof(int));
    memset(loop_set_prefix_id, -1, size * sizeof(int));
    memset(next, -1, max_prefix_num * sizeof(int));
    memset(break_size, -1, max_prefix_num * sizeof(int));
    memset(restrict_last, -1, size * sizeof(int));
    memset(restrict_next, -1, max_prefix_num * sizeof(int));
    memset(restrict_index, -1, max_prefix_num * sizeof(int));

    total_prefix_num = 0;
    total_restrict_num = 0;
    in_exclusion_optimize_num = 0;
    is_vertex_induced = true;

    if (!check_connectivity())
        throw std::runtime_error("pattern is not connected");

    build_loop_invariant();

    set_in_exclusion_optimize_redundancy();

    // print_schedule();
}

Schedule_IEP::~Schedule_IEP() {
    delete[] adj_mat;
    delete[] father_prefix_id;
    delete[] last;
    delete[] next;
    delete[] loop_set_prefix_id;
    delete[] prefix;
    delete[] restrict_last;
    delete[] restrict_next;
    delete[] restrict_index;
}

// note: this function no longer takes `order` as a parameter, instead, it uses `this->adj_mat` directly

// 性能预测模型
double Schedule_IEP::new_estimate_schedule_restrict(const std::vector<std::pair<int, int>> &pairs, int v_cnt, e_index_t e_cnt, long long tri_cnt) {
    int max_degree = get_max_degree();

    double p_size[max_degree];
    double pp_size[max_degree];

    double p0 = e_cnt * 1.0 / v_cnt / v_cnt;           // 边存在的概率
    double p1 = tri_cnt * 1.0 * v_cnt / e_cnt / e_cnt; // 三角形数量乘以顶点数除以楔形的数量

    // 访问到第i层时，第i层的大小
    p_size[0] = v_cnt;
    for (int i = 1; i < max_degree; ++i) {
        p_size[i] = p_size[i - 1] * p0;
    }
    // 相关的每一层的三角形的数量？
    pp_size[0] = 1;
    for (int i = 1; i < max_degree; ++i) {
        pp_size[i] = pp_size[i - 1] * p1;
    }

    // 初始化
    std::vector<std::pair<int, int>> restricts = pairs;
    int restricts_size = restricts.size();
    std::sort(restricts.begin(), restricts.end());

    // 结果
    double sum[restricts_size];
    for (int i = 0; i < restricts_size; ++i)
        sum[i] = 0;

    // 更新sum数组
    int tmp[size];
    for (int i = 0; i < size; ++i)
        tmp[i] = i;
    do {
        for (int i = 0; i < restricts_size; ++i)
            if (tmp[restricts[i].first] > tmp[restricts[i].second]) {
                sum[i] += 1;
            } else
                break;
    } while (std::next_permutation(tmp, tmp + size));

    // 计算sum
    double total = 1;
    for (int i = 2; i <= size; ++i)
        total *= i;
    for (int i = 0; i < restricts_size; ++i)
        sum[i] = sum[i] / total;
    for (int i = restricts_size - 1; i > 0; --i)
        sum[i] /= sum[i - 1];

    // 初始化invariant的大小
    std::vector<int> invariant_size[size];
    for (int i = 0; i < size; ++i)
        invariant_size[i].clear();

    // 性能预测
    double val = 1;

    // 从有相互依赖的顶点开始处理
    for (int i = size - in_exclusion_optimize_num - 1; i >= 0; --i) { // 每一层往上走

        // 与该顶点相关的边数
        int cnt_forward = 0;
        for (int j = 0; j < i; ++j)
            if (adj_mat[INDEX(j, i, size)])
                ++cnt_forward;

        // 记录每条边它是第几个
        int c = cnt_forward;
        for (int j = i - 1; j >= 0; --j)
            if (adj_mat[INDEX(j, i, size)])
                invariant_size[j].push_back(c--);

        // 计算每条边的时间消耗
        for (int j = 0; j < invariant_size[i].size(); ++j)
            if (invariant_size[i][j] > 1)
                val += p_size[1] * pp_size[invariant_size[i][j] - 2] * std::log2(p_size[1]); // 边存在概率 * 边相关的三角形数量 * 概率存在的对数
        // val += p_size[1] * pp_size[invariant_size[i][j] - 2] + p_size[1];
        // val += 1;

        // 乘以约束减小的时间消耗
        for (int j = 0; j < restricts_size; ++j)
            if (restricts[j].second == i)
                val *= sum[j];

        if (i) {
            val *= p_size[1] * pp_size[cnt_forward - 1];
        } else {
            val *= p_size[0];
        }
    }

    return val;
}

void Schedule_IEP::copy_adj_mat_from(const std::vector<int> &vec, const int *src_adj_mat) {
    int rank[size];
    for (int i = 0; i < size; ++i)
        rank[vec[i]] = i;
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            adj_mat[INDEX(rank[i], rank[j], size)] = src_adj_mat[INDEX(i, j, size)];
}

void Schedule_IEP::update_loop_invariant_for_fsm() {
    int max_prefix_num = size * (size - 1) / 2;
    memset(father_prefix_id, -1, max_prefix_num * sizeof(int));
    memset(last, -1, size * sizeof(int));
    memset(next, -1, max_prefix_num * sizeof(int));
    total_prefix_num = 0;
    prefix_target = new int[max_prefix_num];
    memset(prefix_target, -1, max_prefix_num * sizeof(int));

    loop_set_prefix_id[0] = -1;
    for (int i = 1; i < size;
         ++i) // 为每个loop_set独立地建立prefix，因为每个点的label不同（暂不考虑相同的情况），所以不同点的loop_set的prefix没有关系
    {
        int last_prefix = -1;
        for (int j = 0; j < i; ++j)
            if (adj_mat[INDEX(i, j, size)]) {
                int father = last_prefix;
                father_prefix_id[total_prefix_num] = father;
                next[total_prefix_num] = last[j];
                last[j] = total_prefix_num;
                prefix_target[total_prefix_num] = i;
                last_prefix = total_prefix_num++;
            }
        loop_set_prefix_id[i] = last_prefix;
    }
    assert(total_prefix_num <= max_prefix_num);
}

void Schedule_IEP::setup_optimization_info(bool use_in_exclusion_optimize) {
    int max_prefix_num = size * (size - 1) / 2;
    memset(father_prefix_id, -1, max_prefix_num * sizeof(int));
    memset(last, -1, size * sizeof(int));
    memset(next, -1, max_prefix_num * sizeof(int));
    memset(break_size, -1, max_prefix_num * sizeof(int));
    memset(restrict_last, -1, size * sizeof(int));
    memset(restrict_next, -1, max_prefix_num * sizeof(int));

    total_prefix_num = 0;
    total_restrict_num = 0;

    if (use_in_exclusion_optimize) {
        std::vector<int> I;
        I.clear();
        for (int i = 0; i < size; ++i)
            I.push_back(i);
        in_exclusion_optimize_num = get_vec_optimize_num(I);
        if (in_exclusion_optimize_num <= 1) {
            // printf("Can not use in_exclusion_optimize with this schedule\n");
            in_exclusion_optimize_num = 0;
        } else {
            // printf("use in_exclusion_optimize with size %d\n", in_exclusion_optimize_num);
            init_in_exclusion_optimize();
        }
    } else {
        in_exclusion_optimize_num = 0;
    }

    build_loop_invariant(in_exclusion_optimize_num);
}

bool Schedule_IEP::check_connectivity() const {
    // The I-th vertex must connect with at least one vertex from 0 to i-1.
    for (int i = 1; i < size; ++i) {
        bool valid = false;
        for (int j = 0; j < i; ++j) {
            if (adj_mat[INDEX(i, j, size)]) {
                valid = true;
                break;
            }
        }
        if (!valid)
            return false;
    }
    return true;
}

int Schedule_IEP::get_max_degree() const {
    int mx = 0;
    for (int i = 0; i < size; ++i) {
        int cnt = 0;
        for (int j = 0; j < size; ++j)
            cnt += adj_mat[INDEX(i, j, size)];
        if (cnt > mx)
            mx = cnt;
    }
    return mx;
}

void Schedule_IEP::build_loop_invariant(int in_exclusion_optimize_num) {

    // 为每个模式顶点记录前缀顶点
    int *tmp_data = new int[size];
    loop_set_prefix_id[0] = -1;
    for (int i = 1; i < size; ++i) {
        int data_size = 0;
        for (int j = 0; j < i; ++j)
            if (adj_mat[INDEX(i, j, size)])
                tmp_data[data_size++] = j; // 记录所有与当前顶点i相关的顶点
        loop_set_prefix_id[i] = find_father_prefix(data_size, tmp_data);
    }
    basic_prefix_num = total_prefix_num;

    for (int prefix_id = 0; prefix_id < basic_prefix_num; ++prefix_id) {
        const int *data = prefix[prefix_id].get_data_ptr();
        int data_size = prefix[prefix_id].get_size();

        int full_connect_cnt = 0;
        for (int node = 0; node < data[data_size - 1]; ++node) {
            bool in_data = false;
            bool not_connect = false;
            for (int i = 0; i < data_size; ++i) {
                if (node == data[i]) {
                    in_data = true;
                    break;
                }
                if (!adj_mat[INDEX(node, data[i], size)]) {
                    not_connect = true;
                    break;
                }
            }
            if (in_data == false && not_connect == false)
                ++full_connect_cnt;
        }
        break_size[prefix_id] = full_connect_cnt;
    }

    if (in_exclusion_optimize_num > 0) {
        // printf("begin to build IEP loop invariant, basic prefix num = %d\n", basic_prefix_num);
        // IEP loop invariant
        in_exclusion_optimize_vertex_id.clear();
        in_exclusion_optimize_coef.clear();
        in_exclusion_optimize_flag.clear();

        for (int optimize_rank = 0; optimize_rank < in_exclusion_optimize_group.size(); ++optimize_rank) {
            const std::vector<std::vector<int>> &cur_graph = in_exclusion_optimize_group[optimize_rank];
            long long val = in_exclusion_optimize_val[optimize_rank];

            for (int cur_graph_rank = 0; cur_graph_rank < cur_graph.size(); ++cur_graph_rank) {

                int data_size = 0;
                for (int node = 0; node < size - in_exclusion_optimize_num; ++node) {
                    for (int i = 0; i < cur_graph[cur_graph_rank].size(); ++i)
                        if (adj_mat[INDEX(node, cur_graph[cur_graph_rank][i] + size - in_exclusion_optimize_num, size)]) {
                            tmp_data[data_size++] = node;
                            break;
                        }
                }

                int my_id = find_father_prefix(data_size, tmp_data);

                int equal = -1;
                for (int i = 0; i < in_exclusion_optimize_vertex_id.size(); ++i)
                    if (in_exclusion_optimize_vertex_id[i] == my_id) {
                        equal = i;
                        break;
                    }
                if (equal == -1) {
                    in_exclusion_optimize_vertex_id.push_back(my_id);
                    equal = in_exclusion_optimize_vertex_id.size() - 1;

                    int in_set = 0;
                    int full_connect = 0;
                    for (int i = 0; i < size - in_exclusion_optimize_num; ++i) {
                        int cnt = 0;
                        bool hit = false;
                        for (int j = 0; j < data_size; ++j) {
                            if (tmp_data[j] == i) {
                                hit = true;
                                break;
                            }
                            if (adj_mat[INDEX(i, tmp_data[j], size)])
                                ++cnt;
                        }
                        if (hit)
                            ++in_set;
                        else {
                            if (cnt == data_size)
                                ++full_connect;
                        }
                    }
                    if (in_set + full_connect == size - in_exclusion_optimize_num) {
                        in_exclusion_optimize_vertex_flag.push_back(true);
                        in_exclusion_optimize_vertex_coef.push_back(full_connect);
                    } else {
                        in_exclusion_optimize_vertex_flag.push_back(false);
                        in_exclusion_optimize_vertex_coef.push_back(full_connect);
                    }
                }

                in_exclusion_optimize_ans_pos.push_back(equal);

                if (cur_graph_rank == cur_graph.size() - 1) {
                    in_exclusion_optimize_coef.push_back(val);
                    in_exclusion_optimize_flag.push_back(true);
                } else {
                    in_exclusion_optimize_coef.push_back(0);
                    in_exclusion_optimize_flag.push_back(false);
                }
            }
        }

        for (int i = 0; i < in_exclusion_optimize_vertex_flag.size(); ++i) {
            if (in_exclusion_optimize_vertex_flag[i] && i < in_exclusion_optimize_vertex_id.size()) {
                int prefix_id = in_exclusion_optimize_vertex_id[i];
                // printf("prefix_id=%d i=%d in_exclusion_optimize_vertex_id.size()=%ld\n", prefix_id, i, in_exclusion_optimize_vertex_id.size());
                if (prefix[prefix_id].get_has_child() == false)
                    prefix[prefix_id].set_only_need_size(true);
            }

            // printf("total prefix num = %d\n", total_prefix_num);
        }
    }

    // last和next都用于记录一系列的顶点
    for (int i = 0; i < size; ++i)
        if (last[i] != -1) {
            int *buf = new int[size];
            int cnt = 0;

            for (int id = last[i]; id != -1; id = next[id]) {
                buf[cnt++] = id;
            }

            last[i] = buf[cnt - 1];
            for (int id = cnt - 1; id > 0; --id) {
                next[buf[id]] = buf[id - 1];
            }
            next[buf[0]] = -1;

            delete[] buf;
        }

    delete[] tmp_data;
}
// 输入的是某个模式顶点的边表，模式图以有向图的形式存储，因此存在顶点没有边表
int Schedule_IEP::find_father_prefix(int data_size, const int *data) {
    if (data_size == 0)
        return -1;
    int num = data[data_size - 1];                                                // 相关的最后一个顶点
    for (int prefix_id = last[num]; prefix_id != -1; prefix_id = next[prefix_id]) // 对于这个顶点所有的前缀顶点进行处理
        if (prefix[prefix_id].equal(data_size, data))                             // 如果这个前缀顶点与当前顶点相同，那么直接放回。
            return prefix_id;

    // not found, create new prefix and find its father prefix id recursively
    int father = find_father_prefix(data_size - 1, data);
    prefix[father].set_has_child(true);
    father_prefix_id[total_prefix_num] = father;
    next[total_prefix_num] = last[num];
    last[num] = total_prefix_num;
    prefix[total_prefix_num].init(data_size, data);

    ++total_prefix_num;
    return total_prefix_num - 1;
}

// 记录约束
void Schedule_IEP::add_restrict(const std::vector<std::pair<int, int>> &restricts) {
    restrict_pair = restricts;
    for (unsigned int i = 0; i < restrict_pair.size();) {
        bool tag = true;
        for (unsigned int j = 0; j < restrict_pair.size(); ++j) {
            if (i != j && restrict_pair[j].first == restrict_pair[i].first)
                for (unsigned int k = 0; k < restrict_pair.size(); ++k)
                    if (i != k && j != k && restrict_pair[k].second == restrict_pair[i].second && restrict_pair[j].second == restrict_pair[k].first) {
                        tag = false;
                        break;
                    }
            if (tag == false)
                break;
        }
        if (tag == false) {
            restrict_pair.erase(restrict_pair.begin() + i);
        } else
            ++i;
    }

    int max_prefix_num = size * (size - 1) / 2;
    memset(restrict_last, -1, size * sizeof(int));
    memset(restrict_next, -1, max_prefix_num * sizeof(int));
    total_restrict_num = 0;
    for (const auto &p : restrict_pair) {
        // p.first must be greater than p.second
        restrict_index[total_restrict_num] = p.first; // 第一个位置记录限制顶点
        restrict_next[total_restrict_num] = restrict_last[p.second];
        restrict_last[p.second] =
            total_restrict_num; // 互相记录，last记录的是顶点的约束位置，next记录的也是约束位置（从last开始链式记录）， index记录的是当前顶点的顶点
        ++total_restrict_num;
    }
}

int Schedule_IEP::get_multiplicity() const {
    std::vector<std::vector<int>> isomorphism_vec = get_isomorphism_vec();
    return isomorphism_vec.size();
}

void Schedule_IEP::aggressive_optimize(std::vector<std::pair<int, int>> &ordered_pairs) const {
    std::vector<std::vector<int>> isomorphism_vec = get_isomorphism_vec();

    std::vector<std::vector<std::vector<int>>> permutation_groups;
    permutation_groups.clear();
    for (const std::vector<int> &v : isomorphism_vec)
        permutation_groups.push_back(calc_permutation_group(v, size));

    ordered_pairs.clear();

    // delete permutation group which contains 1 permutation with 2 elements and some permutation with 1 elements,
    // and record the corresponding restriction.
    for (unsigned int i = 0; i < permutation_groups.size();) {
        int two_element_number = 0;
        std::pair<int, int> found_pair;
        for (const std::vector<int> &v : permutation_groups[i])
            if (v.size() == 2) {
                ++two_element_number;
                found_pair = std::pair<int, int>(v[0], v[1]);
            } else if (v.size() != 1) {
                two_element_number = -1;
                break;
            }
        if (two_element_number == 1) {
            permutation_groups.erase(permutation_groups.begin() + i);
            isomorphism_vec.erase(isomorphism_vec.begin() + i);
            ordered_pairs.push_back(found_pair);
            assert(found_pair.first < found_pair.second);
        } else
            ++i;
    }

    Pattern base_dag(size);
    for (const std::pair<int, int> &pair : ordered_pairs)
        base_dag.add_ordered_edge(pair.first, pair.second);

    bool changed = true;
    while (changed && isomorphism_vec.size() != 1) {
        // use restrictions to delete other isomophism
        for (unsigned int i = 0; i < isomorphism_vec.size();) {
            Pattern test_dag(base_dag);
            const std::vector<int> &iso = isomorphism_vec[i];
            for (const std::pair<int, int> &pair : ordered_pairs)
                test_dag.add_ordered_edge(iso[pair.first], iso[pair.second]);
            if (test_dag.is_dag() == false) // is not dag means conflict
            {
                permutation_groups.erase(permutation_groups.begin() + i);
                isomorphism_vec.erase(isomorphism_vec.begin() + i);
            } else
                ++i;
        }

        changed = false;
        std::pair<int, int> found_pair;
        for (unsigned int i = 0; i < permutation_groups.size();) {
            int two_element_number = 0;
            for (const std::vector<int> &v : permutation_groups[i])
                if (v.size() == 2) {
                    ++two_element_number;
                    found_pair = std::pair<int, int>(v[0], v[1]);
                    break;
                }
            if (two_element_number >= 1) {
                permutation_groups.erase(permutation_groups.begin() + i);
                isomorphism_vec.erase(isomorphism_vec.begin() + i);
                assert(found_pair.first < found_pair.second);
                ordered_pairs.push_back(found_pair);
                base_dag.add_ordered_edge(found_pair.first, found_pair.second);
                changed = true;
                break;
            } else
                ++i;
        }
    }
    assert(isomorphism_vec.size() == 1);
}

// Schedule_IEP::aggressive_optimize(...) can only get one valid restrictions
// but in this function, we try our best to find more restrictions
// WARNING: the restrictions in ordered_pairs_vector may NOT CORRECT

// 在前面不是已经获取了一个可行的调度方案集合？已经确定同构了？

// 处理所有与当前模式同构的调度方案，与前面的操作似乎不一定有交集
void Schedule_IEP::aggressive_optimize_get_all_pairs(std::vector<std::vector<std::pair<int, int>>> &ordered_pairs_vector) {

    // 获取所有与当前处理模式同构的调度方案
    std::vector<std::vector<int>> isomorphism_vec = get_isomorphism_vec();

    std::vector<std::vector<std::vector<int>>> permutation_groups;
    permutation_groups.clear();
    for (const std::vector<int> &v : isomorphism_vec)                  // 处理所有的同构的调度方案
        permutation_groups.push_back(calc_permutation_group(v, size)); //

    // 初始化
    ordered_pairs_vector.clear();
    std::vector<std::pair<int, int>> ordered_pairs;
    ordered_pairs.clear();

    // delete permutation group which contains 1 permutation with 2 elements and some permutation with 1 elements,
    // and record the corresponding restriction.

    // 删除所有可能产生同构的排列组合
    for (unsigned int i = 0; i < permutation_groups.size();) {
        int two_element_number = 0;
        std::pair<int, int> found_pair;
        for (const std::vector<int> &v : permutation_groups[i])
            if (v.size() == 2) {
                ++two_element_number;
                found_pair = std::pair<int, int>(v[0], v[1]);
            } else if (v.size() != 1) {
                two_element_number = -1;
                break;
            }
        if (two_element_number == 1) {
            permutation_groups.erase(permutation_groups.begin() + i);
            isomorphism_vec.erase(isomorphism_vec.begin() + i);
            ordered_pairs.push_back(found_pair);
            assert(found_pair.first < found_pair.second);
        } else
            ++i;
    }

    Pattern base_dag(size);
    for (const std::pair<int, int> &pair : ordered_pairs)
        base_dag.add_ordered_edge(pair.first, pair.second);

    aggressive_optimize_dfs(base_dag, isomorphism_vec, permutation_groups, ordered_pairs, ordered_pairs_vector);
}

// 这又是干嘛？
void Schedule_IEP::aggressive_optimize_dfs(Pattern base_dag, std::vector<std::vector<int>> isomorphism_vec,
                                           std::vector<std::vector<std::vector<int>>> permutation_groups,
                                           std::vector<std::pair<int, int>> ordered_pairs,
                                           std::vector<std::vector<std::pair<int, int>>> &ordered_pairs_vector) {

    for (unsigned int i = 0; i < isomorphism_vec.size();) {
        Pattern test_dag(base_dag);
        const std::vector<int> &iso = isomorphism_vec[i];
        for (const std::pair<int, int> &pair : ordered_pairs)
            test_dag.add_ordered_edge(iso[pair.first], iso[pair.second]);
        if (test_dag.is_dag() == false) // is not dag means conflict
        {
            permutation_groups.erase(permutation_groups.begin() + i);
            isomorphism_vec.erase(isomorphism_vec.begin() + i);
        } else
            ++i;
    }

    if (isomorphism_vec.size() == 1) {
        ordered_pairs_vector.push_back(ordered_pairs);
        return;
    }

    std::pair<int, int> found_pair;
    for (unsigned int i = 0; i < permutation_groups.size();) {
        int two_element_number = 0;
        for (const std::vector<int> &v : permutation_groups[i])
            if (v.size() == 2) {
                ++two_element_number;
                found_pair = std::pair<int, int>(v[0], v[1]);
                std::vector<std::vector<std::vector<int>>> next_permutation_groups = permutation_groups;
                std::vector<std::vector<int>> next_isomorphism_vec = isomorphism_vec;
                std::vector<std::pair<int, int>> next_ordered_pairs = ordered_pairs;
                Pattern next_base_dag = base_dag;

                next_permutation_groups.erase(next_permutation_groups.begin() + i);
                next_isomorphism_vec.erase(next_isomorphism_vec.begin() + i);
                assert(found_pair.first < found_pair.second);
                next_ordered_pairs.push_back(found_pair);
                next_base_dag.add_ordered_edge(found_pair.first, found_pair.second);

                aggressive_optimize_dfs(next_base_dag, next_isomorphism_vec, next_permutation_groups, next_ordered_pairs, ordered_pairs_vector);
            }
        if (two_element_number >= 1) {
            break;
        } else {
            ++i;
        }
    }
}

void Schedule_IEP::GraphZero_aggressive_optimize(std::vector<std::pair<int, int>> &ordered_pairs) const {
    std::vector<std::vector<int>> Aut;
    GraphZero_get_automorphisms(Aut);

    std::vector<std::pair<int, int>> L;
    L.clear();

    for (int v = 0; v < size; ++v) { // iterate all elements in schedule
        std::vector<std::vector<int>> stabilized_aut;
        stabilized_aut.clear();

        for (int i = 0; i < Aut.size(); ++i) {
            std::vector<int> &x = Aut[i];
            if (x[v] == v) {
                stabilized_aut.push_back(x);
            } else {
                int x1 = v, x2 = x[v];
                if (x1 > x2) {
                    int tmp = x1;
                    x1 = x2;
                    x2 = tmp;
                }
                bool tag = true;
                std::pair<int, int> cur = std::make_pair(x1, x2);
                for (int j = 0; j < L.size(); ++j)
                    if (L[j] == cur) {
                        tag = false;
                        break;
                    }
                if (tag) {
                    L.push_back(cur);
                }
            }
        }
        Aut = stabilized_aut;
    }

    ordered_pairs.clear(); // In GraphZero paper, this vector's name is 'L'

    for (int i = 0; i < L.size(); ++i) {
        bool tag = true;
        for (int j = 0; j < ordered_pairs.size(); ++j)
            if (L[i].second == ordered_pairs[j].second) {
                tag = false;
                if (L[i].first > ordered_pairs[j].first)
                    ordered_pairs[j].first = L[i].first;
                break;
            }
        if (tag)
            ordered_pairs.push_back(L[i]);
    }
}

void Schedule_IEP::GraphZero_get_automorphisms(std::vector<std::vector<int>> &Aut) const {
    int p[size];
    Aut.clear();
    for (int i = 0; i < size; ++i)
        p[i] = i;
    do {
        bool tag = true;
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j)
                if (adj_mat[INDEX(i, j, size)] != adj_mat[INDEX(p[i], p[j], size)]) {
                    tag = false;
                    break;
                }
            if (!tag)
                break;
        }
        if (tag) {
            std::vector<int> tmp;
            tmp.clear();
            for (int i = 0; i < size; ++i)
                tmp.push_back(p[i]);
            Aut.push_back(tmp);
        }
    } while (std::next_permutation(p, p + size));
}

// 获取所有与当前模式同构的调度方案，与添加约束打破对称性有关
std::vector<std::vector<int>> Schedule_IEP::get_isomorphism_vec() const {
    // 获取全部的排列
    unsigned int pow = 1;
    for (int i = 2; i <= size; ++i)
        pow *= i;
    std::vector<std::vector<int>> vec;
    vec.clear();
    bool use[size];
    for (int i = 0; i < size; ++i)
        use[i] = false;
    std::vector<int> tmp_vec;
    get_full_permutation(vec, use, tmp_vec, 0);
    assert(vec.size() == pow);

    // 获取所有的不同构的调度方案
    std::vector<std::vector<int>> isomorphism_vec;
    isomorphism_vec.clear();
    for (const std::vector<int> &v : vec) {
        bool flag = true;
        for (int i = 0; i < size; ++i)
            for (int j = i + 1; j < size; ++j)
                if (adj_mat[INDEX(i, j, size)] != 0)
                    if (adj_mat[INDEX(v[i], v[j], size)] == 0) // not isomorphism
                    {
                        flag = false;
                        break;
                    }
        if (flag == true)
            isomorphism_vec.push_back(v);
    }
    return isomorphism_vec;
}

void Schedule_IEP::get_full_permutation(std::vector<std::vector<int>> &vec, bool use[], std::vector<int> tmp_vec, int depth) const {
    if (depth == size) {
        vec.push_back(tmp_vec);
        return;
    }
    for (int i = 0; i < size; ++i)
        if (use[i] == false) {
            use[i] = true;
            tmp_vec.push_back(i);
            get_full_permutation(vec, use, tmp_vec, depth + 1);
            tmp_vec.pop_back();
            use[i] = false;
        }
}

// 计算随机排列的组合。输入一个调度方案以及大小，输出？
std::vector<std::vector<int>> Schedule_IEP::calc_permutation_group(const std::vector<int> vec, int size) {
    bool use[size];
    for (int i = 0; i < size; ++i)
        use[i] = false;
    std::vector<std::vector<int>> res;
    res.clear();

    //
    for (unsigned int i = 0; i < vec.size(); ++i)
        if (use[i] == false) { // 如果当前顶点没有被处理过
            std::vector<int> tmp_vec;
            tmp_vec.clear();
            tmp_vec.push_back(i); // 记录当前顶点的索引
            use[i] = true;
            int x = vec[i];           // 对于当前顶点进行处理
            while (use[x] == false) { // 如果当前顶点作为索引时，里面的顶点没有处理过
                use[x] = true;
                tmp_vec.push_back(x); // 记录当前顶点的索引
                x = vec[x];
            }
            res.push_back(tmp_vec);
        }
    return res;
}

void Schedule_IEP::performance_modeling(int *best_order, std::vector<std::vector<int>> &candidates, int v_cnt, e_index_t e_cnt) {
    int *order;
    int *rank;

    double *p_size;
    int max_degree = get_max_degree();
    p_size = new double[max_degree];

    double p = e_cnt * 1.0 / v_cnt / v_cnt;

    p_size[0] = v_cnt;
    for (int i = 1; i < max_degree; ++i) {
        p_size[i] = p_size[i - 1] * p;
    }

    order = new int[size];
    rank = new int[size];

    double min_val = 1e18;
    bool have_best = false;
    std::vector<int> invariant_size[size];
    for (const std::vector<int> &vec : candidates) {
        for (int i = 0; i < size; ++i)
            order[i] = vec[i];
        // check whether it is valid schedule
        bool is_valid = true;
        for (int i = 1; i < size; ++i) {
            bool have_edge = false;
            for (int j = 0; j < i; ++j)
                if (adj_mat[INDEX(order[i], order[j], size)]) {
                    have_edge = true;
                    break;
                }
            if (have_edge == false) {
                is_valid = false;
                break;
            }
        }
        if (is_valid == false)
            continue;

        for (int i = 0; i < size; ++i)
            rank[order[i]] = i;
        int *cur_adj_mat;
        cur_adj_mat = new int[size * size];
        for (int i = 0; i < size; ++i)
            for (int j = 0; j < size; ++j)
                cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

        std::vector<std::pair<int, int>> restricts;
        // TODO BUG!!!!!
        GraphZero_aggressive_optimize(restricts);
        int restricts_size = restricts.size();
        std::sort(restricts.begin(), restricts.end());
        double *sum;
        sum = new double[restricts_size];
        for (int i = 0; i < restricts_size; ++i)
            sum[i] = 0;
        int *tmp;
        tmp = new int[size];
        for (int i = 0; i < size; ++i)
            tmp[i] = i;
        do {
            for (int i = 0; i < restricts_size; ++i)
                if (tmp[restricts[i].first] > tmp[restricts[i].second]) {
                    sum[i] += 1;
                } else
                    break;
        } while (std::next_permutation(tmp, tmp + size));
        double total = 1;
        for (int i = 2; i <= size; ++i)
            total *= i;
        for (int i = 0; i < restricts_size; ++i)
            sum[i] = sum[i] / total;
        for (int i = restricts_size - 1; i > 0; --i)
            sum[i] /= sum[i - 1];

        double val = 1;
        for (int i = 0; i < size; ++i)
            invariant_size[i].clear();
        for (int i = size - 1; i >= 0; --i) {
            int cnt_forward = 0;
            int cnt_backward = 0;
            for (int j = 0; j < i; ++j)
                if (cur_adj_mat[INDEX(j, i, size)])
                    ++cnt_forward;
            for (int j = i + 1; j < size; ++j)
                if (cur_adj_mat[INDEX(j, i, size)])
                    ++cnt_backward;

            int c = cnt_forward;
            for (int j = i - 1; j >= 0; --j)
                if (cur_adj_mat[INDEX(j, i, size)])
                    invariant_size[j].push_back(c--);

            for (int j = 0; j < invariant_size[i].size(); ++j)
                if (invariant_size[i][j] > 1)
                    val += p_size[invariant_size[i][j] - 1] + p_size[1];
            for (int j = 0; j < restricts_size; ++j)
                if (restricts[j].second == i)
                    val *= sum[j];
            val *= p_size[cnt_forward];
        }
        if (have_best == false || val < min_val) {
            have_best = true;
            for (int i = 0; i < size; ++i)
                best_order[i] = order[i];
            min_val = val;
        }
        delete[] sum;
        delete[] tmp;
        delete[] cur_adj_mat;
    }

    delete[] order;
    delete[] rank;
    delete[] p_size;
}

void Schedule_IEP::bug_performance_modeling(int *best_order, std::vector<std::vector<int>> &candidates, int v_cnt, e_index_t e_cnt) {
    int *order;
    int *rank;

    double *p_size;
    int max_degree = get_max_degree();
    p_size = new double[max_degree];

    double p = e_cnt * 1.0 / v_cnt / v_cnt;

    p_size[0] = v_cnt;
    for (int i = 1; i < max_degree; ++i) {
        p_size[i] = p_size[i - 1] * p;
    }

    order = new int[size];
    rank = new int[size];

    double min_val;
    bool have_best = false;
    std::vector<int> invariant_size[size];
    for (const std::vector<int> &vec : candidates) {
        for (int i = 0; i < size; ++i)
            order[i] = vec[i];
        // check whether it is valid schedule
        bool is_valid = true;
        for (int i = 1; i < size; ++i) {
            bool have_edge = false;
            for (int j = 0; j < i; ++j)
                if (adj_mat[INDEX(order[i], order[j], size)]) {
                    have_edge = true;
                    break;
                }
            if (have_edge == false) {
                is_valid = false;
                break;
            }
        }
        if (is_valid == false)
            continue;

        for (int i = 0; i < size; ++i)
            rank[order[i]] = i;
        int *cur_adj_mat;
        cur_adj_mat = new int[size * size];
        for (int i = 0; i < size; ++i)
            for (int j = 0; j < size; ++j)
                cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

        std::vector<std::vector<std::pair<int, int>>> restricts_vector;
        restricts_generate(cur_adj_mat, restricts_vector);
        for (int restricts_rank = 0; restricts_rank < restricts_vector.size(); ++restricts_rank) {
            std::vector<std::pair<int, int>> &restricts = restricts_vector[restricts_rank];
            int restricts_size = restricts.size();
            std::sort(restricts.begin(), restricts.end());
            double *sum;
            sum = new double[restricts_size];
            for (int i = 0; i < restricts_size; ++i)
                sum[i] = 0;
            int *tmp;
            tmp = new int[size];
            for (int i = 0; i < size; ++i)
                tmp[i] = i;
            do {
                for (int i = 0; i < restricts_size; ++i)
                    if (tmp[restricts[i].first] > tmp[restricts[i].second]) {
                        sum[i] += 1;
                    } else
                        break;
            } while (std::next_permutation(tmp, tmp + size));
            double total = 1;
            for (int i = 2; i <= size; ++i)
                total *= i;
            for (int i = 0; i < restricts_size; ++i)
                sum[i] = sum[i] / total;
            for (int i = restricts_size - 1; i > 0; --i)
                sum[i] /= sum[i - 1];

            double val = 1;
            for (int i = 0; i < size; ++i)
                invariant_size[i].clear();
            for (int i = size - 1; i >= 0; --i) {
                int cnt_forward = 0;
                int cnt_backward = 0;
                for (int j = 0; j < i; ++j)
                    if (cur_adj_mat[INDEX(j, i, size)])
                        ++cnt_forward;
                for (int j = i + 1; j < size; ++j)
                    if (cur_adj_mat[INDEX(j, i, size)])
                        ++cnt_backward;

                int c = cnt_forward;
                for (int j = i - 1; j >= 0; --j)
                    if (cur_adj_mat[INDEX(j, i, size)])
                        invariant_size[j].push_back(c--);

                for (int j = 0; j < invariant_size[i].size(); ++j)
                    if (invariant_size[i][j] > 1)
                        val += p_size[invariant_size[i][j] - 1] + p_size[1];
                for (int j = 0; j < restricts_size; ++j)
                    if (restricts[j].second == i)
                        val *= sum[j];
                val *= p_size[cnt_forward];
            }
            if (have_best == false || val < min_val) {
                have_best = true;
                for (int i = 0; i < size; ++i)
                    best_order[i] = order[i];
                min_val = val;
            }
            delete[] sum;
            delete[] tmp;
        }
        delete[] cur_adj_mat;
    }

    delete[] order;
    delete[] rank;
    delete[] p_size;
}

void Schedule_IEP::new_performance_modeling(int *best_order, std::vector<std::vector<int>> &candidates, int v_cnt, e_index_t e_cnt,
                                            long long tri_cnt) {
    int *order;
    int *rank;

    double *p_size;
    double *pp_size;
    int max_degree = get_max_degree();
    p_size = new double[max_degree];
    pp_size = new double[max_degree];

    double p0 = e_cnt * 1.0 / v_cnt / v_cnt;
    double p1 = tri_cnt * 1.0 * v_cnt / e_cnt / e_cnt;

    p_size[0] = v_cnt;
    for (int i = 1; i < max_degree; ++i) {
        p_size[i] = p_size[i - 1] * p0;
    }
    pp_size[0] = 1;
    for (int i = 1; i < max_degree; ++i) {
        pp_size[i] = pp_size[i - 1] * p1;
    }

    order = new int[size];
    rank = new int[size];

    double min_val;
    bool have_best = false;
    std::vector<int> invariant_size[size];
    for (const std::vector<int> &vec : candidates) {
        for (int i = 0; i < size; ++i)
            order[i] = vec[i];
        // check whether it is valid schedule
        bool is_valid = true;
        for (int i = 1; i < size; ++i) {
            bool have_edge = false;
            for (int j = 0; j < i; ++j)
                if (adj_mat[INDEX(order[i], order[j], size)]) {
                    have_edge = true;
                    break;
                }
            if (have_edge == false) {
                is_valid = false;
                break;
            }
        }
        if (is_valid == false)
            continue;

        for (int i = 0; i < size; ++i)
            rank[order[i]] = i;
        int *cur_adj_mat;
        cur_adj_mat = new int[size * size];
        for (int i = 0; i < size; ++i)
            for (int j = 0; j < size; ++j)
                cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

        std::vector<std::vector<std::pair<int, int>>> restricts_vector;
        restricts_generate(cur_adj_mat, restricts_vector);
        for (int restricts_rank = 0; restricts_rank < restricts_vector.size(); ++restricts_rank) {
            std::vector<std::pair<int, int>> &restricts = restricts_vector[restricts_rank];
            int restricts_size = restricts.size();
            std::sort(restricts.begin(), restricts.end());
            double *sum;
            sum = new double[restricts_size];
            for (int i = 0; i < restricts_size; ++i)
                sum[i] = 0;
            int *tmp;
            tmp = new int[size];
            for (int i = 0; i < size; ++i)
                tmp[i] = i;
            do {
                for (int i = 0; i < restricts_size; ++i)
                    if (tmp[restricts[i].first] > tmp[restricts[i].second]) {
                        sum[i] += 1;
                    } else
                        break;
            } while (std::next_permutation(tmp, tmp + size));
            double total = 1;
            for (int i = 2; i <= size; ++i)
                total *= i;
            for (int i = 0; i < restricts_size; ++i)
                sum[i] = sum[i] / total;
            for (int i = restricts_size - 1; i > 0; --i)
                sum[i] /= sum[i - 1];

            double val = 1;
            for (int i = 0; i < size; ++i)
                invariant_size[i].clear();
            for (int i = size - 1; i >= 0; --i) {
                int cnt_forward = 0;
                int cnt_backward = 0;
                for (int j = 0; j < i; ++j)
                    if (cur_adj_mat[INDEX(j, i, size)])
                        ++cnt_forward;
                for (int j = i + 1; j < size; ++j)
                    if (cur_adj_mat[INDEX(j, i, size)])
                        ++cnt_backward;

                int c = cnt_forward;
                for (int j = i - 1; j >= 0; --j)
                    if (cur_adj_mat[INDEX(j, i, size)])
                        invariant_size[j].push_back(c--);

                for (int j = 0; j < invariant_size[i].size(); ++j)
                    if (invariant_size[i][j] > 1)
                        val += p_size[1] * pp_size[invariant_size[i][j] - 2] + p_size[1];
                val += 1;
                for (int j = 0; j < restricts_size; ++j)
                    if (restricts[j].second == i)
                        val *= sum[j];
                val *= p_size[1] * pp_size[cnt_forward - 1];
            }
            if (have_best == false || val < min_val) {
                have_best = true;
                for (int i = 0; i < size; ++i)
                    best_order[i] = order[i];
                min_val = val;
            }
            delete[] sum;
            delete[] tmp;
        }
        delete[] cur_adj_mat;
    }

    delete[] order;
    delete[] rank;
    delete[] p_size;
    delete[] pp_size;
}

void Schedule_IEP::init_in_exclusion_optimize() {
    int optimize_num = in_exclusion_optimize_num;

    assert(in_exclusion_optimize_num > 1);

    int *id;
    id = new int[optimize_num];

    int *in_exclusion_val;
    in_exclusion_val = new int[optimize_num * 2];

    for (int n = 1; n <= optimize_num; ++n) {
        DisjointSetUnion dsu(n);
        int m = n * (n - 1) / 2;

        in_exclusion_val[2 * n - 2] = 0;
        in_exclusion_val[2 * n - 1] = 0;

        if (n == 1) {
            ++in_exclusion_val[0];
            continue;
        }

        std::pair<int, int> edge[m];
        e_index_t e_cnt = 0;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < i; ++j)
                edge[e_cnt++] = std::make_pair(i, j);

        for (int s = 0; s < (1 << m); ++s) {
            dsu.init();
            int bit_cnt = 0;
            for (int i = 0; i < m; ++i)
                if (s & (1 << i)) {
                    ++bit_cnt;
                    dsu.merge(edge[i].first, edge[i].second);
                }
            if (dsu.get_set_size() == 1) {
                if (bit_cnt & 1)
                    ++in_exclusion_val[2 * n - 1];
                else
                    ++in_exclusion_val[2 * n - 2];
            }
        }
    }

    in_exclusion_optimize_group.clear();
    in_exclusion_optimize_val.clear();

    get_in_exclusion_optimize_group(0, id, 0, in_exclusion_val);

    delete[] id;
    delete[] in_exclusion_val;
}

void Schedule_IEP::get_in_exclusion_optimize_group(int depth, int *id, int id_cnt, int *in_exclusion_val) {
    if (depth == in_exclusion_optimize_num) {
        int *size = new int[id_cnt];
        for (int i = 0; i < id_cnt; ++i)
            size[i] = 0;
        for (int i = 0; i < in_exclusion_optimize_num; ++i)
            size[id[i]]++;
        int val[2];
        val[0] = in_exclusion_val[size[0] * 2 - 2];
        val[1] = in_exclusion_val[size[0] * 2 - 1];
        for (int i = 1; i < id_cnt; ++i) {
            int tmp0 = val[0];
            int tmp1 = val[1];

            val[0] = tmp0 * in_exclusion_val[size[i] * 2 - 2] + tmp1 * in_exclusion_val[size[i] * 2 - 1];
            val[1] = tmp0 * in_exclusion_val[size[i] * 2 - 1] + tmp1 * in_exclusion_val[size[i] * 2 - 2];
        }

        std::vector<std::vector<int>> group;
        group.clear();
        for (int i = 0; i < id_cnt; ++i) {
            std::vector<int> cur;
            cur.clear();
            for (int j = 0; j < in_exclusion_optimize_num; ++j)
                if (id[j] == i)
                    cur.push_back(j);
            group.push_back(cur);
        }

        in_exclusion_optimize_group.push_back(group);
        in_exclusion_optimize_val.push_back(val[0] - val[1]);

        delete[] size;
        return;
    }

    id[depth] = id_cnt;

    get_in_exclusion_optimize_group(depth + 1, id, id_cnt + 1, in_exclusion_val);

    for (int i = 0; i < id_cnt; ++i) {
        id[depth] = i;
        get_in_exclusion_optimize_group(depth + 1, id, id_cnt, in_exclusion_val);
    }
}

void Schedule_IEP::print_schedule() const {
    printf("Schedule_IEP:\n");
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j)
            printf("%d", adj_mat[INDEX(i, j, size)]);
        puts("");
    }
    printf("loop_set_prefix_ids:");
    for (int i = 0; i < size; ++i)
        printf(" %d", loop_set_prefix_id[i]);
    printf("\nlast:");
    for (int i = 0; i < size; ++i)
        printf(" %d", last[i]);
    printf("\nnext:");
    for (int i = 0; i < size * (size - 1) / 2; ++i)
        printf(" %d", next[i]);
    printf("\nfather_prefix_id:");
    for (int i = 0; i < size; ++i)
        printf(" %d", father_prefix_id[i]);
    printf("\nrestrictions:");
    for (auto &pair : restrict_pair)
        printf(" (%d, %d)", pair.first, pair.second);
    printf("\nbreak_size:");
    for (int i = 0; i < size * (size - 1) / 2; i++)
        printf(" %d", break_size[i]);
    printf("\n");
}

void Schedule_IEP::GraphZero_performance_modeling(int *best_order, int v_cnt, e_index_t e_cnt) {
    int *order;
    int *rank;

    double *p_size;
    double *anti_p;
    p_size = new double[size];
    anti_p = new double[size];

    double p = e_cnt * 2.0 / v_cnt / v_cnt;

    p_size[0] = v_cnt;
    for (int i = 1; i < size; ++i) {
        p_size[i] = p_size[i - 1] * p;
        printf("p %d %.6lf\n", i, p_size[i]);
    }
    anti_p[0] = 1;
    for (int i = 1; i < size; ++i) {
        anti_p[i] = anti_p[i - 1] * (1 - p);
        printf("anti %d %.6lf\n", i, anti_p[i]);
    }

    order = new int[size];
    rank = new int[size];

    for (int i = 0; i < size; ++i)
        order[i] = i;
    double min_val = 1e18;
    bool have_best = false;
    do {
        // check whether it is valid schedule
        bool is_valid = true;
        for (int i = 1; i < size; ++i) {
            bool have_edge = false;
            for (int j = 0; j < i; ++j)
                if (adj_mat[INDEX(order[i], order[j], size)]) {
                    have_edge = true;
                    break;
                }
            if (have_edge == false) {
                is_valid = false;
                break;
            }
        }
        if (is_valid == false)
            continue;

        for (int i = 0; i < size; ++i)
            rank[order[i]] = i;
        int *cur_adj_mat;
        cur_adj_mat = new int[size * size];
        for (int i = 0; i < size; ++i)
            for (int j = 0; j < size; ++j)
                cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

        std::vector<std::pair<int, int>> restricts;
        GraphZero_aggressive_optimize(restricts);
        int restricts_size = restricts.size();
        std::sort(restricts.begin(), restricts.end());
        double *sum;
        sum = new double[restricts_size];
        for (int i = 0; i < restricts_size; ++i)
            sum[i] = 0;
        int *tmp;
        tmp = new int[size];
        for (int i = 0; i < size; ++i)
            tmp[i] = i;
        do {
            for (int i = 0; i < restricts_size; ++i)
                if (tmp[restricts[i].first] > tmp[restricts[i].second]) {
                    sum[i] += 1;
                } else
                    break;
        } while (std::next_permutation(tmp, tmp + size));
        double total = 1;
        for (int i = 2; i <= size; ++i)
            total *= i;
        for (int i = 0; i < restricts_size; ++i)
            sum[i] = sum[i] / total;
        for (int i = restricts_size - 1; i > 0; --i)
            sum[i] /= sum[i - 1];

        double val = 1e18;
        for (int i = size - 1; i >= 0; --i) {
            int cnt_forward = 0;
            for (int j = 0; j < i; ++j)
                if (cur_adj_mat[INDEX(j, i, size)])
                    ++cnt_forward;

            for (int j = 0; j < restricts_size; ++j)
                if (restricts[j].second == i)
                    val *= sum[j];
            //      val *= p_size[cnt_forward + 1] * anti_p[i - cnt_forward];
            val *= v_cnt;
            for (int j = 0; j < i - cnt_forward; ++j)
                val *= (1 - p);
            for (int j = 0; j < cnt_forward; ++j)
                val *= p;
        }
        if (have_best == false || val <= min_val) {
            have_best = true;
            for (int i = 0; i < size; ++i)
                best_order[i] = order[i];
            min_val = val;
            printf("gz upd %.10lf\n", val);
        }

        delete[] cur_adj_mat;
        delete[] sum;
        delete[] tmp;

    } while (std::next_permutation(order, order + size));

    delete[] order;
    delete[] rank;
    delete[] p_size;
    delete[] anti_p;
}

void Schedule_IEP::restrict_selection(int v_cnt, e_index_t e_cnt, long long tri_cnt,
                                      std::vector<std::vector<std::pair<int, int>>> ordered_pairs_vector,
                                      std::vector<std::pair<int, int>> &best_restricts) const {

    double *p_size;
    double *pp_size;
    int max_degree = get_max_degree();
    p_size = new double[max_degree];
    pp_size = new double[max_degree];

    double p0 = e_cnt * 1.0 / v_cnt / v_cnt;
    double p1 = tri_cnt * 1.0 * v_cnt / e_cnt / e_cnt;

    p_size[0] = v_cnt;
    for (int i = 1; i < max_degree; ++i) {
        p_size[i] = p_size[i - 1] * p0;
    }

    pp_size[0] = 1;
    for (int i = 1; i < max_degree; ++i) {
        pp_size[i] = pp_size[i - 1] * p1;
    }

    double min_val;
    bool have_best = false;
    std::vector<int> invariant_size[size];

    for (int cur_restricts = 0; cur_restricts < ordered_pairs_vector.size(); ++cur_restricts) {
        std::vector<std::pair<int, int>> &restricts = ordered_pairs_vector[cur_restricts];
        int restricts_size = restricts.size();
        std::sort(restricts.begin(), restricts.end());

        double *sum;
        sum = new double[restricts_size];
        for (int i = 0; i < restricts_size; ++i)
            sum[i] = 0;

        int *tmp;
        tmp = new int[size];
        for (int i = 0; i < size; ++i)
            tmp[i] = i;

        do {
            for (int i = 0; i < restricts_size; ++i)
                if (tmp[restricts[i].first] > tmp[restricts[i].second]) {
                    sum[i] += 1;
                } else
                    break;
        } while (std::next_permutation(tmp, tmp + size));

        double total = 1;
        for (int i = 2; i <= size; ++i)
            total *= i;

        for (int i = 0; i < restricts_size; ++i)
            sum[i] = sum[i] / total;
        for (int i = restricts_size - 1; i > 0; --i)
            sum[i] /= sum[i - 1];

        double val = 1;
        for (int i = 0; i < size; ++i)
            invariant_size[i].clear();

        for (int i = size - 1; i >= 0; --i) {
            int cnt_forward = 0;
            int cnt_backward = 0;
            for (int j = 0; j < i; ++j)
                if (adj_mat[INDEX(j, i, size)])
                    ++cnt_forward;
            for (int j = i + 1; j < size; ++j)
                if (adj_mat[INDEX(j, i, size)])
                    ++cnt_backward;

            int c = cnt_forward;
            for (int j = i - 1; j >= 0; --j)
                if (adj_mat[INDEX(j, i, size)])
                    invariant_size[j].push_back(c--);

            for (int j = 0; j < invariant_size[i].size(); ++j)
                if (invariant_size[i][j] > 1)
                    val += p_size[1] * pp_size[invariant_size[i][j] - 2] + p_size[1];
            val += 1;
            for (int j = 0; j < restricts_size; ++j)
                if (restricts[j].second == i)
                    val *= sum[j];
            val *= p_size[1] * pp_size[cnt_forward - 1];
        }
        if (have_best == false || val < min_val) {
            have_best = true;
            best_restricts = restricts;
            min_val = val;
        }

        delete[] sum;
        delete[] tmp;
    }

    delete[] p_size;
    delete[] pp_size;
    assert(have_best);
}

// 根据模式生成约束  不同的处理方式应该有不同的约束
void Schedule_IEP::restricts_generate(const int *cur_adj_mat, std::vector<std::vector<std::pair<int, int>>> &restricts) {
    Schedule_IEP schedule(cur_adj_mat, get_size());
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

// 计算调度方案中底层不相交顶点的数量，是调度方案是否优越的一个衡量。
int Schedule_IEP::get_vec_optimize_num(const std::vector<int> &vec) {
    // 首先判断是否后续处理的顶点与前面处理的顶点至少存在一条边，前面已经实现过这个功能？
    bool is_valid = true;
    for (int i = 1; i < size; ++i) {
        bool have_edge = false;
        for (int j = 0; j < i; ++j)
            if (adj_mat[INDEX(vec[i], vec[j], size)]) {
                have_edge = true;
                break;
            }
        if (have_edge == false) {
            is_valid = false;
            break;
        }
    }
    if (!is_valid)
        return -1;

    // 根据能够产生多少个不相关顶点进行计数？
    for (int k = 2; k <= size - 2; ++k) { // 之所以k最大可能为size - 2，是因为第一次取了一条边而不是一个点，所以前两个点都不可能在IEP内
        bool flag = true;
        for (int i = size - k + 1; i < size; ++i) // 从当前顶点处理到最后一个顶点，查看当前顶点是否与前一个顶点之间存在边
            if (adj_mat[INDEX(vec[size - k], vec[i], size)]) {
                flag = false;
                break;
            }
        if (flag == false)
            return k - 1;
    }
    return size - 2;
}

double Schedule_IEP::our_estimate_schedule_restrict(const std::vector<int> &order, const std::vector<std::pair<int, int>> &pairs, int v_cnt,
                                                    e_index_t e_cnt, long long tri_cnt) {
    int max_degree = get_max_degree();

    double p_size[max_degree];
    double pp_size[max_degree];

    double p0 = e_cnt * 1.0 / v_cnt / v_cnt;
    double p1 = tri_cnt * 1.0 * v_cnt / e_cnt / e_cnt;

    p_size[0] = v_cnt;
    for (int i = 1; i < max_degree; ++i) {
        p_size[i] = p_size[i - 1] * p0;
    }
    pp_size[0] = 1;
    for (int i = 1; i < max_degree; ++i) {
        pp_size[i] = pp_size[i - 1] * p1;
    }

    int rank[size];
    for (int i = 0; i < size; ++i)
        rank[order[i]] = i;

    int *cur_adj_mat;
    cur_adj_mat = new int[size * size];
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

    std::vector<std::pair<int, int>> restricts = pairs;
    int restricts_size = restricts.size();
    std::sort(restricts.begin(), restricts.end());

    double sum[restricts_size];
    for (int i = 0; i < restricts_size; ++i)
        sum[i] = 0;

    int tmp[size];
    for (int i = 0; i < size; ++i)
        tmp[i] = i;
    do {
        for (int i = 0; i < restricts_size; ++i)
            if (tmp[restricts[i].first] > tmp[restricts[i].second]) {
                sum[i] += 1;
            } else
                break;
    } while (std::next_permutation(tmp, tmp + size));

    double total = 1;
    for (int i = 2; i <= size; ++i)
        total *= i;
    for (int i = 0; i < restricts_size; ++i)
        sum[i] = sum[i] / total;
    for (int i = restricts_size - 1; i > 0; --i)
        sum[i] /= sum[i - 1];

    std::vector<int> invariant_size[size];
    for (int i = 0; i < size; ++i)
        invariant_size[i].clear();

    double val = 1;
    for (int i = size - 1; i >= 0; --i) {
        int cnt_forward = 0;
        int cnt_backward = 0;
        for (int j = 0; j < i; ++j)
            if (cur_adj_mat[INDEX(j, i, size)])
                ++cnt_forward;
        for (int j = i + 1; j < size; ++j)
            if (cur_adj_mat[INDEX(j, i, size)])
                ++cnt_backward;

        int c = cnt_forward;
        for (int j = i - 1; j >= 0; --j)
            if (cur_adj_mat[INDEX(j, i, size)])
                invariant_size[j].push_back(c--);

        for (int j = 0; j < invariant_size[i].size(); ++j)
            if (invariant_size[i][j] > 1)
                val += p_size[1] * pp_size[invariant_size[i][j] - 2] + p_size[1];
        val += 1;
        for (int j = 0; j < restricts_size; ++j)
            if (restricts[j].second == i)
                val *= sum[j];
        if (i) {
            val *= p_size[1] * pp_size[cnt_forward - 1];
        } else {
            val *= p_size[0];
        }
    }
    delete[] cur_adj_mat;

    return val;
}

double Schedule_IEP::GraphZero_estimate_schedule_restrict(const std::vector<int> &order, const std::vector<std::pair<int, int>> &pairs, int v_cnt,
                                                          e_index_t e_cnt) {
    int max_degree = get_max_degree();

    double p_size[max_degree];
    double p = e_cnt * 1.0 / v_cnt / v_cnt;

    p_size[0] = v_cnt;
    for (int i = 1; i < max_degree; ++i) {
        p_size[i] = p_size[i - 1] * p;
    }

    int rank[size];
    for (int i = 0; i < size; ++i)
        rank[order[i]] = i;

    int *cur_adj_mat;
    cur_adj_mat = new int[size * size];
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

    std::vector<std::pair<int, int>> restricts = pairs;
    int restricts_size = restricts.size();
    std::sort(restricts.begin(), restricts.end());

    double sum[restricts_size];
    for (int i = 0; i < restricts_size; ++i)
        sum[i] = 0;

    int tmp[size];
    for (int i = 0; i < size; ++i)
        tmp[i] = i;
    do {
        for (int i = 0; i < restricts_size; ++i)
            if (tmp[restricts[i].first] > tmp[restricts[i].second]) {
                sum[i] += 1;
            } else
                break;
    } while (std::next_permutation(tmp, tmp + size));

    double total = 1;
    for (int i = 2; i <= size; ++i)
        total *= i;

    for (int i = 0; i < restricts_size; ++i)
        sum[i] = sum[i] / total;
    for (int i = restricts_size - 1; i > 0; --i)
        sum[i] /= sum[i - 1];

    std::vector<int> invariant_size[size];
    for (int i = 0; i < size; ++i)
        invariant_size[i].clear();

    double val = 1;
    for (int i = size - 1; i >= 0; --i) {
        int cnt_forward = 0;
        int cnt_backward = 0;
        for (int j = 0; j < i; ++j)
            if (cur_adj_mat[INDEX(j, i, size)])
                ++cnt_forward;
        for (int j = i + 1; j < size; ++j)
            if (cur_adj_mat[INDEX(j, i, size)])
                ++cnt_backward;

        int c = cnt_forward;
        for (int j = i - 1; j >= 0; --j)
            if (cur_adj_mat[INDEX(j, i, size)])
                invariant_size[j].push_back(c--);

        for (int j = 0; j < invariant_size[i].size(); ++j)
            if (invariant_size[i][j] > 1)
                val += p_size[invariant_size[i][j] - 1] + p_size[1];
        for (int j = 0; j < restricts_size; ++j)
            if (restricts[j].second == i)
                val *= sum[j];
        val *= p_size[cnt_forward];
    }

    delete[] cur_adj_mat;

    return val;
}

double Schedule_IEP::Naive_estimate_schedule_restrict(const std::vector<int> &order, const std::vector<std::pair<int, int>> &pairs, int v_cnt,
                                                      e_index_t e_cnt) {

    double p = e_cnt * 2.0 / v_cnt / v_cnt;

    int rank[size];
    for (int i = 0; i < size; ++i)
        rank[order[i]] = i;

    int *cur_adj_mat;
    cur_adj_mat = new int[size * size];
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

    std::vector<std::pair<int, int>> restricts = pairs;
    int restricts_size = restricts.size();
    std::sort(restricts.begin(), restricts.end());

    double sum[restricts_size];
    for (int i = 0; i < restricts_size; ++i)
        sum[i] = 0;
    int tmp[size];

    for (int i = 0; i < size; ++i)
        tmp[i] = i;
    do {
        for (int i = 0; i < restricts_size; ++i)
            if (tmp[restricts[i].first] > tmp[restricts[i].second]) {
                sum[i] += 1;
            } else
                break;
    } while (std::next_permutation(tmp, tmp + size));

    double total = 1;
    for (int i = 2; i <= size; ++i)
        total *= i;

    for (int i = 0; i < restricts_size; ++i)
        sum[i] = sum[i] / total;
    for (int i = restricts_size - 1; i > 0; --i)
        sum[i] /= sum[i - 1];

    double val = 1;
    for (int i = size - 1; i >= 0; --i) {
        int cnt_forward = 0;
        for (int j = 0; j < i; ++j)
            if (cur_adj_mat[INDEX(j, i, size)])
                ++cnt_forward;

        for (int j = 0; j < restricts_size; ++j)
            if (restricts[j].second == i)
                val *= sum[j];
        val *= v_cnt;
        for (int j = 0; j < i - cnt_forward; ++j)
            val *= (1 - p);
        for (int j = 0; j < cnt_forward; ++j)
            val *= p;
    }

    delete[] cur_adj_mat;

    return val;
}

// 移除不合法的顶点处理顺序，模式中没有边存在的会被移除
void Schedule_IEP::remove_invalid_permutation(std::vector<std::vector<int>> &candidate_permutations) {
    for (unsigned int i = 0; i < candidate_permutations.size();) {
        const auto &vec = candidate_permutations[i]; // 处理每一个调度方案
        bool tag = true;
        for (int x = 1; x < size; ++x) {
            bool have_edge = false;
            for (int y = 0; y < x; ++y)
                if (adj_mat[INDEX(vec[x], vec[y], size)]) {
                    have_edge = true;
                    break;
                }
            if (!have_edge) {
                tag = false;
                break;
            }
        }
        if (tag) {
            ++i;
        } else {
            candidate_permutations.erase(candidate_permutations.begin() + i);
        }
    }
}

int Schedule_IEP::get_in_exclusion_optimize_num_when_not_optimize() {
    std::vector<int> I;
    for (int i = 0; i < size; ++i)
        I.push_back(i);
    return get_vec_optimize_num(I);
}

void Schedule_IEP::set_in_exclusion_optimize_redundancy() {
    int tmp = get_in_exclusion_optimize_num();
    if (tmp <= 1) {
        in_exclusion_optimize_redundancy = 1;
    } else {
        Graph *complete;
        DataLoader *D = new DataLoader();
        assert(D->load_complete(complete, get_size()));
        delete D;
        in_exclusion_optimize_redundancy = 1;
        long long ans = complete->pattern_matching(*this, 1);
        set_in_exclusion_optimize_num(0);
        long long true_ans = complete->pattern_matching(*this, 1);
        set_in_exclusion_optimize_num(tmp);
        delete complete;
        in_exclusion_optimize_redundancy = ans / true_ans;
    }
}
