#pragma once

#include "disjoint_set_union.h"
#include "pattern.h"
#include "prefix.h"
#include "types.h"

#include <cstdint>
#include <vector>

class Schedule_IEP {
  public:
    // TODO : more kinds of constructors to construct different Schedules from one Pattern
    Schedule_IEP(const Pattern &pattern, bool &is_pattern_valid, int performance_modeling_type, int restricts_type, bool use_in_exclusion_optimize,
                 int v_cnt, int64_t e_cnt, long long tri_cnt = 0, bool vertex_induced = false);
    // performance_modeling type = 0 : not use modeling
    //                      type = 1 : use our modeling
    //                      type = 2 : use GraphZero's modeling
    //                      type = 3 : use naive modeling
    // restricts_type = 0 : not use restricts
    //                = 1 : use our restricts
    //                = 2 : use GraphZero's restricts
    Schedule_IEP(const int *_adj_mat, int _size);
    Schedule_IEP(const Schedule_IEP &s) = delete;
    Schedule_IEP &operator=(const Schedule_IEP &s) = delete;
    ~Schedule_IEP();
    inline int get_total_prefix_num() const { return total_prefix_num; }
    inline int get_basic_prefix_num() const { return basic_prefix_num; }
    inline int get_father_prefix_id(int prefix_id) const { return father_prefix_id[prefix_id]; }
    inline int *get_father_prefix_id_ptr() const { return father_prefix_id; }
    inline int get_loop_set_prefix_id(int loop) const { return loop_set_prefix_id[loop]; }
    inline int *get_loop_set_prefix_id_ptr() const { return loop_set_prefix_id; }
    inline bool get_prefix_only_need_size(int prefix_id) const { return prefix[prefix_id].get_only_need_size(); }
    inline int get_size() const { return size; }
    inline int get_last(int i) const { return last[i]; }
    inline int *get_last_ptr() const { return last; }
    inline int get_next(int i) const { return next[i]; }
    inline int *get_next_ptr() const { return next; }
    inline int *get_break_size_ptr() const { return break_size; }
    inline int get_prefix_target(int i) const { return prefix_target[i]; }
    inline int *get_prefix_target_ptr() const { return prefix_target; }
    inline int get_in_exclusion_optimize_num() const { return in_exclusion_optimize_num; }
    int get_in_exclusion_optimize_num_when_not_optimize();
    void add_restrict(const std::vector<std::pair<int, int>> &restricts);
    inline int get_total_restrict_num() const { return total_restrict_num; }
    inline int get_restrict_last(int i) const { return restrict_last[i]; }
    inline int *get_restrict_last_ptr() const { return restrict_last; }
    inline int get_restrict_next(int i) const { return restrict_next[i]; }
    inline int *get_restrict_next_ptr() const { return restrict_next; }
    inline int get_restrict_index(int i) const { return restrict_index[i]; }
    inline int *get_restrict_index_ptr() const { return restrict_index; }
    inline int get_k_val() const { return k_val; } // see below (the k_val's definition line) before using this function
    int get_max_degree() const;
    int get_multiplicity() const;
    void aggressive_optimize(std::vector<std::pair<int, int>> &ordered_pairs) const;
    void aggressive_optimize_get_all_pairs(std::vector<std::vector<std::pair<int, int>>> &ordered_pairs_vector);
    void aggressive_optimize_dfs(Pattern base_dag, std::vector<std::vector<int>> isomorphism_vec,
                                 std::vector<std::vector<std::vector<int>>> permutation_groups, std::vector<std::pair<int, int>> ordered_pairs,
                                 std::vector<std::vector<std::pair<int, int>>> &ordered_pairs_vector);
    void restrict_selection(int v_cnt, e_index_t e_cnt, long long tri_cnt, std::vector<std::vector<std::pair<int, int>>> ordered_pairs_vector,
                            std::vector<std::pair<int, int>> &best_restricts) const;
    void restricts_generate(const int *cur_adj_mat, std::vector<std::vector<std::pair<int, int>>> &restricts);

    void GraphZero_aggressive_optimize(std::vector<std::pair<int, int>> &ordered_pairs) const;
    void GraphZero_get_automorphisms(std::vector<std::vector<int>> &Aut) const;

    std::vector<std::vector<int>> get_isomorphism_vec() const;
    static std::vector<std::vector<int>> calc_permutation_group(const std::vector<int> vec, int size);
    inline const int *get_adj_mat_ptr() const { return adj_mat; }

    inline void set_in_exclusion_optimize_redundancy(long long redundancy) { in_exclusion_optimize_redundancy = redundancy; }
    inline long long get_in_exclusion_optimize_redundancy() const { return in_exclusion_optimize_redundancy; }

    void print_schedule() const;

    void update_loop_invariant_for_fsm();

    std::vector<std::pair<int, int>> restrict_pair;

    std::vector<int> in_exclusion_optimize_vertex_id;
    std::vector<bool> in_exclusion_optimize_vertex_flag;
    std::vector<int> in_exclusion_optimize_vertex_coef;

    std::vector<int> in_exclusion_optimize_coef;
    std::vector<bool> in_exclusion_optimize_flag;
    std::vector<int> in_exclusion_optimize_ans_pos;

    int *break_size = nullptr;

    bool is_vertex_induced;

  private:
    int *adj_mat = nullptr;
    int *father_prefix_id = nullptr;
    int *last = nullptr;
    int *next = nullptr;
    int *loop_set_prefix_id = nullptr;
    int *prefix_target =
        nullptr; // 这个是给带label时使用的，因为带label时，需要提前知道一个prefix最终是为了给哪个点作为循环集合来确定prefix中点的label，比如这个prefix经过几次求交后，得到的集合要给pattern中的第4个点作为循环集合，那么target就是4
    Prefix *prefix = nullptr;
    int *restrict_last = nullptr;
    int *restrict_next = nullptr;
    int *restrict_index = nullptr;
    int size;
    int total_prefix_num;
    int basic_prefix_num;
    int total_restrict_num;
    int in_exclusion_optimize_num;
    int k_val; // inner k loop, WARNING: this val not always meaningful @TODO here
               // only when performance_modeling_type == 1 , this val will be calculated.
    long long in_exclusion_optimize_redundancy;

    std::vector<std::vector<std::vector<int>>> in_exclusion_optimize_group;
    std::vector<int> in_exclusion_optimize_val;

    bool check_connectivity() const;
    void setup_optimization_info(bool use_in_exclusion_optimize = true);
    void copy_adj_mat_from(const std::vector<int> &vec, const int *src_adj_mat);
    void build_loop_invariant(int in_exclusion_optimize_num = 0);
    int find_father_prefix(int data_size, const int *data);
    void get_full_permutation(std::vector<std::vector<int>> &vec, bool use[], std::vector<int> tmp_vec, int depth) const;
    void performance_modeling(int *best_order, std::vector<std::vector<int>> &candidates, int v_cnt, e_index_t e_cnt);
    void bug_performance_modeling(int *best_order, std::vector<std::vector<int>> &candidates, int v_cnt, e_index_t e_cnt);
    void new_performance_modeling(int *best_order, std::vector<std::vector<int>> &candidates, int v_cnt, e_index_t e_cnt, long long tri_cnt);
    void GraphZero_performance_modeling(int *best_order, int v_cnt, e_index_t e_cnt);

    double new_estimate_schedule_restrict(const std::vector<std::pair<int, int>> &restrictions, int v_cnt, e_index_t e_cnt, long long tri_cnt);
    double our_estimate_schedule_restrict(const std::vector<int> &order, const std::vector<std::pair<int, int>> &pairs, int v_cnt, e_index_t e_cnt,
                                          long long tri_cnt);
    double GraphZero_estimate_schedule_restrict(const std::vector<int> &order, const std::vector<std::pair<int, int>> &pairs, int v_cnt,
                                                e_index_t e_cnt);
    double Naive_estimate_schedule_restrict(const std::vector<int> &order, const std::vector<std::pair<int, int>> &paris, int v_cnt, e_index_t e_cnt);

    void get_in_exclusion_optimize_group(int depth, int *id, int id_cnt, int *in_exclusion_val);
    // use principle of inclusion-exclusion to optimize
    void init_in_exclusion_optimize();

    int get_vec_optimize_num(const std::vector<int> &vec);

    void remove_invalid_permutation(std::vector<std::vector<int>> &candidate_permutations);

    inline void set_in_exclusion_optimize_num(int num) { in_exclusion_optimize_num = num; }

    void set_in_exclusion_optimize_redundancy();
};
