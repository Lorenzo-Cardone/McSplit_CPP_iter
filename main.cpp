#include "graph.h"

#include <algorithm>
#include <functional>
#include <numeric>
#include <chrono>
#include <iostream>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <mutex>
#include <thread>
#include <condition_variable>
#include <atomic>
#include <map>
#include <list>
#include <cassert>
#include <cstdlib>
#include <argp.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>


using std::vector;
using std::cout;
using std::endl;

using std::chrono::steady_clock;
using std::chrono::milliseconds;
using std::chrono::duration_cast;

static void fail(std::string msg) {
    std::cerr << msg << std::endl;
    exit(1);
}

enum Heuristic { min_max, min_product };

/*******************************************************************************
                             Command-line arguments
*******************************************************************************/

static char doc[] = "Find a maximum clique in a graph in DIMACS format\vHEURISTIC can be min_max or min_product";
static char args_doc[] = "HEURISTIC FILENAME1 FILENAME2";
static struct argp_option options[] = {
    {"quiet", 'q', 0, 0, "Quiet output"},
    {"verbose", 'v', 0, 0, "Verbose output"},
    {"dimacs", 'd', 0, 0, "Read DIMACS format"},
    {"lad", 'l', 0, 0, "Read LAD format"},
    {"connected", 'c', 0, 0, "Solve max common CONNECTED subgraph problem"},
    {"directed", 'i', 0, 0, "Use directed graphs"},
    {"labelled", 'a', 0, 0, "Use edge and vertex labels"},
    {"vertex-labelled-only", 'x', 0, 0, "Use vertex labels, but not edge labels"},
    {"big-first", 'b', 0, 0, "First try to find an induced subgraph isomorphism, then decrement the target size"},
    {"timeout", 't', "timeout", 0, "Specify a timeout (seconds)"},
    {"threads", 'T', "threads", 0, "Specify how many threads to use"},
    {"randomize_seed", 'r', "randomize_seed", 0, "Randomize the order of the nodes (default=0 for no randomization)"},
    {"new_solver", 'n', 0, 0, "Use the new solver implementation"},
    { 0 }
};

static struct {
    bool quiet;
    bool verbose;
    bool dimacs;
    bool lad;
    bool connected;
    bool directed;
    bool edge_labelled;
    bool vertex_labelled;
    bool big_first;
    bool new_solver;
    Heuristic heuristic;
    char *filename1;
    char *filename2;
    int timeout;
    int threads;
    int arg_num;
    size_t random_seed = 0;
} arguments;

static std::atomic<bool> abort_due_to_timeout;
std::atomic<size_t> global_position {0};

void set_default_arguments() {
    arguments.quiet = false;
    arguments.verbose = false;
    arguments.dimacs = false;
    arguments.lad = false;
    arguments.connected = false;
    arguments.directed = false;
    arguments.edge_labelled = false;
    arguments.vertex_labelled = false;
    arguments.big_first = false;
    arguments.filename1 = NULL;
    arguments.filename2 = NULL;
    arguments.timeout = 0;
    arguments.threads = std::thread::hardware_concurrency();
    arguments.arg_num = 0;
    arguments.new_solver = false;
}

static error_t parse_opt (int key, char *arg, struct argp_state *state) {
    switch (key) {
        case 'd':
            if (arguments.lad)
                fail("The -d and -l options cannot be used together.\n");
            arguments.dimacs = true;
            break;
        case 'l':
            if (arguments.dimacs)
                fail("The -d and -l options cannot be used together.\n");
            arguments.lad = true;
            break;
        case 'q':
            arguments.quiet = true;
            break;
        case 'v':
            arguments.verbose = true;
            break;
        case 'c':
            //if (arguments.directed)
            //    fail("The connected and directed options can't be used together.");
            arguments.connected = true;
            break;
        case 'i':
            //if (arguments.connected)
            //    fail("The connected and directed options can't be used together.");
            arguments.directed = true;
            break;
        case 'a':
            if (arguments.vertex_labelled)
                fail("The -a and -x options can't be used together.");
            arguments.edge_labelled = true;
            arguments.vertex_labelled = true;
            break;
        case 'x':
            if (arguments.edge_labelled)
                fail("The -a and -x options can't be used together.");
            arguments.vertex_labelled = true;
            break;
        case 'b':
            arguments.big_first = true;
            break;
        case 't':
            arguments.timeout = std::stoi(arg);
            break;
        case 'T':
            arguments.threads = std::stoi(arg);
            break;
        case 'r':
            arguments.random_seed = std::stoul(arg);
            break;
        case 'n':
            arguments.new_solver = true;
            break;
        case ARGP_KEY_ARG:
            if (arguments.arg_num == 0) {
                if (std::string(arg) == "min_max")
                    arguments.heuristic = min_max;
                else if (std::string(arg) == "min_product")
                    arguments.heuristic = min_product;
                else
                    fail("Unknown heuristic (try min_max or min_product)");
            } else if (arguments.arg_num == 1) {
                arguments.filename1 = arg;
            } else if (arguments.arg_num == 2) {
                arguments.filename2 = arg;
            } else {
                argp_usage(state);
            }
            arguments.arg_num++;
            break;
        case ARGP_KEY_END:
            if (arguments.arg_num == 0)
                argp_usage(state);
            break;
        default: return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc };

/*******************************************************************************
                                 MCS functions
*******************************************************************************/

struct VtxPair {
    int v;
    int w;
    VtxPair(int v, int w): v(v), w(w) {}
};

struct Bidomain {
    int l,        r;        // start indices of left and right sets
    int left_len, right_len;
    bool is_adjacent;
    Bidomain(int l, int r, int left_len, int right_len, bool is_adjacent):
            l(l),
            r(r),
            left_len (left_len),
            right_len (right_len),
            is_adjacent (is_adjacent) { };
};

struct AtomicIncumbent{
    std::atomic<unsigned> value;
    uint64_t pairing_tested = 0;
    struct timespec s;
    std::mutex incumbent_mutex;

    AtomicIncumbent()
    {
        value.store(0, std::memory_order_seq_cst);
    }

    bool update(unsigned v, uint64_t currently_tested_pairings)
    {
        while (true) {
            unsigned cur_v = value.load(std::memory_order_seq_cst);
            if (v > cur_v) {
                std::unique_lock<std::mutex> guard(incumbent_mutex);
                if (value.compare_exchange_strong(cur_v, v, std::memory_order_seq_cst))
                {
                    pairing_tested = currently_tested_pairings;
                    clock_gettime(CLOCK_MONOTONIC, &s);
                    return true;
                }
            }
            else
                return false;
        }
    }
};

using PerThreadIncumbents = std::map<std::thread::id, vector<VtxPair> >;
struct PerThreadDataStruct {
    vector<VtxPair> best_sol;                   // best sol found
    uint64_t best_sol_nodes;                    // number of tested pairs when best sol was found (a bit useless when multithreading)
    struct timespec best_sol_time;              // time when best sol was first found (useless if best sol is not global best)
    vector<VtxPair> first_backtrack_sol;        // first sol after which we had to backtrack (a bit useless when multithreading)
    uint64_t first_backtrack_sol_nodes;         // number of tested pairs when we first had to backtrack (a bit useless when multithreading)
    struct timespec first_backtrack_sol_time;   // time when first backtrack took place (a bit useless when multithreading)
};
using PerThreadData = std::map<std::thread::id, PerThreadDataStruct>;

const constexpr int split_levels = 4;

struct Position
{
    std::array<unsigned, split_levels + 1> values;
    unsigned depth;

    Position()
    {
        std::fill(values.begin(), values.end(), 0);
        depth = 0;
    }

    bool operator< (const Position & other) const
    {
        if (depth < other.depth)
            return true;
        else if (depth > other.depth)
            return false;

        for (unsigned p = 0 ; p < split_levels + 1 ; ++p)
            if (values.at(p) < other.values.at(p))
                return true;
            else if (values.at(p) > other.values.at(p))
                return false;

        return false;
    }

    void add(unsigned d, unsigned v)
    {
        depth = d;
        if (d <= split_levels)
            values[d] = v;
    }
};
std::atomic<unsigned> global_position_index {0};

struct HelpMe
{
    struct Task
    {
        const std::function<void (unsigned long long &)> * func;
        int pending;
    };

    std::mutex general_mutex;
    std::condition_variable cv;
    std::map<Position, Task> tasks;
    std::atomic<bool> finish;

    vector<std::thread> threads;

    std::list<milliseconds> times;
    std::list<unsigned long long> nodes;

    HelpMe(int n_threads) :
        finish(false)
    {
        for (int t = 0 ; t < n_threads ; ++t)
            threads.emplace_back([this, n_threads, t] {
                    milliseconds total_work_time = milliseconds::zero();
                    unsigned long long this_thread_nodes = 0;
                    while (! finish.load()) {
                        std::unique_lock<std::mutex> guard(general_mutex);
                        bool did_something = false;
                        for (std::map<Position, Task>::iterator task = tasks.begin() ; task != tasks.end() ; ++task) {
                            // std::cout<< task->first.depth << " - ";
                             //for(int m = 0; m < task->first.values.size(); m++)
                                //std::cout << task->first.values[m] << " ";
                            //std::cout << std::endl;
                            if (task->second.func) { // whait for a function to be associated to this task by help_me_with()
                                auto f = task->second.func;
                                ++task->second.pending;
                                guard.unlock(); // so now other threads can stars executing this same function? why

                                auto start_work_time = steady_clock::now(); // local start time

                                (*f)(this_thread_nodes);

                                auto work_time = duration_cast<milliseconds>(steady_clock::now() - start_work_time);
                                total_work_time += work_time;

                                guard.lock();
                                task->second.func = nullptr; // why only now the function reference is removed?
                                if (0 == --task->second.pending) //decrement pending functions for this task
                                    cv.notify_all();            // notify all threads waiting on CV that this thread has finished his task

                                did_something = true;
                                break;
                            }
                        }

                        if ((! did_something) && (! finish.load()))
                            cv.wait(guard); // if nothing has be done in the previous if scope, just wait for it
                    }

                    std::unique_lock<std::mutex> guard(general_mutex);
                    times.push_back(total_work_time);
                    nodes.push_back(this_thread_nodes);
                    });
    }
    auto kill_workers() -> void{
        {
            std::unique_lock<std::mutex> guard(general_mutex);
            finish.store(true);
            cv.notify_all();
        }

        for (std::thread & t : threads)
            t.join();

        threads.clear();

        if (! times.empty()) {
            cout << "Thread work times";
            for (auto & t : times)
                cout << " " << t.count();
            cout << endl;
            times.clear();
        }
    }
    ~HelpMe(){
        kill_workers();
    }
    HelpMe(const HelpMe &) = delete;
    void get_help_with(
            const Position & position,
            const std::function<void (unsigned long long &)> & main_func,
            const std::function<void (unsigned long long &)> & thread_func,
            unsigned long long & main_nodes){
        std::map<Position, HelpMe::Task>::iterator task;

        {
            std::unique_lock<std::mutex> guard(general_mutex);
            auto r = tasks.emplace(position, HelpMe::Task{ &thread_func, 0 });
            assert(r.second);
            task = r.first;
            cv.notify_all();
        }

        main_func(main_nodes);

        {
            std::unique_lock<std::mutex> guard(general_mutex);
            while (0 != task->second.pending)
                cv.wait(guard);
            tasks.erase(task);
        }
    }
};

struct HelpMeNew
{
    struct Task
    {
        const std::function<void (unsigned long long &, std::vector<VtxPair>, std::vector<int>, std::vector<int>)> * func;
        std::vector<VtxPair> * current_sol;
        std::vector<int> * left;
        std::vector<int> * right;
        int pending;
    };

    std::mutex general_mutex;
    std::condition_variable cv;
    std::map<Position, HelpMeNew::Task> tasks;
    std::atomic<bool> finish;

    vector<std::thread> threads;

    std::list<milliseconds> times;
    std::list<unsigned long long> nodes;

    HelpMeNew(int n_threads) :
        finish(false)
    {
        for (int t = 0 ; t < n_threads ; ++t)
            threads.emplace_back([this, n_threads, t] {
                    milliseconds total_work_time = milliseconds::zero();
                    unsigned long long this_thread_nodes = 0;
                    while (! finish.load()) {
                        std::unique_lock<std::mutex> guard(general_mutex);
                        bool did_something = false;
                        for (std::map<Position, HelpMeNew::Task>::iterator task = tasks.begin() ; task != tasks.end() ; ++task) {
                            // std::cout<< task->first.depth << " - ";
                             //for(int m = 0; m < task->first.values.size(); m++)
                                //std::cout << task->first.values[m] << " ";
                            //std::cout << std::endl;
                            if (task->second.func) { // whait for a function to be associated to this task by help_me_with()
                                auto f = task->second.func;
                                ++task->second.pending;
                                guard.unlock(); // so now other threads can stars executing this same function? why

                                auto start_work_time = steady_clock::now(); // local start time

                                (*f)(this_thread_nodes, *(task->second.current_sol), *(task->second.left), *(task->second.right));

                                auto work_time = duration_cast<milliseconds>(steady_clock::now() - start_work_time);
                                total_work_time += work_time;

                                guard.lock();
                                task->second.func = nullptr; // why only now the function reference is removed?
                                if (0 == --task->second.pending) //decrement pending functions for this task
                                    cv.notify_all();            // notify all threads waiting on CV that this thread has finished his task

                                did_something = true;
                                break;
                            }
                        }

                        if ((! did_something) && (! finish.load()))
                            cv.wait(guard); // if nothing has be done in the previous if scope, just wait for it
                    }

                    std::unique_lock<std::mutex> guard(general_mutex);
                    times.push_back(total_work_time);
                    nodes.push_back(this_thread_nodes);
                    });
    }
    auto kill_workers() -> void{
        {
            std::unique_lock<std::mutex> guard(general_mutex);
            finish.store(true);
            cv.notify_all();
        }

        for (std::thread & t : threads)
            t.join();

        threads.clear();

        if (! times.empty()) {
            cout << "Thread work times";
            for (auto & t : times)
                cout << " " << t.count();
            cout << endl;
            times.clear();
        }
    }
    ~HelpMeNew(){
        kill_workers();
    }
    HelpMeNew(const HelpMeNew &) = delete;
    void get_help_with(
            const Position & position,
            const std::function<void ()> & main_func,
            std::function<void (unsigned long long &, std::vector<VtxPair>, std::vector<int>, std::vector<int>)> & thread_func,
            std::vector<VtxPair> &current_sol,
            std::vector<int> &left,
            std::vector<int> &right
        ){
        std::map<Position, HelpMeNew::Task>::iterator task;

        {
            std::unique_lock<std::mutex> guard(general_mutex);
            auto r = tasks.emplace(position, HelpMeNew::Task{ &thread_func, &current_sol, &left, &right, 0 });
            assert(r.second);
            task = r.first;
            cv.notify_all();
        }

        main_func();

        {
            std::unique_lock<std::mutex> guard(general_mutex);
            while (0 != task->second.pending)
                cv.wait(guard);
            tasks.erase(task);
        }
    }
};

bool check_sol(const Graph & g0, const Graph & g1 , const vector<VtxPair> & solution) {
    vector<bool> used_left(g0.n, false);
    vector<bool> used_right(g1.n, false);
    for (unsigned int i=0; i<solution.size(); i++) {
        struct VtxPair p0 = solution[i];
        if (used_left[p0.v] || used_right[p0.w])
            return false;
        used_left[p0.v] = true;
        used_right[p0.w] = true;
        if (g0.label[p0.v] != g1.label[p0.w])
            return false;
        for (unsigned int j=i+1; j<solution.size(); j++) {
            struct VtxPair p1 = solution[j];
            if ( (g0.adjset[p0.v].contains(p1.v) ? g0.adjset[p0.v].at(p1.v) : 0) != (g1.adjset[p0.w].contains(p1.w) ? g1.adjset[p0.w].at(p1.w) : 0) ) {
                return false;
            }
        }
    }
    return true;
}

int calc_bound(const vector<Bidomain>& domains) {
    int bound = 0;
    for (const Bidomain &bd : domains) {
        bound += std::min(bd.left_len, bd.right_len);
    }
    return bound;
}

int find_min_value(const vector<int>& arr, int start_idx, int len) {
    int min_v = INT_MAX;
    for (int i=0; i<len; i++)
        if (arr[start_idx + i] < min_v)
            min_v = arr[start_idx + i];
    return min_v;
}

int select_bidomain(const vector<Bidomain>& domains, const vector<int> & left,
        int current_matching_size){
    // Select the bidomain with the smallest max(leftsize, rightsize), breaking
    // ties on the smallest vertex index in the left set
    int min_size = INT_MAX;
    int min_tie_breaker = INT_MAX;
    int best = -1;
    for (unsigned int i=0; i<domains.size(); i++) {
        const Bidomain &bd = domains[i];
        if (arguments.connected && current_matching_size>0 && !bd.is_adjacent) continue;
        int len = arguments.heuristic == min_max ?
                std::max(bd.left_len, bd.right_len) :
                bd.left_len * bd.right_len;
        if (len < min_size) {
            min_size = len;
            min_tie_breaker = find_min_value(left, bd.l, bd.left_len);
            best = i;
        } else if (len == min_size) {
            int tie_breaker = find_min_value(left, bd.l, bd.left_len);
            if (tie_breaker < min_tie_breaker) {
                min_tie_breaker = tie_breaker;
                best = i;
            }
        }
    }
    return best;
}

bool select_bidomain_no_idx(vector<Bidomain>& domains, const vector<int> & left,
        int current_matching_size){
    // Select the bidomain with the smallest max(leftsize, rightsize), breaking
    // ties on the smallest vertex index in the left set
    int min_size = INT_MAX;
    int min_tie_breaker = INT_MAX;
    int best = -1;
    for (unsigned int i=0; i<domains.size(); i++) {
        const Bidomain &bd = domains[i];
        if (arguments.connected && current_matching_size>0 && !bd.is_adjacent) continue;
        int len = arguments.heuristic == min_max ?
                std::max(bd.left_len, bd.right_len) :
                bd.left_len * bd.right_len;
        if (len < min_size) {
            min_size = len;
            min_tie_breaker = find_min_value(left, bd.l, bd.left_len);
            best = i;
        } else if (len == min_size) {
            int tie_breaker = find_min_value(left, bd.l, bd.left_len);
            if (tie_breaker < min_tie_breaker) {
                min_tie_breaker = tie_breaker;
                best = i;
            }
        }
    }
    if (best == -1) {
        return false;
    }
    std::swap(domains[best], domains[domains.size() - 1]);
    return true;
}

// Returns length of left half of array
int partition(vector<int>& all_vv, int start, int len, const std::unordered_map<size_t, unsigned int> & adjrow) {
    int i=0;
    for (int j=0; j<len; j++) {
        if (adjrow.contains(all_vv[start+j])) {
            std::swap(all_vv[start+i], all_vv[start+j]);
            i++;
        }
    }
    return i;
}

// multiway is for directed and/or labelled graphs
vector<Bidomain> filter_domains(const vector<Bidomain> & d, vector<int> & left,
        vector<int> & right, const Graph & g0, const Graph & g1, int v, int w,
        bool multiway){
    vector<Bidomain> new_d;
    new_d.reserve(d.size());
    for (const Bidomain &old_bd : d) {
        int l = old_bd.l;
        int r = old_bd.r;
        // After these two partitions, left_len and right_len are the lengths of the
        // arrays of vertices with edges from v or w (int the directed case, edges
        // either from or to v or w)
        int left_len = partition(left, l, old_bd.left_len, g0.adjset[v]);
        int right_len = partition(right, r, old_bd.right_len, g1.adjset[w]);
        int left_len_noedge = old_bd.left_len - left_len;
        int right_len_noedge = old_bd.right_len - right_len;
        if (left_len_noedge && right_len_noedge)
            new_d.push_back({l+left_len, r+right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent});
        if (multiway && left_len && right_len) {
            auto& adjrow_v = g0.adjset[v];
            auto& adjrow_w = g1.adjset[w];
            auto l_begin = std::begin(left) + l;
            auto r_begin = std::begin(right) + r;
            std::sort(l_begin, l_begin+left_len, [&](int a, int b)
                    { return (adjrow_v.contains(a) ? adjrow_v.at(a) : 0) < (adjrow_v.contains(b) ? adjrow_v.at(b) : 0); });
            std::sort(r_begin, r_begin+right_len, [&](int a, int b)
                    { return (adjrow_w.contains(a) ? adjrow_w.at(a) : 0) < (adjrow_w.contains(b) ? adjrow_w.at(b) : 0); });
            int l_top = l + left_len;
            int r_top = r + right_len;
            while (l<l_top && r<r_top) {
                unsigned int left_label = (adjrow_v.contains(left[l]) ? adjrow_v.at(left[l]) : 0);
                unsigned int right_label = (adjrow_w.contains(right[r]) ? adjrow_w.at(right[r]) : 0);
                if (left_label < right_label) {
                    l++;
                } else if (left_label > right_label) {
                    r++;
                } else {
                    int lmin = l;
                    int rmin = r;
                    do { l++; } while (l<l_top && ((adjrow_v.contains(left[l]) ? adjrow_v.at(left[l]) : 0) == left_label));
                    do { r++; } while (r<r_top && ((adjrow_w.contains(right[r]) ? adjrow_w.at(right[r]) : 0)==left_label));
                    new_d.push_back({lmin, rmin, l-lmin, r-rmin, true});
                }
            }
        } else if (left_len && right_len) {
            new_d.push_back({l, r, left_len, right_len, true});
        }
    }
    return new_d;
}

// returns the index of the smallest value in arr that is >w.
// Assumption: such a value exists
// Assumption: arr contains no duplicates
// Assumption: arr has no values==INT_MAX
int index_of_next_smallest(const vector<int>& arr, int start_idx, int len, int w) {
    int idx = -1;
    int smallest = INT_MAX;
    for (int i=0; i<len; i++) {
        if (arr[start_idx + i]>w && arr[start_idx + i]<smallest) {
            smallest = arr[start_idx + i];
            idx = i;
        }
    }
    return idx;
}

void remove_vtx_from_left_domain(vector<int>& left, Bidomain& bd, int v){
    int i = 0;
    while(left[bd.l + i] != v) i++;
    std::swap(left[bd.l+i], left[bd.l+bd.left_len-1]);
    bd.left_len--;
}

void remove_bidomain(vector<Bidomain>& domains, int idx) {
    domains[idx] = domains[domains.size()-1];
    domains.pop_back();
}

uint find_smallest_and_move_to_back (vector<int> &nodes, uint start, uint end, int larger_that)
{
    uint smallest = UINT_MAX;
    uint idx_smallest = UINT_MAX;
    for (uint idx = start; idx < end; idx++) {
        if (smallest > (uint)nodes[idx] && nodes[idx] > larger_that) {
            smallest = (uint)nodes[idx];
            idx_smallest = idx;
        }
    }
    if (idx_smallest != UINT_MAX) {
        nodes[idx_smallest] = nodes [end - 1];
        nodes[end - 1] = (int)smallest;
    }
    return smallest;
}

uint solve_first_graph (vector<int> &nodes, Bidomain &bd)
{
    //println!("v: {} {}", bd.left_start, bd.left_len);
    uint vtx = find_smallest_and_move_to_back(nodes, bd.l, bd.l + bd.left_len, -1);
    bd.left_len -= 1;
    return vtx;
}

uint solve_second_graph (vector<int> &nodes, Bidomain &bd, int larger_that)
{
    //println!("w: {} {}", bd.right_start, bd.right_len);
    uint vtx = find_smallest_and_move_to_back(nodes, bd.r, bd.r + bd.right_len, larger_that);
    bd.right_len -= 1;
    return vtx;
}

void print_solution (vector<VtxPair> sol) {
	for (const auto &val : sol) {
		cout << "(" << val.v << " - " << val.w << ") ";
	}
	cout << endl;
	return;
}

void new_solve (const Graph & g0, const Graph & g1,
        vector<VtxPair> & best_sol,
        uint64_t & best_sol_nodes,
        struct timespec & best_sol_time,
        vector<VtxPair> & first_backtrack_sol,
        uint64_t & first_backtrack_sol_nodes,
        struct timespec & first_backtrack_sol_time,
        vector<Bidomain> & starting_bidomain,
        vector<int> & left, vector<int> & right,
        unsigned long long &global_nodes) 
{
    uint max_size = std::max(g0.n, g1.n) * 2;
    int depth = 0;
    vector<VtxPair> current_sol;
    vector<uint> current_bidomain = vector<uint>(max_size, 0);
    vector<vector<Bidomain>> bidomains = vector<vector<Bidomain>>();

    bidomains.emplace_back(starting_bidomain);

    //let mut bound : Vec<usize> = vec![0; max_size];
    //bound[0] = calc_bound(&bidomains[0]);

    int v = INT_MAX;
    int w = -1;

    uint bound = 0;

    //Bidomain *bd = nullptr;
    while (depth >= 0) {
        if ((depth % 2) == 0) {
            if (abort_due_to_timeout) {
                break;
            }
            //print_solution(current_sol);
            global_nodes += 1;
            //show(&current_sol, &bidomains[(depth/2) as usize], &left, &right);
            bound = current_sol.size() + calc_bound(bidomains[depth/2]);
            
            if (bound <= best_sol.size()) {
                depth -= 1;
                if (first_backtrack_sol.size() == 0) {
                    first_backtrack_sol = current_sol;
                    first_backtrack_sol_nodes = global_nodes;
                    first_backtrack_sol_time = best_sol_time;
                }
                if (depth < 0) {
                    continue;
                }
                v = current_sol.back().v;
                w = current_sol.back().w;
                current_sol.pop_back();
                bidomains.pop_back();
                bidomains[depth/2][current_bidomain[depth/2]].right_len += 1;
                continue;
            }
            w = -1;
            current_bidomain[depth/2] = select_bidomain(bidomains[depth/2], left, current_sol.size());
            if (current_bidomain[depth/2] == UINT_MAX) {
                depth -= 1;
                if (first_backtrack_sol.size() == 0) {
                    first_backtrack_sol = current_sol;
                    first_backtrack_sol_nodes = global_nodes;
                    first_backtrack_sol_time = best_sol_time;
                }
                v = current_sol.back().v;
                w = current_sol.back().w;
                current_sol.pop_back();
                bidomains.pop_back();
                bidomains[depth/2][current_bidomain[depth/2]].right_len += 1;
                //bd = &bidomains[depth/2][current_bidomain[depth/2]];
                continue;
            }
            //bd = &bidomains[depth/2][current_bidomain[depth/2]];
            v = solve_first_graph(left, bidomains[depth/2][current_bidomain[depth/2]]);
            depth += 1;
        }
        else {
            w = solve_second_graph(right, bidomains[depth/2][current_bidomain[depth/2]], w);
            if (w != -1) { 
                current_sol.emplace_back(VtxPair(v, w));

                if (current_sol.size() > best_sol.size()) {
                    //print_sol(&current_sol);
                    best_sol = current_sol;
                    best_sol_nodes = global_nodes;
                    clock_gettime(CLOCK_MONOTONIC, &best_sol_time);
                }

                //{
                //    for (int i = 0; i < current_sol.size(); i++) {
                //        cout << "|" << current_sol[i].v << " " << current_sol[i].w << "|" << " ";
                //    }
                //    cout << current_sol.size() << "/" << bound << endl;
                //}
                
                bidomains.emplace_back(filter_domains(bidomains[depth/2], left, right, g0, g1, v, w, arguments.directed || arguments.edge_labelled));

                depth += 1;
            }
            else {
                bidomains[depth/2][current_bidomain[depth/2]].right_len += 1;
                depth -= 1;
                if (first_backtrack_sol.size() == 0) {
                    first_backtrack_sol = current_sol;
                    first_backtrack_sol_nodes = global_nodes;
                    first_backtrack_sol_time = best_sol_time;
                }

                if (bidomains[depth/2][current_bidomain[depth/2]].left_len == 0) {
                    remove_bidomain(bidomains[depth/2], current_bidomain[depth/2]);
                }
            }
        }
    }
    cout << "counter: " << global_nodes << endl;
}

void new_solve_par (const Graph & g0, const Graph & g1,
        AtomicIncumbent & global_incumbent,
        PerThreadData & per_thread_data,
        vector<int> & left, vector<int> & right,
        unsigned long long &global_nodes,
        int starting_depth,
        vector<VtxPair>& current_sol,
        vector<vector<Bidomain>>& bidomains,
        HelpMeNew &help_me
    ) 
{
    int depth = starting_depth;
    PerThreadDataStruct &my_data = per_thread_data[std::this_thread::get_id()];

    int v = INT_MAX;
    int w = -1;

    uint bound = 0;

    while (depth >= starting_depth) {
        if ((depth % 2) == 0) {

            if (current_sol.size() > my_data.best_sol.size()) {
                my_data.best_sol = current_sol;
                my_data.best_sol_nodes = global_nodes;
                clock_gettime(CLOCK_MONOTONIC, &my_data.best_sol_time);
            }
            
            if (abort_due_to_timeout) {
                break;
            }
            global_nodes += 1;
            bound = current_sol.size() + calc_bound(bidomains[depth/2]);
            
            if (bound <= my_data.best_sol.size()) {
                depth -= 1;
                if (my_data.first_backtrack_sol.size() == 0) {
                    my_data.first_backtrack_sol = current_sol;
                    my_data.first_backtrack_sol_nodes = global_nodes;
                    my_data.first_backtrack_sol_time = my_data.best_sol_time;
                }
                if (depth < 0) {
                    continue;
                }
                v = current_sol.back().v;
                w = current_sol.back().w;
                current_sol.pop_back();
                bidomains.pop_back();
                bidomains[depth/2].back().right_len += 1;
                continue;
            }
            w = -1;
            bool found = select_bidomain_no_idx(bidomains[depth/2], left, current_sol.size());
            if (!found) {
                depth -= 1;
                if (my_data.first_backtrack_sol.size() == 0) {
                    my_data.first_backtrack_sol = current_sol;
                    my_data.first_backtrack_sol_nodes = global_nodes;
                    my_data.first_backtrack_sol_time = my_data.best_sol_time;
                }
                v = current_sol.back().v;
                w = current_sol.back().w;
                current_sol.pop_back();
                bidomains.pop_back();
                bidomains[depth/2].back().right_len += 1;
                continue;
            }
            v = solve_first_graph(left, bidomains[depth/2].back());
            depth += 1;
        }
        else {
            // decide if there are too many "waiting tasks" in help_me
            // if so, just proceed sequentially
            // otherwise, offload the work to help_me (might be interesting to share only if the w has a different "best match")
            if (help_me.tasks.size() >= 2*arguments.threads || bidomains[depth/2].back().right_len <= 1) {
                w = solve_second_graph(right, bidomains[depth/2].back(), w);
                if (w != -1) { 
                    current_sol.emplace_back(VtxPair(v, w));
                    
                    bidomains.emplace_back(filter_domains(bidomains[depth/2], left, right, g0, g1, v, w, arguments.directed || arguments.edge_labelled));

                    depth += 1;
                }
                else {
                    bidomains[depth/2].back().right_len += 1;
                    depth -= 1;
                    if (my_data.first_backtrack_sol.size() == 0) {
                        my_data.first_backtrack_sol = current_sol;
                        my_data.first_backtrack_sol_nodes = global_nodes;
                        my_data.first_backtrack_sol_time = my_data.best_sol_time;
                    }

                    if (bidomains[depth/2].back().left_len == 0) {
                        // remove bidomain
                        bidomains[depth/2].pop_back();
                    }
                }
            }
            else 
            {
                std::atomic<int> shared_i{ 0 };
                const int i_end = bidomains[depth/2].back().right_len + 2; /* including the null */
                vector<Bidomain> domains_to_share = bidomains[depth/2];

                std::function<void (unsigned long long &, std::vector<VtxPair>, std::vector<int>, std::vector<int>)> helper_function = [&shared_i, &g0, &g1, &global_incumbent, &per_thread_data, depth,
                                    i_end, &help_me, &domains_to_share, current_sol, left, right] (unsigned long long &global_nodes_helper, std::vector<VtxPair> help_cur_sol, std::vector<int> help_left, std::vector<int> help_right) {
                    
                    int which_i_should_i_run_next = shared_i++;

                    if (which_i_should_i_run_next >= i_end)
                        return; /* don't waste time recomputing */

                    PerThreadDataStruct &my_data = per_thread_data[std::this_thread::get_id()];
                    std::vector<std::vector<Bidomain>> help_bidomains;
                    help_bidomains.resize(depth/2);
                    help_bidomains.back() = domains_to_share;

                    int help_v = help_left[help_bidomains[depth/2].back().l + help_bidomains[depth/2].back().left_len];
                    int help_w = -1;

                    for (int i = 0; (i < i_end) && (which_i_should_i_run_next < i_end); i++) {

                        help_w = solve_second_graph(help_right, help_bidomains[depth/2].back(), help_w);

                        if (i != which_i_should_i_run_next) {
                            continue;
                        }

                        if (help_w != -1) { 
                            help_cur_sol.emplace_back(VtxPair(help_v, help_w));
                            
                            help_bidomains.emplace_back(filter_domains(help_bidomains[depth/2], help_left, help_right, g0, g1, help_v, help_w, arguments.directed || arguments.edge_labelled));

                            int help_depth = depth + 1;
                            // recursive call
                            new_solve_par(g0, g1, global_incumbent, per_thread_data, help_left, help_right, global_nodes_helper, help_depth, help_cur_sol, help_bidomains, help_me);
                        }
                        else {
                            return;
                        }

                        which_i_should_i_run_next = shared_i++;
                    }
                    
                };

                std::function<void ()> main_function = [&]() {
                    int which_i_should_i_run_next = shared_i++;

                    if (which_i_should_i_run_next >= i_end)
                        return; /* don't waste time recomputing */

                    for (int i = 0; i < i_end; i++) {
                        w = solve_second_graph(right, bidomains[depth/2].back(), w);

                        if (i != which_i_should_i_run_next) {
                            continue;
                        }

                        if (w != -1) { 
                            current_sol.emplace_back(VtxPair(v, w));
                            
                            bidomains.emplace_back(filter_domains(bidomains[depth/2], left, right, g0, g1, v, w, arguments.directed || arguments.edge_labelled));

                            depth += 1;
                            // recursive call
                            new_solve_par(g0, g1, global_incumbent, per_thread_data, left, right, global_nodes, depth, current_sol, bidomains, help_me);
                        }
                        else {
                            return;
                        }
                    }
                };

                Position position;
                position.add(0, global_position_index++);
                help_me.get_help_with(position, main_function, helper_function, current_sol, left, right);
            }
        }
    }
    cout << "counter: " << global_nodes << endl;
}

void solve_nopar(const unsigned depth, const Graph & g0, const Graph & g1,
        AtomicIncumbent & global_incumbent,
        vector<VtxPair> & my_incumbent,
        vector<VtxPair> & current, vector<Bidomain> & domains,
        vector<int> & left, vector<int> & right, const unsigned int matching_size_goal,
        unsigned long long & my_thread_nodes){
    if (abort_due_to_timeout)
        return;


    my_thread_nodes++;

    if (my_incumbent.size() < current.size()) {
        my_incumbent = current;
        global_incumbent.update(current.size(), my_thread_nodes);
    }

    unsigned int bound = current.size() + calc_bound(domains);
    if (bound <= global_incumbent.value || bound < matching_size_goal)
        return;

    if (arguments.big_first && global_incumbent.value == matching_size_goal)
        return;

    int bd_idx = select_bidomain(domains, left, current.size());
    if (bd_idx == -1)   // In the MCCS case, there may be nothing we can branch on
        return;
    Bidomain &bd = domains[bd_idx];

    bd.right_len--;
    std::atomic<int> shared_i{ 0 };

    //std::cout << "solve_nopar - " << std::this_thread::get_id() << std::endl;


    int v = find_min_value(left, bd.l, bd.left_len);
    remove_vtx_from_left_domain(left, domains[bd_idx], v);
    int w = -1;
    const int i_end = bd.right_len + 2; /* including the null */

    for (int i = 0 ; i < i_end /* not != */ ; i++) {
        if (i != i_end - 1) {
            int idx = index_of_next_smallest(right, bd.r, bd.right_len+1, w);
            w = right[bd.r + idx];

            // swap w to the end of its colour class
            right[bd.r + idx] = right[bd.r + bd.right_len];
            right[bd.r + bd.right_len] = w;

            auto new_domains = filter_domains(domains, left, right, g0, g1, v, w,
                    arguments.directed || arguments.edge_labelled);
            current.push_back(VtxPair(v, w));
            solve_nopar(depth + 1, g0, g1, global_incumbent, my_incumbent, current, new_domains, left, right, matching_size_goal, my_thread_nodes);
            current.pop_back();
        }
        else {
            // Last assign is null. Keep it in the loop to simplify parallelism.
            bd.right_len++;
            if (bd.left_len == 0)
                remove_bidomain(domains, bd_idx);

            solve_nopar(depth + 1, g0, g1, global_incumbent, my_incumbent, current, domains, left, right, matching_size_goal, my_thread_nodes);
        }
    }
}

void solve(const unsigned depth, const Graph & g0, const Graph & g1,
                AtomicIncumbent & global_incumbent,
                PerThreadIncumbents & per_thread_incumbents,
                vector<VtxPair> & current, vector<Bidomain> & domains,
                vector<int> & left, vector<int> & right, const unsigned int matching_size_goal,
                const Position & position, HelpMe & help_me, unsigned long long & my_thread_nodes){

    if (abort_due_to_timeout)
        return;
    my_thread_nodes++;
    if (per_thread_incumbents.find(std::this_thread::get_id())->second.size() < current.size()) {
        per_thread_incumbents.find(std::this_thread::get_id())->second = current;
        global_incumbent.update(current.size(), my_thread_nodes);
    }

    unsigned int bound = current.size() + calc_bound(domains);
    if (bound <= global_incumbent.value || bound < matching_size_goal)
        return;
    if (arguments.big_first && global_incumbent.value == matching_size_goal)
        return;

    int bd_idx = select_bidomain(domains, left, current.size());
    if (bd_idx == -1)   // In the MCCS case, there may be nothing we can branch on
        return;
    Bidomain &bd = domains[bd_idx];

    bd.right_len--;
    std::atomic<int> shared_i{ 0 };
    const int i_end = bd.right_len + 2; /* including the null */

    // std::cout << "solve - " << std::hash<std::thread::id>{}(std::this_thread::get_id())%100000 << std::endl;
    //usleep(std::hash<std::thread::id>{}(std::this_thread::get_id())%100000+1000000);
    // Version of the loop used by helpers
    std::function<void (unsigned long long &)> helper_function = [&shared_i, &g0, &g1, &global_incumbent, &per_thread_incumbents, &position, &depth,
        i_end, matching_size_goal, current, domains, left, right, &help_me] (unsigned long long & help_thread_nodes) {
        int which_i_should_i_run_next = shared_i++;

      // std::cout << "helper funtion - " << std::this_thread::get_id() << std::endl;

        if (which_i_should_i_run_next >= i_end)
            return; /* don't waste time recomputing */

        /* recalculate to this point */
        vector<VtxPair> help_current = current;
        vector<Bidomain> help_domains = domains;
        vector<int> help_left = left, help_right = right;

        /* rerun important stuff from before the loop */
        int help_bd_idx = select_bidomain(help_domains, help_left, help_current.size());
        if (help_bd_idx == -1)   // In the MCCS case, there may be nothing we can branch on
            return;
        Bidomain &help_bd = help_domains[help_bd_idx];

        int help_v = find_min_value(help_left, help_bd.l, help_bd.left_len);
        remove_vtx_from_left_domain(help_left, help_domains[help_bd_idx], help_v);

        int help_w = -1;

        for (int i = 0 ; i < i_end /* not != */ ; i++) {
            if (i != i_end - 1) {
                int idx = index_of_next_smallest(help_right, help_bd.r, help_bd.right_len+1, help_w);
                help_w = help_right[help_bd.r + idx];

                // swap w to the end of its colour class
                help_right[help_bd.r + idx] = help_right[help_bd.r + help_bd.right_len];
                help_right[help_bd.r + help_bd.right_len] = help_w;
                // don't perform recursions over vertices already considered by other threads
                if (i == which_i_should_i_run_next) {
                    which_i_should_i_run_next = shared_i++;
                    auto new_domains = filter_domains(help_domains, help_left, help_right, g0, g1, help_v, help_w,
                            arguments.directed || arguments.edge_labelled);
                    help_current.push_back(VtxPair(help_v, help_w));
                    if (depth > split_levels) {
                        solve_nopar(depth + 1, g0, g1, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, help_current, new_domains, help_left, help_right, matching_size_goal, help_thread_nodes);
                    }
                    else {
                        auto new_position = position;
                        new_position.add(depth, i + 1);
                        solve(depth + 1, g0, g1, global_incumbent, per_thread_incumbents, help_current, new_domains, help_left, help_right, matching_size_goal, new_position, help_me, help_thread_nodes);
                    }
                    help_current.pop_back();
                }
            }
            else {
                // Last assign is null. Keep it in the loop to simplify parallelism.
                help_bd.right_len++;
                if (help_bd.left_len == 0)
                    remove_bidomain(help_domains, help_bd_idx);

                if (i == which_i_should_i_run_next) {
                    which_i_should_i_run_next = shared_i++;
                    if (depth > split_levels) {
                        solve_nopar(depth + 1, g0, g1, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, help_current, help_domains, help_left, help_right, matching_size_goal, help_thread_nodes);
                    }
                    else {
                        auto new_position = position;
                        new_position.add(depth, i + 1);
                        solve(depth + 1, g0, g1, global_incumbent, per_thread_incumbents, help_current, help_domains, help_left, help_right, matching_size_goal, new_position, help_me, help_thread_nodes);
                    }
                }
            }
        }
    };

    // Grab this first, before advertising that we can get help
    int which_i_should_i_run_next = shared_i++;

    // Version of the loop used by the main thread
    std::function<void (unsigned long long &)> main_function = [&] (unsigned long long & main_thread_nodes) {
        int v = find_min_value(left, bd.l, bd.left_len);
        remove_vtx_from_left_domain(left, domains[bd_idx], v);
        int w = -1;
        //std::cout << "main funtion - " << std::this_thread::get_id() << std::endl;

        for (int i = 0 ; i < i_end /* not != */ ; i++) {
            if (i != i_end - 1) {
                int idx = index_of_next_smallest(right, bd.r, bd.right_len+1, w);
                w = right[bd.r + idx];

                // swap w to the end of its colour class
                right[bd.r + idx] = right[bd.r + bd.right_len];
                right[bd.r + bd.right_len] = w;

                if (i == which_i_should_i_run_next) {
                    which_i_should_i_run_next = shared_i++;
                    auto new_domains = filter_domains(domains, left, right, g0, g1, v, w,
                            arguments.directed || arguments.edge_labelled);
                    current.push_back(VtxPair(v, w));
                    if (depth > split_levels) {
                        solve_nopar(depth + 1, g0, g1, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, current, new_domains, left, right, matching_size_goal, main_thread_nodes);
                    }
                    else {
                        auto new_position = position;
                        new_position.add(depth, i + 1);
                        solve(depth + 1, g0, g1, global_incumbent, per_thread_incumbents, current, new_domains, left, right, matching_size_goal, new_position, help_me, main_thread_nodes);
                    }
                    current.pop_back();
                }
            }
            else {
                // Last assign is null. Keep it in the loop to simplify parallelism.
                bd.right_len++;
                if (bd.left_len == 0)
                    remove_bidomain(domains, bd_idx);

                if (i == which_i_should_i_run_next) {
                    which_i_should_i_run_next = shared_i++;
                    if (depth > split_levels) {
                        solve_nopar(depth + 1, g0, g1, global_incumbent, per_thread_incumbents.find(std::this_thread::get_id())->second, current, domains, left, right, matching_size_goal, main_thread_nodes);
                    }
                    else {
                        auto new_position = position;
                        new_position.add(depth, i + 1);
                        solve(depth + 1, g0, g1, global_incumbent, per_thread_incumbents, current, domains, left, right, matching_size_goal, new_position, help_me, main_thread_nodes);
                    }
                }
            }
        }
    };

    if (depth <= split_levels)
        help_me.get_help_with(position, main_function, helper_function, my_thread_nodes);
    else
        main_function(my_thread_nodes);
}

struct SolInfo {
    vector<VtxPair> sol;
    uint64_t nodes;
    struct timespec time;
};

std::pair<vector<VtxPair>, unsigned long long> mcs(const Graph & g0, const Graph & g1, SolInfo &best_sol_info, SolInfo &first_backtrack_sol_info) {
    vector<int> left;  // the buffer of vertex indices for the left partitions
    vector<int> right;  // the buffer of vertex indices for the right partitions
    //std::cout << "mcs - " << std::this_thread::get_id() << std::endl;

    srand(time(nullptr));
    auto domains = vector<Bidomain> {};

    std::set<unsigned int> left_labels;
    std::set<unsigned int> right_labels;
    for (unsigned int label : g0.label) left_labels.insert(label);
    for (unsigned int label : g1.label) right_labels.insert(label);
    std::set<unsigned int> labels;  // labels that appear in both graphs
    std::set_intersection(std::begin(left_labels),
                          std::end(left_labels),
                          std::begin(right_labels),
                          std::end(right_labels),
                          std::inserter(labels, std::begin(labels)));

    // Create a bidomain for each label that appears in both graphs
    for (unsigned int label : labels) {
        int start_l = left.size();
        int start_r = right.size();

        for (size_t i=0; i<g0.n; i++)
            if (g0.label[i]==label)
                left.push_back(i);
        for (size_t i=0; i<g1.n; i++)
            if (g1.label[i]==label)
                right.push_back(i);

        int left_len = left.size() - start_l;
        int right_len = right.size() - start_r;
        domains.push_back({start_l, start_r, left_len, right_len, false});
    }

    AtomicIncumbent global_incumbent;
    vector<VtxPair> incumbent;
    unsigned long long global_nodes = 0;
    std::atomic<unsigned long long> atomic_global_nodes{0};

    if (arguments.big_first) {
        for (size_t k=0; k<g0.n; k++) {
            unsigned int goal = g0.n - k;
            auto left_copy = left;
            auto right_copy = right;
            auto domains_copy = domains;
            vector<VtxPair> current;
            PerThreadIncumbents per_thread_incumbents;
            per_thread_incumbents.emplace(std::this_thread::get_id(), vector<VtxPair>());
            PerThreadData per_thread_data;
            per_thread_data.emplace(std::this_thread::get_id(), PerThreadDataStruct());
            Position position;
            if (!arguments.new_solver) {
                HelpMe help_me(arguments.threads - 1);
                for (auto & t : help_me.threads) {
                    per_thread_incumbents.emplace(t.get_id(), vector<VtxPair>());
                }
                new_solve(g0, g1, best_sol_info.sol, best_sol_info.nodes, best_sol_info.time, first_backtrack_sol_info.sol, first_backtrack_sol_info.nodes, first_backtrack_sol_info.time, domains_copy, left_copy, right_copy, global_nodes);
                help_me.kill_workers();
                for (auto & n : help_me.nodes) {
                    global_nodes += n;
                }
                for (auto & i : per_thread_incumbents) {
                    if (i.second.size() > incumbent.size()) {
                        incumbent = i.second;
                    }
                }
            }
            else {
                HelpMeNew help_me(arguments.threads - 1);
                for (auto & t : help_me.threads) {
                    per_thread_data.emplace(t.get_id(), PerThreadDataStruct());
                }
                vector<vector<Bidomain>> bidomains;
                bidomains.emplace_back(domains);
                new_solve_par(g0, g1, global_incumbent, per_thread_data, left, right, global_nodes, 0, current, bidomains, help_me);
                help_me.kill_workers();
                for (auto & n : help_me.nodes) {
                    global_nodes += n;
                }
                for (auto & data_pair : per_thread_data) {
                    auto & data = data_pair.second;
                    if (data.best_sol.size() > incumbent.size()) {
                        incumbent = data.best_sol;
                    }
                }
            }
            //solve(0, g0, g1, global_incumbent, per_thread_incumbents, current, domains_copy, left_copy, right_copy, goal, position, help_me, global_nodes);
            //for (auto & i : per_thread_incumbents)
            //    if (i.second.size() > incumbent.size())
            //        incumbent = i.second;
            incumbent = best_sol_info.sol;
            if (global_incumbent.value == goal || abort_due_to_timeout) break;
            if (!arguments.quiet) cout << "Upper bound: " << goal-1 << std::endl;
        }

    } else {
        vector<VtxPair> current;
        PerThreadIncumbents per_thread_incumbents;
        per_thread_incumbents.emplace(std::this_thread::get_id(), vector<VtxPair>());
        PerThreadData per_thread_data;
        per_thread_data.emplace(std::this_thread::get_id(), PerThreadDataStruct());
        Position position;
        if (!arguments.new_solver) {
            HelpMe help_me(arguments.threads - 1);
            for (auto & t : help_me.threads) {
                per_thread_incumbents.emplace(t.get_id(), vector<VtxPair>());
            }
            new_solve(g0, g1, best_sol_info.sol, best_sol_info.nodes, best_sol_info.time, first_backtrack_sol_info.sol, first_backtrack_sol_info.nodes, first_backtrack_sol_info.time, domains, left, right, global_nodes);
            help_me.kill_workers();
            for (auto & n : help_me.nodes) {
                global_nodes += n;
            }
            for (auto & i : per_thread_incumbents) {
                if (i.second.size() > incumbent.size()) {
                    incumbent = i.second;
                }
            }
        }
        else {
            HelpMeNew help_me(arguments.threads - 1);
            for (auto & t : help_me.threads) {
                per_thread_data.emplace(t.get_id(), PerThreadDataStruct());
            }
            vector<vector<Bidomain>> bidomains;
            bidomains.emplace_back(domains);
            new_solve_par(g0, g1, global_incumbent, per_thread_data, left, right, global_nodes, 0, current, bidomains, help_me);
            help_me.kill_workers();
            for (auto & n : help_me.nodes) {
                global_nodes += n;
            }
            for (auto & data_pair : per_thread_data) {
                auto & data = data_pair.second;
                if (data.best_sol.size() > incumbent.size()) {
                    incumbent = data.best_sol;
                }
            }
        }
        //solve(0, g0, g1, global_incumbent, per_thread_incumbents, current, domains, left, right, 1, position, help_me, global_nodes);
        //for (auto & i : per_thread_incumbents)
        //    if (i.second.size() > incumbent.size())
        //        incumbent = i.second;
        incumbent = best_sol_info.sol;
    }

    return { incumbent, global_nodes };
}

vector<int> calculate_degrees(const Graph & g) {
    vector<int> degree(g.n, 0);
    for (size_t v=0; v<g.n; v++) {
        for (size_t w=0; w<g.n; w++) {
            unsigned int mask = 0xFFFFu;
            if (g.adjset[v].contains(w) & mask) degree[v]++;
            if (g.adjset[v].contains(w) & ~mask) degree[v]++;  // inward edge, in directed case
        }
    }
    return degree;
}

int sum(const vector<int> & vec) {
    return std::accumulate(std::begin(vec), std::end(vec), 0);
}

int main(int argc, char** argv) {
    set_default_arguments();
    argp_parse(&argp, argc, argv, 0, 0, 0);

    if (arguments.random_seed != 0) {
        srand(arguments.random_seed);
    }

    char format = arguments.dimacs ? 'D' : arguments.lad ? 'L' : 'B';
    //std::vector<std::vector<unsigned int>> mat0 = {{0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,1,1},{0,0,0,0,0,1,0,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0},{0,0,0,0,0,0,1,0,0,1,0,0,0,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0},{1,0,0,0,0,0,1,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1},{0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0},{0,0,0,1,1,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0},{1,1,0,0,0,0,1,0,1,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0},{0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,1},{1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0},{0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,1},{0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,0},{0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0},{0,1,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0},{0,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0},{1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},{1,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,1,0},{1,0,1,1,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0},{0,0,0,0,0,0,1,0,0,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,0,0,0,0,0,0,0,0},{0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,0},{0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1},{1,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0},{1,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0}};
    //std::vector<std::vector<unsigned int>> mat1 = {{0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,1,0,0},{0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,1,0,0},{0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,1,1,0,0,0,1,0,0,1,0,0,0,0},{0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,1,0,0,0,1,0,1,0,0,0,0,0,0},{1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,1},{0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0},{0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0},{0,0,0,1,1,1,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0},{0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,1,0,0,1,0,0,0},{1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0},{0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,1,0},{0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0},{0,1,1,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0},{0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,1,1,0,0,1,0,0,0,0},{1,0,1,1,0,1,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0},{1,0,0,1,0,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0},{0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0},{1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0},{0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,1,0,1,0,0,0,0},{0,0,0,0,0,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,1},{0,0,1,1,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,1,0,0,1,0,1,0,0,0,0,0},{0,0,0,0,0,0,1,0,0,0,1,1,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0},{1,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0}};
    //struct Graph g0 = graphFromMtx(mat0);
    //struct Graph g1 = graphFromMtx(mat1);
    struct Graph g0 = readGraph(arguments.filename1, format, arguments.directed, arguments.edge_labelled, arguments.vertex_labelled);
    struct Graph g1 = readGraph(arguments.filename2, format, arguments.directed, arguments.edge_labelled, arguments.vertex_labelled);

    std::thread timeout_thread;
    std::mutex timeout_mutex;
    std::condition_variable timeout_cv;
    abort_due_to_timeout.store(false);
    bool aborted = false;

    if (0 != arguments.timeout) {
        timeout_thread = std::thread([&] {
                auto abort_time = steady_clock::now() + std::chrono::seconds(arguments.timeout);
                {
                    /* Sleep until either we've reached the time limit,
                     * or we've finished all the work. */
                    std::unique_lock<std::mutex> guard(timeout_mutex);
                    while (! abort_due_to_timeout.load()) {
                        if (std::cv_status::timeout == timeout_cv.wait_until(guard, abort_time)) {
                            /* We've woken up, and it's due to a timeout. */
                            aborted = true;
                            break;
                        }
                    }
                }
                abort_due_to_timeout.store(true);
                });
    }
struct timespec s, finish;
	float time_elapsed;
	clock_gettime(CLOCK_MONOTONIC, &s);
    //auto start = steady_clock::now();

    vector<int> g0_deg = calculate_degrees(g0);
    vector<int> g1_deg = calculate_degrees(g1);

    // As implemented here, g1_dense and g0_dense are false for all instances
    // in the Experimental Evaluation section of the paper.  Thus,
    // we always sort the vertices in descending order of degree (or total degree,
    // in the case of directed graphs.  Improvements could be made here: it would
    // be nice if the program explored exactly the same search tree if both
    // input graphs were complemented.
    vector<int> vv0(g0.n);
    std::iota(std::begin(vv0), std::end(vv0), 0);
    if (arguments.random_seed != 0) {
        std::random_shuffle(std::begin(vv0), std::end(vv0));
    }
    bool g1_dense = (size_t)sum(g1_deg) > g1.n*(g1.n-1);
    std::stable_sort(std::begin(vv0), std::end(vv0), [&](int a, int b) {
        return g1_dense ? (g0_deg[a]<g0_deg[b]) : (g0_deg[a]>g0_deg[b]);
    });
    vector<int> vv1(g1.n);
    std::iota(std::begin(vv1), std::end(vv1), 0);
    if (arguments.random_seed != 0) {
        std::random_shuffle(std::begin(vv1), std::end(vv1));
    }
    bool g0_dense = (size_t)sum(g0_deg) > g0.n*(g0.n-1);
    std::stable_sort(std::begin(vv1), std::end(vv1), [&](int a, int b) {
        return g0_dense ? (g1_deg[a]<g1_deg[b]) : (g1_deg[a]>g1_deg[b]);
    });

    struct Graph g0_sorted = induced_subgraph(g0, vv0);
    struct Graph g1_sorted = induced_subgraph(g1, vv1);

    SolInfo best_sol_info, first_backtrack_sol_info;
    std::pair<vector<VtxPair>, unsigned long long> solution = mcs(g0_sorted, g1_sorted, best_sol_info, first_backtrack_sol_info);

    // Convert to indices from original, unsorted graphs
    for (auto& vtx_pair : solution.first) {
        vtx_pair.v = vv0[vtx_pair.v];
        vtx_pair.w = vv1[vtx_pair.w];
    }
	clock_gettime(CLOCK_MONOTONIC, &finish);
	time_elapsed = (finish.tv_sec - s.tv_sec);
	time_elapsed += (finish.tv_nsec - s.tv_nsec) / 1000000000.0;



    /* Clean up the timeout thread */
    if (timeout_thread.joinable()) {
        {
            std::unique_lock<std::mutex> guard(timeout_mutex);
            abort_due_to_timeout.store(true);
            timeout_cv.notify_all();
        }
        timeout_thread.join();
    }

    if (!check_sol(g0, g1, solution.first))
        fail("*** Error: Invalid solution\n");

    cout << "Solution size " << solution.first.size() << std::endl;
    for (size_t i=0; i<g0.n; i++)
        for (size_t j=0; j<solution.first.size(); j++)
            if ((size_t)solution.first[j].v == i)
                cout << "(" << solution.first[j].v << " -> " << solution.first[j].w << ") ";
    cout << std::endl;

    cout << "Nodes:                      " << solution.second << endl;
    cout << "CPU time (ms):              " << time_elapsed << endl;



	fprintf(stdout, ">>> %ld - %015.010f\n", solution.first.size(), (double)(time_elapsed));
    //cout << ">>> " << solution.first.size() << " - " << (double)(end-begin)/CLOCKS_PER_SEC << endl;

    cout << "### " << best_sol_info.sol.size() << " " << best_sol_info.nodes << " " << (double)((best_sol_info.time.tv_sec - s.tv_sec) + (best_sol_info.time.tv_nsec - s.tv_nsec) / 1000000000.0) << endl;
    cout << "--- " << first_backtrack_sol_info.sol.size() << " " << first_backtrack_sol_info.nodes << " " << (double)((first_backtrack_sol_info.time.tv_sec - s.tv_sec) + (first_backtrack_sol_info.time.tv_nsec - s.tv_nsec) / 1000000000.0) << endl;

    if (aborted)
        cout << "TIMEOUT" << endl;

    //g0.printGraphMtx();
    //cout << endl;
    //g1.printGraphMtx();
}
