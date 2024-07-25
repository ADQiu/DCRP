#ifndef DCRP_TRAVELING_SALESMAN_PROBLEM_HH
#define DCRP_TRAVELING_SALESMAN_PROBLEM_HH
#include <vector>
#include <list>
#include <set>
#include <limits>
#include "point.hh"
#include "simplegraph.hh"
#include "random.hh"
namespace dcrp
{
    struct TSPNode
    {
        int name;
        int adj[2];
        TSPNode(int i): name(i){adj[0] = adj[1] = i;}
        TSPNode(const TSPNode& node)
        {
            name = node.name;
            adj[0] = node.adj[0];
            adj[1] = node.adj[1];
        }
        ~TSPNode(){}
        void exchange() { std::swap(adj[0], adj[1]); }
    };

    class DoubleTSPNodeList
    {
    public:
        DoubleTSPNodeList() : n_(0) {}
        DoubleTSPNodeList(int n);
        DoubleTSPNodeList(const DoubleTSPNodeList& );
        DoubleTSPNodeList(const std::vector<int>& arr);
        ~DoubleTSPNodeList() = default;

        void             load      (const std::vector<int>& arr);
        std::vector<int> dump      () const;
        inline int       size      () const        { return n_;}
        inline TSPNode&  operator[](int i)         { return nodes_[i];}
        inline const TSPNode& operator[](int i) const { return nodes_[i];}
        inline int       next      (int i) const   { return nodes_[i].adj[1]; }
        inline int&      next      (int i)         { return nodes_[i].adj[1]; }
        inline int       prev      (int i) const   { return nodes_[i].adj[0]; }
        inline int&      prev      (int i)         { return nodes_[i].adj[0]; }
        void             flip      (int x, int y);
        void             insert    (int x, int mid, int y);
        void             drop      (int x, int mid, int y);
        void             toggle(int u1, int u2, int v1, int v2);
        std::vector<int> tour_within(int i) const;
        int              tour_len_within(int i) const;
        void             reverse_subtour(int i);
        bool             sequential(int i, int j, int k) const;
        bool             same_tour(int i, int j) const;
        std::vector<std::vector<int>> dump_tours() const;
        inline void      bind      (int i, int j) { nodes_[i].adj[1] = j; nodes_[j].adj[0] = i; }
        void             set_single_tour(bool flag) {single_tour_ = flag;}

    protected:
        std::vector<TSPNode> nodes_;
        int n_;
        bool single_tour_ = true;
    };

    bool adj_valid(const DoubleTSPNodeList& tours);
    bool has_subtour(const DoubleTSPNodeList& tours);
    bool has_outer_into(const DoubleTSPNodeList& tours);

    struct Edgelook
    {
        int other;
        double diff;
        int over;
        int mm;
        Edgelook() {diff=std::numeric_limits<double>::infinity();}
        Edgelook(int other, double diff, int over) : 
            other(other), diff(diff), over(over) {}
    };

    enum class LKH_act_type
    {
        FLIP,
        SEPARATE
    };

    struct LKH_action
    {
        int t1,t2,t3,t4;
        LKH_act_type type;
        LKH_action(int t1, int t2, int t3, int t4, LKH_act_type type) : 
            t1(t1), t2(t2), t3(t3), t4(t4), type(type) {}
    };

    class LKHAgent
    {
    public:
        LKHAgent(
            DoubleTSPNodeList &sol, 
            const std::vector<std::vector<double>> &dist, 
            const std::vector<std::vector<int>> &sorted_adj,
            int depth=4)
            : sol_(sol), dist_(dist), sorted_adj_(sorted_adj), depth_(depth) {}
        ~LKHAgent() {}
        void FLIP(int t1, int t2, int t3, int t4)
        {
            sol_.flip(t2, t3);
            actions_.push_back(LKH_action(t1,t2,t3,t4,LKH_act_type::FLIP));
        }
        void UNFLIP(int t1, int t2, int t3, int t4)
        {
            sol_.flip(t3, t2);
            actions_.pop_back();
        }
        double FLIP_gain(int t1, int t2, int t3, int t4) const
        {
            return dist_[t1][t2] + dist_[t3][t4] - dist_[t1][t3] - dist_[t2][t4];
        }
        virtual void TOGGLE(int t1, int t2, int t3, int t4)
        {
            sol_.toggle(t1, t2, t3, t4);
            actions_.push_back(LKH_action(t1,t2,t3,t4,LKH_act_type::SEPARATE));
        }
        virtual void UNTOGGLE(int t1, int t2, int t3, int t4)
        {
            sol_.toggle(t1, t4, t3, t2);
            actions_.pop_back();
        }
        double TOGGLE_gain(int t1, int t2, int t3, int t4) const
        {
            return dist_[t1][t2] + dist_[t3][t4] - dist_[t1][t4] - dist_[t3][t2];
        }
        
    protected:
        int inline key(int x, int y) const
        {
            return std::min(x, y) * sol_.size() + std::max(x, y);
        }
        void inline mark_deleted(int x, int y)
        {
            edge_deleted_.insert(key(x, y));
        }
        void inline mark_added(int x, int y)
        {
            edge_added_.insert(key(x, y));
        }
        void inline unmark_deleted(int x, int y)
        {
            edge_deleted_.erase(key(x, y));
        }
        void inline unmark_added(int x, int y)
        {
            edge_added_.erase(key(x, y));
        }
        bool inline has_deleted(int x, int y) const
        {
            return edge_deleted_.find(key(x, y)) != edge_deleted_.end();
        }
        bool inline has_added(int x, int y) const
        {
            return edge_added_.find(key(x, y)) != edge_added_.end();
        }

        void look_ahead(double gain, int first, int last, int level, std::list<Edgelook> &edgelist);
        double step(double gain, int first, int last, int level, double &Gstar);
        double improve_tour(int t1);

    public:
        double lin_kernighan(int root);
        virtual void recover();
        void clear_cache()
        {
            edge_deleted_.clear();
            edge_added_.clear();
            actions_.clear();
            queue_.clear();
        }
        
        DoubleTSPNodeList &sol_;
        const std::vector<std::vector<double>> &dist_;
        const std::vector<std::vector<int>> &sorted_adj_;

    protected:
        std::list<LKH_action> actions_;
        std::list<int> queue_;
        std::set<int> edge_deleted_, edge_added_;
        int depth_;
    };

    class TSP_Annealing
    {
    public:
        TSP_Annealing(const std::vector<std::vector<double>> &distance);
        ~TSP_Annealing();

        void solve();
        
        void set_init_tempe(double tempe) { INIT_TEMPE_ = tempe; }
        void set_max_iteration(size_t max_iter) { MAX_ITER_ = max_iter; }

        DoubleTSPNodeList get_best_solution() const { return sol_best_; }
        double get_best_objective() const { return obj_best_; }
        void print_stats();

    protected:
        double INIT_TEMPE_ = 10.0;
        double tempe_cur_;
        size_t MAX_ITER_ = 100;
        size_t cur_iter_ = 0;
        size_t no_best_update_ = 0;
        size_t no_update_ = 0;
        Random rd_;
        double obj_best_;
        DoubleTSPNodeList sol_best_;

        unsigned update_stats_ = 0;
        unsigned best_update_stats_ = 0;
        const std::vector<std::vector<double>> &distance_;
        std::vector<std::vector<int>> sorted_adj_;

        double _evaluate(const DoubleTSPNodeList &sol) const;
        DoubleTSPNodeList _init_solution();
        DoubleTSPNodeList _generate_new_solution(const DoubleTSPNodeList &sol);
        double _accept_probability(double newone, double oldone) const;

        virtual double _temperature_update(double tempe) { return tempe * 0.7; }
        virtual bool _should_temperature_update() { return (double)no_best_update_ > MAX_ITER_ / 10; }

        inline bool _has_cross(int u1, int u2, int v1, int v2) const;
        void _cross_update(DoubleTSPNodeList &sol);

        inline double weight(int i, int j) const { return distance_[i][j]; }
    };

    double TSP_greedy(
        const std::vector<std::vector<double>> &dist,
        std::vector<int> &best_tour);
    double TSP_greedy(
        const std::vector<std::vector<double>> &dist,
        const std::vector<int> &candidates,
        std::vector<int> &best_tour);
}
#endif //DCRP_TRAVELING_SALESMAN_PROBLEM_HH