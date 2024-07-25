#ifndef DCRP_CCP_HH
#define DCRP_CCP_HH
#include <vector>
#include <map>
#include "tsp.hh"
#include "point.hh"


namespace dcrp
{
    
    enum NeighborType
    {
        INTURN_MERGE_REROUTE_SPLIT=0,
        SUBNEAR_MERGE_REROUTE_SPLIT,
        NEARESTS_REROUTE_SPLIT
    };

    class VDLNSCCP
    {
    public:
        VDLNSCCP(const std::vector<std::vector<double>>& distance, const std::vector<int>& lengths);
        ~VDLNSCCP() {}

        void lns_solve();
        void vdlns_solve(int maxdep=3);

        void CS_simple_solve();
        
        std::vector<std::vector<int>> get_solution() const;
        double get_best_objective() const;

        void set_outer_max_iter(int num) {MAX_ITER_ = num;}

        NeighborType neighbor = INTURN_MERGE_REROUTE_SPLIT;
        bool outputflag=false;

    protected:
        std::vector<std::vector<double>> dist_;
        std::vector<int> lengths_; 
        std::vector<std::vector<int>> sorted_adj_;
        DoubleTSPNodeList sol_best_;
        Random rd_;
        int n_,m_;
        double obj_best_;
        int MAX_ITER_ = 20000;
        int iter_cur_ = 0;
        double MAX_TIME = 600;

        std::set<int> edge_deleted_, edge_added_;
        std::list<int> lkhqueue_, cross_queue_;

        DoubleTSPNodeList _init_solution();
        void _destroy_repair(DoubleTSPNodeList& sol, NeighborType t);
        double _evaluate(const DoubleTSPNodeList&) const;
        double _evaluate_subtour(const DoubleTSPNodeList&, int root) const;
        
        // 邻域算子
        void _nearest_merge_reroute_split_update(DoubleTSPNodeList& sol);
        void _tours_near_merge_reroute_split_update(DoubleTSPNodeList& sol);
        void _mild_subnear_tour_merge_reroute_split_update(DoubleTSPNodeList& sol);

        // 辅助函数
        double _reroute_subtour(DoubleTSPNodeList& sol, int root);
        void _reroute_subtour(DoubleTSPNodeList& sol, const std::vector<int>& candidates);
        double _split(DoubleTSPNodeList& sol, int root, const std::vector<int>& sizes);
        double _close_eval(int u0, int u1, int v0, int v1) const;

        double _VOPT_update(DoubleTSPNodeList&, int root);
        double _VOPT_update(DoubleTSPNodeList&);
    };
} 
#endif // DCRP_CCP_HH
