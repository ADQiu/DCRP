#ifndef DCRP_CYCLECOVERINGGRB_HH
#define DCRP_CYCLECOVERINGGRB_HH
#include <vector>
#include <map>
#include "point.hh"
#include "gurobi_c++.h"

namespace dcrp
{
    class CCPGRB
    {
    public: 
        CCPGRB(const std::vector<std::vector<double>>& dist, const std::vector<int>& lengths);
        ~CCPGRB();

        void solve();
        
        std::vector<std::vector<int>> get_solution() const {return best_sol;}
        double get_best_objective() const {return best_objective;}
    
    private:
        int n,m;
        std::vector<std::vector<double>> dist_;
        std::vector<int> lengths_;
        GRBEnv* env_=nullptr;
        GRBModel* model_=nullptr;
        std::vector<std::vector<int>> best_sol;
        double best_objective=0.0;

        void create_model();
        std::vector<std::vector<std::vector<GRBVar>>> vars;
        void find_tours();
    };

} 
#endif  