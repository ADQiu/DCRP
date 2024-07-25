#include "cyclecovering.hh"
#include "argheap.hpp"
#include<map>
#include<iterator>
#include<algorithm>
#include<thread>
#include<iostream>
#include <chrono>
using namespace std;
using namespace dcrp;

static const double POSINF = numeric_limits<double>::infinity();
static const double NEGINF = numeric_limits<double>::lowest();


VDLNSCCP::VDLNSCCP(const vector<vector<double>> &distance, const vector<int> &lengths)
{
    lengths_ = lengths;
    sort(lengths_.begin(),lengths_.end());
    for (const auto &arr : distance)
    {
        dist_.push_back(arr);
    }
    n_ = dist_.size();
    m_ = lengths_.size();
    for(int i=0;i<n_;i++)
    {
        sorted_adj_.push_back({});
        for(int j=0;j<n_;j++)
        {
            if(j==i) continue;
            sorted_adj_.back().push_back(j);
        }
        sort(sorted_adj_.back().begin(), sorted_adj_.back().end(), [&](int k1, int k2) -> bool
        {
            return dist_[i][k1] < dist_[i][k2];
        });
    }
}

vector<vector<int>> VDLNSCCP::get_solution() const 
{
    return sol_best_.dump_tours();
}

double VDLNSCCP::get_best_objective() const 
{
    return obj_best_;
}

void VDLNSCCP::vdlns_solve(int maxdepth)
{
    auto starttime = chrono::steady_clock::now();
    double duration=0.0;
    auto sol = _init_solution();
    double obj = _evaluate(sol);
    sol_best_ = sol;
    obj_best_ = obj;
    if(outputflag)
        cout<<"Initial obj : "<<obj<<endl;
    int besthold_iter = 0;
    int MAX_BEST_HOLD=max(n_,200);
    double Temperature=1.0, alpha=0.98;
    
    for(int outiter=0;;outiter++)
    {
        DoubleTSPNodeList new_sol(sol);
        auto existing_tours = new_sol.dump_tours();
        int k0 = rd_.randint(m_);
        vector<int> node_assign(n_, -1);
        for (int k = 0; k < m_; k++)
        {
            for (int i : existing_tours[k])
                node_assign[i] = k;
        }
        vector<int> sizes{(int)existing_tours[k0].size()};
        set<int> involved_nodes;
        for(auto ii : existing_tours[k0])
            involved_nodes.insert(ii);
        int q=1;
        bool cur_update=false, best_update=false;
        for(;q<=maxdepth;q++)
        {
            set<int> to_add_tours;
            for(int i: involved_nodes)
            {
                for(int j:sorted_adj_[i])
                {
                    if(involved_nodes.find(j) == involved_nodes.end())
                    {
                        to_add_tours.insert(node_assign[j]);
                        break;
                    }
                }
            }

            if(to_add_tours.size()==0)
                break;
            
            for(int k:to_add_tours)
            {
                sizes.push_back((int)existing_tours[k].size());
                involved_nodes.insert(existing_tours[k].begin(), existing_tours[k].end());
            }

            rd_.shuffle(sizes);
            _reroute_subtour(new_sol, vector<int>(involved_nodes.begin(), involved_nodes.end()));
            _split(new_sol, existing_tours[k0][0], sizes);
            double new_obj = _evaluate(new_sol);
            if(new_obj < obj_best_-1e-6)
            {
                sol_best_ = new_sol;
                obj_best_ = new_obj;
                if(outputflag)
                    cout<<"BestUpdate "<<new_obj<< " at iteration " << outiter <<endl;
                best_update=true;
            }
            double t=new_obj-obj;
            if(t<0 || ((std::exp(-t/Temperature-0.1)>rd_.random())))
            {
                sol = new_sol;
                obj = new_obj;
                cur_update=true;
                break;
            }
        }

        if(best_update)
            besthold_iter=0;
        else 
            besthold_iter++;

        if(besthold_iter>=MAX_BEST_HOLD)
            break;

        auto endtime = chrono::steady_clock::now();
        duration = chrono::duration_cast<chrono::milliseconds>(endtime-starttime).count()/1000.0;
        if(duration>MAX_TIME)
            break;
        
        Temperature*=alpha;
    }
}

void VDLNSCCP::lns_solve()
{
    auto starttime = chrono::steady_clock::now();
    double duration=0.0;
    auto sol = _init_solution();
    double obj = _evaluate(sol);
    sol_best_ = sol;
    obj_best_ = obj;
    if(outputflag)
        cout<<"Initial obj : "<<obj<<endl;
    int besthold_iter = 0;
    int MAX_BEST_HOLD=max(n_,200);
    for(int outiter=0;;outiter++)
    {
        DoubleTSPNodeList new_sol(sol);
        _destroy_repair(new_sol,neighbor);

        double new_obj = _evaluate(new_sol);
        if(new_obj<obj_best_-1e-6)
        {
            sol_best_ = new_sol;
            obj_best_ = new_obj;
            if(outputflag)
                cout<<"BestUpdate "<<new_obj<< " at iteration " << outiter <<endl;
            besthold_iter = 0;
        }
        else 
            besthold_iter++;
        if(besthold_iter>=MAX_BEST_HOLD)
            break;
        auto endtime = chrono::steady_clock::now();
        
        duration = chrono::duration_cast<chrono::milliseconds>(endtime-starttime).count()/1000.0;
        if(duration>MAX_TIME)
            break;

        double t = new_obj - obj;
        if(t<0 || std::exp(-t)>rd_.random())
        {
            sol = new_sol;
            obj = new_obj;
        }
    }
}

void VDLNSCCP::CS_simple_solve()
{
    auto starttime = chrono::steady_clock::now();
    double duration=0.0;
    auto sol = _init_solution();
    double obj = _evaluate(sol);
    sol_best_ = sol;
    obj_best_ = obj;
}

double VDLNSCCP::_reroute_subtour(DoubleTSPNodeList& sol, int root)
{
    double old_val = _evaluate_subtour(sol,root);
    vector<int> nodes;
    int u=root;
    do 
    {
        nodes.push_back(u);
        u=sol.next(u);
    } while(u!=root);
    
    vector<int> newtour;
    TSP_greedy(dist_, nodes, newtour);
    for(int i=0;i<newtour.size();i++)
    {
        sol.bind(newtour[i], newtour[(i+1)%newtour.size()]);
    }
    _VOPT_update(sol,root);
    double new_val = _evaluate_subtour(sol,root);
    return old_val-new_val;
}

void VDLNSCCP::_reroute_subtour(DoubleTSPNodeList& sol, const vector<int>& candidates)
{
    vector<int> newtour;
    TSP_greedy(dist_, candidates, newtour);
    for(int i=0;i<newtour.size();i++)
    {
        sol.bind(newtour[i], newtour[(i+1)%newtour.size()]);
    }
    _VOPT_update(sol,candidates[0]);

}

void VDLNSCCP::_tours_near_merge_reroute_split_update(DoubleTSPNodeList& sol)
{
    auto existing_tours = sol.dump_tours();
    vector<int> node_assign(n_,-1);
    for(int k=0;k<m_;k++)
    {
        for(int i: existing_tours[k])
            node_assign[i] = k;
    }

    vector<double> old_avg_costs;
    for(auto& tour: existing_tours)
    {
        old_avg_costs.push_back(_evaluate_subtour(sol,tour[0])/tour.size());
    }
    double min_avg_cost = *min_element(old_avg_costs.begin(),old_avg_costs.end());

    vector<int> candidate_tour_indices;
    for(int k=0;k<m_;k++)
    {
        if(old_avg_costs[k] > min_avg_cost+1e-6)
            candidate_tour_indices.push_back(k);
    }

    if(candidate_tour_indices.size() == 0)
        return ;

    int k0 = candidate_tour_indices[rd_.randint(candidate_tour_indices.size())];
    vector<int> sizes{(int)existing_tours[k0].size()};
    vector<int> involved_nodes(existing_tours[k0]);
    set<int> involved_tours{k0};
    vector<int> last_added_nodes(existing_tours[k0]);
    set<int> to_add_tours;
    double old_cost = _evaluate(sol);
    double new_cost = old_cost;
    int retry=0;

TRYAGAIN:
    retry++;
    to_add_tours.clear();
    for(int i: last_added_nodes)
    {
        for(int j: sorted_adj_[i])
        {
            if(involved_tours.find(node_assign[j]) == involved_tours.end())
            {
                to_add_tours.insert(node_assign[j]);
                break;
            }
        }
    }

    if(to_add_tours.size() == 0)
    {
        return;
    }    

    last_added_nodes.clear();

    for(int k: to_add_tours)
    {
        sizes.push_back(existing_tours[k].size());
        last_added_nodes.insert(last_added_nodes.end(),existing_tours[k].begin(),existing_tours[k].end());
        involved_nodes.insert(involved_nodes.end(),existing_tours[k].begin(),existing_tours[k].end());
        involved_tours.insert(k);
    }
    
    rd_.shuffle(sizes);
    _reroute_subtour(sol,involved_nodes);
    _split(sol, existing_tours[k0][0], sizes);
    new_cost = _evaluate(sol);
    if(retry > 2)
        return;
    if(new_cost >= old_cost)
        goto TRYAGAIN;
}   

/**
 * @brief 
 *  每次只从k_0的近邻添加新近邻tour
 * @param sol 
 */
void VDLNSCCP::_mild_subnear_tour_merge_reroute_split_update(DoubleTSPNodeList& sol)
{
    auto existing_tours = sol.dump_tours();
    vector<int> node_assign(n_,-1);
    for(int k=0;k<m_;k++)
    {
        for(int i: existing_tours[k])
            node_assign[i] = k;
    }

    vector<double> old_avg_costs;
    for(auto& tour: existing_tours)
    {
        old_avg_costs.push_back(_evaluate_subtour(sol,tour[0])/tour.size());
    }
    double min_avg_cost = *min_element(old_avg_costs.begin(),old_avg_costs.end());

    vector<int> candidate_tour_indices;
    for(int k=0;k<m_;k++)
    {
        if(old_avg_costs[k] > min_avg_cost+1e-6)
            candidate_tour_indices.push_back(k);
    }

    if(candidate_tour_indices.size() == 0)
        return ;

    int k0 = candidate_tour_indices[rd_.randint(candidate_tour_indices.size())];
    int nk0 = existing_tours[k0].size();
    vector<int> sizes{(int)existing_tours[k0].size()};
    vector<int> involved_nodes(existing_tours[k0]);
    set<int> involved_tours{k0};
    set<int> to_add_tours;
    double old_cost = _evaluate(sol);
    double new_cost = old_cost;
    int retry=0;
    vector<vector<int>::iterator> adj_iters;
    for(int i=0;i<nk0;i++)
    {
        adj_iters.push_back(sorted_adj_[existing_tours[k0][i]].begin());
    }

TRYAGAIN:
    retry++;
    to_add_tours.clear();
    for(int ng=0;ng<nk0;ng++)
    {
        while(adj_iters[ng]!=sorted_adj_[existing_tours[k0][ng]].end())
        {
            int j = *adj_iters[ng];
            if(involved_tours.find(node_assign[j]) == involved_tours.end())
            {
                to_add_tours.insert(node_assign[j]);
                adj_iters[ng]++;
                break;   
            }
            adj_iters[ng]++;
        }
    }

    if(to_add_tours.size() == 0)
    {
        return;
    }

    for(int k: to_add_tours)
    {
        sizes.push_back(existing_tours[k].size());
        involved_nodes.insert(involved_nodes.end(),existing_tours[k].begin(),existing_tours[k].end());
        involved_tours.insert(k);
    }

    rd_.shuffle(sizes);
    _reroute_subtour(sol,involved_nodes);
    _split(sol, existing_tours[k0][0], sizes);
    new_cost = _evaluate(sol);
    if(retry > 2)
        return;
    if(new_cost >= old_cost)
        goto TRYAGAIN;
}

void VDLNSCCP::_nearest_merge_reroute_split_update(DoubleTSPNodeList& sol)
{
    int u;

    cross_queue_.clear();
    auto order = rd_.randints(n_,n_,false);
    for(int i=0;i<n_;i++)
    {
        cross_queue_.push_back(order[i]);
    }
    while (cross_queue_.size() > 0)
    {
        u = cross_queue_.front();
        cross_queue_.pop_front();

        int unext = sol.next(u);
        double w = dist_[u][unext];
        vector<int> candidates;
        for(int v:sorted_adj_[u])
        {
            if (dist_[u][v] > w + 1e-6)
                break;
            if (sol.same_tour(u, v))
                continue;
            candidates.push_back(v);
        }
        
        if(candidates.size()==0) continue;

        vector<int> sizes{sol.tour_len_within(u)};
        vector<int> involved(sol.tour_within(u));
        double old_cost = _evaluate(sol);
        double new_cost = old_cost;
        while(
            old_cost - new_cost <1e-6 && 
            candidates.size()>0&&
            sizes.size()<m_&&
            sizes.size()<5)
        {
            int k = rd_.randint(0, candidates.size() - 1);
            int v = candidates[k];
            candidates.erase(candidates.begin() + k);
            for(int k=candidates.size()-1;k>=0;k--)
            {
                if(sol.same_tour(v,candidates[k]))
                {
                    candidates.erase(candidates.begin()+k);
                }
            }
            auto vtour = sol.tour_within(v);
            sizes.push_back(vtour.size());
            involved.insert(involved.end(),vtour.begin(),vtour.end());

            _reroute_subtour(sol, involved);
            _split(sol, u, sizes);
            new_cost = _evaluate(sol);
        }
        // }
    }
}

static void descending_insert_with_key(vector<int>& v, const vector<double>& key)
{
    v.clear();
    int n=key.size();
    for(int i=0;i<n;i++)
    {
        if(v.empty())
        {
            v.push_back(i);
            continue;
        }
        int l=0;
        int r=v.size()-1;

        if(key[v[l]]<=key[i])
        {
            v.insert(v.begin(),i);
            continue;
        }
        if(key[v[r]]>=key[i])
        {
            v.push_back(i);
            continue;
        }
        // descending insert
        while(l<r)
        {
            int mid = (l+r)/2;
            if(key[v[mid]]>=key[i])
            {
                r = mid;
            }
            else
            {
                l = mid+1;
            }
        }
        v.insert(v.begin()+l,i);
    }
}

double VDLNSCCP::_split(DoubleTSPNodeList& sol, int root, const vector<int>& sizes)
{
    int start = root;
    double gain = 0.0;
    LKHAgent agent(sol, dist_, sorted_adj_,2);
    for(int i=0;i+1<sizes.size();i++)
    {
        int t1 = start;
        int t2 = sol.next(t1);
        int t3 = t1;
        for(int count = 0;count<sizes[i];count++)
        {
            t3 = sol.next(t3);
        }
        int t4 = sol.next(t3);
        int s1,s2,s3,s4;
        double best_split=NEGINF;
        do 
        {
            if(agent.TOGGLE_gain(t1,t2,t3,t4)>best_split)
            {
                best_split = agent.TOGGLE_gain(t1,t2,t3,t4);
                s1 = t1; s2 = t2; s3 = t3; s4 = t4;
            }
            t1 = sol.next(t1);
            t2 = sol.next(t2);
            t3 = sol.next(t3);
            t4 = sol.next(t4);
        } while (t1!=start);
        gain += agent.TOGGLE_gain(s1,s2,s3,s4);
        agent.TOGGLE(s1,s2,s3,s4);
        gain += agent.lin_kernighan(s1);
        gain += agent.lin_kernighan(s3);
        start = s1;
    }
    return gain;
}

DoubleTSPNodeList VDLNSCCP::_init_solution()
{
    TSP_Annealing tsp_solver(dist_);
    tsp_solver.solve();
    DoubleTSPNodeList out = tsp_solver.get_best_solution();
    if(m_>1) out.set_single_tour(false);
    
    if(outputflag)
        cout<<"Begin with a tsp tour of length "<<_evaluate(out)<<endl;
    _split(out,0,lengths_);
    return out;
}

void VDLNSCCP::_destroy_repair(DoubleTSPNodeList& sol, NeighborType t)
{
    switch (t)
    {
    case INTURN_MERGE_REROUTE_SPLIT:
        _tours_near_merge_reroute_split_update(sol);
        break;
    case SUBNEAR_MERGE_REROUTE_SPLIT:
        _mild_subnear_tour_merge_reroute_split_update(sol);
        break;
    case NEARESTS_REROUTE_SPLIT:
        _nearest_merge_reroute_split_update(sol);
        break;
    default:
        _mild_subnear_tour_merge_reroute_split_update(sol);
        break;
    }
}

double VDLNSCCP::_evaluate(const DoubleTSPNodeList& sol) const 
{
    // set<int> unvisited;
    // for(int i=0;i<n_;i++)
    // {
    //     unvisited.insert(i);
    // }
    // double out = 0.0;
    // while(!unvisited.empty())
    // {
    //     int root = *unvisited.begin();
    //     int u = root;
    //     int l=0;
    //     do 
    //     {  
    //         unvisited.erase(u);
    //         l++;
    //         out += dist_[u][sol.next(u)];
    //         u=sol.next(u);
    //     } while (u!=root);
    // }
    // return out;
    double out=0.0;
    for(int i=0;i<n_;i++)
    {
        out += dist_[i][sol.next(i)];
    }
    return out;
}

double VDLNSCCP::_evaluate_subtour(const DoubleTSPNodeList& sol, int root) const 
{
    int u = root;
    double out=0.0;
    do 
    {
        out += dist_[u][sol.next(u)];
        u=sol.next(u);
    } while (u!=root);
    return out;
}

double VDLNSCCP::_close_eval(int u0, int u1, int v0, int v1) const
{
    return (dist_[u0][v1] + dist_[v0][u1]) - (dist_[u0][u1] + dist_[v0][v1]);
}

double VDLNSCCP::_VOPT_update(DoubleTSPNodeList& sol, int root)
{
    LKHAgent lkh(sol, dist_, sorted_adj_,2);
    double val = lkh.lin_kernighan(root);
    return val;
}

double VDLNSCCP::_VOPT_update(DoubleTSPNodeList& sol)
{
    auto tours = sol.dump_tours();
    double gain=0.0;
    for(auto& t: tours)
    {
        gain+=_VOPT_update(sol,t[0]);
    }
    return gain;
}