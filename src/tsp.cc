#include "tsp.hh"
// #include"graph.hh"
// #include"shortestpath.hh"
// #include"amopt/algorithms/minimum_spanning_tree.h"
#include "random.hh"
#include <limits>
#include <algorithm>
#include <map>
#include <thread>
#include <stack>
#include <list>
#define MAX(X,Y) ((X)<(Y)?(Y):(X))
#define MIN(X,Y) ((X)>(Y)?(Y):(X))

namespace dcrp
{
    using namespace std;
    
    DoubleTSPNodeList::DoubleTSPNodeList(int n) : n_(n)
    {
        for(int i=0;i<n;i++)
        {
            nodes_.push_back(TSPNode(i));
        }
    }

    DoubleTSPNodeList::DoubleTSPNodeList(const DoubleTSPNodeList& arr)
    {
        nodes_ = arr.nodes_;
        n_ = arr.n_;
        single_tour_ = arr.single_tour_;
    }

    DoubleTSPNodeList::DoubleTSPNodeList(const vector<int>& arr)
    {
        load(arr);
    }

    void DoubleTSPNodeList::load(const vector<int>& arr)
    {
        n_ = arr.size();
        nodes_.clear();
        for(int i=0;i<n_;i++)
        {
            nodes_.push_back(TSPNode(i));
        }
        for(int i=0;i<n_;i++)
        {
            nodes_[arr[i]].adj[0] = arr[(i-1+n_)%n_];
            nodes_[arr[i]].adj[1] = arr[(i+1)%n_];
        }
    }

    vector<int> DoubleTSPNodeList::dump() const 
    {
        return tour_within(0);
    }

    vector<vector<int>> DoubleTSPNodeList::dump_tours() const 
    {
        set<int> unvisited;
        for(int i=0;i<n_;i++)
        {
            unvisited.insert(i);
        }
        vector<vector<int>> out;
        while(!unvisited.empty())
        {
            int i = *unvisited.begin();
            vector<int> tour = tour_within(i);
            for(int j:tour)
            {
                unvisited.erase(j);
            }
            out.push_back(tour);
        }
        return out;
    }

    void DoubleTSPNodeList::flip(int x, int y)
    {
        int u = x;
        int xprev = prev(x);
        int ynext = next(y);
        while (u!=y)
        {
            nodes_[u].exchange();
            u = nodes_[u].adj[0];
        }
        swap(nodes_[u].adj[0],nodes_[u].adj[1]);
        swap(nodes_[x].adj[1],nodes_[y].adj[0]);
        swap(nodes_[xprev].adj[1],nodes_[ynext].adj[0]);
    }
    
    void DoubleTSPNodeList::insert(int x, int mid, int y)
    {
        nodes_[x].adj[1] = nodes_[y].adj[0] = mid;
        nodes_[mid].adj[0] = x;
        nodes_[mid].adj[1] = y;
    }
    
    void DoubleTSPNodeList::toggle(int u1, int u2, int v1, int v2)
    {
        single_tour_=false;
        bind(u1,v2);
        bind(v1,u2);
    }

    void DoubleTSPNodeList::drop(int x, int mid, int y)
    {
        nodes_[x].adj[1] = y;
        nodes_[y].adj[0] = x;
        nodes_[mid].adj[0] = nodes_[mid].adj[1] = mid;
    }

    vector<int> DoubleTSPNodeList::tour_within(int i) const
    {
        int u = i;
        vector<int> out;
        do 
        {
            out.push_back(u);
            u = next(u);
        } while (u!=i);
        return out;
    }

    int DoubleTSPNodeList::tour_len_within(int i) const 
    {
        int len = 0;
        int u = i;
        do 
        {
            ++len;
            u = next(u);
        } while(u!=i);
        return len;
    }

    void DoubleTSPNodeList::reverse_subtour(int root)
    {
        int u = root;
        do 
        {
            nodes_[u].exchange();
            u = next(u);
        } while (u!=root);
    }

    bool DoubleTSPNodeList::sequential(int i, int j, int k) const 
    {
        int u = i;
        do 
        {
            if(u==j) return true;
            else if(u==k) return false;
            u = next(u);
        } while (u!=i);
        return false;
    }

    bool DoubleTSPNodeList::same_tour(int i, int j) const
    {
        if(single_tour_) 
            return true;
        
        int u=i;
        do 
        {
            if(u==j) return true;
            u = next(u);
        } while (u != i);
        return false;
    } 

    bool adj_valid(const DoubleTSPNodeList& tours)
    {
        for(int i=0;i<tours.size();++i)
        {
            int u = tours.next(i);
            if(tours.prev(u)!=i) return false;
            u = tours.prev(i);
            if(tours.next(u)!=i) return false;
        }
        return true;
    }

    bool has_subtour(const DoubleTSPNodeList& tours)
    {
        vector<bool> visited(tours.size(), false);
        int u = 0;
        int edges_visited=0;

        visited[u]=true;
        while(edges_visited<tours.size())
        {
            u = tours.next(u);
            ++edges_visited;
            if(visited[u] && edges_visited<tours.size())
            {
                return true;
            }
            visited[u]=true;
        }
        return false;
    }

    bool has_outer_into(const DoubleTSPNodeList& tours)
    {
        vector<bool> visited(tours.size(), false);
        int edges_visited = 0;
        while(edges_visited<tours.size())
        {
            int start = -1;
            for(int i=0;i<tours.size();++i)
            {
                if(visited[i]) continue;
                start = i;
                break;
            }
            if(start==-1) return false;

            int u = start;
            visited[u]=true;
            do
            {
                u=tours.next(u);
                edges_visited++;
                if(u==start) break;
                if(visited[u]) return true;
                visited[u]=true;
            } while (u!=start);
        }
        return false;
    }
    
//////////////////////////////////////////
    #define MAX_LEVEL 4

    const static double POSINF = std::numeric_limits<double>::infinity();

    static void add_tour_to_queue(list<int>& q, const DoubleTSPNodeList& sol, int root)
    {
        int x = root;
        do 
        {
            q.push_back(x);
            x = sol.next(x);
        } while (x!=root);
    }

    static void ascend_replace(list<Edgelook>& edgelist, const Edgelook& newone)
    {
        auto it = edgelist.begin();
        while(it!=edgelist.end() && it->diff<newone.diff)
        {
            it++;
        }
        edgelist.insert(it, newone);
        edgelist.pop_back();
    }

    void LKHAgent::look_ahead(double gain, int first, int last, int level, std::list<Edgelook> &edgelist)
    {
        int ahead;
        if (level == 0)
            ahead = 4;
        else if (level == 1)
            ahead = 3;
        else if (level == 2 || level == 3)
            ahead = 2;
        else
            ahead = 1;

        edgelist.assign(ahead, Edgelook());
        int lastnext = sol_.next(last);
        for (int critic : sorted_adj_[last])
        {
            if (dist_[last][critic] > gain)
                break;
            if (critic == first || critic == lastnext || has_deleted(last, critic))
                continue;
            if (!sol_.same_tour(last, critic))
                continue;
            int prev = sol_.prev(critic);
            if (has_added(prev, critic))
                continue;
            double val = dist_[last][critic] - dist_[prev][critic];
            if (edgelist.back().diff > val + 1e-6)
            {
                ascend_replace(edgelist, Edgelook(critic, val, prev));
            }
        }

        while (!edgelist.empty() && (edgelist.back().diff == POSINF))
            edgelist.pop_back();
    }

    double LKHAgent::step(double gain, int first, int last, int level, double &Gstar)
    {
        if (level >= MAX_LEVEL || level>= depth_)
        {
            return 0;
        }

        int critic, newlast, hit = 0;
        double oldG = gain;
        double val;

        list<Edgelook> edgelist;
        look_ahead(gain, first, last, level, edgelist);
        for (auto &e : edgelist)
        {
            critic = e.other;
            newlast = e.over;

            gain = oldG - e.diff;
            val = gain - dist_[newlast][first];
            if (val > Gstar + 1e-6)
            {
                Gstar = val;
                hit++;
            }
            FLIP(first, last, newlast, critic);
            if (level < MAX_LEVEL)
            {
                mark_added(last, critic);
                mark_deleted(newlast, critic);
                hit += step(gain, first, newlast, level + 1, Gstar);
                unmark_added(last, critic);
                unmark_deleted(newlast, critic);
            }
            if (!hit)
                UNFLIP(first, last, newlast, critic);
            else
            {
                queue_.push_back(critic);
                queue_.push_back(newlast);
                return 1;
            }
        }
        return 0;
    }

    double LKHAgent::improve_tour(int t1)
    {
        int t2 = sol_.next(t1);
        double gain = dist_[t1][t2];
        double Gstar = 0.0;
        mark_deleted(t1, t2);
        step(gain, t1, t2, 0, Gstar);
        unmark_deleted(t1, t2);
        if (Gstar > 1e-6)
        {
            queue_.push_back(t1);
            queue_.push_back(t2);
        }
        return Gstar;
    }

    double LKHAgent::lin_kernighan(int root)
    {
        double totalwin = 0.0;
        queue_.clear();
        add_tour_to_queue(queue_, sol_, root);
        while (!queue_.empty())
        {
            int t1 = queue_.front();
            queue_.pop_front();
            totalwin += improve_tour(t1);
        }
        return totalwin;
    }

    void LKHAgent::recover()
    {
        for (auto it = actions_.rbegin(); it != actions_.rend(); it++)
        {
            if(it->type == LKH_act_type::FLIP)
                sol_.flip(it->t3, it->t2);
            else if(it->type == LKH_act_type::SEPARATE)
                sol_.toggle(it->t1, it->t4, it->t3, it->t2);
        }
        clear_cache();
    }

/////////////////////////////////////////

    static double oneway_greedy(
        const vector<vector<double>> &dist_,
        const vector<int> &candidates,
        const map<int, vector<int>> &sorted_adj,
        int root,
        list<int> &tour)
    {
        tour.clear();
        double total = 0.0;
        map<int, bool> visited;
        for (int x : candidates)
        {
            visited[x] = false;
        }
        int m = candidates.size();
        tour.push_back(root);
        visited[root] = true;
        while (tour.size() < m)
        {
            for (int j : sorted_adj.at(tour.back()))
            {
                if (!visited[j])
                {
                    total += dist_[tour.back()][j];
                    tour.push_back(j);
                    visited[j] = true;
                    break;
                }
            }
        }
        total += dist_[tour.back()][root];
        return total;
    }

    double TSP_greedy(
        const vector<vector<double>> &dist,
        const vector<int>& candidates,
        vector<int> &best_tour)
    {
        map<int,vector<int>> sorted_adj;
        for(int x:candidates)
        {
            sorted_adj[x] = {};
            for(int y:candidates)
            {
                if(x==y) continue;
                sorted_adj[x].push_back(y);
            }
            sort(sorted_adj[x].begin(), sorted_adj[x].end(), [&](int i, int j) -> bool
            {
                return dist[x][i] < dist[x][j];
            });
        }
        list<int> tour;
        Random rd;
        int root = rd.randint(candidates.size());
        root = candidates[root];


        double obj = oneway_greedy(dist, candidates, sorted_adj, root, tour);
        best_tour.assign(tour.begin(), tour.end());
        return obj;
    }

    double TSP_greedy(
        const vector<vector<double>> &dist,
        vector<int> &best_tour)
    {
        vector<int> candidates;
        for(int i=0;i<dist.size();i++)
        {
            candidates.push_back(i);
        }
        return TSP_greedy(dist,candidates,best_tour);
    }

// /////////////////////////////////////////

//     TSPAnnealingAgent::TSPAnnealingAgent(
//         DoubleTSPNodeList& sol,
//         const vector<vector<double>> &dist, 
//         const vector<vector<int>> &sorted_adj,
//         int root):
//         sol_(sol),dist_(dist),n_(dist.size()),root_(root),sorted_adj_(sorted_adj)
//     {
//         m_ = sol_.tour_len_within(root);
//     }

//     double TSPAnnealingAgent::process(int MAX_ITER)
//     {
//         double delta=0.0;
//         {
//             LKHAgent lkh_first(sol_, dist_, sorted_adj_);
//             delta+=lkh_first.lin_kernighan(root_);
//         }
//         DoubleTSPNodeList sol_cur(sol_);
//         double obj_best = evaluate(sol_cur);
//         double obj_cur = obj_best;
//         double cur_temp = 10.0;
//         int cur_iter_ = 0;
//         int no_best_update_ = 0, no_update_ = 0;

//         for(cur_iter_=0;cur_iter_<MAX_ITER;cur_iter_++)
//         {
//             LKHAgent lkh_cur(sol_cur,dist_,sorted_adj_);
//             double gain = _random_swap(lkh_cur);
//             gain += lkh_cur.lin_kernighan(root_);
//             if(obj_cur - gain < obj_best-1e-6)
//             {
//                 sol_ = sol_cur;
//                 delta += obj_best - obj_cur + gain;
//                 obj_best = obj_cur - gain;
//                 no_best_update_ = 0;
//             }
//             else 
//             {
//                 no_best_update_++;
//             }

//             if(gain>0 || rd_.bernoulli(exp(gain/cur_temp)))
//             {
//                 obj_cur-=gain;
//                 no_update_ = 0;
//             }
//             else 
//             {
//                 no_update_++;
//                 lkh_cur.recover();
//             }

//             if((double)no_best_update_>MAX_ITER/10.0)
//             {
//                 cur_temp*=0.9;
//                 no_best_update_ = 0;
//             }
//         }
//         return delta;
//     }

//     double TSPAnnealingAgent::evaluate(DoubleTSPNodeList& sol) const
//     {
//         int u = root_;
//         double val = 0.0;
//         do 
//         {
//             val += dist_[u][sol.next(u)];
//             u = sol.next(u);
//         } while (u!=root_);
//         return val;
//     }

//     double TSPAnnealingAgent::_random_swap(LKHAgent& agent)
//     {
//         int t1 = root_;

//         int n0 = rd_.randint(m_);
//         int n1 = rd_.randint(2,m_-5);
//         int n2 = rd_.randint(2,m_-n1-3);
        
//         for(int i=0;i<n0;i++)
//             t1 = agent.sol_.next(t1);
//         int t3 = t1;
//         for(int i=0;i<n1;i++)
//             t3 = agent.sol_.next(t3);
//         int t5 = t3;
//         for(int i=0;i<n2;i++)
//             t5 = agent.sol_.next(t5);
//         int t2 = agent.sol_.next(t1);
//         int t4 = agent.sol_.next(t3);
//         int t6 = agent.sol_.next(t5);

//         double val = 0.0;
//         val+=agent.FLIP_gain(t1,t2,t3,t4);
//         agent.FLIP(t1,t2,t3,t4);
//         val+=agent.FLIP_gain(t2,t4,t5,t6);
//         agent.FLIP(t2,t4,t5,t6);
//         val+=agent.FLIP_gain(t4,t6,t1,t3);
//         agent.FLIP(t4,t6,t1,t3);
//         return val;
//     }

// /////////////////////////////////////////

    TSP_Annealing::TSP_Annealing(const vector<vector<double>>& dist):distance_(dist)
    {
        INIT_TEMPE_=4.0;
        int n = dist.size();
        sorted_adj_.assign(n,vector<int>());
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
                if(i!=j)
                {
                    sorted_adj_[i].push_back(j);
                }
            }
            std::sort(sorted_adj_[i].begin(),sorted_adj_[i].end(),[&](int it1, int it2)->bool
            {
                return dist[i][it1] < dist[i][it2];
            });
        }
    }

    TSP_Annealing::~TSP_Annealing() 
    {
    }

    double TSP_Annealing::_evaluate(const DoubleTSPNodeList& sol) const 
    {
        size_t n = sol.size();
        if(n==0) return -1;
        double obj = 0.0;
        for(size_t i = 0; i < n; i++)
        {
            obj+= distance_[i][sol.next(i)];
        }
        return obj;
    }

    DoubleTSPNodeList TSP_Annealing::_init_solution()
    {
        if (sol_best_.size() == 0)
        {
            vector<int> tour;
            TSP_greedy(distance_, tour);
            sol_best_.load(tour);
            LKHAgent(sol_best_, distance_, sorted_adj_).lin_kernighan(0);
            obj_best_ = _evaluate(sol_best_);
        }
        return sol_best_;
    }

    DoubleTSPNodeList TSP_Annealing::_generate_new_solution(const DoubleTSPNodeList& sol)
    {
        size_t n = sol.size();
        DoubleTSPNodeList sol_new(sol);

        int n0 = rd_.randint(n);
        int n1 = rd_.randint(2,n-6);
        int n2 = rd_.randint(2,n-4-n1);
        int n3 = rd_.randint(2,n-2-n1-n2);
        int t1 = 0;
        for(int i=0;i<n0;i++) t1 = sol_new.next(t1);
        int t2 = sol_new.next(t1);
        int t3 = t1;
        for(int i=0;i<n1;i++) t3 = sol_new.next(t3);
        int t4 = sol_new.next(t3);
        int t5 = t3;
        for(int i=0;i<n2;i++) t5 = sol_new.next(t5);
        int t6 = sol_new.next(t5);
        int t7 = t5;
        for(int i=0;i<n3;i++) t7 = sol_new.next(t7);
        int t8 = sol_new.next(t7);

        sol_new.flip(t2,t3);
        sol_new.flip(t4,t5);
        sol_new.flip(t6,t7);
        sol_new.flip(t8,t1);

        _cross_update(sol_new);
        return sol_new;


        // vecInt selects = rd_.randints(n,3,false);
        // int t1,t2,t3,t4,t5,t6;
        // t1 = selects[0]; t3 = selects[1]; t5 = selects[2];
        // t2 = sol_new.next(t1); t4 = sol_new.next(t3); t6 = sol_new.next(t5);
        // while (t1 == t4 || t1 == t6 ||
        //        t3 == t2 || t3 == t6 ||
        //        t5 == t2 || t5 == t4)
        // {
        //     selects = rd_.randints(n,3,false);
        //     t1 = selects[0]; t3 = selects[1]; t5 = selects[2];
        //     t2 = sol_new.next(t1); t4 = sol_new.next(t3); t6 = sol_new.next(t5);
        // }

        // if(!sol_new.sequential(t1,t3,t5)) 
        // {
        //     swap(t3,t5);
        //     swap(t4,t6);
        // }
        // sol_new.flip(t2,t3);
        // sol_new.flip(t4,t5);
        // sol_new.flip(t6,t1);

        // _cross_update(sol_new);
        // return sol_new;
    }

    double TSP_Annealing::_accept_probability(double newone, double oldone) const
    {
        return 1/(1+std::exp(newone/oldone/tempe_cur_));
    }

    bool TSP_Annealing::_has_cross(int u1, int u2, int v1, int v2) const 
    {
        // u1 -> u2 -> .... -> v1 -> v2
        // try: u1 -> v1 -> ... -> u2 -> v2
        return distance_[u1][v1] + distance_[u2][v2] < distance_[u1][u2] + distance_[v1][v2]-1e-7;
    }

    void TSP_Annealing::_cross_update(DoubleTSPNodeList& sol)
    {
        LKHAgent lkh(sol,distance_,sorted_adj_,4);
        lkh.lin_kernighan(0);
    }

    void TSP_Annealing::solve()
    {
        double delta;
        DoubleTSPNodeList sol_cur = this->_init_solution();
        double obj_cur = this->_evaluate(sol_cur), obj_new;
        tempe_cur_ = INIT_TEMPE_;

        cur_iter_ = 0;
        while (cur_iter_ < MAX_ITER_)
        {
            // 产生新解
            DoubleTSPNodeList sol_new = this->_generate_new_solution(sol_cur);
            obj_new = this->_evaluate(sol_new);

            // 是否优于已发现最优解
            delta = obj_new - obj_best_;
            if (delta < 0)
            {
                sol_best_ = sol_new;
                obj_best_ = obj_new;
                no_best_update_ = 0;
                best_update_stats_++;
            }
            else
                no_best_update_++;

            // 是否更新当前解
            delta = obj_new - obj_cur;
            if (delta <= 0 || rd_.bernoulli(this->_accept_probability(obj_new, obj_cur)))
            {
                sol_cur = sol_new;
                obj_cur = obj_new;
                no_update_ = 0;
                update_stats_++;
            }
            else
                no_update_++;

            // 是否更新温度
            if (this->_should_temperature_update())
                tempe_cur_ = this->_temperature_update(tempe_cur_);

            cur_iter_++;
        }
    }

    void TSP_Annealing::print_stats()
    {
        std::cout << "MAX_ITER = " << MAX_ITER_;
        std::cout << ", update " << update_stats_ << " times, best solution update " << best_update_stats_ << " times.\n";
    }
} 
