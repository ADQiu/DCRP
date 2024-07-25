#include "shortestpath.hh"
#include "argheap.hpp"
#include <queue>
#include <tuple>
#include <algorithm>
#include <limits>

using namespace std;

namespace dcrp
{
    void MinBendShortestPath::reset()
    {
        root_=g_.num_vertex();
        predecessors_.assign(g_.num_vertex(),vector<size_t>());
        cb_.assign(g_.num_vertex(),CostBend(numeric_limits<double>::infinity(),numeric_limits<int>::infinity()));
    }

    size_t MinBendShortestPath::get_predecessor_num(size_t v) const
    {
        if(v>=predecessors_.size()) return 0;
        return predecessors_[v].size();
    }

    size_t MinBendShortestPath::predecessor(size_t v, size_t n_direc) const
    {
        if(v>=predecessors_.size()) return numeric_limits<size_t>::infinity();
        if(n_direc >= predecessors_[v].size()) return numeric_limits<size_t>::infinity();
        return predecessors_[v][n_direc];
    }

    double MinBendShortestPath::distance(size_t v) const
    {
        if(v>=cb_.size()) return numeric_limits<double>::infinity();
        return cb_[v].first;
    }

    int MinBendShortestPath::num_bend(size_t v) const 
    {
        if(v>=cb_.size()) return numeric_limits<int>::infinity();
        return cb_[v].second;
    }

    bool MinBendShortestPath::get_path(size_t v, vector<size_t>& pathvec, size_t n_direc) const
    {
        pathvec.clear();
        if (predecessors_.size() <= v) return false;
        if (predecessors_[v].size() == 0) return false;

        if(n_direc >= predecessors_[v].size()) n_direc = 0;

        double REL_ERR = g_.REL_ERR();
        double WEAK_PARA_ERR = g_.WEAK_PARALLEL_ERR();
        size_t q = predecessors_[v][n_direc];
        pathvec.push_back(v);
        pathvec.push_back(q);
        size_t u = q;
        Point d = (g_.vertex(v)-g_.vertex(q)).normalized();
        while (predecessors_[u].size() > 0)
        {
            size_t itr = 0; // 同向前继的序号
            Point up = g_.vertex(u);
            for(size_t itr2 = 1; itr2 < predecessors_[u].size();itr2++)
            {
                Point q1 = g_.vertex(predecessors_[u][itr]);
                Point q2 = g_.vertex(predecessors_[u][itr2]);
                if((!d.IsWeakParallel(up-q1,REL_ERR,WEAK_PARA_ERR))&&
                    d.IsWeakParallel(up-q2,REL_ERR,WEAK_PARA_ERR))
                    itr = itr2;
            }
            u = predecessors_[u][itr];
            d = (up - g_.vertex(u)).normalized();
            pathvec.push_back(u);
        }
        reverse(pathvec.begin(), pathvec.end());
        return true;
    }

    vector<size_t> MinBendShortestPath::get_path(size_t v, size_t n_direc) const
    {
        vector<size_t> out;
        get_path(v,out,n_direc);
        return out;
    }

    void MinBendShortestPath::solve(size_t root, size_t expected_end)
    {
        reset();
        root_ = root;
        size_t n = g_.num_vertex();
        vector<bool> visited(n,false);
        double REL_ERR = g_.REL_ERR();
        double ABS_ERR = g_.REL_ERR();
        double WEAK_PARA_ERR = g_.WEAK_PARALLEL_ERR();
        using CB = CostBend;
        const CB MAX_CB(numeric_limits<double>::infinity(), numeric_limits<int>::max(), ABS_ERR);
        for(auto& cb:cb_) cb.err = ABS_ERR;

        ArgHeap<CB> h(vector<CB>(n, MAX_CB));
        size_t v;
        Point d;

        h.update(root, CB(0.0,-1,ABS_ERR));

        while (h.size() > 0)
        {
            v = h.pop();
            visited[v] = true;
            cb_[v] = h.get(v);
            if (v == expected_end)
                break;
            if (h.get(v) == MAX_CB)
                break;
            if ( (expected_end < n) && h.get(expected_end)<h.get(v))
                continue;
            for(auto nb : g_.reachable_neighbors(v))
            {
                size_t u = nb.first;
                double w = nb.second;
                if (visited[u]) continue;
                Point d = g_.vertex(u) - g_.vertex(v);
                CB pont(cb_[v].first+w,cb_[v].second+1,ABS_ERR);
                for(auto x: predecessors_[v])
                {
                    Point direc0 = g_.vertex(v)-g_.vertex(x); 
                    if(d.IsWeakParallel(direc0,REL_ERR,WEAK_PARA_ERR))
                    {
                        pont.second -= 1;
                        break;
                    }
                }
                if(pont < h.get(u))
                {
                    h.update(u,pont);
                    predecessors_[u] = {v};
                }
                else if(!(h.get(u)<pont))
                {
                    predecessors_[u].push_back(v);
                }
            }
        }
    }

    vector<size_t> MinBendShortestPath::predecessors(size_t v) const 
    {
        if(v>=predecessors_.size()) return {};
        return predecessors_[v];
    }
}