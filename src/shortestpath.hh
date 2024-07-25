#ifndef DCRP_SHORTEST_PATH_HH
#define DCRP_SHORTEST_PATH_HH
#include "simplegraph.hh"
#include <limits>

namespace dcrp
{
    class CostBend
    {
        public:
        double first;
        int second;
        double err=1e-6;
        CostBend(double cost=0.0, int bend=1) :first(cost),second(bend) {}
        CostBend(double cost, int bend, double e) :first(cost),second(bend),err(e) {}
        CostBend(const CostBend& r) :first(r.first),second(r.second),err(r.err) {}
        ~CostBend() {}

        bool operator<(const CostBend & right) const
        {
            if (first + err < right.first) return true;
            if (right.first + err < first) return false;
            return second < right.second;
        }
        bool operator==(const CostBend& r) const
        {
            return (first == r.first)&&(second == r.second);
        }

        CostBend& operator=(const CostBend& r)
        {
            first = r.first;
            second = r.second;
            err = r.err;
            return *this;
        }

        CostBend& operator+=(const CostBend& r)
        {
            first += r.first;
            second += r.second;
            return *this;
        }

        CostBend operator+(const CostBend& r) const
        {
            CostBend res(*this);
            res += r;
            return res;
        }
    };

    const CostBend MAX_COSTBEND(std::numeric_limits<double>::max(), std::numeric_limits<int>::max());

    class MinBendShortestPath
    {
    public:
        MinBendShortestPath(const SimpleGeometricGraph& g):g_(g) {}
        ~MinBendShortestPath() {}

        size_t get_predecessor_num(size_t v) const;
        size_t predecessor(size_t v, size_t n_pred) const;
        double distance(size_t v) const;
        int num_bend(size_t v) const;
        bool get_path(size_t v, std::vector<size_t>& pathvec, size_t n_direc=0) const;
        std::vector<size_t> get_path(size_t v, size_t n_direc=0) const;

        void solve(size_t root, size_t expected_end_node = std::numeric_limits<size_t>::max());

        std::vector<size_t> predecessors(size_t v) const;

    private:
        const SimpleGeometricGraph &g_;
        size_t root_;
        std::vector<std::vector<size_t>> predecessors_;
        std::vector<CostBend> cb_;
        void reset();
    };
}
#endif //DCRP_SHORTEST_PATH_HH