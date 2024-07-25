#include "simplegraph.hh"

namespace dcrp
{
    using namespace std;

    SimpleGraph::SimpleGraph(): num_vertex_(0), num_edge_(0) {}

    
    SimpleGraph::SimpleGraph(size_t num_vertex)
    {
        adjoint_.resize(num_vertex);
        num_vertex_ = num_vertex;
    }

    
    SimpleGraph::SimpleGraph(const SimpleGraph &other)
    {
        adjoint_.resize(other.num_vertex_);
        for (size_t i = 0; i < other.adjoint_.size(); i++)
            adjoint_[i] = other.adjoint_[i];
        num_vertex_ = other.num_vertex_;
    }

    
    size_t SimpleGraph::add_edge(size_t v1, size_t v2, double e)
    {
        if (v1 >= num_vertex() || v2 >= num_vertex()) set_vertex_num(((v1>v2)?v1:v2)+1);
        bool adding = adjoint_[v1].find(v2) == adjoint_[v1].end();
        if (v1!=v2)
        {
            this->adjoint_[v1][v2] = e;
            this->adjoint_[v2][v1] = e;
        }
        if(adding)
            num_edge_++;
        return num_edge_-1;
    }

    
    bool SimpleGraph::remove_edge(size_t v1, size_t v2)
    {
        if (this->adjoint_[v1].find(v2) != this->adjoint_[v1].end())
        {
            this->adjoint_[v1].erase(v2);
            this->adjoint_[v2].erase(v1);
            num_edge_--;
            return true;
        }
        return false;
    }

    map<size_t, double> SimpleGraph::reachable_neighbors(size_t i) const
    {
        if(i >= num_vertex_ || i<0) return std::map<size_t, double>();
        return adjoint_[i];
    }

    void SimpleGraph::set_vertex_num(size_t num_vertex)
    {
        if (num_vertex < num_vertex_)
        {
            for (size_t i = 0; i < num_vertex && i < adjoint_.size(); i++)
            {
                vector<size_t> tmp;
                for (auto it : this->adjoint_[i])
                {
                    if (it.first >= num_vertex)
                        tmp.push_back(it.first);
                }
                for (auto it : tmp)
                    this->adjoint_[i].erase(it);
            }
        }
        this->adjoint_.resize(num_vertex);
        this->num_vertex_ = num_vertex;
    }

    
    double SimpleGraph::get_weight(size_t v1, size_t v2) const
    {
        return adjoint_[v1].at(v2);
    }

    bool SimpleGraph::has_edge(size_t v1, size_t v2) const
    {
        return adjoint_[v1].find(v2) != adjoint_[v1].end();
    }
    
    void SimpleGraph::set_edge_weight(size_t v1, size_t v2, double e)
    {
        bool adding = adjoint_[v1].find(v2) == adjoint_[v1].end();
        if (v1!=v2)
        {
            this->adjoint_[v1][v2] = e;
            this->adjoint_[v2][v1] = e;
        }
        if(adding)
            num_edge_++;
    }

    
    SimpleGeometricGraph::SimpleGeometricGraph() {}

    
    SimpleGeometricGraph::SimpleGeometricGraph(const SimpleGeometricGraph &other)
        : SimpleGraph(other)
    {
        vertex_.clear();
        vertex_ = other.vertex_;
        REL_ERR_ = other.REL_ERR_;
        ABS_ERR_ = other.ABS_ERR_;
        WEAK_PARALLEL_ERR_ = other.WEAK_PARALLEL_ERR_;
    }

    
    size_t SimpleGeometricGraph::add_vertex(const Point &v)
    {
        size_t t = find_vertex(v);
        if (t >= this->num_vertex_)
        {
            vertex_.push_back(v);
            this->num_vertex_++;
            this->adjoint_.push_back(map<size_t, double>());
            return this->num_vertex_ - 1;
        }
        else
            return t;
    }

    
    bool SimpleGeometricGraph::edge_detect(const Point& p, size_t& v1, size_t& v2) const
    {
        for(size_t i = 0; i < this->num_vertex_; i++)
        {
            for(auto& e: this->adjoint_[i])
            {
                if(e.first <= i)
                    continue;
                double t;
                if(LineContainPoint(vertex_[i], vertex_[e.first]-vertex_[i],LineType::SEGMENT,p,t))
                {
                    v1 = i, v2 = e.first;
                    return true;
                }
            }
        }
        return false;
    }
    
    Point SimpleGeometricGraph::get_direction(size_t v1, size_t v2, bool normalize) const
    {
        Point dir = vertex_[v2] - vertex_[v1];
        if (normalize)
            return dir.normalized();
        return dir;
    }

    
    bool SimpleGeometricGraph::is_point_on_edge(size_t v1, size_t v2, const Point &p) const
    {
        Point dir = vertex_[v2] - vertex_[v1];
        Point dir_ = p - vertex_[v1];
        return dir.dot(dir_) < ABS_ERR_;
    }

    size_t SimpleGeometricGraph::find_vertex(const Point &p) const
    {
        for (size_t i = 0; i < vertex_.size(); i++)
        {
            if (vertex_[i].distance(p) < ABS_ERR_)
                return i;
        }
        return vertex_.size();
    }
} 
