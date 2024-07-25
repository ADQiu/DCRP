#ifndef DCRP_SIMPLEGRAPH_HH
#define DCRP_SIMPLEGRAPH_HH
#include "point.hh"
#include <vector>
#include <tuple>
#include <map>
#include <functional>
////////////////////////////////////////////
// 简单图的最大特点是不需要一个线性表维护边集
// 只需要根据邻接关系即可确定边集
//
// 不使用线性表维护边集的最大好处是边集的删除操作较快
////////////////////////////////////////////
namespace dcrp
{
    using std::map;
    using std::vector;

    class SimpleGraph
    {
    protected:
        vector<map<size_t, double>> adjoint_;
        size_t num_vertex_ = 0;
        size_t num_edge_ = 0;

    public:
        SimpleGraph();
        SimpleGraph(size_t num_vertex);
        SimpleGraph(const SimpleGraph &other);
        ~SimpleGraph() {}

        void set_vertex_num(size_t num_vertex);
        bool has_edge(size_t v1, size_t v2) const;
        std::map<size_t, double> reachable_neighbors(size_t i) const; 
        double get_weight(size_t v1, size_t v2) const;

        const vector<map<size_t, double>> &adjoint() const { return adjoint_; }
        const map<size_t, double> &adjoint(size_t v) const { return adjoint_[v]; }

        size_t add_edge(size_t v1, size_t v2, double weight = 1.0);
        bool remove_edge(size_t v1, size_t v2);
        void set_edge_weight(size_t v1, size_t v2, double weight);

        size_t num_vertex() const 
        { 
            return num_vertex_; 
        }

        size_t num_edge() const 
        { 
            return num_edge_; 
        }

        virtual std::map<size_t, double> neighbors(size_t v) const {return reachable_neighbors(v);}
    };

    class SimpleGeometricGraph : public SimpleGraph
    {
    protected:
        std::vector<Point> vertex_;
        double REL_ERR_ = 0.001;
        double ABS_ERR_ = 0.001;
        double WEAK_PARALLEL_ERR_ = 0.3;

    public:
        SimpleGeometricGraph();
        SimpleGeometricGraph(const SimpleGeometricGraph &other);
        ~SimpleGeometricGraph() {}

        size_t add_vertex(const Point &p);
        bool edge_detect(const Point& p,size_t& v1, size_t& v2) const;

        Point get_direction(size_t v1, size_t v2, bool normalize = false) const;
        bool is_point_on_edge(size_t v1, size_t v2, const Point &p) const;

        /**
         * @brief 相对误差
         * 
         * @return double 
         */
        double REL_ERR() const { return REL_ERR_; }

        /**
         * @brief 绝对误差
         * 
         * @return double 
         */
        double ABS_ERR() const { return ABS_ERR_; }

        /**
         * @brief 弱平行误差
         * 
         * @return double 
         */
        double WEAK_PARALLEL_ERR() const { return WEAK_PARALLEL_ERR_; }

        /**
         * @brief 设置相对误差
         * 
         * @param err 
         */
        void set_REL_ERR(double err) { REL_ERR_ = err; }

        /**
         * @brief 设置绝对误差
         * 
         * @param err 
         */
        void set_ABS_ERR(double err) { ABS_ERR_ = err; }

        /**
         * @brief 设置弱平行误差
         * 
         * @param err 
         */
        void set_WEAK_PARALLEL_ERR(double err) { WEAK_PARALLEL_ERR_ = err; }

        /**
         * @brief 寻找坐标点p的结点序号
         * 若结点集中不包含该坐标点，则返回结点集大小
         * @param p 坐标点
         * @return size_t 
         */
        size_t find_vertex(const Point &p) const;
        
        Point vertex(size_t v) const 
        { 
            return vertex_[v]; 
        }
    };
} 
#endif // DCRP_SIMPLEGRAPH_HH