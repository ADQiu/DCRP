
/*****************************
 * 输入数据格式
 ***************************** 
 * 连线对象数量
 * 矩形障碍物数量
 * 分组数量
 * 分组大小1 分组大小2 分组大小3 ...
 * 连线对象1.x 连线对象1.y
 * 连线对象2.x 连线对象2.y
 * ...
 * 障碍物1.xmin 障碍物1.xmax 障碍物1.ymin 障碍物1.ymax
 * 障碍物2.xmin 障碍物2.xmax 障碍物2.ymin 障碍物2.ymax
 * ...
 * 
 ******************************
 * 输出数据格式
 *****************************  
 * 分组数量
 * 分组大小1 分组大小2 分组大小3 ...
 * 第1组对象1.x 第1组对象1.y;第1组对象2.x 第1组对象2.y;...
 * 第2组对象1.x 第2组对象1.y;第2组对象2.x 第2组对象2.y;...、
 * ...
 ******************************/

#pragma once
#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include "cyclecovering.hh"
#include "point.hh"
#include "shortestpath.hh"
#include <chrono>
#include <list>

using namespace std;

namespace dcrp
{
    template<typename T>
    struct Rectangle
    {
        T x0, x1, y0, y1;
        Rectangle(T x0, T x1, T y0, T y1) : x0(x0), x1(x1), y0(y0), y1(y1) {}
        bool cover(T x, T y) const { return (x0 < x && x < x1 && y0 < y && y < y1); }
    };

    class DCRPAbstract
    {
    public:
        DCRPAbstract(const string &filename);
        ~DCRPAbstract() {}

        void build_graph();
        void find_paths();
        void solve();
        void write_result(const string &filename);

        void virtual ccp() = 0;
        double objective() { return objective_value; }

        int n_module, n_barrier;
        vector<Point> modules;
        vector<Rectangle<double>> barriers;
        vector<int> lengths;
        vector<vector<double>> dists;
        vector<vector<int>> sol;
        double objective_value;
        double pure_ccp_time = 0.0;

        SimpleGeometricGraph graph;
        vector<int> mod2ver;
        map<int, int> ver2mod;

    private:
        void extract_xy(list<double> &xlist, list<double> &ylist, double tol = 1);
    };

}