#include"capacitated_tree_grb.hh"
#include"cyclecovering.hh"
#include"tsp.hh"
#include<map>
#include<iterator>
#include<algorithm>
#include<thread>
#include<iostream>
#include <chrono>
using namespace std;
using namespace dcrp;


TreeGRB::TreeGRB(const vector<vector<double>>& dist, const vector<int>& lengths):dist_(dist),lengths_(lengths),n(dist.size()),m(lengths.size())
{
    env_ = new GRBEnv(true);
    env_->start();
    env_->set(GRB_IntParam_OutputFlag, 1);
    env_->set(GRB_DoubleParam_TimeLimit, 3600);
    // env_->set(GRB_IntParam_LazyConstraints, 1);

    x_.reserve(n*n);
    z_.reserve(n);
}

TreeGRB::~TreeGRB()
{
    if(env_ != nullptr)
        delete env_;
    if(model_!=nullptr)
        delete model_;
}

void TreeGRB::find_tours()
{
    best_sol.clear();
    best_objective = 0.0;
    
    for(int i=0;i<n;i++) 
    {
        double z = z_[i].get(GRB_DoubleAttr_X);
        if (z > 0.5)
        {
            vector<int> related;
            related.push_back(i);
            for(int j=0;j<n;j++)
            {
                if(x_[j][i].get(GRB_DoubleAttr_X) > 0.5)
                    related.push_back(j);
            }

            int nk = related.size();
            vector<vector<double>> subdist(nk, vector<double>(nk,0.0));
            for(int ik=0;ik<nk;ik++)
            {
                for(int jk=0;jk<nk;jk++)
                {
                    subdist[ik][jk] = dist_[related[ik]][related[jk]];
                }
            }
            TSP_Annealing tspsolver(subdist);
            tspsolver.solve();
            auto tspnodelist = tspsolver.get_best_solution();
            auto tourk = tspnodelist.dump_tours()[0];
            vector<int> tour;
            for (int ik:tourk)
                tour.push_back(related[ik]);
            best_sol.push_back(tour);

            for(int ik=0;ik<nk;ik++)
            {
                best_objective += dist_[tour[ik]][tour[(ik+1)%nk]];
            }
        }
    }
}


void TreeGRB::create_model()
{
    int length_min = lengths_[0];
    int length_max = lengths_[0];
    for(int i=0;i<m;i++)
    {
        length_min = min(length_min, lengths_[i]);
        length_max = max(length_max, lengths_[i]);
    }
    
    model_ = new GRBModel(env_);
    for (int i = 0; i < n; i++)
    {
        x_.push_back({});
        for (int j = 0; j < n; j++)
        {
            x_[i].push_back(model_->addVar(0, 1, dist_[i][j], GRB_BINARY));
        }
        z_.push_back(model_->addVar(0, 1, 0.0, GRB_BINARY));
    }

    // diagonal vanishing
    {
        GRBLinExpr expr = 0.0;
        for(int i=0;i<n;i++)
        {
            expr+=x_[i][i];
        }
        model_->addConstr(expr == 0);
    }

    // output 1 degree 
    {
        for (int i = 0; i < n; i++)
        {
            GRBLinExpr expr = 0.0;
            for (int j = 0; j < n; j++)
            {
                expr += x_[i][j];
            }
            expr += z_[i];
            model_->addConstr(expr == 1);
        }
    }
    {
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
                model_->addConstr(z_[i]>=x_[j][i]);
            }
        }
    }

    {
        for (int i = 0; i < n; i++)
        {
            GRBLinExpr expr = 0.0;
            expr += length_min * z_[i];
            for (int j = 0; j < n; j++)
            {
                expr -= x_[j][i];
            }
            
            model_->addConstr(expr <= 1);
        }
    }
    // {
    //     for (int i = 0; i < n; i++)
    //     {
    //         GRBLinExpr expr = 0.0;
    //         expr += length_max * z_[i];
    //         for (int j = 0; j < n; j++)
    //         {
    //             expr -= x_[j][i];
    //         }
            
    //         model_->addConstr(expr >= 1);
    //     }
    // }

    // partition size
    {
        GRBLinExpr expr = 0.0;
        for (int i = 0; i < n; i++)
        {
            expr += z_[i];
        }
        model_->addConstr(expr == m);
    }

    // objective
    {
        GRBLinExpr expr = 0.0;
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
                expr += dist_[i][j]*x_[i][j];
        }
        model_->setObjective(expr, GRB_MINIMIZE);
    }
}

void TreeGRB::solve()
{
    create_model();
    model_->optimize();
    find_tours();
}
