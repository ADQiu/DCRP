#include"cyclecovering_grb.hh"
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

void findsubtour(const vector<vector<double>>& sol, vector<int>& tour);

class subtourelim: public GRBCallback
{
public:
    vector<vector<vector<GRBVar>>>& vars;
    int n;
    int m;
    vector<int>& lengths;
    subtourelim(vector<vector<vector<GRBVar>>>& xvars, int xn, int xm,vector<int>& xlengths)
    : vars(xvars), n(xn), m(xm), lengths(xlengths)
    {}
protected:
    void callback()
    {
        try
        {
            if( where == GRB_CB_MIPSOL)
            {
                for(int k=0;k<m;k++)
                {
                    vector<vector<double>> xvals(n,vector<double>(n,0));
                    vector<int> tour(n,0);
                    int i,j,len;
                    for(i=0;i<n;i++)
                    {
                        for(j=0;j<n;j++)
                            xvals[i][j] = getSolution(vars[k][i][j]);
                    }

                    findsubtour(xvals,tour);
                    len = tour.size();
                    
                    for(int k2=0;k2<m;k2++)
                    {
                        if(lengths[k2]<= len)
                            continue;
                        GRBLinExpr expr = 0;
                        for(i=0;i<len;i++)
                        {
                            for(j=0;j<len;j++)
                                expr+=vars[k2][tour[i]][tour[j]];
                        }
                        addLazy(expr <= len-1);
                    }
                }
            }
        }
        catch (GRBException e)
        {
            cout << "Error number: " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        }
        catch (...)
        {
            cout << "Error during callback" << endl;
        }
    }
};

CCPGRB::CCPGRB(const vector<vector<double>>& dist, const vector<int>& lengths):dist_(dist),lengths_(lengths),n(dist.size()),m(lengths.size())
{
    env_ = new GRBEnv(true);
    env_->start();
    env_->set(GRB_IntParam_OutputFlag, 1);
    env_->set(GRB_DoubleParam_TimeLimit, 600);
    env_->set(GRB_IntParam_LazyConstraints, 1);

    vars.reserve(n*n*m);
}

CCPGRB::~CCPGRB()
{
    if(env_ != nullptr)
        delete env_;
    if(model_!=nullptr)
        delete model_;
}

void CCPGRB::find_tours()
{
    best_sol.clear();
    best_objective = model_->get(GRB_DoubleAttr_ObjVal);
    for(int k=0;k<m;k++)
    {
        vector<vector<double>> ksol(n,vector<double>(n,0));
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
                ksol[i][j] = vars[k][i][j].get(GRB_DoubleAttr_X);
            }
        }
        vector<int> tour(n,0);
        int len;

        findsubtour(ksol, tour);
        len = tour.size();
        
        best_sol.push_back({});
        for(int i=0;i<len;i++)
        {
            best_sol.back().push_back(tour[i]);
        }
    }
}


void CCPGRB::create_model()
{
    model_ = new GRBModel(env_);
    for(int k=0;k<m;k++)
    {
        vars.push_back({});
        for (int i = 0; i < n; i++)
        {
            vars[k].push_back({});
            for (int j = 0; j < n; j++)
            {
                vars[k][i].push_back(model_->addVar(0, 1, dist_[i][j], GRB_BINARY));
            }
        }
    }

    // diagonal vanishing
    {
        GRBLinExpr expr = 0.0;
        for(int i=0;i<n;i++)
        {
            for(int k=0;k<m;k++)
                expr+=vars[k][i][i];
        }
        model_->addConstr(expr == 0);
    }

    // consistency
    for(int i=0;i<n;i++)
    {
        for(int k=0;k<m;k++)
        {
            GRBLinExpr expr = 0.0;
            for(int j=0;j<n;j++)
            {
                expr += vars[k][i][j];
                expr -= vars[k][j][i];
            }
            model_->addConstr(expr == 0);
        }
    }

    // size
    for(int k=0;k<m;k++)
    {
        GRBLinExpr expr = 0.0;
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
                expr += vars[k][i][j];
            }
        }
        model_->addConstr(expr == lengths_[k]);
    }

    // partition
    for(int i=0;i<n;i++)
    {
        GRBLinExpr expr = 0.0;
        for(int j=0;j<n;j++)
        {
            for(int k=0;k<m;k++)
                expr += vars[k][i][j];
        }
        model_->addConstr(expr == 1);
    }

    // objective
    {
        GRBLinExpr expr = 0.0;
        for(int k=0;k<m;k++)
        {
            for(int i=0;i<n;i++)
            {
                for(int j=0;j<n;j++)
                    expr += dist_[i][j]*vars[k][i][j];
            }
        }
        model_->setObjective(expr, GRB_MINIMIZE);
    }
}

void CCPGRB::solve()
{
    create_model();
    subtourelim cb(vars,n,m,lengths_);
    model_->setCallback(&cb);
    model_->optimize();
    find_tours();
}

void findsubtour(const vector<vector<double>>& sol, vector<int>& tour)
{
    int n = sol.size();
    vector<bool> seen(n,false);
    map<int,vector<int>> neighbor;
    int i,j, next, cur;

    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            if(sol[i][j]>0.5)
            {
                neighbor[i].push_back(j);
                neighbor[j].push_back(i);
            }
        }
    }

    tour.clear();
    cur = neighbor.begin()->first;
    while(true)
    {
        seen[cur]=true;
        tour.push_back(cur);
        next = neighbor[cur][0];
        if(seen[next] && neighbor[cur].size()>1)
            next = neighbor[cur][1];
        cur = next;
        if(seen[cur])
            break;
    }
}
