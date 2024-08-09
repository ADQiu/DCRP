#include"dcrp.hh"
#include "gurobi_c++.h"
#include"capacitated_tree_grb.hh"
#include "string.h"

using namespace std;
using namespace dcrp;
class DCRPGRB: public DCRPAbstract
{
public:
    DCRPGRB(const string& filename):DCRPAbstract(filename) {}
    ~DCRPGRB() {}

    void ccp() override
    {
        TreeGRB solver(dists, lengths);
        solver.solve();
        sol = solver.get_solution();
        objective_value = solver.get_best_objective();
    }
};

int main(int argc, char** argv)
{
    bool stat_flag=false;

    if(argc<3)
    {
        cout<<"Bad number of parameters"<<endl;
        cout<<"Usage: "<<argv[0]<<" <inputfile> <outputfile> [-s <statistic file>] [-d <max depth>]"<<endl;
        return 1;
    }
    string filename(argv[1]);
    string outputfile(argv[2]);

    // string filename("toy1.dcrp");
    // string outputfile("toy1.solu");
    
    string statistic;
    int maxdep=3;

    for (int i = 3; i < argc; ++i)
    {
        if (strcmp(argv[i], "-s") == 0)
        {
            statistic = argv[i + 1];
            stat_flag = true;
        }
        else if (strcmp(argv[i], "-d") == 0)
        {
            maxdep = atoi(argv[i + 1]);
        }
    }

    if(stat_flag)
    {
        ifstream stat0(statistic, ios::in);
        bool status_stat = stat0.good();
        stat0.close();
        if (!status_stat)
        {
            ofstream stat1(statistic, ios::out);
            stat1 << "data,time,output,objective\n";
            stat1.close();
        }
    }

    ofstream stat(statistic,ios::app);
    DCRPGRB dcrp(filename);
    auto start = chrono::steady_clock::now();
    dcrp.solve();
    auto end = chrono::steady_clock::now();
    auto dura = chrono::duration_cast<chrono::milliseconds>(end-start).count()/1000.0;
    cout<< "solving " << filename <<" using time: "<<dura<<" s"<<endl;
    dcrp.write_result(outputfile);

    if(stat_flag)
    {
        stat << filename << ", " << dcrp.pure_ccp_time << "s, ";
        stat << outputfile << ", ";
        stat << dcrp.objective() << endl;
        stat.close();
    }

    return 0;
}