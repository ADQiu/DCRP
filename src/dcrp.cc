#include "dcrp.hh"

using namespace dcrp;

static bool covered_by_one(const vector<Rectangle<int>> &tuples, int x, int y)
{
    for (const auto &t : tuples)
    {
        if (t.cover(x, y))
            return true;
    }
    return false;
}

static list<double>::const_iterator find_ascend_pos(const list<double> &l, double x, double err = 1)
{
    list<double>::const_iterator it = l.begin();
    while (it != l.end() && *it + err < x)
        it++;
    return it;
}

static void ascend_insert(list<double> &l, double x, double err = 1)
{
    list<double>::const_iterator it = find_ascend_pos(l, x, err);
    if (it == l.end() || x + err < *it)
        l.insert(it, x);
}

static int find_pos(const vector<double> &arr, double val, int tol = 1)
{
    int i = 0;
    while (i < arr.size() && arr[i] + tol < val)
        i++;
    return i;
}

DCRPAbstract::DCRPAbstract(const string &filename)
{
    ifstream in(filename);
    int n_lengths;
    in >> n_module >> n_barrier >> n_lengths;

    for (int i = 0; i < n_lengths; i++)
    {
        int num;
        in >> num;
        lengths.push_back(num);
    }
    for (int i = 0; i < n_module; i++)
    {
        double x, y;
        in >> x >> y;
        modules.push_back(Point(x, y, 0.0));
    }
    for (int i = 0; i < n_barrier; i++)
    {
        double x0, x1, y0, y1;
        in >> x0 >> x1 >> y0 >> y1;
        barriers.push_back(Rectangle<double>(x0, x1, y0, y1));
    }
    in.close();
}

void DCRPAbstract::extract_xy(list<double> &xlist, list<double> &ylist, double tol)
{

    for (auto &p : modules)
    {
        ascend_insert(xlist, p.x(), tol);
        ascend_insert(ylist, p.y(), tol);
    }
    for (auto &b : barriers)
    {
        ascend_insert(xlist, b.x0, tol);
        ascend_insert(xlist, b.x1, tol);
        ascend_insert(ylist, b.y0, tol);
        ascend_insert(ylist, b.y1, tol);
    }
}

void DCRPAbstract::build_graph()
{
    mod2ver.clear();
    double tol = 1;
    list<double> xlist, ylist;
    extract_xy(xlist, ylist, tol);

    vector<double> xs(xlist.begin(), xlist.end()), ys(ylist.begin(), ylist.end());

    // find range of barriers
    vector<Rectangle<int>> barange;
    for (auto &b : barriers)
    {
        barange.push_back(Rectangle<int>(find_pos(xs, b.x0, tol),
                                    find_pos(xs, b.x1, tol),
                                    find_pos(ys, b.y0, tol),
                                    find_pos(ys, b.y1, tol)));
    }

    int xsn = xs.size(), ysn = ys.size();
    for (int i = 0; i < xsn; i++)
    {
        for (int j = 0; j < ysn; j++)
        {
            graph.add_vertex(Point(xs[i], ys[j], 0.0));
            int n_cur = i * ysn + j;
            if (covered_by_one(barange, i, j))
                continue;
            if (i > 0)
            {
                int n_prev = i * ysn + j - ysn;
                if (!covered_by_one(barange, i - 1, j))
                {
                    graph.add_edge(n_cur, n_prev, xs[i] - xs[i - 1]);
                }
            }
            if (j > 0)
            {
                int n_prev = i * ysn + j - 1;
                if (!covered_by_one(barange, i, j - 1))
                {
                    graph.add_edge(n_cur, n_prev, ys[j] - ys[j - 1]);
                }
            }
        }
    }

    // find xnum and ynum for each module
    vector<pair<int, int>> modnums;
    for (auto &p : modules)
    {
        modnums.push_back(make_pair(find_pos(xs, p.x(), tol),
                                    find_pos(ys, p.y(), tol)));
    }

    for (auto &pa : modnums)
    {
        mod2ver.push_back(pa.first * ysn + pa.second);
        ver2mod[pa.first * ysn + pa.second] = (int)(mod2ver.size() - 1);
    }
}

void DCRPAbstract::find_paths()
{
    dists.assign(n_module, vector<double>(n_module, numeric_limits<double>::max()));
    MinBendShortestPath solver(graph);
    for (int i = 0; i < n_module - 1; i++)
    {
        solver.solve(mod2ver[i]);
        for (int j = i + 1; j < n_module; j++)
            dists[i][j] = dists[j][i] = solver.distance(mod2ver[j]);
    }
}

void DCRPAbstract::write_result(const string &filename)
{
    ofstream fout(filename);
    fout << lengths.size() << endl;
    for (auto &s : sol)
        fout << s.size() << " ";
    fout << endl;
    for (auto &s : sol)
    {
        for (auto &i : s)
        {
            auto &pnt = modules[i];
            fout << pnt.x() << " " << pnt.y() << ";";
        }
        fout << endl;
    }
    fout.close();
}

void DCRPAbstract::solve()
{
    build_graph();
    find_paths();
    auto start = chrono::steady_clock::now();
    ccp();
    auto end = chrono::steady_clock::now();
    pure_ccp_time = chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000.0;
}