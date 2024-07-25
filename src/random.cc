#include "random.hh"
#include <stdexcept>
#define ASSERT_POSITIVE(x) {if(x<=0) throw std::runtime_error("Random: non-positive value");}
#define ASSERT_NONNEGATIVE(x) {if(x<0) throw std::runtime_error("Random: negative value");}
namespace dcrp
{
    using std::vector;
    Random::Random()
    {
        std::random_device rd;
        mtgen_ = std::mt19937(rd());
    }

    Random::Random(int seed)
    {
        mtgen_ = std::mt19937(seed);
    }

    Random::~Random() {}

    double Random::random()
    {
        return (double) mtgen_() / (double)std::mt19937::max();
    }

    double Random::random(double a, double b)
    {
        return a + (double) mtgen_() / (double)std::mt19937::max()*(b-a);
    }

    std::vector<double> Random::random(int num)
    {
        vector<double> out;
        for(int i = 0; i < num;i++) out.push_back(random());
        return out;
    }

    std::vector<double> Random::random(double a, double b, int num)
    {
        vector<double> out;
        for(int i = 0; i < num;i++) out.push_back(random(a,b));
        return out;
    }

    int Random::randint(int b)
    {
        ASSERT_POSITIVE(b);
        return mtgen_() % b;
    }

    int Random::randint(int a, int b)
    {
        if(a >= b-1)
            return a;
        return a + randint(b - a);
    }

    int Random::randint(int b, const std::vector<double>& p)
    {
        ASSERT_POSITIVE(b);
        double psum=0.0;
        for(double p0:p)
            psum += MAX(p0,0.0);

        double pr = random();
        double s=0.0;
        for(int i = 0; i < b ; i++)
        {
            s += ((p[i]>0)?p[i]:0)/(psum + 1e-12);
            if(s >= pr)
                return i;
        }
        return b-1;
    }

    vector<int> Random::randints(int b, int num, bool replace)
    {
        ASSERT_POSITIVE(b);
        ASSERT_NONNEGATIVE(num);
        if(!replace)
            ASSERT_NONNEGATIVE(b-num);

        vector<int> out(num);
        int n = 0;
        int* source = new int[b];
        for (int i = 0; i < b; i++)
            source[i] = i;

        int tail = b;
        int j;
        while (n < num)
        {
            j = randint(tail);
            out[n++] = source[j];
            if (!replace)
            {
                std::swap(source[tail-1],source[j]);
                tail--;
            }
            if(tail <= 0)
                break;
        }  
        delete[] source;
        return out;
    }

    vector<int> Random::randints(int b, int num, const std::vector<double>& p, bool replace)
    {
        ASSERT_POSITIVE(b);
        ASSERT_NONNEGATIVE(num);
        if(!replace)
            ASSERT_NONNEGATIVE(num-b);

        vector<int> out(num);
        int n=0;
        int* source = new int[b];
        double denom=0.0;
        for(double p0: p) denom += MAX(p0,0.0);
        for(int i=0;i<b;i++) source[i] = i;

        int tail = b;
        int j;
        double spl,cursum;
        int q;
        while (n < num)
        {
            spl = random();
            cursum = 0.0;
            for(j=0;j<tail;j++)
            {
                q = source[j];
                cursum += MAX(p[q],0.0)/(denom+1e-12);
                if (spl <= cursum)
                    break;
            }
            if(j == tail)
                j--;
            q = source[j];
            out[n++] = q; 
            if(!replace)
            {
                std::swap(source[j],source[tail-1]);
                tail--;
                denom -= MAX(p[q],0.0);
                if(tail <= 0)
                    break;
            }
        }
        delete[] source;
        return out;
    }

    bool Random::bernoulli(double p)
    {
        return random() < p;
    }

} 
