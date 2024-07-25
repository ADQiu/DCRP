#ifndef DCRP_RANDOM_HH
#define DCRP_RANDOM_HH
#include<random>
#include<vector>
#include<algorithm>
namespace dcrp
{
    class Random
    {
    public:
        Random();
        Random(int seed);
        ~Random();

        /**
         * @brief 从[0,1]区间上的均匀分布采样
         * 
         * @return double 
         */
        double random();

        /**
         * @brief 从[a,b]区间上的均匀分布采样
         * 
         * @param a 区间左端点
         * @param b 区间右端点
         * @return double 
         */
        double random(double a, double b);

        /**
         * @brief 从[0,1]区间上的均匀分布多次采样
         * 
         * @param num 采样次数
         * @return std::vector<double> 
         */
        std::vector<double> random(int num);

        /**
         * @brief 从[a,b]区间上的均匀分布多次采样
         * 
         * @param a 区间左端点
         * @param b 区间右端点
         * @param num 采样次数
         * @return std::vector<double> 
         */
        std::vector<double> random(double a, double b, int num);

        /**
         * @brief 从[0,b)的区间上采样整数
         * 
         * @param b 区间右端点
         * @return int 
         */
        int randint(int b = 2);

        /**
         * @brief 从[a,b)的区间上采样整数
         * 
         * @param a 区间左端点
         * @param b 区间右端点
         * @return int 
         */
        int randint(int a, int b);

        /**
         * @brief 从[0,b)的区间上以概率向量p采样整数
         * 
         * @param b 区间右端点
         * @param p 概率向量
         * @return int 
         */
        int randint(int b, const std::vector<double>& p);

        /**
         * @brief 从[0,b)的区间上多次采样整数
         * 
         * @param b 区间右端点
         * @param num 采样次数
         * @param replace 是否允许多次采样同一个整数
         * @return std::vector<int> 
         */
        std::vector<int> randints(int b, int num, bool replace = true);

        /**
         * @brief 从[0,b)的区间上以概率向量p多次采样整数
         * 
         * @param b 区间右端点
         * @param num 采样次数
         * @param p 概率向量
         * @param replace 是否允许多次采样同一个整数
         * @return std::vector<int> 
         */
        std::vector<int> randints(int b, int num, const std::vector<double>& p, bool replace = true);

        /**
         * @brief 从概率为p的伯努利分布采样
         * 
         * @param p 概率参数
         * @return true 
         * @return false 
         */
        bool bernoulli(double p);

        /**
         * @brief 从向量arr中等概率地选择一个元素
         * 
         * @param arr 备选元素构成的向量
         * @return T 
         */
        template <typename T>
        T choice(const std::vector<T> &arr)
        {
            int num = arr.size();
            return arr[randint(num)];
        }

        /**
         * @brief 从向量arr中以概率p选择一个元素
         * 
         * @param arr  备选元素构成的向量
         * @param p 概率向量
         * @return T 
         */
        template <typename T>
        T choice(const std::vector<T> &arr, const std::vector<double> &p)
        {
            int i = randint(arr.size(), p);
            return arr[i];
        }

        /**
         * @brief 从向量arr中等概率地选择num个元素
         * 
         * @param arr 备选元素构成的向量
         * @param num 采样次数
         * @param replace 是否允许多次采样同一个元素
         * @return std::vector<T> 
         */
        template <typename T>
        std::vector<T> choice(const std::vector<T> &arr, int num, bool replace = true)
        {
            std::vector<int> ind = randints(arr.size(), num, replace);
            std::vector<T> out(ind.size());
            for(int i = 0; i < ind.size(); ++i)
                out[i] = arr[ind[i]];
            return out;
        }

        /**
         * @brief 从向量arr中以概率p选择num个元素
         * 
         * @param arr 备选元素构成的向量
         * @param num 采样次数
         * @param p 概率向量
         * @param replace 是否允许多次采样同一个元素
         * @return std::vector<T> 
         */
        template <typename T>
        std::vector<T> choice(const std::vector<T> &arr, int num, const std::vector<double> &p, bool replace = true)
        {
            std::vector<int> ind = randints(arr.size(), num, p, replace);
            std::vector<T> out(ind.size());
            for(int i = 0; i < ind.size(); ++i)
                out[i] = arr[ind[i]];
            return out;
        }

        template <typename T>
        void shuffle(std::vector<T> &arr)
        {
            std::shuffle(arr.begin(), arr.end(), mtgen_);
        }

    protected:
        std::mt19937 mtgen_;

        template<typename T>
        inline T MAX(T a, T b){ return (a)<(b)?(b):(a); }
    };
    
} 


#endif //DCRP_RANDOM_HH