#ifndef DCRP_ARGHEAP_HPP
#define DCRP_ARGHEAP_HPP
#include<vector>
#include<functional>
#include<limits>
namespace dcrp
{
    template <typename T>
    class ArgHeap 
    {
    public:
        /**
         * @brief Construct a new Arg Heap object
         * 
         * @param values 
         * 注意：数据类型T的任意两个变量t1和t2，t1<t2 有定义
         */
        ArgHeap(const std::vector<T> &values);
        ~ArgHeap() 
        {
            delete[] vals_;
            delete[] indices_;
            delete[] arg_indices_;
        }
        size_t pop();
        void update(size_t i_val, const T& new_val);
        size_t size() const { return size_; }
        bool empty() const { return size_==0; }
        T get(size_t i) const { return vals_[i]; }

    private:
        T *vals_;
        size_t *indices_;
        size_t *arg_indices_;
        size_t size_, capacity_;
        void InitHeap();
        void Heapify(size_t i);
        void Swap(size_t i, size_t j);
    };

    template <typename T>
    ArgHeap<T>::ArgHeap(const std::vector<T> &values)
    {
        size_ = values.size();
        capacity_ = values.size();
        vals_ = new T[capacity_];
        indices_ = new size_t[capacity_];
        arg_indices_ = new size_t[capacity_];
        
        for(size_t i=0;i<capacity_;++i)
        {
            vals_[i] = values[i];
            indices_[i] = i;
            arg_indices_[i] = i;
        }

        
        InitHeap();
    }

    template <typename T>
    void ArgHeap<T>::InitHeap()
    {
        if(size_==0) return;
        for (size_t i = (size_-1) / 2+1; i > 0; --i)
        {
            Heapify(i-1);
        }
    }

    template <typename T>
    void ArgHeap<T>::Heapify(size_t i)
    {
        size_t l = 2 * i + 1;
        size_t r = 2 * i + 2;
        size_t smallest = i;
        if (l < size_ && (vals_[indices_[l]]< vals_[indices_[i]]))
        {
            smallest = l;
        }
        if (r < size_ && (vals_[indices_[r]]< vals_[indices_[smallest]]))
        {
            smallest = r;
        }
        if (smallest != i)
        {
            Swap(i, smallest);
            Heapify(smallest);
        }
    }

    template <typename T>
    void ArgHeap<T>::Swap(size_t i, size_t j)
    {
        std::swap(indices_[i], indices_[j]);
        arg_indices_[indices_[i]] = i;
        arg_indices_[indices_[j]] = j;
    }

    template <typename T>
    size_t ArgHeap<T>::pop()
    {
        if(size_ == 0) return std::numeric_limits<size_t>::max();
        
        size_t ret = indices_[0];
        indices_[0] = indices_[size_-1];
        indices_[size_-1] = ret;
        --size_;
        arg_indices_[indices_[0]] = 0;
        arg_indices_[ret] = size_;
        Heapify(0);
        return ret;
    }

    template <typename T>
    void ArgHeap<T>::update(size_t i_val, const T& new_val)
    {
        vals_[i_val] = new_val;
        size_t k = arg_indices_[i_val];
        if(k>= size_) return;
        while(k>0&&(vals_[indices_[k]] < vals_[indices_[(k-1)/2]]))
        {
            Swap(k, (k-1)/2);
            k = (k-1)/2;
        }
        Heapify(k);
    }
}
#endif