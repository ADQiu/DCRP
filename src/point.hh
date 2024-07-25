#ifndef DCRP_POINT_HH
#define DCRP_POINT_HH
#include <iostream>
#include <string>
#include <vector>
#include <tuple>

namespace dcrp
{
	/**
	 * @brief 三维点类，同时也是向量
	 * 
	 */
	class Point
	{
	private:
		double x_;
		double y_;
		double z_;

	public:
		Point(double x_val=0.0, double y_val=0.0, double z_val=0.0);
		Point(const Point &p);
		~Point();
		inline double x() const { return x_; }
		inline double y() const { return y_; }
		inline double z() const { return z_; }
		
		/**
		 * @brief 计算向量的欧氏范数
		 * 
		 * @return double 
		 */
		double norm() const;
		
		/**
		 * @brief 计算两点之间的欧氏距离
		 * 
		 * @param v 另一点
		 * @return double 
		 */
		double distance(const Point &v) const;

		/**
		 * @brief 计算该向量与另一向量v的叉乘（向量积）
		 * 
		 * @param v 
		 * @return Point 
		 */
		Point cross(const Point &v) const;

        /**
		 * @brief 计算该向量与另一向量v的内积
		 * 
		 * @param v 
		 * @return value 
		 */
        double dot(const Point &v) const;

		/**
		 * @brief 以z方向为轴将向量【顺时针】旋转angle角度
		 * 
		 * @param angle 角度（注意不是弧度）
		 * @return Point 
		 */
		Point rotate(const int &angle) const;

		/**
		 * @brief 计算当前向量的单位向量
		 * 
		 * @return Point 
		 */
		Point normalized() const;

		Point operator+(const Point &v) const;
		Point operator-(const Point &v) const;
		Point operator*(const double &t) const;
		friend Point operator*(double t, const Point& p) {return p*t;}
		double operator*(const Point &v) const;
		Point operator/(const double &t) const;
		double operator/(const Point &v) const;
		double operator[](const int &i) const;
		Point& operator=(const Point &v);
		bool operator==(const Point& v) const;
		bool operator!=(const Point& v) const;

		/**
		 * @brief 判断当前向量与另一向量是否平行
		 * 
		 * @param v 向量
		 * @param REL_ERR 精度
		 * @return true 
		 * @return false 
		 */
		bool IsParallel(const Point &v, double REL_ERR = 0.01) const;

		/**
		 * @brief 判断当前向量与另一向量是否同向
		 * 
		 * @param v 向量
		 * @param REL_ERR 精度
		 * @return true 
		 * @return false 
		 */
		bool IsSameDirection(const Point &v, double REL_ERR = 0.01) const;

		/**
		 * @brief 判断当前向量与另一向量是否反向
		 * 
		 * @param v 向量
		 * @param REL_ERR 精度 
		 * @return true 
		 * @return false 
		 */
		bool IsReverseDirection(const Point &v, double REL_ERR = 0.01) const;
		
		/**
		 * @brief 判断当前向量与另一向量是否弱平行
		 * 弱平行即 差不多平行 ，比如(0,0,1)与(0,0,1.1)弱平行
		 * @param v 另一向量
		 * @param REL_ERR 精度
		 * @param WEAK_PARALLEL_ERR 弱平行误差界（需要小于1/sqrt(2)，约0.7）
		 * @return true 
		 * @return false 
		 */
		bool IsWeakParallel(const Point &v, double REL_ERR = 0.01,double WEAK_PARALLEL_ERR = 0.3) const;
		
		/**
		 * @brief 判断当前向量与另一向量是否弱同向
		 * 弱同向即 差不多同向 ，比如(0,0,1)与(0,0,1.1)弱同向
		 * @param v 另一向量
		 * @param REL_ERR 精度
		 * @param WEAK_PARALLEL_ERR 弱同向误差界（需要小于1/sqrt(2)，约0.7）
		 * @return true 
		 * @return false 
		 */
		bool IsWeakSameDirection(const Point &v, double REL_ERR = 0.01,double WEAK_PARALLEL_ERR = 0.3) const;
		
		/**
		 * @brief 判断当前向量与另一向量是否弱反向
		 * 弱反向即 差不多反向 ，比如(0,0,1)与(0,0,-1.1)弱反向
		 * @param v 另一向量
		 * @param REL_ERR 精度
		 * @param WEAK_PARALLEL_ERR 弱反向误差界（需要小于1/sqrt(2)，约0.7）
		 * @return true 
		 * @return false 
		 */
		bool IsWeakReverseDirection(const Point &v, double REL_ERR = 0.01,double WEAK_PARALLEL_ERR = 0.3) const;

		friend std::ostream &operator<<(std::ostream &out, const Point &p)
		{
			out << "[" << p.x_ << "," << p.y_ << "," << p.z_ << "]";
			return out;
		}
		
		/**
		 * @brief 打印当前向量
		 * 
		 * @return std::string 
		 */
		std::string ToStr() const
		{
			std::string ss;
			ss = "[" + std::to_string(x_) + "," + std::to_string(y_) + "," + std::to_string(z_) + "]";
			return ss;
		}
		double l1distance(const Point &v) const;
	};
	

	/**
	 * @brief 直线种类
	 * 
	 */
	enum class LineType
	{
		SEGMENT,
		RAY,
		LINE
	};

	bool LineContainPoint(const Point &start, const Point & direct, LineType type,const Point &point, double& t);
}
#endif // DCRP_POINT_HH