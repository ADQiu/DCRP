#include "point.hh"
#include <cmath>

constexpr auto ERR = 1.0e-2;
template <typename L, typename R>
auto Min(const L &l, const R &r) -> decltype(l < r ? l : r) { return l < r ? l : r; }
template <typename L, typename R>
auto Max(const L &l, const R &r) -> decltype(l > r ? l : r) { return l > r ? l : r; }
#define PI 3.14159265358979323846

using namespace std;

namespace dcrp
{
	Point::Point(double dx, double dy, double dz)
	{
		x_ = dx;
		y_ = dy;
		z_ = dz;
	}

	Point::Point(const Point &p) : x_(p.x_), y_(p.y_), z_(p.z_) {}

	Point::~Point() {}

	double Point::norm() const
	{
		return sqrt(x_ * x_ + y_ * y_ + z_ * z_);
	}
	double Point::distance(const Point &v) const
	{
		return sqrt((x_ - v.x_) * (x_ - v.x_) + (y_ - v.y_) * (y_ - v.y_) + (z_ - v.z_) * (z_ - v.z_));
	}

	Point Point::normalized() const
	{
		return *this / this->norm();
	}

	Point Point::cross(const Point &v) const
	{
		double x1 = y_ * v.z_ - z_ * v.y_;
		double x2 = z_ * v.x_ - x_ * v.z_;
		double x3 = x_ * v.y_ - y_ * v.x_;
		return Point(x1, x2, x3);
	}
    
    double Point::dot(const Point &v) const{
        return (this->x_*v.x_+this->y_*v.y_+this->z_*v.z_);
    }

	Point Point::rotate(const int &angle) const {
		if (angle == 90)
			return Point(this->y_, -this->x_, this->z_);
		else if (angle == -90)
			return Point(-this->y_, this->x_, this->z_); 
		
		double theta = angle * PI / 180;
		double cos_theta = cos(theta);
		double sin_theta = sin(theta);
		double x1 = x_ * cos_theta + y_ * sin_theta;
		double y1 = -x_ * sin_theta + y_ * cos_theta;
		return Point(x1, y1, z_);
	}

	Point Point::operator+(const Point &v) const
	{
		return Point(x_ + v.x_, y_ + v.y_, z_ + v.z_);
	}

	Point Point::operator-(const Point &v) const
	{
		return Point(x_ - v.x_, y_ - v.y_, z_ - v.z_);
	}

	Point Point::operator*(const double &t) const
	{
		return Point(t * x_, t * y_, t * z_);
	}

	double Point::operator*(const Point &v) const
	{
		return x_ * v.x_ + y_ * v.y_ + z_ * v.z_;
	}

	Point Point::operator/(const double &t) const
	{
		return Point(x_ / t, y_ / t, z_ / t);
	}

	double Point::operator/(const Point &v) const
	{
		return (*this) * v / v.norm() / v.norm();
	}

	Point& Point::operator=(const Point &v)
	{
		x_ = v.x_;
		y_ = v.y_;
		z_ = v.z_;
		return *this;
	}

	double Point::operator[](const int &i) const
	{
		if (i == 0 || i == -3)
			return x_;
		else if (i == 1 || i == -2)
			return y_;
		else if (i == 2 || i == -1)
			return z_;
		else
			return numeric_limits<double>::quiet_NaN();
	}

	bool Point::operator==(const Point& v) const {
		if (((*this) - v).norm() < ERR)
			return true;
		else
			return false;
		return false;
	}

	bool Point::operator!=(const Point& v) const {
		if (((*this) - v).norm() >= ERR)
			return true;
		else
			return false;
		return false;
	}

	bool Point::IsParallel(const Point &v, double REL_ERR) const
	{
		if (this->norm() < REL_ERR || v.norm() < REL_ERR)
		{
			return true;
		}

		return (this->normalized().distance(v.normalized()) < REL_ERR) || ((this->normalized() + v.normalized()).norm() < REL_ERR);
	}

	bool Point::IsSameDirection(const Point &v, double REL_ERR) const
	{
		return IsParallel(v,REL_ERR) && ((*this) / v) > 0;
	}

	bool Point::IsReverseDirection(const Point &v, double REL_ERR) const
	{
		return IsParallel(v,REL_ERR) && ((*this) * v) < 0;
	}
	bool Point::IsWeakParallel(const Point &v, double REL_ERR, double WEAK_PARALLEL_ERR) const
	{
		double norm0 = norm();
		double normv = v.norm();

		if ((norm0 < REL_ERR) || (normv < REL_ERR))
			return true;

		double normal_dist = 0.0;
		double neg_normal_dist = 0.0;

		normal_dist += (x_ / norm0 - v.x_ / normv) * (x_ / norm0 - v.x_ / normv);
		normal_dist += (y_ / norm0 - v.y_ / normv) * (y_ / norm0 - v.y_ / normv);
		normal_dist += (z_ / norm0 - v.z_ / normv) * (z_ / norm0 - v.z_ / normv);
		normal_dist = std::sqrt(normal_dist);

		neg_normal_dist += (x_ / norm0 + v.x_ / normv) * (x_ / norm0 + v.x_ / normv);
		neg_normal_dist += (y_ / norm0 + v.y_ / normv) * (y_ / norm0 + v.y_ / normv);
		neg_normal_dist += (z_ / norm0 + v.z_ / normv) * (z_ / norm0 + v.z_ / normv);
		neg_normal_dist = std::sqrt(neg_normal_dist);
		return (normal_dist < WEAK_PARALLEL_ERR) || (neg_normal_dist < WEAK_PARALLEL_ERR);
	}

	bool Point::IsWeakSameDirection(const Point &v, double REL_ERR, double WEAK_PARALLEL_ERR) const
	{
		return IsWeakParallel(v, REL_ERR, WEAK_PARALLEL_ERR) && ((*this) * v) > 0;
	}
	bool Point::IsWeakReverseDirection(const Point &v, double REL_ERR, double WEAK_PARALLEL_ERR) const
	{
		return IsWeakParallel(v, REL_ERR, WEAK_PARALLEL_ERR) && ((*this) * v) < 0;
	}

	double Point::l1distance(const Point &v) const
	{
		return fabs(x_ - v.x_) + fabs(y_ - v.y_) + fabs(z_ - v.z_);
	}

	bool LineContainPoint(const Point &start, const Point& direc, LineType type, const Point &pnt, double& t)
	{
		Point p = pnt - start;
		if ( !p.IsParallel(direc, 1e-6))
			return false;
		t = p/direc;
		return !(((type != LineType::LINE) && (t < - 1e-6)) || 
				 ((type == LineType::SEGMENT) && (t> 1+1e-6)));
	}

}