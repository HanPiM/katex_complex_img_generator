#pragma once

#include <iostream>
#include <iomanip>

#define _IC inline constexpr

class kunit
{
public:

	_IC kunit()noexcept :_val(0) {}
	_IC kunit(double x)noexcept :_val(x) {}
	_IC kunit(const kunit& x) : _val(x._val) {}
	~kunit() = default;

	_IC bool operator<(kunit x)const { return _val < x._val; }
	_IC bool operator>(kunit x)const { return _val > x._val; }
	_IC bool operator<=(kunit x)const { return _val <= x._val; }
	_IC bool operator>=(kunit x)const { return _val >= x._val; }
	_IC bool operator==(kunit x)const { return _eq(_val, x._val); }
	_IC bool operator!=(kunit x)const { return !_eq(_val, x._val); }
	_IC const kunit& operator+=(const kunit& x) { _val += x._val; return x; }
	_IC const kunit& operator-=(const kunit& x) { _val -= x._val; return x; }
	_IC const kunit& operator*=(double x) { _val *= x; return *this; }
	_IC const kunit& operator/=(double x) { _val /= x; return *this; }

	_IC kunit operator+(kunit x)const { return kunit(_val + x._val); }
	_IC kunit operator-(kunit x)const { return kunit(_val - x._val); }
	_IC kunit operator*(double x)const { return kunit(_val * x); }
	_IC kunit operator/(double x)const { return kunit(_val / x); }

	_IC const kunit& operator=(const kunit& x) { _val = x._val; return x; }

	_IC kunit operator-() { return kunit(-_val); }

	_IC double val()const { return _val; }
	_IC explicit operator double() { return _val; }

	friend std::ostream& operator << (std::ostream& os, const kunit& x)
	{
		int tmp = int(x._val * 100) % 100;
		int w = 2;
		if (tmp == 0)w = 0;
		else if (tmp % 10 == 0)w = 1;
		auto old = os.precision();
		os << std::fixed << std::setprecision(w);

		if (x != 0 && x.abs() < 0.01)os << x._val * 65536 << "sp";
		else os << x._val << "pt";

		os.precision(old);
		return os;
	}

	_IC kunit abs()const { return kunit(_abs(_val)); }

private:
	static _IC double _abs(double x) { return x > 0 ? x : -x; }
	static _IC double _sp(double x) { return x / 65536.0; }
	_IC bool _divisible(double x, double eps = 1)const
	{
		double tmp = _abs(_val);
		return tmp - int(tmp / x) * x <= eps;
	}
	_IC double _as_unit(double base)const
	{
		return _val / base;
	}

	static _IC bool _eq(double a, double b)
	{
		constexpr double eps = double(_sp(0.01));
		return a > b ? (a - b <= eps) : (b - a <= eps);
	}

	double _val;
};

#undef _IC

#define _MAKE(name,expr)\
inline constexpr kunit name(double x){return kunit(expr);}\
inline constexpr kunit operator ""##_##name (unsigned long long _raw)\
{return name(double(_raw));}\
inline constexpr kunit operator ""##_##name (long double _raw)\
{return name(double(_raw));}

_MAKE(pt, x) _MAKE(pc, 12 * x) _MAKE(dd, 1238.0 * x / 1157.0)
_MAKE(cc, 14856.0 * x / 1157.0) _MAKE(nd, 685.0 * x / 642.0)
_MAKE(nc, 1370.0 * x / 107.0) _MAKE(sp, x / 65536.0)
// 1pt = 1/72.27 in => 1in = 72.27pt
_MAKE(in, 72.27 * x)
// 1bp = 1/70 in => 72.27/70 pt
_MAKE(bp, 72.27 * x / 70.0)
// 1px = 1bp ???
_MAKE(px, 72.27 * x / 70.0)
// 1in = 2.54cm, 1in = 72.27pt => 1cm = 7227/254 pt
_MAKE(cm, 7227 * x / 254.0)
// 1mm = 1/100 cm = 7227/25400 pt
_MAKE(mm, 7227 * x / 25400.0)