#pragma once

#include <vector>
#include <functional>

#include <variant>

#include "katex_def.h"

#define _IC inline constexpr

enum class border_style
{
	none,
	expend_inward,
	expend_middle,
	expend_outward
};

class serializer
{
private:
	template <typename _t>
	static _IC const _t& _min(const _t& a, const _t& b)
		noexcept(noexcept(b < a)) {
		return a < b ? a : b;
	}
	template <typename _t>
	static _IC const _t& _max(const _t& a, const _t& b)
		noexcept(noexcept(b > a)) {
		return a > b ? a : b;
	}
	template <typename _t>
	static _IC _t _abs(const _t& x)
		noexcept(noexcept(x > 0, -x)) {
		return x > 0 ? x : -x;
	}

	static _IC bool _deq(double a, double b, double eps = 1e-7)
	{
		return _abs(a - b) <= eps;
	}

public:

	struct mode
	{
		enum
		{
			no_def = 1,
			raw_unit = 1 << 1,
			for_discuss = (1 << 2) | no_def
		};
	};

	class range
	{
	public:
		_IC range()noexcept {}
		_IC range(double a, double b)noexcept :
			_s(_min(a, b)), _e(_max(a, b)) {}
		_IC double begin()const noexcept { return _s; }
		_IC double end()const noexcept { return _e; }
		_IC double size()const noexcept { return _e - _s; }

		_IC bool has(double v)const noexcept
		{
			return (_s <= v && v <= _e) || _deq(v, _s) || _deq(v, _e);
		}
		_IC bool empty()const noexcept { return size() < 1e-8; }

		inline void assign(double s, double e)noexcept
		{
			_s = _min(s, e);
			_e = _max(s, e);
		}

		inline friend std::ostream& operator<<(std::ostream& o, const range& r)
		{
			o << '[' << r.begin() << ',' << r.end() << ']';
			return o;
		}

	private:
		double _s = 0, _e = 0;
	};

	struct rect;

	struct point
	{
		_IC point() noexcept {}
		_IC point(double x, double y)noexcept :x(x), y(y) {}
		_IC point(const std::initializer_list<double> vals)
			: x(*vals.begin()), y(*(vals.begin() + 1)) {}

		_IC bool in_rect(const rect& r)const noexcept
		{
			return r.xrange().has(x) && r.yrange().has(y);
		}

		_IC bool operator ==(const point& p)const
		{
			return _deq(x, p.x) && _deq(y, p.y);
		}

		constexpr void cut_into(const rect& r)noexcept
		{
			if (x < r.sx())x = r.sx();
			if (x > r.ex())x = r.ex();
			if (y < r.sy())y = r.sy();
			if (y > r.ey())y = r.ey();
		}
		constexpr point get_cut_into(const rect& r)const noexcept
		{
			point res = *this;
			res.cut_into(r);
			return res;
		}

		inline void assign(double xx, double yy)noexcept
		{
			x = xx, y = yy;
		}
		inline void assign(const point& p)noexcept
		{
			x = p.x, y = p.y;
		}

		inline friend std::ostream& operator<<(std::ostream& o, const point& p)
		{
			o << '(' << p.x << ',' << p.y << ')';
			return o;
		}

		double x = 0, y = 0;
	};

	struct rect
	{
		_IC rect() noexcept {}
		_IC rect(double x1, double y1, double x2, double y2) noexcept :
			xr(x1, x2), yr(y1, y2) {}
		_IC rect(const point& s, const point& e) :
			xr(s.x, e.x), yr(s.y, e.y) {}

		_IC double width()const noexcept { return xr.size(); }
		_IC double height()const noexcept { return yr.size(); }

		_IC const range& xrange()const noexcept { return xr; }
		_IC const range& yrange()const noexcept { return yr; }

		_IC point spos()const noexcept { return point(xr.begin(), yr.begin()); }
		_IC point epos()const noexcept { return point(xr.end(), yr.end()); }

		_IC double sx()const noexcept { return xr.begin(); }
		_IC double ex()const noexcept { return xr.end(); }
		_IC double sy()const noexcept { return yr.begin(); }
		_IC double ey()const noexcept { return yr.end(); }

		range xr, yr;
	};

	static constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
	static constexpr kunit min_zoom_x = 4_pt;
	static constexpr kunit min_zoom_y = 1_pt;

	inline int get_mode()const noexcept { return _mode; }
	inline void set_mode(int x)noexcept { _mode = x; }
	inline void add_mode(int x)noexcept { _mode |= x; }

	inline void add(const std::string& s)
	{
		_do_last_op();
		_o << s;
	}

	inline void finish() { _do_last_op(); }

	inline void color(unsigned char r, unsigned char g, unsigned char b)
	{
		color((r << 16) | (g << 8) | b);
	}
	
	void color(unsigned int hex);
	void kern(kunit x);
	void rule(kunit w, kunit h, kunit offset = 0);

	void border(kunit w, kunit h,
		kunit offset = 0,
		kunit topw = 0.1, kunit bottomw = 0.1,
		kunit leftw = 0.1, kunit rightw = 0.1,
		border_style style = border_style::none
	);

	inline void raisebox(kunit x)
	{
		_do_last_op();
		if (x != 0)_o << "\\raisebox{" << x << '}';
	}

	void basic_line(kunit w, kunit h = 0_pt);

	using unary_func = std::function<double(double)>;

	kunit draw(
		unary_func f, rect draw_range = rect(),
		kunit w = 0, kunit h = 0,
		unary_func next_pos = [](double x) {return x + 1; }
	);

	kunit draw(
		unary_func fx, unary_func fy,
		unary_func next_t, range t_range,
		rect draw_range,
		kunit w = 0, kunit h = 0
	);

	kunit draw(
		const std::vector<point>& points, rect draw_range,
		kunit w = 0, kunit h = 0
	);

	inline kunit draw_axes(const rect& rc, kunit w = 0, kunit h = 0)
	{
		return draw_axes(rc.sx(), rc.ex(), rc.sy(), rc.ey(), w, h);
	}
	kunit draw_axes(
		double x1, double x2, double y1, double y2,
		kunit w = 0, kunit h = 0
	);

private:

	enum
	{
		_op_color,

		_op_kern

	};
	 
	int _mode = 0;

	int _last_opt = -1;
	std::variant <
		unsigned int,
		kunit
	> _last_op_data;
	std::function<void(const void*)> _last_op_func = [](const void*) {};

	void _do_last_op();

	std::ostream& _o = std::cout;
};
