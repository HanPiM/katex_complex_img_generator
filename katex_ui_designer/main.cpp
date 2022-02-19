#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <map>
#include <functional>

#include <exception>

#include <string>

#include <vector>

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

public:

	struct mode
	{
		enum
		{
			no_def = 1,
			raw_unit = 1 << 1,
			for_discuss = no_def
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
		_IC rect() noexcept{}
		_IC rect(double x1, double y1, double x2, double y2) noexcept :
			xr(x1, x2), yr(y1, y2) {}
		_IC rect(const point& s, const point& e) :
			xr(s.x, e.x), yr(s.y, e.y) {}

		_IC double width()const noexcept{ return xr.size(); }
		_IC double height()const noexcept{ return yr.size(); }

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

	int get_mode()const { return _mode; }
	void set_mode(int x) { _mode = x; }
	void add_mode(int x) { _mode |= x; }

	void add(const std::string& s) { _o << s; }

	void color(unsigned char r, unsigned char g, unsigned char b)
	{
		color((r << 16) | (g << 8) | b);
	}
	void color(unsigned int hex)
	{
		auto oldch = _o.fill(); _o.fill('0');
		_o << "\\color{#" << std::setw(6) << std::hex << hex << std::dec << "}";
		_o.fill(oldch);
	}
	
	void kern(kunit x)
	{
		if (x != 0)_o << "\\kern{" << x << '}';
	}
	void rule(kunit w, kunit h, kunit offset = 0)
	{
		if (w == 0 && h == 0)return;
		_o << "\\rule";
		if (offset != 0)_o << '[' << offset << ']';
		_o << '{' << w << "}{" << h << '}';
	}
	void border(kunit w, kunit h,
		kunit offset = 0,
		kunit topw = 0.1, kunit bottomw = 0.1,
		kunit leftw = 0.1, kunit rightw = 0.1,
		border_style style = border_style::none
	) {
		if (bottomw.abs() <= 0.1_pt)bottomw = 0;
		if (topw.abs() <= 0.1_pt)topw = 0;
		
		if (style != border_style::none)
		{
			switch (style)
			{
			case border_style::expend_inward:break;
			case border_style::expend_middle:
				kern(-leftw / 2); offset -= bottomw / 2;
				w += (leftw + rightw) / 2;
				h += (topw + bottomw) / 2;
				break;
			case border_style::expend_outward:
				kern(-leftw); offset -= bottomw;
				w += leftw + rightw;
				h += topw + bottomw;
				break;
			default:
				return;
				break;
			}
		}

		rule(leftw, h, offset); kern(-leftw);
		rule(w, bottomw, offset); kern(-w);
		rule(w, topw, h - topw + offset); kern(-rightw);
		rule(rightw, h, offset);

		if (style != border_style::none)
		{
			switch (style)
			{
			case border_style::expend_middle:
				kern(-rightw / 2);
				break;
			case border_style::expend_outward:
				kern(-rightw);
				break;
			default:
				break;
			}
		}
	}

	void raisebox(kunit x)
	{
		if (x != 0)_o << "\\raisebox{" << x << '}';
	}

#define _MAKE_DEF(name,def) 

	void basic_line(kunit w, kunit h = 0_pt)
	{
		if (h.abs() <= 0.5_pt)
		{
			rule(w, 0_pt, 0);
			return;
		}
		kern(2_pt);
		_o << '\\';
		if (w < 0)
		{
			w = -w;
			h = -h;
		}
		if (h < 0)
		{
			h = -h;
			_o << 'b';
		}
		_o << "cancel{\\phantom{";
		rule(w, h, 0);
		kern(-4_pt);
		_o << "}}";
		kern((_mode & mode::for_discuss) ? 4_pt : 2_pt);
	}

	kunit draw(
		std::function<double(double)> f, const rect& draw_range = rect(),
		kunit w = 0, kunit h = 0,
		std::function<double(double)> next_pos = [](double x) {return x + 1; }
	)
	{
		std::vector<point> points;
		double x = draw_range.sx();
		double y = f(x);
		points.emplace_back(x, y);

		for (x = next_pos(x); x <= draw_range.ex() || _deq(x, draw_range.ex()); x = next_pos(x))
		{
			y = f(x);
			points.emplace_back(x, y);
		}

		return draw(points, draw_range, w, h);
	}

	kunit draw(
		const std::vector<point>& points, rect draw_range,
		kunit w = 0, kunit h = 0
	) {
		if (points.size() < 2)return 0;

		const point& p0 = points[0];
		double minx = p0.x, maxx = p0.x;
		double miny = p0.y, maxy = p0.y;
		double mindeltax = _abs(points[1].x - p0.x);
		for (size_t i = 1; i < points.size(); ++i)
		{
			auto& p = points[i];
			auto& lastp = points[i - 1];
			minx = _min(minx, p.x);
			maxx = _max(maxx, p.x);
			miny = _min(miny, p.y);
			maxy = _max(maxy, p.y);
			if (!_deq(p.x, lastp.x))
				mindeltax = _min(mindeltax, _abs(p.x - lastp.x));
		}
		if (draw_range.xrange().empty())draw_range.xr = range(minx, maxx);
		if (draw_range.yrange().empty())draw_range.yr = range(miny, maxy);

		range& xr = draw_range.xr, & yr = draw_range.yr;

		if (xr.end() < minx || xr.begin() > maxx)return 0;
		if (yr.end() < miny || yr.begin() > maxy)return 0;

		kunit zoom_x = w / xr.size();
		kunit zoom_y = h / yr.size();
		if (zoom_x * mindeltax < min_zoom_x)zoom_x = min_zoom_x / mindeltax;
		if (zoom_y == 0)zoom_y = zoom_x;
		if (zoom_y < min_zoom_y)zoom_y = min_zoom_y;
		//printf("(%lf,%lf,%lf,%lf)\n", minx, miny, maxx, maxy);

		//std::cout << "zx" << zoom_x << " zy" << zoom_y << "\n";

		/*
		* 如果 p 在范围内或 p.x == lastp.x 则结果会与 p 有极大关系
		* 此时传入的 lastp, p 的位置调换会有不同结果
		* 反之则无影响
		*/
		auto get_drawpos = [&draw_range, &xr, &yr](const point& lastp, const point& p, point& drawpos)
		{
			if (p.in_rect(draw_range))drawpos = p;
			else
			{
				if (!_deq(lastp.x, p.x))
				{
					double k = (lastp.y - p.y) / (lastp.x - p.x);
					double b = p.y - k * p.x;

					double sxy = k * xr.begin() + b;
					double exy = k * xr.end() + b;

					//printf("\n__%lf x + %lf\n", k, b);

					const range& tmp = range(lastp.x, p.x);

					//std::cout << xr << yr << '\n' << tmp << '\n';


					// 左右两边
					if (tmp.has(xr.begin()) && yr.has(sxy))
						drawpos.assign(xr.begin(), sxy);
					else if (tmp.has(xr.end()) && yr.has(exy))
						drawpos.assign(xr.end(), exy);
					else
					{
						// 上下两边
						double ans1 = (yr.begin() - b) / k;
						double ans2 = (yr.end() - b) / k;

						//printf("\n__%lf %lf ans1:%lf ans2:%lf\n",yr.begin(),yr.end(), ans1, ans2);

						if (tmp.has(ans1))drawpos.assign(ans1, yr.begin());
						else if (tmp.has(ans2))drawpos.assign(ans2, yr.end());
						else throw std::logic_error("reached the impossible branch");
						// 可以证明不可能到达这里

						//std::cout << drawpos;
					}
				}
				else // 直线
				{
					drawpos.assign(
						p.x,
						range(lastp.y, p.y).has(yr.begin())
						? yr.begin() : yr.end()
					);
				}
			}

			drawpos.x -= xr.begin();
			drawpos.y -= yr.begin();

		};

		double deltax, deltay;
		point drawpos, lastdrawpos;

		size_t beg_idx = 1;
		double offset = 0;
		while (
			beg_idx < points.size()
			&& (
				std::isnan(points[beg_idx].y)
				||!points[beg_idx].in_rect(draw_range)
			)
		) {
			offset += points[beg_idx].x - points[beg_idx - 1].x;
			beg_idx++;
		}

		//printf("beg_idx %d offset %lf\n", beg_idx, offset);

		kern(zoom_x * offset);

		bool has_lastdrawpos = 0;
		if (!std::isnan(points[beg_idx - 1].y))
		{
			get_drawpos(points[beg_idx], points[beg_idx - 1], lastdrawpos);
			has_lastdrawpos = !std::isnan(points[beg_idx].y);
		}
		
//		printf("(%lf,%lf)->(%lf,%lf)\n",
//	lastdrawpos.x+xr.begin(), lastdrawpos.y+yr.begin(),
//	drawpos.x+xr.begin(), drawpos.y+yr.begin()
//);
		for (size_t i = beg_idx; i < points.size(); ++i)
		{
			const point& p = points[i];
			const point& lastp = points[i - 1];

			//std::cout << "\n[cur]:" << p << "\n";

			bool no_draw = false;
			if (std::isnan(p.y) || std::isnan(lastp.y))
			{
				no_draw = true;
				has_lastdrawpos = !std::isnan(p.y);
			}
			if (!lastp.in_rect(draw_range) && !p.in_rect(draw_range))
			{
				no_draw = true;
				has_lastdrawpos = p.in_rect(draw_range);
			}
			if (no_draw)
			{
				// 留空
				kern(zoom_x* (p.x - lastp.x));
				if (has_lastdrawpos)
				{
					lastdrawpos.assign(p.x - xr.begin(), p.y - yr.begin());
				}
				continue;
			}

			//std::cout << "\np:" << p << '\n';

			get_drawpos(lastp, p, drawpos);

			if (has_lastdrawpos)
			{
				if (drawpos == lastdrawpos)
				{
					continue;
				}

				//printf("\n__(%lf,%lf)->(%lf,%lf) ", lastp.x, lastp.y, p.x, p.y);
				//printf("(%lf,%lf)->(%lf,%lf)\n\n",
				//	lastdrawpos.x + xr.begin(), lastdrawpos.y + yr.begin(),
				//	drawpos.x + xr.begin(), drawpos.y + yr.begin()
				//);

				deltax = drawpos.x - lastdrawpos.x;
				deltay = drawpos.y - lastdrawpos.y;

				if (deltax < 0)kern(zoom_x * deltax);
				auto tmp = zoom_y * (lastdrawpos.y + (deltay < 0 ? deltay : 0));
				if (tmp != 0)
				{
					raisebox(tmp);
					_o << "{$";
				}
				basic_line(zoom_x * deltax, zoom_y * deltay);
				if (tmp != 0)_o << "$}";
			}
			else
			{
				//std::cout << "lastp:" << lastp;
				//puts("___");
			}

			has_lastdrawpos = true;
			lastdrawpos = drawpos;
		}
		return zoom_x * xr.size();
	}

	inline kunit draw_axes(const rect& rc, kunit w = 0, kunit h = 0)
	{
		return draw_axes(rc.sx(), rc.ex(), rc.sy(), rc.ey(), w, h);
	}
	kunit draw_axes(
		double x1, double x2, double y1, double y2,
		kunit w = 0, kunit h = 0
	) {
		double sizx = _abs(x1 - x2), sizy = _abs(y1 - y2);
		kunit zoom_x = w / sizx;
		kunit zoom_y = h / sizy;
		if (zoom_x < min_zoom_x)zoom_x = min_zoom_x;
		if (zoom_y == 0)zoom_y = zoom_x;
		if (zoom_y < min_zoom_y)zoom_y = min_zoom_y;

		rule(zoom_x * sizx, 0_sp, zoom_y * _abs(y1));
		double tmp = _abs(x1);
		kern(-zoom_x * sizx + zoom_x * tmp);
		rule(1_sp, zoom_y * sizy);
		return zoom_x * tmp;
	}

	//void draw_func(
	//	std::function<double(double)> f, const rect& show_range,
	//	double xunit = 0, double yunit = 0,
	//)
	//{
	//	double last = f(s), offset;
	//	double tmp;
	//	double miny = last, maxy = last;
	//	std::vector<double> points;
	//	for (double x = s + xunit_length; x <= t; x += xunit_length)
	//	{
	//		tmp = f(x);
	//		miny = std::min(miny, tmp);
	//		maxy = std::max(maxy, tmp);
	//		points.push_back(tmp);
	//	}
	//
	//	if (_deq(yunit_length, 0))yunit_length = xunit_length;
	//
	//	kunit xkunit_length = w / ((t - s) / xunit_length - 1);
	//	kunit ykunit_length = h / ((maxy - miny) / yunit_length);
	//
	//	if (xkunit_length == 0)xkunit_length = 5_pt;
	//	if (ykunit_length == 0)ykunit_length = xkunit_length * (yunit_length / xunit_length);
	//
	//	auto x_tokunit = [&xunit_length, &xkunit_length](double x)
	//	{
	//		return xkunit_length * (x / xunit_length);
	//	};
	//	auto y_tokunit = [&yunit_length, &ykunit_length](double y)
	//	{
	//		return ykunit_length * (y / yunit_length);
	//	};
	//
	//	for (double y : points)
	//	{
	//		offset = y - last;
	//		tmp = last + (offset < 0 ? offset : 0);
	//		if (!_deq(tmp, 0))_o << "\\raisebox{" << y_tokunit(tmp) << "}{$";
	//		line(xkunit_length, y_tokunit(offset));
	//		if (!_deq(tmp, 0))_o << "$}";
	//		last = y;
	//	}
	//}

private:

	static _IC bool _deq(double a, double b, double eps = 1e-7)
	{
		return _abs(a - b) <= eps;
	}

	using _func_t = std::function<void(serializer&,int)>;

	int _mode = 0;

	std::map<std::string, _func_t> _func_map;

	std::ostream& _o = std::cout;
};

//class lexer
//{
//public:
//	using int_type = int;
//	using size_type = size_t;
//	using char_type = unsigned char;
//	using getchar_func = std::function<char_type()>;
//
//	using string_type = std::basic_string<
//		char_type, std::char_traits<char_type>,
//		std::allocator<char_type>
//	>;
//	enum class token_type
//	{
//		unknown,
//
//		number,
//		command,
//		text,
//		string,
//
//		l_paren,	// (
//		r_paren,	// )
//		l_bracket,	// [
//		r_bracket,	// ]
//		l_brace,	// {
//		r_brace,	// }
//
//		add,		// +
//		sub,		// -
//		mul,		// *
//		div,		// /
//
//		add_eq,		// +=
//		sub_eq,		// -=
//		mul_eq,		// *=
//		div_eq,		// /=
//
//		assign,		// =
//
//		logic_and,	// &&
//		logic_or,	// ||
//		logic_not,	// !
//
//		eq,			// ==
//		neq,		// !=
//		greate,		// >
//		greate_eq,	// >=
//		less,		// <
//		less_eq,	// <=
//
//		comma,		// ,
//
//		question,	// ?
//		colon		// :
//	};
//	struct token
//	{
//		token_type type = token_type::unknown;
//		string_type data;
//		inline const char* c_str() { return (const char*)data.c_str(); }
//		const char* type_name()
//		{
//			switch (type)
//			{
//			case lexer::token_type::unknown: return "unknown";
//			case lexer::token_type::number: return "number";
//			case lexer::token_type::command: return "command";
//			case lexer::token_type::text: return "text";
//			case lexer::token_type::string: return "string";
//			case lexer::token_type::l_paren: return "l_paren";
//			case lexer::token_type::r_paren: return "r_paren";
//			case lexer::token_type::l_bracket: return "l_bracket";
//			case lexer::token_type::r_bracket: return "r_bracket";
//			case lexer::token_type::l_brace: return "l_brace";
//			case lexer::token_type::r_brace: return "r_brace";
//			case lexer::token_type::add: return "add";
//			case lexer::token_type::sub: return "sub";
//			case lexer::token_type::mul: return "mul";
//			case lexer::token_type::div: return "div";
//			case lexer::token_type::add_eq: return "add_eq";
//			case lexer::token_type::sub_eq: return "sub_eq";
//			case lexer::token_type::mul_eq: return "mul_eq";
//			case lexer::token_type::div_eq: return "div_eq";
//			case lexer::token_type::assign: return "assign";
//			case lexer::token_type::logic_and: return "logic_and";
//			case lexer::token_type::logic_or: return "logic_or";
//			case lexer::token_type::logic_not: return "logic_not";
//			case lexer::token_type::eq: return "eq";
//			case lexer::token_type::neq: return "neq";
//			case lexer::token_type::greate: return "greate";
//			case lexer::token_type::greate_eq: return "greate_eq";
//			case lexer::token_type::less: return "less";
//			case lexer::token_type::less_eq: return "less_eq";
//			case lexer::token_type::comma: return "comma";
//			case lexer::token_type::question: return "question";
//			case lexer::token_type::colon: return "colon";
//			default:
//				break;
//			}
//			return "unknown";
//		}
//	};
//
//	lexer(std::istream& is)
//		:_getch([&is]() {char ch = 0; is.get(ch); return char_type(ch); }),
//		_cur_ch(0), _next_ch(_getch()) {getnext_char();}
//	lexer(getchar_func f) :_getch(f), _cur_ch(0), _next_ch(_getch()) { getnext_char(); }
//
//	/**
//	 * @brief 将下一个字符读入到当前字符
//	 * @return 更新后的当前字符
//	*/
//	char_type getnext_char()
//	{
//		_cur_ch = _next_ch;
//		_next_ch = _getch();
//		_cur_pos++;
//		_check_endl();
//		return _cur_ch;
//	}
//	/**
//	 * @brief 
//	 * @return 当前字符
//	*/
//	inline char_type getcur_char()const { return _cur_ch; }
//	/**
//	 * @brief 查看下一个字符但不读入
//	 * @return 下一个字符
//	*/
//	inline char_type looknext_char()const { return _next_ch; }
//	/**
//	 * @brief 尝试匹配字符，如是期望值则读入
//	 * @param ch 要进行匹配的字符
//	 * @return 是否为期望值
//	*/
//	bool match_char(char_type ch)
//	{
//		if (looknext_char() == ch)
//		{
//			getnext_char();
//			return true;
//		}
//		return false;
//	}
//	/**
//	 * @brief 跳过空白内容
//	*/
//	void skip_blanks()
//	{
//		while (isspace(_cur_ch))
//		{
//			getnext_char();
//		}
//	}
//	/**
//	 * @brief 跳过一行
//	*/
//	void skip_aline()
//	{
//		getnext_char();
//		while (_cur_ch != '\0')
//		{
//			if (_cur_ch == '\n')
//			{
//				getnext_char();
//				break;
//			}
//			getnext_char();
//		}
//	}
//	/**
//	 * @brief 跳过注释
//	*/
//	inline void skip_comment() { skip_aline(); }
//
//	token& getcur_token() { return _cur_tok; }
//
//	inline bool is_end() { return _cur_ch == 0; }
//
//	void getnext_token()
//	{
//#define _CASE1(ch,_type) case ch: _cur_tok.type = token_type::##_type; break
//#define _CASE2(ch,_type1,_type2) \
//		case ch: _cur_tok.type = match_char('=') \
//			? token_type::##_type1 \
//			: token_type::##_type2;break
//
//		skip_blanks();
//		_cur_tok.type = token_type::unknown;
//		_cur_tok.data.clear();
//
//		while (1)
//		{
//			switch (_cur_ch)
//			{
//				_CASE1('(', l_paren);
//				_CASE1(')', r_paren);
//				_CASE1('[', l_bracket);
//				_CASE1(']', r_bracket);
//				_CASE1('{', l_brace);
//				_CASE1('}', r_brace);
//			case '%':skip_comment(); break;
//			case '@':
//				_is_expr_mode = !_is_expr_mode;
//				getnext_char(); skip_blanks();
//				continue;
//			default:
//				constexpr auto sepchars = "?,:+-*/><!=|&\\";
//				if (_is_expr_mode && strchr(sepchars, _cur_ch) != NULL)
//				{
//					switch (_cur_ch)
//					{
//						_CASE1('?', question);
//						_CASE1(',', comma);
//						_CASE1(':', colon);
//
//						_CASE2('+', add_eq, add);
//						_CASE2('-', sub_eq, sub);
//						_CASE2('*', mul_eq, mul);
//						_CASE2('/', div_eq, div);
//						_CASE2('>', greate_eq, greate);
//						_CASE2('<', less_eq, less);
//						_CASE2('!', neq, logic_not);
//						_CASE2('=', eq, assign);
//
//					case '|':
//						if (match_char('|'))_cur_tok.type = token_type::logic_or;
//						else
//						{
//							_err("bit_or is disabled");
//							return;
//						}
//						break;
//					case '&':
//						if (match_char('&'))_cur_tok.type = token_type::logic_or;
//						else
//						{
//							_err("bit_and is disabled");
//							return;
//						}
//						break;
//					case '\\':
//					{
//						getnext_char();
//						_cur_tok.type = token_type::command;
//						while (!is_end() && isalpha(_cur_ch) && strchr(sepchars, _cur_ch) == NULL)
//						{
//							_cur_tok.data.push_back(_cur_ch);
//							getnext_char();
//						}
//						return;
//					}
//					default:
//						break;
//					}
//					break;
//				}
//				_cur_tok.type = token_type::text;
//				while (!is_end() && strchr("()[]{}%@", _cur_ch) == NULL && !isspace(_cur_ch))
//				{
//					_cur_tok.data.push_back(_cur_ch);
//					getnext_char();
//				}
//				return;
//			}
//			getnext_char();
//			break;
//		}
//		
//#undef _CASE1
//#undef _CASE2
//	}
//
//private:
//
//	void _warn(const std::string& msg)
//	{
//
//	}
//	void _err(const std::string& msg)
//	{
//
//	}
//
//	static bool _is_idchar(char_type ch, bool first = false)
//	{
//		return isalpha(ch) || ch == '_';
//	}
//
//	bool _check_endl()
//	{
//		if (_cur_ch == '\n')
//		{
//			_cur_line++;
//			_cur_pos = 0;
//			return true;
//		}
//		return false;
//	}
//
//	void _get_string()
//	{
//		_cur_tok.type = token_type::string;
//		_cur_tok.data.clear();
//		while (true)
//		{
//			getnext_char();
//			if (_check_endl())
//			{
//				_err("string is terminated by \\n");
//				return;
//			}
//			if (_cur_ch == '\0')
//			{
//				_err("unterminated string");
//				return;
//			}
//			if (_cur_ch == '"')break;
//			if (_cur_ch == '\\')
//			{
//				getnext_char();
//				char_type ch = 0;
//				switch (_cur_ch)
//				{
//				case 'n':ch = '\n'; break;
//				case 't':ch = '\t'; break;
//				case '\\':
//				case '"':ch = _cur_ch; break;
//				default:
//					ch = _cur_ch;
//					_warn("unsupport escape \\" + std::to_string(ch));
//					break;
//				}
//			}
//			else _cur_tok.data.push_back(_cur_ch);
//		}
//	}
//
//	getchar_func _getch;
//	char_type _cur_ch;
//	char_type _next_ch;
//	size_type _cur_line = 0;
//	size_type _cur_pos = 0;
//
//	size_type _tab_size = 4;
//
//	bool _is_expr_mode = 0;
//
//	token _cur_tok;
//};
//
//void expand_macro(std::istream& is, std::ostream& os)
//{
//	typename std::istream::char_type ch;
//	std::map<std::string, std::function<void()> > macros;
//	std::string name;
//	while (is >> ch)
//	{
//		if (ch == '%')
//		{
//			while (is >> ch)
//			{
//				if (ch == '\n')break;
//			}
//		}
//		else if (ch == '\\')
//		{
//			name.clear();
//			while (is >> ch)
//			{
//				if (!isalpha(ch))break;
//				name.push_back(ch);
//			}
//			auto it = macros.find(name);
//			if (it != macros.end())it->second();
//			else
//			{
//				os << '\\' << name;
//			}
//		}
//		else os << ch;
//	}
//}

int main()
{
//	constexpr kunit a = 1145_sp;
//	constexpr double b = a.val() + 1;
//
//#define out(x) std::cout<<#x<<": "<<1##_##x<<'\n'
//	
//	out(pt);
//	out(in);
//	out(bp);
//	out(pc);
//	out(dd);
//	out(cc);
//	out(nc);
//	out(sp);
//	out(cm);
//	out(mm);

	//std::cout << a;

	std::cout << "$$\n";

	kunit w = 200_pt;
	kunit h = 0_pt;

	serializer s;

	//s.color(255, 0, 0);

	//s.rule(1_sp, 400_pt);
	//s.kern(-1_sp);
	//s.rule(200_pt, 0);

	//s.kern(-200_pt);

	//s.color(0);

	auto draw_range = serializer::rect(-18, -18, 18, 14);

	s.add_mode(serializer::mode::for_discuss);

	s.color(0xcccccc);
	s.kern(-s.draw_axes(draw_range, 288_pt));

	//s.color(0);

	s.add("\\color{red}");

	auto t = [](double x)
	{
		double s = x < 0 ? -1 : 1;
		return asin(s * pow(abs(x) / 16.0, 1.0 / 3.0));
	};

	s.kern(
		-s.draw(
			[&t](double x)
			{
				if (abs(abs(x) - 16) < 1e-8)x = x < 0 ? -16 : 16;
				return 13 * cos(t(x)) - 5 * cos(2 * t(x)) - 2 * cos(3 * t(x)) - cos(4 * t(x));
			},
			draw_range, 0, 0, [](double x) {return x + 1; }
		)
	);
	s.add("\\color{green}");

	s.kern(
		-s.draw(
			[&t](double x)
			{
				if (abs(abs(x) - 16) < 1e-8)x = x < 0 ? -16 : 16;
				return -13 * cos(t(x)) - 5 * cos(2 * t(x)) + 2 * cos(3 * t(x)) - cos(4 * t(x));
			},
			draw_range, 0, 0, [](double x) {return x + 1; }
				)
	);
	//s.draw_func(
	//	[](double x) {
	//		return abs(x * x - 5 * x + 3);
	//	},
	//	0, 5,
	//	0.1, 0,
	//	w, h
	//);

	//s.kern(-w);
	//s.color(0, 255, 0);

	//s.draw_func(
	//	[](double x) {return pow(x, -x); },
	//	0, 2,
	//	0.1, 0,
	//	w, h
	//);

	//s.kern(-w);
	//s.color(0, 0, 255);

	//s.draw_func(
	//	[](double x) {return pow(x, 2); },
	//	0, 2,
	//	0.1, 0,
	//	w, h
	//);

	std::cout << "\n$$";
	return 0;
}
