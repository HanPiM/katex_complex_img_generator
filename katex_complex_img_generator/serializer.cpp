#include "serializer.h"

void serializer::_do_last_op()
{
	const void* p = NULL;
	switch (_last_opt)
	{
#define _CASE(_op) case _op: p = &std::get<_op>(_last_op_data); break

		_CASE(_op_kern);
		_CASE(_op_color);
		
	default:
		break;
#undef _CASE
	}
	//puts("_");
	_last_op_func(p);
	_last_op_func = [](const void*) {};
	_last_opt = -1;
	//puts("_____");
}

void serializer::color(unsigned int hex)
{
	if (_last_opt == _op_color)
	{
		std::get<_op_color>(_last_op_data) = hex;
	}
	else
	{
		_do_last_op();
		_last_opt = _op_color;
		_last_op_data = hex;
		_last_op_func = [this](const void* p)
		{
			unsigned int hex = *(unsigned int*)p;
			auto oldch = _o.fill(); _o.fill('0');
			_o << "\\color{#" << std::hex;
			if (
				(hex & 0xf) == ((hex >> 4) & 0xf)
				&& ((hex >> 8) & 0xf) == ((hex >> 12) & 0xf)
				&& ((hex >> 16) & 0xf) == ((hex >> 20) & 0xf)
				) {
				_o << ((hex >> 16) & 0xf) << ((hex >> 8) & 0xf) << (hex & 0xf);
			}
			else _o << std::setw(6) << hex;
			_o << std::dec << "}";
			_o.fill(oldch);
		};
	}
}

void serializer::kern(kunit x)
{
	if (_last_opt == _op_kern)
	{
		std::get<kunit>(_last_op_data) += x;
		return;
	}
	else
	{
		_do_last_op();
		_last_opt = _op_kern;
		_last_op_data = x;
		_last_op_func = [this](const void* p)
		{
			kunit x = *(kunit*)p;
			if (x != 0)_o << "\\kern{" << x << '}';
		};
	}
}

void serializer::rule(kunit w, kunit h, kunit offset)
{
	if (w == 0 && h == 0)return;
	_do_last_op();
	_o << "\\rule";
	if (offset != 0)_o << '[' << offset << ']';
	_o << '{' << w << "}{" << h << '}';
}

void serializer::border(
	kunit w, kunit h, kunit offset,
	kunit topw, kunit bottomw, kunit leftw, kunit rightw,
	border_style style
)
{
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

void serializer::basic_line(kunit w, kunit h)
{
	if (h.abs() <= 0.5_pt)
	{
		rule(w, 0_pt, 0);
		return;
	}
	kern(2_pt);
	_do_last_op();
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
	_do_last_op();
	_o << "}}";
	kern((_mode & mode::for_discuss) ? 4_pt : 2_pt);
}

kunit serializer::draw(
	unary_func f, rect draw_range,
	kunit w, kunit h, unary_func next_pos
) {
	double pi = acos(-1.0);
	if (draw_range.xr.empty())draw_range.xr.assign(-2.0 * pi, 2.0 * pi);

	return draw(
		[](double t) {return t; },
		[&f](double t) {return f(t); },
		next_pos, draw_range.xrange(),
		draw_range, w, h
	);
	//std::vector<point> points;
	//double x = draw_range.sx();
	//double y = f(x);
	//points.emplace_back(x, y);

	//for (x = next_pos(x); x <= draw_range.ex() || _deq(x, draw_range.ex()); x = next_pos(x))
	//{
	//	y = f(x);
	//	points.emplace_back(x, y);
	//}

	//return draw(points, draw_range, w, h);
}

kunit serializer::draw(
	unary_func fx, unary_func fy,
	unary_func next_t, range t_range,
	rect draw_range, kunit w, kunit h
) {
	std::vector<point> points;

	double t = t_range.begin();
	double x = fx(t);
	double y = fy(t);
	points.emplace_back(x, y);

	for (
		t = next_t(t);
		t <= t_range.end() || _deq(t, t_range.end());
		t = next_t(t)
	) {
		x = fx(t);
		y = fy(t);
		points.emplace_back(x, y);
	}

	return draw(points, draw_range, w, h);
}

kunit serializer::draw(
	const std::vector<point>& points, rect draw_range,
	kunit w, kunit h
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
			if (_deq(lastp.x, p.x))
			{
				drawpos.assign(
					p.x,
					range(lastp.y, p.y).has(yr.begin())
					? yr.begin() : yr.end()
				);
			}
			else // 不平行则需要计算是与哪一边相交
			{
				double k = (lastp.y - p.y) / (lastp.x - p.x);
				double b = p.y - k * p.x;

				double sxy = k * xr.begin() + b;
				double exy = k * xr.end() + b;

				const range& lineseg_xr = range(lastp.x, p.x);

				// 左右两边
				if (lineseg_xr.has(xr.begin()) && yr.has(sxy))
					drawpos.assign(xr.begin(), sxy);
				else if (lineseg_xr.has(xr.end()) && yr.has(exy))
					drawpos.assign(xr.end(), exy);
				else
				{
					// 上下两边
					double ans1 = (yr.begin() - b) / k;
					double ans2 = (yr.end() - b) / k;

					if (lineseg_xr.has(ans1))drawpos.assign(ans1, yr.begin());
					else if (lineseg_xr.has(ans2))drawpos.assign(ans2, yr.end());
					else throw std::logic_error("reached the impossible branch");
					// 可以证明不可能到达这里
				}
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
			|| !points[beg_idx].in_rect(draw_range)
			)
		) {
		offset += points[beg_idx].x - points[beg_idx - 1].x;
		beg_idx++;
	} // 跳过开头那些不符合要求的
	kern(zoom_x * offset);

	bool has_lastdrawpos = 0;
	if (!std::isnan(points[beg_idx - 1].y))
	{
		get_drawpos(points[beg_idx], points[beg_idx - 1], lastdrawpos);
		has_lastdrawpos = !std::isnan(points[beg_idx].y);
	}

	for (size_t i = beg_idx; i < points.size(); ++i)
	{
		const point& p = points[i];
		const point& lastp = points[i - 1];

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
			// 留空不画
			kern(zoom_x * (p.x - lastp.x));
			if (has_lastdrawpos)
			{
				lastdrawpos.assign(p.x - xr.begin(), p.y - yr.begin());
			}
			continue;
		}

		get_drawpos(lastp, p, drawpos);

		if (has_lastdrawpos)
		{
			if (drawpos == lastdrawpos)continue;

			deltax = drawpos.x - lastdrawpos.x;
			deltay = drawpos.y - lastdrawpos.y;

			if (deltax < 0)kern(zoom_x * deltax);
			auto tmp = zoom_y * (lastdrawpos.y + (deltay < 0 ? deltay : 0));
			if (tmp != 0)
			{
				raisebox(tmp);
				_do_last_op();
				_o << "{$";
			}
			basic_line(zoom_x * deltax, zoom_y * deltay);
			if (tmp != 0)
			{
				_do_last_op();
				_o << "$}";
			}
		}

		has_lastdrawpos = true;
		lastdrawpos = drawpos;
	}
	return zoom_x * xr.size();
}

kunit serializer::draw_axes(
	double x1, double x2, double y1, double y2, kunit w, kunit h
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