#include "serializer.h"

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <map>

#include <exception>

#include <string>

int main()
{
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

	//s.add_mode(serializer::mode::for_discuss);

	s.color(0x00cccc);
	s.kern(-s.draw_axes(draw_range, 144_pt));

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

	s.kern(
		-s.draw(
			[&t](double x)
			{
				return x * x;
			},
			draw_range, 0, 0, [](double x) {return x + 1; }
				)
	);
	s.finish();

	std::cout << "\n$$";
	return 0;
}
