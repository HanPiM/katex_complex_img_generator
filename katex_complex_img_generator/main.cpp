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

	kunit h = 0_pt;

	serializer s;

	auto draw_range = serializer::rect(-18, -18, 18, 14);

	kunit w = 288_pt;

	//s.add_mode(serializer::mode::for_discuss);

	s.color(0xcccccc);
	s.kern(-s.draw_axes(draw_range, w));

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
			draw_range, 0, 0, [](double x) {return x + 0.5; }
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
			draw_range, 0, 0, [](double x) {return x + 0.5; }
				)
	);

	//s.color(0xff0000);

	//s.kern(
	//	-s.draw(
	//		[](double x)
	//		{
	//			return 15.0 / x;
	//		},
	//		draw_range, 0, 0, [](double x) {return x + 1; }
	//			)
	//);

	//s.kern(
	//	-s.draw(
	//		[](double x)
	//		{
	//			return pow(x, sin(x));
	//		},
	//		draw_range
	//			)
	//);


	s.kern(w);

	s.finish();

	std::cout << "\n$$";
	return 0;
}
