//test_rt.cpp
//

#include <utility>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <print>

import ts;
import rt;

using namespace ts;

suite test_tests {
	{"Test1", [](auto env) {
		env.check(true);
	}},
	{"Test2", [](auto env) {
		env.check(false, "hi");
	}},
	scenario("Test3") = [](auto env) {
		env.check(false);
	},
	expect_that("tests should work") = [](auto env) {
		env.check(false);
		env.check(true);
		env.check(false);
	},
};

auto default_world() 
{
	using namespace rt; 
	point_light light {point(-10, 10, -10), color(1, 1, 1)};
	sphere s1;
	s1.surface = {
			.color = color(0.8, 1.0, 0.6),
			.diffuse = 0.7,
			.specular = 0.2};
	sphere s2;
	s2.transform = scaling(0.5, 0.5, 0.5);
	world w;
	w.objects.push_back(s1);
	w.objects.push_back(s2);
	w.light = light;
	return w;
}

suite rt_tests {
	scenario("A tuple with w=1.0 is a point") = [](auto check) {
		auto a = rt::tuple {4.3, -4.2, 3.1, 1.0};
		check(a.x == 4.3);
		check(a[0] == 4.3);
		check(a.y == -4.2);
		check(a[1] == -4.2);
		check(a.z == 3.1);
		check(a[2] == 3.1);
		check(a.w == 1.0);
		check(a[3] == 1.0);
		check(a.is_point());
		check(!a.is_vector());
		a.x = 4.0;
		check(a[0] == 4.0);
		a[0] = 3.0;
		check(a.x == 3.0);
	},
	scenario("A tuple with w=0 is a vector") = [](auto check) {
		auto a = rt::tuple {4.3, -4.2, 3.1, 0.0};
		check(a.x == 4.3);
		check(a.y == -4.2);
		check(a.z == 3.1);
		check(a.w == 0.0);
		check(!a.is_point());
		check(a.is_vector());
	},
	scenario("point() creates tuples with w=1") = [](auto check) {
		auto p = rt::point(4, -4, 3);
		check(p == rt::tuple{4, -4, 3, 1});
	},
	scenario("vector() creates tupels with w=0") = [](auto check) {
		auto v = rt::vector(4, -4, 3);
		check(v == rt::tuple{4, -4, 3, 0});
	},
	scenario("adding two tuples") = [](auto check) {
		using namespace rt;
		auto a1 = tuple{3, -2, 5, 1};
		auto a2 = tuple{-2, 3, 1, 0};
		check(a1 + a2 == tuple{1, 1, 6, 1});
	},
	scenario("subtracting a vector from a point") = [](auto check) {
		using namespace rt;
		auto p = point(3, 2, 1);
		auto v = vector(5, 6, 7);
		check(p - v == point(-2, -4, -6));
	},
	scenario("subtracting two vectors") = [](auto check) {
		using namespace rt;
		auto v1 = vector(3, 2, 1);
		auto v2 = vector(5, 6, 7);
		check(v1 - v2 == vector(-2, -4, -6));
	},
	scenario("subtracting a vector from a zero vector") = [](auto check) {
		using namespace rt;
		auto zero = vector(0, 0, 0);
		auto v = vector(1, -2, 3);
		check(zero - v == vector(-1, 2, -3));
	},
	scenario("negating a tuple") = [](auto check) {
		using namespace rt;
		auto a = tuple{1, -2, 3, -4};
		check(-a == tuple{-1, 2, -3, 4});
	},
	scenario("Multiplying a tuple by a scalar") = [](auto check) {
		using namespace rt;
		auto a = tuple{1, -2, 3, -4};
		check(a * 3.5 == tuple{3.5, -7, 10.5, -14});
	},
	scenario("dividing a tuple by a scalar") = [](auto check) {
		using namespace rt;
		auto a = tuple{1, -2, 3, -4};
		check(a / 2.0 == tuple{0.5, -1, 1.5, -2});
	},
	scenario("compute magnitudes") = [](auto check) {
		using namespace rt;
		{
			auto v = vector(1, 0, 0);
			check(magnitude(v) == 1.0);
		}{
			auto v = vector(0, 1, 0);
			check(magnitude(v) == 1.0);
		}{
			auto v = vector(0, 0, 1);
			check(magnitude(v) == 1.0);
		}{
			auto v = vector(1, 2, 3);
			check(magnitude(v) == std::sqrt(14.0));
		}{
			auto v = vector(-1, -2, -3);
			check(magnitude(v) == std::sqrt(14.0));
		}
	},
	scenario("normalizing vector") = [](auto check) {
		using namespace rt;
		{
			auto v = vector(4, 0, 0);
			check.that(normalize(v) == vector(1, 0, 0));
		} {
			auto v = vector(1, 2, 3);
			check.that(normalize(v) == vector(0.26726, 0.53452, 0.80178));
		}
	},
	scenario("magnitude of normal vect") = [](auto check) {
		using namespace rt;
		auto v = vector(1, 2, 3);
		auto norm = normalize(v);
		check(magnitude(norm) == 1.0);
	},
	scenario("dot product of two tuples") = [](auto check) {
		using namespace rt;
		auto a = vector(1, 2, 3);
		auto b = vector(2, 3, 4);
		check(dot(a, b) == 20.0);
	},
	scenario("cross product of two vectors") = [](auto check) {
		using namespace rt;
		auto a = vector(1, 2, 3);
		auto b = vector(2, 3, 4);
		check(cross(a,b) == vector(-1, 2, -1));
		check(cross(b,a) == vector(1, -2, 1));
	},
	scenario("colors are red, green blue tuples") = [](auto check) {
		using namespace rt;
		auto c = color(-0.5, 0.4, 1.7);
		check(red(c) == -0.5);
		check(green(c) == 0.4);
		check(blue(c) == 1.7);
	},
	scenario("adding colors") = [](auto check) {
		using namespace rt;
		auto c1 = color(0.9, 0.6, 0.75);
		auto c2 = color(0.7, 0.1, 0.25);
		check(c1 + c2 == color(1.6, 0.7, 1.0));
	},
	scenario("subtracting colors") = [](auto check) {
		using namespace rt;
		auto c1 = color(0.9, 0.6, 0.75);
		auto c2 = color(0.7, 0.1, 0.25);
		check(c1 - c2 == color(0.2, 0.5, 0.5));
	},
	scenario("multiplying a color by a scalar") = [](auto check) {
		using namespace rt;
		auto c = color(0.2, 0.3, 0.4);
		check(c * 2 == color(0.4, 0.6, 0.8));
	},
	scenario("multiplying a color by a scalar") = [](auto check) {
		using namespace rt;
		auto c1 = color(1, 0.2, 0.4);
		auto c2 = color(0.9, 1, 0.1);
		check(c1 * c2 == color(0.9, 0.2, 0.04));
	},
	scenario("creating a canvas") = [](auto check) {
		using namespace rt;
		auto c = canvas{10, 20};
		check(c.width() == 10);
		check(c.height() == 20);
		check(c[0,0] == color(0, 0, 0));
	},
	scenario("writing pixels") = [](auto check) {
		using namespace rt;
		auto c = canvas{10, 20};
		auto red = color(1, 0, 0);
		c[2,3] = red;
		check(c[2,3] == red);
	},
	scenario("constructing ppm") = [](auto check) {
		using namespace rt;
		auto c = canvas(5, 3);
		c[0,0] = color(1.5, 0, 0);
		c[2,1] = color(0, 0.5, 0);
		c[4,2] = color(-0.5, 0, 1);
		auto ppm = c.to_ppm();
		check(ppm == 
				"P3\n5 3\n255\n"
				"255 0 0\n0 0 0\n0 0 0\n0 0 0\n0 0 0\n"
				"0 0 0\n0 0 0\n0 128 0\n0 0 0\n0 0 0\n"
				"0 0 0\n0 0 0\n0 0 0\n0 0 0\n0 0 255\n"
				);
	},
	scenario("constructing and inspecting a 4x4 matrix") = [](auto check) {
		using namespace rt;
		auto m = matrix<4>{1, 2, 3, 4,
			5.5, 6.5, 7.5, 8.5,
			9, 10, 11, 12,
			13.5, 14.5, 15.5, 16.5};
		check(m[0,0] == 1);
		check(m[0,3] == 4);
		check(m[1,0] == 5.5);
		check(m[1,2] == 7.5);
		check(m[2,2] == 11);
		check(m[3,0] == 13.5);
		check(m[3,2] == 15.5);
	},
	scenario("matrix equality") = [](auto check) {
		using namespace rt;
		matrix a {1, 2, 3, 4, 5, 6, 7, 8,
			9, 8, 7, 6, 5, 4, 3, 2};
		matrix b {1, 2, 3, 4, 5, 6, 7, 8,
			9, 8, 7, 6, 5, 4, 3, 2};
		check(a == b);
	},
	scenario("matrix inequality") = [](auto check) {
		using namespace rt;
		matrix a {1, 2, 3, 4, 5, 6, 7, 8,
			9, 8, 7, 6, 5, 4, 3, 2};
		matrix b {2, 3, 4, 5, 6, 7, 8, 9,
			8, 7, 6, 5, 4, 3, 2, 1};
		check(a != b);
	},
	scenario("multiply two matrices") = [](auto check) {
		using namespace rt;
		matrix a {1, 2, 3, 4, 5, 6, 7, 8,
			9, 8, 7, 6, 5, 4, 3, 2};
		matrix b {-2, 1, 2, 3,
			3, 2, 1, -1,
			4, 3, 6, 5,
			1, 2, 7, 8};
		check(a * b == matrix{20, 22, 50, 48,
				44, 54, 114, 108,
				40, 58, 110, 102,
				16, 26, 46, 42});
	},
	scenario("matrix multiplied by a tuple") = [](auto check) {
		using namespace rt;
		matrix a {1, 2, 3, 4, 
			2, 4, 4, 2,
			8, 6, 4, 1,
			0, 0, 0, 1};
		tuple b {1, 2, 3, 1};

		check(a * b == tuple{18, 24, 33, 1});
	},
	scenario("multiplying a matrix by the identity matrix") = [](auto check) {
		using namespace rt;
		matrix a {0, 1, 2, 4,
			1, 2, 4, 8,
			2, 4, 8, 16,
			4, 8, 16, 32};
		check(a * identity == a);
	},
	scenario("transposing a matrix") = [](auto check) {
		using namespace rt;
		matrix a {0, 9, 3, 0,
			9, 8, 0, 8,
			1, 8, 5, 3,
			0, 0, 5, 8};
		check(transpose(a) == matrix{
				0, 9, 1, 0,
				9, 8, 8, 0,
				3, 0, 5, 5,
				0, 8, 3, 8});
	},
	scenario("calculating the determinant of a 2x2 matrix") = [](auto check) {
		using namespace rt;
		matrix<2> a {1, 5, -3, 2};
		check(determinant(a) == 17);
	},
	scenario("a submatrix of a 3x3 matrix is a 2x2 matrix") = [](auto check) {
		using namespace rt;
		matrix<3> a {1, 5, 0,
			-3, 2, 7,
			0, 6, -3};
		check(a.submatrix(0,2) == matrix<2>{-3, 2, 0, 6});
	},
	scenario("a submatrix of a 4x4 matrix is a 3x3 matrix") = [](auto check) {
		using namespace rt;
		matrix<4> a = {
			-6, 1, 1, 6,
			-8, 5, 8, 6,
			-1, 0, 8, 2,
			-7, 1, -1, 1};
		check(a.submatrix(2,1) == matrix<3>{
				-6, 1, 6,
				-8, 8, 6,
				-7, -1, 1});
	},
	scenario("calculating a minor of a 3x3 matrix") = [](auto check) {
		using namespace rt;
		matrix<3> a {
			3, 5, 0,
			2, -1, -7,
			6, -1, 5};
		auto b = a.submatrix(1,0);
		check(determinant(b) == 25);
		check(a.minor(1,0) == 25);
	},
	scenario("calculating a cofactor of a 3x3 matrix") = [](auto check) {
		using namespace rt;
		matrix<3> a {3, 5, 0,
			2, -1, -7,
			6, -1, 5};
		check(a.minor(0,0) == -12);
		check(a.cofactor(0,0) == -12);
		check(a.minor(1,0) == 25);
		check(a.cofactor(1,0) == -25);
	},
	scenario("calculating the determinant of a 3x3 matrix") = [](auto check) {
		using namespace rt;
		matrix<3> a {1, 2, 6, -5, 8, -4, 2, 6, 4};
		check(a.cofactor(0, 0) == 56);
		check(a.cofactor(0, 1) == 12);
		check(a.cofactor(0, 2) == -46);
		check(determinant(a) == -196);
	},
	scenario("calculating the determinant of a 4x4 matrix") = [](auto check) {
		using namespace rt;
		matrix a {
			-2, -8, 3, 5,
			-3, 1, 7, 3,
			1, 2, -9, 6,
			-6, 7, 7, -9};
		check(a.cofactor(0, 0) == 690);
		check(a.cofactor(0, 1) == 447);
		check(a.cofactor(0, 2) == 210);
		check(a.cofactor(0, 3) == 51);
		check(determinant(a) == -4071);
	},
	scenario("testing an invertible matrix for invertibility") = [](auto check) {
		using namespace rt;
		matrix a {
			6, 4, 4, 4,
			5, 5, 7, 6,
			4, -9, 3, -7,
			9, 1, 7, -6};
		check(determinant(a) == -2120);
		check(is_invertible(a));
	},
	scenario("testing a noninvertible matrix for invertibility") = [](auto check) {
		using namespace rt;
		matrix a {
			-4, 2, -2, -3,
			9, 6, 2, 6,
			0, -5, 1, -5,
			0, 0, 0, 0};
		check(determinant(a) == 0);
		check(!is_invertible(a));
	},
	scenario("calculating the inverse of a matrix") = [](auto check) {
		using namespace rt;
		matrix a {
			-5, 2, 6, -8,
			1, -5, 1, 8,
			7, 7, -6, -7,
			1, -3, 7, 4};
		auto b = inverse(a);
		check(determinant(a) == 532);
		check(a.cofactor(2, 3) == -160);
		check(b[3, 2] == -160.0/532);
		check(a.cofactor(3, 2) == 105);
		check(b[2, 3] == 105.0/532);
		check(b == matrix {
				0.21805, 0.45113, 0.24060, -0.04511,
				-0.80827, -1.45677, -0.44361, 0.52068,
				-0.07895, -0.22368, -0.05263, 0.19737,
				-0.52256, -0.81391, -0.30075, 0.30639});
	},
	scenario("calculating the inverse of another matrix") = [](auto check) {
		using namespace rt;
		matrix a {
			8, -5, 9, 2,
			7, 5, 6, 1,
			-6, 0, 9, 6,
			-3, 0, -9, -4};
		check(inverse(a) == matrix{
				-0.15385, -0.15385, -0.28205, -0.53846,
				-0.07692, 0.12308, 0.02564, 0.03077,
				0.35897, 0.35897, 0.43590, 0.92308,
				-0.69231, -0.69231, -0.76923, -1.92308});
	},
	scenario("calculating the inverse of a third matrix") = [](auto check) {
		using namespace rt;
		matrix a {
			9, 3, 0, 9,
			-5, -2, -6, -3,
			-4, 9, 6, 4,
			-7, 6, 6, 2 };
		check(inverse(a) == matrix{
				-0.04074, -0.07778, 0.14444, -0.22222,
				-0.07778, 0.03333, 0.36667, -0.33333,
				-0.02901, -0.14630, -0.10926, 0.12963,
				0.17778, 0.06667, -0.26667, 0.33333});
	},
	scenario("multiplying a product by its inverse") = [](auto check) {
		using namespace rt;
		matrix a {
			3, -9, 7, 3,
			3, -8, 2, -9,
			-4, 4, 4, 1,
			-6, 5, -1, 1};
		matrix b {
			8, 2, 2, 2,
			3, -1, 7, 0,
			7, 0, 5, 4,
			6, -2, 0, 5};
		matrix c = a * b;
		check(c * inverse(b) == a);
	},
	scenario("multiplying by a translation matrix") = [](auto check) {
		using namespace rt;
		auto transform = translation(5, -3, 2);
		auto p = point(-3, 4, 5);
		check(transform * p == point(2, 1, 7));
	},
	scenario("multiplying by the inverse of a translation matrix") = [](auto check) {
		using namespace rt;
		auto transform = translation(5, -3, 2);
		auto inv = inverse(transform);
		auto p = point(-3, 4, 5);
		check(inv * p == point(-8, 7, 3));
	},
	scenario("translation does not affect vectors") = [](auto check) {
		using namespace rt;
		auto transform = translation(5, -3, 2);
		auto v = vector(-3, 4, 5);
		check(transform * v == v);
	},
	scenario("a scaling matrix applied to a point") = [](auto check) {
		using namespace rt;
		auto transform = scaling(2, 3, 4);
		auto p = point(-4, 6, 8);
		check(transform * p == point(-8, 18, 32));
	},
	scenario("a scaling matrix applied to a vector") = [](auto check) {
		using namespace rt;
		auto transform = scaling(2, 3, 4);
		auto v = vector(-4, 6, 8);
		check(transform * v == vector(-8, 18, 32));
	},
	scenario("multiplying by the inverse of a scaling matrix") = [](auto check) {
		using namespace rt;
		auto transform = scaling(2, 3, 4);
		auto inv = inverse(transform);
		auto v = vector(-4, 6, 8);
		check(inv * v == vector(-2, 2, 2));
	},
	scenario("reflection is scaling by a negative value") = [](auto check) {
		using namespace rt;
		auto transform = scaling(-1, 1, 1);
		auto p = point(2, 3, 4);
		check(transform * p == point(-2, 3, 4));
	},
	scenario("rotating a point around the x axis") = [](auto check) {
		using namespace rt;
		using namespace std::numbers;
		auto p = point(0, 1, 0);
		auto half_q = rotation_x(pi / 4);
		auto full_q = rotation_x(pi / 2);
		check(half_q * p == point(0, sqrt2/2, sqrt2/2));
		check(full_q * p == point(0, 0, 1));
	},
	scenario("rotating a point around the y axis") = [](auto check) {
		using namespace rt;
		using namespace std::numbers;
		auto p = point(0, 0, 1);
		auto half_q = rotation_y(pi / 4);
		auto full_q = rotation_y(pi / 2);
		check(half_q * p == point(sqrt2/2, 0, sqrt2/2));
		check(full_q * p == point(1, 0, 0));
	},
	scenario("rotating a point around the z axis") = [](auto check) {
		using namespace rt;
		using namespace std::numbers;
		auto p = point(0, 1, 0);
		auto half_q = rotation_z(pi / 4);
		auto full_q = rotation_z(pi / 2);
		check(half_q * p == point(-sqrt2/2, sqrt2/2, 0));
		check(full_q * p == point(-1, 0, 0));
	},
	scenario("a shearing transformation moves x in proportion to y") = [](auto check) {
		using namespace rt;
		auto transform = shearing(1, 0, 0, 0, 0, 0);
		auto p = point(2, 3, 4);
		check(transform * p == point(5, 3, 4));
	},
	scenario("a shearing transformation moves x in proportion to z") = [](auto check) {
		using namespace rt;
		auto transform = shearing(0, 1, 0, 0, 0, 0);
		auto p = point(2, 3, 4);
		check(transform * p == point(6, 3, 4));
	},
	scenario("a shearing transformation moves y in proportion to x") = [](auto check) {
		using namespace rt;
		auto transform = shearing(0, 0, 1, 0, 0, 0);
		auto p = point(2, 3, 4);
		check(transform * p == point(2, 5, 4));
	},
	scenario("a shearing transformation moves y in proportion to z") = [](auto check) {
		using namespace rt;
		auto transform = shearing(0, 0, 0, 1, 0, 0);
		auto p = point(2, 3, 4);
		check(transform * p == point(2, 7, 4));
	},
	scenario("a shearing transformation moves y in proportion to x") = [](auto check) {
		using namespace rt;
		auto transform = shearing(0, 0, 0, 0, 1, 0);
		auto p = point(2, 3, 4);
		check(transform * p == point(2, 3, 6));
	},
	scenario("a shearing transformation moves z in proportion to y") = [](auto check) {
		using namespace rt;
		auto transform = shearing(0, 0, 0, 0, 0, 1);
		auto p = point(2, 3, 4);
		check(transform * p == point(2, 3, 7));
	},
	scenario("individual transformation are applied in sequence") = [](auto check) {
		using namespace rt;
		using namespace std::numbers;
		auto p = point(1, 0, 1);
		auto a = rotation_x(pi / 2);
		auto b = scaling(5, 5, 5);
		auto c = translation(10, 5, 7);
		auto p2 = a * p;
		check(p2 == point(1, -1, 0));
		auto p3 = b * p2;
		check(p3 == point(5, -5, 0));
		auto p4 = c * p3;
		check(p4 == point(15, 0, 7));
	},
	scenario("individual transformation are applied in sequence") = [](auto check) {
		using namespace rt;
		using namespace std::numbers;
		auto p = point(1, 0, 1);
		auto a = rotation_x(pi / 2);
		auto b = scaling(5, 5, 5);
		auto c = translation(10, 5, 7);
		auto t = c * b * a;
		check(t * p == point(15, 0, 7));
	},
	scenario("creating and querying a ray") = [](auto check) {
		using namespace rt;
		auto origin = point(1, 2, 3);
		auto direction = vector(4, 5, 6);
		ray r {origin, direction};
		check(r.origin == origin);
		check(r.direction == direction);
	},
	scenario("computing a point from a distance") = [](auto check) {
		using namespace rt;
		ray r {point(2, 3, 4), vector(1, 0, 0)};
		check(position(r, 0) == point(2, 3, 4));
		check(position(r, 1) == point(3, 3, 4));
		check(position(r, -1) == point(1, 3, 4));
		check(position(r, 2.5) == point(4.5, 3, 4));
	},
	scenario("a ray intersects a sphere at two points") = [](auto check) {
		using namespace rt;
		ray r {point(0, 0, -5), vector(0, 0, 1)};
		sphere s {};
		auto xs = intersect(s, r);
		check(xs.size() == 2);
		check(xs[0].t == 4.0);
		check(xs[1].t == 6.0);
	},
	scenario("a ray intersects a sphere at a tangent") = [](auto check) {
		using namespace rt;
		ray r {point(0, 1, -5), vector(0, 0, 1)};
		sphere s {};
		auto xs = intersect(s, r);
		check(xs.size() == 2);
		check(xs[0].t == 5.0);
		check(xs[1].t == 5.0);
	},
	scenario("a ray misses a sphere") = [](auto check) {
		using namespace rt;
		ray r {point(0, 2, -5), vector(0, 0, 1)};
		sphere s {};
		auto xs = intersect(s, r);
		check(xs.size() == 0);
	},
	scenario("a ray originates inside a pehre") = [](auto check) {
		using namespace rt;
		ray r {point(0, 0, 0), vector(0, 0, 1)};
		sphere s {};
		auto xs = intersect(s, r);
		check(xs.size() == 2);
		check(xs[0].t == -1.0);
		check(xs[1].t == 1.0);
	},
	scenario("a sphere is behind a ray") = [](auto check) {
		using namespace rt;
		ray r {point(0, 0, 5), vector(0, 0, 1)};
		sphere s {};
		auto xs = intersect(s, r);
		check(xs.size() == 2);
		check(xs[0].t == -6.0);
		check(xs[1].t == -4.0);
	},
	scenario("an intersection encapsulates t and object") = [](auto check) {
		using namespace rt;
		sphere s {};
		intersection i {3.5, s};
		check(i.t == 3.5);
		check(i.object == s);
	},
	scenario("aggregating intersections") = [](auto check) {
		using namespace rt;
		sphere s {};
		intersection i1 {1, s};
		intersection i2 {2, s};
		auto xs = intersections {i1, i2};
		check(xs.size() == 2);
		check(xs[0].t == 1);
		check(xs[1].t == 2);
	},
	scenario("intersect sets the object on the intersection") = [](auto check) {
		using namespace rt;
		ray r {point(0, 0, -5), vector(0, 0, 1)};
		sphere s {};
		auto xs = intersect(s,r);
		check(xs.size() == 2);
		check(xs[0].object == s);
		check(xs[1].object == s);
	},
	scenario("the hit, when all intersections have a positive t") = [](auto check) {
		using namespace rt;
		sphere s {};
		intersection i1 {1, s};
		intersection i2 {2, s};
		intersections xs = {i2, i1};
		auto i = hit(xs);
		check(i.value() == i1);
	},
	scenario("the hit, when some intersections have a negative t") = [](auto check) {
		using namespace rt;
		sphere s {};
		intersection i1 {-1, s};
		intersection i2 {1, s};
		intersections xs = {i2, i1};
		auto i = hit(xs);
		check(i.value() == i2);
	},
	scenario("the hit, when all intersections have a negative t") = [](auto check) {
		using namespace rt;
		sphere s {};
		intersection i1 {-2, s};
		intersection i2 {-1, s};
		intersections xs = {i2, i1};
		auto i = hit(xs);
		check(!i);
	},
	scenario("the hit is always the lowest nonnegative itnersection") = [](auto check) {
		using namespace rt;
		sphere s {};
		intersection i1 {5, s};
		intersection i2 {7, s};
		intersection i3 {-3, s};
		intersection i4 {2, s};
		intersections xs = {i1, i2, i3, i4};
		auto i = hit(xs);
		check(*i == i4);
	},
	scenario("translating a ray") = [](auto check) {
		using namespace rt;
		ray r {point(1,2,3), vector(0,1,0)};
		auto m {translation(3,4,5)};
		auto r2 {transform(r, m)};
		check(r2.origin == point(4, 6, 8));
		check(r2.direction == vector(0, 1, 0));
	},
	scenario("the hit is always the lowest nonnegative itnersection") = [](auto check) {
		using namespace rt;
		ray r {point(1, 2, 3), vector(0, 1, 0)};
		auto m {scaling(2, 3, 4)};
		auto r2 {transform(r, m)};
		check(r2.origin == point(2, 6, 12));
		check(r2.direction == vector(0, 3, 0));
	},
	scenario("a sphere's default transformation") = [](auto check) {
		using namespace rt;
		sphere s {};
		check(s.transform == identity);
	},
	scenario("changing a sphere's transformation") = [](auto check) {
		using namespace rt;
		sphere s {};
		auto t {translation(2, 3, 4)};
		s.transform = t;
		check(s.transform == t);
	},
	scenario("intersecting a scaled sphere with a ray") = [](auto check) {
		using namespace rt;
		ray r {point(0,0, -5), vector(0, 0, 1)};
		sphere s {};
		s.transform = scaling(2, 2, 2);
		auto xs {intersect(s,r)};
		check(xs.size() == 2);
		check(xs[0].t == 3);
		check(xs[1].t == 7);
	},
	scenario("normal on a sphere at a point on the x axis") = [](auto check) {
		using namespace rt;
		sphere s {};
		auto n {normal_at(s, point(1, 0, 0))};
		check(n == vector(1, 0, 0));
	},
	scenario("normal on a sphere at a point on the y axis") = [](auto check) {
		using namespace rt;
		sphere s {};
		auto n {normal_at(s, point(0, 1, 0))};
		check(n == vector(0, 1, 0));
	},
	scenario("normal on a sphere at a point on the z axis") = [](auto check) {
		using namespace rt;
		sphere s {};
		auto n {normal_at(s, point(0, 0, 1))};
		check(n == vector(0, 0, 1));
	},
	scenario("normal on a sphere at a nonaxial point") = [](auto check) {
		using namespace rt;
		using namespace std::numbers;
		sphere s {};
		auto n {normal_at(s, point(sqrt3/3, sqrt3/3, sqrt3/3))};
		check(n == vector(sqrt3/3, sqrt3/3, sqrt3/3));
	},
	scenario("normal is a normalized vector") = [](auto check) {
		using namespace rt;
		using namespace std::numbers;
		sphere s {};
		auto n {normal_at(s, point(sqrt3/3, sqrt3/3, sqrt3/3))};
		check(n == normalize(n));
	},
	scenario("computing normal on a translated sphere") = [](auto check) {
		using namespace rt;
		sphere s {};
		s.transform = translation(0, 1, 0);
		auto n {normal_at(s, point(0, 1.70711, -0.70711))};
		check(n == vector(0, 0.70711, -0.70711));
	},
	scenario("computing the normal on a transformed sphere") = [](auto check) {
		using namespace rt;
		using namespace std::numbers;
		sphere s {};
		auto m = scaling(1, 0.5, 1) * rotation_z(pi/5);
		s.transform = m;
		//s.transform = m;
		auto n = normal_at(s, point(0, sqrt2/2, -sqrt2/2));
		check(n == vector(0, 0.97014, -0.24254));
	},
	scenario("reflecting a vector approaching 45 degrees") = [](auto check) {
		using namespace rt;
		auto v = vector(1, -1, 0);
		auto n = vector(0, 1, 0);
		auto r = reflect(v, n);
		check(r == vector(1, 1, 0));
	},
	scenario("reflecting a vector off a slanted surface") = [](auto check) {
		using namespace rt;
		using namespace std::numbers;
		auto v = vector(0, -1, 0);
		auto n = vector(sqrt2/2, sqrt2/2, 0);
		auto r = reflect(v, n);
		check(r == vector(1, 0, 0));
	},
	scenario("a point light has a position and intensity") = [](auto check) {
		using namespace rt;
		auto intensity = color(1, 1, 1);
		auto position = point(0, 0, 0);
		point_light light {position, intensity};
		check(light.position == position);
		check(light.intensity == intensity);
	},
	scenario("the default material") = [](auto check) {
		using namespace rt;
		material m {};
		check(m.color == color(1, 1, 1));
		check(m.ambient == 0.1);
		check(m.diffuse == 0.9);
		check(m.specular == 0.9);
		check(m.shininess == 200.0);
	},
	scenario("a sphere has a default material") = [](auto check) {
		using namespace rt;
		sphere s {};
		auto m = s.surface;
		check(m == material {});
	},
	scenario("a sphere may be assigned a material") = [](auto check) {
		using namespace rt;
		sphere s {};
		material m {};
		m.ambient = 1;
		s.surface = m;
		check(s.surface == m);
	},
	scenario("lighting with the eye between the light and the surface") = [](auto check) {
		using namespace rt;
		material m {};
		auto position = point(0, 0, 0);
		auto eyev = vector(0, 0, -1);
		auto normalv = vector(0, 0, -1);
		point_light light {point(0, 0, -10), color(1, 1, 1)};
		auto result = lighting(m, light, position, eyev, normalv);
		check(result == color(1.9, 1.9, 1.9));
	},
	scenario("lighting with the eye between the light and the surface, eye offset 45") = [](auto check) {
		using namespace rt;
		using namespace std::numbers;
		material m {};
		auto position = point(0, 0, 0);

		auto eyev = vector(0, sqrt2/2, -sqrt2/2);
		auto normalv = vector(0, 0, -1);
		point_light light {point(0, 0, -10), color(1, 1, 1)};
		auto result = lighting(m, light, position, eyev, normalv);
		check(result == color(1.0, 1.0, 1.0));
	},
	scenario("lighting with the eye opposite surface, light offset 45") = [](auto check) {
		using namespace rt;
		material m {};
		auto position = point(0, 0, 0);

		auto eyev = vector(0, 0, -1);
		auto normalv = vector(0, 0, -1);
		point_light light {point(0, 10, -10), color(1, 1, 1)};
		auto result = lighting(m, light, position, eyev, normalv);
		check(result == color(0.7364, 0.7364, 0.7364));
	},
	scenario("lightin with eye in the pathe of the reflection vector") = [](auto check) {
		using namespace rt;
		using namespace std::numbers;
		material m {};
		auto position = point(0, 0, 0);

		auto eyev = vector(0, -sqrt2/2, -sqrt2/2);
		auto normalv = vector(0, 0, -1);
		point_light light {point(0, 10, -10), color(1, 1, 1)};
		auto result = lighting(m, light, position, eyev, normalv);
		check(result == color(1.6364, 1.6364, 1.6364));
	},
	scenario("lighting with the light behind the surface") = [](auto check) {
		using namespace rt;
		using namespace std::numbers;
		material m {};
		auto position = point(0, 0, 0);

		auto eyev = vector(0, 0, -1);
		auto normalv = vector(0, 0, -1);
		point_light light {point(0, 0, 10), color(1, 1, 1)};
		auto result = lighting(m, light, position, eyev, normalv);
		check(result == color(0.1, 0.1, 0.1));
	},
	scenario("creating a world") = [](auto check) {
		using namespace rt;
		world w {};
		check(w.objects.size() == 0);
		check(!w.light);
	},
	scenario("the default world") = [](auto check) {
		using namespace rt;
		point_light light {point(-10, 10, -10), color(1, 1, 1)};
		sphere s1;
		s1.surface = {
				.color = color(0.8, 1.0, 0.6),
				.diffuse = 0.7,
				.specular = 0.2};
		sphere s2;
		s2.transform = scaling(0.5, 0.5, 0.5);
		world w = default_world();
		check(w.objects[0] == s1);
		check(w.objects[1] == s2);
		check(w.light && *w.light == light);
	},
	scenario("intersect a world with a ray") = [](auto check) {
		using namespace rt;
		auto w = default_world();
		ray r {point(0, 0, -5), vector(0, 0, 1)};
		auto xs = intersect_world(w, r);
		check(xs.size() == 4);
		check(xs[0].t == 4);
		check(xs[1].t == 4.5);
		check(xs[2].t == 5.5);
		check(xs[3].t == 6);
	},
	scenario("precomputing the state of an intersection") = [](auto check) {
		using namespace rt;
		ray r {point(0, 0, -5), vector(0, 0, 1)};
		sphere shape {};
		auto i = intersection(4, shape);
		auto comps = prepare_computations(i, r);
		check(comps.t == i.t);
		check(comps.object == i.object);
		check(comps.point == point(0, 0, -1));
		check(comps.eyev == vector(0, 0, -1));
		check(comps.normalv == vector(0, 0, -1));
	},
	scenario("the hit, when an intersection occurs on the outside") = [](auto check) {
		using namespace rt;
		ray r {point(0, 0, -5), vector(0, 0, 1)};
		sphere shape {};
		auto i = intersection(4, shape);
		auto comps = prepare_computations(i, r);
		check(comps.inside == false);
	},
	scenario("the hit, when an intersection occurs on the inside") = [](auto check) {
		using namespace rt;
		ray r {point(0, 0, 0), vector(0, 0, 1)};
		sphere shape {};
		auto i = intersection(1, shape);
		auto comps = prepare_computations(i, r);
		check(comps.point == point(0, 0, 1));
		check(comps.eyev == vector(0, 0, -1));
		check(comps.inside == true);
		check(comps.normalv == vector(0, 0, -1));
	},
	scenario("shading an intersection") = [](auto check) {
		using namespace rt;
		auto w = default_world();
		ray r {point(0, 0, -5), vector(0, 0, 1)};
		auto& shape = w.objects[0];
		auto i = intersection(4, shape);
		auto comps = prepare_computations(i, r);
		auto c = shade_hit(w, comps);
		check(c == color(0.38066, 0.47583, 0.2855));
	},
	scenario("shading an intersection from the inside") = [](auto check) {
		using namespace rt;
		auto w = default_world();
		w.light = point_light(point(0, 0.25, 0), color(1, 1, 1));
		ray r {point(0, 0, 0), vector(0, 0, 1)};
		auto& shape = w.objects[1];
		auto i = intersection(0.5, shape);
		auto comps = prepare_computations(i, r);
		auto c = shade_hit(w, comps);
		check(c == color(0.90498, 0.90498, 0.90498));
	},
	scenario("the color when a ray misses") = [](auto check) {
		using namespace rt;
		auto w = default_world();
		ray r {point(0, 0, -5), vector(0, 1, 0)};
		auto c = color_at(w, r);
		check(c == color(0, 0, 0));
	},
	scenario("the color when a ray hits") = [](auto check) {
		using namespace rt;
		auto w = default_world();
		ray r {point(0, 0, -5), vector(0, 0, 1)};
		auto c = color_at(w, r);
		check(c == color(0.38066, 0.47583, 0.2855));
	},
	scenario("the color with an interection behind the ray") = [](auto check) {
		using namespace rt;
		auto w = default_world();
		auto& outer = w.objects[0];
		outer.surface.ambient = 1;
		auto& inner = w.objects[1];
		inner.surface.ambient = 1;
		ray r {point(0, 0, 0.75), vector(0, 0, -1)};
		auto c = color_at(w, r);
		check(c == inner.surface.color);
	},
	scenario("the transformation matrix for the default orientation") = [](auto check) {
		using namespace rt;
		auto from = point(0, 0, 0);
		auto to = point(0, 0, -1);
		auto up = vector(0, 1, 0);
		auto t = view_transform(from, to, up);
		check(t == identity);
	},
	scenario("a view transformation looking in positive z") = [](auto check) {
		using namespace rt;
		auto from = point(0, 0, 0);
		auto to = point(0, 0, 1);
		auto up = vector(0, 1, 0);
		auto t = view_transform(from, to, up);
		check(t == scaling(-1, 1, -1));
	},
	scenario("view transformation moves the world") = [](auto check) {
		using namespace rt;
		auto from = point(0, 0, 8);
		auto to = point(0, 0, 0);
		auto up = vector(0, 1, 0);
		auto t = view_transform(from, to, up);
		check(t == translation(0, 0, -8));
	},
	scenario("an arbitrary view transformation") = [](auto check) {
		using namespace rt;
		auto from = point(1, 3, 2);
		auto to = point(4, -2, 8);
		auto up = vector(1, 1, 0);
		auto t = view_transform(from, to, up);
		check(t == matrix {
				-0.50709, 0.50709, 0.67612, -2.36643,
				0.76772, 0.60609, 0.12122, -2.82843,
				-0.35857, 0.59761, -0.71714, 0.00000,
				0.00000, 0.00000, 0.00000, 1.00000});
	},
	scenario("constructing a camera") = [](auto check) {
		using namespace rt;
		using namespace std::numbers;
		camera c {
			.hsize = 160,
			.vsize = 120,
			.field_of_view = pi/2 };
		check(c.transform == identity);
	},

};

int main() {
	rt_tests.run();
}
