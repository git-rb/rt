//rt.cpp - ray tracer module
//
//

module;

#include <stdexcept>
#include <format>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

export module rt;

export namespace rt {
	constexpr double eps = 0.0001;
	bool appx_eq(double a, double b) 
	{
		return std::abs(a-b) < eps;	
	}

	struct tuple {
		double x,y,z,w;
		template <typename T>
		decltype(auto) operator[](this T&& self, int idx) 
		{
			auto&& fwd = std::forward<T>(self);
			switch (idx) {
				case 0: return (fwd.x);
				case 1: return (fwd.y);
				case 2: return (fwd.z);
				case 3: return (fwd.w);
				default: throw std::out_of_range {"tuple index operator out of range"};
			}
		}
		bool is_point() const { return w == 1.0; }
		bool is_vector() const { return w == 0; }

		bool operator==(tuple const& other) const
		{
			return appx_eq(x, other.x) &&
				appx_eq(y, other.y) &&
				appx_eq(z, other.z) &&
				appx_eq(w, other.w);
		}
	};

	std::ostream& operator<<(std::ostream& out, tuple const& t) 
	{
		return out << std::format("[{}, {}, {}, {}]", t.x, t.y, t.z, t.w);
	}

	tuple point(double x, double y, double z) { return tuple{x, y, z, 1.0}; }
	tuple vector(double x, double y, double z) { return tuple{x, y, z, 0.0}; }
	tuple color(double r, double g, double b) { return tuple{r, g, b, 0.0}; }

	double red(tuple const& t) { return t.x; }
	double green(tuple const& t) { return t.y; }
	double blue(tuple const& t) { return t.z; }


	tuple operator+(tuple const& a, tuple const& b)
	{
		return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
	}
	tuple operator-(tuple const& a, tuple const& b)
	{
		return {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
	}
	tuple operator-(tuple const& a)
	{
		return {-a.x, -a.y, -a.z, -a.w};
	}
	tuple operator*(tuple const& t, double f)
	{
		return {t.x * f, t.y * f, t.z * f, t.w * f};
	}
	tuple operator*(tuple const& a, tuple const& b)
	{
		return {a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w};
	}
	tuple operator/(tuple const& t, double f)
	{
		return {t.x / f, t.y / f, t.z / f, t.w / f};
	}
	double magnitude(tuple const& t)
	{
		return std::sqrt(t.x * t.x + t.y * t.y + t.z * t.z + t.w * t.w);
	}
	tuple normalize(tuple const& t)
	{
		return t / magnitude(t);
	}
	double dot(tuple const& a, tuple const& b)
	{
		return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
	}
	tuple cross(tuple const& a, tuple const& b)
	{
		return vector(a.y * b.z - a.z * b.y,
				a.z * b.x - a.x * b.z,
				a.x * b.y - a.y * b.x);
	}


	struct canvas {
	private:
		int width_{};
		int height_ {};
		std::vector<tuple> data_ {static_cast<size_t>(width_ * height_), color(0,0,0)};

	public:
		canvas(int w, int h) : width_{w}, height_{h} {}

		int width() const { return width_; }
		int height() const { return height_; }
		decltype(auto) operator[](this auto&& self, int x, int y)
		{
			return std::forward<decltype(self)>(self).data_[y * self.width_ + x];
		}

		std::string to_ppm()
		{
			std::string ppm = std::format("P3\n{} {}\n255\n", width_, height_);
			auto out = std::back_inserter(ppm);

			for (int y {}; y < height_; ++y)
				for (int x {}; x < width_; ++x) {
					auto c = operator[](x,y);
					std::format_to(out, "{} {} {}\n", 
							std::clamp<int>(red(c) * 256, 0, 255),
							std::clamp<int>(green(c) * 256, 0, 255),
							std::clamp<int>(blue(c) * 256, 0, 255));
				}
			return ppm;
		}
	};

	struct undefined_tag {};
	constexpr undefined_tag undefined;

	template <int size = 4>
	struct matrix {
	private:
		std::array<double, size*size> data_;
	
	public:
		matrix() 
			{ std::ranges::fill(data_.begin(), data_.end(), 0.0); }

		matrix(undefined_tag) { }

		matrix(std::initializer_list<double> init)
		{
			auto n = std::min(init.size(), data_.size());
			auto [_,out_it] = std::ranges::copy_n(init.begin(), n, data_.begin());
			std::ranges::fill(out_it, data_.end(), 0.0);
		}

		template <typename Self>
		decltype(auto) operator[](this Self&& self, int row, int col)
		{
			if (row < 0 || size <= row || col < 0 || size <= col) throw std::out_of_range{"row or column out of range"};
			return std::forward<Self>(self).data_[row * size + col];
		}

		auto data() 
			{ return data_.data(); }
		auto const data() const
			{ return data_.data(); }

		bool operator==(matrix<size> const& o) const
			{ return std::ranges::equal(data_, o.data_, appx_eq); }

		matrix<size-1> submatrix(int row, int col) const
		{
			matrix<size-1> result {undefined};
			auto out_it = result.data();
			for (int r {}; r < size; ++r) {
				if (r == row) continue;
				for (int c {}; c < size; ++c) {
					if (c == col) continue;
					*out_it++ = operator[](r,c);
				}
			}
			return result;
		}

		double minor(int row, int col) const
			{ return determinant(submatrix(row,col)); }

		double cofactor(int row, int col) const
		{ return determinant(submatrix(row, col)) * ((row + col) % 2 ? -1 : 1); }
	};

	template<int size>
	matrix<size> operator*(matrix<size> const& a, matrix<size> const& b)
	{
		matrix result {};
		for (int r {}; r < size; ++r)
			for (int c {}; c < size; ++c)
				for (int x {}; x < size; ++x)
					result[r,c] += a[r,x] * b[x,c];
		return result;
	}

	tuple operator*(matrix<4> const& a, tuple const& b)
	{
		tuple result {};
		for (int r {}; r < 4; ++r)
			for (int c {}; c < 4; ++c)
				result[r] += a[r,c] * b[c];
		return result;
	}
	tuple operator*(tuple const& a, matrix<4> const& b)
		{ return b * a; }

	const matrix identity = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};

	template <int size>
	matrix<size> transpose(matrix<size> const& m)
	{
		matrix result {undefined};
		for (int r {}; r < size; ++r)
			for (int c {}; c < size; ++c)
				result[c,r] = m[r,c];
		return result;
	}

	double determinant(matrix<2> const& m)
		{ return m[0,0] * m[1,1] - m[0,1] * m[1,0]; }

	template <int size>
	double determinant(matrix<size> const& m)
	{
		double sum {};
		for (int col {}; col < size; ++col)
			sum += m[0, col] * m.cofactor(0, col);
		return sum;
	}

	template <int size>
	bool is_invertible(matrix<size> const& m)
		{ return determinant(m) != 0; }

	struct not_invertible : std::runtime_error {
		not_invertible() : runtime_error {"attempted to invert a noninvertible matrix"} {}
	};

	template <int size>
	matrix<size> inverse(matrix<size> const& m)
	{
		matrix<size> ret {undefined};
		auto det = determinant(m);
		if (det == 0) throw not_invertible {};
		for (int r {}; r < size; ++r)
			for (int c {}; c < size; ++c)
				ret[c, r] = m.cofactor(r, c) / det;
		return ret;
	}



	auto translation(double x, double y, double z)
	{
		return matrix {
			1, 0, 0, x,
			0, 1, 0, y,
			0, 0, 1, z,
			0, 0, 0, 1 };
	}
	
	auto scaling(double x, double y, double z)
	{
		return matrix {
			x, 0, 0, 0,
			0, y, 0, 0,
			0, 0, z, 0,
			0, 0, 0, 1 };
	}

	auto rotation_x(double r)
	{
		auto cosr = std::cos(r);
		auto sinr = std::sin(r);
		return matrix {
			1, 0, 0, 0,
			0, cosr, -sinr, 0,
			0, sinr, cosr, 0,
			0, 0, 0, 1 };
	}

	auto rotation_y(double r)
	{
		auto cosr = std::cos(r);
		auto sinr = std::sin(r);
		return matrix {
			cosr, 0, sinr, 0,
			0, 1, 0, 0,
			-sinr, 0, cosr, 0,
			0, 0, 0, 1 };
	}

	auto rotation_z(double r)
	{
		auto cosr = std::cos(r);
		auto sinr = std::sin(r);
		return matrix {
			cosr, -sinr, 0, 0,
			sinr, cosr, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1 };
	}

	auto shearing(double xy, double xz, double yx, double yz, double zx, double zy)
	{
		return matrix {
			1, xy, xz, 0,
			yx, 1, yz, 0,
			zx, zy, 1, 0,
			0, 0, 0, 1 };
	}


	struct ray {
		tuple origin;
		tuple direction;
	};

	tuple position(ray const& r, double t)
		{ return r.origin + r.direction * t; }

	ray transform(ray const& r, matrix<4> const& transformation) 
		{ return {r.origin * transformation, r.direction * transformation}; }

	const tuple black = color(0, 0, 0);
	const tuple white = color(1, 1, 1);

	struct material {
		tuple color = white;
		double ambient = 0.1;
		double diffuse = 0.9;
		double specular = 0.9;
		double shininess = 200.0;

		bool operator==(material const&) const = default;
	};

	struct sphere {
		matrix<4> transform = identity;
		material surface {};

		bool operator==(sphere const&) const = default;
	};

	struct intersection {
		double t;
		sphere object;

		bool operator==(intersection const&) const = default;
	};

	using intersections = std::vector<intersection>;

	intersections intersect(sphere const& s, ray const& ray)
	{
		auto r {transform(ray, inverse(s.transform))};
		auto sphere_to_ray = r.origin - point(0, 0, 0);
		auto a = dot(r.direction, r.direction);
		auto b = 2 * dot(r.direction, sphere_to_ray);
		auto c = dot(sphere_to_ray, sphere_to_ray) - 1;
		auto discriminant = b * b - 4 * a * c;
		if (discriminant < 0) return {};
		auto t1 = (-b - std::sqrt(discriminant)) / (2 * a);
		auto t2 = (-b + std::sqrt(discriminant)) / (2 * a);
		return {{t1, s}, {t2, s}};
	}

	auto normal_at(sphere const& s, tuple const& p)
	{
		auto object_point = inverse(s.transform) * p;
		auto object_nl = object_point - point(0, 0, 0);
		auto world_nl = transpose(inverse(s.transform)) * object_nl;
		world_nl.w = 0;
		return normalize(world_nl);
	}

	std::optional<intersection> hit(intersections is)
	{
		std::ranges::sort(is, [](auto&& a, auto&& b){return a.t < b.t; });
		auto it = is.begin();
		while (it != is.end() && it->t < 0) ++it;
		if (it == is.end()) return {};
		return {*it};
	}

	auto reflect(tuple const& in, tuple const& normal)
		{ return in - normal * 2 * dot(in , normal); }

	struct point_light {
		tuple position = point(0, 0, 0);
		tuple intensity = color(1, 1, 1);
		bool operator==(point_light const&) const = default;
	};

	auto lighting(material const& mat, point_light const& light, tuple const& point,
			tuple const& eyev, tuple const& normalv)
	{
		auto effective_color = mat.color * light.intensity;
		auto lightv = normalize(light.position - point);
		auto ambient = effective_color * mat.ambient;
		auto light_dot_normal = dot(lightv, normalv);
		tuple diffuse;
		tuple specular;

		if (light_dot_normal < 0) {
			diffuse = black;
			specular = black;
		} else
			diffuse = effective_color * mat.diffuse * light_dot_normal;

		auto reflectv = reflect(-lightv, normalv);
		auto reflect_dot_eye = dot(reflectv, eyev);

		if (reflect_dot_eye <= 0)
			specular = black;
		else {
			auto factor = std::pow(reflect_dot_eye, mat.shininess);
			specular = light.intensity * mat.specular * factor;
		}
		return ambient + diffuse + specular;
	}

	struct world {
		std::vector<sphere> objects {};
		std::optional<point_light> light {};
	};

	auto intersect_world(world const& w, ray const& r)
	{
		std::vector<intersection> is {};
		for (auto const& obj : w.objects)
			for (auto&& i : intersect(obj, r))
				is.push_back(std::move(i));
		std::ranges::sort(is, [](auto&& a, auto&& b){ return a.t < b.t; });
		return is;
	}

	struct computations {
		double t {};
		sphere object {};
		tuple point {};
		tuple eyev {};
		tuple normalv {};
		bool inside {};
	};

	auto prepare_computations(intersection const& inter, ray const& r)
	{
		computations comps {};
		comps.t = inter.t;
		comps.object = inter.object;
		comps.point = position(r, comps.t);
		comps.eyev = -r.direction;
		comps.normalv = normal_at(comps.object, comps.point);
		if(dot(comps.normalv, comps.eyev) < 0) {
			comps.inside = true;
			comps.normalv = -comps.normalv;
		} else {
			comps.inside = false;
		}
		return comps;
	}

	auto shade_hit(world const& w, computations const& comps)
	{
		if (!w.light) return black;
		return lighting(comps.object.surface, *w.light, comps.point,
				comps.eyev, comps.normalv);
	}

	auto color_at(world const& w, ray const& r)
	{
		auto op_hit = hit(intersect_world(w,r));
		if (!op_hit) return black;
		return shade_hit(w, prepare_computations(*op_hit, r));
	}

	auto view_transform(tuple const& from, tuple const& to, tuple const& up)
	{
		auto forward = normalize(to - from);
		auto upn = normalize(up);
		auto left = cross(forward, upn);
		auto true_up = cross(left, forward);
		matrix orientation {
			left.x, left.y, left.z, 0,
			true_up.x, true_up.y, true_up.z, 0,
			-forward.x, -forward.y, -forward.z, 0,
			0, 0, 0, 1};
		return orientation * translation(-from.x, -from.y, -from.z);
	}

	struct camera {
		double hsize {};
		double vsize {};
		double field_of_view {};
		matrix<4> transform = identity;
	};

}
