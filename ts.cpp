//ts.cpp - testing suite

module;

#include <string_view>
#include <memory>
#include <source_location>
#include <print>
#include <vector>

export module ts;

export namespace ts {

	namespace printer {
		inline void fail_loc(std::source_location loc)
		{
			std::print("{}:{}:{}:", loc.file_name(), loc.line(), loc.column());
		}
		inline void failure(std::string_view desc, std::source_location loc, std::string_view check_desc)
		{
			fail_loc(loc);
			std::print(" assertion failed in test \"{}\"", desc);
			if (!empty(check_desc))
				std::print(": {}", check_desc);
			std::println();
		}
		inline void summary(auto test_n, auto check_n, auto fail_n)
		{
			if (fail_n <= 0) std::print("** PASS ** ");
			else std::print("!! FAIL !! ");
			std::println("{} checks in {} tests passed. {} checks failed.", check_n, test_n, fail_n);
		}
	}

	struct test_stats {
		int checked{};
		int failed{};
	};

	struct test_env {
		test_env(std::string_view desc, test_stats& stats)
			:m {desc}, m_stats {stats} {}

		struct assertion {
			assertion(bool bool_val, std::source_location location = std::source_location::current()) 
				: val {bool_val}, loc {location} {}

			operator bool() const { return val; }
			std::source_location location() const { return loc; }

			private:
				bool val;
				std::source_location loc;
		};

		void operator()(assertion val, std::string_view check_desc = "")
		{
			++m_stats.checked;
			if (val) return;
			++m_stats.failed;
			printer::failure(m.test_desc, val.location(), check_desc);
		}
		void check(assertion val, std::string_view check_desc = "")
		{
			operator()(val, check_desc);
		}
		void that(assertion val, std::string_view check_desc = "") 
		{ 
			operator()(val, check_desc); 
		}

	private:
		struct {
			std::string_view test_desc;		
		} m;
		test_stats& m_stats;
	};

	//
	// type-erased test class
	//
	struct test {/*{{{*/

		struct erased {/*{{{*/
			virtual void operator()(test_env env) const = 0;
			virtual std::unique_ptr<erased> clone() const = 0;
			virtual ~erased() {}
		};

		template <typename F>
		struct concrete : erased {
			concrete(F func) : func {func} {}

			void operator()(test_env env) const override { func(env); }

			std::unique_ptr<erased> clone() const override 
			{ 
				return std::unique_ptr<erased> {new concrete<F>{*this}}; 
			}
			F func;
		};/*}}}*/

		std::string_view desc;
		std::unique_ptr<erased> pfunc;

		template <typename F>
		test(std::string_view description, F test_func)
			: desc {description},
				pfunc {new concrete<F>{test_func}} 
		{ }

		test(test const& t)
			: desc {t.desc}, 
				pfunc {t.pfunc->clone()}
		{ }

		test(test&& t) = default;

		std::string_view description() const { return desc; }

		void operator()(test_env env) const { if (pfunc) (*pfunc)(env); }
	};/*}}}*/

	struct test_case { //{{{
		std::string_view desc;

		test_case(std::string_view description) 
			: desc {description}
		{ }

		template <typename F>
		test operator=(F&& func) const 
		{ 
			return test {desc, std::forward<F>(func)}; 
		}
	}; //}}}

	using scenario = test_case;
	using expect_that = test_case;

	struct suite {/*{{{*/
		std::vector<test> tests;

		suite(std::initializer_list<test> init)
			: tests {init} {}

		void run() const 
		{
			int test_n {};
			test_stats stats;
			for (auto&& test : tests) {
				test_env env {test.description(), stats};
				test(env);
				++test_n;
			}
			printer::summary(test_n, stats.checked, stats.failed);
		}

	};/*}}}*/

}
