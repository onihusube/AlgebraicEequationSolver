#pragma once

#include "QuadraticEquation.hpp"

namespace AESolver {

	template<typename T>
	auto SolveCubicEquation(const T p, const T q) -> typename traits::complexs_tuple<3, T>::type {
		//p/3
		const auto p_d3 = p / T(3.0);
		//q/2
		const auto q_half = q / T(2.0);
		//D = (q/2)^2 + (p/3)^3
		const auto D = q_half * q_half + p_d3 * p_d3 * p_d3;

		T x1{};
		std::complex<T> x2{}, x3{};

		if (std::fabs(D) < 1.0E-8) {
			// D == 0 重解

			auto tmp = std::cbrt(q_half);
			x1 = T(-2.0) * tmp;
			x2 = x3 = tmp;
		}
		else if (0.0 < D) {
			// 0 < D １つの実数解と二つの共役な複素数解

			//√D
			const auto D_root = std::sqrt(D);
			//u = v = ∛(-q/2 ± √D)
			const auto u = std::cbrt(-q_half + D_root);
			const auto v = std::cbrt(-q_half - D_root);

			x1 = u + v;
			//-(u + v)/2
			const auto real = T(-0.5) * x1;
			//√3
			constexpr T root3 = T(1.7320508075688772935274463415059);
			//√3(u - v)/2
			const auto imag = root3 * (u - v) * T(0.5);

			x2 = { real,  imag };
			x3 = { real, -imag };
		}
		else {
			// D < 0 異なる３つの実数解

			//2√(-p/3)
			const auto sqrt_p_d3 = T(2.0) * std::sqrt(-p_d3);
			//θ/3
			const auto arg = std::arg(std::complex<T>{ -q_half, std::sqrt(-D) }) / T(3.0);
			//2π/3
			constexpr T pi2d3 = T(2.0 * 3.1415926535897932384626433832795 / 3.0);

			x1 = sqrt_p_d3 * std::cos(arg);
			x2 = sqrt_p_d3 * std::cos(arg + pi2d3);
			x3 = sqrt_p_d3 * std::cos(arg + pi2d3 + pi2d3);
		}

		return std::make_tuple(x1, x2, x3);
	}

	template<typename T>
	auto SolveCubicEquation(const T A, const T B, const T C) -> typename traits::complexs_tuple<3, T>::type {
		//A/3
		auto A_div3 = A / T(3.0);

		const T p = B - A_div3 * A;

		//A^3/27
		auto q_t = A_div3 * A_div3 * A_div3;
		//2A^3/27
		q_t += q_t;

		const T q = q_t - A_div3 * B + C;

		auto [x1, x2, x3] = SolveCubicEquation(p, q);

		//xn = yn - A/3
		return std::make_tuple(x1 - A_div3, x2 - A_div3, x3 - A_div3);
	}


	template<typename T>
	auto SolveCubicEquation(const T a, const T b, const T c, const T d) -> typename traits::complexs_tuple<3, T>::type {
		constexpr std::tuple<std::complex<T>> zero{};

		if (a == 0.0) {
			//bx^2 + cx + d = 0
			//x = (-b ± √(b^2 - 4ac))/2a
			auto x = SolveQuadraticEquation(b, c, d);

			return std::tuple_cat(std::move(x), zero);
		}

		return SolveCubicEquation(b / a, c / a, d / a);
	}
}