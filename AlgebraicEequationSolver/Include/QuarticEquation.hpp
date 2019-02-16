#pragma once

#include "CubicEquation.hpp"

namespace AESolver {

	template<typename T>
	auto BiquadraticEquation(const T a, const T c, const T e) -> typename traits::complexs_tuple<4, T>::type {
		auto [x1, x2] = SolveQuadraticEquation(a, c, e);

		auto sqrt_x1 = std::sqrt(x1);
		auto sqrt_x2 = std::sqrt(x2);

		return std::make_tuple(sqrt_x1, -sqrt_x1, sqrt_x2, -sqrt_x2);
	}

	template<typename T>
	auto SolveQuarticEquation(const T p, const T q, const T r) -> typename traits::complexs_tuple<4, T>::type {
		if (q == T(0.0)) {
			//複二次式
			return BiquadraticEquation(T(1.0), p, r);
		}
		else if (r == T(0.0)) {
			constexpr std::tuple<std::complex<T>> zero{};

			auto x = SolveCubicEquation(p, q);

			return std::tuple_cat(std::move(x), zero);
		}

		//tの式の各係数を求める
		//t^3 -pt^2-4rt+(4pr-q^2) = 0

		const auto r4 = T(4.0)*r;
		const auto c =r4*p - q * q;

		std::complex<T> t_c{};

		//tのうち実数のものを選択
		std::tie(t_c, std::ignore, std::ignore) = SolveCubicEquation(-p, -r4, c);
		
		const auto t = t_c.real();

		//m, nを求める

		//t-p
		const auto t_p = t - p;

		if (T(0.0) <= t_p) {
			//m = √(t-p)
			const auto m = std::sqrt(t_p);

			//n = q/2√(t-p)
			const auto n = q / (T(2.0)*m);

			//t/2
			const auto half_t = T(0.5)*t;

			//2つの二次式を解く！
			auto&& y_12 = SolveQuadraticEquation(T(1.0),  m, half_t - n);
			auto&& y_34 = SolveQuadraticEquation(T(1.0), -m, half_t + n);

			return std::tuple_cat(std::move(y_12), std::move(y_34));
		}
		else {
			//負の数の平方根を求めなければならない場合

			//m = √(t-p)
			const std::complex<T> m = { 0.0, std::sqrt(-t_p) };

			//n = q/2√(t-p)
			const auto n = q / (T(2.0)*m);

			//t/2
			const auto half_t = T(0.5)*t;

			//複素係数二次方程式
			auto&& y_12 = SolveQuadraticEquation({ T(1.0) },  m, half_t - n);
			auto&& y_34 = SolveQuadraticEquation({ T(1.0) }, -m, half_t + n);

			return std::tuple_cat(std::move(y_12), std::move(y_34));
		}
	}

	template<typename T>
	auto SolveQuarticEquation(const T A, const T B, const T C, const T D) -> typename traits::complexs_tuple<4, T>::type {
		//(A/4)
		const auto A_d4 = A / T(4.0);
		//(A/4)^2
		const auto A_d4_sq = A_d4 * A_d4;
		//(A/4)^3
		const auto A_d4_cu = A_d4 * A_d4_sq;

		const auto p = T(-6.0) * A_d4_sq + B;
		const auto q = T(2.0) * (T(4.0) * A_d4_cu - B * A_d4) + C;
		const auto r = (T(-3.0) * A_d4_cu - C) * A_d4 + B * A_d4_sq + D;

		auto [y1, y2, y3, y4] = SolveQuarticEquation(p, q, r);

		return std::make_tuple(y1 - A_d4, y2 - A_d4, y3 - A_d4, y4 - A_d4);
	}

	template<typename T>
	auto SolveQuarticEquation(const T a, const T b, const T c, const T d, const T e) -> typename traits::complexs_tuple<4, T>::type {
		if (a == T(0.0)) {
			constexpr std::tuple<std::complex<T>> zero{};

			auto x = SolveCubicEquation(b, c, d, e);
			return std::tuple_cat(std::move(x), zero);
		}
		else if (b == T(0.0) && d == T(0.0)) {
			//複二次式
			return BiquadraticEquation(a, c, e);
		}

		return SolveQuarticEquation(b / a, c / a, d / a, e / a);
	}

}