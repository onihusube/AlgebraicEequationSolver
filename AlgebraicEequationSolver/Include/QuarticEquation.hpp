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
		//t^3 -(1/2)pt^2-rt-(1/8)q^2+(1/2)pr = 0
		//t^3 + at^2 + bt + c = 0として名前付け

		const auto half_p = T(0.5) * p;
		const auto a = -half_p;
		const auto b = -r;
		const auto c = T(0.125) * q * q + half_p * r;

		std::complex<T> t_c{};

		//tのうち実数のものを選択
		std::tie(t_c, std::ignore, std::ignore) = SolveCubicEquation(a, b, c);

		const auto t = t_c.real();

		//m, nを求める

		//2t-p
		const auto t2mp = t + t - p;

		//m = √(2t-p)
		const auto m = std::sqrt(T(2.0) * t - p);

		//n = -q(√(2t-p))/2(2t-p)
		const auto n = -q * m / (T(2.0) * t2mp);

		//2つの二次式を解く！

		auto&& y_12 = SolveQuadraticEquation(T(1.0), m, t + n);
		auto&& y_34 = SolveQuadraticEquation(T(1.0), -m, t - n);

		return std::tuple_cat(std::move(y_12), std::move(y_34));
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