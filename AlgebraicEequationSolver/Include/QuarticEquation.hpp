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


	}

	template<typename T>
	auto SolveQuarticEquation(const T A, const T B, const T C, const T D) -> typename traits::complexs_tuple<4, T>::type {
		//(A/4)
		const auto A_d4 = A / T(4.0);
		//(A/4)^2
		const auto A_d4_sq = A_d4 * A_d4;
		//(A/4)^3
		const auto A_d4_cu = A_d4 * A_d4_sq;

		const auto P = T(-6.0) * A_d4 + B;
		const auto q = T(2.0) * (T(4.0) * A_d4_cu - B * A_d4) + C;
		const auto r = T(-3.0) * A_d4_cu + B * A_d4_sq - C * A_d4 + D;

		auto [y1, y2, y3, y4] = SolveQuarticEquation(p, q, r);

		return std::make_tuple(y1 - A_d4, y2 - A_d4, y3 - A_d4, y4 - A_d4)
	}

	template<typename T>
	auto SolveQuarticEquation(const T a, const T b, const T c, const T d, const T e) -> typename traits::complexs_tuple<4, T>::type {
		if (a == T(0.0)) {
			constexpr std::tuple<std::complex<T>> zero{};

			auto x = SolveCubicEquation(b, c, d, e);
			return std::tuple_cat(std::move(x), zero);
		}
		else if (b == T(0.0) && d = T(0.0)) {
			//複二次式
			return BiquadraticEquation(a, c, e);
		}

		return SolveQuarticEquation(b / a, c / a, d / a, e / a);
	}

}