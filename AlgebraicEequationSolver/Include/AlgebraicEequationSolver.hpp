#pragma once

#include <tuple>
#include <complex>
#include <cmath>

#include "QuadraticEquation.hpp"
#include "CubicEquation.hpp"
#include "QuarticEquation.hpp"

namespace AESolver {

	/**
	* @brief 各n次方程式solver実装へのディスパッチ
	*/
	namespace detail {

		/**
		* @brief ax = 0 → x = 0
		*/
		template<typename T>
		constexpr auto solve_nth_degree_equation_impl(const T a) -> typename traits::complexs_tuple<1, T>::type {
			return {};
		}

		template<typename T>
		constexpr auto solve_nth_degree_equation_impl(const T a, const T b) {
			return SolveLinearEquation(a, b);
		}

		template<typename T>
		auto solve_nth_degree_equation_impl(const T a, const T b, const T c) {
			return SolveQuadraticEquation(a, b, c);
		}

		template<typename T>
		auto solve_nth_degree_equation_impl(const T a, const T b, const T c, const T d) {
			return SolveCubicEquation(a, b, c, d);
		}

		template<typename T>
		auto solve_nth_degree_equation_impl(const T a, const T b, const T c, const T d, const T e) {
			return SolveQuarticEquation(a, b, c, d, e);
		}

		/**
		* @brief 入力がおかしい時
		* @detail 4 < n
		* @detail 異なる型が混じっている etc...
		*/
		auto solve_nth_degree_equation_impl(...) -> typename traits::complexs_tuple<1, double>::type = delete;
	}

	/**
	* @brief n次方程式を解く（現在n=4まで）
	* @tparam T 値型、double等を想定、すべて同じ型である必要がある
	* @param a_n 各項の係数（from n to 0、最高次の係数から順に）
	* @return std::tuple<std::complex<T>...> n個の解を格納したtuple（実数解も複素数型として返す）
	*/
	template<typename... T>
	constexpr auto solve_nth_degree_equation(T... a_n) {
		return detail::solve_nth_degree_equation_impl(a_n...);
	}
}