#pragma once

#include <tuple>
#include <complex>
#include <cmath>

namespace AESolver::traits {

	template<std::size_t N, typename T, typename...Complexs>
	struct complexs_tuple_impl : complexs_tuple_impl<N - 1, T, std::complex<T>, Complexs...>
	{};

	template<typename T, typename...Complexs>
	struct complexs_tuple_impl<0, T, Complexs...> {
		using type = std::tuple<Complexs...>;
	};

	template<std::size_t N, typename T>
	struct complexs_tuple : complexs_tuple_impl<N, T> {};

	template<typename T>
	struct complexs_tuple<0, T>;
}


namespace AESolver {

	template<typename T>
	auto SolveLinearEquation(const T a, const T b) -> std::tuple<T> {
		//ax + b = 0
		//x = -b / a
		return { -a / b };
	}

	template<typename T>
	auto SolveQuadraticEquation(const T a, const T b, const T c) -> typename traits::complexs_tuple<2, T>::type {
		if (a == T(0.0)) {
			constexpr std::tuple<std::complex<T>> zero{};

			auto x = SolveLinearEquation(b, c);

			return std::tuple_cat(std::move(x), zero);
		}

		//D = b^2 - 4ac
		const auto D = b * b - T(4.0) * a * c;

		std::complex<T> x1{}, x2{};

		if (std::fabs(D) < 1.0E-12) {
			//重解
			x1 = x2 = -b / (a + a);
		}
		else if (D < 0.0) {
			//互いに共役な複素数解
			const auto a2 = a + a;
			const auto real = -b / a2;
			const auto imag = std::sqrt(D) / a2;
			x1 = { real,  imag };
			x2 = { real, -imag };
		}
		else {
			//異なる実数解

			//-b - √D
			const auto tmp = -b - std::sqrt(D);
			//2c/(-b - √D)
			x1 = (c + c) / tmp;
			//(-b - √D)/2a
			x2 = tmp / (a + a);
		}

		return std::make_tuple(x1, x2);
	}

}