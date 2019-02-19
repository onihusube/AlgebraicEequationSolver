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
	constexpr auto SolveLinearEquation(const T a, const T b) -> typename traits::complexs_tuple<1, T>::type {
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
			const auto imag = std::sqrt(-D) / a2;
			x1 = { real,  imag };
			x2 = { real, -imag };
		}
		else {
			//異なる実数解

			//桁落ち対策
			if (T(0.0) < b) {
				//-b - √D
				const auto tmp = -b - std::sqrt(D);
				//2c/(-b - √D)
				x1 = (c + c) / tmp;
			}
			else {
				//(-b + √D)
				const auto tmp = -b + std::sqrt(D);
				//(-b + √D)/2a
				x1 = tmp / (a + a);
			}

			//解と係数の関係（x1・x2 = c/a）より
			x2 = (c / a) / x1;
		}

		return std::make_tuple(x1, x2);
	}

	template<typename T>
	auto SolveQuadraticEquation(const std::complex<T> a, const std::complex<T> b, const std::complex<T> c) -> typename traits::complexs_tuple<2, T>::type {
		//D = b^2 - 4ac
		const auto D = b * b - T(4.0)*a*c;
		//√D
		const auto root_D = std::sqrt(D);

		const auto a2 = a + a;

		constexpr T zero(0.0);

		//桁落ち対策
		if (zero < b.real() && zero < b.imag()) {
			if (root_D.real() <= zero && root_D.imag() <= zero) {
				//-b - √D
				const auto tmp = -b - root_D;
				//2c/(-b - √D)
				auto x1 = (c + c) / tmp;
				auto x2 = (c / a) / x1;
				return std::make_tuple(std::move(x1), std::move(x2));
			}
		}
		else if(b.real() <= zero && b.imag() <= zero) {
			if (zero < root_D.real() && zero < root_D.imag()) {
				const auto tmp = -b + root_D;
				//(-b + √D)/2a
				auto x1 = tmp / a2;
				auto x2 = (c / a) / x1;
				return std::make_tuple(std::move(x1), std::move(x2));
			}
		}

		//桁落ち対策を諦めた

		//(-b + √D)/2a
		auto x1 = (-b + root_D) / a2;
		//(-b - √D)/2a
		auto x2 = (-b - root_D) / a2;

		return std::make_tuple(std::move(x1), std::move(x2));
	}
}