#include "pch.h"
#include "CppUnitTest.h"
#include "Include/AlgebraicEequationSolver.hpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace AlgebraicEequationSolver
{
	TEST_CLASS(AlgebraicEequationSolverTest)
	{
	public:

		TEST_METHOD(QuadraticEquationTest)
		{
			using namespace std::complex_literals;

			using zd = std::complex<double>;

			{
				//複素係数二次方程式
				auto[x1, x2] = AESolver::SolveQuadraticEquation(zd{ 1.0,1.0 }, -zd{ 2,-1.0 }, zd{ -1.0, 2.0 });

				Assert::IsTrue( 1.0 == x1.real());
				Assert::IsTrue(-2.0 == x1.imag());
				Assert::IsTrue(-0.5 == x2.real());
				Assert::IsTrue( 0.5 == x2.imag());
			}

			{
				//桁落ち対策のチェック
				auto[x1, x2] = AESolver::solve_nth_degree_equation(0.2876, -53.14, 9.872);

				Assert::AreEqual(184.584554016612735, x1.real(), 1.0E-13);
				Assert::AreEqual(0.185960587003398633, x2.real(), 1.0E-16);
			}
		}
		
		TEST_METHOD(CubicEquationTest)
		{
			{
				auto [x1, x2, x3] = AESolver::SolveCubicEquation(1.0f, 6.0f, 9.0f, 4.0f);
				Assert::IsTrue(-4.0f == x1.real());
				Assert::IsTrue(-1.0f == x2.real());
				Assert::IsTrue(-1.0f == x3.real());
			}

			{
				auto [x1, x2, x3] = AESolver::SolveCubicEquation(2.0, -3.0, -5.0, 6.0);
				Assert::IsTrue(2.0 == x1.real());
				Assert::AreEqual(-1.5, x2.real(), 1.0E-8);
				Assert::AreEqual(1.0, x3.real(), 1.0E-8);
			}

			{
				auto [x1, x2, x3] = AESolver::SolveCubicEquation(1.0, -4.0, 5.0, 6.0);
				Assert::AreEqual(-0.7161886589931, x1.real(), 1.0E-8);

				Assert::AreEqual(2.3580943294966, x2.real(), 1.0E-8);
				Assert::AreEqual(1.6784135260534, x2.imag(), 1.0E-8);

				Assert::AreEqual(2.3580943294966, x3.real(), 1.0E-8);
				Assert::AreEqual(-1.6784135260534, x3.imag(), 1.0E-8);
			}

			{
				auto [x1, x2, x3] = AESolver::solve_nth_degree_equation(0.0, 3.0, 14.0, 15.0);
				Assert::AreEqual(-1.6666666666667, x1.real(), 1.0E-8);
				Assert::AreEqual(-3.0, x2.real());
				Assert::AreEqual(0.0, x3.real());
			}
		}


		TEST_METHOD(QuarticEquationTest) {

			{
				//複二次式
				auto [x1, x2, x3, x4] = AESolver::SolveQuarticEquation(1.0, -12.0, 49.0, -78.0, 40.0);

				Assert::IsTrue(5.0 == x1.real());
				Assert::IsTrue(1.0 == x2.real());
				Assert::IsTrue(4.0 == x3.real());
				Assert::IsTrue(2.0 == x4.real());
			}

			{
				auto[x1, x2, x3, x4] = AESolver::SolveQuarticEquation(1.0, 4.0, -21.0, 10.0, -1.0);

				Assert::AreEqual(0.1400549446403, x1.real(), 1.0E-8);
				Assert::AreEqual(-7.1400549446403, x2.real(), 1.0E-8);
				Assert::AreEqual(2.6180339887499, x3.real(), 1.0E-8);
				Assert::AreEqual(0.3819660112501, x4.real(), 1.0E-8);
			}

			{
				auto[x1, x2, x3, x4] = AESolver::SolveQuarticEquation(19.0, 4.0, 58.0, 32.0, 7.0);

				Assert::AreEqual(-0.2760718664709, x1.real(), 1.0E-8);
				Assert::AreEqual(-0.2760718664709, x2.real(), 1.0E-8);
				Assert::AreEqual(0.1708087085762, x3.real(), 1.0E-8);
				Assert::AreEqual(0.1708087085762, x4.real(), 1.0E-8);

				Assert::AreEqual(0.204312122584, x1.imag(), 1.0E-8);
				Assert::AreEqual(-0.204312122584, x2.imag(), 1.0E-8);
				Assert::AreEqual(1.759010733644, x3.imag(), 1.0E-8);
				Assert::AreEqual(-1.759010733644, x4.imag(), 1.0E-8);
			}

			{
				auto[x1, x2, x3, x4] = AESolver::solve_nth_degree_equation(1.0, 1.0, 1.0, 1.0, 0.0);

				Assert::AreEqual(0.0, x1.real(), 1.0E-8);
				Assert::AreEqual(-1.0, x2.real(), 1.0E-8);
				Assert::AreEqual(0.0, x3.real(), 1.0E-8);
				Assert::AreEqual(0.0, x4.real(), 1.0E-8);

				Assert::AreEqual(1.0, x3.imag(), 1.0E-8);
				Assert::AreEqual(-1.0, x4.imag(), 1.0E-8);
			}
		}
	};
}
