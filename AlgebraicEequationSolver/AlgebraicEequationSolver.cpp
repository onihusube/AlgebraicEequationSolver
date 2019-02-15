#include "pch.h"
#include "CppUnitTest.h"
#include "Include/AlgebraicEequationSolver.hpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace AlgebraicEequationSolver
{
	TEST_CLASS(AlgebraicEequationSolverTest)
	{
	public:
		
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
				auto [x1, x2, x3] = AESolver::SolveCubicEquation(0.0, 3.0, 14.0, 15.0);
				Assert::AreEqual(-1.6666666666667, x1.real(), 1.0E-8);
				Assert::AreEqual(-3.0, x2.real());
				Assert::AreEqual(0.0, x3.real());
			}
		}
	};
}
