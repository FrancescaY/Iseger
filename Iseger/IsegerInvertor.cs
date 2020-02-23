using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace Iseger
{
	/// <summary>
	/// Source: Peter den Iseger. Numerical Transform Inversion Using Gaussian Quadrature.
	/// Probability in the Engineering and Informational Sciences, 2006, v. 20, pp 1-44 doi:10.1017/S0269964806060013.
	/// </summary>
	/// <seealso cref="PikkaMath.IntegralTransforms.Laplace.IComplexInvertor" />
	public class IsegerInvertor
	{
		#region Private Data		
		/// <summary>
		/// Iseger's coefficients.
		/// Borrowed from: https://www.cs.hs-rm.de/~weber/lapinv/deniseger/deniseger.htm
		/// </summary>
		private static Dictionary<int, (double Alpha, double Lambda)[]> IsegerCoefficients	= new Dictionary<int, (double Alpha, double Lambda)[]>()
		{
			{
				16,
				new (double Alpha, double Lambda)[]
				{
					(1.00000000000000, 0),
					(1.00000000000004, 6.28318530717958),
					(1.00000015116847, 12.5663706962589),
					(1.00081841700481, 18.8502914166954),
					(1.09580332705189, 25.2872172156717),
					(2.00687652338724, 34.2969716635260),
					(5.94277512934943, 56.1725527716607),
					(54.9537264520382, 170.533131190126),
				}
			},

			{
				32,
				new (double Alpha, double Lambda)[]
				{
					 (1.00000000000000, 0),
					 (1.00000000000000, 6.28318530717958),
					 (1.00000000000000, 12.5663706143592),
					 (1.00000000000000, 18.8495559215388),
					 (1.00000000000000, 25.1327412287184),
					 (1.00000000000895, 31.4159265359035),
					 (1.00000004815464, 37.6991118820067),
					 (1.00003440685547, 43.9823334683971),
					 (1.00420404867308, 50.2716029125234),
					 (1.09319461846681, 56.7584358919044),
					 (1.51528642466058, 64.7269529917882),
					 (2.41320766467140, 76.7783110023797),
					 (4.16688127092229, 96.7780294888711),
					 (8.37770013129610, 133.997553190014),
					 (23.6054680083019, 222.527562038705),
					 (213.824023377988, 669.650134867713),
				}
			},

			{
				48,
				new (double Alpha, double Lambda)[]
				{
					(1.00000000000000,  0),
					(1.00000000000000,  6.28318530717957),
					(1.00000000000000,  12.5663706143592),
					(1.00000000000000,  18.8495559215388),
					(1.00000000000000,  25.1327412287183),
					(1.00000000000000,  31.4159265358979),
					(1.00000000000000,  37.6991118430775),
					(1.00000000000000,  43.9822971502571),
					(1.00000000000000,  50.2654824574367),
					(1.00000000000234,  56.5486677646182),
					(1.00000000319553,  62.8318530747628),
					(1.00000128757818,  69.1150398188909),
					(1.00016604436873,  75.3984537709689),
					(1.00682731991922,  81.6938697567735),
					(1.08409730759702,  88.1889420301504),
					(1.36319173228680,  95.7546784637379),
					(1.85773538601497,  105.767553649199),
					(2.59022367414073,  119.58751936774),
					(3.73141804564276,  139.158762677521),
					(5.69232680539143,  168.156165377339),
					(9.54600616545647,  214.521886792255),
					(18.8912132110256,  298.972429369901),
					(52.7884611477405,  497.542914576338),
					(476.4483318696360, 1494.71066227687),
				}
			}
		};
		#endregion

		/// <summary>
		/// Computes the value of the Laplace original for a given Laplace image and value of the (time) argument.
		/// </summary>
		/// <param name="LaplaceImage">The Laplace image.</param>
		/// <param name="t">The time value.</param>
		/// <param name="parameters">Optional parameters of an overloading inversion method.</param>
		/// <returns>
		/// The numerical value of the original for the argument.
		/// </returns>
		/// <exception cref="NotImplementedException"></exception>
		public double GetValue(Func<Complex, Complex> LaplaceImage, double t, params object[] parameters)
		{
			throw new NotImplementedException();
		}

		/// <summary>
		/// Computes the values of the Laplace original for a given Laplace image and a sequence of values of the argument with given step.
		/// </summary>
		/// <param name="LaplaceImage">The Laplace image.</param>
		/// <param name="deltaT">The argument step.</param>
		/// <param name="numberOfOutputValues">The number of output values of the result, should be a power of 2.</param>
		/// <param name="parameters">
		///		If the number of parameters is either 1 or 2, 
		///		- the first parameter is interpreted as the critical abscissa for the Laplace image 
		///			which is the value of an abscissa in the complex plane such that c is greater than the real part 
		///			of all singularities of the  Laplace image and it is bounded on the line (https://en.wikipedia.org/wiki/Inverse_Laplace_transform);
		///		- the second one as the degree of Gauss quadrature (supported values are 16, 32, 48).
		///	</param>
		/// <returns>
		///		Array of values of the original function in points k * deltaT, k = 0, ..., m.
		///	</returns>
		public double[] GetValues
									(
										Func<Complex, Complex>	LaplaceImage, 
										double					deltaT, 
										int						numberOfOutputValues, 
										params					object[] parameters
									)
		{
			if (deltaT <= 0)
			{
				throw new ArgumentException($"The value of delta is invalid: {deltaT}.");
			}

			if (numberOfOutputValues < 2)
			{
				throw new ArgumentException($"The number of output values is invalid: {numberOfOutputValues} (must be >= 2).");
			}

			double criticalAbscissa = 0;

			if (parameters.Length >= 1 && parameters[0] is double)
			{
				criticalAbscissa	= (double)parameters[0];
			}

			int						quadratureDegree = 16;

			if (parameters.Length >= 2 && parameters[1] is int)
			{
				quadratureDegree	= (int)parameters[1];
			}

			if (!IsegerCoefficients.ContainsKey(quadratureDegree))
			{
				throw new KeyNotFoundException($"The number of quadrature nodes {quadratureDegree} is not supported. Must be 16, 32, or 48.");
			}

			int m					= numberOfOutputValues;
			int mm					= 2;

			while (mm <= numberOfOutputValues)
			{
				mm					*= 2;
			}

			if (mm < numberOfOutputValues)
			{
				mm					*= 2;
			}

			numberOfOutputValues	= mm;
			int m2					= 8 * numberOfOutputValues;
			double b				= 44.0 / m2;

			(double Beta, double Lambda)[] coeffs	= IsegerCoefficients[quadratureDegree];

			double[] y	= new double[m + 1];

			for (int k = 1; k <= m; k++)
			{
				y[k]				= 0;
			}

			Complex[] imageValues	= new Complex[m2 + 1];

			for (int k = 0; k <= m2; k++)
			{
				double sum			= 0;

				for (int j = 0; j < quadratureDegree / 2; j++)
				{
					Complex z		= b + Complex.ImaginaryOne * (coeffs[j].Lambda + 2.0 * Math.PI * (double)k / m2);
					z				= criticalAbscissa + z / deltaT;
					sum				+= coeffs[j].Beta * LaplaceImage(z).Real;
				}

				imageValues[k]		= 2.0 * sum / deltaT;
			}

			imageValues[0]			= (imageValues[0] + imageValues[m2]) / 2.0;

			double[] inverseFft		= FastFourierTransform.InverseTransform(imageValues.Take(m2).ToArray());

			int m4					= m2 / 4;

			double[] result	= new double[numberOfOutputValues];

			for (int j = 0; j < result.Length; j++)
			{
				double expArg		= b * j;

				if (criticalAbscissa > 0)
				{
					expArg			+= criticalAbscissa *(j * deltaT);
				}

				result[j]			= inverseFft[j] * Math.Exp(expArg) / m4;
			}

			return result;
		}
	}
}
