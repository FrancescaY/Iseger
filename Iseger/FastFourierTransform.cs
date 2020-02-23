using System;
using System.Numerics;

/// Classes describing intergal transforms.
namespace Iseger
{
	/// <summary>
	/// Provides an abstraction for a Fast Fourier Transformation (FFT).
	/// This is to understand as a functionality to perform direct 
	/// and inverse Fourier transformations of real functions.
	/// </summary>
	public static class FastFourierTransform
	{
		#region Public Enums
		/// <summary>
		/// How to pad the tail of the function values if the number of the probes is not a power of two:
		/// </summary>
		public enum PaddingMode
		{
			/// <summary>
			/// Pad with zeros
			/// </summary>
			Zeros = 0,

			/// <summary>
			/// Pad with the last value of the series
			/// </summary>
			LastValue = 1
		}
		#endregion

		#region Public Functionality
		/// <summary>
		/// Direct FFT according to Cooley-Tukey routine on double data.
		/// The spectrum is an array of complex numbers with central symmetry: 
		/// \f$ s_i	= s_{N-i} \f$, where \f$ N \f$ is the actual number of telements in the spectrum.
		/// \todo This is a redundancy that could be eliminated in later versions. So long it exists, 
		/// it is essential to copy the values to the right part of the spectrum
		/// after such actions as filtering, so as to make the spectrum symmetrical again.
		/// </summary>
		/// <param name="data">Array of double probe data.</param>
		/// <param name="paddingMode">Mode to pad trailng values.</param>
		/// <returns>Calculated spectrum.</returns>
		public static Complex[] DirectTransform(double[] data, PaddingMode paddingMode = PaddingMode.Zeros)
		{
			Complex[] spectrum = PrepareData(data, paddingMode);

			return CooleyTukey(spectrum, false);
		}

		/// <summary>
		/// Added for convenience.
		/// Direct FFT according to Cooley-Tukey routine on byte data.
		/// </summary>
		/// <param name="data">Array of byte probe data.</param>
		/// <param name="paddingMode">Mode to pad trailng values.</param>
		/// <returns>Calculated spectrum.</returns>
		public static Complex[] DirectTransform(byte[] data, PaddingMode paddingMode = PaddingMode.Zeros)
		{
			Complex[] spectrum = PrepareData(data, paddingMode);

			return CooleyTukey(spectrum, false);
		}

		/// <summary>
		/// Inverse FFT according to Cooley-Tukey routine.
		/// </summary>
		/// <param name="spectrum">Spectrum to transform from.</param>
		/// <returns>Array of data.</returns>
		public static double[] InverseTransform(Complex[] spectrum)
		{
			Complex[] dataComplex	= CooleyTukey(spectrum, true);
			double[] data			= new double[dataComplex.Length];

			for (int i = 0; i < dataComplex.Length; i ++)
			{
				data[i]			= dataComplex[i].Real;
			}

			return data;
		}
		#endregion

		#region Filters
		/// <summary>
		/// 
		/// </summary>
		/// <param name="spectrum"></param>
		/// <param name="from"></param>
		/// <param name="to"></param>
		internal static void Window(ref Complex[] spectrum, int from, int to)
		{
			for (int i = 0; i < from; i++)
			{
				spectrum[i]							= Complex.Zero;
				spectrum[spectrum.Length - i - 1]	= Complex.Zero;
			}

			for (int i = to + 1; i <= spectrum.Length / 2; i++)
			{
				spectrum[i]							= Complex.Zero;
				spectrum[spectrum.Length - i - 1]	= Complex.Zero;
			}
		}

		/// <summary>
		/// Low pass filter. Passes only harmonics with indexes from 0 to \e upperLimit.
		/// </summary>
		/// <param name="spectrum">Array of data to filter.</param>
		/// <param name="upperLimit">Number of the upper harmonic.</param>
		internal static void LowPass(ref Complex[] spectrum, int upperLimit)
		{
			Window(ref spectrum, 0, upperLimit);
		}

		/// <summary>
		/// 
		/// </summary>
		/// <param name="spectrum"></param>
		/// <param name="lowerLimit"></param>
		internal static void HighPass(ref Complex[] spectrum, int lowerLimit)
		{
			Window(ref spectrum, lowerLimit, spectrum.Length - 1);
		}
		#endregion

		#region Private Auxiliary
		/// <summary>
		/// Calculates the next number being a power of two such as it is greater or equal to nSize.
		/// </summary>
		/// <param name="size">Number to find the next power of 2</param>
		/// <returns>The next power of 2 gor the number.</returns>
		private static int NextPowerOf2(int size)
		{
			int figure = 1;

			if (size < 1)
			{
				return 0;	// error
			}

			if (size == 1)
			{
				return 1;
			}

			int nResult = size;

			while (nResult > 1)
			{
				nResult /= 2;
				figure *= 2;
			}

			if (size > figure)
			{
				figure *= 2;
			}

			return figure;
		}

		/// <summary>
		/// Preparation of spectrum according to given double data.
		/// The size of spectrum will be incresed to the next power of two
		/// and the real parts of the spectrum will be initialized with \e arData.
		/// </summary>
		/// <param name="data">Initial data</param>
		/// <param name="paddingMode">Mode to pad trailing elements.</param>
		/// <returns>Prepared complex array.</returns>
		private static Complex[] PrepareData(double[] data, PaddingMode paddingMode)
		{
			int spectrumSize = NextPowerOf2(data.Length);

			Complex[] spectrum = new Complex[spectrumSize];

			Complex trail = Complex.Zero;

			switch (paddingMode)
			{
				case PaddingMode.Zeros:
					trail = new Complex();
					break;

				case PaddingMode.LastValue:
					trail = new Complex(data[data.Length - 1], 0);
					break;
			}

			for (int i = 0; i < spectrum.Length; i++)
			{
				if (i < data.Length)
				{
					spectrum[i] = new Complex(data[i], 0);
				}
				else
				{
					spectrum[i] = trail;
				}
			}

			return spectrum;
		}

		/// <summary>
		/// Added for convenience.
		/// Preparation of spectrum according to given byte data.
		/// The size of spectrum will be incresed to the next power of two and the real parts of the spectrum will be initialized with \e arData.
		/// </summary>
		/// <param name="data">Initial data</param>
		/// <param name="paddingMode">Mode to pad trailing elements.</param>
		/// <returns>Prepared complex array.</returns>
		private static Complex[] PrepareData(byte[] data, PaddingMode paddingMode)
		{
			int spectrumSize = NextPowerOf2(data.Length);

			Complex[] spectrum = new Complex[spectrumSize];

			Complex trail = Complex.Zero;

			switch (paddingMode)
			{
				case PaddingMode.Zeros:
					trail = new Complex();
					break;

				case PaddingMode.LastValue:
					trail = new Complex(data[data.Length - 1], 0);
					break;
			}

			for (int i = 0; i < spectrum.Length; i++)
			{
				if (i < data.Length)
				{
					spectrum[i] = new Complex(data[i], 0);
				}
				else
				{
					spectrum[i] = trail;
				}
			}

			return spectrum;
		}

		/// <summary>
		/// The Cooley-Tukey algorithm.
		/// </summary>
		/// <param name="spectrum">Initially prepared sprectrum.</param>
		/// <param name="isInverseTransform">True if inverse FFT is to be calculated</param>
		/// <returns>Calculated spectrum of the function.</returns>
		private static Complex[] CooleyTukey(Complex[] spectrum, bool isInverseTransform)
		{
			Complex[] arTemp = new Complex[spectrum.Length];

			int powerOf2 = IntLog2(spectrum.Length);

			for (int i = 0; i < spectrum.Length; i++)
			{
				int p = BitReversion(i, powerOf2);

				arTemp[p] = spectrum[i];
			}

			for (int i = 0; i < spectrum.Length; i++)
			{
				spectrum[i] = arTemp[i];
			}

			for (int i = 0; i < spectrum.Length; i += 2)
			{
				Complex sum			= spectrum[i] + spectrum[i + 1];
				Complex difference	= spectrum[i] - spectrum[i + 1];

				spectrum[i]			= sum;
				spectrum[i + 1]		= difference;
			}

			int span				= 4;
			int numberOfSections	= spectrum.Length / span;

			Complex[] w				= new Complex[spectrum.Length];

			while (span <= spectrum.Length)
			{
				double argument		= -2.0 * Math.PI / span;

				if (isInverseTransform)
				{
					argument		= -argument;
				}

				Complex z = Complex.FromPolarCoordinates(1.0, argument);

				w[0] = Complex.One;

				for (int j = 1; j < span; j++)
				{
					w[j] = w[j - 1] * z;
				}

				int spanPrevious	= span / 2;

				for (int s = 0; s < numberOfSections; s++)
				{
					for (int k = 0; k < spanPrevious; k++)
					{
						int targetIndex							= s * span + k;
						int evenIndex							= targetIndex;
						int oddIndex							= targetIndex + spanPrevious;

						Complex even							= spectrum[evenIndex];
						Complex odd								= spectrum[oddIndex];

						spectrum[targetIndex]					= even + w[k] * odd;
						spectrum[targetIndex + spanPrevious]	= even - w[k] * odd;
					}
				}

				span											*= 2;
				numberOfSections								= spectrum.Length / span;
			}

			double factor = (isInverseTransform) ? 0.5 : 2.0 / spectrum.Length;

			for (int i = 0; i < spectrum.Length; i++)
			{
				spectrum[i] *= factor;
			}

			return spectrum;
		}

		/// <summary>
		/// Index with reversed bitd for thr realization of the FFT butterfly.
		/// </summary>
		/// <param name="j">Index for which to calculate the reversed index.</param>
		/// <param name="power">Power of 2 within which the transformation takes place.</param>
		/// <returns>Index with reversed bits.</returns>
		private static int BitReversion(int j, int power)
		{
			int i, j1, j2, k;

			j1	= j;
			k	= 0;

			for (i = 1; i <= power; i++)
			{
				j2	= j1 / 2;
				k	= k * 2 + (j1 - 2 * j2);
				j1	= j2;
			}

			return k;
		}

		/// <summary>
		/// Integer binary logarithm for an integer number.
		/// </summary>
		/// <param name="n">Number for which to calculate logarithm.</param>
		/// <returns>Log_2(n)</returns>
		private static int IntLog2(int n)
		{
			int count = 0;

			while (n > 1)
			{
				n /= 2;
				count++;
			}

			return count;
		}
		#endregion
	}
}