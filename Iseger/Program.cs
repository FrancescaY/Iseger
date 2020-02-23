using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace Iseger
{
	class Program
	{
		static IsegerInvertor _invertor	= new IsegerInvertor();

		static void Main(string[] args)
		{
			TestIseger();
		}

		static void TestIseger()
		{
			Func<Complex, Complex> F	= ((Complex p) => {return 1 / (1 + p);});
			Func<double, double> f		= ((double t) => {return Math.Exp(-t);});

			double delta		= 0.1;
			int numberOfValues	= 20;

			double[] values	= _invertor.GetValues(F, delta, numberOfValues);

			for (int i = 0; i < values.Length; i++)
			{
				Console.WriteLine($"t = {i * delta}\tIseger = {values[i]}\tExact = {f(delta * i)}");
			}

			Console.ReadKey();
		}
	}
}
