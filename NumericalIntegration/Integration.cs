using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace NumericalIntegration
{


    class Integration
    {
        static void Main(string[] args)
        {
            double x1 = (double)0.0;
            double x2 = (double)1.0;
            double n = 4;
            double tolerance = (double)10.0e-12;

            Func<double, double> f = fFunc;

            double integralValue;

            // COS function evaluation
            integralValue = ConvergenceMachine(f, tolerance, x1, x2, n);
            double result = (double)Math.Exp((double)1.0) + integralValue;
            Console.WriteLine("g({0}) = {1}", 1.0, Math.Round(result, 12).ToString());

            Console.WriteLine();

            //Natural Distribution func evaluation
            Console.WriteLine("Calculating N(0.1)");
            x2 = 0.1;
            f = NDistFunc;
            integralValue = ConvergenceMachine(f, tolerance, x1, x2, n);
            result = (double)0.5 + (1.0 / (Math.Sqrt(2.0 * Math.PI))) * integralValue;
            Console.WriteLine("N({0}) = {1}", 0.1, Math.Round(result, 12).ToString());

            Console.WriteLine();

            Console.WriteLine("Calculating N(0.5)");
            result = 0;
            x2 = 0.5;
            Func<double, double> func = NDistFunc;
            integralValue = ConvergenceMachine(func, tolerance, x1, x2, n);
            result = (double)0.5 + (1.0 / (Math.Sqrt(2.0 * Math.PI))) * integralValue;
            Console.WriteLine("N({0}) = {1}", 0.5, Math.Round(result, 12).ToString());

            Console.WriteLine();

            Console.WriteLine("Calculating N(1.0)");
            x2 = 1.0;
            f = NDistFunc;
            integralValue = ConvergenceMachine(f, tolerance, x1, x2, n);
            result = (double)0.5 + (1.0 / (Math.Sqrt(2.0 * Math.PI))) * integralValue;
            Console.WriteLine("N({0}) = {1}", 1.0, Math.Round(result, 12).ToString());


         
            Console.WriteLine();
            Console.WriteLine("Problem #9: ");
            Console.WriteLine("(i)");
            f = ZeroCurve;
            x1 = 0;
            x2 = 0.5;
            tolerance = 10.0e-6;
            integralValue = ConvergenceMachine(f, tolerance, x1, x2, n);
            result = Math.Exp(-(integralValue));
            Console.WriteLine("6 month discount factor: " + Math.Round(result, 6));

            Console.WriteLine();
            f = ZeroCurve;
            x1 = 0;
            x2 = 1.0;
            tolerance = 10.0e-6;
            integralValue = ConvergenceMachine(f, tolerance, x1, x2, n);
            result = Math.Exp(-(integralValue));
            Console.WriteLine("1 year discount factor: " + Math.Round(result, 6));

            Console.WriteLine();
            f = ZeroCurve;
            x1 = 0;
            x2 = 1.5;
            tolerance = 10.0e-6;
            integralValue = ConvergenceMachine(f, tolerance, x1, x2, n);
            result = Math.Exp(-(integralValue));
            Console.WriteLine("18 months discount factor: " + Math.Round(result, 6));

            Console.WriteLine();
            f = ZeroCurve;
            x1 = 0;
            x2 = 2.0;
            tolerance = 10.0e-8;
            integralValue = ConvergenceMachine(f, tolerance, x1, x2, n);
            result = Math.Exp(-(integralValue));
            Console.WriteLine("2 year discount factor :" + Math.Round(result, 8));


            Console.WriteLine();
            Console.WriteLine("(ii)");

            
            Console.WriteLine();
        

        }



        static double Integrate(Func<double, double> f, double x1, double x2, double n)
        {

            double h = (x2 - x1) / (double)n;
            double integrationSum;
            integrationSum = (f(x1) + f(x2)) / (double)6.0;

            for (int i = 1; i < n; i++)
                integrationSum += f(x1 + i * h) / (double)3.0;

            for (int i = 1; i <= (n); i++)
                integrationSum += (double)2.0 * f(x1 + ((double)i - (double)0.5) * h) / (double)3.0;

            integrationSum *= h;
            return integrationSum;

        }

        static double ConvergenceMachine(Func<double, double> f, double t, double x1, double x2, double n)
        {
            double valPrev = Integrate(f, x1, x2, n);
            Console.WriteLine(n.ToString() + " " + Math.Round(valPrev, 12).ToString());
            n *= 2;
            double valCurrent = Integrate(f, x1, x2, n);
            Console.WriteLine(n.ToString() + " " + Math.Round(valCurrent, 12).ToString());

            while (Math.Abs(valCurrent - valPrev) > t)
            {
                valPrev = valCurrent;
                n *= 2;
                valCurrent = Integrate(f, x1, x2, n);
                Console.WriteLine(n.ToString() + "  " + Math.Round(valCurrent, 12).ToString());
            }

            return valCurrent;
        }

        static double fFunc(double x)
        {
            double r = (double)Math.Sqrt(Math.Cos(Math.Abs((double)x))) / (double)Math.PI;
            return r;
        }

        static double NDistFunc(double x)
        {
            double r = (double)Math.Exp(-(x * x) * 0.5);//  also /2
            return r;
        }

        static double ZeroCurve(double t)
        {
            return (double)0.05 / (1.0 + 2 * Math.Exp(-((1 + t)*(1 + t))));
        }


        //passing n = number of cash flows
        // t = vector of cash flow dates
        // v = vector of cash flows
        // Func<double, double> f =zeroRate 
        //output B = bond price
       /* static double BondCalculator(double n, Math. ) 
        {
            double bondValue = 0;
            for(int i =1; i<(x/2); int++)
            {
                    bondValue +=(r/frequency) *fV*Math.Exp(-(f(x))*x);
            }
            bondValue += (fV+r/frequency) *fV*Math.Exp(-(f(x))*x);
            
        }*/

        static double ZeroRateProblem6(double x)
        {
              return 0.05+0.005*Math.Sqrt(1+x);
        }

     /*   static double[] DiscountFactorCalculator(Func<double, double> f, double[] t, int m)
        {

            double [] cashFlows = new double[m];
            for(int i = 0; i<m; i++)
                cashFlows[i] = Math.Exp(-(f(t[i]))*t[i]);
            return cashFlows;         
        }

        static double[] GetCashFlows(double r, int m, int frequency) {
            double[] cashFlows = new double[m];
            for(int i =0; i<m; i++)
                cashFlows[i] = r/frequency*100;
            cashFlows[m] = 100 + (r / frequency * 100);
            return cashFlows;
        }*/

    }
}
