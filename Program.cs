using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Threading.Tasks;
using System.Windows.Forms;

using MathNet.Numerics.Interpolation;
using MathNet.Numerics;
using System.IO;
using System.Globalization;
using System.Text.RegularExpressions;

using MathNet.Numerics.Distributions;


namespace IGRF_Module
{
    internal static class Program
    {
        /// <summary>
        /// Главная точка входа для приложения.
        /// </summary>
        [STAThread]
        static void Main()
        {

            // Примеры использования функций и получения с них значений
            //var tt = legendre_poly(1, 0.42);
            //var pp = PnmNorm(2, 3, 0.42);
            //var p = pp.Item1;

            //var gh = MyReadCOEF(new DateTime(2023, 7, 23, 23, 59, 0));
            //var igrf = igrfCalculate(2300, 4900, 3500, new DateTime(2023, 1, 1, 0, 0, 0));

            
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            Application.Run(new Form1());
        }

    }
}
