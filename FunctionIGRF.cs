using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IGRF_Module
{
    internal class FunctionIGRF
    {



        /// <summary>
        /// Функция перевода декартовых координат [X, Y, Z] в сферические [склонение, наклонение, горизонтальная_составляющая, тотальная_сила] 
        /// </summary>
        /// <param name="x">в км</param>
        /// <param name="y">в км</param>
        /// <param name="z">в км</param>
        /// <returns>{ dec [рад], inc [рад], hoz, eff [км] }</returns>
        public static double[] xyz2dhif(double x, double y, double z)
        {
            //double FACT = 180 / Math.PI;

            double hsq = x * x + y * y;
            double hoz = Math.Sqrt(hsq);
            double eff = Math.Sqrt(hsq + z * z);
            double dec = Math.Atan2(y, x);
            double inc = Math.Atan2(z, hoz);
            //return new double[] { dec * FACT, inc * FACT, hoz, eff };
            return new double[] { dec, inc, hoz, eff };
        }
        /// <summary>
        /// Функция вычисления присоединенных полиномов Лежандра нормированных по Шмидту и их производных до Pnm. Вводится n и вычисляется до n=m.
        /// Также вводится угол в радианах.
        /// Пример: при n=2, будут вычислены полиномы и их производные n=1 m=0, n=1 m=1, n=2 m=0, n=2 m=1, n=2 m=2
        /// </summary>
        /// <param name="nmax">Порядок полинома</param>
        /// <param name="theta">Угол в радианах (предварительно вычисленное значение cos(alpha))</param>
        /// <returns>Трехмерный массив с полиномами и их производными</returns>
        public static double[,,] legendre_poly(int nmax, double theta)
        {
            double costh = theta;
            double sinth = Math.Sqrt(1 - costh * costh);

            double[,,] Pnm = new double[nmax + 1, nmax + 2, 1];
            Pnm[0, 0, 0] = 1;
            Pnm[1, 1, 0] = sinth;

            double[] rootn = new double[2 * nmax * nmax + 1];
            for (int i = 0; i < rootn.Length; i++)
            {
                rootn[i] = Math.Sqrt(i);
            }

            for (int m = 0; m < nmax; m++)
            {
                double Pnm_tmp = rootn[m + m + 1] * Pnm[m, m, 0];
                Pnm[m + 1, m, 0] = costh * Pnm_tmp;

                if (m > 0)
                {
                    Pnm[m + 1, m + 1, 0] = sinth * Pnm_tmp / rootn[m + m + 2];
                }

                for (int n = m + 2; n < nmax + 1; n++)
                {
                    int d = n * n - m * m;
                    int e = n + n - 1;
                    Pnm[n, m, 0] = ((e * costh * Pnm[n - 1, m, 0] - rootn[d - e] * Pnm[n - 2, m, 0]) / rootn[d]);
                }
            }

            Pnm[0, 2, 0] = -Pnm[1, 1, 0];
            Pnm[1, 2, 0] = Pnm[1, 0, 0];
            for (int n = 2; n < nmax + 1; n++)
            {
                Pnm[0, n + 1, 0] = -Math.Sqrt((n * n + n) / 2) * Pnm[n, 1, 0];
                Pnm[1, n + 1, 0] = ((Math.Sqrt(2 * (n * n + n)) * Pnm[n, 0, 0] - Math.Sqrt((n * n + n - 2)) * Pnm[n, 2, 0]) / 2);

                for (int m = 2; m < n; m++)
                {
                    Pnm[m, n + 1, 0] = (0.5 * (Math.Sqrt((n + m) * (n - m + 1)) * Pnm[n, m - 1, 0] - Math.Sqrt((n + m + 1) * (n - m)) * Pnm[n, m + 1, 0]));
                }

                Pnm[n, n + 1, 0] = Math.Sqrt(2 * n) * Pnm[n, n - 1, 0] / 2;
            }

            return Pnm;

        }
        /// <summary>
        /// Функция обертка для legendre_poly, чтобы брать именно последний член вычисления.
        /// Возвращает полином и производную для заданного m и n.
        /// </summary>
        /// <param name="m">степень полинома</param>
        /// <param name="n">порядок полинома</param>
        /// <param name="x">Угол в радианах (предварительно вычисленное значение cos(alpha))</param>
        /// <returns>Pnm, dPnm</returns>
        public static (double, double) PnmNorm(int m, int n, double x)
        {
            double Pnm = 0;
            double dPnm = 0;

            var val = legendre_poly(n, x);
            if (m == 0)
            {
                Pnm = val[n, 0, 0];
                dPnm = val[0, n + 1, 0];
            }
            else
            {
                Pnm = val[n, m, 0];
                dPnm = val[m, n + 1, 0];
            }

            return (Pnm, dPnm);
        }
        /// <summary>
        /// Считываем файл коэффициентов IGRF13 с именем path и экстраполируем значения на дату date.
        /// !!!Внимание. Предварительно нужно стереть в последнем столбце тире и оставить последний год, то
        /// есть чтобы было не 2020-2025, а просто 2025.
        /// Возвращаем 2 матрицы для коэффициентов g и h.
        /// </summary>
        /// <param name="date">дата в которую нужно узнать значение поля по модели igrf</param>
        /// <param name="path">Название файла коэффициентов модели в папке с exe исполняемым файлом</param>
        /// <returns>2 матрицы для коэффициентов g и h</returns>
        public static (double[,], double[,]) MyReadCOEF(DateTime date, string path, int nmax=13)
        {
            int sz = nmax + 1;
            double[,] g = new double[sz, sz];
            double[,] h = new double[sz, sz];
            bool flag = true;

            //string path = "IGRF13.txt";

            string[] lines = File.ReadAllLines(path);

            string[] head_name;
            string[] head_type;
            double[] head_years = { 0 };


            // Создаем культуру, в которой точка будет являться разделителем десятичной части
            CultureInfo culture = (CultureInfo)CultureInfo.CurrentCulture.Clone();
            culture.NumberFormat.NumberDecimalSeparator = ".";
            //var t = Convert.ToDouble("1900.0", culture);

            foreach (string line in lines)
            {
                string[] s = line.Trim().Split().Where(x => !string.IsNullOrWhiteSpace(x)).ToArray();

                // Беру заголовки
                if (flag == true)
                {
                    if (line[0] == '#')
                    {
                        continue;
                    }

                    if (line.Substring(0, 2) == "cs")
                    {
                        head_name = s.ToArray();
                        continue;
                    }

                    if (line.Substring(0, 2) == "gh")
                    {
                        head_type = s.ToArray();

                        head_years = s.Skip(3).Select(j => Convert.ToDouble(j, culture)).ToArray();
                        flag = false;
                        continue;
                    }
                }

                string k = s[0];

                int n = int.Parse(s[1]);
                int m = int.Parse(s[2]);

                double[] coeffs = s.Skip(3).Select(j => Convert.ToDouble(j, culture)).ToArray();
                coeffs[coeffs.Length - 1] = coeffs[coeffs.Length - 2] + 5 * coeffs[coeffs.Length - 1];

                var f = Interpolate.CubicSpline(head_years, coeffs);
                double dateFloat = YearInFloat(date);
                double coef = f.Interpolate(dateFloat);

                if(nmax >= n) { // Если захочется поэксперементировать с количеством N
                    if (k == "g")
                    {
                        g[n, m] = coef;
                    }
                    else if (k == "h")
                    {
                        h[n, m] = coef;
                    }
                } 
            }

            return (g, h);
        }
        /// <summary>
        /// Функция для MyReadCOEF. Переводит значение даты в год с плавающей точкой. Это нужно для экстраполяции.
        /// Например, половина года 2023, будет 2022.5
        /// </summary>
        /// <param name="date">Дата</param>
        /// <returns>Значение года в десятичном виде</returns>
        public static double YearInFloat(DateTime date)
        {
            DateTime startOfYear = new DateTime(date.Year, 1, 1, 0, 0, 0); // Начало года
            //int daysPassed = (date - startOfYear).Days; // Дней прошло с начала года
            int daysAlls = (new DateTime(date.Year + 1, 1, 1, 0, 0, 0)
                         - startOfYear).Days; // Всего дней в году

            TimeSpan diffResult = date.Subtract(startOfYear); // Дней дробных прошло с начала года
            //double dayssss = diffResult.TotalDays; // Дней дробных прошло с начала года

            double result = diffResult.TotalDays / daysAlls;

            double Year = date.Year + result;

            return Year;
        }
        /// <summary>
        /// Функция вычисления значения компонент магнитного поля заданной точке и в указанную дату. 
        /// </summary>
        /// <param name="X">Точка по X в которой находится КА</param>
        /// <param name="Y">Точка по Y в которой находится КА</param>
        /// <param name="Z">Точка по Z в которой находится КА</param>
        /// <param name="date">Дата измерения поля</param>
        /// <param name="N">Выпуск модели IGRF по умолчанию 13</param>
        /// <param name="Rz">Радиус Земли в км по умолчанию 6371.2</param>
        /// <param name="path">имя файла содержащего коэффициенты модели по умолчанию "IGRF13.txt"</param>
        /// <returns>Bx, By, Bz и полное значение поля F</returns>
        public static (double, double, double, double, string) igrfCalculate(double X, double Y, double Z, DateTime date,
            int N = 13, double Rz = 6371.2, string path = "IGRF13.txt")
        {
            double FACT = 180 / Math.PI; // Деление на константу переводит в  радианы, умножение в градусы

            var coordSphere = xyz2dhif(X, Y, Z);
            var r = coordSphere[3]; // Радиус от центра Земли
            var theta = Math.PI / 2 - coordSphere[1]; // Широта
            var lambda = coordSphere[0]; // Долгота

            //r = 6470; theta = (90-25) / FACT; lambda = 50 / FACT; // Тестовое рабочее


            var gh = MyReadCOEF(date, path, N);
            var g = gh.Item1;
            var h = gh.Item2;


            var x = Math.Cos(theta);

            double resultUr = 0.0;
            double resultUt = 0.0;
            double resultUl = 0.0;

            // Получим Ur
            double result0 = 0;
            double result1 = 0;
            double result2 = 0;

            double radius = 6371.2 / r;
            double koef_start1 = 0.0;
            double koef_start2 = 0.0;
            for (int n = 1; n <= N; n++)
            {
                result0 = 0;
                result1 = 0;
                result2 = 0;

                koef_start1 = Math.Pow(Rz / r, n + 2);
                koef_start2 = -1 / (Math.Sin(theta));
                for (int m = 0; m <= n; m++)
                {
                    var Pnm = PnmNorm(m, n, x);
                    result0 += Pnm.Item1 * (g[n, m] * Math.Cos(m * lambda) + h[n, m] * Math.Sin(m * lambda));
                    result1 += Pnm.Item2 * (g[n, m] * Math.Cos(m * lambda) + h[n, m] * Math.Sin(m * lambda));
                    result2 += Pnm.Item1 * (m * (-g[n, m] * Math.Sin(m * lambda) + h[n, m] * Math.Cos(m * lambda)));
                }
                resultUr += Math.Pow(radius, n + 2) * (n + 1) * result0;
                resultUt += koef_start1 * result1;
                resultUl += koef_start1 * result2;
            }
            resultUt = -resultUt;
            resultUl = resultUl * koef_start2;

            //Br = resultUr
            //Bp = resultUl
            //Bt = resultUt

            var Bx = -resultUt;
            var By = resultUl;
            var Bz = -resultUr;

            //var ans = xyz2dhif(Bx, By, Bz);



            //return (ans[0]*FACT, ans[1]*FACT, ans[2]);
            return (Bx, By, Bz, Math.Sqrt(Bx * Bx + By * By + Bz * Bz), date.ToString("yyyy-M-d-H-m-s"));

        }

    }


}
