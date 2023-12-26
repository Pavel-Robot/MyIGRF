
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using static IGRF_Module.FunctionIGRF;
using static System.Windows.Forms.VisualStyles.VisualStyleElement;


namespace IGRF_Module
{
    public partial class Form1 : Form
    {
        public string pathInputFile;
        public string pathOutputFile;

        public int valN;
        public double valRz;
        public string valNameCoef;

        
        

        public Form1()
        {
            InitializeComponent();
            textBox1.Text = "Выберите входной файл name.txt со строками X Y Z DATE(год-мес-день-час-мес-сек). Пример:1.1 2.2 3.3 2023-1-2-3-4-5";
            
            textBox2.Text = "13";
            textBox3.Text = "6371.2";
            textBox4.Text = "IGRF13.txt";
        }
        /// <summary>
        /// Функция считывает по пути заданному сразу все входные данные
        /// Потом подставляет их построчно в функцию igrfCalculate
        /// Считает для всех значение
        /// Сохраняет по новому пути файл
        /// </summary>
        public static void igrfVector(string pathInputFile, string pathOutputFile, int valN, double valRz, string valNameCoef)
        {
            //var igrf = igrfCalculate(2300, 4900, 3500, new DateTime(2023, 1, 1, 0, 0, 0));
            //valN = int.Parse(TextBox2.text);
            //double valRz = Convert.ToDouble(TextBox3, culture);
            //string valNameCoef = TextBox4.text;

            string[] lines = File.ReadAllLines(pathInputFile);

            // Создаем культуру, в которой точка будет являться разделителем десятичной части
            CultureInfo culture = (CultureInfo)CultureInfo.CurrentCulture.Clone();
            culture.NumberFormat.NumberDecimalSeparator = ".";

            string[] s;
            string[] d;
            List<string> Output = new List<string>(); // Объект уведомлений

            foreach (string line in lines)
            {
                s = line.Split(' '); // Разбиваем строку по пробелам и получаем массив строк

                d = s[3].Split('-');
                
                var igrf = igrfCalculate(Convert.ToDouble(s[0], culture),
                    Convert.ToDouble(s[1], culture),
                    Convert.ToDouble(s[2], culture),
                    new DateTime(int.Parse(d[0]), int.Parse(d[1]), int.Parse(d[2]), int.Parse(d[3]), int.Parse(d[4]), int.Parse(d[5])),
                    valN,
                    valRz,
                    valNameCoef);

                Output.Add(Convert.ToString(igrf.Item1, culture) + " " +
                           Convert.ToString(igrf.Item2, culture) + " " +
                           Convert.ToString(igrf.Item3, culture) + " " +
                           Convert.ToString(igrf.Item4, culture) + " " +
                           igrf.Item5
                           );

                
            }

            //Тут сохранение в файл
            using (StreamWriter writer = new StreamWriter(pathOutputFile))
            {
                foreach (string item in Output)
                {
                    writer.WriteLine(item);
                }
            }
        }

        private void button1_Click(object sender, EventArgs e)
        {
            OpenFileDialog openFile = new OpenFileDialog();
            openFile.Filter = "Text files(*.txt)|*.txt|All files(*.*)|*.*";

            if (openFile.ShowDialog() == DialogResult.OK)
            {
                textBox1.Text = openFile.FileName;
                pathInputFile = openFile.FileName;
            }
        }

        private void textBox1_TextChanged(object sender, EventArgs e)
        {

        }

        private void MagneticField_Load(object sender, EventArgs e)
        {

        }

        private void label3_Click(object sender, EventArgs e)
        {

        }

        private void textBox4_TextChanged(object sender, EventArgs e)
        {

        }

        private void label2_Click(object sender, EventArgs e)
        {

        }

        private void label1_Click(object sender, EventArgs e)
        {

        }

        private void textBox5_TextChanged(object sender, EventArgs e)
        {

        }

        private void button3_Click(object sender, EventArgs e)
        {
            // Создаем культуру, в которой точка будет являться разделителем десятичной части
            CultureInfo culture = (CultureInfo)CultureInfo.CurrentCulture.Clone();
            culture.NumberFormat.NumberDecimalSeparator = ".";

            
            valN = int.Parse(textBox2.Text);
            valRz = Convert.ToDouble(textBox3.Text, culture);
            valNameCoef = textBox4.Text;

            igrfVector(pathInputFile, pathOutputFile, valN, valRz, valNameCoef);

            //если выполнено
            MessageBox.Show("Программа выполнена");

        }

        private void button2_Click(object sender, EventArgs e)
        {
            SaveFileDialog FBD = new SaveFileDialog();
            FBD.Filter = "Text files(*.txt)|*.txt|All files(*.*)|*.*";

            if (FBD.ShowDialog() == DialogResult.OK)
            {
                textBox5.Text = FBD.FileName;
                pathOutputFile = FBD.FileName;

            }
        }
    }
}
