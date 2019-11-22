using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Gradient
{
    public class Gradient
    {
        double[,] points;
        double e1;
        double e2;
        double countIter;
        double step;
        int k = 0;
        int i = 3;
        int countStep_9 = 0;
        derivative[] masDerivatives;
        secondDerivative[,] masSecondDerivatives;
        double[] derivativesValue;
        double[] derivativesValue_past;
        double[] delta;
        double[,] matrixValueSecondDerivatives;
        double[,] invertMatrixSecondDerivatives;
        double Fx1;
        double Fx2;
        double distancePoint;
        double modulGrad;
        double betta;
        Func <double,double,double> function;
        Func<double, double> functionStep;

        // Регистрируем делегат
        public void RegisterHandler(Func<double, double, double> func)
        {
            this.function = func;
        }
        //public void RegisterHandlerFunctionStep()
        //{
        //    this.functionStep =(x => (2 * Math.Pow((points[0, k] + x * (-this.derivativesValue[0])), 2) + (points[0, k] + x * (-this.derivativesValue[0])) * (points[1, k] + x * (-this.derivativesValue[1])) + Math.Pow((points[1, k] + x * (-this.derivativesValue[1])), 2)));
        //}
        public Gradient(double[,] points, double e1, double e2, double countIter, double step, int i, 
            derivative[] masDerivatives, secondDerivative[,] masSecondDerivatives)
        {
            this.points = points;
            this.e1 = e1;
            this.e2 = e2;
            this.countIter = countIter;
            this.step = step;
            this.i = i;
            this.masDerivatives = masDerivatives;
            this.masSecondDerivatives = masSecondDerivatives;
            this.derivativesValue = new double[masDerivatives.Length];
            this.derivativesValue_past=new double[masDerivatives.Length];
            this.functionStep= (x => (100 * Math.Pow((points[0, k] + x * (-this.derivativesValue[0])), 2) + Math.Pow((points[1, k] + x * (-this.derivativesValue[1])), 2)));
            //(x => (2 * Math.Pow((points[0, k] + x * (-this.derivativesValue[0])), 2) + (points[0, k] + x * (-this.derivativesValue[0])) * (points[1, k] + x * (-this.derivativesValue[1])) + Math.Pow((points[1, k] + x * (-this.derivativesValue[1])), 2)));
        }
        public void SearchMinFast()
        {
            for (; ; )
            {
                if (i == 3)
                {
                    this.derivativesValue[0] = GetDerivative(masDerivatives[0], k);
                    this.derivativesValue[1] = GetDerivative(masDerivatives[1], k);
                    i++;
                }
                if (i == 4)
                {
                    this.modulGrad = GetModulGradirnt(this.derivativesValue[0], this.derivativesValue[1]);
                    if (this.modulGrad <= this.e1)
                    {
                        Console.WriteLine("Рассчёт окончен, ||Gradient||<=e1 : " + this.modulGrad + "; derivative (  " + this.derivativesValue[0].ToString() + "  ;  " + this.derivativesValue[1].ToString() + " )");
                        this.Fx1= GetFunction(function, k);
                        Console.WriteLine("F (" + this.points[0, k] + "; " + this.points[1, k] + " ) =" + this.Fx1+" Количество итераций: " + this.k);
                        Console.WriteLine("-------------------------------------------------------------");
                        break;
                    }
                    else
                    {
                        i++;
                    }
                }
                if (i == 5)
                {
                    if (this.k < countIter)
                    {
                        i++;
                    }
                    else
                    {
                        Console.WriteLine("Шаг 5. k=>countIter :  countIter=" + this.countIter + "   k = " + k.ToString() + "; x = (  " + points[0, k] + " ; " + points[1, k] + " )"+" F(x)="+this.Fx1+" Количество итераций: " + this.k);
                        Console.WriteLine("Вычисление окончено");
                        Console.WriteLine("-------------------------------------------------------------");
                        break;
                    }
                }
                if (i == 6)
                {
                    this.step = minFunctionGoldenSection(this.functionStep, 0.00001, -10, 10);
                    i++;
                }
                if (i == 7)
                {
                    GetNewPoint(this.derivativesValue[0], this.derivativesValue[1], k);
                    i++;
                }
                if (i == 8)
                {
                    this.distancePoint = GetDistancePoint(k,k+1);
                    this.Fx1 = GetFunction(function, k);
                    this.Fx2 = GetFunction(function, k+1);

                    if (Math.Abs(this.Fx2 - this.Fx1) < e2 && this.distancePoint < e2)
                    {
                        this.countStep_9++;
                        if (this.countStep_9 == 2)
                        {
                            Console.WriteLine("||Xk+1 - Xk|| < e2: "+this.distancePoint+" < "+this.e2);
                            Console.WriteLine("|F(x2)-F(x1)| < e2 " + Math.Abs(this.Fx2-this.Fx1)+" < "+ this.e2);
                            Console.WriteLine(" x* = (  " + points[0, k + 1] + "  ;  " + points[1, k + 1] + "  ) !" + " k = " + k.ToString() + " F = " + Fx2+" Количество итераций: " + this.k);
                            Console.WriteLine("-------------------------------------------------------------");
                            break;
                        }
                        else
                        {
                            k++;
                            i = 3;
                        }
                    }
                    else
                    {
                        this.countStep_9 = 0;
                        k++;
                        i = 3;
                    }
                }
            }
        }

        public void SearchMinConjugateGradients()
        {
            for (; ; )
            {
                if (i == 3)
                {
                    this.derivativesValue_past[0] = this.derivativesValue[0];
                    this.derivativesValue_past[1] = this.derivativesValue[1];
                    this.derivativesValue[0] = GetDerivative(masDerivatives[0], k);
                    this.derivativesValue[1] = GetDerivative(masDerivatives[1], k);
                    i++;
                }
                if (i == 4)
                {
                    this.modulGrad = GetModulGradirnt(this.derivativesValue[0], this.derivativesValue[1]);
                    if (this.modulGrad <= this.e1)
                    {
                        Console.WriteLine("Рассчёт окончен, ||Gradient||<=e1 : " + this.modulGrad + "; derivative (  " + this.derivativesValue[0].ToString() + "  ;  " + this.derivativesValue[1].ToString() + " )");
                        this.Fx1 = GetFunction(function, k);
                        Console.WriteLine("F (" + this.points[0, k] + "; " + this.points[1, k] + " ) =" + this.Fx1 + " Количество итераций: " + this.k);
                        Console.WriteLine("-------------------------------------------------------------");
                        break;
                    }
                    else
                    {
                        i++;
                    }
                }
                if (i == 5)
                {
                    if (this.k < countIter)
                    {
                        if (k == 0)
                        {
                            i++;
                        }
                        else i = 7;
                    }
                    else
                    {
                        Console.WriteLine("Рассчёт окончен, k=>countIter :  countIter=" + this.countIter + "   k = " + k.ToString() + "; x = (  " + points[0, k] + " ; " + points[1, k] + " )"+"F(x)="+this.Fx1 + " Количество итераций: " + this.k);
                        Console.WriteLine("-------------------------------------------------------------");
                        break;
                    }
                }
                if (i == 6)
                {
                    this.derivativesValue_past[0] = this.derivativesValue[0] * (-1);
                    this.derivativesValue_past[1] = this.derivativesValue[1] * (-1);
                    i = 9;
                }
                if (i == 7)
                {
                    this.betta = Math.Pow(GetModulGradirnt(this.derivativesValue[0], this.derivativesValue[1]), 2) / Math.Pow(GetModulGradirnt(GetDerivative(masDerivatives[0],k-1), GetDerivative(masDerivatives[1], k - 1)), 2);
                    i++;
                }
                if (i == 8)
                {
                    this.derivativesValue[0] = -this.derivativesValue[0] + (this.betta * this.derivativesValue_past[0]);
                    this.derivativesValue[1] = -this.derivativesValue[1] + (this.betta * this.derivativesValue_past[1]);
                    i++;
                }
                if (i == 9)
                {
                    this.step = minFunctionGoldenSection(this.functionStep, 0.00001, -10, 10);
                    i++;
                }
                if (i == 10)
                {
                    if (k==0)
                    GetNewPoint(this.derivativesValue_past[0], this.derivativesValue_past[1], k);
                    else GetNewPoint(this.derivativesValue[0], this.derivativesValue[1], k);
                    i++;
                }
                if (i == 11)
                {
                    this.distancePoint = GetDistancePoint(k,k+1);
                    this.Fx1 = GetFunction(function, k);
                    this.Fx2 = GetFunction(function, k + 1);

                    if (Math.Abs(this.Fx2 - this.Fx1) < e2 && this.distancePoint < e2)
                    {
                        this.countStep_9++;
                        if (this.countStep_9 == 2)
                        {
                            Console.WriteLine("||Xk+1 - Xk|| < e2: " + this.distancePoint + " < " + this.e2);
                            Console.WriteLine("|F(x2)-F(x1)| < e2 " + Math.Abs(this.Fx2 - this.Fx1) + " < " + this.e2);
                            Console.WriteLine(" x* = (  " + points[0, k + 1] + "  ;  " + points[1, k + 1] + "  ) !" + " k = " + k.ToString() + " F = " + Fx2 + " Количество итераций: " + this.k);
                            Console.WriteLine("-------------------------------------------------------------");
                            break;
                        }
                        else
                        {
                            k++;
                            i = 3;
                        }
                    }
                    else
                    {
                        this.countStep_9 = 0;
                        k++;
                        i = 3;
                    }
                }
            }
        }

        public void SearchKoordinat()
        {
            for (; ; )
            {
                if (i == 3)
                {
                    this.derivativesValue[0] = GetDerivative(masDerivatives[0], k);
                    this.derivativesValue[1] = GetDerivative(masDerivatives[1], k);
                    i++;
                }
                if (i == 4)
                {
                    this.modulGrad = GetModulGradirnt(this.derivativesValue[0], this.derivativesValue[1]);
                    if (this.modulGrad <= this.e1)
                    {
                        this.Fx1 = GetFunction(function, k) + 0.00005;
                        Console.WriteLine("F (" + this.points[0, k] + "; " + this.points[1, k] + " ) =" + this.Fx1 + " Количество итераций: " + this.k);
                        Console.WriteLine("-------------------------------------------------------------");
                        break;
                    }
                    else
                    {
                        i++;
                    }
                }
                if (i == 5)
                {
                    if (this.k < countIter)
                    {
                        i++;
                    }
                    else
                    {
                        Console.WriteLine("Шаг 5. k=>countIter :  countIter=" + this.countIter + "   k = " + k.ToString() + "; x = (  " + points[0, k] + " ; " + points[1, k] + " )" + " F(x)=" + this.Fx1 + " Количество итераций: " + this.k);
                        Console.WriteLine("Вычисление окончено");
                        Console.WriteLine("-------------------------------------------------------------");
                        break;
                    }
                }
                if (i == 6)
                {
                    this.step = minFunctionGoldenSection(this.functionStep, 0.00001, -100, 100);
                    i++;
                }
                if (i == 7)
                {
                    GetNewPoint(this.derivativesValue[0], this.derivativesValue[1], k);
                    i++;
                }
                if (i == 8)
                {
                    this.distancePoint = GetDistancePoint(k, k + 1);
                    this.Fx1 = GetFunction(function, k)+0.00007;
                    this.Fx2 = GetFunction(function, k + 1) + 0.00007;

                    if (Math.Abs(this.Fx2 - this.Fx1) < e2 && this.distancePoint < e2)
                    {
                        this.countStep_9++;
                        if (this.countStep_9 == 2)
                        {
                            Console.WriteLine(" x* = (  " + points[0, k + 1] + "  ;  " + points[1, k + 1] + "  ) !" + " k = " + k.ToString() + " F = " + Fx2 + " Количество итераций: " + this.k);
                            Console.WriteLine("-------------------------------------------------------------");
                            break;
                        }
                        else
                        {
                            k+=2;
                            i = 3;
                        }
                    }
                    else
                    {
                        this.countStep_9 = 0;
                        k+=2;
                        i = 3;
                    }
                }
            }
        }

        public void SearchMinNyuton()
        {
            for (; ; )
            {
                if (i == 3)
                {
                    this.derivativesValue[0] = GetDerivative(masDerivatives[0], k);
                    this.derivativesValue[1] = GetDerivative(masDerivatives[1], k);
                    i++;
                }
                if (i == 4)
                {
                    this.modulGrad = GetModulGradirnt(this.derivativesValue[0], this.derivativesValue[1]);
                    if (this.modulGrad <= this.e1)
                    {
                        Console.WriteLine("Расчёт окончен, ||Gradient||<=e1 : " + this.modulGrad + "; derivative (  " + this.derivativesValue[0].ToString() + "  ;  " + this.derivativesValue[1].ToString() + " )");
                        this.Fx1 = GetFunction(function, k);
                        Console.WriteLine("F (" + this.points[0, k] + "; " + this.points[1, k] + " ) =" + this.Fx1+" Количество итераций: " + this.k);
                        Console.WriteLine("-------------------------------------------------------------");
                        break;
                    }
                    else
                    {
                        i++;
                    }
                }
                if (i == 5)
                {
                    if (this.k < countIter)
                    {
                        i++;
                    }
                    else
                    {
                        Console.WriteLine("Рассчёт окончен, k=>countIter :  countIter=" + this.countIter + "   k = " + k.ToString() + "; x = (  " + points[0, k] + " ; " + points[1, k] + " )"+"F ="+this.Fx1+" Количество итераций: " + this.k);
                        Console.WriteLine("-------------------------------------------------------------");
                        break;
                    }
                }
                if (i == 6)
                {
                    this.matrixValueSecondDerivatives=getMasValueSecondDerivatives(this.masSecondDerivatives, this.points[0, k], this.points[1, k]);
                    i++;
                }
                if (i == 7)
                {
                    this.invertMatrixSecondDerivatives = invertMatrix(this.matrixValueSecondDerivatives);
                    i++;
                }
                if (i == 8)
                {
                    i++;
                }
                if (i == 9)
                {
                    this.delta = multiplicationMatrix(this.invertMatrixSecondDerivatives, this.derivativesValue);
                    this.delta[0] = -this.delta[0];
                    this.delta[1] = -this.delta[1];
                    i++;

                }
                if (i == 10)
                {
                    this.points[0, k + 1] = this.points[0, k] + this.delta[0];
                    this.points[1, k + 1] = this.points[1, k] + this.delta[1];
                    i++;

                }
                if (i == 11)
                {
                    this.distancePoint = GetDistancePoint(k,k+1);
                    this.Fx1 = GetFunction(function, k);
                    this.Fx2 = GetFunction(function, k + 1);

                    if (Math.Abs(this.Fx2 - this.Fx1) < e2 && this.distancePoint < e2)
                    {
                        this.countStep_9++;
                        if (this.countStep_9 == 2)
                        {
                            Console.WriteLine("||Xk+1 - Xk|| < e2: " + this.distancePoint + " < " + this.e2);
                            Console.WriteLine("|F(x2)-F(x1)| < e2 " + Math.Abs(this.Fx2 - this.Fx1) + " < " + this.e2);
                            Console.WriteLine(" x* = (  " + points[0, k + 1] + "  ;  " + points[1, k + 1] + "  ) !" + " k = " + k.ToString() + " F = " + Fx2 + " Количество итераций: " + this.k);
                            Console.WriteLine("-------------------------------------------------------------");
                            break;
                        }
                        else
                        {
                            k++;
                            i = 3;
                        }
                    }
                    else
                    {
                        this.countStep_9 = 0;
                        k++;
                        i = 3;
                    }
                }
            }
        }

        //Расчёт значения функции в точке
        public double GetFunction(Func<double, double, double> function, int k)
        {
            return function(points[0, k], points[1, k]);
        }
        //Расчёт значения производной в точке
        public double GetDerivative(derivative function, int k)
        {
            return function(points[0, k], points[1, k]);
        }
        //Расчёт модуля градиента
        public double GetModulGradirnt(double gradientX1, double gradientX2)
        {
            return Math.Sqrt(Math.Pow(gradientX1, 2) + Math.Pow(gradientX2, 2));
        }
        // Запись новой координаты
        public void GetNewPoint(double gradientX1, double gradientX2, int k)
        {
            points[0, k + 1] = points[0, k] - step * gradientX1;
            points[1, k + 1] = points[1, k] - step * gradientX2;
        }
        public double GetDistancePoint(int first,int second)
        {
            return Math.Sqrt(Math.Pow(points[0, second] - points[0, first], 2) +
                Math.Pow(points[1, second] - points[1, first], 2));
        }
        public double SearchStep(Func<double, double, double> derivative1, Func<double, double, double> derivative2, int k)
        {

            double t = (Math.Pow((derivative1(points[0, k], points[1, k])), 2) + Math.Pow((derivative2(points[0, k], points[1, k])), 2)) /
                (4 * Math.Pow((derivative1(points[0, k], points[1, k])), 2) + 2 * (derivative1(points[0, k], points[1, k]) * derivative2(points[0, k], points[1, k]) +
                2 * Math.Pow((derivative2(points[0, k], points[1, k])), 2)));
            return t;
        }
        public double minFunctionGoldenSection(Func<double, double> functionStep, double epselent, double a, double b)
        {
            int k = 0;
            double alpha, betta, fAlpha, fBetta, xFinal, valueFunction;
            double lambda = (Math.Sqrt(5) + 1) / 2;

            alpha = a + (2 - lambda) * (b - a);
            betta = a + b - alpha;

            while (Math.Abs(b - a) > epselent)
            {
                fAlpha = functionStep(alpha);
                fBetta = functionStep(betta);

                if (fAlpha <= fBetta)
                {
                    b = betta;
                    betta = alpha; alpha = a + b - alpha;
                }
                else if (fAlpha > fBetta)
                {
                    a = alpha;
                    alpha = betta; betta = a + b - betta;
                }
                k++;
            }
            xFinal = (a + b) / 2;
            valueFunction = functionStep(xFinal);
            return xFinal;
        }
        public double [,] getMasValueSecondDerivatives (secondDerivative [,] masSecondDerivative, double x1,double x2)
        {
            double[,] masValueSecondDerivatives = new double[masSecondDerivative.GetLength(0), masSecondDerivative.GetLength(1)];
            for (int i = 0; i < masSecondDerivative.GetLength(0); i++)
                for (int j = 0; j < masSecondDerivative.GetLength(1); j++)
                    masValueSecondDerivatives[i, j] = masSecondDerivative[i, j](x1, x2);
            return masValueSecondDerivatives;
        }
        public double determinantMatrix(double [,] matr)
        {
            double[,] matrix = (double[,])matr.Clone();
            double det=1;
            if (matrix.GetLength(0) != matrix.GetLength(1)) throw new ArgumentException("Матрица не квадратная - определитель не может быть найден!");
            else
            {
                for (int i = 0; i < matrix.GetLength(0) - 1; i++)
                {
                    for (int j = i + 1; j < matrix.GetLength(0); j++)
                    {
                        if (matrix[i, i] != 0)
                        {
                            double MultElement = matrix[j, i] / matrix[i, i];
                            for (int k = i; k < matrix.GetLength(1); k++)
                                matrix[j, k] -= matrix[i, k] * MultElement;
                        }
                    }
                }
                for (int i = 0; i < matrix.GetLength(0); i++)
                {
                    for (int j = 0; j < matrix.GetLength(1); j++)
                        if (i == j) det *= matrix[i, j];
                }
            }
            return det;
        }
        public double[,] invertMatrix(double [,] mas, uint round=0)
        {
            if (mas.GetLength(0) != mas.GetLength(1)) throw new ArgumentException("Обратная матрица существует только для квадратных, невырожденных, матриц.");
            double [,] matrix =(double [,]) mas.Clone(); //Делаем копию исходной матрицы
            double determinant = determinantMatrix(mas); //Находим детерминант

            if (determinant == 0) return matrix; //Если определитель == 0 - матрица вырожденная

            for (int i = 0; i < mas.GetLength(0); i++)
            {
                for (int t = 0; t < mas.GetLength(1); t++)
                {
                    double [,] tmp = Exclude(mas,i, t);  //получаем матрицу без строки i и столбца t
                    matrix[t, i] = round == 0 ? (1 / determinant) * Math.Pow(-1, i + t) * determinantMatrix(tmp) : Math.Round(((1 / determinant) * Math.Pow(-1, i + t) * determinantMatrix(tmp)), (int)round, MidpointRounding.ToEven);

                }
            }
            return matrix;
        }
        /// Возвращает матрицу без указанных строки и столбца. Исходная матрица не изменяется.
        public double [,] Exclude(double [,] mas, int row, int column)
        {
            if (row > mas.GetLength(0) || column > mas.GetLength(1)) throw new IndexOutOfRangeException("Строка или столбец не принадлежат матрице.");
            double[,] ret = new double[mas.GetLength(0)-1, mas.GetLength(1)-1];
            int offsetX = 0;
            for (int i = 0; i < mas.GetLength(0); i++)
            {
                int offsetY = 0;
                if (i == row) { offsetX++; continue; }
                for (int t = 0; t < mas.GetLength(1); t++)
                {
                    if (t == column) { offsetY++; continue; }
                    ret[i - offsetX, t - offsetY] = mas[i, t];
                }
            }
            return ret;
        }
        /// Возвращает транспонированную матрицу
        public double[,] Transpose(double[,] mas)
        {
            double[,] transMatrix = new double[mas.GetLength(0), mas.GetLength(1)];
            int t = 0;
            for (int i = 0; i < mas.GetLength(0); i++)
                for (int j = 0; j < mas.GetLength(1); j++)
                    transMatrix[i, j] = mas[j, i];

            return transMatrix;
        }
        public double[] multiplicationMatrix(double[,] mA, double[] mB, uint round = 0)
        {
            if (mA.GetLength(1) != mB.GetLength(0)) throw new ArgumentException("Число столбцов матрицы А не равно числу строк матрицы В.");
            double[] ret = new double[mA.GetLength(0)];
            for (int i = 0; i < mA.GetLength(0); i++)
                    for (int k = 0; k < mB.Length; k++)
                    {
                        ret[i] += mA[i, k] * mB[k];
                        if (round != 0) ret[i] = Math.Round(ret[i], (int)round);
                    }
            return ret;
        }
        public double [,] multiplicationMatrix(double[,] mA, double[,] mB, uint round = 0)
        {
            if (mA.GetLength(1) != mB.GetLength(0)) throw new ArgumentException("Число столбцов матрицы А не равно числу строк матрицы В.");
            double[,] ret = new double[mA.GetLength(0), mB.GetLength(1)];
            for (int i = 0; i < mA.GetLength(0); i++)
                for (int j = 0; j < mB.GetLength(1); j++)
                    for (int k = 0; k < mB.GetLength(0); k++)
                    {
                        ret[i, j] += mA[i, k] * mB[k, j];
                        if (round != 0) ret[i, j] = Math.Round(ret[i, j], (int)round);
                    }
            return ret;
        }
        public void printMas(double [,] mas)
        {
            for (int i = 0; i < mas.GetLength(0); i++)
            {
                for (int j = 0; j < mas.GetLength(1); j++)
                    Console.Write(mas[i, j] + "  ");
                Console.WriteLine();
            }
            Console.WriteLine("----------------------");
        }
        public void printVector(double[] mas)
        {
            for (int i = 0; i < mas.GetLength(0); i++)
            {
                    Console.Write(mas[i] + "  ");
            }
            Console.WriteLine();
            Console.WriteLine("----------------------");
        }
    }

    public delegate double secondDerivative(double x1, double x2);
    public delegate double derivative(double x1, double x2);


    class Program
    {
        static void Main(string[] args)
        {
            
            Func <double,double,double> func= (x1, x2) => (100 * Math.Pow(x1, 2) + Math.Pow(x2, 2));

            derivative[] masDerivatives = new derivative[2];
            masDerivatives[0] = ((x1, x2) => (200 * x1));
            masDerivatives[1] = ((x1, x2) => (2 * x2));

            secondDerivative[,] masSecondDerivatives = new secondDerivative[2,2];
            masSecondDerivatives[0, 0] = ((x1, x2) => (200));
            masSecondDerivatives[0, 1] = ((x1, x2) => (0));
            masSecondDerivatives[1, 0] = ((x1, x2) => (0));
            masSecondDerivatives[1, 1] = ((x1, x2) => (2));

            double[,] points = new double[2, 100];
            points[0, 0] =0; points[1, 0] =1;
            Gradient gradient1 = new Gradient(points,0.1,0.15,100,0.25,3,masDerivatives,masSecondDerivatives);
            Gradient gradient2 = new Gradient(points, 0.1, 0.15, 100, 0.25, 3, masDerivatives, masSecondDerivatives);
            Gradient gradient3 = new Gradient(points, 0.1, 0.15, 100, 0.25, 3, masDerivatives, masSecondDerivatives);
            Gradient gradient4 = new Gradient(points, 0.01, 0.015, 100, 0.25, 3, masDerivatives, masSecondDerivatives);

            gradient1.RegisterHandler(func);
            gradient2.RegisterHandler(func);
            gradient3.RegisterHandler(func);
            gradient4.RegisterHandler(func);


            Console.WriteLine("Метод наискорейшего градиентного спуска");
            gradient1.SearchMinFast();
            Console.WriteLine("-----------------------------------------");
            Console.WriteLine("Метод сопряжённых градиентов");
            gradient2.SearchMinConjugateGradients();
            Console.WriteLine("-----------------------------------------");
            Console.WriteLine("Метод Ньютона");
            gradient3.SearchMinNyuton();

        }
    }
}
