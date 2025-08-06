import java.util.Arrays;

public class Main {

    public static double f1(double x) {
        return Math.cos(x);
    }

    public static double f2(double x) {
        return Math.abs(x) - 1;
    }

    // Создание равномерных узлов
    public static double[] uniformNodes(double a, double b, int n) {
        double[] nodes = new double[n + 1];
        double step = (b - a) / n;

        for (int i = 0; i <= n; i++) {
            nodes[i] = a + i * step;
        }
        return nodes;
    }

    // Создание узлов Чебышева
    public static double[] chebyshevNodes(double a, double b, int n) {
        double[] nodes = new double[n + 1];

        for (int i = 0; i <= n; i++) {
            double cosValue = Math.cos((2 * i + 1) * Math.PI / (2 * (n + 1)));
            nodes[i] = (a + b) / 2 + (b - a) / 2 * cosValue;
        }
        return nodes;
    }

    // Интерполяция Ньютона
    public static double newtonInterpolation(double[] xNodes, double[] yNodes, double x) {
        int n = xNodes.length;
        double[][] dividedDifferences = new double[n][n];

        // Заполнение таблицы разделённых разностей
        for (int i = 0; i < n; i++) {
            dividedDifferences[i][0] = yNodes[i];
        }
        for (int j = 1; j < n; j++) {
            for (int i = 0; i < n - j; i++) {
                dividedDifferences[i][j] = (dividedDifferences[i + 1][j - 1] - dividedDifferences[i][j - 1])
                        / (xNodes[i + j] - xNodes[i]);
            }
        }

        // Вычисление многочлена Ньютона
        double result = dividedDifferences[0][0];
        for (int i = 1; i < n; i++) {
            double term = dividedDifferences[0][i];
            for (int j = 0; j < i; j++) {
                term *= (x - xNodes[j]);
            }
            result += term;
        }
        return result;
    }

    public static void main(String[] args) {
        double a = -3, b = 3;
        int[] degrees = {3, 5, 7, 10, 15, 20}; // Степени интерполяции

        System.out.println("\nМетод Ньютона:");
        interpolateFunction(a, b, degrees, "newton", Main::f1, "cos(x)");
        interpolateFunction(a, b, degrees, "newton", Main::f2, "|x| - 1");
    }

    // Выполнение интерполяции для заданной функции
    public static void interpolateFunction(double a, double b, int[] degrees, String method,
                                           Function<Double, Double> func, String funcName) {
        System.out.println("Интерполяция функции " + funcName + " методом " + method + ":");

        for (int n : degrees) {
            double[] xNodes = uniformNodes(a, b, n);
            double[] yNodes = Arrays.stream(xNodes).map(func::apply).toArray();

            double[] xChebyshev = chebyshevNodes(a, b, n);
            double[] yChebyshev = Arrays.stream(xChebyshev).map(func::apply).toArray();

            System.out.println("\nСтепень интерполяции n = " + n);

            System.out.println("Равномерные узлы:");
            for (int i = 0; i < xNodes.length; i++) {
                double interpolatedValue = newtonInterpolation(xNodes, yNodes, xNodes[i]);
                System.out.printf("x = %.4f, f(x) = %.4f, Интерполяция = %.4f%n",
                        xNodes[i], yNodes[i], interpolatedValue);
            }

            System.out.println("Узлы Чебышева:");
            for (int i = 0; i < xChebyshev.length; i++) {
                double interpolatedValue = newtonInterpolation(xNodes, yNodes, xNodes[i]);
                System.out.printf("x = %.4f, f(x) = %.4f, Интерполяция = %.4f%n",
                        xChebyshev[i], yChebyshev[i], interpolatedValue);
            }
        }
    }

    // Функциональный интерфейс для передачи функций
    @FunctionalInterface
    interface Function<T, R> {
        R apply(T t);
    }
}