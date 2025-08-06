import java.util.Arrays;

public class Main {

    // 1. Исходная функция
    public static double f(double x) {
        return x * Math.sin(x);
    }

    // 2. Функция для вычисления конечных разностей
    public static double[][] finiteDifferences(double[] x, double[] y) {
        int n = x.length;
        double[][] F = new double[n][n];

        // Заполнение первого столбца значениями функции
        for (int i = 0; i < n; i++) {
            F[i][0] = y[i];
        }

        // Вычисление конечных разностей
        for (int j = 1; j < n; j++) {
            for (int i = 0; i < n - j; i++) {
                F[i][j] = F[i + 1][j - 1] - F[i][j - 1];
            }
        }

        return F;
    }

    // 3. Интерполяция Ньютона
    public static double newtonInterpolation(double[] x, double[] y, double xNew) {
        int n = x.length;
        double[][] F = finiteDifferences(x, y);

        double h = x[1] - x[0]; // Шаг
        double t = (xNew - x[0]) / h; // Параметр t

        double result = y[0];
        double product = 1;

        // Вычисление многочлена Ньютона
        for (int i = 1; i < n; i++) {
            product *= (t - (i - 1)) / i;
            result += F[0][i] * product;
        }

        return result;
    }

    public static void main(String[] args) {
        // 4. Параметры задачи
        double a = 0, b = 2;   // Интервал [a, b]
        double h = 0.1;        // Шаг
        double newH = 0.05;    // Новый шаг
        int m = 4;             // Количество точек для интерполяции

        // 5. Исходная таблица значений
        int originalSize = (int) ((b - a) / h) + 1;
        double[] xOriginal = new double[originalSize];
        double[] yOriginal = new double[originalSize];
        for (int i = 0; i < originalSize; i++) {
            xOriginal[i] = a + i * h;
            yOriginal[i] = f(xOriginal[i]);
        }

        // 6. Новая таблица значений
        int newSize = (int) ((b - a) / newH) + 1;
        double[] xNew = new double[newSize];
        double[] yNew = new double[newSize];
        for (int i = 0; i < newSize; i++) {
            xNew[i] = a + i * newH;
        }

        // 7. Интерполяция для каждой новой точки
        for (int i = 0; i < newSize; i++) {
            // Находим ближайший узел
            int idx = Arrays.binarySearch(xOriginal, xNew[i]);
            if (idx < 0) {
                idx = -(idx + 1) - 1;
            }

            // Определяем диапазон для интерполяции
            int start = Math.max(0, idx - m / 2 + 1);
            int end = Math.min(start + m, xOriginal.length);
            start = Math.max(0, end - m);

            // Выбираем точки для интерполяции
            double[] xInterp = Arrays.copyOfRange(xOriginal, start, end);
            double[] yInterp = Arrays.copyOfRange(yOriginal, start, end);

            // Вычисляем интерполированное значение
            yNew[i] = newtonInterpolation(xInterp, yInterp, xNew[i]);
        }

        // 8. Округление результатов
        for (int i = 0; i < newSize; i++) {
            yNew[i] = Math.round(yNew[i] * 100000.0) / 100000.0;
        }

        // 9. Печать таблиц
        System.out.println("Исходная таблица значений (шаг h=0.1):");
        for (int i = 0; i < originalSize; i++) {
            System.out.printf("x=%.1f, f(x)=%.4f%n", xOriginal[i], yOriginal[i]);
        }

        System.out.println("\nНовая таблица значений (шаг h=0.05):");
        for (int i = 0; i < newSize; i++) {
            System.out.printf("x=%.2f, f(x)=%.4f%n", xNew[i], yNew[i]);
        }
    }
}