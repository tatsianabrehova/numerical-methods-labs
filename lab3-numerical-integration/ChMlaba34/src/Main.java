import java.util.function.DoubleUnaryOperator;

public class Main {

    // Подынтегральная функция: f(x) = cos((x + 1)^2)
    public static DoubleUnaryOperator f = x -> Math.cos(Math.pow(x + 1, 2));

    // Квадратурная формула Эрмита
    public static Result integrate(DoubleUnaryOperator func, double a, double b, double epsilon) {
        int n = 1;
        double I_last = Double.POSITIVE_INFINITY;
        double I_new = 0;

        while (Math.abs(I_new - I_last) >= epsilon) {
            I_last = I_new;
            I_new = 0;

            for (int k = 0; k <= n; k++) {
                double t_k = Math.cos((2 * k + 1) * Math.PI / (2 * n + 2)); // Узлы на [-1, 1]
                double A_k = Math.PI / (n + 1); // Веса

                // Преобразуем t_k из [-1, 1] в [a, b]
                double x_transformed = ((b - a) / 2.0) * t_k + (a + b) / 2.0;

                // Значение подынтегральной функции
                double fx = func.applyAsDouble(x_transformed);

                I_new += A_k * fx;
            }

            I_new *= (b - a) / 2.0; // Корректировка длины отрезка
            n++;
        }

        return new Result(I_new, n - 1);
    }

    // Класс для хранения результата
    public static class Result {
        public final double integral;
        public final int nodes;

        public Result(double integral, int nodes) {
            this.integral = integral;
            this.nodes = nodes;
        }
    }

    // Точка входа
    public static void main(String[] args) {
        double a = -1;
        double b = 1;
        double epsilon = 1e-6;

        Result result = integrate(f, a, b, epsilon);

        System.out.printf("Приближенное значение интеграла: %.10f%n", result.integral);
        System.out.printf("Число узлов: %d%n", result.nodes);
    }
}