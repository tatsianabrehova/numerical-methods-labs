import java.util.function.DoubleUnaryOperator;

public class Main {

    // === Метод правых прямоугольников ===
    public static double rightRectangles(DoubleUnaryOperator f, double a, double b, double h) {
        int n = (int) ((b - a) / h);
        double sum = 0;
        for (int k = 1; k <= n; k++) {
            sum += f.applyAsDouble(a + k * h);}
        return h * sum;}

    // === Метод трапеций ===
    public static double trapezoids(DoubleUnaryOperator f, double a, double b, double h) {
        int n = (int) ((b - a) / h);
        double sum = 0.5 * (f.applyAsDouble(a) + f.applyAsDouble(b));
        for (int k = 1; k < n; k++) {
            sum += f.applyAsDouble(a + k * h);}
        return h * sum;}

    // === Правило Рунге ===
    public static class RungeResult {
        public final double integral;
        public final int iterations;

        public RungeResult(double integral, int iterations) {
            this.integral = integral;
            this.iterations = iterations;
        }
    }

    public static RungeResult rungeRule(int m, java.util.function.BiFunction<DoubleUnaryOperator, Double, Double> method, DoubleUnaryOperator f, double a, double b, double epsilon) {
        double h = b - a;
        double I_prev = method.apply(f, h);
        h /= 2;
        double I_curr = method.apply(f, h);

        double R = (I_curr - I_prev) / (1 - Math.pow(0.5, m));
        int iterations = 1;

        while (Math.abs(R) >= epsilon) {
            h /= 2;
            I_prev = I_curr;
            I_curr = method.apply(f, h);
            R = (I_curr - I_prev) / (1 - Math.pow(0.5, m));
            iterations++;
        }

        return new RungeResult(I_curr, iterations);
    }

    // === Основной метод main ===
    public static void main(String[] args) {
        DoubleUnaryOperator f = x -> Math.sin(x * x) / (x + 1); // Подынтегральная функция
        double a = 1, b = 3;
        double epsilon = 1e-3;

        System.out.println("=== Метод правых прямоугольников ===");
        test((f1, h) -> rightRectangles(f1, a, b, h), f, epsilon, 1);

        System.out.println("\n=== Метод трапеций ===");
        test((f1, h) -> trapezoids(f1, a, b, h), f, epsilon, 2);

    }

    private static void test(java.util.function.BiFunction<DoubleUnaryOperator, Double, Double> method, DoubleUnaryOperator f, double epsilon, int m) {
        long startTime = System.nanoTime();
        RungeResult result = rungeRule(m, method, f, 1, 3, epsilon);
        long endTime = System.nanoTime();

        double durationSeconds = (endTime - startTime) / 1_000_000_000.0;

        System.out.printf("Интеграл: %.12f%n", result.integral);
        System.out.printf("Число разбиений: %d%n", result.iterations);
        System.out.printf("Время выполнения: %.3f секунд%n", durationSeconds);
    }
}