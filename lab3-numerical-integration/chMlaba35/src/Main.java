import java.util.ArrayList;
import java.util.List;

public class Main {

    // Точность вычисления интеграла
    private static final double EPSILON_INTEGRAL = 1e-9;

    // Точность для метода Ньютона
    private static final double EPSILON_NEWTON = 1e-4;

    // Подынтегральная функция g(t) = e^(t^2) * (1 + t^3)
    public static double g(double t) {
        return Math.exp(Math.pow(t, 2)) * (1 + Math.pow(t, 3));
    }

    // Формула Симпсона
    public static double simpson(DoubleFunction f, double a, double b, double h) {
        int n = (int) ((b - a) / h);
        if (n % 2 != 0) n--; // n должно быть четным

        double result = f.apply(a);
        double x = a;
        for (int i = 1; i < n; i++) {
            x += h;
            if (i % 2 == 1) {
                result += 4 * f.apply(x);
            } else {
                result += 2 * f.apply(x);
            }
        }
        result += f.apply(b);
        return result * h / 3;
    }

    // Правило Рунге для автоматического выбора шага
    public static double rungeRuleSimpson(DoubleFunction f, double a, double b, double epsilon) {
        double h = b - a;
        double I_prev = simpson(f, a, b, h);
        h /= 2;
        double I_curr = simpson(f, a, b, h);

        double R = Math.abs(I_curr - I_prev);
        while (R >= epsilon) {
            h /= 2;
            I_prev = I_curr;
            I_curr = simpson(f, a, b, h);
            R = Math.abs(I_curr - I_prev);
        }
        return I_curr;
    }

    // Функция f(x) = ∫₁ˣ g(t) dt − 1.7
    public static double f(double x) {
        return rungeRuleSimpson(Main::g, 1, x, EPSILON_INTEGRAL) - 1.7;
    }

    // Метод бисекции (для уточнения отрезка)
    public static void bisect(DoubleFunction f, double a, double b, double epsilon) {
        System.out.println("=== Метод бисекции ===");
        List<String> table = new ArrayList<>();
        table.add(String.format("%-6s %-6s %-6s %-12s %-12s %-6s", "a", "b", "b-a", "f(a)", "f(b)", "(a+b)/2"));
        table.add("-------------------------------------------------------------");

        double fa = f.apply(a);
        double fb = f.apply(b);
        table.add(String.format("%.2f   %.2f   %.2f   %.6f     %.6f     %.2f", a, b, b - a, fa, fb, (a + b) / 2));

        while (b - a > epsilon) {
            double c = (a + b) / 2;
            double fc = f.apply(c);
            if (fc * fa >= 0) {
                a = c;
                fa = fc;
            } else {
                b = c;
                fb = fc;
            }
            table.add(String.format("%.2f   %.2f   %.2f   %.6f     %.6f     %.2f", a, b, b - a, fa, fb, (a + b) / 2));
        }

        for (String row : table) {
            System.out.println(row);
        }
        System.out.println();
    }

    // Вторая производная g'(x) = e^(x²)*x*(2+3x+2x³)
    public static double derivativeG(double x) {
        return Math.exp(x * x) * x * (2 + 3 * x + 2 * Math.pow(x, 3));
    }

    // Метод Ньютона
    public static void newtonMethod(double x0) {
        System.out.println("=== Метод Ньютона ===");
        System.out.println("x₀ = " + x0);

        double h0 = -f(x0) / g(x0);
        System.out.println("h₀ = " + h0);

        double s0Start = x0;
        double s0End = x0 + 2 * h0;
        System.out.println("s₀ = [" + s0Start + "; " + s0End + "]");

        // Проверка условия f(x)*f'(x) ≠ 0 на концах отрезка
        System.out.println("f(s₀[0]) * g(s₀[0]) = " + f(s0Start) * g(s0Start));
        System.out.println("f(s₀[-1]) * g(s₀[-1]) = " + f(s0End) * g(s0End));

        // Вычисление M = max |g’(x)| на [s₀]
        double M = 0;
        for (double x = s0Start; x <= s0End; x += 0.0001) {
            double dg = derivativeG(x);
            if (Math.abs(dg) > M) {
                M = Math.abs(dg);
            }
        }

        boolean condition = 2 * Math.abs(h0) * M <= Math.abs(g(x0));
        System.out.println("Условие 2|h₀|M ≤ |g(x₀)| выполнено? " + condition);

        // Итерации Ньютона
        System.out.println("\nИтерационный процесс:");
        System.out.printf("%-10s %-20s%n", "x_k", "|x_k+1 - x_k|");
        double xk = x0;
        double xk1 = phi(xk);
        System.out.printf("%.8f   %.2e%n", xk1, Math.abs(xk1 - xk));

        while (Math.abs(xk1 - xk) >= EPSILON_NEWTON) {
            xk = xk1;
            xk1 = phi(xk);
            System.out.printf("%.8f   %.2e%n", xk1, Math.abs(xk1 - xk));
        }

        System.out.println("\nПолученный корень: " + xk1);
        System.out.println("Значение f(x_k1): " + f(xk1));
    }

    // Функция φ(x) = x - f(x)/f'(x)
    public static double phi(double x) {
        return x - f(x) / g(x);
    }

    // Интерфейс для функции
    @FunctionalInterface
    interface DoubleFunction {
        double apply(double x);
    }

    // Точка входа
    public static void main(String[] args) {
        // Отделение корня
        System.out.println("=== Значения функции f(x) ===");
        System.out.println("f(1.2) = " + f(1.2));
        System.out.println("f(2) = " + f(2));

        // Метод бисекции
        bisect(Main::f, 1.2, 2, 0.1);

        // Метод Ньютона
        newtonMethod(1.2);
    }
}