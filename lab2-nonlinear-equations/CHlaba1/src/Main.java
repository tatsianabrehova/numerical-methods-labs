
public class Main {

    public static void main(String[] args) {
        double epsilon = 1e-6;
        double x0 = 3.25; // Начальное приближение

        // Метод Ньютона с постоянной производной
        double rootNewton = newtonMethod(x0, epsilon);
        System.out.println("Метод Ньютона: " + rootNewton);

        // Метод секущих
        double rootSecant = secantMethod(x0, x0 + 0.1, epsilon);
        System.out.println("Метод секущих: " + rootSecant);

        // Метод Стеффенсона
        double rootSteffensen = steffensenMethod(x0, epsilon);
        System.out.println("Метод Стеффенсона: " + rootSteffensen);

        // Метод Чебышева третьего порядка
        double rootChebyshev = chebyshevMethod(x0, epsilon);
        System.out.println("Метод Чебышева: " + rootChebyshev);
    }

    // Функция f(x) = 2sin(3x) - (x^2 - 4x + 3)
    public static double f(double x) {
        return 2 * Math.sin(3 * x) - (x * x - 4 * x + 3);
    }

    // Производная f'(x) = 6cos(3x) - 2x + 4
    public static double f_derivative(double x) {
        return 6 * Math.cos(3 * x) - 2 * x + 4;
    }

    // Вторая производная f''(x) = -18sin(3x) - 2
    public static double f_secondDerivative(double x) {
        return -18 * Math.sin(3 * x) - 2;
    }

    // Метод Ньютона с постоянной производной
    public static double newtonMethod(double x0, double epsilon) {
        double xk = x0;
        double xk1 = xk - f(xk) / f_derivative(xk);
        int k = 1;
        System.out.printf("Iteration %d: x = %.6f%n", k, xk1);
        while (Math.abs(xk1 - xk) >= epsilon) {
            xk = xk1;
            xk1 = xk - f(xk) / f_derivative(xk);
            k++;
            System.out.printf("Iteration %d: x = %.6f%n", k, xk1);
        }
        System.out.println("Iterations for Newton: " + k);
        return xk1;
    }

    // Метод секущих
    public static double secantMethod(double x0, double x1, double epsilon) {
        double xk = x0;
        double xk1 = x1;
        double xk2 = xk1 - f(xk1) * (xk1 - xk) / (f(xk1) - f(xk));
        int k = 2;
        System.out.printf("Iteration %d: x = %.6f%n", k, xk1);
        while (Math.abs(xk2 - xk1) >= epsilon) {
            xk = xk1;
            xk1 = xk2;
            xk2 = xk1 - f(xk1) * (xk1 - xk) / (f(xk1) - f(xk));
            k++;
            System.out.printf("Iteration %d: x = %.6f%n", k, xk1);
        }
        System.out.println("Iterations for Secant: " + k);
        return xk2;
    }

    // Метод Стеффенсона
    public static double steffensenMethod(double x0, double epsilon) {
        double xk = x0;
        double xk1 = xk - Math.pow(f(xk), 2) / (f(xk + f(xk)) - f(xk));
        int k = 1;
        System.out.printf("Iteration %d: x = %.6f%n", k, xk1);
        while (Math.abs(xk1 - xk) >= epsilon) {
            xk = xk1;
            xk1 = xk - Math.pow(f(xk), 2) / (f(xk + f(xk)) - f(xk));
            k++;
            System.out.printf("Iteration %d: x = %.6f%n", k, xk1);
        }
        System.out.println("Iterations for Steffensen: " + k);
        return xk1;
    }

    // Метод Чебышева третьего порядка
    public static double chebyshevMethod(double x0, double epsilon) {
        double xk = x0;
        double xk1 = xk - (f(xk) / f_derivative(xk)) - (Math.pow(f(xk), 2) * f_secondDerivative(xk)) / (2 * Math.pow(f_derivative(xk), 3));
        int k = 1;
        System.out.printf("Iteration %d: x = %.6f%n", k, xk1);
        while (Math.abs(xk1 - xk) >= epsilon) {
            xk = xk1;
            xk1 = xk - (f(xk) / f_derivative(xk)) - (Math.pow(f(xk), 2) * f_secondDerivative(xk)) / (2 * Math.pow(f_derivative(xk), 3));
            k++;
            System.out.printf("Iteration %d: x = %.6f%n", k, xk1);
        }
        System.out.println("Iterations for Chebyshev: " + k);
        return xk1;
    }
}