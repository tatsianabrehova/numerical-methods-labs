import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

public class Main {
    // Сделали класс статическим
    class PolynomialFunction {
        private final double[] coefficients;

        public PolynomialFunction(double[] coefficients) {
            int i = coefficients.length - 1;
            while (i > 0 && coefficients[i] == 0) {
                i--;
            }
            this.coefficients = Arrays.copyOf(coefficients, i + 1);
        }

        public double value(double x) {
            double result = 0.0;
            for (int i = 0; i < coefficients.length; i++) {
                result += coefficients[i] * Math.pow(x, i);
            }
            return result;
        }

        public PolynomialFunction derivative() {
            if (coefficients.length <= 1) {
                return new PolynomialFunction(new double[]{0});
            }
            double[] derivativeCoefficients = new double[coefficients.length - 1];
            for (int i = 1; i < coefficients.length; i++) {
                derivativeCoefficients[i - 1] = coefficients[i] * i;
            }
            return new PolynomialFunction(derivativeCoefficients);
        }

        public double derivativeValue(double x) {
            double result = 0.0;
            for (int i = 1; i < coefficients.length; i++) {
                result += coefficients[i] * i * Math.pow(x, i - 1);
            }
            return result;
        }

        @Override
        public String toString() {
            if (coefficients.length == 0 || (coefficients.length == 1 && coefficients[0] == 0)) {
                return "0";
            }
            StringBuilder sb = new StringBuilder();
            for (int i = coefficients.length - 1; i >= 0; i--) {
                if (coefficients[i] == 0) {
                    continue;
                }
                if (i != coefficients.length - 1) {
                    sb.append(coefficients[i] >= 0 ? " + " : " - ");
                } else {
                    if (coefficients[i] < 0) {
                        sb.append("-");
                    }
                }
                double absCoeff = Math.abs(coefficients[i]);
                if (i == 0 || absCoeff != 1) {
                    sb.append(absCoeff);
                }
                if (i > 0) {
                    sb.append("x");
                    if (i > 1) {
                        sb.append("^").append(i);
                    }
                }
            }
            return sb.toString();
        }
    }


    public static void main(String[] args) {
        // Define the matrix A and vector f
        double[][] AData = {
                {15, -2, 3.5, 1, -0.1},
                {1, -8, -3.1, 1, 2.3},
                {-1, 3, 30, 2, 4.6},
                {0, 0.1, 2, 15, 2},
                {1, -2, 0.4, 3.2, -17}
        };
        double[] fData = {1, 0, 20, -3, 5};

        // Solve the system using Gaussian elimination
        double[] coef = gaussianElimination(AData, fData);
        System.out.println("Coefficients of the polynomial: " + Arrays.toString(coef));

        // Define the polynomial function P(x)
        PolynomialFunction P = new Main().new PolynomialFunction(coef);

        // Find the interval where the root is located using the bisection method
        double a = -3.0;
        double b = 0.0;
        double epsilon = 1e-1;
        List<double[]> tableList = new ArrayList<>();
        tableList.add(new double[]{a, b, b - a, P.value(a), P.value(b), (a + b) / 2});
        double c = 0.0;
        while (b - a > epsilon) {
            c = (a + b) / 2;
            if (P.value(c) * P.value(a) >= 0) {
                a = c;
            } else {
                b = c;
            }
            tableList.add(new double[]{a, b, b - a, P.value(a), P.value(b), (a + b) / 2});
        }

        // Print the table
        System.out.println("Bisection Method Table:");
        System.out.printf("%10s %10s %10s %10s %10s %10s%n", "a", "b", "b-a", "f(a)", "f(b)", "(a+b)/2");
        for (double[] row : tableList) {
            System.out.printf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f%n", row[0], row[1], row[2], row[3], row[4], row[5]);
        }

        // Now, refine the root using Newton's method
        double x0 = (a + b) / 2;
        double root = newtonMethod(P, x0, 1e-6);
        System.out.println("Approximate root: " + root);
    }

    // Gaussian elimination with partial pivoting
    public static double[] gaussianElimination(double[][] A, double[] f) {
        int n = A.length;
        double[][] M = new double[n][n];
        for (int i = 0; i < n; i++) {
            System.arraycopy(A[i], 0, M[i], 0, n);
        }
        double[] b = f.clone();

        for (int k = 0; k < n; k++) {
            // Partial pivoting
            int maxIndex = k;
            for (int i = k + 1; i < n; i++) {
                if (Math.abs(M[i][k]) > Math.abs(M[maxIndex][k])) {
                    maxIndex = i;
                }
            }

            // Swap rows
            if (maxIndex != k) {
                double[] tempRow = M[k];
                M[k] = M[maxIndex];
                M[maxIndex] = tempRow;

                double temp = b[k];
                b[k] = b[maxIndex];
                b[maxIndex] = temp;
            }

            // Eliminate
            for (int i = k + 1; i < n; i++) {
                double factor = M[i][k] / M[k][k];
                for (int j = k; j < n; j++) {
                    M[i][j] -= factor * M[k][j];
                }
                b[i] -= factor * b[k];
            }
        }

        // Back substitution
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < n; j++) {
                sum += M[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / M[i][i];
            System.out.printf("Root %d: x = %.6f%n", i, x[i]);
        }

        return x;
    }

    // Newton's method for finding the root
    public static double newtonMethod(PolynomialFunction P, double x0, double epsilon) {
        double xk = x0;
        System.out.println(P);
        double xk1 = xk - P.value(xk) / P.derivativeValue(xk);
        int k = 1;
        while (Math.abs(xk1 - xk) >= epsilon) {
            xk = xk1;
            xk1 = xk - P.value(xk) / P.derivativeValue(xk);
            k++;
        }
        System.out.println("Iterations for Newton: " + k);
        return xk1;
    }
}