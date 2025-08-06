public class Main {

    // Метод для вычисления нормы вектора
    public static double vectorNorm(double[] vector) {
        double sum = 0.0;
        for (double v : vector) {
            sum += v * v;
        }
        return Math.sqrt(sum);
    }

    // Метод для умножения матрицы на вектор
    public static double[] matrixVectorMultiply(double[][] matrix, double[] vector) {
        int n = matrix.length;
        double[] result = new double[n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                result[i] += matrix[i][j] * vector[j];
            }
        }
        return result;
    }

    // Метод для вычисления определителя матрицы
    public static double determinant(double[][] matrix) {
        int n = matrix.length;
        if (n == 1) return matrix[0][0];
        if (n == 2) {
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        }
        double det = 0.0;
        for (int k = 0; k < n; k++) {
            det += Math.pow(-1, k) * matrix[0][k] * determinant(minor(matrix, 0, k));
        }
        return det;
    }

    // Метод для вычисления минора матрицы
    public static double[][] minor(double[][] matrix, int row, int col) {
        int n = matrix.length;
        double[][] result = new double[n - 1][n - 1];
        int r = 0;
        for (int i = 0; i < n; i++) {
            if (i == row) continue;
            int c = 0;
            for (int j = 0; j < n; j++) {
                if (j == col) continue;
                result[r][c] = matrix[i][j];
                c++;
            }
            r++;
        }
        return result;
    }

    // Метод для вычисления обратной матрицы
    public static double[][] inverse(double[][] matrix) {
        int n = matrix.length;
        double det = determinant(matrix);
        if (det == 0) throw new IllegalArgumentException("Matrix is singular, cannot be inverted.");
        double[][] adjugate = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                adjugate[j][i] = Math.pow(-1, i + j) * determinant(minor(matrix, i, j));
            }
        }
        double[][] inverse = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                inverse[i][j] = adjugate[i][j] / det;
            }
        }
        return inverse;
    }

    // Метод для выполнения градиентного спуска
    public static double[] gradientDescent(double[][] A, double[] b, double epsilon) {
        int n = b.length;
        double[] x_k = b.clone(); // Начальное приближение x_0 = g
        double[] r_k = matrixVectorMultiply(A, x_k); // Невязка
        for (int i = 0; i < n; i++) {
            r_k[i] -= b[i];
        }

        double[] x_k1 = new double[n];
        int k = 0;

        while (true) {
            k++;
            double numerator = 0.0;
            double denominator = 0.0;

            for (int i = 0; i < n; i++) {
                numerator += r_k[i] * r_k[i];
            }

            double[] Ar_k = matrixVectorMultiply(A, r_k);
            for (int i = 0; i < n; i++) {
                denominator += Ar_k[i] * r_k[i];
            }

            double tau = numerator / denominator; // Коэффициент спуска
            for (int i = 0; i < n; i++) {
                x_k1[i] = x_k[i] - tau * r_k[i];
            }

            if (vectorNorm(x_k1) - vectorNorm(x_k) < epsilon) {
                break; // Условие выхода
            }

            x_k = x_k1.clone();
            r_k = matrixVectorMultiply(A, x_k);
            for (int i = 0; i < n; i++) {
                r_k[i] -= b[i];
            }
        }

        return x_k1;
    }

    public static void main(String[] args) {
        // Матрица A
        double[][] A = {
                {4, -1, 0, 0},
                {-1, 4, -1, 0},
                {0, -1, 4, -1},
                {0, 0, -1, 3}
        };

        // Вектор b
        double[] b = {15, 10, 10, 10};

        // Точность
        double epsilon = 1e-6;

        // Градиентный спуск
        double[] x = gradientDescent(A, b, epsilon);
        System.out.println("Solution x:");
        for (double v : x) {
            System.out.printf("%.6f ", v);
        }

        // Обратная матрица
        double[][] A_inv = inverse(A);
        System.out.println("\nInverse of A:");
        for (double[] row : A_inv) {
            for (double v : row) {
                System.out.printf("%.6f ", v);
            }
            System.out.println();
        }

        // Определитель
        double det = determinant(A);
        System.out.println("Determinant of A: " + det);

        // Норма
        double norm = vectorNorm(matrixVectorMultiply(A, x));
        System.out.println("Norm of A: " + norm);

        // Невязка
        double[] residual = matrixVectorMultiply(A, x);
        for (int i = 0; i < b.length; i++) {
            residual[i] -= b[i];
        }
        System.out.println("Residual (Ax - b):");
        for (double v : residual) {
            System.out.printf("%.6f ", v);
        }
    }
}
