import java.util.Arrays;  // Import the Arrays class

public class Main {

    public static void main(String[] args) {
        // Initialize a 4x4 matrix A
        double[][] A = {
                {4, 1, 2, 3},
                {1, 5, 6, 7},
                {2, 6, 8, 9},
                {3, 7, 9, 10}
        };

        int N = 4;
        double[] SX = new double[N];
        double ALPHA, BETA, G;

        for (int I = 0; I < N - 1; I++) {
            // Initialize SX to 0
            for (int K = I; K < N; K++) {
                SX[K] = A[N - 1][I] * A[N - 1][K];
            }
            for (int J = N - 2; J >= I; J--) {
                for (int K = I; K < N; K++) {
                    SX[K] = SX[K] + A[J][I] * A[J][K];
                }
            }

            ALPHA = Math.sqrt(SX[I]);
            if (A[I][I] != 0) {
                BETA = 1.0 / ALPHA;
                for (int J = I + 1; J < N; J++) {
                    A[J][I] = A[J][I] * BETA;
                }
                SX[I] = A[I][I] * BETA + signum(A[I][I]);
                A[I][I] = ALPHA;
                G = 1.0 / Math.abs(SX[I]); // 1/gamma

                for (int K = I + 1; K < N; K++) {
                    SX[K] = SX[K] * BETA * G + signum(A[I][K], SX[I]);
                    A[I][K] = A[I][K] + SX[K] * SX[I];
                    for (int J = I + 1; J < N; J++) {
                        A[J][K] = A[J][K] - A[J][I] * SX[K];
                    }
                }
            } else {
                if (ALPHA != 0) {
                    BETA = 1.0 / ALPHA;
                    for (int J = I + 1; J < N; J++) {
                        A[J][I] = A[J][I] * BETA;
                    }
                    SX[I] = -1;
                    A[I][I] = ALPHA;
                    G = 1.0; // 1/gamma
                    for (int K = I + 1; K < N; K++) {
                        SX[K] = SX[K] * BETA * G + signum(A[I][K], SX[I]);
                        A[I][K] = A[I][K] + SX[K] * SX[I];
                        for (int J = I + 1; J < N; J++) {
                            A[J][K] = A[J][K] - A[J][I] * SX[K];
                        }
                    }
                } else {
                    SX[I] = 1;
                    G = 2.0;
                    for (int K = I + 1; K < N; K++) {
                        SX[K] = 2.0;
                        A[I][K] = A[I][K] - SX[K];
                    }
                }
            }
        }

        // Print the transformed matrix A
        System.out.println("Transformed Matrix A:");
        printMatrix(A);

        // Step 1: Calculate the determinant
        double det = determinant(A);
        System.out.println("Determinant of A: " + det);

        // Step 2: Calculate the norm (Frobenius norm)
        double norm = frobeniusNorm(A);
        System.out.println("Frobenius Norm of A: " + norm);

        // Step 3: Find the inverse of the matrix
        double[][] inverseA = inverse(A);
        System.out.println("Inverse of A:");
        printMatrix(inverseA);

        // Step 4: Solving the system A * x = b
        double[] b = {1, 2, 3, 4}; // Right-hand side column vector
        double[] x = solveSystem(A, b);
        System.out.println("Solution x:");
        System.out.println(Arrays.toString(x));
    }

    // Custom signum function to mimic Fortran's SIGN(X, Y)
    public static double signum(double x) {
        return x < 0 ? -1 : 1;
    }

    // Custom signum function to mimic Fortran's SIGN(X, Y)
    public static double signum(double x, double y) {
        return y < 0 ? -Math.abs(x) : Math.abs(x);
    }

    // Method to print a matrix
    public static void printMatrix(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                System.out.print(matrix[i][j] + " ");
            }
            System.out.println();
        }
    }

    // Method to calculate the determinant of a matrix using Gaussian elimination
    public static double determinant(double[][] matrix) {
        int N = matrix.length;
        double[][] A = new double[N][N];
        for (int i = 0; i < N; i++) {
            System.arraycopy(matrix[i], 0, A[i], 0, N);
        }

        double det = 1;
        for (int i = 0; i < N; i++) {
            int maxRow = i;
            for (int j = i + 1; j < N; j++) {
                if (Math.abs(A[j][i]) > Math.abs(A[maxRow][i])) {
                    maxRow = j;
                }
            }
            if (i != maxRow) {
                double[] temp = A[i];
                A[i] = A[maxRow];
                A[maxRow] = temp;
                det *= -1;
            }

            for (int j = i + 1; j < N; j++) {
                double factor = A[j][i] / A[i][i];
                for (int k = i; k < N; k++) {
                    A[j][k] -= factor * A[i][k];
                }
            }
            det *= A[i][i];
        }
        return det;
    }

    // Method to compute the Frobenius norm of a matrix
    public static double frobeniusNorm(double[][] matrix) {
        double norm = 0;
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                norm += matrix[i][j] * matrix[i][j];
            }
        }
        return Math.sqrt(norm);
    }

    // Method to compute the inverse of a matrix using Gaussian elimination
    public static double[][] inverse(double[][] matrix) {
        int N = matrix.length;
        double[][] A = new double[N][N];
        double[][] I = new double[N][N];

        // Create an identity matrix I
        for (int i = 0; i < N; i++) {
            I[i][i] = 1;
        }

        // Copy the original matrix to A
        for (int i = 0; i < N; i++) {
            System.arraycopy(matrix[i], 0, A[i], 0, N);
        }

        for (int i = 0; i < N; i++) {
            double pivot = A[i][i];
            if (pivot == 0) {
                throw new ArithmeticException("Matrix is singular and cannot be inverted.");
            }

            for (int j = 0; j < N; j++) {
                A[i][j] /= pivot;
                I[i][j] /= pivot;
            }

            for (int j = 0; j < N; j++) {
                if (i != j) {
                    double factor = A[j][i];
                    for (int k = 0; k < N; k++) {
                        A[j][k] -= factor * A[i][k];
                        I[j][k] -= factor * I[i][k];
                    }
                }
            }
        }

        return I;
    }

    // Method to solve the system A * x = b
    public static double[] solveSystem(double[][] A, double[] b) {
        int N = A.length;
        double[][] inverseA = inverse(A);
        double[] x = new double[N];

        // Multiply the inverse of A by b to get x
        for (int i = 0; i < N; i++) {
            x[i] = 0;
            for (int j = 0; j < N; j++) {
                x[i] += inverseA[i][j] * b[j];
            }
        }
        return x;
    }
}
