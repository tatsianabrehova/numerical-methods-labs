import java.util.Arrays;

public class Main {

    public static Result gaussian(double[][] matrix_, double[] inhomogeneity) {
        int n = matrix_.length;
        double[][] matrix = new double[n][n];
        for (int i = 0; i < n; i++) {
            System.arraycopy(matrix_[i], 0, matrix[i], 0, n);
        }
        double[] b = Arrays.copyOf(inhomogeneity, inhomogeneity.length);

        int[] indexes = new int[n];
        for (int i = 0; i < n; i++) {
            indexes[i] = i;
        }

        int insertions = 0;

        // Forward elimination with partial pivoting
        for (int k = 0; k < n; k++) {
            int leadingRow = k;
            for (int i = k + 1; i < n; i++) {
                if (Math.abs(matrix[i][k]) > Math.abs(matrix[leadingRow][k])) {
                    leadingRow = i;
                }
            }

            // Swap rows in matrix and b vector
            if (leadingRow != k) {
                double[] tempRow = matrix[k];
                matrix[k] = matrix[leadingRow];
                matrix[leadingRow] = tempRow;

                double tempB = b[k];
                b[k] = b[leadingRow];
                b[leadingRow] = tempB;

                insertions++;
            }

            for (int i = k + 1; i < n; i++) {
                double factor = matrix[i][k] / matrix[k][k];
                for (int j = k; j < n; j++) {
                    matrix[i][j] -= factor * matrix[k][j];
                }
                b[i] -= factor * b[k];
            }
        }

        // Back substitution
        double[] results = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0;
            for (int j = i + 1; j < n; j++) {
                sum += matrix[i][j] * results[j];
            }
            results[i] = (b[i] - sum) / matrix[i][i];
        }

        // Compute determinant
        double determinant = (insertions % 2 == 0 ? 1 : -1);
        for (int i = 0; i < n; i++) {
            determinant *= matrix[i][i];
        }

        // Compute discrepancy vector
        double[] discrepancyVector = new double[n];
        for (int i = 0; i < n; i++) {
            double ax = 0;
            for (int j = 0; j < n; j++) {
                ax += matrix_[i][j] * results[j];
            }
            discrepancyVector[i] = ax - inhomogeneity[i];
        }

        // Compute inverse matrix
        double[][] inverseMatrix = computeInverse(matrix_);

        // Compute discrepancy of inverse matrix
        double[][] identity = identityMatrix(n);
        double[][] inverseDiscrepancy = matrixSubtract(matrixMultiply(matrix_, inverseMatrix), identity);

        // Compute condition number
        double normA = matrixNorm(matrix_);
        double normAInv = matrixNorm(inverseMatrix);
        double conditionNumber = normA * normAInv;

        return new Result(results, discrepancyVector, determinant, inverseMatrix, inverseDiscrepancy, conditionNumber);
    }

    private static double[][] computeInverse(double[][] matrix) {
        int n = matrix.length;
        double[][] augmented = new double[n][2 * n];
        for (int i = 0; i < n; i++) {
            System.arraycopy(matrix[i], 0, augmented[i], 0, n);
            augmented[i][n + i] = 1.0;
        }

        for (int k = 0; k < n; k++) {
            for (int i = 0; i < n; i++) {
                if (i != k) {
                    double factor = augmented[i][k] / augmented[k][k];
                    for (int j = 0; j < 2 * n; j++) {
                        augmented[i][j] -= factor * augmented[k][j];
                    }
                }
            }
            double diagFactor = augmented[k][k];
            for (int j = 0; j < 2 * n; j++) {
                augmented[k][j] /= diagFactor;
            }
        }

        double[][] inverse = new double[n][n];
        for (int i = 0; i < n; i++) {
            System.arraycopy(augmented[i], n, inverse[i], 0, n);
        }
        return inverse;
    }

    private static double[][] matrixMultiply(double[][] a, double[][] b) {
        int rows = a.length;
        int cols = b[0].length;
        int sharedDim = a[0].length;
        double[][] result = new double[rows][cols];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                for (int k = 0; k < sharedDim; k++) {
                    result[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        return result;
    }

    private static double[][] matrixSubtract(double[][] a, double[][] b) {
        int rows = a.length;
        int cols = a[0].length;
        double[][] result = new double[rows][cols];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[i][j] = a[i][j] - b[i][j];
            }
        }
        return result;
    }

    private static double[][] identityMatrix(int size) {
        double[][] identity = new double[size][size];
        for (int i = 0; i < size; i++) {
            identity[i][i] = 1.0;
        }
        return identity;
    }

    private static double matrixNorm(double[][] matrix) {
        double maxSum = 0;
        for (double[] row : matrix) {
            double rowSum = 0;
            for (double value : row) {
                rowSum += Math.abs(value);
            }
            maxSum = Math.max(maxSum, rowSum);
        }
        return maxSum;
    }

    public static class Result {
        public final double[] solutions;
        public final double[] discrepancy;
        public final double determinant;
        public final double[][] inverseMatrix;
        public final double[][] inverseDiscrepancy;
        public final double conditionNumber;

        public Result(double[] solutions, double[] discrepancy, double determinant, double[][] inverseMatrix, double[][] inverseDiscrepancy, double conditionNumber) {
            this.solutions = solutions;
            this.discrepancy = discrepancy;
            this.determinant = determinant;
            this.inverseMatrix = inverseMatrix;
            this.inverseDiscrepancy = inverseDiscrepancy;
            this.conditionNumber = conditionNumber;
        }

        @Override
        public String toString() {
            return "Solutions: " + Arrays.toString(solutions) + ", Discrepancy: " + Arrays.toString(discrepancy) + ", Determinant: " + determinant +
                    ", Condition Number: " + conditionNumber + ", Inverse Matrix: " + Arrays.deepToString(inverseMatrix) + ", Inverse Discrepancy: " + Arrays.deepToString(inverseDiscrepancy);
        }
    }

    public static void main(String[] args) {
        double[][] matrix = {
                {2, 5, 4, 1},
                {1, 3, 2, 1},
                {2, 10, 9, 7},
                {3, 8, 9, 2}
        };

        double[] inhomogeneity = {20, 11, 40, 37};

        Result result = gaussian(matrix, inhomogeneity);
        System.out.println(result);
    }
}
