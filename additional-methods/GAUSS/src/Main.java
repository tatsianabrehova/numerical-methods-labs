import java.util.Arrays;
import java.util.Scanner;
public class Main {

    // Округление до 3
    public static double roundToThreeDecimals(double value) {
        return Math.round(value * 1000.0) / 1000.0;
    }

    // округление всех элементов матрицы до 3 знаков после запятой
    public static void roundMatrix(double[][] matrix) {
        int rows = matrix.length;
        int cols = matrix[0].length;

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                matrix[i][j] = roundToThreeDecimals(matrix[i][j]);
            }
        }
    }

    // метод Гаусса
    public static double gaussElimination(double[][] matrix) {
        int rows = matrix.length;
        int cols = matrix[0].length;
        double determinant = 1.0; // Начальное значение определителя

        // Прямой ход
        for (int i = 0; i < rows; i++) {
            // деление текущей строки на ведущий элемент
            double mainElement = matrix[i][i];
            determinant *= mainElement;
            for (int j = i; j < cols; j++) {
                matrix[i][j] /= mainElement;
                matrix[i][j] = roundToThreeDecimals(matrix[i][j]);
            }

            // зануление элементов
            for (int k = i + 1; k < rows; k++) {
                double minus = matrix[k][i];
                for (int j = i; j < cols; j++) {
                    matrix[k][j] -= minus * matrix[i][j];
                    matrix[k][j] = roundToThreeDecimals(matrix[k][j]);
                }
            }
        }
        return roundToThreeDecimals(determinant);
    }

    // Функция нахождения решения после метода Гаусса
    public static double[] findSolution(double[][] matrix) {
        int rows = matrix.length;
        double[] solution = new double[rows];

        // Обратный ход для нахождения решений
        for (int i = rows - 1; i >= 0; i--) {
            solution[i] = matrix[i][rows];  // Последний столбец содержит b
            for (int j = i + 1; j < rows; j++) {
                solution[i] -= matrix[i][j] * solution[j];
            }
            solution[i] = roundToThreeDecimals(solution[i]);
        }

        return solution;
    }

    // Функция для вычисления невязки r = b - Ax
    public static double[] calculateResidual(double[][] originalMatrix, double[] solution, double[] rhs) {
        int n = originalMatrix.length;
        double[] residual = new double[n];

        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = 0; j < n; j++) {
                sum += originalMatrix[i][j] * solution[j];
            }
            residual[i] = roundToThreeDecimals(rhs[i] - sum); // r = b - Ax
        }

        return residual;
    }





    // Вывод
    public static void printMatrix(double[][] matrix) {
        for (double[] row : matrix) {
            for (double value : row) {
                System.out.printf("%8.3f", value);
            }
            System.out.println();
        }
    }
    // Вывод вектора
    public static void printVector(double[] vector) {
        for (double value : vector) {
            System.out.printf("%8.3f", value);
            System.out.println();
        }
    }

    // Функция для вычисления обратной матрицы
    public static double[][] inverseMatrix(double[][] matrix) {
        int n = matrix.length;
        double[][] augmentedMatrix = new double[n][2 * n]; // Расширенная матрица

        // Создание расширенной матрицы: исходная матрица слева, единичная матрица справа
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                augmentedMatrix[i][j] = matrix[i][j];
            }
            augmentedMatrix[i][i + n] = 1.0; // Единичная матрица
        }

        // Прямой ход метода Гаусса
        for (int i = 0; i < n; i++) {
            double mainElement = augmentedMatrix[i][i];

            // Делим строку на ведущий элемент
            for (int j = 0; j < 2 * n; j++) {
                augmentedMatrix[i][j] /= mainElement;
                augmentedMatrix[i][j] = roundToThreeDecimals(augmentedMatrix[i][j]);
            }

            // Зануление элементов ниже диагонали
            for (int k = i + 1; k < n; k++) {
                double factor = augmentedMatrix[k][i];
                for (int j = 0; j < 2 * n; j++) {
                    augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
                    augmentedMatrix[k][j] = roundToThreeDecimals(augmentedMatrix[k][j]);
                }
            }
        }

        // Обратный ход для зануления элементов выше диагонали
        for (int i = n - 1; i >= 0; i--) {
            for (int k = i - 1; k >= 0; k--) {
                double factor = augmentedMatrix[k][i];
                for (int j = 0; j < 2 * n; j++) {
                    augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
                    augmentedMatrix[k][j] = roundToThreeDecimals(augmentedMatrix[k][j]);
                }
            }
        }

        // Извлекаем правую часть матрицы (обратная матрица)
        double[][] inverseMatrix = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                inverseMatrix[i][j] = augmentedMatrix[i][j + n];
            }
        }

        return inverseMatrix;
    }


    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        System.out.print("Введите размер (n): ");
        int n = scanner.nextInt();
        // Создание и заполнение матрицы по формуле 1 / (i + j)
        double[][] matrix = new double[n][n+1];
        double[][] originalMatrix = new double[n][n]; // Для сохранения исходной матрицы
        double[] rhs = new double[n]; // Вектор правых частей
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                matrix[i][j] = 1.0 / (i + j + 1); // +1, чтобы избежать деления на 0
                originalMatrix[i][j] = matrix[i][j]; // Сохраняем оригинальную матрицу
            }
        }
        matrix[0][n]= 1;

        for (int i = 0; i < n; i++) {
            rhs[i] = matrix[i][n]; // Вектор правых частей
        }


        System.out.println("Ввод:");
        printMatrix(matrix);

        double determinant = gaussElimination(matrix);

        System.out.println("\nВывод:");
        printMatrix(matrix);

        // Нахождение решения
        double[] solution = findSolution(matrix);

        System.out.println("\nРешение:");
        printVector(solution);

        // Вычисление невязки
        double[] residual = calculateResidual(originalMatrix, solution, rhs);

        System.out.println("\nНевязка:");
        printVector(residual);

        System.out.println("\nОпределитель: " + determinant);

        double[][] inverse = inverseMatrix(matrix);

        System.out.println("\nОбратная матрица:");
        printMatrix(inverse);

    }
}
