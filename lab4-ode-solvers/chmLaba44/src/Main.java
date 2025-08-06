import java.util.ArrayList;
import java.util.Arrays;

public class Main {

    // Точность вычисления
    private static final double EPSILON = 1e-4;

    // Начальные условия
    private static final double[] U0 = {0.2, 0.0};

    // Временной интервал
    private static final double T_START = 0.0;
    private static final double T_END = 1.0;

    // Метод main
    public static void main(String[] args) {
        System.out.println("=== Явный метод трапеций ===");
        runAlgorithm(Main::trapezoidMethod, 2);

        System.out.println("\n=== Метод Рунге–Кутты 3-го порядка ===");
        runAlgorithm(Main::rungeKuttaThirdOrder, 3);

        System.out.println("\n=== Метод Адамса 4-го порядка (с адаптивным шагом) ===");
        runAdamsFourthOrder();
    }

    // Обертка для запуска методов
    private static void runAlgorithm(Method method, int m) {
        Result result = integrate(method, U0, T_START, T_END, EPSILON, m);
        double[][] solution = result.solution;
        double[] nodes = result.nodes;

        System.out.println("Число узлов: " + nodes.length);
        System.out.println("Максимальная погрешность:");
        double maxError = maxError(solution, nodes);
        System.out.println(maxError);

        // Можно добавить вывод графика через JavaFX или сохранение в CSV
    }

    // Правая часть системы ОДУ
    private static double[] model(double[] u, double t) {
        double u1 = u[0];
        double u2 = u[1];

        double du1dt = -10 * u1 * u1 + 4 * u1 + 3 * u2 * u2 + 0.875;
        double du2dt = 11 * u2 * u2 + (-20 * u1 + 4) * u2 - 7.25;

        return new double[]{du1dt, du2dt};
    }

    // Аналитическое решение
    private static double[] realSolution(double t) {
        double y1 = 0.2 - 4.2 * Math.tanh(-7 * t) + 3.75 * Math.tanh(-7.5 * t);
        double y2 = -7.0 * Math.tanh(-7 * t) + 7.5 * Math.tanh(-7.5 * t);
        return new double[]{y1, y2};
    }

    // Формула явного метода трапеций
    private static double[] trapezoidMethod(double[] y, double t, double h, Function f) {
        double[] k1 = f.apply(y, t);
        double[] yTemp = new double[y.length];
        for (int i = 0; i < y.length; i++) {
            yTemp[i] = y[i] + h * k1[i];
        }
        double[] k2 = f.apply(yTemp, t + h);
        for (int i = 0; i < y.length; i++) {
            y[i] = y[i] + h / 2 * (k1[i] + k2[i]);
        }
        return y;
    }

    // Метод Рунге–Кутты 3-го порядка
    private static double[] rungeKuttaThirdOrder(double[] y, double t, double h, Function f) {
        double[] k0 = scale(f.apply(y, t), h);
        double[] y1 = add(y, scale(k0, 1.0 / 3));
        double[] k1 = scale(f.apply(y1, t + h / 3), h);
        double[] y2 = add(y, scale(k1, 2.0 / 3));
        double[] k2 = scale(f.apply(y2, t + 2 * h / 3), h);
        return add(y, add(scale(k0, 1.0 / 4), scale(k2, 3.0 / 4)));
    }

    // Адаптивный интегратор по правилу Рунге
    private static Result integrate(Method method, double[] y0, double tStart, double tEnd, double epsilon, int m) {
        ArrayList<double[]> solutionList = new ArrayList<>();
        ArrayList<Double> timeList = new ArrayList<>();

        double[] y = Arrays.copyOf(y0, y0.length);
        double t = tStart;

        solutionList.add(Arrays.copyOf(y, y.length));
        timeList.add(t);

        double h = Math.pow(epsilon, 1.0 / m); // начальный шаг
        double k = 1.05;

        while (t < tEnd) {
            boolean stepAccepted = false;
            while (!stepAccepted) {
                double tNext = t + h;
                if (tNext > tEnd) tNext = tEnd;

                double[] yFullStep = method.apply(copy(y), t, h, Main::model);
                double[] yHalfStep = method.apply(copy(y), t, h / 2, Main::model);
                yHalfStep = method.apply(yHalfStep, t + h / 2, h / 2, Main::model);

                double error = max(abs(subtract(yHalfStep, yFullStep)));

                if (error < epsilon) {
                    y = yHalfStep;
                    t = tNext;
                    solutionList.add(Arrays.copyOf(y, y.length));
                    timeList.add(t);
                    h *= k;
                    stepAccepted = true;
                } else {
                    h /= 2;
                }
            }
        }

        double[][] solutionArray = new double[solutionList.size()][];
        for (int i = 0; i < solutionArray.length; i++) {
            solutionArray[i] = solutionList.get(i);
        }

        double[] timeArray = new double[timeList.size()];
        for (int i = 0; i < timeArray.length; i++) {
            timeArray[i] = timeList.get(i);
        }

        return new Result(solutionArray, timeArray);
    }

    // --- Реализация метода Адамса 4-го порядка ---
    private static void runAdamsFourthOrder() {
        // Используем Рунге–Кутта 4-го порядка для первых 3 точек
        double[] y0 = {0.2, 0.0};
        double tStart = 0.0;
        double tEnd = 1.0;
        double h = 1e-4;

        double[][] rkPoints = new double[3][];
        double[] times = new double[3];

        double t = tStart;
        double[] y = Arrays.copyOf(y0, y0.length);
        rkPoints[0] = copy(y);
        times[0] = t;

        for (int i = 1; i < 3; i++) {
            y = rungeKuttaFourthOrder(y, t, h, Main::model);
            t += h;
            rkPoints[i] = copy(y);
            times[i] = t;
        }

        // Запуск метода Адамса
        ArrayList<double[]> adamsSolution = new ArrayList<>();
        ArrayList<Double> adamsTime = new ArrayList<>();

        for (double[] pt : rkPoints) adamsSolution.add(copy(pt));
        for (double tm : times) adamsTime.add(tm);

        double[] lastY = copy(rkPoints[2]);
        double lastT = times[2];

        double currentH = h;

        while (lastT < tEnd) {
            double[] nextY = adamsFourthOrder(lastY, lastT, currentH, adamsSolution, adamsTime, Main::model);
            lastT += currentH;
            adamsSolution.add(copy(nextY));
            adamsTime.add(lastT);
        }

        double[][] sol = new double[adamsSolution.size()][];
        for (int i = 0; i < sol.length; i++) {
            sol[i] = adamsSolution.get(i);
        }

        double[] nodes = new double[adamsTime.size()];
        for (int i = 0; i < nodes.length; i++) {
            nodes[i] = adamsTime.get(i);
        }

        System.out.println("Число узлов (метод Адамса): " + nodes.length);
        System.out.println("Максимальная погрешность (метод Адамса):");
        System.out.println(maxError(sol, nodes));
    }

    private static double[] adamsFourthOrder(double[] y, double t, double h, ArrayList<double[]> prevY, ArrayList<Double> prevT, Function f) {
        int n = prevY.size();
        double[] f0 = f.apply(prevY.get(n - 1), prevT.get(n - 1));
        double[] f1 = f.apply(prevY.get(n - 2), prevT.get(n - 2));
        double[] f2 = f.apply(prevY.get(n - 3), prevT.get(n - 3));

        double[] dy = new double[y.length];
        for (int i = 0; i < y.length; i++) {
            dy[i] = h / 24 * (9 * f0[i] + 19 * f1[i] - 5 * f2[i]);
        }

        double[] newY = add(y, dy);
        return newY;
    }

    // --- Вспомогательные функции ---

    @FunctionalInterface
    interface Function {
        double[] apply(double[] y, double t);
    }

    @FunctionalInterface
    interface Method {
        double[] apply(double[] y, double t, double h, Function f);
    }

    private static class Result {
        double[][] solution;
        double[] nodes;

        public Result(double[][] solution, double[] nodes) {
            this.solution = solution;
            this.nodes = nodes;
        }
    }

    private static double[] add(double[] a, double[] b) {
        double[] res = new double[a.length];
        for (int i = 0; i < a.length; i++) res[i] = a[i] + b[i];
        return res;
    }

    private static double[] subtract(double[] a, double[] b) {
        double[] res = new double[a.length];
        for (int i = 0; i < a.length; i++) res[i] = a[i] - b[i];
        return res;
    }

    private static double[] scale(double[] a, double s) {
        double[] res = new double[a.length];
        for (int i = 0; i < a.length; i++) res[i] = a[i] * s;
        return res;
    }

    private static double[] copy(double[] a) {
        return Arrays.copyOf(a, a.length);
    }

    private static double[] abs(double[] a) {
        double[] res = new double[a.length];
        for (int i = 0; i < a.length; i++) res[i] = Math.abs(a[i]);
        return res;
    }

    private static double max(double[] a) {
        double max = 0;
        for (double v : a) max = Math.max(max, Math.abs(v));
        return max;
    }

    private static double maxError(double[][] solution, double[] nodes) {
        double maxErr = 0;
        for (int i = 0; i < solution.length; i++) {
            double[] real = realSolution(nodes[i]);
            double[] diff = abs(subtract(real, solution[i]));
            maxErr = Math.max(maxErr, max(diff));
        }
        return maxErr;
    }

    private static double[] rungeKuttaFourthOrder(double[] y, double t, double h, Function f) {
        double[] k1 = scale(f.apply(y, t), h);
        double[] k2 = scale(f.apply(add(y, scale(k1, 0.5)), t + h / 2), h);
        double[] k3 = scale(f.apply(add(y, scale(k2, 0.5)), t + h / 2), h);
        double[] k4 = scale(f.apply(add(y, k3), t + h), h);
        return add(y, add(add(scale(k1, 1.0 / 6), scale(k2, 1.0 / 3)),
                add(scale(k3, 1.0 / 3), scale(k4, 1.0 / 6))));
    }
}