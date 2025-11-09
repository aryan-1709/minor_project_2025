package pso_package;

import java.util.Random;

public class SimulatedAnnealingScheduler {

    private static final Random rand = new Random();

    private static final double T_START = 500.0;
    private static final double T_END = 0.01;
    private static final double COOL = 0.98;
    private static final int STEPS = 500;

    private static double[][] execMatrix;
    private static double[][] commMatrix;

    private static double computeMakespan(int[] sol) {
        double[] vmFinish = new double[Constants.NO_OF_DATA_CENTERS];

        for (int task = 0; task < Constants.NO_OF_TASKS; task++) {
            int vm = sol[task];
            vmFinish[vm] += execMatrix[task][vm] + commMatrix[task][vm];
        }

        double max = 0;
        for (double t : vmFinish) if (t > max) max = t;
        return max;
    }

    private static void randomNeighbor(int[] sol) {
        int task = rand.nextInt(Constants.NO_OF_TASKS);
        int oldVM = sol[task], newVM = oldVM;
        while (newVM == oldVM) newVM = rand.nextInt(Constants.NO_OF_DATA_CENTERS);
        sol[task] = newVM;
    }

    public static int[] runSA(double[][] execM, double[][] commM) {

        execMatrix = execM;
        commMatrix = commM;

        int[] curr = new int[Constants.NO_OF_TASKS];
        int[] best = new int[Constants.NO_OF_TASKS];

        for (int i = 0; i < Constants.NO_OF_TASKS; i++) {
            curr[i] = rand.nextInt(Constants.NO_OF_DATA_CENTERS);
            best[i] = curr[i];
        }

        double currCost = computeMakespan(curr);
        double bestCost = currCost;
        double T = T_START;

        while (T > T_END) {
            for (int s = 0; s < STEPS; s++) {
                int[] next = curr.clone();
                randomNeighbor(next);

                double nextCost = computeMakespan(next);
                if (nextCost < currCost || Math.exp(-(nextCost - currCost) / T) > rand.nextDouble()) {
                    curr = next;
                    currCost = nextCost;
                }

                if (currCost < bestCost) {
                    best = curr.clone();
                    bestCost = currCost;
                }
            }
            T *= COOL;
        }

        System.out.println("SA Finished. Best Estimated Makespan = " + bestCost);
        return best;
    }
}
