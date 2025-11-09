package pso_package;

import java.util.*;

public class AntColonyScheduler {

    // Compute makespan of a permutation (same as SA)
    static int computeMakespan(List<Integer> perm, double[][] execMatrix) {
        int m = execMatrix[0].length;
        PriorityQueue<int[]> pq = new PriorityQueue<>(Comparator.comparingInt(a -> a[0]));
        for (int i = 0; i < m; i++) pq.add(new int[]{0, i});

        for (int task : perm) {
            int[] top = pq.poll();
            int start = top[0];
            int finish = (int) (start + execMatrix[task][top[1]]);
            pq.add(new int[]{finish, top[1]});
        }

        int makespan = 0;
        while (!pq.isEmpty()) makespan = Math.max(makespan, pq.poll()[0]);
        return makespan;
    }

    // Main ACO logic
    public static int[] runACO(double[][] execMatrix, double[][] commMatrix) {
        int n = execMatrix.length; // tasks
        int m = execMatrix[0].length; // machines

        int numAnts = 20;
        int iterations = 100;
        double alpha = 1.0;
        double beta = 2.0;
        double evaporation = 0.5;
        double Q = 500.0;

        double[][] pheromone = new double[n][m];
        for (int i = 0; i < n; i++) Arrays.fill(pheromone[i], 1.0);

        int[] bestSol = new int[n];
        int bestCost = Integer.MAX_VALUE;

        Random rand = new Random();

        for (int it = 0; it < iterations; it++) {
            List<int[]> antSolutions = new ArrayList<>();
            List<Integer> antCosts = new ArrayList<>();

            for (int ant = 0; ant < numAnts; ant++) {
                int[] assign = new int[n];
                for (int task = 0; task < n; task++) {
                    double[] probs = new double[m];
                    double sum = 0.0;
                    for (int mach = 0; mach < m; mach++) {
                        double eta = 1.0 / (execMatrix[task][mach] + 1e-6);
                        probs[mach] = Math.pow(pheromone[task][mach], alpha) * Math.pow(eta, beta);
                        sum += probs[mach];
                    }

                    double r = rand.nextDouble() * sum;
                    double cumulative = 0.0;
                    int chosen = 0;
                    for (int mach = 0; mach < m; mach++) {
                        cumulative += probs[mach];
                        if (r <= cumulative) {
                            chosen = mach;
                            break;
                        }
                    }
                    assign[task] = chosen;
                }

                // Convert assignment → permutation for makespan calc
                List<Integer> perm = new ArrayList<>();
                for (int i = 0; i < n; i++) perm.add(i);
                int cost = computeMakespan(perm, execMatrix);
                antSolutions.add(assign);
                antCosts.add(cost);

                if (cost < bestCost) {
                    bestCost = cost;
                    bestSol = Arrays.copyOf(assign, n);
                }
            }

            // Evaporate pheromone
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    pheromone[i][j] *= (1 - evaporation);
                }
            }

            // Deposit pheromone
            for (int a = 0; a < numAnts; a++) {
                int[] sol = antSolutions.get(a);
                double contrib = Q / antCosts.get(a);
                for (int t = 0; t < n; t++) {
                    pheromone[t][sol[t]] += contrib;
                }
            }

            System.out.println("Iteration " + it + " Best Makespan: " + bestCost);
        }

        System.out.println("\nACO Optimization Completed. Best Makespan = " + bestCost);
        return bestSol;
    }
}
