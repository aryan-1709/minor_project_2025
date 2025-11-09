package pso_package;
import java.util.*;

public class SimulatedAnnealingScheduler {

    // compute makespan of a given permutation
    static int computeMakespan(List<Integer> perm, List<Integer> p, int m) {
        // priority queue storing (finish_time, machine_id)
        PriorityQueue<int[]> pq = new PriorityQueue<>((a, b) -> a[0] - b[0]);
        for(int i = 0; i < m; i++) {
            pq.offer(new int[]{0, i});
        }

        for(int task : perm) {
            int[] top = pq.poll();
            int start = top[0];
            int finish = start + p.get(task);
            pq.offer(new int[]{finish, top[1]});
        }

        int makespan = 0;
        while(!pq.isEmpty()) {
            makespan = Math.max(makespan, pq.poll()[0]);
        }
        return makespan;
    }

    // generate neighbor by swapping two random tasks
    static void randomNeighbor(List<Integer> perm) {
        Random r = new Random();
        int n = perm.size();
        int i = r.nextInt(n), j = r.nextInt(n);
        while(j == i) j = r.nextInt(n);
        Collections.swap(perm, i, j);
    }

    // simulated annealing optimization
    static List<Integer> simulatedAnnealing(List<Integer> p, int m) {
        Random rand = new Random();
        int n = p.size();

        List<Integer> curr = new ArrayList<>();
        for(int i = 0; i < n; i++) curr.add(i);

        Collections.shuffle(curr);

        List<Integer> best = new ArrayList<>(curr);

        int currCost = computeMakespan(curr, p, m);
        int bestCost = currCost;

        double T = 500.0;
        double T_min = 0.01;
        double alpha = 0.98;
        int steps = 1000;

        while(T > T_min) {
            for(int s = 0; s < steps; s++) {
                List<Integer> next = new ArrayList<>(curr);
                randomNeighbor(next);
                int nextCost = computeMakespan(next, p, m);

                int diff = nextCost - currCost;

                if(diff < 0 || Math.exp(-diff / T) > rand.nextDouble()) {
                    curr = next;
                    currCost = nextCost;
                }

                if(currCost < bestCost) {
                    best = new ArrayList<>(curr);
                    bestCost = currCost;
                }
            }
            T *= alpha;
        }
        return best;
    }

    public static void main(String[] args) {
        Random rand = new Random();

        int n = 20; // number of tasks
        int m = 4;  // number of machines

        List<Integer> p = new ArrayList<>();
        for(int i = 0; i < n; i++) p.add(rand.nextInt(100) + 1);

        List<Integer> best = simulatedAnnealing(p, m);
        int bestCost = computeMakespan(best, p, m);

        System.out.println("Task processing times: " + p);
        System.out.println("Best Makespan: " + bestCost);
        System.out.println("Best Permutation: " + best);
    }
}
