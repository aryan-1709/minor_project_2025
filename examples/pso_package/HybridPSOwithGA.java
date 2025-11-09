package pso_package;

import org.cloudbus.cloudsim.*;
import org.cloudbus.cloudsim.core.CloudSim;

import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Hybrid PSO + GA scheduler for CloudSim task assignment.
 * - PSO stage uses Local Best Murmuration PSO (PSOMbest) velocity update (Twumasi et al., 2024).
 * - GA stage follows the research article "Multivalent Optimizer-Based Hybrid Genetic Algorithm
 *   for Task Scheduling in Cloud Applications" (Malik et al., 2025) for GA operators.
 *
 * Prints:
 *  - Iteration progress (estimated best makespan)
 *  - Final cloudlet -> VM mapping
 *  - Actual simulation makespan from CloudSim
 */
public class HybridPSOwithGA {

    // ---------- Problem data ----------
    private static double[][] commMatrix;
    private static double[][] execMatrix;

    // ---------- CloudSim entities ----------
    private static List<Cloudlet> cloudletList;
    private static List<Vm> vmList;
    private static Datacenter[] datacenter;

    // ---------- Swarm / state ----------
    private static NewParticle[] swarm;        // Use your NewParticle or Particle type
    private static int[] gBestPosition;
    private static double gBestFitness;

    // ---------- PSO (PSOMbest) params (Twumasi et al., 2024) ----------
    private static final int NUM_PARTICLES = 150;
    private static final int NUM_CLUSTERS  = 5;
    private static final int MAX_ITERATIONS = 1000;   // you can lower to 500 for speed
    private static double W = 1.0;
    private static final double WDAMP = 0.99;
    private static final double C1 = 1.5;
    private static final double C2 = 2.0;
    private static final double C3 = 1.0;

    // ---------- GA params (Malik et al., 2025) ----------
    private static final double PCROSS = 0.8;     // two-point crossover probability
    private static final double PMUT   = 0.2;     // mutation probability (two mutation points)
    private static final int    GA_PERIOD = 50;   // apply GA every 50 PSO iterations
    private static final int    NO_IMPROVE_LIMIT = 100; // early stop if no improvement

    private static final Random rand = new Random();


    // =========================================================================================
    // Entry point
    // =========================================================================================
    public static void main(String[] args) {
        Log.printLine("Starting Hybrid PSO + GA Scheduler...");

        // 1) Build matrices
        new GenerateMatrices();
        execMatrix = GenerateMatrices.getExecMatrix();
        commMatrix = GenerateMatrices.getCommMatrix();

        // 2) Optimize by Hybrid PSO + GA
        int[] best = runHybrid();

        // 3) Simulate in CloudSim using best mapping
        runCloudSim(best);
    }


    // =========================================================================================
    // Hybrid PSO (PSOMbest) + GA
    // =========================================================================================
    private static int[] runHybrid() {
        Log.printLine("Running PSO(Mbest) + GA ...");

        // Initialize swarm
        swarm = new NewParticle[NUM_PARTICLES];
        gBestFitness = Double.MAX_VALUE;
        gBestPosition = new int[Constants.NO_OF_TASKS];

        for (int i = 0; i < NUM_PARTICLES; i++) {
            swarm[i] = new NewParticle(Constants.NO_OF_TASKS, Constants.NO_OF_DATA_CENTERS);
            swarm[i].evaluateFitness(commMatrix, execMatrix);
            if (swarm[i].getFitness() < gBestFitness) {
                gBestFitness = swarm[i].getFitness();
                gBestPosition = swarm[i].getPosition().clone();
            }
        }

        int noImprove = 0;

        for (int iter = 0; iter < MAX_ITERATIONS; iter++) {

            // ------ K-means clustering to get murmuration groups ------
            Map<Integer, List<NewParticle>> clusters = kMeansCluster(swarm, NUM_CLUSTERS);

            // ------ find local best (Mbest) per cluster ------
            Map<Integer, NewParticle> localBest = new HashMap<>();
            for (Map.Entry<Integer, List<NewParticle>> e : clusters.entrySet()) {
                NewParticle best = e.getValue().stream()
                        .min(Comparator.comparingDouble(NewParticle::getFitness))
                        .orElse(null);
                if (best != null) localBest.put(e.getKey(), best);
            }

            // ------ PSO (PSOMbest) velocity + position update ------
            for (NewParticle p : swarm) {
                int cid = assignCluster(p, clusters);
                NewParticle mBest = localBest.get(cid);

                double[] v = p.getVelocity();
                int[] x = p.getPosition();
                int[] pBest = p.getPersonalBest();

                for (int d = 0; d < x.length; d++) {
                    double r1 = rand.nextDouble();
                    double r2 = rand.nextDouble();
                    double r3 = rand.nextDouble();

                    double newV = W * v[d]
                            + C1 * r1 * (pBest[d]       - x[d])
                            + C2 * r2 * (gBestPosition[d] - x[d])
                            + C3 * r3 * (mBest.getPosition()[d] - x[d]);

                    v[d] = newV;

                    int newX = x[d] + (int) Math.round(newV);
                    if (newX < 0) newX = 0;
                    if (newX >= Constants.NO_OF_DATA_CENTERS)
                        newX = Constants.NO_OF_DATA_CENTERS - 1;

                    x[d] = newX;
                }

                p.setVelocity(v);
                p.setPosition(x);
                p.evaluateFitness(commMatrix, execMatrix);

                if (p.getFitness() < p.getPersonalBestFitness()) {
                    p.setPersonalBest(x.clone());
                    p.setPersonalBestFitness(p.getFitness());
                }
                if (p.getFitness() < gBestFitness) {
                    gBestFitness = p.getFitness();
                    gBestPosition = x.clone();
                    noImprove = 0;
                }
            }

            // Inertia damping
            W *= WDAMP;

            // ------ Periodic GA injection (Malik et al., 2025) ------
            if (iter > 0 && iter % GA_PERIOD == 0) {
                applyGAtoSwarm();
                // Refresh global best after GA
                for (NewParticle p : swarm) {
                    p.evaluateFitness(commMatrix, execMatrix);
                    if (p.getFitness() < p.getPersonalBestFitness()) {
                        p.setPersonalBest(p.getPosition().clone());
                        p.setPersonalBestFitness(p.getFitness());
                    }
                    if (p.getFitness() < gBestFitness) {
                        gBestFitness = p.getFitness();
                        gBestPosition = p.getPosition().clone();
                        noImprove = 0;
                    }
                }
                Log.printLine("Applied GA at iteration " + iter + " | Best Makespan = " + gBestFitness);
            }

            // Logging
            if (iter % 50 == 0) {
                Log.printLine("Iteration " + iter + " | Best Makespan = " + gBestFitness);
            }

            // Early stopping if no improvement
            noImprove++;
            if (noImprove >= NO_IMPROVE_LIMIT) {
                Log.printLine("Early stop: no improvement for " + NO_IMPROVE_LIMIT + " iterations.");
                break;
            }
        }

        Log.printLine("Hybrid completed. Estimated Best Makespan: " + gBestFitness);
        return gBestPosition;
    }


    // =========================================================================================
    // GA: selection (top-50%), two-point crossover (Pc=0.8), mutation (Pm=0.2 with two points)
    // =========================================================================================
    private static void applyGAtoSwarm() {
        // Sort by fitness ascending
        Arrays.sort(swarm, Comparator.comparingDouble(NewParticle::getFitness));

        int half = swarm.length / 2;
        List<NewParticle> parents = new ArrayList<>(Arrays.asList(swarm).subList(0, half));
        List<NewParticle> offspring = new ArrayList<>();

        while (offspring.size() < half) {
            // Select two random parents from top half
            NewParticle p1 = parents.get(rand.nextInt(parents.size()));
            NewParticle p2 = parents.get(rand.nextInt(parents.size()));

            int[] c1 = p1.getPosition().clone();
            int[] c2 = p2.getPosition().clone();

            if (rand.nextDouble() < PCROSS) {
                // two-point crossover
                int n = c1.length;
                int a = rand.nextInt(n);
                int b = rand.nextInt(n);
                if (a > b) { int tmp = a; a = b; b = tmp; }
                for (int i = a; i <= b; i++) {
                    int tmp = c1[i];
                    c1[i] = c2[i];
                    c2[i] = tmp;
                }
            }

            // mutation with two mutation points (Mp1, Mp2)
            mutateTwoPoints(c1);
            mutateTwoPoints(c2);

            // build children particles (need NewParticle(int[] pos, int numVMs))
            NewParticle ch1 = new NewParticle(c1, Constants.NO_OF_DATA_CENTERS);
            NewParticle ch2 = new NewParticle(c2, Constants.NO_OF_DATA_CENTERS);
            ch1.evaluateFitness(commMatrix, execMatrix);
            ch2.evaluateFitness(commMatrix, execMatrix);
            offspring.add(ch1);
            if (offspring.size() < half) offspring.add(ch2);
        }

        // Replace worst half with offspring
        for (int i = 0; i < half; i++) {
            swarm[swarm.length - 1 - i] = offspring.get(i);
        }
    }

    private static void mutateTwoPoints(int[] chrom) {
        if (rand.nextDouble() < PMUT) {
            int n = chrom.length;
            int mp1 = rand.nextInt(n);
            int mp2 = rand.nextInt(n);
            // reassign task->VM to random VM at two points
            chrom[mp1] = rand.nextInt(Constants.NO_OF_DATA_CENTERS);
            chrom[mp2] = rand.nextInt(Constants.NO_OF_DATA_CENTERS);
        }
    }


    // =========================================================================================
    // CloudSim flow (same style/printing as your PSOScheduler)
    // =========================================================================================
    private static void runCloudSim(int[] mapping) {
        Log.printLine("Initializing CloudSim with Hybrid schedule...");
        try {
            CloudSim.init(1, Calendar.getInstance(), false);

            // Datacenters
            datacenter = new Datacenter[Constants.NO_OF_DATA_CENTERS];
            for (int i = 0; i < Constants.NO_OF_DATA_CENTERS; i++) {
                datacenter[i] = DatacenterCreator.createDatacenter("Datacenter_" + i);
            }

            // Broker
            DatacenterBroker broker = new DatacenterBroker("Hybrid_Broker");
            int brokerId = broker.getId();

            // VMs and Cloudlets
            vmList = createVM(brokerId, Constants.NO_OF_DATA_CENTERS);
            cloudletList = createCloudlet(brokerId, Constants.NO_OF_TASKS, 0, mapping);

            broker.submitVmList(vmList);
            broker.submitCloudletList(cloudletList);

            CloudSim.startSimulation();
            List<Cloudlet> received = broker.getCloudletReceivedList();
            CloudSim.stopSimulation();

            printCloudletList(received);

        } catch (Exception e) {
            e.printStackTrace();
            Log.printLine("Simulation terminated due to an unexpected error.");
        }
    }

    private static List<Vm> createVM(int userId, int vms) {
        List<Vm> list = new ArrayList<>();
        long size = 10000;
        int ram = 512;
        int mips = 250;
        long bw = 1000;
        int pesNumber = 1;
        String vmm = "Xen";

        for (int i = 0; i < vms; i++) {
            Vm vm = new Vm(i, userId, mips, pesNumber, ram, bw, size, vmm,
                    new CloudletSchedulerSpaceShared());
            list.add(vm);
        }
        return list;
    }

    private static List<Cloudlet> createCloudlet(int userId, int cloudlets, int idShift, int[] mapping) {
        List<Cloudlet> list = new ArrayList<>();
        long fileSize = 300, outputSize = 300;
        int pesNumber = 1;
        UtilizationModel um = new UtilizationModelFull();

        for (int i = 0; i < cloudlets; i++) {
            int vmIndex = mapping[i];
            int vmId = vmList.get(vmIndex).getId();
            long length = (long) (1e3 * (execMatrix[i][vmIndex] + commMatrix[i][vmIndex]));
            Cloudlet c = new Cloudlet(idShift + i, length, pesNumber, fileSize, outputSize, um, um, um);
            c.setUserId(userId);
            c.setVmId(vmId);
            list.add(c);
        }
        return list;
    }

    private static void printCloudletList(List<Cloudlet> list) {
        String indent = "    ";
        DecimalFormat dft = new DecimalFormat("###.##");
        Log.printLine();
        Log.printLine("========== HYBRID PSO+GA OUTPUT ==========");
        Log.printLine("Cloudlet ID" + indent + "STATUS" + indent +
                "Data center ID" + indent + "VM ID" + indent + "Time" +
                indent + "Start Time" + indent + "Finish Time" + indent + "Waiting Time");

        for (Cloudlet cl : list) {
            if (cl.getCloudletStatus() == Cloudlet.SUCCESS) {
                Log.print(indent + cl.getCloudletId() + indent + "SUCCESS" + indent +
                        cl.getResourceId() + indent + indent + cl.getVmId() + indent + indent +
                        dft.format(cl.getActualCPUTime()) + indent +
                        dft.format(cl.getExecStartTime()) + indent +
                        dft.format(cl.getFinishTime()) + indent +
                        dft.format(cl.getWaitingTime()));
                Log.printLine();
            }
        }

        // Actual makespan from simulation
        double makespan = 0;
        for (Cloudlet c : list) makespan = Math.max(makespan, c.getFinishTime());
        Log.printLine("Actual Makespan (Simulation): " + makespan);

        // Also print final mapping for clarity
        Log.printLine("\n=== Final Cloudlet to VM Mapping (gBest) ===");
        for (int i = 0; i < gBestPosition.length; i++) {
            Log.printLine("Cloudlet " + i + " -> VM " + gBestPosition[i]);
        }
        Log.printLine("Estimated Best Makespan (from optimizer): " + gBestFitness);
    }


    // =========================================================================================
    // K-means helpers (murmuration clustering for PSOMbest)
    // =========================================================================================
    private static Map<Integer, List<NewParticle>> kMeansCluster(NewParticle[] swarm, int k) {
        Map<Integer, List<NewParticle>> clusters = new HashMap<>();
        List<NewParticle> centroids = new ArrayList<>();
        for (int i = 0; i < k; i++) centroids.add(swarm[rand.nextInt(swarm.length)]);

        boolean changed = true;
        int iter = 0, maxIter = 100;

        while (changed && iter < maxIter) {
            changed = false;
            clusters.clear();
            for (int i = 0; i < k; i++) clusters.put(i, new ArrayList<>());

            for (NewParticle p : swarm) {
                int nearest = getNearestCentroid(p, centroids);
                clusters.get(nearest).add(p);
            }

            List<NewParticle> newCentroids = new ArrayList<>();
            for (int i = 0; i < k; i++) {
                NewParticle best = clusters.get(i).stream()
                        .min(Comparator.comparingDouble(NewParticle::getFitness))
                        .orElse(centroids.get(i));
                newCentroids.add(best);
            }

            for (int i = 0; i < k; i++) {
                if (newCentroids.get(i) != centroids.get(i)) {
                    changed = true; break;
                }
            }
            centroids = newCentroids;
            iter++;
        }
        return clusters;
    }

    private static int getNearestCentroid(NewParticle p, List<NewParticle> centroids) {
        double minDist = Double.MAX_VALUE;
        int idx = 0;
        for (int i = 0; i < centroids.size(); i++) {
            double d = euclideanDistance(p.getPosition(), centroids.get(i).getPosition());
            if (d < minDist) { minDist = d; idx = i; }
        }
        return idx;
    }

    private static int assignCluster(NewParticle p, Map<Integer, List<NewParticle>> clusters) {
        int bestCluster = 0;
        double minDist = Double.MAX_VALUE;
        for (Map.Entry<Integer, List<NewParticle>> e : clusters.entrySet()) {
            for (NewParticle c : e.getValue()) {
                double d = euclideanDistance(p.getPosition(), c.getPosition());
                if (d < minDist) { minDist = d; bestCluster = e.getKey(); }
            }
        }
        return bestCluster;
    }

    private static double euclideanDistance(int[] a, int[] b) {
        double sum = 0.0;
        for (int i = 0; i < a.length; i++) {
            double diff = a[i] - b[i];
            sum += diff * diff;
        }
        return Math.sqrt(sum);
    }
}
