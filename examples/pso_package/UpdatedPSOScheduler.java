//package pso_package;
//
//import org.cloudbus.cloudsim.*;
//import java.util.*;
//import java.util.stream.Collectors;
//
///**
// * Improved PSO Scheduler implementing the Local Best Murmuration Particle approach (PSOMbest)
// * strictly following:
// * "A novel improvement of particle swarm optimization using an improved velocity update
// * function based on local best murmuration particle" — Twumasi et al. (2024)
// */
//
//public class UpdatedPSOScheduler {
//	
//	/**public static void main(String[] args) {
//	    System.out.println("Initializing Local Best Murmuration PSO (PSOMbest) Test...");
//
//	    // Example dummy matrices for testing (replace with real values later)
//	    int numTasks = Constants.NO_OF_TASKS;
//	    int numVMs = Constants.NO_OF_DATA_CENTERS;
//
//	    double[][] execMatrix = new double[numTasks][numVMs];
//	    double[][] commMatrix = new double[numTasks][numVMs];
//
//	    Random rand = new Random();
//	    for (int i = 0; i < numTasks; i++) {
//	        for (int j = 0; j < numVMs; j++) {
//	            execMatrix[i][j] = 10 + rand.nextDouble() * 10;
//	            commMatrix[i][j] = rand.nextDouble();
//	        }
//	    }
//
//	    // Dummy lists for CloudSim (not used in this direct test)
//	    List<Cloudlet> cloudlets = new ArrayList<>();
//	    List<Vm> vms = new ArrayList<>();
//	    Datacenter[] datacenters = new Datacenter[1];
//
//	    runPSOMbest(cloudlets, vms, datacenters, commMatrix, execMatrix);
//	}*/
//	public static void main(String[] args) {
//	    System.out.println("=== Running Local Best Murmuration PSO (PSOMbest) Scheduler ===");
//
//	    // -------------------------------------------------------------
//	    // 1️⃣ Define simulation parameters (from Constants)
//	    // -------------------------------------------------------------
////	    Constants.NO_OF_TASKS = 10;         // Example: 10 Cloudlets
////	    Constants.NO_OF_DATA_CENTERS = 4;   // Example: 4 VMs
//
//	    int numTasks = Constants.NO_OF_TASKS;
//	    int numVMs = Constants.NO_OF_DATA_CENTERS;
//
//	    // -------------------------------------------------------------
//	    // 2️⃣ Create dummy execution and communication matrices
//	    // -------------------------------------------------------------
//	    double[][] execMatrix = new double[numTasks][numVMs];
//	    double[][] commMatrix = new double[numTasks][numVMs];
//
//	    Random rand = new Random();
//	    for (int i = 0; i < numTasks; i++) {
//	        for (int j = 0; j < numVMs; j++) {
//	            execMatrix[i][j] = 10 + rand.nextDouble() * 5;  // execution time between 10–15
//	            commMatrix[i][j] = rand.nextDouble() * 2;       // communication time between 0–2
//	        }
//	    }
//
//	    // -------------------------------------------------------------
//	    // 3️⃣ Dummy CloudSim entities (not used here, placeholders)
//	    // -------------------------------------------------------------
//	    List<Cloudlet> cloudlets = new ArrayList<>();
//	    List<Vm> vms = new ArrayList<>();
//	    Datacenter[] datacenters = new Datacenter[1];
//
//	    // -------------------------------------------------------------
//	    // 4️⃣ Run the PSOMbest algorithm
//	    // -------------------------------------------------------------
//	    runPSOMbest(cloudlets, vms, datacenters, commMatrix, execMatrix);
//
//	    // -------------------------------------------------------------
//	    // 5️⃣ Display final mapping of Cloudlets → VMs
//	    // -------------------------------------------------------------
//	    System.out.println("\n=== Final Cloudlet to VM Mapping (Global Best Position) ===");
//	    for (int i = 0; i < gBestPosition.length; i++) {
//	        System.out.printf("Cloudlet %d → VM %d%n", i, gBestPosition[i]);
//	    }
//
//	    // -------------------------------------------------------------
//	    // 6️⃣ Compute and print final makespan using best mapping
//	    // -------------------------------------------------------------
//	    double finalMakespan = calculateMakespan(gBestPosition, execMatrix, commMatrix);
//	    System.out.println("\nFinal Makespan (Best Schedule): " + finalMakespan);
//	}
//	
//	// --------------------------------------------------------------------
//	// Overloaded Makespan Calculation — with explicit matrices
//	// --------------------------------------------------------------------
//	private static double calculateMakespan(int[] solution, double[][] execMatrix, double[][] commMatrix) {
//	    double makespan = 0.0;
//	    double[] vmFinishTime = new double[Constants.NO_OF_DATA_CENTERS];
//
//	    for (int i = 0; i < Constants.NO_OF_TASKS; i++) {
//	        int vmId = solution[i];
//	        double execTime = execMatrix[i][vmId];
//	        double commTime = commMatrix[i][vmId];
//	        double cost = execTime + commTime;
//	        vmFinishTime[vmId] += cost;
//	    }
//
//	    for (double finishTime : vmFinishTime) {
//	        if (finishTime > makespan) {
//	            makespan = finishTime;
//	        }
//	    }
//
//	    return makespan;
//	}
//
//
//
//    private static List<Cloudlet> cloudletList;
//    private static List<Vm> vmList;
//    private static Datacenter[] datacenter;
//    private static double[][] commMatrix;
//    private static double[][] execMatrix;
//
//    private static NewParticle[] swarm;
//    private static int[] gBestPosition;
//    private static double gBestFitness;
//
//    private static final Random rand = new Random();
//
//    // --- PSOMbest Parameters (From Table 3 of the Paper) ---
//    private static final int MAX_ITERATIONS = 1000;
//    private static final int NUM_PARTICLES = 150;
//    private static final int NUM_CLUSTERS = 5;
//
//    private static double W = 1.0;
//    private static final double WDAMP = 0.99;
//    private static final double C1 = 1.5;
//    private static final double C2 = 2.0;
//    private static final double C3 = 1.0;
//
//    public static void runPSOMbest(List<Cloudlet> cloudlets, List<Vm> vms,
//                                   Datacenter[] datacenters,
//                                   double[][] commMat, double[][] execMat) {
//
//        cloudletList = cloudlets;
//        vmList = vms;
//        datacenter = datacenters;
//        commMatrix = commMat;
//        execMatrix = execMat;
//
//        swarm = new NewParticle[NUM_PARTICLES];
//        gBestFitness = Double.MAX_VALUE;
//
//        // --- Initialize swarm ---
//        for (int i = 0; i < NUM_PARTICLES; i++) {
//            swarm[i] = new NewParticle(Constants.NO_OF_TASKS, Constants.NO_OF_DATA_CENTERS);
//            swarm[i].evaluateFitness(commMatrix, execMatrix);
//            if (swarm[i].getFitness() < gBestFitness) {
//                gBestFitness = swarm[i].getFitness();
//                gBestPosition = swarm[i].getPosition().clone();
//            }
//        }
//
//        // --- Main PSO Loop ---
//        for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
//
//            // Step 1: Cluster particles using K-means
//            Map<Integer, List<NewParticle>> clusters = kMeansCluster(swarm, NUM_CLUSTERS);
//
//            // Step 2: Find local best (Mbest) for each cluster
//            Map<Integer, NewParticle> localBestMap = new HashMap<>();
//            for (Map.Entry<Integer, List<NewParticle>> entry : clusters.entrySet()) {
//            	NewParticle bestInCluster = entry.getValue().stream()
//                        .min(Comparator.comparingDouble(NewParticle::getFitness))
//                        .orElse(null);
//                if (bestInCluster != null)
//                    localBestMap.put(entry.getKey(), bestInCluster);
//            }
//
//            // Step 3: Update each particle’s velocity and position
//            for (int i = 0; i < NUM_PARTICLES; i++) {
//            	NewParticle p = swarm[i];
//                int clusterId = assignCluster(p, clusters);
//                NewParticle mBest = localBestMap.get(clusterId);
//
//                double[] velocity = p.getVelocity();
//                int[] position = p.getPosition();
//                int[] pBest = p.getPersonalBest();
//
//                for (int d = 0; d < position.length; d++) {
//                    double r1 = rand.nextDouble();
//                    double r2 = rand.nextDouble();
//                    double r3 = rand.nextDouble();
//
//                    double newVelocity = W * velocity[d]
//                            + C1 * r1 * (pBest[d] - position[d])
//                            + C2 * r2 * (gBestPosition[d] - position[d])
//                            + C3 * r3 * (mBest.getPosition()[d] - position[d]);
//                    velocity[d] = newVelocity;
//
//                    int newPos = position[d] + (int) Math.round(newVelocity);
//                    // Ensure position bounds
//                    if (newPos < 0) newPos = 0;
//                    if (newPos >= Constants.NO_OF_DATA_CENTERS)
//                        newPos = Constants.NO_OF_DATA_CENTERS - 1;
//
//                    position[d] = newPos;
//                }
//
//                // Update position and velocity
//                p.setVelocity(velocity);
//                p.setPosition(position);
//
//                // Evaluate new fitness
//                p.evaluateFitness(commMatrix, execMatrix);
//
//                // Update personal best
//                if (p.getFitness() < p.getPersonalBestFitness()) {
//                    p.setPersonalBest(position.clone());
//                    p.setPersonalBestFitness(p.getFitness());
//                }
//
//                // Update global best
//                if (p.getFitness() < gBestFitness) {
//                    gBestFitness = p.getFitness();
//                    gBestPosition = position.clone();
//                }
//            }
//
//            // Step 4: Apply inertia damping
//            W *= WDAMP;
//
//            // Optionally print progress
//            if (iter % 50 == 0)
//                System.out.println("Iteration " + iter + " | Best Fitness: " + gBestFitness);
//        }
//
//        System.out.println("Final Best Makespan: " + gBestFitness);
//    }
//
//    // --------------------------------------------------------------------
//    // --- Helper: K-Means clustering for murmuration groups ---
//    // --------------------------------------------------------------------
//    private static Map<Integer, List<NewParticle>> kMeansCluster(NewParticle[] swarm, int k) {
//        Map<Integer, List<NewParticle>> clusters = new HashMap<>();
//        // Initialize centroids randomly
//        List<NewParticle> centroids = new ArrayList<>();
//        for (int i = 0; i < k; i++) {
//            centroids.add(swarm[rand.nextInt(swarm.length)]);
//        }
//
//        boolean changed = true;
//        int maxIter = 100;
//        int iter = 0;
//
//        while (changed && iter < maxIter) {
//            changed = false;
//            clusters.clear();
//            for (int i = 0; i < k; i++) clusters.put(i, new ArrayList<>());
//
//            // Assign each particle to nearest centroid
//            for (NewParticle p : swarm) {
//                int nearest = getNearestCentroid(p, centroids);
//                clusters.get(nearest).add(p);
//            }
//
//            // Update centroids
//            List<NewParticle> newCentroids = new ArrayList<>();
//            for (int i = 0; i < k; i++) {
//                NewParticle best = clusters.get(i).stream()
//                        .min(Comparator.comparingDouble(NewParticle::getFitness))
//                        .orElse(centroids.get(i));
//                newCentroids.add(best);
//            }
//
//            // Check if centroids changed
//            for (int i = 0; i < k; i++) {
//                if (newCentroids.get(i) != centroids.get(i)) {
//                    changed = true;
//                    break;
//                }
//            }
//            centroids = newCentroids;
//            iter++;
//        }
//        return clusters;
//    }
//
//    private static int getNearestCentroid(NewParticle p, List<NewParticle> centroids) {
//        double minDist = Double.MAX_VALUE;
//        int idx = 0;
//        for (int i = 0; i < centroids.size(); i++) {
//            double dist = euclideanDistance(p.getPosition(), centroids.get(i).getPosition());
//            if (dist < minDist) {
//                minDist = dist;
//                idx = i;
//            }
//        }
//        return idx;
//    }
//
//    private static int assignCluster(NewParticle p, Map<Integer, List<NewParticle>> clusters) {
//        int bestCluster = 0;
//        double minDist = Double.MAX_VALUE;
//        for (Map.Entry<Integer, List<NewParticle>> entry : clusters.entrySet()) {
//            for (NewParticle c : entry.getValue()) {
//                double dist = euclideanDistance(p.getPosition(), c.getPosition());
//                if (dist < minDist) {
//                    minDist = dist;
//                    bestCluster = entry.getKey();
//                }
//            }
//        }
//        return bestCluster;
//    }
//
//    private static double euclideanDistance(int[] a, int[] b) {
//        double sum = 0.0;
//        for (int i = 0; i < a.length; i++) {
//            double diff = a[i] - b[i];
//            sum += diff * diff;
//        }
//        return Math.sqrt(sum);
//    }
//
//}


package pso_package;

import org.cloudbus.cloudsim.*;
import org.cloudbus.cloudsim.core.CloudSim;
import java.text.DecimalFormat;
import java.util.*;

/**
 * Improved PSO Scheduler implementing Local Best Murmuration PSO (PSOMbest)
 * based strictly on:
 * "A novel improvement of particle swarm optimization using an improved velocity update
 * function based on local best murmuration particle" — Twumasi et al. (2024)
 */
public class UpdatedPSOScheduler {

    private static List<Cloudlet> cloudletList;
    private static List<Vm> vmList;
    private static Datacenter[] datacenter;
    private static double[][] commMatrix;
    private static double[][] execMatrix;

    private static NewParticle[] swarm;
    private static int[] gBestPosition;
    private static double gBestFitness;
    private static final Random rand = new Random();

    // --- PSOMbest Parameters (Table 3 of the Paper) ---
    private static final int MAX_ITERATIONS = 500;
    private static final int NUM_PARTICLES = 150;
    private static final int NUM_CLUSTERS = 5;
    private static double W = 1.0;
    private static final double WDAMP = 0.99;
    private static final double C1 = 1.5;
    private static final double C2 = 2.0;
    private static final double C3 = 1.0;

    // -------------------------------------------------------------
    // 1️⃣ Main Entry Point
    // -------------------------------------------------------------
    public static void main(String[] args) {
        Log.printLine("Starting Local Best Murmuration PSO Scheduler...");

        // Step 1: Initialize execution and communication matrices
        new GenerateMatrices();
        execMatrix = GenerateMatrices.getExecMatrix();
        commMatrix = GenerateMatrices.getCommMatrix();

        // Step 2: Run PSOMbest to find optimal cloudlet→VM mapping
        int[] psoSolution = runPSOMbest();

        // Step 3: Simulate in CloudSim using PSOMbest schedule
        runCloudSim(psoSolution);
    }

    // -------------------------------------------------------------
    // 2️⃣ PSOMbest Algorithm Implementation
    // -------------------------------------------------------------
    private static int[] runPSOMbest() {
        Log.printLine("Running Local Best Murmuration PSO (PSOMbest)...");

        swarm = new NewParticle[NUM_PARTICLES];
        gBestFitness = Double.MAX_VALUE;
        gBestPosition = new int[Constants.NO_OF_TASKS];

        // --- Initialize swarm ---
        for (int i = 0; i < NUM_PARTICLES; i++) {
            swarm[i] = new NewParticle(Constants.NO_OF_TASKS, Constants.NO_OF_DATA_CENTERS);
            swarm[i].evaluateFitness(commMatrix, execMatrix);
            if (swarm[i].getFitness() < gBestFitness) {
                gBestFitness = swarm[i].getFitness();
                gBestPosition = swarm[i].getPosition().clone();
            }
        }

        // --- Main Iterative Optimization ---
        for (int iter = 0; iter < MAX_ITERATIONS; iter++) {

            // Step 1: Cluster particles into murmuration groups (K-means)
            Map<Integer, List<NewParticle>> clusters = kMeansCluster(swarm, NUM_CLUSTERS);

            // Step 2: Find local best (Mbest) per cluster
            Map<Integer, NewParticle> localBestMap = new HashMap<>();
            for (Map.Entry<Integer, List<NewParticle>> entry : clusters.entrySet()) {
                NewParticle best = entry.getValue().stream()
                        .min(Comparator.comparingDouble(NewParticle::getFitness))
                        .orElse(null);
                if (best != null)
                    localBestMap.put(entry.getKey(), best);
            }

            // Step 3: Update each particle’s velocity and position
            for (NewParticle p : swarm) {
                int clusterId = assignCluster(p, clusters);
                NewParticle mBest = localBestMap.get(clusterId);

                double[] velocity = p.getVelocity();
                int[] position = p.getPosition();
                int[] pBest = p.getPersonalBest();

                for (int d = 0; d < position.length; d++) {
                    double r1 = rand.nextDouble();
                    double r2 = rand.nextDouble();
                    double r3 = rand.nextDouble();

                    double newVelocity = W * velocity[d]
                            + C1 * r1 * (pBest[d] - position[d])
                            + C2 * r2 * (gBestPosition[d] - position[d])
                            + C3 * r3 * (mBest.getPosition()[d] - position[d]);

                    velocity[d] = newVelocity;

                    int newPos = position[d] + (int) Math.round(newVelocity);
                    if (newPos < 0) newPos = 0;
                    if (newPos >= Constants.NO_OF_DATA_CENTERS)
                        newPos = Constants.NO_OF_DATA_CENTERS - 1;
                    position[d] = newPos;
                }

                // Update particle state
                p.setVelocity(velocity);
                p.setPosition(position);
                p.evaluateFitness(commMatrix, execMatrix);

                // Update personal and global bests
                if (p.getFitness() < p.getPersonalBestFitness()) {
                    p.setPersonalBest(position.clone());
                    p.setPersonalBestFitness(p.getFitness());
                }
                if (p.getFitness() < gBestFitness) {
                    gBestFitness = p.getFitness();
                    gBestPosition = position.clone();
                }
            }

            // Step 4: Dampen inertia weight
            W *= WDAMP;

            // Step 5: Print iteration summary
            if (iter % 50 == 0)
                Log.printLine("Iteration " + iter + ": Best Makespan = " + gBestFitness);
        }

        Log.printLine("PSOMbest completed. Estimated Best Makespan: " + gBestFitness);
        return gBestPosition;
    }

    // -------------------------------------------------------------
    // 3️⃣ CloudSim Execution using PSOMbest Mapping
    // -------------------------------------------------------------
    private static void runCloudSim(int[] psoSolution) {
        Log.printLine("Initializing CloudSim Simulation with PSOMbest Solution...");

        try {
            int num_user = 1;
            Calendar calendar = Calendar.getInstance();
            boolean trace_flag = false;

            CloudSim.init(num_user, calendar, trace_flag);

            // Create datacenters
            datacenter = new Datacenter[Constants.NO_OF_DATA_CENTERS];
            for (int i = 0; i < Constants.NO_OF_DATA_CENTERS; i++) {
                datacenter[i] = DatacenterCreator.createDatacenter("Datacenter_" + i);
            }

            // Create Broker
            DatacenterBroker broker = new DatacenterBroker("PSO_Broker");
            int brokerId = broker.getId();

            // Create VMs and Cloudlets
            vmList = createVM(brokerId, Constants.NO_OF_DATA_CENTERS);
            cloudletList = createCloudlet(brokerId, Constants.NO_OF_TASKS, 0, psoSolution);

            broker.submitVmList(vmList);
            broker.submitCloudletList(cloudletList);

            // Start Simulation
            CloudSim.startSimulation();

            List<Cloudlet> newList = broker.getCloudletReceivedList();
            CloudSim.stopSimulation();

            // Print results
            printCloudletList(newList);

        } catch (Exception e) {
            e.printStackTrace();
            Log.printLine("Simulation terminated unexpectedly.");
        }
    }

    // -------------------------------------------------------------
    // Helper Methods for CloudSim Setup
    // -------------------------------------------------------------
    private static List<Vm> createVM(int userId, int vms) {
        List<Vm> list = new ArrayList<>();
        long size = 10000;
        int ram = 512;
        int mips = 250;
        long bw = 1000;
        int pesNumber = 1;
        String vmm = "Xen";

        for (int i = 0; i < vms; i++) {
            Vm vm = new Vm(i, userId, mips, pesNumber, ram, bw, size, vmm, new CloudletSchedulerSpaceShared());
            list.add(vm);
        }
        return list;
    }

    private static List<Cloudlet> createCloudlet(int userId, int cloudlets, int idShift, int[] psoSolution) {
        List<Cloudlet> list = new ArrayList<>();
        long fileSize = 300;
        long outputSize = 300;
        int pesNumber = 1;
        UtilizationModel utilizationModel = new UtilizationModelFull();

        for (int i = 0; i < cloudlets; i++) {
            int vmIndex = psoSolution[i];
            int vmId = vmList.get(vmIndex).getId();
            long length = (long) (1e3 * (execMatrix[i][vmIndex] + commMatrix[i][vmIndex]));
            Cloudlet cloudlet = new Cloudlet(idShift + i, length, pesNumber, fileSize, outputSize,
                    utilizationModel, utilizationModel, utilizationModel);
            cloudlet.setUserId(userId);
            cloudlet.setVmId(vmId);
            list.add(cloudlet);
        }
        return list;
    }

    private static void printCloudletList(List<Cloudlet> list) {
        String indent = "    ";
        DecimalFormat dft = new DecimalFormat("###.##");
        Log.printLine();
        Log.printLine("========== PSOMbest OUTPUT ==========");
        Log.printLine("Cloudlet ID" + indent + "STATUS" + indent + "Datacenter ID" +
                indent + "VM ID" + indent + "Time" + indent + "Start Time" +
                indent + "Finish Time" + indent + "Waiting Time");

        for (Cloudlet cloudlet : list) {
            if (cloudlet.getCloudletStatus() == Cloudlet.SUCCESS) {
                Log.print(indent + cloudlet.getCloudletId() + indent + "SUCCESS" + indent +
                        cloudlet.getResourceId() + indent + indent +
                        cloudlet.getVmId() + indent + indent +
                        dft.format(cloudlet.getActualCPUTime()) + indent +
                        dft.format(cloudlet.getExecStartTime()) + indent +
                        dft.format(cloudlet.getFinishTime()) + indent +
                        dft.format(cloudlet.getWaitingTime()));
                Log.printLine();
            }
        }

        // Calculate actual makespan
        double makespan = calcActualMakespan(list);
        Log.printLine("Actual Makespan (Simulation): " + makespan);
    }

    private static double calcActualMakespan(List<Cloudlet> list) {
        double makespan = 0;
        for (Cloudlet c : list) {
            if (c.getFinishTime() > makespan) makespan = c.getFinishTime();
        }
        return makespan;
    }

    // -------------------------------------------------------------
    // Helper Methods for K-Means and Distance
    // -------------------------------------------------------------
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
                    changed = true;
                    break;
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
            double dist = euclideanDistance(p.getPosition(), centroids.get(i).getPosition());
            if (dist < minDist) {
                minDist = dist;
                idx = i;
            }
        }
        return idx;
    }

    private static int assignCluster(NewParticle p, Map<Integer, List<NewParticle>> clusters) {
        int bestCluster = 0;
        double minDist = Double.MAX_VALUE;
        for (Map.Entry<Integer, List<NewParticle>> entry : clusters.entrySet()) {
            for (NewParticle c : entry.getValue()) {
                double dist = euclideanDistance(p.getPosition(), c.getPosition());
                if (dist < minDist) {
                    minDist = dist;
                    bestCluster = entry.getKey();
                }
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
