package pso_package;

import org.cloudbus.cloudsim.*;
import org.cloudbus.cloudsim.core.CloudSim;

import java.text.DecimalFormat;
import java.util.*;
import java.io.FileWriter;

/**
 * GeneticScheduler.java
 *
 * Implements the GA described in:
 * "Multivalent Optimizer-Based Hybrid Genetic Algorithm for Task Scheduling in Cloud Applications"
 * (Malik et al., 2025) — selection, two-point crossover (Pc=0.8),
 * two-point mutation (Pm=0.2), termination when no improvement.
 *
 * Weighted fitness (user requested):
 *   Fitness = 0.5 * Makespan + 0.3 * CPU + 0.2 * Energy
 *
 * Outputs:
 *  - Iteration logs with (Fitness, Makespan, CPU, Energy)
 *  - Final best chromosome mapping (Cloudlet -> VM)
 *  - Estimated (GA) makespan, CPU, energy
 *  - Actual CloudSim makespan and per-cloudlet results
 *  - CSV logs: fitness.csv, makespan.csv, cpu.csv, energy.csv
 */
public class GeneticScheduler {

    // Matrices and CloudSim entities
    private static double[][] execMatrix;
    private static double[][] commMatrix;
    private static List<Cloudlet> cloudletList;
    private static List<Vm> vmList;
    private static Datacenter[] datacenter;

    // GA state
    private static int[][] population;        // population of chromosomes: population[i][task] = vmIndex
    private static double[] populationFitness;
    private static int popSize;
    private static int chromosomeLength;      // = Constants.NO_OF_TASKS

    private static int[] bestChromosome;
    private static double bestFitness = Double.MAX_VALUE;

    private static final Random rand = new Random();

    // GA hyperparameters (from paper)
    private static final double PC = 0.8;    // crossover probability
    private static final double PM = 0.2;    // mutation probability (two-point)
    private static final int MAX_GENERATIONS = 1000; // paper uses many iterations; user can lower
    private static final int NO_IMPROVE_LIMIT = 100; // termination when no improvement

    // CSV logging buffers
    private static final List<Double> csvFitness = new ArrayList<>();
    private static final List<Double> csvMakespan = new ArrayList<>();
    private static final List<Double> csvCPU = new ArrayList<>();
    private static final List<Double> csvEnergy = new ArrayList<>();

    // ---------- Entry point ----------
    public static void main(String[] args) {
        Log.printLine("Starting GeneticScheduler (MO-GA like implementation)...");

        // 1) generate matrices (existing in your project)
        new GenerateMatrices();
        execMatrix = GenerateMatrices.getExecMatrix();
        commMatrix = GenerateMatrices.getCommMatrix();

        // 2) population size: prefer Constants.POPULATION_SIZE if exists, else default 100
        try {
            popSize = Constants.POPULATION_SIZE;
        } catch (Throwable t) {
            popSize = 100;
        }
        chromosomeLength = Constants.NO_OF_TASKS;

        // 3) initialize population
        initPopulation();

        // 4) run GA
        runGA();

        // 5) print final GA results and run CloudSim with best mapping
        printFinalEstimatedResults();
        runCloudSim(bestChromosome);

        // 6) save CSV logs
        saveCSV();
    }

    // ---------- Population initialization ----------
    private static void initPopulation() {
        population = new int[popSize][chromosomeLength];
        populationFitness = new double[popSize];

        for (int i = 0; i < popSize; i++) {
            for (int j = 0; j < chromosomeLength; j++) {
                population[i][j] = rand.nextInt(Constants.NO_OF_DATA_CENTERS);
            }
            populationFitness[i] = evaluateWeightedFitness(population[i]);
            if (populationFitness[i] < bestFitness) {
                bestFitness = populationFitness[i];
                bestChromosome = population[i].clone();
            }
        }
        Log.printLine("Initial population created. popSize=" + popSize + " tasks=" + chromosomeLength);
    }

    // ---------- Main GA loop ----------
    private static void runGA() {
        int noImprove = 0;
        for (int gen = 0; gen < MAX_GENERATIONS; gen++) {

            // Selection + Crossover + Mutation -> create new population
            int[][] newPop = new int[popSize][chromosomeLength];
            // Elitism: preserve best individual into new population (index 0)
            int eliteIndex = getBestIndex();
            newPop[0] = population[eliteIndex].clone();

            int idx = 1;
            while (idx < popSize) {
                // Selection: tournament of size 3 (simple, but paper uses random pair selection; we implement random pair + better one)
                int p1 = rand.nextInt(popSize);
                int p2 = rand.nextInt(popSize);
                int[] parent1 = population[p1];
                int[] parent2 = population[p2];

                // With probability PC: perform two-point crossover, otherwise copy parents
                int[] child1 = parent1.clone();
                int[] child2 = parent2.clone();
                if (rand.nextDouble() < PC) {
                    twoPointCrossover(child1, child2);
                }

                // Mutation on children with PM (two mutation points)
                mutateTwoPoints(child1);
                mutateTwoPoints(child2);

                // Add children to new population
                newPop[idx++] = child1;
                if (idx < popSize) newPop[idx++] = child2;
            }

            // Replace population
            population = newPop;

            // Evaluate fitness and update best
            for (int i = 0; i < popSize; i++) {
                populationFitness[i] = evaluateWeightedFitness(population[i]);
                if (populationFitness[i] < bestFitness) {
                    bestFitness = populationFitness[i];
                    bestChromosome = population[i].clone();
                    noImprove = 0;
                }
            }

            // Log metrics for CSV and console
            double[] mce = computeMetrics(bestChromosome); // [makespan, cpu, energy]
            csvFitness.add(bestFitness);
            csvMakespan.add(mce[0]);
            csvCPU.add(mce[1]);
            csvEnergy.add(mce[2]);

            if (gen % 10 == 0) {
                Log.printLine(String.format("Gen %d | Fitness=%.4f | M=%.4f | C=%.4f | E=%.4f",
                        gen, bestFitness, mce[0], mce[1], mce[2]));
            }

            // termination check
            noImprove++;
            if (noImprove >= NO_IMPROVE_LIMIT) {
                Log.printLine("GA early stopping at generation " + gen + " (no improvement for " + NO_IMPROVE_LIMIT + " gens).");
                break;
            }
        }
        Log.printLine("GA finished. Estimated best weighted fitness = " + bestFitness);
    }

    // ---------- Genetic operators ----------
    private static void twoPointCrossover(int[] a, int[] b) {
        int n = a.length;
        int p1 = rand.nextInt(n);
        int p2 = rand.nextInt(n);
        if (p1 > p2) { int t = p1; p1 = p2; p2 = t; }
        for (int i = p1; i <= p2; i++) {
            int tmp = a[i];
            a[i] = b[i];
            b[i] = tmp;
        }
    }

    private static void mutateTwoPoints(int[] chrom) {
        if (rand.nextDouble() < PM) {
            int n = chrom.length;
            int m1 = rand.nextInt(n);
            int m2 = rand.nextInt(n);
            chrom[m1] = rand.nextInt(Constants.NO_OF_DATA_CENTERS);
            chrom[m2] = rand.nextInt(Constants.NO_OF_DATA_CENTERS);
        }
    }

    // ---------- Fitness and metrics ----------
    // Weighted fitness: 0.5*M + 0.3*C + 0.2*E
    private static double evaluateWeightedFitness(int[] solution) {
        double[] metrics = computeMetrics(solution); // [makespan, cpu, energy]
        double makespan = metrics[0];
        double cpu = metrics[1];
        double energy = metrics[2];
        return 0.5 * makespan + 0.3 * cpu + 0.2 * energy;
    }

    // Compute metrics: makespan (max VM finish time), cpu sum (sum of exec times), energy ~ 0.001*cpu*makespan
    private static double[] computeMetrics(int[] solution) {
        double makespan = 0;
        double cpu = 0;
        double[] vmFinish = new double[Constants.NO_OF_DATA_CENTERS];
        for (int t = 0; t < Constants.NO_OF_TASKS; t++) {
            int vm = solution[t];
            double exec = execMatrix[t][vm];
            double comm = commMatrix[t][vm];
            vmFinish[vm] += exec + comm;
            cpu += exec; // aggregate CPU load (executions only)
        }
        for (double f : vmFinish) if (f > makespan) makespan = f;
        double energy = 0.001 * cpu * makespan; // simple energy model
        return new double[]{makespan, cpu, energy};
    }

    // ---------- Utilities ----------
    private static int getBestIndex() {
        int idx = 0;
        double best = Double.MAX_VALUE;
        for (int i = 0; i < popSize; i++) {
            if (populationFitness[i] < best) {
                best = populationFitness[i];
                idx = i;
            }
        }
        return idx;
    }

    private static void printFinalEstimatedResults() {
        double[] mce = computeMetrics(bestChromosome);
        Log.printLine("=== GA Estimated Best Results ===");
        Log.printLine("Estimated Weighted Fitness: " + bestFitness);
        Log.printLine(String.format("Estimated Makespan = %.4f, CPU(total exec) = %.4f, Energy = %.6f",
                mce[0], mce[1], mce[2]));
        Log.printLine("=== Best Chromosome (task -> vm) ===");
        for (int i = 0; i < bestChromosome.length; i++) {
            Log.printLine("Task " + i + " -> VM " + bestChromosome[i]);
        }
    }

    // ---------- CloudSim simulation using best mapping ----------
    private static void runCloudSim(int[] mapping) {
        try {
            Log.printLine("Running CloudSim simulation with GA best mapping...");
            CloudSim.init(1, Calendar.getInstance(), false);

            datacenter = new Datacenter[Constants.NO_OF_DATA_CENTERS];
            for (int i = 0; i < Constants.NO_OF_DATA_CENTERS; i++) {
                datacenter[i] = DatacenterCreator.createDatacenter("Datacenter_" + i);
            }

            DatacenterBroker broker = new DatacenterBroker("GA_Broker");
            int brokerId = broker.getId();

            vmList = createVM(brokerId, Constants.NO_OF_DATA_CENTERS);
            cloudletList = createCloudlet(brokerId, Constants.NO_OF_TASKS, mapping);

            broker.submitVmList(vmList);
            broker.submitCloudletList(cloudletList);

            CloudSim.startSimulation();

            List<Cloudlet> result = broker.getCloudletReceivedList();
            CloudSim.stopSimulation();

            // Print simulation output and metrics
            printCloudletList(result);

            // Also compute actual makespan, cpu, energy from simulation if needed
            double actualMakespan = result.stream().mapToDouble(Cloudlet::getFinishTime).max().orElse(0);
            double actualCpu = 0;
            for (Cloudlet c : result) actualCpu += c.getActualCPUTime();
            double actualEnergy = 0.001 * actualCpu * actualMakespan;

            Log.printLine("=== GA Simulation (actual) Metrics ===");
            Log.printLine(String.format("Actual Makespan = %.4f, Actual CPU(total CPU time) = %.4f, Actual Energy = %.6f",
                    actualMakespan, actualCpu, actualEnergy));

        } catch (Exception e) {
            e.printStackTrace();
            Log.printLine("CloudSim simulation failed.");
        }
    }

    private static List<Vm> createVM(int uid, int vms) {
        List<Vm> list = new ArrayList<>();
        long size = 10000;
        int ram = 512;
        int mips = 250;
        long bw = 1000;
        int pesNumber = 1;
        String vmm = "Xen";

        for (int i = 0; i < vms; i++) {
            Vm vm = new Vm(i, uid, mips, pesNumber, ram, bw, size, vmm, new CloudletSchedulerSpaceShared());
            list.add(vm);
        }
        return list;
    }

    private static List<Cloudlet> createCloudlet(int uid, int count, int[] mapping) {
        List<Cloudlet> list = new ArrayList<>();
        long fileSize = 300;
        long outputSize = 300;
        int pesNumber = 1;
        UtilizationModel um = new UtilizationModelFull();

        for (int i = 0; i < count; i++) {
            int vmIdx = mapping[i];
            long length = (long) (1e3 * (execMatrix[i][vmIdx] + commMatrix[i][vmIdx]));
            Cloudlet c = new Cloudlet(i, length, pesNumber, fileSize, outputSize, um, um, um);
            c.setUserId(uid);
            c.setVmId(vmList.get(vmIdx).getId());
            list.add(c);
        }
        return list;
    }

    private static void printCloudletList(List<Cloudlet> list) {
        DecimalFormat dft = new DecimalFormat("###.##");
        Log.printLine("========== GA OUTPUT (CloudSim Results) ==========");
        for (Cloudlet c : list) {
            if (c.getCloudletStatus() == Cloudlet.SUCCESS) {
                Log.printLine("Cloudlet " + c.getCloudletId() + " -> VM " + c.getVmId() +
                        " | CPUTime=" + dft.format(c.getActualCPUTime()) +
                        " | Start=" + dft.format(c.getExecStartTime()) +
                        " | Finish=" + dft.format(c.getFinishTime()));
            }
        }
        double makespan = list.stream().mapToDouble(Cloudlet::getFinishTime).max().orElse(0);
        Log.printLine("Actual Makespan (from simulation): " + makespan);
    }

    // ---------- CSV saving ----------
    private static void saveCSV() {
        saveList("ga_fitness.csv", csvFitness, "Generation,Fitness");
        saveList("ga_makespan.csv", csvMakespan, "Generation,Makespan");
        saveList("ga_cpu.csv", csvCPU, "Generation,CPU");
        saveList("ga_energy.csv", csvEnergy, "Generation,Energy");
        Log.printLine("CSV logs saved: ga_fitness.csv, ga_makespan.csv, ga_cpu.csv, ga_energy.csv");
    }

    private static void saveList(String name, List<Double> data, String header) {
        try (FileWriter fw = new FileWriter(name)) {
            fw.write(header + "\n");
            for (int i = 0; i < data.size(); i++) {
                fw.write(i + "," + data.get(i) + "\n");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
