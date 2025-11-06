package pso_package;

import org.cloudbus.cloudsim.*;
import org.cloudbus.cloudsim.core.CloudSim;
import org.cloudbus.cloudsim.provisioners.BwProvisionerSimple;
import org.cloudbus.cloudsim.provisioners.PeProvisionerSimple;
import org.cloudbus.cloudsim.provisioners.RamProvisionerSimple;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

/**
 * Main class to run a Particle Swarm Optimization (PSO) scheduler for task assignment
 * in a CloudSim environment.
 *
 * Replaces the SJF_Scheduler logic.
 */
public class PSOScheduler {

    // --- Simulation-related Class Members ---
    private static List<Cloudlet> cloudletList;
    private static List<Vm> vmList;
    private static Datacenter[] datacenter;
    private static double[][] commMatrix; // Communication Time Matrix
    private static double[][] execMatrix; // Execution Time Matrix

    // --- PSO-related Class Members ---
    private static Particle[] swarm;            // The population of particles
    private static int[] gBestPosition;         // The global best solution
    private static double gBestFitness;         // The fitness of the global best solution
    private static final Random rand = new Random();

    // --- PSO Configuration Parameters ---
    private static final int MAX_ITERATIONS = 100;    // Number of iterations to run PSO
    private static final double W = 0.729;          // Inertia weight
    private static final double C1 = 1.49445;       // Cognitive (personal) coefficient
    private static final double C2 = 1.49445;       // Social (global) coefficient

    /**
     * Main entry point.
     * 1. Runs the PSO algorithm to find an optimal schedule.
     * 2. Runs the CloudSim simulation using that schedule.
     */
    public static void main(String[] args) {
        Log.printLine("Starting PSO Scheduler...");

        // 1. Initialize Communication and Execution Matrices
        new GenerateMatrices();
        execMatrix = GenerateMatrices.getExecMatrix();
        commMatrix = GenerateMatrices.getCommMatrix();

        // 2. Run the PSO algorithm to find the best solution
        int[] psoSolution = runPSO();

        // 3. Run the CloudSim simulation with the solution found by PSO
        runCloudSim(psoSolution);
    }

    /**
     * Executes the Particle Swarm Optimization algorithm.
     * @return The best solution (task-to-VM mapping) found.
     */
    private static int[] runPSO() {
        Log.printLine("Running PSO Algorithm...");
        initializeSwarm();
        gBestFitness = Double.MAX_VALUE;
        gBestPosition = new int[Constants.NO_OF_TASKS];

        for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
            for (Particle p : swarm) {
                // 1. Calculate the fitness (makespan) for the particle's current position
                double currentFitness = calculateMakespan(p.position);

                // 2. Update pBest
                if (currentFitness < p.pBestFitness) {
                    p.pBestFitness = currentFitness;
                    System.arraycopy(p.position, 0, p.pBestPosition, 0, Constants.NO_OF_TASKS);
                }

                // 3. Update gBest
                if (currentFitness < gBestFitness) {
                    gBestFitness = currentFitness;
                    System.arraycopy(p.position, 0, gBestPosition, 0, Constants.NO_OF_TASKS);
                }
            }

            // 4. Update velocity and position for all particles
            for (Particle p : swarm) {
                updateVelocity(p);
                updatePosition(p);
            }
            
            Log.printLine("Iteration " + iter + ": Best Makespan = " + gBestFitness);
        }

        Log.printLine("PSO finished. Final best makespan (estimated): " + gBestFitness);
        return gBestPosition;
    }

    /**
     * Initializes the swarm of particles.
     */
    private static void initializeSwarm() {
        swarm = new Particle[Constants.POPULATION_SIZE];
        for (int i = 0; i < Constants.POPULATION_SIZE; i++) {
            swarm[i] = new Particle();
        }
    }

    /**
     * This is the PSO fitness function.
     * It estimates the makespan for a given solution (task-to-VM mapping)
     * based on the exec/comm time matrices.
     *
     * @param solution An array where solution[i] = VM_ID for Task_i.
     * @return The estimated makespan (total finish time).
     */
    private static double calculateMakespan(int[] solution) {
        double makespan = 0.0;
        // vmFinishTime[j] = the time when VM 'j' will be free.
        double[] vmFinishTime = new double[Constants.NO_OF_DATA_CENTERS];

        for (int i = 0; i < Constants.NO_OF_TASKS; i++) {
            int vmId = solution[i]; // The VM assigned to this task
            
            // Get the cost of this task (i) on this VM (vmId)
            double execTime = execMatrix[i][vmId];
            double commTime = commMatrix[i][vmId];
            double cost = execTime + commTime;
            
            // The task 'i' can only start after the VM 'vmId' is free.
            // So, its finish time is the VM's current finish time + its own cost.
            vmFinishTime[vmId] += cost;
        }

        // The makespan is the maximum finish time among all VMs.
        for (double finishTime : vmFinishTime) {
            if (finishTime > makespan) {
                makespan = finishTime;
            }
        }
        return makespan;
    }

    /**
     * Updates a particle's velocity based on its pBest and the gBest.
     */
    private static void updateVelocity(Particle p) {
        for (int i = 0; i < Constants.NO_OF_TASKS; i++) {
            double r1 = rand.nextDouble();
            double r2 = rand.nextDouble();

            // v = w*v + c1*r1*(pBest - x) + c2*r2*(gBest - x)
            double cognitive = C1 * r1 * (p.pBestPosition[i] - p.position[i]);
            double social = C2 * r2 * (gBestPosition[i] - p.position[i]);
            
            p.velocity[i] = (W * p.velocity[i]) + cognitive + social;
        }
    }

    /**
     * Updates a particle's position (solution) based on its new velocity.
     * This is a "discrete" PSO, so we map the new continuous position
     * back to a valid VM ID (0 to NO_OF_DATA_CENTERS - 1).
     */
    private static void updatePosition(Particle p) {
        for (int i = 0; i < Constants.NO_OF_TASKS; i++) {
            // x = x + v
            double newPosition = p.position[i] + p.velocity[i];

            // Map the continuous position back to a discrete VM ID
            // We use a simple modulo operator.
            // Math.round() helps discretize, Math.abs() handles negative velocities.
            int newVmId = (int) (Math.abs(Math.round(newPosition)) % Constants.NO_OF_DATA_CENTERS);
            
            p.position[i] = newVmId;
        }
    }


    /**
     * Initializes and runs the CloudSim simulation using the provided PSO solution.
     * This method is adapted from your SJF_Scheduler.main() method.
     *
     * @param psoSolution The task-to-VM mapping from runPSO().
     */
    private static void runCloudSim(int[] psoSolution) {
        Log.printLine("Initializing CloudSim...");
        try {
            int num_user = 1;
            Calendar calendar = Calendar.getInstance();
            boolean trace_flag = false;

            CloudSim.init(num_user, calendar, trace_flag);

            // Second step: Create Datacenters
            datacenter = new Datacenter[Constants.NO_OF_DATA_CENTERS];
            for (int i = 0; i < Constants.NO_OF_DATA_CENTERS; i++) {
                datacenter[i] = DatacenterCreator.createDatacenter("Datacenter_" + i);
            }

            // Third step: Create Broker
            // We use the *standard* DatacenterBroker, not the SJF one.
            DatacenterBroker broker = createBroker("Broker_0");
            int brokerId = broker.getId();

            // Fourth step: Create VMs and Cloudlets
            vmList = createVM(brokerId, Constants.NO_OF_DATA_CENTERS);
            
            // Create cloudlets using the PSO solution for binding
            cloudletList = createCloudlet(brokerId, Constants.NO_OF_TASKS, 0, psoSolution);

            broker.submitVmList(vmList);
            broker.submitCloudletList(cloudletList);

            // Fifth step: Starts the simulation
            CloudSim.startSimulation();

            // Final step: Print results
            List<Cloudlet> newList = broker.getCloudletReceivedList();
            CloudSim.stopSimulation();

            printCloudletList(newList);

            Log.printLine(PSOScheduler.class.getName() + " finished!");
        } catch (Exception e) {
            e.printStackTrace();
            Log.printLine("The simulation has been terminated due to an unexpected error");
        }
    }

    /**
     * Creates a standard DatacenterBroker.
     */
    private static DatacenterBroker createBroker(String name) throws Exception {
        return new DatacenterBroker(name);
    }

    /**
     * Creates the list of Virtual Machines.
     * Creates one VM per Datacenter.
     * VM IDs will be 0, 1, 2, 3, 4.
     */
    private static List<Vm> createVM(int userId, int vms) {
        LinkedList<Vm> list = new LinkedList<Vm>();

        long size = 10000;
        int ram = 512;
        int mips = 250;
        long bw = 1000;
        int pesNumber = 1;
        String vmm = "Xen";

        Vm[] vm = new Vm[vms];
        for (int i = 0; i < vms; i++) {
            // Give VMs simple IDs (0, 1, 2...)
            vm[i] = new Vm(i, userId, mips, pesNumber, ram, bw, size, vmm, new CloudletSchedulerSpaceShared());
            list.add(vm[i]);
        }
        return list;
    }

    /**
     * Creates the list of Cloudlets (Tasks).
     * Binds each cloudlet to a VM according to the psoSolution.
     *
     * @param psoSolution The mapping array from PSO (solution[i] = vm_index).
     */
    private static List<Cloudlet> createCloudlet(int userId, int cloudlets, int idShift, int[] psoSolution) {
        LinkedList<Cloudlet> list = new LinkedList<Cloudlet>();

        long fileSize = 300;
        long outputSize = 300;
        int pesNumber = 1;
        UtilizationModel utilizationModel = new UtilizationModelFull();

        Cloudlet[] cloudlet = new Cloudlet[cloudlets];

        for (int i = 0; i < cloudlets; i++) {
            // Get the VM *index* (0-4) from the PSO solution
            int vmIndex = psoSolution[i];
            
            // Get the *actual* VM ID from the vmList (which will also be 0-4)
            int vmId = vmList.get(vmIndex).getId();

            // Calculate task length based on exec and comm time for that (task, vm) pair
            long length = (long) (1e3 * (commMatrix[i][vmIndex] + execMatrix[i][vmIndex]));
            
            cloudlet[i] = new Cloudlet(idShift + i, length, pesNumber, fileSize, outputSize, utilizationModel, utilizationModel, utilizationModel);
            cloudlet[i].setUserId(userId);
            
            // *** This is the key step: Assign the Cloudlet to the VM chosen by PSO ***
            cloudlet[i].setVmId(vmId);
            
            list.add(cloudlet[i]);
        }
        return list;
    }

    /**
     * Prints the simulation output.
     */
    private static void printCloudletList(List<Cloudlet> list) {
        int size = list.size();
        Cloudlet cloudlet;

        String indent = "    ";
        Log.printLine();
        Log.printLine("========== OUTPUT ==========");
        Log.printLine("Cloudlet ID" + indent + "STATUS" +
                indent + "Data center ID" +
                indent + "VM ID" +
                indent + indent + "Time" +
                indent + "Start Time" +
                indent + "Finish Time" +
                indent + "Waiting Time");

        DecimalFormat dft = new DecimalFormat("###.##");
        dft.setMinimumIntegerDigits(2);
        for (int i = 0; i < size; i++) {
            cloudlet = list.get(i);
            Log.print(indent + dft.format(cloudlet.getCloudletId()) + indent + indent);

            if (cloudlet.getCloudletStatus() == Cloudlet.SUCCESS) {
                Log.print("SUCCESS");

                Log.printLine(indent + indent + dft.format(cloudlet.getResourceId()) +
                        indent + indent + indent + dft.format(cloudlet.getVmId()) +
                        indent + indent + dft.format(cloudlet.getActualCPUTime()) +
                        indent + indent + dft.format(cloudlet.getExecStartTime()) +
                        indent + indent + indent + dft.format(cloudlet.getFinishTime()) +
                        indent + indent + indent + dft.format(cloudlet.getWaitingTime()));
            }
        }
        
        // Calculate the *actual* makespan from the simulation results
        double makespan = calcActualMakespan(list);
        Log.printLine("Actual Makespan (from simulation): " + makespan);
    }

    /**
     * Calculates the *actual* makespan from the simulation results
     * by finding the maximum finish time of any cloudlet.
     */
    private static double calcActualMakespan(List<Cloudlet> list) {
        double makespan = 0;
        for (Cloudlet cloudlet : list) {
            if (cloudlet.getFinishTime() > makespan) {
                makespan = cloudlet.getFinishTime();
            }
        }
        return makespan;
    }
}