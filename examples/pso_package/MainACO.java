package pso_package;

public class MainACO {

    public static void main(String[] args) {

        System.out.println("\n--- Running Ant Colony Optimization Scheduler ---\n");

        // Step 1: Generate matrices
        new GenerateMatrices();
        double[][] execMatrix = GenerateMatrices.getExecMatrix();
        double[][] commMatrix = GenerateMatrices.getCommMatrix();

        // Step 2: COPY MATRICES INTO PSOScheduler (if needed for CloudSim)
        PSOScheduler.execMatrix = execMatrix;
        PSOScheduler.commMatrix = commMatrix;

        // Step 3: Run ACO
        int[] acoSolution = AntColonyScheduler.runACO(execMatrix, commMatrix);

        // Step 4: Run CloudSim using ACO solution
        PSOScheduler.runCloudSim(acoSolution, "ACO");


        System.out.println("\n--- Simulation Finished (ACO Based Scheduling) ---\n");
    }
}
