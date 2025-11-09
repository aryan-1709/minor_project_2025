package pso_package;

public class MainSA {

    public static void main(String[] args) {

        System.out.println("\n--- Running Simulated Annealing Scheduler ---\n");

        // Step 1: Generate matrices
        new GenerateMatrices();
        double[][] execMatrix = GenerateMatrices.getExecMatrix();
        double[][] commMatrix = GenerateMatrices.getCommMatrix();

        // ✅ Step 2: COPY MATRICES INTO PSOScheduler (VERY IMPORTANT)
        PSOScheduler.execMatrix = execMatrix;
        PSOScheduler.commMatrix = commMatrix;

        // Step 3: Run SA
        int[] saSolution = SimulatedAnnealingScheduler.runSA(execMatrix, commMatrix);

        // Step 4: Run CloudSim using SA solution
        PSOScheduler.runCloudSim(saSolution);

        System.out.println("\n--- Simulation Finished (SA Based Scheduling) ---\n");
    }
}
