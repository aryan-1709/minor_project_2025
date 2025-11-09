package pso_package;

public class MainCompare {

    public static void main(String[] args) {

        System.out.println("\n--- Running Comparative Study: SA vs ACO ---\n");

        // Step 1: Generate same input matrices once
        new GenerateMatrices();
        double[][] execMatrix = GenerateMatrices.getExecMatrix();
        double[][] commMatrix = GenerateMatrices.getCommMatrix();

        // Step 2: Copy to PSOScheduler
        PSOScheduler.execMatrix = execMatrix;
        PSOScheduler.commMatrix = commMatrix;

        // Step 3: Run Simulated Annealing
        long startSA = System.currentTimeMillis();
        int[] saSolution = SimulatedAnnealingScheduler.runSA(execMatrix, commMatrix);
        long endSA = System.currentTimeMillis();
        double saTime = (endSA - startSA) / 1000.0;

        System.out.println("\n>>> Running CloudSim with SA Solution...");
        PSOScheduler.runCloudSim(saSolution, "SA");

        // Step 4: Run Ant Colony Optimization
        long startACO = System.currentTimeMillis();
        int[] acoSolution = AntColonyScheduler.runACO(execMatrix, commMatrix);
        long endACO = System.currentTimeMillis();
        double acoTime = (endACO - startACO) / 1000.0;

        System.out.println("\n>>> Running CloudSim with ACO Solution...");
        PSOScheduler.runCloudSim(acoSolution, "ACO");

        // Step 5: Print comparative stats
        System.out.println("\n================= COMPARISON REPORT =================");
        System.out.printf("Algorithm\t|\tRuntime (s)\n");
        System.out.println("------------------------------------------------------");
        System.out.printf("Simulated Annealing\t|\t%.2f\n", saTime);
        System.out.printf("Ant Colony Optimization\t|\t%.2f\n", acoTime);
        System.out.println("======================================================\n");

        System.out.println("--- Comparative Study Finished ---");
    }
}
