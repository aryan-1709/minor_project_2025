package pso_package;

import java.util.Random;

/**
 * Represents a single particle in the PSO swarm.
 * Each particle holds a potential solution (position), its velocity,
 * and its personal best (pBest) solution found so far.
 */
public class Particle {

    // --- Public fields for easy access from the main scheduler ---

    // The particle's current position (the solution).
    // The index is the Task ID (0 to 29).
    // The value is the VM ID (0 to 4) it's assigned to.
    public int[] position;

    // The particle's current velocity.
    public double[] velocity;

    // The particle's personal best position found so far.
    public int[] pBestPosition;

    // The fitness (makespan) of the pBestPosition.
    public double pBestFitness;

    private static final Random rand = new Random();

    /**
     * Constructor for a new Particle.
     * Initializes the particle with a random position (solution)
     * and a zero velocity.
     */
    public Particle() {
        this.position = new int[Constants.NO_OF_TASKS];
        this.velocity = new double[Constants.NO_OF_TASKS];
        this.pBestPosition = new int[Constants.NO_OF_TASKS];
        this.pBestFitness = Double.MAX_VALUE;

        // Initialize with a random solution
        for (int i = 0; i < Constants.NO_OF_TASKS; i++) {
            // Assign Task 'i' to a random VM (0 to NO_OF_DATA_CENTERS - 1)
            this.position[i] = rand.nextInt(Constants.NO_OF_DATA_CENTERS);
            this.velocity[i] = 0.0; // Start with zero velocity
        }

        // Copy the initial position to pBest
        System.arraycopy(this.position, 0, this.pBestPosition, 0, Constants.NO_OF_TASKS);
    }
}