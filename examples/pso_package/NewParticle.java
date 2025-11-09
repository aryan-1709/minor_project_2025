package pso_package;

import java.util.Random;

/**
 * Particle class for Local Best Murmuration PSO (PSOMbest)
 * Implements the structure described in:
 * "A novel improvement of particle swarm optimization using an improved velocity update
 * function based on local best murmuration particle" — Twumasi et al. (2024)
 */
public class NewParticle {

    private int[] position;            // Current position (solution)
    private double[] velocity;         // Velocity vector
    private int[] personalBest;        // Personal best position
    private double fitness;            // Current fitness value
    private double personalBestFitness;// Best fitness found by this particle

    private static final Random rand = new Random();

    // --------------------------------------------------------------------
    // Constructor: Initialize particle with random position and velocity
    // --------------------------------------------------------------------
    public NewParticle(int numTasks, int numVMs) {
        position = new int[numTasks];
        velocity = new double[numTasks];
        personalBest = new int[numTasks];

        // Random initialization
        for (int i = 0; i < numTasks; i++) {
            position[i] = rand.nextInt(numVMs);        // assign task to a random VM
            velocity[i] = rand.nextDouble() * 2 - 1;   // random velocity in [-1, 1]
            personalBest[i] = position[i];
        }

        fitness = Double.MAX_VALUE;
        personalBestFitness = Double.MAX_VALUE;
    }
    
    public NewParticle(int[] givenPosition, int numVMs) {
        this.position = givenPosition.clone();
        this.velocity = new double[givenPosition.length];
        this.personalBest = givenPosition.clone();
        this.fitness = Double.MAX_VALUE;
        this.personalBestFitness = Double.MAX_VALUE;
    }

    // --------------------------------------------------------------------
    // Evaluate fitness using makespan calculation (as defined in PSOScheduler)
    // --------------------------------------------------------------------
    public void evaluateFitness(double[][] commMatrix, double[][] execMatrix) {
        this.fitness = calculateMakespan(position, commMatrix, execMatrix);
        if (this.fitness < this.personalBestFitness) {
            this.personalBestFitness = this.fitness;
            this.personalBest = position.clone();
        }
    }

    // --------------------------------------------------------------------
    // Makespan calculation (same as in PSOScheduler)
    // --------------------------------------------------------------------
    private static double calculateMakespan(int[] solution, double[][] commMatrix, double[][] execMatrix) {
        double makespan = 0.0;
        double[] vmFinishTime = new double[Constants.NO_OF_DATA_CENTERS];

        for (int i = 0; i < Constants.NO_OF_TASKS; i++) {
            int vmId = solution[i];
            double execTime = execMatrix[i][vmId];
            double commTime = commMatrix[i][vmId];
            double cost = execTime + commTime;
            vmFinishTime[vmId] += cost;
        }

        for (double finishTime : vmFinishTime) {
            if (finishTime > makespan) {
                makespan = finishTime;
            }
        }
        return makespan;
    }

    // --------------------------------------------------------------------
    // Getters and Setters
    // --------------------------------------------------------------------
    public int[] getPosition() {
        return position;
    }

    public void setPosition(int[] position) {
        this.position = position;
    }

    public double[] getVelocity() {
        return velocity;
    }

    public void setVelocity(double[] velocity) {
        this.velocity = velocity;
    }

    public int[] getPersonalBest() {
        return personalBest;
    }

    public void setPersonalBest(int[] personalBest) {
        this.personalBest = personalBest;
    }

    public double getFitness() {
        return fitness;
    }

    public double getPersonalBestFitness() {
        return personalBestFitness;
    }

    public void setPersonalBestFitness(double value) {
        this.personalBestFitness = value;
    }

    // --------------------------------------------------------------------
    // Utility: Clone particle
    // --------------------------------------------------------------------
    public NewParticle clone() {
        NewParticle copy = new NewParticle(position.length, Constants.NO_OF_DATA_CENTERS);
        copy.position = position.clone();
        copy.velocity = velocity.clone();
        copy.personalBest = personalBest.clone();
        copy.fitness = fitness;
        copy.personalBestFitness = personalBestFitness;
        return copy;
    }
}
