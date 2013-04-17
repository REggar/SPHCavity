import java.util.ArrayList;
import java.util.List;

public class Grid {
	
	// This class contains information of the grid

	public int width; // Width of the grid in cells
	public int height; // Height of the grid in cells
	
	Cell[][] grid = new Cell[width][height]; // Initialise the grid
	
	public void clean(){
		// Clear all contents of grid.
		for (int i = 0; i < grid.length-1; i++) {
			for (int j = 0; j < grid[i].length-1; j++) {
				grid[i][j].clean();
			}
		}
	}
	
	public void addToCell(Particle particle){
		// Add particle to which cell of grid
		grid[particle.gridx][particle.gridy].add(particle);
	}
	
	public List<Particle> particles(Particle particle){
		// List all particles in range of particle of interest for linked list algorithm
		ArrayList<Particle> particles = new ArrayList<Particle>();
		for (int i = -1; i <= 1; i++) {
			if(particle.gridx+i<0 || particle.gridx+i>width) continue;
			for (int j = -1; j <= 1; j++) {
				if(particle.gridy+j<0 || particle.gridy+j>height) continue;
				particles.addAll(grid[particle.gridx+i][particle.gridy+j].contents());
			}
		}
		return particles;
	}
	
	public Grid(List<Particle> particleList){
		// Initialise grid
		this.width = (int) Math.ceil(((GenerateParticles.BOX_WIDTH/GenerateParticles.PARTICLE_SPACING)+6)/(2*SPHCavityProblem.H));
		this.height = (int) Math.ceil(((GenerateParticles.BOX_HEIGHT/GenerateParticles.PARTICLE_SPACING)+6)/(2*SPHCavityProblem.H));
		
		this.grid = new Cell[width+1][height+1];
		
		for (int i = 0; i < width+1; i++) {
			for (int j = 0; j < height+1; j++) {
				this.grid[i][j] = new Cell();
			}
		}
	}
}