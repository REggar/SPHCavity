import java.util.ArrayList;
import java.util.List;

public class GenerateParticles {

	// This class generates the initial box of particles

	public static double PARTICLE_SPACING = 1; // Distance between particles
	public static int BOX_WIDTH = 50; // Width of box
	public static int BOX_HEIGHT = 50; // Height of box
	public static int DUMMY_WIDTH = 2; // How many dummy particles at edges

	public List<Particle> particleList; 
	
	// Function to generate particles
	
	public List<Particle> generate() {
		particleList = new ArrayList<Particle>();
		int n=0;
		
		for(double i = 0; i < BOX_WIDTH+2*DUMMY_WIDTH; i++){
			for(double j = 0; j < BOX_HEIGHT+2*DUMMY_WIDTH; j++){
				// Adding all the particles to the list of particles
				if(i>=(BOX_WIDTH+DUMMY_WIDTH) || i< DUMMY_WIDTH) particleList.add(new Particle(1f, i*PARTICLE_SPACING , j*PARTICLE_SPACING , 0f, 0f, true, false, false, 0));
				else if(j >=(BOX_HEIGHT+DUMMY_WIDTH)) particleList.add(new Particle(1f, i*PARTICLE_SPACING , j*PARTICLE_SPACING , 0f, 0f, true, false, false, 0));
				else if(j < DUMMY_WIDTH){
					n++;
					particleList.add(new Particle(1f, i*PARTICLE_SPACING , j*PARTICLE_SPACING , SPHCavityProblem.PLATE_VELOCITY, 0f, true, false, true, n));
				}
				else {
					n++;
					particleList.add(new Particle(1f, i*PARTICLE_SPACING, j*PARTICLE_SPACING, 0f, 0f, false, false, false, n));		
				}
			}
		}
		return particleList;
	}

	public GenerateParticles() {
	}
}