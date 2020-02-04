package org.opensourcephysics.stp.percolation;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.frames.*;

public class fgrApp extends AbstractSimulation {
	/* This subroutine implements the algorithm for fission gas release fraction 
	 * The Forsberg-Massih algorithm is adopted as descripbed in Ref.1
	 * Ref.1  P.C. Millett et al. /Journal of Nuclear Materials 424 (2012) 176–182
	 */
	  Scalar2DFrame grid = new Scalar2DFrame("FGR percolation algorithm");
	  FGClusters lattice;
	  double tDisplay,t1,t2,dt;
	  double t=0.0; // current time
	  public void initialize() {
		    int L = control.getInt("Lattice X-axis size L");
		    int M = control.getInt("Lattice Y-axis size M");
		    grid.resizeGrid(L, M);
		    lattice = new FGClusters(L, M);
		    tDisplay = control.getDouble("Display lattice at this time");
		    t1 = control.getDouble("Beginning time for plots");
		    t2 = control.getDouble("Ending time for plots");
		    dt = control.getDouble("dt for plots");
		    // adds sites to new cluster
		    lattice.newLattice();
		    displayLattice();//paint the lattice
	  }
	  
	  public void doStep() {
		    control.clearMessages();
		    t+=dt;
		    control.println("Trial "+t);
		    grid.setMessage("t = "+tDisplay);
		    displayLattice();//paint the lattice
	  }
	  
	  private void displayLattice() {
		    double display[] = new double[lattice.N];//color value of each site
		    for(int s = 0;s<lattice.N;s++) {
		      //display[s] = lattice.getClusterSize(s);
		      display[s] = s;// for test
		    }
		    grid.setAll(display);
		  }
	  
	  public void reset() {
		  control.setValue("Lattice X-axis size L", 128);
		   control.setValue("Lattice Y-axis size M", 128);
		   control.setValue("Beginning time for plots", 0.0);
		   control.setValue("Ending time for plots", 100000);
		   control.setValue("dt for plots", 1000);
		   control.setValue("Display lattice at this time", 50000);
	  }
	  
	public static void main(String args[]) {
		SimulationControl.createApp(new fgrApp());
		//fgrTest grain = new fgrTest();
	    //grain.fgrCalculate();
	  }
}
