package org.opensourcephysics.stp.percolation;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.display.Dataset;
import org.opensourcephysics.frames.*;

public class fgrApp extends AbstractSimulation {
	/* This subroutine implements the algorithm for fission gas release fraction 
	 * The Forsberg-Massih algorithm is adopted as descripbed in Ref.1
	 * Ref.1  P.C. Millett et al. /Journal of Nuclear Materials 424 (2012) 176–182
	 */
	  Scalar2DFrame grid = new Scalar2DFrame("FGR percolation algorithm");
	  Scalar2DFrame grid2 = new Scalar2DFrame("GB status");
	  PlotFrame plot1 = new PlotFrame("Radial node No.", "FGR Fraction", "Fission Gas Release Fraction");
	  FGClusters lattice;
	  double tDisplay,t1,t2,dt;
	  double t; // current time
	  int L,M;
	  public void initialize() {
		    L = control.getInt("Lattice X-axis size L");
		    M = control.getInt("Lattice Y-axis size M");
		    tDisplay = 86400*control.getDouble("Display lattice at this time, days");// convert time from days to seconds
		    t1 = 86400*control.getDouble("Beginning time for plots, days");// convert time from days to seconds
		    t2 = 86400*control.getDouble("Ending time for plots, days");// convert time from days to seconds
		    dt = 86400*control.getDouble("dt for plots, days");// convert time from days to seconds
		    t=0.0;
		    grid.resizeGrid(L, M);
		    grid2.resizeGrid(L, M);
		    lattice = new FGClusters(L, M, dt);
		    // adds sites to new cluster
		    lattice.newLattice();
		    lattice.dt= dt;
		    displayLattice();//paint the lattice
	  }
	  
	  public void doStep() {
		    control.clearMessages();
		    t+=dt;
		    double tDays = t/86400;
		    control.println("Current Time(days)="+tDays);
		    grid.setMessage("t = "+tDays);
		    grid2.setMessage("t = "+tDays);
	        lattice.updateSite();
		    //for(int i = 0;i<lattice.N;i++) {

		        //meanClusterSize[i] += lattice.getMeanClusterSize();
		        //P_infinity[i] += (double) lattice.getSpanningClusterSize()/lattice.numSitesOccupied;
		        //P_span[i] += (lattice.getSpanningClusterSize()==0 ? 0 : 1);
		    //  }
		    displayLattice();//paint the lattice
		    lattice.updateVenting();//
		    displayLattice();//paint the lattice
		    plotAverages();
	  }
	  private void plotAverages() {
		  plot1.clearData();
		  plot1.setConnected(1, true);
		  plot1.setMarkerShape(1, Dataset.NO_MARKER);
		  //plot1.showDataTable(true);
		  for (int x=0;x<L;x++) {
			  plot1.append(0, x, lattice.getFGRP(x));
			  plot1.append(1, x, lattice.getFGRD(x));
		  }
	  }
	  private void displayLattice() {
		    double display[] = new double[lattice.N];//color value of each site
		    for(int s = 0;s<lattice.N;s++) {
		      display[s] = lattice.getClusterSize(s);
		      //display[s] = s;// for test
		    }
		    grid.setAll(display);
		    for(int s = 0;s<lattice.N;s++) {
			      display[s] = lattice.gbState[s];
			      //display[s] = s;// for test
			    }
			grid2.setAll(display);
		  }
	  
	  public void reset() {
		   control.setValue("Lattice X-axis size L", 10);
		   control.setValue("Lattice Y-axis size M", 10);
		   control.setValue("Beginning time for plots, days", 0.0);
		   control.setValue("Ending time for plots, days", 1000);
		   control.setValue("dt for plots, days", 10);//1e6 s = 11.57 days
		   control.setValue("Display lattice at this time, days", 500);
	  }
	  
	public static void main(String args[]) {
		SimulationControl.createApp(new fgrApp());
		//fgrTest grain = new fgrTest();
	    //grain.fgrCalculate();
	  }
}
