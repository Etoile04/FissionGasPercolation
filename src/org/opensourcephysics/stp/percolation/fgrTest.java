package org.opensourcephysics.stp.percolation;
import org.opensourcephysics.frames.PlotFrame;
import org.opensourcephysics.display.Dataset;
//import javax.swing.JFrame.*;

public class fgrTest {
	/* This subroutine implements the algorithm for fission gas release fraction 
	 * The Forsberg-Massih algorithm is adopted as descripbed in Ref.1
	 * Ref.1  P.C. Millett et al. /Journal of Nuclear Materials 424 (2012) 176–182
	 */
	int noGB;      // number of the calculated grain boundary in the bond network
	int xGB;       // x position of the calculated grain boundary in the bond network
	int yGB;       // y position of the calculated grain boundary in the bond network
	int xGB_max;   // maximum number of the bond network in x axis
	// why is there no yGB_max?because T is dependent on X position but not Y position
	int GBstate = 0; // the calculated grain boundary state index- 0, closed;1, satured; 2,vented; 
	double T;      // local temperature of the calculated grain boundary
	//double beta_g; // local gas production rate, unit: atoms/m^3/s
	//double Dg;     // gas atom diffusion coefficient, unit: m^2/s
	//double Nt;     // the number of gas atoms per grain boundary area at time t, unit: atoms/m^2
	//double Nsat;   // the saturation density of grain boundary, unit: atoms/m^2
	double bgb = 1.0e-5;    // the resolution rate, unit: 1/s
	double theta = 50; // bubble’s contact angle, unit: degree
	double a = 1.0e-5; // the radius of grain, unit: m
	double rb = 5e-7;  // the radius of curvature of the grain boundary bubble, unit: m
	double Xgb_c = 0.5; //threshold fraction of grain boundary gas bubble coverage
	double Pext = 0.0;  // external applied hydro-static pressure, unit: Pa
	double gamma = 0.5; // the bubble surface energy, unit: J/m^2
	double delta = 1.0e-8; //the resolution depth, unit: m
	double delta_t = 1000; // time step, unit:second
	double Ni, Nt; //the previous and current values of Nt
	double C_fgr=0.0;// the amount of fission gas released after saturation,unit:mol
	final double Navgd=6.022e23;//Avogadro constant
	final double kB = 1.380649e-23; // Boltzmann constant, unit:J/K
//
	public double T(int xGB, int xGB_max) {
		double Tmin = 800; //unit: K
		double Tmax = 1200;//unit: K
		double r = xGB;
		double r0= xGB_max;
		double T = Tmax - (Tmax-Tmin)*(r/r0)*(r/r0); //Eq.(10) in Ref.1
		return T;
	}
//
	public double Dg(double T) {
		double Dg;
		if (T>1650)
			Dg = 1.09e-17*Math.exp(-6614/T);//Eq.(12) in Ref.1
		else if (T>1381) 
			Dg = 2.14e-13*Math.exp(-22884/T);//Eq.(12) in Ref.1
		else 
			Dg = 1.51e-17*Math.exp(-9508/T);//Eq.(12) in Ref.1
		return Dg;
	}
//	
	public double beta_g(double T) {
		double beta_g = -2.218e18+ 3.854e15*T;
		return beta_g;
	}
//
	public double f_theta (double theta){
		double theta_r = Math.toRadians(theta); // unit conversion from degree to radians
		double f_theta = 1 - 1.5* Math.cos(theta_r) + 0.5* Math.pow(Math.cos(theta_r), 3); 
		return f_theta;
	}
	public double Nsat (double T){
		double theta_r = Math.toRadians(theta); // unit conversion from degree to radians
		//double f_theta = 1 - 1.5* Math.cos(theta_r) + 0.5* Math.pow(Math.cos(theta_r), 3); 
		double Nsat = 4*rb*f_theta(theta)/3/kB/T/Math.sin(theta_r)/Math.sin(theta_r)*Xgb_c*(2*gamma/rb+Pext); // Eq.(9) in Ref.1,
		return Nsat;
	}
//
	public double f0(double T,double t) {
		double f0;
		double omega=Dg(T)*t/a/a;
		double pi_2=1/Math.PI/Math.PI;
		if (omega>pi_2)
			f0= beta_g(T)*a/3*(1-0.608*Math.exp(-Math.PI*Math.PI*Dg(T)*t/a/a));//Eq.(7) in Ref.1, long times;
		else
			f0 = beta_g(T)*a/3*(6/a*Math.sqrt(Dg(T)*t/Math.PI))-3*Dg(T)*t/a/a;//Eq.(6) in Ref.1, short times;
		return f0;
	}
//
	public double dNdt (double T,double t){
		//double Agb = 4.0* Math.PI *a*a;//the area of the grain boundary of a particular grain, unit: m^2
		//t is time,s
		//T is temperature, K
		//double f0;//The Booth flux
		//double t1=t+delta_t;
		double JFlux;
		if(t>0)
			JFlux= 2*f0(T, t)*(1-bgb*delta*Ni/2/Dg(T)/beta_g(T)/t);//Eq.(8) in Ref.1 solved with the Euler scheme;
		else
			JFlux= 2*f0(T, t+delta_t);
		return JFlux;
	}
//
	public double Nt (double T, double t){
		//t is time,s
		//T is temperature, K
		//Ni is the Nt in previous step
		//double f0;//The Booth flux
		double Np,Nc, Nt;// predicted, corrected values of Nt in improved Euler's scheme
		double t1=t+delta_t;
		if(t>0) {
			Np = 2*delta_t*f0(T, t) + (1-bgb*delta*f0(T, t)*delta_t/Dg(T)/beta_g(T)/t)*Ni;
			Nc = 2*delta_t*f0(T, t1) + Ni-bgb*delta*f0(T, t1)*delta_t/Dg(T)/beta_g(T)/t1*Np;
			Nt=0.5*Np+0.5*Nc;
				}
		else {
			Nt = 2*delta_t*f0(T, t);	
			}
		if(Nt>Nsat(T)) {
			//when Nt exceeds Nsat for the first time, GB state is set to 1, and GB is overpressurized. 
			//Nt=Nsat(T);//test 
			GBstate=1;
		}
		//if the grain boundary is vented in previous 
		if(GBstate==2) {
			C_fgr+=Nt*Math.PI*a*a/Navgd;//gas released amount is accumulated, unit:mol
			//Ni=0.0;// gas density is reset to zero after venting
			Nt=0.0;// gas density is reset to zero after venting
			GBstate=0;//grain boundary state is reset to 'closed'
		}
		return Nt;
	}
//	
	public void fgrCalculate(){
		double Ngas[]=new double[200];
		Ni=0;
		double T= this.T(xGB, xGB_max);
		double t=0;
		//double deltaN;
		int stepNo=120;
		Ngas[0]=0.0;
		int j;
		//plot setting
		PlotFrame frame= new PlotFrame("time(s)","Number density at Grain Boudnary","FGR test");
		frame.setConnected(1, true);
		frame.setMarkerShape(1, Dataset.NO_MARKER);
		//
		for(int i=0;i<=stepNo;i++) {
			j=i+1;
			//Ngas[j]= Ngas[i];
			//deltaN = dNdt(T,t)*delta_t;
			t=t+delta_t;
			Ngas[j]= Nt(T,t);
			Ni=Ngas[j];
			frame.append(0, t, Ni);
			frame.append(1, t, Nsat(T));
			//System.out.println(j);
			//System.out.println(Ni);
			//System.out.println(deltaN);
			//System.out.println(Nsat(T));
		}
		frame.setVisible(true);
		frame.setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);
		frame.setXPointsLinked(true);
		frame.setXYColumnNames(0, "time(s)", "N_boundary");
		frame.setXYColumnNames(1, "time(s)", "Nsat");
		frame.setRowNumberVisible(true);
	}
//
	public static void main(String args[]) {
		fgrTest grain = new fgrTest();
	    grain.xGB = 100;
	    grain.xGB_max = 800;
	    //int xGB = grain.xGB;
	    //int xGB_max=grain.xGB_max;
	    //double T = grain.T(xGB, xGB_max);
	    //double beta_g =beta_g(T);
	    //double Dg = Dg(T);
	    //double Nsat= Nsat(T);
	    //System.out.println(T);
	    //System.out.println(grain.beta_g(T));
	    //System.out.println(grain.Dg(T));
	    //System.out.println(grain.Nsat(T));
	    grain.delta_t = 500000;
	    grain.fgrCalculate();
	  }
}
