package org.opensourcephysics.stp.percolation;
//import org.opensourcephysics.frames.PlotFrame;
//import org.opensourcephysics.display.Dataset;
//import javax.swing.JFrame.*;

public class FGR {
	/* This subroutine implements the algorithm for fission gas release fraction 
	 * The Forsberg-Massih algorithm is adopted as descripbed in Ref.1
	 * Ref.1  P.C. Millett et al. /Journal of Nuclear Materials 424 (2012) 176–182
	 */
	int noGB;      // number of the calculated grain boundary in the bond network
	int xGB;       // x position of the calculated grain boundary in the bond network
	int yGB;       // y position of the calculated grain boundary in the bond network
	int xGB_max;   // maximum number of the bond network in x axis
	// why is there no yGB_max?because T is dependent on X position but not Y position
	int GBstate; // the calculated grain boundary state index- 0, closed;1, satured; 2,overpressurized;
	double t;       // time since irradiation
	double local_t; //time since last venting, used for Booth flux 
	double T;      // local temperature of the calculated grain boundary
	//double beta_g; // local gas production rate, unit: atoms/m^3/s
	//double Dg;     // gas atom diffusion coefficient, unit: m^2/s
	//double Nt;     // the number of gas atoms per grain boundary area at time t, unit: atoms/m^2
	//double Nsat;   // the saturation density of grain boundary, unit: atoms/m^2
	double bgb = 1.0e-5;    // the resolution rate, unit: 1/s
	double theta; // bubble’s contact angle, unit: degree
	double a = 1.0e-5; // the radius of grain, unit: m
	double rb = 5e-7;  // the radius of curvature of the grain boundary bubble, unit: m
	double Xgb_c = 0.5; //threshold fraction of grain boundary gas bubble coverage
	double Pext = 0.0;  // external applied hydro-static pressure, unit: Pa
	double gamma = 0.5; // the bubble surface energy, unit: J/m^2
	double delta = 1.0e-8; //the resolution depth, unit: m
	double dt; // time step, unit:second
	double Ni, Nt; //the previous and current values of Nt, unit:atoms/m^2
	double C_fgr=0.0;// the amount of fission gas released after saturation,unit:mol
	double C_fgp; // fission gas produced at time t, unit:mol
	double fgrPercolation=0.0;// fission gas release fraction according to percolation
	double fgrDiffusion=0.0;  // fission gas release fraction according to diffusion theory
	final double Navgd=6.022e23;//Avogadro constant
	final double kB = 1.380649e-23; // Boltzmann constant, unit:J/K
	final double PI = Math.PI;
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
	public double Dg(double T) {//Eq.(12) in Ref.1
		double Dg;
		if (T>1650)
			Dg = 1.09e-17*Math.exp(-6614/T);
		else if (T>1381) 
			Dg = 2.14e-13*Math.exp(-22884/T);
		else 
			Dg = 1.51e-17*Math.exp(-9508/T);
		return Dg;
	}
	public double DTurbull(double T) {//Eq.(12) in Ref.1
		double D1=7.6e-10*Math.exp(-35000/T);
		double D2=1.77e-15*Math.exp(-13800/T);
		double D3=2e-21;
		double Dg=1/(D1+D2+D3)+1e15;
		Dg=1.0/Dg;
		return Dg;
	}
//	
	public double beta_g(double T) {
		double beta_g = -2.218e18+ 3.854e15*T;
		return beta_g;
	}
//
	public double f_theta (double theta){
		//double theta_r = Math.toRadians(theta); // unit conversion from degree to radians
		double f_theta = 1 - 1.5* Math.cos(theta) + 0.5* Math.pow(Math.cos(theta), 3); 
		return f_theta;
	}
	public double Nsat (double T){
		//double theta_r = Math.toRadians(theta); // unit conversion from degree to radians
		//double f_theta = 1 - 1.5* Math.cos(theta_r) + 0.5* Math.pow(Math.cos(theta_r), 3); 
		double Nsat = 4*rb*f_theta(theta)/3/kB/T/Math.sin(theta)/Math.sin(theta)*Xgb_c*(2*gamma/rb+Pext); // Eq.(9) in Ref.1,
		return Nsat;
	}
//
	public double f0(double T,double t) {
		//t is time in previous step,s
		double f0;
		double beta_g = beta_g(T);
		double Dg=Dg(T);
		double omega=Dg*t/a/a;
		double pi_2=1/PI/PI;
		if (omega>pi_2)
			f0= beta_g*a/3*(1-0.608*Math.exp(-PI*PI*Dg*t/a/a));//Eq.(7) in Ref.1, long times;
		else
			f0 = beta_g*a/3*(6/a*Math.sqrt(Dg*t/PI))-3*Dg*t/a/a;//Eq.(6) in Ref.1, short times;
		return f0;
	}
//
	public double dNdt (double T,double t){
		//double Agb = 4.0* PI *a*a;//the area of the grain boundary of a particular grain, unit: m^2
		//t is time in previous step,s
		//T is temperature, K
		//double f0;//The Booth flux
		//double t1=t+dt;
		double JFlux;
		if(t>0)
			JFlux= 2*f0(T, t)*(1-bgb*delta*Ni/2/Dg(T)/beta_g(T)/t);//Eq.(8) in Ref.1 solved with the Euler scheme;
		else
			JFlux= 2*f0(T, t+dt);
		return JFlux;
	}
//
	public double Nt (double T, double Ni){
		//t is time in previous step,s
		//T is temperature, K
		//Ni is the Nt in previous step
		//double f0;//The Booth flux
		double Np,Nc, Nt;// predicted, corrected values of Nt in improved Euler's scheme
		double t1=local_t+dt;
		double beta_g = beta_g(T);
		double Dg=Dg(T);
		if(local_t>0) {
			Np = 2*dt*f0(T, local_t) + (1-bgb*delta*f0(T, local_t)*dt/Dg/beta_g/local_t)*Ni;
			Nc = 2*dt*f0(T, t1) + Ni-bgb*delta*f0(T, t1)*dt/Dg/beta_g/t1*Np;
			Nt=0.5*Np+0.5*Nc;
				}
		else 
			Nt = 2*dt*f0(T, local_t);	
		//when Nt exceeds Nsat for the first time, GB state is set to 1, and GB is overpressurized.	
		if(Nt>Nsat(T)) {
 			//Nt=Nsat(T);//test 
			GBstate = 1;
		}
		return Nt;
	}
//
	public void C_fgp(double dt,double T) {
		double beta_g = beta_g(T);
		C_fgp += beta_g*4/3*PI*a*a*a/Navgd*dt;		//update total fission gas produced,unit:mol
	}
	public void venting() {
		C_fgr+=Nt*4*PI*a*a/Navgd;//gas released amount is accumulated, unit:mol
		fgrPercolation = C_fgr/C_fgp;//update FGR fraction in the percolation model
		Nt = 0.0;// gas density is reset to zero after venting
		GBstate = 0;//grain boundary state is reset to 'closed'
		local_t = 0.0;// reset the local time following venting
	}
	public double fgrFM(double t,double T) {
		//double omega = DTurbull(T)*t/a/a;
		double omega = Dg(T)*t/a/a;
		double fgrFM;
		if(omega>0.1)
			fgrFM = 1-0.0662*omega*(1-0.9239*Math.exp(-PI*PI*omega));//Eq.(4) in Ref.1
		else if(omega>0&&omega<=0.1)
			fgrFM = 4*Math.sqrt(omega/PI)-1.5*omega; //Eq.(3) in Ref.1
		else
			fgrFM = 0.0;
		return fgrFM;
	}
	public FGR(double dt,int xGB,int xGB_max,double Ni ){
		this.dt = dt;
		this.xGB = xGB_max;
	    this.xGB_max = xGB_max;
	    this.Ni = Ni;   //Grain boudaries gas density at previous time step
	    this.theta = Math.toRadians(60+Math.random()*20);//bubble contact angle that is randomly chosen between 40° and 80°
	    T = this.T(xGB, xGB_max);
	    this.GBstate = 0;
	    this.t = 0.0;
	    this.local_t = 0.0;
	}
//	
	public int fgrCalculate(double Ni){
		this.Ni = Ni;   //Grain boudaries gas density at previous time step		
//		double deltaN;
//			deltaN = dNdt(T,t)*dt;
			t+=dt;
			local_t+=dt;
			C_fgp(dt,T);
			fgrDiffusion=fgrFM(t,T);
			this.Nt = Nt(T,Ni);// use the local time, which restarts following venting
		if(this.Nt>=Nsat(T)) 
			return 1;
		else
			return 0; 
	}
}
