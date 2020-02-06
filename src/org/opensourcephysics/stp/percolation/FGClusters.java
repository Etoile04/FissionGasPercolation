/*
 * Open Source Physics software is free software as described near the bottom of this code file.
 *
 * For additional information and documentation on Open Source Physics please see: 
 * <http://www.opensourcephysics.org/>
 */

package org.opensourcephysics.stp.percolation;

/**
 *  Clusters implements the Newman-Ziff algorithm for identifying clusters.
 *
 *  @author Jan Tobochnik, Wolfgang Christian, Harvey Gould
 *  @version 1.0 06/15/05
 */
public class FGClusters {
  static private final int EMPTY = Integer.MIN_VALUE; // most negative integer
  public int L;                                       // X-axis dimension of lattice, no. of columns
  public int M;                                       // Y-axis dimension of lattice, no. of rows
  public int N;                                       // N = L*M, L is num. of lines, M is num. of rows
  public double dt;
  public int numSitesOccupied;                        // number of occupied lattice sites
  public int[] numClusters;                           // number of clusters of size s, n_s
  public int numOpenGBs;                              // number of open grain boundaries at specific P
  public double accumulatedFGR;						  // amount of fission gas released, unit: mols
  private FGR grainBoundaries[];                       // 
  private double[] Ni;                                // 
  // secondClusterMoment stores sum{s^2 n_s}, where sum is over all clusters (not counting spanning cluster)
  // first cluster moment, sum{s n_s} equals  numSitesOccupied.
  // mean cluster size S is defined as S = secondClusterMoment/numSitesOccupied

  private int secondClusterMoment;
  // spanningClusterSize, number of sites in a spanning cluster; 0 if it doesn't exist
  // assume at most one spanning cluster
  private int spanningClusterSize;
  // gbState[n] gives index of nth site. - 0, closed;1, saturated; 2,vented

  private int[] gbState;
  // parent[] array serves three purposes: stores cluster size when site
  // is root. Otherwise, it stores index of the site's "parent" or is
  // EMPTY. The root is found from an occupied site by recursively following the
  // parent array. Recursion terminates when we encounter a negative value in the
  // parent array, which indicates we have found the unique cluster root.
  // if (parent[s] >= 0) parent[s] is parent site index
  // if (0 > parent[s] > EMPTY) s is root of size -parent[s]
  // if (parent[s] == EMPTY) site s is empty (unoccupied)

  private int[] parent;
  // A spanning cluster touches both left and right boundaries of lattice.
  // As clusters are merged, we maintain this information in following arrays
  // at roots. For example, if root of a cluster is at site 7 and this cluster
  // touches the left side, then touchesLeft[7] == true.

  private boolean[] touchesRght,touchesDown;

  public FGClusters(int L, int M, double dt) {
    this.L = L;
    this.M = M;
    this.dt = dt;
    N = L*M;
    numClusters = new int[N+1];
    gbState = new int[N];
    parent = new int[N];
    Ni =new double [N];
    touchesRght = new boolean[N];
    touchesDown = new boolean[N];
    grainBoundaries = new FGR[N];
  }

  public void newLattice() {
    //setOccupationOrder(); // choose order in which sites are occupied
    // initially all sites are empty, and there are no clusters
    numSitesOccupied = secondClusterMoment = spanningClusterSize = 0;
    accumulatedFGR = 0.0;
    for(int s = 0;s<N;s++) {
      numClusters[s] = 0;
      parent[s] = EMPTY;
      gbState[s] = 0;
      Ni[s] = 0.0;
      //initialize FGR[s] according to grainboudaries properties
      int xGB = s%L;
      int xGB_max = L;
      grainBoundaries[s] = new FGR(dt, xGB, xGB_max, Ni[s]);
    }
    // initially left boundary touchesLeft, right boundary touchesRight
    for(int s = 0;s<N;s++) {
      //touchesLeft[s] = (s%L==0);
      touchesRght[s] = (s%L==L-1);
      touchesDown[s] = (s/L==0);//modified
    }
  }

  // update site state according to FGR diffusion model
  public void updateSite() {
	  for(int s = 0;s<N;s++) {
		  gbState[s] = grainBoundaries[s].fgrCalculate(Ni[s]);
		  // if site s is newly saturated	  
		  if(gbState[s]==1) {
			  // store new cluster's size in parent[]; negative sign distinguishes
			  // newSite as a root, with a size value. Positive values correspond
			  // to non-root sites with index pointers.
			  parent[s] = -1;
			  numClusters[1]++;
			  numSitesOccupied++;
			  secondClusterMoment++;
			  // merge newSite with occupied neighbors.  root is index of merged cluster root at each step
			  int root = s;
			  for(int j = 0;j<4;j++) {
			      // neighborSite is jth site neighboring newly added site newSite
			      int neighborSite = getNeighbor(s, j);
			      if((neighborSite!=EMPTY)&&(parent[neighborSite]!=EMPTY)) {
			        root = mergeRoots(root, findRoot(neighborSite));
			        //if neighboring site is vented to right or down side, so is the newly added site newSite
			        touchesRght[s]=touchesRght[neighborSite];
			        touchesDown[s]=touchesDown[neighborSite];
			      }
			  }
		  }
		  Ni[s] = grainBoundaries[s].Nt;//update Ni[s] for next step 
	  }
	  // assign site venting status according to touchesRght and touchesDown
	  for(int s = 0;s<N;s++) {
		  if (touchesRght[s]==true||touchesDown[s]==true) {
			  gbState[s] = grainBoundaries[s].GBstate = 2;// vented
			  accumulatedFGR += grainBoundaries[s].C_fgr; 
			  Ni[s] =  0.0; // update Ni[s] for next step
			  grainBoundaries[s].venting();
		  }
		  else
			  grainBoundaries[s].Ni = Ni[s]; // update Ni[s] for next step
	  }
  }  
	  /*   // if all sites are occupied, we can't add anymore
    if(numSitesOccupied==N) {
      return;
    }

  }
/*
  public int getOpenGBs() {
	// sum open GBs
	  int j=0;
	  numOpenGBs = 0;
	  int k=0;
	  while(j<N) {
		  if (touchesRght[j]==true||touchesDown[j]==true) {
			  k=findRoot(j);
			  numOpenGBs-= parent[k];
			  j=j-parent[k];
		  }
		  else {
			  j++;
		  }
	  }
	return numOpenGBs;
}
*/

// gets size of  cluster to which site s belongs.
  public int getClusterSize(int s) {
    return(parent[s]==EMPTY) ? 0 : -parent[findRoot(s)];
  }

  // returns size of spanning cluster if it exists, otherwise 0
  public int getSpanningClusterSize() {
    return spanningClusterSize;
  }

  // returns S (mean cluster size); sites belonging to spanning cluster not counted in cluster moments
  public double getMeanClusterSize() {
    int spanSize = getSpanningClusterSize();
    // subtract sites in spanning cluster
    double correctedSecondMoment = secondClusterMoment-spanSize*spanSize;
    double correctedFirstMoment = numSitesOccupied-spanSize;
    if(correctedFirstMoment>0) {
      return correctedSecondMoment/correctedFirstMoment;
    }
	return 0;
  }

  // given a site index s, returns site index representing the root of cluster to which s belongs.
  private int findRoot(int s) {
    if(parent[s]<0) {
      return s; // root site (with size -parent[s])
    }
	// first link parent[s] to the cluster's root to improve performance
      // (path compression); then return this value.
      parent[s] = findRoot(parent[s]);
    return parent[s];
  }

  // returns jth neighbor of site s; j can be 0 (left), 1 (right), 2 (down), or 3
  // (above). If no neighbor exists because of boundary, return value EMPTY.
  // change this method for periodic boundary conditions.
  // compare this method to same method in PercolationApp
  private int getNeighbor(int s, int j) {
    switch(j) {
    case 0 :
      return(s%L==0) ? EMPTY : s-1;   // left
    case 1 :
      return(s%L==L-1) ? EMPTY : s+1; // right
    case 2 :
      return(s/L==0) ? EMPTY : s-L;   // down
    case 3 :
      return(s/L==M-1) ? EMPTY : s+L; // above
    default :
      return EMPTY;
    }
  }
/*
  // fills order[] array with random permutation of site indices. First order[]
  // is set to the identity permutation. Then for values of i in {1...N-1}, swap
  // values of order[i] with order[r], where r is a random index in {i+1 ...N}.
  private void setOccupationOrder() {
    for(int s = 0;s<N;s++) {
      order[s] = s;
    }
    for(int s = 0;s<N-1;s++) {
      int r = s+(int) (Math.random()*(N-s));
      int temp = order[s];
      order[s] = order[r];
      order[r] = temp;
    }
  }
*/
  // utility method to square an integer
  private int sqr(int x) {
    return x*x;
  }

  // merges two root sites into one to represent cluster merging.
  // use heuristic that root of smaller cluster points to root of larger
  // cluster to improve performance.
  // parent[root] stores negative cluster size.
  private int mergeRoots(int r1, int r2) {
    // clusters are uniquely identified by their root sites. If they
    // are the same, clusters are already merged, and we need do nothing
    if(r1==r2) {
      return r1;
      // if r1 has smaller cluster size than r2, reverse (r1,r2) labels
    } else if(-parent[r1]<-parent[r2]) {
      return mergeRoots(r2, r1);//switch the order so that the first cluster is always bigger than the second
    } else { // (-parent[r1] > -parent[r2])
      // update cluster count, and second cluster moment to account for
      // loss of two small clusters and gain of one bigger cluster
      numClusters[-parent[r1]]--;
      numClusters[-parent[r2]]--;
      numClusters[-parent[r1]-parent[r2]]++;
      secondClusterMoment += sqr(parent[r1]+parent[r2])-sqr(parent[r1])-sqr(parent[r2]);
      // cluster at r1 now includes sites of old cluster at r2
      parent[r1] += parent[r2];
      // make r1 new parent of r2
      parent[r2] = r1;
      // if r2 touched down or right, then so does merged cluster r1
      //touchesLeft[r1] |= touchesLeft[r2];
      touchesRght[r1] |= touchesRght[r2];
      touchesDown[r1] |= touchesDown[r2];
      // if cluster at r1 touches Down or Right boundaries, then find the biggest size
      if(touchesDown[r1]||touchesRght[r1]) {
      // if cluster at r1 spans lattice, then remember its size
      /*
       if(touchesLeft[r1]&&touchesRght[r1]) {
       spanningClusterSize = -parent[r1]; 
       */
    	if(spanningClusterSize<-parent[r1]) {
    		spanningClusterSize=-parent[r1];
    	}  
      }
      // return new root site r1
      return r1;
    }
  }
}

/* 
 * Open Source Physics software is free software; you can redistribute
 * it and/or modify it under the terms of the GNU General Public License (GPL) as
 * published by the Free Software Foundation; either version 2 of the License,
 * or(at your option) any later version.

 * Code that uses any portion of the code in the org.opensourcephysics package
 * or any subpackage (subdirectory) of this package must must also be be released
 * under the GNU GPL license.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
 * or view the license online at http://www.gnu.org/copyleft/gpl.html
 *
 * Copyright (c) 2007  The Open Source Physics project
 *                     http://www.opensourcephysics.org
 */
