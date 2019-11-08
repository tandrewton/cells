/*

	Methods file for 2D active particles

*/


// include file
#include "deformableParticles2D.h"
#include "cellPacking2D.h"

// namespace
using namespace std;

// constants
const double PI = 4*atan(1);


// FUNCTION TO INITALIZE PARTICLE POSITIONS AS IF THEY WERE ACTIVE SOFT PARTICLES
// WITH DIAMETER SIGMA
// 	** In Pipe flow geometry, to simulate ABP flow for zebrafish study
// 
// 	NOTE: radii are input in units of sigma, the mean particle diameter
// 
// 	** ALSO, use fireMinimizeSP and spAttractiveForces, two functions NOT in this file
void cellPacking2D::initializeActiveStickySP(vector<double>& radii, int NV, double phiDisk, double sizeDispersion, double Lscale){
	// local variables
	int ci, vi, d, nvtmp;
	double r1, r2, g1, radsum;
	double xpos, ypos;
	double xmin, xmax, ymin, ymax;
	double calA;
	double rtmp, l0tmp, a0tmp;
	double delval = 1.0;

	// minimum number of vertices
	const int nvmin = 12;

	// check inputs to stick SP initialization
	if (radii.size() < NCELLS){
		cout << "	** ERROR: in initializing sticky SP, input radii vector size = " << radii.size() << ", which is != NCELLS (= " << NCELLS << "). ending." << endl;
		exit(1);
	}

	// output to console
	cout << "		-- In active stickySP initialization, initializing active SP particles" << endl;

	// initialize length scales as gaussian random variables (will becomes area square roots)
	radsum = 0.0;
	for (ci=0; ci<NCELLS; ci++){
		// generate random numbers
		r1 = drand48();
		r2 = drand48();

		// calculate gaussian random variable using Box-Muller transform
		g1 = sqrt(-2.0*log(r1))*cos(2*PI*r2);

		// get radius
		radii.at(ci) = 0.5*(g1*sizeDispersion + 1.0);

		// add to lenscales sum for boundary size
		radsum += radii.at(ci)*radii.at(ci);
	}

	// determine box length from particle sizes and input packing fraction
	L.at(0) = sqrt(Lscale*PI*radsum/phiDisk);
	L.at(1) = sqrt(PI*radsum/(Lscale*phiDisk));

	// set phi to input
	phi = phiDisk;

	// reseed rng
	srand48(56835698*seed);

	// initialize cell information
	cout << "		-- Ininitializing cell objects" << endl;
	for (ci=0; ci<NCELLS; ci++){
		// boundary information
		for (d=0; d<NDIM; d++)
			cell(ci).setL(d,L.at(d));

		// x-direction is periodic, y-direction is not
		cell(ci).setpbc(0,1);
		cell(ci).setpbc(1,0);

		// number of vertices ( SIGMA SETS # OF VERTS )
		nvtmp = round(2.0*radii.at(ci)*NV);
		if (nvtmp > nvmin)
 			cell(ci).setNV(nvtmp);
		else
			cell(ci).setNV(nvmin);

		// array information
		cell(ci).initializeVertices();
		cell(ci).initializeCell();

		// initial length of polygon side
		l0tmp = 2.0*radii.at(ci)*sin(PI/nvtmp);

		// use rtmp slightly smaller than lenscale, so no overlaps at end
		rtmp = radii.at(ci) - 0.25*delval*l0tmp;

		// calculate a0 and l0 based on fact that they are regular polygons
		a0tmp = 0.5*nvtmp*pow(rtmp,2.0)*sin(2.0*PI/nvtmp);
		l0tmp = 2.0*rtmp*sin(PI/nvtmp);

		// set preferred area and length 
		cell(ci).seta0(a0tmp);
		cell(ci).setl0(l0tmp);
		cell(ci).setdel(delval);
	}

	// initialize particle positions
	cout << "		-- Ininitializing cell positions" << endl;
	for (ci=0; ci<NCELLS; ci++){
		// set min and max values of positions
		xmin = radii.at(ci);
		xmax = L.at(0) - radii.at(ci);
		ymin = radii.at(ci);
		ymax = L.at(1) - radii.at(ci);

		// get random location in pipe
		xpos = (xmax-xmin)*drand48() + xmin;
		ypos = (ymax-ymin)*drand48() + ymin;

		// set as initial position of com
		cell(ci).setCPos(0,xpos);
		cell(ci).setCPos(1,ypos);

		// initialize vertices as a regular polygon
		cell(ci).regularPolygon();
	}

	// initial time scales (t_0^2 = mass*sigma/f_0, mass = 0.25*PI*sigma^2)
	cout << "		-- Ininitializing time scale" << endl;
	dt = 0.05*sqrt(0.25*PI);
	dt0 = dt;

	// use FIRE in hopper geometry to relax overlaps
	cout << "		-- Using FIRE to relax overlaps..." << endl;
	fireMinimizeSP(radii,0.0);
}

void cellPacking2D::spActivePipeForces(vector<double>& radii){
	// local variables
	int ci, cj, vi, d;
	double l1, l2;
	double meanAttract, meanDiam;
	double energyScale;
	double overlap = 0.0;
	double uv = 0.0;
	double ftmp, utmp;
	vector<double> distanceVec(NDIM,0.0);
	double contactDistance = 0.0; 
	double centerDistance = 0.0;

	// fraction of a for l1
	const double l1Scale = 0.5;

	// vector to store number of attractive contacts per cell
	vector<int> nac(NCELLS,0);

	// reset virial stresses to 0
	sigmaXX = 0.0;
	sigmaXY = 0.0;
	sigmaYX = 0.0;
	sigmaYY = 0.0;

	// reset forces
	for (ci=0; ci<NCELLS; ci++){
		for (d=0; d<NDIM; d++)
			cell(ci).setCForce(d,0.0);
	}

	// reset contacts
	resetContacts();

	// get number of attractive contacts
	// loop over cell pairs
	for (ci=0; ci<NCELLS; ci++){
		for (cj=0; cj<NCELLS; cj++){

			// determine if ci and cj engage in attractive bond
			meanAttract = 0.5*(cell(ci).geta()*radii.at(ci) + cell(cj).geta()*radii.at(cj));
			meanDiam = radii.at(ci) + radii.at(cj);

			// get cell-cell distance
			centerDistance = 0.0;
			for (d=0; d<NDIM; d++)
				centerDistance += pow(cell(ci).cellDistance(cell(cj),d),2.0);
			centerDistance = sqrt(centerDistance);

			if (centerDistance < meanAttract + meanDiam)
				addContact(ci,cj);
		}
	}

	// determine number of attractive contacts on each particle
	for (ci=0; ci<NCELLS; ci++){
		for (cj=ci+1; cj<NCELLS; cj++){
			nac.at(ci) += contacts(ci,cj);
			nac.at(cj) += contacts(ci,cj);
		}
	}

	// loop over cells and cell pairs, calculate shape and interaction forces
	for (ci=0; ci<NCELLS; ci++){

		// loop over pairs, add info to contact matrix
		for (cj=ci+1; cj<NCELLS; cj++){
			
			// contact distance
			contactDistance = radii.at(ci) + radii.at(cj);

			// attractive shell
			meanAttract = 0.5*(cell(ci).geta() + cell(cj).geta());

			// center-to-center distance
			centerDistance = 0.0;
			for (d=0; d<NDIM; d++){
				// vectorial quantity
				distanceVec.at(d) = cell(ci).cellDistance(cell(cj),d);

				// add to distance
				centerDistance += pow(distanceVec.at(d),2);
			}
			centerDistance = sqrt(centerDistance);

			// attraction quantities
			l2 = meanAttract;
			l1 = 0.5*l1Scale*l2*(1.0/nac.at(cj) + 1.0/nac.at(ci));

			// check if within interaction zone
			if (centerDistance < meanAttract + contactDistance){
				// overlap scale
				overlap = centerDistance/contactDistance;

				// energy scale
				energyScale = contactDistance;

				// IF in zone to use repulsive force (and, if attractiveParam > 0, bottom of attractive well)
				if (centerDistance < contactDistance*(1 + l1)){
					// interaction potential
					utmp = 0.5 * energyScale * (pow(1 - overlap,2) - l1*l2);

					// scalar part of force
					ftmp = 1.0 - overlap;
				}
				// IF attractiveParam > 0, in attractive well
				else if (centerDistance >= contactDistance*(1 + l1) && centerDistance < contactDistance*(1 + l2) && meanAttract > 0.0){
					// interaction potential
					utmp =  -(0.5*energyScale*l1/(l2 - l1)) * pow(overlap - 1 - l2,2);

					// scalar part of force
					ftmp = (l1/(l2 - l1)) * (overlap - 1.0 - l2);
				}

				// add to u and f based on utmp and ftmp

				// potential energy from utmp
				for (vi=0; vi<cell(ci).getNV(); vi++)
					cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp/cell(ci).getNV());
				for (vi=0; vi<cell(cj).getNV(); vi++)
					cell(cj).setUInt(vi,cell(cj).uInt(vi) + utmp/cell(cj).getNV());

				// forces from ftmp
				for (d=0; d<NDIM; d++){
					// unit vector
					uv = distanceVec.at(d)/centerDistance;

					// add to forces (MIND FORCE DIRECTION; rij points from i -> j, so need extra minus sign)
					cell(ci).setCForce(d,cell(ci).cforce(d) - ftmp*uv);
					cell(cj).setCForce(d,cell(cj).cforce(d) + ftmp*uv);

					// add to virial stresses
					if (d == 0){
						sigmaXX += ftmp*uv*distanceVec.at(0);
						sigmaXY += ftmp*uv*distanceVec.at(1);
					}
					else{
						sigmaYX += ftmp*uv*distanceVec.at(0);
						sigmaYY += ftmp*uv*distanceVec.at(1);
					}
				}
			} 
		}
	}

	// normalize virial stresses by the area
	sigmaXX /= L.at(0)*L.at(1);
	sigmaXY /= L.at(0)*L.at(1);
	sigmaYX /= L.at(0)*L.at(1);
	sigmaYY /= L.at(0)*L.at(1);
}

void cellPacking2D::spActivePipeWallForces(vector<double>& radii){
	// local variables
	int ci, vi;
	double overlap, sigma, x, y, lwy;
	double ftmp, utmp;

	// loop over cells
	for (ci=0; ci<NCELLS; ci++){
		// get sigma (2*radius)
		sigma = 2*radii.at(ci);

		// get particle positions
		x = cell(ci).cpos(0);
		y = cell(ci).cpos(1);

		// if true, interacting with bottom wall
		if (y < 0.5*sigma){
			// vector from wall to particle
			lwy = y;

			// overlap with wall
			overlap = 2.0*lwy/sigma;

			// add to y force ONLY (points in positive y direction)
			ftmp = 1 - overlap;
			cell(ci).setCForce(1,cell(ci).cforce(1) + ftmp);

			// add to energies
			utmp = 0.25*sigma*pow(1 - overlap,2);
			for (vi=0; vi<cell(ci).getNV(); vi++)
				cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp/cell(ci).getNV());

			// virial stress in YY direction
			sigmaYY += ftmp*lwy;
		}

		// if true, interacting with top wall
		if (y > L.at(1) - 0.5*sigma){
			// vector from particle to wall
			lwy = L.at(1) - y;

			// overlap with wall
			overlap = 2.0*lwy/sigma;

			// add to y force ONLY (points in negative y direction)
			ftmp = 1 - overlap;
			cell(ci).setCForce(1,cell(ci).cforce(1) - ftmp);

			// add to energies
			utmp = 0.25*sigma*pow(1 - overlap,2);
			for (vi=0; vi<cell(ci).getNV(); vi++)
				cell(ci).setUInt(vi,cell(ci).uInt(vi) + utmp/cell(ci).getNV());

			// virial stress in YY direction 
			sigmaYY -= ftmp*lwy;
		}
	}
}

void cellPacking2D::spActivePipeNVE(vector<double>& radii, double T0){
	// local variables
	int t, ci, d;
	double Pvirial, K;
	int closed = 1;

	// check that NT has been set 
	if (NT <= 0){
		cout << "	** ERROR: in active pipe SP NVE, sim length NT = " << NT << ", which is <= 0. ending." << endl;
		exit(1);
	}

	// initialize velocities using Gaussian random variables
	vector<double> pmean(NDIM,0.0);
	double r1, r2, grv, vscale, mtmp;

	// loop over velocities, give them initial conditions
	for (ci=0; ci<NCELLS; ci++){
		for (d=0; d<NDIM; d++){
			// draw uniform random variables
			r1 = drand48();
			r2 = drand48();

			// use Box-Muller trnsfrm to get GRV
			grv = sqrt(-2.0*log(r1))*cos(2*PI*r2);

			// add to cell velocity
			cell(ci).setCVel(d,grv);		
		}
	}

	// get system momentum
	for (d=0; d<NDIM; d++)
		pmean.at(d) /= NCELLS;

	// subtract off mean
	K = 0.0;
	for (ci=0; ci<NCELLS; ci++){
		for (d=0; d<NDIM; d++){
			// particle mass
			mtmp = PI*pow(radii.at(ci),2);

			// subtract of com motion
			cell(ci).setCVel(d,cell(ci).cvel(d) - pmean.at(d)/mtmp);

			// calc ek
			K += 0.5*mtmp*pow(cell(ci).cvel(d),2);
		}
	}

	// scale velocities so K to start is T0
	vscale = sqrt(T0/K);
	for (ci=0; ci<NCELLS; ci++){
    	for (d=0; d<NDIM; d++)
        	cell(ci).setCVel(d,cell(ci).cvel(d)*vscale);
    }

	// loop over time, run NVE dynamics
	for (t=0; t<NT; t++){
		// update virial pressure
		Pvirial = 0.5*(sigmaXX + sigmaYY);

		// update kinetic energy based on com velocity
		K = 0.0;
		for (ci=0; ci<NCELLS; ci++)
			K += 0.5*(PI*pow(radii.at(ci),2))*sqrt(pow(cell(ci).cvel(0),2) + pow(cell(ci).cvel(1),2));

		// output some information to console
		if (t % NPRINT == 0){
			cout << "===================================================" << endl << endl;
			cout << " 		Soft, active Particle NVE, t = " << t << endl << endl;
			cout << "===================================================" << endl;

			// print if object has been opened already
			if (packingPrintObject.is_open()){
				cout << "	* Printing SP center positions to file" << endl;
				printPositionsStickySP(radii);
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing SP energy to file" << endl;
				printEnergyStickySP();
			}

			if (statPrintObject.is_open()){
				cout << "	* Printing SP energy to file" << endl;
				printSystemContacts();
			}
			
			cout << "	* Run data:" << endl;
			cout << "	* K 		= " << K << endl;
			cout << "	* virial P 	= " << Pvirial << endl;
			cout << "	* phi 		= " << phi << endl;
			cout << "	* dt 		= " << dt << endl;
			cout << endl << endl;
		}

		// verlet position update
		spPosVerlet();

		// reset contacts before force calculation
		resetContacts();

		// calculate forces between disks
		// spActivePipeForces(radii);
		spAttractiveForces(radii,0.1);
		spActivePipeWallForces(radii);

		// verlet velocity update
		spVelVerlet(radii);
	}
}
/*
void cellPacking2D::spActivePipeFlow(vector<double>& radii, double v0){
	// local variables
	int t, ci, vi;

	// loop over time, integrate overdamped eqn of motion with active motility
	for (t=0; t<NT; t++){
		// update virial pressure
		Pvirial = 0.5*(sigmaXX + sigmaYY);

		// update kinetic energy based on com velocity
		K = 0.0;
		for (ci=0; ci<NCELLS; ci++)
			K += 0.5*(PI*pow(radii.at(ci),2))*sqrt(pow(cell(ci).cvel(0),2) + pow(cell(ci).cvel(1),2));

		// output some information to console
		if (t % NPRINT == 0){
			cout << "===================================================" << endl << endl;
			cout << " 		Soft Particle FLOW, g = " << g << ", t = " << t << endl << endl;
			cout << "===================================================" << endl;

			// print if object has been opened already
			if (packingPrintObject.is_open()){
				cout << "	* Printing SP center positions to file" << endl;
				printHopperSP(radii,w0,w,th,0.0);
			}
			
			if (energyPrintObject.is_open()){
				cout << "	* Printing SP energy to file" << endl;
				printSystemEnergy(t,Pvirial,K);
			}
			
			cout << "	* Run data:" << endl;
			cout << "	* K 		= " << K << endl;
			cout << "	* virial P 	= " << Pvirial << endl;
			cout << "	* phi 		= " << phi << endl;
			cout << "	* dt 		= " << dt << endl;
			cout << "	* g 		= " << g << endl;
			cout << endl << endl;
		}

		// update positions based on forces (EULER)
		for (ci=0; ci<NCELLS; ci++){
			// loop over dimensions, update positions and reset forces for next time
			for (d=0; d<NDIM; d++){
				// get velocities (= forces in overdamped regime)
				veltmp = cell(ci).cvel(d);

				// if new position in outflow region, place back in hopper
				postmp = cell(ci).cpos(d) + dt*veltmp;

				// update positions (EULER STEP)
				cell(ci).setCPos(d,postmp);

				// update velocities
				cell(ci).setCVel(d,veltmp);

				// reset forces
				cell(ci).setCForce(d,0.0);
			}

			// set interaction energy to 0
			for (vi=0; vi<cell(ci).getNV(); vi++)
				cell(ci).setUInt(vi,0.0);
		}

		// reset contacts before force calculation
		resetContacts();

		// calculate forces between disks
		// spActivePipeForces(radii);
		spAttractiveForces(radii,0.1);
		spActivePipeWallForces(radii);

		// update velocities based on forces
		for (ci=0; ci<NCELLS; ci++){
			for (d=0; d<NDIM; d++)
				cell(ci).setCVel(d,cell(ci).cforce(d));
		}
	}
}
*/















