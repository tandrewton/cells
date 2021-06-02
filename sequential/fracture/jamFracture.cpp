/* 

	MAIN FILE FOR 2D ATTRACTIVE JAMMING CODE VIA ISOTROPIC COMPRESSION
            OF NCELLS USING LINKED-LIST SPEED UP

        CELLS ARE BIDISPERSE, WITH BENDING ENERGY (KB) AND VERTEX-VERTEX
            ATTRACTION OPTIONS

		-- Compresses packing to jamming
		-- Prints position / shape info every dphiPrint
		-- Supports bending energy (kb) and vertex-vertex attraction

        Andrew Ton 
        06/01/2021, majority of algorithm courtesy of Jack Treado

*/

// preprocessor directives
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>

#define NDIM 2
#define NNN 4

// namespace
using namespace std;

// GLOBAL CONSTANTS
const double PI = 4 * atan(1);
const int w = 10;
const int wnum = 25;
const int pnum = 14;

// simulation constants
const int nvmin = 12;
const double timeStepMag = 0.005;
const double phiJMin = 0.6;
const double dphiGrow = 0.001;
const double dphiPrint = 0.01;
const double sizeRatio = 1.4;
const double sizeFraction = 0.5;

// FIRE constants for initial minimizations (SP + DP)
const double alpha0 = 0.2;
const double finc = 1.1;
const double fdec = 0.5;
const double falpha = 0.99;
const double Ftol = 1e-8;

const int NSKIP = 2e4;
const int NMIN = 10;
const int NNEGMAX = 1000;
const int NDELAY = 10;
const int itmax = 5e7;

// DP force constants
const double ka = 1.0;   // area spring (should be = 1)
const double eint = 0.5; // baseline interaction energy
const double del = 1.0;  // radius of vertices in units of l0

// FUNCTION PROTOTYPES

// indexing
int gindex(int ci, int vi, vector<int> &szList);
void cindices(int &ci, int &vi, int gi, int NCELLS, vector<int> &szList);

// particle geometry
double area(vector<double> &vpos, int ci, vector<double> &L, vector<int> &nv, vector<int> &szList);
double perimeter(vector<double> &vpos, int ci, vector<double> &L, vector<int> &nv, vector<int> &szList);

// print to file
void printPos(ofstream &posout, vector<double> &vpos, vector<double> &vrad, vector<double> &a0, vector<double> &calA0, vector<double> &L, vector<int> &cij, vector<int> &nv, vector<int> &szList, double phi, int NCELLS);

// retire the last cell's degree of freedom information from all vectors
void deleteLastCell(int &smallN, int &largeN, int &NCELLS, int &NVTOT, int &cellDOF, int &vertDOF, vector<int> &szList,
                    vector<int> &nv, vector<int> &list, vector<double> &vvel, vector<double> &vpos, vector<double> &vF, vector<int> &im1,
                    vector<int> &ip1, int &vim1, int &vip1, double phi, vector<double> a0, vector<double> l0, vector<double> L, int largeNV, int smallNV);

// MAIN
int main(int argc, char const *argv[])
{
    // variables for indexing loops
    int i, ci, cj, vi, vj, gi, gj, d;

    // parameters to be read in
    int NCELLS, largeNV, smallNV, smallN, largeN, NVTOT, cellDOF, vertDOF, seed, NT;
    double calA0Input, phi, Ptol, phiMin, kl, kb, att, NT_dbl;

    // read in parameters from command line input
    // test: g++ -O3 -std=c++11 sequential/fracture/jamFracture.cpp -o frac.o
    // test: ./frac.o 12 20 1.08 0.4 1e-7 1.0 0 0.01 1 4e5 pos.test shape.test energy.test
    //
    // PARAMETERS:
    // 1. NCELLS 		= # of dpm particles
    // 2. NV 			= # of vertices on smaller cells, larger cells = round(NV*1.4)
    // 3. calA0 		= shape parameter of all cells, initial a0 = 1
    // 4. phiMin 		= min packing fraction (start)
    // 5. Ptol          = pressure tolerance for jamming loop
    // 6. kl 			= perimeter spring constant
    // 7. kb 			= bending spring constant
    // 8. att 			= short-range attraction parameter
    // 9. seed 			= random number generation seed
    // 10. NT           = number of timesteps to simulate
    // 11. positionFile	= position data file string
    // 12. shapeFile	= shape data file string
    // 13. energyFile   = energy data file string

    string NCELLS_str = argv[1];
    string NV_str = argv[2];
    string calA0_str = argv[3];
    string phiMin_str = argv[4];
    string Ptol_str = argv[5];
    string kl_str = argv[6];
    string kb_str = argv[7];
    string att_str = argv[8];
    string seed_str = argv[9];
    string NT_str = argv[10];
    string positionFile = argv[11];
    string shapeFile = argv[12];
    string energyFile = argv[13];

    stringstream NCELLSss(NCELLS_str);
    stringstream NVss(NV_str);
    stringstream calA0ss(calA0_str);
    stringstream phiMinss(phiMin_str);
    stringstream klss(kl_str);
    stringstream kbss(kb_str);
    stringstream Ptolss(Ptol_str);
    stringstream attss(att_str);
    stringstream seedss(seed_str);
    stringstream NTss(NT_str);

    NCELLSss >> NCELLS;
    NVss >> smallNV;
    calA0ss >> calA0Input;
    phiMinss >> phiMin;
    Ptolss >> Ptol;
    klss >> kl;
    kbss >> kb;
    attss >> att;
    seedss >> seed;
    NTss >> NT_dbl;

    //attraction value for jamming
    double att_nve = att;
    att = 0.0;

    // cast input NT_dbl to integer
    NT = (int)NT_dbl;

    // seed random number generator
    srand48(seed);

    // open xyz file
    ofstream posout;
    posout.open(positionFile.c_str());
    if (!posout.is_open())
    {
        cout << "	** ERROR: position file " << positionFile << " could not be opened, ending." << endl;
        return 1;
    }

    ofstream shapeout;
    shapeout.open(shapeFile.c_str());
    if (!shapeout.is_open())
    {
        cout << "	** ERROR: shape file " << shapeFile << " could not be opened, ending." << endl;
        return 1;
    }

    // open jam file, which holds jammed configuration (and whatever subsequent NVE and defect stuff performed on it)
    ofstream jamout;
    string jamFile = positionFile + ".jam";
    jamout.open(jamFile.c_str());
    if (!jamout.is_open())
    {
        cout << "   ** ERROR: jam file " << jamFile << "could not be opened, ending." << endl;
        return 1;
    }

    // open jam file, which holds jammed configuration (and whatever subsequent NVE and defect stuff performed on it)
    ofstream enout;
    enout.open(energyFile.c_str());
    if (!enout.is_open())
    {
        cout << "	** ERROR: energy file " << energyFile << " could not be opened, ending." << endl;
        return 1;
    }

    // total number of vertices
    largeNV = round(sizeRatio * smallNV);

    // total number of vertices
    smallN = round(sizeFraction * NCELLS);
    largeN = NCELLS - smallN;
    NVTOT = smallNV * smallN + largeNV * largeN;

    // szList and nv (keep track of global vertex indices)
    // nv is the number of vertices of cell # ci
    // szList is the cumulative number of vertices from cell # 0 to cell # ci
    cout << "makin vectors" << endl;
    vector<int> szList(NCELLS, 0);
    vector<int> nv(NCELLS, 0);
    nv.at(0) = smallNV;
    for (ci = 1; ci < NCELLS; ci++)
    {
        if (ci < smallN)
        {
            nv.at(ci) = smallNV;
            szList.at(ci) = szList.at(ci - 1) + nv.at(ci - 1);
        }
        else
        {
            nv.at(ci) = largeNV;
            szList.at(ci) = szList.at(ci - 1) + nv.at(ci - 1);
        }
    }

    // degree of freedom counts
    cellDOF = NDIM * NCELLS;
    vertDOF = NDIM * NVTOT;

    // save list of adjacent vertices
    vector<int> im1(NVTOT, 0);
    vector<int> ip1(NVTOT, 0);
    int vim1, vip1;
    for (ci = 0; ci < NCELLS; ci++)
    {
        for (vi = 0; vi < nv.at(ci); vi++)
        {
            // wrap local indices
            vim1 = (vi - 1 + nv.at(ci)) % nv.at(ci);
            vip1 = (vi + 1) % nv.at(ci);

            // get global wrapped indices
            gi = gindex(ci, vi, szList);
            im1.at(gi) = gindex(ci, vim1, szList);
            ip1.at(gi) = gindex(ci, vip1, szList);
        }
    }

    // fundamental MD time units
    double dtMD, dt0, dt;

    dtMD = 1.0;
    dt0 = timeStepMag * dtMD;
    dt = dt0;

    // output opening statement to console
    cout << "=======================================================" << endl
         << endl;

    cout << "		confluence.cpp 								" << endl;
    cout << "		Jack Treado, 2021   							" << endl;

    cout << "		NCELLS 		= " << NCELLS << "					" << endl;
    cout << "		NV (small)	= " << smallNV << "						" << endl;
    cout << "		NV (large)	= " << largeNV << "						" << endl;
    cout << "		NVTOT 		= " << NVTOT << "					" << endl;
    cout << endl;

    cout << "		calA0 		= " << calA0Input << "				" << endl;

    cout << "		phiMin 		= " << phiMin << "					" << endl
         << endl;

    cout << "		kl 			= " << kl << "						" << endl;
    cout << "		kb 			= " << kb << "						" << endl;
    cout << "		att 		= " << att << " 					" << endl;
    cout << "		seed 		= " << seed << "					" << endl
         << endl;

    cout << "		pos file 	= " << positionFile << "			" << endl;
    cout << "		shape file 	= " << shapeFile << "				" << endl
         << endl;

    cout << "=======================================================" << endl
         << endl;

    /* * * * * * * * * * * * * * * * * * 

				INITIALIZE

					PARTICLES

	 * * * * * * * * * * * * * * * * * */

    // initialization variables
    int nvtmp;
    double a0tmp, lenscale, calA0tmp, areaSum = 0.0;

    // initialize vectors for storing coordinates, shape information
    vector<double> vrad(NVTOT, 1.0);
    vector<double> drad(NCELLS, 1.0);

    vector<double> vpos(vertDOF, 0.0);
    vector<double> dpos(cellDOF, 0.0);

    vector<double> a0(NCELLS, 1.0);
    vector<double> l0(NCELLS, 1.0);
    vector<double> calA0(NCELLS, 1.0);

    // shape parameters
    double smallCalA0 = calA0Input * smallNV * tan(PI / smallNV) / PI;
    double largeCalA0 = calA0Input * largeNV * tan(PI / largeNV) / PI;

    // initialize effective disk radius (for minimization), and l0 parameter
    for (ci = 0; ci < NCELLS; ci++)
    {
        // set initial area
        if (ci < smallN)
        {
            lenscale = 1.0;
            a0tmp = 1.0;
            calA0tmp = smallCalA0;
            nvtmp = smallNV;
        }
        else
        {
            lenscale = sizeRatio;
            a0tmp = sizeRatio * sizeRatio;
            calA0tmp = largeCalA0;
            nvtmp = largeNV;
        }
        calA0.at(ci) = calA0tmp;

        // store preferred area
        a0.at(ci) = a0tmp;

        // set disk radius
        drad.at(ci) = 1.1 * sqrt((2.0 * a0tmp) / (nvtmp * sin(2.0 * PI / nvtmp)));

        // set l0, vector radius
        l0.at(ci) = 2.0 * lenscale * sqrt(PI * calA0tmp) / nvtmp;
        gi = szList.at(ci);
        for (vi = 0; vi < nvtmp; vi++)
            vrad.at(gi + vi) = 0.5 * l0.at(ci) * del;

        // add to sum of particle areas (including contribution from vertices)
        areaSum += a0tmp + 0.25 * PI * pow(l0.at(ci) * del, 2.0) * (0.5 * nvtmp - 1);
        cout << "drad = " << drad.at(ci) << ", disk area = " << PI * pow(drad.at(ci), 2.0) << ", l0 = " << l0.at(ci) << ", a0 = " << a0tmp << ", vrad = " << vrad.at(gi) << endl;
    }

    // determine box lengths from particle sizes and input packing fraction
    vector<double> L(NDIM, 1.0);
    for (d = 0; d < NDIM; d++)
        L.at(d) = sqrt(areaSum / phiMin);
    phi = phiMin;

    // initialize cell centers randomly
    for (ci = 0; ci < cellDOF; ci += 2)
        dpos.at(ci) = L[ci % 2] * drand48();
    for (ci = cellDOF - 1; ci > 0; ci -= 2)
        dpos.at(ci) = L[ci % 2] * drand48();

    // initialize contact network
    int NCTCS = 0.5 * NCELLS * (NCELLS - 1);
    vector<int> cij(NCTCS, 0);

    /* * * * * * * * * * * * * * * * * * 

			INITIALIZE

				BOX-LINKED-LIST

	 * * * * * * * * * * * * * * * * * */

    // Cell-linked-list variables
    double boxLengthScale = 3.0;
    double llscale = l0.at(NCELLS - 1);

    // box lengths in each direction
    vector<int> sb(NDIM, 0);
    vector<double> lb(NDIM, 0.0);
    int NBX = 1;
    for (d = 0; d < NDIM; d++)
    {
        // determine number of cells along given dimension by rmax
        sb[d] = round(L[d] / (boxLengthScale * llscale));

        // just in case, if < 3, change to 3 so box neighbor checking will work
        if (sb[d] < 3)
            sb[d] = 3;

        // determine box length by number of cells
        lb[d] = L[d] / sb[d];

        // count total number of cells
        NBX *= sb[d];
    }

    // neighboring boxes for each box (4 neighbors / box)
    int nntmp;
    int scx = sb[0];
    vector<vector<int>> nn;
    nn.resize(NBX);

    // loop over cells, save forward neighbors for each box
    for (i = 0; i < NBX; i++)
    {
        // reshape entry
        nn[i].resize(NNN);

        // neighbors
        nn[i][0] = (i + 1) % NBX;        // right neighbor (i+1)
        nn[i][1] = (i + scx) % NBX;      // top neighbor (j+1)
        nntmp = (i + NBX - scx) % NBX;   // bottom neighbor (j-1)
        nn[i][2] = (nn[i][1] + 1) % NBX; // top-right neighbor (i+1, j+1)
        nn[i][3] = nntmp + 1;            // bottom-right neighbor (i+1, j-1)

        // right-hand bc (periodic)
        if ((i + 1) % scx == 0)
        {
            nn[i][0] = i - scx + 1;
            nn[i][2] = nn[i][1] - scx + 1;
            nn[i][3] = nntmp - scx + 1;
        }
    }

    // linked-list variables
    vector<int> head(NBX, 0);
    vector<int> last(NBX, 0);
    vector<int> list(NVTOT + 1, 0);

    /* * * * * * * * * * * * * * * * * * 

			INITIAL FIRE

					MINIMZATION

	 * * * * * * * * * * * * * * * * * */

    // ----------------------------
    //
    // S P  M I N I M I Z A T I O N
    //
    // ----------------------------

    // initialize disk velocity and force vectors
    vector<double> dv(cellDOF, 0.0);
    vector<double> dF(cellDOF, 0.0);
    vector<double> dFold(cellDOF, 0.0);

    // FIRE VARIABLES
    double P = 0;
    double fnorm = 0;
    double vnorm = 0;
    double alpha = alpha0;

    double dtmax = 10 * dt0;
    double dtmin = 1e-8 * dt0;

    int npPos = 0;
    int npNeg = 0;
    int npPMin = 0;

    int fireit = 0;
    double fcheck = 10 * Ftol;

    // interaction variables
    double xi, yi, xj, yj, dx, dy, fx, fy, rij, sij, ftmp;

    // loop until force relaxes
    while ((fcheck > Ftol || npPMin < NMIN) && fireit < itmax)
    {
        // VV POSITION UPDATE
        for (i = 0; i < cellDOF; i++)
        {
            dpos[i] += dt * dv[i] + 0.5 * dt * dt * dF[i];
            dF[i] = 0;
        }

        // FORCE UPDATE
        for (ci = 0; ci < NCELLS; ci++)
        {
            xi = dpos[NDIM * ci];
            yi = dpos[NDIM * ci + 1];
            for (cj = ci + 1; cj < NCELLS; cj++)
            {
                xj = dpos[NDIM * cj];
                yj = dpos[NDIM * cj + 1];

                // contact distance
                sij = drad[ci] + drad[cj];

                // true distance
                dx = xj - xi;
                dx = dx - L[0] * round(dx / L[0]);
                if (dx < sij)
                {
                    dy = yj - yi;
                    dy = dy - L[1] * round(dy / L[1]);
                    if (dy < sij)
                    {
                        rij = sqrt(dx * dx + dy * dy);
                        if (rij < sij)
                        {
                            ftmp = eint * (1.0 - (rij / sij)) / sij;
                            fx = ftmp * (dx / rij);
                            fy = ftmp * (dy / rij);

                            dF[NDIM * ci] -= fx;
                            dF[NDIM * ci + 1] -= fy;

                            dF[NDIM * cj] += fx;
                            dF[NDIM * cj + 1] += fy;
                        }
                    }
                }
            }
        }

        // VV VELOCITY UPDATE
        for (i = 0; i < cellDOF; i++)
        {
            dv[i] += 0.5 * (dF[i] + dFold[i]) * dt;
            dFold[i] = dF[i];
        }

        // FIRE UPDATE
        // compute fnorm, vnorm and P
        fnorm = 0.0;
        vnorm = 0.0;
        P = 0.0;
        for (i = 0; i < cellDOF; i++)
        {
            fnorm += dF[i] * dF[i];
            vnorm += dv[i] * dv[i];
            P += dv[i] * dF[i];
        }
        fnorm = sqrt(fnorm);
        vnorm = sqrt(vnorm);

        // update fcheck based on fnorm (= force per degree of freedom)
        fcheck = fnorm / (NDIM * NCELLS);

        // update npPMin
        if (fcheck < Ftol && fireit > NDELAY)
            npPMin++;
        else
            npPMin = 0;

        // print to console
        if (fireit % NSKIP == 0)
        {
            cout << endl
                 << endl;
            cout << "===========================================" << endl;
            cout << "		I N I T I A L 				" << endl;
            cout << " 	F I R E 						" << endl;
            cout << "		M I N I M I Z A T I O N 	" << endl;
            cout << "===========================================" << endl;
            cout << endl;
            cout << "	** fireit = " << fireit << endl;
            cout << "	** fcheck = " << fcheck << endl;
            cout << "	** vnorm = " << vnorm << endl;
            cout << "	** dt = " << dt << endl;
            cout << "	** P = " << P << endl;
            cout << "	** Pdir = " << P / (fnorm * vnorm) << endl;
            cout << "	** alpha = " << alpha << endl;
        }

        // Step 1. adjust simulation based on net motion of degrees of freedom
        if (P > 0)
        {
            // increase positive counter
            npPos++;

            // reset negative counter
            npNeg = 0;

            // alter simulation if enough positive steps have been taken
            if (npPos > NMIN)
            {
                // change time step
                if (dt * finc < dtmax)
                    dt *= finc;

                // decrease alpha
                alpha *= falpha;
            }
        }
        else
        {
            // reset positive counter
            npPos = 0;

            // increase negative counter
            npNeg++;

            // check if simulation is stuck
            if (npNeg > NNEGMAX)
            {
                cout << "	** ERROR: During initial FIRE minimization, P < 0 for too long, so ending." << endl;
                return 1;
            }

            // take half step backwards, reset velocities
            for (i = 0; i < cellDOF; i++)
            {
                // take half step backwards
                dpos[i] -= 0.5 * dt * dv[i] + 0.25 * dt * dt * dF[i];

                // reset velocities
                dv[i] = 0.0;
            }

            // decrease time step if past initial delay
            if (fireit > NDELAY)
            {
                // decrease time step
                if (dt * fdec > dtmin)
                    dt *= fdec;

                // reset alpha
                alpha = alpha0;
            }
        }

        // update velocities (s.d. vs inertial dynamics) only if forces are acting
        if (fnorm > 0)
        {
            for (i = 0; i < cellDOF; i++)
                dv[i] = (1 - alpha) * dv[i] + alpha * (vnorm / fnorm) * dF[i];
        }

        // update iterator
        fireit++;
    }
    // check if FIRE converged
    if (fireit == itmax)
    {
        cout << "	** FIRE minimization did not converge, fireit = " << fireit << ", itmax = " << itmax << "; ending." << endl;
        return 1;
    }
    else
    {
        cout << endl
             << endl;
        cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
        cout << "===========================================" << endl;
        cout << " 	F I R E 						" << endl;
        cout << "		M I N I M I Z A T I O N 	" << endl;
        cout << "	C O N V E R G E D! 				" << endl
             << endl;

        cout << "	(for initial disk minimization) " << endl;
        cout << "===========================================" << endl;
        cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
        cout << endl;
        cout << "	** fireit = " << fireit << endl;
        cout << "	** fcheck = " << fcheck << endl;
        cout << "	** vnorm = " << vnorm << endl;
        cout << "	** dt = " << dt << endl;
        cout << "	** P = " << P << endl;
        cout << "	** Pdir = " << P / (fnorm * vnorm) << endl;
        cout << "	** alpha = " << alpha << endl;
    }

    // initialize vertex positions based on cell centers
    for (ci = 0; ci < NCELLS; ci++)
    {
        for (vi = 0; vi < nv.at(ci); vi++)
        {
            // get global vertex index
            gi = gindex(ci, vi, szList);

            // get distance from cell center to vertex if reg poly
            lenscale = sqrt((2.0 * a0.at(ci)) / (nv.at(ci) * sin((2.0 * PI) / nv.at(ci))));

            // set vertex positions
            vpos.at(NDIM * gi) = lenscale * cos((2.0 * PI * vi) / nv.at(ci)) + dpos.at(NDIM * ci) + 1e-2 * l0[ci] * drand48();
            vpos.at(NDIM * gi + 1) = lenscale * sin((2.0 * PI * vi) / nv.at(ci)) + dpos.at(NDIM * ci + 1) + 1e-2 * l0[ci] * drand48();
        }
    }

    /* * * * * * * * * * * * * * * * * * 

		COMPRESS TO

			JAMMING ONSET

	 * * * * * * * * * * * * * * * * * */

    // initialize disk velocity and force vectors
    vector<double> vvel(vertDOF, 0.0);
    vector<double> vF(vertDOF, 0.0);
    vector<double> vFold(vertDOF, 0.0);

    // jamming check variables
    bool undercompressed, overcompressed, jammed;
    int k, kmax, xind, yind;
    double rH, r0, rL, drgrow, drshrink, scaleFactor;
    double pcheck;

    // max number of jamming iteractions
    k = 0;
    kmax = 1e4;

    ///// compute scale factors
    ///double scaleFactor = sqrt((phiMin + dphiGrow) / phiMin);

    // initial boolean variables
    undercompressed = 0;
    overcompressed = 0;
    jammed = 0;

    // jamming bounds
    rH = -1;
    rL = -1;

    // compute scale factors
    drgrow = sqrt((phi + dphiGrow) / phi);
    drshrink = sqrt((phi - 0.5 * dphiGrow) / phi);

    // saved state
    vector<double> vposSave(vertDOF, 0.0);
    vector<double> vradSave(NVTOT, 0.0);
    vector<double> a0Save(NCELLS, 0.0);
    vector<double> l0Save(NCELLS, 0.0);

    // save initial state
    vposSave = vpos;
    vradSave = vrad;
    a0Save = a0;
    l0Save = l0;

    // linked list variables
    int boxid, bi, bj, pi, pj, sbtmp;
    int d0, dend;

    // temporary tolerance (to speed initial compression)
    double Ftoltmp;
    if (Ftol < 0.01 * Ptol)
        Ftoltmp = 0.01 * Ptol;

    // total potential energy
    double U = 0.0;

    // length unit variable
    double rho0 = 0.0;

    // shape force variables
    double fa, fl, fb, l0tmp, atmp, li, lim1, cx, cy;
    double da, dli, dlim1;
    double lim2x, lim2y, lim1x, lim1y, lix, liy, lip1x, lip1y;
    double rim2x, rim2y, rim1x, rim1y, rix, riy, rip1x, rip1y, rip2x, rip2y;
    double cutij, shellij;
    double lastPrintPhi = phi;

    // print frame info
    int printF = 0;

    // compress to jamming, relax U and F using FIRE
    while (!jammed && k < kmax)
    {
        /*//to test energy conservation, break out of loop at lower density and see if velocity verlet implementation is good
        if (phi > 0.5)
            break;*/
        // update iterator
        k++;

        // check to see if cell linked-list needs to be redrawn
        if ((l0[NCELLS - 1] / llscale) > 0.9 * boxLengthScale)
        {
            // print to console
            cout << "\t ** Resetting linked-list: old llscale = " << llscale << ", new llscale = " << l0[NCELLS - 1] << endl;
            cout << "\t ** Resetting linked-list: old NBX = " << NBX;

            // reset linked list for larger particles
            llscale = l0[NCELLS - 1];

            // reset box length info
            fill(sb.begin(), sb.end(), 0);
            fill(lb.begin(), lb.end(), 0.0);
            for (i = 0; i < NBX; i++)
                nn[i].clear();

            NBX = 1;
            for (d = 0; d < NDIM; d++)
            {
                // determine number of cells along given dimension by rmax
                sb[d] = round(L[d] / (boxLengthScale * llscale));

                // just in case, if < 3, change to 3 so box neighbor checking will work
                if (sb[d] < 3)
                    sb[d] = 3;

                // determine box length by number of cells
                lb[d] = L[d] / sb[d];

                // count total number of cells
                NBX *= sb[d];
            }
            cout << ", new NBX = " << NBX << endl;

            // neighboring boxes for each box (4 neighbors / box)
            scx = sb[0];
            nn.clear();
            nn.resize(NBX);

            // loop over cells, save forward neighbors for each box
            for (i = 0; i < NBX; i++)
            {
                // reshape entry
                nn[i].resize(NNN);

                // neighbors
                nn[i][0] = (i + 1) % NBX;        // right neighbor (i+1)
                nn[i][1] = (i + scx) % NBX;      // top neighbor (j+1)
                nntmp = (i + NBX - scx) % NBX;   // bottom neighbor (j-1)
                nn[i][2] = (nn[i][1] + 1) % NBX; // top-right neighbor (i+1, j+1)
                nn[i][3] = nntmp + 1;            // bottom-right neighbor (i+1, j-1)

                // right-hand bc (periodic)
                if ((i + 1) % scx == 0)
                {
                    nn[i][0] = i - scx + 1;
                    nn[i][2] = nn[i][1] - scx + 1;
                    nn[i][3] = nntmp - scx + 1;
                }
            }

            // reset head/last/list
            head.clear();
            head.resize(NBX);

            last.clear();
            last.resize(NBX);

            fill(list.begin(), list.end(), 0);
        }

        // update tolerance
        if (phi > 0.6)
            Ftoltmp = Ftol;

        // RESET FIRE VARIABLES
        P = 0;
        fnorm = 0;
        vnorm = 0;
        alpha = alpha0;

        dtmax = 10.0 * dt0;
        dtmin = 1e-2 * dt0;
        dt = dt0;

        npPos = 0;
        npNeg = 0;
        npPMin = 0;

        fireit = 0;
        fcheck = 10 * Ftoltmp;

        // reset forces
        for (i = 0; i < vertDOF; i++)
        {
            vF[i] = 0.0;
            vFold[i] = 0.0;
            vvel[i] = 0.0;
        }

        // set length constant (variable due to particle growth)
        rho0 = sqrt(a0.at(0));

        // RELAX FORCES USING FIRE
        while (fcheck > Ftoltmp && fireit < itmax)
        {
            // VV POSITION UPDATE
            for (i = 0; i < vertDOF; i++)
            {
                // update position
                vpos[i] += dt * vvel[i] + 0.5 * dt * dt * vF[i];

                // recenter in box
                if (vpos[i] > L[i % NDIM])
                    vpos[i] -= L[i % NDIM];
                else if (vpos[i] < 0)
                    vpos[i] += L[i % NDIM];

                // reset forces
                vF[i] = 0;
            }

            // reset linked list variables
            fill(list.begin(), list.end(), 0);
            fill(head.begin(), head.end(), 0);
            fill(last.begin(), last.end(), 0);

            // sort vertices into linked list
            for (gi = 0; gi < NVTOT; gi++)
            {
                // 1. get cell id of current particle position
                boxid = 0;
                sbtmp = 1;
                for (d = 0; d < NDIM; d++)
                {
                    // add d index to 1d list
                    boxid += floor(vpos[NDIM * gi + d] / lb[d]) * sbtmp;

                    // increment dimensional factor
                    sbtmp *= sb[d];
                }

                // 2. add to head list or link within list
                // NOTE: particle ids are labelled starting from 1, setting to 0 means end of linked list
                if (head[boxid] == 0)
                {
                    head[boxid] = gi + 1;
                    last[boxid] = gi + 1;
                }
                else
                {
                    list[last[boxid]] = gi + 1;
                    last[boxid] = gi + 1;
                }
            }

            // reset contact networks
            fill(cij.begin(), cij.end(), 0);

            // FORCE UPDATE

            // interaction forces (USE BOX LINKED LIST)
            pcheck = 0.0;
            U = 0.0;
            for (bi = 0; bi < NBX; bi++)
            {

                // get start of list of particles
                pi = head[bi];

                // loop over linked list
                while (pi > 0)
                {
                    // real particle index
                    gi = pi - 1;

                    // next particle in list
                    pj = list[pi];

                    // loop down neighbors of pi in same cell
                    while (pj > 0)
                    {
                        // real index of pj
                        gj = pj - 1;

                        if (gj == ip1[gi] || gj == im1[gi])
                        {
                            pj = list[pj];
                            continue;
                        }

                        // contact distance
                        sij = vrad[gi] + vrad[gj];

                        // attraction distances
                        shellij = (1.0 + att) * sij;
                        cutij = (1.0 + 0.5 * att) * sij;

                        // particle distance
                        dx = vpos[NDIM * gj] - vpos[NDIM * gi];
                        dx -= L[0] * round(dx / L[0]);
                        if (dx < shellij)
                        {
                            dy = vpos[NDIM * gj + 1] - vpos[NDIM * gi + 1];
                            dy -= L[1] * round(dy / L[1]);
                            if (dy < shellij)
                            {
                                rij = sqrt(dx * dx + dy * dy);
                                // if two tumor cells, check for attraction
                                if (rij < shellij && rij > cutij && att > 0.0)
                                {
                                    // force scale
                                    ftmp = eint * ((rij / sij) - 1.0 - att) / sij;
                                    fx = ftmp * (dx / rij);
                                    fy = ftmp * (dy / rij);

                                    // add to forces
                                    vF[NDIM * gi] -= fx;
                                    vF[NDIM * gi + 1] -= fy;

                                    vF[NDIM * gj] += fx;
                                    vF[NDIM * gj + 1] += fy;

                                    // add to virial expression for pressure
                                    pcheck += dx * fx + dy * fy;

                                    // add to potential energy
                                    U -= 0.5 * eint * pow((rij / sij) - 1.0 - att, 2.0);

                                    // add to contacts
                                    cindices(ci, vi, gi, NCELLS, szList);
                                    cindices(cj, vj, gj, NCELLS, szList);

                                    if (ci > cj)
                                        cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2]++;
                                    else if (ci < cj)
                                        cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2]++;
                                }
                                else if (rij < cutij && rij > sij && att > 0.0)
                                {
                                    // force scale
                                    ftmp = eint * (1 - (rij / sij)) / sij;
                                    fx = ftmp * (dx / rij);
                                    fy = ftmp * (dy / rij);

                                    // add to forces
                                    vF[NDIM * gi] -= fx;
                                    vF[NDIM * gi + 1] -= fy;

                                    vF[NDIM * gj] += fx;
                                    vF[NDIM * gj + 1] += fy;

                                    // add to virial expression for pressure
                                    pcheck += dx * fx + dy * fy;

                                    // add to potential energy
                                    U += 0.5 * eint * (pow(1.0 - (rij / sij), 2.0) - 0.5 * att * att);

                                    // add to contacts
                                    cindices(ci, vi, gi, NCELLS, szList);
                                    cindices(cj, vj, gj, NCELLS, szList);

                                    if (ci > cj)
                                        cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2]++;
                                    else if (ci < cj)
                                        cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2]++;
                                }
                                // otherwise, check for two overlapping purely-repulsive vertices
                                else if (rij < sij)
                                {
                                    // force scale
                                    ftmp = eint * (1 - (rij / sij)) / sij;
                                    fx = ftmp * (dx / rij);
                                    fy = ftmp * (dy / rij);

                                    // add to forces
                                    vF[NDIM * gi] -= fx;
                                    vF[NDIM * gi + 1] -= fy;

                                    vF[NDIM * gj] += fx;
                                    vF[NDIM * gj + 1] += fy;

                                    // add to virial expression for pressure
                                    pcheck += dx * fx + dy * fy;

                                    // add to potential energy
                                    U += 0.5 * eint * (pow(1.0 - (rij / sij), 2.0) - 0.5 * att * att);

                                    // add to contacts
                                    cindices(ci, vi, gi, NCELLS, szList);
                                    cindices(cj, vj, gj, NCELLS, szList);

                                    if (ci > cj)
                                        cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2]++;
                                    else if (ci < cj)
                                        cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2]++;
                                }
                            }
                        }

                        // update pj
                        pj = list[pj];
                    }

                    // test overlaps with forward neighboring cells
                    for (bj = 0; bj < NNN; bj++)
                    {
                        // get first particle in neighboring cell
                        pj = head[nn[bi][bj]];

                        // loop down neighbors of pi in same cell
                        while (pj > 0)
                        {
                            // real index of pj
                            gj = pj - 1;

                            if (gj == ip1[gi] || gj == im1[gi])
                            {
                                pj = list[pj];
                                continue;
                            }
                            // contact distance
                            sij = vrad[gi] + vrad[gj];

                            // attraction distances
                            shellij = (1.0 + att) * sij;
                            cutij = (1.0 + 0.5 * att) * sij;

                            // particle distance
                            dx = vpos[NDIM * gj] - vpos[NDIM * gi];
                            dx -= L[0] * round(dx / L[0]);
                            if (dx < shellij)
                            {
                                dy = vpos[NDIM * gj + 1] - vpos[NDIM * gi + 1];
                                dy -= L[1] * round(dy / L[1]);
                                if (dy < shellij)
                                {
                                    rij = sqrt(dx * dx + dy * dy);
                                    // if two tumor cells, check for attraction
                                    if (rij < shellij && rij > cutij && att > 0.0)
                                    {
                                        // force scale
                                        ftmp = eint * ((rij / sij) - 1.0 - att) / sij;
                                        fx = ftmp * (dx / rij);
                                        fy = ftmp * (dy / rij);

                                        // add to forces
                                        vF[NDIM * gi] -= fx;
                                        vF[NDIM * gi + 1] -= fy;

                                        vF[NDIM * gj] += fx;
                                        vF[NDIM * gj + 1] += fy;

                                        // add to virial expression for pressure
                                        pcheck += dx * fx + dy * fy;

                                        // add to potential energy
                                        U -= 0.5 * eint * pow((rij / sij) - 1.0 - att, 2.0);

                                        // add to contacts
                                        cindices(ci, vi, gi, NCELLS, szList);
                                        cindices(cj, vj, gj, NCELLS, szList);

                                        if (ci > cj)
                                            cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2]++;
                                        else if (ci < cj)
                                            cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2]++;
                                    }
                                    else if (rij < cutij && rij > sij && att > 0.0)
                                    {
                                        // force scale
                                        ftmp = eint * (1 - (rij / sij)) / sij;
                                        fx = ftmp * (dx / rij);
                                        fy = ftmp * (dy / rij);

                                        // add to forces
                                        vF[NDIM * gi] -= fx;
                                        vF[NDIM * gi + 1] -= fy;

                                        vF[NDIM * gj] += fx;
                                        vF[NDIM * gj + 1] += fy;

                                        // add to virial expression for pressure
                                        pcheck += dx * fx + dy * fy;

                                        // add to potential energy
                                        U += 0.5 * eint * (pow(1.0 - (rij / sij), 2.0) - 0.5 * att * att);

                                        // add to contacts
                                        cindices(ci, vi, gi, NCELLS, szList);
                                        cindices(cj, vj, gj, NCELLS, szList);

                                        if (ci > cj)
                                            cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2]++;
                                        else if (ci < cj)
                                            cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2]++;
                                    }
                                    // otherwise, check for two overlapping purely-repulsive vertices
                                    else if (rij < sij)
                                    {
                                        // force scale
                                        ftmp = eint * (1 - (rij / sij)) / sij;
                                        fx = ftmp * (dx / rij);
                                        fy = ftmp * (dy / rij);

                                        // add to forces
                                        vF[NDIM * gi] -= fx;
                                        vF[NDIM * gi + 1] -= fy;

                                        vF[NDIM * gj] += fx;
                                        vF[NDIM * gj + 1] += fy;

                                        // add to virial expression for pressure
                                        pcheck += dx * fx + dy * fy;

                                        // add to potential energy
                                        U += 0.5 * eint * (pow(1.0 - (rij / sij), 2.0) - 0.5 * att * att);

                                        // add to contacts
                                        cindices(ci, vi, gi, NCELLS, szList);
                                        cindices(cj, vj, gj, NCELLS, szList);

                                        if (ci > cj)
                                            cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2]++;
                                        else if (ci < cj)
                                            cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2]++;
                                    }
                                }
                            }

                            // update pj
                            pj = list[pj];
                        }
                    }

                    // update pi index to be next
                    pi = list[pi];
                }
            }

            // normalize pressure by box area (make dimensionless with extra factor of rho)
            pcheck *= (rho0 / (2.0 * L[0] * L[1]));

            // shape forces (loop over global vertex labels)
            ci = 0;
            for (gi = 0; gi < NVTOT; gi++)
            {

                // -- Area force (and get cell index ci)
                if (ci < NCELLS)
                {
                    if (gi == szList[ci])
                    {
                        // compute shape parameter
                        nvtmp = nv[ci];
                        a0tmp = a0[ci];
                        l0tmp = l0[ci];

                        // compute area deviation
                        atmp = area(vpos, ci, L, nv, szList);
                        da = (atmp / a0tmp) - 1.0;

                        // add to potential energy
                        U += 0.5 * pow(da, 2.0);

                        // shape force parameters (kl and kl are unitless energy ratios)
                        fa = da * (rho0 / a0tmp); // derivation from the fact that rho0^2 does not necessarily cancel a0tmp
                        fl = kl * (rho0 / l0tmp);
                        fb = kb * (rho0 / (l0tmp * l0tmp));

                        // compute cell center of mass
                        xi = vpos[NDIM * gi];
                        yi = vpos[NDIM * gi + 1];
                        cx = xi;
                        cy = yi;
                        for (vi = 1; vi < nvtmp; vi++)
                        {
                            dx = vpos.at(NDIM * (gi + vi)) - xi;
                            dx -= L[0] * round(dx / L[0]);

                            dy = vpos.at(NDIM * (gi + vi) + 1) - yi;
                            dy -= L[1] * round(dy / L[1]);

                            xi += dx;
                            yi += dy;

                            cx += xi;
                            cy += yi;
                        }
                        cx /= nvtmp;
                        cy /= nvtmp;

                        // get coordinates relative to center of mass
                        rix = vpos[NDIM * gi] - cx;
                        riy = vpos[NDIM * gi + 1] - cy;

                        // get (prior) adjacent vertices
                        rim1x = vpos[NDIM * im1[gi]] - cx;
                        rim1x -= L[0] * round(rim1x / L[0]);

                        rim1y = vpos[NDIM * im1[gi] + 1] - cy;
                        rim1y -= L[1] * round(rim1y / L[1]);

                        rim2x = vpos[NDIM * im1[im1[gi]]] - cx;
                        rim2x -= L[0] * round(rim2x / L[0]);

                        rim2y = vpos[NDIM * im1[im1[gi]] + 1] - cy;
                        rim2y -= L[1] * round(rim2y / L[1]);

                        // increment cell index
                        ci++;
                    }
                }

                // get next adjacent vertices
                rip1x = vpos.at(NDIM * ip1[gi]) - cx;
                rip1x -= L[0] * round(rip1x / L[0]);

                rip1y = vpos.at(NDIM * ip1[gi] + 1) - cy;
                rip1y -= L[1] * round(rip1y / L[1]);

                // -- Area force
                vF[NDIM * gi] += 0.5 * fa * (rim1y - rip1y);
                vF[NDIM * gi + 1] += 0.5 * fa * (rip1x - rim1x);

                // -- Perimeter force

                // segment vector elements
                lim1x = rix - rim1x;
                lim1y = riy - rim1y;

                lix = rip1x - rix;
                liy = rip1y - riy;

                // segment lengths
                lim1 = sqrt(lim1x * lim1x + lim1y * lim1y);
                li = sqrt(lix * lix + liy * liy);

                // segment deviations
                dlim1 = (lim1 / l0tmp) - 1.0;
                dli = (li / l0tmp) - 1.0;

                // add to potential energy
                U += 0.5 * kl * pow(dli, 2.0);

                // add to forces
                vF[NDIM * gi] += fl * (dli * (lix / li) - dlim1 * (lim1x / lim1));
                vF[NDIM * gi + 1] += fl * (dli * (liy / li) - dlim1 * (lim1y / lim1));

                // -- Bending force
                if (kb > 0)
                {
                    // segment vectors for ip2
                    rip2x = vpos[NDIM * ip1[ip1[gi]]] - cx;
                    rip2x -= L[0] * round(rip2x / L[0]);

                    rip2y = vpos[NDIM * ip1[ip1[gi]] + 1] - cy;
                    rip2y -= L[1] * round(rip2y / L[1]);

                    lip1x = rip2x - rip1x;
                    lip1y = rip2y - rip1y;

                    lim2x = rim1x - rim2x;
                    lim2y = rim1y - rim2y;

                    // add to potential energy
                    U += 0.5 * kb * (pow(lix - lim1x, 2.0) + pow(liy - lim1y, 2.0)) / (l0tmp * l0tmp);

                    // add to force
                    vF[NDIM * gi] += fb * (3.0 * (lix - lim1x) + lim2x - lip1x);
                    vF[NDIM * gi + 1] += fb * (3.0 * (liy - lim1y) + lim2y - lip1y);
                }

                // update old coordinates
                rim2x = rim1x;
                rim1x = rix;
                rix = rip1x;

                rim2y = rim1y;
                rim1y = riy;
                riy = rip1y;
            }

            // VV VELOCITY UPDATE
            for (i = 0; i < vertDOF; i++)
            {
                vvel[i] += 0.5 * (vF[i] + vFold[i]) * dt;
                vFold[i] = vF[i];
            }

            // FIRE UPDATE
            // compute fnorm, vnorm and P
            fnorm = 0.0;
            vnorm = 0.0;
            P = 0.0;
            for (i = 0; i < vertDOF; i++)
            {
                fnorm += vF[i] * vF[i];
                vnorm += vvel[i] * vvel[i];
                P += vvel[i] * vF[i];
            }
            fnorm = sqrt(fnorm);
            vnorm = sqrt(vnorm);

            // update fcheck based on fnorm (= force per degree of freedom)
            fcheck = fnorm / sqrt(NCELLS);

            // update npPMin
            if (fcheck < Ftol)
                npPMin++;
            else
                npPMin = 0;

            // print to console
            if (fireit % NSKIP == 0)
            {
                cout << endl
                     << endl;
                cout << "===========================================" << endl;
                cout << " 	F I R E 						" << endl;
                cout << "		M I N I M I Z A T I O N 	" << endl;
                cout << "===========================================" << endl;
                cout << endl;
                cout << "	** fireit 	= " << fireit << endl;
                cout << "	** fcheck 	= " << fcheck << endl;
                cout << "	** pcheck 	= " << pcheck << endl;
                cout << "	** U 		= " << U << endl;
                cout << "   ** Ftoltmp     = " << Ftoltmp << endl;
                cout << "   ** Ptol     = " << Ptol << endl;
                cout << "	** vnorm 	= " << vnorm << endl;
                cout << "	** dt 		= " << dt << endl;
                cout << "	** P 		= " << P << endl;
                cout << "	** Pdir 	= " << P / (fnorm * vnorm) << endl;
                cout << "	** alpha 	= " << alpha << endl;
            }

            // Step 1. adjust simulation based on net motion of degrees of freedom
            if (P > 0)
            {
                // increase positive counter
                npPos++;

                // reset negative counter
                npNeg = 0;

                // alter simulation if enough positive steps have been taken
                if (npPos > NDELAY)
                {
                    // change time step
                    if (dt * finc < dtmax)
                        dt *= finc;

                    // decrease alpha
                    alpha *= falpha;
                }
            }
            else
            {
                // reset positive counter
                npPos = 0;

                // increase negative counter
                npNeg++;

                // check if simulation is stuck
                if (npNeg > NNEGMAX)
                {
                    cout << "	** ERROR: During initial FIRE minimization, P < 0 for too long, so ending." << endl;
                    return 1;
                }

                // take half step backwards, reset velocities
                for (i = 0; i < vertDOF; i++)
                {
                    // take half step backwards
                    vpos[i] -= 0.5 * dt * vvel[i] + 0.25 * dt * dt * vF[i];

                    // reset vertex velocities
                    vvel[i] = 0.0;
                }

                // decrease time step if past initial delay
                if (fireit > NDELAY)
                {
                    // decrease time step
                    if (dt * fdec > dtmin)
                        dt *= fdec;

                    // reset alpha
                    alpha = alpha0;
                }
            }

            // update velocities (s.d. vs inertial dynamics) only if forces are acting
            if (fnorm > 0)
            {
                for (i = 0; i < vertDOF; i++)
                    vvel[i] = (1 - alpha) * vvel[i] + alpha * (vnorm / fnorm) * vF[i];
            }

            // update iterator
            fireit++;
        }
        // check if FIRE converged
        if (fireit == itmax)
        {
            cout << "	** FIRE minimization did not converge, fireit = " << fireit << ", itmax = " << itmax << "; ending." << endl;
            return 1;
        }
        else
        {
            cout << endl;
            cout << "===========================================" << endl;
            cout << " 	F I R E 						" << endl;
            cout << "		M I N I M I Z A T I O N 	" << endl;
            cout << "	C O N V E R G E D! 				" << endl;
            cout << "	** at k = " << k << " 			" << endl;
            cout << "   ** line 1436 in code " << endl;
            cout << "===========================================" << endl;
            cout << endl;
            cout << "	** fireit 	= " << fireit << endl;
            cout << "	** fcheck 	= " << fcheck << endl;
            cout << "	** pcheck 	= " << pcheck << endl;
            cout << "	** U 		= " << U << endl;

            cout << "	** vnorm 	= " << vnorm << endl;
            cout << "	** dt 		= " << dt << endl;
            cout << "	** P 		= " << P << endl;
            cout << "	** Pdir 	= " << P / (fnorm * vnorm) << endl;
            cout << "	** alpha 	= " << alpha << endl;
            cout << endl
                 << endl;
        }

        // boolean check for jamming
        undercompressed = ((pcheck < 2.0 * Ptol && rH < 0) || (pcheck < Ptol && rH > 0));
        overcompressed = (pcheck > 2.0 * Ptol && phi > phiJMin);
        jammed = (pcheck < 2.0 * Ptol && pcheck > Ptol && rH > 0 && rL > 0);

        // output to console
        cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
        cout << "===============================================" << endl
             << endl;
        cout << " 	Q U A S I S T A T I C  				" << endl;
        cout << " 	  	I S O T R O P I C 				" << endl;
        cout << "			C O M P R E S S I O N 		" << endl
             << endl;
        cout << "===============================================" << endl;
        cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
        cout << endl;
        cout << "	* k 			= " << k << endl;
        cout << "	* scaleFactor 	= " << scaleFactor << endl;
        cout << "	* phi 			= " << phi << endl;
        cout << "	* r0 			= " << r0 << endl;
        cout << "	* rH 			= " << rH << endl;
        cout << "	* rL 			= " << rL << endl;
        cout << "	* fcheck 		= " << fcheck << endl;
        cout << "	* pcheck 		= " << pcheck << endl;
        cout << "	* U 		 	= " << U << endl;
        cout << "	* undercompressed = " << undercompressed << endl;
        cout << "	* overcompressed = " << overcompressed << endl;
        cout << "	* jammed = " << jammed << endl
             << endl;
        cout << endl;

        // update particle sizes based on target check
        if (rH < 0)
        {
            // if still undercompressed, then grow until overcompressed found
            if (undercompressed)
            {
                r0 = rho0;
                scaleFactor = sqrt((phi + dphiGrow) / phi);
            }
            // if first overcompressed, decompress by dphi/2 until unjamming
            else if (overcompressed)
            {
                // current = upper bound length scale r
                rH = rho0;

                // save first overcompressed state
                r0 = rH;
                vposSave = vpos;
                vradSave = vrad;
                a0Save = a0;
                l0Save = l0;

                // compute new scale factor
                scaleFactor = drshrink;

                // print to console
                cout << "	-- -- overcompressed for the first time, scaleFactor = " << scaleFactor << endl;
            }
            else
            {
                std::cout << "error: undercompressed and overcompressed flags are not set.\n";
                return 1;
            }
        }
        else
        {
            if (rL < 0)
            {
                // if first undercompressed, save last overcompressed state, beging root search
                if (undercompressed)
                {
                    // current = new lower bound length scale r
                    rL = rho0;

                    // load state
                    vpos = vposSave;
                    vrad = vradSave;
                    a0 = a0Save;
                    l0 = l0Save;

                    // compute new scale factor by root search
                    scaleFactor = 0.5 * (rH + rL) / r0;

                    // print to console
                    cout << "	-- -- undercompressed for the first time, scaleFactor = " << scaleFactor << endl;
                    cout << "	-- -- BEGINNING ROOT SEARCH IN ENTHALPY MIN PROTOCOL..." << endl;
                }
                // if still overcompressed, decrement again
                else if (overcompressed)
                {
                    // current = upper bound length scale r
                    rH = rho0;

                    // save overcompressed state
                    r0 = rH;
                    vposSave = vpos;
                    vradSave = vrad;
                    a0Save = a0;
                    l0Save = l0;

                    // keep shrinking at same rate until unjamming
                    scaleFactor = drshrink;

                    // print to console
                    cout << "	-- -- overcompressed, still no unjamming, scaleFactor = " << scaleFactor << endl;
                }
                else
                {
                    std::cout << "error: undercompressed and overcompressed flags are not set.\n";
                    return 1;
                }
            }
            else
            {
                // if found undercompressed state, go to state between undercompressed and last overcompressed states (from saved state)
                if (undercompressed)
                {
                    // current = new lower bound length scale r
                    rL = rho0;

                    // load state
                    vpos = vposSave;
                    vrad = vradSave;
                    a0 = a0Save;
                    l0 = l0Save;

                    // compute new scale factor
                    scaleFactor = 0.5 * (rH + rL) / r0;

                    // print to console
                    cout << "	-- -- undercompressed, scaleFactor = " << scaleFactor << endl;
                }
                else if (overcompressed)
                {
                    // current = new upper bound length scale r
                    rH = rho0;

                    // load state
                    vpos = vposSave;
                    vrad = vradSave;
                    a0 = a0Save;
                    l0 = l0Save;

                    // compute new scale factor
                    scaleFactor = 0.5 * (rH + rL) / r0;

                    // print to console
                    cout << "	-- -- overcompressed, scaleFactor = " << scaleFactor << endl;
                }
                else if (jammed)
                {
                    cout << "	** At k = " << k << ", target pressure found!" << endl;
                    cout << "	** fcheck = " << fcheck << endl;
                    cout << "	** pcheck = " << pcheck << endl;
                    cout << "	** U = " << U << endl;
                    cout << " NOT WRITING ENTHALPY-MINIMIZED CONFIG TO .jam FILE" << endl;
                    cout << " ENDING COMPRESSION SIMULATION" << endl;
                    //printPos(jamout, vpos, vrad, a0, calA0, L, cij, nv, szList, phi, NCELLS);
                    break;
                }
                else
                {
                    std::cout << "error: neither under nor overcompressed nor jammed\n";
                    return 1;
                }
            }
        }

        // print position, energy + shape information
        if (abs(lastPrintPhi - phi) > dphiPrint)
        {
            // print shape info
            cout << "\t** PRINTING SHAPE INFO TO FILE... " << endl;
            shapeout << setw(6) << printF;
            shapeout << setw(15) << fireit;
            shapeout << setw(15) << phi;
            shapeout << setw(15) << pcheck;
            shapeout << setw(15) << U;
            for (ci = 0; ci < NCELLS; ci++)
            {
                shapeout << setw(15) << perimeter(vpos, ci, L, nv, szList);
                shapeout << setw(15) << area(vpos, ci, L, nv, szList);
            }
            shapeout << endl;

            // print position info
            cout << "\t** PRINTING POSITIONS TO FILE... " << endl;
            printPos(posout, vpos, vrad, a0, calA0, L, cij, nv, szList, phi, NCELLS);

            // update last phi when printed, for next time
            lastPrintPhi = phi;
            printF++;
        }

        // grow or shrink particles by scale factor
        phi = 0.0;
        for (ci = 0; ci < NCELLS; ci++)
        {
            // scale preferred lengths
            l0[ci] *= scaleFactor;
            a0[ci] *= scaleFactor * scaleFactor;

            // first global index for ci
            gi = szList.at(ci);

            // compute cell center of mass
            xi = vpos[NDIM * gi];
            yi = vpos[NDIM * gi + 1];
            cx = xi;
            cy = yi;
            for (vi = 1; vi < nv.at(ci); vi++)
            {
                dx = vpos.at(NDIM * (gi + vi)) - xi;
                dx -= L[0] * round(dx / L[0]);

                dy = vpos.at(NDIM * (gi + vi) + 1) - yi;
                dy -= L[1] * round(dy / L[1]);

                xi += dx;
                yi += dy;

                cx += xi;
                cy += yi;
            }
            cx /= nv.at(ci);
            cy /= nv.at(ci);

            for (vi = 0; vi < nv.at(ci); vi++)
            {
                // x and y inds
                xind = NDIM * (gi + vi);
                yind = xind + 1;

                // closest relative position
                dx = vpos[xind] - cx;
                dx -= L[0] * round(dx / L[0]);

                dy = vpos[yind] - cy;
                dy -= L[1] * round(dy / L[1]);

                // update vertex positions
                vpos[xind] += (scaleFactor - 1.0) * dx;
                vpos[yind] += (scaleFactor - 1.0) * dy;

                // scale vertex radii
                vrad[gi + vi] *= scaleFactor;
            }

            // update packing fraction
            phi += a0[ci] + 0.25 * PI * pow(l0[ci] * del, 2.0) * (0.5 * nv[ci] - 1);
        }
        phi /= L[0] * L[1];
        scaleFactor = sqrt((phi + dphiGrow) / phi);
    }
    if (k == kmax)
    {
        cout << "	** ERROR: IN 2d cell jamming finding, k reached kmax without finding jamming. Ending." << endl;
        return 1;
    }

    // print shape info
    cout << "\t** PRINTING SHAPE INFO TO FILE... " << endl;
    shapeout << setw(6) << printF;
    shapeout << setw(15) << fireit;
    shapeout << setw(15) << phi;
    shapeout << setw(15) << pcheck;
    shapeout << setw(15) << U;
    for (ci = 0; ci < NCELLS; ci++)
    {
        shapeout << setw(15) << perimeter(vpos, ci, L, nv, szList);
        shapeout << setw(15) << area(vpos, ci, L, nv, szList);
    }
    shapeout << endl;

    // print position info
    cout << "\t** PRINTING POSITIONS TO FILE... " << endl;
    printPos(posout, vpos, vrad, a0, calA0, L, cij, nv, szList, phi, NCELLS);

    /* * * * * * * * * * * * * * * * * * 

			NVE

				DYNAMICS

	 * * * * * * * * * * * * * * * * * */

    // DEBUG MENU
    cout << endl
         << endl;
    cout << "===========================================" << endl;
    cout << "			N V E " << endl;
    cout << "===========================================" << endl;
    cout << endl;
    cout << "	** tt = n/a" << endl;
    cout << "	** smallN = " << smallN << endl;
    cout << "	** largeN = " << largeN << endl;
    cout << "	** NCELLS = " << NCELLS << endl;
    cout << "	** NVTOT = " << NVTOT << endl;
    cout << "	** cellDOF = " << cellDOF << endl;
    cout << "	** vertDOF = " << vertDOF << endl;
    cout << "	** dim szList = " << szList.size() << endl;
    cout << "	** dim nv = " << nv.size() << endl;
    cout << "	** dim list = " << list.size() << endl;
    cout << "	** dim vvel = " << vvel.size() << endl;
    cout << "	** dim vpos = " << vpos.size() << endl;
    cout << "	** dim vF = " << vF.size() << endl;

    // NVE VARIABLES
    int tt;
    double K;
    att = att_nve;

    // reset for NVE
    dt = dt0;

    //initialize velocities according to temperature T0 and zero total linear momentum.
    double T0 = 0.0001;
    double v_cm_x = 0.0;
    double v_cm_y = 0.0;
    for (i = 0; i < vertDOF; i++)
    {
        // box muller transform
        double r1 = drand48();
        double r2 = drand48();
        double grv = sqrt(-2.0 * log(r1)) * cos(2.0 * PI * r2);

        // draw random velocity
        vvel[i] = sqrt(T0) * grv;
        if (i % NDIM == 0)
        {
            v_cm_x += vvel[i];
        }
        else if (i % NDIM == 1)
        {
            v_cm_y += vvel[i];
        }
        else
        {
            cout << "Error with number of dimensions in velocity initialization, exiting...\n";
            return 1;
        }
    }
    for (i = 0; i < vertDOF; i++)
    {
        if (i % NDIM == 0)
        {
            vvel[i] -= v_cm_x / NVTOT;
        }
        else
        {
            vvel[i] -= v_cm_y / NVTOT;
        }
    }
    cout << "\t** Looping over time for NT = " << NT << " time steps with dt = " << dt << endl;
    for (tt = 0; tt < NT; tt++)
    {
        //after equilibrating, delete the last cell
        if (tt == 2e5 && NT >= 2e5)
        {
            cout << "Deleting particle!\n";
            deleteLastCell(smallN, largeN, NCELLS, NVTOT, cellDOF, vertDOF, szList,
                           nv, list, vvel, vpos, vF, im1, ip1, vim1, vip1, phi, a0, l0, L, largeNV, smallNV);
        }
        // VV POSITION UPDATE
        for (i = 0; i < vertDOF; i++)
        {
            // update position
            vpos[i] += dt * vvel[i] + 0.5 * dt * dt * vF[i];

            // recenter in box
            if (vpos[i] > L[i % NDIM])
                vpos[i] -= L[i % NDIM];
            else if (vpos[i] < 0)
                vpos[i] += L[i % NDIM];

            // reset forces
            vF[i] = 0;
        }

        // reset linked list
        for (gi = 0; gi < NVTOT + 1; gi++)
            list[gi] = 0;

        // reset linked list head
        for (i = 0; i < NBX; i++)
        {
            head[i] = 0;
            last[i] = 0;
        }

        // sort vertices into linked list
        for (gi = 0; gi < NVTOT; gi++)
        {
            // 1. get cell id of current particle position
            boxid = 0;
            sbtmp = 1;
            for (d = 0; d < NDIM; d++)
            {
                // add d index to 1d list
                boxid += floor(vpos[NDIM * gi + d] / lb[d]) * sbtmp;

                // increment dimensional factor
                sbtmp *= sb[d];
            }

            // 2. add to head list or link within list
            // NOTE: particle ids are labelled starting from 1, setting to 0 means end of linked list
            if (head[boxid] == 0)
            {
                head[boxid] = gi + 1;
                last[boxid] = gi + 1;
            }
            else
            {
                list[last[boxid]] = gi + 1;
                last[boxid] = gi + 1;
            }
        }

        // FORCE UPDATE

        // interaction forces (USE BOX LINKED LIST)
        pcheck = 0.0;
        U = 0.0;
        for (bi = 0; bi < NBX; bi++)
        {

            // get start of list of particles
            pi = head[bi];

            // loop over linked list
            while (pi > 0)
            {
                // real particle index
                gi = pi - 1;

                // next particle in list
                pj = list[pi];

                // loop down neighbors of pi in same cell
                while (pj > 0)
                {
                    // real index of pj
                    gj = pj - 1;

                    if (gj == ip1[gi] || gj == im1[gi])
                    {
                        pj = list[pj];
                        continue;
                    }

                    // contact distance
                    sij = vrad[gi] + vrad[gj];

                    // attraction distances
                    shellij = (1.0 + att) * sij;
                    cutij = (1.0 + 0.5 * att) * sij;

                    // particle distance
                    dx = vpos[NDIM * gj] - vpos[NDIM * gi];
                    dx -= L[0] * round(dx / L[0]);
                    if (dx < shellij)
                    {
                        dy = vpos[NDIM * gj + 1] - vpos[NDIM * gi + 1];
                        dy -= L[1] * round(dy / L[1]);
                        if (dy < shellij)
                        {
                            rij = sqrt(dx * dx + dy * dy);
                            // if two tumor cells, check for attraction
                            if (rij < shellij && rij > cutij && att > 0.0)
                            {
                                // force scale
                                ftmp = eint * ((rij / sij) - 1.0 - att) / sij;
                                fx = ftmp * (dx / rij);
                                fy = ftmp * (dy / rij);

                                // add to forces
                                vF[NDIM * gi] -= fx;
                                vF[NDIM * gi + 1] -= fy;

                                vF[NDIM * gj] += fx;
                                vF[NDIM * gj + 1] += fy;

                                // add to virial expression for pressure
                                pcheck += dx * fx + dy * fy;

                                // add to potential energy
                                U -= 0.5 * eint * pow((rij / sij) - 1.0 - att, 2.0);

                                // add to contacts
                                cindices(ci, vi, gi, NCELLS, szList);
                                cindices(cj, vj, gj, NCELLS, szList);

                                if (ci > cj)
                                    cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2]++;
                                else if (ci < cj)
                                    cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2]++;
                            }
                            else if (rij < cutij && rij > sij && att > 0.0)
                            {
                                // force scale
                                ftmp = eint * (1 - (rij / sij)) / sij;
                                fx = ftmp * (dx / rij);
                                fy = ftmp * (dy / rij);

                                // add to forces
                                vF[NDIM * gi] -= fx;
                                vF[NDIM * gi + 1] -= fy;

                                vF[NDIM * gj] += fx;
                                vF[NDIM * gj + 1] += fy;

                                // add to virial expression for pressure
                                pcheck += dx * fx + dy * fy;

                                // add to potential energy
                                U += 0.5 * eint * (pow(1.0 - (rij / sij), 2.0) - 0.5 * att * att);

                                // add to contacts
                                cindices(ci, vi, gi, NCELLS, szList);
                                cindices(cj, vj, gj, NCELLS, szList);

                                if (ci > cj)
                                    cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2]++;
                                else if (ci < cj)
                                    cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2]++;
                            }
                            // otherwise, check for two overlapping purely-repulsive vertices
                            else if (rij < sij)
                            {
                                // force scale
                                ftmp = eint * (1 - (rij / sij)) / sij;
                                fx = ftmp * (dx / rij);
                                fy = ftmp * (dy / rij);

                                // add to forces
                                vF[NDIM * gi] -= fx;
                                vF[NDIM * gi + 1] -= fy;

                                vF[NDIM * gj] += fx;
                                vF[NDIM * gj + 1] += fy;

                                // add to virial expression for pressure
                                pcheck += dx * fx + dy * fy;

                                // add to potential energy
                                U += 0.5 * eint * (pow(1.0 - (rij / sij), 2.0) - 0.5 * att * att);

                                // add to contacts
                                cindices(ci, vi, gi, NCELLS, szList);
                                cindices(cj, vj, gj, NCELLS, szList);

                                if (ci > cj)
                                    cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2]++;
                                else if (ci < cj)
                                    cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2]++;
                            }
                        }
                    }

                    // update pj
                    pj = list[pj];
                }

                // test overlaps with forward neighboring cells
                for (bj = 0; bj < NNN; bj++)
                {
                    // get first particle in neighboring cell
                    pj = head[nn[bi][bj]];

                    // loop down neighbors of pi in same cell
                    while (pj > 0)
                    {
                        // real index of pj
                        gj = pj - 1;

                        if (gj == ip1[gi] || gj == im1[gi])
                        {
                            pj = list[pj];
                            continue;
                        }
                        // contact distance
                        sij = vrad[gi] + vrad[gj];

                        // attraction distances
                        shellij = (1.0 + att) * sij;
                        cutij = (1.0 + 0.5 * att) * sij;

                        // particle distance
                        dx = vpos[NDIM * gj] - vpos[NDIM * gi];
                        dx -= L[0] * round(dx / L[0]);
                        if (dx < shellij)
                        {
                            dy = vpos[NDIM * gj + 1] - vpos[NDIM * gi + 1];
                            dy -= L[1] * round(dy / L[1]);
                            if (dy < shellij)
                            {
                                rij = sqrt(dx * dx + dy * dy);
                                // if two tumor cells, check for attraction
                                if (rij < shellij && rij > cutij && att > 0.0)
                                {
                                    // force scale
                                    ftmp = eint * ((rij / sij) - 1.0 - att) / sij;
                                    fx = ftmp * (dx / rij);
                                    fy = ftmp * (dy / rij);

                                    // add to forces
                                    vF[NDIM * gi] -= fx;
                                    vF[NDIM * gi + 1] -= fy;

                                    vF[NDIM * gj] += fx;
                                    vF[NDIM * gj + 1] += fy;

                                    // add to virial expression for pressure
                                    pcheck += dx * fx + dy * fy;

                                    // add to potential energy
                                    U -= 0.5 * eint * pow((rij / sij) - 1.0 - att, 2.0);

                                    // add to contacts
                                    cindices(ci, vi, gi, NCELLS, szList);
                                    cindices(cj, vj, gj, NCELLS, szList);

                                    if (ci > cj)
                                        cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2]++;
                                    else if (ci < cj)
                                        cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2]++;
                                }
                                else if (rij < cutij && rij > sij && att > 0.0)
                                {
                                    // force scale
                                    ftmp = eint * (1 - (rij / sij)) / sij;
                                    fx = ftmp * (dx / rij);
                                    fy = ftmp * (dy / rij);

                                    // add to forces
                                    vF[NDIM * gi] -= fx;
                                    vF[NDIM * gi + 1] -= fy;

                                    vF[NDIM * gj] += fx;
                                    vF[NDIM * gj + 1] += fy;

                                    // add to virial expression for pressure
                                    pcheck += dx * fx + dy * fy;

                                    // add to potential energy
                                    U += 0.5 * eint * (pow(1.0 - (rij / sij), 2.0) - 0.5 * att * att);

                                    // add to contacts
                                    cindices(ci, vi, gi, NCELLS, szList);
                                    cindices(cj, vj, gj, NCELLS, szList);

                                    if (ci > cj)
                                        cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2]++;
                                    else if (ci < cj)
                                        cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2]++;
                                }
                                // otherwise, check for two overlapping purely-repulsive vertices
                                else if (rij < sij)
                                {
                                    // force scale
                                    ftmp = eint * (1 - (rij / sij)) / sij;
                                    fx = ftmp * (dx / rij);
                                    fy = ftmp * (dy / rij);

                                    // add to forces
                                    vF[NDIM * gi] -= fx;
                                    vF[NDIM * gi + 1] -= fy;

                                    vF[NDIM * gj] += fx;
                                    vF[NDIM * gj + 1] += fy;

                                    // add to virial expression for pressure
                                    pcheck += dx * fx + dy * fy;

                                    // add to potential energy
                                    U += 0.5 * eint * (pow(1.0 - (rij / sij), 2.0) - 0.5 * att * att);

                                    // add to contacts
                                    cindices(ci, vi, gi, NCELLS, szList);
                                    cindices(cj, vj, gj, NCELLS, szList);

                                    if (ci > cj)
                                        cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2]++;
                                    else if (ci < cj)
                                        cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2]++;
                                }
                            }
                        }

                        // update pj
                        pj = list[pj];
                    }
                }

                // update pi index to be next
                pi = list[pi];
            }
        }

        // normalize pressure by box area (make dimensionless with extra factor of rho)
        pcheck *= (rho0 / (2.0 * L[0] * L[1]));

        // shape forces (loop over global vertex labels)
        ci = 0;
        for (gi = 0; gi < NVTOT; gi++)
        {

            // -- Area force (and get cell index ci)
            if (ci < NCELLS)
            {
                if (gi == szList[ci])
                {
                    // compute shape parameter
                    nvtmp = nv[ci];
                    a0tmp = a0[ci];
                    l0tmp = l0[ci];

                    // compute area deviation
                    atmp = area(vpos, ci, L, nv, szList);
                    da = (atmp / a0tmp) - 1.0;

                    // add to potential energy
                    U += 0.5 * pow(da, 2.0);

                    // shape force parameters (kl and kl are unitless energy ratios)
                    fa = da * (rho0 / a0tmp); // derivation from the fact that rho0^2 does not necessarily cancel a0tmp
                    fl = kl * (rho0 / l0tmp);
                    fb = kb * (rho0 / (l0tmp * l0tmp));

                    // compute cell center of mass
                    xi = vpos[NDIM * gi];
                    yi = vpos[NDIM * gi + 1];
                    cx = xi;
                    cy = yi;
                    for (vi = 1; vi < nvtmp; vi++)
                    {
                        dx = vpos.at(NDIM * (gi + vi)) - xi;
                        dx -= L[0] * round(dx / L[0]);

                        dy = vpos.at(NDIM * (gi + vi) + 1) - yi;
                        dy -= L[1] * round(dy / L[1]);

                        xi += dx;
                        yi += dy;

                        cx += xi;
                        cy += yi;
                    }
                    cx /= nvtmp;
                    cy /= nvtmp;

                    // get coordinates relative to center of mass
                    rix = vpos[NDIM * gi] - cx;
                    riy = vpos[NDIM * gi + 1] - cy;

                    // get (prior) adjacent vertices
                    rim1x = vpos[NDIM * im1[gi]] - cx;
                    rim1x -= L[0] * round(rim1x / L[0]);

                    rim1y = vpos[NDIM * im1[gi] + 1] - cy;
                    rim1y -= L[1] * round(rim1y / L[1]);

                    rim2x = vpos[NDIM * im1[im1[gi]]] - cx;
                    rim2x -= L[0] * round(rim2x / L[0]);

                    rim2y = vpos[NDIM * im1[im1[gi]] + 1] - cy;
                    rim2y -= L[1] * round(rim2y / L[1]);

                    // increment cell index
                    ci++;
                }
            }

            // get next adjacent vertices
            rip1x = vpos.at(NDIM * ip1[gi]) - cx;
            rip1x -= L[0] * round(rip1x / L[0]);

            rip1y = vpos.at(NDIM * ip1[gi] + 1) - cy;
            rip1y -= L[1] * round(rip1y / L[1]);

            // -- Area force
            vF[NDIM * gi] += 0.5 * fa * (rim1y - rip1y);
            vF[NDIM * gi + 1] += 0.5 * fa * (rip1x - rim1x);

            // -- Perimeter force

            // segment vector elements
            lim1x = rix - rim1x;
            lim1y = riy - rim1y;

            lix = rip1x - rix;
            liy = rip1y - riy;

            // segment lengths
            lim1 = sqrt(lim1x * lim1x + lim1y * lim1y);
            li = sqrt(lix * lix + liy * liy);

            // segment deviations
            dlim1 = (lim1 / l0tmp) - 1.0;
            dli = (li / l0tmp) - 1.0;

            // add to potential energy
            U += 0.5 * kl * pow(dli, 2.0);

            // add to forces
            vF[NDIM * gi] += fl * (dli * (lix / li) - dlim1 * (lim1x / lim1));
            vF[NDIM * gi + 1] += fl * (dli * (liy / li) - dlim1 * (lim1y / lim1));

            // -- Bending force
            if (kb > 0)
            {
                // segment vectors for ip2
                rip2x = vpos[NDIM * ip1[ip1[gi]]] - cx;
                rip2x -= L[0] * round(rip2x / L[0]);

                rip2y = vpos[NDIM * ip1[ip1[gi]] + 1] - cy;
                rip2y -= L[1] * round(rip2y / L[1]);

                lip1x = rip2x - rip1x;
                lip1y = rip2y - rip1y;

                lim2x = rim1x - rim2x;
                lim2y = rim1y - rim2y;

                // add to potential energy
                U += 0.5 * kb * (pow(lix - lim1x, 2.0) + pow(liy - lim1y, 2.0)) / (l0tmp * l0tmp);

                // add to force
                vF[NDIM * gi] += fb * (3.0 * (lix - lim1x) + lim2x - lip1x);
                vF[NDIM * gi + 1] += fb * (3.0 * (liy - lim1y) + lim2y - lip1y);
            }

            // update old coordinates
            rim2x = rim1x;
            rim1x = rix;
            rix = rip1x;

            rim2y = rim1y;
            rim1y = riy;
            riy = rip1y;
        }

        // VV VELOCITY UPDATE
        for (i = 0; i < vertDOF; i++)
        {
            vvel[i] += 0.5 * (vF[i] + vFold[i]) * dt;
            vFold[i] = vF[i];
        }

        // print message console, print energy to file
        if (tt % int(NSKIP / 100.) == 0)
        {
            // compute kinetic energy
            K = 0;
            for (i = 0; i < vertDOF; i++)
            {
                K += 0.5 * pow(vvel[i], 2.0);
            }

            // print energies
            cout << "\t** PRINTING ENERGIES TO FILE... " << endl;
            enout << tt << "  " << K << "  " << U << "  " << K + U << "  " << endl;
        }
        // print debug to console, print vertex positions
        if (tt % NSKIP == 0)
        {

            cout << endl
                 << endl;
            cout << "===========================================" << endl;
            cout << "			N V E 							" << endl;
            cout << "===========================================" << endl;
            cout << endl;
            cout << "	** tt = " << tt << endl;
            cout << "	** U = " << U << endl;
            cout << "	** K = " << K << endl;
            cout << "	** E = " << U + K << endl;
            cout << "	** p = " << pcheck << endl;
            cout << "	** smallN = " << smallN << endl;
            cout << "	** largeN = " << largeN << endl;
            cout << "	** NCELLS = " << NCELLS << endl;
            cout << "	** NVTOT = " << NVTOT << endl;
            cout << "	** cellDOF = " << cellDOF << endl;
            cout << "	** vertDOF = " << vertDOF << endl;
            cout << "	** dim szList = " << szList.size() << endl;
            cout << "	** dim nv = " << nv.size() << endl;
            cout << "	** dim list = " << list.size() << endl;
            cout << "	** dim vvel = " << vvel.size() << endl;
            cout << "	** dim vpos = " << vpos.size() << endl;
            cout << "	** dim vF = " << vF.size() << endl;

            // print vertex positions to check placement
            cout << "\t** PRINTING POSITIONS TO FILE... " << endl;
            printPos(jamout, vpos, vrad, a0, calA0, L, cij, nv, szList, phi, NCELLS);
        }
    }

    // close open objects
    shapeout.close();
    posout.close();
    cout << "	** FINISHED MAIN FOR cellNVE.cpp, ENDING." << endl
         << endl
         << endl;
    return 0;
}

/* 

	&&&&&&&&&&&&&&&&&&&&&&& FUNCTION DEFINITIONS &&&&&&&&&&&&&&&&&&&&&

	FUNCTIONS DEFINED

	gindex 			: returns global vertex index (gi) given cell (ci) and local vertex index (vi)
	cindex 			: returns cell index (ci) given global vertex index (gi)

	area 			: returns area of cell ci
	perimeter 		: returns perimeter of cell ci

	printPos 		: output vertex positions to .pos file for processing and visualization

    deleteLastCell  : delete last cell, reindex all relevant vectors to account for loss of a cell, its vertices, and the corresponding degrees of freedom.

	&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 

*/

// -- INDEXING

// get global vertex index gi given input cell index ci and vertex index vi
int gindex(int ci, int vi, vector<int> &szList)
{
    return szList[ci] + vi;
}

// get cell index ci and vertex index
void cindices(int &ci, int &vi, int gi, int NCELLS, vector<int> &szList)
{
    if (gi >= szList[NCELLS - 1])
    {
        ci = NCELLS - 1;
        vi = gi - szList[NCELLS - 1];
    }
    else
    {
        for (int i = 1; i < NCELLS; i++)
        {
            if (szList[i] > gi)
            {
                ci = i - 1;
                vi = gi - szList[i - 1];
                break;
            }
        }
    }
}

// -- CELL SHAPE

// get cell area
double area(vector<double> &vpos, int ci, vector<double> &L, vector<int> &nv, vector<int> &szList)
{
    // local variables
    int vi, vip1, gi, gip1;
    double dx, dy, xi, yi, xip1, yip1, areaVal = 0.0;

    // initial position: vi = 0
    gi = gindex(ci, 0, szList);
    xi = vpos[NDIM * gi];
    yi = vpos[NDIM * gi + 1];

    // loop over vertices of cell ci, get area by shoe-string method
    for (vi = 0; vi < nv.at(ci); vi++)
    {
        // next vertex
        vip1 = (vi + 1) % nv.at(ci);
        gip1 = gindex(ci, vip1, szList);

        // get positions (check minimum images)
        dx = vpos[NDIM * gip1] - xi;
        dx -= L[0] * round(dx / L[0]);
        xip1 = xi + dx;

        dy = vpos[NDIM * gip1 + 1] - yi;
        dy -= L[1] * round(dy / L[1]);
        yip1 = yi + dy;

        // increment area
        areaVal += xi * yip1 - xip1 * yi;

        // set next coordinates
        xi = xip1;
        yi = yip1;
    }
    areaVal *= 0.5;

    return abs(areaVal);
}

// get cell perimeter
double perimeter(vector<double> &vpos, int ci, vector<double> &L, vector<int> &nv, vector<int> &szList)
{
    // local variables
    int vi, vip1, gi, gip1;
    double dx, dy, xi, yi, xip1, yip1, l, perimVal = 0.0;

    // initial position: vi = 0
    gi = gindex(ci, 0, szList);
    xi = vpos[NDIM * gi];
    yi = vpos[NDIM * gi + 1];

    // loop over vertices of cell ci, get perimeter
    for (vi = 0; vi < nv.at(ci); vi++)
    {
        // next vertex
        vip1 = (vi + 1) % nv.at(ci);
        gip1 = gindex(ci, vip1, szList);

        // get positions (check minimum images)
        dx = vpos[NDIM * gip1] - xi;
        dx -= L[0] * round(dx / L[0]);
        xip1 = xi + dx;

        dy = vpos[NDIM * gip1 + 1] - yi;
        dy -= L[1] * round(dy / L[1]);
        yip1 = yi + dy;

        // compute segment length
        l = sqrt(dx * dx + dy * dy);

        // add to perimeter
        perimVal += l;

        // update coordinates
        xi = xip1;
        yi = yip1;
    }

    // return perimeter
    return perimVal;
}

// -- PRINT TO FILE

// print cell positions
void printPos(ofstream &posout, vector<double> &vpos, vector<double> &vrad, vector<double> &a0, vector<double> &calA0, vector<double> &L, vector<int> &cij, vector<int> &nv, vector<int> &szList, double phi, int NCELLS)
{
    // local variables
    int ci, cj, vi, gi, ctmp, zc, zv;
    double xi, yi, dx, dy, Lx, Ly;

    // save box sizes
    Lx = L.at(0);
    Ly = L.at(1);

    // print information starting information
    posout << setw(w) << left << "NEWFR"
           << " " << endl;
    posout << setw(w) << left << "NUMCL" << setw(w) << right << NCELLS << endl;
    posout << setw(w) << left << "PACKF" << setw(wnum) << setprecision(pnum) << right << phi << endl;

    // print box sizes
    posout << setw(w) << left << "BOXSZ";
    posout << setw(wnum) << setprecision(pnum) << right << Lx;
    posout << setw(wnum) << setprecision(pnum) << right << Ly;
    posout << endl;

    // print coordinate for rest of the cells
    for (ci = 0; ci < NCELLS; ci++)
    {

        // get cell contact data
        zc = 0;
        zv = 0;
        for (cj = 0; cj < NCELLS; cj++)
        {
            if (ci != cj)
            {
                // grab contact info from entry ci, cj
                ctmp = 0;
                if (ci > cj)
                    ctmp = cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2];
                else
                    ctmp = cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2];

                // add to contact information
                zv += ctmp;
                if (ctmp > 0)
                    zc++;
            }
        }

        // cell information
        posout << setw(w) << left << "CINFO";
        posout << setw(w) << right << nv.at(ci);
        posout << setw(w) << right << zc;
        posout << setw(w) << right << zv;
        posout << setw(wnum) << right << a0.at(ci);
        posout << setw(wnum) << right << calA0.at(ci);
        posout << endl;

        // get initial vertex positions
        gi = gindex(ci, 0, szList);
        xi = vpos.at(NDIM * gi);
        yi = vpos.at(NDIM * gi + 1);

        // place back in box center
        xi = fmod(xi, Lx);
        yi = fmod(yi, Ly);

        posout << setw(w) << left << "VINFO";
        posout << setw(w) << left << ci;
        posout << setw(w) << left << 0;

        // output initial vertex information
        posout << setw(wnum) << setprecision(pnum) << right << vrad.at(gi);
        posout << setw(wnum) << setprecision(pnum) << right << xi;
        posout << setw(wnum) << setprecision(pnum) << right << yi;
        posout << endl;

        // vertex information for next vertices
        for (vi = 1; vi < nv.at(ci); vi++)
        {
            // get global vertex index for next vertex
            gi++;

            // get next vertex positions (use MIC)
            dx = vpos.at(NDIM * gi) - xi;
            dx -= Lx * round(dx / Lx);
            xi += dx;

            dy = vpos.at(NDIM * gi + 1) - yi;
            dy -= Ly * round(dy / Ly);
            yi += dy;

            // Print indexing information
            posout << setw(w) << left << "VINFO";
            posout << setw(w) << left << ci;
            posout << setw(w) << left << vi;

            // output vertex information
            posout << setw(wnum) << setprecision(pnum) << right << vrad.at(gi);
            posout << setw(wnum) << setprecision(pnum) << right << xi;
            posout << setw(wnum) << setprecision(pnum) << right << yi;
            posout << endl;
        }
    }

    // print end frame
    posout << setw(w) << left << "ENDFR"
           << " " << endl;
}

void deleteLastCell(int &smallN, int &largeN, int &NCELLS, int &NVTOT, int &cellDOF, int &vertDOF, vector<int> &szList,
                    vector<int> &nv, vector<int> &list, vector<double> &vvel, vector<double> &vpos, vector<double> &vF, vector<int> &im1,
                    vector<int> &ip1, int &vim1, int &vip1, double phi, vector<double> a0, vector<double> l0, vector<double> L, int largeNV, int smallNV)
{

    //delete a particle, re-index all N-dependent vectors to account for this.

    //numbers of cells
    smallN = smallN;
    largeN -= 1;
    NCELLS -= 1;

    // total number of vertices
    NVTOT = smallNV * smallN + largeNV * largeN;

    // degree of freedom counts
    cellDOF = NDIM * NCELLS;
    vertDOF = NDIM * NVTOT;

    //adjust szList and nv, which keep track of global vertex indices
    szList.pop_back();
    nv.pop_back();
    //remove one largeNV's worth of indices from vectors
    for (int i = 0; i < largeNV; i++)
    {
        list.pop_back();
        for (int j = 0; j < NDIM; j++)
        {
            vvel.pop_back();
            vpos.pop_back();
            vF.pop_back();
        }
    }

    // save list of adjacent vertices
    im1 = vector<int>(NVTOT, 0);
    ip1 = vector<int>(NVTOT, 0);
    for (int ci = 0; ci < NCELLS; ci++)
    {
        for (int vi = 0; vi < nv.at(ci); vi++)
        {
            // wrap local indices
            vim1 = (vi - 1 + nv.at(ci)) % nv.at(ci);
            vip1 = (vi + 1) % nv.at(ci);

            // get global wrapped indices
            int gi = gindex(ci, vi, szList);
            im1.at(gi) = gindex(ci, vim1, szList);
            ip1.at(gi) = gindex(ci, vip1, szList);
        }
        // update packing fraction
        phi += a0[ci] + 0.25 * PI * pow(l0[ci] * del, 2.0) * (0.5 * nv[ci] - 1);
    }
    phi /= L[0] * L[1];
}