//=====================================//
//          CorsikaConverter           //
//     --------------------------      //
//                                     //
//  For the NOvA Collaboration         //
//                                     //
//  This code  produces  a root and a  //
//  text file from CORSIKA's  binary.  //
//  The text file output is in HEPEVT  //
//  format and is  used to  integrate  //
//  CORSIKA   with  GEANT  in  NOvA's  //
//  framework.                         //
//                                     //
//     Written by                      //
//        Stefano Tognini              //
//        stognini@fnal.gov            //
//                                     //
//        HEP Group                    //
//        Federal University of Goias  //
//        Brazil                       //
//=====================================//


//_____________________uboone Parameters used_____________________


const int seed = 3;								// seed for the random number generator
// EAS are randomly distributed inside a rectangular area at the top of the atmosphere
// The area is (lowerRandomLimitX < X < upperRandomLimitX) * (lowerRandomLimitZ < Z < upperRandomLimitZ) cm^2
const double boxsize=7500/2; //extend box beyond TPC this much (cm)
const double lowerRandomLimitX = 0-boxsize;			// in cm
const double upperRandomLimitX = 256.36+boxsize;			// in cm
const double lowerRandomLimitZ = 0-boxsize;				// in cm
const double upperRandomLimitZ = 1036.8+boxsize;			// in cm

// NOvA Far Detector's time window
const double spillTime = 6.4e-3;				// in s

// Constant from the cosmic ray primaries' flux
//double timeFitParameter = 18000;			// I(E) = timeFitParameter * E^-gamma nucleons/m^2.s.sr [Cosmic Rays Review. PDG, Phys. Rev. D86, 010001. 2012]

// For the shower randomization in a spill time window
// It is important to have a negative lowerSpillRandomLimit due to the fact that every EAS needs a minimum time to reach surface.
// If lowerSpillRandomLimit = 0, then all spills will only have cosmic particles hitting the detector after a while, missing
// particles that could reach it around t = 0

const double lowerSpillRandomLimit = -3.2e6;	// in ns
const double upperSpillRandomLimit = 3.2e6;	// in ns


// Limits of the box for particles at surface. Particles outside these limits will not be considered in the GEANT4 simulation
//const double lowerBoxLimitX = 0-boxsize;			// in cm
//const double upperBoxLimitX = 256.36+boxsize;			// in cm
//const double lowerBoxLimitZ = 0-boxsize;			// in cm
//const double upperBoxLimitZ = 1036.8+boxsize;		// in cm


// Far Detector's surface height
// Value taken by using SurfaceY() in the ART Framework (from Geometry.h)
// This is the start y value
const double FDHalfHeight = 116.5;				// in cm, this is the top of the TPC
//After checking for cryo intersection with start position, project back to this height for output
const double ProjectToHeight = 1800;				// in cm, this is just above the top of LArTF

//_________________________________________________________



//_____________________NOVA Parameters used_____________________


/*
// EAS are randomly distributed inside a rectangular area at the top of the atmosphere
// The area is (lowerRandomLimitX < X < upperRandomLimitX) * (lowerRandomLimitZ < Z < upperRandomLimitZ) cm^2

const int seed = 3;								// seed for the random number generator
const double lowerRandomLimitX = -750;			// in cm
const double upperRandomLimitX = +750;			// in cm
const double lowerRandomLimitZ = 0;				// in cm
const double upperRandomLimitZ = +6000;			// in cm


// NOvA Far Detector's time window

const double spillTime = 0.0005;				// in s


// Constant from the cosmic ray primaries' flux

const double timeFitParameter = 18000;			// I(E) = timeFitParameter * E^-gamma nucleons/m^2.s.sr [Cosmic Rays Review. PDG, Phys. Rev. D86, 010001. 2012]


// For the shower randomization in a spill time window
// It is important to have a negative lowerSpillRandomLimit due to the fact that every EAS needs a minimum time to reach surface.
// If lowerSpillRandomLimit = 0, then all spills will only have cosmic particles hitting the detector after a while, missing
// particles that could reach it around t = 0

const double lowerSpillRandomLimit = -300000;	// in ns
const double upperSpillRandomLimit = +200000;	// in ns


// Limits of the box for particles at surface. Particles outside these limits will not be considered in the GEANT4 simulation
// The NOvA FD is at -750cm < X < 750cm and 0cm < Z < 6000cm, obeying the NuMI Global Coordinate System

const double lowerBoxLimitX = -5000;			// in cm
const double upperBoxLimitX = +5000;			// in cm
const double lowerBoxLimitZ = -2000;			// in cm
const double upperBoxLimitZ = +8000;			// in cm


// Far Detector's surface height
// Value taken by using SurfaceY() in the ART Framework (from Geometry.h)

const double FDHalfHeight = 994;				// in cm

//_________________________________________________________
*/







#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <strstream>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TRandom.h"
#include "TH1F.h"


using namespace std;
using namespace TMath;

const float radToDeg = TMath::RadToDeg();
const float pi = TMath::Pi();

float block[273]; 			// all data blocks have this size
float garbage;


//tpc,cryobounds: boundaries of tpc (AV),cryostat
double tpcbounds[6] = {0.0,256.35,-116.5,116.5,0.0,1036.8};
double bbox=300.0; //amount to go out in each dimension away from TPCAV = 2.5 m + 0.5 m to account for difference in cryostat boundaries from AV
double cryobounds[6] = {tpcbounds[0]-bbox,tpcbounds[1]+bbox,tpcbounds[2]-bbox,tpcbounds[3]+bbox,tpcbounds[4]-bbox,tpcbounds[5]+bbox};

//Function copied from LArSoft CRYHelper module to project positions back to wordvol
void ProjectToBoxEdge(	const double xyz[], const double dxyz[], double xyzout[]){
    // make the world box slightly smaller so that the projection to 
    // the edge avoids possible rounding errors later on with Geant4
    double fBoxDelta=1.e-5;
    double xlo = -74100 + fBoxDelta;
    double xhi = 74100 - fBoxDelta;
    double ylo = -53000 + fBoxDelta;
    double yhi = ProjectToHeight - fBoxDelta;
    double zlo = -74100 + fBoxDelta;
    double zhi = 74100 - fBoxDelta;

    // Assume we're inside the box!
    /*if(xyz[0] < xlo || xyz[0] > xhi ||
       xyz[1] < ylo || xyz[1] > yhi || 
       xyz[2] < zlo || xyz[2] > zhi)
      throw cet::exception("CRYHelper") << "Projection to edge is outside"
                                        << " bounds of world box:\n "
                                        << "\tx: " << xyz[0] << " ("
                                        << xlo << "," << xhi << ")\n"
                                        << "\ty: " << xyz[1] << " ("
                                        << ylo << "," << yhi << ")\n"
                                        << "\tz: " << xyz[2] << " ("
                                        << zlo << "," << zhi << ")";*/

    // Compute the distances to the x/y/z walls
    double dx = 99.E99;
    double dy = 99.E99;
    double dz = 99.E99;
    if      (dxyz[0] > 0.0) { dx = (xhi-xyz[0])/dxyz[0]; }
    else if (dxyz[0] < 0.0) { dx = (xlo-xyz[0])/dxyz[0]; }
    if      (dxyz[1] > 0.0) { dy = (yhi-xyz[1])/dxyz[1]; }
    else if (dxyz[1] < 0.0) { dy = (ylo-xyz[1])/dxyz[1]; }
    if      (dxyz[2] > 0.0) { dz = (zhi-xyz[2])/dxyz[2]; }
    else if (dxyz[2] < 0.0) { dz = (zlo-xyz[2])/dxyz[2]; }
    
    // Choose the shortest distance
    double d = 0.0;
    if      (dx < dy && dx < dz) d = dx;
    else if (dy < dz && dy < dx) d = dy;
    else if (dz < dx && dz < dy) d = dz;
    
    // Make the step
    for (int i = 0; i < 3; ++i) {
      xyzout[i] = xyz[i] + dxyz[i]*d;
    }
  }

bool CheckTPCIntersection( const double xyz[], const double dxyz[]){
	//taken from LArSoft code//
	//determine if this track will intersect the tpc
	//calculate the intersection point with each cryostat surface
	bool intersects_cryo = false;
	for (int bnd=0; bnd!=6; ++bnd) {
	  if (bnd<2) {
		double p2[3] = {cryobounds[bnd],  xyz[1] + (dxyz[1]/dxyz[0])*(cryobounds[bnd] - xyz[0]), xyz[2] + (dxyz[2]/dxyz[0])*(cryobounds[bnd] - xyz[0])};
		if ( p2[1] >= cryobounds[2] && p2[1] <= cryobounds[3] && p2[2] >= cryobounds[4] && p2[2] <= cryobounds[5] ) {
		 intersects_cryo = true;
		 //cout<<"Intersected bound: "<<bnd<<" at "<<p2[0]<<","<<p2[1]<<","<<p2[2]<<" ;started at "<<xyz[0]<<","<<xyz[1]<<","<<xyz[2]<<";pointed "<<dxyz[0]<<","<<dxyz[1]<<","<<dxyz[2]<<endl;
		 break;
		}
	  }
	  else if (bnd>=2 && bnd<4) {
		double p2[3] = {xyz[0] + (dxyz[0]/dxyz[1])*(cryobounds[bnd] - xyz[1]), cryobounds[bnd], xyz[2] + (dxyz[2]/dxyz[1])*(cryobounds[bnd] - xyz[1])};
		if ( p2[0] >= cryobounds[0] && p2[0] <= cryobounds[1] && p2[2] >= cryobounds[4] && p2[2] <= cryobounds[5] ) {
		  intersects_cryo = true;
		  //cout<<"Intersected bound: "<<bnd<<" at "<<p2[0]<<","<<p2[1]<<","<<p2[2]<<" ;started at "<<xyz[0]<<","<<xyz[1]<<","<<xyz[2]<<";pointed "<<dxyz[0]<<","<<dxyz[1]<<","<<dxyz[2]<<endl;
		  break;
		}
	  }
	  else if (bnd>=4) {
		double p2[3] = {xyz[0] + (dxyz[0]/dxyz[2])*(cryobounds[bnd] - xyz[2]), xyz[1] + (dxyz[1]/dxyz[2])*(cryobounds[bnd] - xyz[2]), cryobounds[bnd]};
		if ( p2[0] >= cryobounds[0] && p2[0] <= cryobounds[1] && p2[1] >= cryobounds[2] && p2[1] <= cryobounds[3] ) {
		  intersects_cryo = true;
		  //cout<<"Intersected bound: "<<bnd<<" at "<<p2[0]<<","<<p2[1]<<","<<p2[2]<<" ;started at "<<xyz[0]<<","<<xyz[1]<<","<<xyz[2]<<";pointed "<<dxyz[0]<<","<<dxyz[1]<<","<<dxyz[2]<<endl;
		  break;
		}
	  }
	}
	
	return intersects_cryo;
	//return true;
}

double wrapvar( double var, double low, double high){
  //wrap variable using modulus so that it's always between low and high
  return (var - (high - low) * floor(var/(high-low))) + low;
}
double getParticleMass(int corsikaID) {
	
	switch (corsikaID) {
		case 1:	return 0.;			// gamma
		case 2: return 0.00511;		// e+
		case 3:return 0.00511;		// e-
		case 5: return 0.10566;		// mu+
		case 6: return 0.10566;		// mu-
		case 7: return 0.13497;		// pi0
		case 8: return 0.13957;		// pi+
		case 9: return 0.13957;		// pi-
		case 10: return 0.49762;	// K0_L
		case 11: return 0.49368;	// K+
		case 12: return 0.49368;	// K-
		case 13: return 0.93956;	// n
		case 14: return 0.93827;	// p
		case 15: return 0.93827;	// pbar
		case 16: return 0.49762;	// K0_S
		case 17: return 0.54786;	// eta
		case 18: return 1.11568;	// Lambda
		case 19: return 1.18937;	// Sigma+
		case 20: return 1.19264;	// Sigma0
		case 21: return 1.19745;	// Sigma-
		case 22: return 1.31486;	// Cascade0
		case 23: return 1.32171;	// Cascade-
		case 24: return 1.67245;	// Omega-
		case 25: return 0.93956;	// nbar
		case 26: return 1.11568;	// Lambdabar
		case 27: return 1.19745;	// Sigma-bar
		case 28: return 1.19264;	// Sigma0bar
		case 29: return 1.18937;	// Sigma+bar
		case 30: return 1.31486;	// Cascade0bar
		case 31: return 1.32171;	// Cascade+bar
		case 32: return 1.67245;	// Omega+bar
			
		case 50: return 0.78265;	// omega
		case 51: return 0.77526;	// rho0
		case 52: return 0.77511;	// rho+
		case 53: return 0.77511;	// rho-
		case 54: return 1.23055;	// Delta++
		case 55: return 1.2349;		// Delta+
		case 56: return 1.2313;		// Delta0
		case 57: return 1.2349;		// Delta-
		case 58: return 1.23055;	// Delta--bar
		case 59: return 1.2349;		// Delta-bar
		case 60: return 1.2313;		// Delta0bar
		case 61: return 1.2349;		// Delta+bar
		case 62: return 0.89581;	// K*0
		case 63: return 0.89547;	// K*+
		case 64: return 0.89547;	// K*-
		case 65: return 0.89581;	// K*0bar
		case 66: return 0.;			// nu_e
		case 67: return 0.;			// nu_ebar
		case 68: return 0.;			// nu_mu
		case 69: return 0.;			// nu_mubar
		case 116: return 1.86486;	// D0
		case 117: return 1.86962;	// D+
		case 118: return 1.86962;	// D-bar
		case 119: return 1.86486;	// D0bar
		case 120: return 1.96850;	// D+_s
		case 121: return 1.96850;	// D-_sbar
		case 122: return 2.9837;	// eta_c
		case 123: return 2.00699;	// D*0
		case 124: return 2.01029;	// D*+
		case 125: return 2.01029;	// D*-bar
		case 126: return 2.00699;	// D*0bar
		case 127: return 2.1123;	// D*+_s
		case 128: return 2.1123;	// D*-_s
			
		case 130: return 3.09692;	// J/Psi
		case 131: return 1.77682;	// tau+
		case 132: return 1.77682;	// tau-
		case 133: return 0.;		// nu_tau
		case 134: return 0.;		// nu_taubar
			
		case 137: return 2.28646;	// Lambda+_c
		case 138: return 2.4678;	// Cascade+_c
		case 139: return 2.47088;	// Cascade0_c
		case 140: return 2.45398;	// Sigma++_c
		case 141: return 2.4529;	// Sigma+_c
		case 142: return 2.45374;	// Sigma0_c
		case 143: return 2.5756;	// Cascade'+_c
		case 144: return 2.5779;	// Cascade'0_c
		case 145: return 2.6952;	// Omega0_c
		case 149: return 2.28646;	// Lambda-_cbar
		case 150: return 2.4678;	// Cascade-_cbar
		case 151: return 2.47088;	// Cascade0_cbar
		case 152: return 2.45398;	// Sigma--_cbar
		case 153: return 2.4529;	// Sigma-_cbar
		case 154: return 2.45374;	// Sigma0_cbar
		case 155: return 2.5756;	// Cascade'-_cbar
		case 156: return 2.5779;	// Cascade'0_cbar
		case 157: return 2.6952;	// Omega0_cbar
		case 161: return 2.45398;	// Sigma*++_c
		case 162: return 2.4529;	// Sigma*+_c
		case 163: return 2.45374;	// Sigma*0_c
			
		case 171: return 2.45398;	// Sigma*--_cbar
		case 172: return 2.4529;	// Sigma*-_cbar
		case 173: return 2.45374;	// Sigma*0_cbar
		case 176: return 5.27958;	// B0
		case 177: return 5.27926;	// B+
		case 178: return 5.27926;	// B-bar
		case 179: return 5.27958;	// B0bar
		case 180: return 5.36677;	// B0_s
		case 181: return 5.36677;	// B0_sbar
		case 182: return 6.2745;	// B+_c
		case 183: return 6.2745;	// B-_cbar
		case 184: return 5.6194;	// Lambda0_b
		case 185: return 5.1855;	// Sigma-_b
		case 186: return 5.8113;	// Sigma+_b
		case 187: return 5.7878;	// Cascade0_b
		case 188: return 5.7911;	// Cascade-_b
		case 189: return 6.071;		// Omega-_b
		case 190: return 5.6194;	// Lambda0_bbar
		case 191: return 5.8113;	// Sigma+_bbar
		case 192: return 5.1855;	// Sigma-_bbar
		case 193: return 5.7878;	// Cascade0_bbar
		case 194: return 5.7911;	// Cascade+_bbar
		case 195: return 6.071;		// Omega+_bbar
	}
	return 0;
}





int corsikaToHepevtID(int corsikaID) {
	
	switch (corsikaID) {
			
		case 1: return 22;			// gamma
		case 2:	return  -11;		// e+
		case 3:	return 11;			// e-
		case 5: return -13;			// mu+
		case 6: return 13;			// mu-
		case 7: return 111;			// pi0
		case 8:	return 211;			// pi+
		case 9: return -211;		// pi-
		case 10: return 130;		// K0_L
		case 11: return 321;		// K+
		case 12: return -321;		// K-
		case 13: return 2112;		// n
		case 14: return 2212;		// p
		case 15: return -2212;		// pbar
		case 16: return 310;		// K0_S
		case 17: return 221;		// eta
		case 18: return 3122;		// Lambda
		case 19: return 3222;		// Sigma+
		case 20: return 3212;		// Sigma0
		case 21: return 3112;		// Sigma-
		case 22: return 3322;		// Cascade0
		case 23: return 3312;		// Cascade-
		case 24: return 3334;		// Omega-
		case 25: return -2112;		// nbar
		case 26: return -3122;		// Lambdabar
		case 27: return -3112;		// Sigma-bar
		case 28: return -3212;		// Sigma0bar
		case 29: return -3222;		// Sigma+bar
		case 30: return -3322;		// Cascade0bar
		case 31: return -3312;		// Cascade+bar
		case 32: return -3334;		// Omega+bar
			
		case 50: return 223;		// omega
		case 51: return 113;		// rho0
		case 52: return 213;		// rho+
		case 53: return -213;		// rho-
		case 54: return 2224;		// Delta++
		case 55: return 2214;		// Delta+
		case 56: return 2114;		// Delta0
		case 57: return 1114;		// Delta-
		case 58: return -2224;		// Delta--bar
		case 59: return -2214;		// Delta-bar
		case 60: return -2114;		// Delta0bar
		case 61: return -1114;		// Delta+bar
		case 62: return 10311;		// K*0
		case 63: return 10321;		// K*+
		case 64: return -10321;		// K*-
		case 65: return -10311;		// K*0bar
		case 66: return 12;			// nu_e
		case 67: return -12;		// nu_ebar
		case 68: return 14;			// nu_mu
		case 69: return -14;		// nu_mubar
			
		case 116: return 421;		// D0
		case 117: return 411;		// D+
		case 118: return -411;		// D-bar
		case 119: return -421;		// D0bar
		case 120: return 431;		// D+_s
		case 121: return -431;		// D-_sbar
		case 122: return 441;		// eta_c
		case 123: return 423;		// D*0
		case 124: return 413;		// D*+
		case 125: return -413;		// D*-bar
		case 126: return -423;		// D*0bar
		case 127: return 433;		// D*+_s
		case 128: return -433;		// D*-_s
			
		case 130: return 443;		// J/Psi
		case 131: return -15;		// tau+
		case 132: return 15;		// tau-
		case 133: return 16;		// nu_tau
		case 134: return -16;		// nu_taubar
			
		case 137: return 4122;		// Lambda+_c
		case 138: return 4232;		// Cascade+_c
		case 139: return 4132;		// Cascade0_c
		case 140: return 4222;		// Sigma++_c
		case 141: return 4212;		// Sigma+_c
		case 142: return 4112;		// Sigma0_c
		case 143: return 4322;		// Cascade'+_c
		case 144: return 4312;		// Cascade'0_c
		case 145: return 4332;		// Omega0_c
		case 149: return -4122;		// Lambda-_cbar
		case 150: return -4232;		// Cascade-_cbar
		case 151: return -4132;		// Cascade0_cbar
		case 152: return -4222;		// Sigma--_cbar
		case 153: return -4212;		// Sigma-_cbar
		case 154: return -4112;		// Sigma0_cbar
		case 155: return -4322;		// Cascade'-_cbar
		case 156: return -4312;		// Cascade'0_cbar
		case 157: return -4332;		// Omega0_cbar
		case 161: return 4224;		// Sigma*++_c
		case 162: return 1214;		// Sigma*+_c
		case 163: return 4114;		// Sigma*0_c
			
		case 171: return -4224;		// Sigma*--_cbar
		case 172: return -1214;		// Sigma*-_cbar
		case 173: return -4114;		// Sigma*0_cbar
		case 176: return 511;		// B0
		case 177: return 521;		// B+
		case 178: return -521;		// B-bar
		case 179: return -511;		// B0bar
		case 180: return 531;		// B0_s
		case 181: return -531;		// B0_sbar
		case 182: return 541;		// B+_c
		case 183: return -541;		// B-_cbar
		case 184: return 5122;		// Lambda0_b
		case 185: return 5112;		// Sigma-_b
		case 186: return 5222;		// Sigma+_b
		case 187: return 5232;		// Cascade0_b
		case 188: return 5132;		// Cascade-_b
		case 189: return 5332;		// Omega-_b
		case 190: return -5112;		// Lambda0_bbar
		case 191: return -5222;		// Sigma+_bbar
		case 192: return -5112;		// Sigma-_bbar
		case 193: return -5232;		// Cascade0_bbar
		case 194: return -5132;		// Cascade+_bbar
		case 195: return -5332;		// Omega+_bbar
	}
	
	return 0;
}






//###########################//
//                           //
// The main code starts here //
//                           //
//###########################//




int main(int argc, char *argv[]) {
	
	if (argc == 1) {
		cerr << endl;
		cerr << "Usage: ./corsikaConverter K DAT000001 DAT000002 ... DATnnnnnn" << "\a" << endl;
		cerr << "For more information, type ./corsikaConverter -h" << endl;
		cerr << endl;
		
		return 0;
	}
	
	
	
	
	
	char *helpCommand = argv[1];
	
	if ( (argc == 2) && (strcmp(helpCommand,"-h") == 0) ) {
		cout << endl;
		cout << " -------------------------------------" << endl;
		cout << " |          CorsikaConverter         |" << endl;
		cout << " |                                   |" << endl;
		cout << " | For the NOvA Collaboration        |" << endl;
		cout << " |                                   |" << endl;
		cout << " | Written by                        |" << endl;
		cout << " |       Stefano Tognini             |" << endl;
		cout << " |       stognini@fnal.gov           |" << endl;
		cout << " |                                   |" << endl;
		cout << " |       HEP Group                   |" << endl;
		cout << " |       Federal University of Goias |" << endl;
		cout << " |       Brazil                      |" << endl;
		cout << " -------------------------------------" << endl;
		cout << endl;
		cout << "This software produces 3 files for each CORSIKA's simulation binary file. ";
		cout << "The first is a ROOT file that considers events as cosmic rays air showers, ";
		cout << "the second is another ROOT file that considers one event as a 500 miroseconds spill ";
		cout << "and the third is text output file can be used to integrate CORSIKA with GEANT. ";
		cout << "It obeys the PDG particle numbering scheme and the HEPEVT format. For further details see" << endl;
		cout << "http://pdg.lbl.gov/2012/reviews/rpp2012-rev-monte-carlo-numbering.pdf" << endl;
		cout << "http://cepa.fnal.gov/psm/simulation/mcgen/lund/pythia_manual/pythia6.3/pythia6301/node39.html" << endl;
		cout << endl;
		cout << "Usage: ./corsikaConverter K DAT000001 DAT000002 ... DATnnnnnn" << endl;
		cout << endl;
		
		return 0;
	}
	
	
  //constant for computing flux of simulated primary type
	double K = atof(argv[1]);
	
	//______Reading all data files...
	
	for (int ifile = 1; ifile < argc - 1; ifile++) {
		
		int i = 1;							// block number
		bool running = true;				// end run condition
		bool validParticleSubBlock = true;	// so the code does not read a LONG block as a Particle Data Sub-Block
		int loops = 0;						// number of loops, necessary to skip the garbage data between blocks
		
		
		if (argc >= 2) {
			cerr << endl;
		}
		
		
		cerr << "For " << argv[ifile+1] << ":" << endl;
		cerr << "	Checking source file...";
		
		FILE *in;
		in = fopen(argv[ifile+1],"r");
		
		if (in == NULL) {
			cerr << " FAILED!" << endl;
			cerr << "	" << argv[ifile+1] << " does not exist" << endl;
			cerr << endl;
			
			return 0;
		}
		
		fread(&garbage,4,1,in);		// skipping the first (and useless) information
		
		
		
		//______Creating the root file
		
		
		char outputRootShowersFilename[256];
		strcpy(outputRootShowersFilename, argv[ifile+1]);
		strcat(outputRootShowersFilename,"_showers.root");
		
		
		TFile *rootShowersOutput = new TFile(outputRootShowersFilename,"RECREATE");
		
		
		
		
		
		
		//______Finding out the size of the branches
		
		int totalNumberOfShowers = 0;
		int particlesInCurrentEvent = 0;
		int numberOfParticlesPerEvent = 0;
		int currentEvent = 0;
		bool corruptedFile = true;
		
		while (running && !feof(in)) {
			
			fread(block, sizeof(block), 1, in);	// read blocks with 273 lines each
			char blockName[5];
			memcpy(blockName, block, 4);
			blockName[4] = '\0';
			
			loops++;
			
			if (loops % 21 == 0) {				// 2 garbage data floats at every 21 blocks
				fread(&garbage,4,1,in);
				fread(&garbage,4,1,in);
			}
			
			if (strcmp(blockName,"EVTE") == 0) {
				particlesInCurrentEvent = block[2] + block[3] + block[4] + block[5];
				
				if (particlesInCurrentEvent > numberOfParticlesPerEvent) {
					
					numberOfParticlesPerEvent = particlesInCurrentEvent;
				}
				currentEvent++;
			}
			
			else if (strcmp(blockName,"RUNE") == 0) {
				totalNumberOfShowers = block[2];
				running = false;		// end run condition
				corruptedFile = false;	// if there is no RUNE, the file ended unexpectedly
				
				cerr << " Done!" << endl;
			}
		}
		
		
		
		if (corruptedFile == true) {
			
			cerr << " FAILED!" << endl;
			
			if (currentEvent == 0) {
				
				cerr << "	No CORSIKA structures found" << endl;
				cerr << "	Please check file format" << endl;
				
			}
			
			else {
				
				cerr << "	Unexpected end of file on event " << currentEvent + 1 << endl;
				cerr << "	" << argv[ifile+1] << " may be corrupted" << endl;
			}
			
			cerr << endl;
			
			remove(outputRootShowersFilename);
			
			return 0;
		}
		
		
		
		
		TRandom rand(seed);
		double randomX = 0;
		double randomZ = 0;
		
		
		
		//______Now indexing the branches and creating the root file
		
		
		int runNumber;
		int NSHOW;
		double dateOfRun;
		double corsikaVersion;
		int lowEnergyModelFlag;
		int highEnergyModelFlag;
		double energySlope;
		double upperLimitOfEnergyRange;
		double lowerLimitOfEnergyRange;
		double energyCutOffHadron;
		double energyCutOffMuon;
		double energyCutOffElectron;
		double energyCutOffPhoton;
		
		
		int *ParticleID;
		ParticleID = new int[numberOfParticlesPerEvent];
		double *MuonAdditionalInfo;
		MuonAdditionalInfo = new double[numberOfParticlesPerEvent];
		double *MuonBirthAltitude;
		MuonBirthAltitude = new double[numberOfParticlesPerEvent];
		double *ParticlePx;
		ParticlePx = new double[numberOfParticlesPerEvent];
		double *ParticlePy;
		ParticlePy = new double[numberOfParticlesPerEvent];
		double *ParticlePz;
		ParticlePz = new double[numberOfParticlesPerEvent];
		double *ParticleX;
		ParticleX = new double[numberOfParticlesPerEvent];
		double *ParticleY;
		ParticleY = new double[numberOfParticlesPerEvent];
		double *ParticleZ;
		ParticleZ = new double[numberOfParticlesPerEvent];
		double *ParticleEnergy;
		ParticleEnergy = new double[numberOfParticlesPerEvent];
    double *ParticleMass;
		ParticleMass = new double[numberOfParticlesPerEvent];
		double *ParticleTheta;
		ParticleTheta = new double[numberOfParticlesPerEvent];
		double *ParticleCosTheta;
		ParticleCosTheta = new double[numberOfParticlesPerEvent];
		double *ParticlePhi;
		ParticlePhi = new double[numberOfParticlesPerEvent];
		double *ParticleTime;
		ParticleTime = new double[numberOfParticlesPerEvent];
		
		
		
		//		cout << endl;
		//		fprintf(stderr, "Valor de ParticlePx: %p",(void*)ParticlePx);
		//		cout << endl;
		
		
		int Primary;
		double N_gamma;
		double N_e;
		double N_had;
		double N_mu;
		double FirstInteraction;
		double PrimaryTheta;
		double PrimaryCosTheta;
		double PrimaryPhi;
		double EventEnergy;
		double PrimaryPx;
		double PrimaryPy;
		double PrimaryPz;
		double PrimPtot;
		
		
		
		
		TTree *input = new TTree("Input", "Input");
		TTree *particles = new TTree("Particles","Particles");
		TTree *events = new TTree("Events","Events");
		
		
		input -> Branch("RUNNR",&runNumber, "RUNNR/I");
		input -> Branch("Date",&dateOfRun, "DATE/D");
		input -> Branch("Version",&corsikaVersion, "VERSION/D");
		input -> Branch("NSHOW",&NSHOW, "NSHOW/I");
		input -> Branch("Model_High",&highEnergyModelFlag, "Model_High/I");
		input -> Branch("Model_Low",&lowEnergyModelFlag, "Model_Low/I");
		input -> Branch("ESLOPE",&energySlope, "ESLOPE/D");
		input -> Branch("ERANGE_High",&upperLimitOfEnergyRange, "ERANGE_High/D");
		input -> Branch("ERANGE_Low",&lowerLimitOfEnergyRange, "ERANGE_Low/D");
		input -> Branch("ECUTS_Hadron",&energyCutOffHadron, "ECUTS_Hadron/D");
		input -> Branch("ECUTS_Muon",&energyCutOffMuon, "ECUTS_Muon/D");
		input -> Branch("ECUTS_Electron",&energyCutOffElectron, "ECUTS_Electron/D");
		input -> Branch("ECUTS_Photon",&energyCutOffPhoton, "ECUTS_Photon/D");
		
		
		particles -> Branch("ParticleID",ParticleID, Form("ParticleID[%i]/I", numberOfParticlesPerEvent));
		particles -> Branch("MuonAdditionalInfo",MuonAdditionalInfo, Form("MuonAdditionalInfo[%i]/D", numberOfParticlesPerEvent));
		particles -> Branch("MuonBirthAltitude",MuonBirthAltitude, Form("MuonBirthAltitude[%i]/D", numberOfParticlesPerEvent));
		particles -> Branch("ParticlePx",ParticlePx, Form("ParticlePx[%i]/D", numberOfParticlesPerEvent));
		particles -> Branch("ParticlePy",ParticlePy, Form("ParticlePy[%i]/D", numberOfParticlesPerEvent));
		particles -> Branch("ParticlePz",ParticlePz, Form("ParticlePz[%i]/D", numberOfParticlesPerEvent));
		particles -> Branch("ParticleX",ParticleX, Form("ParticleX[%i]/D", numberOfParticlesPerEvent));
		particles -> Branch("ParticleZ",ParticleZ, Form("ParticleZ[%i]/D", numberOfParticlesPerEvent));
		particles -> Branch("ParticleTime",ParticleTime, Form("ParticleTime[%i]/D", numberOfParticlesPerEvent));
		particles -> Branch("ParticleEnergy",ParticleEnergy, Form("ParticleEnergy[%i]/D", numberOfParticlesPerEvent));
		particles -> Branch("ParticleTheta",ParticleTheta, Form("ParticleTheta[%i]/D", numberOfParticlesPerEvent));
		particles -> Branch("ParticleCosTheta",ParticleCosTheta, Form("ParticleCosTheta[%i]/D", numberOfParticlesPerEvent));
		particles -> Branch("ParticlePhi",ParticlePhi, Form("ParticlePhi[%i]/D", numberOfParticlesPerEvent));
		
		
		events -> Branch("Primary",&Primary,"Primary/I");
		events -> Branch("EventEnergy",&EventEnergy,"EventEnergy/D");
		events -> Branch("N_gamma",&N_gamma,"N_gamma/D");
		events -> Branch("N_e",&N_e,"N_e/D");
		events -> Branch("N_had",&N_had,"N_had/D");
		events -> Branch("N_mu",&N_mu,"N_mu/D");
		events -> Branch("FirstInteraction",&FirstInteraction,"FirstInteraction/D");
		events -> Branch("PrimaryPx",&PrimaryPx,"PrimaryPx/D");
		events -> Branch("PrimaryPy",&PrimaryPy,"PrimaryPy/D");
		events -> Branch("PrimaryPz",&PrimaryPz,"PrimaryPz/D");
		events -> Branch("PrimaryTheta",&PrimaryTheta, "PrimaryTheta/D");
		events -> Branch("PrimaryCosTheta",&PrimaryCosTheta, "PrimaryCosTheta/D");
		events -> Branch("PrimaryPhi",&PrimaryPhi,"PrimaryPhi/D");
		
		for (int zeros=0; zeros < numberOfParticlesPerEvent; zeros++)	{
			
			ParticleID[zeros]			= 0;
			MuonAdditionalInfo[zeros]	= 0;
			MuonBirthAltitude[zeros]	= 0;
			ParticlePx[zeros] 			= 0;
			ParticlePy[zeros] 			= 0;
			ParticlePz[zeros] 			= 0;
			ParticleX[zeros]			= 0;
			ParticleZ[zeros]			= 0;
			ParticleTime[zeros]			= 0;
			ParticleEnergy[zeros] 		= 0;
			ParticleTheta[zeros]		= 0;
			ParticleCosTheta[zeros]		= 0;
			ParticlePhi[zeros] 			= 0;
		}
		
		
		//______For the LONG blocks. The number 55 refers to the total number of steps in the LONG blocks.
		//______It may vary according to the longitudinal step size and the atmosphere's parameterization.
		
		//		for (int zeros=0; zeros < 55; zeros++) {
		//
		//			longBlockVariable[depthStepCount] = 0;
		//		}
		
		
		
		in = fopen(argv[ifile+1],"r");
		fread(&garbage,4,1,in);	// skipping the first (and useless) information
		
		i = 1;					// block number
		running = true;			// end run condition
		loops = 0;				// number of loops, necessary to skip the garbage data between blocks
		
		
		int depthStepCount = 0;
		int particleDataSubBlockCounter = 0;
		int currentShowerNumber = 0;
		
		int currentPercent = 0;
		int printedPercent = 0;
		
		
		
		
		while (running && !feof(in)) {
			
			fread(block,sizeof(block),1,in);	// blocks have 273 informations
			char blockName[5];
			memcpy(blockName,block,4);
			blockName[4] = '\0';
			
			loops++;
			
			if (loops % 21 == 0) {				// 2 garbage data floats between every 21 blocks
				
				fread(&garbage,4,1,in);
				fread(&garbage,4,1,in);
			}
			
			
			currentPercent = (int)(100 * currentShowerNumber) / totalNumberOfShowers;
			
			if (currentPercent > printedPercent) {
				
				printedPercent = currentPercent;
				cerr << "\r	Creating " << outputRootShowersFilename << "... " << printedPercent << "%";
				
			}
			
			
			
			
			//__________RUN HEADER
			
			if (strcmp(blockName,"RUNH")==0) {
				
				runNumber = block[1];
				dateOfRun = block[2];
				corsikaVersion = block[3];
				NSHOW = block[92];
				energySlope = block[15];
				upperLimitOfEnergyRange = block[17];
				lowerLimitOfEnergyRange = block[16];
				energyCutOffHadron = block[20];
				energyCutOffMuon = block[21];
				energyCutOffElectron = block[22];
				energyCutOffPhoton = block[23];
				
				i++;
			}
			
			
			
			//__________EVENT HEADER
			
			else if (strcmp(blockName,"EVTH")==0) {
				
				currentShowerNumber++;
				validParticleSubBlock = true;
				
				
				Primary = block[2];
				EventEnergy = block[3];
				FirstInteraction = sqrt(block[6]*block[6]);
				
				
				randomX = rand.Uniform(lowerRandomLimitX, upperRandomLimitX);
				randomZ = rand.Uniform(lowerRandomLimitZ, upperRandomLimitZ);
				
				
				
				
				// NOvA FD Global Coordinate System (same as NuMI)
				// X is horizontal
				// Y is vertical, pointing to the sky
				// Z is longitudinal, passing through the FD center.
				// (0,0,0) is in the front face of the FD, at its nominal midpoint
				
				
				PrimaryPx = -block[7];
				PrimaryPy = -block[9];
				PrimaryPz = block[8];
				
				PrimPtot = sqrt(PrimaryPx*PrimaryPx + PrimaryPy*PrimaryPy + PrimaryPz*PrimaryPz);
				PrimaryTheta = acos(PrimaryPy / PrimPtot)*radToDeg;
				PrimaryTheta = 180 - PrimaryTheta;			// This makes PrimaryTheta = Zenith
				PrimaryCosTheta = cos ( PrimaryTheta/radToDeg );
				PrimaryPhi = atan(PrimaryPz/PrimaryPx)*radToDeg;
				
				
				if (PrimaryPx < 0 && PrimaryPz > 0) PrimaryPhi = PrimaryPhi + 180;
				if (PrimaryPx > 0 && PrimaryPz < 0) PrimaryPhi = PrimaryPhi + 360;
				if (PrimaryPx < 0 && PrimaryPz < 0) PrimaryPhi = PrimaryPhi + 180;
				
				
				if (i == 2) {
					
					lowEnergyModelFlag = block[74];
					highEnergyModelFlag = block[75];
					
					input->Fill();
				}
				
				i++;
			}
			
			
			
			//__________EVENT END
			
			else if (strcmp(blockName,"EVTE")==0) {
				
				N_gamma = block[2];
				N_e = block[3];
				N_had = block[4];
				N_mu = block[5];
				
				for (int zeros = particleDataSubBlockCounter; zeros < numberOfParticlesPerEvent; zeros++) {
					
					ParticleID[zeros]			= 0;
					MuonAdditionalInfo[zeros]	= 0;
					MuonBirthAltitude[zeros]	= 0;
					ParticlePx[zeros] 			= 0;
					ParticlePy[zeros] 			= 0;
					ParticlePz[zeros] 			= 0;
					ParticleX[zeros]  			= 0;
					ParticleZ[zeros]   			= 0;
					ParticleTime[zeros]			= 0;
					ParticleEnergy[zeros]		= 0;
					ParticleTheta[zeros]		= 0;
					ParticleCosTheta[zeros]		= 0;
					ParticlePhi[zeros]  		= 0;
				}
				
				//				for (int zeros = depthStepCount; zeros < 55; zeros++) {
				//
				// 					longBlockVariable[zeros] = 0;
				// 				}
				
				
				events->Fill();
				particles->Fill();
				//				longblocks->Fill();
				
				
				
				particleDataSubBlockCounter = 0;
				
				
				depthStepCount = 0;
				
				i++;
			}
			
			
			
			
			
			//__________LONGITUDINAL SUB-BLOCK (not used now. It's here just in case)
			
			else if (strcmp(blockName,"LONG")==0) {
				/*
				 float totallongblocks = block[4]/10;
				 totallongblocks = totallongblocks - (int)totallongblocks;
				 totallongblocks = totallongblocks*10;
				 
				 for (int n=1; n<27; n++) {
				 
				 if (block[10*n+3] == 0 && block[5] > 1) {
				 
				 // nothing goes here. Seriously, don't write anything here.
				 }
				 
				 else {
				 // Ex.: longBlockVariable[depthStepCount] = block[10*n + 1];
				 depthStepCount++;
				 }
				 
				 i++;
				 
				 }
				 */
			}
			
			
			
			
			
			
			//__________END OF RUN
			
			else if (strcmp(blockName,"RUNE")==0) {
				
				cerr << "\r	Creating " << outputRootShowersFilename << "... Done!";
				
				running = false;		// end run condition
			}
			
			
			
			//__________PARTICLE DATA SUB-BLOCKS
			
			else {
				
				
				
				if ( validParticleSubBlock == true ) {
					
					
					for (int l=1; l<40; l++) { // each sub-block has up to 39 particles. If a sub-block is not fulfilled, trailing zeros are added
						
						int k = 7*(l-1);
						double id = block[k];
						
						//cout << "	" << k << " " << block[k] << " " << block[k+1] << " " << block[k+2] << " " << block[k+3] << " " << block[k+4] << " " << block[k+5] << " " << block[k+6] << endl;	// for verification purposes only
						
						
						if ( (int)(id / 1000) != 0 ) {
							
							
							
							if (id > 75000) {
								
								int idDivision = (int)(id / 1000);
								int idRemainder = (int)id % 1000;
								
								if ( idDivision == 75 && idRemainder/100 == 5 ) MuonAdditionalInfo[particleDataSubBlockCounter] = 8;	// pi+
								if ( idDivision == 75 && idRemainder/100 == 0 ) MuonAdditionalInfo[particleDataSubBlockCounter] = 11;	// K+
								if ( idDivision == 76 && idRemainder/100 == 5 ) MuonAdditionalInfo[particleDataSubBlockCounter] = 9;	// pi-
								if ( idDivision == 76 && idRemainder/100 == 0 ) MuonAdditionalInfo[particleDataSubBlockCounter] = 12;	// K-
								
								MuonBirthAltitude[particleDataSubBlockCounter] = block[k+6] / 100;
							}
							
							else if ( id < 75000 ) {
								
								id = (int)(id/1000);
								
								
								// NOvA FD Global Coordinate System (Same as NuMI)
								double P_x = -block[k+1];
								double P_y = -block[k+3];
								double P_z = block[k+2];
								
								double mass = getParticleMass(id);
								double Energy = sqrt(P_x*P_x + P_y*P_y + P_z*P_z + mass*mass);
								
								double P_t = sqrt(P_x*P_x + P_y*P_y + P_z*P_z);
								double ThetaPart = acos(P_y / P_t)*radToDeg;
								double PhiPart = atan(P_z/P_x)*radToDeg;
								
								ThetaPart = 180 - ThetaPart;	// This makes ParticleTheta = Zenith_mu
								
								if (P_x < 0 && P_z > 0) PhiPart = PhiPart + 180;
								if (P_x > 0 && P_z < 0) PhiPart = PhiPart + 360;
								if (P_x < 0 && P_z < 0) PhiPart = PhiPart + 180;
								
								
								
								
								ParticleID[particleDataSubBlockCounter]			= (int)id;
								ParticlePx[particleDataSubBlockCounter]			= P_x;
								ParticlePy[particleDataSubBlockCounter]			= P_y;
								ParticlePz[particleDataSubBlockCounter]			= P_z;
								//ParticleX[particleDataSubBlockCounter]		= -block[k+4] + randomX;	// in cm
								//ParticleZ[particleDataSubBlockCounter]		= block[k+5] + randomZ;		// in cm
								//ParticleX[particleDataSubBlockCounter]			= rand.Uniform(lowerBoxLimitX, upperBoxLimitX);		// in cm
								//ParticleZ[particleDataSubBlockCounter]			= rand.Uniform(lowerBoxLimitZ, upperBoxLimitZ);		// in cm
								//wrap particle positions around detector box using fmod
								//ParticleX[particleDataSubBlockCounter]		= wrapvar(-block[k+4],lowerRandomLimitX,upperRandomLimitX);	// in cm
								//ParticleZ[particleDataSubBlockCounter]		= wrapvar(block[k+5],lowerRandomLimitZ,upperRandomLimitZ);		// in cm
								//changed to use raw positions, wraping is done for _spills output
								ParticleX[particleDataSubBlockCounter]		= -block[k+4];
								ParticleZ[particleDataSubBlockCounter]		= block[k+5];
								ParticleTime[particleDataSubBlockCounter]		= block[k+6]; 				// in ns
								ParticleEnergy[particleDataSubBlockCounter]		= Energy;
								ParticleTheta[particleDataSubBlockCounter]		= ThetaPart;
								ParticleCosTheta[particleDataSubBlockCounter]	= cos( ThetaPart / radToDeg );
								ParticlePhi[particleDataSubBlockCounter]		= PhiPart;
								
								particleDataSubBlockCounter++;
								
							}
						}
						
						else {
							
							validParticleSubBlock = false;
							
						}
					}
				}
				
				else { /* not a Particle Data Sub-Block. Don't do anything. */ }
			}
		}
		
		cerr << endl;
		
		
		input->Write();
		particles->Write();
		events->Write();
		//		longblocks->Write();
		
		
		rootShowersOutput->Close();
		
		
		delete ParticleID;
		delete MuonAdditionalInfo;
		delete MuonBirthAltitude;
		delete ParticlePx;
		delete ParticlePy;
		delete ParticlePz;
		delete ParticleX;
		delete ParticleZ;
		delete ParticleEnergy;
		delete ParticleTheta;
		delete ParticleCosTheta;
		delete ParticlePhi;
		delete ParticleTime;
		
		currentPercent = 0;
		printedPercent = 0;
		
		
		
		
		
		
		//______Time distribution
		
		
		cerr << "	Defining the number of spills... ";
		
		
		char outputSpillsFilename[256];
		strcpy(outputSpillsFilename, argv[ifile+1]);
		strcat(outputSpillsFilename,"_spills.txt");
		
		char outputRootSpillsFilename[256];
		strcpy(outputRootSpillsFilename, argv[ifile+1]);
		strcat(outputRootSpillsFilename,"_spills.root");
		
		ofstream outputSpills;
		outputSpills.open(outputSpillsFilename, ios::out);
		
		
		TFile *rootfile;
		TTree *TreeParticles;
		TTree *TreeEvents;
		TTree *TreeInput;
		
		
		rootfile = new TFile(outputRootShowersFilename);
		TreeParticles = (TTree*)rootfile->Get("Particles");
		TreeEvents = (TTree*)rootfile->Get("Events");
		TreeInput = (TTree*)rootfile->Get("Input");
		
		
		
		TreeInput->GetEntry();
		
		
		double totalTimeWindow = 0;			// in s
		double steradians = pi;				// MAYBE IS NECESSARY TO FIX THIS.
		double showersArea = (upperRandomLimitX/100 - lowerRandomLimitX/100) * (upperRandomLimitZ/100 - lowerRandomLimitZ/100); //in m^2
		int showersPerSpill;
		int numberOfSpills;
		int currentSpill;
		int numberOfParticlesPerSpill;
		int currentShowerInSpill;
		int particleInShower;
		double *randomTimes;
		
		
		upperLimitOfEnergyRange = (TreeInput->GetLeaf("ERANGE_High"))->GetValue();
		lowerLimitOfEnergyRange = (TreeInput->GetLeaf("ERANGE_Low"))->GetValue();
		energySlope = (TreeInput->GetLeaf("ESLOPE"))->GetValue();
		
		double oneMinusGamma = 1 + energySlope;
		
    //Get primary code and compute number of nucleons
		//int PrimaryCode=(TreeEvents->GetLeaf("Primary"))->GetValue(0);
    //int NumNucleons=Primary/100;
    //if (NumNucleons==0) NumNucleons=1;
    //cout<<"Found primary type of "<<Primary<<" which has "<<NumNucleons<<" nucleons";

		double EiToOneMinusGamma = Power(lowerLimitOfEnergyRange, oneMinusGamma );
		double EfToOneMinusGamma = Power(upperLimitOfEnergyRange, oneMinusGamma );		
		
    //Compute total time window based on I(E/nucleon)=1.8e4 (E/1GeV)^(-alpha) [nucleons/m^2/s/sr/GeV] 
    //replaced 1.8e4 by K which is passed in via the command line
		totalTimeWindow = NSHOW / ( steradians * showersArea * K * (EfToOneMinusGamma - EiToOneMinusGamma) / oneMinusGamma );
    printf("\nTotal Time Window (NSHOW=%d/(steradians=%f * showersArea=%f * K=%f *(EfToOneMinusGamma=%f - EiToOneMinusGamma=%f)/oneMinusGamma=%f)))=%f\n",
      NSHOW, steradians, showersArea, K, EfToOneMinusGamma, EiToOneMinusGamma, oneMinusGamma , totalTimeWindow);		
		numberOfSpills = (int)(totalTimeWindow / spillTime);
		printf("\nNumber of spills (%f/%f)=%d\n", totalTimeWindow, spillTime, numberOfSpills);
		
		
		
		if (numberOfSpills < 2) {
			
			showersPerSpill = totalNumberOfShowers;
			numberOfSpills = 1;
			
			cerr << "\r	Now creating " << outputRootSpillsFilename << " and " << outputSpillsFilename << "... ";
			
		}
		
		
		
		
		//showersPerSpill = (int)(totalNumberOfShowers / numberOfSpills);
    showersPerSpill = (int)((float)NSHOW / (totalTimeWindow / spillTime));
		printf("\nShowers per spill (%d/%f)=%d\n", NSHOW, (totalTimeWindow / spillTime), showersPerSpill);
    printf("Wasting %d showers...\n", NSHOW-showersPerSpill*numberOfSpills);
		randomTimes = new double[showersPerSpill];
		
		
		
		
		//______Finding out the spill with the biggest number of particles from the first root file to define the branches of the new root file
		
		int biggestNumberOfParticlesPerSpill = 0;
		
		for (currentSpill = 0; currentSpill < numberOfSpills; currentSpill++) {
			
			numberOfParticlesPerSpill = 0;
			
			for (currentShowerInSpill = 0; currentShowerInSpill < showersPerSpill; currentShowerInSpill++) {
				
				TreeEvents->GetEntry(currentSpill*showersPerSpill + currentShowerInSpill);
				TreeParticles->GetEntry(currentSpill*showersPerSpill + currentShowerInSpill);
				
				int particlesInCurrentShower = (TreeEvents->GetLeaf("N_gamma"))->GetValue() + (TreeEvents->GetLeaf("N_e"))->GetValue() + (TreeEvents->GetLeaf("N_had"))->GetValue() + (TreeEvents->GetLeaf("N_mu"))->GetValue();
				
				
				for (particleInShower = 0; particleInShower < particlesInCurrentShower; particleInShower++) {
					double outputParticleX	= wrapvar((TreeParticles->GetLeaf("ParticleX"))->GetValue(particleInShower),lowerRandomLimitX,upperRandomLimitX);
					double outputParticleZ	= wrapvar((TreeParticles->GetLeaf("ParticleZ"))->GetValue(particleInShower),lowerRandomLimitZ,upperRandomLimitZ);
					
					double dxyz[3] = {(TreeParticles->GetLeaf("ParticlePx"))->GetValue(particleInShower), (TreeParticles->GetLeaf("ParticlePy"))->GetValue(particleInShower), (TreeParticles->GetLeaf("ParticlePz"))->GetValue(particleInShower)};
					double xyz[3] ={outputParticleX,FDHalfHeight,outputParticleZ};
					if(CheckTPCIntersection(xyz,dxyz)) numberOfParticlesPerSpill++;
				}
				
			}
			
			if (numberOfParticlesPerSpill > biggestNumberOfParticlesPerSpill) {
				
				biggestNumberOfParticlesPerSpill = numberOfParticlesPerSpill;
			}
		}
		
		
		
		
		TFile *rootSpillsOutput = new TFile(outputRootSpillsFilename,"RECREATE");
		
		
		input = new TTree("Input", "Input");
		particles = new TTree("Particles", "Particles");
		TTree *spills = new TTree("Spills", "Spills");
		TTree *summary = new TTree("Summary","Summary");
		
		
		
		ParticleID = new int[biggestNumberOfParticlesPerSpill];
		ParticlePx = new double[biggestNumberOfParticlesPerSpill];
		ParticlePy = new double[biggestNumberOfParticlesPerSpill];
		ParticlePz = new double[biggestNumberOfParticlesPerSpill];
		ParticleX = new double[biggestNumberOfParticlesPerSpill];
		ParticleY = new double[biggestNumberOfParticlesPerSpill];
    ParticleZ = new double[biggestNumberOfParticlesPerSpill];
		ParticleEnergy = new double[biggestNumberOfParticlesPerSpill];
    ParticleMass = new double[biggestNumberOfParticlesPerSpill];
		ParticleTheta = new double[biggestNumberOfParticlesPerSpill];
		ParticleCosTheta = new double[biggestNumberOfParticlesPerSpill];
		ParticlePhi = new double[biggestNumberOfParticlesPerSpill];
		ParticleTime = new double[biggestNumberOfParticlesPerSpill];
		
		int NumberOfParticlesInSpill;
		int NumberOfMuonsInSpill;
		int NumberOfElectronsInSpill;
		int NumberOfGammasInSpill;
		
		
		
		input -> Branch("RUNNR",&runNumber, "RUNNR/I");
		input -> Branch("Date",&dateOfRun, "DATE/D");
		input -> Branch("Version",&corsikaVersion, "VERSION/D");
		input -> Branch("NSHOW",&NSHOW, "NSHOW/I");
		input -> Branch("Model_High",&highEnergyModelFlag, "Model_High/I");
		input -> Branch("Model_Low",&lowEnergyModelFlag, "Model_Low/I");
		input -> Branch("ESLOPE",&energySlope, "ESLOPE/D");
		input -> Branch("ERANGE_High",&upperLimitOfEnergyRange, "ERANGE_High/D");
		input -> Branch("ERANGE_Low",&lowerLimitOfEnergyRange, "ERANGE_Low/D");
		input -> Branch("ECUTS_Hadron",&energyCutOffHadron, "ECUTS_Hadron/D");
		input -> Branch("ECUTS_Muon",&energyCutOffMuon, "ECUTS_Muon/D");
		input -> Branch("ECUTS_Electron",&energyCutOffElectron, "ECUTS_Electron/D");
		input -> Branch("ECUTS_Photon",&energyCutOffPhoton, "ECUTS_Photon/D");
		
		particles -> Branch("ParticleID",ParticleID, Form("ParticleID[%i]/I", biggestNumberOfParticlesPerSpill));
		particles -> Branch("ParticlePx",ParticlePx, Form("ParticlePx[%i]/D", biggestNumberOfParticlesPerSpill));
		particles -> Branch("ParticlePy",ParticlePy, Form("ParticlePy[%i]/D", biggestNumberOfParticlesPerSpill));
		particles -> Branch("ParticlePz",ParticlePz, Form("ParticlePz[%i]/D", biggestNumberOfParticlesPerSpill));
		particles -> Branch("ParticleX",ParticleX, Form("ParticleX[%i]/D", biggestNumberOfParticlesPerSpill));
    particles -> Branch("ParticleY",ParticleY, Form("ParticleY[%i]/D", biggestNumberOfParticlesPerSpill));
		particles -> Branch("ParticleZ",ParticleZ, Form("ParticleZ[%i]/D", biggestNumberOfParticlesPerSpill));
		particles -> Branch("ParticleTime",ParticleTime, Form("ParticleTime[%i]/D", biggestNumberOfParticlesPerSpill));
		particles -> Branch("ParticleEnergy",ParticleEnergy, Form("ParticleEnergy[%i]/D", biggestNumberOfParticlesPerSpill));
    particles -> Branch("ParticleMass",ParticleMass, Form("ParticleMass[%i]/D", biggestNumberOfParticlesPerSpill));
		particles -> Branch("ParticleTheta",ParticleTheta, Form("ParticleTheta[%i]/D", biggestNumberOfParticlesPerSpill));
		particles -> Branch("ParticleCosTheta",ParticleCosTheta, Form("ParticleCosTheta[%i]/D", biggestNumberOfParticlesPerSpill));
		particles -> Branch("ParticlePhi",ParticlePhi, Form("ParticlePhi[%i]/D", biggestNumberOfParticlesPerSpill));
		
		spills -> Branch("NumberOfParticlesInSpill", &NumberOfParticlesInSpill, "NumberOfParticlesInSpill/I");
		spills -> Branch("NumberOfMuonsInSpill", &NumberOfMuonsInSpill, "NumberOfMuonsInSpill/I");
		spills -> Branch("NumberOfElectronsInSpill", &NumberOfElectronsInSpill, "NumberOfElectronsInSpill/I");
		spills -> Branch("NumberOfGammasInSpill", &NumberOfGammasInSpill, "NumberOfGammasInSpill/I");
		
		
		
		summary -> Branch("NumberOfSpills", &numberOfSpills, "NumberOfSpills/D");
		summary -> Branch("ShowersPerSpill", &showersPerSpill, "ShowersPerSpill/I");
		summary -> Branch("TotalTimeWindow", &totalTimeWindow, "TotalTimeWindow/D");
		summary -> Branch("AreaAtSurface", &showersArea, "AreaAtSurface/D");
		summary -> Branch("BiggestNumberOfParticlesPerSpill", &biggestNumberOfParticlesPerSpill, "BiggestNumberOfParticlesPerSpill/I");
		
		
		input->Fill();
		summary->Fill();
		
		
		for (currentSpill = 0; currentSpill < numberOfSpills; currentSpill++) {
			
			currentPercent = (int)(100 * currentSpill) / numberOfSpills;
			
			if (currentPercent > printedPercent) {
				
				printedPercent = currentPercent;
				cerr << "\r	Now creating " << outputRootSpillsFilename << " and " << outputSpillsFilename << "... " << printedPercent << "%";
				
			}
			
			for (int zeros = 0; zeros < biggestNumberOfParticlesPerSpill; zeros++) {
				
				ParticleID[zeros]		= 0;
				ParticlePx[zeros] 		= 0;
				ParticlePy[zeros] 		= 0;
				ParticlePz[zeros] 		= 0;
				ParticleX[zeros]  		= 0;
        ParticleY[zeros]  		= 0;
				ParticleZ[zeros]   		= 0;
				ParticleTime[zeros]		= 0;
				ParticleEnergy[zeros]	= 0;
        ParticleMass[zeros]	= 0;
				ParticleTheta[zeros]	= 0;
				ParticleCosTheta[zeros]	= 0;
				ParticlePhi[zeros]  	= 0;
			}
			
			
			
			NumberOfParticlesInSpill = 0;
			numberOfParticlesPerSpill = 0;
			
			int numberOfMuonsInSpill = 0;
			int numberOfElectronsInSpill = 0;
			int numberOfGammasInSpill = 0;
			
			
			for (currentShowerInSpill = 0; currentShowerInSpill < showersPerSpill; currentShowerInSpill++) {
				
				TreeEvents->GetEntry(currentSpill*showersPerSpill + currentShowerInSpill);
				TreeParticles->GetEntry(currentSpill*showersPerSpill + currentShowerInSpill);
				
				randomTimes[currentShowerInSpill] = rand.Uniform(lowerSpillRandomLimit, upperSpillRandomLimit);
				
				
				int particlesInCurrentShower = (TreeEvents->GetLeaf("N_gamma"))->GetValue() + (TreeEvents->GetLeaf("N_e"))->GetValue() + (TreeEvents->GetLeaf("N_had"))->GetValue() + (TreeEvents->GetLeaf("N_mu"))->GetValue();
				
				for (particleInShower = 0; particleInShower < particlesInCurrentShower; particleInShower++) {
					
					int particleIDinSpill	= (TreeParticles->GetLeaf("ParticleID"))->GetValue(particleInShower);
					double outputParticleX	= (TreeParticles->GetLeaf("ParticleX"))->GetValue(particleInShower);
					double outputParticleZ	= (TreeParticles->GetLeaf("ParticleZ"))->GetValue(particleInShower);
					
					//if ( outputParticleX > lowerBoxLimitX && outputParticleX < upperBoxLimitX    &&    outputParticleZ > lowerBoxLimitZ && outputParticleZ < upperBoxLimitZ) {
						
					double dxyz[3] = {(TreeParticles->GetLeaf("ParticlePx"))->GetValue(particleInShower), (TreeParticles->GetLeaf("ParticlePy"))->GetValue(particleInShower), (TreeParticles->GetLeaf("ParticlePz"))->GetValue(particleInShower)};
					double xyz[3] ={outputParticleX,FDHalfHeight,outputParticleZ};
					if(CheckTPCIntersection(xyz,dxyz)){
						NumberOfParticlesInSpill++;
						if (particleIDinSpill == 5 || particleIDinSpill == 6) numberOfMuonsInSpill++;
						if (particleIDinSpill == 2 || particleIDinSpill == 3) numberOfElectronsInSpill++;
						if (particleIDinSpill == 1) numberOfGammasInSpill++;
					}
					
						
					//}
				}
			}
			
			
			
			// For the text file. This is the first line of each event
			outputSpills << currentSpill << " " << NumberOfParticlesInSpill << endl;
			
			numberOfParticlesPerSpill = NumberOfParticlesInSpill;
			NumberOfMuonsInSpill = numberOfMuonsInSpill;
			NumberOfGammasInSpill = numberOfGammasInSpill;
			NumberOfElectronsInSpill = numberOfElectronsInSpill;
			
			spills->Fill();
			
			
			
			
			int particleCountedInSpill = 0;
			
			for (currentShowerInSpill = 0; currentShowerInSpill < showersPerSpill; currentShowerInSpill++) {
				
				
				TreeEvents->GetEntry(currentSpill*showersPerSpill + currentShowerInSpill);
				TreeParticles->GetEntry(currentSpill*showersPerSpill + currentShowerInSpill);
				
				
				int particlesInCurrentShower = (TreeEvents->GetLeaf("N_gamma"))->GetValue() + (TreeEvents->GetLeaf("N_e"))->GetValue() + (TreeEvents->GetLeaf("N_had"))->GetValue() + (TreeEvents->GetLeaf("N_mu"))->GetValue();
				
				for (particleInShower = 0; particleInShower < particlesInCurrentShower; particleInShower++) {
					
					
					double outputParticleX		= (TreeParticles->GetLeaf("ParticleX"))->GetValue(particleInShower);
					double outputParticleZ		= (TreeParticles->GetLeaf("ParticleZ"))->GetValue(particleInShower);
					
					
					double dxyz[3] = {(TreeParticles->GetLeaf("ParticlePx"))->GetValue(particleInShower), (TreeParticles->GetLeaf("ParticlePy"))->GetValue(particleInShower), (TreeParticles->GetLeaf("ParticlePz"))->GetValue(particleInShower)};
					double xyz[3] ={outputParticleX,FDHalfHeight,outputParticleZ};
					if(CheckTPCIntersection(xyz,dxyz)){	
						
						double outputParticleID		= (TreeParticles->GetLeaf("ParticleID"))->GetValue(particleInShower);
						double outputParticlePx		= (TreeParticles->GetLeaf("ParticlePx"))->GetValue(particleInShower);
						double outputParticlePy		= (TreeParticles->GetLeaf("ParticlePy"))->GetValue(particleInShower);
						double outputParticlePz		= (TreeParticles->GetLeaf("ParticlePz"))->GetValue(particleInShower);
						double outputParticleEnergy	= (TreeParticles->GetLeaf("ParticleEnergy"))->GetValue(particleInShower);
						double outputparticleTime	= (TreeParticles->GetLeaf("ParticleTime"))->GetValue(particleInShower) + randomTimes[currentShowerInSpill];
						
						//Project back to world volume
						double proj_dxyz[3] = {-(TreeParticles->GetLeaf("ParticlePx"))->GetValue(particleInShower), -(TreeParticles->GetLeaf("ParticlePy"))->GetValue(particleInShower), -(TreeParticles->GetLeaf("ParticlePz"))->GetValue(particleInShower)}; //reverse the direction
						double proj_xyz[3];
						ProjectToBoxEdge(xyz,proj_dxyz,proj_xyz);
						
						// THE FILE IS WRITTEN IN THE FOLLOWING ORDER:
						// status code << PDG code << first mother << second mother << first daughter << second daughter << Px << Py << Pz << energy << mass << x << y << z << time of particle production
						
						outputSpills << "1 " << corsikaToHepevtID( outputParticleID ) << " 0 0 0 0 " << outputParticlePx << " " << outputParticlePy << " " << outputParticlePz << " " <<
						outputParticleEnergy << " " << getParticleMass( outputParticleID ) << " " << proj_xyz[0] << " " << proj_xyz[1] << " " << proj_xyz[2] << " " << outputparticleTime << endl;
						
						ParticleID[particleCountedInSpill]		= corsikaToHepevtID( outputParticleID );
						ParticlePx[particleCountedInSpill]		= outputParticlePx;
						ParticlePy[particleCountedInSpill]		= outputParticlePy;
						ParticlePz[particleCountedInSpill]		= outputParticlePz;
						ParticleEnergy[particleCountedInSpill]	= outputParticleEnergy;
            ParticleMass[particleCountedInSpill]	= getParticleMass( outputParticleID );
						ParticleX[particleCountedInSpill]		= proj_xyz[0];//outputParticleX;
            ParticleY[particleCountedInSpill]		= proj_xyz[1];
						ParticleZ[particleCountedInSpill]		= proj_xyz[2];//outputParticleZ;
						ParticleTime[particleCountedInSpill]	= outputparticleTime;
						ParticleTheta[particleCountedInSpill]	= (TreeParticles->GetLeaf("ParticleTheta"))->GetValue(particleInShower);
						ParticleCosTheta[particleCountedInSpill]= cos( (TreeParticles->GetLeaf("ParticleTheta"))->GetValue(particleInShower) / radToDeg );
						ParticlePhi[particleCountedInSpill]		= (TreeParticles->GetLeaf("ParticlePhi"))->GetValue(particleInShower);
						
						
						particleCountedInSpill++;
					}
				//	}
				}
			}
			
			particles->Fill();
			
		}
		
		
		
		
		input->Write();
		particles->Write();
		spills->Write();
		summary->Write();
		
		outputSpills.close();
		rootfile->Close();
		rootSpillsOutput->Close();
		
		cerr << "\r	Now creating " << outputRootSpillsFilename << " and " << outputSpillsFilename << "... Done!" << endl;
		
		
		delete randomTimes;
		delete ParticleID;
		delete ParticlePx;
		delete ParticlePy;
		delete ParticlePz;
		delete ParticleX;
		delete ParticleZ;
		delete ParticleEnergy;
		delete ParticleTheta;
		delete ParticleCosTheta;
		delete ParticlePhi;
		delete ParticleTime;
		
	}
	
	cerr << endl;
	
}



