#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>

#include <AMReX_TagBox.H>
#include <AMReX_EBFArrayBox.H>

#include <AMReX_PROB_AMR_F.H>
#include <AMReX_ParmParse.H>

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>

#include "CNS_index_macros.H"
#include "CNS_parm.H"
#include "cns_prob_parm.H"

#include <cmath>
#include <algorithm>

// Our include
#include "../EBGeometry/EBGeometry.hpp"

#define PI 3.141592653589793238
#define PREATM 101325

using namespace amrex;

AMREX_GPU_DEVICE
inline
void
cns_initdata (int i, int j, int k, amrex::Array4<amrex::Real> const& state,
              amrex::GeometryData const& geomdata, Parm const& parm, ProbParm const& prob_parm)
{
    using amrex::Real;

    const Real* prob_lo = geomdata.ProbLo();
    const Real* prob_hi = geomdata.ProbHi();
    const Real* dx = geomdata.CellSize();

    // Calculate cell center coordinates
    Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
    Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
    Real z = prob_lo[2] + (k + Real(0.5)) * dx[2];
    
    // Lengths
    Real lz = prob_hi[2] - prob_lo[2];
    Real ly = prob_hi[1] - prob_lo[1];
    Real lx = prob_hi[0] - prob_lo[0];

    // Cylinder center (on the left x-face) and radius
    Real y_center = prob_parm.yCenter;
    Real z_center = prob_parm.zCenter;
    Real rad_cyl = std::sqrt((y - y_center) * (y - y_center) + (z - z_center) * (z - z_center));

    // Sphere center and radius (inside the cylinder)
    const Real sphere_center_x =  prob_lo[0] + prob_parm.xsh;
    const Real sphere_center_y = y_center; 
    const Real sphere_center_z = z_center;
    const Real sphere_radius = prob_parm.sphere_radius;

                     
    // Define inlet mach number
    Real M = std::sqrt(prob_parm.od_factor) * prob_parm.Mcj;
    Real od = std::sqrt(prob_parm.od_factor);
    Real a0 = std::sqrt(parm.eos_gamma * parm.Rsp * prob_parm.T0);
    Real u0 = prob_parm.Mobj * a0;

    // Atmospheric conditions
    Real pAtm = 101325.0;
    Real rhoAtm = 1.225;
    Real TAtm = 300;
    
    Real Pt, Tt, rhot, uxt, uyt, uzt, Yt;

    // Initialize based on position
    //if (rad_cyl <= prob_parm.rCyl && x <= prob_parm.LCyl){
        	
            // Define centers and radius for the three spherical hotspots
	    Real radius = prob_parm.sphere_radius;

	    // Coordinates for three hotspots
	    Real x_center_1 = sphere_center_x;
	    Real y_center_1 = sphere_center_y;
	    Real z_center_1 = sphere_center_z + 3.0 * radius;  // Adjust spacing

	    Real x_center_2 = sphere_center_x;
	    Real y_center_2 = sphere_center_y;
	    Real z_center_2 = sphere_center_z - 3.0 * radius;  // Adjust spacing
	    
	    Real x_center_3 = sphere_center_x;
	    Real y_center_3 = sphere_center_y + 3.0 * radius;  // Adjust spacing
	    Real z_center_3 = sphere_center_z;
	    
	    Real x_center_4 = sphere_center_x;
	    Real y_center_4 = sphere_center_y - 3.0 * radius;  // Adjust spacing
	    Real z_center_4 = sphere_center_z;

	    // Check distances to each hotspot
	    Real dist_sq_1 = std::sqrt((x - x_center_1) * (x - x_center_1) +
                           (y - y_center_1) * (y - y_center_1)
				#if (AMREX_SPACEDIM == 3)
				+ (z - z_center_1) * (z - z_center_1)
				#endif
				);

	    Real dist_sq_2 = std::sqrt((x - x_center_2) * (x - x_center_2) +
                           (y - y_center_2) * (y - y_center_2)
				#if (AMREX_SPACEDIM == 3)
				+ (z - z_center_2) * (z - z_center_2)
				#endif
				);
				
	   Real dist_sq_3 = std::sqrt((x - x_center_3) * (x - x_center_3) +
                           (y - y_center_3) * (y - y_center_3)
				#if (AMREX_SPACEDIM == 3)
				+ (z - z_center_3) * (z - z_center_3)
				#endif
				);
				
	   Real dist_sq_4 = std::sqrt((x - x_center_4) * (x - x_center_4) +
                           (y - y_center_4) * (y - y_center_4)
				#if (AMREX_SPACEDIM == 3)
				+ (z - z_center_4) * (z - z_center_4)
				#endif
				);

	    // Check if point lies in any of the hotspots
	    //if (dist_sq_1 <= radius || dist_sq_2 <= radius || dist_sq_3 <= radius || dist_sq_4 <= radius){
		     	// A ball of partially reacted hot material
		     	// Inside one of the hotspots
		//	Pt = 4.5e6;
		//	Real Temp = 2500;
		//	rhot = Pt / (parm.Rsp * Temp);
		//	uxt = Real(0.0);
		//	uyt = Real(0.0);
		//#if (AMREX_SPACEDIM == 3)
		//	uzt = Real(0.0);
		//#endif
		//	Yt = Real(0.5);  // Initialize hotspot specific state
		//     }
		if(i >= prob_parm.ilo && i <= prob_parm.ihi){
		    // Initialize with the ZND profile
		    // "ind" gives the array index (ihi is the index of shock location / von-Neumann state)
		    int ind = prob_parm.ihi - i;
		    // Obtain fuel mass fraction
		    Yt   = Real(1.0) - prob_parm.rpznd[ind]; 

		    Real c1 = Real(1.0) + (parm.eos_gamma * M * M) - (std::sqrt(Yt)*(M*M - Real(1.0)));
		    Real d0d1 = c1 / ((parm.eos_gamma + Real(1.0))*M*M);
		    // Obtain mass density
		    rhot = prob_parm.rho0 / d0d1;

		    c1 = Real(1.0) + (parm.eos_gamma*M*M) + (parm.eos_gamma*(M*M - Real(1.0))*std::sqrt(Yt));
		    Real p1p0 = c1 / (parm.eos_gamma + Real(1.0));
		    // Obtain pressure
		    Pt = prob_parm.p0 * p1p0;
		
		    Real vel = od * prob_parm.Dcj * (Real(1.0) - d0d1) + (u0*d0d1);
		    // vel gives the velocity in lab-ccordinate frame
		    uxt = vel;
		    uyt = -0.01*uxt*sin(2.0*PI*y/ly);
		    #if (AMREX_SPACEDIM == 3)
		    	uzt = -0.01*uxt*sin(2.0*PI*z/lz);
		    #endif	
		}
		else{
		     	// Inside the channel but outside ZND
			Pt = prob_parm.p0;
			Real Temp = prob_parm.T0;
            		rhot = Pt / (parm.Rsp * Temp);
			uxt = Real(0.0);
			uyt = Real(0.0);
			#if (AMREX_SPACEDIM == 3)
				    Real uzt = Real(0.0);
			#endif
			Yt = Real(1.0);  // Reactants present
		}   
        
        //}
           
    	//else {
	//	// Outside the cylinder (non-reactive region)
	//	Pt = pAtm;
	//	rhot = rhoAtm;
	//	Tt = TAtm;
	//	uxt = Real(0.0);
	//	uyt = Real(0.0);
	//	#if (AMREX_SPACEDIM == 3)
	//		    Real uzt = Real(0.0);
	//	#endif
	//	Yt = Real(0.0); // No reactants
    	//}
    
	// Set state variables
	state(i, j, k, URHO) = rhot;
	state(i, j, k, UMX) = rhot * uxt;
	state(i, j, k, UMY) = rhot * uyt;
	#if (AMREX_SPACEDIM == 3)
		state(i, j, k, UMZ) = rhot * uzt;
	#endif
	state(i, j, k, URHOY) = rhot * Yt;
	Real et = Pt / (parm.eos_gamma - Real(1.0));
	state(i, j, k, UEINT) = et;
	state(i, j, k, UEDEN) = et + Real(0.5) * (rhot * (uxt * uxt + uyt * uyt + uzt * uzt));
	state(i, j, k, UTEMP) = Pt / (parm.Rsp * rhot);
	state(i,j,k,SFOIL) = Real(0.0);
}


AMREX_GPU_HOST
AMREX_FORCE_INLINE
void
init_probparams (amrex::GeometryData const& geomdata, Parm const& parm, ProbParm& prob_parm, ProbParm& dprob_parm)
{
    using amrex::Real;

    prob_parm.rho0 = prob_parm.p0 / (parm.Rsp * prob_parm.T0);
    Real a0 = std::sqrt(parm.eos_gamma * parm.Rsp * prob_parm.T0);
    
    Real t1 = Real(0.5) * parm.q_dim * prob_parm.rho0 * (parm.eos_gamma*parm.eos_gamma-1.0) 
            / (prob_parm.p0 * parm.eos_gamma);
    prob_parm.Dcj = a0 * (std::sqrt(1.+t1) + std::sqrt(t1));
    prob_parm.Mcj = prob_parm.Dcj / a0;

    Real od = std::sqrt(prob_parm.od_factor);
    Real mach = od * prob_parm.Mcj;

    // Get the von-Neumann state variables (use normal shock relations)
    Real eta = 1.0 / (mach * mach);

    Real dvnd0 = (parm.eos_gamma+1.0) / (parm.eos_gamma - 1.0 + 2.0*eta);
    prob_parm.rhovn = prob_parm.rho0 * dvnd0;
    
    Real pvnp0 = 1. + (((2.0*parm.eos_gamma)/(parm.eos_gamma+1.0))*((1.0/eta) - 1.0));
    prob_parm.pvn = prob_parm.p0 * pvnp0;

    prob_parm.Tvn = prob_parm.pvn / (parm.Rsp * prob_parm.rhovn);

    Real u0 = prob_parm.Mobj * a0;
    prob_parm.upvn = (od * prob_parm.Dcj * (1. - (1.0/dvnd0))) + (u0 / dvnd0);
    prob_parm.uwvn = (od*prob_parm.Dcj) - prob_parm.upvn;

    // Get the lower and upper indices (global domain indices) for the ZND profile
    // Smaller index corresponds to CJ state and larger index corresponds to von-Neumann state
    // In the wave fixed frame, a detonation propagating from left to right has negative reactant and
    // product velocity
    // xsh is the location of the shock (von-Neumann state)
    const Real* prob_lo = geomdata.ProbLo();
    const Real* prob_hi = geomdata.ProbHi();
    const Real* dx      = geomdata.CellSize();
    
    Real hi_i = ((prob_parm.xsh - prob_lo[0])/dx[0]) - Real(0.5);
    prob_parm.ihi = static_cast<int>(hi_i);
    Real lo_i = ((prob_parm.xsh - 8.0*prob_parm.xd  - prob_lo[0])/dx[0]) - Real(0.5);
    prob_parm.ilo = static_cast<int>(lo_i);
    int asize = prob_parm.ihi - prob_parm.ilo + 1;
    Real tol = 1.e-5;

    prob_parm.rpznd   = (Real*)The_Arena()->alloc(sizeof(Real)*asize);
    dprob_parm.rpznd  = (Real*)The_Arena()->alloc(sizeof(Real)*asize);
    prob_parm.rpznd[0]  = Real(0.0);
    
    for(int i = 0; i < asize-1; ++i){
        // Get the density at point i
        Real c1 = Real(1.0) + (parm.eos_gamma * mach * mach) 
                - (std::sqrt(Real(1.0) - prob_parm.rpznd[i])*(mach*mach - Real(1.0)));
        Real d0d1 = c1 / ((parm.eos_gamma + Real(1.0))*(mach*mach));
        Real den = prob_parm.rho0 / d0d1;

        // Get the pressure at point i
        c1 = Real(1.0) + (parm.eos_gamma*mach*mach)
           + (parm.eos_gamma*(mach*mach - Real(1.0))*std::sqrt(Real(1.0) - prob_parm.rpznd[i]));
        Real p1p0 = c1 / (parm.eos_gamma + Real(1.0));
        Real pre = prob_parm.p0 * p1p0;

        Real Temp = pre / (parm.Rsp * den);
        Real rrate = parm.pre_exp * den * (1. - prob_parm.rpznd[i]) 
                   * std::exp(-parm.Ea_dim / (parm.Ru * Temp));

        Real vel = od*prob_parm.Dcj * (1 - d0d1) + (u0 * d0d1);
        Real dr =  (od*prob_parm.Dcj) - vel;
        if(fabs(rrate) <= tol && fabs(dr) <= tol)
            prob_parm.rpznd[i+1] = 1.0;
        else
            prob_parm.rpznd[i+1] = prob_parm.rpznd[i] + (dx[0] * rrate / dr);

        prob_parm.rpznd[i+1] = amrex::min(prob_parm.rpznd[i+1],Real(1.0));
     }

    Real d0d1 = (Real(1.0) + (parm.eos_gamma*mach*mach)) / ((parm.eos_gamma + Real(1.0))*(mach*mach));
    prob_parm.rhocj = prob_parm.rho0 / d0d1;

    Real p1p0 = (Real(1.0) + (parm.eos_gamma*mach*mach)) / (parm.eos_gamma + Real(1.0));
    prob_parm.pcj = prob_parm.p0 * p1p0;

    prob_parm.Tcj = prob_parm.pcj / (parm.Rsp * prob_parm.rhocj);
    prob_parm.upcj = od * prob_parm.Dcj * (1. - d0d1) + (d0d1 * u0);

    prob_parm.shloc = 0.0; prob_parm.flameloc = 0.0;
    dprob_parm.shloc = 0.0; dprob_parm.flameloc = 0.0; 
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
cns_probspecific_bc(
  const amrex::Real x[AMREX_SPACEDIM],
  const amrex::Real s_int[NGROW][NUM_STATE],
  amrex::Real s_ext[NUM_STATE],
  const int idir,
  const int i,
  const int j,
#if AMREX_SPACEDIM==3
  const int k,
#endif
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  Parm const& parm,
  amrex::Array4<amrex::Real> const& dat)
{  

}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
cns_tag_probspecific(int i, int j, int k, amrex::Array4<char> const& tag, 
                    amrex::Array4<amrex::Real > const& sarr, 
                    amrex::GeometryData const& geomdata, 
                    char tagval, 
                    amrex::Array4<amrex::EBCellFlag const> const& flag,
                    Parm const& parm, ProbParm const& prob_parm, amrex::Real time,
                    int level)
{
    using amrex::Real;

    if(flag(i,j,k).isRegular() && flag(i+1,j,k).isRegular() && flag(i-1,j,k).isRegular()
        && flag(i,j+1,k).isRegular() && flag(i,j-1,k).isRegular() && flag(i,j,k+1).isRegular() && flag(i,j,k-1).isRegular() && level < prob_parm.refuptolev){

        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx      = geomdata.CellSize();
        Real x = prob_lo[0] + (i+Real(0.5))*dx[0];
        Real y = prob_lo[1] + (j+Real(0.5))*dx[1];
#if AMREX_SPACEDIM==3
        Real z = prob_lo[2] + (k+Real(0.5))*dx[2];
#endif

        Real axY = 0.0, axT = 0.0, axP = 0.0;
        Real ayY = 0.0, ayT = 0.0, ayP = 0.0;
        Real azY = 0.0, azT = 0.0, azP = 0.0;
        
        if(flag(i,j,k).isConnected(1,0,0)){
            axY = amrex::Math::abs( (sarr(i+1,j,k,URHOY)/sarr(i+1,j,k,URHO)) 
                - (sarr(i,j,k,URHOY)/sarr(i,j,k,URHO)) );    
            axT = amrex::Math::abs( sarr(i+1,j,k,UTEMP) - sarr(i,j,k,UTEMP) );
            axP = (parm.eos_gamma-1.)*amrex::Math::abs( sarr(i+1,j,k,UEINT) - sarr(i,j,k,UEINT) );
        }
    
        if(flag(i,j,k).isConnected(0,1,0)){
            ayY = amrex::Math::abs( (sarr(i,j+1,k,URHOY)/sarr(i,j+1,k,URHO)) 
                - (sarr(i,j,k,URHOY)/sarr(i,j,k,URHO)) );
            ayT = amrex::Math::abs( sarr(i,j+1,k,UTEMP) - sarr(i,j,k,UTEMP) );
            ayP = (parm.eos_gamma-1.)*amrex::Math::abs( sarr(i,j+1,k,UEINT) - sarr(i,j,k,UEINT) );
        }
        
        if(flag(i,j,k).isConnected(0,0,1)){
            ayY = amrex::Math::abs( (sarr(i,j,k+1,URHOY)/sarr(i,j,k+1,URHO)) 
                - (sarr(i,j,k,URHOY)/sarr(i,j,k,URHO)) );
            ayT = amrex::Math::abs( sarr(i,j,k+1,UTEMP) - sarr(i,j,k,UTEMP) );
            ayP = (parm.eos_gamma-1.)*amrex::Math::abs( sarr(i,j,k+1,UEINT) - sarr(i,j,k,UEINT) );
        }
    
        if(flag(i,j,k).isConnected(-1,0,0)){
            axY = amrex::max(axY, 
            amrex::Math::abs((sarr(i,j,k,URHOY)/sarr(i,j,k,URHO))-(sarr(i-1,j,k,URHOY)/sarr(i-1,j,k,URHO))));

            axT = amrex::max(axT,amrex::Math::abs(sarr(i,j,k,UTEMP)-sarr(i-1,j,k,UTEMP)));

            axP = amrex::max(axP,
                (parm.eos_gamma-1.)*amrex::Math::abs(sarr(i,j,k,UEINT)-sarr(i-1,j,k,UEINT)));
        }

        if(flag(i,j,k).isConnected(0,-1,0)){
            ayY = amrex::max(ayY, 
            amrex::Math::abs((sarr(i,j,k,URHOY)/sarr(i,j,k,URHO))-(sarr(i,j-1,k,URHOY)/sarr(i,j-1,k,URHO))));

            ayT = amrex::max(ayT, amrex::Math::abs(sarr(i,j,k,UTEMP)-sarr(i,j-1,k,UTEMP)));

            ayP = amrex::max(ayP, (parm.eos_gamma-1.) 
                * amrex::Math::abs(sarr(i,j,k,UEINT)-sarr(i,j-1,k,UEINT)));
        }
        
        if(flag(i,j,k).isConnected(0,0,-1)){
            ayY = amrex::max(ayY, 
            amrex::Math::abs((sarr(i,j,k,URHOY)/sarr(i,j,k,URHO))-(sarr(i,j,k-1,URHOY)/sarr(i,j,k-1,URHO))));

            ayT = amrex::max(ayT, amrex::Math::abs(sarr(i,j,k,UTEMP)-sarr(i,j,k-1,UTEMP)));

            ayP = amrex::max(ayP, (parm.eos_gamma-1.) 
                * amrex::Math::abs(sarr(i,j,k,UEINT)-sarr(i,j,k-1,UEINT)));
        }

        Real gymag = std::sqrt(axY*axY + ayY*ayY + azY*azY);
        Real gTmag = std::sqrt(axT*axT + ayT*ayT + azT*azT);
        Real gPmag = std::sqrt(axP*axP + ayP*ayP + azP*azP);

        Real minloc = prob_parm.flameloc - prob_parm.xreflo;
        Real maxloc = prob_parm.flameloc + prob_parm.xrefhi;
    
        if((!flag(i,j,k).isCovered())){
            if(gymag >= prob_parm.deltaY 
                || gTmag >= prob_parm.deltaT
                || gPmag >= prob_parm.deltaP*PREATM
                || (x >= minloc && x <= minloc && y >= prob_parm.yreflo && y <= prob_parm.yrefhi
#if AMREX_SPACEDIM==3
                && z >= prob_parm.zreflo && z <= prob_parm.zrefhi
#endif
                ) 
                ){
                    tag(i,j,k) = tagval;
            } 
        }
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real
get_flame_location (int i, int j, int k, 
    amrex::GeometryData const& geomdata, 
    amrex::Array4<amrex::Real const> const& state,
    Parm const& parm) noexcept      
{
    using amrex::Real; 

    Real flamelocate = geomdata.ProbLo(0);    
    // FLAME LOCATION IS THE LARGEST VALUE OF x FOR
    // WHICH REACTANT MASS FRACTION DROPS BELOW 0.5
    if(state(i,j,k,URHOY) / state(i,j,k,URHO) < 0.5)
       flamelocate = geomdata.ProbLo(0) + (i+0.5)*geomdata.CellSize(0);

    return flamelocate;
} 

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real
get_shock_location (int i, int j, int k, 
    amrex::GeometryData const& geomdata, 
    amrex::Array4<Real const> const& state,
    Parm const& parm) noexcept      
{
    Real shloc = geomdata.ProbLo(0);
    // SHOCK LOCATION IS THE LARGEST VALUE OF x FOR
    // WHICH PRESSURE IS GREATER THAN 1.5 atm
    if((parm.eos_gamma-1.)*state(i,j,k,UEINT) > 1.5*PREATM)
        shloc = geomdata.ProbLo(0) + (i+0.5)*geomdata.CellSize(0);
    return shloc;
} 


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real
global_energy_release_rate (int i, int j, int k, 
    amrex::GeometryData const& geomdata, 
    amrex::Array4<Real const> const& state,
    Parm const& parm) noexcept      
{
    Real omegarhoq = state(i,j,k,URHO) * parm.pre_exp * state(i,j,k,URHOY) 
                          * std::exp(-parm.Ea_dim / (parm.Ru * state(i,j,k,UTEMP)))
                          * parm.q_dim;
    return omegarhoq;
}   

// THIS FUNCTION IS CALLED BY THE CPU AT THE END OF EVERY TIMESTEP
// (IE) TIMESTEP AT LEVEL ZERO
// THIS FUNCTION CAN BE USED TO COMPUTE THE LOCAL SHOCK AND REACTION FRONT 
// LOCATION, REACTION FRONT SPEED ETC. THESE QUANTITIES MUST BE DEFINED 
// IN prob_parm
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void
cns_probspecific_func (amrex::MultiFab& State, 
                       amrex::GeometryData const& geomdata, 
                       int write_data, Parm const& parm, 
                       ProbParm& hprob_parm, 
                       ProbParm& dprob_parm, 
                       amrex::Real time, amrex::Real dt
                       , int level)
{

}

/*!
  @brief This is an AMReX-capable version of the EBGeometry BVH accelerator.
  @tparam T The floating-point precision (float or double)
  @tparam Meta The metadata type for triangles
  @tparam BV The bounding volume type
  @tparam K The tree degree (4 for quadtree, 8 for octree)
*/
template <class T, class Meta, class BV, size_t K>
class AMReXSDF
{
public:
  /*!
    @brief Full constructor
    @param[in] a_filename STL file to load
    @param[in] a_use_bvh Whether to use BVH acceleration
  */
  AMReXSDF(const std::string& a_filename, const bool a_use_bvh)
  {
    if (a_use_bvh) {
      // Read into a BVH-accelerated implicit function
      m_sdf = EBGeometry::Parser::readIntoLinearBVH<T, Meta, BV, K>(a_filename);
    }
    else {
      // Read into a non-accelerated mesh
      m_sdf = EBGeometry::Parser::readIntoMesh<T, Meta>(a_filename);
    }
  }

  /*!
    @brief AMReX's implicit function interface - 3D coordinate version
  */
  Real operator()(AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
  {
    return m_sdf->value(EBGeometry::Vec3T<T>(AMREX_D_DECL(x, y, z)));
  }
  
  // Add a getter method to access the protected member
  std::shared_ptr<EBGeometry::ImplicitFunction<T>> getSdf() const {
    return m_sdf;
  }

  /*!
    @brief AMReX's implicit function interface - RealArray version
  */
  inline Real operator()(const RealArray& p) const noexcept
  {
    return this->operator()(AMREX_D_DECL(p[0], p[1], p[2]));
  }

protected:
  /*!
    @brief The EBGeometry SDF (signed distance function)
  */
  std::shared_ptr<EBGeometry::ImplicitFunction<T>> m_sdf;
};


AMREX_GPU_HOST
AMREX_FORCE_INLINE
void
init_eb_geometry(const amrex::Geometry& geom, const int max_coarsening_level)
{
 
}
#endif
