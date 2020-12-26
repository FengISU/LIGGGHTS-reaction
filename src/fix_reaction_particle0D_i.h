/* ----------------------------------------------------------------------
    This is the 
    Reaction extention for

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    Author: Fenglei Qi, Iowa State University

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:

    Fenglei QI

------------------------------------------------------------------------- */

//#include "fix_reaction_particle0D.h"
#ifndef LMP_FIX_REACTION_PARTICLE0D_I
#define LMP_FIX_REACTION_PARTICLE0D_I

using namespace LAMMPS_NS;
//for using odeint

//This is the functor called by odeint
void FixReactionParticle0D::operator()( const odeint_state_type &Y, odeint_state_type &dY, const double /* t */)
{
    double* r = new double [numReaction+1];
    double* M = new double [numComp+1];
    
    M[1] = 162.0;
    M[2] = 132.0;
    M[3] = 258.0;
    M[4] = 436.0;
    M[5] = 422.0;
    M[6] = 162.0;
    M[7] = 132.0;
    M[8] = 132.0;
    M[9] = 258.0;
    M[10] = 378.0;
    M[11] = 208.0;
    M[12] = 44.0;
    M[13] = 28.0;
    M[14] = 30.0;
    M[15] = 2.0;
    M[16] = 12.0;
    M[17] = 18.0;
    M[18] = 60.0;
    M[19] = 58.0;
    M[20] = 58.0;
    M[21] = 72.0;
    M[22] = 126.0;
    M[23] = 162.0;
    M[24] = 132.0;
    M[25] = 150.0;
    M[26] = 94.0;
    M[27] = 208.0;
    M[28] = 30.0;
    M[29] = 32.0;
    M[30] = 44.0;
    M[31] = 46.0;
    M[32] = 18.0;
    M[33] = 2.0;
    M[34] = 28.0;
    M[35] = 44.0;
    M[36] = 16.0;
    M[37] = 28.0;
    // reaction rate
    r[1] = k[1] * Y[1];
    r[2] = k[2] * Y[1];
    r[3] = k[3] * Y[6];
    r[4] = k[4] * temperature * Y[6];
    r[5] = k[5] * Y[2];
    r[6] = k[6] * Y[7];
    r[7] = k[7] * temperature * Y[7];
    r[8] = k[8] * Y[8];
    r[9] = k[9] * Y[3];
    r[10] = k[10] * Y[4];
    r[11] = k[11] * Y[5];
    r[12] = k[12] * Y[9];
    r[13] = k[13] * Y[10];
    r[14] = k[14] * temperature * Y[11];
    r[15] = k[15] * Y[11];
    r[16] = k[16] * Y[12];
    r[17] = k[17] * Y[13];
    r[18] = k[18] * Y[14];
    r[19] = k[19] * Y[15];
    r[20] = k[20] * Y[17];

    //mass balance if mass_fraction is indicated
    //solid
    dY[1] = -r[1] - r[2];
    dY[2] = -r[5];
    dY[3] = -r[9];
    dY[4] = -r[10];
    dY[5] = -r[11];
    dY[6] = r[1] -r[3] - r[4];
    dY[7] = 0.4 * r[5] - r[6] -r[7];
    dY[8] = 0.6 * r[5] - r[8];
    dY[9] = 0.35 * M[9]/M[3] *r[9] - r[12];
    dY[10] = M[10]/M[4] * r[10] + M[10]/M[5] * r[11] -r[13];
    dY[11] = M[11]/M[10] * r[13] - r[14] -r[15];
    dY[12] = 0.8 * M[12]/M[8] * r[8] - r[16];
    dY[13] = 0.8 * M[13]/M[9] * r[12] + 1.4 * M[13]/M[10] * r[13];
    dY[13] += M[13]/M[11] * r[15] - r[17];
    dY[14] = 0.8 * M[14]/M[8] * r[8] + M[14]/M[3] * r[9] + M[14]/M[9] * r[12];
    dY[14] += 0.6 * M[14]/M[10] * r[13] + 0.5 * M[14]/M[11] * r[15] - r[18];
    dY[15] = 0.75 * M[15]/M[7] * r[6] + 0.1 * M[15]/M[10] * r[13] - r[19];
    dY[16] = 6.0 * M[16]/M[1] * r[2] + 0.61 * M[16]/M[6] * r[3] + 0.675* M[16]/M[7] * r[6];
    dY[16] += M[16]/M[8] * r[8] + 5.735 * M[16]/M[3] * r[9] + 6.4 * M[16]/M[9] * r[12];
    dY[16] += 4.15 * M[16]/M[10] * r[13] + 5.5 * M[16]/M[11] * r[15];

    //liquid
    dY[17] = -1*r[20];

    //vapor
    dY[18] = 0.95 * M[18]/M[6] * r[3];
    dY[19] = 0.25 * M[19]/M[6] * r[3];
    dY[20] = 0.2 * M[20]/M[6] * r[3] + M[20]/M[4] * r[10] + 0.2 * M[20]/M[11] * r[15];
    dY[21] = 0.35 * M[21]/M[9] * r[12];
    dY[22] = 0.25 * M[22]/M[6] * r[3];
    dY[23] = r[4];
    dY[24] = r[7];
    dY[25] = 0.1 * M[25]/M[3] * r[9] + 0.3 * M[25]/M[9] * r[12];
    dY[26] = 0.08 * M[26]/M[3] * r[9] + 0.2 * M[26]/M[9] * r[12];
    dY[27] = r[14];
    dY[28] = 0.5 * M[28]/M[7] * r[6] + 0.7 * M[28]/M[8] * r[8] + 0.2 * M[28]/M[11] * r[15];
    dY[29] = 0.25 * M[29]/M[7] * r[6] + 0.25* M[29]/M[8] * r[8];
    dY[29] += M[29]/M[10] * r[13] + 0.4 * M[29]/M[11] * r[15];
    dY[30] = 0.2 * M[30]/M[6] * r[3] + 0.2 * M[30]/M[11] * r[15];
    dY[31] = 0.125 * M[31]/M[7] * r[6] + 0.125 * M[31]/M[8] * r[8];
    dY[32] = 5.0 * M[32]/M[1] * r[2] + 0.9 * M[32]/M[6] * r[3] + 0.125 * M[32]/M[7] * r[6];
    dY[32] += 0.125 * M[32]/M[8] * r[8] + M[32]/M[3] * r[9] + 0.7 * M[32]/M[9] * r[12];
    dY[32] += M[32]/M[10] * r[13] + M[32]/M[11] * r[15] + r[20];

    //gas
    dY[33] = M[33]/M[14] * r[18] + r[19];
    dY[34] = 0.23 * M[34]/M[6] * r[3] + 1.4 * M[34]/M[7] * r[6] + 0.32 * M[34]/M[3] * r[9];
    dY[34] += 0.5 * M[34]/M[11] * r[15] + r[17] + M[35]/M[14] * r[18];
    dY[35] = 0.16 * M[35]/M[6] * r[3] + 0.8 * M[35]/M[7] * r[6] + 0.2 * M[35]/M[8] * r[8];
    dY[35] += M[35]/M[5] * r[11] + r[16];
    dY[36] = 0.1 * M[36]/M[6] * r[3] + 0.625 * M[36]/M[7] * r[6] + 0.5 * M[36]/M[8] * r[8];
    dY[36] += 0.495 * M[36]/M[3] * r[9] + 0.65 * M[36]/M[9] * r[12] + 0.45 * M[36]/M[10] * r[13];
    dY[36] += 0.6 * M[36]/M[11] * r[15];
    dY[37] = 0.25 * M[37]/M[7] * r[6] + 0.25 * M[37]/M[8] * r[8] + 0.41 * M[37]/M[3] * r[9];
    dY[37] += 0.6 * M[37]/M[9] * r[12] + 0.2 * M[37]/M[10] * r[13] + 0.65 * M[37]/M[11] * r[15];

    delete [] r;
    delete [] M;

}

//This is the wrapper for calling odeint integration
void FixReactionParticle0D::ODEsolver( int i)
{
    double R = 8.31446; //unit J/mol/K
    // Defining the odeint stepper
    double abs_err = 1.0e-10;
    double rel_err = 1.0e-6;

    temperature = tempOld[i];

    for(int j = 0; j < numReaction; j++)
    {
        k[j+1] = frequencyFactor[j] * exp(-1.0*activation[j]/(R*tempOld[i]));
    }
    if(temperature < 353.15)
        k[20] = 0.0;
     //heatsource
    heatscheme g =  &FixReactionParticle0D::heat_reaction;   
    double *reaction_heat;
    double *rmass = atom->rmass;
    double  capacity = fix_capacity->vector_atom[i];
    double  high_heat_value = 0.;
    double  dt_dem = update->dt;
 
    reaction_heat = (this->*g)(k,Y,i);

    // finding the largest temperature changes due to reaction heat
    for(int j = 0; j<numReaction; j++)
    {
        if (
            reaction_heat[j+1]>high_heat_value ||
            reaction_heat[j+1]< -1.0*high_heat_value
           )
            high_heat_value = reaction_heat[j+1];
    } 

    double temp= fabs(high_heat_value*dt_dem/(rmass[i]*capacity));
    int n_ode = 1;
    double dt_ode = dt_;
    
    //reducing time step if temperature change caused by reaction heat is too large
    if ( 
         temp > 1e-4
       )
    {
        n_ode = (int)(temp*1e4)+1;
        n_ode *= nevery_;
        dt_ode /= double(n_ode);
    }
    
    double ptemp = tempOld[i];
    double total_source = 0.;
    double flux = fix_flux->vector_atom[i];
    for (int m = 0; m<n_ode; m++)
    {
        for(int j = 0; j < numReaction; j++)
        {
            k[j+1] = frequencyFactor[j] * exp(-1.0*activation[j]/(R*ptemp));
        }
        if(ptemp < 353.15)
            k[20] = 0.0;
        //updating Y at dt_ode interval
        int steps = integrate_adaptive(make_controlled< error_stepper_type >( abs_err , rel_err ), boost::ref(*this),Y,0.0, dt_ode, 0.5*dt_ode);
        //updating particle teperature
        if (n_ode == 1) break; //only one step, then get out of loop

        reaction_heat = (this->*g)(k,Y,i);
        double heat_source=0.0;
        for (int j = 0; j < numReaction;j++)
            heat_source += reaction_heat[j+1];

        total_source += heat_source*dt_ode; //summation
        
        if(fabs(capacity) > 1e-8)
            ptemp += (
                    flux +
                    heat_source
                    ) * dt_ode / (rmass[i]*capacity);
        
    }

    if(fix_heatSource)
    {
        heatSource = fix_heatSource->vector_atom;
        heatSource[i] = 0.0;
        if(n_ode == 1)
        {
            for(int j = 0; j < numReaction; j++)
            {
                heatSource[i] += reaction_heat[j+1];
            }
        }
        else
        {
            heatSource[i] = total_source/dt_;
        }
        //no interaction with other particles, no need to do forward or reverse communication
    }
}

//parameter: initial composition, reaction_rate, delta_H, rmass0; return H
//delta_H and rmass0 can be used to introduce the heat of reaction
inline double* FixReactionParticle0D::heat_reaction(const vector<double>& k, const vector<double>& Y, int i)
{
    H[1] = -k[1]*Y[1]*delta_H[1]*rmass0[i];
    H[2] = -k[2]*Y[1]*delta_H[2]*rmass0[i];
    H[3] = -k[3]*Y[6]*delta_H[3]*rmass0[i];
    H[4] = -k[4] * tempOld[i] * Y[6]*delta_H[4]*rmass0[i];
    H[5] = -k[5]*Y[2]*delta_H[5]*rmass0[i];
    H[6] = -k[6]*Y[7]*delta_H[6]*rmass0[i];
    H[7] = -k[7] * tempOld[i] * Y[7]*delta_H[7]*rmass0[i];
    H[8] = -k[8]*Y[8]*delta_H[8]*rmass0[i];
    H[9] = -k[9]*Y[3]*delta_H[9]*rmass0[i];
    H[10] = -k[10]*Y[4]*delta_H[10]*rmass0[i];
    H[11] = -k[11]*Y[5]*delta_H[11]*rmass0[i];
    H[12] = -k[12]*Y[9]*delta_H[12]*rmass0[i];
    H[13] = -k[13]*Y[10]*delta_H[13]*rmass0[i];
    H[14] = -k[14] * tempOld[i] * Y[11]*delta_H[14]*rmass0[i];
    H[15] = -k[15]*Y[11]*delta_H[15]*rmass0[i];
    H[16] = -k[16]*Y[12]*delta_H[16]*rmass0[i];
    H[17] = -k[17]*Y[13]*delta_H[17]*rmass0[i];
    H[18] = -k[18]*Y[14]*delta_H[18]*rmass0[i];
    H[19] = -k[19]*Y[15]*delta_H[19]*rmass0[i];
    H[20] = -k[20]*Y[17]*delta_H[20]*rmass0[i];
    return H;
}


#endif
