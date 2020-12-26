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
#ifndef LMP_FIX_PROPERTY_ATOM_REACTION_I
#define LMP_FIX_PROPERTY_ATOM_REACTION_I


using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- 
   unit for temperature is K.
   return a0 + a1*T + a2*T^2
-----------------------------------------------------------------------*/
inline double FixPropertyAtomReaction::property_T(int type, double temp)
{
    return temp_coeff[type][0] + temp_coeff[type][1]*temp + temp_coeff[type][2]*temp*temp;

}

/* ---------------------------------------------------------------------- 
   unit for temperature is K.
   return a0 + a1*T + a2*T^2
-----------------------------------------------------------------------*/

inline double FixPropertyAtomReaction::species_T(int type, double temp)
{
    return temp_s_coeff[type][0] + temp_s_coeff[type][1]*temp + temp_s_coeff[type][2]*temp*temp;
}

/* ---------------------------------------------------------------------- 
  Users can define their own functions here.
  User has access to variables: composition, Temp, numComp, numSolid, numLiquid, numGas.
  Variable 'composition' stores the fraction of each solid, liquid and gas species
  with dimensions nlocal by numComp.
  'Temp' varible stores the temperature for each particle with dimension of [nlocal].

  numComp -- number of composition defined in fix_reaction
  numSolid -- number of solid species
  numLiquid -- number of liquid species
  numGas -- number of gas species

  Property_s stores  the property of species with dimension of numcomp.

-----------------------------------------------------------------------*/
inline double FixPropertyAtomReaction::property_species( int iparticle,const double* property_s)
{
    double property = 0.0;
    if(strcmp(variablename,"thermalCapacity")==0)
    {
        for( int i = 0; i < numSolid+numLiquid; i++)
        {
            property += composition[iparticle][i]*property_s[i];
        }
    }
    else if (strcmp(variablename,"thermalConductivity") == 0)
    {
        double  temp = 0.0;
        for(int i = 0; i < numSolid + numLiquid; i++)
        {
            temp += composition[iparticle][i];
        }

        for (int i = 0; i < numSolid + numLiquid; i++)
        {
            property += composition[iparticle][i] / temp * property_s[i];
        }
    }
    else
        error->all(FLERR, "Property_species function is not defined for this property");
        /*
        for( int i = 0; i < numSolid+numLiquid; i++)
        {
            property += composition[iparticle][i]*property_s[i];
        }
        */

    return property;
}

#endif
