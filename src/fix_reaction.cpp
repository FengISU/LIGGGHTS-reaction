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
#include "fix_reaction.h"

#include <cstring>
#include "atom.h"
#include "fix_property_atom.h"
#include "fix_property_global.h"
#include "force.h"
#include "pair_gran.h"
#include "compute_pair_gran_local.h"
#include "modify.h"
#include "error.h"
#include "stdlib.h"
#include "comm.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixReaction::FixReaction(LAMMPS *lmp, int narg, char **arg):
Fix(lmp,narg,arg), 
frequencyFactor(0), 
activation(0), 
composition0(0),
temp(0),
tempOld(0),
heatReaction(0)
{
    if (!atom->sphere_flag)
        error->all(FLERR,"Fix reaction onlys work for granular particles");
    
    if (narg < 13)
        error->fix_error(FLERR,this,"not enough arguments");

    int iarg = 3;
    bool hasargs = true;
     
    while (iarg < narg && hasargs)
    {
        
        hasargs = false;
        if (strcmp(arg[iarg],"type") == 0)
        {
            if (narg < iarg + 2)
                error->fix_error(FLERR,this,"not enough arguments for 'type'");
            iarg++;
            ptype = atoi(arg[iarg++]);
            if ( !(ptype > 0) )
                error->fix_error(FLERR,this,"Num of type < 0 required");
            else if (
                        ptype >
                        atom->get_properties()->max_type()
                    )
                error->fix_error(FLERR,this,"Num of type > max_type required");

            hasargs = true;
        }
        else if (strcmp(arg[iarg],"reaction") == 0)
        {
            if (narg < iarg + 2)
                error->fix_error(FLERR,this,"not enough arguments for 'reaction'");
            iarg++;
            numReaction = atoi(arg[iarg++]);
            if ( !(numReaction > 0) )
                error->fix_error(FLERR,this,"Num of reaction > 0 required");
            hasargs = true;
        }
        else if (strcmp(arg[iarg],"composition") == 0)
        {       
            if (narg < iarg + 2)
                error->fix_error(FLERR,this,"not enough arguments for 'composition'");
            iarg++;
            fractionType = arg[iarg++];

            if (strcmp(fractionType,"default") == 0)
                fractionType = "mass_fraction";
            else if( strcmp(fractionType,"mass_fraction") == 0)
                fractionType = "mass_fraction";
            else                
                error->fix_error(FLERR,this,"Only component mass fraction type is supported");
            
            hasargs = true;
        }
        else if (strcmp(arg[iarg],"solid") == 0)
        {
            if( narg < iarg + 2)
                error->fix_error(FLERR,this,"not enough arguments for solid");
            iarg++;
            numSolid = atoi(arg[iarg++]);
            if ( !(numSolid > 0) )
                error->fix_error(FLERR,this,"Num of solids > 0 required");
            hasargs = true;
        }
        else if (strcmp(arg[iarg],"liquid") == 0)
        {
            if( narg < iarg + 2)
                error->fix_error(FLERR,this,"not enough arguments for liquid");
            iarg++;
            numLiquid = atoi(arg[iarg++]);
            if ( numLiquid < 0 )
                error->fix_error(FLERR,this,"Negative values assigned for numLiquid");
            hasargs = true;

        }
        else if (strcmp(arg[iarg],"gas") == 0)
        {
            if( narg < iarg + 2)
                error->fix_error(FLERR,this,"not enough arguments for gas");
            iarg++;
            numGas = atoi(arg[iarg++]);
            if ( numGas < 0 )
                error->fix_error(FLERR,this,"Negative values assigned for numGas");
            hasargs = true;
        }
    }

    numComp = numSolid + numLiquid + numGas;
    fix_frequencyFactor = NULL;
    fix_activation = NULL;
    fix_composition0 = NULL;
    fix_tempold = NULL;
    fix_temp = NULL;
   
    peratom_flag = 1;      
    size_peratom_cols = 0; 
    peratom_freq = 1;

    scalar_flag = 1;
    global_freq = 1;
    cpl = NULL;
    rad_mass_vary_flag = 1;

    if (comm->me == 0)
    {
        printf("********************\n");
        printf("check input for reaction\n");
        printf("numReaction=%d\n",numReaction);
        printf("numSolid=%d\n",numSolid);
        printf("numLiquid=%d\n",numLiquid);
        printf("numGas=%d\n",numGas);
        printf("numComp=%d\n",numComp);
        printf("********************\n");
    }
}

FixReaction::~FixReaction()
{
    if (frequencyFactor)
        delete [] frequencyFactor; 
    if (activation)
        delete [] activation;
    if (composition0)
        delete [] composition0;
    if(heatReaction)
        delete [] heatReaction;
}

/* ---------------------------------------------------------------------- */
void FixReaction::post_create()
{
    //register a copy of old temperature for particles
    fix_tempold = 
    static_cast<FixPropertyAtom*>
    (
        modify->find_fix_property("tempOld","property/atom","scalar",0,0,style,false)
    );

    if(!fix_tempold)
    {
        const char* fixarg[9];
        fixarg[0] = "tempOld";
        fixarg[1] = "all";
        fixarg[2] = "property/atom";
        fixarg[3] = "tempOld";
        fixarg[4] = "scalar";
        fixarg[5] = "yes";
        fixarg[6] = "yes";
        fixarg[7] = "no";
        fixarg[8] = "273.15";
        modify->add_fix(9,const_cast<char**>(fixarg));
        fix_tempold=static_cast<FixPropertyAtom*>
        (
            modify->find_fix_property("tempOld","property/atom","scalar",0,0,style)
        );
    }

}

/* ---------------------------------------------------------------------- */

void FixReaction::updatePtrs(){

  temp = fix_temp->vector_atom;
  tempOld = fix_tempold->vector_atom;
}

/* ---------------------------------------------------------------------- */

void FixReaction::init()
{
    int n = modify->n_fixes_style("multisphere");
    if (n > 0)
        error->fix_error(FLERR,this,"may not be used together with fix multisphere");
    
    if (!atom->radius_flag || !atom->rmass_flag)
        error->fix_error(FLERR,this,"must use a granular atom styple");

    //check if a fix of  this style already exists
    if(modify->n_fixes_style(style)>1)
        error->fix_error(FLERR,this,"Cannot have more than one fix of this style");

    if (!force->pair_match("gran",0))
        error->fix_error(FLERR,this,"needs a granular pair styple to be used");

    pair_gran = static_cast<PairGran*> (force->pair_match("gran",0));
    history_flag = pair_gran->is_history();

    
    if (frequencyFactor) delete [] frequencyFactor;
    if (activation) delete [] activation;
    if (composition0) delete [] composition0;
    if (heatReaction) delete [] heatReaction;

    frequencyFactor = new double[numReaction];
    activation = new double [numReaction];
    composition0 = new double [numComp];
    heatReaction = new double [numReaction];

    fix_frequencyFactor = 
    static_cast<FixPropertyGlobal*>
    (
        modify->find_fix_property("frequencyFactor","property/global","vector",numReaction,0,style)
    );

    for (int i = 0; i < numReaction; i++)
        frequencyFactor[i] = fix_frequencyFactor->compute_vector(i);

    fix_activation = 
    static_cast<FixPropertyGlobal*>
    (
        modify->find_fix_property("activation","property/global","vector",numReaction,0,style)
    );

    for (int i = 0; i < numReaction; i++)
        activation[i] = fix_activation->compute_vector(i);

    fix_composition0=
    static_cast<FixPropertyGlobal*>
    (
        modify->find_fix_property("composition0","property/global","vector",numComp,0,style)
    );

    for (int i = 0; i < numComp; i++)
        composition0[i] = fix_composition0->compute_vector(i);
    
    fix_heatofReaction = 
    static_cast<FixPropertyGlobal*>
    (
        modify->find_fix_property("heatReaction","property/global","vector",numReaction,0,style)
    );

    for (int i = 0; i < numReaction; i++)
        heatReaction[i] = fix_heatofReaction->compute_vector(i);

    fix_tempold = 
    static_cast<FixPropertyAtom*>
    (
        modify->find_fix_property("tempOld","property/atom","scalar",0,0,style)
    );

    fix_temp = 
    static_cast<FixPropertyAtom*>
    (
        modify->find_fix_property("Temp","property/atom","scalar",0,0,style)
    );

    updatePtrs();
}

/* ---------------------------------------------------------------------- */

int FixReaction::setmask()
{
    int mask = 0;
    mask |= INITIAL_INTEGRATE;
    mask |= PRE_FORCE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixReaction::initial_integrate(int vflag)
{
    updatePtrs();
    //store the particles' temperature before updating temperature
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
    {
        tempOld[i] = temp[i];
    }

    //update ghosts
    fix_tempold->do_forward_comm();
}
/* ---------------------------------------------------------------------- */

void FixReaction::pre_force(int vflag)
{
    fix_tempold->do_forward_comm();
}
/* ---------------------------------------------------------------------- */

void FixReaction::cpl_evaluate(class ComputePairGranLocal * cpl){

  char *mystyle = style;
  char *emsg = new char[100];
  sprintf(emsg, "Fix %s does not implement cpl_evaluate().\n", mystyle);
  error->all(FLERR, emsg);

}

/* ---------------------------------------------------------------------- */

void FixReaction::register_compute_pair_local(class ComputePairGranLocal *ptr){

  char *mystyle = style;
  char *emsg = new char[100];
  sprintf(emsg, "Fix %s does not implement register_compute_pair_local().\n", mystyle);
  error->all(FLERR, emsg);

}

/* ---------------------------------------------------------------------- */

void FixReaction::unregister_compute_pair_local(class ComputePairGranLocal *ptr){

  char *mystyle = style;
  char *emsg = new char[100];
  sprintf(emsg, "Fix %s does not implement unregister_compute_pair_local().\n", mystyle);
  error->all(FLERR, emsg);

}