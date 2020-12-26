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
#include "fix_reaction_particle0D.h"

#include "stdlib.h"
#include "math.h"
#include "atom.h"
#include "error.h" 
#include "compute_pair_gran_local.h"
#include "fix_property_atom.h"
#include "fix_property_global.h"
#include "fix_property_atom_reaction.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "modify.h"
#include "update.h"
#include "math_extra_liggghts.h"

using namespace LAMMPS_NS;
using namespace FixConst;

//for using odeint
using namespace std;
using namespace boost::numeric::odeint;

/* ---------------------------------------------------------------------- */

FixReactionParticle0D::FixReactionParticle0D(class LAMMPS *lmp, int narg, char **arg)
:
FixReaction(lmp, narg, arg),
ngrid(0)
{
    particleResolvedDim = 0;
    nevery_= 1;
    iarg_ = 15;
    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
        hasargs = false;
        if(strcmp(arg[iarg_++],"every") == 0)
        {
            if(iarg_+1 > narg)
                error->fix_error(FLERR,this,"not enough arguments for keyword 'every'");
            nevery_ = atoi (arg[iarg_++]);
            hasargs = true;
        }
    }

    composition = NULL;
    H = NULL;
    delta_H = NULL;
    dt_ = 0.0;
    restart_flag=false;
    numrun == 0;
    fix_capacity = NULL;
}

/* ---------------------------------------------------------------------- */

FixReactionParticle0D::~FixReactionParticle0D()
{
    //LAMMPS_NS::FixReaction::~FixReaction();
    if (H) delete [] H;
    if (delta_H) delete [] delta_H;
}

/* ---------------------------------------------------------------------- */

void FixReactionParticle0D::post_create()
{
    FixReaction::post_create();

    fix_rmass0 = 
    static_cast<FixPropertyAtom*>
    (
        modify->find_fix_property("rmass0","property/atom","scalar",0,0,this->style,false)
    );

    if(fix_rmass0 == NULL)
    {
        const char* fixarg[9];
        fixarg[0] = "rmass0";
        fixarg[1] = "all";
        fixarg[2] = "property/atom";
        fixarg[3] = "rmass0";
        fixarg[4] = "scalar";
        fixarg[5] = "yes";
        fixarg[6] = "yes";
        fixarg[7] = "no";
        fixarg[8] = "0.0";
        fix_rmass0 = static_cast<FixPropertyAtom*>
        (
            modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style)
        );
    }

    fix_radius0 = 
    static_cast<FixPropertyAtom*>
    (
        modify->find_fix_property("radius0","property/atom","scalar",0,0,this->style,false)
    );

    if(fix_radius0 == NULL)
    {
        const char* fixarg[9];
        fixarg[0] = "radius0";
        fixarg[1] = "all";
        fixarg[2] = "property/atom";
        fixarg[3] = "radius0";
        fixarg[4] = "scalar";
        fixarg[5] = "yes";
        fixarg[6] = "yes";
        fixarg[7] = "no";
        fixarg[8] = "0.0";
        fix_radius0 = static_cast<FixPropertyAtom*>
        (
            modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style)
        );
    }

    fix_density0 = 
    static_cast<FixPropertyAtom*>
    (
        modify->find_fix_property("density0","property/atom","scalar",0,0,this->style,false)
    );

    if(fix_density0 == NULL)
    {
        const char* fixarg[9];
        fixarg[0] = "density0";
        fixarg[1] = "all";
        fixarg[2] = "property/atom";
        fixarg[3] = "density0";
        fixarg[4] = "scalar";
        fixarg[5] = "yes";
        fixarg[6] = "yes";
        fixarg[7] = "no";
        fixarg[8] = "0.0";
        fix_density0 = static_cast<FixPropertyAtom*>
        (
            modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style)
        );
    }

    //register composition fraction
    if (particleResolvedDim == 0)
    {
        fix_composition = 
        static_cast<FixPropertyAtom*> 
        (
            modify->find_fix_property("compositionFraction","property/atom","vector",numComp,0,this->style,false)
        );
        if (!fix_composition)
        {
            char** fixarg = new char* [8+numComp]; //don't forget to release the space
            fixarg[0] = "compositionFraction";
            fixarg[1] = "all";
            fixarg[2] = "property/atom";
            fixarg[3] = "compositionFraction";
            fixarg[4] = "vector";
            fixarg[5] = "yes";
            fixarg[6] = "yes";
            fixarg[7] = "no";

            int iter;
            for (iter=0;iter<numComp;iter++)
            {
                fixarg[8+iter] = "0.0";
            }
            fix_composition = static_cast<FixPropertyAtom*>
            (
                modify->add_fix_property_atom(8+numComp,const_cast<char**>(fixarg),style)
            );
            
            delete [] fixarg;
        }
    }

}

/* ---------------------------------------------------------------------- */

void FixReactionParticle0D::pre_delete(bool unfixflag)
{
    //inform cpl that this fix is deleted;
    if(cpl && unfixflag) cpl->reference_deleted();
}

/* ---------------------------------------------------------------------- */

int FixReactionParticle0D::setmask()
{
  int mask = FixReaction::setmask();
  mask |= INITIAL_INTEGRATE;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixReactionParticle0D::init()
{
    //initialize base class
    FixReaction::init();
    //initialization of rmass0, density0, and radius0;
    fix_rmass0 = static_cast<FixPropertyAtom*>
        (
            modify->find_fix_property("rmass0","property/atom","scalar",0,0,this->style)
        );

    rmass0 = fix_rmass0->vector_atom;

    fix_density0 = 
    static_cast<FixPropertyAtom*>
        (
            modify->find_fix_property("density0","property/atom","scalar",0,0,this->style)
        );

    density0 = fix_density0->vector_atom;

    fix_radius0 = static_cast<FixPropertyAtom*>
        (
            modify->find_fix_property("radius0","property/atom","scalar",0,0,this->style)
        );

    radius0 = fix_radius0->vector_atom;

    fix_composition = 
    static_cast<FixPropertyAtom*> 
        (
            modify->find_fix_property("compositionFraction","property/atom","vector",numComp,0,this->style)
        );

    composition = fix_composition->array_atom; 
    fix_capacity =  static_cast<FixPropertyAtomReaction*>(modify->find_fix_property("thermalCapacity","property/atom/reaction","scalar",1,0,style));
    fix_flux=static_cast<FixPropertyAtom*>(modify->find_fix_property("heatFlux","property/atom","scalar",0,0,style));
    //create_attribute is set to 0, we don't allow default initialization of these properties
    fix_rmass0->create_attribute = 0;
    fix_density0->create_attribute = 0;
    fix_radius0->create_attribute = 0;
    fix_composition->create_attribute = 0;
    //consider restart, does restart process save first_update? no, --> program goes through init() again
    //only initialize in the first run
    if((!fix_rmass0->recent_restart) && (numrun == 0) )
    {

        int nlocal = atom->nlocal;
        double* rmass = atom->rmass;
        double* radius = atom->radius;
        double* density  = atom->density;
        int *mask = atom->mask;
        int *type = atom->type;

        for(int iparticle = 0; iparticle < nlocal; iparticle++)
        {
            if(type[iparticle] == ptype)
            {
                rmass0[iparticle] = rmass[iparticle];
                radius0[iparticle] = radius[iparticle];
                density0[iparticle] = density[iparticle];
            }
        }

        //Comment: Here we are going to initialize the composition of existing particles
        //in the system. How about the new particles that are to be fed in through
        // comment: this problem is delt with in ParticleToInsert class in the function insert()

        for (int i = 0; i < nlocal; i++)
        {
            if(type[i] == ptype)
            {
                for(int j = 0; j < numComp; j++)
                {
                    composition[i][j] = composition0[j];                     
                }
            }
        
        }

    }

    numrun++;
    //forward info
    fix_rmass0->do_forward_comm();
    fix_radius0->do_forward_comm();
    fix_density0->do_forward_comm();
    fix_composition->do_forward_comm();


    if(H) delete [] H;
    H = new double [numReaction+1];
    if(delta_H) delete [] delta_H;
    delta_H = new double [numReaction+1];
    for (int i = 0; i < numReaction; i++)
        delta_H[i+1] = heatReaction[i];
    //reaction time interval
    double dt = update->dt;
    dt_ = dt * nevery_;

    //odeint wrapper variable
    Y.resize(numComp+1);
    Yold.resize(numComp+1);
    k.resize(numReaction+1);

    fix_heatSource = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatSource","property/atom","scalar",0,0,style));

    if(fix_heatSource)
        heatSource = fix_heatSource->vector_atom;
    
    if(fix_composition->recent_restart)
        restart_flag = true;

    #if 0
    printf("********************\n");
    printf("Print all varibales after init function \n");
    printf("input for fix reaction:\n");
    printf("numComp = %d, numSolid=%d, numLiquid=%d, numGas=%d,numReaction=%d\n",numComp,numSolid,numLiquid,numGas,numReaction);

    printf("Reaction parameters \n");   
    for (int i =0; i < numReaction; i++)
    {
        printf("A=%.5f, E=%.5f\n",frequencyFactor[i],activation[i]);
    }
    
    printf("Initial composition: ");
    for (int i =0; i < numComp; i++)
    {
        printf("%.4f, ", composition0[i]);
    }
    printf("\n");

    int nlocal = atom->nlocal;
    
    printf("particle initial temperature: \n");
    for (int i =0; i < nlocal; i++)
    {
        printf("%d: %.5f \n",atom->tag[i],tempOld[i]);
    }

    printf("particle initial mass: \n");
    for (int i =0; i < nlocal; i++)
    {
        printf("%d: %.5f \n",atom->tag[i],rmass0[i]);
    }
    printf("fix_composition address: %d \n", fix_rmass0->vector_atom);


    printf("particle initial density: \n");
    for (int i =0; i < nlocal; i++)
    {
        printf("%d: %.5f \n",atom->tag[i],density0[i]);
    }    
    printf("fix_composition address: %d \n", fix_density0->vector_atom);

    printf("particle initial radius: \n");
    for (int i =0; i < nlocal; i++)
    {
        printf("%d: %.5f \n",atom->tag[i],radius0[i]);
    }

    int *mask = atom->mask;
    for (int i = 0; i < nlocal; i++)
    {
        if(mask[i] & groupbit)
        {
            for(int j = 0; j < numComp; j++)
            {
                printf("********************\n");           
                printf("particle %d: %.4f, \n", atom->tag[i],fix_composition->array_atom[i][j]);
                printf("********************\n");                      
            }
        }
        
    }

    printf("fix_composition address: %d \n", fix_composition->array_atom);
    printf("just_created:%d, create_attribute:%d",fix_rmass0->just_created,fix_rmass0->create_attribute);
    printf("********************\n");

    #endif
    //printf("This is for debugging\n");
}

void FixReactionParticle0D::initial_integrate(int vflag)
{
    FixReaction::initial_integrate(vflag);

    //array or vector size may change w.r.t. particle number nmax
    // memory reallocation does not gurantee the same address of these storage variables
    //fix_rmass0->do_forward_comm(); //value does not change every step
    //fix_radius0->do_forward_comm();
    //fix_density0->do_forward_comm();
    rmass0 = fix_rmass0->vector_atom;
    density0 = fix_density0->vector_atom;
    radius0 = fix_radius0->vector_atom;
    composition = fix_composition->array_atom; 
    
    if (restart_flag)
    {
        fix_composition->do_forward_comm();
        restart_flag=false;
    }
//    else if (update->ntimestep%nevery_ == 0)
//        fix_composition->do_forward_comm();

}

/* ---------------------------------------------------------------------- */
void FixReactionParticle0D::pre_force(int vflag)
{
    if(neighbor->ago == 0)
    {
        FixReaction::pre_force(vflag);
        
        fix_rmass0->do_forward_comm();
        fix_radius0->do_forward_comm();
        fix_density0->do_forward_comm();
        fix_composition->do_forward_comm();

        //fix_composition->do_forward_comm();
    }   
}
/* ---------------------------------------------------------------------- */

void FixReactionParticle0D::post_force(int vflag)
{
    if(update->ntimestep % nevery_) 
        return;

    //printf("time step in post_force: %d\n",update->ntimestep);
    //caculating reaction rates
    double R = 8.31446; //unit J/mol/K
    int nlocal = atom->nlocal;
    int *mask = atom->mask;
    int *type = atom->type;

    double fractionSum = 0.;

    //reaction time
    double dt = update->dt;
    dt_ = dt * nevery_;
    
    FixReaction::updatePtrs();
    rmass0 = fix_rmass0->vector_atom;
    density0 = fix_density0->vector_atom;
    radius0 = fix_radius0->vector_atom;
    composition = fix_composition->array_atom; 

    for (int iparticle = 0; iparticle < nlocal; iparticle++)
    {
        if (type[iparticle] == ptype)
        {
            //copy 
            for(int j = 0; j < numComp; j++)
            {
                Y[j+1] = composition[iparticle][j];
            }

            Yold.assign(Y.begin(),Y.end()); //vector copying

            //update Y
            FixReactionParticle0D::ODEsolver( iparticle);        
           
            fractionSum = 0.0;
            //check composition fraction bounding  to [0,1]
            for(int j = 0; j < numComp; j++)
            {
                if(Y[j+1] < 1e-8)
                    Y[j+1] = 0.0;
                if(Y[j+1] > 1.0)
                    Y[j+1] = 1.0-1e-8;

                fractionSum += Y[j+1];
            }

            //update composition
            if(fractionSum!=0.0)
            {
                for(int j = 0; j < numComp; j++)
                {
                    composition[iparticle][j] = Y[j+1] / fractionSum; //make sure fraction summation is 1

                }
            }
        }

    } 

   // fix_composition->do_forward_comm();

}

/* ---------------------------------------------------------------------- */

void FixReactionParticle0D::end_of_step()
{
    if(update->ntimestep % nevery_) 
        return;

    //printf("time step in end_of_step: %d\n",update->ntimestep);           
    double* rmass = atom->rmass;
    double* density = atom->density;
    double* radius = atom->radius;
    int nlocal = atom->nlocal;
    int *mask = atom->mask;
    int *type = atom->type;
    double temp;

    for(int i = 0; i < nlocal; i++)
    {
        if(type[i] == ptype)
        {
            //fraction of solid and liquid phases that account towards to particle mass
            temp=0.;
            for(int j = 0; j < numSolid + numLiquid; j++)
            {
                temp += composition[i][j];
            }

            rmass[i] = MathExtraLiggghts::max(temp, 1e-6) * rmass0[i];
            density[i] = MathExtraLiggghts::max(temp, 1e-6) * density0[i]; //assuming particle volume does not change.
        }
    } 
}

/* ---------------------------------------------------------------------- */

void FixReactionParticle0D::cpl_evaluate(ComputePairGranLocal *caller)
{
    if(caller != cpl) error->all(FLERR,"Illegal situation in FixReactionParticle0D::cpl_evaluate");
}

/* ---------------------------------------------------------------------- */

void FixReactionParticle0D::register_compute_pair_local(ComputePairGranLocal *ptr)
{
   
   if(cpl != NULL)
      error->all(FLERR,"Fix reaction/particle0 allows only one compute of type pair/local");
   cpl = ptr;
}

/* ---------------------------------------------------------------------- */

void FixReactionParticle0D::unregister_compute_pair_local(ComputePairGranLocal *ptr)
{
   
   if(cpl != ptr)
       error->all(FLERR,"Illegal situation in FixReactionParticle0D::unregister_compute_pair_local");
   cpl = NULL;
}
