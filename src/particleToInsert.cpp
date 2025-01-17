
/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

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
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include "particleToInsert.h"
#include "math.h"
#include "error.h"
#include "update.h"
#include "domain.h"
#include "atom.h"
#include "atom_vec.h"
#include "fix_property_atom.h"
#include "fix_property_atom_reaction.h"
#include "vector_liggghts.h"
#include "math_extra_liggghts.h"
#include "modify.h"

#include "fix_reaction_particle0D.h"
#include "fix_heat_gran.h"
using namespace LAMMPS_NS;

ParticleToInsert::ParticleToInsert(LAMMPS* lmp,int ns) : Pointers(lmp)
{
        groupbit = 0;

        distorder = -1;

        nspheres = ns;

        memory->create(x_ins,nspheres,3,"x_ins");
        radius_ins = new double[nspheres];

        atom_type_vector = new int[nspheres];
        atom_type_vector_flag = false;

        fix_property = 0;
        fix_property_value = 0.;
}

/* ---------------------------------------------------------------------- */

ParticleToInsert::~ParticleToInsert()
{
        memory->destroy(x_ins);
        delete []radius_ins;
        delete []atom_type_vector;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsert::insert()
{
    // perform the actual insertion
    // add particles, set coordinate and radius
    // set group mask to "all" plus fix groups

    int inserted = 0;
    int nfix = modify->nfix;
    Fix **fix = modify->fix;
        //add by Fenglei Qi, 2016/8/15
    FixPropertyAtomReaction* fix_conductivity_=NULL;
    FixPropertyAtomReaction* fix_capacity_ = NULL;
    FixReactionParticle0D* fix_reaction = static_cast<FixReactionParticle0D*>(modify->find_fix_style("reaction", 0));
    FixHeatGran* fix_heat = static_cast<FixHeatGran*>(modify->find_fix_style("heat/gran", 0));
    char style[] = "ParitcleToInsert";
    if (fix_reaction || fix_heat)
    {
    	fix_conductivity_ =
        static_cast<FixPropertyAtomReaction*>
        (
            modify->find_fix_property("thermalConductivity","property/atom/reaction","scalar",1,0,style)
        );
    	fix_capacity_ =
        static_cast<FixPropertyAtomReaction*>
        (
            modify->find_fix_property("thermalCapacity","property/atom/reaction","scalar",1,0,style)
        );
    }

    for(int i = 0; i < nspheres; i++)
    {
        
        //if (domain->is_in_extended_subdomain(x_ins[i]))
        //{
                
                inserted++;
                if(atom_type_vector_flag)
                    atom->avec->create_atom(atom_type_vector[i],x_ins[i]);
                else
                    atom->avec->create_atom(atom_type,x_ins[i]);
                int m = atom->nlocal - 1;
                atom->mask[m] = 1 | groupbit;
                vectorCopy3D(v_ins,atom->v[m]);
                vectorCopy3D(omega_ins,atom->omega[m]);
                atom->radius[m] = radius_ins[i];
                atom->density[m] = density_ins;
                
                atom->rmass[m] = (1==nspheres)? (mass_ins) : (4.18879020479/*4//3*pi*/*radius_ins[i]*radius_ins[i]*radius_ins[i]*density_ins);

                //pre_set_arrays() called above
                for (int j = 0; j < nfix; j++)
                   if (fix[j]->create_attribute) fix[j]->set_arrays(m);

                // apply fix property setting coming from fix insert
                // this overrides the set_arrays call above
                if(fix_property)
                    fix_property->vector_atom[m] = fix_property_value;
        
                //fix_property->do_forward_comm(); //do forward communication, added by Fenglei Qi
           
             //add by Fenglei Qi, 2016/8/15
            if(fix_reaction)
            {
                if(atom->type[m] == fix_reaction->ptype)
                {
                    fix_reaction->fix_rmass0->vector_atom[m] = atom->rmass[m];
                    fix_reaction->fix_radius0->vector_atom[m] = atom->radius[m];
                    fix_reaction->fix_density0->vector_atom[m] = atom->density[m];

                    for(int j = 0; j < fix_reaction->numComp; j++)
                    {
                        fix_reaction->fix_composition->array_atom[m][j]= fix_reaction->composition0[j];
                    }
                }
                else
                {
                    fix_reaction->fix_rmass0->vector_atom[m] = 0.0;
                    fix_reaction->fix_radius0->vector_atom[m] = 0.0;
                    fix_reaction->fix_density0->vector_atom[m] = 0.0;

                    for(int j = 0; j < fix_reaction->numComp; j++)
                    {
                        fix_reaction->fix_composition->array_atom[m][j]= 0.0;
                    }                    
                } 

                //fix_reaction->fix_composition->do_forward_comm(); 
            } 

            if(fix_conductivity_)
            {
                fix_conductivity_->vector_atom[m] = fix_conductivity_->defaultvalues[atom->type[m]-1];
                //fix_conductivity_->do_forward_comm(); //will reneighbor after insertion
            }
            if(fix_capacity_)
            {
                fix_capacity_->vector_atom[m] = fix_capacity_->defaultvalues[atom->type[m]-1];
                //fix_capacity_->do_forward_comm(); //will reneighbor after insertion

            }
        //}

    }

    #if 0
    printf("*********************\n");
    double** temp = fix_reaction->fix_composition->array_atom;
    int m = atom->nlocal - 1;
    //for (int i =0; i<atom->nlocal;i++)
    //    printf("atom(type %d) %d: %.5f, %.5f, %.5f\n",atom->type[i], atom->tag[i],temp[i][0],temp[i][1],temp[][2]);
    printf("atom(type %d) %d: %.5f, %.5f, %.5f\n",atom->type[m], atom->tag[m],temp[m][0],temp[m][1],temp[m][2]);

    printf("*********************\n");
    #endif
    
    return inserted;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsert::check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, double **xnear, int &nnear)
{
    if(nspheres > 1)
        return check_near_set_x_v_omega_ms(x,v, omega,quat,xnear,nnear);

    // check sphere against all others in xnear
    // if no overlap add to xnear
    double del[3], rsq, radsum;

    vectorCopy3D(x,x_ins[0]);

    for(int i = 0; i < nnear; i++)
    {
        vectorSubtract3D(x_ins[0],xnear[i],del);
        rsq = vectorMag3DSquared(del);
        
        radsum = radius_ins[0] + xnear[i][3];

        // no success in overlap
        if (rsq <= radsum*radsum) return 0;
    }

    // no overlap with any other - success

    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    // add to xnear
    vectorCopy3D(x_ins[0],xnear[nnear]);
    xnear[nnear][3] = radius_ins[0];
    nnear++;

    return 1;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsert::check_near_set_x_v_omega_ms(double *x,double *v, double *omega, double *quat, double **xnear, int &nnear)
{
    // x is position where insertion should take place
    // v and omega are the velocity and omega for the newly inserted particles
    double rel[3],xins_j_try[3];
    double del[3], rsq, radsum;

    // check insertion position, take quat into account
    // relative position of spheres to each other already stored at this point
    // check sphere against all others in xnear
    for(int j = 0; j < nspheres; j++)
    {
        // take orientation into account; x_bound_ins is in the global coordinate system
        // calculate xins_j_try for every sphere and check if would work
        vectorSubtract3D(x_ins[j],x_bound_ins,rel);
        MathExtraLiggghts::vec_quat_rotate(rel,quat);
        vectorAdd3D(rel,x,xins_j_try);

        for(int i = 0; i < nnear; i++)
        {
           vectorSubtract3D(xins_j_try,xnear[i],del);
           rsq = vectorMag3DSquared(del);
           radsum = radius_ins[j] + xnear[i][3];

           // no success in overlap
           if (rsq <= radsum*radsum)
            return 0;
        }
    }

    // no overlap with any other - success
    // set x_ins, v_ins and omega_ins
    for(int j = 0; j < nspheres; j++)
    {
        vectorSubtract3D(x_ins[j],x_bound_ins,rel);
        MathExtraLiggghts::vec_quat_rotate(rel,quat);
        vectorAdd3D(rel,x,x_ins[j]);
    }
    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    // add to xnear for future checks
    for(int j = 0; j < nspheres; j++)
    {
        vectorCopy3D(x_ins[j],xnear[nnear]);
        xnear[nnear][3] = radius_ins[j];
        nnear++;
    }

    return nspheres;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsert::set_x_v_omega(double *x, double *v, double *omega, double *quat)
{
    double rel[3];

    // x is position where insertion should take place
    // v and omega are the velocity and omega for the newly inserted particles

    // add insertion position
    // relative position of spheres to each other already stored at this point
    // also take quat into account
    for(int j = 0; j < nspheres; j++)
    {
        // if only one sphere, then x_bound = x_ins and there is
        // no relevant orientation
        if(1 == nspheres)
            vectorAdd3D(x_ins[j],x,x_ins[j]);
        // if > 1 sphere, take orientation into account
        // x_bound_ins is in the global coordinate system
        else
        {
            vectorSubtract3D(x_ins[j],x_bound_ins,rel);
            MathExtraLiggghts::vec_quat_rotate(rel,quat);
            vectorAdd3D(rel,x,x_ins[j]);
        }
    }

    // set velocity and omega
    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    return nspheres;
}

/* ---------------------------------------------------------------------- */

void ParticleToInsert::scale_pti(double r_scale)
{
    double r_scale3 = r_scale*r_scale*r_scale;

    for(int i = 0; i < nspheres; i++)
    {
        radius_ins[i] *= r_scale;
        vectorScalarMult3D(x_ins[i],r_scale);
    }

    volume_ins *= r_scale3;
    mass_ins *= r_scale3;

    r_bound_ins *= r_scale;

    vectorScalarMult3D(x_bound_ins,r_scale);
}
