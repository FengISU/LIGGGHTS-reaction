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

#ifdef FIX_CLASS

FixStyle(reaction/particle0D, FixReactionParticle0D)
FixStyle(reaction, FixReactionParticle0D)

#else

#ifndef LMP_FIX_REACTION_PARTICLE0D_H
#define LMP_FIX_REACTION_PARTICLE0D_H

#include "fix_reaction.h"
#include "particleToInsert.h"
#include "atom.h"
#include "fix_property_atom_reaction.h"
#include "update.h"
//for using odeint ODE solver
#include <vector>
#include "boost/numeric/odeint.hpp"
#include "boost/ref.hpp"

using namespace std;
using namespace boost::numeric::odeint;

namespace LAMMPS_NS {

    class FixReactionParticle0D : public FixReaction {
    public:
        friend class LAMMPS_NS::ParticleToInsert;
        FixReactionParticle0D(class LAMMPS *, int, char **);
        ~FixReactionParticle0D();
        virtual void post_create();
        virtual void pre_delete(bool);
        virtual void init();
        virtual int setmask();
        virtual void initial_integrate(int vflag);
        virtual void pre_force(int);
        virtual void post_force(int);
        virtual void end_of_step();

        //refering to the code of fix_heat_gran
        // per default these three methods throw errors.
        virtual void cpl_evaluate(class ComputePairGranLocal *);
        virtual void register_compute_pair_local(class ComputePairGranLocal *);
        virtual void unregister_compute_pair_local(class ComputePairGranLocal *);
        
        int nevery_; //integrate only this many time steps
        bool restart_flag;

    protected:
        int iarg_;
        int ngrid; //intra-particle resolution
        double dt_; //integration time interval
        class FixPropertyAtom *fix_composition;
        class FixPropertyAtom *fix_rmass0; //record particle initial mass
        class FixPropertyAtom *fix_density0; //record particle initial density
        class FixPropertyAtom *fix_radius0; //record particle initial radius
        class FixPropertyAtomReaction* fix_capacity;
        class FixPropertyAtom *fix_flux;
        double **composition; //pointers for accessing fix property arrays
        double *rmass0; //pointers for accessing fix property vectors
        double *density0;
        double *radius0;
        
        double *H;
        double *delta_H;

        class FixPropertyAtom* fix_heatSource;
        double *heatSource;

        int numrun; //in case of multiple runs in an input script
    public:
        //For odeint wrapper
        typedef std::vector<double> odeint_state_type;
        odeint_state_type Y; //composition
        odeint_state_type Yold;
        odeint_state_type k; //reaction rate
        double temperature; // particle temperature;

        typedef runge_kutta_dopri5<odeint_state_type> error_stepper_type;
        typedef controlled_runge_kutta<error_stepper_type> controlled_stepper_type;
        controlled_stepper_type controlled_stepper;
        //overload operator () in which the ode functor is defined for odeint to call
        inline void operator()(const odeint_state_type &Y, odeint_state_type &dY, const double /* t */);
        
    private:
        inline void ODEsolver( int i ); 
        inline double* heat_reaction(const vector<double> &, const vector<double> &, int);
        typedef double* (FixReactionParticle0D::*heatscheme) (const vector<double>&,const vector<double>&, int i);
        
    };
    #include "fix_reaction_particle0D_i.h"
}
#endif
#endif
