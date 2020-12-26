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

#ifndef LMP_FIX_REACTION_H
#define LMP_FIX_REACTION_H

#include "fix.h"

namespace LAMMPS_NS {

    class FixReaction : public Fix {
    public:
        FixReaction(class LAMMPS *, int, char**);
        ~FixReaction();
        virtual void post_create();
        virtual void pre_delete (bool unfixflag){ UNUSED(unfixflag);}
        virtual void init();
        virtual int setmask();
        virtual void initial_integrate(int vflag);
        virtual void pre_force(int);
        virtual void post_force(int)=0;
        //virtual void end_of_step();

        //refering to the code of fix_heat_gran
        // per default these three methods throw errors.
        virtual void cpl_evaluate(class ComputePairGranLocal *);
        virtual void register_compute_pair_local(class ComputePairGranLocal *);
        virtual void unregister_compute_pair_local(class ComputePairGranLocal *);
        void updatePtrs();
        int check_numComp() { return numComp; }
        int check_numSolid() { return numSolid; }
        int check_numLiquid(){ return numLiquid; }
        int check_numGas() { return numGas; }
        int check_reactive_particle() { return ptype; }
    protected:
        int particleResolvedDim; // 0/1/2/3
        int ptype;
        int numComp;
        int numSolid;
        int numLiquid;
        int numGas;
        int numReaction;
        char *fractionType;

        class ComputePairGranLocal *cpl;
        class FixPropertyGlobal *fix_frequencyFactor;
        class FixPropertyGlobal *fix_activation;
        class FixPropertyGlobal *fix_composition0;
        class FixPropertyGlobal* fix_heatofReaction;
        class FixPropertyAtom *fix_tempold;
        class FixPropertyAtom *fix_temp;

        double *frequencyFactor;
        double *activation;
        double *composition0;
        double *temp;
        double *tempOld;
        double *heatReaction;

        class PairGran *pair_gran;
        int history_flag;


    };
}

#endif