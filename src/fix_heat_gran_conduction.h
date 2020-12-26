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

#ifdef FIX_CLASS

FixStyle(heat/gran/conduction,FixHeatGranCond)
FixStyle(heat/gran,FixHeatGranCond)

#else

#ifndef LMP_FIX_HEATGRAN_CONDUCTION_H
#define LMP_FIX_HEATGRAN_CONDUCTION_H

#include "fix_heat_gran.h"

namespace LAMMPS_NS {

  class FixHeatGranCond : public FixHeatGran {
  public:
    FixHeatGranCond(class LAMMPS *, int, char **);
    ~FixHeatGranCond();
    virtual void post_create();
    void pre_delete(bool);

    int setmask();
    void init();
    virtual void post_force(int);

    void cpl_evaluate(class ComputePairGranLocal *);
    void register_compute_pair_local(ComputePairGranLocal *);
    void unregister_compute_pair_local(ComputePairGranLocal *);

  protected:
    int iarg_;

  private:
    template <int,int> void post_force_eval(int,int);


    class FixPropertyAtomReaction* fix_conductivity_;
    double *conductivity_;

    //considering the properties of interstitial fluids, Fenglei
    class FixPropertyGlobal *fix_fluid_conductivity_; //consider property change with temp
    double fluid_conductivity_;

    //considering the heat radiation, Fenglei QI
    class FixPropertyGlobal *fix_radiation_emissivity;
    double *emissivity_;

    // model for contact area calculation
    int area_calculation_mode_;

    double fixed_contact_area_;

    // for heat transfer area correction
    int area_correction_flag_;
    double const* const* deltan_ratio_;

    //adding struct data type
    struct pfp_param //particle-fluid-particle pathsway(pfp)
    {
    	double lowBound,upBound;
    	double condi,condj, condf;
    	double radi,radj;
    	double rad_eff;
    	double half_gap;

    	pfp_param():
    	             lowBound(0.0),upBound(0.0),
    	             condi(0.0),condj(0.0),condf(0.0),
    	             radi (0.0),radj(0.0),
    	             rad_eff(0.0),
    	             half_gap(0.0) {}
    };
    //adding functions
    inline double pfp_heat_transfer(const FixHeatGranCond::pfp_param &pp, const double &rr);
    typedef double (FixHeatGranCond::*hfun)(const FixHeatGranCond::pfp_param &pp, const double &rr);
    double num_integrate(FixHeatGranCond::hfun fun,const FixHeatGranCond::pfp_param &pfp,int numIntervals);
    inline double viewfactor_correlation_uncontacted(const FixHeatGranCond::pfp_param &pp); //Fenglei QI, 2016/2/21
    inline double viewfactor_correlation_contacted(const FixHeatGranCond::pfp_param &pp); //Fenglei QI, 2016/2/21
  };

}

#endif
#endif

