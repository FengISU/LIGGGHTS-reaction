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

#include "fix_heat_gran_conduction.h"

#include "atom.h"
#include "compute_pair_gran_local.h"
#include "fix_property_atom.h"
#include "fix_property_atom_reaction.h"
#include "fix_property_global.h"
#include "force.h"
#include "math_extra.h"
#include "math_extra_liggghts.h"
#include "properties.h"
#include "modify.h"
#include "neigh_list.h"
#include "pair_gran.h"
#include "update.h" //added by Fenglei Qi

#include "comm.h"

using namespace LAMMPS_NS;
using namespace FixConst;
// These definitions are used in particle-fluid-particle heat transfer caculation
#define CONFIG_ENABLE_DEBUG_LOG 0
#define EPSILON 0.44
#define INTERVALS 8 //Fenglei QI, 2016/1/27
#define DMAX_FACTOR 0.5 //Fenglei QI, 2016/1/27
#define RADIATION 1  //Fenglei Qi, 2016/2/20
//#define EMISSIVITY 0.9 //Fenglei Qi, 2016/2/20
// modes for conduction contact area calaculation
// same as in fix_wall_gran.cpp

enum{ CONDUCTION_CONTACT_AREA_OVERLAP,
      CONDUCTION_CONTACT_AREA_CONSTANT,
      CONDUCTION_CONTACT_AREA_PROJECTION};

/* ---------------------------------------------------------------------- */

FixHeatGranCond::FixHeatGranCond(class LAMMPS *lmp, int narg, char **arg) :
  FixHeatGran(lmp, narg, arg),
  fix_conductivity_(0),
  fix_fluid_conductivity_(0),
  fix_radiation_emissivity(0),
  conductivity_(0),
  fluid_conductivity_(0),
  emissivity_(0),
  area_calculation_mode_(CONDUCTION_CONTACT_AREA_OVERLAP),
  fixed_contact_area_(0.),
  area_correction_flag_(0),
  deltan_ratio_(0)
{
  iarg_ = 5;

  bool hasargs = true;
  while(iarg_ < narg && hasargs)
  {
    hasargs = false;
    if(strcmp(arg[iarg_],"contact_area") == 0) {

      if(strcmp(arg[iarg_+1],"overlap") == 0)
        area_calculation_mode_ =  CONDUCTION_CONTACT_AREA_OVERLAP;
      else if(strcmp(arg[iarg_+1],"projection") == 0)
        area_calculation_mode_ =  CONDUCTION_CONTACT_AREA_PROJECTION;
      else if(strcmp(arg[iarg_+1],"constant") == 0)
      {
        if (iarg_+3 > narg)
            error->fix_error(FLERR,this,"not enough arguments for keyword 'contact_area constant'");
        area_calculation_mode_ =  CONDUCTION_CONTACT_AREA_CONSTANT;
        fixed_contact_area_ = force->numeric(FLERR,arg[iarg_+2]);
        if (fixed_contact_area_ <= 0.)
            error->fix_error(FLERR,this,"'contact_area constant' value must be > 0");
        iarg_++;
      }
      else error->fix_error(FLERR,this,"expecting 'overlap', 'projection' or 'constant' after 'contact_area'");
      iarg_ += 2;
      hasargs = true;
    } else if(strcmp(arg[iarg_],"area_correction") == 0) {
      if (iarg_+2 > narg) error->fix_error(FLERR,this,"not enough arguments for keyword 'area_correction'");
      if(strcmp(arg[iarg_+1],"yes") == 0)
        area_correction_flag_ = 1;
      else if(strcmp(arg[iarg_+1],"no") == 0)
        area_correction_flag_ = 0;
      else error->fix_error(FLERR,this,"expecting 'yes' otr 'no' after 'area_correction'");
      iarg_ += 2;
      hasargs = true;
    } else if(strcmp(style,"heat/gran/conduction") == 0)
        error->fix_error(FLERR,this,"unknown keyword");
  }

  if(CONDUCTION_CONTACT_AREA_OVERLAP != area_calculation_mode_ && 1 == area_correction_flag_)
    error->fix_error(FLERR,this,"can use 'area_correction' only for 'contact_area = overlap'");
}

/* ---------------------------------------------------------------------- */

FixHeatGranCond::~FixHeatGranCond()
{

  //if (conductivity_)
  //  delete conductivity_;
  if (emissivity_)
    delete [] emissivity_;
}

/* ---------------------------------------------------------------------- */

void FixHeatGranCond::post_create()
{
  FixHeatGran::post_create();
}

/* ---------------------------------------------------------------------- */

void FixHeatGranCond::pre_delete(bool unfixflag)
{

  // tell cpl that this fix is deleted
  if(cpl && unfixflag) cpl->reference_deleted();

}

/* ---------------------------------------------------------------------- */

int FixHeatGranCond::setmask()
{
  int mask = FixHeatGran::setmask();
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixHeatGranCond::init()
{
  // initialize base class
  FixHeatGran::init();
  const double *Y, *nu, *Y_orig;
  double expo, Yeff_ij, Yeff_orig_ij, ratio;
  int max_type = atom->get_properties()->max_type();
  //if (conductivity_) delete []conductivity_;
  if (emissivity_) delete [] emissivity_; //addding by Fenglei Qi, 3/18/2016
  //conductivity_ = new double[max_type];
  emissivity_ =new double[max_type]; //addding by Fenglei Qi, 3/18/2016
  fix_conductivity_ =
  static_cast<FixPropertyAtomReaction*>
    (
      modify->find_fix_property("thermalConductivity","property/atom/reaction","scalar",1,0,style)
    );
  //adding fluid properties-conductivity
  fix_fluid_conductivity_ =
    static_cast<FixPropertyGlobal*> (modify->find_fix_property("fluidThermalConductivity","property/global","scalar",1,0,style));
  // pre-calculate conductivity for possible contact material combinations
  //for(int i=1;i< max_type+1; i++)
  //    for(int j=1;j<max_type+1;j++)
  //    {
  //        conductivity_[i-1] = fix_conductivity_->compute_vector(i-1);
  //        if(conductivity_[i-1] < 0.)
  //          error->all(FLERR,"Fix heat/gran/conduction: Thermal conductivity must not be < 0");
  //    }
  conductivity_ = fix_conductivity_->vector_atom;
  //adding fluid properties-conductivity
  fluid_conductivity_ = fix_fluid_conductivity_->compute_scalar();
  // calculate heat transfer correction
  //adding radiation emissivity, Fenglei Qi
  #if RADIATION
    fix_radiation_emissivity =
        static_cast<FixPropertyGlobal*>(modify->find_fix_property("emissivity","property/global","peratomtype",max_type,0,style));  
    for(int i=1;i<max_type+1;i++)
      {
          emissivity_[i-1] = fix_radiation_emissivity->compute_vector(i-1);
          if(emissivity_[i-1] < 0.)
            error->all(FLERR,"Fix heat/gran/conduction: Thermal radiation emissivity must not be < 0");
      }
  #endif
  if(area_correction_flag_)
  {
    if(!force->pair_match("gran",0))
        error->fix_error(FLERR,this,"area correction only works with using granular pair styles");

    //expo = 1./pair_gran->stressStrainExponent();
      expo =0.2;
    Y = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0,style))->get_values();
    nu = static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0,style))->get_values();
    Y_orig = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,style))->get_values();

    // allocate a new array within youngsModulusOriginal
    static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,style))->new_array(max_type,max_type);

    // feed deltan_ratio into this array
    for(int i = 1; i < max_type+1; i++)
    {
      for(int j = 1; j < max_type+1; j++)
      {
        Yeff_ij      = 1./((1.-pow(nu[i-1],2.))/Y[i-1]     +(1.-pow(nu[j-1],2.))/Y[j-1]);
        Yeff_orig_ij = 1./((1.-pow(nu[i-1],2.))/Y_orig[i-1]+(1.-pow(nu[j-1],2.))/Y_orig[j-1]);
        ratio = pow(Yeff_ij/Yeff_orig_ij,expo);
        
        static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,style))->array_modify(i-1,j-1,ratio);
      }
    }

    // get reference to deltan_ratio
    deltan_ratio_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,style))->get_array_modified();
  }

  updatePtrs();

  // error checks on coarsegraining
  
}

/* ---------------------------------------------------------------------- */

void FixHeatGranCond::post_force(int vflag)
{

  if(history_flag == 0 && CONDUCTION_CONTACT_AREA_OVERLAP == area_calculation_mode_)
    post_force_eval<0,CONDUCTION_CONTACT_AREA_OVERLAP>(vflag,0);
  if(history_flag == 1 && CONDUCTION_CONTACT_AREA_OVERLAP == area_calculation_mode_)
    post_force_eval<1,CONDUCTION_CONTACT_AREA_OVERLAP>(vflag,0);

  if(history_flag == 0 && CONDUCTION_CONTACT_AREA_CONSTANT == area_calculation_mode_)
    post_force_eval<0,CONDUCTION_CONTACT_AREA_CONSTANT>(vflag,0);
  if(history_flag == 1 && CONDUCTION_CONTACT_AREA_CONSTANT == area_calculation_mode_)
    post_force_eval<1,CONDUCTION_CONTACT_AREA_CONSTANT>(vflag,0);

  if(history_flag == 0 && CONDUCTION_CONTACT_AREA_PROJECTION == area_calculation_mode_)
    post_force_eval<0,CONDUCTION_CONTACT_AREA_PROJECTION>(vflag,0);
  if(history_flag == 1 && CONDUCTION_CONTACT_AREA_PROJECTION == area_calculation_mode_)
    post_force_eval<1,CONDUCTION_CONTACT_AREA_PROJECTION>(vflag,0);
}

/* ---------------------------------------------------------------------- */

void FixHeatGranCond::cpl_evaluate(ComputePairGranLocal *caller)
{
  if(caller != cpl) error->all(FLERR,"Illegal situation in FixHeatGranCond::cpl_evaluate");

  if(history_flag == 0 && CONDUCTION_CONTACT_AREA_OVERLAP == area_calculation_mode_)
    post_force_eval<0,CONDUCTION_CONTACT_AREA_OVERLAP>(0,1);
  if(history_flag == 1 && CONDUCTION_CONTACT_AREA_OVERLAP == area_calculation_mode_)
    post_force_eval<1,CONDUCTION_CONTACT_AREA_OVERLAP>(0,1);

  if(history_flag == 0 && CONDUCTION_CONTACT_AREA_CONSTANT == area_calculation_mode_)
    post_force_eval<0,CONDUCTION_CONTACT_AREA_CONSTANT>(0,1);
  if(history_flag == 1 && CONDUCTION_CONTACT_AREA_CONSTANT == area_calculation_mode_)
    post_force_eval<1,CONDUCTION_CONTACT_AREA_CONSTANT>(0,1);

  if(history_flag == 0 && CONDUCTION_CONTACT_AREA_PROJECTION == area_calculation_mode_)
    post_force_eval<0,CONDUCTION_CONTACT_AREA_PROJECTION>(0,1);
  if(history_flag == 1 && CONDUCTION_CONTACT_AREA_PROJECTION == area_calculation_mode_)
    post_force_eval<1,CONDUCTION_CONTACT_AREA_PROJECTION>(0,1);
}



/* ---------------------------------------------------------------------- */
//Adding new heat transfer model
// This heat transfer model includes particle-fluid-particle model (uncontacted/contacted)
// and particle-particle conduction heat transfer (contacted)
// The mathematical model refers to 
//*******************
// Qinfu Hou, Jieqing Gan, Zongyan Zhou, and Aibing Yu (2015).Particle scale study of heat
// transfer in packed and fluidized beds. Chapter four in
//Advances in Chemical Eingeering, Volume 46
//*******************
//define two parameters:
// solid fraction, epsilon=0.5;
// In the uncontacted scenario, the largest distance between particle considered is
// dmax=2*sqrt[2]*r (only applied to monodisperse system)

template <int HISTFLAG,int CONTACTAREA>
void FixHeatGranCond::post_force_eval(int vflag,int cpl_flag)
{
  //debugging
  //if(comm->me == 0)
  //printf("time step into fix heat gran: %d\n", update->ntimestep);
  double hc,hc_pfp,hnc_pfp,contactArea,contactRad,tmpcontactRad,delta_n,flux,dirFlux[3]; // add contactRad,hc_pfp variable
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double radi,radj,radsum,rsq,r,tcoi,tcoj;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *contact_flag,**first_contact_flag;
  double tmp; //add temp variable for use

  int newton_pair = force->newton_pair;

  if (strcmp(force->pair_style,"hybrid")==0)
    error->warning(FLERR,"Fix heat/gran/conduction implementation may not be valid for pair style hybrid");
  if (strcmp(force->pair_style,"hybrid/overlay")==0)
    error->warning(FLERR,"Fix heat/gran/conduction implementation may not be valid for pair style hybrid/overlay");

  inum = pair_gran->list->inum;
  ilist = pair_gran->list->ilist;
  numneigh = pair_gran->list->numneigh;
  firstneigh = pair_gran->list->firstneigh;
  if(HISTFLAG) first_contact_flag = pair_gran->listgranhistory->firstneigh;

  double *radius = atom->radius;
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  updatePtrs();
  conductivity_ = fix_conductivity_->vector_atom;
  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    if(HISTFLAG) contact_flag = first_contact_flag[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      if (!(mask[i] & groupbit) && !(mask[j] & groupbit)) continue;

      //if(!HISTFLAG)
      //{
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        radj = radius[j];
        radsum = radi + radj;
      //}
      
      //debugging 
      #if CONFIG_ENABLE_DEBUG_LOG
          if (update->ntimestep %10000 ==0)
          {
            FixPropertyAtomReaction* fix_capacity =  static_cast<FixPropertyAtomReaction*>(modify->find_fix_property("thermalCapacity","property/atom/reaction","scalar",1,0,style));
            char message[500];           
            sprintf(message,"Parameters in heat transfer: i: %d, j: %d, nlocal: %d (%d cpu), condi: %.4e, condj: %.4e, cpi:%.4e,cpj:%.4e\n",
              i,j,atom->nlocal,comm->me,conductivity_[i],conductivity_[j], fix_capacity->vector_atom[i],fix_capacity->vector_atom[j]);
            error->message(FLERR,message,1);
          }
      #endif


      //if ((HISTFLAG && contact_flag[jj]) || (!HISTFLAG && (rsq < radsum*radsum))) {  //contact
      if (rsq<radsum*radsum) { //contact 
        //if(HISTFLAG)
        //{
        //  delx = xtmp - x[j][0];
        //  dely = ytmp - x[j][1];
        //  delz = ztmp - x[j][2];
        //  rsq = delx*delx + dely*dely + delz*delz;
        //  radj = radius[j];
        //  radsum = radi + radj;
        //  if(rsq >= radsum*radsum) continue;
        //}

        r = sqrt(rsq);
        //delta_n = radsum - r;//added by Fenglei Qi, 2015.12.18
        if(CONTACTAREA == CONDUCTION_CONTACT_AREA_OVERLAP)
        {
            tmp = radi*radi - radj*radj + rsq;
            tmp *= 0.5 / r; //r should not be zero
            contactRad = sqrt(radi*radi-tmp*tmp);
            tmpcontactRad = contactRad;
            delta_n = (radsum - r)/2.0;//added by Fenglei Qi, 2015.12.18
            if(area_correction_flag_)
            {
              contactRad *= deltan_ratio_[type[i]-1][type[j]-1];
              double minrad = MathExtraLiggghts::min(radi,radj);
              delta_n = minrad - sqrt(minrad*minrad - contactRad*contactRad);
            }
            #if 0
            if(area_correction_flag_)
            {
              //delta_n = radsum - r;//remove this line
              delta_n *= deltan_ratio_[type[i]-1][type[j]-1];
              r = radsum - delta_n;
            }

            //contact area of the two spheres
            //contactArea = - M_PI/4 * ( (r-radi-radj)*(r+radi-radj)*(r-radi+radj)*(r+radi+radj) )/(r*r);
            // Using Contact radius
            tmp = radi*radi - radj*radj + r*r;
            tmp *= 1.0 / 2.0 / r; //r should not be zero
            contactRad = sqrt(radi*radi-tmp*tmp);
            #endif
        }
        else if (CONTACTAREA == CONDUCTION_CONTACT_AREA_CONSTANT)
        {
            contactArea = fixed_contact_area_;
            contactRad	= sqrt(contactArea/M_PI); //Using contact radius for consistency
            
        }
        else if (CONTACTAREA == CONDUCTION_CONTACT_AREA_PROJECTION)
        {
            double rmax = MathExtraLiggghts::max(radi,radj);
            contactArea = M_PI*rmax*rmax;
            contactRad	= rmax; //Using contact radius for consistency
        }
        tcoi = conductivity_[i];
        tcoj = conductivity_[j];

        if (tcoi < SMALL || tcoj < SMALL) hc = 0.;
        //else hc = 4.*tcoi*tcoj/(tcoi+tcoj)*sqrt(contactArea);
        else hc = 4.*tcoi*tcoj/(tcoi+tcoj)*contactRad;
        flux = (Temp[j]-Temp[i])*hc;

        heatFluxDistribution[i][0] += flux;
        if (newton_pair || j < nlocal)
        {
        	heatFluxDistribution[j][0] -= flux;
        }

        //Heat transfer through particle-fluid-particle pattern in particle contacting scenarioF
        struct pfp_param contacting_pfp;
        
        contacting_pfp.radi = radi;
        contacting_pfp.radj = radj;
        contacting_pfp.condi = tcoi;
        contacting_pfp.condj = tcoj;
        contacting_pfp.condf = fluid_conductivity_;
        double minrad = MathExtraLiggghts::min(radi,radj);
        contacting_pfp.rad_eff = 0.560 * minrad * pow(1.0-EPSILON,-1.0/3.0);
        contacting_pfp.half_gap = -1.0*delta_n; //only valid with CONDUCTION_CONTACT_AREA_OVERLAP model
        //contacting_pfp.half_gap =(sqrt(rsq)-radi-radj)/2.0;
        contacting_pfp.lowBound = contactRad;
        contacting_pfp.upBound = minrad * contacting_pfp.rad_eff;
        contacting_pfp.upBound /= sqrt(contacting_pfp.rad_eff*contacting_pfp.rad_eff+(minrad+contacting_pfp.half_gap)*(minrad+contacting_pfp.half_gap)); //Fenglei QI, 2016/01/27
        
        hc_pfp = num_integrate(&FixHeatGranCond::pfp_heat_transfer,contacting_pfp,INTERVALS);
        #if CONFIG_ENABLE_DEBUG_LOG
          if (update->ntimestep %1000 ==0)
          {

            char message[500];
            sprintf(message,"Particle-particle interaction: scaled dist between particle %d and %d: dist: %.6e, contactrad: %.6e (orig: %.6e); corrected gap: %.6e\n",
            atom->tag[i],atom->tag[j],radsum,contactRad,tmpcontactRad,delta_n);
            error->message(FLERR,message,1);
            
            sprintf(message,"Parameters in heat transfer: condi: %.4e, condj: %.4e, condg: %.4e,lowB:%.4e, upB: %.4e, radeff: %.4e\n",
                    tcoi,tcoj,fluid_conductivity_,contacting_pfp.lowBound,contacting_pfp.upBound,contacting_pfp.rad_eff);
            error->message(FLERR,message,1);
            
            sprintf(message,"Parameters in heat transfer: hc: %.4e, hc_pfp: %.4e", hc,hc_pfp);
            error->message(FLERR,message,1);
          }
        #endif
        tmp = (Temp[j]-Temp[i])*hc_pfp;
        flux += tmp; //separately storing???

        heatFluxDistribution[i][1] += tmp;
        if (newton_pair || j < nlocal)
        {
        	heatFluxDistribution[j][1] -= tmp;
        }

        //Heat transfer through particle-surface radiation in particle contacting scenario
        if(RADIATION) //Fenglei QI, 2016/2/21
        {
          double emissivityi = emissivity_[type[i]-1];
          double emissivityj = emissivity_[type[j]-1];
          double capheight = minrad - sqrt(minrad*minrad-contacting_pfp.upBound*contacting_pfp.upBound);
          double radiationarea = 2.0*M_PI*minrad*(capheight-delta_n);
          double viewfactor = viewfactor_correlation_contacted(contacting_pfp);
          double coefficent = 5.670367e-8*radiationarea/((1.0-emissivityi)/emissivityi+(1.0-emissivityj)/emissivityj+2.0/(1.0+viewfactor)); //Fenglei QI, 3/20/2016
          double tempi=Temp[i]*Temp[i];
          double tempj=Temp[j]*Temp[j];
          double radheat = coefficent*(tempj*tempj - tempi*tempi);
          flux += radheat;

          heatFluxDistribution[i][3] += radheat;
          if (newton_pair || j < nlocal)
          {
            heatFluxDistribution[j][3] -= radheat;
          } 
        }

        //debugging
         if (!(flux== flux) || (flux >1.e8) ) // NAN or infi
        {
            printf("Error with flux caculation \n");
            printf("Particle in %d node and with id %d \n", comm->me, atom->tag[i]);
            char message[500];
            sprintf(message,"Parameters in pfp heat transfer: condi: %.4e, condj: %.4e, condg: %.4e,lowB:%.4e, upB: %.4e, radeff: %.4e\n",
                    tcoi,tcoj,fluid_conductivity_,contacting_pfp.lowBound,contacting_pfp.upBound,contacting_pfp.rad_eff);
            error->message(FLERR,message,1);
            printf("pp: %.4e, pfp1: %.4e, pfp2: %.4e, radiation: %.4e\n",heatFluxDistribution[i][0],heatFluxDistribution[i][1],heatFluxDistribution[i][2],heatFluxDistribution[i][3]);
            printf("particle %d, x %.4e, y %.4e,z %.4e \n", atom->tag[i],x[i][0],x[i][1],x[i][1]);
            printf("particle %d, x %.4e, y %.4e,z %.4e \n", atom->tag[j],x[j][0],x[j][1],x[j][1]);
            printf("hc %.4e, hc_pfp %.4e, tempi %.4e, tempj %.4e", hc,hc_pfp,Temp[i],Temp[j]);
            printf("\n");
        }

        dirFlux[0] = flux*delx;
        dirFlux[1] = flux*dely;
        dirFlux[2] = flux*delz;

        if(!cpl_flag)
        {
          //Add half of the flux (located at the contact) to each particle in contact
          heatFlux[i] += flux;
          directionalHeatFlux[i][0] += 0.50 * dirFlux[0];
          directionalHeatFlux[i][1] += 0.50 * dirFlux[1];
          directionalHeatFlux[i][2] += 0.50 * dirFlux[2];
          if (newton_pair || j < nlocal)
          {
            heatFlux[j] -= flux;
            directionalHeatFlux[j][0] += 0.50 * dirFlux[0];
            directionalHeatFlux[j][1] += 0.50 * dirFlux[1];
            directionalHeatFlux[j][2] += 0.50 * dirFlux[2];
          }
        
        }

        if(cpl_flag && cpl) cpl->add_heat(i,j,flux); //what is the purpose
      }
      else //uncontacted
      {
      	//Heat transfer through particle-fluid-particle pattern in particle noncontacting scenario
        //Justify if the particles are considerring to conduct heat transfer
        double half_gap = (sqrt(rsq) - radsum)/2.0;
        double minrad=MathExtraLiggghts::min(radi,radj);
        if (half_gap/minrad > DMAX_FACTOR) continue;

        tcoi = conductivity_[i];
        tcoj = conductivity_[j];

        struct pfp_param contacting_pfp;
        contacting_pfp.radi = radi;
        contacting_pfp.radj = radj;
        contacting_pfp.condi = tcoi;
        contacting_pfp.condj = tcoj;
        contacting_pfp.condf = fluid_conductivity_;
        contacting_pfp.rad_eff = 0.560 * minrad * pow(1.0-EPSILON,-1.0/3.0);
        contacting_pfp.half_gap =half_gap;
        contacting_pfp.lowBound = 0.0;
        contacting_pfp.upBound = minrad * contacting_pfp.rad_eff;
        contacting_pfp.upBound /= sqrt(contacting_pfp.rad_eff*contacting_pfp.rad_eff+(minrad+contacting_pfp.half_gap)*(minrad+contacting_pfp.half_gap));
        
        hnc_pfp = num_integrate(&FixHeatGranCond::pfp_heat_transfer,contacting_pfp,INTERVALS);

        flux = (Temp[j]-Temp[i])*hnc_pfp; //separately storing???

        heatFluxDistribution[i][2] += flux;
        if (newton_pair || j < nlocal)
        {
        	heatFluxDistribution[j][2] -= flux;
        }
        
        //Heat transfer through particle-surface radiation in particle uncontacting scenario
        if(RADIATION) //Fenglei QI, 2016/2/21
        {
          double emissivityi = emissivity_[type[i]-1];
          double emissivityj = emissivity_[type[j]-1];
          double capheight = minrad - sqrt(minrad*minrad-contacting_pfp.upBound*contacting_pfp.upBound);
          double radiationarea = 2.0*M_PI*minrad*capheight;
          double viewfactor = viewfactor_correlation_uncontacted(contacting_pfp);
          double coefficent = 5.670367e-8*radiationarea/((1.0-emissivityi)/emissivityi+(1.0-emissivityj)/emissivityj+2.0/(1.0+viewfactor)); //Fenglei QI, 3/20/2016
          double tempi=Temp[i]*Temp[i];
          double tempj=Temp[j]*Temp[j];  
          double radheat = coefficent*(tempj*tempj - tempi*tempi);
          flux += radheat;

          heatFluxDistribution[i][3] += radheat;
          if (newton_pair || j < nlocal)
          {
            heatFluxDistribution[j][3] -= radheat;
          }
        }

        //debugging
         if (!(flux== flux) || (flux >1.e8) ) // NAN or infi
        {
            printf("Error with flux caculation \n");
            printf("Particle in %d node and with id %d \n", comm->me, atom->tag[i]);
            printf("pp: %.4e, pfp1: %.4e, pfp2: %.4e, radiation: %.4e\n",heatFluxDistribution[i][0],heatFluxDistribution[i][1],heatFluxDistribution[i][2],heatFluxDistribution[i][3]);
            printf("particle %d, x %.4e, y %.4e,z %.4e \n", atom->tag[i],x[i][0],x[i][1],x[i][1]);
            printf("particle %d, x %.4e, y %.4e,z %.4e \n", atom->tag[j],x[j][0],x[j][1],x[j][1]);
            printf(" hnc_pfp %.4e, tempi %.4e, tempj %.4e", hnc_pfp,Temp[i],Temp[j]);
            printf("\n");
        }


        dirFlux[0] = flux*delx;
        dirFlux[1] = flux*dely;
        dirFlux[2] = flux*delz;
        if(!cpl_flag)
        {
          //Add half of the flux (located at the contact) to each particle in contact
          heatFlux[i] += flux;
          directionalHeatFlux[i][0] += 0.50 * dirFlux[0];
          directionalHeatFlux[i][1] += 0.50 * dirFlux[1];
          directionalHeatFlux[i][2] += 0.50 * dirFlux[2];
          if (newton_pair || j < nlocal)
          {
            heatFlux[j] -= flux;
            directionalHeatFlux[j][0] += 0.50 * dirFlux[0];
            directionalHeatFlux[j][1] += 0.50 * dirFlux[1];
            directionalHeatFlux[j][2] += 0.50 * dirFlux[2];
          }

        }

        if(cpl_flag && cpl) cpl->add_heat(i,j,flux);
      }
    }
  }

  if(newton_pair) fix_heatFluxDistribution->do_reverse_comm(); //adding reverse comm, Fenglei
  if(newton_pair) fix_heatFlux->do_reverse_comm();
  if(newton_pair) fix_directionalHeatFlux->do_reverse_comm();
}


/* ---------------------------------------------------------------------- */
#if 0
template <int HISTFLAG,int CONTACTAREA>
void FixHeatGranCond::post_force_eval(int vflag,int cpl_flag)
{
  double hc,contactArea,delta_n,flux,dirFlux[3];
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double radi,radj,radsum,rsq,r,tcoi,tcoj;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *contact_flag,**first_contact_flag;

  int newton_pair = force->newton_pair;

  if (strcmp(force->pair_style,"hybrid")==0)
    error->warning(FLERR,"Fix heat/gran/conduction implementation may not be valid for pair style hybrid");
  if (strcmp(force->pair_style,"hybrid/overlay")==0)
    error->warning(FLERR,"Fix heat/gran/conduction implementation may not be valid for pair style hybrid/overlay");

  inum = pair_gran->list->inum;
  ilist = pair_gran->list->ilist;
  numneigh = pair_gran->list->numneigh;
  firstneigh = pair_gran->list->firstneigh;
  if(HISTFLAG) first_contact_flag = pair_gran->listgranhistory->firstneigh;

  double *radius = atom->radius;
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  updatePtrs(); //post force update ptrs, it might be changed during  initial integration

  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    if(HISTFLAG) contact_flag = first_contact_flag[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      if (!(mask[i] & groupbit) && !(mask[j] & groupbit)) continue;

      if(!HISTFLAG)
      {
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        radj = radius[j];
        radsum = radi + radj;
      }

      if ((HISTFLAG && contact_flag[jj]) || (!HISTFLAG && (rsq < radsum*radsum))) {  //contact
        
        if(HISTFLAG)
        {
          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx*delx + dely*dely + delz*delz;
          radj = radius[j];
          radsum = radi + radj;
          if(rsq >= radsum*radsum) continue;
        }

        r = sqrt(rsq);

        if(CONTACTAREA == CONDUCTION_CONTACT_AREA_OVERLAP)
        {
            
            if(area_correction_flag_)
            {
              delta_n = radsum - r;
              delta_n *= deltan_ratio_[type[i]-1][type[j]-1];
              r = radsum - delta_n;
            }

            //contact area of the two spheres
            contactArea = - M_PI/4 * ( (r-radi-radj)*(r+radi-radj)*(r-radi+radj)*(r+radi+radj) )/(r*r);
        }
        else if (CONTACTAREA == CONDUCTION_CONTACT_AREA_CONSTANT)
            contactArea = fixed_contact_area_;
        else if (CONTACTAREA == CONDUCTION_CONTACT_AREA_PROJECTION)
        {
            double rmax = MathExtraLiggghts::max(radi,radj);
            contactArea = M_PI*rmax*rmax;
        }

        tcoi = conductivity_[type[i]-1];
        tcoj = conductivity_[type[j]-1];
        if (tcoi < SMALL || tcoj < SMALL) hc = 0.;
        else hc = 4.*tcoi*tcoj/(tcoi+tcoj)*sqrt(contactArea);

        flux = (Temp[j]-Temp[i])*hc;

        dirFlux[0] = flux*delx;
        dirFlux[1] = flux*dely;
        dirFlux[2] = flux*delz;
        if(!cpl_flag)
        {
          //Add half of the flux (located at the contact) to each particle in contact
          heatFlux[i] += flux;
          directionalHeatFlux[i][0] += 0.50 * dirFlux[0];
          directionalHeatFlux[i][1] += 0.50 * dirFlux[1];
          directionalHeatFlux[i][2] += 0.50 * dirFlux[2];
          if (newton_pair || j < nlocal)
          {
            heatFlux[j] -= flux;
            directionalHeatFlux[j][0] += 0.50 * dirFlux[0];
            directionalHeatFlux[j][1] += 0.50 * dirFlux[1];
            directionalHeatFlux[j][2] += 0.50 * dirFlux[2];
          }

        }

        if(cpl_flag && cpl) cpl->add_heat(i,j,flux);
      }
    }
  }

  if(newton_pair) fix_heatFlux->do_reverse_comm();
  if(newton_pair) fix_directionalHeatFlux->do_reverse_comm();
}

#endif
//End Original Code

/* ----------------------------------------------------------------------
   register and unregister callback to compute
------------------------------------------------------------------------- */

void FixHeatGranCond::register_compute_pair_local(ComputePairGranLocal *ptr)
{
   
   if(cpl != NULL)
      error->all(FLERR,"Fix heat/gran/conduction allows only one compute of type pair/local");
   cpl = ptr;
}

void FixHeatGranCond::unregister_compute_pair_local(ComputePairGranLocal *ptr)
{
   
   if(cpl != ptr)
       error->all(FLERR,"Illegal situation in FixHeatGranCond::unregister_compute_pair_local");
   cpl = NULL;
}

/* ----------------------------------------------------------------------
   Additional functions
------------------------------------------------------------------------- */
double FixHeatGranCond::pfp_heat_transfer(const FixHeatGranCond::pfp_param &pp, const double &rr)
{
  double rad = MathExtraLiggghts::min(pp.radi,pp.radj);
  double temp = sqrt(rad*rad-rr*rr); //Fenglei Qi, 2016/1/27
  double tmp = (temp-rr*(rad+pp.half_gap)/pp.rad_eff)*(1.0/pp.condi+ 1.0/pp.condj);
       tmp+=2.0*(rad+pp.half_gap-temp)/pp.condf;

  return 2.0*M_PI*rr/tmp;
}

double FixHeatGranCond::num_integrate(FixHeatGranCond::hfun fun,const FixHeatGranCond::pfp_param &pfp, int numIntervals)
{
  
  double intervals=(pfp.upBound-pfp.lowBound)/(double)numIntervals;
  double result=0.0,tmp;

    for (int i=0;i<numIntervals;i++)
    {
      tmp=(this->*fun)(pfp,pfp.lowBound+1.0*i*intervals)+(this->*fun)(pfp,pfp.lowBound+1.0*(i+1)*intervals);
      result+=0.5*tmp*intervals;
    }
    return result;
}

double FixHeatGranCond::viewfactor_correlation_contacted(const FixHeatGranCond::pfp_param &pfp)
{
  double rad = MathExtraLiggghts::min(pfp.radi,pfp.radj);
  double x = 2.0*(rad+pfp.half_gap)/rad;
  double y = 0.4225*x-0.3371;
  return y;
}
double FixHeatGranCond::viewfactor_correlation_uncontacted(const FixHeatGranCond::pfp_param &pfp)
{
  double rad = MathExtraLiggghts::min(pfp.radi,pfp.radj);
  double x = 2.0*(rad+pfp.half_gap)/rad;
  double temp = pow(0.6794*0.6794+0.25*x*x,7.364);
  double y = 0.06233*pow(x,7.051)/temp;
  return y;
}
