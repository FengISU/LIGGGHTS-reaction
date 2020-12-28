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

    Christoph Kloss (DCS Computing GmbH, Linz, JKU Linz)
    Richard Berger (JKU Linz)
    Philippe Seil (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_wall_gran.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "pair_gran.h"
#include "fix_rigid.h"
#include "fix_mesh.h"
#include "fix_contact_history_mesh.h"
#include "modify.h"
#include "respa.h"
#include "memory.h"
#include "comm.h"
#include "error.h"
#include "fix_property_atom.h"
#include "fix_property_atom_reaction.h"
#include "fix_contact_property_atom_wall.h"
#include "math_extra.h"
#include "math_extra_liggghts.h"
#include "compute_pair_gran_local.h"
#include "fix_neighlist_mesh.h"
#include "fix_mesh_surface_stress.h"
#include "tri_mesh.h"
#include "primitive_wall.h"
#include "primitive_wall_definitions.h"
#include "mpi_liggghts.h"
#include "neighbor.h"
#include "contact_interface.h"
#include "fix_property_global.h"
#include <vector>
#include "granular_wall.h"
#ifdef SUPERQUADRIC_ACTIVE_FLAG
  #include "math_extra_liggghts_superquadric.h"
#endif

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace LAMMPS_NS::PRIMITIVE_WALL_DEFINITIONS;
using namespace LIGGGHTS::Walls;
using namespace LIGGGHTS::ContactModels;

const double SMALL = 1e-12;

//heat transfer parameters, added by Fenglei Qi
#define CONFIG_ENABLE_DEBUG_LOG 0
#define EPSILON 0.44
#define INTERVALS 8 //Fenglei QI,2016/01/16
#define DMAX_FACTOR 0.5 //sqrt(2)-1
#define RADIATION 1 //Fenglei QI, 2016/2/21
//#define EMISSIVITY 0.9 //Fenglei QI, 2016/2/21
#define LARGE_DISTANCE 100000 // corresponding to the LARGE_TRIMESH in tri_mesh_I.h
  // modes for conduction contact area calaculation
  // same as in fix_heat_gran_conduction.cpp

  enum{ CONDUCTION_CONTACT_AREA_OVERLAP,
        CONDUCTION_CONTACT_AREA_CONSTANT,
        CONDUCTION_CONTACT_AREA_PROJECTION};

/* ---------------------------------------------------------------------- */

FixWallGran::FixWallGran(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
    // wall/gran requires gran properties
    // sph not
    if (strncmp(style,"wall/gran",9) == 0 && (!atom->radius_flag || !atom->omega_flag || !atom->torque_flag))
        error->fix_error(FLERR,this,"requires atom attributes radius, omega, torque");

    // defaults
    store_force_ = false;
    store_force_contact_ = false;
    stress_flag_ = false;
    n_FixMesh_ = 0;
    dnum_ = 0;
    skinDistance_ = 0.0;

    r0_ = 0.;

    shear_ = 0;
    shearDim_ = shearAxis_ = -1;
    vectorZeroize3D(shearAxisVec_);

    atom_type_wall_ = 1; // will be overwritten during execution, but other fixes require a value here

    // initializations
    fix_wallforce_ = 0;
    fix_wallforce_contact_ = 0;
    fix_rigid_ = NULL;
    heattransfer_flag_ = NULL; //Modifed by Fenglei Qi, 2016/3/20

    FixMesh_list_ = NULL;
    primitiveWall_ = NULL;
    fix_history_primitive_ = NULL;

    rebuildPrimitiveNeighlist_ = false;

    addflag_ = 0;
    cwl_ = NULL;

    computeflag_ = 1;

    meshwall_ = -1;

    track_energy_ = false;

    Temp_wall = -1.;
    fixed_contact_area_ = 0.;
    Q = Q_add = 0.;

    area_calculation_mode_ = CONDUCTION_CONTACT_AREA_OVERLAP;

    // parse args
    //style = new char[strlen(arg[2])+2];
    //strcpy(style,arg[2]);

    iarg_ = 3;
    narg_ = narg;

    int nremaining = narg - 3;
    char ** remaining_args = &arg[3];

    int64_t variant = Factory::instance().selectVariant("gran", nremaining, remaining_args);
    impl = Factory::instance().create("gran", variant, lmp, this);
   
    if(!impl && 0 == strncmp(style,"wall/gran",9)) 
    {
        printf("ERROR: Detected problem with model '%s' (and possible subsequent arguments). \n", arg[4]);
        printf("ERROR: The user must ensure that the contact model for this fix exists in the list of variants. \n");
        printf("ERROR: For the list of variants see the 'style_tangential_model.h' and 'style_normal_model.h' file in the src directory of your LIGGGHTS installation. \n");
        error->fix_error(FLERR,this,"unknown contact model");
    }

    iarg_ = narg - nremaining;

    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
        hasargs = false;

        if (strcmp(arg[iarg_],"primitive") == 0) {
           iarg_++;
           meshwall_ = 0;

           if (meshwall_ == 1)
             error->fix_error(FLERR,this,"'mesh' and 'primitive' are incompatible, choose either of them");

           if (strcmp(arg[iarg_++],"type"))
             error->fix_error(FLERR,this,"expecting keyword 'type'");
           atom_type_wall_ = force->inumeric(FLERR,arg[iarg_++]);
           if (atom_type_wall_ < 1 || atom_type_wall_ > atom->ntypes)
             error->fix_error(FLERR,this,"1 <= type <= max type as defined in create_box'");

           char *wallstyle = arg[iarg_++];
           int nPrimitiveArgs = PRIMITIVE_WALL_DEFINITIONS::numArgsPrimitiveWall(wallstyle);
           
           if(narg-iarg_ < nPrimitiveArgs)
            error->fix_error(FLERR,this,"not enough arguments for primitive wall");

           double * argVec = new double[nPrimitiveArgs];
           for(int i=0;i<nPrimitiveArgs;i++)
           {
             
             argVec[i] = force->numeric(FLERR,arg[iarg_++]);
           }

           bool setflag = false;
           for(int w=0;w<(int)PRIMITIVE_WALL_DEFINITIONS::NUM_WTYPE;w++)
           {
             
             if(strcmp(wallstyle,PRIMITIVE_WALL_DEFINITIONS::wallString[w]) == 0)
             {
               primitiveWall_ = new PrimitiveWall(lmp,(PRIMITIVE_WALL_DEFINITIONS::WallType)w,nPrimitiveArgs,argVec);
               setflag = true;
               break;
             }
           }
           if(!setflag) error->fix_error(FLERR,this,"unknown primitive wall style");
           hasargs = true;
           delete[] argVec;
        } else if (strcmp(arg[iarg_],"mesh") == 0) {
           hasargs = true;
           meshwall_ = 1;
           iarg_ += 1;
        } else if (strcmp(arg[iarg_],"track_energy") == 0) {
           hasargs = true;
           track_energy_ = true;
           iarg_ += 1;
        } else if (strcmp(arg[iarg_],"store_force") == 0) {
           if (iarg_+2 > narg)
              error->fix_error(FLERR,this," not enough arguments");
           if (strcmp(arg[iarg_+1],"yes") == 0) store_force_ = true;
           else if (strcmp(arg[iarg_+1],"no") == 0) store_force_ = false;
           else error->fix_error(FLERR,this,"expecting 'yes' or 'no' after keyword 'store_force'");
           hasargs = true;
           iarg_ += 2;
        } else if (strcmp(arg[iarg_],"store_force_contact") == 0) {
           if (iarg_+2 > narg)
              error->fix_error(FLERR,this," not enough arguments");
           if (strcmp(arg[iarg_+1],"yes") == 0) store_force_contact_ = true;
           else if (strcmp(arg[iarg_+1],"no") == 0) store_force_contact_ = false;
           else error->fix_error(FLERR,this,"expecting 'yes' or 'no' after keyword 'store_force_contact_'");
           hasargs = true;
           iarg_ += 2;
        } else if (strcmp(arg[iarg_],"n_meshes") == 0) {
          if (meshwall_ != 1)
             error->fix_error(FLERR,this,"have to use keyword 'mesh' before using 'n_meshes'");
          if (iarg_+2 > narg)
             error->fix_error(FLERR,this,"not enough arguments");
          n_FixMesh_ = atoi(arg[iarg_+1]);
          if(n_FixMesh_ < 1)
              error->fix_error(FLERR,this,"'n_meshes' > 0 required");
          hasargs = true;
          iarg_ += 2;
        } else if (strcmp(arg[iarg_],"meshes") == 0) {
          if (meshwall_ != 1)
             error->fix_error(FLERR,this,"have to use keyword 'mesh' before using 'meshes'");
          if(n_FixMesh_ == 0)
              error->fix_error(FLERR,this,"have to define 'n_meshes' before 'meshes'");
          if (narg < iarg_+1+n_FixMesh_)
              error->fix_error(FLERR,this,"not enough arguments");

          FixMesh_list_ = new FixMeshSurface*[n_FixMesh_];
          for(int i = 1; i <= n_FixMesh_; i++)
          {
              int f_i = modify->find_fix(arg[iarg_+i]);
              if (f_i == -1)
                  error->fix_error(FLERR,this,"could not find fix mesh id you provided");
              if (strncmp(modify->fix[f_i]->style,"mesh/surface",12))
                  error->fix_error(FLERR,this,"the fix belonging to the id you provided is not of type mesh");
              FixMesh_list_[i-1] = static_cast<FixMeshSurface*>(modify->fix[f_i]);

              if(FixMesh_list_[i-1]->trackStress())
                stress_flag_ = true;
              
          }
          hasargs = true;
          iarg_ += 1+n_FixMesh_;
        } else if (strcmp(arg[iarg_],"shear") == 0) {
          if (iarg_+3 > narg)
            error->fix_error(FLERR,this,"not enough arguments for 'shear'");
          if(!primitiveWall_)
            error->fix_error(FLERR,this,"have to define primitive wall before 'shear'. For mehs walls, please use fix move/mesh");

          if (strcmp(arg[iarg_+1],"x") == 0) shearDim_ = 0;
          else if (strcmp(arg[iarg_+1],"y") == 0) shearDim_ = 1;
          else if (strcmp(arg[iarg_+1],"z") == 0) shearDim_ = 2;
          else error->fix_error(FLERR,this,"illegal 'shear' dim");
          vshear_ = force->numeric(FLERR,arg[iarg_+2]);
          shear_ = 1;

          // update axis for cylinder etc if needed
          if(shearDim_ != primitiveWall_->axis())
          {
            shearAxis_ = primitiveWall_->axis();
            shearAxisVec_[shearAxis_] = vshear_;
          }

          hasargs = true;
          iarg_ += 3;
        } else if (strcmp(arg[iarg_],"temperature") == 0) {
            if (iarg_+1 >= narg)
              error->fix_error(FLERR,this,"not enough arguments for 'temperature'");
            Temp_wall = force->numeric(FLERR,arg[iarg_+1]);
            hasargs = true;
            iarg_ += 2;
        } else if(strcmp(arg[iarg_],"contact_area") == 0) {

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
        }
    }

    if(impl)
      impl->settings(narg - iarg_, &arg[iarg_]);

    // error checks

    if(meshwall_ == -1 && primitiveWall_ == 0)
        error->fix_error(FLERR,this,"Need to use define style 'mesh' or 'primitive'");

    if(meshwall_ == 1 && !FixMesh_list_)
        error->fix_error(FLERR,this,"Need to provide the number and a list of meshes by using 'n_meshes' and 'meshes'");
}

/* ---------------------------------------------------------------------- */

void FixWallGran::post_create()
{
    if(strncmp(style,"wall/gran",9) != 0)
    {
      // case non-granular (sph)
      dnum_ = 0;
    }

    // register storage for wall force if required
    if(store_force_)
    {
          char *wallforce_name = new char[strlen(style)+1+6];
          strcpy(wallforce_name,"force_");
          strcat(wallforce_name,id);
          char **fixarg = new char*[11];
          fixarg[0] = wallforce_name;
          fixarg[1] = (char *) "all";
          fixarg[2] = (char *) "property/atom";
          fixarg[3] = wallforce_name;
          fixarg[4] = (char *) "vector";
          fixarg[5] = (char *) "no";    // restart
          fixarg[6] = (char *) "no";    // communicate ghost
          fixarg[7] = (char *) "no";    // communicate rev
          fixarg[8] = (char *) "0.";
          fixarg[9] = (char *) "0.";
          fixarg[10] = (char *) "0.";
          modify->add_fix(11,fixarg);
          fix_wallforce_ =
              static_cast<FixPropertyAtom*>(modify->find_fix_property(wallforce_name,"property/atom","vector",3,0,style));
          delete []fixarg;
          delete []wallforce_name;
   }

   if(store_force_contact_ && 0 == meshwall_)
   {
        char **fixarg = new char*[19];
        char fixid[200],ownid[200];
        sprintf(fixid,"contactforces_%s",id);
        sprintf(ownid,"%s",id);
        fixarg[0]=fixid;
        fixarg[1]=(char *) "all";
        fixarg[2]=(char *) "contactproperty/atom/wall";
        fixarg[3]=fixid;
        fixarg[4]=(char *) "6";
        fixarg[5]=(char *) "fx";
        fixarg[6]=(char *) "0";
        fixarg[7]=(char *) "fy";
        fixarg[8]=(char *) "0";
        fixarg[9]=(char *) "fz";
        fixarg[10]=(char *) "0";
        fixarg[11]=(char *) "tx";
        fixarg[12]=(char *) "0";
        fixarg[13]=(char *) "ty";
        fixarg[14]=(char *) "0";
        fixarg[15]=(char *) "tz";
        fixarg[16]=(char *) "0";
        fixarg[17]=(char *) "primitive";
        fixarg[18]=ownid;
        modify->add_fix(19,fixarg);
        fix_wallforce_contact_ = static_cast<FixContactPropertyAtomWall*>(modify->find_fix_id(fixid));
        delete []fixarg;
   }

   // create neighbor list for each mesh
   
   for(int i=0;i<n_FixMesh_;i++)
   {
       
       FixMesh_list_[i]->createWallNeighList(igroup);
       FixMesh_list_[i]->createContactHistory(dnum());

       if(store_force_contact_)
         FixMesh_list_[i]->createMeshforceContact();
   }

   // contact history for primitive wall
   if(meshwall_ == 0 && dnum_ > 0)
   {
          char *hist_name = new char[strlen(id)+1+10];
          strcpy(hist_name,"history_");
          strcat(hist_name,id);
          char **fixarg = new char*[8+dnum_];
          fixarg[0] = hist_name;
          fixarg[1] = (char *) "all";
          fixarg[2] = (char *) "property/atom";
          fixarg[3] = hist_name;
          fixarg[4] = (char *) "vector";
          fixarg[5] = (char *) "yes";    // restart
          fixarg[6] = (char *) "no";    // communicate ghost
          fixarg[7] = (char *) "no";    // communicate rev
          for(int i = 8; i < 8+dnum_; i++)
              fixarg[i] = (char *) "0.";
          modify->add_fix(8+dnum_,fixarg);
          fix_history_primitive_ =
              static_cast<FixPropertyAtom*>(modify->find_fix_property(hist_name,"property/atom","vector",dnum_,0,style));
          delete []fixarg;
          delete []hist_name;
   }
}

/* ---------------------------------------------------------------------- */

void FixWallGran::pre_delete(bool unfixflag)
{
    if(unfixflag && store_force_)
        modify->delete_fix(fix_wallforce_->id);
    if(unfixflag && fix_history_primitive_)
        modify->delete_fix(fix_history_primitive_->id);

    if(unfixflag && store_force_contact_)
        modify->delete_fix(fix_wallforce_contact_->id);

    if(unfixflag && cwl_)
       error->fix_error(FLERR,this,"need to uncompute the active compute wall/gran/local before unfixing the wall");

    if(unfixflag)
    {
       for(int i=0;i<n_FixMesh_;i++)
       {
           
           FixMesh_list_[i]->deleteWallNeighList();
           FixMesh_list_[i]->deleteContactHistory();
       }
    }
}

/* ---------------------------------------------------------------------- */

FixWallGran::~FixWallGran()
{
    if(primitiveWall_ != 0) delete primitiveWall_;
    if(FixMesh_list_) delete [] FixMesh_list_;
    if(heattransfer_flag_) delete [] heattransfer_flag_;
    delete impl;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::setmask()
{
    int mask = 0;
    mask |= PRE_NEIGHBOR;
    mask |= PRE_FORCE;
    mask |= POST_FORCE;
    mask |= POST_FORCE_RESPA;
    return mask;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::min_type()
{
    return atom_type_wall_;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::max_type()
{
    return atom_type_wall_;
}

/* ---------------------------------------------------------------------- */

PrimitiveWall* FixWallGran::primitiveWall()
{ return primitiveWall_; }

/* ---------------------------------------------------------------------- */

void FixWallGran::init()
{
    dt_ = update->dt;

    // case granular
    if(strncmp(style,"wall/gran",9) == 0)
    {
        // check if a fix rigid is registered - important for damp
        fix_rigid_ = static_cast<FixRigid*>(modify->find_fix_style_strict("rigid",0));

        if (strcmp(update->integrate_style,"respa") == 0)
          nlevels_respa_ = ((Respa *) update->integrate)->nlevels;

        if(impl)
          impl->init_granular();
        else
        {
          // init for derived classes
          init_granular();
        }

        // disallow more than one wall of non-rimitive style
        
        if(is_mesh_wall())
        {
            int nfix = modify->n_fixes_style("wall/gran");
            for (int ifix = 0; ifix < nfix; ifix++)
            {
                FixWallGran *fwg = static_cast<FixWallGran*>(modify->find_fix_style("wall/gran",ifix));
                if (fwg == this) continue;
                if (fwg->is_mesh_wall())
                    error->fix_error(FLERR,this,"More than one wall of type 'mesh' is not supported");
            }
        }
    }
    
}

/* ---------------------------------------------------------------------- */

void FixWallGran::setup(int vflag)
{
    init_heattransfer(); //modifed by Fenglei QI, 2o16/3/20
    if (strstr(update->integrate_style,"verlet"))
    {
      pre_neighbor();
      pre_force(vflag);
      post_force(vflag);
    }
    else
    {
      ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa_-1);
      post_force_respa(vflag,nlevels_respa_-1,0);
      ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa_-1);
    }

    //init_heattransfer();
}

/* ----------------------------------------------------------------------
   neighbor list via fix wall/gran is only relevant for primitive walls
------------------------------------------------------------------------- */

void FixWallGran::pre_neighbor()
{
    rebuildPrimitiveNeighlist_ = (primitiveWall_ != 0);
}

void FixWallGran::pre_force(int vflag)
{
    double halfskin = neighbor->skin*0.5;
    int nlocal = atom->nlocal;

    x_ = atom->x;
    radius_ = atom->radius;
    cutneighmax_ = neighbor->cutneighmax;
#ifdef SUPERQUADRIC_ACTIVE_FLAG
    quat_ = atom->quaternion;
    shape_ = atom->shape;
    roundness_ = atom->roundness;
#endif

    // build neighlist for primitive walls
    
    if(rebuildPrimitiveNeighlist_)
      primitiveWall_->buildNeighList(radius_ ? halfskin:(r0_+halfskin),x_,radius_,nlocal);

    rebuildPrimitiveNeighlist_ = false;
}

/* ----------------------------------------------------------------------
   force on each atom calculated via post_force
   called via verlet
------------------------------------------------------------------------- */

void FixWallGran::post_force(int vflag)
{
    computeflag_ = 1;
    shearupdate_ = 1;
    if (update->setupflag) shearupdate_ = 0;
    addflag_ = 0;

    post_force_wall(vflag);
}

/* ----------------------------------------------------------------------
   force on each atom calculated via post_force
   called via compute wall/gran
------------------------------------------------------------------------- */

void FixWallGran::post_force_pgl()
{
    computeflag_ = 0;
    shearupdate_ = 0;
    addflag_ = 1;

    post_force_wall(0);
}

/* ----------------------------------------------------------------------
   post_force
------------------------------------------------------------------------- */

void FixWallGran::post_force_wall(int vflag)
{
  // set pointers and values appropriately
  nlocal_ = atom->nlocal;
  x_ = atom->x;
  f_ = atom->f;
  radius_ = atom->radius;
  rmass_ = atom->rmass;

#ifdef SUPERQUADRIC_ACTIVE_FLAG
  quat_ = atom->quaternion;
  shape_ = atom->shape;
  roundness_ = atom->roundness;
#endif

  if(fix_rigid_)
  {
      body_ = fix_rigid_->body;
      masstotal_ = fix_rigid_->masstotal;
  }

  if(fix_wallforce_)
    wallforce_ = fix_wallforce_->array_atom;

  cutneighmax_ = neighbor->cutneighmax;

  if(nlocal_ && !radius_ && r0_ == 0.)
    error->fix_error(FLERR,this,"need either per-atom radius or r0_ being set");

    if(store_force_)
    {
        for(int i = 0; i < nlocal_; i++)
        {
            vectorZeroize3D(wallforce_[i]);
        }
    }

  if(meshwall_ == 1)
    post_force_mesh(vflag);
  else
    post_force_primitive(vflag);

  if(meshwall_ == 0 && store_force_contact_)
    fix_wallforce_contact_->do_forward_comm();

  if(meshwall_ == 1 && store_force_contact_)
  {
    for(int imesh = 0; imesh < n_FixMesh_; imesh++)
        FixMesh_list_[imesh]->meshforceContact()->do_forward_comm();
  }
}

/* ---------------------------------------------------------------------- */

void FixWallGran::post_force_respa(int vflag, int ilevel, int iloop)
{
    if (ilevel == nlevels_respa_-1) post_force(vflag);
}

/* ----------------------------------------------------------------------
   post_force for mesh wall
------------------------------------------------------------------------- */

void FixWallGran::post_force_mesh(int vflag)
{
    
    // contact properties
    double v_wall[3],bary[3];
    double delta[3],deltan;
    MultiVectorContainer<double,3,3> *vMeshC;
    double ***vMesh;
    int nlocal = atom->nlocal;
    int nTriAll;

    SurfacesIntersectData sidata;
    sidata.is_wall = true;

    for(int iMesh = 0; iMesh < n_FixMesh_; iMesh++)
    {
      TriMesh *mesh = FixMesh_list_[iMesh]->triMesh();
      nTriAll = mesh->sizeLocal() + mesh->sizeGhost();
      FixContactHistoryMesh *fix_contact = FixMesh_list_[iMesh]->contactHistory();

      // mark all contacts for delettion at this point
      
      if(fix_contact) fix_contact->markAllContacts();

      if(store_force_contact_)
        fix_wallforce_contact_ = FixMesh_list_[iMesh]->meshforceContact();

      // get neighborList and numNeigh
      FixNeighlistMesh * meshNeighlist = FixMesh_list_[iMesh]->meshNeighlist();

      vectorZeroize3D(v_wall);
      vMeshC = mesh->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v");

      atom_type_wall_ = FixMesh_list_[iMesh]->atomTypeWall();

      // moving mesh
      if(vMeshC)
      {
        vMesh = vMeshC->begin();

        // loop owned and ghost triangles
        for(int iTri = 0; iTri < nTriAll; iTri++)
        {
          const std::vector<int> & neighborList = meshNeighlist->get_contact_list(iTri);
          const int numneigh = neighborList.size();
          for(int iCont = 0; iCont < numneigh; iCont++)
          {
            
            const int iPart = neighborList[iCont];

            // do not need to handle ghost particles
            if(iPart >= nlocal) continue;

            int idTri = mesh->id(iTri);

            #ifdef SUPERQUADRIC_ACTIVE_FLAG
                if(atom->superquadric_flag) {
                  sidata.pos_i = x_[iPart];
                  sidata.quat_i = quat_[iPart];
                  sidata.shape_i = shape_[iPart];
                  sidata.roundness_i = roundness_[iPart];
                  Superquadric particle(sidata.pos_i, sidata.quat_i, sidata.shape_i, sidata.roundness_i);
                  if(mesh->sphereTriangleIntersection(iTri, sidata.pos_i, radius_[iPart])) //check for Bounding Sphere-triangle intersection
                    deltan = mesh->resolveTriSuperquadricContactBary(iTri, radius_ ? radius_[iPart]:r0_, sidata.pos_i, delta, sidata.contact_point, particle, bary);
                  else
                    deltan = LARGE_TRIMESH;
                  sidata.is_non_spherical = true; //by default it is false
                } else
                  deltan = mesh->resolveTriSphereContactBary(iPart,iTri,radius_ ? radius_[iPart]:r0_ ,x_[iPart],delta,bary);
            #else
                deltan = mesh->resolveTriSphereContactBary(iPart,iTri,radius_ ? radius_[iPart]:r0_ ,x_[iPart],delta,bary);
            #endif
            
            if(deltan > skinDistance_) //allow force calculation away from the wall
            {
                if (deltan > LARGE_DISTANCE) continue; // Indicating the edge or corner in trimesh is not active,in this case delta is not caculated
                if(heattransfer_flag_[iMesh]) //heat transfer in noncontacting scenario, added by Fenglei Qi
                {
                  SurfacesCloseData scdata;
                  scdata.is_wall = true;
                  scdata.i = iPart;
                  scdata.j = iTri;
                  scdata.radi = radius_[iPart];
                  scdata.radsum = radius_[iPart];
                  scdata.delta[0] = delta[0];
                  scdata.delta[1] = delta[1];
                  scdata.delta[2] = delta[2];
                  scdata.rsq =  delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
                  #if 0 //CONFIG_ENABLE_DEBUG_LOG
                  if (update->ntimestep %1000 ==0)
                  {
                    char message[500];
                    sprintf(message,"Particles interaction with wall: dist between particle %d and wall: pradius: %.4e;distance: %.4e; gap: %.4e\n",
                    	    atom->tag[iPart],radius_[iPart],sqrt(scdata.rsq),deltan);
                    error->message(FLERR,message,1);
                  }
                  #endif
                  //here add a error remind if sqrt(rsq)!= radius_[iPart]+deltan, Fenglei Qi 2015/12/18
                  addNoncontactingHeatFlux(mesh,iPart,& scdata,1.);
                }
            }
            else
            {
              
              if(fix_contact && ! fix_contact->handleContact(iPart,idTri,sidata.contact_history)) continue;

              for(int i = 0; i < 3; i++)
                v_wall[i] = (bary[0]*vMesh[iTri][0][i] + bary[1]*vMesh[iTri][1][i] + bary[2]*vMesh[iTri][2][i]);

              sidata.i = iPart;
              sidata.is_wall = true; //added by Fenglei QI, March 4 2016 
              sidata.j = iTri; //added by Fenglei QI, March 4 2016
              sidata.deltan = -deltan;
              sidata.delta[0] = -delta[0];
              sidata.delta[1] = -delta[1];
              sidata.delta[2] = -delta[2];
              post_force_eval_contact(sidata, v_wall,iMesh,FixMesh_list_[iMesh],mesh,iTri);
            }

          }
        }
      }
      // non-moving mesh - do not calculate v_wall, use standard distance function
      else
      {
        // loop owned and ghost particles
        for(int iTri = 0; iTri < nTriAll; iTri++)
        {
          const std::vector<int> & neighborList = meshNeighlist->get_contact_list(iTri);
          const int numneigh = neighborList.size();
          
          for(int iCont = 0; iCont < numneigh; iCont++)
          {
            const int iPart = neighborList[iCont];

            // do not need to handle ghost particles
            if(iPart >= nlocal) continue;

            int idTri = mesh->id(iTri);
            #ifdef SUPERQUADRIC_ACTIVE_FLAG
                if(atom->superquadric_flag) {
                  sidata.pos_i = x_[iPart];
                  sidata.quat_i = quat_[iPart];
                  sidata.shape_i = shape_[iPart];
                  sidata.roundness_i = roundness_[iPart];
                  Superquadric particle(sidata.pos_i, sidata.quat_i, sidata.shape_i, sidata.roundness_i);
                  if(mesh->sphereTriangleIntersection(iTri, sidata.pos_i, radius_[iPart])) //check for Bounding Sphere-triangle intersection
                    deltan = mesh->resolveTriSuperquadricContact(iTri, radius_ ? radius_[iPart]:r0_, sidata.pos_i, delta, sidata.contact_point, particle);
                  else
                    deltan = LARGE_TRIMESH;
                  sidata.is_non_spherical = true; //by default it is false
                } else
                  deltan = mesh->resolveTriSphereContact(iPart,iTri,radius_ ? radius_[iPart]:r0_,x_[iPart],delta);
            #else
                deltan = mesh->resolveTriSphereContact(iPart,iTri,radius_ ? radius_[iPart]:r0_,x_[iPart],delta);
            #endif
            
            if(deltan > skinDistance_) //allow force calculation away from the wall
            {
              if (deltan > LARGE_DISTANCE) continue; // Indicating the edge or corner in trimesh is not active, in this case delta is not caculated
              if(heattransfer_flag_[iMesh]) //heat transfer in noncontacting scenario, added by Fenglei Qi
              {
                SurfacesCloseData scdata;
                scdata.is_wall = true; 
                scdata.i = iPart;
                scdata.j = iTri; //added by Fenglei Qi, March 04, 2016
                scdata.radi = radius_[iPart];
                scdata.radsum = radius_[iPart];
                scdata.delta[0] = delta[0];
                scdata.delta[1] = delta[1];
                scdata.delta[2] = delta[2];
                scdata.rsq =  delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
                addNoncontactingHeatFlux(mesh,iPart,& scdata,1.);
              }
              
            }
             
            else
            {
              
              if(fix_contact && ! fix_contact->handleContact(iPart,idTri,sidata.contact_history)) continue;
              
              sidata.i = iPart;
              sidata.is_wall = true;//aaded by Fenglei Qi, March 04,2016
              sidata.j = iTri; //added by Fenglei Qi, March 04,2016
              sidata.deltan = -deltan;
              sidata.delta[0] = -delta[0];
              sidata.delta[1] = -delta[1];
              sidata.delta[2] = -delta[2];
              post_force_eval_contact(sidata, v_wall,iMesh,FixMesh_list_[iMesh],mesh,iTri);
            }
          }
        }
      }

      // clean-up contacts
      
      if(fix_contact) fix_contact->cleanUpContacts();
    }

}

/* ----------------------------------------------------------------------
   post_force for primitive wall
------------------------------------------------------------------------- */

void FixWallGran::post_force_primitive(int vflag)
{
  int *mask = atom->mask;

  SurfacesIntersectData sidata;
  sidata.is_wall = true;

  // contact properties
  double delta[3]={},deltan,rdist[3];
  double v_wall[] = {0.,0.,0.};
  double **c_history = 0;

  if(dnum() > 0)
    c_history = fix_history_primitive_->array_atom;

  // if shear, set velocity accordingly
  if (shear_) v_wall[shearDim_] = vshear_;

  // loop neighbor list
  int *neighborList;
  int nNeigh = primitiveWall_->getNeighbors(neighborList);

  for (int iCont = 0; iCont < nNeigh ; iCont++, neighborList++)
  {
    int iPart = *neighborList;

    if(!(mask[iPart] & groupbit)) continue;

    deltan = primitiveWall_->resolveContact(x_[iPart],radius_?radius_[iPart]:r0_,delta);

    if(deltan > skinDistance_) //allow force calculation away from the wall
    {
      if(c_history) vectorZeroizeN(c_history[iPart],dnum_);
      if(heattransfer_flag_[0]) //heat transfer in noncontacting scenario, added by Fenglei Qi
      {
        SurfacesCloseData scdata;
        scdata.is_wall = true;
        scdata.i = iPart;
        scdata.radi = radius_[iPart];
        scdata.radsum = radius_[iPart];
        scdata.delta[0] = delta[0];
        scdata.delta[1] = delta[1];
        scdata.delta[2] = delta[2];
        scdata.rsq =  delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
        #if 0 //CONFIG_ENABLE_DEBUG_LOG
          if (update->ntimestep %1000 ==0)
            {
              char message[500];
              sprintf(message,"Particles interaction with wall: dist between particle %d and wall: pradius: %.4e;distance: %.4e; gap: %.4e\n",
              atom->tag[iPart],radius_[iPart],sqrt(scdata.rsq),deltan);
              error->message(FLERR,message,1);
            }
        #endif
        addNoncontactingHeatFlux(NULL,iPart,& scdata,1.);
      }      
    }
    else
    {
      bool particle_wall_intersection = true; //always true for spheres
      #ifdef SUPERQUADRIC_ACTIVE_FLAG
          if(atom->superquadric_flag) {
            double sphere_contact_point[3];
            vectorAdd3D(x_[iPart], delta, sphere_contact_point);
            sidata.pos_i = x_[iPart];
            sidata.quat_i = quat_[iPart];
            sidata.shape_i = shape_[iPart];
            sidata.roundness_i = roundness_[iPart];
            double closestPoint[3], closestPointProjection[3], point_of_lowest_potential[3];
            Superquadric particle(sidata.pos_i, sidata.quat_i, sidata.shape_i, sidata.roundness_i);
            particle_wall_intersection = particle.plane_intersection(delta, sphere_contact_point, closestPoint, point_of_lowest_potential);
            deltan = -MathExtraLiggghtsSuperquadric::point_wall_projection(delta, sphere_contact_point, closestPoint, closestPointProjection);
            vectorSubtract3D(closestPoint, closestPointProjection, delta);
            vectorCopy3D(closestPoint, sidata.contact_point);
            sidata.is_non_spherical = true; //by default it is false
          }
      #endif
      if(particle_wall_intersection) { //always true for spheres
          if(shear_ && shearAxis_ >= 0)
          {
              primitiveWall_->calcRadialDistance(x_[iPart],rdist);
              vectorCross3D(shearAxisVec_,rdist,v_wall);
              
          }
          sidata.i = iPart;
          sidata.contact_history = c_history ? c_history[iPart] : NULL;
          sidata.deltan = -deltan;
          sidata.delta[0] = -delta[0];
          sidata.delta[1] = -delta[1];
          sidata.delta[2] = -delta[2];
          post_force_eval_contact(sidata,v_wall);
      } else {
        if(c_history) vectorZeroizeN(c_history[iPart],dnum_);
      }
    }
  }
}

void FixWallGran::compute_force(SurfacesIntersectData &, double *)
{
    
}

/* ----------------------------------------------------------------------
   actually calculate force, called for both mesh and primitive
------------------------------------------------------------------------- */

inline void FixWallGran::post_force_eval_contact(SurfacesIntersectData & sidata, double * v_wall, int iMesh, FixMeshSurface *fix_mesh, TriMesh *mesh, int iTri)
{
  const int iPart = sidata.i;

  // deltan > 0 in compute_force
  // but negative in distance algorithm
  sidata.r = (radius_ ? radius_[iPart] : r0_) - sidata.deltan; // sign of corrected, because negative value is passed
  sidata.rsq = sidata.r*sidata.r;
  sidata.meff = rmass_ ? rmass_[iPart] : atom->mass[atom->type[iPart]];
  sidata.area_ratio = 1.;

  sidata.computeflag = computeflag_;
  sidata.shearupdate = shearupdate_;
  sidata.jtype = atom_type_wall_;

  double force_old[3]={}, f_pw[3];

  // if force should be stored - remember old force
  if(store_force_ || stress_flag_)
    vectorCopy3D(f_[iPart],force_old);

  // add to cwl
  if(cwl_ && addflag_)
  {
      double contactPoint[3];
      vectorSubtract3D(x_[sidata.i],sidata.delta,contactPoint);
      cwl_->add_wall_1(iMesh,mesh->id(iTri),iPart,contactPoint,v_wall);
  }

  if(impl)
    impl->compute_force(this, sidata, v_wall,mesh,iTri);
  else
    compute_force(sidata, v_wall); // LEGACY CODE (SPH)

  // if force should be stored or evaluated
  if(store_force_ || stress_flag_)
  {
    vectorSubtract3D(f_[iPart],force_old,f_pw);

    if(store_force_)
        vectorAdd3D (wallforce_[iPart], f_pw, wallforce_[iPart]);

    if(stress_flag_ && fix_mesh->trackStress())
    {
        double delta[3];
        delta[0] = -sidata.delta[0];
        delta[1] = -sidata.delta[1];
        delta[2] = -sidata.delta[2];
        static_cast<FixMeshSurfaceStress*>(fix_mesh)->add_particle_contribution
        (
           iPart,f_pw,delta,iTri,v_wall
        );
    }
  }

  // add heat flux
  if(mesh && heattransfer_flag_[iMesh])
    addHeatFlux(mesh,iPart,sidata,1.);// modifed by Fenglei Qi, 2016.3.5
  else if(!mesh && heattransfer_flag_[0])
    addHeatFlux(mesh,iPart,sidata,1.);// modifed by Fenglei Qi, 2016.3.20
}

/* ---------------------------------------------------------------------- */

int FixWallGran::is_moving()
{
    if(is_mesh_wall())
    {
        for(int i = 0; i < n_FixMesh_; i++) {
            if(FixMesh_list_[i]->mesh()->isMoving())
               return 1;
        }
        return 0;
    }
    return shear_;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::n_contacts_local()
{
    if (!is_mesh_wall() || dnum() == 0) return 0;

    int ncontacts = 0;
    for(int i = 0; i < n_FixMesh_; i++)
        ncontacts += FixMesh_list_[i]->contactHistory()->n_contacts();

    return ncontacts;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::n_contacts_all()
{
    int ncontacts = n_contacts_local();
    MPI_Sum_Scalar(ncontacts,world);
    return ncontacts;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::n_contacts_local(int contact_groupbit)
{
    if (!is_mesh_wall() || dnum() == 0) return 0;

    int ncontacts = 0;
    for(int i = 0; i < n_FixMesh_; i++)
        ncontacts += FixMesh_list_[i]->contactHistory()->n_contacts(contact_groupbit);

    return ncontacts;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::n_contacts_all(int contact_groupbit)
{
    int ncontacts = n_contacts_local(contact_groupbit);
    MPI_Sum_Scalar(ncontacts,world);
    return ncontacts;
}

/* ----------------------------------------------------------------------
   register and unregister callback to compute
------------------------------------------------------------------------- */

void FixWallGran::register_compute_wall_local(ComputePairGranLocal *ptr,int &dnum_compute)
{
   if(cwl_ != NULL)
     error->fix_error(FLERR,this,"Fix wall/gran allows only one compute of type wall/gran/local");
   cwl_ = ptr;
   dnum_compute = dnum_; //history values
}

void FixWallGran::unregister_compute_wall_local(ComputePairGranLocal *ptr)
{
   if(cwl_ != ptr)
     error->fix_error(FLERR,this,"Illegal situation in FixWallGran::unregister_compute_wall_local");
   cwl_ = NULL;
}

/* ---------------------------------------------------------------------- */

void FixWallGran::init_heattransfer()
{
    fppa_T = NULL;
    fppa_hf = NULL;
    deltan_ratio = NULL;
    //added by Fenglei Qi, 2016/03/20
    if(primitiveWall_)
    {
      heattransfer_flag_ =new bool[1]; 
    }
    else if(is_mesh_wall())
    {
      heattransfer_flag_ =new bool[n_FixMesh_];
    }
    
    // decide if heat transfer is to be calculated
    
    if (!is_mesh_wall() && Temp_wall < 0.) { heattransfer_flag_[0]=false; return;}
    else if(!is_mesh_wall())
    {
      heattransfer_flag_[0] = true;
    }
    else if (is_mesh_wall())
    {
        int heatflag = 0;
        for(int imesh = 0; imesh < n_meshes(); imesh++)
        {
            heattransfer_flag_[imesh] =  mesh_list()[imesh]->mesh()->prop().getGlobalProperty<ScalarContainer<double> >("Temp") != NULL;
            heatflag = heatflag || heattransfer_flag_[imesh];
        }

        if(!heatflag) return;
    }

    // heat transfer is to be calculated - continue with initializations

    // set flag so addHeatFlux function is called
    //heattransfer_flag_ = true;

    // if(screen && comm->me == 0) fprintf(screen,"Initializing wall/gran heat transfer model\n");
    fppa_T = static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp","property/atom","scalar",1,0,style));
    fppa_hf = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatFlux","property/atom","scalar",1,0,style));
    fix_conductivity_ = 
    static_cast<FixPropertyAtomReaction*>
    (
      modify->find_fix_property("thermalConductivity","property/atom/reaction","scalar",1,0,style)
    );
    th_cond = fix_conductivity_->vector_atom;
    #if RADIATION
    emissivity = static_cast<FixPropertyGlobal*>(modify->find_fix_property("emissivity","property/global","peratomtype",0,0,style))->get_values();
    #endif

    //extra heat transfer parameters, added by Fenglei Qi
    fppa_hf_dist = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatFluxDistribution","property/atom","vector",5,0,style));
    th_condf = static_cast<FixPropertyGlobal*> (modify->find_fix_property("fluidThermalConductivity","property/global","scalar",1,0,style))->get_values();
    // if youngsModulusOriginal defined, get deltan_ratio
    Fix* ymo_fix = modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",0,0,style,false);
    // deltan_ratio is defined by heat transfer fix, see if there is one
    int n_htf = modify->n_fixes_style("heat/gran/conduction");

    // get deltan_ratio set by the heat transfer fix
    if(ymo_fix && n_htf) deltan_ratio = static_cast<FixPropertyGlobal*>(ymo_fix)->get_array_modified();
}

/* ---------------------------------------------------------------------- */

void FixWallGran::addHeatFlux(TriMesh *mesh,int ip, SurfacesIntersectData & sidata, double area_ratio)
{
    //r is the distance between the sphere center and wall
    double delta_n = sidata.deltan;
    double tcop, tcowall, hc, Acont=0.0, r, rcont, temprcont;//adding rcont by Fenglei Qi
    double reff_wall = atom->radius[ip];
    int itype = atom->type[ip];
    double ri = atom->radius[ip];
    double triArea, surfEffRad; //added by Fenglei QI, March 4, 2016 
    double temp;
    th_cond = fix_conductivity_->vector_atom;   //added by Fenglei Qi, Oct,4,2016
    if(meshwall_)
    {
      //char message[1000];
      //sprintf(message,"using mesh wall");
      //error->message(FLERR,message,1);
      triArea = mesh->areaElem(sidata.j);
      surfEffRad = sqrt(triArea/M_PI);
      //surfEffRad = ri; //temprary
      Temp_wall = (*mesh->prop().getGlobalProperty< ScalarContainer<double> >("Temp"))(0);
    }
    else
      surfEffRad = LARGE_DISTANCE;

//    if(mesh)
//        Temp_wall = (*mesh->prop().getGlobalProperty< ScalarContainer<double> >("Temp"))(0);

    double *Temp_p = fppa_T->vector_atom;
    double *heatflux = fppa_hf->vector_atom;
    double **heatflux_dist = fppa_hf_dist->array_atom; 

    if(CONDUCTION_CONTACT_AREA_OVERLAP == area_calculation_mode_)
    {
        r = ri - delta_n; //r is not being correctly if area correction is set up, comment by Fenglei QI, 2016,2,23
        rcont = sqrt(ri*ri - r*r);
        temprcont = rcont;
        if(deltan_ratio)
        {
          rcont *= deltan_ratio[itype-1][atom_type_wall_-1];
          delta_n = ri - sqrt(ri*ri - rcont*rcont);
        }
        
        #if 0
        if(deltan_ratio)
           delta_n *= deltan_ratio[itype-1][atom_type_wall_-1];

        r = ri - delta_n;

        Acont = (reff_wall*reff_wall-r*r)*M_PI*area_ratio; //contact area sphere-wall
        rcont = sqrt(ri*ri-r*r); //caculating contact surface radius, added by Fenglei Qi
        #endif
    }
    else if (CONDUCTION_CONTACT_AREA_CONSTANT == area_calculation_mode_)
    {
        Acont = fixed_contact_area_;
        rcont = sqrt(Acont / M_PI); // keep consistency, added by Fenglei Qi
    }
    else if (CONDUCTION_CONTACT_AREA_PROJECTION == area_calculation_mode_)
    {
        Acont = M_PI*ri*ri;
        rcont = sqrt(Acont /M_PI); // keep consistency, added by Fenglei Qi
    }

    tcop = th_cond[ip]; //types start at 1, array at 0
    tcowall = fix_conductivity_->output_default(atom_type_wall_);

    if ((fabs(tcop) < SMALL) || (fabs(tcowall) < SMALL)) hc = 0.;
    //else hc = 4.*tcop*tcowall/(tcop+tcowall)*sqrt(Acont);
    else hc = 4.*tcop*tcowall/(tcop+tcowall)*rcont;
    
    if(computeflag_)
    {
        temp = (Temp_wall-Temp_p[ip]) * hc;
        heatflux[ip] += temp;
        heatflux_dist[ip][0] += temp; //heat transfer through contacting surface, Fenglei Qi
        Q_add += temp * update->dt;
    }
    if(cwl_ && addflag_)
        cwl_->add_heat_wall(ip,(Temp_wall-Temp_p[ip]) * hc);
    
    #if 0
        if (update->ntimestep %1000 ==0)
        {
            char message[500];
            sprintf(message,"Parameters in heat transfer: particle: %d, radius: %.4e, contact radius, %.4e, deltan: %.4e\n",atom->tag[ip],ri,rcont,delta_n);
            error->message(FLERR,message,1);
        }
    #endif
    //heat transfer through particle-fluid-particle in contacting scenario, added by Fenglei Qi
    struct pfp_param wcontacting_pfp;

    wcontacting_pfp.radi = ri;
    wcontacting_pfp.condi = tcop;
    wcontacting_pfp.condw = tcowall;
    wcontacting_pfp.condf = th_condf[0];
    wcontacting_pfp.rad_eff = 0.560 * ri * pow(1.0-EPSILON,-1.0/3.0);
    wcontacting_pfp.gap = -delta_n;
    wcontacting_pfp.lowBound = rcont;
    double upB;
    upB = ri * wcontacting_pfp.rad_eff;//Fenglei QI, 2016/1/26, Modified, 2016/3/4
    upB /=sqrt(wcontacting_pfp.rad_eff*wcontacting_pfp.rad_eff+(ri+wcontacting_pfp.gap)*(ri+wcontacting_pfp.gap));
    wcontacting_pfp.upBound = MathExtraLiggghts::min(surfEffRad,upB);;
        
    double hc_pfp = num_integrate(&FixWallGran::pfp_heat_transfer,wcontacting_pfp,INTERVALS);
    if(computeflag_)
      {
        temp = (Temp_wall-Temp_p[ip])*hc_pfp;
        heatflux[ip] += temp; //separately storing??? 
        heatflux_dist[ip][1] += temp;
        Q_add += temp* update->dt;
      }
    #if CONFIG_ENABLE_DEBUG_LOG
        if (update->ntimestep %1000 ==0)
          {
            char message[500];
            sprintf(message,"Parameters in heat transfer: particle: %d, radius: %.6e, distant:%.6e,contact radius, %.6e (org: %.6e)",atom->tag[ip],ri,r,rcont,temprcont);
            error->message(FLERR,message,1);
            
            sprintf(message,"Particles contacting information: kp: %.4e, kw: %.4e, kg,: %.4e, scaled delta_n: %.6e,lowBound, %.4e, upbound,%.4e ", tcop,tcowall,th_condf[0],delta_n,rcont,wcontacting_pfp.upBound);
            error->message(FLERR,message,1);
            
            sprintf(message,"Heat transfer information:contacting pp: %.4e, contacting pfp: %.4e",hc,hc_pfp);  
            error->message(FLERR,message,1);

          }
    #endif
    //Heat transfer through particle-surface radiation in particle contacting scenario, Fenglei QI,2016/2/22
    if(RADIATION)
      { 
        double emissivityp = emissivity[itype-1]; //types start at 1, array at 0
        double emissivityw = emissivity[atom_type_wall_-1]; 
        double capheight = ri - sqrt(ri*ri-wcontacting_pfp.upBound*wcontacting_pfp.upBound);
        double radiationarea_sphere = 2.0*M_PI*ri*(capheight-delta_n);
        double radiationarea_wall = M_PI*(wcontacting_pfp.upBound+wcontacting_pfp.lowBound)*(wcontacting_pfp.upBound-wcontacting_pfp.lowBound);
        double viewfactor_sphere = viewfactor_correlation_contacted_sphere(wcontacting_pfp);
        double viewfactor_wall = viewfactor_sphere*radiationarea_sphere/radiationarea_wall;
        double resistence = 1.0/(1.0/radiationarea_sphere/(1.0-viewfactor_sphere)+1.0/radiationarea_wall/(1.0-viewfactor_wall));
        double coefficent = (1.0-emissivityp)/emissivityp/radiationarea_sphere+(1.0-emissivityw)/emissivityw/radiationarea_wall + 1.0/(radiationarea_sphere*viewfactor_sphere+resistence);
        double tempi=Temp_p[ip]*Temp_p[ip];
        double tempwall=Temp_wall*Temp_wall;
        double radheat = 5.6696e-8*(tempwall*tempwall - tempi*tempi)/coefficent;  
        #if 0
        if (update->ntimestep %10000 ==0)
          {
            char message[1000];
            sprintf(message,"Parameters in heat transfer: particle: %d, radius: %.6e,contact radius, %.6e (org: %.6e)",atom->tag[ip],ri,rcont,temprcont);
            error->message(FLERR,message,1);
            sprintf(message,"Particles contacting information: kp: %.4e, kw: %.4e, kg,: %.4e, scaled delta_n: %.6e,lowBound, %.4e, upbound,%.4e ", tcop,tcowall,th_condf[0],delta_n,rcont,wcontacting_pfp.upBound);
            error->message(FLERR,message,1);
            sprintf(message,"Parameters in radiation heat transfer: particle: %d, distant:%.6e, viewfactor: %.6e (wall:%.6e), radheat: %.6e\n",atom->tag[ip],ri+wcontacting_pfp.gap,viewfactor_sphere,viewfactor_wall,radheat);
            error->message(FLERR,message,1);
            sprintf(message,"cap height(%.6e),radiation surface area: sphere cap (%.6e), wall(%.6e),resistence (%.6e),coefficent(%.6e), Temp(wall,%.6e,pi:%.6e)\n",capheight,radiationarea_sphere,radiationarea_wall,resistence,coefficent,Temp_wall,Temp_p[ip]);
            error->message(FLERR,message,1);

          } 
        #endif

        if(computeflag_)
        {
          heatflux[ip] += radheat; //separately storing??? 
          heatflux_dist[ip][3] += radheat;
          Q_add += radheat * update->dt;
        }
      }

    //debugging
    if (!(heatflux[ip]== heatflux[ip]) || (heatflux[ip] >1.e8) ) // NAN or infi
      {
        printf("Error with wall flux caculation \n");
        printf("Particle in %d node and with id %d \n", comm->me, atom->tag[ip]);
        printf("pp: %.4e, pfp1: %.4e, pfp2: %.4e, radiation: %.4e\n",heatflux_dist[ip][0],heatflux_dist[ip][1],heatflux_dist[ip][2],heatflux_dist[ip][3]);
        printf("particle %d, x %.4e, y %.4e,z %.4e \n", atom->tag[ip], atom->x[ip][0], atom->x[ip][1], atom->x[ip][1]);
        printf(" hc %.4e, hc_pfp %.4e, tempi %.4e, tempj %.4e", hc, hc_pfp,Temp_p[ip],Temp_wall);
        printf("\n");
      }
    //if(cwl_ && addflag_)
    //  cwl_->add_heat_wall(ip,(Temp_wall-Temp_p[ip]) * hc);
    
}


/*-----------------------------------------------------------------------*/

void FixWallGran::addNoncontactingHeatFlux(TriMesh *mesh,int ip, SurfacesCloseData const *scdata, double area_ratio)
{ 
    //heat transfer through particle-fluid-particle in contacting scenario, added by Fenglei Qi
    double ri = atom->radius[ip];
    double delta_n = sqrt(scdata->rsq) - ri;
    if(delta_n/ri > DMAX_FACTOR ) return ;
    double triArea,center[3],surfNorm[3], surfEffRad; //added by Fenglei QI, March 4, 2016   
    th_cond = fix_conductivity_->vector_atom;   //added by Fenglei Qi, Oct,4,2016 
    if(meshwall_)
    {
      triArea = mesh->areaElem(scdata->j);
      mesh->surfaceNorm(scdata->j,surfNorm);
      mesh->center(scdata->j,center);      
      surfEffRad = sqrt(triArea/M_PI);
      //surfEffRad=ri;//temperary
    }
    else
    {
      surfEffRad = LARGE_DISTANCE;
    }

    
    struct pfp_param wcontacting_pfp;
    double tcop, tcowall;
    int itype = atom->type[ip];
    tcop = th_cond[ip]; //types start at 1, array at 0
    tcowall = fix_conductivity_->output_default(atom_type_wall_);

    wcontacting_pfp.radi = ri;
    wcontacting_pfp.condi = tcop;
    wcontacting_pfp.condw = tcowall;
    wcontacting_pfp.condf = th_condf[0];
    wcontacting_pfp.rad_eff = 0.560 * ri * pow(1.0-EPSILON,-1.0/3.0);
    wcontacting_pfp.gap = delta_n;
    wcontacting_pfp.lowBound = 0.;
    double upB;
    upB = ri * wcontacting_pfp.rad_eff;
    upB /= sqrt(wcontacting_pfp.rad_eff*wcontacting_pfp.rad_eff+(ri+delta_n)*(ri+delta_n)); ////Fenglei QI,2016/1/27 modifed 2016/3/4
    wcontacting_pfp.upBound = MathExtraLiggghts::min(surfEffRad,upB);

    if(mesh)
        Temp_wall = (*mesh->prop().getGlobalProperty< ScalarContainer<double> >("Temp"))(0);

    double *Temp_p = fppa_T->vector_atom;
    double *heatflux = fppa_hf->vector_atom;
    double **heatflux_dist = fppa_hf_dist->array_atom; 

    double hnc_pfp = num_integrate(&FixWallGran::pfp_heat_transfer,wcontacting_pfp,INTERVALS);
    if(computeflag_)
    {
        heatflux[ip] += (Temp_wall-Temp_p[ip]) * hnc_pfp;
        heatflux_dist[ip][2] += (Temp_wall-Temp_p[ip]) * hnc_pfp; //heat transfer through contacting surface, Fenglei Qi
        Q_add += (Temp_wall-Temp_p[ip]) * hnc_pfp * update->dt;
    }
    //heat transfer through radiation in uncontacting scenario, Fenglei QI, 2016/2/22
    if(RADIATION) 
    {
      double emissivityp = emissivity[itype-1]; //types start at 1, array at 0
      double emissivityw = emissivity[atom_type_wall_-1]; 
      double capheight = ri - sqrt(ri*ri-wcontacting_pfp.upBound*wcontacting_pfp.upBound);
      double radiationarea_sphere = 2.0*M_PI*ri*capheight;
      double radiationarea_wall = M_PI*wcontacting_pfp.upBound*wcontacting_pfp.upBound;
      double viewfactor_sphere = viewfactor_correlation_uncontacted_sphere(wcontacting_pfp);
      double viewfactor_wall = viewfactor_sphere*radiationarea_sphere/radiationarea_wall;
      double resistence = 1.0/(1.0/radiationarea_sphere/(1.0-viewfactor_sphere)+1.0/radiationarea_wall/(1.0-viewfactor_wall));
      double coefficent = (1.0-emissivityp)/emissivityp/radiationarea_sphere+(1.0-emissivityw)/emissivityw/radiationarea_wall + 1.0/(radiationarea_sphere*viewfactor_sphere+resistence);
      double tempi=Temp_p[ip]*Temp_p[ip];
      double tempwall=Temp_wall*Temp_wall;
      double radheat = 5.670367e-8*(tempwall*tempwall - tempi*tempi)/coefficent;   
      #if 0       
      if (update->ntimestep %10000 ==0)
          {
            char message[1000];
            sprintf(message,"Parameters in radiation heat transfer: particle: %d, distant:%.6e, viewfactor: %.6e (wall:%.6e), radheat: %.6e\n",atom->tag[ip],ri+wcontacting_pfp.gap,viewfactor_sphere,viewfactor_wall,radheat);
            error->message(FLERR,message,1);
          } 
      #endif
      if(computeflag_)
        {
          heatflux[ip] += radheat; //separately storing??? 
          heatflux_dist[ip][3] += radheat;
          Q_add += radheat * update->dt;
        }
    }

    //debugging
    if (!(heatflux[ip]== heatflux[ip]) || (heatflux[ip] >1.e8) ) // NAN or infi
      {
        printf("Error with wall flux caculation \n");
        printf("Particle in %d node and with id %d \n", comm->me, atom->tag[ip]);
        printf("pp: %.4e, pfp1: %.4e, pfp2: %.4e, radiation: %.4e\n",heatflux_dist[ip][0],heatflux_dist[ip][1],heatflux_dist[ip][2],heatflux_dist[ip][3]);
        printf("particle %d, x %.4e, y %.4e,z %.4e \n", atom->tag[ip], atom->x[ip][0], atom->x[ip][1], atom->x[ip][1]);
        printf("hnc_pfp %.4e, tempi %.4e, tempj %.4e", hnc_pfp,Temp_p[ip],Temp_wall);
        printf("\n");
      }
    //if(cwl_ && addflag_)
    //  cwl_->add_heat_wall(ip,(Temp_wall-Temp_p[ip]) * hc);
}

/* ----------------------------------------------------------------------
   Additional functions
------------------------------------------------------------------------- */
double FixWallGran::pfp_heat_transfer(const FixWallGran::pfp_param &pp, const double &rr)
{
  double rad = pp.radi;
  double temp = sqrt(rad*rad-rr*rr);//Fenglei QI,2016/1/27
  double tmp = (temp-rr*(rad+pp.gap)/pp.rad_eff)*(1.0/pp.condi);
  tmp+=(rad+pp.gap-temp)/pp.condf;

  return 2.0*M_PI*rr/tmp;
}

double FixWallGran::num_integrate(FixWallGran::hfun fun,const FixWallGran::pfp_param &pfp, int numIntervals)
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

double FixWallGran::viewfactor_correlation_contacted_sphere(const FixWallGran::pfp_param &pfp)
{
  double rad = pfp.radi;
  double gap = pfp.gap;
  double x = (rad+gap)/rad;
  double y =0.5795*x+0.1134;
  return y;

}
double FixWallGran::viewfactor_correlation_contacted_wall(const FixWallGran::pfp_param &pfp)
{
  double rad = pfp.radi;
  double gap = pfp.gap;
  double x = (rad+gap)/rad;
  double y = 0.0*x;
  return y;
}
double FixWallGran::viewfactor_correlation_uncontacted_sphere(const FixWallGran::pfp_param &pfp)
{
  double rad = pfp.radi;
  double gap = pfp.gap;
  double x = (rad+gap)/rad;
  double temp = pow(0.6794*0.6794+0.25*x*x,4.801);
  double y = 0.1346*pow(x,1.825)/temp;
  return y;
}
double FixWallGran::viewfactor_correlation_uncontacted_wall(const FixWallGran::pfp_param &pfp)
{
  double rad = pfp.radi;
  double gap = pfp.gap;
  double x = (rad+gap)/rad;
  double y = 0.0*x;
  return y;
}
