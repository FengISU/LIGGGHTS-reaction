/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   CFD-DEM Coupling Stuff
------------------------------------------------------------------------- */

#include <mpi.h>
#include <string.h>
#include <stdio.h>
#include "library_cfd_coupling.h"
#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "fix_cfd_coupling.h"
#include "fix_multisphere.h"
#include "cfd_regionmodel.h"
#include "memory.h"
#include "error.h"
#include "comm.h"
#include "variable.h"
#include "cfd_datacoupling.h"
#include "cfd_datacoupling_one2one.h"

using namespace LAMMPS_NS;

#define LMP_GROW_DELTA 11000
/*NL*/ #define LMP_OF_DEBUG false

/* ---------------------------------------------------------------------- */

int liggghts_get_maxtag(void *ptr)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  return lmp->atom->tag_max();
}

/* ---------------------------------------------------------------------- */

int liggghts_get_maxtag_ms(void *ptr)
{
  // currently no possibility to delete multisphere bodies
  // so just return # of bodies

  LAMMPS *lmp = (LAMMPS *) ptr;
  FixMultisphere *fix_ms = static_cast<FixMultisphere*>(lmp->modify->find_fix_style_strict("multisphere",0));
  if(!fix_ms) return 0;
  return fix_ms->tag_max_body();
}

/* ---------------------------------------------------------------------- */

int liggghts_get_ntypes_ms(void *ptr)
{
  // currently no possibility to delete multisphere bodies
  // so just return # of bodies

  LAMMPS *lmp = (LAMMPS *) ptr;
  FixMultisphere *fix_ms = static_cast<FixMultisphere*>(lmp->modify->find_fix_style_strict("multisphere",0));
  if(!fix_ms) return 0;
  return fix_ms->ntypes();
}

/* ---------------------------------------------------------------------- */

double* liggghts_get_vclump_ms(void *ptr)
{
  // currently no possibility to delete multisphere bodies
  // so just return # of bodies

  LAMMPS *lmp = (LAMMPS *) ptr;
  FixMultisphere *fix_ms = static_cast<FixMultisphere*>(lmp->modify->find_fix_style_strict("multisphere",0));
  if(!fix_ms) return 0;
  return fix_ms->vclump();
}

/* ---------------------------------------------------------------------- */

double liggghts_get_variable(void *ptr, const char *variablename)
{
    LAMMPS *lmp = (LAMMPS *) ptr;
    int size = strlen(variablename);
    char* vn = new char[size+1];
    strcpy(vn,variablename); // Fengei Qi, making compatible with LIGGGHTs321version
    int var = lmp->input->variable->find(vn); 
    delete [] vn;
    if (var < 0)
    {
      char error_message[128] = {0};
      snprintf(error_message, 127,"Failed to find DEM variable '%s' requested by OF, aborting.\n",variablename);
      lmp->error->all(FLERR,error_message);
    }
    return lmp->input->variable->compute_equal(var);
}

/* ---------------------------------------------------------------------- */

void* locate_coupling_fix(void *ptr)
{
    LAMMPS *lmp = (LAMMPS *) ptr;
    int ifix = -1;
    for(int i=0;i<lmp->modify->nfix;i++)
      if(strcmp(lmp->modify->fix[i]->style,"couple/cfd") == 0)
        ifix = i;

    if(ifix ==-1) lmp->error->all(FLERR,"No fix of style 'couple/cfd' found, aborting.");

    return ((void*)lmp->modify->fix[ifix]);
}

/* ---------------------------------------------------------------------- */

void data_liggghts_to_of(const char *name, const char *type, void *ptr, void *&data, const char *datatype)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->push(name,type,data,datatype);
}

/* ---------------------------------------------------------------------- */

void data_of_to_liggghts(const char *name,const char *type,void *ptr,void *data,const char* datatype)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->pull(name,type,data,datatype);
}

/* ---------------------------------------------------------------------- */

//NP update region model
void update_rm(void *ptr)
{
    LAMMPS *lmp = (LAMMPS *) ptr;
    //FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    locate_coupling_fix(ptr);
    //CfdRegionmodel *rm = fcfd->rm;

    //NP call region model
    //if(rm) rm->rm_update();
    lmp->error->all(FLERR,"Region model update not implemented aborting.");
}

/* ---------------------------------------------------------------------- */

void allocate_external_int(int    **&data, int len2,int len1,int    initvalue,void *ptr)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->allocate_external(data,len2,len1,initvalue);
}
/* ---------------------------------------------------------------------- */

void allocate_external_int(int    **&data, int len2,const char *keyword,int    initvalue,void *ptr)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->allocate_external(data,len2,keyword,initvalue);
}

/* ---------------------------------------------------------------------- */

void allocate_external_double(double **&data, int len2,int len1,double initvalue,void *ptr)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->allocate_external(data,len2,len1,initvalue);
}

/* ---------------------------------------------------------------------- */

void allocate_external_double(double **&data, int len2,const char* keyword,double initvalue,void *ptr)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->allocate_external(data,len2,keyword,initvalue);
}

/* ---------------------------------------------------------------------- */

//NP check if all requested quantities have been communicated
//NP    since last call of this function
void check_datatransfer(void *ptr)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->check_datatransfer();
}

/* ---------------------------------------------------------------------- */

double** o2o_liggghts_get_boundingbox(void *ptr)
{
    LAMMPS *lmp = (LAMMPS *) ptr;
    double** bbox = new double*[2];
    bbox[0] = lmp->domain->sublo;
    bbox[1] = lmp->domain->subhi;
    return bbox;
}

/* ---------------------------------------------------------------------- */

void o2o_data_of_to_liggghts
(
    const char *name,
    const char *type,
    void *ptr,
    void *data,
    const char* datatype,
    const int* ids,
    const int ncollected
)
{
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    CfdDatacouplingOne2One* dc = static_cast<CfdDatacouplingOne2One*>(fcfd->get_dc());

    if(strcmp(datatype,"double") == 0)
        dc->pull_mpi<double>(name, type, data, ids, ncollected);
    else if(strcmp(datatype,"int") == 0)
        dc->pull_mpi<int>(name, type, data, ids, ncollected);
 //   else error->one(FLERR,"Illegal call to CfdDatacouplingOne2One::pull, valid datatypes are 'int' and double'");

}

