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

FixStyle(property/atom/reaction,FixPropertyAtomReaction)

#else

#ifndef LMP_FIX_PROPERTY_ATOM_REACTION_H
#define LMP_FIX_PROPERTY_ATOM_REACTION_H

#include <cstring>
#include "fix.h"
#include "input.h"
#include "fix_property_atom.h"
#include "error.h"

namespace LAMMPS_NS {

class FixPropertyAtomReaction : public FixPropertyAtom {
 friend class Set;
 friend class FixPropertyAtomUpdateFix;
 friend class FixPropertyAtomRandom;
 friend class ParticleToInsert;
 public:
  FixPropertyAtomReaction(class LAMMPS *, int, char **,bool parse = true);
  ~FixPropertyAtomReaction();
  virtual void init();
  virtual int setmask();
  virtual void initial_integrate(int vflag);
  virtual void pre_force(int vflag);
  virtual void set_arrays(int);
  virtual void restart(char *);
  virtual double compute_vector(int n);
  virtual void mark_tracers(int ilo, int ihi) { UNUSED(ilo); UNUSED(ihi); }
  virtual double output_default(int i);

  virtual Fix* check_fix(const char *varname,const char *svmstyle,int len1,int len2,const char *caller,bool errflag);

  void do_forward_comm();
  void do_reverse_comm();
  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int,int);
  void pre_set_arrays();

  void set_all(double value);

  void write_restart(FILE *);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();
  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
 protected:
  void parse_args(int narg, char **arg);
  void readinparameters();
  void parameter_allocation();
  void parameter_deallocation();
  inline double property_T(int, double);
  inline double species_T(int,double);
  inline double property_species(int, const double*);

  FixPropertyAtom* fix_temp;
  FixPropertyAtom* fix_composition;
  double *Temp; // access to fix_composition vector_atom
  double **composition; // access to fix_composition array_atom
  int nevery_; //update properties every this step
  bool restart_flag;
 private:
  char *variablename;   // name of the variable (used for identification by other fixes)
  int data_style;       // 0 if a scalar is registered, 1 if vector
  int commGhost;        // 1 if communicated to ghost particles (via pack_comm/unpack_comm), 0 if not
  int commGhostRev;     // 1 if rev communicated from ghost particles (via pack_comm_rev/unpack_comm_rev), 0 if not
  int nvalues;
  double *defaultvalues; // default values at particle creation

  int ntypes;
  int numComp;
  int numSolid;
  int numLiquid;
  int numGas;
  int ndefault;

  int* f_temp; //particle property dependency on temperature
  double** temp_coeff; // coefficents, only quadratic coefficients are supported
  int* f_species; //particle property dependency on species
  int* f_temp_s; //species dependency on temperature
  double** temp_s_coeff; //coefficents for species, only quadratic coefficents are supported
  
  bool infileflag;
  char infile[256];         // infile

  // in case of initialization from property - name of property
  char *propertyname;
  double *property;
}; //end class

#include "fix_property_atom_reaction_i.h"

}
#endif
#endif
