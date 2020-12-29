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
#include "stdio.h"
#include "unistd.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_property_atom.h"
#include "fix_property_atom_reaction.h"
#include "atom.h"
#include "memory.h"
#include "error.h"

#include "pair_gran.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "group.h"
#include "timer.h"
#include "neighbor.h"

#include "fix_reaction_particle0D.h"
#include "fix_reaction.h"
#include "cJSON.h"

#include "mpi_liggghts.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define EPSILON 0.001

/* ---------------------------------------------------------------------- */

FixPropertyAtomReaction::FixPropertyAtomReaction(LAMMPS *lmp, int narg, char **arg, bool parse) :
  FixPropertyAtom(lmp,narg,arg,false),
  propertyname(0),
  property(0),
  restart_flag(false),
  fix_temp(NULL),
  fix_composition(NULL),
  Temp(NULL),
  composition(NULL),
  f_temp(NULL),
  temp_coeff(NULL),
  f_species(NULL),
  f_temp_s(NULL),
  temp_s_coeff(NULL)
{
    
    if(parse) parse_args(narg,arg);
}

//parse the input commends
void FixPropertyAtomReaction::parse_args(int narg, char **arg)
{
    // Check args
    if (narg < 9) error->all(FLERR,"Illegal fix property/atom command, not enough arguments");
    if (narg > 29) error->warning(FLERR,"Vector length in fix property/atom larger than 20. Are you sure you want that?");

    // Read args
    
    int n = strlen(arg[3]) + 1;
    variablename = new char[n ];
    strcpy(variablename,arg[3]);

    //Only support scalar type
    if (strcmp(arg[4],"scalar") == 0) data_style = FIXPROPERTY_ATOM_SCALAR;
    else if (strcmp(arg[4],"vector") == 0) error->all(FLERR,"Only scalar type is support for defining particle property");
    else error->all(FLERR,"Unknown style for fix property/atom. Valid styles are 'scalar' or 'vector'");

    nvalues = 1;

    if (strcmp(arg[5],"yes") == 0)
    {
            restart_peratom = 1;
            restart_global = 1;
    }
    else if (strcmp(arg[5],"no") == 0)
    {
         restart_peratom = 0;
         restart_global = 0;
    }
    else error->all(FLERR,"Unknown restart style for fix property/atom. Valid styles are 'yes' or 'no'");

    if (strcmp(arg[6],"yes") == 0) commGhost = 1;
    else if (strcmp(arg[6],"no") == 0) commGhost = 0;
    else error->all(FLERR,"Unknown communicate_ghost style for fix property/atom. Valid styles are 'yes' or 'no'");

    if (strcmp(arg[7],"yes") == 0) commGhostRev = 1;
    else if (strcmp(arg[7],"no") == 0) commGhostRev = 0;
    else error->all(FLERR,"Unknown communicate_reverse_ghost style for fix property/atom. Valid styles are 'yes' or 'no'");
    
    infileflag = false;
    if(strcmp(arg[narg-2],"file") ==0)
    {
      strcpy(infile,arg[narg-1]);
      infileflag = true;
      ndefault = narg - 10;
    }
    else
      ndefault = narg-8;  

    ntypes = atom->ntypes;

    if (ndefault != ntypes)
      error->all(FLERR,"The number of default values is not equal to the number of atom types");

    defaultvalues = new double[ndefault];

    // fix handles properties that need to be initialized at particle creation
    create_attribute = 1;

    //Initilize the default array

    for (int j = 0; j < ndefault; j++)
      {
          if(strcmp(arg[8+j],"none") == 0)
          {
              create_attribute = 0;
              continue;
          }
          defaultvalues[j] = force->numeric(FLERR,arg[8+j]);
      }

    size_peratom_cols = 0;

    peratom_flag=1; 
    peratom_freq=1;
    extvector=0; 

    if (commGhost) comm_forward = nvalues;
    if (commGhostRev) comm_reverse = nvalues;

    // perform initial allocation of atom-based array
    // register with Atom class
    vector_atom = NULL; array_atom = NULL;
    grow_arrays(atom->nmax); 
    atom->add_callback(0); 
    if (restart_peratom) atom->add_callback(1); 

    // init all arrays since dump may access it on timestep 0
    // or a variable may access it before first run
    
    int nlocal = atom->nlocal;
    int* type = atom->type;
    if(create_attribute)
    {

      for (int i = 0; i < nlocal; i++)
      {
        vector_atom[i] = defaultvalues[type[i]-1];
      }

    }

    // check if there is already a fix that tries to register a property with the same name
    
    for (int ifix = 0; ifix < modify->nfix; ifix++)
        if ((modify->fix[ifix]) && (strcmp(modify->fix[ifix]->style,style) == 0) && (strcmp(((FixPropertyAtomReaction*)(modify->fix[ifix]))->variablename,variablename)==0) )
            error->fix_error(FLERR,this,"there is already a fix that registers a variable of the same name");

    // flags for vector output
    //vector_flag = 1;
    size_vector = nvalues;
    global_freq = 1;
    extvector = 1;
    
    ntypes = 0;
    numComp = 0;
    numSolid = 0;
    numLiquid = 0;
    numGas = 0;

    nevery_ = 1;
   //printf("file name:%s\n",infile);
    //printf("default values: %f, %f\n",defaultvalues[0],defaultvalues[1]);
}

/* ---------------------------------------------------------------------- */

void FixPropertyAtomReaction::parameter_allocation()
{
    ntypes = atom->ntypes;
    f_temp = new int[ntypes];
    f_species = new int[ntypes];
    temp_coeff = new double*[ntypes];
    for (int i = 0; i < ntypes; i++)
    {
      temp_coeff[i] = new double[3];
    }
    //printf("ntypes:%d\n",ntypes);    
    
    FixReaction* fix_reaction = NULL;
    
    //fix_reaction = static_cast<FixReaction*>(modify->find_fix_style("reaction", 0));
    for(int ifix = 0; ifix < modify->nfix; ifix++)
    {
      if(modify->fix[ifix] && strncmp(modify->fix[ifix]->style,"reaction",8) == 0)
      {
        fix_reaction = static_cast<FixReaction*> (modify->fix[ifix]);
        break;
      }
    }
    
    //printf("fix_style is %s\n",fix_reaction->style);

    if (fix_reaction)
    {
      numComp = fix_reaction->check_numComp();
      numSolid = fix_reaction->check_numSolid();
      numLiquid = fix_reaction->check_numLiquid();
      numGas = fix_reaction->check_numGas();
    }

    if (numComp > 0)
    {
      f_temp_s = new int[numComp];
      temp_s_coeff = new double*[numComp];
      for(int i = 0; i < numComp; i++)
      {
        temp_s_coeff[i] = new double[3];
      }
    }
    else
    {
      f_temp_s = NULL;
      temp_s_coeff = NULL;
    }

    //printf("numComp:%d\n",numComp);
}

/* ---------------------------------------------------------------------- 
Function reads in the property definition in a specified file. Two scenarios
are supported.
(1) The particle property is only function of temperature. In this case, 
0 value is given to the "composition". The "temperature" may or may 
not be defined in the file.
(2) The particle property is function of species (and temperature). In this 
case, the actual number of species is given to the "composition". The 
"temperature" must be defined in this case. 

For the relationship between property and temperature, only "constant" 
and "quadratic" are supported.

When inputing the factors for a quadratic expression, e.g. a0 +a1*T+a2*T^2,
the "factor" array is [a0, a1, a2].

Property file example for a binary mixture system:

{
  "particle 1":
  {
    "type": 1,
    "composition": 0,
    "temperature":
    {  
        "composition0":
      {
        "ID": 0,
        "type":"quadratic",
        "factor":[0.0,0.1,0.2]
      }
    }
  },
  "particle 2":
  {
    "type": 2,
    "composition": 3,
    "temperature":
    {
      "composition1":
      {
        "ID": 1,
        "type":"quadratic",
        "factor":[1.0,0.1,0.2]
      },
      "composition2":
      {
        "ID": 2,
        "type":"quadratic",
        "factor":[2.0,0.1,0.2]
      },
      "composition3":
      {
        "ID": 3,
        "type":"quadratic",
        "factor":[3.0,0.1,0.2]
      }
    }
  }
}

The composition ID should be  the same as defined in fix_reaction.
-----------------------------------------------------------------------*/

void FixPropertyAtomReaction::readinparameters()
{
   char cwd[2048];
   if (getcwd(cwd, sizeof(cwd)) != NULL)
    {
       strcat(cwd,"/");
       strcat(cwd,infile);
    }
   else
       error->all(FLERR,"Can't get property file path:");

   //parse the file
    FILE *f;
    int len;
    char *data;

    f=fopen(cwd,"rb");
    if(!f)
      error->all(FLERR,"Can't open property file");

    fseek(f,0,SEEK_END);
    len=ftell(f);
    fseek(f,0,SEEK_SET);
    data = new char[len+1];
    fread(data,1,len,f);
    fclose(f);
    cJSON *json = cJSON_Parse(data);
      
    char* out;
    out = cJSON_Print(json);
    if(comm->me == 0)
    {
      printf("The property definition in %s is:\n",infile);
      printf("%s\n",out);
    }
    free(out);
   
    //check on the number of particle type    
    //char *out = cJSON_Print(json);
    //printf("%s\n",out);
    //error->all(FLERR,"Done test");
    int numtype=0;
    if(json)
    {
      cJSON* particle=json->child;
      if(particle) numtype++;
      cJSON* temp=particle;
      while(temp->next)
      {
        temp=temp->next;
        numtype++;
      }
    }

    if (numtype != atom->ntypes)
       error->all(FLERR,"Input of particle properties is not equal to defined particle type");
    
    //seperate particle properites and obtain particle property objects
    cJSON** particles = new cJSON*[numtype];
    cJSON* temp = json->child;
    
    int i=0;
    particles[0] = temp;
    while(temp->next)
    {
      particles[++i] = temp->next;
      temp=temp->next;
    }
   
    //parsing properties for each particle--particle layer
    int type;
    int composition;
    cJSON* temperature;
    for(int j = 0; j < numtype; j++)
    {
      type = cJSON_GetObjectItem(particles[j],"type")->valueint;
      if(type < 1 || type > atom->ntypes)
        error->all(FLERR, "The particle type definition is wrong in the property file");
      
      composition = cJSON_GetObjectItem(particles[j],"composition")->valueint;
      if (composition < 0 || composition > numComp)
        error->all(FLERR, "The composition definition is wrong in the property file");
      else
        f_species[type-1] = composition;

      temperature = cJSON_GetObjectItem(particles[j],"temperature");

      if(temperature)
        f_temp[type-1] = 1;
      else
        f_temp[type-1] = 0;
    }

    i=0;
    for( int j = 0; j < numtype; j++)
    {
      if(f_species[j] > 0 ) i++;
    }

    if(i>1) error->all(FLERR,"Only one type particle is allowed to be reactive");

    //parsing properties for each particle--temperature child layer
    
    FixReaction* fix_reaction;
    for(int ifix = 0; ifix < modify->nfix; ifix++)
    {
      if(modify->fix[ifix] && strncmp(modify->fix[ifix]->style,"reaction",8) == 0)
      {
        fix_reaction = static_cast<FixReaction*> (modify->fix[ifix]);
        break;
      }
    }

    for(int j = 0; j < numtype; j++)
    {
      int type = cJSON_GetObjectItem(particles[j],"type")->valueint;
      
      if (f_temp[type-1] == 0 && f_species[type-1] > 0)
          error->all(FLERR, "Composition is specified but the dependency on the temperature is not defined");
      else if (f_temp[type-1] == 0 && f_species[type-1] == 0)
          continue;

      temperature = cJSON_GetObjectItem(particles[j],"temperature");
      cJSON* child = temperature->child;
      int sibling = 0;
      if(child) sibling++;
      else error->all(FLERR, "No relationship between property and tempture is defined");
      while(child->next)
      {
        sibling++;
        child = child->next;
      }
      if(f_species[type-1] == 0 && sibling > 1)
        error->all(FLERR, "Definition in composition and temperature does not match");
      else if (f_species[type-1] > 0 && sibling != numComp)
        error->all(FLERR, "Definition in composition and temperature does not match");

      if(f_species[type-1] > 0 && type != fix_reaction->check_reactive_particle())
        error->all(FLERR, "Defined reative particle does not match the definition in fix_reaction");
    }

    for(int j = 0; j < numtype; j++)
    {
        int type = cJSON_GetObjectItem(particles[j],"type")->valueint;
        
        if (f_temp[type-1] == 0 && f_species[type-1] > 0)
          error->all(FLERR, "Composition is specified but the dependency on the temperature is not defined");
        else if (f_temp[type-1] == 0 && f_species[type-1] == 0)
          continue;

        temperature = cJSON_GetObjectItem(particles[j],"temperature");
        int ID;
        char* relationtype;
        cJSON* factor; 
        cJSON* child = temperature->child;
        cJSON* temp;
        if(f_species[type-1] == 0 )
        {
          ID = cJSON_GetObjectItem(child,"ID")->valueint;
          relationtype = cJSON_GetObjectItem(child,"type")->valuestring;
          factor = cJSON_GetObjectItem(child,"factor");
          if(!strcmp(relationtype,"constant") && !strcmp(relationtype,"quadratic"))
            error->all(FLERR,"Wrong property-temperature relationship");
          temp = factor->child;
          for( int k = 0; k < 3; k++) //quadratic, index < 3
          {
            temp_coeff[type-1][k] = temp->valuedouble;
            temp = temp->next;
          }
        }
        else
        {
          for(int k = 0; k < numComp; k++)
          {
            ID = cJSON_GetObjectItem(child,"ID")->valueint;
            f_temp_s[ID-1] =1;
            relationtype = cJSON_GetObjectItem(child,"type")->valuestring;
            factor = cJSON_GetObjectItem(child,"factor");
            if(!strcmp(relationtype,"constant") && !strcmp(relationtype,"quadratic"))
              error->all(FLERR,"Wrong property-temperature relationship");
            temp = factor->child;
            for( int k = 0; k < 3; k++) //quadratic, index < 3
            {
              temp_s_coeff[ID-1][k] = temp->valuedouble;
              temp = temp->next;
            }
            child = child->next;
          }
        }
    }

    cJSON_Delete(json);

    delete [] data;
    delete [] particles;
}

/* ---------------------------------------------------------------------- */

void FixPropertyAtomReaction::parameter_deallocation()
{
    delete [] f_temp;
    delete [] f_species;

    for(int i = 0; i < ntypes; i++)
      delete [] temp_coeff[i];
    delete [] temp_coeff;
    
    if(numComp)
    {
      delete [] f_temp_s;
      for(int i = 0; i < numComp; i++)
        delete [] temp_s_coeff[i];
      delete [] temp_s_coeff;
    }
}
/* ---------------------------------------------------------------------- */

FixPropertyAtomReaction::~FixPropertyAtomReaction()
{
  // unregister callbacks to this fix from Atom class
  //atom->delete_callback(id,0);
  //if (restart_peratom) atom->delete_callback(id,1);

  // delete locally stored arrays
  delete [] variablename;
  delete [] defaultvalues;
  if(propertyname) delete [] propertyname;

  //if (data_style) memory->destroy(array_atom);
  //else memory->destroy(vector_atom);

  parameter_deallocation();
}

/* ---------------------------------------------------------------------- */

void FixPropertyAtomReaction::init()
{
    //allocation for reaction parameters
    parameter_allocation();
    //read in particle thermo properties
    if(infileflag)
      readinparameters();
    else
      return;
    
    bool temp_flag = false;
    bool reaction_flag = false;
    
    //check if we need to call fix_temp
    for(int i = 0 ; i < ntypes; i++ )
    {
      if (f_temp[i] > 0 && infileflag )
        {
          temp_flag = true;
          break;
        }
    }

    //check if we need to call fix_composition
    for(int i = 0 ; i < ntypes; i++ )
    {
      if (f_species[i] > 0 && infileflag )
        {
          reaction_flag = true;
          break;
        }
    }

    if (temp_flag)
    {
      fix_temp = static_cast<FixPropertyAtom*>
      (
        modify->find_fix_property("Temp","property/atom","scalar",0,0,style)
      );
      if(!fix_temp)
        error->all(FLERR, "Can't find fix_temperature property");
      Temp = fix_temp->vector_atom;
    }

    if(reaction_flag)
    {
      fix_composition = static_cast<FixPropertyAtom*> 
      (
        modify->find_fix_property("compositionFraction","property/atom","vector",numComp,0,this->style)
      );

      if(!fix_composition)
        error->all(FLERR,"Can't find fix_reaction");
      composition = fix_composition->array_atom;
      FixReactionParticle0D* fix_reaction = static_cast<FixReactionParticle0D*>(modify->find_fix_style("reaction", 0));
      nevery_ = fix_reaction->nevery_;
    }

    if (recent_restart)
    {
      restart_flag = true;
    }

}

/* ---------------------------------------------------------------------- */

int FixPropertyAtomReaction::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= INITIAL_INTEGRATE;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- 
Update the atom property nevery steps before doing heat transfer or 
reaction calculation. The nevery is 1 or nevery from reaction
-----------------------------------------------------------------------*/
void FixPropertyAtomReaction::initial_integrate(int vflag)
{
  if (restart_flag)
  {
    restart_flag = false;
    this->do_forward_comm();
  }
  else if(update->ntimestep % nevery_) 
    return;

  if (!infileflag) {
      int nlocal = atom->nlocal;
      int *type = atom->type;
      for (int i = 0; i < nlocal; i++)
      {
        vector_atom[i] = defaultvalues[type[i]-1];
      }
      this->do_forward_comm();
      return;
  }

  if(fix_temp)
    Temp = fix_temp->vector_atom;
  if(fix_composition)
    composition = fix_composition->array_atom;

  int nlocal = atom->nlocal;
  int *type = atom->type;
  int ptype;
  double *property_s;
  if (numComp) property_s = new double[numComp];
  
  for(int i = 0; i < nlocal; i++)
  {
    ptype = type[i]-1;
    if (f_temp[ptype] && f_species[ptype])
    {
      for (int j = 0; j < numComp; j++)
      {
        property_s[j] = species_T(j,Temp[i]);
      }
      vector_atom[i] = property_species(i, property_s);
    }
    else if (f_temp[ptype])
    {
        vector_atom[i]= property_T(ptype,Temp[i]);
    }
  }

  this->do_forward_comm();

  if (numComp) delete [] property_s;
}

void FixPropertyAtomReaction::pre_force(int vflag)
{
    if(neighbor->ago == 0)
      this->do_forward_comm();
}
/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixPropertyAtomReaction::set_arrays(int i)
{
  vector_atom[i] = property ? property[i] : defaultvalues[0];
}

/* ---------------------------------------------------------------------- */

Fix* FixPropertyAtomReaction::check_fix(const char *varname,const char *svmstyle,int len1,int len2,const char *caller,bool errflag)
{
    char errmsg[400];
    if(strcmp(varname,variablename) == 0)
    {
        if(strcmp(svmstyle,"scalar") == 0) len1 = 1;

        // check variable style
        if(
            (strcmp(svmstyle,"scalar") == 0 && data_style != FIXPROPERTY_ATOM_SCALAR) ||
            (strcmp(svmstyle,"vector") == 0 && data_style != FIXPROPERTY_ATOM_VECTOR)
        )
        {
            if(errflag)
            {
                sprintf(errmsg,"%s style required for fix property/atom variable %s for usage with caller %s",
                        svmstyle,varname,caller);
                error->all(FLERR,errmsg);
            }
            else return NULL;
        }

        // check length
        if(len1 > nvalues)
        {
            if(errflag)
            {
                sprintf(errmsg,"Fix property/atom variable %s has wrong length (length is %d but length %d expected) for usage with caller %s",
                        varname,nvalues,len1,caller);
                error->all(FLERR,errmsg);
            }
            else return NULL;
        }

        // success
        return static_cast<Fix*>(this);
    }

    return NULL;
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixPropertyAtomReaction::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;
  int nvalues_re;

  nvalues_re = static_cast<int> (list[n++]);
  if(nvalues_re != nvalues)
    error->fix_error(FLERR,this,"restarted fix has incompatible data size");
}

/* ----------------------------------------------------------------------
   return components of property sum on fix group, n = 0..nvalues-1
------------------------------------------------------------------------- */

double FixPropertyAtomReaction::compute_vector(int n)
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  double value = 0.;

  for (int i = 0; i < nlocal; i++)
  {
      if (mask[i] & groupbit)
      {
        value += vector_atom[i];
      }
  }

  MPI_Sum_Scalar(value,world);
  return value;
}


double FixPropertyAtomReaction::output_default(int i) //i is the particle type
{
  if(i < 1 || i > ndefault)
    error->fix_error(FLERR,this,"the index is out of bound");
  else 
    return defaultvalues[i-1]; 
}

/* ----------------------------------------------------------------------
   forward and backward comm to be used by other fixes as needed
------------------------------------------------------------------------- */

void FixPropertyAtomReaction::do_forward_comm()
{
    timer->stamp();
    if (commGhost) comm->forward_comm_fix(this);
    else error->all(FLERR,"FixPropertyAtom: Faulty implementation - forward_comm invoked, but not registered");
    timer->stamp(TIME_COMM);
}

void FixPropertyAtomReaction::do_reverse_comm()
{
   timer->stamp();
   if (commGhostRev)  comm->reverse_comm_fix(this);
   else error->all(FLERR,"FixPropertyAtom: Faulty implementation - reverse_comm invoked, but not registered");
   timer->stamp(TIME_COMM);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixPropertyAtomReaction::memory_usage()
{
    int nmax = atom->nmax;
    double bytes = nmax * nvalues * sizeof(double);
    return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixPropertyAtomReaction::grow_arrays(int nmax)
{
    if (data_style) memory->grow(array_atom,nmax,nvalues,"FixPropertyAtom:array_atom");
    else memory->grow(vector_atom, nmax, "FixPropertyAtom:vector_atom");
    
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixPropertyAtomReaction::copy_arrays(int i, int j, int delflag)
{
    if (data_style) for(int k=0;k<nvalues;k++) array_atom[j][k] = array_atom[i][k];
    else vector_atom[j]=vector_atom[i];
}

/* ----------------------------------------------------------------------
   called before set_arrays is called for each atom
------------------------------------------------------------------------- */

void FixPropertyAtomReaction::pre_set_arrays()
{
    
    property = 0;
    if(propertyname)
    {
        
        int len1,len2;
        
        property = (double*) atom->get_properties()->find_property(propertyname,"scalar-atom",len1,len2);
        if(!property)
        {
            char errstr[200];
            sprintf(errstr,"Property %s not found",propertyname);
            
            error->fix_error(FLERR,this,errstr);
        }

    }
}

/* ----------------------------------------------------------------------
   set all atoms values
------------------------------------------------------------------------- */

void FixPropertyAtomReaction::set_all(double value)
{
    
    int nlocal = atom->nlocal;
    if (data_style)
    {
        for(int i = 0; i < nlocal; i++)
        {
            for(int k=0;k<nvalues;k++)
                array_atom[i][k] = value;
        }
    }
    else
    {
        for(int i = 0; i < nlocal; i++)
            vector_atom[i] = value;
    }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixPropertyAtomReaction::pack_exchange(int i, double *buf)
{
    if (data_style) for(int k=0;k<nvalues;k++) buf[k] = array_atom[i][k];
    else buf[0] = vector_atom[i];
    return nvalues;
}

/* ----------------------------------------------------------------------
   unpack values into local atom-based arrays after exchange
------------------------------------------------------------------------- */

int FixPropertyAtomReaction::unpack_exchange(int nlocal, double *buf)
{
    if (data_style) for(int k=0;k<nvalues;k++) array_atom[nlocal][k] = buf[k];
    else vector_atom[nlocal]=buf[0];
    return nvalues;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixPropertyAtomReaction::pack_restart(int i, double *buf)
{
  buf[0] = static_cast<double>(nvalues+1);
  if (data_style) for(int k=0;k<nvalues;k++) buf[k+1] = array_atom[i][k];
  else buf[1] = vector_atom[i];

  return (nvalues+1);
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixPropertyAtomReaction::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  if (data_style) for(int k=0;k<nvalues;k++) array_atom[nlocal][k] = extra[nlocal][m++];
  else vector_atom[nlocal] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixPropertyAtomReaction::maxsize_restart()
{
  return nvalues+1;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixPropertyAtomReaction::size_restart(int nlocal)
{
  return nvalues+1;
}

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

int FixPropertyAtomReaction::pack_comm(int n, int *list, double *buf,
                             int pbc_flag, int *pbc)
{
    int i,j;
    //we dont need to account for pbc here
    int m = 0;
    for (i = 0; i < n; i++) {
      j = list[i];
      if (data_style) for(int k=0;k<nvalues;k++) buf[m++] = array_atom[j][k];
      else buf[m++] = vector_atom[j];
    }
    return nvalues;
}

/* ---------------------------------------------------------------------- */

void FixPropertyAtomReaction::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
      if (data_style) for(int k=0;k<nvalues;k++) array_atom[i][k]=buf[m++];
      else vector_atom[i]=buf[m++];
  }

}

/* ---------------------------------------------------------------------- */

int FixPropertyAtomReaction::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (data_style) for(int k=0;k<nvalues;k++) buf[m++] = array_atom[i][k];
    else buf[m++] = vector_atom[i];
  }
  return nvalues;
}

/* ---------------------------------------------------------------------- */

void FixPropertyAtomReaction::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (data_style) for(int k=0;k<nvalues;k++) array_atom[j][k]+=buf[m++];
    else vector_atom[j]+=buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void FixPropertyAtomReaction::write_restart(FILE *fp)
{
  int n = 0;
  double list[1];
  list[n++] = static_cast<double>(nvalues);

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}
