/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include "fix_conc_berendsen.h"
#include "fix_adapt_fep.h"
#include "atom.h"
#include "atom_vec.h"
#include "molecule.h"
#include "force.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "random_park.h"
#include "kspace.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "group.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "compute_conc.h"
#include "error.h"
#include "region.h"
#include "memory.h"
#include "neighbor.h"
#include "math_extra.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define NUM_TO_MOL 1661.129568 // Actually it is NUM/angstrom^3 -> mol/L

enum{NOBIAS,BIAS};
enum{CONSTANT,EQUAL};
enum{PAIR,KSPACE,ATOM,MOLECULE};
enum{CHARGE};

/* ---------------------------------------------------------------------- */

FixConcBerendsen::FixConcBerendsen(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  cstr(NULL), id_temp(NULL), cflag(0), random_equal(NULL), fixshake(NULL), idshake(NULL)
{
  if (narg < 10) error->all(FLERR,"Illegal fix temp/berendsen command");

  // Berendsen thermostat should be applied every step

  dynamic_group_allow = 1;
  nevery = 1;
  scalar_flag = 1;
  global_freq = nevery;
  extscalar = 1;
  shakeflag = 0;

  cstr = NULL;
  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    cstr = new char[n];
    strcpy(cstr,&arg[3][2]);
    cstyle = EQUAL;
  } else {
    c_start = force->numeric(FLERR,arg[3]);
    c_target = c_start;
    cstyle = CONSTANT;
  }

  c_stop = force->numeric(FLERR,arg[4]);
  c_period = force->numeric(FLERR,arg[5]);

  concfactor = force->inumeric(FLERR,arg[7]);//0:Using averaged concentration, 1:Using Cation concentration, 2:Using Anion concentration

  typeCat_num = 0;
  typeAni_num = 0;
  typeCat_frac_num = 0;
  typeAni_frac_num = 0;

  Cat_ratio = 1;
  Ani_ratio = 1;
  Cat_charge = 1.0;
  Ani_charge = -1.0;

  if (strcmp(arg[8],"atom")==0) {
    mode = ATOM;
    typeCat = force->inumeric(FLERR,arg[9]);//type of Cation
    typeAni = force->inumeric(FLERR,arg[10]);//type of Anion
    typeCat_frac = typeCat+2;//type of fractional Cation
    typeAni_frac = typeAni+2;//type of fractional Anion
    typeCat_num = 1;
    typeAni_num = 1;
    typeCat_frac_num = 1;
    typeAni_frac_num = 1;
    int iarg = 11;
    while (iarg < narg) {
	if (strcmp(arg[iarg],"ratio") == 0) {
	    Cat_ratio = force->inumeric(FLERR,arg[iarg+1]);
	    Ani_ratio = force->inumeric(FLERR,arg[iarg+2]);
	    iarg += 3;
	}
	else if (strcmp(arg[iarg],"charge") == 0) {
	    Cat_charge = force->numeric(FLERR,arg[iarg+1]);
	    Ani_charge = force->numeric(FLERR,arg[iarg+2]);
	    iarg += 3;
	}
	else error->all(FLERR,"Illegal fix conc/berendsen command");
    }
  }
  else if (strcmp(arg[8],"mol")==0) {
    mode = MOLECULE;
    int molbase = atom->find_molecule(arg[9]);
    if (molbase == -1) error->all(FLERR,"Molecule template ID for fix conc/berendsen does not exist");
    onemols = atom->molecules;
    if (onemols[molbase]->nset != 4) error->all(FLERR,"Fix conc/berendsen molecule template must contain exactly four molecules in cation, anion, cation_frac, and anion_frac order");
    molCat = molbase;
    molAni = molbase+1;
    molCat_frac = molbase+2;
    molAni_frac = molbase+3;

    // Assign represent atom type(first atom) & number
    // type* : min atom type, type*_max : max atom_type
    int typetemp;
    typeCat = onemols[molCat]->type[0];
    typeCat_max = typeCat;
    for ( int i = 0; i < onemols[molCat]->natoms; i++ ) {
      typetemp = onemols[molCat]->type[i];
      if (typetemp < typeCat) typeCat = typetemp;
      if (typetemp > typeCat_max) typeCat_max = typetemp;
    }
    for ( int i = 0; i < onemols[molCat]->natoms; i++ ) {
      if (typeCat == onemols[molCat]->type[i]) typeCat_num++;
    }
    typeAni = onemols[molAni]->type[0];
    typeAni_max = typeAni;
    for ( int i = 0; i < onemols[molAni]->natoms; i++ ) {
      typetemp = onemols[molAni]->type[i];
      if (typetemp < typeAni) typeAni = typetemp;
      if (typetemp > typeAni_max) typeAni_max = typetemp;
    }
    for ( int i = 0; i < onemols[molAni]->natoms; i++ ) {
      if (typeAni == onemols[molAni]->type[i]) typeAni_num++;
    }
    typeCat_frac = onemols[molCat_frac]->type[0];
    typeCat_frac_max = typeCat_frac;
    for ( int i = 0; i < onemols[molCat_frac]->natoms; i++ ) {
      typetemp = onemols[molCat_frac]->type[i];
      if (typetemp < typeCat_frac) typeCat_frac = typetemp;
      if (typetemp > typeCat_frac_max) typeCat_frac_max = typetemp;
    }
    for ( int i = 0; i < onemols[molCat_frac]->natoms; i++ ) {
      if (typeCat_frac == onemols[molCat_frac]->type[i]) typeCat_frac_num++;
    }
    typeAni_frac = onemols[molAni_frac]->type[0];
    typeAni_frac_max = typeAni_frac;
    for ( int i = 0; i < onemols[molAni_frac]->natoms; i++ ) {
      typetemp = onemols[molAni_frac]->type[i];
      if (typetemp < typeAni_frac) typeAni_frac = typetemp;
      if (typetemp > typeAni_frac_max) typeAni_frac_max = typetemp;
    }
    for ( int i = 0; i < onemols[molAni_frac]->natoms; i++ ) {
      if (typeAni_frac == onemols[molAni_frac]->type[i]) typeAni_frac_num++;
    }

    int iarg = 10;
    while (iarg < narg) {
	if (strcmp(arg[iarg],"ratio") == 0) {
	    Cat_ratio = force->inumeric(FLERR,arg[iarg+1]);
	    Ani_ratio = force->inumeric(FLERR,arg[iarg+2]);
	    iarg += 3;
	}
        else if (strcmp(arg[iarg],"shake") == 0) {
            if (iarg+2 > narg) error->all(FLERR,"Illegal fix conc/berendsen command");
            int n = strlen(arg[iarg+1]) + 1;
            delete [] idshake;
            idshake = new char[n];
            strcpy(idshake,arg[iarg+1]);
            shakeflag = 1;
            iarg += 2;
        }
	else error->all(FLERR,"Illegal fix conc/berendsen command");
    }
  }
  else error->all(FLERR,"Illegal fix conc/berendsen command"); 

  // error checks

  if (c_period <= 0.0)
    error->all(FLERR,"Fix temp/berendsen period must be > 0.0");

  if (Cat_ratio >= Ani_ratio) {
    if (Ani_ratio != 1) error->all(FLERR,"Cation Anion ratio should be 1:x or x:1");
  }
  else {
    if (Cat_ratio != 1) error->all(FLERR,"Cation Anion ratio should be 1:x or x:1");
  }

  // create a new compute conc style
  // id = fix-ID + conc, compute group = fix group

  int n = strlen(id) + 6;
  id_temp = new char[n];
  strcpy(id_temp,id);
  strcat(id_temp,"_conc");

  iregion = domain->find_region(arg[6]);
  int x = strlen(arg[6]) + 1;
  idregion = new char[x];
  strcpy(idregion,arg[6]);
  

  char **newarg = new char*[4];
  newarg[0] = id_temp;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "conc";
  newarg[3] = idregion;
  modify->add_compute(4,newarg);
  delete [] newarg;
  delete [] idregion;
  cflag = 1;

  energy = 0;

  // Start fix_conc_berendsen settings
  // adapt[].var is ignored
  adapt = new Adapt[3];//pair,kspace & atom or molecule
  nadapt = 3;
  //pair
  adapt[0].which = PAIR;
  n = strlen("lj/charmm/coul/long/soft") + 1;
  adapt[0].pstyle = new char[n];
  strcpy(adapt[0].pstyle,"lj/charmm/coul/long/soft");
  n = strlen("lambda") + 1;
  adapt[0].pparam = new char[n];
  strcpy(adapt[0].pparam,"lambda");
  if (mode == ATOM) {
    adapt[0].ilo = typeCat-2;//Water solvent assumed, mino
    adapt[0].ihi = typeAni;
    adapt[0].jlo = typeCat_frac;
    adapt[0].jhi = typeAni_frac;
  }
  if (mode == MOLECULE) {
    adapt[0].ilo = typeCat-2;//Water solvent assumed, mino
    adapt[0].ihi = typeAni_max;
    adapt[0].jlo = typeCat_frac;
    adapt[0].jhi = typeAni_frac_max;
  }
  //printf("%d\t%d\t%d\t%d\n",adapt[0].ilo,adapt[0].ihi,adapt[0].jlo,adapt[0].jhi);
  //kspace
  adapt[1].which = KSPACE;
  if (mode == ATOM) {
    //atom
    adapt[2].which = ATOM;
    adapt[2].aparam = CHARGE;
    chgflag = 1;
    adapt[2].ilo = typeCat_frac;
    adapt[2].ihi = typeAni_frac;
    // For the case where addtional atom types exist between typeAni_max & typeCat_frac
    adapt[2].jlo = typeCat;
    adapt[2].jhi = typeAni;
  }
  else if (mode == MOLECULE) {
    //molecule
    adapt[2].which = MOLECULE;
    adapt[2].aparam = CHARGE;
    chgflag = 1;
    adapt[2].ilo = typeCat_frac;
    adapt[2].ihi = typeAni_frac_max;
    // For the case where addtional atom types exist between typeAni_max & typeCat_frac
    adapt[2].jlo = typeCat;
    adapt[2].jhi = typeAni_max;
  }
  // optional keywords for fix_adapt_fep
  // reset=no, scale=no, after=no
  resetflag = 0;
  scaleflag = 0;
  afterflag = 1;

  // allocate pair style arrays
  n = atom->ntypes;
  for (int m = 0; m < 3; m++)
    if (adapt[m].which == PAIR)
      memory->create(adapt[m].array_orig,n+1,n+1,"adapt:array_orig");

  id_fix_diam = id_fix_chg = NULL;

  ngroups = 0;
  int ngroupsmax = 0;
  groupstrings = NULL;
  ngrouptypes = 0;
  int ngrouptypesmax = 0;
  grouptypestrings = NULL;
  grouptypes = NULL;
  grouptypebits = NULL;
  if (ngroups >= ngroupsmax) {
    ngroupsmax = ngroups+1;
    groupstrings = (char **)
      memory->srealloc(groupstrings,
                       ngroupsmax*sizeof(char *),
                       "fix_gcmc:groupstrings");
  }
  n = strlen(arg[1]) + 1;
  groupstrings[ngroups] = new char[n];
  strcpy(groupstrings[ngroups],arg[1]);
  ngroups++;
  
  // set up reneighboring

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep;
  
  // set up seed
  int seed = 259348;
  random_equal = new RanPark(lmp,seed);
}

/* ---------------------------------------------------------------------- */

FixConcBerendsen::~FixConcBerendsen()
{
  delete [] cstr;

  // delete concentration if fix created it

  if (cflag) modify->delete_compute(id_temp);
  delete [] id_temp;

  //For fix_adapt_fep
  for (int m = 0; m < nadapt; m++) {
    if (adapt[m].which == PAIR) {
      delete [] adapt[m].pstyle;
      delete [] adapt[m].pparam;
      memory->destroy(adapt[m].array_orig);
    }
  }
  delete [] adapt;

  // check nfix in case all fixes have already been deleted

  if (id_fix_diam && modify->nfix) modify->delete_fix(id_fix_diam);
  if (id_fix_chg && modify->nfix) modify->delete_fix(id_fix_chg);
  delete [] id_fix_diam;
  delete [] id_fix_chg;
  if (ngroups > 0) {
    for (int igroup = 0; igroup < ngroups; igroup++)
      delete [] groupstrings[igroup];
    memory->sfree(groupstrings);
  }

  if (ngrouptypes > 0) {
    memory->destroy(grouptypes);
    memory->destroy(grouptypebits);
    for (int igroup = 0; igroup < ngrouptypes; igroup++)
      delete [] grouptypestrings[igroup];
    memory->sfree(grouptypestrings);
  }
  delete random_equal;
  delete [] idshake;
}

/* ---------------------------------------------------------------------- */

int FixConcBerendsen::setmask()
{
  int mask = 0;
  //mask |= END_OF_STEP;
  mask |= PRE_EXCHANGE;
  mask |= THERMO_ENERGY;
  mask |= POST_RUN;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixConcBerendsen::init()
{
  // check variable

  if (cstr) {
    cvar = input->variable->find(cstr);
    if (cvar < 0)
      error->all(FLERR,"Variable name for fix temp/berendsen does not exist");
    if (input->variable->equalstyle(cvar)) cstyle = EQUAL;
    else error->all(FLERR,"Variable for fix temp/berendsen is invalid style");
  }

  int icompute = modify->find_compute(id_temp);
  if (icompute < 0)
    error->all(FLERR,"Temperature ID for fix temp/berendsen does not exist");
  concentration = modify->compute[icompute];

  // if shakeflag defined, check for SHAKE fix
  // its molecule template must be same as this one

  // From here, mino
  // For the fix shake mol tag should be used.
  // We can define molecule with several txt files
  // Then the test should be performed.
  fixshake = NULL;
  if (shakeflag && mode == MOLECULE) {
    int ifix = modify->find_fix(idshake);
    if (ifix < 0) error->all(FLERR,"Fix conc/berendsen shake fix does not exist");
    fixshake = modify->fix[ifix];
    int tmp;
    Molecule **shake_molecules = (Molecule **) fixshake->extract("onemol",tmp);
    if (shake_molecules == NULL) {
      error->all(FLERR,"Fix shake must use the mol keyword for fix conc/berendsen");
    }
    if (shake_molecules[0]->nset != 4) {
      error->all(FLERR,"Cation, anion, fractional cation, and anion should be defined as one molecule template ID in fix shake for fix conc/berendsen");
    }
    if (shake_molecules != &onemols[molCat]) {
      error->all(FLERR,"Fix shake and fix conc/berendsen must have same molecule template set");
    }
  }

  if (modify->check_rigid_group_overlap(groupbit))
    error->warning(FLERR,"Cannot thermostat atoms in rigid bodies");

  // init for fix_adapt_fep
  int i,j;

  // allow a dynamic group only if ATOM attribute not used

//  if (group->dynamic[igroup])
//    for (int i = 0; i < nadapt; i++)
//      if (adapt[i].which == ATOM)
//        error->all(FLERR,"Cannot use dynamic group with fix adapt/fep atom");

  // setup and error checks
  anypair = 0;

  for (int m = 0; m < nadapt; m++) {
    Adapt *ad = &adapt[m];

    if (ad->which == PAIR) {
      anypair = 1;
      Pair *pair = NULL;

      if (lmp->suffix_enable) {
        char psuffix[128];
        strcpy(psuffix,ad->pstyle);
        strcat(psuffix,"/");
        strcat(psuffix,lmp->suffix);
        pair = force->pair_match(psuffix,1);
      }
      if (pair == NULL) pair = force->pair_match(ad->pstyle,1);
      if (pair == NULL)
        error->all(FLERR, "Fix adapt/fep pair style does not exist");
      ad->pdim = 2;
      void *ptr = pair->extract(ad->pparam,ad->pdim);
      if (ptr == NULL)
        error->all(FLERR,"Fix adapt/fep pair style param not supported");

      ad->pdim = 2;
      if (ad->pdim == 0) ad->scalar = (double *) ptr;
      if (ad->pdim == 2) ad->array = (double **) ptr;

      // if pair hybrid, test that ilo,ihi,jlo,jhi are valid for sub-style

      if (ad->pdim == 2 && (strcmp(force->pair_style,"hybrid") == 0 ||
                            strcmp(force->pair_style,"hybrid/overlay") == 0)) {
        PairHybrid *pair = (PairHybrid *) force->pair;
        for (i = ad->ilo; i <= ad->ihi; i++) {
          for (j = MAX(ad->jlo,i); j <= ad->jhi; j++) {
            if (!pair->check_ijtype(i,j,ad->pstyle))
              error->all(FLERR,"Fix adapt/fep type pair range is not valid for "
                         "pair hybrid sub-style");
	  }
	}
      }

    } else if (ad->which == KSPACE) {
      if (force->kspace == NULL)
        error->all(FLERR,"Fix adapt/fep kspace style does not exist");
      kspace_scale = (double *) force->kspace->extract("scale");

    } else if (ad->which == ATOM) {
      if (ad->aparam == CHARGE) {
        if (!atom->q_flag)
          error->all(FLERR,"Fix adapt/fep requires atom attribute charge");
      }
    } else if (ad->which == MOLECULE) {
      if (ad->aparam == CHARGE) {
        if (!atom->q_flag)
          error->all(FLERR,"Fix adapt/fep requires atom attribute charge");
      }
    }

  }

  // make copy of original pair array values

  for (int m = 0; m < nadapt; m++) {
    Adapt *ad = &adapt[m];
    if (ad->which == PAIR && ad->pdim == 2) {
      for (i = ad->ilo; i <= ad->ihi; i++)
        for (j = MAX(ad->jlo,i); j <= ad->jhi; j++)
          ad->array_orig[i][j] = ad->array[i][j];
    }
  }

  // fixes that store initial per-atom values

  if (id_fix_chg) {
    int ifix = modify->find_fix(id_fix_chg);
    if (ifix < 0) error->all(FLERR,"Could not find fix adapt storage fix ID");
    fix_chg = (FixStore *) modify->fix[ifix];
  }

  // construct group bitmask for all new atoms
  // aggregated over all group keywords

  groupbitall = 1 | groupbit;
  for (int igroup = 0; igroup < ngroups; igroup++) {
    int jgroup = group->find(groupstrings[igroup]);
    if (jgroup == -1)
      error->all(FLERR,"Could not find specified fix gcmc group ID");
    groupbitall |= group->bitmask[jgroup];
  }

  // construct group type bitmasks
  // not aggregated over all group keywords

  if (ngrouptypes > 0) {
    memory->create(grouptypebits,ngrouptypes,"fix_gcmc:grouptypebits");
    for (int igroup = 0; igroup < ngrouptypes; igroup++) {
      int jgroup = group->find(grouptypestrings[igroup]);
      if (jgroup == -1)
        error->all(FLERR,"Could not find specified fix gcmc group ID");
      grouptypebits[igroup] = group->bitmask[jgroup];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixConcBerendsen::pre_exchange()
{
  double c_current,c_frac,n_frac;//M unit
  double volume = concentration->vol_assign();//region volume
  double **soft_current;
  int itmp=2;
  soft_current = (double **) force->pair->extract("lambda",itmp);
  concentration->compute_vector();
  if (concfactor == 0) {//Constant averaged concentration mode
      c_current = 0.5 * (concentration->vector[typeCat-1]/double(typeCat_num)/double(Cat_ratio) + concentration->vector[typeAni-1]/double(typeAni_num)/double(Ani_ratio));
      n_frac = 0.5 * (soft_current[typeCat_frac][typeCat_frac] + soft_current[typeAni_frac][typeAni_frac]);
      c_frac = n_frac * NUM_TO_MOL / volume;
      c_current += c_frac;
  }
  else if (concfactor == 1) {//Constant cation concentration mode
      c_current=concentration->vector[typeCat-1]/double(typeCat_num)/double(Cat_ratio);
      n_frac = soft_current[typeCat_frac][typeCat_frac];
      c_frac = n_frac * NUM_TO_MOL / volume;
      c_current += c_frac;
  }
  else if (concfactor == 2) {//Constant anion concentration mode
      c_current=concentration->vector[typeAni-1]/double(typeAni_num)/double(Ani_ratio);
      n_frac = soft_current[typeAni_frac][typeAni_frac];
      c_frac = n_frac * NUM_TO_MOL / volume;
      c_current += c_frac;
  }
  else error->all(FLERR, "Invalid conc_factor value");
  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;
  int del_t;

  // For the case when you read restart & reset timestep -> next_reneighbor is not well-defined at the initial.
  // If you set the unrelated timestep at the initial, it will cause error.
  if (next_reneighbor < update->ntimestep) del_t = update->ntimestep - next_reneighbor;
  else del_t = update->ntimestep;

  // set current c_target
  // if variable conc, evaluate variable, wrap with clear/add

  if (cstyle == CONSTANT) {
    c_target = c_start + delta * (c_stop-c_start);
  }
  else {
    modify->clearstep_compute();
    c_target = input->variable->compute_equal(cvar);
    if (c_target < 0.0)
      error->one(FLERR,
                 "Fix temp/berendsen variable returned negative concentration");
    modify->addstep_compute(update->ntimestep + nevery);
  }

  // mol/L to pts again to evaluate the paritcle insertion/deletion
  double n_current = c_current * volume / NUM_TO_MOL;
  double n_target = c_target * volume / NUM_TO_MOL;
  
  // Calculate del_n
  double del_n = (n_target - n_current) * update->dt * del_t / c_period;
  double n_next = n_current + del_n;
  double adapt_next = n_frac + del_n;

  if (del_n >= 1.0 || del_n <= -1.0)
      error->one(FLERR,
                 "Too large change for Fix conc/berendsen(>1 pts/timestep), Please increase the damping parameter");
  double region_xlo = domain->regions[iregion]->extent_xlo;
  double region_xhi = domain->regions[iregion]->extent_xhi;
  double region_ylo = domain->regions[iregion]->extent_ylo;
  double region_yhi = domain->regions[iregion]->extent_yhi;
  double region_zlo = domain->regions[iregion]->extent_zlo;
  double region_zhi = domain->regions[iregion]->extent_zhi;
  
  if (adapt_next >= 0 && adapt_next <= 1) {//just parameter change
      change_param(adapt_next);
  }
  else if (adapt_next > 1) {//particle insertion
      double adapt_insertion;
      adapt_insertion = adapt_next - 1.0;
      change_param(1.0);//change parameter to 1 && change atom type
      // create new atom pair
      double coord[3];
      int i;
      domain->pbc();
      comm->exchange();
      atom->nghost = 0;
      comm->borders();
      if (mode == ATOM) {
	// Create Cations to match the ratio
	for (i = 0; i < Cat_ratio; i++) {
          coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
          coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
          coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
          if (coord[0] >= domain->sublo[0] && coord[0] < domain->subhi[0] &&
              coord[1] >= domain->sublo[1] && coord[1] < domain->subhi[1] &&
              coord[2] >= domain->sublo[2] && coord[2] < domain->subhi[2]) {
                atom->avec->create_atom(typeCat_frac,coord);
	        int m=atom->nlocal-1;
  	        atom->mask[m] = groupbitall;
	        atom->v[m][0] = 0.0;
	        atom->v[m][1] = 0.0;
	        atom->v[m][2] = 0.0;
	        atom->q[m] = Cat_charge;
	        modify->create_attribute(m);
	        printf("Cation created\n");
	  }
          bigint nblocal = atom->nlocal;
          MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
          if (atom->tag_enable) {
              atom->tag_extend();
              if (atom->map_style) {
                atom->map_init();
              }
          }
          atom->nghost = 0;
          comm->borders();
	}
      }
      if (mode == MOLECULE) {
	// Create Cations to match the ratio
	for (int j = 0; j < Cat_ratio; j++) {
          coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
          coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
          coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
  	  tagint maxmol = 0;
	  int i;
          for (i = 0; i < atom->nlocal; i++) maxmol = MAX(maxmol,atom->molecule[i]);
          tagint maxmol_all;
          MPI_Allreduce(&maxmol,&maxmol_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
          maxmol_all++;
          int insertion_molecule = maxmol_all;
          tagint maxtag = 0;
          for (i = 0; i < atom->nlocal; i++) maxtag = MAX(maxtag,atom->tag[i]);
          tagint maxtag_all;
          MPI_Allreduce(&maxtag,&maxtag_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
	  onemols = atom->molecules;
          onemols[molCat_frac]->compute_center();
	  double r[3],rotmat[3][3],quat[4];
          double rsq = 1.1;
          while (rsq > 1.0) {
            r[0] = 2.0*random_equal->uniform() - 1.0;
            r[1] = 2.0*random_equal->uniform() - 1.0;
            r[2] = 2.0*random_equal->uniform() - 1.0;
            rsq = MathExtra::dot3(r, r);
          }
          double theta = random_equal->uniform() * MY_2PI;
          MathExtra::norm3(r);
          MathExtra::axisangle_to_quat(r,theta,quat);
          MathExtra::quat_to_mat(quat,rotmat);

	  int nlocalprev = atom->nlocal;

          if (coord[0] >= domain->sublo[0] && coord[0] < domain->subhi[0] &&
              coord[1] >= domain->sublo[1] && coord[1] < domain->subhi[1] &&
              coord[2] >= domain->sublo[2] && coord[2] < domain->subhi[2]) {
	    for (i = 0; i < onemols[molCat_frac]->natoms; i++) {
              double xtmp[3];
	      MathExtra::matvec(rotmat,onemols[molCat_frac]->dx[i],xtmp);
	      xtmp[0] += coord[0];
	      xtmp[1] += coord[1];
	      xtmp[2] += coord[2];

	      // need to adjust image flags in remap()
              imagezero = ((imageint) IMGMAX << IMG2BITS) |
                 ((imageint) IMGMAX << IMGBITS) | IMGMAX;
	       
	      imageint imagetmp = imagezero;

	      domain->remap(xtmp,imagetmp);
	      if (!domain->inside(xtmp))
	          error->one(FLERR,"Fix gcmc put atom outside box");
            
              atom->avec->create_atom(onemols[molCat_frac]->type[i],xtmp);
              int m = atom->nlocal - 1;
          
              // add to groups
              // optionally add to type-based groups
          
              atom->image[m] = imagetmp;
              atom->molecule[m] = insertion_molecule;
              atom->tag[m] = maxtag_all + i + 1;
              atom->v[m][0] = 0.0;
              atom->v[m][1] = 0.0;
              atom->v[m][2] = 0.0;
          
              atom->add_molecule_atom(onemols[molCat_frac],i,m,maxtag_all);
              modify->create_attribute(m);
	    }
	    printf("Cation created\n");
	  }
	  // shake
	  if (shakeflag) fixshake->set_molecule(nlocalprev,maxtag_all,molCat_frac-molCat,coord,NULL,quat);

          atom->natoms += onemols[molCat_frac]->natoms;
          atom->nbonds += onemols[molCat_frac]->nbonds;
          atom->nangles += onemols[molCat_frac]->nangles;
          atom->ndihedrals += onemols[molCat_frac]->ndihedrals;
          atom->nimpropers += onemols[molCat_frac]->nimpropers;
          if (atom->tag_enable) {
              if (atom->map_style) {
                 atom->map_init();
              } 
          }
          atom->nghost = 0;
          comm->borders();
        }
      }

      // Anion part
      if (mode == ATOM) {
        // Create Anions to match the ratio
	for (int j = 0; j < Ani_ratio; j++) {
          coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
          coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
          coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
 
          if (coord[0] >= domain->sublo[0] && coord[0] < domain->subhi[0] &&
              coord[1] >= domain->sublo[1] && coord[1] < domain->subhi[1] &&
              coord[2] >= domain->sublo[2] && coord[2] < domain->subhi[2]) {
  	      atom->avec->create_atom(typeAni_frac,coord);
	      int m=atom->nlocal-1;
	      atom->mask[m] = groupbitall;
	      atom->v[m][0] = 0.0;
	      atom->v[m][1] = 0.0;
	      atom->v[m][2] = 0.0;
	      atom->q[m] = Ani_charge;
	      modify->create_attribute(m);
	      printf("Anion created\n");
	  }
          bigint nblocal = atom->nlocal;
          MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
          if (atom->tag_enable) {
              atom->tag_extend();
              if (atom->map_style) {
                atom->map_init();
              }
          }
          atom->nghost = 0;
          comm->borders();
	}
      }
      if (mode == MOLECULE) {
        // Create Anions to match the ratio
	for (int j = 0; j < Ani_ratio; j++) {
          coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
          coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
          coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
          tagint maxmol = 0;
	  int i;
          for (i = 0; i < atom->nlocal; i++) maxmol = MAX(maxmol,atom->molecule[i]);
          tagint maxmol_all;
          MPI_Allreduce(&maxmol,&maxmol_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
          maxmol_all++;
          int insertion_molecule = maxmol_all;
          tagint maxtag = 0;
          for (i = 0; i < atom->nlocal; i++) maxtag = MAX(maxtag,atom->tag[i]);
          tagint maxtag_all;
          MPI_Allreduce(&maxtag,&maxtag_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
	  onemols = atom->molecules;
          onemols[molAni_frac]->compute_center();
	  double r[3],rotmat[3][3],quat[4];
          double rsq = 1.1;
          while (rsq > 1.0) {
            r[0] = 2.0*random_equal->uniform() - 1.0;
            r[1] = 2.0*random_equal->uniform() - 1.0;
            r[2] = 2.0*random_equal->uniform() - 1.0;
            rsq = MathExtra::dot3(r, r);
          }
          double theta = random_equal->uniform() * MY_2PI;
          MathExtra::norm3(r);
          MathExtra::axisangle_to_quat(r,theta,quat);
          MathExtra::quat_to_mat(quat,rotmat);
	  
	  int nlocalprev = atom->nlocal;

          if (coord[0] >= domain->sublo[0] && coord[0] < domain->subhi[0] &&
              coord[1] >= domain->sublo[1] && coord[1] < domain->subhi[1] &&
              coord[2] >= domain->sublo[2] && coord[2] < domain->subhi[2]) {
	    for (i = 0; i <  onemols[molAni_frac]->natoms; i++) {
              double xtmp[3];
	      MathExtra::matvec(rotmat,onemols[molAni_frac]->dx[i],xtmp);
	      xtmp[0] += coord[0];
	      xtmp[1] += coord[1];
	      xtmp[2] += coord[2];

	      // need to adjust image flags in remap()
              imagezero = ((imageint) IMGMAX << IMG2BITS) |
               ((imageint) IMGMAX << IMGBITS) | IMGMAX;
	      imageint imagetmp = imagezero;
	      domain->remap(xtmp,imagetmp);
	      if (!domain->inside(xtmp))
	          error->one(FLERR,"Fix gcmc put atom outside box");
            
              atom->avec->create_atom(onemols[molAni_frac]->type[i],xtmp);
              int m = atom->nlocal - 1;
          
              // add to groups
              // optionally add to type-based groups
          
              atom->image[m] = imagetmp;
              atom->molecule[m] = insertion_molecule;
              atom->tag[m] = maxtag_all + i + 1;
              atom->v[m][0] = 0.0;
              atom->v[m][1] = 0.0;
              atom->v[m][2] = 0.0;
          
              atom->add_molecule_atom(onemols[molAni_frac],i,m,maxtag_all);
              modify->create_attribute(m);
	    }
	  printf("Anion created\n");
	  }
	  // shake
	  if (shakeflag) fixshake->set_molecule(nlocalprev,maxtag_all,molAni_frac-molCat,coord,NULL,quat);

          atom->natoms += onemols[molAni_frac]->natoms;
          atom->nbonds += onemols[molAni_frac]->nbonds;
          atom->nangles += onemols[molAni_frac]->nangles;
          atom->ndihedrals += onemols[molAni_frac]->ndihedrals;
          atom->nimpropers += onemols[molAni_frac]->nimpropers;
          if (atom->tag_enable) {
              if (atom->map_style) {
                 atom->map_init();
              }
          }
          atom->nghost = 0;
          comm->borders();
        }
      }
      domain->pbc();
      comm->exchange();
      atom->nghost = 0;
      comm->borders();
      change_param(adapt_insertion);
      if (force->kspace) force->kspace->qsum_qsq();
  }
  else {//particle deletion
    double adapt_deletion = adapt_next + 1.0;
    change_param(0.0);//change_param
    int i;
    if (mode == ATOM) {
      //delete_atom_pair
      int tot_deletion_number,count,count_all;
      tot_deletion_number = Cat_ratio + Ani_ratio;
      count = 0;
      count_all = 0;
      while (count_all < tot_deletion_number) {
        for (i=0;i<atom->nlocal;i++) {
            if (atom->type[i] == typeCat_frac) {
	       atom->avec->copy(atom->nlocal-1,i,1);
               atom->nlocal--;
	       printf("Cation deleted\n");
	       count++;
            }
            if (atom->type[i] == typeAni_frac) {
	       atom->avec->copy(atom->nlocal-1,i,1);
               atom->nlocal--;
	       printf("Anion deleted\n");
	       count++;
            }
        }
	MPI_Allreduce(&count,&count_all,1,MPI_INT,MPI_SUM,world);
        bigint nblocal = atom->nlocal;
        MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world); 
        if (atom->map_style) atom->map_init();
      }
      if (count_all > tot_deletion_number) error->all(FLERR,"Fractional ions are too many");
      //Choose Cat, Ani in the region -> change atom type
      int i_temp,i_Cattemp[Cat_ratio],i_Anitemp[Ani_ratio];
      int i_Cattemp_all[Cat_ratio],i_Anitemp_all[Ani_ratio];
      for (int j = 0; j < Cat_ratio; j++) {
	  i_Cattemp[j] = -1;
          for (i=0;i<atom->nlocal;i++) {
	    if (atom->type[i]==typeCat) {
	        if (domain->regions[iregion]->match(atom->x[i][0],atom->x[i][1],atom->x[i][2])) {
                  if (j == 0) i_Cattemp[j] = atom->tag[i];
                  else {
		      for (int k = 0; k < j; k++) {
		        if (i_Cattemp[k] == atom->tag[i]) continue;
			i_Cattemp[j] = atom->tag[i];
		      }
		  }
	        }
	    }
          }
          i_Cattemp_all[j]=-1;
          MPI_Allreduce(&i_Cattemp[j],&i_Cattemp_all[j],1,MPI_INT,MPI_MAX,world);
          if (i_Cattemp_all[j] == -1) error->all(FLERR,"No cation exists in the bulk region.\n");
      }
      for (int j = 0; j < Ani_ratio; j++) {
	  i_Anitemp[j] = -1;
          for (i=0;i<atom->nlocal;i++) {
	    if (atom->type[i]==typeAni) {
	        if (domain->regions[iregion]->match(atom->x[i][0],atom->x[i][1],atom->x[i][2])) {
	       	  if (j == 0) i_Anitemp[j] = atom->tag[i];
		  else {
		      for (int k = 0; k < j; k++) {
			  if (i_Anitemp[k] == atom->tag[i]) continue;
			  i_Anitemp[j] = atom->tag[i];
		      }
		  }
	        }
	    }
          }
          i_Anitemp_all[j]=-1;
          MPI_Allreduce(&i_Anitemp[j],&i_Anitemp_all[j],1,MPI_INT,MPI_MAX,world);
          if (i_Anitemp_all[j] == -1) error->all(FLERR,"No anion exists in the bulk region.\n");
      }
     
      // Change Cat, Ani to fractional ions
      for (int j = 0; j < Cat_ratio; j++) {
        for (i=0;i<atom->nlocal;i++) {
	    if (atom->tag[i]==i_Cattemp_all[j] && atom->type[i]==typeCat) {
	        atom->type[i] = typeCat_frac;
	        atom->tag[i] = atom->natoms-Ani_ratio-j;
	    }
	}
      }
      for (int j = 0; j < Ani_ratio; j++) {
	  for (i=0;i<atom->nlocal;i++) {
	    if (atom->tag[i]==i_Anitemp_all[j] && atom->type[i]==typeAni) {
	        atom->type[i]=typeAni_frac;
	        atom->tag[i] = atom->natoms-j;
	    }
	  }
      }
      // Change atom id for last two atoms
      for (int j = 0; j < Cat_ratio; j++) {
        for (i=0;i<atom->nlocal;i++) {
	    if (atom->tag[i]==atom->natoms-Ani_ratio-j && atom->type[i]==typeCat) atom->tag[i]=i_Cattemp_all[j];
        }
      }
      for (int j = 0; j < Ani_ratio; j++) {
        for (i=0;i<atom->nlocal;i++) {
	    if (atom->tag[i]==atom->natoms-j && atom->type[i]==typeAni) atom->tag[i]=i_Anitemp_all[j];
        }
      }
    }
    if (mode == MOLECULE) {
      onemols = atom->molecules;
      // Delete fractional molecule
      int tot_deletion_number,count,count_all;
      tot_deletion_number = Cat_ratio * onemols[molCat]->natoms + Ani_ratio * onemols[molAni]->natoms;
      count = 0;
      count_all = 0;
      while (count_all < tot_deletion_number) {
        for (i=0;i<atom->nlocal;i++) {
            if ((atom->type[i] >= typeCat_frac) && (atom->type[i] <= typeCat_frac_max)) {
	       atom->avec->copy(atom->nlocal-1,i,1);
               atom->nlocal--;
	       count++;
            }
            if ((atom->type[i] >= typeAni_frac) && (atom->type[i] <= typeAni_frac_max)) {
	       atom->avec->copy(atom->nlocal-1,i,1);
               atom->nlocal--;
	       count++;
            }
        }
	MPI_Allreduce(&count,&count_all,1,MPI_INT,MPI_SUM,world);
      }
      if (count_all > tot_deletion_number) {
	      printf("Deleted %d Frational ions %d\n",count_all, tot_deletion_number);
	      error->all(FLERR,"Fractional ions are too many");
      }
      if (comm->me ==0) printf("Cations & anions are deleted\n");
      bigint nblocal = atom->nlocal;
      MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world); 
      atom->nbonds -= onemols[molCat_frac]->nbonds * Cat_ratio;
      atom->nangles -= onemols[molCat_frac]->nangles * Cat_ratio;
      atom->ndihedrals -= onemols[molCat_frac]->ndihedrals * Cat_ratio;
      atom->nimpropers -= onemols[molCat_frac]->nimpropers * Cat_ratio;
      atom->nbonds -= onemols[molAni_frac]->nbonds * Ani_ratio;
      atom->nangles -= onemols[molAni_frac]->nangles * Ani_ratio;
      atom->ndihedrals -= onemols[molAni_frac]->ndihedrals * Ani_ratio;
      atom->nimpropers -= onemols[molAni_frac]->nimpropers * Ani_ratio;

      if (atom->map_style) atom->map_init();

      // Count number of Cations & Anions: For the case of charged system
      int count_cat, count_cat_all, count_ani, count_ani_all;
      count_cat = count_ani = 0;
      for (i=0;i<atom->nlocal;i++) {
          if ((atom->type[i] >= typeAni) && (atom->type[i] <= typeAni_max)){
             count_ani++;
          }
          if ((atom->type[i] >= typeCat) && (atom->type[i] <= typeCat_max)){
             count_cat++;
          }
      }
      MPI_Allreduce(&count_cat,&count_cat_all,1,MPI_INT,MPI_SUM,world);
      MPI_Allreduce(&count_ani,&count_ani_all,1,MPI_INT,MPI_SUM,world);
      int num_cat=count_cat_all/onemols[molCat]->natoms;
      int num_ani=count_ani_all/onemols[molAni]->natoms;
      int charging_flag; // 2 if excess anion exist, 1 if excess cation exist, 0 if PZC
      int num_ion_neu; // total ion number - ion number for charging
      num_ion_neu = count_cat_all + count_ani_all; // PZC case
      if (comm->me == 0) printf("num_ion_neu before %d\n",num_ion_neu);
      if (num_cat%Cat_ratio != 0 || num_ani%Ani_ratio != 0) {
	      if (num_cat%Cat_ratio != 0) {
		      charging_flag = 1;
		      num_ion_neu = count_ani_all + num_ani/Ani_ratio*Cat_ratio*onemols[molCat]->natoms;
	      }
	      if (num_ani%Ani_ratio != 0) {
		      charging_flag = 2;
		      num_ion_neu = count_cat_all + num_cat/Cat_ratio*Ani_ratio*onemols[molAni]->natoms;
	      }
      }
      else {
	      if (num_cat/Cat_ratio > num_ani/Ani_ratio) {
		      charging_flag = 1;
		      num_ion_neu = count_ani_all + num_ani/Ani_ratio*Cat_ratio*onemols[molCat]->natoms;
	      }
	      else if (num_cat/Cat_ratio < num_ani/Ani_ratio) {
		      charging_flag = 2;
		      num_ion_neu = count_cat_all + num_cat/Cat_ratio*Ani_ratio*onemols[molAni]->natoms;
	      }
	      else charging_flag = 0;
      }
      if (comm->me == 0) printf("num_ion_neu %d\n",num_ion_neu);
      // Choose Cat, Ani in the region(minimum atom tag) -> charge atom type
      int i_temp,i_Cattemp[Cat_ratio],i_Anitemp[Ani_ratio],i_Cattemp_all[Cat_ratio],i_Anitemp_all[Ani_ratio];
      int max_cat_i, max_ani_i, max_cat_i_all, max_ani_i_all;
      int max_cat_mol, max_ani_mol, max_cat_mol_all, max_ani_mol_all;
      int cat_order, ani_order;
      max_cat_i = max_ani_i = 0;
      for (int j = 0; j < Cat_ratio; j++) {
        i_Cattemp[j] = -1;
        for (i=0;i<atom->nlocal;i++) {
	    if (atom->type[i]>=typeCat && atom->type[i]<=typeCat_max) {
	        if (domain->regions[iregion]->match(atom->x[i][0],atom->x[i][1],atom->x[i][2])) {
	  	    if (j == 0) {
		      i_Cattemp[j] = atom->tag[i];
		      // cat_order: In cation, natoms-1, natoms-2, ..., 0 order
		      // change every cation id into first atom of each molecule
		      if (charging_flag == 2 || charging_flag == 0) {
		        cat_order = (((atom->natoms-i_Cattemp[j])%(onemols[molCat]->natoms*Cat_ratio+onemols[molAni]->natoms*Ani_ratio))-(onemols[molAni]->natoms*Ani_ratio))%(onemols[molCat]->natoms);
		        if (cat_order != (onemols[molCat]->natoms - 1)) i_Cattemp[j] -= (onemols[molCat]->natoms - 1) - cat_order;
		      }
		      else {
			if ((atom->natoms-i_Cattemp[j]) >= num_ion_neu) {// cation for charging is selected
			  cat_order = (atom->natoms-i_Cattemp[j]-num_ion_neu)%(onemols[molCat]->natoms);
			  if (cat_order != (onemols[molCat]->natoms - 1)) i_Cattemp[j] -= (onemols[molCat]->natoms - 1) - cat_order;
			}
			else {
			  cat_order = (((atom->natoms-i_Cattemp[j])%(onemols[molCat]->natoms*Cat_ratio+onemols[molAni]->natoms*Ani_ratio))-(onemols[molAni]->natoms*Ani_ratio))%(onemols[molCat]->natoms);
			  if (cat_order != (onemols[molCat]->natoms - 1)) i_Cattemp[j] -= (onemols[molCat]->natoms - 1) - cat_order;
			}
		      }
  		    }
		    else {
		      i_temp = atom->tag[i];
		      // cat_order: In cation, natoms-1, natoms-2, ..., 0 order
		      // change every cation id into first atom of each molecule
		      if (charging_flag == 2 || charging_flag == 0) {
		        cat_order = (((atom->natoms-i_temp)%(onemols[molCat]->natoms*Cat_ratio+onemols[molAni]->natoms*Ani_ratio))-(onemols[molAni]->natoms*Ani_ratio))%(onemols[molCat]->natoms);
		        if (cat_order != (onemols[molCat]->natoms - 1)) i_temp -= (onemols[molCat]->natoms - 1) - cat_order;
		      }
		      else {
			if ((atom->natoms-i_temp) >= num_ion_neu) {// cation for charging is selected
			  cat_order = (atom->natoms-i_temp-num_ion_neu)%(onemols[molCat]->natoms);
			  if (cat_order != (onemols[molCat]->natoms - 1)) i_temp -= (onemols[molCat]->natoms - 1) - cat_order;
			}
			else {
			  cat_order = (((atom->natoms-i_temp)%(onemols[molCat]->natoms*Cat_ratio+onemols[molAni]->natoms*Ani_ratio))-(onemols[molAni]->natoms*Ani_ratio))%(onemols[molCat]->natoms);
			  if (cat_order != (onemols[molCat]->natoms - 1)) i_temp -= (onemols[molCat]->natoms - 1) - cat_order;
			}
		      }
                      int check_repetition=0;
		      for (int k = 0; k < j; k++) {
 		        if (i_Cattemp_all[k] == i_temp) check_repetition++;
	              }
                      if (check_repetition == 0) i_Cattemp[j] = i_temp;
		    }
	        }
            }
        }
	//printf("%d\t%d\n",j,i_Cattemp[j]);
        i_Cattemp_all[j]=-1;
        MPI_Allreduce(&i_Cattemp[j],&i_Cattemp_all[j],1,MPI_INT,MPI_MAX,world);
        if (i_Cattemp_all[j] == -1) error->all(FLERR,"Cation for fractionalization is not exist.\n");
      }
      // adjust the i_Cattemp_all[j] when the atom id is already in the frac ids -> make it unchanged for atom ids
      for (int j = 0; j < Cat_ratio; j++) {
	      // i_Cattemp_all[j] is in the frac ids
	      if (i_Cattemp_all[j] >= atom->natoms-onemols[molCat]->natoms*Cat_ratio-onemols[molAni]->natoms*Ani_ratio+1 && i_Cattemp_all[j] <= atom->natoms-onemols[molAni]->natoms*Ani_ratio) {
		      // Prevent the atom id changing
		      if ( ((atom->natoms - (i_Cattemp_all[j])) - (onemols[molCat]->natoms*(j+1) + onemols[molAni]->natoms*Ani_ratio - 1)) != 0) {
			      for (int k = 0; k < Cat_ratio; k++) {
				      if (j==k) continue;
				      if (i_Cattemp_all[j] == (i_Cattemp_all[k]+(atom->natoms - i_Cattemp_all[k]) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1))) {
					      int switch_ids;
					      switch_ids = i_Cattemp_all[j];
					      i_Cattemp_all[j] = i_Cattemp_all[k];
					      i_Cattemp_all[k] = switch_ids;
				      }
			      }
		      }
	      }
      }
      // adjust the i_Anitemp_all[j] when the atom id is already in the frac ids -> make it unchanged for atom ids
      for (int j = 0; j < Ani_ratio; j++) {
	      // i_Anitemp_all[j] is in the frac ids
	      if (i_Anitemp_all[j] >= atom->natoms-onemols[molAni]->natoms*Ani_ratio+1 && i_Anitemp_all[j] <= atom->natoms) {
		      // Prevent the atom id changing
		      if ((atom->natoms - (i_Anitemp_all[j]+onemols[molAni]->natoms*(j+1)-1)) != 0) {
			      for (int k = 0; k < Cat_ratio; k++) {
				      if (j==k) continue;
				      if (i_Anitemp_all[j] == (i_Anitemp_all[k]+(atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1)))) {
					      int switch_ids;
					      switch_ids = i_Anitemp_all[j];
					      i_Anitemp_all[j] = i_Anitemp_all[k];
					      i_Anitemp_all[k] = switch_ids;
				      }
			      }
		      }
	      }
      }

      for (int j = 0; j < Ani_ratio; j++) {
        i_Anitemp[j] = -1;
        for (i=0;i<atom->nlocal;i++) {
	    if (atom->type[i]>=typeAni && atom->type[i]<=typeAni_max) {
	        if (domain->regions[iregion]->match(atom->x[i][0],atom->x[i][1],atom->x[i][2])) {
	  	    if (j == 0) {
		      i_Anitemp[j] = atom->tag[i];
		      // ani_order: In anion, natoms-1, natoms-2, ..., 0 order
		      // change every cation id into first atom of each molecule
		      if (charging_flag == 1 || charging_flag == 0) {
			ani_order = (atom->natoms-i_Anitemp[j])%(onemols[molCat]->natoms*Cat_ratio+onemols[molAni]->natoms*Ani_ratio)%(onemols[molAni]->natoms);
			if (ani_order != (onemols[molAni]->natoms - 1)) i_Anitemp[j] -= (onemols[molAni]->natoms - 1) - ani_order;
		      }
		      else {
			if ((atom->natoms-i_Anitemp[j]) >= num_ion_neu) {// anion for charging is selected
			  ani_order = (atom->natoms-i_Anitemp[j]-num_ion_neu)%(onemols[molAni]->natoms);
			  if (ani_order != (onemols[molAni]->natoms - 1)) i_Anitemp[j] -= (onemols[molAni]->natoms - 1) - ani_order;
			}
			else {
			  ani_order = (atom->natoms-i_Anitemp[j])%(onemols[molCat]->natoms*Cat_ratio+onemols[molAni]->natoms*Ani_ratio)%(onemols[molAni]->natoms);
			  if (ani_order != (onemols[molAni]->natoms - 1)) i_Anitemp[j] -= (onemols[molAni]->natoms - 1) - ani_order;
			}
		      }
                    }
		    else {
		      i_temp = atom->tag[i];
		      // ani_order: In anion, natoms-1, natoms-2, ..., 0 order
		      // change every cation id into first atom of each molecule
		      if (charging_flag == 1 || charging_flag == 0) {
			ani_order = (atom->natoms-i_temp)%(onemols[molCat]->natoms*Cat_ratio+onemols[molAni]->natoms*Ani_ratio)%(onemols[molAni]->natoms);
			if (ani_order != (onemols[molAni]->natoms - 1)) i_temp -= (onemols[molAni]->natoms - 1) - ani_order;
		      }
		      else {
			if ((atom->natoms-i_temp) >= num_ion_neu) {// anion for charging is selected
			  ani_order = (atom->natoms-i_temp-num_ion_neu)%(onemols[molAni]->natoms);
			  if (ani_order != (onemols[molAni]->natoms - 1)) i_temp -= (onemols[molAni]->natoms - 1) - ani_order;
			}
			else {
			  ani_order = (atom->natoms-i_temp)%(onemols[molCat]->natoms*Cat_ratio+onemols[molAni]->natoms*Ani_ratio)%(onemols[molAni]->natoms);
			  if (ani_order != (onemols[molAni]->natoms - 1)) i_temp -= (onemols[molAni]->natoms - 1) - ani_order;
			}
		      }
                      int check_repetition=0;
		      for (int k = 0; k < j; k++) {
 		        if (i_Anitemp_all[k] == i_temp) check_repetition++;
	              }
                      if (check_repetition == 0) i_Anitemp[j] = i_temp;
		    }
	        }
            }
        }
        i_Anitemp_all[j]=-1;
        MPI_Allreduce(&i_Anitemp[j],&i_Anitemp_all[j],1,MPI_INT,MPI_MAX,world);
        if (i_Anitemp_all[j] == -1) error->all(FLERR,"Anion for fractionalization is not exist.\n");
      }

      int *num_bond = atom->num_bond;
      int *num_angle = atom->num_angle;
      int *num_dihedral = atom->num_dihedral;
      int *num_improper = atom->num_improper;
      tagint **bond_atom = atom->bond_atom;
      tagint **angle_atom1 = atom->angle_atom1;
      tagint **angle_atom2 = atom->angle_atom2;
      tagint **angle_atom3 = atom->angle_atom3;
      tagint **dihedral_atom1 = atom->dihedral_atom1;
      tagint **dihedral_atom2 = atom->dihedral_atom2;
      tagint **dihedral_atom3 = atom->dihedral_atom3;
      tagint **dihedral_atom4 = atom->dihedral_atom4;
      tagint **improper_atom1 = atom->improper_atom1;
      tagint **improper_atom2 = atom->improper_atom2;
      tagint **improper_atom3 = atom->improper_atom3;
      tagint **improper_atom4 = atom->improper_atom4;
      int **nspecial = atom->nspecial;
      tagint **special = atom->special;

      // Change atom type to frac & change atom id to the last
      if (comm->me==0) {
        printf("natoms %d\n",atom->natoms);
        printf("ids %d\t%d\t%d\n",i_Cattemp_all[0],i_Cattemp_all[1],i_Anitemp_all[0]);
      }
      // Change bond, angle, dihedral, improper information
      for (int k = 0; k < Cat_ratio; k++) {
        for (i=0;i<atom->nlocal;i++) {
  	    if (onemols[molCat]->bondflag) {
	        for (int j = 0; j < num_bond[i]; j++) {
		    if (bond_atom[i][j] >= i_Cattemp_all[k] && bond_atom[i][j]<=i_Cattemp_all[k]+onemols[molCat]->natoms-1) {
		        bond_atom[i][j] += (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		    }
	        }
	    }
	    if (onemols[molCat]->angleflag) {
                for (int j = 0; j < num_angle[i]; j++) {
                    if (angle_atom1[i][j] >= i_Cattemp_all[k] && angle_atom1[i][j]<=i_Cattemp_all[k]+onemols[molCat]->natoms-1) {
                        angle_atom1[i][j] += (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		    }
                    if (angle_atom2[i][j] >= i_Cattemp_all[k] && angle_atom2[i][j]<=i_Cattemp_all[k]+onemols[molCat]->natoms-1) {
                        angle_atom2[i][j] += (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		    }
                    if (angle_atom3[i][j] >= i_Cattemp_all[k] && angle_atom3[i][j]<=i_Cattemp_all[k]+onemols[molCat]->natoms-1) {
                        angle_atom3[i][j] += (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		    }
	        }
	    }
	    if (onemols[molCat]->dihedralflag) {
	        for (int j = 0; j < num_dihedral[i]; j++) {
		    if (dihedral_atom1[i][j] >= i_Cattemp_all[k] && dihedral_atom1[i][j]<=i_Cattemp_all[k]+onemols[molCat]->natoms-1) {
		        dihedral_atom1[i][j] += (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		    }
		    if (dihedral_atom2[i][j] >= i_Cattemp_all[k] && dihedral_atom2[i][j]<=i_Cattemp_all[k]+onemols[molCat]->natoms-1) {
		        dihedral_atom2[i][j] += (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		    }
		    if (dihedral_atom3[i][j] >= i_Cattemp_all[k] && dihedral_atom3[i][j]<=i_Cattemp_all[k]+onemols[molCat]->natoms-1) {
		        dihedral_atom3[i][j] += (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		    }
		    if (dihedral_atom4[i][j] >= i_Cattemp_all[k] && dihedral_atom4[i][j]<=i_Cattemp_all[k]+onemols[molCat]->natoms-1) {
		        dihedral_atom4[i][j] += (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		    }
	        }
	    }
	    if (onemols[molCat]->improperflag) {
	        for (int j = 0; j < num_improper[i]; j++) {
		    if (improper_atom1[i][j] >= i_Cattemp_all[k] && improper_atom1[i][j]<=i_Cattemp_all[k]+onemols[molCat]->natoms-1) {
		        improper_atom1[i][j] += (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		    }
		    if (improper_atom2[i][j] >= i_Cattemp_all[k] && improper_atom2[i][j]<=i_Cattemp_all[k]+onemols[molCat]->natoms-1) {
		        improper_atom2[i][j] += (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		    }
		    if (improper_atom3[i][j] >= i_Cattemp_all[k] && improper_atom3[i][j]<=i_Cattemp_all[k]+onemols[molCat]->natoms-1) {
		        improper_atom3[i][j] += (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		    }
		    if (improper_atom4[i][j] >= i_Cattemp_all[k] && improper_atom4[i][j]<=i_Cattemp_all[k]+onemols[molCat]->natoms-1) {
		        improper_atom4[i][j] += (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		    }
	        }
	    }
	    if (onemols[molCat]->specialflag) {
	        for (int j = 0; j < nspecial[i][2]; j++) {
		    if (special[i][j] >= i_Cattemp_all[k] && special[i][j]<=i_Cattemp_all[k]+onemols[molCat]->natoms-1) {
		        special[i][j] += (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		    }
	        }
	    }
        }
      }
      for (int k = 0; k < Ani_ratio; k++) {
        for (i=0;i<atom->nlocal;i++) {
	  if (onemols[molAni]->bondflag) {
	      for (int j = 0; j < num_bond[i]; j++) {
		  if (bond_atom[i][j] >= i_Anitemp_all[k] && bond_atom[i][j]<=i_Anitemp_all[k]+onemols[molAni]->natoms-1) {
		      bond_atom[i][j] += (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
	      }
	  }
	  if (onemols[molAni]->angleflag) {
              for (int j = 0; j < num_angle[i]; j++) {
                  if (angle_atom1[i][j] >= i_Anitemp_all[k] && angle_atom1[i][j]<=i_Anitemp_all[k]+onemols[molAni]->natoms-1) {
                      angle_atom1[i][j] += (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
                  if (angle_atom2[i][j] >= i_Anitemp_all[k] && angle_atom2[i][j]<=i_Anitemp_all[k]+onemols[molAni]->natoms-1) {
                      angle_atom2[i][j] += (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
                  if (angle_atom3[i][j] >= i_Anitemp_all[k] && angle_atom3[i][j]<=i_Anitemp_all[k]+onemols[molAni]->natoms-1) {
                      angle_atom3[i][j] += (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
	      }
	  }
	  if (onemols[molAni]->dihedralflag) {
              for (int j = 0; j < num_dihedral[i]; j++) {
                  if (dihedral_atom1[i][j] >= i_Anitemp_all[k] && dihedral_atom1[i][j]<=i_Anitemp_all[k]+onemols[molAni]->natoms-1) {
                      dihedral_atom1[i][j] += (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
                  if (dihedral_atom2[i][j] >= i_Anitemp_all[k] && dihedral_atom2[i][j]<=i_Anitemp_all[k]+onemols[molAni]->natoms-1) {
                      dihedral_atom2[i][j] += (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
                  if (dihedral_atom3[i][j] >= i_Anitemp_all[k] && dihedral_atom3[i][j]<=i_Anitemp_all[k]+onemols[molAni]->natoms-1) {
                      dihedral_atom3[i][j] += (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
                  if (dihedral_atom4[i][j] >= i_Anitemp_all[k] && dihedral_atom4[i][j]<=i_Anitemp_all[k]+onemols[molAni]->natoms-1) {
                      dihedral_atom4[i][j] += (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
	      }
	  }
	  if (onemols[molAni]->improperflag) {
              for (int j = 0; j < num_improper[i]; j++) {
                  if (improper_atom1[i][j] >= i_Anitemp_all[k] && improper_atom1[i][j]<=i_Anitemp_all[k]+onemols[molAni]->natoms-1) {
                      improper_atom1[i][j] += (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
                  if (improper_atom2[i][j] >= i_Anitemp_all[k] && improper_atom2[i][j]<=i_Anitemp_all[k]+onemols[molAni]->natoms-1) {
                      improper_atom2[i][j] += (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
                  if (improper_atom3[i][j] >= i_Anitemp_all[k] && improper_atom3[i][j]<=i_Anitemp_all[k]+onemols[molAni]->natoms-1) {
                      improper_atom3[i][j] += (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
                  if (improper_atom4[i][j] >= i_Anitemp_all[k] && improper_atom4[i][j]<=i_Anitemp_all[k]+onemols[molAni]->natoms-1) {
                      improper_atom4[i][j] += (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
	      }
	  }
	  if (onemols[molAni]->specialflag) {
	      for (int j = 0; j < nspecial[i][2]; j++) {
		  if (special[i][j] >= i_Anitemp_all[k] && special[i][j]<=i_Anitemp_all[k]+onemols[molAni]->natoms-1) {
		      special[i][j] += (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
	      }
	  }
        }
      }

      // Change atoms to fractional atoms
      for (int j = 0; j < Cat_ratio; j++) {
        for (i=0;i<atom->nlocal;i++) {
	    if (atom->type[i]>=typeCat && atom->type[i]<=typeCat_max && atom->tag[i]>=i_Cattemp_all[j] && atom->tag[i]<=i_Cattemp_all[j]+onemols[molCat]->natoms-1) {
	        printf("cat to catfrac before %d\t%d\n",atom->tag[i],atom->type[i]);
	        atom->type[i] = typeCat_frac_max - (typeCat_max - atom->type[i]);
	        atom->tag[i] += ((atom->natoms - (i_Cattemp_all[j])) - (onemols[molCat]->natoms*(j+1) + onemols[molAni]->natoms*Ani_ratio - 1));
		if (shakeflag) fixshake->update_arrays(i,((atom->natoms - (i_Cattemp_all[j])) - (onemols[molCat]->natoms*(j+1) + onemols[molAni]->natoms*Ani_ratio - 1)));
                printf("cat to catfrac after %d\t%d\n",atom->tag[i],atom->type[i]);
	    }
        }
      }
      for (int j = 0; j < Ani_ratio; j++) {
        for (i=0;i<atom->nlocal;i++) {
	    if (atom->type[i]>=typeAni && atom->type[i]<=typeAni_max && atom->tag[i]>=i_Anitemp_all[j] && atom->tag[i]<=i_Anitemp_all[j]+onemols[molAni]->natoms-1) {
	        printf("ani to anifrac before %d\t%d\n",atom->tag[i],atom->type[i]);
	        atom->type[i] = typeAni_frac_max - (typeAni_max - atom->type[i]);
	        atom->tag[i] += ((atom->natoms - (i_Anitemp_all[j]+onemols[molAni]->natoms*(j+1)-1)));
		if (shakeflag) fixshake->update_arrays(i,((atom->natoms - (i_Anitemp_all[j]+onemols[molAni]->natoms*(j+1)-1))));
                printf("ani to anifrac after %d\t%d\n",atom->tag[i],atom->type[i]);
	    }
        }
      }

      // Change bond, angle, dihedral, improper information
      for (int k = 0; k < Cat_ratio; k++) {
      for (i=0;i<atom->nlocal;i++) {
	  if (onemols[molCat]->bondflag) {
	      for (int j = 0; j < num_bond[i]; j++) {
		  if (atom->type[i]>=typeCat && atom->type[i]<=typeCat_max && bond_atom[i][j]>=atom->natoms-onemols[molCat]->natoms*Cat_ratio-onemols[molAni]->natoms*Ani_ratio+1 && bond_atom[i][j]<=atom->natoms-onemols[molAni]->natoms*Ani_ratio) {
		      bond_atom[i][j] -= (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		  }
	      }
	  }
	  if (onemols[molCat]->angleflag) {
              for (int j = 0; j < num_angle[i]; j++) {
		  if (atom->type[i]>=typeCat && atom->type[i]<=typeCat_max && angle_atom1[i][j]>=atom->natoms-onemols[molCat]->natoms*Cat_ratio-onemols[molAni]->natoms*Ani_ratio+1 && angle_atom1[i][j]<=atom->natoms-onemols[molAni]->natoms*Ani_ratio) {
                      angle_atom1[i][j] -= (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		  }
		  if (atom->type[i]>=typeCat && atom->type[i]<=typeCat_max && angle_atom2[i][j]>=atom->natoms-onemols[molCat]->natoms*Cat_ratio-onemols[molAni]->natoms*Ani_ratio+1 && angle_atom2[i][j]<=atom->natoms-onemols[molAni]->natoms*Ani_ratio) {
                      angle_atom2[i][j] -= (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		  }
		  if (atom->type[i]>=typeCat && atom->type[i]<=typeCat_max && angle_atom3[i][j]>=atom->natoms-onemols[molCat]->natoms*Cat_ratio-onemols[molAni]->natoms*Ani_ratio+1 && angle_atom3[i][j]<=atom->natoms-onemols[molAni]->natoms*Ani_ratio) {
                      angle_atom3[i][j] -= (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		  }
	      }
	  }
	  if (onemols[molCat]->dihedralflag) {
	      for (int j = 0; j < num_dihedral[i]; j++) {
		  if (atom->type[i]>=typeCat && atom->type[i]<=typeCat_max && dihedral_atom1[i][j]>=atom->natoms-onemols[molCat]->natoms*Cat_ratio-onemols[molAni]->natoms*Ani_ratio+1 && dihedral_atom1[i][j]<=atom->natoms-onemols[molAni]->natoms*Ani_ratio) {
		      dihedral_atom1[i][j] -= (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		  }
		  if (atom->type[i]>=typeCat && atom->type[i]<=typeCat_max && dihedral_atom2[i][j]>=atom->natoms-onemols[molCat]->natoms*Cat_ratio-onemols[molAni]->natoms*Ani_ratio+1 && dihedral_atom2[i][j]<=atom->natoms-onemols[molAni]->natoms*Ani_ratio) {
		      dihedral_atom2[i][j] -= (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		  }
		  if (atom->type[i]>=typeCat && atom->type[i]<=typeCat_max && dihedral_atom3[i][j]>=atom->natoms-onemols[molCat]->natoms*Cat_ratio-onemols[molAni]->natoms*Ani_ratio+1 && dihedral_atom3[i][j]<=atom->natoms-onemols[molAni]->natoms*Ani_ratio) {
		      dihedral_atom3[i][j] -= (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		  }
		  if (atom->type[i]>=typeCat && atom->type[i]<=typeCat_max && dihedral_atom4[i][j]>=atom->natoms-onemols[molCat]->natoms*Cat_ratio-onemols[molAni]->natoms*Ani_ratio+1 && dihedral_atom4[i][j]<=atom->natoms-onemols[molAni]->natoms*Ani_ratio) {
		      dihedral_atom4[i][j] -= (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		  }
	      }
	  }
	  if (onemols[molCat]->improperflag) {
	      for (int j = 0; j < num_improper[i]; j++) {
		  if (atom->type[i]>=typeCat && atom->type[i]<=typeCat_max && improper_atom1[i][j]>=atom->natoms-onemols[molCat]->natoms*Cat_ratio-onemols[molAni]->natoms*Ani_ratio+1 && improper_atom1[i][j]<=atom->natoms-onemols[molAni]->natoms*Ani_ratio) {
		      improper_atom1[i][j] -= (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		  }
		  if (atom->type[i]>=typeCat && atom->type[i]<=typeCat_max && improper_atom2[i][j]>=atom->natoms-onemols[molCat]->natoms*Cat_ratio-onemols[molAni]->natoms*Ani_ratio+1 && improper_atom2[i][j]<=atom->natoms-onemols[molAni]->natoms*Ani_ratio) {
		      improper_atom2[i][j] -= (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		  }
		  if (atom->type[i]>=typeCat && atom->type[i]<=typeCat_max && improper_atom3[i][j]>=atom->natoms-onemols[molCat]->natoms*Cat_ratio-onemols[molAni]->natoms*Ani_ratio+1 && improper_atom3[i][j]<=atom->natoms-onemols[molAni]->natoms*Ani_ratio) {
		      improper_atom3[i][j] -= (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		  }
		  if (atom->type[i]>=typeCat && atom->type[i]<=typeCat_max && improper_atom4[i][j]>=atom->natoms-onemols[molCat]->natoms*Cat_ratio-onemols[molAni]->natoms*Ani_ratio+1 && improper_atom4[i][j]<=atom->natoms-onemols[molAni]->natoms*Ani_ratio) {
		      improper_atom4[i][j] -= (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		  }
	      }
	  }
	  if (onemols[molCat]->specialflag) {
	      for (int j = 0; j < nspecial[i][2]; j++) {
		  if (atom->type[i]>=typeCat && atom->type[i]<=typeCat_max && special[i][j]>=atom->natoms-onemols[molCat]->natoms*Cat_ratio-onemols[molAni]->natoms*Ani_ratio+1 && special[i][j]<=atom->natoms-onemols[molAni]->natoms*Ani_ratio) {
		      special[i][j] -= (atom->natoms - (i_Cattemp_all[k])) - (onemols[molCat]->natoms*(k+1) + onemols[molAni]->natoms*Ani_ratio - 1);
		  }
	      }
	  }
      }
      }

      for (int k = 0; k < Ani_ratio; k++) {
      for (i=0;i<atom->nlocal;i++) {
	  if (onemols[molAni]->bondflag) {
	      for (int j = 0; j < num_bond[i]; j++) {
		  if (atom->type[i]>=typeAni && atom->type[i]<=typeAni_max && bond_atom[i][j]>=atom->natoms-onemols[molAni]->natoms*Ani_ratio+1 && bond_atom[i][j]<=atom->natoms) {
		      bond_atom[i][j] -= (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
	      }
	  }
	  if (onemols[molAni]->angleflag) {
              for (int j = 0; j < num_angle[i]; j++) {
		  if (atom->type[i]>=typeAni && atom->type[i]<=typeAni_max && angle_atom1[i][j]>=atom->natoms-onemols[molAni]->natoms*Ani_ratio+1 && angle_atom1[i][j]<=atom->natoms) {
                      angle_atom1[i][j] -= (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
		  if (atom->type[i]>=typeAni && atom->type[i]<=typeAni_max && angle_atom2[i][j]>=atom->natoms-onemols[molAni]->natoms*Ani_ratio+1 && angle_atom2[i][j]<=atom->natoms) {
                      angle_atom2[i][j] -= (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
		  if (atom->type[i]>=typeAni && atom->type[i]<=typeAni_max && angle_atom3[i][j]>=atom->natoms-onemols[molAni]->natoms*Ani_ratio+1 && angle_atom3[i][j]<=atom->natoms) {
                      angle_atom3[i][j] -= (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
	      }
	  }
	  if (onemols[molAni]->dihedralflag) {
              for (int j = 0; j < num_dihedral[i]; j++) {
		  if (atom->type[i]>=typeAni && atom->type[i]<=typeAni_max && dihedral_atom1[i][j]>=atom->natoms-onemols[molAni]->natoms*Ani_ratio+1 && dihedral_atom1[i][j]<=atom->natoms) {
                      dihedral_atom1[i][j] -= (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
		  if (atom->type[i]>=typeAni && atom->type[i]<=typeAni_max && dihedral_atom2[i][j]>=atom->natoms-onemols[molAni]->natoms*Ani_ratio+1 && dihedral_atom2[i][j]<=atom->natoms) {
                      dihedral_atom2[i][j] -= (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
		  if (atom->type[i]>=typeAni && atom->type[i]<=typeAni_max && dihedral_atom3[i][j]>=atom->natoms-onemols[molAni]->natoms*Ani_ratio+1 && dihedral_atom3[i][j]<=atom->natoms) {
                      dihedral_atom3[i][j] -= (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
		  if (atom->type[i]>=typeAni && atom->type[i]<=typeAni_max && dihedral_atom4[i][j]>=atom->natoms-onemols[molAni]->natoms*Ani_ratio+1 && dihedral_atom4[i][j]<=atom->natoms) {
                      dihedral_atom4[i][j] -= (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
	      }
	  }
	  if (onemols[molAni]->improperflag) {
              for (int j = 0; j < num_improper[i]; j++) {
		  if (atom->type[i]>=typeAni && atom->type[i]<=typeAni_max && improper_atom1[i][j]>=atom->natoms-onemols[molAni]->natoms*Ani_ratio+1 && improper_atom1[i][j]<=atom->natoms) {
                      improper_atom1[i][j] -= (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
		  if (atom->type[i]>=typeAni && atom->type[i]<=typeAni_max && improper_atom2[i][j]>=atom->natoms-onemols[molAni]->natoms*Ani_ratio+1 && improper_atom2[i][j]<=atom->natoms) {
                      improper_atom2[i][j] -= (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
		  if (atom->type[i]>=typeAni && atom->type[i]<=typeAni_max && improper_atom3[i][j]>=atom->natoms-onemols[molAni]->natoms*Ani_ratio+1 && improper_atom3[i][j]<=atom->natoms) {
                      improper_atom3[i][j] -= (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
		  if (atom->type[i]>=typeAni && atom->type[i]<=typeAni_max && improper_atom4[i][j]>=atom->natoms-onemols[molAni]->natoms*Ani_ratio+1 && improper_atom4[i][j]<=atom->natoms) {
                      improper_atom4[i][j] -= (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
	      }
	  }
	  if (onemols[molAni]->specialflag) {
	      for (int j = 0; j < nspecial[i][2]; j++) {
		  if (atom->type[i]>=typeAni && atom->type[i]<=typeAni_max && special[i][j]>=atom->natoms-onemols[molAni]->natoms*Ani_ratio+1 && special[i][j]<=atom->natoms) {
		      special[i][j] -= (atom->natoms - (i_Anitemp_all[k]+onemols[molAni]->natoms*(k+1)-1));
		  }
	      }
	  }
      }
      }

      for (int j = 0; j < Cat_ratio; j++) {
        for (i=0;i<atom->nlocal;i++) {
            if (atom->type[i]>=typeCat && atom->type[i]<=typeCat_max && atom->tag[i]>=atom->natoms-onemols[molCat]->natoms*Cat_ratio-onemols[molAni]->natoms*Ani_ratio+1+onemols[molCat]->natoms*(Cat_ratio-j-1) && atom->tag[i]<=atom->natoms-onemols[molAni]->natoms*Ani_ratio-onemols[molCat]->natoms*j) {
                printf("catfrac to cat before %d\t%d\n",atom->tag[i],atom->type[i]);
	        atom->tag[i] -= ((atom->natoms - (i_Cattemp_all[j])) - (onemols[molCat]->natoms*(j+1) + onemols[molAni]->natoms*Ani_ratio - 1));
		if (shakeflag) fixshake->update_arrays(i,-((atom->natoms - (i_Cattemp_all[j])) - (onemols[molCat]->natoms*(j+1) + onemols[molAni]->natoms*Ani_ratio - 1)));
                printf("catfrac to cat after %d\t%d\n",atom->tag[i],atom->type[i]);
	    }
        }
      }
      for (int j = 0; j < Ani_ratio; j++) {
        for (i=0;i<atom->nlocal;i++) {
	    if (atom->type[i]>=typeAni && atom->type[i]<=typeAni_max && atom->tag[i]>=atom->natoms-onemols[molAni]->natoms*Ani_ratio+1+onemols[molAni]->natoms*(Ani_ratio-j-1) && atom->tag[i]<=atom->natoms-onemols[molAni]->natoms*j) {
                printf("anifrac to ani before %d\t%d\n",atom->tag[i],atom->type[i]);
	        atom->tag[i] -= ((atom->natoms - (i_Anitemp_all[j]+onemols[molAni]->natoms*(j+1)-1)));
		if (shakeflag) fixshake->update_arrays(i,-((atom->natoms - (i_Anitemp_all[j]+onemols[molAni]->natoms*(j+1)-1))));
                printf("anifrac to ani after %d\t%d\n",atom->tag[i],atom->type[i]);
	    }
        }
      }
    }
    if (atom->map_style) atom->map_init();
    atom->nghost = 0;
    comm->borders();
    for (i=0; i<atom->nlocal; i++) {
	//printf("final %d\t%d\n",atom->tag[i],atom->type[i]);
    }
    change_param(adapt_deletion);
  }
  next_reneighbor = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void FixConcBerendsen::post_run()
{
  double c_frac,n_frac;//M unit
  double **soft_current;
  int itmp=2;
  soft_current = (double **) force->pair->extract("lambda",itmp);
  if (concfactor == 0) {//Constant averaged concentration mode
      n_frac = 0.5 * (soft_current[typeCat_frac][typeCat_frac] + soft_current[typeAni_frac][typeAni_frac]);
  }
  else if (concfactor == 1) {//Constant cation concentration mode
      n_frac = soft_current[typeCat_frac][typeCat_frac];
  }
  else if (concfactor == 2) {//Constant anion concentration mode
      n_frac = soft_current[typeAni_frac][typeAni_frac];
  }
  else error->all(FLERR, "Invalid conc_factor value");
 
  if (comm->me == 0) { 
    printf("#######\n");
    printf("Chemostat: final soft lambda value is %lf\n",n_frac);
  }
}

/* ---------------------------------------------------------------------- */

void FixConcBerendsen::change_param(double value)
{
  // change parameter with adapt_next
  // variable evaluation may invoke computes so wrap with clear/add
  int i,j;
  modify->clearstep_compute();

  for (int m = 0; m < nadapt; m++) {
    Adapt *ad = &adapt[m];

    // set global scalar or type pair array values

    if (ad->which == PAIR) {
      if (ad->pdim == 0) {
        if (scaleflag) *ad->scalar = value * ad->scalar_orig;
        else *ad->scalar = value;
      } else if (ad->pdim == 2) {
        if (scaleflag)
          for (i = ad->ilo; i <= ad->ihi; i++) {
            for (j = MAX(ad->jlo,i); j <= ad->jhi; j++) {
              ad->array[i][j] = value*ad->array_orig[i][j];
	    }
	  }
        else {
	  // set lambda between atom type 1~typeAni_max & typeCat_frac ~ typeAni_frac
          for (i = ad->ilo; i <= ad->ihi; i++) {
            for (j = MAX(ad->jlo,i); j <= ad->jhi; j++) {
              ad->array[i][j] = value;
	      //printf("%d\t%d\t%lf\n",i,j,value);
	    }
	  }
	  // set lambda between frac atoms
          for (i = ad->jlo; i <= ad->jhi; i++) {
            for (j = ad->jlo; j <= ad->jhi; j++) {
              ad->array[i][j] = value;
	      //printf("%d\t%d\t%lf\n",i,j,value);
	    }
	  }
	}
      }
      //ad->array[typeCat_frac][typeAni_frac]=1.0;

    // set kspace scale factor

    }
    if (ad->which == KSPACE) {
    //  *kspace_scale = value; //delete due to unstable issue
    }
    // set per atom values, also make changes for ghost atoms
    if (ad->which == ATOM) {
      if (ad->aparam == CHARGE) {
        int *atype = atom->type;
        double *q = atom->q;
        int *mask = atom->mask;
        int nlocal = atom->nlocal;
        int nall = nlocal + atom->nghost;

        for (i = 0; i < nall; i++) {
          if (atype[i] >= ad->ilo && atype[i] <= ad->ihi) {
            if (mask[i] & groupbit) {
		if (atype[i]==typeCat_frac) {
		    q[i] = value; 
		    if (value == 1.0) atom->type[i]=typeCat;
		}
		if (atype[i]==typeAni_frac) {
		    q[i] = -value;
		    if (value == 1.0) atom->type[i]=typeAni;
		}
	    }
	  }
	}
      }
    }
    if (ad->which == MOLECULE) {
      if (ad->aparam == CHARGE) {
        int *atype = atom->type;
	int *atag = atom->tag;
        double *q = atom->q;
        int *mask = atom->mask;
        int nlocal = atom->nlocal;
        int nall = nlocal + atom->nghost;
        
	onemols = atom->molecules;
        int cat_frac,ani_frac,cat_frac_max,ani_frac_max;
        int count_cat, count_ani;
	cat_frac = ani_frac = 0;
	count_cat = count_ani = 0;
        for (i = 0; i < nall; i++) {
          if (atype[i] >= ad->ilo && atype[i] <= ad->ihi) {
            if (mask[i] & groupbit) {
              if (atype[i]>=typeCat_frac && atype[i]<=typeCat_frac_max) {
		  if (count_cat == 0) cat_frac = atag[i];
                  else if (cat_frac <= atag[i]) cat_frac = atag[i];
		  count_cat++;
	      }
              if (atype[i]>=typeAni_frac && atype[i]<=typeAni_frac_max) {
		  if (count_ani == 0) ani_frac = atag[i];
                  else if (ani_frac <= atag[i]) ani_frac = atag[i];
		  count_ani++;
	      }
	    }
	  }
	}
        MPI_Allreduce(&cat_frac,&cat_frac_max,1,MPI_INT,MPI_MAX,world);
        MPI_Allreduce(&ani_frac,&ani_frac_max,1,MPI_INT,MPI_MAX,world);
        // Scale molecular charge with lambda
        for (i = 0; i < nall; i++) {
          if (atype[i] >= ad->ilo && atype[i] <= ad->ihi) {
            if (mask[i] & groupbit) {
              if (atag[i] >= cat_frac_max-onemols[molCat_frac]->natoms*Cat_ratio+1 && atag[i] <= cat_frac_max) {
		  atom->q[i] = onemols[molCat_frac]->q[(atag[i]-cat_frac_max+onemols[molCat_frac]->natoms*Cat_ratio-1)%onemols[molCat_frac]->natoms] * value;
	      }
              if (atag[i] >= ani_frac_max-onemols[molAni_frac]->natoms*Ani_ratio+1 && atag[i] <= ani_frac_max) {
		  atom->q[i] = onemols[molAni_frac]->q[(atag[i]-ani_frac_max+onemols[molAni_frac]->natoms*Ani_ratio-1)%onemols[molAni_frac]->natoms] * value;
	      }
	    }
	  }
	}
        for (i = 0; i < nall; i++) {
          if (atype[i] >= ad->ilo && atype[i] <= ad->ihi) {
            if (mask[i] & groupbit) {
              if (value == 1.0) atom->type[i] -= ad->ihi - ad->jhi; // typeAni_frac_max - typeAni_max
	    }
	  }
	}
      }
    }
  }

  modify->addstep_compute(update->ntimestep + nevery);

  // re-initialize pair styles if any PAIR settings were changed
  // this resets other coeffs that may depend on changed values,
  // and also offset and tail corrections

  if (anypair) force->pair->reinit();

  // reset KSpace charges if charges have changed

  if (chgflag && force->kspace) force->kspace->qsum_qsq();
}

/* ---------------------------------------------------------------------- */

int FixConcBerendsen::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"temp") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    if (cflag) {
      modify->delete_compute(id_temp);
      cflag = 0;
    }
    delete [] id_temp;
    int n = strlen(arg[1]) + 1;
    id_temp = new char[n];
    strcpy(id_temp,arg[1]);

    int icompute = modify->find_compute(id_temp);
    if (icompute < 0)
      error->all(FLERR,"Could not find fix_modify concentration ID");
    concentration = modify->compute[icompute];

    if (concentration->tempflag == 0)
      error->all(FLERR,
                 "Fix_modify concentration ID does not compute concentration");
    if (concentration->igroup != igroup && comm->me == 0)
      error->warning(FLERR,"Group for fix_modify temp != fix group");
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

void FixConcBerendsen::reset_target(double c_new)
{
  c_target = c_start = c_stop = c_new;
}

/* ---------------------------------------------------------------------- */

double FixConcBerendsen::compute_scalar()
{
  return energy;
}

/* ----------------------------------------------------------------------
   extract thermostat properties
------------------------------------------------------------------------- */

void *FixConcBerendsen::extract(const char *str, int &dim)
{
  dim=0;
  if (strcmp(str,"c_target") == 0) {
    return &c_target;
  }
  return NULL;
}
