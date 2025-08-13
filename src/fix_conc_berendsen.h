/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(conc/berendsen,FixConcBerendsen)

#else

#ifndef LMP_FIX_CONC_BERENDSEN_H
#define LMP_FIX_CONC_BERENDSEN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixConcBerendsen : public Fix {
 public:
  FixConcBerendsen(class LAMMPS *, int, char **);
  ~FixConcBerendsen();
  int setmask();
  void init();
  //void end_of_step();
  void pre_exchange();
  void post_run();
  int modify_param(int, char **);
  void reset_target(double);
  double compute_scalar();
  virtual void *extract(const char *, int &);

 private:
  int which;
  int typeCat,typeCat_frac,typeAni,typeAni_frac,concfactor;
  int typeCat_num, typeCat_frac_num, typeAni_num, typeAni_frac_num;
  int typeCat_max, typeCat_frac_max, typeAni_max, typeAni_frac_max;
  int molCat, molAni, molCat_frac, molAni_frac;
  int Cat_ratio, Ani_ratio;
  double Cat_charge, Ani_charge;
  int mode; //ATOM or MOLECULE
  double c_start,c_stop,c_period,c_target;
  double energy;
  int cstyle,cvar;
  char *cstr;
  char *idregion;
  void change_param(double value);
  
  int groupbitall;          // group bitmask for inserted atoms
  int ngroups;              // number of group-ids for inserted atoms
  char** groupstrings;      // list of group-ids for inserted atoms
  int ngrouptypes;          // number of type-based group-ids for inserted atoms
  char** grouptypestrings;  // list of type-based group-ids for inserted atoms
  int* grouptypebits;       // list of type-based group bitmasks
  int* grouptypes;          // list of type-based group types

  int iregion;
  char *id_temp;
  //class Compute *temperature;
  class Compute *concentration;//2022.06.17 mino, change variable name
  int cflag;

  int nadapt,anypair,chgflag,resetflag,scaleflag,afterflag;
  char *id_fix_diam,*id_fix_chg;
  //copying from fix_adapt
  class FixStore *fix_diam,*fix_chg;
  struct Adapt {
    int which,ivar;
    char *var;
    char *pstyle,*pparam;
    int ilo,ihi,jlo,jhi;
    int pdim;
    double *scalar,scalar_orig;
    double **array,**array_orig;
    int aparam;
  };

  Adapt *adapt;
  double *kspace_scale;
  class RanPark *random_equal;

  // molecule
  class Molecule **onemols;
  imageint imagezero;
};

}

#endif
#endif
