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

/* ----------------------------------------------------------------------
   Contributing author: Paul Crozier (SNL)
   Modification : Minho M. Kim (KAIST)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_spring_bound.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "domain.h"
#include "force.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define SMALL 1.0e-10

enum{TETHER};

/* ---------------------------------------------------------------------- */

FixSpringBound::FixSpringBound(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  group2(NULL)
{
  if (narg < 9) error->all(FLERR,"Illegal fix spring bound command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 4;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  dynamic_group_allow = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  if (strcmp(arg[3],"tether") == 0) {
    if (narg != 9) error->all(FLERR,"Illegal fix spring bound command");
    styleflag = TETHER;
    k_spring = force->numeric(FLERR,arg[4]);
    xflag = yflag = zflag = 1;
    if (strcmp(arg[5],"NULL") == 0) xflag = 0;
    else xc = force->numeric(FLERR,arg[5]);
    if (strcmp(arg[6],"NULL") == 0) yflag = 0;
    else yc = force->numeric(FLERR,arg[6]);
    if (strcmp(arg[7],"NULL") == 0) zflag = 0;
    else zc = force->numeric(FLERR,arg[7]);
    r0 = force->numeric(FLERR,arg[8]);
    if (r0 < 0) error->all(FLERR,"R0 < 0 for fix spring bound command");
  } else error->all(FLERR,"Illegal fix spring bound command");

  ftotal[0] = ftotal[1] = ftotal[2] = ftotal[3] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixSpringBound::~FixSpringBound()
{
  delete [] group2;
}

/* ---------------------------------------------------------------------- */

int FixSpringBound::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSpringBound::init()
{
  // recheck that group 2 has not been deleted

  if (group2) {
    igroup2 = group->find(group2);
    if (igroup2 == -1)
      error->all(FLERR,"Fix spring couple group ID does not exist");
    group2bit = group->bitmask[igroup2];
  }

  masstotal = group->mass(igroup);

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixSpringBound::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixSpringBound::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSpringBound::post_force(int vflag)
{
  if (styleflag == TETHER) spring_tether();
}

/* ---------------------------------------------------------------------- */

void FixSpringBound::spring_tether()
{
  double xcm[3];

  if (group->dynamic[igroup])
    masstotal = group->mass(igroup);

  group->xcm(igroup,masstotal,xcm);

  // fx,fy,fz = components of k * (r-r0) / masstotal

  double dx,dy,dz,fx,fy,fz,r,dr;

  dx = xcm[0] - xc;
  dy = xcm[1] - yc;
  dz = xcm[2] - zc;
  if (!xflag) dx = 0.0;
  if (!yflag) dy = 0.0;
  if (!zflag) dz = 0.0;
  r = sqrt(dx*dx + dy*dy + dz*dz);
  r = MAX(r,SMALL);
  dr = r - r0;

  // Applying force only when the dr >= 0.0 -> Smoothly bound the particles
  if (dr >= 0.0) {
    fx = k_spring*dx*dr/r;
    fy = k_spring*dy*dr/r;
    fz = k_spring*dz*dr/r;
    ftotal[0] = -fx;
    ftotal[1] = -fy;
    ftotal[2] = -fz;
    ftotal[3] = sqrt(fx*fx + fy*fy + fz*fz);
    if (dr < 0.0) ftotal[3] = -ftotal[3];
    espring = 0.5*k_spring * dr*dr;

    if (masstotal > 0.0) {
      fx /= masstotal;
      fy /= masstotal;
      fz /= masstotal;
    }
  }
  else {
    fx = 0.0;
    fy = 0.0;
    fz = 0.0;
    ftotal[0] = -fx;
    ftotal[1] = -fy;
    ftotal[2] = -fz;
    ftotal[3] = sqrt(fx*fx + fy*fy + fz*fz);
  }
  // apply restoring force to atoms in group

  double **f = atom->f;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double massone;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massone = rmass[i];
        f[i][0] -= fx*massone;
        f[i][1] -= fy*massone;
        f[i][2] -= fz*massone;
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massone = mass[type[i]];
        f[i][0] -= fx*massone;
        f[i][1] -= fy*massone;
        f[i][2] -= fz*massone;
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixSpringBound::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSpringBound::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy of stretched spring
------------------------------------------------------------------------- */

double FixSpringBound::compute_scalar()
{
  return espring;
}

/* ----------------------------------------------------------------------
   return components of total spring force on fix group
------------------------------------------------------------------------- */

double FixSpringBound::compute_vector(int n)
{
  return ftotal[n];
}
