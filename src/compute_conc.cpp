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

#include <mpi.h>
#include <string.h>
#include "compute_conc.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "comm.h"
#include "group.h"
#include "error.h"
#include "region.h"

using namespace LAMMPS_NS;

#define NUM_TO_MOL 1661.129568 

/* ---------------------------------------------------------------------- */

ComputeConc::ComputeConc(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal compute conc command");

  scalar_flag = vector_flag = 1;
  //size_vector = 2;
  size_vector = atom->ntypes;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;

  //iregion = domain->find_region(arg[3]);

  vector = new double[size_vector];

  int n = strlen(arg[3]) + 1;
  id_temp = new char[n];
  strcpy(id_temp,arg[3]);
  iregion = domain->find_region(id_temp);
  delete [] id_temp;

  region_xlo = domain->regions[iregion]->extent_xlo;
  region_xhi = domain->regions[iregion]->extent_xhi;
  region_ylo = domain->regions[iregion]->extent_ylo;
  region_yhi = domain->regions[iregion]->extent_yhi;
  region_zlo = domain->regions[iregion]->extent_zlo;
  region_zhi = domain->regions[iregion]->extent_zhi;

}

/* ---------------------------------------------------------------------- */

ComputeConc::~ComputeConc()
{
  if (!copymode)
    delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeConc::setup()
{
  dynamic = 0;
  if (dynamic_user || group->dynamic[igroup]) dynamic = 1;
}

/* ---------------------------------------------------------------------- */

double ComputeConc::vol_assign()
{
  vol = ( region_xhi - region_xlo ) * ( region_yhi - region_ylo) * ( region_zhi - region_zlo );
  return vol; 
}

/* ---------------------------------------------------------------------- */

void ComputeConc::compute_vector()
{
  int i;

  invoked_vector = update->ntimestep;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double t[atom->ntypes];
  double **soft_current;
  int itmp=2;
  soft_current = (double **) force->pair->extract("lambda",itmp);
  for (i = 0; i < atom->ntypes; i++) t[i] = 0;
  for (i = 0; i < nlocal; i++){
    if ( mask[i] & groupbit  && domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]) ) {
       t[type[i]-1] += soft_current[atom->type[i]][atom->type[i]];
    }
  }

  MPI_Allreduce(t,vector,atom->ntypes,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < atom->ntypes; i++) vector[i] *= ( NUM_TO_MOL / vol_assign() ) ;
}
