/*
 *  This file is part of libsharp.
 *
 *  libsharp is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libsharp is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libsharp; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libsharp is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file sharp_internal.h
 *  Internally used functionality for the spherical transform library.
 *
 *  Copyright (C) 2006-2012 Max-Planck-Society
 *  \author Martin Reinecke \author Dag Sverre Seljebotn
 */

#ifndef PLANCK_SHARP_INTERNAL_H
#define PLANCK_SHARP_INTERNAL_H

#ifdef __cplusplus
#error This header file cannot be included from C++, only from C
#endif

#include "sharp.h"

#define SHARP_MAXTRANS 1

typedef struct
  {
  hpsharp_jobtype type;
  int spin;
  int add_output;
  int nmaps, nalm;
  int flags;
  void **map;
  void **alm;
  complex double *phase;
  double *norm_l;
  complex double *almtmp;
  const hpsharp_geom_info *ginfo;
  const hpsharp_alm_info *ainfo;
  double time;
  int ntrans;
  unsigned long long opcnt;
  } hpsharp_job;

int hpsharp_get_nv_max (void);
int hpsharp_nv_oracle (hpsharp_jobtype type, int spin, int ntrans);

#endif
