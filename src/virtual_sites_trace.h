#ifndef _VIRTUAL_SITES_TRACE_H
#define _VIRTUAL_SITES_TRACE_H

#include "config.h"
#include "particle_data.h"

#ifdef _VIRTUAL_SITES_TRACE

//Update Position ~ Euler/Runge-Kutta
void update_mol_pos_particle(Particle *);
//Update Velocity ~ Get interpolated velocity of LB
void update_mol_vel_particle(Particle *);
//Since no 'real' particles are involved, this function will stay empty
void distribute_mol_force();

#endif

#endif