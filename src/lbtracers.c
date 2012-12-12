#include "lbtracers.h"
#include "lb.h"

#ifdef LBTRACERS

extern double time_step;

//Update Position ~ Euler
void update_mol_pos_particle(Particle *p) {
	//Do Euler for particle p; assume velocity has already been calculated 
	// & is stored in particle data
	int j;
	for(j=0;j<3;j++) {
	    //velocity scaled by *time_step -> simple addition
		p->r.p[j] = p->r.p[j] + p->m.v[j];
	}
}

//Update Velocity ~ Get interpolated velocity of LB
void update_mol_vel_particle(Particle *p) {
	int j;
	double v_int[3];// v_intg[3];
	double p_temp[3];
	
	for(j=0;j<3;j++) {
	    p_temp[j] = p->r.p[j];
	}
	
	//Get interpolated velocity from LB
	//lb_lbfluid_get_interpolated_velocity_global(p_temp, v_int);
	lb_lbfluid_get_interpolated_velocity(p_temp,v_int);
	//lb_lbfluid_get_interpolated_velocity_global(p_temp, v_intg);
	
	//if(p->p.identity == 0) {
	  //fprintf(stderr, "(%lf %lf %lf) - (%lf %lf %lf)\n", v_int[0], v_int[1], v_int[2], v_intg[0], v_intg[1], v_intg[2]);
	//}
	
	//rescale velocities on LB-level to MD-level (see viscous coupling)
	for(j=0;j<3;j++) {
		v_int[j] = v_int[j] * time_step;
	}
	
	for(j=0;j<3;j++) {
		p->m.v[j] = v_int[j];
	}
}

//Distribute forces
void distribute_mol_force() {
	//All forces these sites are subject to are influencing the LB fluid, not other
	//particles, therefore no forces need to be distributed here. => Do nothing
}


#endif