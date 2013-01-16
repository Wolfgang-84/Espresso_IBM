#include "lbtracers.h"
#include "lb.h"
#include "integrate.h"

#ifdef LBTRACERS

extern double time_step;
//extern double skin (?)

//Update Position ~ Euler
void update_mol_pos_particle(Particle *p) {
	//Do Euler for particle p; assume velocity has already been calculated 
	// & is stored in particle data
	int j;
	double skin2 = SQR(0.5 * skin);
	for(j=0;j<3;j++) {
	    //velocity scaled by *time_step -> simple addition
	    // Check if a particle might have crossed a box border (Verlet criterium); 
	    //if possible resort_particles = 1
		p->r.p[j] = p->r.p[j] + p->m.v[j];
		if(distance2(p[j].r.p,p[j].l.p_old) > skin2 ) resort_particles = 1;
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
	lb_lbfluid_get_interpolated_velocity_lbtrace(p_temp,v_int);
	
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