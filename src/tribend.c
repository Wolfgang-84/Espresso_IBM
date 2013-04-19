#include "tribend.h"

#ifdef TRIBEND
#include "communication.h"


int tribend_set_params(int bond_type, int ind1, int ind2, int ind3, int ind4, int boo, double max, double kb) {
    Particle p1, p2, p3, p4;
    double n1l[3], n2l[3], n1[3], n2[3];
    double dx1[3], dx2[3], dx3[3];
    double theta0, desc;
    double tmp[3];
    double sc;
    
    if(bond_type<0) {
      return ES_ERROR;
    }
    
    make_bond_type_exist(bond_type);
    
    get_particle_data(ind1, &p1);
    get_particle_data(ind2, &p2);
    get_particle_data(ind3, &p3);
    get_particle_data(ind4, &p4);

    
    //Get vectors of triangles
    get_mi_vector(dx1, p1.r.p, p3.r.p);
    get_mi_vector(dx2, p2.r.p, p3.r.p);
    get_mi_vector(dx3, p4.r.p, p3.r.p);

    
    //Get normals on triangle; pointing outwards by definition of indices sequence
    vector_product(dx1, dx2, n1l);
    vector_product(dx1, dx3, n2l);
    
    
    if(boo == 0) {
	n2l[0]=-1*n2l[0]; n2l[1]=-1*n2l[1]; n2l[2]=-1*n2l[2];
    }
    
    
    unit_vector(n1l,n1);
    unit_vector(n2l,n2);
    
    
    //calculate theta by taking the acos of the scalar n1*n2
    sc = scalar(n1,n2);
  
    if(sc>1.0) {
      sc = 1.0;
    }
    
    theta0 = acos(sc);
    vector_product(n1,n2,tmp);
    
    desc = scalar(dx1,tmp);
    
    if(desc<0) {
	theta0 = -theta0;
    }
    
    
    //effective springconstant = bending resistance^(1/3); Krueger2012
    bonded_ia_params[bond_type].p.tribend.kb = sqrt(3)*kb;
    bonded_ia_params[bond_type].p.tribend.theta0 = theta0;
    bonded_ia_params[bond_type].p.tribend.boo = boo;
    bonded_ia_params[bond_type].p.tribend.max = max;
    
    bonded_ia_params[bond_type].type = TRIBEND_IA;
    bonded_ia_params[bond_type].num = 3;
    
    mpi_bcast_ia_params(bond_type, -1);
    
    //free allocated paticles
    free_particle(&p1);
    free_particle(&p2);
    free_particle(&p3);
    free_particle(&p4);
    
    return ES_OK;
    
}

int tribend_reset_params(int bond_type, double bood, double theta0, double kb, double max) {
    
    int boo = bood;
  
    if(bond_type<0) {
      return ES_ERROR;
    }
    
    make_bond_type_exist(bond_type);
    
    bonded_ia_params[bond_type].p.tribend.boo = boo;
    bonded_ia_params[bond_type].p.tribend.theta0 = theta0;
    bonded_ia_params[bond_type].p.tribend.kb = kb;
    bonded_ia_params[bond_type].p.tribend.max = max;
    
    bonded_ia_params[bond_type].type = TRIBEND_IA;
    bonded_ia_params[bond_type].num = 3;
    
    mpi_bcast_ia_params(bond_type, -1);
    
    return ES_OK;
}

#endif
