#ifndef TRIBEND_H
#define TRIBEND_H

#include "utils.h"
#include "interaction_data.h"
#include "particle_data.h"
#include "grid.h"


#ifdef TRIBEND

int tribend_set_params(int bond_type, int ind1, int ind2, int ind3, int ind4, int boo, double max, double kb);
int tribend_reset_params(int bond_type, double bood, double theta0, double kb, double max);

MDINLINE int calc_tribend_force(Particle *p1, Particle *p2, Particle *p3, Particle *p4, Bonded_ia_parameters *iaparams) { 
  double theta, Ai, Aj;
  double dx1[3], dx2[3], dx3[3], n1[3], n2[3];
  double Pre, sc, len;
  double v1l[3], v2l[3], v1[3], v2[3], tmp[3], tmp2[3], term1[3], term2[3];
  double direc[3];
  double desc, DTh;
  int i;
  
  
  //Get vectors making up the two triangles
  get_mi_vector(dx1, p1->r.p, p3->r.p);
  get_mi_vector(dx2, p2->r.p, p3->r.p);
  get_mi_vector(dx3, p4->r.p, p3->r.p);

  
  //Get normals on triangle; pointing outwards by definition of indices sequence
  vector_product(dx1, dx2, n1);
  vector_product(dx1, dx3, n2);
    
  if(iaparams->p.tribend.boo == 0) {
      n2[0]=-1*n2[0]; n2[1]=-1*n2[1]; n2[2]=-1*n2[2];
  }
   
  //Get 2*area of triangles out of the magnitude of the resulting normals and make the latter unity 
  Ai = sqrt(n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2]);
  n1[0] = n1[0]/Ai; n1[1]=n1[1]/Ai; n1[2]=n1[2]/Ai; 
  
  Aj = sqrt(n2[0]*n2[0] + n2[1]*n2[1] + n2[2]*n2[2]);
  n2[0] = n2[0]/Aj; n2[1]=n2[1]/Aj; n2[2]=n2[2]/Aj; 
  
  //Get the prefactor for the force term; due to numerical errors the scalar product of two normalized paralell vectors
  // can become slightly above one, which would wreck the acos function > Check that
  sc = scalar(n1,n2);
  
  if(sc>1.0) {
      sc = 1.0;
  }
  
  //Get theta as angle between normals
  theta = acos(sc);

  // Find out if angle is positive or negative
  vector_product(n1,n2,direc);
  desc = scalar(dx1,direc);
  

  if(desc<0) {
      theta = -1.0*theta;
  }
  
  
  //Calculate delta between current and equilibrium angle
  DTh = theta-iaparams->p.tribend.theta0;  
  
  
  if(theta>0) {
    Pre = 1.0*iaparams->p.tribend.kb * sin(DTh);  
  } else {
    Pre = -1.0*iaparams->p.tribend.kb * sin(DTh); 
  }

  
  //Calculate Force on each point. Since the sum of all forces is not zero, the forces will be added right here
  //instead of at the end of the calc_force function
  
  for(i=0; i<3; i++) {
      v1l[i] = n2[i]-sc*n1[i];
      v2l[i] = n1[i]-sc*n2[i];
  }

  
  len = sqrt(sqrlen(v1l));

  
  if(len>0) {
      for(i=0;i <3; i++)
	v1[i]=v1l[i]/len;
  }
  
  len = sqrt(sqrlen(v2l));

  
  if(len>0) {
      for(i=0;i <3; i++)
	v2[i]=v2l[i]/len;
  }
  

  
  //Force for particle 1:
  get_mi_vector(tmp,p2->r.p,p3->r.p); get_mi_vector(tmp2, p3->r.p, p4->r.p);
  vector_product(tmp,v1, term1); vector_product(tmp2,v2, term2);

  
  for(i=0;i<3;i++) {
      p1->f.f[i] += Pre*(term1[i]/Ai + term2[i]/Aj);
  }
  
  //Force for particle 2:
  get_mi_vector(tmp,p3->r.p,p1->r.p);
  vector_product(tmp,v1, term1);

  
  for(i=0;i<3;i++) {
      p2->f.f[i] += Pre*(term1[i]/Ai);
  }
  
  //Force for Particle 3:
  get_mi_vector(tmp,p1->r.p,p2->r.p); get_mi_vector(tmp2, p4->r.p, p1->r.p);
  vector_product(tmp,v1, term1); vector_product(tmp2,v2, term2);
  
  //printf("\np3f: ");
  
  for(i=0;i<3;i++) {
      p3->f.f[i] += Pre*(term1[i]/Ai + term2[i]/Aj);
  }
  
  //Force for Particle 4:
  get_mi_vector(tmp,p1->r.p,p3->r.p);
  vector_product(tmp,v2, term1);
 
  
  for(i=0;i<3;i++) {
      p4->f.f[i] += Pre*(term1[i]/Aj);
  }
 
  
  return 0;
  
}

#endif
#endif
