#ifndef TRIEL_H
#define TRIEL_H

#include "utils.h"
#include "interaction_data.h"
#include "particle_data.h"
#include "grid.h"

#ifdef TRIELASTIC

int triel_set_params(int bond_type, int ind1, int ind2, int ind3, double max, double ks, double ka);
int triel_reset_params(int bond_type, double lo, double lpo, double cospo, double sinpo, double Area0, double max, double ks, double ka);

//Use knowledge that the x-axis in rotates system is parallel to r(p1->p2) in original system;
//To find the corresponding unit vector to y in the rotated system, construct vector perpendicular to r(p1->p2);
MDINLINE void RotateForces(double f1_rot[2], double f2_rot[2], double f1[3], double f2[3], double r1[3], double r2[3]) {
	double xu[3];
	double y[3], yu[3];
	double sca;
	int i;
	unit_vector(r1,xu);
	sca = scalar(r2,xu);
	for(i=0; i<3; i++) {
		y[i] = r2[i] - sca * xu[i];
	}
	
	unit_vector(y,yu);

	f1[0] = f1_rot[0] * xu[0] + f1_rot[1] * yu[0]; f1[1] = f1_rot[0] * xu[1] + f1_rot[1] * yu[1]; f1[2] = f1_rot[0] * xu[2] + f1_rot[1] * yu[2]; 
        f2[0] = f2_rot[0] * xu[0] + f2_rot[1] * yu[0]; f2[1] = f2_rot[0] * xu[1] + f2_rot[1] * yu[1]; f2[2] = f2_rot[0] * xu[2] + f2_rot[1] * yu[2];
    
}

MDINLINE int calc_triel_force(Particle *p_ind1, Particle *p_ind2, Particle *p_ind3,
			      Bonded_ia_parameters *iaparams, double force1[3], double force2[3]) 
{
    double dxy, dxx, dyy;
    double gxy, gyx, gxx, gyy;
    double e1, e2, i1, i2;
    double a1, a2, b1, b2;
    double l, lp, sinp, cosp;
    double vec1[3], vec2[3], vecpro[3];
    double f1_rot[2], f2_rot[2];
    
    //Calculate the current shape of the triangle (l,lp,cos(phi),sin(phi));
    //l = length between 1 and 3
    get_mi_vector(vec2, p_ind3->r.p, p_ind1->r.p);
    //vecsub(p_ind3->r.p,p_ind1->r.p,vec2);
    l = sqrt (sqrlen(vec2));
    //lp = lenght between 1 and 2
    get_mi_vector(vec1, p_ind2->r.p, p_ind1->r.p);
	//vecsub(p_ind2->r.p,p_ind1->r.p,vec1);
    lp = sqrt (sqrlen(vec1));
	//cosp / sinp angle functions between these vectors; calculated directly via the producs
	cosp = scalar(vec1,vec2)/(lp*l);
	vector_product(vec1, vec2, vecpro);
	sinp = sqrt(sqrlen(vecpro))/(l*lp);
    
    if( (lp-iaparams->p.triel.lpo > iaparams->p.triel.maxdist) ||  (l-iaparams->p.triel.lo > iaparams->p.triel.maxdist)) {
    	return 1;
    }
    
    
    //Calculate forces in common plane (after assumed triangle rotation/translation in xy-plane);
    //Note that certain geometries and parameters (e.g. a3=0) can be used to speed the code up.
    //For now it will be played safe and done in detail. 
	dxx = lp/iaparams->p.triel.lpo;
	dxy = ((l*cosp/iaparams->p.triel.lo) - (lp*iaparams->p.triel.cospo/iaparams->p.triel.lpo)) / iaparams->p.triel.sinpo;
	dyy = (l*sinp)/(iaparams->p.triel.lo * iaparams->p.triel.sinpo);
	gxx = SQR(dxx);
	gxy = dxx*dxy;
	gyx = dxx*dxy;
	gyy = SQR(dxy) + SQR(dyy);
	i1 = gxx + gyy - 2;
	i2 = gxx * gyy - gxy * gyx - 1;
	e1 = iaparams->p.triel.ks*(i1+1)/6.0;
        e2 = (-1)*iaparams->p.triel.ks/6.0 + iaparams->p.triel.ka*i2/6.0;
    
    //For sake of better readability shorten the call for the triangle's constants:
    a1 = iaparams->p.triel.a1; a2 = iaparams->p.triel.a2;
	b1 = iaparams->p.triel.b1; b2 = iaparams->p.triel.b2;
	
    f1_rot[0] = iaparams->p.triel.Area0*((-1)*e1*((2*a1*dxx)+(2*b1*dxy))+ (-1)*e2*((gyy*2*a1*dxx)+(-2.0*gyx*(a1*dxy+b1*dxx))+(gxx*2*b1*dxy)));
    f1_rot[1] = iaparams->p.triel.Area0*((-1)*e1*((2*b1*dyy))+ (-1)*e2*((-2.0*gyx*a1*dyy)+(gxx*2*b1*dyy)));
    f2_rot[0] = iaparams->p.triel.Area0*((-1)*e1*((2*a2*dxx)+(2*b2*dxy))+ (-1)*e2*((gyy*2*a2*dxx)+(-2.0*gyx*(a2*dxy+b2*dxx))+(gxx*2*b2*dxy)));
    f2_rot[1] = iaparams->p.triel.Area0*((-1)*e1*((2*b2*dyy))+ (-1)*e2*((-2.0*gyx*a2*dyy)+(gxx*2*b2*dyy)));
    
    //fprintf(stderr, "Forces 2D: (%lf %lf) & (%lf %lf)\n", f1_rot[0], f1_rot[1], f2_rot[0], f2_rot[1]);
    //Rotate forces back into original position of triangle
    RotateForces(f1_rot,f2_rot,force1,force2, vec1, vec2); 
    //fprintf(stderr, "Forces Rot: (%lf %lf %lf) & (%lf %lf %lf\n", force1[0], force1[1], force1[2], force2[0], force2[1], force2[2]); 
     return 0;
}

#endif
#endif
