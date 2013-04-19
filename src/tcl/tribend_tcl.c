#include "tribend_tcl.h"
#include "tribend.h"

#ifdef TRIBEND
#include "communication.h"

int tclcommand_inter_parse_tribend(Tcl_Interp *interp, int bond_type, int argc, char **argv) {
  
    int i1, i2, i3, i4, boo;
    double kb, max, theta0, bood;
    
    if(argc != 8 && argc != 5) {
      Tcl_AppendResult(interp, "The bending triangle initially needs 7 parameters: "
		     "<ind1> <ind2> <ind3> <ind4> <bool> <kb> <maxdist> ", (char *) NULL);
    } else if(argc == 8) {
    
      if (! ARG_IS_I(1, i1)) {
	Tcl_AppendResult(interp, "Triangle index must be int: "
		     "<ind1> ", (char *) NULL);
	return TCL_ERROR;
      }

      if (! ARG_IS_I(2, i2)) {
	Tcl_AppendResult(interp, "Triangle index must be int: "
		     "<ind2> ", (char *) NULL);
	return TCL_ERROR;
      }
  
      if (! ARG_IS_I(3, i3)) {
	Tcl_AppendResult(interp, "Triangle index must be int: "
		     "<ind3> ", (char *) NULL);
	return TCL_ERROR;
      }
      
      if (! ARG_IS_I(4, i4)) {
	Tcl_AppendResult(interp, "Triangle index must be int: "
		     "<ind4> ", (char *) NULL);
	return TCL_ERROR;
      }
      
      if (! ARG_IS_I(5, boo)) {
	Tcl_AppendResult(interp, "Triangle bool must be int: "
		     "<bool> ", (char *) NULL);
	return TCL_ERROR;
      } 
    
      if (! ARG_IS_D(6, kb)) {
	Tcl_AppendResult(interp, "Triangle needs bending modulus (double) : "
		     "<kb> ", (char *) NULL);
	return TCL_ERROR;
      }
  
      if (! ARG_IS_D(7, max)) {
	Tcl_AppendResult(interp, "Triangle needs maximal distance (double) : "
		     "<max> ", (char *) NULL);
	return TCL_ERROR;
      }
  
      CHECK_VALUE(tribend_set_params(bond_type, i1, i2, i3, i4, boo, max,kb), "bond type must be nonnegative");
      
    } else if(argc == 5) {
      
      if (! ARG_IS_D(1, bood)) {
	Tcl_AppendResult(interp, "Fatal Error Resurrecting tribend"
		     "<bool> ", (char *) NULL);
	return TCL_ERROR;
      }
      
      if (! ARG_IS_D(2, theta0)) {
	Tcl_AppendResult(interp, "Fatal Error Resurrecting tribend"
		     "<theta0> ", (char *) NULL);
	return TCL_ERROR;
      }
      
      if (! ARG_IS_D(3, kb)) {
	Tcl_AppendResult(interp, "Fatal Error Resurrecting tribend"
		     "<kb> ", (char *) NULL);
	return TCL_ERROR;
      }
      
      if (! ARG_IS_D(4, max)) {
	Tcl_AppendResult(interp, "Fatal Error Resurrecting tribend"
		     "<max> ", (char *) NULL);
	return TCL_ERROR;
      }
      
       CHECK_VALUE(tribend_reset_params(bond_type, bood, theta0, kb, max), "bond type must be nonnegative");
      
    }
}

int tclprint_to_result_tribendIA(Tcl_Interp *interp, Bonded_ia_parameters *params)
{
  char buffer[TCL_DOUBLE_SPACE];

   Tcl_PrintDouble(interp, params->p.tribend.boo, buffer);
  Tcl_AppendResult(interp, "tribend ", buffer," ", (char *) NULL);
   Tcl_PrintDouble(interp, params->p.tribend.theta0, buffer);
  Tcl_AppendResult(interp," ",buffer, " ", (char *) NULL);
   Tcl_PrintDouble(interp, params->p.tribend.kb, buffer);
  Tcl_AppendResult(interp," ", buffer, " ", (char *) NULL);
   Tcl_PrintDouble(interp, params->p.tribend.max, buffer);
  Tcl_AppendResult(interp," ", buffer, " ", (char *) NULL);
  
  return TCL_OK;
}

#endif