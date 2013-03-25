#include "vvolume_tcl.h"
#include "vvolume.h"
#include "communication.h"
#include "parser.h"

int tclcallback_vescnum(Tcl_Interp *interp, void *_data) {
  int data = *(int *)_data;
  
  if (data < 0 || data > 200) {
    Tcl_AppendResult(interp, "vescnum must be positive and smaller than 200.", (char *) NULL);
    return (TCL_ERROR);
  }
  vescnum = data;
  mpi_bcast_parameter(FIELD_VESCNUM);
  return (TCL_OK);
}

int tclcallback_vvolo(Tcl_Interp *interp, void *_data) {
	double *data = _data;
	int i;
	
	for(i=0; i<200; i++) {
		
		if(data[i]<0.0) {
			Tcl_AppendResult(interp, "illegal value, Volume must be positive", (char *) NULL);
  		    return (TCL_ERROR);
		}
		
		VVolo[i]=data[i];
	}
	
	mpi_bcast_parameter(FIELD_VVOLO);
	
	return (TCL_OK);
}