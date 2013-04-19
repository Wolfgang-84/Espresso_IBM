#ifndef TRIBEND_TCL_H
#define TRIBEND_TCL_H

#include "parser.h"
#include "interaction_data.h"

#ifdef TRIBEND

int tclcommand_inter_parse_tribend(Tcl_Interp *interp, int bond_type, int argc, char **argv);
int tclprint_to_result_tribendIA(Tcl_Interp *interp, Bonded_ia_parameters *params);

#endif
#endif
