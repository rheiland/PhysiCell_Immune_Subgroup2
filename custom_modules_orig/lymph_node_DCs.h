#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

#include "./submodel_data_structures.h" 

#ifndef __lymph_node_DCs__
#define __lymph_node_DCs__
	
extern Submodel_Information lymph_node_DCs_info; 
	
void lymph_node_DCs_model_setup( void );

void lymph_node_DCs_model( double dt );

#endif 