#include "./lymph_node_DCs.h" 

using namespace PhysiCell; 

std::string lymph_node_DCs_version = "0.4.0"; 

Submodel_Information lymph_node_DCs_info; 

void lymph_node_DCs_model_setup( void )
{
		// set version
	lymph_node_DCs_info.name = "lymph node DCs"; 
	lymph_node_DCs_info.version = lymph_node_DCs_version; 
		// set functions 
	lymph_node_DCs_info.main_function = lymph_node_DCs_model; 
	lymph_node_DCs_info.phenotype_function = NULL; 
	lymph_node_DCs_info.mechanics_function = NULL; 
		// what microenvironment variables do I need? 

		// what custom data do I need? 
	//external_immune_info.parameters.doubles.push_back( "DM" );
	//external_immune_info.parameters.doubles.push_back( "TC" );

	// submodel_registry.register_model( internal_viral_dynamics_info ); 
	lymph_node_DCs_info.register_model();
	
	return; 
}

void lymph_node_DCs_model( double dt )
{
	// bookkeeping -- find microenvironment variables we need

	extern double DM;
	
	std::cout<<DM<<std::endl;
	
	//DM=;
	
	return; 
}
	
