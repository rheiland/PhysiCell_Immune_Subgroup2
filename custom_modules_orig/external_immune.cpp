#include "./external_immune.h" 

using namespace PhysiCell; 

std::string external_immune_version = "0.4.0"; 

Submodel_Information external_immune_info; 


void external_immune_model_setup( void )
{
		// set version
	external_immune_info.name = "external immune"; 
	external_immune_info.version = external_immune_version; 
		// set functions 
	external_immune_info.main_function = external_immune_model; 
	external_immune_info.phenotype_function = NULL; 
	external_immune_info.mechanics_function = NULL; 
		// what microenvironment variables do I need? 

		// what custom data do I need? 
	//external_immune_info.parameters.doubles.push_back( "DM" );
	//external_immune_info.parameters.doubles.push_back( "TC" );

	// submodel_registry.register_model( internal_viral_dynamics_info ); 
	external_immune_info.register_model();
	
	return; 
}

void external_immune_model( double dt )
{
	extern double DM;
	DM -= DM*parameters.doubles("DM_clearance_rate")*dt;
	
	return; 
}
	
