#include "./internal_viral_response.h" 

using namespace PhysiCell; 

std::string internal_virus_response_version = "0.4.0"; 

Submodel_Information internal_virus_response_model_info; 

void simple_internal_virus_response_model_setup( void )
{
	// set up the model 
		// set version info 
	internal_virus_response_model_info.name = "internal viral response"; 
	internal_virus_response_model_info.version = internal_virus_response_version; 
		// set functions 
	internal_virus_response_model_info.main_function = NULL; 
	internal_virus_response_model_info.phenotype_function = simple_internal_virus_response_model; 
	internal_virus_response_model_info.mechanics_function = NULL; 
	
	// what microenvironment variables do you need 
	internal_virus_response_model_info.microenvironment_variables.push_back( "virion" ); 
	internal_virus_response_model_info.microenvironment_variables.push_back( "VTEST" ); 		
	internal_virus_response_model_info.microenvironment_variables.push_back( "interferon 1" ); 	
	internal_virus_response_model_info.microenvironment_variables.push_back( "pro-inflammatory cytokine" ); 	
	internal_virus_response_model_info.microenvironment_variables.push_back( "debris" ); 	
	internal_virus_response_model_info.microenvironment_variables.push_back( "chemokine" ); 			
		// what microenvironment variables do you expect? 
		// what custom data do I need? 
	internal_virus_response_model_info.cell_variables.push_back( "r_max" ); 
	internal_virus_response_model_info.cell_variables.push_back( "max_apoptosis_half_max" ); 
	internal_virus_response_model_info.cell_variables.push_back( "apoptosis_hill_power" ); 	
	internal_virus_response_model_info.cell_variables.push_back( "infected_cell_chemokine_secretion_activated" );
	internal_virus_response_model_info.cell_variables.push_back( "infected_cell_chemokine_secretion_rate" );
	internal_virus_response_model_info.cell_variables.push_back( "activated_cytokine_secretion_rate" );
	
	
	internal_virus_response_model_info.cell_variables.push_back( "VEn" ); 	
	internal_virus_response_model_info.cell_variables.push_back( "Vnuc" ); 
	internal_virus_response_model_info.cell_variables.push_back( "VRel" ); 
	
		// register the submodel  
	internal_virus_response_model_info.register_model();	
		// set functions for the corresponding cell definition 
		
	
	return; 
}

void simple_internal_virus_response_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	
	// if amount of intracellular virions is less than 1 virion, don't do anything
	double Vnuc = pCell->custom_data["Vnuc"];
	if(pCell->custom_data["Vnuc"]<1e-6)
	{return;}
	
	static int lung_epithelial_type = get_cell_definition( "lung epithelium" ).type; 
	// if not lung epithelium, do nothing 
	if( pCell->type != lung_epithelial_type )
	{ return; } 
	
			
	static int vtest_external = microenvironment.find_density_index( "VTEST" ); 	
	static int chemokine_index = microenvironment.find_density_index( "chemokine" );
	static int IFN_index = microenvironment.find_density_index( "interferon 1" );
	static int proinflammatory_cytokine_index = microenvironment.find_density_index( "pro-inflammatory cytokine");
	
	if(pCell->custom_data["antiviral_state"]>0.5)
{		//cell is in antiviral state so should not be producing virus or signalling to immune cells
		pCell->phenotype.secretion.secretion_rates[proinflammatory_cytokine_index] = 0;
		pCell->phenotype.secretion.secretion_rates[chemokine_index] = 0;
		pCell->custom_data["infected_cell_chemokine_secretion_activated"] = 0.5; 
		pCell->phenotype.secretion.secretion_rates[vtest_external] = 0;
		pCell->phenotype.molecular.internalized_total_substrates[vtest_external] = 0;
		pCell->custom_data["Vnuc"] = 0;
		return;
	}
	
	double Vvoxel = microenvironment.mesh.voxels[1].volume;
		
	simple_viral_secretion_model( pCell, phenotype, dt );
			
				
	if( Vnuc > parameters.doubles("v_rep")/Vvoxel - 1e-16 ) 
	{
		pCell->custom_data["infected_cell_chemokine_secretion_activated"] = 1.0; 
	}

	if( pCell->custom_data["infected_cell_chemokine_secretion_activated"] > 0.5 && phenotype.death.dead == false )
	{
		phenotype.secretion.saturation_densities[chemokine_index] = 1.0;
		phenotype.secretion.secretion_rates[IFN_index]=parameters.doubles("IFN_secretion_rate");

		// (Adrianne) adding pro-inflammatory cytokine secretion by infected cells
		if(pCell->custom_data["antiviral_state"]<0.5)
		{
			pCell->phenotype.secretion.secretion_rates[proinflammatory_cytokine_index] = pCell->custom_data["activated_cytokine_secretion_rate"];
			pCell->phenotype.secretion.secretion_rates[chemokine_index] = pCell->custom_data["infected_cell_chemokine_secretion_rate"];//rate;
		}
		
		//phenotype.secretion.secretion_rates[vtest_external]  = 1;
	}
	
	if(pCell->custom_data["antiviral_state"]>0.5)
	{
			pCell->phenotype.secretion.secretion_rates[proinflammatory_cytokine_index] = 0;
			pCell->phenotype.secretion.secretion_rates[chemokine_index] = 0;
			pCell->custom_data["infected_cell_chemokine_secretion_activated"] = 0.5; 
	}
		
	return; 	
}

void simple_viral_secretion_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	
	static int vtest_external = microenvironment.find_density_index( "VTEST" ); 
	static int proinflammatory_cytokine_index = microenvironment.find_density_index( "pro-inflammatory cytokine");
			
	double Vvoxel = microenvironment.mesh.voxels[1].volume;
		
	if(pCell->phenotype.molecular.internalized_total_substrates[vtest_external]*Vvoxel>8e3 && PhysiCell_globals.current_time>pCell->custom_data["eclipse_time"])
	{
		pCell->phenotype.secretion.secretion_rates[vtest_external] = parameters.doubles("kRel");
	}
	else
	{pCell->phenotype.secretion.secretion_rates[vtest_external] = 0;}	

	if(pCell->custom_data["antiviral_state"]>0.5) // cell is in an antiviral state
	{			
		pCell->phenotype.secretion.secretion_rates[vtest_external] = 0;
		pCell->phenotype.secretion.secretion_rates[proinflammatory_cytokine_index] = 0;
		pCell->phenotype.molecular.internalized_total_substrates[vtest_external] = 0;
		pCell->custom_data["Vnuc"] = 0;
	}
		
	return;
}




