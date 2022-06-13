#include "./epithelium_submodel.h" 

using namespace PhysiCell; 

std::string epithelium_submodel_version = "0.4.0"; 

Submodel_Information epithelium_submodel_info; 

void epithelium_contact_function( Cell* pC1, Phenotype& p1, Cell* pC2, Phenotype& p2, double dt )
{
	// elastic adhesions 
	standard_elastic_contact_function( pC1,p1, pC2, p2, dt );
	
	return; 
}

void epithelium_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	
	static int debris_index = microenvironment.find_density_index( "debris");
	static int chemokine_index = microenvironment.find_density_index( "chemokine" );
	static int proinflammatory_cytokine_index = microenvironment.find_density_index( "pro-inflammatory cytokine");
	static int vtest_external = microenvironment.find_density_index( "VTEST" ); 
		
		
	static int apoptosis_index = pCell->phenotype.death.find_death_model_index( "apoptosis" ); 
	
	// viral dynamics model 
	internal_viral_dynamics_info.phenotype_function(pCell,phenotype,dt); 
	// internal_virus_model(pCell,phenotype,dt);
	
	// viral response model 
	internal_virus_response_model_info.phenotype_function(pCell,phenotype,dt); 
	// internal_virus_response_model(pCell,phenotype,dt);	
	
	// T-cell based death
	TCell_induced_apoptosis(pCell, phenotype, dt ); 
	
	// (Adrianne) ROS induced cell death model
	ROS_induced_apoptosis(pCell, phenotype, dt);
			
	// if I am dead, remove all adhesions 
	if( phenotype.death.dead == true )
	{
		// detach all attached cells 
		// remove_all_adhesions( pCell ); 
		phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"]; 
		phenotype.secretion.secretion_rates[vtest_external] = 0; 
	}

	static int IFN_index = microenvironment.find_density_index( "interferon 1" );
	double IFN_internal = pCell->nearest_density_vector()[IFN_index];
	double IC_50_IFN = parameters.doubles("IC_50_IFN");
	double Vvoxel = microenvironment.mesh.voxels[1].volume;
	
	double IFN_prob = IFN_internal/(IC_50_IFN+IFN_internal);
		
	double prob_prob = UniformRandom();
		
	if( pCell->custom_data["antiviral_state_time"] <1e5 && pCell->custom_data["antiviral_state_time"] >=PhysiCell_globals.current_time)
	{ // do nothin
	}
	else if (pCell->custom_data["antiviral_state_time"] < PhysiCell_globals.current_time)
	{
		//cell in antirival state
		pCell->custom_data["antiviral_state"] = 1;
		pCell->phenotype.secretion.secretion_rates[proinflammatory_cytokine_index] = 0;
		pCell->phenotype.secretion.secretion_rates[chemokine_index] = 0;
		pCell->phenotype.secretion.secretion_rates[vtest_external] = 0;
		pCell->phenotype.molecular.internalized_total_substrates[vtest_external] = 0;
		pCell->custom_data["Vnuc"] = 0;	
		
		if( pCell->custom_data["antiviral_state_timer"]<PhysiCell_globals.current_time )// antiviral state timer has expired or hasn't started 
		{
			pCell->custom_data["antiviral_state_timer"] = PhysiCell_globals.current_time+parameters.doubles("tau_IFN");
		}// else it hasn't finished it's antiviral state timer
	}	
	else if(prob_prob<IFN_prob && pCell->custom_data["Vnuc"]<parameters.doubles("Infection_detection_threshold")/Vvoxel) // if stimulation is sufficient cell is antiviral
	{
		// cell enters antiviral state
		//pCell->custom_data["antiviral_state"] = 1;
		pCell->custom_data["antiviral_state_time"] = PhysiCell_globals.current_time+parameters.doubles("IFN_delay");	
	    
		//pCell->phenotype.secretion.secretion_rates[proinflammatory_cytokine_index] = 0;
		//pCell->phenotype.secretion.secretion_rates[chemokine_index] = 0;
		//pCell->phenotype.secretion.secretion_rates[vtest_external] = 0;
		//pCell->phenotype.molecular.internalized_total_substrates[vtest_external] = 0;
		//pCell->custom_data["Vnuc"] = 0;	
		
		//if( pCell->custom_data["antiviral_state_timer"]<PhysiCell_globals.current_time )// antiviral state timer has expired or hasn't started 
		//{
		//	pCell->custom_data["antiviral_state_timer"] = PhysiCell_globals.current_time+parameters.doubles("tau_IFN");
		//}// else it hasn't finished it's antiviral state timer
	}
	else //antiviral stimulation wasn't sufficient
	{
		//check if antiviral state time has expired or if it's not in an antiviral state
		if( pCell->custom_data["antiviral_state_timer"]<PhysiCell_globals.current_time )// it's not meant to stay in antiviral state
		{			
			pCell->custom_data["antiviral_state"] = 0;
		}
	}
	/*
	if(pCell->custom_data["antiviral_state_timer"]<PhysiCell_globals.current_time && prob_prob>IFN_prob )
	{
		pCell->custom_data["antiviral_state"] = 0;
	}
	else if( prob_prob<IFN_prob && pCell->custom_data["antiviral_state"]<1)
	{
		// cell enters antiviral state
		pCell->custom_data["antiviral_state"] = 1;
		
		// start the counter for the antiviral state - only lasts for tau_IFN hours
		pCell->custom_data["antiviral_state_timer"] = PhysiCell_globals.current_time+parameters.doubles("tau_IFN");
		pCell->phenotype.secretion.secretion_rates[proinflammatory_cytokine_index] = 0;
		pCell->phenotype.secretion.secretion_rates[chemokine_index] = 0;
		pCell->phenotype.secretion.secretion_rates[vtest_external] = 0;
		pCell->phenotype.molecular.internalized_total_substrates[vtest_external] = 0;
		pCell->custom_data["Vnuc"] = 0;	
	}*/
		
	// if I am dead, don't bother executing this function again 
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL; 
	}
	
	
	if(pCell->custom_data["antiviral_state"]>0.5)
	{
		pCell->phenotype.secretion.secretion_rates[vtest_external] = 0;
		pCell->phenotype.molecular.internalized_total_substrates[vtest_external] = 0;
		pCell->custom_data["Vnuc"] = 0;	
	}
	
	return; 
}

void epithelium_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int debris_index = microenvironment.find_density_index( "debris");
	static int proinflammatory_cytokine_index = microenvironment.find_density_index( "pro-inflammatory cytokine");
	
	pCell->is_movable = false; 
	
	// if I'm dead, don't bother 
	if( phenotype.death.dead == true )
	{
		// the cell death functions don't automatically turn off custom functions, 
		// since those are part of mechanics. 
		// remove_all_adhesions( pCell ); 
		
		// Let's just fully disable now. 
		pCell->functions.custom_cell_rule = NULL; 
		pCell->functions.contact_function = NULL; 

		phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"]; 
		phenotype.secretion.secretion_rates[proinflammatory_cytokine_index] = 0; 
		return; 
	}	
	
	// this is now part of contact_function 
	/*
	// if I'm adhered to something ... 
	if( pCell->state.neighbors.size() > 0 )
	{
		// add the elastic forces 
		extra_elastic_attachment_mechanics( pCell, phenotype, dt );
	}
	*/
	return; 
}

void epithelium_submodel_setup( void )
{
	Cell_Definition* pCD;
	
	// set up any submodels you need 
	// viral replication 
	
	// receptor trafficking 
	simple_receptor_dynamics_model_setup(); // done 
	// viral replication 
	simple_internal_virus_model_setup();	
	// single-cell response  
	simple_internal_virus_response_model_setup(); 
 	
	// set up epithelial cells
		// set version info 
	epithelium_submodel_info.name = "epithelium model"; 
	epithelium_submodel_info.version = epithelium_submodel_version; 
		// set functions 
	epithelium_submodel_info.main_function = NULL; 
	epithelium_submodel_info.phenotype_function = epithelium_phenotype; 
	epithelium_submodel_info.mechanics_function = epithelium_mechanics; 
	
		// what microenvironment variables do you expect? 
	epithelium_submodel_info.microenvironment_variables.push_back( "virion" ); 
	epithelium_submodel_info.microenvironment_variables.push_back( "VTEST" ); 
	epithelium_submodel_info.microenvironment_variables.push_back( "interferon 1" ); 
	epithelium_submodel_info.microenvironment_variables.push_back( "pro-inflammatory cytokine" ); 
	epithelium_submodel_info.microenvironment_variables.push_back( "chemokine" ); 
		// what custom data do I need? 
	//epithelium_submodel_info.cell_variables.push_back( "something" ); 
		// register the submodel  
	epithelium_submodel_info.register_model();	
		// set functions for the corresponding cell definition 
	pCD = find_cell_definition( "lung epithelium" ); 
	pCD->functions.update_phenotype = epithelium_submodel_info.phenotype_function;
	pCD->functions.custom_cell_rule = epithelium_submodel_info.mechanics_function;
	pCD->functions.contact_function = epithelium_contact_function; 
	
	return;
}

void TCell_induced_apoptosis( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int apoptosis_index = phenotype.death.find_death_model_index( "Apoptosis" ); 
	static int debris_index = microenvironment.find_density_index( "debris" ); 
	static int proinflammatory_cytokine_index = microenvironment.find_density_index("pro-inflammatory cytokine");
	static int virion_index = microenvironment.find_density_index("virion");
	static int vtest_index = microenvironment.find_density_index("VTEST");
	
	
	
	if( pCell->custom_data["TCell_contact_time"] > pCell->custom_data["TCell_contact_death_threshold"] )
	{
		// make sure to get rid of all adhesions! 
		// detach all attached cells 
		// remove_all_adhesions( pCell ); 
		
		#pragma omp critical
		{
		//std::cout << "\t\t\t\t" << pCell << " (of type " << pCell->type_name <<  ") died from T cell contact" << std::endl; 
		}
		
		// induce death 
		pCell->start_death( apoptosis_index ); 
		
		pCell->phenotype.secretion.secretion_rates[proinflammatory_cytokine_index] = 0; 
		pCell->phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"]; 
		pCell->phenotype.molecular.fraction_released_at_death[virion_index] = parameters.doubles("virus_fraction_released_after_apoptosis");
		pCell->phenotype.molecular.fraction_released_at_death[vtest_index] = parameters.doubles("virus_fraction_released_after_apoptosis");
		
		pCell->functions.update_phenotype = NULL; 
	}
	
	return; 
}


void ROS_induced_apoptosis( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int apoptosis_index = phenotype.death.find_death_model_index( "Apoptosis" ); 
	static int ROS_index = microenvironment.find_density_index( "ROS" ); 
	double ROS_amount = pCell->nearest_density_vector()[ROS_index];
	static int debris_index = microenvironment.find_density_index( "debris" ); 
	static int proinflammatory_cytokine_index = microenvironment.find_density_index("pro-inflammatory cytokine");
	
	double epsilon_ROS = parameters.doubles("epsilon_ROS");
	
	double prob_apoptosis = ROS_amount/(ROS_amount+epsilon_ROS);
	
	if( UniformRandom() < prob_apoptosis )
	{
		// make sure to get rid of all adhesions! 
		// detach all attached cells 
		// remove_all_adhesions( pCell ); 
		
		#pragma omp critical
		{
		//std::cout << "\t\t\t\t" << pCell << " (of type " << pCell->type_name <<  ") died from ROS" << std::endl; 
		}
		
		// induce death 
		pCell->start_death( apoptosis_index ); 
		
		pCell->phenotype.secretion.secretion_rates[proinflammatory_cytokine_index] = 0; 
		pCell->phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"]; 
		
		pCell->functions.update_phenotype = NULL; 
	}
	
	return; 
}
