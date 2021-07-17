#include "./internal_viral_dynamics.h" 

using namespace PhysiCell; 

std::string internal_virus_replication_version = "0.4.0"; 

Submodel_Information internal_viral_dynamics_info; 

void simple_internal_virus_model_setup(void)
{
	// set version
	internal_viral_dynamics_info.name = "internal viral replication dynamics"; 
	internal_viral_dynamics_info.version = internal_virus_replication_version; 
		// set functions 
	internal_viral_dynamics_info.main_function = NULL; 
	internal_viral_dynamics_info.phenotype_function = simple_internal_virus_model; 
	internal_viral_dynamics_info.mechanics_function = NULL; 
	
	// what microenvironment variables do you need 
	internal_viral_dynamics_info.microenvironment_variables.push_back( "virion" ); 	
	internal_viral_dynamics_info.microenvironment_variables.push_back( "interferon 1" ); 		

	// what custom data do I need? 
	internal_viral_dynamics_info.cell_variables.push_back( "VEn" ); 		
	internal_viral_dynamics_info.cell_variables.push_back( "Vnuc" ); 
		
	// submodel_registry.register_model( internal_viral_dynamics_info ); 
	internal_viral_dynamics_info.register_model();
	return;
}


void simple_internal_virus_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	if( phenotype.death.dead == true ) // if cell is dead don't replicate virus
	{ return; } 

	static int lung_epithelial_type = get_cell_definition( "lung epithelium" ).type; 
		// if not lung epithelium, do nothing 
	if( pCell->type != lung_epithelial_type )
	{ return; } 

	// Simple virus replication model
	static int IFN_index = microenvironment.find_density_index( "interferon 1" );
	simple_intracellular_replication_model(pCell, phenotype, dt);
			
	return;
}
void simple_intracellular_replication_model(  Cell* pCell, Phenotype& phenotype, double dt )
{	

	static int vtest_external = microenvironment.find_density_index( "VTEST" ); 
	static int proinflammatory_cytokine_index = microenvironment.find_density_index( "pro-inflammatory cytokine");
	static int apoptosis_index = pCell->phenotype.death.find_death_model_index( "apoptosis" ); 
	
	double Vconc = pCell->phenotype.molecular.internalized_total_substrates[vtest_external];
	double Vvoxel = microenvironment.mesh.voxels[1].volume;
	static double gamnuc = parameters.doubles("gamnuc");
	static double alpha = parameters.doubles("alpha")/Vvoxel;
	double v_rep = parameters.doubles("v_rep");
	double tau_rel = parameters.doubles("tau_rel");
	
	// cell isn't in an antiviral state and has enough intracellular virus
	if(pCell->custom_data["antiviral_state"]<0.5&&pCell->phenotype.molecular.internalized_total_substrates[vtest_external]*Vvoxel>v_rep) 
	{	
		pCell->phenotype.molecular.internalized_total_substrates[vtest_external] += gamnuc*Vconc*(1-Vconc/alpha)*dt;
		if(pCell->custom_data["eclipse_time"]<1)
		{pCell->custom_data["eclipse_time"] = PhysiCell_globals.current_time+tau_rel;}
		pCell->phenotype.death.rates[apoptosis_index] = parameters.doubles("infected_cell_death_rate");
	
	}
	else if(pCell->custom_data["antiviral_state"]>0.5) //Cell is in an antiviral state stop secreting cytokine and virions
	{
		pCell->phenotype.molecular.internalized_total_substrates[vtest_external] = 0;
		phenotype.secretion.secretion_rates[proinflammatory_cytokine_index] = 0;
	}
	// cell not in antiviral state but not sufficiently infected
	else if(pCell->phenotype.molecular.internalized_total_substrates[vtest_external]*Vvoxel<11) 
	{
			double prob_infection_recognition = 0.3; // probability at low MOI a cell realises it's infected
			if(UniformRandom()<prob_infection_recognition)
			{pCell->phenotype.molecular.internalized_total_substrates[vtest_external] = 0;}
			
	}
	pCell->custom_data["Vnuc"] = pCell->phenotype.molecular.internalized_total_substrates[vtest_external];	
	
	return;	
}
