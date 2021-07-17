/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL;  
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	int virion_index = microenvironment.find_density_index( "virion" ); 
	int VTEST_index = microenvironment.find_density_index( "VTEST" ); 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	// register the submodels 
	// (which ensures that the cells have all the internal variables they need) 
	
	Cell_Definition* pCD = find_cell_definition( "lung epithelium" ); 
	pCD->phenotype.molecular.fraction_released_at_death[virion_index] = 
		parameters.doubles("virus_fraction_released_at_death"); 
	pCD->phenotype.molecular.fraction_released_at_death[VTEST_index] = 
		parameters.doubles("virus_fraction_released_at_death"); 

	immune_submodels_setup();
	// receptor_dynamics_model_setup(); 
	// internal_virus_model_setup();
	// internal_virus_response_model_setup();
	epithelium_submodel_setup(); 

	submodel_registry.display( std::cout ); 
		
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	static int nV = microenvironment.find_density_index( "virion" ); 
	
	static int vtest_external = microenvironment.find_density_index( "VTEST" ); 
	
	choose_initialized_voxels();
	
	// create some cells near the origin
	
	Cell* pC;
	
	// hexagonal cell packing 
	Cell_Definition* pCD = find_cell_definition("lung epithelium"); 
	
	double cell_radius = pCD->phenotype.geometry.radius; 
	double spacing = 0.95 * cell_radius * 2.0; 
	
	double x_min = microenvironment.mesh.bounding_box[0] + cell_radius; 
	double x_max = microenvironment.mesh.bounding_box[3] - cell_radius; 

	double y_min = microenvironment.mesh.bounding_box[1] + cell_radius; 
	double y_max = microenvironment.mesh.bounding_box[4] - cell_radius; 
	
	double x = x_min; 
	double y = y_min; 
	
	double center_x = 0.5*( x_min + x_max ); 
	double center_y = 0.5*( y_min + y_max ); 
	
	double triangle_stagger = sqrt(3.0) * spacing * 0.5; 
	
	// find hte cell nearest to the center 
	double nearest_distance_squared = 9e99; 
	Cell* pNearestCell = NULL; 

	int density_virions_4 = 10;
	
	int large_bron = parameters.ints("large_bron");
	int small_bron = parameters.ints("small_bron");
	double density_virions = parameters.doubles("density_virion");
	
	int bronchcentre1_x = 1500-4000;
	int bronchcentre1_y = 2500-2500;	
	int bronchcentre2_x = 2500-4000;
	int bronchcentre2_y = 3500-2500;
	int bronchcentre3_x = 4000-4000;
	int bronchcentre3_y = 4000-2500;
	int bronchcentre4_x = 5500-4000;
	int bronchcentre4_y = 3000-2500;
	int bronchcentre5_x = 6500-4000;
	int bronchcentre5_y = 2200-2500;
	int bronchcentre6_x = 4422-4000;
	int bronchcentre6_y = 3811-2500;
	int bronchcentre7_x = 5953-4000;
	int bronchcentre7_y = 3439-2500;
	int bronchcentre8_x = 2298-4000;
	int bronchcentre8_y = 1588-2500;
	int bronchcentre9_x = 4869-4000;
	int bronchcentre9_y = 601-2500;
	int bronchcentre10_x = 5181-4000;
	int bronchcentre10_y = 2241-2500;
	int bronchcentre11_x = 1552-4000;
	int bronchcentre11_y = 999-2500;
	int bronchcentre12_x = 3935-4000;
	int bronchcentre12_y = 2692-2500;
	int bronchcentre13_x = 909-4000;
	int bronchcentre13_y = 1963-2500;
	int bronchcentre14_x = 3185-4000;
	int bronchcentre14_y = 1215-2500;
	int bronchcentre15_x = 6671-4000;
	int bronchcentre15_y = 4276-2500;
	int bronchcentre16_x = 7318-4000;
	int bronchcentre16_y = 2992-2500;
	int bronchcentre17_x = 1845-4000;
	int bronchcentre17_y = 4067-2500;		
		
	if( parameters.bools( "initial_condition_large_tissue_bronchiole") == true )
	{
		// first bronchiole radius 120
		std::vector<double> position = {0,0,0};
				
		int n = 0; 
		while( y < y_max )
		{
			while( x < x_max)
			{
				if((x-bronchcentre1_x)*(x-bronchcentre1_x)+(y-bronchcentre1_y)*(y-bronchcentre1_y)>large_bron*large_bron &&
				(x-bronchcentre2_x)*(x-bronchcentre2_x)+(y-bronchcentre2_y)*(y-bronchcentre2_y)>large_bron*large_bron &&
				(x-bronchcentre3_x)*(x-bronchcentre3_x)+(y-bronchcentre3_y)*(y-bronchcentre3_y)>large_bron*large_bron &&
				(x-bronchcentre4_x)*(x-bronchcentre4_x)+(y-bronchcentre4_y)*(y-bronchcentre4_y)>large_bron*large_bron&&
				(x-bronchcentre5_x)*(x-bronchcentre5_x)+(y-bronchcentre5_y)*(y-bronchcentre5_y)>large_bron*large_bron &&
				(x-bronchcentre6_x)*(x-bronchcentre6_x)+(y-bronchcentre6_y)*(y-bronchcentre6_y)>small_bron*small_bron&&
				(x-bronchcentre7_x)*(x-bronchcentre7_x)+(y-bronchcentre7_y)*(y-bronchcentre7_y)>small_bron*small_bron&&
				(x-bronchcentre8_x)*(x-bronchcentre8_x)+(y-bronchcentre8_y)*(y-bronchcentre8_y)>small_bron*small_bron&&
				(x-bronchcentre9_x)*(x-bronchcentre9_x)+(y-bronchcentre9_y)*(y-bronchcentre9_y)>small_bron*small_bron&&
				(x-bronchcentre10_x)*(x-bronchcentre10_x)+(y-bronchcentre10_y)*(y-bronchcentre10_y)>small_bron*small_bron&&
				(x-bronchcentre11_x)*(x-bronchcentre11_x)+(y-bronchcentre11_y)*(y-bronchcentre11_y)>small_bron*small_bron&&
				(x-bronchcentre12_x)*(x-bronchcentre12_x)+(y-bronchcentre12_y)*(y-bronchcentre12_y)>small_bron*small_bron&&
				(x-bronchcentre13_x)*(x-bronchcentre13_x)+(y-bronchcentre13_y)*(y-bronchcentre13_y)>small_bron*small_bron&&
				(x-bronchcentre14_x)*(x-bronchcentre14_x)+(y-bronchcentre14_y)*(y-bronchcentre14_y)>small_bron*small_bron&&
				(x-bronchcentre15_x)*(x-bronchcentre15_x)+(y-bronchcentre15_y)*(y-bronchcentre15_y)>small_bron*small_bron&&
				(x-bronchcentre16_x)*(x-bronchcentre16_x)+(y-bronchcentre16_y)*(y-bronchcentre16_y)>small_bron*small_bron&&
				(x-bronchcentre17_x)*(x-bronchcentre17_x)+(y-bronchcentre17_y)*(y-bronchcentre17_y)>small_bron*small_bron)
				{
					pC = create_cell( get_cell_definition("lung epithelium" ) ); 
					pC->assign_position( x,y, 0.0 );
				
					double dx = x - center_x;
					double dy = y - center_y; 
					
					double temp = dx*dx + dy*dy; 
					if( temp < nearest_distance_squared )
					{
						nearest_distance_squared = temp;
						pNearestCell = pC; 
					}
				}
				x += spacing; 
			}
			x = x_min; 
			
			n++; 
			y += triangle_stagger; 
			// in odd rows, shift 
			if( n % 2 == 1 )
			{
				x += 0.5 * spacing; 
			}
		}		
	}
	else
		{
		
		int n = 0; 
		while( y < y_max )
		{
			while( x < x_max)
			{
				pC = create_cell( get_cell_definition("lung epithelium" ) ); 
				pC->assign_position( x,y, 0.0 );
				
				double dx = x - center_x;
				double dy = y - center_y; 
				
				double temp = dx*dx + dy*dy; 
				if( temp < nearest_distance_squared )
				{
					nearest_distance_squared = temp;
					pNearestCell = pC; 
				}
				x += spacing; 
			}
			x = x_min; 
			
			n++; 
			y += triangle_stagger; 
			// in odd rows, shift 
			if( n % 2 == 1 )
			{
				x += 0.5 * spacing; 
			}
		}
	}
	
	
	if( parameters.bools( "initial_condition_large_tissue_bronchiole") == true )
	{
				
		std::vector<double> position = {0,0,0};
		
		for( int n=0 ; n < microenvironment.mesh.voxels.size() ; n++ )
		{
			std::vector<double> Vectpos = microenvironment.mesh.voxels[n].center;
			
			if((Vectpos[0]-bronchcentre1_x)*(Vectpos[0]-bronchcentre1_x)+(Vectpos[1]-bronchcentre1_y)*(Vectpos[1]-bronchcentre1_y)<large_bron*large_bron)//in bronchiole 1
			{microenvironment(n)[vtest_external] += density_virions;}
			else if((Vectpos[0]-bronchcentre2_x)*(Vectpos[0]-bronchcentre2_x)+(Vectpos[1]-bronchcentre2_y)*(Vectpos[1]-bronchcentre2_y)<large_bron*large_bron)//in bronchiole 2
			{microenvironment(n)[vtest_external] += density_virions;}
			else if((Vectpos[0]-bronchcentre3_x)*(Vectpos[0]-bronchcentre3_x)+(Vectpos[1]-bronchcentre3_y)*(Vectpos[1]-bronchcentre3_y)<large_bron*large_bron)//in bronchiole 3
			{microenvironment(n)[vtest_external] += density_virions;}
			else if((Vectpos[0]-bronchcentre4_x)*(Vectpos[0]-bronchcentre4_x)+(Vectpos[1]-bronchcentre4_y)*(Vectpos[1]-bronchcentre4_y)<large_bron*large_bron)//in bronchiole 4
			{microenvironment(n)[vtest_external] += density_virions;}
			else if((Vectpos[0]-bronchcentre5_x)*(Vectpos[0]-bronchcentre5_x)+(Vectpos[1]-bronchcentre5_y)*(Vectpos[1]-bronchcentre5_y)<large_bron*large_bron)//in bronchiole 5
			{microenvironment(n)[vtest_external] += density_virions;}
			else if((Vectpos[0]-bronchcentre6_x)*(Vectpos[0]-bronchcentre6_x)+(Vectpos[1]-bronchcentre6_y)*(Vectpos[1]-bronchcentre6_y)<small_bron*small_bron)//in bronchiole 6
			{microenvironment(n)[vtest_external] += density_virions;}
			else if((Vectpos[0]-bronchcentre7_x)*(Vectpos[0]-bronchcentre7_x)+(Vectpos[1]-bronchcentre7_y)*(Vectpos[1]-bronchcentre7_y)<small_bron*small_bron)//in bronchiole 7
			{microenvironment(n)[vtest_external] += density_virions;}
			else if((Vectpos[0]-bronchcentre8_x)*(Vectpos[0]-bronchcentre8_x)+(Vectpos[1]-bronchcentre8_y)*(Vectpos[1]-bronchcentre8_y)<small_bron*small_bron)//in bronchiole 8
			{microenvironment(n)[vtest_external] += density_virions;}
			else if((Vectpos[0]-bronchcentre9_x)*(Vectpos[0]-bronchcentre9_x)+(Vectpos[1]-bronchcentre9_y)*(Vectpos[1]-bronchcentre9_y)<small_bron*small_bron)//in bronchiole 9
			{microenvironment(n)[vtest_external] += density_virions;}
			else if((Vectpos[0]-bronchcentre10_x)*(Vectpos[0]-bronchcentre10_x)+(Vectpos[1]-bronchcentre10_y)*(Vectpos[1]-bronchcentre10_y)<small_bron*small_bron)//in bronchiole 10
			{microenvironment(n)[vtest_external] += density_virions;}
			else if((Vectpos[0]-bronchcentre11_x)*(Vectpos[0]-bronchcentre11_x)+(Vectpos[1]-bronchcentre11_y)*(Vectpos[1]-bronchcentre11_y)<small_bron*small_bron)//in bronchiole 11
			{microenvironment(n)[vtest_external] += density_virions;}
			else if((Vectpos[0]-bronchcentre12_x)*(Vectpos[0]-bronchcentre12_x)+(Vectpos[1]-bronchcentre12_y)*(Vectpos[1]-bronchcentre12_y)<small_bron*small_bron)//in bronchiole 12
			{microenvironment(n)[vtest_external] += density_virions;}
			else if((Vectpos[0]-bronchcentre13_x)*(Vectpos[0]-bronchcentre13_x)+(Vectpos[1]-bronchcentre13_y)*(Vectpos[1]-bronchcentre13_y)<small_bron*small_bron)//in bronchiole 13
			{microenvironment(n)[vtest_external] += density_virions;}
			else if((Vectpos[0]-bronchcentre14_x)*(Vectpos[0]-bronchcentre14_x)+(Vectpos[1]-bronchcentre14_y)*(Vectpos[1]-bronchcentre14_y)<small_bron*small_bron)//in bronchiole 14
			{microenvironment(n)[vtest_external] += density_virions;}
			else if((Vectpos[0]-bronchcentre15_x)*(Vectpos[0]-bronchcentre15_x)+(Vectpos[1]-bronchcentre15_y)*(Vectpos[1]-bronchcentre15_y)<small_bron*small_bron)//in bronchiole 15
			{microenvironment(n)[vtest_external] += density_virions;}
			else if((Vectpos[0]-bronchcentre16_x)*(Vectpos[0]-bronchcentre16_x)+(Vectpos[1]-bronchcentre16_y)*(Vectpos[1]-bronchcentre16_y)<small_bron*small_bron)//in bronchiole 16
			{microenvironment(n)[vtest_external] += density_virions;}
			else if((Vectpos[0]-bronchcentre17_x)*(Vectpos[0]-bronchcentre17_x)+(Vectpos[1]-bronchcentre17_y)*(Vectpos[1]-bronchcentre17_y)<small_bron*small_bron)//in bronchiole 17
			{microenvironment(n)[vtest_external] += density_virions;}
			
		}
		
	}
		
	// now place immune cells 
	
	initial_immune_cell_placement();
	
	return; 
}

std::vector<std::string> epithelium_coloring_function( Cell* pCell )
{
	std::vector<std::string> output( 4, "black" ); 
	
		double Vvoxel = microenvironment.mesh.voxels[1].volume;
		
	if( pCell->phenotype.death.dead == false )
	{
		double Vnuc = pCell->custom_data["Vnuc" ]*Vvoxel;
				
		double interpolation = 0; 
		if( Vnuc < 1 )
		{ interpolation = 0; } 
		if( Vnuc >= 1.0 && Vnuc < 10 )
		{ interpolation = 0.25;} 
		if( Vnuc >= 10.0 && Vnuc < 100 )
		{ interpolation = 0.5; } 
		if( Vnuc >= 100.0 && Vnuc < 1000 )
		{ interpolation = 0.75; } 
		if( Vnuc >= 1000.0 )
		{ interpolation = 1.0;}    	
		
		int red = (int) floor( 255.0 * interpolation ) ; 
		int green = red; 
		int blue = 255 - red; 

		char color [1024]; 
		sprintf( color, "rgb(%u,%u,%u)" , red,green,blue ); 

		output[0] = color; 
		output[2] = color; 
		output[3] = color; 
	}
	
	return output; 
}

std::vector<std::string> tissue_coloring_function( Cell* pCell )
{
	static int lung_epithelial_type = get_cell_definition( "lung epithelium" ).type; 
	
	static int CD8_Tcell_type = get_cell_definition( "CD8 Tcell" ).type; 
	static int Macrophage_type = get_cell_definition( "macrophage" ).type; 
	static int Neutrophil_type = get_cell_definition( "neutrophil" ).type; 
	static int DC_type = get_cell_definition( "DC" ).type; 
	static int CD4_Tcell_type = get_cell_definition( "CD4 Tcell" ).type; 
	
	// start with white 
	
	std::vector<std::string> output = {"white", "black", "white" , "white" };	
	
	if( pCell->phenotype.death.dead == true )
	{
		if( pCell->type != lung_epithelial_type )
		{
			output[0] = parameters.strings("apoptotic_immune_color");		
			output[2] = output[0]; 		
			output[3] = output[0]; 	
			return output; 
		}

		output[0] = parameters.strings("apoptotic_epithelium_color");	
		output[2] = output[0];
		output[3] = output[0];
		return output; 
	}

	if( pCell->phenotype.death.dead == false && pCell->type == lung_epithelial_type )
	{
		if(pCell->custom_data["antiviral_state"]>0.5)
		{
			output[0] = "blue";	
			output[2] = "blue";
			output[3] = "blue";	
			return output; 
		}
		// color by virion 
		output = epithelium_coloring_function(pCell); 
		return output; 
	}
	
	if( pCell->phenotype.death.dead == false && pCell->type == CD8_Tcell_type )
	{
		output[0] = parameters.strings("CD8_Tcell_color");  
		output[2] = output[0];
		output[3] = output[0];
		return output; 
	}
	
	// (Adrianne) adding CD4 T cell colouring
	if( pCell->phenotype.death.dead == false && pCell->type == CD4_Tcell_type )
	{
		output[0] = parameters.strings("CD4_Tcell_color");  
		output[2] = output[0];
		output[3] = output[0];
		return output; 
	}

	if( pCell->phenotype.death.dead == false && pCell->type == Macrophage_type )
	{
		std::string color = parameters.strings("Macrophage_color");  
		if( pCell->custom_data["activated_immune_cell" ] > 0.5 )
		{ color = parameters.strings("activated_macrophage_color"); }
		
		// (Adrianne) added colours to show when macrophages are exhausted and when they are hyperactivated
		if( pCell->phenotype.volume.total> pCell->custom_data["threshold_macrophage_volume"] )// macrophage exhausted
		{ color = parameters.strings("exhausted_macrophage_color"); }
		else if( pCell->custom_data["ability_to_phagocytose_infected_cell"] == 1)// macrophage has been activated to kill infected cells by T cell
		{ color = parameters.strings("hyperactivated_macrophage_color"); }
		
		output[0] = color; 
		output[2] = output[0];
		output[3] = output[0];
		return output; 
	}

	if( pCell->phenotype.death.dead == false && pCell->type == Neutrophil_type )
	{
		output[0] = parameters.strings("Neutrophil_color");  
		output[2] = output[0];
		output[3] = output[0];
		return output; 
	}
	
	//(Adrianne) adding colour for DCs
	if( pCell->phenotype.death.dead == false && pCell->type == DC_type )
	{
		std::string color = parameters.strings("DC_color");  
		if( pCell->custom_data["activated_immune_cell" ] > 0.5 )
		{ color = parameters.strings("activated_DC_color"); }
	
		output[0] = color; 
		output[2] = output[0];
		output[3] = output[0];
		return output; 
	}

	return output; 
}

bool Write_SVG_circle_opacity( std::ostream& os, double center_x, double center_y, double radius, double stroke_size, 
                       std::string stroke_color , std::string fill_color , double opacity )
{
 os << "  <circle cx=\"" << center_x << "\" cy=\"" << center_y << "\" r=\"" << radius << "\" stroke-width=\"" << stroke_size 
    << "\" stroke=\"" << stroke_color << "\" fill=\"" << fill_color 
	<< "\" fill-opacity=\"" << opacity << "\"/>" << std::endl; 
 return true; 
}


//
void SVG_plot_virus( std::string filename , Microenvironment& M, double z_slice , double time, std::vector<std::string> (*cell_coloring_function)(Cell*) )
{
	static double X_lower = M.mesh.bounding_box[0];
	static double X_upper = M.mesh.bounding_box[3];
 
	static double Y_lower = M.mesh.bounding_box[1]; 
	static double Y_upper = M.mesh.bounding_box[4]; 

	static double plot_width = X_upper - X_lower; 
	static double plot_height = Y_upper - Y_lower; 

	static double font_size = 0.025 * plot_height; // PhysiCell_SVG_options.font_size; 
	static double top_margin = font_size*(.2+1+.2+.9+.5 ); 
	
	static double epithelial_opacity = parameters.doubles("epithelial_opacity");
	static double non_epithelial_opacity = parameters.doubles("non_epithelial_opacity"); 

	// open the file, write a basic "header"
	std::ofstream os( filename , std::ios::out );
	if( os.fail() )
	{ 
		std::cout << std::endl << "Error: Failed to open " << filename << " for SVG writing." << std::endl << std::endl; 

		std::cout << std::endl << "Error: We're not writing data like we expect. " << std::endl
		<< "Check to make sure your save directory exists. " << std::endl << std::endl
		<< "I'm going to exit with a crash code of -1 now until " << std::endl 
		<< "you fix your directory. Sorry!" << std::endl << std::endl; 
		exit(-1); 
	} 
	
	Write_SVG_start( os, plot_width , plot_height + top_margin );

	// draw the background 
	Write_SVG_rect( os , 0 , 0 , plot_width, plot_height + top_margin , 0.002 * plot_height , "white", "white" );

	// write the simulation time to the top of the plot
 
	char* szString; 
	szString = new char [1024]; 
 
	int total_cell_count = all_cells->size(); 
 
	double temp_time = time; 

	std::string time_label = formatted_minutes_to_DDHHMM( temp_time ); 
 
	sprintf( szString , "Current time: %s, z = %3.2f %s", time_label.c_str(), 
		z_slice , PhysiCell_SVG_options.simulation_space_units.c_str() ); 
	Write_SVG_text( os, szString, font_size*0.5,  font_size*(.2+1), 
		font_size, PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() );
	sprintf( szString , "%u agents" , total_cell_count ); 
	Write_SVG_text( os, szString, font_size*0.5,  font_size*(.2+1+.2+.9), 
		0.95*font_size, PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() );
	
	delete [] szString; 


	// add an outer "g" for coordinate transforms 
	
	os << " <g id=\"tissue\" " << std::endl 
	   << "    transform=\"translate(0," << plot_height+top_margin << ") scale(1,-1)\">" << std::endl; 
	   
	// prepare to do mesh-based plot (later)
	
	double dx_stroma = M.mesh.dx; 
	double dy_stroma = M.mesh.dy; 
	
	os << "  <g id=\"ECM\">" << std::endl; 
  
	int ratio = 1; 
	double voxel_size = dx_stroma / (double) ratio ; 
  
	double half_voxel_size = voxel_size / 2.0; 
	double normalizer = 78.539816339744831 / (voxel_size*voxel_size*voxel_size); 

	os << "  </g>" << std::endl; 
	
	static Cell_Definition* pEpithelial = find_cell_definition( "lung epithelium" ); 
 
	// plot intersecting epithelial cells 
	os << "  <g id=\"cells\">" << std::endl; 
	for( int i=0 ; i < total_cell_count ; i++ )
	{
		Cell* pC = (*all_cells)[i]; // global_cell_list[i]; 
  
		static std::vector<std::string> Colors; 
		if( fabs( (pC->position)[2] - z_slice ) < pC->phenotype.geometry.radius 
			&& pC->type == pEpithelial->type )
		{
			double r = pC->phenotype.geometry.radius ; 
			double rn = pC->phenotype.geometry.nuclear_radius ; 
			double z = fabs( (pC->position)[2] - z_slice) ; 
   
			Colors = cell_coloring_function( pC ); 

			os << "   <g id=\"cell" << pC->ID << "\">" << std::endl; 
  
			// figure out how much of the cell intersects with z = 0 
   
			double plot_radius = sqrt( r*r - z*z ); 

			Write_SVG_circle_opacity( os, (pC->position)[0]-X_lower, (pC->position)[1]-Y_lower, 
				plot_radius , 0.5, Colors[1], Colors[0] , epithelial_opacity ); 
			os << "   </g>" << std::endl;
		}
	}
	
	// plot intersecting non=epithelial cells 
	for( int i=0 ; i < total_cell_count ; i++ )
	{
		Cell* pC = (*all_cells)[i]; // global_cell_list[i]; 
  
		static std::vector<std::string> Colors; 
		if( fabs( (pC->position)[2] - z_slice ) < pC->phenotype.geometry.radius 
			&& pC->type != pEpithelial->type )
		{
			double r = pC->phenotype.geometry.radius ; 
			double rn = pC->phenotype.geometry.nuclear_radius ; 
			double z = fabs( (pC->position)[2] - z_slice) ; 
   
			Colors = cell_coloring_function( pC ); 

			os << "   <g id=\"cell" << pC->ID << "\">" << std::endl; 
  
			// figure out how much of the cell intersects with z = 0 
   
			double plot_radius = sqrt( r*r - z*z ); 

			Write_SVG_circle_opacity( os, (pC->position)[0]-X_lower, (pC->position)[1]-Y_lower, 
				plot_radius , 0.5, Colors[1], Colors[0] , non_epithelial_opacity ); 
			os << "   </g>" << std::endl;
		}
		
	}
	
	os << "  </g>" << std::endl; 
	
	os << " </g>" << std::endl; 
 
 
	double bar_margin = 0.025 * plot_height; 
	double bar_height = 0.01 * plot_height; 
	double bar_width = PhysiCell_SVG_options.length_bar; 
	double bar_stroke_width = 0.001 * plot_height; 
	
	std::string bar_units = PhysiCell_SVG_options.simulation_space_units; 
	// convert from micron to mm
	double temp = bar_width;  

	if( temp > 999 && std::strstr( bar_units.c_str() , PhysiCell_SVG_options.mu.c_str() )   )
	{
		temp /= 1000;
		bar_units = "mm";
	}
	// convert from mm to cm 
	if( temp > 9 && std::strcmp( bar_units.c_str() , "mm" ) == 0 )
	{
		temp /= 10; 
		bar_units = "cm";
	}
	
	szString = new char [1024];
	sprintf( szString , "%u %s" , (int) round( temp ) , bar_units.c_str() );
 
	Write_SVG_rect( os , plot_width - bar_margin - bar_width  , plot_height + top_margin - bar_margin - bar_height , 
		bar_width , bar_height , 0.002 * plot_height , "rgb(255,255,255)", "rgb(0,0,0)" );
	Write_SVG_text( os, szString , plot_width - bar_margin - bar_width + 0.25*font_size , 
		plot_height + top_margin - bar_margin - bar_height - 0.25*font_size , 
		font_size , PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() ); 
	
	delete [] szString; 

	// plot runtime 
	szString = new char [1024]; 
	RUNTIME_TOC(); 
	std::string formatted_stopwatch_value = format_stopwatch_value( runtime_stopwatch_value() );
	Write_SVG_text( os, formatted_stopwatch_value.c_str() , bar_margin , top_margin + plot_height - bar_margin , 0.75 * font_size , 
		PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() );
	delete [] szString; 

	// draw a box around the plot window
	Write_SVG_rect( os , 0 , top_margin, plot_width, plot_height , 0.002 * plot_height , "rgb(0,0,0)", "none" );
	
	// close the svg tag, close the file
	Write_SVG_end( os ); 
	os.close();
 
	return; 
}


