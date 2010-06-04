/* 
   Copyright (C) 1997 1998 1999 2000 Dr. Preston B. Moore

Dr. Preston B. Moore
Associate Director, Center for Molecular Modeling (CMM)
University of Pennsylvania, Department of Chemistry, Box 188 
231 S. 34th St. Philadelphia, PA 19104-6323 USA
EMAIL: moore@cmm.chem.upenn.edu  
WWW: http://www.cmm.upenn.edu/~moore

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

/* routine to set up the keywords for the simulation inputs */

#include "md.h"

#define NUM_DICT_SIM (84)

void set_sim_keyword(int *num_dict,WORD **dict)
{
  *num_dict = NUM_DICT_SIM;
  *dict = (WORD *)cmalloc(NUM_DICT_SIM*sizeof(WORD));
  
  strcpy((*dict)[0],"calculation");   
  strcpy((*dict)[1],"ensemble");
  strcpy((*dict)[2],"ntime");         
  strcpy((*dict)[3],"istart");
  strcpy((*dict)[4],"iperd");         
  strcpy((*dict)[5],"dt");
  strcpy((*dict)[6],"t_ext");         
  strcpy((*dict)[7],"p_ext");
  strcpy((*dict)[8],"nlen");              
  strcpy((*dict)[9],"neighbor_list");

  strcpy((*dict)[10],"iwr_screen");   
  strcpy((*dict)[11],"iwr_dump");
  strcpy((*dict)[12],"iwr_inst");     
  strcpy((*dict)[13],"iwr_confp");
  strcpy((*dict)[14],"iwr_confv");    
  strcpy((*dict)[15],"resamvel");
  strcpy((*dict)[16],"ivec_init");    
  strcpy((*dict)[17],"ivec_update");
  strcpy((*dict)[18],"iwr_system_vectors");
  strcpy((*dict)[19],"npoints");

  strcpy((*dict)[20],"sim_name");     
  strcpy((*dict)[21],"dpo_file");
  strcpy((*dict)[22],"dpi_file");     
  strcpy((*dict)[23],"ham_file");
  strcpy((*dict)[24],"pos_file");     
  strcpy((*dict)[25],"vel_file");
  strcpy((*dict)[26],"eng_file");     
  strcpy((*dict)[27],"mol_set_file");
  strcpy((*dict)[28],"inter_file");   
  strcpy((*dict)[29],"bond_file");

  strcpy((*dict)[30],"bend_file");    
  strcpy((*dict)[31],"torsion_file");
  strcpy((*dict)[32],"onefour_file"); 
  strcpy((*dict)[33],"nma_file");
  strcpy((*dict)[34],"vec_file");     
  strcpy((*dict)[35],"vec_conf");
  strcpy((*dict)[36],"sys_vec_file");             
  strcpy((*dict)[37],"mol_vec_file");
  strcpy((*dict)[38],"spectrum_file");             
  strcpy((*dict)[39],"force_file");

  strcpy((*dict)[40],"nres_lrf");     
  strcpy((*dict)[41],"nres_intra");
  strcpy((*dict)[42],"len_nhc");      
  strcpy((*dict)[43],"nbar");
  strcpy((*dict)[44],"rheal_res");    
  strcpy((*dict)[45],"skin");
  strcpy((*dict)[46],"tau_nhc");      
  strcpy((*dict)[47],"tau_vol");
  strcpy((*dict)[48],"alp_ewald");    
  strcpy((*dict)[49],"kmax");

  strcpy((*dict)[50],"iseed");        
  strcpy((*dict)[51],"shift");
  strcpy((*dict)[52],"units");        
  strcpy((*dict)[53],"ntable");
  strcpy((*dict)[54],"scale_14");     
  strcpy((*dict)[55],"scale_14_e");
  strcpy((*dict)[56],"rcute_max");    
  strcpy((*dict)[57],"rcute_min");
  strcpy((*dict)[58],"rcute_resp");   
  strcpy((*dict)[59],"nstates");      

  strcpy((*dict)[60],"num_real");     
  strcpy((*dict)[61],"num_imag");
  strcpy((*dict)[62],"freq_min");     
  strcpy((*dict)[63],"freq_max");
  strcpy((*dict)[64],"num_subst");    
  strcpy((*dict)[65],"num_vecsub");
  strcpy((*dict)[66],"rescalevel");
  strcpy((*dict)[67],"scale_subspace");
  strcpy((*dict)[68],"ncell_div");
  strcpy((*dict)[69],"long_rc");

  strcpy((*dict)[70],"iwr_conff");
  strcpy((*dict)[71],"min_type");
  strcpy((*dict)[72],"linkcell");
  strcpy((*dict)[73],"max_exclude");
  strcpy((*dict)[74],"extern_field");
  strcpy((*dict)[75],"extern_file");
  strcpy((*dict)[76],"min_tolerance");
  strcpy((*dict)[77],"wall_clock");
  strcpy((*dict)[78],"eval_file");
  strcpy((*dict)[79],"nres_tors");
  strcpy((*dict)[80],"nres_nhc");
  strcpy((*dict)[81],"yoshidas");
  strcpy((*dict)[82],"scaleeps");
  strcpy((*dict)[83],"scalecharge");
}
