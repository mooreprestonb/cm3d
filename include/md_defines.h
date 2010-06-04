/* 
   Copyright (C) 1997-2001 Dr. Preston B. Moore

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

#define MAXWORDLEN 256
#define MAXLINELEN 512

#define META_KEY_CHAR '~'
#define KEY_WORD_CHAR '\\'
#define NEWLINE       '\n'
#define NILL          '\0'
#define WILDCARD      "X"

#define META_KEY_INTER  "inter_parm"
#define META_KEY_BOND   "bond_parm"
#define META_KEY_BONDX  "bondx_parm"
#define META_KEY_BEND   "bend_parm"
#define META_KEY_BENDQUART   "bendquart_parm"
#define META_KEY_TORS   "torsion_parm"
#define META_KEY_ONFO   "onefour_parm"
#define MAXPOWER_TORS   7

#define POT_TYPE_NULL   "null"
#define POT_TYPE_LJ     "lennard-jones"
#define POT_TYPE_LJ64   "lennard-jones64"
#define POT_TYPE_LJ96   "lennard-jones96"
#define POT_TYPE_LJ124  "lennard-jones124"
#define POT_TYPE_WILL   "williams"
#define POT_TYPE_AZIZ   "aziz"
#define POT_TYPE_HYDBND "hyd-bond12_10"
#define POT_TYPE_HARM   "harm"
#define POT_TYPE_POWER  "power"
#define POT_TYPE_BONDX  "bond_cross"
#define POT_TYPE_MORSE  "morse"
#define POT_TYPE_COSINE "cosine"
#define POT_TYPE_QUARTIC "quartic"
#define POT_TYPE_TAB    "tabulated"
#define POT_TYPE_DZUG   "dzugutov"

#define MOL_DEF "mol_def"

#define IMODES 9
#define ITMAX 3
#define NP_BOND 5 /* degree of polynomial for bond +1 ie degree 4 will have 5*/
#define NP_BEND 5 /* degree of polynomial for bend +1  */
