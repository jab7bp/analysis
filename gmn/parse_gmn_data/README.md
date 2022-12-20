#parse_gmn_data

This script parses down the raw data files for GMn down into one single file after applying a wide, global cut. 

The data/ROOTfiles are organized by:

*pass: 
  0 or 1 (as of Dec. 20, 2022). This specifies which sub-folder the data is in under /work/halla/sbs/sbs-gmn/pass#

*kine: 
  SBS kinematic setting: 4, 7, 8, 9, 11, 14
  Defines the kinematic setting --> used in locating the data/ROOTfile in: /work/halla/sbs/sbs-gmn/pass#/SBS#
  
*sbsfieldscale:
  SBS field in terms of percentage from 0 to 100%.
  This one is used a bit more generally:
    - Defines the run numbers for that particular kinematic that match the SBS field setting (looked up in beam_variables.h)
    - Defines the output file/sub-grouping for the parsed file.
    - Is used to know whether or not the dxdy plot should have a separated p or n spot, or just a single spot. 
  
*run_target:
  Defines the target used in the run/setting.
    - Used in locating the data/ROOTfile in: /work/halla/sbs/sbs-gmn/pass#/SBS#/TARGET
    - Used in determining fits for p and n in dxdy.
    
These are the primary inputs for the script. 

If these are defined then most everything else is pulled from beam_variables.h
For instance, runnum_vec ( a vector to hold the run numbers to be analyzed ) will be filled with all of the run numbers corresponding to kine, sbsfield, and run_target.
These numbers are stored in lookup_parsed_runnums in beam_variables.h

Another primary input is "parser_cut_vec".

*parser_cut_vec:
  Defines the cuts to be applied on the files. This is the parsing cut.
  the first 6 cuts are basic and can apply to ALL kinematic settings.
  The cuts that follow "bb.tr.n==1" are the cuts that are unique to each setting. 
  
  Before parsing the files you should run a script that plots the variables that you are interested in cutting. 
  For instance: 
    Pre-Shower: bb.ps.e
    E/p: (bb.sh.e + bb.ps.e)/(bb.tr.p[0])
    HCal Cluster Energy: sbs.hcal.e
    Total Energer (Shower + Pre-Shower): bb.sh.e + bb.ps.e
    W2 (from tree variable): e.kine.W2
    
  From your plots you can select the cuts you want to apply and enter them here.
  Keep them broad and wide enough that you are not losing real data!
  
  A parser_cut_vec is defined per kinematic. This is because each one gets a unique set of cuts and should be considered unique. 
 
The output file will be named following this pattern:
  gmn_parsed_{TARGET}_SBS{KINEMATIC}_mag{FIELD SETTING}.root
  
  
    
