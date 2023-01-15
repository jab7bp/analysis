This directory is for the scripts that will generate a Neighbor Ratio plot for all APVs in a run.

The script relies on the prime GMn replay. No strip variables are need in output file. The strip comparison is performed on the replay-level.
The only criteria is that there be two sets of files: XTALK_ON and XTALK_OFF.

Another interesting and annoying aspect of this is that due to a bug in one of the codes, a middle step requires analyzing the replay files into ALL_MODULES AND LAST_MODULE. For some reason, in the ALL_MODULES version the very last module is not properly analyzed. So, we must run a separate analysis where we only look at the final, LAST, module.

#1. MAIN REPLAY

  ##submit-xtalk-jobs.sh

    Description: This produces a controlled set of replays to be used down the line for the Crosstalk Ratio analysis. Technically, you could use any set of output replay files but this one pulls in from the DB_XTALK dir and therefore can pull certain criteria from that database.

  ##run-xtalk.sh

   Called by submit-xtalk-jobs.sh
   
  -------------------------------------------------------------------
  ##Inputs:
  
  ##submit-xtalk-jobs runnum maxsegments
  
  **Output files: e1209019_fullreplay_$RUNNUM$_stream0_seg#_#_xtalk_replay_Histogramming_WITH_corr.root)
