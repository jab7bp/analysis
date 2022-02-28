from ROOT import *
from beam_variables import *
import numpy as np

def plot_Ep(runnum, rootFiles):
	global hEp, hEp_stack, hEp_stack_can

	hEp = []
	hEp_stack_title = "E/p plots"
	hEp_stack_can = TCanvas("hEp_stack_can", hEp_stack_title)
	hEp_stack = THStack("hEp_stack", hEp_stack_title)

	bins = 200
	min_bin = 0.5 
	max_bin = 1.25

	lineColor = 1

	#Run info to put in plot
	run_info = TPaveText(0.55,0.75,0.89,0.89,"NDC")
	#run_info.AddText("Momentum transfer, Q^{2} = " + str(np.round(q_squared, 3)) + " (GeV/c)^{2}")
	run_info.AddText("Target: " + TargetFromRun(int(runnum[0])))
	run_info.AddText("Energy: " + BeamEnergyFromRun(int(runnum[0])))
	run_info.SetTextSize(.025)
	run_info.SetFillColor(0)


	for run in runnum:

		index = runnum.index(run)

		tCh_Ep = TChain("T")
		tCh_Ep.Add(rootFiles[index])

		hist_Ep = TH1F('hist_Ep_' + run , 'Energy over momentum', bins, min_bin, max_bin)
		
		cut_ps = TCut("bb.ps.e>0.1")

		eq = "((bb.sh.e + bb.ps.e)/(bb.tr.p))>> hist_Ep_" + run

		tCh_Ep.Draw(eq, cut_ps)

		#PLOT DISPLAY STUFF#
		if(lineColor == 5):
			lineColor+=1
		hist_Ep.SetLineColor(lineColor)
		lineColor+=1

		hEp_stack.Add(hist_Ep)


	hEp_stack.Draw("nostack")
	hEp_stack_can.Update()

	return 1