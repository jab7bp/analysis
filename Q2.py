from ROOT import *
from beam_variables import *
import numpy as np

def plot_Q2(runnum, rootFiles, save):

	global hist_Q2, hQ2_stack, Q2_stack_can, hQ2_mean

	Q, Q2, hQ2, hQ2_mean, hQ2_max, num_events = [], [], [], [], [], []

	Me = (0.0005)*0.938

	bins = 200

	#Canvas and Information to put on plot
	legend = TLegend(0.1, 0.6, 0.38, 0.9)
	legend.SetHeader("Beam current and # of Events in Peak Region", "C")
	legend.SetTextSize(0.02)

	#Run info to put in plot
	run_info = TPaveText(0.55,0.75,0.89,0.89,"NDC")
	run_info.AddText("Target: " + TargetFromRun(int(runnum[0])))
	run_info.AddText("Required hits on track: 3+ ")
	run_info.AddText("Shower & Pre-shower cuts: " + str(CutFromRun(int(runnum[0]))[0]) + " < E/p < " + str(CutFromRun(int(runnum[0]))[1]))
	run_info.SetTextSize(.025)
	run_info.SetFillColor(0)

	stack_title = "Momentum transfer, Q-squared"
	Q2_stack_can = TCanvas("Q2_stacked_can", stack_title)
	hQ2_stack = THStack("hQ2_stack", stack_title)


	lineColor = 1	
	
###############################################################################################
###******* WHERE THE WORK IS DONE -------- CALCULATING Q-Squared Momentum Transfer**********###
###############################################################################################

	for run in runnum:
		E = BeamEnergyFromRun(int(run))

		index = runnum.index(run)
		min_bin = qBinFromRun(int(run))[0]
		max_bin = qBinFromRun(int(run))[1]

		cut_min = CutFromRun(int(run))[0]
		cut_max = CutFromRun(int(run))[1]


		tCh_Q2 = TChain("T")
		tCh_Q2.Add(rootFiles[index])

		hist_Q2 = TH1F('hist_Q2_' + run, 'Q-Squared Transfer', bins, min_bin, max_bin)

		cut1 = TCut("bb.gem.track.nhits[0] > 3")
		cut2 = TCut("(bb.sh.e + bb.ps.e)/bb.tr.p > " + str(cut_min) + " && (bb.ps.e + bb.ps.e)/bb.tr.p < " + str(cut_max))

		all_cuts = cut1 + cut2

		###Variables: E = Beam Energy; bb.tr.p[0] = scattered electron energy; bb.tr.pz[0] = scattered electron energy in "z"
		###Mp = mass of proton; angle = (bb.tr.pz[0]/bb.tr.p[0])

		eq_Q2 = "2*%f*bb.tr.p[0]*(1-cos(bb.tr.pz[0]/bb.tr.p[0]))>> hist_Q2_" + run


		tCh_Q2.Draw(eq_Q2 % E, all_cuts)

		num_events.append(hist_Q2.Integral(hist_Q2.GetXaxis().FindBin(min_bin), hist_Q2.GetXaxis().FindBin(max_bin)))

		#PLOT DISPLAY STUFF#
		if(lineColor == 5):
			lineColor+=1
		hist_Q2.SetLineColor(lineColor)
		lineColor+=1

		legend.AddEntry(hist_Q2, "Current: " + str(CurrentFromRun(int(run))) + " \mu A; " + str(int(num_events[index])) + ' ' + "events" )

		##Store calculated variables in arrays##
		hQ2.append(hist_Q2)
		hQ2_stack.Add(hist_Q2)
		hQ2_mean.append(hist_Q2.GetMean(1))
		hQ2_max.append(hQ2[index].GetMaximum())


	
	hQ2_stack.Draw("nostack")

	#Line info to show peak cuts
	q2_mean_line = TLine(np.mean(hQ2_mean), 0, np.mean(hQ2_mean), max(hQ2_max)*1.05)
	q2_mean_line.SetLineStyle(6)
	legend.AddEntry(q2_mean_line, "Average Q-Squared = " + str(np.round(np.mean(hQ2_mean) , 3)) + " (GeV/c)^{2}", "L")
	legend.SetTextSize(0.02)

	
	hQ2_stack.Draw("nostack")
	hQ2_stack.GetXaxis().SetLabelOffset(0.02)
	legend.Draw("same")
	q2_mean_line.Draw("same")
	run_info.Draw()
	hQ2_stack.GetXaxis().SetTitle("Momentum Transfer, Q^{2} [(GeV/c)^{2}]")
	hQ2_stack.GetXaxis().SetTitleOffset(1.3)
	hQ2_stack.GetYaxis().SetTitle("Counts [N]")

	Q2_stack_can.Update()

	if(save):
		Q2_stack_can.Print("./momentum-transfers_plot.pdf")
		
	return 1