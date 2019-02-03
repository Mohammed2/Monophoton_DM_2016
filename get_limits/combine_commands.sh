# -------------------------------
# --- Installation (one-time) ---
# -------------------------------
# First, go to desired base directory for installation.
# Then run these steps:
export SCRAM_ARCH=slc6_amd64_gcc530
cmsrel CMSSW_8_1_0
cd CMSSW_8_1_0/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v7.0.7
scramv1 b clean
scramv1 b
#  Install the portion of CombineHarvester needed to run impacts
#   Starting from CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit
cd ../..
bash <(curl -s https://raw.githubusercontent.com/cms-analysis/CombineHarvester/master/CombineTools/scripts/sparse-checkout-ssh.sh)
scram b -j 8


# ----------------------------------------------
# --- How to run everything in HiggsCombine ---
# ----------------------------------------------

# Go to CMSSW_8_1_0/src/HiggsCombine
#  Set up working area here
cmsenv

# Import needed histograms and transfer factors
# Replace /nfs_scratch/jjbuchanan/CMSSW_8_0_26_patch1/src/ControlRegionStudy_clean
# with your working directory
./copy_histos.sh

## Make workspace(s) ##
#  Choose which workspace(s) to make
#  Set connectWZ to true or false, based on whether the Wgamma and Zgamma estimates should be connected
#  VERY IMPORTANT: Record the text output at the end --- this is used to scale limits obtained with AsymptoticLimits
root -l -q -b createWorkspaces_Pt_HaloFit.C
root -l -q -b createWorkspaces_Pt_HaloFit_DMinterpolation.C

## Make data card ##
combineCards.py -S SA=datacard_signal_above0p5_yesZW_halo.txt SB=datacard_signal_below0p5_halo.txt CR_ME=datacard_monoele.txt CR_MM=datacard_monomu.txt CR_DM=datacard_dimu.txt CR_DE=datacard_diele.txt > comb_card.txt
combineCards.py -S SA=datacard_signal_above0p5_yesZW_halo_DMinterpolation.txt SB=datacard_signal_below0p5_halo_DMinterpolation.txt CR_ME=datacard_monoele.txt CR_MM=datacard_monomu.txt CR_DM=datacard_dimu.txt CR_DE=datacard_diele.txt > comb_card_DMinterpolation.txt

## Get limits and plots for one signal sample ##
#  Open comb_card.txt and replace workspace.root with a specific workspace name
#  Remove -t -1 for unblinded limits
combine -M AsymptoticLimits comb_card.txt -t -1
combine -M AsymptoticLimits comb_card_DMinterpolation.txt -t -1
# -- Run fit diagnostics --
text2workspace.py comb_card.txt
combine comb_card.root -M FitDiagnostics --saveShapes --saveWithUncertainties  --saveNormalizations
# -- Get total pre- and post-fit yields --
python ../test/mlfitNormsToText.py fitDiagnostics.root --uncertainties >& pre_and_post_fit_yields.txt
# -- Get magnitude of nuisance shifts --
python ../test/diffNuisances.py fitDiagnostics.root -g absoluteNuisanceDifferences.root
# -- Get pulls --
python ../test/diffNuisances.py fitDiagnostics.root -g pulls.root --pullDef relDiffAsymErrs
# -- Make phoPt plots --
root -l -q -b prefit_znng_SA_plotter.C
root -l -q -b prefit_znng_SB_plotter.C
root -l -q -b prefit_weng_plotter.C
root -l -q -b prefit_wmng_plotter.C
root -l -q -b prefit_zeeg_plotter.C
root -l -q -b prefit_zmmg_plotter.C
#  * To make background-only postfit plots, make sure shapes_fit_b (and not shapes_fit_s) is used everywhere
#    and make sure "SA_b_" (and not "SA_s_") is added to plotTitle.
#    For signal+background fits, do the opposite.
root -l -q -b postfit_znng_SA_plotter.C
root -l -q -b postfit_znng_SB_plotter.C
root -l -q -b postfit_weng_plotter.C
root -l -q -b postfit_wmng_plotter.C
root -l -q -b postfit_zeeg_plotter.C
root -l -q -b postfit_zmmg_plotter.C
# -- Impacts --
#  * CombineHarvester needs to be installed, according to the instructions above
#  * Remove -t -1 for unblinded impacts
text2workspace.py comb_card.txt -m 125 -o comb_card.root
combineTool.py -M Impacts -d comb_card.root -m 125 --robustFit 1 --expectSignal=0 -t -1 --doInitialFit
combineTool.py -M Impacts -d comb_card.root -m 125 --robustFit 1 --expectSignal=0 -t -1 --doFits --parallel 24
combineTool.py -M Impacts -d comb_card.root -m 125 -o impacts.json_comb_card
plotImpacts.py -i impacts.json_comb_card -o impacts_comb_card
#  * Signal significance:
#  * Must create workspaces with actual dark matter expected yield to be able to get actual significance
#  * Remove -t -1 for unblinded significance
combine -M Significance comb_card.txt -t -1 --expectSignal=1


## Running over every signal sample ##
#  Be sure template_card_Pt_HaloFit.txt is in the working directory and up to date
#   template_card_Pt_HaloFit.txt should be a duplicate of comb_card.txt with workspace.root as the workspace
#  The points to run over should be uncommented in createWorkspaces_Pt_HaloFit.C
root -l -q -b createWorkspaces_Pt_HaloFit.C
  #  VERY IMPORTANT: Copy the final text output to createDatacards_Pt_HaloFit.py
  #  as well as print_true_limits_Pt.py
  #   This contains the scaling factors needed for the final limits
python createDatacards_Pt_HaloFit.py
  #  The next commands could take awhile to finish
nohup ./limit_commands_LO_Pt.sh &
nohup ./limit_commands_NLO_Pt.sh &
nohup ./limit_commands_ADD_Pt.sh &
nohup ./limit_commands_EWK_Pt.sh &
nohup ./limit_commands_NLOfromInterpolation_Pt.sh &
  #  Gathers the limit results in a plotter-friendly format
python print_true_limits_Pt.py
  #  Copy resulting .txt files to wherever 2D limit plots are made
