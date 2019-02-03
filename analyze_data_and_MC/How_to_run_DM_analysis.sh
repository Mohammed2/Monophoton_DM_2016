## Set up the basics ##
# Go to CMSSW 8_0_26_patch1/src
# Set up working directory here
# Make sure rootcom and MakeCondorFiles.csh are in the working directory
cmsenv
voms-proxy-init --voms cms --valid 168:00


## Finding acceptances ##
# Run over dark matter samples
python analyze_central_signals.py
python analyze_DM_LO.py
python analyze_DM_NLO.py
# Run over data and MC background samples
# The commands given here for submit_*.sh are indicative of how they can be run,
# but it is often preferable to run partial fragments of the full code.
./submit_data.sh
./submit_pdfscale.sh
./submit_JESPES.sh


## Putting it all together ##
# Everything in "Finding acceptances" must be finished
hadd -f ZnnG_data_above0p5_all.root ZnnG_data*000*_above0p5.root
hadd -f ZnnG_data_below0p5_all.root ZnnG_data*000*_below0p5.root
hadd -f ZeeG_data_all.root ZeeG_data*000*.root
hadd -f ZmmG_data_all.root ZmmG_data*000*.root
hadd -f WenG_data_all.root WenG_data*000*.root
hadd -f WmnG_data_all.root WmnG_data*000*.root
hadd -f ZnnG_wenu_above0p5_all.root ZnnG_wenu*000*_above0p5.root
hadd -f ZnnG_wenu_below0p5_all.root ZnnG_wenu*000*_below0p5.root
hadd -f WenG_wenu_all.root WenG_wenu*000*.root
hadd -f ZnnG_bhalo_above0p5_all.root ZnnG_bhalo*000*_above0p5.root
hadd -f ZnnG_bhalo_below0p5_all.root ZnnG_bhalo*000*_below0p5.root
hadd -f ZnnG_qcd_above0p5_all.root ZnnG_qcd*000*_above0p5.root
hadd -f ZnnG_qcd_below0p5_all.root ZnnG_qcd*000*_below0p5.root
hadd -f ZeeG_qcd_all.root ZeeG_qcd*000*.root
hadd -f ZmmG_qcd_all.root ZmmG_qcd*000*.root
hadd -f WenG_qcd_all.root WenG_qcd*000*.root
hadd -f WmnG_qcd_all.root WmnG_qcd*000*.root
# Once these are all successfully merged, delete the separate files
rm ZnnG_data*000*_above0p5.root
rm ZnnG_data*000*_below0p5.root
rm ZeeG_data*000*.root
rm ZmmG_data*000*.root
rm WenG_data*000*.root
rm WmnG_data*000*.root
rm ZnnG_wenu*000*_above0p5.root
rm ZnnG_wenu*000*_below0p5.root
rm WenG_wenu*000*.root
rm ZnnG_bhalo*000*_above0p5.root
rm ZnnG_bhalo*000*_below0p5.root
rm ZnnG_qcd*000*_above0p5.root
rm ZnnG_qcd*000*_below0p5.root
rm ZeeG_qcd*000*.root
rm ZmmG_qcd*000*.root
rm WenG_qcd*000*.root
rm WmnG_qcd*000*.root


# Evaluate pdf and scale systematics
root -l -q -b pdfscale_plotter.C

# Get prefit histograms
# Make sure znng_*0p5_plotter.C has been updated based on print_plotDM_commands.py (see "Finding DM cross sections")
root -l -q -b znng_above0p5_plotter.C
root -l -q -b znng_below0p5_plotter.C
root -l -q -b zeeg_plotter.C
root -l -q -b zmmg_plotter.C
root -l -q -b weng_plotter.C
root -l -q -b wmng_plotter.C

# Compute transfer factors
root -l -q -b transfer_histo_collator.C


## Finding DM cross sections ##
# Shouldn't need to do this unless the desired DM samples have changed
# Make sure samples_central.txt, das_client.py, and ana.py are in working area
# Make sure samples_central.txt is up to date (from multicrab.py)
# Make sure analyze_DM_NLO.py and analyze_DM_LO.py have been run and the .out files are in working area
python print_files_lister_commands.py > files_lister_commands.sh
chmod 777 files_lister_commands.sh
./files_lister_commands.sh
python print_ana_commands.py > ana_commands.sh
chmod 777 ana_commands.sh
# Issue: All commands in ana_commands want to run at once.
# Maybe run a few at a time, so as not to fill up processes on the login machine.
./ana_commands.sh
# First run ana_commands.sh to get xsecs,
# and also run analyze_central_signals.py to get nevents
python print_plotDM_commands.py
# Copy and paste the output into znng_above0p5_plotter.C and znng_below0p5_plotter.C
