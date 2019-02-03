import os
import re

ntuple_folders = []
ntuple_folders.extend(os.listdir('/hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET_DM_EWK/'))
ntuple_folders.extend(os.listdir('/hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_NLODM_pta130_central/'))
ntuple_folders.extend(os.listdir('/hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET_ADD/'))

output_path = '/nfs_scratch/jjbuchanan/CMSSW_8_0_26_patch1/src/ControlRegionStudy_clean/'
fileNames = os.listdir(output_path)

sample_names = []
xsecs = []
nevents_totals = []

for folderName in ntuple_folders:
  if folderName[:len('DarkMatter_')] == 'DarkMatter_' or folderName[:len('ADDmonoPhoton')] == 'ADDmonoPhoton':
    sample_name = 'SAMPLE_NAME NOT SET'
    outfile_name = 'OUTFILE_NAME NOT SET'
    txtfile_name = 'TXTFILE_NAME NOT SET'
    
    folderName_parts = folderName.split('_')
    if folderName_parts[2] == 'EWK':
      sample_name = 'DM_EWK_'+folderName_parts[4]
      outfile_name = 'ZnnG_JESPES_'+sample_name
      txtfile_name = 'xsecinfo_'+sample_name
    elif folderName_parts[2] == 'NLO':
      sample_name = 'DM_NLO_'+folderName_parts[3]+'_'+folderName_parts[4]+'_'+folderName_parts[5]+'_pta130'
      outfile_name = 'ZnnG_JESPES_'+sample_name
      txtfile_name = 'xsecinfo_DM'+folderName_parts[3]+folderName_parts[4]+folderName_parts[5]+'_NLO_pta130'
    elif folderName_parts[0] == 'ADDmonoPhoton':
      sample_name = 'ADDmonoPhoton_'+folderName_parts[1]+'_'+folderName_parts[2]
      # outfile_name = 'ZnnG_JESPES_'+sample_name
      outfile_name = 'ZnnG_JESPES_ADDmonoPhoton_'+folderName_parts[1]+folderName_parts[2]
      txtfile_name = 'xsecinfo_ADD_'+folderName_parts[1]+'_'+folderName_parts[2]
    
    # Set nevents_total
    nevents_total = 'NEVENTS_NOT_SET'
    regex_outfile = re.compile(".*("+outfile_name+"\.out)")
    results = [m.group(0) for fileName in fileNames for m in [regex_outfile.search(fileName)] if m]
    if len(results) == 1:
      with open(results[0], 'r') as outfile_text:
        prefix = 'Number of events inspected (minus negative gen. weights): '
        for line in outfile_text:
          line = line.strip()
          if len(line) > len(prefix) and line[:len(prefix)] == prefix:
            nevents_total = (line.split())[-1]
    else:
      print 'No .out file found for sample_name='+sample_name
      print 'results='+str(results)
      continue
    
    # Set xsec
    xsec = 'XSEC_NOT_SET'
    regex_txtfile = re.compile("("+txtfile_name+"\.txt)")
    results = [m.group(0) for fileName in fileNames for m in [regex_txtfile.search(fileName)] if m]
    if len(results) == 1:
      with open(results[0], 'r') as outfile_text:
        prefix = 'After filter: final cross section = '
        for line in outfile_text:
          line = line.strip()
          if len(line) > len(prefix) and line[:len(prefix)] == prefix:
            xsec = (line.split())[6]
    else:
      print 'No .txt file found for txtfile_name='+txtfile_name
      print 'results='+str(results)
      continue
    
    sample_names.append(sample_name)
    xsecs.append(xsec)
    nevents_totals.append(nevents_total)

for i in range(0, len(sample_names)):
  print '  plot_DM("Photon_Et_range", "'+sample_names[i]+'", '+xsecs[i]+', '+nevents_totals[i]+');'
print ''
for i in range(0, len(sample_names)):
  print '  plot_DM("pfMET", "'+sample_names[i]+'", '+xsecs[i]+', '+nevents_totals[i]+');'
print ''
for i in range(0, len(sample_names)):
  print '  plot_DM("Mt", "'+sample_names[i]+'", '+xsecs[i]+', '+nevents_totals[i]+');'
