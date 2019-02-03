import os
import subprocess

ntuple_folders = []
ntuple_folders.extend(os.listdir('/hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET_DM_EWK/'))
ntuple_folders.extend(os.listdir('/hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_NLODM_pta130_central/'))
ntuple_folders.extend(os.listdir('/hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET_ADD/'))

# DEBUG
# print 'ntuple_folders:'
# print ntuple_folders

for folderName in ntuple_folders:
  # DEBUG
  # print 'folderName: %s' % folderName
  
  if folderName[:len('DarkMatter_')] == 'DarkMatter_' or folderName[:len('ADDmonoPhoton')] == 'ADDmonoPhoton':
    # DEBUG
    # print 'good folderName'
    
    folderName_parts = folderName.split('_')
    folder_path = 'FOLDER_PATH NOT SET'
    sample_name = 'SAMPLE_NAME NOT SET'
    if folderName_parts[2] == 'EWK':
      folder_path = '/hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET_DM_EWK/'+folderName+'/'
      sample_name = 'DM_EWK_'+folderName_parts[4]
    elif folderName_parts[2] == 'NLO':
      folder_path = '/hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_NLODM_pta130_central/'+folderName+'/'
      sample_name = 'DM_NLO_'+folderName_parts[3]+'_'+folderName_parts[4]+'_'+folderName_parts[5]+'_pta130'
    elif folderName_parts[0] == 'ADDmonoPhoton':
      folder_path = '/hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET_ADD/'+folderName+'/'
      sample_name = 'ADDmonoPhoton_'+folderName_parts[1]+folderName_parts[2]
    
    # DEBUG
    print 'sample_name: %s' % sample_name
    
    still_expanding_folder = True
    
    while(still_expanding_folder):
      fileNames = os.listdir(folder_path)
      
      foundRootFile = False
      for fileName in fileNames:
        if len(fileName) > 5 and fileName[-len('.root'):] == '.root':
          foundRootFile = True
      
      # If this directory contains ROOT files, analyze them,
      # then move to the next folder in ntuple_folders.
      if foundRootFile:
        command_args = []
        command_args.append('./MakeCondorFiles.csh')
        command_args.append('ZnnG_mc_JESPES')
        command_args.append(folder_path+' ')
        command_args.append('ZnnG_JESPES_'+sample_name+'.root')
        command_args.append('-1')
        command_args.append('1000')
        command_args.append('ZnnG_JESPES_'+sample_name)
        # DEBUG
        print ' '.join(command_args)
        subprocess.call(command_args)
        still_expanding_folder = False
      
      # If this directory doesn't contain a ROOT file, open what it does contain and look inside.
      elif len(fileNames) == 1:
        folder_path += fileNames[0]+'/'
      
      else:
        print 'Error: Unable to find a ROOT file or open a folder'
        print 'folder_path = '+folder_path
        print 'fileNames = '+str(fileNames)
        still_expanding_folder = False
