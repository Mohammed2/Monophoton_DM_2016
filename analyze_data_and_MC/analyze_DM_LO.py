import os
import subprocess

nlodm_folders = os.listdir('/hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/')

for folderName in nlodm_folders:
  if folderName[:len('DarkMatter_')] == 'DarkMatter_':
    folder_path = '/hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/'+folderName+'/'
    folderName_parts = folderName.split('_')
    sample_name = 'DM_LO_'+folderName_parts[2]+'_'+folderName_parts[3]+'_'+folderName_parts[4]
    
    expanding_folder = True
    
    while(expanding_folder):
      fileNames = os.listdir(folder_path)
      
      foundRootFile = False
      for fileName in fileNames:
        if len(fileName) > 5 and fileName[-len('.root'):] == '.root':
          foundRootFile = True
      
      # If this directory contains ROOT files, analyze them,
      # then move to the next folder in nlodm_folders.
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
        # print ' '.join(command_args)
        subprocess.call(command_args)
        expanding_folder = False
      
      # If this directory doesn't contain a ROOT file, open what it does contain and look inside.
      elif len(fileNames) == 1:
        folder_path += fileNames[0]+'/'
      
      else:
        print 'Error: Unable to find a ROOT file or open a folder'
        print 'folder_path = '+folder_path
        print 'fileNames = '+str(fileNames)
        expanding_folder = False
