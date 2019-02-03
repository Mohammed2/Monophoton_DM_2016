./rootcom postAnalyzer_ZnnG_mc_znng_JESPES ZnnG_mc_znng_JESPES
./rootcom postAnalyzer_ZnnG_mc_zllg_genPtOver300_JESPES ZnnG_mc_zllg_genPtOver300_JESPES
./rootcom postAnalyzer_ZnnG_mc_zllg_genPtUnder300_JESPES ZnnG_mc_zllg_genPtUnder300_JESPES
./rootcom postAnalyzer_ZnnG_mc_wg_JESPES ZnnG_mc_wg_JESPES
./rootcom postAnalyzer_ZnnG_mc_JESPES ZnnG_mc_JESPES
./rootcom postAnalyzer_ZnnG_mc_WZ_JESPES ZnnG_mc_WZ_JESPES

./rootcom postAnalyzer_ZeeG_mc_zllg_genPtOver300_JESPES ZeeG_mc_zllg_genPtOver300_JESPES
./rootcom postAnalyzer_ZeeG_mc_zllg_genPtUnder300_JESPES ZeeG_mc_zllg_genPtUnder300_JESPES
./rootcom postAnalyzer_ZeeG_mc_WZ_JESPES ZeeG_mc_WZ_JESPES
./rootcom postAnalyzer_ZeeG_mc_JESPES ZeeG_mc_JESPES

./rootcom postAnalyzer_ZmmG_mc_zllg_genPtOver300_JESPES ZmmG_mc_zllg_genPtOver300_JESPES
./rootcom postAnalyzer_ZmmG_mc_zllg_genPtUnder300_JESPES ZmmG_mc_zllg_genPtUnder300_JESPES
./rootcom postAnalyzer_ZmmG_mc_WZ_JESPES ZmmG_mc_WZ_JESPES
./rootcom postAnalyzer_ZmmG_mc_JESPES ZmmG_mc_JESPES

./rootcom postAnalyzer_WenG_mc_zllg_genPtOver300_JESPES WenG_mc_zllg_genPtOver300_JESPES
./rootcom postAnalyzer_WenG_mc_zllg_genPtUnder300_JESPES WenG_mc_zllg_genPtUnder300_JESPES
./rootcom postAnalyzer_WenG_mc_wg_JESPES WenG_mc_wg_JESPES
./rootcom postAnalyzer_WenG_mc_WZ_JESPES WenG_mc_WZ_JESPES
./rootcom postAnalyzer_WenG_mc_JESPES WenG_mc_JESPES

./rootcom postAnalyzer_WmnG_mc_zllg_genPtOver300_JESPES WmnG_mc_zllg_genPtOver300_JESPES
./rootcom postAnalyzer_WmnG_mc_zllg_genPtUnder300_JESPES WmnG_mc_zllg_genPtUnder300_JESPES
./rootcom postAnalyzer_WmnG_mc_wg_JESPES WmnG_mc_wg_JESPES
./rootcom postAnalyzer_WmnG_mc_WZ_JESPES WmnG_mc_WZ_JESPES
./rootcom postAnalyzer_WmnG_mc_JESPES WmnG_mc_JESPES



nohup ./ZnnG_mc_wg_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET_mitExt/CRAB_UserFiles/crab_WGJets_LO_mitExt/171005_190151/0000/ ZnnG_JESPES_WGJets_mitExt.root -1 1000 >& ZnnG_JESPES_WGJets_mitExt.txt &
nohup ./WenG_mc_wg_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET_mitExt/CRAB_UserFiles/crab_WGJets_LO_mitExt/171005_190151/0000/ WenG_JESPES_WGJets_mitExt.root -1 1000 >& WenG_JESPES_WGJets_mitExt.txt &
nohup ./WmnG_mc_wg_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET_mitExt/CRAB_UserFiles/crab_WGJets_LO_mitExt/171005_190151/0000/ WmnG_JESPES_WGJets_mitExt.root -1 1000 >& WmnG_JESPES_WGJets_mitExt.txt &
nohup ./ZnnG_mc_znng_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET_others/ZNuNuGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZNuNuGJets_LO_ext/170828_151900/0000/ ZnnG_JESPES_ZNuNuGJets_ext.root -1 1000 >& ZnnG_JESPES_ZNuNuGJets_ext.txt &
nohup ./ZnnG_mc_wg_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET_others/WGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_WGJets_LO_ext/170828_152007/0000/ ZnnG_JESPES_WGJets_ext.root -1 1000 >& ZnnG_JESPES_WGJets_ext.txt &
nohup ./WenG_mc_wg_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET_others/WGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_WGJets_LO_ext/170828_152007/0000/ WenG_JESPES_WGJets_ext.root -1 1000 >& WenG_JESPES_WGJets_ext.txt &
nohup ./WmnG_mc_wg_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET_others/WGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_WGJets_LO_ext/170828_152007/0000/ WmnG_JESPES_WGJets_ext.root -1 1000 >& WmnG_JESPES_WGJets_ext.txt &

nohup ./ZnnG_mc_zllg_genPtOver300_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET/ZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets_LO/170715_001458/0000/ ZnnG_JESPES_ZLLGJets_130_over300.root -1 1000 >& ZnnG_JESPES_ZLLGJets_genPtOver300.txt &
nohup ./ZnnG_mc_zllg_genPtUnder300_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET/ZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets_LO/170715_001458/0000/ ZnnG_JESPES_ZLLGJets_130_under300.root -1 1000 >& ZnnG_JESPES_ZLLGJets_genPtUnder300.txt &
nohup ./ZnnG_mc_zllg_genPtOver300_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET/ZLLGJets_MonoPhoton_PtG-300_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets_LO_highPtext/170715_001525/0000/ ZnnG_JESPES_ZLLGJets_300.root -1 1000 >& ZnnG_JESPES_ZLLGJets_300.txt &
nohup ./ZeeG_mc_zllg_genPtOver300_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET/ZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets_LO/170715_001458/0000/ ZeeG_JESPES_ZLLGJets_130_over300.root -1 1000 >& ZeeG_JESPES_ZLLGJets_genPtOver300.txt &
nohup ./ZeeG_mc_zllg_genPtUnder300_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET/ZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets_LO/170715_001458/0000/ ZeeG_JESPES_ZLLGJets_130_under300.root -1 1000 >& ZeeG_JESPES_ZLLGJets_genPtUnder300.txt &
nohup ./ZeeG_mc_zllg_genPtOver300_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET/ZLLGJets_MonoPhoton_PtG-300_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets_LO_highPtext/170715_001525/0000/ ZeeG_JESPES_ZLLGJets_300.root -1 1000 >& ZeeG_JESPES_ZLLGJets_300.txt &
nohup ./ZmmG_mc_zllg_genPtOver300_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET/ZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets_LO/170715_001458/0000/ ZmmG_JESPES_ZLLGJets_130_over300.root -1 1000 >& ZmmG_JESPES_ZLLGJets_genPtOver300.txt &
nohup ./ZmmG_mc_zllg_genPtUnder300_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET/ZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets_LO/170715_001458/0000/ ZmmG_JESPES_ZLLGJets_130_under300.root -1 1000 >& ZmmG_JESPES_ZLLGJets_genPtUnder300.txt &
nohup ./ZmmG_mc_zllg_genPtOver300_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET/ZLLGJets_MonoPhoton_PtG-300_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets_LO_highPtext/170715_001525/0000/ ZmmG_JESPES_ZLLGJets_300.root -1 1000 >& ZmmG_JESPES_ZLLGJets_300.txt &
nohup ./WenG_mc_zllg_genPtOver300_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET/ZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets_LO/170715_001458/0000/ WenG_JESPES_ZLLGJets_130_over300.root -1 1000 >& WenG_JESPES_ZLLGJets_genPtOver300.txt &
nohup ./WenG_mc_zllg_genPtUnder300_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET/ZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets_LO/170715_001458/0000/ WenG_JESPES_ZLLGJets_130_under300.root -1 1000 >& WenG_JESPES_ZLLGJets_genPtUnder300.txt &
nohup ./WenG_mc_zllg_genPtOver300_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET/ZLLGJets_MonoPhoton_PtG-300_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets_LO_highPtext/170715_001525/0000/ WenG_JESPES_ZLLGJets_300.root -1 1000 >& WenG_JESPES_ZLLGJets_300.txt &
nohup ./WmnG_mc_zllg_genPtOver300_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET/ZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets_LO/170715_001458/0000/ WmnG_JESPES_ZLLGJets_130_over300.root -1 1000 >& WmnG_JESPES_ZLLGJets_genPtOver300.txt &
nohup ./WmnG_mc_zllg_genPtUnder300_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET/ZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets_LO/170715_001458/0000/ WmnG_JESPES_ZLLGJets_130_under300.root -1 1000 >& WmnG_JESPES_ZLLGJets_genPtUnder300.txt &
nohup ./WmnG_mc_zllg_genPtOver300_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET/ZLLGJets_MonoPhoton_PtG-300_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets_LO_highPtext/170715_001525/0000/ WmnG_JESPES_ZLLGJets_300.root -1 1000 >& WmnG_JESPES_ZLLGJets_300.txt &

nohup ./ZnnG_mc_znng_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET_others/ZNuNuGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZNuNuGJets_LO/170828_151828/0000/ ZnnG_JESPES_ZNuNuGJets.root -1 1000 >& ZnnG_JESPES_ZNuNuGJets.txt &
nohup ./ZnnG_mc_wg_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET_others/WGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_WGJets_LO/170828_151934/0000/ ZnnG_JESPES_WGJets.root -1 1000 >& ZnnG_JESPES_WGJets.txt &
./MakeCondorFiles.csh ZnnG_mc_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-40To100/170224_162416/0000/  ZnnG_JESPES_GJets_HT-40To100.root -1 1000 ZnnG_JESPES_GJets_HT-40To100
./MakeCondorFiles.csh ZnnG_mc_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-100To200/170224_162449/0000/  ZnnG_JESPES_GJets_HT-100To200.root -1 1000 ZnnG_JESPES_GJets_HT-100To200
./MakeCondorFiles.csh ZnnG_mc_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-200To400/170224_162521/0000/  ZnnG_JESPES_GJets_HT-200To400.root -1 1000 ZnnG_JESPES_GJets_HT-200To400
./MakeCondorFiles.csh ZnnG_mc_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-400To600/170224_162554/0000/  ZnnG_JESPES_GJets_HT-400To600.root -1 1000 ZnnG_JESPES_GJets_HT-400To600
./MakeCondorFiles.csh ZnnG_mc_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-600ToInf/170224_162630/0000/  ZnnG_JESPES_GJets_HT-600ToInf.root -1 1000 ZnnG_JESPES_GJets_HT-600ToInf
./MakeCondorFiles.csh ZnnG_mc_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/WToMuNu_M-100_TuneCUETP8M1_13TeV-pythia8/crab_WToMuNu/170224_162311/0000/  ZnnG_JESPES_WToMuNu.root -1 1000 ZnnG_JESPES_WToMuNu
./MakeCondorFiles.csh ZnnG_mc_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/WToTauNu_M-100_TuneCUETP8M1_13TeV-pythia8-tauola/crab_WToTauNu/170224_162343/0000/  ZnnG_JESPES_WToTauNu.root -1 1000 ZnnG_JESPES_WToTauNu
./MakeCondorFiles.csh ZnnG_mc_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTGJets/170224_162704/0000/  ZnnG_JESPES_TTGJets.root -1 1000 ZnnG_JESPES_TTGJets
./MakeCondorFiles.csh ZnnG_mc_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/DiPhotonJets_MGG-80toInf_13TeV_amcatnloFXFX_pythia8/crab_DiPhotonJets/170224_162923/0000/  ZnnG_JESPES_Diphoton.root -1 1000 ZnnG_JESPES_Diphoton
./MakeCondorFiles.csh ZnnG_mc_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8/crab_TG/170224_161901/0000/  ZnnG_JESPES_TGJets.root -1 1000 ZnnG_JESPES_TGJets
./MakeCondorFiles.csh ZnnG_mc_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/ZZ_TuneCUETP8M1_13TeV-pythia8/crab_ZZ/170224_162051/0000/  ZnnG_JESPES_ZZ.root -1 1000 ZnnG_JESPES_ZZ
./MakeCondorFiles.csh ZnnG_mc_WZ_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/WZ_TuneCUETP8M1_13TeV-pythia8/crab_WZ/170224_162012/0000/  ZnnG_JESPES_WZ.root -1 1000 ZnnG_JESPES_WZ
./MakeCondorFiles.csh ZnnG_mc_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/WW_TuneCUETP8M1_13TeV-pythia8/crab_WW/170224_162956/0000/  ZnnG_JESPES_WW.root -1 1000 ZnnG_JESPES_WW

./MakeCondorFiles.csh ZeeG_mc_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTGJets/170224_162704/0000/  ZeeG_JESPES_TTGJets.root -1 1000 ZeeG_JESPES_TTGJets
./MakeCondorFiles.csh ZeeG_mc_WZ_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/WZ_TuneCUETP8M1_13TeV-pythia8/crab_WZ/170224_162012/0000/  ZeeG_JESPES_WZ.root -1 1000 ZeeG_JESPES_WZ

./MakeCondorFiles.csh ZmmG_mc_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTGJets/170224_162704/0000/  ZmmG_JESPES_TTGJets.root -1 1000 ZmmG_JESPES_TTGJets
./MakeCondorFiles.csh ZmmG_mc_WZ_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/WZ_TuneCUETP8M1_13TeV-pythia8/crab_WZ/170224_162012/0000/  ZmmG_JESPES_WZ.root -1 1000 ZmmG_JESPES_WZ

nohup ./WenG_mc_wg_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET_others/WGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_WGJets_LO/170828_151934/0000/ WenG_JESPES_WGJets.root -1 1000 >& WenG_JESPES_WGJets.txt &
./MakeCondorFiles.csh WenG_mc_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTGJets/170224_162704/0000/  WenG_JESPES_TTGJets.root -1 1000 WenG_JESPES_TTGJets
./MakeCondorFiles.csh WenG_mc_WZ_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/WZ_TuneCUETP8M1_13TeV-pythia8/crab_WZ/170224_162012/0000/  WenG_JESPES_WZ.root -1 1000 WenG_JESPES_WZ
./MakeCondorFiles.csh WenG_mc_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8/crab_TG/170224_161901/0000/  WenG_JESPES_TGJets.root -1 1000 WenG_JESPES_TGJets
./MakeCondorFiles.csh WenG_mc_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/DiPhotonJets_MGG-80toInf_13TeV_amcatnloFXFX_pythia8/crab_DiPhotonJets/170224_162923/0000/  WenG_JESPES_Diphoton.root -1 1000 WenG_JESPES_Diphoton
./MakeCondorFiles.csh WenG_mc_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/WW_TuneCUETP8M1_13TeV-pythia8/crab_WW/170224_162956/0000/  WenG_JESPES_WW.root -1 1000 WenG_JESPES_WW

nohup ./WmnG_mc_wg_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17_withGenPhoET_others/WGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_WGJets_LO/170828_151934/0000/ WmnG_JESPES_WGJets.root -1 1000 >& WmnG_JESPES_WGJets.txt &
./MakeCondorFiles.csh WmnG_mc_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTGJets/170224_162704/0000/  WmnG_JESPES_TTGJets.root -1 1000 WmnG_JESPES_TTGJets
./MakeCondorFiles.csh WmnG_mc_WZ_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/WZ_TuneCUETP8M1_13TeV-pythia8/crab_WZ/170224_162012/0000/  WmnG_JESPES_WZ.root -1 1000 WmnG_JESPES_WZ
./MakeCondorFiles.csh WmnG_mc_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8/crab_TG/170224_161901/0000/  WmnG_JESPES_TGJets.root -1 1000 WmnG_JESPES_TGJets
./MakeCondorFiles.csh WmnG_mc_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/DiPhotonJets_MGG-80toInf_13TeV_amcatnloFXFX_pythia8/crab_DiPhotonJets/170224_162923/0000/  WmnG_JESPES_Diphoton.root -1 1000 WmnG_JESPES_Diphoton
./MakeCondorFiles.csh WmnG_mc_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/WW_TuneCUETP8M1_13TeV-pythia8/crab_WW/170224_162956/0000/  WmnG_JESPES_WW.root -1 1000 WmnG_JESPES_WW

# ZeeG and ZmmG: WW = ZZ = TGJets = Diphoton = 0
# WenG and WmnG: ZZ = 0



./MakeCondorFiles.csh ./ZnnG_mc_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/DarkMatter_MonoPhoton_V_Mx-1_Mv-500_gDMgQ1_LO_13TeV-madgraph/crab_DMVMx-1Mv-500/170224_171031/0000/  ZnnG_JESPES_DM_V_Mx-1_Mv-500.root -1 1000 ZnnG_JESPES_DM_V_Mx-1_Mv-500
./MakeCondorFiles.csh ./ZnnG_mc_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/DarkMatter_MonoPhoton_V_Mx-1_Mv-1000_gDMgQ1_LO_13TeV-madgraph/crab_DMVMx-1Mv-1000/170224_170604/0000/  ZnnG_JESPES_DM_V_Mx-1_Mv-1000.root -1 1000 ZnnG_JESPES_DM_V_Mx-1_Mv-1000
./MakeCondorFiles.csh ./ZnnG_mc_JESPES /hdfs/store/user/jjbuch/ntuples_moriond17/RunIISummer16MiniAODv2-PUMoriond17/DarkMatter_MonoPhoton_V_Mx-1_Mv-2000_gDMgQ1_LO_13TeV-madgraph/crab_DMVMx-1Mv-2000/170224_170751/0000/  ZnnG_JESPES_DM_V_Mx-1_Mv-2000.root -1 1000 ZnnG_JESPES_DM_V_Mx-1_Mv-2000
