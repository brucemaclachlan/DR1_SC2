# #!/usr/bin/env bash
python pMHC_analyse_v2.py --ray --omit_mode phenix --mhc_class II --peptide SC2_6 		 	--pdb 8CMD.pdb	--mtz 8CMD.mtz	--chains ABC
python pMHC_analyse_v2.py --ray --omit_mode phenix --mhc_class II --peptide SC2_9 			--pdb 8CME.pdb	--mtz 8CME.mtz	--chains ABC
python pMHC_analyse_v2.py --ray --omit_mode phenix --mhc_class II --peptide SC2_22 			--pdb 8CMF.pdb	--mtz 8CMF.mtz	--chains ABC
python pMHC_analyse_v2.py --ray --omit_mode phenix --mhc_class II --peptide SC2_24 			--pdb 8CMG.pdb	--mtz 8CMG.mtz	--chains ABC
python pMHC_analyse_v2.py --ray --omit_mode phenix --mhc_class II --peptide X2_A2 			--pdb 8CMB.pdb	--mtz 8CMB.mtz	--chains ABC
python pMHC_analyse_v2.py --ray --omit_mode phenix --mhc_class II --peptide X2_A7 		 	--pdb 8CMC.pdb	--mtz 8CMC.mtz	--chains ABC
python pMHC_analyse_v2.py --ray --omit_mode phenix --mhc_class II --peptide X2_A2_omicron		--pdb 8CMH.pdb	--mtz 8CMH.mtz	--chains ABC
python pMHC_analyse_v2.py --ray --omit_mode phenix --mhc_class II --peptide SC2_6_omicron 	--pdb 8CMI.pdb	--mtz 8CMI.mtz	--chains ABC
python pMHC_analyse_v2.py --ray --omit_mode phenix --mhc_class II --peptide SC2_6_omicron 	--pdb 8CMI.pdb	--mtz 8CMI.mtz	--chains DEF
python pMHC_analyse_v2.py --ray --omit_mode phenix --mhc_class II --peptide SC2_6_omicron 	--pdb 8CMI.pdb	--mtz 8CMI.mtz	--chains GHI

python pocket_analysis_v2.py --ray --mhc_class II --peptide SC2_6 		 	--pdb 8CMD.pdb	--chains ABC
python pocket_analysis_v2.py --ray --mhc_class II --peptide SC2_9 			--pdb 8CME.pdb	--chains ABC
python pocket_analysis_v2.py --ray --mhc_class II --peptide SC2_22 			--pdb 8CMF.pdb	--chains ABC
python pocket_analysis_v2.py --ray --mhc_class II --peptide SC2_24 			--pdb 8CMG.pdb	--chains ABC
python pocket_analysis_v2.py --ray --mhc_class II --peptide X2_A2 			--pdb 8CMB.pdb	--chains ABC
python pocket_analysis_v2.py --ray --mhc_class II --peptide X2_A7 		 	--pdb 8CMC.pdb	--chains ABC
python pocket_analysis_v2.py --ray --mhc_class II --peptide X2_A2_omicron		--pdb 8CMH.pdb	--chains ABC
python pocket_analysis_v2.py --ray --mhc_class II --peptide SC2_6_omicron 	--pdb 8CMI.pdb	--chains ABC
python pocket_analysis_v2.py --ray --mhc_class II --peptide SC2_6_omicron 	--pdb 8CMI.pdb	--chains DEF
python pocket_analysis_v2.py --ray --mhc_class II --peptide SC2_6_omicron 	--pdb 8CMI.pdb	--chains GHI