# #!/usr/bin/env bash
python pMHC_analyse_v2.py --do_bsa --ray --do_apbs --mhc_class II --peptide SC2_6 		 	--pdb DR1_SC2_6_refine_36.pdb 				--mtz DR1_SC2_6_refine_36.mtz					--chains ABC
python pMHC_analyse_v2.py --do_bsa --ray --do_apbs --mhc_class II --peptide SC2_9 			--pdb DR1_SC2_9_refine_50.pdb 				--mtz DR1_SC2_9_refine_50.mtz					--chains ABC
python pMHC_analyse_v2.py --do_bsa --ray --do_apbs --mhc_class II --peptide SC2_22 		 	--pdb DR1_SC2_22_refine_41.pdb 				--mtz DR1_SC2_22_refine_41.mtz					--chains ABC
python pMHC_analyse_v2.py --do_bsa --ray --do_apbs --mhc_class II --peptide SC2_24 		 	--pdb DR1_SC2_24_refine_72.pdb 				--mtz DR1_SC2_24_refine_72.mtz					--chains ABC
python pMHC_analyse_v2.py --do_bsa --ray --do_apbs --mhc_class II --peptide X2_A2 			--pdb DR1_SC2_X2_A2_refine_47.pdb 			--mtz DR1_SC2_X2_A2_refine_47.mtz				--chains ABC
python pMHC_analyse_v2.py --do_bsa --ray --do_apbs --mhc_class II --peptide X2_A7 		 	--pdb DR1_SC2_X2_A7_refine_44.pdb 			--mtz DR1_SC2_X2_A7_refine_44.mtz				--chains ABC
python pMHC_analyse_v2.py --do_bsa --ray --do_apbs --mhc_class II --peptide X2_A2_omicron	--pdb DR1_SC2_X2_A2_omicron_refine_23.pdb 	--mtz DR1_SC2_X2_A2_omicron_refine_23.mtz		--chains ABC
python pMHC_analyse_v2.py --do_bsa --ray --do_apbs --mhc_class II --peptide SC2_6_omicron 	--pdb DR1_SC2_6_omicron_refine_46.pdb		--mtz DR1_SC2_6_omicron_refine_46.mtz			--chains ABC
python pMHC_analyse_v2.py --do_bsa --ray --do_apbs --mhc_class II --peptide SC2_6_omicron 	--pdb DR1_SC2_6_omicron_refine_46.pdb		--mtz DR1_SC2_6_omicron_refine_46.mtz			--chains DEF
python pMHC_analyse_v2.py --do_bsa --ray --do_apbs --mhc_class II --peptide SC2_6_omicron 	--pdb DR1_SC2_6_omicron_refine_46.pdb		--mtz DR1_SC2_6_omicron_refine_46.mtz			--chains GHI

python P1_Lys.py --ray --mhc_class II --peptide SC2_6 		 	--pdb DR1_SC2_6_refine_36.pdb 				--mtz DR1_SC2_6_refine_36.mtz					--chains ABC
python P1_Lys.py --ray --mhc_class II --peptide SC2_9 			--pdb DR1_SC2_9_refine_50.pdb 				--mtz DR1_SC2_9_refine_50.mtz					--chains ABC
python P1_Lys.py --ray --mhc_class II --peptide SC2_22 		 	--pdb DR1_SC2_22_refine_41.pdb 				--mtz DR1_SC2_22_refine_41.mtz					--chains ABC
python P1_Lys.py --ray --mhc_class II --peptide SC2_24 		 	--pdb DR1_SC2_24_refine_72.pdb 				--mtz DR1_SC2_24_refine_72.mtz					--chains ABC
python P1_Lys.py --ray --mhc_class II --peptide X2_A2 			--pdb DR1_SC2_X2_A2_refine_47.pdb 			--mtz DR1_SC2_X2_A2_refine_47.mtz				--chains ABC
python P1_Lys.py --ray --mhc_class II --peptide X2_A7 		 	--pdb DR1_SC2_X2_A7_refine_44.pdb 			--mtz DR1_SC2_X2_A7_refine_44.mtz				--chains ABC

python P1_Lys.py --ray --mhc_class II --peptide X2_A2_omicron	--pdb DR1_SC2_X2_A2_omicron_refine_23.pdb 	--mtz DR1_SC2_X2_A2_omicron_refine_23.mtz		--chains ABC

python P1_Lys.py --ray --mhc_class II --peptide SC2_6_omicron 	--pdb DR1_SC2_6_omicron_refine_46.pdb		--mtz DR1_SC2_6_omicron_refine_46.mtz			--chains ABC
python P1_Lys.py --ray --mhc_class II --peptide SC2_6_omicron 	--pdb DR1_SC2_6_omicron_refine_46.pdb		--mtz DR1_SC2_6_omicron_refine_46.mtz			--chains DEF
python P1_Lys.py --ray --mhc_class II --peptide SC2_6_omicron 	--pdb DR1_SC2_6_omicron_refine_46.pdb		--mtz DR1_SC2_6_omicron_refine_46.mtz			--chains GHI


