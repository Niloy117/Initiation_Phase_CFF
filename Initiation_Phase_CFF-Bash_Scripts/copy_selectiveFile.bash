#!/bin/bash

targetDirectry=/home/rniloy/PROJECT_CELL/2D_codes/unit_blck/restructuring_Data_Structure/pulley_point_simulations/WT_pulling_WO_TIIF/CM4PCS/BkUpSimData/

#pulley_point_simulations/WO_pulling_WO_TIIF_Y23_managed/CM1PCS/BkUpSimData/
#WO_pulling_WO_TIIF_Y23_managed/CM2PCS/BkUpSimData/
#WO_pulling_WT_EqualTIIF_Y23_and_CurvtrManaged/CM2PCS/BkUpSimData/

#PULLING_TEST_ON_WO_WT_TIIF/TO_CHK_WHAT_PULL_DOES_ONLY_WITH_BNDRY_ADJSTMNT/BkUpSimData/
#PULLING_TEST_ON_WO_WT_TIIF/BkUpSimData/
#WO_pulling_WT_TIIF_FRM_ZERO_TO_SPCFC_VAL/CM4PCS/BkUpSimData/

#perturbtnOnPCS_withEquvlntPP_TIFF_atLatest_JointNode/CM1PCS/
#perturbtnOnPCS_withEquvlntPP_andNo_TIIF/CM1PCS/

#making_PPbfraft_matchWithExp/prtrbtnTst_CgX0_CgY0/
#remvingVrtclForceAndRetainNodePosition/rescaledVersion/
#remvingVrtclForceAndRetainNodePosition/CM2PCS/
#making_PPbfraft_matchWithExp/CM4PCS/

#Force_LtrlSprNC1_comboTst_EqualLtrlTnsn/CgY=0.75/
#ForceApclLayr_LtSpShrtn_BndryPull/CgY=1.00_l0var_All_pullF_All/
#Force_LtrlSprNC1_comboTst_CgY=1.00/ForceDownLtrlShortnALWAYS_TiltReduced/

#VF_withJointNodeForceAndLtrlShrtn

#VF_withAddedNCPair10
#VF_withAddedNCPair4
#VF_withAddedNCPair3
#VF_EffctOfVF_tst
#VF_CF_WT_Nrt_P_2G5C_m

#ReleasingCntrlN_uptoN2PairCellsApclNode_WT_VFarea_WT_highKS_LtrlIC_ShortnTst
#ReleasingCntrlN_uptoN2PairCellsApclNode_WT_VFarea_WT_highKS
         
#mkdir $targetDirectry

#WO_VF_AreaTst/ForceDown
#ReleasingCntrlN_uptoN2PairCellsApclNode_WT_VFarea
#ReleasingCntrlNN1PairCellsApclNode_WT_VFarea

#CM2ICS/IntrlPltnFrmDirectMeasuredPCSDiagTnsnZero
#CM4ICS/IntrlPltnFrmDirectMeasuredPCSDiagTnsnZero
#CM3ICS/IntrlPltnFrmDirectMeasuredPCSDiagTnsnZero
#CM3ICS/IntrlPltnFrmDirectMeasuredPCS
#CM1PCS/DirectCurvtrMeasrmnt
#CM2PCS/DirectCurvtrMeasrmnt
#CM3PCS/DirectCurvtrMeasrmnt
#CM4PCS/DirectCurvtrMeasrmnt

#CM3PCS/DirectCurvtrMeasrmntTstAftCell
#CM2PCS/DirectCurvtrMeasrmntTst_AddedFlder
#CM2PCS/PressInvPropCurvtrTest
#CM2PCS/l0_eql_l0ExpctdLen_tst

#CM1PCS/k_chge_tst_only
#CM2PCS/k_chge_tst_only
#CM3PCS/k_chge_tst_only

#CM1PCS/after_applying_high_ks_algorithm
#CM2PCS/after_applying_high_ks_algorithm
#CM3PCS/after_applying_high_ks_algorithm
#CM4PCS/after_applying_high_ks_algorithm
#CM1PCS/yNode23ForcePrtrbtn_with_high_ks_algorithm
#CM2PCS/yNode23ForcePrtrbtn_with_high_ks_algorithm
#CM3PCS/yNode23ForcePrtrbtn_with_high_ks_algorithm
#CM4PCS/yNode23ForcePrtrbtn_with_high_ks_algorithm

#NI_DataInitPhase/HighKTest/WT_AR_Fix_BndryCellApBsEquivlncePressAdjst
#NI_DataInitPhase/HighKTest/WT_AR_Fix_BndryCellApBsEquivlnce
#NI_DataInitPhase/HighKTest/WT_AR_Fix
#NI_DataInitPhase/HighKTest/WO_AR_Fix
#CM4PCS
#DataAftRestrctingCellsMeet2/CM1_nodeDist_0pnt06
#ForceWithBsalSideOnly/Remving_Apcl_Side_Force_and_Straightning_CS
#ForceWithApclSideOnly
#SECOND_MEET_withForce
#InitPhaseTrialAftTwoCellMeet/Tr1
#MAKING_prfct_control_state
#DIAGONAL_SPR_TENSION_NOTRELEASE_TEST
#BfrRebuttalCHANGES
#VrticalForceData/CntrlCellApSprShortn00prcnt
#SevenAddedCell
#IfExtrapltnHelps/Cg_30
#propsflnm/Using_step_frm_NonUniStrct
#NI_AddedCellDataUpdtedFeb22/sevenCycl/WithShakingIF
#Extrapolate_NI/woTnsnPresChk

#NI_AddedCelldata/sevenCycl/NonShakingIF_WT3
#WithShakingIF #NonShakingIF_WT3

#NI_data/AftShftdMatchCS1CS4
#NI_data/AvrgingPulleyPropWOOuterLayer/

flNm='NI_modelInitiatn'
#flNm='NI_modelInitiatn_WT_VF'
#flNm='BeginStgAddCell_'

#COUNTMAX=14

COUNT_strt=1
COUNT_fnsh=142
COUNTER=$COUNT_strt

echo $COUNT_strt,$COUNT_fnsh,$COUNTER
echo $targetDirectry

#while [  $COUNTER -le $COUNTMAX ] ; do
while [ "$COUNTER" -ge $COUNT_strt ] && [ "$COUNTER" -le $COUNT_fnsh ]; do
    
    echo The counter is $COUNTER
    
    currID=$COUNTER
    echo $currID
    
    FlnmAW=$(printf "%sAW%04d.dat" "$flNm" "$currID")
    FlnmSW=$(printf "%sSW%04d.dat" "$flNm" "$currID")
    FlnmNW=$(printf "%sNW%04d.dat" "$flNm" "$currID")
    FlnmAP=$(printf "%sAP%04d.dat" "$flNm" "$currID")
    
    echo $FlnmAW $FlnmSW $FlnmNW $FlnmAP
    
    cp $FlnmAW $targetDirectry
    cp $FlnmSW $targetDirectry
    cp $FlnmNW $targetDirectry
    cp $FlnmAP $targetDirectry
    
    let COUNTER=COUNTER+1
    
done

