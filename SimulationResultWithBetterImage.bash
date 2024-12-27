#! /bin/bash

rm Exp0_????.pdf
rm flNm_*.pdf

gnufile1="fill_color_test_powerlaw_stage1E.gnu"
gnufile2="fill_color_test_powerlaw_stage1E_with_addedCell.gnu"
gnufile3="fill_color_test_powerlaw_NIadded.gnu"

gnufileToRUN=$gnufile1 

echo $gnufileToRUN

#gnuplot $gnufileToRUN

#convert -delay 100 Exp0_???.gif OneNodeNI.gif


countMax=84 #21 #111 #187 
#411 #96#367#144#84#42#48#145#85#232#38#46#14#399#44#388#47#224#8#38#27#6#8#16

countStrt=1

count=$countStrt
loopCount=$countStrt

echo variables 
echo parameters values are $cyclNo $NfrmPerCycl $countMax $count $loopCount

while [ $count -le $countMax ] ; do
    
    #echo The counter is $count
    
    prevFlnm=$(printf "Exp0_%04d.eps"  "$count")
    currFlnm=$(printf "Exp0_%04d.pdf"  "$count")
    cropFlnm=$(printf "Exp0_c%04d.pdf" "$count")
    
    echo $prevFlnm $currFlnm
    ps2pdf $prevFlnm $currFlnm
    pdfcrop --margins '40 40 40 40' --clip  $currFlnm $cropFlnm
    
    #ps2pdf -dEPSCrop $prevFlnm
    #ps2pdf -dPDFSETTINGS=/prepress -dEPSCrop $prevFlnm
    
    let count=count+1
    
    flNm=$(printf "flNm_%04d.pdf" "$loopCount")
    echo $flNm
    
    cp $cropFlnm $flNm
    echo $cropFlnm $flNm
    
    let loopCount=loopCount+1
done

pdftk flNm_*.pdf cat output OneNodeNI00.pdf
pdftk OneNodeNI00.pdf cat output OneNodeNI.pdf

BlankPdf=$(printf "Blank.pdf")

finalPdf01=$(printf "stage1E.pdf")
finalPdf02=$(printf "ImageCreatn.pdf")
finalPdf03=$(printf "InitatorCellIncr3.pdf")
finalPdf04=$(printf "ForceIm4.pdf")
finalPdf05=$(printf "StiffArea.pdf")

finalPdf06=$(printf "ApplyingForceVrticallyDown3.pdf") #3 is only aft node added
finalPdf07=$(printf "ApplyingForceVrticallyDown2.pdf") #2 is for only NI,no begin in btwn
finalPdf08=$(printf "ApplyingForceVrticallyDown7NoNode.pdf")

finalPdf09=$(printf "SimlnForRebuttal.pdf")
finalPdf10=$(printf "First_HalfCycle.pdf")
finalPdf11=$(printf "CurvatureReductionAtCellTwoBfrInvgntngCell1.pdf")
finalPdf12=$(printf "lengthMatching.pdf")
finalPdf13=$(printf "Dignl_Tension_Test_Tr2.pdf")
finalPdf14=$(printf "ApclIncr_and_IncrForceToHold1.pdf")      #Frm16-27 [pdf 1]
finalPdf15=$(printf "ReduceTilt_Tr2_Aftpdf2.pdf")             #Frm28-32

finalPdf16=$(printf "Reduce_Force_On_Basal_Boundary7.pdf")
finalPdf17=$(printf "Reduce_Force_On_Apical_Boundary5.pdf")
finalPdf18=$(printf "Force_On_Both_Boundary.pdf")
finalPdf19=$(printf "RestrctAftCM1_withForceOnBothBndry.pdf")

finalPdf20=$(printf "RestrctBfrCM1_withForceDownwards.pdf")
finalPdf21=$(printf "RestrctBfrCM1_withForceDownwards_Tr1.pdf")
finalPdf22=$(printf "RestrctBfrCM1_withForceDownwards_Tr3.pdf")
finalPdf23=$(printf "RestrctBfrCM1_withForceDownwards_CM2ICS_Devlpmnt.pdf")
finalPdf24=$(printf "RestrctBfrCM1_withForceDownwards_CM2PCS_Devlpmnt.pdf")
finalPdf25=$(printf "RestrctBfrCM1_withForceDownwards_CM2PCS_Devlpmnt1.pdf")
finalPdf26=$(printf "RestrctBfrCM1_withForceDownwards_CM3ICS_Devlpmnt.pdf")
finalPdf27=$(printf "RestrctBfrCM1_withForceDownwards_CM3PCS_Devlpmnt.pdf")
finalPdf28=$(printf "RestrctBfrCM1_withForceDownwards_CM4ICS_Devlpmnt.pdf")
finalPdf29=$(printf "RestrctBfrCM1_withForceDownwards_CM4PCS_Devlpmnt2.pdf")

finalPdf30=$(printf "High_ks_test.pdf")
finalPdf31=$(printf "High_ks_test_withEL_ks08.pdf")
finalPdf32=$(printf "High_ks_test_withEL_ks08_AddedCell.pdf")
finalPdf33=$(printf "High_ks_test_withEL_ks08_BndryForceTst.pdf")
finalPdf34=$(printf "High_ks_test_withEL_ks10_A0adjstd.pdf")
finalPdf35=$(printf "High_ks_test_CM2PCS.pdf")

finalPdf35=$(printf "High_ks_test_CM3PCS.pdf")
finalPdf36=$(printf "High_ks_test_CM3PCS_y23Force_Cg40toCgMin15test.pdf")
finalPdf37=$(printf "High_ks_test_CM1PCS.pdf")
finalPdf38=$(printf "High_ks_test_CM1PCS_y23Force_Cg40toCgMin15test_.pdf")
finalPdf39=$(printf "High_ks_test_CM4PCS.pdf")
finalPdf40=$(printf "High_ks_test_CM4PCS_y23Force_Cg40toCgMin15test_tr3.pdf")

finalPdf41=$(printf "High_ks_test_CM_all.pdf")
finalPdf42=$(printf "High_ks_test_CM3PCS_k_chge_tst_only.pdf")
finalPdf43=$(printf "High_ks_test_CM4PCS_k_chge_tst_only.pdf")
finalPdf44=$(printf "High_ks_test_CM2PCS_k_chge_tst_only.pdf")
finalPdf45=$(printf "High_ks_test_CM1PCS_k_chge_tst_only.pdf")

finalPdf46=$(printf "High_ks_test_CM4PCS_frm_CM3PCS_help.pdf")
finalPdf47=$(printf "High_ks_test_CM2PCS_ExpctdL0_tst.pdf")

finalPdf48=$(printf "High_ks_test_CM2PCS_CurvtrWithPD.pdf")
finalPdf49=$(printf "High_ks_test_CM2PCS_CurvtrWithPD1.pdf")
finalPdf50=$(printf "High_ks_test_CM2PCS_CurvtrWithDM.pdf") #Direct_Measurement
finalPdf51=$(printf "High_ks_test_CM3PCS_CurvtrWithDM_Tr43.pdf")

finalPdf52=$(printf "High_ks_test_CM4PCS_CurvtrWithDM_Tr50.pdf")
finalPdf53=$(printf "High_ks_test_CM1PCS_CurvtrWithDM_Tr61.pdf")
finalPdf54=$(printf "High_ks_test_CM2PCS_CurvtrWithDM_Tr01.pdf")

finalPdf55=$(printf "High_ks_test_CM2ICS_UsingHighksPCS_Tr03.pdf")
finalPdf56=$(printf "High_ks_test_CM3ICS_UsingHighksPCS_Tr03.pdf")
finalPdf57=$(printf "High_ks_test_CM4ICS_UsingHighksPCS_Tr01.pdf")

finalPdf58=$(printf "VF.pdf")
finalPdf59=$(printf "VF_AllMembranes.pdf")
finalPdf60=$(printf "VF_AllMembranesAndVF.pdf")
finalPdf61=$(printf "VF16.pdf")

finalPdf62=$(printf "VF_wtSRyp3.pdf")
finalPdf63=$(printf "Tst_PulleyPointForce.pdf")
finalPdf64=$(printf "addedNCPair=3.pdf")
finalPdf65=$(printf "addedNCPair=4.pdf")
finalPdf66=$(printf "addedNCPair=10.pdf")

finalPdf67=$(printf "NCP_CrtclSrfc_11.pdf")
finalPdf68=$(printf "High_ks_test_VrtclForce_Sensitivity.pdf")

finalPdf69=$(printf "CgY_l0Var_pull1.pdf")
finalPdf70=$(printf "CgY_l0Var_pull3.pdf")
finalPdf71=$(printf "CgY_l0Var_pull4.pdf")
finalPdf72=$(printf "CgY_l0Var_pull5.pdf")
finalPdf73=$(printf "CgY_l0Var_pull6.pdf")
finalPdf74=$(printf "CgY_l0Var_pull7.pdf")
finalPdf75=$(printf "CgY_l0Var_pull8.pdf")

finalPdf76=$(printf "ProductionRunInitiatn_aftrFixingCM1PCS.pdf")
finalPdf77=$(printf "Prtrbtn_OnInitPCS.pdf")

finalPdf78=$(printf "chk_CM1PCS.pdf")
finalPdf79=$(printf "CM4PCS_pp.pdf")
finalPdf80=$(printf "CM2PCS_forceRmvAndRetainShapeTst.pdf")
finalPdf81=$(printf "CM1234PCS_forceRmvAndRetainShapeAftRescaling.pdf")
finalPdf82=$(printf "CM1234PCS_CgX0CgY0.pdf")
finalPdf83=$(printf "CM1234PCS_rescaledVrsn_AftrmvVrtclForceAndRetain.pdf")
finalPdf84=$(printf "CM1234PCS_prtrbtnTst_CgX0_CgY0_on_Making_PPbfraft_matchWithExp.pdf")
finalPdf85=$(printf "CM1PCS_pullingTst_aftMakingPP_eqvlntToImg.pdf")
finalPdf86=$(printf "CM1PCS_pullingTst_Equivlnt_TIIF.pdf")

finalPdf87=$(printf "CM4PCS_rndm.pdf")

finalPdf=$finalPdf87

ConnectingPdf=$(printf "Connecting.pdf")
CombinedPdf=$(printf "Combined.pdf")
WOBlank=$(printf "WithOutBlank.pdf")

connectingORnot=0

if [ $connectingORnot == 1 ]; then
    rm $ConnectingPdf
    blankPageNumMinus1=$((countStrt-1))
    echo $blankPageNumMinus1
    pdftk $finalPdf cat 1-$blankPageNumMinus1 output $WOBlank
    cp OneNodeNI.pdf $ConnectingPdf
    pdftk $WOBlank $ConnectingPdf $BlankPdf cat output $CombinedPdf
    mv $CombinedPdf $finalPdf
else
    cp OneNodeNI.pdf $finalPdf
    pdftk $finalPdf $BlankPdf cat output $CombinedPdf
    mv $CombinedPdf $finalPdf
fi

echo $finalPdf

rm Exp0_*.pdf
rm flNm_*.pdf

evince $finalPdf &
