#! /bin/bash

df1=fret_fpp'.dat'
df2=gam'.dat'
df3=inside_minimizer'.dat'
df4=AA0_inPress_strct1'.dat'
df5=AA0_inPress_strct2'.dat'

df6=ll0_inTnsnCmpr_strct1'.dat'
df7=ll0_inTnsnCmpr_strct2'.dat'

df8=ll0Lt_fLt_ftilt'.dat'
df9=Pressure_in_strct1'.dat'
df10=Pressure_in_strct2'.dat'

df11=TnsnCmpr_in_strct1'.dat'
df12=TnsnCmpr_in_strct2'.dat'

df13=adjust_l0'.dat'
df14=InStepVars'.dat'
df15=l0A0kskaCoordntes'.dat'
df16=karea_strct1'.dat'
df17=karea_strct2'.dat'

df18=step'.dat'
df19=grdBendCheck'.dat'
df20=store'.dat'
df21=area_serial'.dat'
df22=grdBendCheck'.dat'
df23=cotSq_and_Phi'.dat'
df24=Phi'.dat'
df25=All_Energy'.dat'
df26=All_PenFunc'.dat'

df27=LftRghtSideCell'.dat'
df28=ApclSpr_BfrAft'.dat'
df29=Step_for_Ift3'.dat'

df30=NodeVars'.dat'
df31=Entered_node.'dat'
df32=check_inside_areaE'.dat' 
df33=PresCrossCheck'.dat' 

df34=ActvRgn'.dat'
df35=reading_for_cycle'.dat'
df36=TimeStepEnergy'.dat'
df37=TimeVsEnergy'.dat'

df38=OriginDis'.dat'
df39=PulleyDis'.dat'
df40=FrameNoAtDiffStrct'.dat'
df41=doubleNodeList'.dat'
df42=read_strct_A0l0Cg'.dat'
df43=read_strct_A0l0Cg_withActvRgn'.dat'
df44=read_strct_A0l0Cg_wwoActvRgn'.dat'

df45=stageConv'dat'

df46=conV_cell_in_Cycls'.dat'
df47=conV_spr_in_Cycls'.dat'
df48=conVData_frm_strctPropsConv'.dat'

df49=Extrapolate'.dat'
df50=Extrapolate_MAE'.dat'
df51=CG_Energy'.dat'

df52=PT_arrays'.dat'
df53=SmthCellSprsArray'.dat'
df54=HmMltpl_Ltrl'.dat'

df55=save_ProgrsnStgProp'.dat'
df56=switchNI_and_save'.dat'
df57=switchNI_and_saveInsd'.dat'

df58=DN_chk'.dat'
df59=DiscrepencyInProgCycl'.dat'
df60=Frame_NI'.dat'
df61=AdjustPressChk'.dat'
df62=reduce_l0_of_Prtclr_cells_Spr'.dat'
df63=reduce_l0_of_Consec_cells_Spr'.dat'
df64=pressChkInLoop'.dat'
df65=PressAdjstmntNIChk'.dat'
df66=NI_switch_and_readStrctProps'.dat'
df67=reduce_l0forprcntChnge'.dat'
df68=rescale'.dat'

df69=maintain_dP'.dat'
df70=OriginDisNI'.dat'
df71=makingSamePrpSegmntd'.dat'
df72=sprCurs_and_sprNxts'.dat'
df73=pressure_adjstmnt_at_bottomCells'.dat'
df74=FrmChkBfrPressAdjstmntBotm'.dat'
df75=get_InsOutsForAC'.dat'

df76=chngA0_0fPairorSingl'.dat'
df77=chngl0_0fPairorSingl'.dat'
df78=area_Manpltn'.dat'

df79=A0l0_prcntChng'.dat'
df80=ManpltdSprArray'.dat'
df81=Area_Changes'.dat'

df82=chkTN_GRD'.dat'
df83=chkNI_GRD'.dat'

df84=diff_A0A_AreaPrpchk'.dat'
df85=diff_A0A_SprsPrpChk'.dat'
df86=keepIncrA0_ToMAtchPress'.dat'
df87=CurvatureComprsn'.dat'

rm $df1 $df2 $df3 $df4 $df5 $df6 $df7 $df8 $df9 $df10 $df11 $df12 $df13 $df14 $df15 $df16 $df17 $df18 $df19 $df20 $df21 $df22 $df23 $df24 $df25 $df26 $df27 $df28 $df29 $df30 $df31 $df32 $df33 $df34 $df35 $df36 $df37 $df38 $df39 $df40 $df41 $df42 $df43 $df44 $df45 $df46 $df47 $df48 $df49 $df50 $df51 $df52 $df53 $df54 $df55 $df56 $df57 $df59 $df60 $df61 $df62 $df63 $df64 $df65 $df67 $df68 $df69 $df70 $df71 $df72 $df73 $df74 $dfd75 $df76 $df77 $df78 $df79 $df80 $df81 $df82 $df83 $df84 $df85 $df86 $df87

project=cff #cff=cephallic furrow formation

f1=input_frm_bash_what_to_build #done
f2=info_modules #done
f3=conversion_for_moving_system #done
f4=diff_stageRoutines #done
f5=generating_system #done
f6=Manipulating_properties #done
f7=transfrm_variables #done
f8=area_to_other_trnsfrm #done
f9=node_to_other_trnsfrm #done
f10=spr_to_other_trnsfrm #done
f11=neighbourInfo_and_trnsfrmInfo_dependentInfo #done

f12=vitelline_fluid_calctn
f13=Energy_tree #done
f14=gradient_calcEnrgy #done
f15=gradient_calcPenF #done
f16=PenF_and_Wfunc #done

f17=mltiple_calls_together #done
f18=Minimizer_Routines_CGMethod #done
f19=PenF_minimizer #done
f20=Time_Stepping #done

f21=Force_calc_switching_off_mech #done
f22=optimizing_shape #done
f23=writing_sbrtns #done

f24=node_insertion_model #done
f25=nodeInsrtd_cntrlStates #done
f26=AddedCell_model #done
f27=ManipltngEquilibrateForInitiatn #done
f28=CS_match
f29=Manipulating_and_Equilibrate #done

f30=redefine_system_122123
f31=calls_for_tests #done

f32=propTrnsfer #done
f33=build_strctures #done
f34=readjust_two_step_trnsfrmtn #done
f35=multiple_cycle #done

f36=Fundamental_routines #done
f37=circle_center #done
f38=complt_simltn #done
f39=main_blcks #done


file1=${f1}'.f08'
file2=${f2}'.f08'
file3=${f3}'.f08'
file4=${f4}'.f08'
file5=${f5}'.f08'
file6=${f6}'.f08'
file7=${f7}'.f08'
file8=${f8}'.f08'
file9=${f9}'.f08'
file10=${f10}'.f08'
file11=${f11}'.f08'
file12=${f12}'.f08'
file13=${f13}'.f08'
file14=${f14}'.f08'
file15=${f15}'.f08'
file16=${f16}'.f08'
file17=${f17}'.f08'
file18=${f18}'.f08'
file19=${f19}'.f08'
file20=${f20}'.f08'
file21=${f21}'.f08'
file22=${f22}'.f08'
file23=${f23}'.f08'
file24=${f24}'.f08'
file25=${f25}'.f08'
file26=${f26}'.f08'
file27=${f27}'.f08'
file28=${f28}'.f08'
file29=${f29}'.f08'
file30=${f30}'.f08'
file31=${f31}'.f08'
file32=${f32}'.f08'
file33=${f33}'.f08'
file34=${f34}'.f08'
file35=${f35}'.f08'
file36=${f36}'.f08'
file37=${f37}'.f08'
file38=${f38}'.f08'
file39=${f39}'.f08'

executableFile=${project}'.out'

compiler='gfortran'

#gflags='-O -Wall -fcheck=all -g -fmax-errors=5 -fbacktrace -ffpe-trap=invalid,zero,overflow'
                  
gflags='-O -fcheck=all -g -fmax-errors=20 -fbacktrace -ffpe-trap=invalid,zero,overflow'

#gflags='-fcheck=all -fmax-errors=5 -fbacktrace -ffpe-trap=invalid,zero,overflow'

rm $executableFile

$compiler $gflags $file1 $file2 $file3 $file4 $file5 $file6 $file7 $file8 $file9 $file10 $file11 $file12 $file13 $file14 $file15 $file16 $file17 $file18 $file19 $file20 $file21 $file22 $file23 $file24 $file25 $file26 $file27 $file28 $file29 $file30 $file31 $file32 $file33 $file34 $file35 $file36 $file37 $file38 $file39 -o $executableFile


Analtcl_or_Numrcl=1

./$executableFile<<EOF


$Analtcl_or_Numrcl

EOF

#command; echo "Process Done" | mail -s "Simulation is done" redowan.a.niloy@ttu.edu
	 
