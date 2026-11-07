module diff_stage
  use system_parameters
  implicit none
  integer :: LastCycl=5
  
  integer :: LC,RC,CLTC,CRTC,CLBC,CRBC,EC !Cell No in Lft,Rght,ClftBot,CrghtBot for S1T1
  integer :: LC_finl,RC_finl,CLBC_finl,CRBC_finl
  
contains
  
  subroutine get_InitiatorCell
    implicit none
    
    if (stageNo==1 .and. stageType==1) then !stage3 may've cntrl prtn
       InitiatorCell = ncl+ncr+1   
    elseif (stageNo==4 .and. stageType==1) then 
       InitiatorCell = N_cell
    else
       write(*,*) "Initiator cell correction "
    endif
    
    !write(*,*) InitiatorCell,"Initiator Cell"
    
  end subroutine get_InitiatorCell

  subroutine S4_S1_conversion(InFn,convS,convA,convN)
    ! This is more like S4T1 to S1T1 conversion, if needed I have to write S4T2 to S1T2 conv
    ! Then I will change this naming, otherwise  it is fine.
    
    implicit none
    integer, intent(in)  :: InFn
    integer, intent(out) :: convS(1:N_spr),convA(1:N_cell),convN(1:(N_node+1))
    
    integer :: i,j
    integer :: nclS4,ncrS4,nccltS4,nccrtS4,ncclbS4,nccrbS4
    integer :: nslS4,nsrS4,nscltS4,nscrtS4,nsclbS4,nscrbS4
    integer :: nvslS4,nvsrS4,nvsclbS4,nvscrbS4
    integer :: nseccltS4,nseccrtS4,nsecclbS4 
    integer :: nneccltS4,nneccrtS4
    integer :: lim1,lim2,lim3,lim4,lim5,lim6
    
    integer :: strtSL,strtSR,cntSL,cntSR
    integer :: strtAL,strtAR,cntAL,cntAR
    integer :: strtNL,strtNR,cntNL,cntNR
    
    integer :: LftRghtInd,res,vagsesh
    
    open(unit=433,file='stageConv.dat',position='append')
    
    if (InFn==1) then ! 1 is for Initial
       nclS4   = 4 ; ncrS4   = 4
       nccltS4 = 1 ; nccrtS4 = 1
       ncclbS4 = 3 ; nccrbS4 = 3
       
    elseif (InFn==2) then ! 2 is for Final
       nclS4   = 3 ; ncrS4   = 3
       nccltS4 = 1 ; nccrtS4 = 1
       ncclbS4 = 4 ; nccrbS4 = 4
       
    endif

    nseccltS4 = 3
    nseccrtS4 = 3
    nsecclbS4 = 3

    nneccltS4 = 3
    nneccrtS4 = 3
    
    nslS4   = nclS4*nsecl
    nsrS4   = ncrS4*nsecr
    nscltS4 = nccltS4*nseccltS4
    nscrtS4 = nccrtS4*nseccrtS4
    nsclbS4 = ncclbS4*nsecclbS4
    nscrbS4 = nccrbS4*nsecclbS4
    
    lim1 = nslS4
    lim2 = lim1 + nsrS4
    lim3 = lim2 + nscltS4
    lim4 = lim3 + nscrtS4
    lim5 = lim4 + nsclbS4 + nscrbS4
    lim6 = N_spr
    
    write(433,*) lim1,lim2,lim3,lim4,lim5,lim6,"lims Spr"
    
    strtSL = 1 ; strtSR = nsl+1 ! not nslS4
    cntSL  = 0 ; cntSR  = 0
    write(433,*) strtSL,strtSR,"strtSL,strtSR"
    write(433,*) " "
    
    LftRghtInd = 1
    
    do i = 1,N_spr
       
       if (i.le.lim1) then
          convS(i) = strtSL + cntSL
          cntSL    = cntSL+1
          
       elseif (i.gt.lim1 .and. i.le.lim2) then
          convS(i) = strtSR + cntSR
          cntSR    = cntSR+1
          
       elseif (i.gt.lim2 .and. i.le.lim3) then
          
          if ((i-lim2)==1) then
             convS(i) = strtSL + cntSL
             
          elseif ((i-lim2)==2) then
             convS(i) = strtSL + (cntSL+2)
             
          elseif ((i-lim2)==3) then
             convS(i) = strtSL + (cntSL+1)
             cntSL = cntSL+nseccltS4
          endif
          
          
       elseif (i.gt.lim3 .and. i.le.lim4) then
          
          if ((i-lim3)==1) then
             convS(i) = strtSR + cntSR
             
          elseif ((i-lim3)==2) then
             convS(i) = strtSR + (cntSR+2)
             
          elseif ((i-lim3)==3) then
             convS(i) = strtSR + (cntSR+1)
             cntSR = cntSR+nseccrtS4
          endif
          
       elseif (i.gt.lim4 .and. i.le.lim5) then
          
          res = (LftRghtInd-1)/nsecclbS4
          vagsesh = mod(res,2)
          LftRghtInd = LftRghtInd+1
          
          if (vagsesh==0) then
             convS(i) = strtSL + cntSL
             cntSL    = cntSL+1
          elseif (vagsesh==1) then
             convS(i) = strtSR + cntSR
             cntSR    = cntSR+1
          endif   
          
       elseif (i.gt.lim5 .and. i.le.lim6) then
          convS(i) = i
       endif
       
       write(433,*) "i=",i, "and its convS=", convS(i)
       
    enddo
    
    write(433,*) " "
    
    lim1 = nclS4
    lim2 = lim1 + ncrS4
    lim3 = lim2 + nccltS4
    lim4 = lim3 + nccrtS4
    lim5 = lim4 + ncclbS4 + nccrbS4
    lim6 = N_cell
    
    write(433,*) lim1,lim2,lim3,lim4,lim5,lim6,"lims Area"
    
    strtAL = 1 ; strtAR = (ncl+1)
    cntAL  = 0 ; cntAR  = 0
    write(433,*) strtAL,strtAR,"strtAL,strtAR"
    write(433,*) " "
    
    do i = 1,N_cell
       
       if (i.le.lim1) then
          convA(i) = strtAL + cntAL
          cntAL    = cntAL+1
          
       elseif (i.gt.lim1 .and. i.le.lim2) then
          convA(i) = strtAR + cntAR
          cntAR    = cntAR+1
          
       elseif (i.gt.lim2 .and. i.le.lim3) then
          convA(i) = strtAL + cntAL
          cntAL    = cntAL+1
          
       elseif (i.gt.lim3 .and. i.le.lim4) then
          convA(i) = strtAR + cntAR
          cntAR    = cntAR+1
          
       elseif (i.gt.lim4 .and. i.le.lim5) then
          
          if (mod(i,2)==1) then
             convA(i) = strtAL + cntAL
             cntAL    = cntAL+1
          elseif (mod(i,2)==0) then
             convA(i) = strtAR + cntAR
             cntAR    = cntAR+1
          endif
          
       elseif (i.gt.lim5 .and. i.le.lim6) then
          convA(i) = i
       endif
       
       write(433,*) "i=",i, "and its convA=", convA(i)
       
    enddo
    
    write(433,*) " "
    
    nvslS4   = nclS4+1
    nvsrS4   = ncrS4+1
    nvsclbS4 = ncclbS4
    nvscrbS4 = nccrbS4
    
    
    lim1 = 2*nvslS4
    lim2 = lim1 + 2*nvsrS4
    lim3 = lim2 + nneccltS4
    lim4 = lim3 + nneccrtS4
    lim5 = lim4 + (2*nvsclbS4 + 2*nvscrbS4)

    write(433,*) lim1,lim2,lim3,lim4,lim5,"lims Node"
    
    
    strtNL = 1 ; strtNR = (2*nvsl) + 1 ! not nvslS4
    cntNL  = 0 ; cntNR  = 0
    write(433,*) strtNL,strtNR,"strtNL,strtNR"
    write(433,*) " "
    
    LftRghtInd = 1
    
    write(433,*) N_node,"N_node in S1T1"
    
    do i = 1,(N_node+1) !
       
       if (i.le.lim1) then
          convN(i) = strtNL + cntNL
          cntNL    = cntNL+1
          
       elseif (i.gt.lim1 .and. i.le. lim2) then
          convN(i) = strtNR + cntNR
          cntNR    = cntNR+1
          
       elseif (i.gt.lim2 .and. i.le. lim3) then
          
          if ((i-lim2)==1) then
             convN(i) = N_node
             
          elseif ((i-lim2)==2) then
             convN(i) = strtNL + cntNL
             cntNL    = cntNL+1
             
          elseif ((i-lim2)==3) then
             convN(i) = strtNL + cntNL
             cntNL    = cntNL+1
          endif
          
       elseif (i.gt.lim3 .and. i.le. lim4) then
          
          if ((i-lim3)==1) then
             convN(i) = N_node
             
          elseif ((i-lim3)==2) then
             convN(i) = strtNR + cntNR
             cntNR    = cntNR+1
             
          elseif ((i-lim3)==3) then
             convN(i) = strtNR + cntNR
             cntNR    = cntNR+1
          endif
          
       elseif (i.gt.lim4 .and. i.le.lim5) then

          res = (LftRghtInd-1)/2
          vagsesh = mod(res,2)
          LftRghtInd = LftRghtInd+1

          if (vagsesh==0) then
             convN(i) = strtNL + cntNL
             cntNL    = cntNL+1
          elseif (vagsesh==1) then
             convN(i) = strtNR + cntNR
             cntNR    = cntNR+1
          endif
          
       endif
       
       write(433,*) "i=",i, "and its convN=", convN(i)
       
    enddo
    
    close(433)
    
  end subroutine S4_S1_conversion
  
  subroutine convert_strctData_from_S4_to_S1
    implicit none
    real*8 :: ks_S1(1:N_spr), l0_S1(1:N_spr)
    real*8 :: ks_S4(1:N_spr), l0_S4(1:N_spr)
    real*8 :: ka_S1(1:N_cell),A0_S1(1:N_cell)
    real*8 :: ka_S4(1:N_cell),A0_S4(1:N_cell)
    real*8 :: CgXNode_S1(1:(N_node+1)),CgXNode_S4(1:(N_node+1))
    
    integer :: convA_SI(1:N_cell),convA_SF(1:N_cell)
    integer :: convS_SI(1:N_spr) ,convS_SF(1:N_spr)
    integer :: conVN_SI(1:(N_node+1)),convN_SF(1:(N_node+1))
    !N_node+1 coz in S1T1 before demolishing cntrlCell N_node=36 but we want space=38
    integer :: InFn
    
    integer :: i,j,jmax
    integer :: strctNum,strctMaxS4T1,N_itm
    integer :: ReadConv
    integer :: convSpr,convArea,convNode
    
    character(len=100) :: flnm
    character(len=100) :: flnmbr
    character(len=100) :: ReadFull_flnm
    character(len=100) :: ConvFull_flnm
    character(len=100) :: WriteFull_flnm
    
    convS_SI = 0 ; convS_SF = 0 ! SF = Structure Final
    convA_SI = 0 ; convA_SF = 0 ! SI = Structure Initial
    convN_SI = 0 ; convN_SF = 0 ! N~node,S~Spr,A~Area
    
    InFn = 1 ; call S4_S1_conversion(InFn,convS_SI,convA_SI,conVN_SI)
    InFn = 2 ; call S4_S1_conversion(InFn,convS_SF,convA_SF,convN_SF)
    
    strctMaxS4T1 = 3
    
    do strctNum = 1,strctMaxS4T1
       
       flnm='strctProps'
       write(flnmbr,'(i1.1,a)') strctNum,'S4T1.dat'
       write(*,*) trim(adjustl(flnm)),trim(adjustl(flnmbr))
       ReadFull_flnm =trim(adjustl(flnm))//trim(adjustl(flnmbr))
       
       flnm='strctPropsConv'
       write(flnmbr,'(i1.1,a)') strctNum,'S1T1.dat'
       write(*,*) trim(adjustl(flnm)),trim(adjustl(flnmbr))
       ConvFull_flnm=trim(adjustl(flnm))//trim(adjustl(flnmbr))

       flnm='WriteStrctPropsConv'
       write(flnmbr,'(i1.1,a)') strctNum,'S1T1.dat'
       write(*,*) trim(adjustl(flnm)),trim(adjustl(flnmbr))
       WriteFull_flnm=trim(adjustl(flnm))//trim(adjustl(flnmbr))
       
       open(unit=21,file=trim(adjustl(ReadFull_flnm)))
       open(unit=22,file=trim(adjustl(ConvFull_flnm)))
       open(unit=23,file=trim(adjustl(WriteFull_flnm)))
       
       N_itm = 3
       
       do ReadConv = 1,2
          
          do i=1,N_itm
             
             if (i==1) jmax=N_spr
             if (i==2) jmax=N_cell
             if (i==3) jmax=N_node+1
             
             
             do j = 1,jmax
                
                if (ReadConv==1) then !reading S4, no need of conversion
                   
                   if (i==1) read(21,*)  ks_S4(j),l0_S4(j)
                   if (i==2) read(21,*)  ka_S4(j),A0_S4(j)
                   if (i==3) read(21,*)  CgXNode_S4(j)
                   
                elseif (ReadConv==2) then !Finding prop for S1, conversion req
                   
                   if (i==1) then
                      
                      if (strctNum==1 .or. strctNum==3) convSpr = convS_SI(j)
                      if (strctNum==2) convSpr = convS_SF(j)
                      
                      ks_S1(convSpr) = ks_S4(j)
                      l0_S1(convSpr) = l0_S4(j)
                      
                   elseif (i==2) then
                      
                      if (strctNum==1 .or. strctNum==3) convArea = convA_SI(j)
                      if (strctNum==2) convArea = convA_SF(j)
                      
                      ka_S1(convArea) = ka_S4(j)
                      A0_S1(convArea) = A0_S4(j)
                      
                   elseif (i==3) then
                      
                      if (strctNum==1 .or. strctNum==3) convNode = convN_SI(j)
                      if (strctNum==2) convNode = convN_SF(j)
                      
                      CgXNode_S1(convNode) = CgXNode_S4(j)
                      
                   endif
                   
                endif
             enddo
             
          enddo
          
       enddo
       
       
       do i = 1,N_itm
          
          if (i==1) jmax=N_spr
          if (i==2) jmax=N_cell
          if (i==3) jmax=N_node
          
          do j = 1,jmax
             
             if (i==1) then
                write(22,*) ks_S1(j),l0_S1(j),j
                write(23,*) "sprNo =",j,"and ks,l0 are", ks_S1(j),l0_S1(j)
                
             elseif (i==2) then
                write(22,*) ka_S1(j),A0_S1(j),j
                write(23,*) "areaNo =",j,"and ka,A0 are", ka_S1(j),A0_S1(j)
                
             elseif (i==3) then
                write(22,*) CgXNode_S1(j)
                write(23,*) "nodeNo =",j,"and CgXNode is", CgXNode_S1(j)
             endif
             
          enddo
          
          write(22,*) " "
          write(23,*) " "
          
       enddo
       
       close(21)
       close(22)
       close(23)
       
    enddo
    
  end subroutine convert_strctData_from_S4_to_S1
  
  
  subroutine read_conVstrctProps_with_adjstmnt(conVstrct)
    implicit none
    integer, intent(in) :: conVstrct
    
    character(len=100) :: flnm1,flnm2
    character(len=100) :: flnmbr1,flnmbr2
    character(len=100) :: full_flnm
    
    integer :: i,j,jmax
    integer :: N_itm
    integer :: sprNo,cellNo
    integer :: ToBeReadCycl
    
    real*8  :: ks_Rff(1:N_spr) ,l0_Rff(1:N_spr) ! Rff = Read from file
    real*8  :: ka_Rff(1:N_cell),A0_Rff(1:N_cell) 
    real*8  :: CgXNode_Rff(1:N_node)
    
    integer :: conV_cell(1:N_cell),conV_spr(1:N_spr)
    
    if (stageNo.ne.1) then
       write(*,*) "check read_conVstrctProps_with_adjstment"
       stop
    endif
    
    if (CyclNo == LastCycl) then
       ToBeReadCycl = CyclNo
    elseif (CyclNo.ne.LastCycl) then
       ToBeReadCycl = CyclNo+1
    endif

    if (CyclNo.ne.LastCycl .and. conVstrct==2) then
       write(*,*) "Not Compatible for conVstrct=2 in other cycle except last cycle"
       stop
    endif
    
    flnm1='strctPropsCycl'
    write(flnmbr1,'(i1.1)') ToBeReadCycl
    flnm2='Conv'
    write(flnmbr2,'(i1.1,a)') conVstrct,'S1T1.dat'
    
    write(*,*)trim(adjustl(flnm1)),trim(adjustl(flnmbr1)),trim(adjustl(flnm2)),trim(adjustl(flnmbr2))
    full_flnm=trim(adjustl(flnm1))//trim(adjustl(flnmbr1))//trim(adjustl(flnm2))//trim(adjustl(flnmbr2))
    
    open(unit=21,file=trim(adjustl(full_flnm)))
    
    N_itm = 3
    
    do i=1,N_itm
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       if (i==3) jmax=N_node
       
       do j = 1,jmax
          if (i==1) read(21,*) ks_Rff(j),l0_Rff(j)
          if (i==2) read(21,*) ka_Rff(j),A0_Rff(j)
          if (i==3) read(21,*) CgXNode_Rff(j)
       enddo
       
    enddo
    
    close(21)
    
    call get_Cyclwise_conversion(conV_cell,conV_spr)
    
    open(unit=22,file='conVData_frm_strctPropsConv.dat',position='append')
    write(22,*) conVstrct,"conVstrct"
    
    do i =1,N_itm
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       if (i==3) jmax=N_node
       
       do j=1,jmax
          
          if (i==1) then
             sprNo    = conV_spr(j)
             k_spr(j) = ks_Rff(sprNo)
             l0(j)    = l0_Rff(sprNo)
             write(22,*) k_spr(j),l0(j),j,"ks,l0,spr"
             
          elseif (i==2) then
             cellNo    = conV_cell(j)
             k_area(j) = ka_Rff(cellNo)
             A0(j)     = A0_Rff(cellNo)
             write(22,*) k_area(j),A0(j),j,"ka,A0,cell"
             
          elseif (i==3) then
             CgXNode(j) = CgXNode_Rff(j)
             write(22,*) CgXNode(j),j,"CgXNode,node"
             
          endif
          
       enddo
       
       write(22,*) " "
    enddo
    
    close(22)
    
  end subroutine read_conVstrctProps_with_adjstmnt
  
  subroutine get_Cyclwise_conversion(conV_cell,conV_spr)
    implicit none
    integer :: CyclDiff
    integer, intent(out) :: conV_cell(1:N_cell),conV_spr(1:N_spr)
    
    !CyclDiff = LastCycl - CyclNo

    if (CyclNo==LastCycl) then
       CyclDiff=0
    elseif (CyclNo.ne.LastCycl) then
       CyclDiff = 1
    endif
    
    !Increasing previously, now 1(except last cycl), as reading from next cycle, not last cylc
    
    LC_finl   = 4 ; RC_finl   = 4 
    CLBC_finl = 3 ; CRBC_finl = 3
    
    LC   = LC_finl   + CyclDiff
    RC   = RC_finl   + CyclDiff
    
    if (CyclNo==1) then
       CLTC = 0
       CRTC = 0
    elseif (CyclNo.gt.1) then
       CLTC = 1
       CRTC = 1
    endif
    
    CLBC = CLBC_finl - CyclDiff
    CRBC = CLBC_finl - CyclDiff
    EC   = 1
    
    if (CLBC.lt.0) CLBC=0
    if (CRBC.lt.0) CRBC=0

    call get_conVcell(CyclDiff,conV_cell)
    call get_conVspr(CyclDiff,conV_spr)
    
  end subroutine get_Cyclwise_conversion
  
  subroutine get_conVcell(CyclDiff,conV_cell)
    implicit none
    integer, intent(in)  :: CyclDiff
    integer, intent(out) :: conV_cell(1:N_cell)
    
    integer :: i

    !This routn is for S1T1, wont work in S4T1 as ncl=4 in S4T1
    
    open(unit=654,file='conV_cell_in_Cycls.dat',position='append')
    write(654,*) ncl,ncr,N_cell,"CELLS"
    write(654,*) " "
    
    do i = 1,N_cell
       
       if (i.le.ncl) then
          
          if ((i-CyclDiff).le.1) then
             conV_cell(i) = 1
          else
             conV_cell(i) = i-CyclDiff
          endif
          
       elseif (i.gt.ncl .and. i.le. (ncl+ncr)) then
          
          if ((i-CyclDiff).le.(ncl+1)) then
             conV_cell(i) = ncl+1
          else
             conV_cell(i) = i-CyclDiff
          endif
          
       elseif (i==N_cell) then
          conV_cell(i) = i
       endif
       
       write(654,*) "i =",i,"and conV_cell(i) =",conV_cell(i)
       
    enddo
    
    close(654)
    
  end subroutine get_conVcell
  
  subroutine get_conVspr(CyclDiff,conV_spr)
    implicit none
    integer, intent(in)  :: CyclDiff
    integer, intent(out) :: conV_spr(1:N_spr)
    
    integer :: i,Initial_SimlSpr
    
    !This routn is for S1T1, wont work in S4T1 as nsl=12 in S4T1
    
    open(unit=654,file='conV_spr_in_Cycls.dat',position='append')
    
    Initial_SimlSpr = CyclDiff*nsecl
    write(654,*) nsl,nsr,N_spr,Initial_SimlSpr,"SPRINGS"

    do i = 1,N_spr
       
       if (i.le.nsl) then   
          
          if ((i-Initial_SimlSpr).le.nsecl) then
             if (mod(i,3)==1) conV_spr(i) = 1
             if (mod(i,3)==2) conV_spr(i) = 2
             if (mod(i,3)==0) conV_spr(i) = 3
          else 
             conV_spr(i) = i-Initial_SimlSpr 
          endif
          
       elseif (i.gt.nsl .and. i.le.(nsl+nsr)) then
          
          if ((i-Initial_SimlSpr).le.(nsl+nsecr)) then
             if (mod(i,3)==1) conV_spr(i) = nsl+1
             if (mod(i,3)==2) conV_spr(i) = nsl+2
             if (mod(i,3)==0) conV_spr(i) = nsl+3
          else
             conV_spr(i) = i-Initial_SimlSpr
          endif
          
       elseif (i==N_spr) then
          conV_spr(i) = i
       endif
       
       write(654,*) "i =",i,"and conV_spr(i) =",conV_spr(i)
    enddo
    
    close(654)
    
  end subroutine get_conVspr
  
end module diff_stage

