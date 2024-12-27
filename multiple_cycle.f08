module mltpl_cycle
  use build_struct !if needed will use readjust_struct
  use nodeInsrtd_cntrlStates
  use Adding_cells
  implicit none
  
  integer,allocatable :: actvCell(:)
  integer,allocatable :: actvSpr(:)
  integer,allocatable :: actvNode(:)
  
  !real*8, allocatable :: A0_stepF(:),ka_stepF(:) !with F=first, from strc 1-3
  !real*8, allocatable :: l0_stepF(:),ks_stepF(:) !if needed will design frm 3-2
  !The four matrices above are now in Manipulate_and_Equilibrate
  
contains
  
  subroutine destroy_and_rebuilt_for_cyclic_trnsfrmtn
    implicit none
    integer :: i,j,jmax,i1
    integer :: PullPushIndicator
    
    if (modelID==2) then
       write(*,*) "routine: get_destroy_and_rebuilt_for_cyclic_trnsfrmtn isn't for modelID=2"
       stop
    endif
    
    Exprmnt = 17
    Frame   = 1
    call destroy_struct
    
    write(*,*) "Check to Adjust strctNo in destroy_and_rebuilt_for_cyclic in multiple_cycle"
    !call sleep(1)
    
    strctNo = 4 !4
    call reinitialize_system_params
    call build_struct_for_cyclic_trnsfrmtn
    call Equilibrate_system
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
    Frame = Frame+1
    call switchto_NI_model_run_and_switchbackto_TN
    
    if (stageNo==4 .and. stageType==1) then
       nclMin=4 ; ncrMin = 4
       N_actvCell = 17 ; N_actvSpr = 49 ; N_actvNode = 4
    elseif (stageNo==4 .and. stageType==2) then
       nclMin=4 ; ncrMin = 4
       N_actvCell = 16 ; N_actvSpr = 48 ; N_actvNode = 4
    elseif (stageNo==1 .and. stageType==1) then
       write(*,*) "adjust in destroy and rebuilt in multiple_cycle.f08"
       stop
       nclMin=8 ; ncrMin = 8
       N_actvCell = 16 ; N_actvSpr = 50 ; N_actvNode = 4
    endif
    
    allocate(actvCell(1:N_actvCell))
    allocate(actvSpr(1:N_actvSpr))
    allocate(actvNode(1:N_actvNode))
    
    deallocate(ks_strctTN,l0_strctTN)
    deallocate(ka_strctTN,A0_strctTN)
    deallocate(CgX_strctTN,CgY_strctTN)
    
    allocate(ks_strctTN(1:N_struct,1:N_spr) ,l0_strctTN(1:N_struct,1:N_spr))
    allocate(ka_strctTN(1:N_struct,1:N_cell),A0_strctTN(1:N_struct,1:N_cell))
    allocate(CgX_strctTN(1:N_struct,1:N_node),CgY_strctTN(1:N_struct,1:N_node))
    
    deallocate(nodeXY_strctTN)
    allocate(nodeXY_strctTN(1:N_struct,1:N_node,1:2))
    
    nodeXY_strctTN = -1.d30
    
    do i=1,N_struct
       
       ks_strctTN(i,1:N_spr)   = k_spr(1:N_spr)
       l0_strctTN(i,1:N_spr)   = l0(1:N_spr)
       ka_strctTN(i,1:N_cell)  = k_area(1:N_cell)
       A0_strctTN(i,1:N_cell)  = A0(1:N_cell)
       CgX_strctTN(i,1:N_node) = 0.0d0
       CgY_strctTN(i,1:N_node) = 0.0d0
       
       if (i==1 .or. i==3) call get_active_region(ncl,i)
       if (i==2) call get_active_region(ncl,i)
       
       call read_strct_A0l0Cg(i)
       PullPushIndicator = 1
       call change_CgNode_forPullingPushing_test(i,PullPushIndicator)
       
       if (i==1.or.i==3) then
          call AreaSprCg_ofStrct_to_actual(i)
          call nodes_to_coordntes(node_xy,coordntes_xy)
          
          call Equilibrate_system
          call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
          Frame = Frame+1
          call switchto_NI_model_run_and_switchbackto_TN
       endif
       
       if (i==3) then
          !call adjust_critical_spr_l0_and_strctl0
          
          open(unit=88,file='PropChkinAlgChng3.dat',position='append')
          
          do i1 = 1,N_spr
             write(88,*) k_spr(i1),l0(i1),i1
          enddo
          write(88,*) " "
          do i1 = 1,N_cell
             write(88,*) k_area(i1),A0(i1),i1
          enddo
          write(88,*) " "
          write(88,*) " "
          close(88)
          
          allocate(A0_I(1:N_cell),ka_I(1:N_cell),l0_I(1:N_spr),ks_I(1:N_spr))
          allocate(A0_F(1:N_cell),ka_F(1:N_cell),l0_F(1:N_spr),ks_F(1:N_spr))
          
          A0_I = 1.0d30 ; ka_I = 1.0d30 ; l0_I = 1.0d30 ; ks_I = 1.0d30
          A0_F = 1.0d30 ; ka_F = 1.0d30 ; l0_F = 1.0d30 ; ks_F = 1.0d30

          A0_I(1:N_cell) = A0_strctTN(1,1:N_cell) ; ka_I(1:N_cell) = ka_strctTN(1,1:N_cell)
          l0_I(1:N_spr)  = l0_strctTN(1,1:N_spr)  ; ks_I(1:N_spr)  = ks_strctTN(1,1:N_spr)
          A0_F(1:N_cell) = A0_strctTN(3,1:N_cell) ; ka_F(1:N_cell) = ka_strctTN(3,1:N_cell)
          l0_F(1:N_spr)  = l0_strctTN(3,1:N_spr)  ; ks_F(1:N_spr)  = ks_strctTN(3,1:N_spr)
          
          
          call Extrapolate_and_adjust_strctProp_ifNeeded_in_MAE(Exprmnt,Frame)
          
          A0_strctTN(3,1:N_cell) = A0_F(1:N_cell) ; ka_strctTN(3,1:N_cell) = ka_F(1:N_cell)
          l0_strctTN(3,1:N_spr)  = l0_F(1:N_spr)  ; ks_strctTN(3,1:N_spr)  = ks_F(1:N_spr)
          
          deallocate(A0_I,ka_I,l0_I,ks_I)
          deallocate(A0_F,ka_F,l0_F,ks_F)
          
       endif
       
       if (i==1) then
          deallocate(coordntesXY_strct1)
          N_mvCoordnte1 = N_mvCoordnte
          allocate(coordntesXY_strct1(1:N_mvCoordnte1))
          coordntesXY_strct1 = coordntes_xy
          
          deallocate(coordntes_strct1)
          NMCWL0A01 = N_mvCoordnte_withl0A0
          allocate(coordntes_strct1(1:NMCWL0A01))
          coordntes_strct1 = coordntes        
       endif
       
       
       if (i==1.or.i==3) then
          do j = 1,N_node
             nodeXY_strctNI(i,j,1:2) = node_xy(j,1:2)
          enddo
       endif
       
    enddo
    
    !stop
    
    !do i = 1,9
       
       !if (i.le.4) jmax=N_spr
       !if ((i.gt.4) .and. (i.le.8)) jmax=N_cell
       !if ((i.gt.8) .and. (i.le.10)) jmax=N_node
       
       !do j = 1,jmax
        !  if (i==5) write(*,*) A0_strct(1,j),j
       !enddo
    
       !write(*,*) " "
    !enddo
    
  end subroutine destroy_and_rebuilt_for_cyclic_trnsfrmtn
  
  
  subroutine build_struct_for_cyclic_trnsfrmtn
    implicit none
    integer :: i
    
    call get_blcks
    call updt_kspr_with_l0
    call updt_karea_with_A0
    
    call get_all_the_transfrms
    call get_all_moving_coordnte_variables
    call get_all_gradient_variables
    call get_all_grdPenF_variables
    
    call print_the_system
    
  end subroutine build_struct_for_cyclic_trnsfrmtn
  
  
  subroutine get_active_region_arrays(iMain)
    implicit none
    integer, intent(in) :: iMain
    
    if (stageNo==1 .and. stageType==1) then
       continue
    else
       write(*,*) "edit get_active_region_arrays for diff stages"
       stop
    endif
    
    if (iMain .le. 1) then
       nclMin=8 ; ncrMin = 8
       N_actvCell = 17 ; N_actvSpr = 50 ; N_actvNode = 4
    elseif (iMain .gt. 1) then
       nclMin=8 ; ncrMin = 8
       N_actvCell = 17 ; N_actvSpr = 49 ; N_actvNode = 4
    endif
    
    allocate(actvCell(1:N_actvCell))
    allocate(actvSpr(1:N_actvSpr))
    allocate(actvNode(1:N_actvNode))
    
  end subroutine get_active_region_arrays
  
  
  subroutine deallocate_active_region_arrays
    implicit none
    
    deallocate(actvCell)
    deallocate(actvSpr)
    deallocate(actvNode)
    
  end subroutine deallocate_active_region_arrays
  
  
  subroutine get_active_region(nclVal,structNo) 
    implicit none
    integer, intent(in) :: nclVal
    integer, intent(in) :: structNo
    
    integer :: diff
    integer :: i,j,jmax
    integer :: cell_num,spr_num
    integer :: spr_cnt,node_cnt
    
    actvCell = 0 ; actvSpr = 0 ; actvNode = 0
    
    diff = nclVal-nclMin
    
    if (structNo==1 .or. structNo==3) then
       
       do i = 1,N_actvCell
          if (i.le.nclMin) then
             actvCell(i) = i+diff
          elseif (i.gt.nclMin.and.i.le.(nclMin+ncrMin)) then
             actvCell(i) = i+2*diff
          elseif ((i.gt.(nclMin+ncrMin)).and.i.le.(N_actvCell)) then
             actvCell(i) = i+2*diff
          endif
       enddo
       
       spr_cnt = 0
       
       do i = 1,N_actvCell
          cell_num = actvCell(i)
          
          do j = 1,3
             spr_num = 3*(cell_num-1) + j
             if (spr_num.gt.N_spr) exit
             spr_cnt = spr_cnt+1
             actvSpr(spr_cnt) = spr_num
          enddo
          
       enddo
       
       do i = 1,N_actvNode
          if (i==1 .or. i==2) then
             actvNode(i) = 2*diff+i
          elseif (i==3 .or. i==4) then
             actvNode(i) = 2*nvsl + i
          endif
       enddo
       
    elseif (structNo==2) then
       
       do i = 1,N_actvCell
          if (i.le.(nclMin-1)) then
             actvCell(i) = i+diff
          elseif(i.gt.(nclMin-1).and.i.le.(nclMin+ncrMin-2))then
             actvCell(i) = i+2*diff
          elseif(i.gt.(nclMin+ncrMin-2).and.i.le.N_actvCell)then
             actvCell(i) = i+2*diff
          endif
       enddo
       
       spr_cnt = 0
       
       do i = 1,N_actvCell
          cell_num = actvCell(i)
          
          do j = 1,3
             spr_num = 3*(cell_num-1) + j
             if (spr_num.gt.N_spr) exit
             spr_cnt = spr_cnt+1
             actvSpr(spr_cnt) = spr_num
          enddo
          
       enddo
       
       do i = 1,N_actvNode
          if (i==1 .or. i==2) then
             actvNode(i) = 2*diff+i
          elseif (i==3 .or. i==4) then
             actvNode(i) = 2*(nvsl-1) + i
          endif
       enddo
    endif
    
    open(unit=35,file='ActvRgn.dat',position='append')
    
    do i = 1,3
       if (i==1) jmax=N_actvCell
       if (i==2) jmax=N_actvSpr
       if (i==3) jmax=N_actvNode
       
       if (i==1) write(35,fmt=*) "Active Cell"," ","Cell Serial"
       if (i==2) write(35,fmt=*) "Active Spr"," ","Spr Serial"
       if (i==3) write(35,fmt=*) "Active Node"," ","Node Serial"
       
       do j=1,jmax
          if (i==1) write(35,fmt=*) actvCell(j),j
          if (i==2) write(35,fmt=*) actvSpr(j),j
          if (i==3) write(35,fmt=*) actvNode(j),j
       enddo
       
       write(35,fmt=*) " "
    enddo
    
    close(35)
    
  end subroutine get_active_region
  
  
  subroutine get_active_region_S1
    implicit none
    integer :: diff
    integer :: i,j,jmax
    integer :: cell_num,spr_num
    integer :: spr_cnt,node_cnt
    
    nclMin=8 ; ncrMin = 8
    
    actvCell = 0 ; actvSpr = 0 ; actvNode = 0
    diff = ncl-nclMin
    
    do i = 1,N_actvCell
       
       if (i.le.nclMin) then
          actvCell(i) = i+diff
       elseif (i.gt.nclMin.and.i.le.(nclMin+ncrMin)) then
          actvCell(i) = i+2*diff
       elseif (i.gt.(nclMin+ncrMin) .and. i.le.N_cell) then
          actvCell(i) = i+2*diff
       endif
       
    enddo
    
    spr_cnt = 0
    
    write(*,*) N_actvSpr,N_spr,"N_actvSpr,N_spr"
    
    do i = 1,N_actvCell
       cell_num = actvCell(i)
       
       do j = 1,3
          spr_num = 3*(cell_num-1) + j
          if (spr_num.gt.N_spr) exit
          spr_cnt = spr_cnt+1
          !write(*,*) spr_cnt,cell_num,i,j,"spr_cnt"
          actvSpr(spr_cnt) = spr_num
       enddo
       
    enddo
    
    do i = 1,N_actvNode
       if (i==1 .or. i==2) then
          actvNode(i) = 2*diff+i
       elseif (i==3 .or. i==4) then
          actvNode(i) = 2*nvsl + 2*diff + (i-2)
       endif
    enddo
    
    open(unit=35,file='ActvRgn.dat',position='append')
    
    do i = 1,3
       if (i==1) jmax=N_actvCell
       if (i==2) jmax=N_actvSpr
       if (i==3) jmax=N_actvNode
       
       if (i==1) write(35,fmt=*) "Active Cell"," ","Cell Serial"
       if (i==2) write(35,fmt=*) "Active Spr"," ","Spr Serial"
       if (i==3) write(35,fmt=*) "Active Node"," ","Node Serial"
       
       do j=1,jmax
          if (i==1) write(35,fmt=*) actvCell(j),j
          if (i==2) write(35,fmt=*) actvSpr(j),j
          if (i==3) write(35,fmt=*) actvNode(j),j
       enddo
       
       write(35,fmt=*) " "
    enddo
    
    close(35)
    
  end subroutine get_active_region_S1
  
  
  subroutine read_strct_A0l0Cg(strctNum)
    implicit none
    integer, intent(in) :: strctNum
    integer :: i,j,jmax
    integer :: N_itm
    integer :: spr_num,cell_num,node_num
    
    character(len=100) :: flnm
    character(len=100) :: flnmbr
    character(len=100) :: full_flnm
    
    real*8 :: CgX_valueRead

    if (modelID==2) then
       write(*,*) "routine: read_strct_A0l0Cg isn't for modelID=2"
       stop
    endif
    
    flnm='strctProps'
    N_itm = 3
    
    
    if (stageNo==1 .and. stageType==1) then
       
       if (strctNum.le.9) then
          write(flnmbr,'(i1.1,a)') strctNum,'S1T1.dat'
       elseif (strctNum.gt.9 .and. strctNum.le.99) then
          write(flnmbr,'(i2.2,a)') strctNum,'S1T1.dat'
       endif
       
    elseif (stageNo==4 .and. stageType==1) then
       write(flnmbr,'(i1.1,a)') strctNum,'S4T1.dat'
    elseif (stageNo==4 .and. stageType==2) then
       write(flnmbr,'(i1.1,a)') strctNum,'S4T2.dat'
    endif
    
    
    write(*,*) trim(adjustl(flnm)),trim(adjustl(flnmbr))
    full_flnm=trim(adjustl(flnm))//trim(adjustl(flnmbr))
    
    open(unit=21,file=trim(adjustl(full_flnm)))
    open(unit=22,file='read_strct_A0l0Cg.dat',position='append')
    
    
    if (strctNum==2) then
       do i = 1,nsecr
          spr_num = 3*(ncl-1)+i
          ks_strctTN(strctNum,spr_num) = k_spr(i)
          l0_strctTN(strctNum,spr_num) = l0(i)
       enddo
       
       ka_strctTN(strctNum,ncl) = k_area(1)
       A0_strctTN(strctNum,ncl) = A0(1)
    endif
    
    do i=1,N_itm
       
       if (i==1) jmax=N_actvSpr
       if (i==2) jmax=N_actvCell
       if (i==3) jmax=N_actvNode
       
       do j = 1,jmax
          if (i==1) then
             spr_num = actvSpr(j)
             read(21,*)  ks_strctTN(strctNum,spr_num),l0_strctTN(strctNum,spr_num)
             write(22,*) ks_strctTN(strctNum,spr_num),l0_strctTN(strctNum,spr_num)
             
          elseif (i==2) then
             cell_num = actvCell(j)
             read(21,*)  ka_strctTN(strctNum,cell_num),A0_strctTN(strctNum,cell_num)
             write(22,*) ka_strctTN(strctNum,cell_num),A0_strctTN(strctNum,cell_num)
             
          elseif (i==3) then
             if (j==1) read(21,*) CgX_valueRead
             node_num = actvNode(j)
             
             if(j.le.(N_actvNode/2)) then
                CgX_strctTN(strctNum,node_num) = CgX_valueRead
             elseif (j.gt.(N_actvNode/2)) then
                CgX_strctTN(strctNum,node_num) = -CgX_valueRead
             endif
             
          endif
       enddo
       write(22,*) " "
    enddo
    
    CgY_strctTN(strctNum,1:N_node) = 0.0d0
    
    close(21)
    close(22)
    
    open(unit=34,file='reading_for_cycle.dat',position='append')
    
    do i = 1,N_itm
       
       if (i==1) then
          jmax=N_spr
          write(34,*) "ks"," ","l0 "
       elseif (i==2) then
          jmax=N_cell
          write(34,*) "ka"," ","A0"
       elseif (i==3) then
          write(34,*) "CgX"
          jmax=N_node
       endif
       
       do j = 1,jmax
          if (i==1) write(34,*) ks_strctTN(strctNum,j),l0_strctTN(strctNum,j),j
          if (i==2) write(34,*) ka_strctTN(strctNum,j),A0_strctTN(strctNum,j),j
          if (i==3) write(34,*) CgX_strctTN(strctNum,j),CgY_strctTN(strctNum,j),j
       enddo
       write(34,*) " "
       
    enddo
    
    close(34)
    
  end subroutine read_strct_A0l0Cg
  
  
  subroutine read_strctPropsA0l0Cg_withActvRgn(strctNum)
    implicit none
    integer, intent(in) :: strctNum
    integer :: i,j,jmax
    integer :: N_itm
    integer :: spr_num,cell_num,node_num
    
    character(len=100) :: flnm
    character(len=100) :: flnmbr
    character(len=100) :: full_flnm
    
    real*8 :: CgX_valueRead,CgY_valueRead
    
    if (modelID==2) then
       write(*,*)"routine:read_strctPropsA0l0Cg_withActvRgn isn't for modelID=2"
       stop
    endif
    
    flnm='strctProps'
    N_itm = 3
    
    
    if (stageNo==1 .and. stageType==1) then
       
       if (strctNum.le.9) then
          write(flnmbr,'(i1.1,a)') strctNum,'S1T1.dat'
       elseif (strctNum.gt.9 .and. strctNum.le.99) then
          write(flnmbr,'(i2.2,a)') strctNum,'S1T1.dat'
       endif
       
    elseif (stageNo==4 .and. stageType==1) then
       write(flnmbr,'(i1.1,a)') strctNum,'S4T1.dat'
    elseif (stageNo==4 .and. stageType==2) then
       write(flnmbr,'(i1.1,a)') strctNum,'S4T2.dat'
       
    endif
    
    write(*,*) trim(adjustl(flnm)),trim(adjustl(flnmbr))
    full_flnm=trim(adjustl(flnm))//trim(adjustl(flnmbr))
    
    open(unit=61,file=trim(adjustl(full_flnm)))
    open(unit=62,file='read_strct_A0l0Cg_withActvRgn.dat',position='append')
    
    ks_strctTN(strctNum,1:N_spr)   = k_spr(1:N_spr)
    l0_strctTN(strctNum,1:N_spr)   = l0(1:N_spr)
    ka_strctTN(strctNum,1:N_cell)  = k_area(1:N_cell)
    A0_strctTN(strctNum,1:N_cell)  = A0(1:N_cell)
    CgX_strctTN(strctNum,1:N_node) = 0.00d0
    CgY_strctTN(strctNum,1:N_node) = 0.00d0
    
    call get_non_actvCell_StrctProps(strctNum)
    
    do i=1,N_itm
       
       if (i==1) jmax=N_actvSpr
       if (i==2) jmax=N_actvCell
       if (i==3) jmax=N_actvNode
       
       if (i==1) write(62,*) strctNum,"strct Num"
       
       do j = 1,jmax
          
          if (i==1) then
             spr_num = actvSpr(j)
             read(61,*)  ks_strctTN(strctNum,spr_num),l0_strctTN(strctNum,spr_num)
             write(62,*) ks_strctTN(strctNum,spr_num),l0_strctTN(strctNum,spr_num),spr_num
             
          elseif (i==2) then
             cell_num = actvCell(j)
             read(61,*)  ka_strctTN(strctNum,cell_num),A0_strctTN(strctNum,cell_num)
             write(62,*) ka_strctTN(strctNum,cell_num),A0_strctTN(strctNum,cell_num),cell_num
             
          elseif (i==3) then
             if (j==1) read(61,*) CgX_valueRead,CgY_valueRead
             node_num = actvNode(j)
             
             if(j.le.(N_actvNode/2)) then
                CgX_strctTN(strctNum,node_num) = CgX_valueRead
             elseif (j.gt.(N_actvNode/2)) then
                CgX_strctTN(strctNum,node_num) = -CgX_valueRead
             endif
             
             write(62,*) CgX_strctTN(strctNum,node_num),CgY_strctTN(strctNum,node_num),node_num
             
          endif
       enddo
       
       write(62,*) " "
       
    enddo
    
    close(61)
    close(62)
    
    open(unit=63,file='read_strct_A0l0Cg_wwoActvRgn.dat',position='append')
    
    do i=1,N_itm
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       if (i==3) jmax=N_node

       if (i==1) write(63,*) strctNum,"strct Num"
       
       do j = 1,jmax
          
          if (i==1) then
             write(63,*) ks_strctTN(strctNum,j),l0_strctTN(strctNum,j),j
          elseif (i==2) then
             write(63,*) ka_strctTN(strctNum,j),A0_strctTN(strctNum,j),j
          elseif (i==3) then
             write(63,*) CgX_strctTN(strctNum,j),CgY_strctTN(strctNum,j),j
          endif
          
       enddo
       
       write(63,*) " "
       
    enddo
    
    close(63)
    
  end subroutine read_strctPropsA0l0Cg_withActvRgn
  
  
  subroutine read_strctPropsA0l0Cg_Updtd(strctNum)
    implicit none
    integer, intent(in) :: strctNum
    
    character(len=100) :: flnmPrf1,flnmPrf2,flnmbr
    character(len=100) :: fullflnm

    integer :: i,j,jmax
    integer :: N_itm
    
    flnmPrf1='strctPropsAddedCell'
    
    if (modelID==1) flnmPrf2='TN'
    if (modelID==2) flnmPrf2='NI'
    
    N_itm = 3
    
    if (strctNum.le.9) then
       write(flnmbr,'(i1.1,a)') strctNum,'S1T1.dat'
    elseif ((strctNum.gt.9) .and. (strctNum.le.99)) then
       write(flnmbr,'(i2.2,a)') strctNum,'S1T1.dat'
    endif
    
    write(*,*)trim(adjustl(flnmPrf1))//trim(adjustl(flnmPrf2))//trim(adjustl(flnmbr))," WritingUPDT"
    fullflnm=trim(adjustl(flnmPrf1))//trim(adjustl(flnmPrf2))//trim(adjustl(flnmbr))
    
    open(unit=43,file=trim(adjustl(fullflnm)))
    
    if (modelID==1) then
       
       ks_strctTN(strctNum,1:N_spr)   = k_spr(1:N_spr)
       l0_strctTN(strctNum,1:N_spr)   = l0(1:N_spr)
       ka_strctTN(strctNum,1:N_cell)  = k_area(1:N_cell)
       A0_strctTN(strctNum,1:N_cell)  = A0(1:N_cell)
       CgX_strctTN(strctNum,1:N_node) = 0.00d0
       CgY_strctTN(strctNum,1:N_node) = 0.00d0
       
    elseif (modelID==2) then
       
       ks_strctNI(strctNum,1:N_spr)   = k_spr(1:N_spr)
       l0_strctNI(strctNum,1:N_spr)   = l0(1:N_spr)
       ka_strctNI(strctNum,1:N_cell)  = k_area(1:N_cell)
       A0_strctNI(strctNum,1:N_cell)  = A0(1:N_cell)
       CgX_strctNI(strctNum,1:N_node) = 0.00d0
       CgY_strctNI(strctNum,1:N_node) = 0.00d0
       
    endif
    
    
    do i = 1,N_itm
       
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
       if (i==3) jmax = N_node
       
       do j = 1,jmax
          
          if (modelID==1) then
             
             if (i==1) read(43,*) ks_strctTN(strctNum,j),l0_strctTN(strctNum,j)
             if (i==2) read(43,*) ka_strctTN(strctNum,j),A0_strctTN(strctNum,j)
             if (i==3) read(43,*) CgX_strctTN(strctNum,j),CgY_strctTN(strctNum,j)
             
          elseif (modelID==2) then
             
             if (i==1) read(43,*) ks_strctNI(strctNum,j),l0_strctNI(strctNum,j)
             if (i==2) read(43,*) ka_strctNI(strctNum,j),A0_strctNI(strctNum,j)
             if (i==3) read(43,*) CgX_strctNI(strctNum,j),CgY_strctNI(strctNum,j)
             
          endif
          
       enddo
       
    enddo
    
    close(43)
    
  end subroutine read_strctPropsA0l0Cg_Updtd
  
  subroutine get_strctProps_readFrmA0l0Cg_Updtd(strctNum)
    implicit none
    integer, intent(in) :: strctNum
    
    if (modelID==1) then
       
       k_spr(1:N_spr)    = ks_strctTN(strctNum,1:N_spr)
       l0(1:N_spr)       = l0_strctTN(strctNum,1:N_spr)
       k_area(1:N_cell)  = ka_strctTN(strctNum,1:N_cell)
       A0(1:N_cell)      = A0_strctTN(strctNum,1:N_cell)
       CgXNode(1:N_node) = CgX_strctTN(strctNum,1:N_node)
       CgYNode(1:N_node) = CgY_strctTN(strctNum,1:N_node)
       
    elseif (modelID==2) then
       
       k_spr(1:N_spr)    = ks_strctNI(strctNum,1:N_spr)
       l0(1:N_spr)       = l0_strctNI(strctNum,1:N_spr)
       k_area(1:N_cell)  = ka_strctNI(strctNum,1:N_cell)
       A0(1:N_cell)      = A0_strctNI(strctNum,1:N_cell)
       CgXNode(1:N_node) = CgX_strctNI(strctNum,1:N_node)
       CgYNode(1:N_node) = CgY_strctNI(strctNum,1:N_node)
       
    endif
    
  end subroutine get_strctProps_readFrmA0l0Cg_Updtd
  
  subroutine read_strctPropsA0l0Cg_FrmStrcBfrMeet(strctBfrMeet,strctAftMeet)
    implicit none
    integer, intent(in) :: strctBfrMeet,strctAftMeet
    character(len=100)  :: flnmPrf1,flnmPrf2
    character(len=100)  :: flnmPrv,flnmAft
    character(len=100)  :: fullflnm_strctBfrMeet,fullflnm_strctAftMeet,saveFlnm
    
    integer :: N_itm
    integer :: cntrlSprNm,nsprsInACell
    integer :: i,j,jmax
    integer :: sprNm,nodeNm
    
    integer :: nsprsBfrMeet,nsprsAftMeet
    integer :: TNnodesBfrMeet,TNnodesAftMeet
    integer :: sprNmBfr,areaNmBfr,nodeNmBfr
    integer :: sprNmAft,areaNmAft,nodeNmAft
    integer :: sgmntdOrNot,whereTo
    
    real*8, allocatable :: ksFrmBfrMeet(:),ksFrmAftMeet(:)
    real*8, allocatable :: l0FrmBfrMeet(:),l0FrmAftMeet(:)
    real*8, allocatable :: kaFrmBfrMeet(:),kaFrmAftMeet(:)
    real*8, allocatable :: A0FrmBfrMeet(:),A0FrmAftMeet(:)
    real*8, allocatable :: CgXFrmBfrMeet(:),CgXFrmAftMeet(:)
    
    if (modelID==1) then
       continue
    elseif (modelID==2) then
       write(*,*) "stop Reason: Beginning of this sbrtn can't be 2"
       write(*,*) "sb:read_strctPropsA0l0Cg_FrmStrcBfrMeet,fl:multiple_cycle.f08"
       stop
    endif
    
    flnmPrf1='strctPropsAddedCell'
    flnmPrf2='TN'
    
    write(flnmPrv,'(i1.1,a)') strctBfrMeet,'S1T1.dat'
    write(flnmAft,'(i1.1,a)') strctAftMeet,'S1T1.dat'
    
    N_itm = 4
    
    fullflnm_strctBfrMeet=trim(adjustl(flnmPrf1))//trim(adjustl(flnmPrf2))//trim(adjustl(flnmPrv))
    fullflnm_strctAftMeet=trim(adjustl(flnmPrf1))//trim(adjustl(flnmPrf2))//trim(adjustl(flnmAft))
    
    write(*,*) trim(adjustl(fullflnm_strctBfrMeet)),"strctBfrMeet name"
    write(*,*) trim(adjustl(fullflnm_strctAftMeet)),"strctAftMeet name"
    
    open(unit=891,file=trim(adjustl(fullflnm_strctBfrMeet)))
    open(unit=892,file=trim(adjustl(fullflnm_strctAftMeet)))
    
    nsprsInACell = 3
    cntrlSprNm   = (N_cell-1)*(nsprsInACell) + 1
    
    nsprsBfrMeet = (N_cell-1)*(nsprsInACell) + (1) + (1) 
    nsprsAftMeet = (N_cell-1)*(nsprsInACell) + (1)
    write(*,*) nsprsBfrMeet,nsprsAftMeet,"nsprs bfr and aft meet"
    write(*,*) cntrlSprNm,"cntrlSprNm"

    TNnodesBfrMeet = (((Hlf_Ncell+1)*2)*2)
    TNnodesAftMeet = (((Hlf_Ncell+1)*2)*2)+1
    write(*,*) TNnodesBfrMeet,TNnodesAftMeet,"TNnodes Bfr-Aft"
    
    allocate(ksFrmBfrMeet(1:nsprsBfrMeet),ksFrmAftMeet(1:nsprsAftMeet))
    allocate(l0FrmBfrMeet(1:nsprsBfrMeet),l0FrmAftMeet(1:nsprsAftMeet))
    
    ksFrmBfrMeet=-1.0d30 ; ksFrmAftMeet=-1.0d30
    l0FrmBfrMeet=-1.0d30 ; l0FrmAftMeet=-1.0d30
    
    allocate(kaFrmBfrMeet(1:N_cell),kaFrmAftMeet(1:N_cell))
    allocate(A0FrmBfrMeet(1:N_cell),A0FrmAftMeet(1:N_cell))
    
    kaFrmBfrMeet=-1.0d30 ; kaFrmAftMeet=-1.0d30
    A0FrmBfrMeet=-1.0d30 ; A0FrmAftMeet=-1.0d30
    
    allocate(CgXFrmBfrMeet(1:TNnodesBfrMeet),CgXFrmAftMeet(1:TNnodesAftMeet))
    CgXFrmBfrMeet=-1.0d30 ; CgXFrmAftMeet=-1.0d30
    CgYNode(1:N_node) = 0.0d0
    
    do i = 1,N_itm
       
       if (i==1) jmax = nsprsBfrMeet
       if (i==2) jmax = N_cell
       if (i==3) jmax = TNnodesBfrMeet
       if (i==4) jmax = N_node
       
       do j = 1,jmax
          
          if (i==1) then
             
             sprNmBfr = j
             
             if (sprNmBfr .lt. cntrlSprNm) then
                
                read(891,*) ksFrmBfrMeet(sprNmBfr),l0FrmBfrMeet(sprNmBfr)
                sprNmAft = sprNmBfr
                
                ksFrmAftMeet(sprNmAft)=ksFrmBfrMeet(sprNmBfr)
                l0FrmAftMeet(sprNmAft)=l0FrmBfrMeet(sprNmBfr)
                
                !write(*,*)  ksFrmBfrMeet(sprNmBfr),l0FrmBfrMeet(sprNmBfr),sprNmBfr,"ks+l0+Bfr"
                !write(*,*)  ksFrmAftMeet(sprNmAft),l0FrmAftMeet(sprNmAft),sprNmAft,"ks+l0+Aft"
                write(892,*) ksFrmAftMeet(sprNmAft),l0FrmAftMeet(sprNmAft)
                
             elseif (sprNmBfr == cntrlSprNm) then
                
                read(891,*) ksFrmBfrMeet(sprNmBfr),l0FrmBfrMeet(sprNmBfr)
                !write(*,*)  ksFrmBfrMeet(sprNmBfr),l0FrmBfrMeet(sprNmBfr),cntrlSprNm,"ks+l0+cntrlSpr"
                
             elseif (sprNmBfr .gt. cntrlSprNm) then
                
                read(891,*) ksFrmBfrMeet(sprNmBfr),l0FrmBfrMeet(sprNmBfr)
                sprNmAft = sprNmBfr-1
                !write(*,*) NAEC_Apcl,NAEC_Bsal,NAEC_Ltrl,"Value of NAEC's when model ID =",modelID
                
                
                ksFrmAftMeet(sprNmAft) = ksFrmBfrMeet(sprNmBfr)
                l0FrmAftMeet(sprNmAft) = l0FrmBfrMeet(sprNmBfr)
                
                !write(*,*)  ksFrmBfrMeet(sprNmBfr),l0FrmBfrMeet(sprNmBfr),sprNmBfr,"ks+l0+Bfr"
                !write(*,*)  ksFrmAftMeet(sprNmAft),l0FrmAftMeet(sprNmAft),sprNmAft,"ks+l0+Aft"
                write(892,*) ksFrmAftMeet(sprNmAft),l0FrmAftMeet(sprNmAft)
                
             endif
             
          elseif (i==2) then
             
             areaNmBfr = j
             
             read(891,*) kaFrmBfrMeet(areaNmBfr),A0FrmBfrMeet(areaNmBfr)
             areaNmAft = areaNmBfr
             kaFrmAftMeet(areaNmAft) = kaFrmBfrMeet(areaNmBfr)
             A0FrmAftMeet(areaNmAft) = A0FrmBfrMeet(areaNmBfr)
             
             !write(*,*)  kaFrmBfrMeet(areaNmBfr),A0FrmBfrMeet(areaNmBfr),areaNmBfr,"ka+A0+Bfr"
             !write(*,*)  kaFrmAftMeet(areaNmAft),A0FrmAftMeet(areaNmAft),areaNmAft,"ka+A0+Aft"
             write(892,*) kaFrmAftMeet(areaNmAft),A0FrmAftMeet(areaNmAft)
             
          elseif (i==3) then
             
             nodeNmBfr = j
             read(891,*) CgXFrmBfrMeet(nodeNmBfr)
             nodeNmAft = nodeNmBfr
             CgXFrmAftMeet(nodeNmAft) = CgXFrmBfrMeet(nodeNmBfr)
             !write(*,*)   CgXFrmBfrMeet(nodeNmBfr),nodeNmBfr,"CgX+Bfr"
             !write(*,*)   CgXFrmAftMeet(nodeNmAft),nodeNmAft,"CgX+Aft"
             write(892,*) CgXFrmAftMeet(nodeNmAft),CgYNode(nodeNmAft)
             
             if (nodeNmBfr==TNnodesBfrMeet) then
                CgXFrmAftMeet(nodeNmAft+1) = CgXFrmBfrMeet(nodeNmBfr)
                !write(*,*)   CgXFrmBfrMeet(nodeNmBfr),   nodeNmBfr,   "CgX+Bfr"
                !write(*,*)   CgXFrmAftMeet(nodeNmAft+1),(nodeNmAft+1),"CgX+Aft+1"
                write(892,*) CgXFrmAftMeet(nodeNmAft+1),CgYNode(nodeNmAft+1)
             endif
             
             
          elseif (i==4) then
             
             if (j==1) then
                write(*,*) N_spr,nsprsAftMeet,"sprs"
                write(*,*) N_node,TNnodesAftMeet,"nodes"
                
                call sleep(2)
                
                k_spr(1:N_spr)    = ksFrmAftMeet(1:N_spr)
                l0(1:N_spr)       = l0FrmAftMeet(1:N_spr)
                k_area(1:N_cell)  = kaFrmAftMeet(1:N_cell)
                A0(1:N_cell)      = A0FrmAftMeet(1:N_cell)
                CgXNode(1:N_node) = CgXFrmAftMeet(1:N_node)
                CgYNode(1:N_node) = 0.0d0
                
                ks_strctTN(strctAftMeet,1:N_spr)   = k_spr(1:N_spr)
                l0_strctTN(strctAftMeet,1:N_spr)   = l0(1:N_spr)
                ka_strctTN(strctAftMeet,1:N_cell)  = k_area(1:N_cell)
                A0_strctTN(strctAftMeet,1:N_cell)  = A0(1:N_cell)
                CgX_strctTN(strctAftMeet,1:N_node) = CgXNode(1:N_node)
                CgY_strctTN(strctAftMeet,1:N_node) = CgYNode(1:N_node)
                
                call Equilibrate_system
                call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
                Frame=Frame+1
             endif
             
             write(892,*) node_xy(j,1:N_dmnsn)
             
             if (j==jmax) then
                flnmPrf2='NI'
                write(flnmAft,'(i1.1,a)') strctAftMeet,'S1T1.dat'
                saveFlnm=trim(adjustl(flnmPrf1))//trim(adjustl(flnmPrf2))//trim(adjustl(flnmAft))
                write(*,*) trim(adjustl(saveFlnm)),"saveFlnm"
                
                sgmntdOrNot=1 ; whereTo=1
                call switchto_NI_model_run_save_switchbackto_TN_seprtlyOrtogthr(&
                     saveFlnm,sgmntdOrNot,whereTo)
                
                ks_strctNI(strctAftMeet,1:N_spr)   = k_spr(1:N_spr)
                l0_strctNI(strctAftMeet,1:N_spr)   = l0(1:N_spr)
                ka_strctNI(strctAftMeet,1:N_cell)  = k_area(1:N_cell)
                A0_strctNI(strctAftMeet,1:N_cell)  = A0(1:N_cell)
                CgX_strctNI(strctAftMeet,1:N_node) = CgXNode(1:N_node)
                CgY_strctNI(strctAftMeet,1:N_node) = CgYNode(1:N_node)
                
                sgmntdOrNot=1 ; whereTo=2
                call switchto_NI_model_run_save_switchbackto_TN_seprtlyOrtogthr(&
                     saveFlnm,sgmntdOrNot,whereTo)
                
             endif
          
          endif
          
       enddo
       
       write(892,*) " "
       
    enddo 
    
    close(891)
    close(892)
    
    
  end subroutine read_strctPropsA0l0Cg_FrmStrcBfrMeet
  
  
  subroutine get_non_actvCell_StrctProps(strctNum)
    implicit none
    integer, intent(in) :: strctNum
    
    integer :: N_nonactvCell
    real*8  :: ka_FAC,A0_FAC ! FAC=First Active Cell
    integer :: FAC,i
    real*8  :: rltv_FvalFAC,rltv_FvalIncrSize,rltv_FvalIncr
    real*8  :: fmax,fval
    integer :: Lcell,Rcell

    if (modelID==2) then
       write(*,*) "routine: get_non_actvCell_StrctProps isn't for modelID=2"
       stop
    endif
    
    if (stageNo==1 .and. stageType==1) then
       continue
    else
       write(*,*) "May need to check the subroutine for nclMin,rltv_FvalFAC,fmax, fl:multiple_cycle, sb:get_non_actvCellProps"
       
    endif
    
    N_nonactvCell = ncl-nclMin
    
    open(unit=54,file='nonActvCell.dat')
    write(54,*) N_nonactvCell,"N_nonactvCell"
    
    if (N_nonactvCell==0) then
       write(*,*) "N_actvCell can't be zero(0) in this subroutine"
       stop
    endif
    
    ka_FAC = 0.0d0 ; A0_FAC = 0.0d0
    FAC    = N_nonactvCell + 1
    
    ka_FAC = k_area(FAC) ; A0_FAC = A0(FAC)
    write(54,*) ka_FAC,A0_FAC,"ka,A0 of FAC"
    
    rltv_FvalFAC  = 0.80d0
    rltv_FvalIncrSize = 0.05d0
    fmax = 1.20d0
    
    do i = 1,N_nonactvCell
       
       rltv_FvalIncr = (N_nonactvCell-(i-1))*rltv_FvalIncrSize
       fval = 1.0d0 + (rltv_FvalFAC+rltv_FvalIncr)*(fmax-1.0d0)
       
       write(54,*) rltv_FvalIncr,fval,i,"rltv_FvalIncr"
       
       Lcell = i ; Rcell = Hlf_Ncell+i
       write(54,*) Lcell,Rcell,"Lcell,Rcell"
       
       A0_strctTN(strctNum,Lcell) = (A0_FAC)*(fval)
       ka_strctTN(strctNum,Lcell) = (ka_FAC)!/(fval)
       write(54,*) ka_strctTN(strctNum,Lcell), A0_strctTN(strctNum,Lcell),"ka,A0 of Lcell"
       
       A0_strctTN(strctNum,Rcell) = (A0_FAC)*(fval)
       ka_strctTN(strctNum,Rcell) = (ka_FAC)!/(fval)
       write(54,*) ka_strctTN(strctNum,Rcell), A0_strctTN(strctNum,Rcell),"ka,A0 of Rcell"
       
    enddo
    
    close(54)
    !stop
    
  end subroutine get_non_actvCell_StrctProps
  
  subroutine AreaSprCg_ofStrct_to_actual(struct_No)
    implicit none
    integer, intent(in) :: struct_No

    if (modelID==2) then
       write(*,*) "routine: AreaSprCg_ofStrct_to_actual isn't for modelID=2"
       stop
    endif
    
    l0(1:N_spr)       = l0_strctTN(struct_No,1:N_spr)
    k_spr(1:N_spr)    = ks_strctTN(struct_No,1:N_spr)
    A0(1:N_cell)      = A0_strctTN(struct_No,1:N_cell)
    k_area(1:N_cell)  = ka_strctTN(struct_No,1:N_cell)
    CgXNode(1:N_node) = CgX_strctTN(struct_No,1:N_node)
    CgYNode(1:N_node) = 0.0d0
    
  end subroutine AreaSprCg_ofStrct_to_actual
  
  subroutine change_CgNode_forPullingPushing_test(structNum,PullPush)
    implicit none
    integer, intent(in) :: structNum
    integer, intent(in) :: PullPush
    integer :: i
    real*8  :: TinyValue
    
    if (modelID==2) then
       write(*,*) "routine: change_CgNode_forPullingPushing_test isn't for modelID=2"
       stop
    endif
    
    TinyValue = 1.0d-15
    write(*,*) structNum,PullPush,"structNum,PullPush"
    
    do i = 1,N_node
       !write(*,*) CgX_strctTN(structNum,i),i,"CgX_strctTN"
       
       if (abs(CgX_strctTN(structNum,i)) .gt. TinyValue) then
          !write(*,*) i,"Entering"
          if (CgX_strctTN(structNum,i) .gt. TinyValue) then!(+)ve val for Pull
             if (PullPush==1)    CgX_strctTN(structNum,i) = 0.15d0
             if (PullPush==(-1)) CgX_strctTN(structNum,i) = -0.15d0 
             
          elseif (CgX_strctTN(structNum,i).lt.TinyValue) then!(-)ve val forPull
             if (PullPush==1)    CgX_strctTN(structNum,i) = -0.15d0
             if (PullPush==(-1)) CgX_strctTN(structNum,i) = 0.15d0
          endif
          
       endif

       !write(*,*) CgX_strctTN(structNum,i),i,"CgX_strctTN"
       
    enddo
    !stop
    
  end subroutine change_CgNode_forPullingPushing_test
  
  subroutine get_critical_spr_tension
    implicit none
    continue
  end subroutine get_critical_spr_tension
  
  subroutine correct_strcture_props
    implicit none
    continue
  end subroutine correct_strcture_props
  
  subroutine adjust_critical_spr_l0_and_strctl0
    implicit none
    
    integer :: CP1(1:2)
    real*8  :: ratio,ratioChng
    integer :: i
    real*8  :: tol_Rdfn
    logical :: lgcl_rdfn

    if (modelID==2) then
       write(*,*) "routine: adjust_critical_spr_l0_and_strctl0 isn't for modelID=2"
       stop
    endif
    
    open(unit=87,file='adjstCritcspr.dat')
    
    tol_Rdfn  = 0.05d0
    lgcl_rdfn = .False.
    
    call get_first_critical_pair_spr(CP1) !Pulley_type=6, edit to use for S1T1 
    write(87,fmt=*) CP1(1:2),"CP1"
    
    ratio = l0(CP1(1))/initial_l0_ofCP
    ratioChng = 0.05d0
    
    write(87,fmt=*) ratio,"ratio bfr"
    write(87,fmt=*) initial_l0_ofCP,"initl0_ofCP"
    
    call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
    write(87,fmt=*) lgcl_rdfn,"lgcl"
    write(87,fmt=*) node_xy(11,1:2),"11"
    write(87,fmt=*) node_xy(23,1:2),"23"
    write(87,fmt=*) node_xy(25,1:2),"25"
    write(87,fmt=*) " "
    
    write(87,fmt=*) tol_Rdfn,"tolRdfn"
    
    if (lgcl_rdfn .eqv. .True.) then
       
       do
          do i=1,2
             l0(CP1(i)) = (ratio-ratioChng)*l0(CP1(i))
          enddo
          
          call Equilibrate_system
          call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
          Frame=Frame+1
          call switchto_NI_model_run_and_switchbackto_TN
          
          call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
          
          if (lgcl_rdfn.eqv..False.) then
             ratio = ratio+ratioChng
             write(87,fmt=*) ratio,"bfr exit"
             write(*,*) "stop in flnm: multiple_cycle,sb: adjust_critical_spr_l0_and_strctl0 "
             stop
             exit
          endif
       enddo
       
       do i = 1,2
          l0(CP1(i)) = ratio*initial_l0_ofCP
       enddo
       
    elseif (lgcl_rdfn .eqv. .False.) then
       do
          ratio = ratio+ratioChng
          
          do i=1,2
             l0(CP1(i)) = ratio*initial_l0_ofCP
          enddo
          
          call Equilibrate_system
          call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
          Frame=Frame+1
          call switchto_NI_model_run_and_switchbackto_TN
          
          call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
          write(87,fmt=*) lgcl_rdfn,"lgclF"
          write(87,fmt=*) tol_Rdfn,"tol_Rdfn"
          write(87,fmt=*) node_xy(11,1:2),"11"
          write(87,fmt=*) node_xy(23,1:2),"23"
          write(87,fmt=*) node_xy(25,1:2),"25"
          write(87,fmt=*) " "
          
          if (lgcl_rdfn.eqv..True.) then
             exit
          endif
       enddo
       
    endif
    
    write(87,fmt=*) ratio,"ratio aft"
    write(87,fmt=*) CP1(1:2),"CP1"
    close(87)
    
    l0_strctTN(3,CP1(1)) = l0(CP1(1))
    l0_strctTN(3,CP1(2)) = l0(CP1(2))
    
  end subroutine adjust_critical_spr_l0_and_strctl0
  
  subroutine Extrapolate_and_adjust_strctProp_ifNeeded
    implicit none
    real*8  :: tol_Rdfn
    logical :: lgcl_rdfn
    real*8  :: incrSize
    
    open(unit=89,file='Extrapolate.dat',position='append')
    
    tol_Rdfn = 0.05d0
    lgcl_rdfn = .False.
    
    call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)

    if (lgcl_rdfn .eqv. .True.) then
       continue
       !write for intrapolate (reverse of extrapolate)
       write(89,fmt=*) lgcl_rdfn,"lgcl_rdfn"
    elseif (lgcl_rdfn .eqv. .False.) then
       write(89,fmt=*) lgcl_rdfn,"lgcl_rdfn"
       
       incrSize = 0.025d0
       
       call allocate_and_initialize_StepVars
       call get_strct_step_in_between_strcts(incrSize)
       
       do
          call Extrapolate_A0l0kaks_for_gettingOverPulley
          
          if (Frame==4) then
             open(unit=87,file='PropChkinAlgChng.dat',position='append')
             
             do i = 1,N_spr
                write(87,*) k_spr(i),l0(i),i
             enddo
             write(87,*) " "
             do i = 1,N_cell
                write(87,*) k_area(i),A0(i),i
             enddo
             write(87,*) " "
             write(87,*) " "
             close(87)
          endif
          
          call Equilibrate_system
          call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
          Frame = Frame+1
          call switchto_NI_model_run_and_switchbackto_TN
          
          call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
          write(89,fmt=*) lgcl_rdfn,(Frame-1),"lgcl_rdfn,Frame"
          
          if (lgcl_rdfn.eqv..True.) then
             call adjst_strctF_Props
             write(*,*) "stop in flnm: multiple_cycle,sb: Extrapolate_and_adjust_strctProp_ifNeeded"
             !stop
             exit
          endif
          
       enddo
       
    endif
    
    close(89)
    
  end subroutine Extrapolate_and_adjust_strctProp_ifNeeded
  
  !subroutine allocate_and_initialize_stepVars
   ! implicit none
    
    !allocate(A0_stepF(1:N_cell),ka_stepF(1:N_cell))
    !allocate(l0_stepF(1:N_spr),ks_stepF(1:N_spr))
    
    !A0_stepF = -1.0d30 ; ka_stepF = -1.0d30
    !l0_stepF = -1.0d30 ; ks_stepF = -1.0d30
    
  !end subroutine allocate_and_initialize_stepVars
  
  !subroutine get_strct_step_from_strc1_to_strc3(incrSize)
   ! implicit none
    !real*8, intent(in) :: incrSize !can also be done for indiviv param incrSize
    !integer :: i,j,jmax
    
    !do i = 1,2
       
       !if (i==1) jmax = N_cell
       !if (i==2) jmax = N_spr
       
       !do j = 1,jmax
        !  if (i==1) then
         !    A0_stepF(j) = incrSize*(A0_strct(3,j)-A0_strct(1,j))  
         !    ka_stepF(j) = incrSize*(ka_strct(3,j)-ka_strct(1,j))
         ! elseif (i==2) then
          !   l0_stepF(j) = incrSize*(l0_strct(3,j)-l0_strct(1,j))  
           !  ks_stepF(j) = incrSize*(ks_strct(3,j)-ks_strct(1,j))
          !endif
       !enddo
       
    !enddo
    
  !end subroutine get_strct_step_from_strc1_to_strc3
  
  ! subroutine Extrapolate_A0l0kaks_for_gettingOverPulley
  !   implicit none
  !   integer :: i,j,jmax
  !   real*8  :: Tiny
    
  !   Tiny = -0.0d0
    
  !   do i = 1,2
  !      if (i==1) jmax=N_cell
  !      if (i==2) jmax=N_spr
       
  !      do j = 1,jmax
  !         if (i==1) then
  !            A0(j)     = A0(j) + A0_stepF(j)
  !            k_area(j) = k_area(j) + ka_stepF(j)
             
  !            if ((A0(j).lt.Tiny).or.(k_area(j).lt.Tiny)) then
  !               write(*,*) "A0/k_area got_negative"
  !               write(*,*) A0(j),k_area(j),j,"A0/k_area/j"
  !               stop
  !            endif
             
  !         elseif (i==2) then
  !            l0(j)     = l0(j) + l0_stepF(j)
  !            k_spr(j)  = k_spr(j) + ks_stepF(j)

  !            if ((l0(j).lt.Tiny).or.(k_spr(j).lt.Tiny)) then
  !               write(*,*) "l0/k_spr got_negative"
  !               write(*,*) l0(j),k_spr(j),j,"l0/k_spr/j"
  !               stop
  !            endif
             
  !         endif
  !      enddo
       
  !   enddo
    
  ! end subroutine Extrapolate_A0l0kaks_for_gettingOverPulley
  
  !subroutine adjst_strct3_Props
   ! implicit none
    !integer :: i,j,jmax
    
    !do i = 1,2
     !  if (i==1) jmax=N_cell
      ! if (i==2) jmax=N_spr
       
       !do j = 1,jmax
         ! if (i==1) then
        !     A0_strct(3,j) = A0(j)
          !   ka_strct(3,j) = k_area(j)
          !elseif (i==2) then
           !  l0_strct(3,j) = l0(j)
           !  ks_strct(3,j) = k_spr(j)
          !endif
       !enddo
       
    !enddo
    
  !end subroutine adjst_strct3_Props
  
  
end module mltpl_cycle

