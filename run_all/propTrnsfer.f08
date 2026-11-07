module PropTrnsferToNxtCell
  use calls_for_tests
  use redefining_system_module
  
  implicit none
  real*8, allocatable  :: ks_PT(:,:),l0_PT(:,:)
  real*8, allocatable  :: ka_PT(:,:),A0_PT(:,:)
  real*8, allocatable  :: CgX_PT(:,:)
  
  integer, allocatable :: cellPT(:),sprPT(:),nodePT(:)

  real*8, allocatable  :: ks_Intr(:) ,l0_Intr(:)
  real*8, allocatable  :: ka_Intr(:),A0_Intr(:)
  integer,allocatable  :: smthSpr(:),smthCell(:)
  
contains
  
  
  subroutine ProgStgPropTrnsfr
    implicit none
    integer :: N_sprPT,N_cellPT,N_nodePT
    integer :: i,j,jmax
    integer :: N_itm,N_strctPT
    
    character(len=100) :: fileInit,fileIntr1,fileIntr2,fileFinl
    character(len=100) :: flnm,flnmInit,flnmIntr1,flnmIntr2,flnmFinl
    
    integer :: fileInitNo,fileIntr1No,fileIntr2No,fileFinlNo
    
    !This routine is for S1T1 and to make more progression cycle
    
    flnm='strctPropsAddedCell'
    
    N_strctPT   = 4
    
    allocate(ks_PT(1:N_strctPT,1:N_spr) ,l0_PT(1:N_strctPT,1:N_spr))
    allocate(ka_PT(1:N_strctPT,1:N_cell),A0_PT(1:N_strctPT,1:N_cell))
    allocate(CgX_PT(1:N_strctPT,1:N_node))
    
    allocate(cellPT(1:N_cell),sprPT(1:N_spr),nodePT(1:N_node))
    
    allocate(ks_Intr(1:N_spr),l0_Intr(1:N_spr))
    allocate(ka_Intr(1:N_cell),A0_Intr(1:N_cell))
    allocate(smthSpr(1:N_spr),smthCell(1:N_cell))
    
    fileInitNo  = 13
    fileIntr1No = 14
    fileIntr2No = 15
    fileFinlNo  = 16
    
    write(fileInit,'(i2.2,a)')  fileInitNo,'S1T1.dat'
    write(fileIntr1,'(i2.2,a)') fileIntr1No,'S1T1.dat'
    write(fileIntr2,'(i2.2,a)') fileIntr2No,'S1T1.dat'
    write(fileFinl,'(i2.2,a)')  fileFinlNo,'S1T1.dat'
    
    flnmInit  =trim(adjustl(flnm))//trim(adjustl(fileInit)) 
    flnmIntr1 =trim(adjustl(flnm))//trim(adjustl(fileIntr1))
    flnmIntr2 =trim(adjustl(flnm))//trim(adjustl(fileIntr2))
    flnmFinl  =trim(adjustl(flnm))//trim(adjustl(fileFinl))
    
    write(*,*) trim(adjustl(flnm))//trim(adjustl(fileInit))
    write(*,*) trim(adjustl(flnm))//trim(adjustl(fileIntr1))
    write(*,*) trim(adjustl(flnm))//trim(adjustl(fileIntr2))
    write(*,*) trim(adjustl(flnm))//trim(adjustl(fileFinl))
    
    call sleep(2)
    
    open(unit=51,file=trim(adjustl(flnmInit)))
    open(unit=52,file=trim(adjustl(flnmIntr1)))
    open(unit=53,file=trim(adjustl(flnmIntr2)))
    open(unit=54,file=trim(adjustl(flnmFinl)))
    
    N_itm = 3
    
    do i=1,N_itm
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       if (i==3) jmax=N_node
       
       do j = 1,jmax
          
          if (i==1) then
             
             read(51,*) ks_PT(1,j),l0_PT(1,j) 
             read(52,*) ks_PT(2,j),l0_PT(2,j) 
             read(53,*) ks_PT(3,j),l0_PT(3,j) 
             read(54,*) ks_PT(4,j),l0_PT(4,j) 
             
          elseif (i==2) then
             
             read(51,*) ka_PT(1,j),A0_PT(1,j)
             read(52,*) ka_PT(2,j),A0_PT(2,j) 
             read(53,*) ka_PT(3,j),A0_PT(3,j)
             read(54,*) ka_PT(4,j),A0_PT(4,j)
             
          elseif (i==3) then
             
             read(51,*) CgX_PT(1,j)
             read(52,*) CgX_PT(2,j)
             read(53,*) CgX_PT(3,j)
             read(54,*) CgX_PT(4,j)     
             
          endif
          
       enddo
       
    enddo
    
    close(51)
    close(52)
    close(53)
    close(54)
    
    
    open(unit=96,file='PT_values1.dat')
    open(unit=97,file='PT_values2.dat')
    open(unit=98,file='PT_values3.dat')
    open(unit=99,file='PT_values4.dat')
    
    do i=1,N_itm
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       if (i==3) jmax=N_node
       
       do j = 1,jmax
          
          if (i==1) then
             
             write(96,*) ks_PT(1,j),l0_PT(1,j) 
             write(97,*) ks_PT(2,j),l0_PT(2,j) 
             write(98,*) ks_PT(3,j),l0_PT(3,j) 
             write(99,*) ks_PT(4,j),l0_PT(4,j) 
             
          elseif (i==2) then
             
             write(96,*) ka_PT(1,j),A0_PT(1,j)
             write(97,*) ka_PT(2,j),A0_PT(2,j) 
             write(98,*) ka_PT(3,j),A0_PT(3,j)
             write(99,*) ka_PT(4,j),A0_PT(4,j)
             
          elseif (i==3) then
             
             write(96,*) CgX_PT(1,j)
             write(97,*) CgX_PT(2,j)
             write(98,*) CgX_PT(3,j)
             write(99,*) CgX_PT(4,j)     
             
          endif
          
       enddo
       
       write(96,*) " "
       write(97,*) " "
       write(98,*) " "
       write(99,*) " "  
       
    enddo
    
    close(96)
    close(97)
    close(98)
    close(99)
    
  end subroutine ProgStgPropTrnsfr
  
  
  subroutine getPT_arrays(ProgCycl,FirstPTCell,BotmCellPExclusn,cellPT,sprPT,nodePT)
    implicit none
    integer, intent(in)  :: ProgCycl,FirstPTcell,BotmCellPExclusn
    integer, intent(out) :: cellPT(1:N_cell),sprPT(1:N_spr),nodePT(1:N_node)
    
    integer :: i,j,jmax,k,N_itm
    integer :: lim1,lim2,lim3,lim4,lim5,lim6
    
    integer :: FirstPTlftCell,FirstPTrghtCell
    integer :: cellPTInpt,cellPTOtpt
    integer :: sprPTInpt(1:3),sprPTOtpt(1:3)
    
    integer :: sprIn,sprOut
    integer :: cntBtmL,cntBtmR
    
    FirstPTlftCell  = FirstPTCell
    Hlf_Ncell       = N_cell/2
    FirstPTrghtCell = Hlf_Ncell + FirstPTCell
    
    lim1 = FirstPTlftCell-1
    lim2 = Hlf_Ncell-1-BotmCellPExclusn
    lim3 = Hlf_Ncell
    lim4 = Hlf_Ncell+FirstPTlftCell-1
    
    if(mod(N_cell,2)==1) then
       lim5 = (N_cell-1)-1-BotmCellPExclusn
       lim6 = N_cell-1
    elseif(mod(N_cell,2)==0) then
       lim5 = N_cell-1+BotmCellPExclusn
       lim6 = N_cell
    endif

    cntBtmL=1 ; cntBtmR=1
    
    do j = 1,N_cell
       
       if (j.le.lim1) then
          cellPT(j) = j
          
          cellPTInpt = j ; cellPTOtpt = cellPT(j)
          call cellPT_to_sprPT(cellPTInpt,cellPTOtpt,sprPTInpt,sprPTOtpt)
          
          do k=1,nsecl
             sprIn = sprPTInpt(k) ; sprOut = sprPTOtpt(k) 
             sprPT(sprIn) = sprOut
          enddo
          
       elseif ((j.gt.lim1).and.(j.le.lim2)) then
          cellPT(j) = j+ProgCycl
          
          cellPTInpt = j ; cellPTOtpt = cellPT(j)
          call cellPT_to_sprPT(cellPTInpt,cellPTOtpt,sprPTInpt,sprPTOtpt)
          
          do k=1,nsecl
             sprIn = sprPTInpt(k) ; sprOut = sprPTOtpt(k) 
             sprPT(sprIn) = sprOut
          enddo
          
       elseif ((j.gt.lim2) .and. (j.le.lim3)) then
          
          if (ProgCycl==1) then
             cellPT(j) = j
             
          elseif (ProgCycl==2) then
             
             if (cntBtmL==1)   cellPT(j) = j + (ProgCycl-1)
             if (cntBtmL.gt.1) cellPT(j) = j
             
             cntBtmL = cntBtmL+1
             
          elseif (ProgCycl==3) then
             
             if (cntBtmL==1)   cellPT(j) = j + (ProgCycl-1)
             if (cntBtmL==2)   cellPT(j) = j + (ProgCycl-2)
             if (cntBtmL.gt.2) cellPT(j) = j
             
             cntBtmL = cntBtmL+1
             
          endif
             
          cellPTInpt = j ; cellPTOtpt = cellPT(j)
          call cellPT_to_sprPT(cellPTInpt,cellPTOtpt,sprPTInpt,sprPTOtpt)
          
          do k=1,nsecl
             sprIn = sprPTInpt(k) ; sprOut = sprPTOtpt(k) 
             sprPT(sprIn) = sprOut
          enddo
          
            
       elseif ((j.gt.lim3) .and. (j.le.lim4)) then
          cellPT(j) = j
          
          cellPTInpt = j ; cellPTOtpt = cellPT(j)
          call cellPT_to_sprPT(cellPTInpt,cellPTOtpt,sprPTInpt,sprPTOtpt)
          
          do k=1,nsecr
             sprIn = sprPTInpt(k) ; sprOut = sprPTOtpt(k) 
             sprPT(sprIn) = sprOut
          enddo
          
       elseif ((j.gt.lim4) .and. (j.le.lim5)) then
          cellPT(j) = j+ProgCycl
          
          cellPTInpt = j ; cellPTOtpt = cellPT(j)
          call cellPT_to_sprPT(cellPTInpt,cellPTOtpt,sprPTInpt,sprPTOtpt)
            
          do k=1,nsecr
             sprIn = sprPTInpt(k) ; sprOut = sprPTOtpt(k) 
             sprPT(sprIn) = sprOut
          enddo
          
       elseif ((j.gt.lim5) .and. (j.le.lim6)) then
          
          if (ProgCycl==1) then
             cellPT(j) = j
             
          elseif (ProgCycl==2) then
             
             if (cntBtmR==1)   cellPT(j) = j + (ProgCycl-1)
             if (cntBtmR.gt.1) cellPT(j) = j
             
             cntBtmR = cntBtmR+1
             
          elseif (ProgCycl==3) then
             
             if (cntBtmR==1)   cellPT(j) = j + (ProgCycl-1)
             if (cntBtmR==2)   cellPT(j) = j + (ProgCycl-2)
             if (cntBtmR.gt.2) cellPT(j) = j
             
             cntBtmR = cntBtmR+1
             
          endif
          
          
          cellPTInpt = j ; cellPTOtpt = cellPT(j)
          call cellPT_to_sprPT(cellPTInpt,cellPTOtpt,sprPTInpt,sprPTOtpt)
          
          do k=1,nsecr
             sprIn = sprPTInpt(k) ; sprOut = sprPTOtpt(k) 
             sprPT(sprIn) = sprOut
          enddo
          
          
       elseif (j.gt.lim6) then
          cellPT(j) = j
          
          sprIn = 3*j-2 ; sprOut = 3*j-2
          sprPT(sprIn) = sprOut
          
       endif
       
    enddo
    
    
    do j=1,N_node
       
       if ((j==(2*FirstPTlftCell-1)) .or. (j==2*FirstPTlftCell)) then !(37/2)
          nodePT(j) = 1
       elseif ((j==((N_node/2)+(2*FirstPTlftCell-1))) .or. (j==((N_node/2)+(2*FirstPTlftCell)))) then
          nodePT(j) = -1
       else
          nodePT(j) = 0
       endif
       
    enddo
    
    open(unit=76,file='PT_arrays.dat',position='append')
    
    N_itm = 3
    
    do i = 1,N_itm
       
       if (i==1) jmax = N_cell
       if (i==2) jmax = N_spr
       if (i==3) jmax = N_node
       
       do j = 1,jmax
          if (i==1) write(76,*) "j=",j,"and cellPT=",cellPT(j)
          if (i==2) write(76,*) "j=",j,"and sprPT=",sprPT(j)
          if (i==3) write(76,*) "j=",j,"and nodePT=",nodePT(j)
       enddo
       
       write(76,*) " "
       write(76,*) lim1,lim2,lim3,lim4,lim5,lim6,"lims for cellPT"
       
    enddo
    
    close(76)
    
  end subroutine getPT_arrays
  
  
  subroutine cellPT_to_sprPT(cellPTInpt,cellPTOtpt,sprPTInpt,sprPTOtpt)
    implicit none
    integer, intent(in)  :: cellPTInpt,cellPTOtpt 
    integer, intent(out) :: sprPTInpt(1:3),sprPTOtpt(1:3)
    
    sprPTInpt(1) = 3*cellPTInpt-2
    sprPTInpt(2) = 3*cellPTInpt-1
    sprPTInpt(3) = 3*cellPTInpt-0
    
    sprPTOtpt(1) = 3*cellPTOtpt-2
    sprPTOtpt(2) = 3*cellPTOtpt-1
    sprPTOtpt(3) = 3*cellPTOtpt-0
    
  end subroutine cellPT_to_sprPT
  
  subroutine get_CgXVal(CgXVal)
    implicit none
    real*8, intent(out) :: CgXVal
    
    real*8  :: ZERO
    integer :: i
    
    ZERO = 0.00d0
    
    do i = 1,N_node
       if (abs(CgXNode(i)) .gt. ZERO) then
          CgXVal = abs(CgXNode(i))
       else
          CgXVal = ZERO
       endif
    enddo
    
  end subroutine get_CgXVal
  
  
  subroutine get_PT_InptOtpt(Prtn,ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F)
    implicit none
    integer, intent(in)  :: Prtn
    real*8,  intent(out) :: ks_I(1:N_spr),ks_F(1:N_spr),l0_I(1:N_spr),l0_F(1:N_spr)
    real*8,  intent(out) :: ka_I(1:N_cell),ka_F(1:N_cell),A0_I(1:N_cell),A0_F(1:N_cell)
    
    integer :: N_itm
    integer :: i,j,jmax
    integer :: sprPTNum,cellPTNum,nodePTNum
    real*8  :: CgXVal
    
    call get_CgXVal(CgXVal)
    
    open(unit=57,file='PT_InptOtpt.dat')
    
    N_itm = 3
    
    do i = 1,N_itm
       
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
       if (i==3) jmax = N_node
       
       do j = 1,jmax
          
          if (i==1) then
             sprPTNum = sprPT(j)
             
             if (prtn==1) then
                ks_I(j) = k_spr(j)
                ks_F(j) = ks_PT(2,sprPTnum)
                l0_I(j) = l0(j)
                l0_F(j) = l0_PT(2,sprPTnum)
                
                write(57,*) ks_I(j),ks_F(j),l0_I(j),l0_F(j),j,"ksI,ksF,l0I,l0F for prtn1"
                
             elseif (prtn==2) then
                ks_I(j) = k_spr(j)
                ks_F(j) = ks_PT(4,sprPTnum)
                l0_I(j) = l0(j)
                l0_F(j) = l0_PT(4,sprPTnum)
                
                write(57,*) ks_I(j),ks_F(j),l0_I(j),l0_F(j),j,"ksI,ksF,l0I,l0F for prtn2"
                
             endif
                
          elseif (i==2) then
             cellPTNum = cellPT(j)
             
             if (prtn==1) then
                ka_I(j)  = k_area(j)
                ka_F(j)  = ka_PT(2,cellPTnum)
                A0_I(j)  = A0(j)
                A0_F(j)  = A0_PT(2,cellPTnum)
                
                write(57,*) ka_I(j),ka_F(j),A0_I(j),A0_F(j),j,"kaI,kaF,A0I,A0F for prtn1"
                
             elseif (prtn==2) then
                ka_I(j)  = k_area(j)
                ka_F(j)  = ka_PT(4,cellPTnum)
                A0_I(j)  = A0(j)
                A0_F(j)  = A0_PT(4,cellPTnum)
                
                write(57,*) ka_I(j),ka_F(j),A0_I(j),A0_F(j),j,"kaI,kaF,A0I,A0F for prtn2"
                
             endif
             
          elseif (i==3) then
             nodePTNum = nodePT(j)
             
             if (nodePTNum==1) then
                CgXNode(j) = CgXVal
             elseif (nodePTNum==(-1)) then
                CgXNode(j) = -CgXVal
             elseif (nodePTNum==0) then
                CgXNode(j) = 0.00d0
             endif

             write(57,*) CgXNode(j),j,"CgXNode,j"
             
          endif
          
       enddo
       
       write(57,*) " "
       
    enddo
    
    
    close(57)
    
  end subroutine get_PT_InptOtpt
  
  
  subroutine get_PT_InptOtpt_WO_Intr(Prtn,ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F)
    implicit none
    integer, intent(in)  :: Prtn
    real*8,  intent(out) :: ks_I(1:N_spr),ks_F(1:N_spr),l0_I(1:N_spr),l0_F(1:N_spr)
    real*8,  intent(out) :: ka_I(1:N_cell),ka_F(1:N_cell),A0_I(1:N_cell),A0_F(1:N_cell)
    
    integer :: N_itm
    integer :: i,j,jmax
    integer :: sprPTnum,cellPTnum,nodePTnum
    real*8  :: CgXVal
    
    call get_CgXVal(CgXVal)
    
    if (prtn==1) then
       call get_smthdIntrProp
    endif
    
    open(unit=57,file='PT_InptOtpt_WoIntr.dat')
    
    N_itm = 3
    
    do i = 1,N_itm
       
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
       if (i==3) jmax = N_node
       
       do j = 1,jmax
          
          if (i==1) then
             sprPTnum = sprPT(j)
             
             if (prtn==1) then
                ks_I(j) = k_spr(j)
                l0_I(j) = l0(j)
                
                if (smthSpr(j)==0) then
                   ks_F(j) = ks_PT(2,sprPTnum)
                   l0_F(j) = l0_PT(2,sprPTnum)
                   
                elseif (smthSpr(j)==1) then
                   ks_F(j) = ks_Intr(j)
                   l0_F(j) = l0_Intr(j)
                   
                endif
                   
                write(57,*) ks_I(j),ks_F(j),l0_I(j),l0_F(j),j,"ksI,ksF,l0I,l0F for prtn1"
                
             elseif (prtn==2) then
                ks_I(j) = k_spr(j)
                ks_F(j) = ks_PT(4,sprPTnum)
                l0_I(j) = l0(j)
                l0_F(j) = l0_PT(4,sprPTnum)
                
                write(57,*) ks_I(j),ks_F(j),l0_I(j),l0_F(j),j,"ksI,ksF,l0I,l0F for prtn2"
                
             endif
                
          elseif (i==2) then
             cellPTnum = cellPT(j)
             
             if (prtn==1) then
                ka_I(j)  = k_area(j)
                A0_I(j)  = A0(j)
                
                if (smthCell(j)==0) then
                   ka_F(j)  = ka_PT(2,cellPTnum)
                   A0_F(j)  = A0_PT(2,cellPTnum)
                   
                elseif (smthCell(j)==1) then
                   ka_F(j)  = ka_Intr(j)
                   A0_F(j)  = A0_Intr(j)
                endif
                   
                write(57,*) ka_I(j),ka_F(j),A0_I(j),A0_F(j),j,"kaI,kaF,A0I,A0F for prtn1"
                
             elseif (prtn==2) then
                ka_I(j)  = k_area(j)
                ka_F(j)  = ka_PT(4,cellPTnum)
                A0_I(j)  = A0(j)
                A0_F(j)  = A0_PT(4,cellPTnum)
                
                write(57,*) ka_I(j),ka_F(j),A0_I(j),A0_F(j),j,"kaI,kaF,A0I,A0F for prtn2"
                
             endif
             
          elseif (i==3) then
             nodePTnum = nodePT(j)
             
             if (nodePTnum==1) then
                CgXNode(j) = CgXVal
             elseif (nodePTnum==(-1)) then
                CgXNode(j) = -CgXVal
             elseif (nodePTnum==0) then
                CgXNode(j) = 0.00d0
             endif
             
             write(57,*) CgXNode(j),j,"CgXNode,j"
             
          endif
          
       enddo
       
       write(57,*) " "
       
    enddo
    
    close(57)
    
    
  end subroutine get_PT_InptOtpt_WO_Intr
  
  subroutine get_smthdIntrProp
    implicit none
    integer :: i,j,jmax
    integer :: N_itm
    real*8  :: initVal,finlVal,intrmdVal
    integer :: sprPTnum,cellPTnum
    
    open(unit=75,file='InitFinlIntrPropVal.dat')
    
    N_itm = 2
    
    do i = 1,N_itm
       
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
       
       do j = 1,jmax
          
          if (i==1) then
             
             sprPTnum = sprPT(j)
             
             initVal = k_spr(j)
             finlVal = ks_PT(4,sprPTnum) 
             
             call get_Intrmd(initVal,finlVal,intrmdVal)
             write(75,*) initVal,finlVal,intrmdVal,j,"ks_val"
             
             ks_Intr(j) = intrmdVal 
             
             initVal = l0(j)
             finlVal = l0_PT(4,sprPTnum)
             
             call get_Intrmd(initVal,finlVal,intrmdVal)
             write(75,*) initVal,finlVal,intrmdVal,j,"l0_val"
             
             l0_Intr(j) = intrmdVal
             
          elseif (i==2) then
             
             cellPTnum = cellPT(j)
             
             initVal = k_area(j)
             finlVal = ka_PT(4,cellPTnum)
             
             call get_Intrmd(initVal,finlVal,intrmdVal)
             write(75,*) initVal,finlVal,intrmdVal,j,"ka_val"
             
             ka_Intr(j) = intrmdVal
             
             initVal = A0(j)
             finlVal = A0_PT(4,cellPTnum)
             
             call get_Intrmd(initVal,finlVal,intrmdVal)
             write(75,*) initVal,finlVal,intrmdVal,j,"A0_val"
             
             A0_Intr(j) = intrmdVal
             
          endif
          
       enddo
       write(75,*) " "
       
    enddo
    
    close(75)
    
    open(unit=74,file='IntrPropVal.dat')
    
    do i = 1,2
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       
       do j=1,jmax
          if (i==1) write(74,*) ks_Intr(j),l0_Intr(j),j,"ks_Intr,l0_Intr,j"
          if (i==2) write(74,*) ka_Intr(j),A0_Intr(j),j,"ka_Intr,A0_Intr,j"
       enddo
       
    enddo
    
    close(74)
    
    
  end subroutine get_smthdIntrProp
  
  
  subroutine get_SmthCellSprsArray(cellNo)
    implicit none
    integer, intent(in) :: cellNo
    integer :: i,j,jmax
    integer :: N_itm
    integer :: spr1,spr2,spr3
    !integer :: cellL,cellR
    integer :: cellL(1:2),cellR(1:2)
    
    N_itm = 2
    
    cellL(1) = cellNo
    cellR(1) = Hlf_Ncell+cellNo
    
    cellL(2) = (cellNo) - 1
    cellR(2) = (Hlf_Ncell+cellNo) - 1
    
    !cellL(3) = (cellNo) + 1
    !cellR(3) = (Hlf_Ncell+cellNo) + 1
    
    open(unit=84,file='SmthCellSprsArray.dat',position='append')
    write(84,*) Hlf_Ncell,cellNo,"Hlf_Ncell,Approaching Cell"
    close(84)
    
    smthCell(1:N_cell) = 1
    smthSpr(1:N_spr)   = 1
    
    do i = 1,N_cell
       
       do j = 1,1
          
          if ((i==cellL(j)).or.(i==cellR(j))) then
             
             smthCell(i)=0
             
             spr1 = 3*i-2
             spr2 = 3*i-1
             spr3 = 3*i-0
             
             smthSpr(spr1)=0
             smthSpr(spr2)=0
             smthSpr(spr3)=0
             
          endif
          
       enddo
       
    enddo
    
    open(unit=73,file='SmthArray.dat')
    
    do i = 1,2
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       
       do j=1,jmax
          if (i==1) write(73,*) smthSpr(j),j,"smthSpr"
          if (i==2) write(73,*) smthCell(j),j,"smthCell"
       enddo
       
    enddo
    
    close(73)
    
  end subroutine get_SmthCellSprsArray
  
  
  
  subroutine get_Intrmd(initVal,finlVal,intrmdVal)
    implicit none
    real*8, intent(in)  :: initVal,finlVal
    real*8, intent(out) :: intrmdVal
    
    real*8 :: N_stp,Hlf_Nstp
    real*8 :: diff,stpSize
    real*8 :: TINY_val
    
    TINY_val = 1.0d-16
    
    N_stp    = 20.0
    Hlf_Nstp = 10.0
    
    diff    = finlVal-initVal
    stpSize = diff/N_stp
    
    if (abs(diff).lt.TINY_val) then
       intrmdVal = initVal
       
    elseif (abs(diff).ge.TINY_val) then
       intrmdVal = initVal + (Hlf_Nstp)*stpSize
    endif
    
  end subroutine get_Intrmd
  
  
end module PropTrnsferToNxtCell
