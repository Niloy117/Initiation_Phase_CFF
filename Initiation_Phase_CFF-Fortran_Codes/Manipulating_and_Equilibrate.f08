module Manipulate_and_equilibrate
  use MAE_initaitn
  use Control_State_Length_Match
  
  implicit none
  
contains
  
  subroutine changeA0_inSteps_andEquilibrate_using_f(Exprmnt,Frame)
    implicit none
    integer, intent(inout) :: Exprmnt,Frame
    real*8  :: InitA0(1:N_cell),FinlA0(1:N_cell)
    integer :: Nstep
    real*8  :: UnitStep,stepSize
    integer :: i
    
    write(*,*) Exprmnt,Frame,"Exp and Frm in changeA0_inSteps routine"
    InitA0=0.0d0 ; FinlA0=0.0d0
    
    call get_InitA0(InitA0)
    call get_FinlA0(FinlA0)
    
    write(*,*) InitA0(1),InitA0(11),InitA0(13),InitA0(16)
    write(*,*) FInlA0(1),FInlA0(11),FInlA0(13),FInlA0(16)
    
    Nstep = 11
    Unitstep = (1.0d0)/(Nstep-1)
    write(*,*) Unitstep,"UnitStep"
    
    do i = 1,Nstep
       
       stepSize = (i-1)*UnitStep
       call stepUpA0(InitA0,FinlA0,stepSize)
       write(*,*) A0(1),A0(11),A0(13),A0(16),"A0s aft stepping up"
       
       !if (strctNo==1.and.i==2) stop
       
       call nodesl0A0_to_coordntes(node_xy,l0,A0,coordntes)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
    enddo
    
    write(*,*) Frame,"Frame at last"
    
  end subroutine changeA0_inSteps_andEquilibrate_using_f
  
  
  subroutine changeA0_inSteps_Equilibrate_and_Optimize(Exprmnt,Frame)
    implicit none
    integer, intent(inout) :: Exprmnt,Frame
    real*8  :: InitA0(1:N_cell),FinlA0(1:N_cell)
    integer :: Nstep
    real*8  :: UnitStep,stepSize
    integer :: i
    
    write(*,*) Exprmnt,Frame,"Exp and Frm in changeA0_inSteps routine"
    InitA0=0.0d0 ; FinlA0=0.0d0
    
    call get_InitA0(InitA0)
    call get_FinlA0(FinlA0)
    
    write(*,*) InitA0(1),InitA0(11),InitA0(13),InitA0(16)
    write(*,*) FInlA0(1),FInlA0(11),FInlA0(13),FInlA0(16)
    
    Nstep = 11
    Unitstep = (1.0d0)/(Nstep-1)
    write(*,*) Unitstep,"UnitStep"
    
    do i = 1,Nstep
       stepSize = (i-1)*UnitStep
       call stepUpA0(InitA0,FinlA0,stepSize)
       write(*,*) A0(1),A0(11),A0(13),A0(16),"A0s aft stepping up"
       call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
       
       call Equilibrate_system
       call coordntes_to_nodes(coordntes_xy,node_xy)
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
    enddo
    
    if (stageNo==1 .and. stageType==1) then
       continue
    elseif (stageNo==4) then
       call Pressure_adjstmnt_decision(node_xy,A0)
    endif
    
    do i = 1,N_cell
       !k_area(i) = (k_area(i))/2.0d0
       !write(*,*) k_area(i),"karea"
    enddo
    
    write(*,*) Frame,"Frame bfr adjust"
    call adjust_spr_length_and_Equilibrate(Exprmnt,Frame)
    write(*,*) Frame,"Frame at last"
    
  end subroutine changeA0_inSteps_Equilibrate_and_Optimize
  
  
  subroutine get_InitA0(InitA0)
    implicit none
    real*8, intent(inout) :: InitA0(1:N_cell)
    integer :: j

    do j = 1,N_cell
       InitA0(j) = A0(j)
    enddo
    
  end subroutine get_InitA0
  
  subroutine get_FinlA0(FinlA0)
    implicit none
    real*8, intent(inout) :: FinlA0(1:N_cell)
    integer :: j
    
    do j = 1,N_cell
       FinlA0(j) = fctr(j)*A0(j)
    enddo
    
  end subroutine get_FinlA0
  
  subroutine stepUpA0(InitA0,FinlA0,stepSize)
    implicit none
    real*8, intent(in) :: InitA0(1:N_cell)
    real*8, intent(in) :: FinlA0(1:N_cell)
    real*8, intent(in) :: stepSize
    
    integer :: j
    
    do j = 1,N_cell
       A0(j) = InitA0(j) + stepSize*(FinlA0(j)-InitA0(j))
    enddo
    
  end subroutine stepUpA0
  
  
  subroutine adjust_spr_length_and_Equilibrate(Exprmnt,Frame)
    implicit none
    integer, intent(inout) :: Exprmnt,Frame
    
    real*8  :: f_Lt(1:N_spr),f_Tilt(1:N_spr)
    integer :: m,mmax,j
    real*8  :: Enr

    Enr = 0.0d0
    mmax = 10
    open(unit=331,file="adjust_l0.dat",position="append")
    open(unit=332,file="ll0Lt_fLt_ftilt.dat",position="append")
    
    do m = 1,mmax
       
       call coordntes_to_nodes(coordntes_xy,node_xy)
       call calc_length(node_xy)
       write(331,fmt=*) l(43),l(46),"l bfr adjustmnt of spr 43,46"
       
       f_Lt = 10.0d5 ; f_Tilt = 10.0d5
       call get_f_Lt(f_Lt)
       write(331,fmt=*) f_Lt(43),f_Lt(46),"f_Lt of spr 43,46"
       call get_f_Tilt(f_Lt,f_Tilt)
       call get_l0s(f_Tilt)
       write(331,fmt=*) f_Tilt(43),f_Tilt(46),"f_Tilt of spr 43,46"
       
       call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
       call Equilibrate_system
       write(331,fmt=*) l(43),l0(43),Lt(43),m,"l,l0,Lt of spr 43 in Iteration No"
       write(331,fmt=*) l(46),l0(46),Lt(46),m,"l,l0,Lt of spr 46 in Iteration No"
       
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       write(331,fmt=*) " "
       
       if (m==10) then
          call coordntes_to_nodes(coordntes_xy,node_xy)
          call calc_length(node_xy)
          
          write(332,fmt=*) "l,l0,Lt,f_Lt,f_tilt"
          
          do j = 1,N_spr
             write(332,fmt=*) l(j),l0(j),Lt(j),f_Lt(j),f_tilt(j)
          enddo
          
          !do j = 1,N_spr
           !  Enr = Enr + k_spr(j)*(l(j)-l0(j))**2
            ! write(332,fmt=*) Enr,j,"Enr at spr_no"
          !enddo
          
          !do j = 1,N_cell
           !  Enr = Enr + k_area(j)*(A(j)-A0(j))**2
           !  write(332,fmt=*) Enr,j,"Enr at cell no"
          !enddo
       endif
       
    enddo
    
    close(331)
    close(332)
    
  end subroutine adjust_spr_length_and_Equilibrate
  
  
  subroutine change_thirdPairA0_inSteps_Equilibrate(Exprmnt,Frame)
    implicit none
    integer, intent(inout) :: Exprmnt,Frame
    real*8  :: InitA0(1:2),FinlA0(1:2)
    integer :: cellPair(1:2)
    real*8  :: prv_fctr(1:2),trgt_fctr(1:2)
    integer :: Nstep
    real*8  :: UnitStep,stepSize
    integer :: i
    
    write(*,*) Exprmnt,Frame,"Exp and Frm in changeA0_inSteps routine"
    
    cellPair(1) = ncl+ncr+3
    cellPair(2) = ncl+ncr+4
    write(*,*) cellPair(1:2)

    InitA0(1) = A0(cellPair(1))
    InitA0(2) = A0(cellPair(2))

    prv_fctr(1) = fctr(cellPair(1))
    prv_fctr(2) = fctr(cellPair(2))

    trgt_fctr(1) = 1.10d0
    trgt_fctr(2) = 1.10d0

    fctr(cellPair(1)) = trgt_fctr(1)/prv_fctr(1)
    fctr(cellPair(2)) = trgt_fctr(2)/prv_fctr(2)

    write(*,*) fctr(cellPair(1)),fctr(cellPair(2))
    
    FinlA0(1) = fctr(cellPair(1))*A0(cellPair(1))
    FinlA0(2) = fctr(cellPair(2))*A0(cellPair(2))
        
    Nstep = 5
    Unitstep = (1.0d0)/(Nstep-1)
    write(*,*) Unitstep,"UnitStep"
    
    do i = 1,Nstep
       
       stepSize = (i-1)*UnitStep
       A0(cellPair(1))=InitA0(1)+stepSize*(FinlA0(1)-InitA0(1))
       A0(cellPair(2))=InitA0(2)+stepSize*(FinlA0(2)-InitA0(2))

       !if (strctNo==1.and.i==2) stop

       call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
    enddo

    do i = 1,N_cell
       write(*,*) A(i),A0(i),fctr(i),i,"A,A0,fctr,cell_no"
    enddo
    
    
    write(*,*) Frame,"Frame at last"
    
    
  end subroutine change_thirdPairA0_inSteps_Equilibrate
  

  subroutine calc_length(dum_Nodes)
    implicit none
    real*8, intent(in) :: dum_Nodes(1:N_node,1:N_dmnsn)
    integer :: spr_nmbr
    integer :: Node1,Node2,Node3
    real*8  :: N1(1:N_dmnsn),N2(1:N_dmnsn),N3(1:N_dmnsn)
    integer :: i
    
    interface
       real*8 function Length(N1,N2,N_dmnsn)
         implicit none
         integer, intent(in) :: N_dmnsn
         real*8 , intent(in) :: N1(1:N_dmnsn),N2(1:N_dmnsn)
       end function Length
       
       real*8 function Length_2Prtn(N1,N2,N3,N_dmnsn)
         implicit none
         integer, intent(in) :: N_dmnsn
         real*8 , intent(in) :: N1(1:N_dmnsn),N2(1:N_dmnsn),N3(1:N_dmnsn)
       end function Length_2Prtn   
    end interface
    
    !write(*,*) N_dmnsn,"N_dmnsn"
    do spr_nmbr = 1,N_spr
       
       if (spr_node(spr_nmbr,0) == 2) then
          !write(*,*) spr_nmbr,"spr_nmbr"
          
          Node1 = spr_node(spr_nmbr,1)
          Node2 = spr_node(spr_nmbr,2)

          
          N1(1:N_dmnsn) = dum_Nodes(Node1,1:N_dmnsn)
          N2(1:N_dmnsn) = dum_Nodes(Node2,1:N_dmnsn)
          
          l(spr_nmbr) = Length(N1,N2,N_dmnsn)
           
       elseif (spr_node(spr_nmbr,0) == 3) then
          !write(*,*) spr_nmbr,"spr_nmbr"
          
          Node1 = spr_node(spr_nmbr,1)
          Node2 = spr_node(spr_nmbr,2)
          Node3 = spr_node(spr_nmbr,3)
          
          N1(1:N_dmnsn) = dum_Nodes(Node1,1:N_dmnsn)
          N2(1:N_dmnsn) = dum_Nodes(Node2,1:N_dmnsn)
          N3(1:N_dmnsn) = dum_Nodes(Node3,1:N_dmnsn)   
          
          l(spr_nmbr) = Length_2Prtn(N1,N2,N3,N_dmnsn)
       
       endif
       
    enddo

  end subroutine calc_length
  
  
  subroutine get_f_Lt(f_Lt)
    implicit none
    real*8, intent(inout) :: f_Lt(1:N_spr) 
    integer :: i
    
    do i = 1,N_spr
       f_Lt(i) = l(i)/Lt(i)
    enddo
    
  end subroutine get_f_Lt
  
  subroutine get_f_Tilt(f_Lt,f_Tilt)
    implicit none
    real*8, intent(in) :: f_Lt(1:N_spr)
    real*8, intent(inout) :: f_Tilt(1:N_spr)
    
    integer :: i
    real*8  :: bita
    
    bita = 0.20d0
    
    do i = 1,N_spr
       f_Tilt(i) = 1.0d0 + bita*(f_Lt(i)-1.0d0)
    enddo

  end subroutine get_f_Tilt
  
  subroutine get_l0s(f_Tilt)
    implicit none
    real*8, intent(in) :: f_Tilt(1:N_spr)
    integer :: i

    do i = 1,N_spr
       l0(i) = (1.0d0/f_Tilt(i))*l0(i)
    enddo
    
  end subroutine get_l0s

  
  subroutine stepwise_shortening_MidSpr_frmtop_and_Equilibrate(Exprmnt,Frame,serial,hmuch)
    implicit none
    integer, intent(in) :: Exprmnt
    integer, intent(inout) :: Frame
    integer, intent(in) :: serial
    real*8 , intent(in) :: hmuch

    integer :: Nstp
    real*8  :: stepS
    integer :: middle_spr(1:2)
    real*8  :: Initl0(1:2)
    integer :: i
    real*8  :: hm

    open(unit=134,file="MidSpr_Shortn.dat")
    
    Nstp  = 5
    stepS = (hmuch)/(Nstp-1)
    
    call get_the_midSpr_frmtop(serial,middle_spr)
    
    Initl0(1) = l0(middle_spr(1))
    Initl0(2) = l0(middle_spr(2))

    write(134,fmt=*) Initl0(1:2),"Initl0"
    
    do i = 1,Nstp
       hm = (i-1)*stepS
       call shortening_sprPair(middle_spr,hm)
       write(134,fmt=*) middle_spr(1),middle_spr(2),"middle sprs"
       write(134,fmt=*) l0(middle_spr(1)),l0(middle_spr(2)),"l0 aft shrtn insd main sbrtn"
       call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame + 1
       call switchto_NI_model_run_and_switchbackto_TN
       
       
       if (i.ne.Nstp) then
          l0(middle_spr(1)) = Initl0(1)
          l0(middle_spr(2)) = Initl0(2)
       endif
    enddo
    
    write(134,fmt=*) Exprmnt,Frame,"Exprmnt,Frame"
    write(134,fmt=*) " "
    close(134)
    
  end subroutine stepwise_shortening_MidSpr_frmtop_and_Equilibrate
  
  
  
  subroutine stepwise_relaxing_DiagSpr_and_Equilibrate(Exprmnt,Frame,hmuch)
    implicit none
    integer, intent(in) :: Exprmnt
    integer, intent(inout) :: Frame
    real*8 , intent(in) :: hmuch
    
    integer :: Nstp
    real*8  :: stepS
    integer :: diag_spr(1:2)
    real*8  :: Initl0(1:2)
    integer :: i
    real*8  :: hm
    
    open(unit=134,file="DiagSpr_Relax.dat")
    
    Nstp  = 5
    stepS = (hmuch)/(Nstp-1)
    
    call get_the_DiagSpr(diag_spr)
    
    Initl0(1) = l0(diag_spr(1))
    Initl0(2) = l0(diag_spr(2))
    
    write(134,fmt=*) Initl0(1:2),"Initl0"
    
    do i = 1,Nstp
       hm = (i-1)*stepS
       call shortening_sprPair(diag_spr,(-hm))
       write(134,fmt=*) diag_spr(1),diag_spr(2),"diag sprs"
       write(134,fmt=*) l0(diag_spr(1)),l0(diag_spr(2)),"l0 aft shrtn insd main sbrtn"
       call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame + 1
       call switchto_NI_model_run_and_switchbackto_TN
       
       if (i.ne.Nstp) then
          l0(diag_spr(1)) = Initl0(1)
          l0(diag_spr(2)) = Initl0(2)
       endif
    enddo
    
    write(134,fmt=*) Exprmnt,Frame,"Exprmnt,Frame"
    write(134,fmt=*) " "
    close(134)
    
  end subroutine stepwise_relaxing_DiagSpr_and_Equilibrate
  
  
  
  subroutine Pressure_adjstmnt_decision(dum_Nodes,dumA0)
    implicit none
    real*8, intent(in) :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8, intent(in) :: dumA0(1:N_cell) 

    real*8  :: Pressure_curr(1:N_cell)
    integer :: cmprRnkin(1:2,1:Hlf_Ncell)
    integer :: i,j
    integer :: cell
    
    interface
       function Pressure(dum_Nodes,dumA0)
         use system_parameters
         use transfrm_info
         implicit none
         real*8  :: Pressure(1:N_cell)
         real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
         real*8  :: dumA0(1:N_cell)
       end function Pressure
    end interface
    
    
    Pressure_curr(1:N_cell) = Pressure(dum_Nodes,dumA0)
    
    open(unit=145,file='AA0_fctr_inPressureAdjst.dat')
    
    write(145,fmt=*) "A","","A0","","fctr","","Pressure","","cell_no"
    
    do i = 1,2
       if(i==1) write(145,fmt=*) "Left Side Cell :"
       if(i==2) write(145,fmt=*) "Right Side Cell :" 

       do j = 1,Hlf_Ncell
          if(i==1) then
             cell = lftSideCell(j)
             write(145,fmt=*) A(cell),A0(cell),fctr(cell),Pressure_curr(cell),cell
          elseif(i==2) then
             cell = rghtSideCell(j)
             write(145,fmt=*) A(cell),A0(cell),fctr(cell),Pressure_curr(cell),cell
          endif
       enddo
       
       write(145,fmt=*) " "
    enddo
    
    close(145)
    
    call get_cmprRnkin(Pressure_curr,cmprRnkin)
    
    do i = 1,2
       do j = 1,Hlf_Ncell
          if (cmprRnkin(i,j).ne.1) then
             write(*,*) "Have to be adjusted"
             write(*,*) cmprRnkin(i,j),i,j,"cmprRnkin,i,j"
             !stop
          endif
          
       enddo
    enddo
    
  end subroutine Pressure_adjstmnt_decision
  
  
  subroutine get_cmprRnkin(Pressure_curr,cmprRnkin)
    implicit none
    real*8, intent(in)  :: Pressure_curr(1:N_cell) 
    integer,intent(out) :: cmprRnkin(1:2,1:Hlf_Ncell)
    
    integer :: PrpsdRnkin(1:2,1:Hlf_Ncell)
    integer :: CurrRnkin(1:2,1:Hlf_Ncell)
    
    
    call get_the_ProposedRanking(PrpsdRnkin)
    call get_the_CurrntRanking(Pressure_curr,CurrRnkin)
    call compare_Proposed_and_CurrntRanking(PrpsdRnkin,CurrRnkin,cmprRnkin)
    
  end subroutine get_cmprRnkin
  
  
  subroutine get_the_ProposedRanking(PrpsdRnkin)
    implicit none
    integer, intent(out) :: PrpsdRnkin(1:2,1:Hlf_Ncell)
    integer :: i
    
    open(unit=189,file='Prpsed_Rnkn.dat')
    
    write(189,fmt=*) strctNo,"strctNo"
    
    do i = 1,2
       
       if (strctNo==1) then
          PrpsdRnkin(i,1) = 3
          PrpsdRnkin(i,2) = 4
          
          PrpsdRnkin(i,3) = 5
          PrpsdRnkin(i,4) = 6
          PrpsdRnkin(i,5) = 8
          PrpsdRnkin(i,6) = 7
          PrpsdRnkin(i,7) = 1
          PrpsdRnkin(i,8) = 2
          
       elseif (strctNo==2) then
          PrpsdRnkin(i,1) = 4
          
          PrpsdRnkin(i,2) = 5
          PrpsdRnkin(i,3) = 6
          PrpsdRnkin(i,4) = 8
          PrpsdRnkin(i,5) = 7
          PrpsdRnkin(i,6) = 1
          PrpsdRnkin(i,7) = 2
          PrpsdRnkin(i,8) = 3
          
       endif
       
       write(189,fmt=*) PrpsdRnkin(i,1:Hlf_Ncell),"for lft/Rght depends on i=",i
       
    enddo
    
    close(189)
    
  end subroutine get_the_ProposedRanking
  
  
  subroutine get_the_CurrntRanking(Pressure_curr,CurrRnkin)
    implicit none
    real*8 , intent(in)  :: Pressure_curr(1:N_cell)
    integer, intent(out) :: CurrRnkin(1:2,Hlf_Ncell)
    
    real*8  :: inpt_array(1:Hlf_Ncell)
    integer :: Ranking(1:Hlf_Ncell)
    integer :: i,j
    
    open(unit=129,file='CurrRanking.dat')
    
    do i = 1,2
       
       do j = 1,Hlf_Ncell
          if (i==1) inpt_array(j) = Pressure_curr(lftSideCell(j))
          if (i==2) inpt_array(j) = Pressure_curr(rghtSideCell(j))
       enddo
       
       call get_Ranking(inpt_array,Hlf_Ncell,Ranking)     
       CurrRnkin(i,1:Hlf_Ncell) = Ranking(1:Hlf_Ncell)
       write(129,fmt=*) Ranking(1:Hlf_Ncell),"Ranking for i =",i
       write(129,fmt=*) CurrRnkin(i,1:Hlf_Ncell),"CurrRnkin for i=",i
    enddo
    
    close(129)
    
  end subroutine get_the_CurrntRanking
  
  
  subroutine compare_Proposed_and_CurrntRanking(PrpsdRnkin,CurrRnkin,cmprRnkin)
    implicit none
    integer, intent(in)  :: PrpsdRnkin(1:2,1:Hlf_Ncell)
    integer, intent(in)  :: CurrRnkin(1:2,1:Hlf_Ncell)
    integer, intent(out) :: cmprRnkin(1:2,1:Hlf_Ncell)
    
    integer :: i,j
    
    open(unit=146,file='cmprRnkin.dat')
    
    write(146,fmt=*) strctNo,"strctNo"
    
    do i = 1,2
       do j = 1,Hlf_Ncell
          if (PrpsdRnkin(i,j)==CurrRnkin(i,j)) then
             cmprRnkin(i,j) = 1
          elseif (PrpsdRnkin(i,j).ne.CurrRnkin(i,j)) then
             cmprRnkin(i,j) = -1
          endif
          write(146,fmt=*) cmprRnkin(i,j)
       enddo
    enddo
    
    close(146)
    
  end subroutine compare_Proposed_and_CurrntRanking
  
  subroutine cmpr_and_decide_to_reduce_tilt(Exprmnt,Frame,thrshAngl,curr_thrshAngl,thrshAnglInfo,thrsh_lgcl)
    implicit none
    integer,intent(in)  :: Exprmnt
    integer,intent(inout) :: Frame
    real*8, intent(in)  :: thrshAngl(1:N_thrshAngl)
    integer,intent(in)  :: thrshAnglInfo(1:N_thrshAngl,1:2)
    real*8, intent(out) :: curr_thrshAngl(1:N_thrshAngl)
    
    logical,intent(inout) :: thrsh_lgcl
    
    real*8  :: hmuch
    
    thrsh_lgcl = .False.
    
    call coordntes_to_nodes(coordntes_xy,node_xy)
    call get_curr_thrshldAngls(coordntes_xy,thrshAnglInfo,curr_thrshAngl)
    call cmpr_thrshAngls(curr_thrshAngl,thrshAngl,thrsh_lgcl)
    
    if (thrsh_lgcl .eqv. .True.) then
       hmuch = 0.15d0
       
       open(unit=908,file='frmcChk.dat')
       write(908,fmt=*) Frame,"Frame no bfr ApclShrtn"
       call StepwiseshortenAllApclSpr_equilibrate_and_cmpr(Exprmnt,Frame,hmuch,thrshAnglInfo,thrshAngl,curr_thrshAngl,thrsh_lgcl)
       write(908,fmt=*) Frame,"Frame no aft ApclShrtn"
       close(908)
    endif
    
  end subroutine cmpr_and_decide_to_reduce_tilt
  
  subroutine StepwiseshortenAllApclSpr_equilibrate_and_cmpr(Exprmnt,Frame,hmuch,thrshAnglInfo,thrshAngl,curr_thrshAngl,thrsh_lgcl)
    implicit none
    
    integer, intent(in)    :: Exprmnt
    integer, intent(inout) :: Frame
    real*8 , intent(in)    :: hmuch
    
    integer , intent(in)   :: thrshAnglInfo(1:N_thrshAngl,1:2)
    real*8 , intent(in)    :: thrshAngl(1:N_thrshAngl)
    real*8 , intent(inout) :: curr_thrshAngl(1:N_thrshAngl)
    logical, intent(inout) :: thrsh_lgcl
    
    integer :: Nstp
    real*8  :: stepS
    integer :: apcl_spr(1:N_cell)
    real*8  :: Initl0(1:N_cell)
    integer :: i,j
    real*8  :: hm
    
    open(unit=135,file="ApclSpr_Short.dat")
    open(unit=136,file="Insd_StepwiseshortenAllApclSpr_equilibrate_and_cmpr.dat")
    
    Nstp  = 5
    stepS = (hmuch)/(Nstp-1)
    
    call get_the_ApclSpr(apcl_spr)
    
    do i = 1,N_cell
       
       if (apcl_spr(i) .ne. (-1)) then
          Initl0(i) = l0(apcl_spr(i))
          write(135,fmt=*) Initl0(i),apcl_spr(i),i,"Initl0,Apcl_spr,cell"
          
       elseif (apcl_spr(i)==(-1)) then
          Initl0(i) = 0.00d0
          write(135,fmt=*) Initl0(i),apcl_spr(i),i,"Initl0,Apcl_spr,cell"
       endif
       
    enddo
    
    
    do i = 1,Nstp
       
       hm = (i-1)*stepS
       call shortening_ApclSpr(apcl_spr,hm)
       call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
       call Equilibrate_system
       
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame + 1
       call switchto_NI_model_run_and_switchbackto_TN
       
       call get_curr_thrshldAngls(coordntes_xy,thrshAnglInfo,curr_thrshAngl)
       call cmpr_thrshAngls(curr_thrshAngl,thrshAngl,thrsh_lgcl)
       
       write(136,fmt=*) thrsh_lgcl,i,"thrsh_lgcl,step"
       
       if (thrsh_lgcl.eqv..False.) then
          write(136,fmt=*) i,"step at which tilt disappear"
          exit
       endif
       
       if (i.ne.Nstp) then
          do j = 1,N_cell
             if (apcl_spr(j) .ne. (-1)) then
                l0(apcl_spr(j)) = Initl0(j)
             elseif (apcl_spr(j) .ne. (-1)) then
                continue
             endif
          enddo
       endif
       
    enddo
    
    write(135,fmt=*) Exprmnt,Frame,"Exprmnt,Frame"
    write(135,fmt=*) " "
    
    close(135)
    close(136)
    
  end subroutine StepwiseshortenAllApclSpr_equilibrate_and_cmpr
  

  subroutine StepWiseChangingl0_AP1AP2_andthen_BP3_andthenLP4(Exprmnt,Frame,fracChangeLP4)
    
    implicit none
    integer, intent(in) :: Exprmnt
    integer, intent(inout) :: Frame
    real*8 , intent(inout) :: fracChangeLP4(1:2)
    
    integer :: AP1(1:2),AP2(1:2),BP3(1:2),LP4(1:2)
    integer :: N_sprUExp !Num of Spr under Expr
    integer :: Nstp
    integer :: i
    
    real*8, allocatable :: l0_curr(:),l0_fin(:),l0_rngeV(:)
    real*8 :: stepSizeV,stepV
    
    open(unit=149,file="StpWChngl0_AP1AP2BP3LP4.dat")
    
    N_sprUExp = 8 ; Nstp = 5
    allocate(l0_curr(1:N_sprUExp),l0_fin(1:N_sprUExp),l0_rngeV(1:N_sprUExp))
    
    call get_first_and_secnd_pair_ApclSpr(AP1,AP2)
    call get_third_pair_BaslSpr(BP3)
    call get_fourth_pair_LtrlSpr(LP4)
    
    write(149,fmt=*) AP1(1:2),"Apcl pair1"
    write(149,fmt=*) AP2(1:2),"Apcl pair2" 
    write(149,fmt=*) BP3(1:2),"Basl pair3"
    write(149,fmt=*) LP4(1:2),"Ltrl pair4"
    
    fracChangeLP4(1) = 1.46d0 !
    fracChangeLP4(2) = 1.46d0 !
    
    call get_the_finVal_ofl0(AP1,AP2,BP3,LP4,N_sprUExp,fracChangeLP4,l0_fin) !prvVal_bfr_trnsfrmtn=Final value here
    
    l0_curr(1) = l0(AP1(1))
    l0_curr(2) = l0(AP1(2))
    l0_curr(3) = l0(AP2(1))
    l0_curr(4) = l0(AP2(2))
    l0_curr(5) = l0(BP3(1))
    l0_curr(6) = l0(BP3(2))
    l0_curr(7) = l0(LP4(1))
    l0_curr(8) = l0(LP4(2))
    
    l0_rngeV = l0_fin - l0_curr
    
    do i = 1,N_sprUExp
       write(149,fmt=*) l0_fin(i),l0_curr(i),l0_rngeV(i),"fin,curr,rngeV"
    enddo
    
    stepSizeV = 1.0d0/(Nstp-1)
    
    do i = 1,Nstp
       stepV = (i-1)*stepSizeV
       
       l0(AP1(1)) = l0_curr(1) + stepV*(l0_rngeV(1))
       l0(AP1(2)) = l0_curr(2) + stepV*(l0_rngeV(2))
       l0(AP2(1)) = l0_curr(3) + stepV*(l0_rngeV(3))
       l0(AP2(2)) = l0_curr(4) + stepV*(l0_rngeV(4))
       
       write(149,fmt=*) l0(AP1(1:2)),l0(AP2(1:2)),i,"l0 AP1,AP2 in step i =",i
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame + 1
       call switchto_NI_model_run_and_switchbackto_TN
       
    enddo
    
    write(149,fmt=*) " "

    stepSizeV = 1.0d0/(Nstp-1)
    
    do i = 1,Nstp
       stepV = (i-1)*stepSizeV
       
       l0(BP3(1)) = l0_curr(5) + stepV*l0_rngeV(5)
       l0(BP3(2)) = l0_curr(6) + stepV*l0_rngeV(6)
       write(149,fmt=*) l0(BP3(1:2)),i,"l0 BP3, in step i =",i
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame + 1
       call switchto_NI_model_run_and_switchbackto_TN
       
    enddo
    
    write(149,fmt=*) " "

    stepSizeV = 1.0d0/(Nstp-1)
    
    do i = 1,Nstp
       stepV = (i-1)*stepSizeV
       
       l0(LP4(1)) = l0_curr(7) + stepV*l0_rngeV(7)
       l0(LP4(2)) = l0_curr(8) + stepV*l0_rngeV(8)
       write(149,fmt=*) l0(LP4(1:2)),i,"l0 LP4, in step i =",i
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame + 1
       call switchto_NI_model_run_and_switchbackto_TN
       
       ! if (i==(Nstp-1)) then
       !    funcChoice=1 ; N_variabls = N_mvCoordnte
       !    call gradient_coordntes(coordntes_xy,grdmv_xy)
       !    grd_mv(1:N_mvCoordnte) = grdmv_xy(1:N_mvCoordnte)
          
       !    Analtcl_or_Numrcl = 2
       !    call gradient_coordntes(coordntes_xy,grdmv_xy)
       !    grd_mv(1:N_mvCoordnte) = grdmv_xy(1:N_mvCoordnte)
          
       !    Analtcl_or_Numrcl = 1
       !    stop
       ! endif
       
    enddo
    
    write(149,fmt=*) " "
    
    close(149)
    
  end subroutine StepWiseChangingl0_AP1AP2_Andthen_BP3_AndthenLP4
  
  subroutine StepWise_unchangingl0_LP4(Exprmnt,Frame,fracChangeLP4)
    implicit none
    integer, intent(in) :: Exprmnt
    integer, intent(inout) :: Frame
    real*8 , intent(in) :: fracChangeLP4(1:2)
    
    integer :: LP4(1:2)
    integer :: N_sprUExp
    integer :: Nstp
    integer :: i
    
    real*8, allocatable :: l0_curr(:),l0_fin(:),l0_rngeV(:) 
    real*8 :: stepSizeV,stepV
    
    open(unit=149,file="StpWUnchngl0_LP4.dat")
    
    N_sprUExp = 2 ; Nstp = 5
    allocate(l0_curr(1:N_sprUExp),l0_fin(1:N_sprUExp),l0_rngeV(1:N_sprUExp))
    
    call get_fourth_pair_LtrlSpr(LP4)
    
    write(149,fmt=*) LP4(1:2),"Ltrl pair4"
    
    call get_the_finVal_for_unchanging(LP4,N_sprUExp,fracChangeLP4,l0_fin)
    
    l0_curr(1) = l0(LP4(1))
    l0_curr(2) = l0(LP4(2))
    
    l0_rngeV = l0_fin - l0_curr
    
    do i = 1,N_sprUExp
       write(149,fmt=*) l0_fin(i),l0_curr(i),l0_rngeV(i),"fin,curr,rngeV"
    enddo
    
    stepSizeV = 1.0d0/(Nstp-1)
    
    do i = 1,Nstp
       stepV = (i-1)*stepSizeV
       
       l0(LP4(1)) = l0_curr(1) + stepV*l0_rngeV(1)
       l0(LP4(2)) = l0_curr(2) + stepV*l0_rngeV(2)
       write(149,fmt=*) l0(LP4(1:2)),i,"l0 LP4, in step i =",i
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)

       if (Frame==68) then
          call print_the_PairNode_Prop(8,16)
          call print_the_PairNode_Prop(19,22)
          call print_the_PairNode_Prop(6,14)
       endif
       
       Frame = Frame + 1
       call switchto_NI_model_run_and_switchbackto_TN
       
    enddo
    
    write(149,fmt=*) " "
    
    close(149)
    
  end subroutine StepWise_unchangingl0_LP4
    
  subroutine get_first_and_secnd_pair_ApclSpr(AP1,AP2)
    implicit none
    integer, intent(out) :: AP1(1:2),AP2(1:2) ! AP=Apcl_Pair1
    integer :: apcl_spr(1:N_cell)
    integer :: cnt1,cnt2
    integer :: i,j
    
    
    call get_the_ApclSpr(apcl_spr)

    cnt1 = 1 ; cnt2 = 1
    
    do i = 1,2
       do j = 1,N_cell

          if (i==1) then
             if (j==1 .or. j==(ncl+1)) then
                AP1(cnt1) = apcl_spr(j)
                cnt1 = cnt1+1
             endif
          elseif (i==2) then
             if (j==2 .or. j==(ncl+2)) then
                AP2(cnt2) = apcl_spr(j)
                cnt2 = cnt2+1
             endif
          endif
          
       enddo
    enddo
    
  end subroutine get_first_and_secnd_pair_ApclSpr

  subroutine get_third_pair_BaslSpr(BP3)
    implicit none
    integer, intent(out) :: BP3(1:2)
    integer :: basl_spr(1:N_cell)
    integer :: cnt
    integer :: j
    
    call get_the_BaslSpr(basl_spr)
    
    cnt = 1
    
    do j = 1,N_cell
       if (j==3 .or. j==(ncl+3)) then
          BP3(cnt) = basl_spr(j)
          cnt = cnt+1
       endif
    enddo

  end subroutine get_third_pair_BaslSpr

  subroutine get_fourth_pair_LtrlSpr(LP4)
    implicit none
    integer, intent(out) :: LP4(1:2)
    integer :: ltrl_spr(1:N_cell)
    integer :: cell_lft,cell_rght
    integer :: cnt
    integer :: j
    
    call get_the_LtrlSpr(ltrl_spr)

    cell_lft  = lftSideCell(4)
    cell_rght = rghtSideCell(4)
    
    cnt = 1
    
    do j = 1,N_cell
       if (j==cell_lft .or. j==cell_rght) then
          LP4(cnt) = ltrl_spr(j)
          cnt = cnt+1
       endif
    enddo
    
  end subroutine get_fourth_pair_LtrlSpr
  
  subroutine get_the_finVal_ofl0(AP1,AP2,BP3,LP4,N_sprUExp,fracChangeLP4,l0_fin)
    implicit none
    integer, intent(in)  :: AP1(1:2),AP2(1:2) !AP=Apical Pair
    integer, intent(in)  :: BP3(1:2) !BP=Basal Pair
    integer, intent(in)  :: LP4(1:2) !LP=Lateral Pair
    integer, intent(in)  :: N_sprUExp
    real*8 , intent(in)  :: fracChangeLP4(1:2)
    real*8 , intent(out) :: l0_fin(1:N_sprUExp)
    integer :: SBI !SBI=strct before Interpolation 
    
    SBI = 1

    if (modelID==2) then
       write(*,*) "routine: get_the_finVal_ofl0 is not for modelID=2"
       stop
    endif
    
    l0_fin(1) = l0_strctTN(SBI,AP1(1)) 
    l0_fin(2) = l0_strctTN(SBI,AP1(2))
    l0_fin(3) = l0_strctTN(SBI,AP2(1))
    l0_fin(4) = l0_strctTN(SBI,AP2(2))
    l0_fin(5) = l0_strctTN(SBI,BP3(1))
    l0_fin(6) = l0_strctTN(SBI,BP3(2))
    
    l0_fin(7) = fracChangeLP4(1) * l0(LP4(1))
    l0_fin(8) = fracChangeLP4(2) * l0(LP4(2))

  end subroutine get_the_finVal_ofl0

  subroutine get_the_finVal_for_unchanging(LP4,N_sprUExp,fracChangeLP4,l0_fin)
    implicit none
    integer, intent(in)  :: LP4(1:2)
    integer, intent(in)  :: N_sprUExp
    real*8 , intent(in)  :: fracChangeLP4(1:2)
    real*8 , intent(out) :: l0_fin(1:N_sprUExp)
    
    l0_fin(1) = l0(LP4(1))/fracChangeLP4(1)
    l0_fin(2) = l0(LP4(2))/fracChangeLP4(2)
    
  end subroutine get_the_finVal_for_unchanging
  
  subroutine print_the_PairNode_Prop(node1,node2)
    implicit none
    integer, intent(in) :: node1,node2
    integer :: spr1,spr2,area1,area2
    integer :: i,imax,imax1,imax2
    
    open(unit=189,file='node_to_cnnctd_elmnt.dat',position='append')
    
    imax1 = node_spr(node1,0)
    imax2 = node_spr(node2,0)
    
    if (imax1.ne.imax2) then
       write(*,*) "PairNodes must have same no of cnnctd spr"
       continue
    endif
    
    imax = imax1
    
    do i = 1,imax
       spr1 = node_spr(node1,i)
       spr2 = node_spr(node2,i)
       write(189,fmt=*) spr1,spr2,"spr1,spr2"
       write(189,fmt=*) k_spr(spr1),k_spr(spr2),l0(spr1),l0(spr2),typ_spr(spr1),typ_spr(spr2)
    enddo
    
    write(189,fmt=*) " "
    
    imax1 = node_area(node1,0)
    imax2 = node_area(node2,0)
    
    if (imax1.ne.imax2) then
       write(*,*) "PairNodes must have same no of cnnctd area"
       continue
    endif
    
    imax = imax1
    
    do i = 1,imax
       area1 = node_area(node1,i)
       area2 = node_area(node2,i)
       write(189,fmt=*) area1,area2,"area1,area2"
       write(189,fmt=*) k_area(area1),k_area(area2),A0(area1),A0(area2),typ_area(area1),typ_area(area2)
    enddo

    write(189,fmt=*) " " 
    
    do i = 1,imax
       area1 = node_area(node1,i)
       area2 = node_area(node2,i)
       write(189,fmt=*) area1,area2,"area1,area2"
       write(189,fmt=*) k_phi(node1,1:4),"k_phi of node1"
       write(189,fmt=*) k_phi(node2,1:4),"k_phi of node2"
    enddo
    
    write(189,fmt=*) " "
    
    write(189,fmt=*) CgXNode(node1),CgXNode(node2),"CgXNode"
    
    write(189,fmt=*) " "
    write(189,fmt=*) " "
    
    close(189)
    
  end subroutine print_the_PairNode_Prop
  
  !stage1E_relatedRoutines
  
  subroutine Stepwise_SprShorteningIC_and_Equilibrate(Exprmnt,Frame,hmMltplL,hmMltplK)
    implicit none
    integer, intent(in) :: Exprmnt
    integer, intent(inout) :: Frame
    real*8 , intent(in) :: hmMltplL(1:3),hmMltplK(1:3) 
    
    real*8  :: hmuch
    integer :: Nstp
    real*8  :: stepS
    integer :: apclSprIC,baslSprIC
    integer :: ltrl1SprIC,ltrl2SprIC
    
    real*8  :: Initl0,Initl01,Initl02
    integer :: i,imax
    real*8  :: hm
    integer :: spr_nm,sprPair(1:2)
    integer :: cnt_ltrl
    real*8  :: TINY_val
    
    TINY_val = 1.0d-15
    
    !First Apical spr shorten,secnd basal spr relax then the lateral pair relax
    
    open(unit=134,file="ApclSprIC_Shortn.dat")
    open(unit=135,file="BaslSprIC_Shortn.dat")
    open(unit=136,file="LtrlPairIC_Shortn.dat")
    
    apclSprIC=0 ; baslSprIC=0 ; ltrl1SprIC=0 ; ltrl2SprIC=0
    cnt_ltrl = 0
    imax = area_spr(InitiatorCell,0)
    
    write(*,*) InitiatorCell,max_spr_area,"InitiatorCell,max_spr_area"
    
    do i = 1,imax
       spr_nm = area_spr(InitiatorCell,i) 
       write(*,*) spr_nm,"spr_nm"
       if (typ_spr(spr_nm)==3) then
          apclSprIC = spr_nm
       elseif (typ_spr(spr_nm)==4) then
          baslSprIC = spr_nm
       elseif (typ_spr(spr_nm)==5) then
          cnt_ltrl = cnt_ltrl+1
          if (cnt_ltrl==1) ltrl1sprIC = spr_nm
          if (cnt_ltrl==2) ltrl2sprIC = spr_nm
       endif
    enddo
    
    !APICAL SPRING SHORTEN
    
    if (apclSprIC.gt.0) then
       
       if (abs(hmMltplL(1)) .gt. TINY_val) then
          hmuch = hmMltplL(1)
          Nstp  = 5
          stepS = (hmuch)/(Nstp-1)
          
          Initl0 = l0(apclSprIC)
          write(134,fmt=*) Initl0,"Initl0"
          
          do i = 1,Nstp
             hm = (i-1)*stepS
             call shortening_singleSpr(apclSprIC,hm)
             write(134,fmt=*) apclSprIC,"apclSprIC"
             write(134,fmt=*) l0(apclSprIC),"l0 aft shrtn insd main sbrtn"
             
             call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
             call Equilibrate_system
             call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
             Frame = Frame + 1
             call switchto_NI_model_run_and_switchbackto_TN
             
             if (i.ne.Nstp) then
                l0(apclSprIC) = Initl0
             endif
          
          enddo
          
          write(134,fmt=*) Exprmnt,Frame,"Exprmnt,Frame"
          write(134,fmt=*) " "
          close(134)
       endif

       if (abs(hmMltplK(1)) .gt. TINY_val) then
          k_spr(apclSprIC) = hmMltplK(1)*k_spr(apclSprIC)
          call Equilibrate_system
          call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
          Frame = Frame+1
          call switchto_NI_model_run_and_switchbackto_TN
       endif
       
    endif
    
    !BASAL SPRING SHORTEN
    
    if (baslSprIC.gt.0) then
       
       if (abs(hmMltplL(2)) .gt. TINY_val) then
          hmuch = hmMltplL(2)
          Nstp  = 5
          stepS = (hmuch)/(Nstp-1)
          
          Initl0 = l0(baslSprIC)
          
          write(135,fmt=*) Initl0,"Initl0"
          
          do i = 1,Nstp
             hm = (i-1)*stepS
             call shortening_singleSpr(baslSprIC,hm)
             write(135,fmt=*) baslSprIC,"baslSprIC"
             write(135,fmt=*) l0(baslSprIC),"l0 aft shrtn insd main sbrtn"
             
             call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
             call Equilibrate_system
             call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
             Frame = Frame + 1
             call switchto_NI_model_run_and_switchbackto_TN
             
             if (i.ne.Nstp) then
                l0(baslSprIC) = Initl0
             endif
             
          enddo
          
          write(135,fmt=*) Exprmnt,Frame,"Exprmnt,Frame"
          write(135,fmt=*) " "
          close(135)
       
       endif
       
       if (abs(hmMltplK(2)) .gt. TINY_val) then
          k_spr(baslSprIC) = hmMltplK(2)*k_spr(baslSprIC)
          call Equilibrate_system
          call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
          Frame = Frame+1
          call switchto_NI_model_run_and_switchbackto_TN
       endif
       
    endif
    
    !LATERAL PAIR SPRING SHORTEN
    
    if ((ltrl1SprIC.gt.0).and.(ltrl2SprIC.gt.0)) then

       if (abs(hmMltplL(3)) .gt. TINY_val) then
          hmuch = hmMltplL(3)
          Nstp  = 5
          stepS = (hmuch)/(Nstp-1)
          
          Initl01 = l0(ltrl1SprIC)
          Initl02 = l0(ltrl2SprIC)
          
          sprPair(1) = ltrl1SprIC ; sprPair(2) = ltrl2SprIC
          
          write(136,fmt=*) Initl01,Initl02,"Initl0s"
          write(136,fmt=*) sprPair(1:2),"sprPair"
          
          do i = 1,Nstp
             
             hm = (i-1)*stepS
             call shortening_sprPair(sprPair,hm)
             write(136,*) ltrl1SprIC,ltrl2SprIC,"ltrlSprIC"
             write(136,*) l0(ltrl1SprIC),l0(ltrl2SprIC),"l0 aft shrtn insd main sbrtn"
             
             call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
             call Equilibrate_system
             call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
             Frame = Frame + 1
             call switchto_NI_model_run_and_switchbackto_TN
             
             if (i.ne.Nstp) then
                l0(ltrl1SprIC) = Initl01
                l0(ltrl2SprIC) = Initl02
             endif
             
          enddo
       
          write(136,fmt=*) Exprmnt,Frame,"Exprmnt,Frame"
          write(136,fmt=*) " "
          close(136)
          
       endif
       
       if (abs(hmMltplK(3)) .gt. TINY_val) then
          
          k_spr(ltrl1SprIC) = hmMltplK(3)*k_spr(ltrl1SprIC)
          k_spr(ltrl2SprIC) = hmMltplK(3)*k_spr(ltrl2SprIC)
          
          call Equilibrate_system
          call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
          Frame = Frame+1
          call switchto_NI_model_run_and_switchbackto_TN
       endif
    endif
    
  end subroutine Stepwise_SprShorteningIC_and_Equilibrate
  
  
  subroutine Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)!single Neighbouring cellPair
    implicit none
    integer, intent(in)    :: Exprmnt
    integer, intent(inout) :: Frame
    integer, intent(in)    :: NCP
    
    real*8 , intent(in) :: hmMltplL(1:3)
    
    integer :: Nstp
    real*8  :: stepS
    
    integer :: apclSprNC1,apclSprNC2
    integer :: baslSprNC1,baslSprNC2
    integer :: ltrlSpr1NC1,ltrlSpr1NC2
    
    real*8  :: Initl0,Initl01,Initl02
    integer :: i,j,jmax
    integer :: NC1,NC2
    real*8  :: hmuch,hm
    integer :: spr_nm,sprPair(1:2)
    integer :: cnt_ltrl
    real*8  :: TINYval
    integer :: spr1,spr2,spr3,spr4
    
    TINYval = 1.0d-15
    
    NC1=(ncl)-(NCP-1)
    NC2=(ncl+ncr)-(NCP-1)
    
    write(*,*) NC1,NC2,max_spr_area,"NC1,NC2,max_spr_area"
    write(*,*) area_spr(NC1,0:max_spr_area),"NC1"
    write(*,*) area_spr(NC2,0:max_spr_area),"NC2"
    
    if (area_spr(NC1,0)==4) then
       spr1 = area_spr(NC1,1) ; spr2 = area_spr(NC1,2)
       spr3 = area_spr(NC1,3) ; spr4 = area_spr(NC1,4)  
       write(*,*) typ_spr(spr1),typ_spr(spr2),typ_spr(spr3),typ_spr(spr4),"typ_spr of NC1"
       
    elseif (area_spr(NC1,0)==3) then
       spr1 = area_spr(NC1,1) ; spr2 = area_spr(NC1,2)
       spr3 = area_spr(NC1,3)   
       write(*,*) typ_spr(spr1),typ_spr(spr2),typ_spr(spr3),"typ_spr of NC1"
       
    endif
    
    if (area_spr(NC2,0)==4) then
       spr1 = area_spr(NC2,1) ; spr2 = area_spr(NC2,2)
       spr3 = area_spr(NC2,3) ; spr4 = area_spr(NC2,4)  
       write(*,*) typ_spr(spr1),typ_spr(spr2),typ_spr(spr3),typ_spr(spr4),"typ_spr of NC2"
       
    elseif (area_spr(NC2,0)==3) then
       spr1 = area_spr(NC2,1) ; spr2 = area_spr(NC2,2)
       spr3 = area_spr(NC2,3) 
       write(*,*) typ_spr(spr1),typ_spr(spr2),typ_spr(spr3),"typ_spr of NC2"
       
    endif
    
    !stop
    
    do i = 1,2
       cnt_ltrl = 0
       
       if (i==1) jmax = area_spr(NC1,0)
       if (i==2) jmax = area_spr(NC2,0)
       
       do j = 1,jmax
          
          if (i==1) spr_nm = area_spr(NC1,j) 
          if (i==2) spr_nm = area_spr(NC2,j)
          
          if (typ_spr(spr_nm)==1.or. typ_spr(spr_nm)==3.or.typ_spr(spr_nm)==6) then
             if (i==1) apclSprNC1 = spr_nm
             if (i==2) apclSprNC2 = spr_nm
             
          elseif (typ_spr(spr_nm)==2 .or. typ_spr(spr_nm)==4) then
             if (i==1) baslSprNC1 = spr_nm
             if (i==2) baslSprNC2 = spr_nm
             
          elseif (typ_spr(spr_nm)==5) then
             
             if (i==1) then
                cnt_ltrl = cnt_ltrl+1
                if (cnt_ltrl==1) ltrlSpr1NC1 = spr_nm
                if (cnt_ltrl==2) continue
             elseif (i==2) then
                cnt_ltrl = cnt_ltrl+1
                if (cnt_ltrl==1) ltrlSpr1NC2 = spr_nm
                if (cnt_ltrl==2) continue
             endif
             
          endif
          
       enddo
       
    enddo
    
    write(*,*) apclSprNC1,apclSprNC2,baslSprNC1,baslSprNC2,"apclNC1 & NC2,baslNC1 & NC2"
    write(*,*) ltrlSpr1NC1,ltrlSpr1NC2,"lt1 NC1 & NC2"
    
    !APICAL SPRING SHORTEN(+hm)/RELAX(-hm)
    
    hmuch = hmMltplL(1)
    
    if (abs(hmuch).gt.TINYval) then
       open(unit=134,file="ApclSprNC_Shortn.dat")
       
       Nstp  = 5
       stepS = (hmuch)/(Nstp-1)
       
       Initl01 = l0(apclSprNC1)
       Initl02 = l0(apclSprNC2)
       
       sprPair(1) = apclSprNC1 ; sprPair(2) = apclSprNC2
       
       write(134,fmt=*) Initl01,Initl02,"Initl0s"
       write(134,fmt=*) sprPair(1:2),"sprPair"
       
       do i = 1,Nstp
       
          hm = (i-1)*stepS
          call shortening_sprPair(sprPair,hm)
          write(134,fmt=*) l0(apclSprNC1),l0(apclSprNC2),"l0 aftShrtn insd main sb"
          
          call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
          call Equilibrate_system
          call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
          Frame = Frame + 1
          call switchto_NI_model_run_and_switchbackto_TN
          
          if (i.ne.Nstp) then
             l0(apclSprNC1) = Initl01
             l0(apclSprNC2) = Initl02
          endif
          
       enddo
       
       write(134,fmt=*) Exprmnt,Frame,"Exprmnt,Frame"
       write(134,fmt=*) " "
       close(134)
       
    endif

    !BASAL SPRING SHORTEN(+hm)/RELAX(-hm)
    
    hmuch = hmMltplL(2)
    
    if (abs(hmuch).gt.TINYval) then
       open(unit=135,file="BaslSprNC_Shortn.dat")
       
       Nstp  = 5
       stepS = (hmuch)/(Nstp-1)
       
       Initl01 = l0(baslSprNC1)
       Initl02 = l0(baslSprNC2)
    
       sprPair(1) = baslSprNC1 ; sprPair(2) = baslSprNC2
       
       write(135,fmt=*) Initl01,Initl02,"Initl0s"
       write(135,fmt=*) sprPair(1:2),"sprPair"
       
       do i = 1,Nstp
          
          hm = (i-1)*stepS
          call shortening_sprPair(sprPair,hm)
          write(135,fmt=*) l0(baslSprNC1),l0(baslSprNC2),"l0 AftShrtn insd main sb"
          
          call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
          call Equilibrate_system
          call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
          Frame = Frame + 1
          call switchto_NI_model_run_and_switchbackto_TN
          
          if (i.ne.Nstp) then
             l0(baslSprNC1) = Initl01
             l0(baslSprNC2) = Initl02
          endif
          
       enddo
    
       write(135,fmt=*) Exprmnt,Frame,"Exprmnt,Frame"
       write(135,fmt=*) " "
       close(135)
       
    endif
       
    !LATERAL SPRING1 SHORTEN(+hm)/RELAX(-hm)
    
    hmuch = hmMltplL(3)
    
    if (abs(hmuch).gt.TINYval) then
       open(unit=136,file="LtrlSpr1NC_Shortn.dat")
       
       Nstp  = 5
       stepS = (hmuch)/(Nstp-1)
       
       Initl01 = l0(ltrlSpr1NC1)
       Initl02 = l0(ltrlSpr1NC2)
       
       sprPair(1) = ltrlSpr1NC1 ; sprPair(2) = ltrlSpr1NC2
       
       write(136,fmt=*) Initl01,Initl02,"Initl0s"
       write(136,fmt=*) sprPair(1:2),"sprPair"
       
       do i = 1,Nstp
          
          hm = (i-1)*stepS
          call shortening_sprPair(sprPair,hm)
          write(136,*) l0(ltrlSpr1NC1),l0(ltrlSpr1NC2),"l0 AftShrtn insd main sb"
          
          call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
          call Equilibrate_system
          call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
          Frame = Frame + 1
          call switchto_NI_model_run_and_switchbackto_TN
          
          if (i.ne.Nstp) then
             l0(ltrlSpr1NC1) = Initl01
             l0(ltrlSpr1NC2) = Initl02
          endif
          
       enddo
       
       write(136,fmt=*) Exprmnt,Frame,"Exprmnt,Frame"
       write(136,fmt=*) " "
       close(136)
       
    endif
    
    
  end subroutine Stepwise_SprShorteningNC_and_Equilibrate
  
  

  subroutine StepwiseSprShorten_MltplNC_Equilibrate(Exprmnt,Frame,NCpairs,decideSpr,hmuch)
    implicit none
    integer, intent(in) :: Exprmnt
    integer, intent(inout) :: Frame
    integer, intent(in) :: NCpairs(1:2)
    integer, intent(in) :: decideSpr
    real*8 , intent(in) :: hmuch
    
    integer :: NCP11,NCP12 !NCP=Neighbouring Cell Pair
    integer :: NCP21,NCP22
    
    integer :: pairNo,i,j,jmax,cnt_ltrl
    integer :: spr_nm,sprPair(1:2)
    real*8  :: hm
    
    integer :: apclSprNCP11,apclSprNCP12
    integer :: apclSprNCP21,apclSprNCP22 
    
    integer :: baslSprNCP11,baslSprNCP12
    integer :: baslSprNCP21,baslSprNCP22
    
    integer :: ltrlSpr2NCP11,ltrlSpr2NCP12
    integer :: ltrlSpr2NCP21,ltrlSpr2NCP22
    
    integer :: Nstp
    real*8  :: stepS
    real*8  :: Initl0P11,Initl0P12,Initl0P21,Initl0P22
    integer :: NCp1,NCp2
    
    !decideSpr = 1 ! 1 for apcl, 2 for basal, 3 for ltrl2
    
    NCp1 = NCpairs(1) ; NCp2 = NCpairs(2)
    
    NCP11 = ncl-(NCp1-1) ; NCP12 = (ncl+ncr)-(NCp1-1)
    NCP21 = ncl-(NCp2-1) ; NCP22 = (ncl+ncr)-(NCp2-1)
    
    write(*,*) ncl,ncr,"ncl,ncr"
    write(*,*) NCP11,NCP12,NCP21,NCP22,"NCPs"
    write(*,*) area_spr(NCP11,0:max_spr_area),"area_spr for NCP11"
    write(*,*) area_spr(NCP12,0:max_spr_area),"area_spr for NCP12"
    write(*,*) area_spr(NCP21,0:max_spr_area),"area_spr_for_NCP21"
    write(*,*) area_spr(NCP22,0:max_spr_area),"area_spr_for_NCP22"
    
    
    do pairNo = 1,2
       
       do i=1,2
          
          if (pairNo==1) then
             if (i==1) jmax = area_spr(NCP11,0) 
             if (i==2) jmax = area_spr(NCP12,0)
          elseif (pairNo==2) then
             if (i==1) jmax = area_spr(NCP21,0) 
             if (i==2) jmax = area_spr(NCP22,0)
          endif
          
          cnt_ltrl = 0
          
          do j=1,jmax
             
             if (pairNo==1) then
                if (i==1) spr_nm = area_spr(NCP11,j) 
                if (i==2) spr_nm = area_spr(NCP12,j)
             elseif (pairNo==2) then
                if (i==1) spr_nm = area_spr(NCP21,j) 
                if (i==2) spr_nm = area_spr(NCP22,j)
             endif
             
             if (typ_spr(spr_nm)==3.or.typ_spr(spr_nm)==1.or.typ_spr(spr_nm)==6) then
                
                if (pairNo==1) then
                   if (i==1) apclSprNCP11 = spr_nm
                   if (i==2) apclSprNCP12 = spr_nm
                elseif (pairNo==2)then
                   if (i==1) apclSprNCP21 = spr_nm
                   if (i==2) apclSprNCP22 = spr_nm
                endif
                
             elseif (typ_spr(spr_nm)==4.or.typ_spr(spr_nm)==2) then
                
                if (pairNo==1) then
                   if (i==1) baslSprNCP11 = spr_nm
                   if (i==2) baslSprNCP12 = spr_nm
                elseif (pairNo==2)then
                   if (i==1) baslSprNCP21 = spr_nm
                   if (i==2) baslSprNCP22 = spr_nm
                endif
                
             elseif (typ_spr(spr_nm)==5) then
                cnt_ltrl = cnt_ltrl+1
                
                if (cnt_ltrl==1) then
                   continue
                elseif (cnt_ltrl==2) then
                   
                   if (pairNo==1) then
                      if (i==1) ltrlSpr2NCP11 = spr_nm
                      if (i==2) ltrlSpr2NCP12 = spr_nm
                   elseif (pairNo==2)then
                      if (i==1) ltrlSpr2NCP21 = spr_nm
                      if (i==2) ltrlSpr2NCP22 = spr_nm
                   endif
                   
                endif
                
             endif
             
          enddo
       enddo
    enddo
    
    open(unit=135,file="SprNCP1NCP2_Shortn.dat")
    write(135,*) ltrlSpr2NCP11,ltrlSpr2NCP12,ltrlSpr2NCP21,ltrlSpr2NCP22,"lt2"
    
    !SPRING SHORTEN
    
    Nstp  = 5
    stepS = (hmuch)/(Nstp-1)
    
    if (decideSpr==1) then
       Initl0P11 = l0(apclSprNCP11)
       Initl0P12 = l0(apclSprNCP12)
       Initl0P21 = l0(apclSprNCP21)
       Initl0P22 = l0(apclSprNCP22)
    elseif (decideSpr==2) then
       Initl0P11 = l0(baslSprNCP11)
       Initl0P12 = l0(baslSprNCP12)
       Initl0P21 = l0(baslSprNCP21)
       Initl0P22 = l0(baslSprNCP22)
    elseif (decideSpr==3) then
       Initl0P11 = l0(ltrlSpr2NCP11)
       Initl0P12 = l0(ltrlSpr2NCP12)
       Initl0P21 = l0(ltrlSpr2NCP21)
       Initl0P22 = l0(ltrlSpr2NCP22)
    endif

    do i = 1,Nstp
       
       hm = (i-1)*stepS

       if (decideSpr==1) then
          sprPair(1) = apclSprNCP11  ; sprPair(2) = apclSprNCP12
       elseif (decideSpr==2) then
          sprPair(1) = baslSprNCP11  ; sprPair(2) = baslSprNCP12
       elseif (decideSpr==3) then
          sprPair(1) = ltrlSpr2NCP11 ; sprPair(2) = ltrlSpr2NCP12
       endif
       
       call shortening_sprPair(sprPair,hm)
       
       if (decideSpr==1) then
          write(135,*) l0(apclSprNCP11),l0(apclSprNCP12),sprPair
       elseif (decideSpr==2) then
          write(135,*) l0(baslSprNCP11),l0(baslSprNCP12),sprPair
       elseif (decideSpr==3) then
          write(135,*) l0(ltrlSpr2NCP11),l0(ltrlSpr2NCP12),sprPair
       endif
       
       
       
       if (decideSpr==1) then
          sprPair(1) = apclSprNCP21  ; sprPair(2) = apclSprNCP22
       elseif (decideSpr==2) then
          sprPair(1) = baslSprNCP21  ; sprPair(2) = baslSprNCP22
       elseif (decideSpr==3) then
          sprPair(1) = ltrlSpr2NCP21 ; sprPair(2) = ltrlSpr2NCP22
       endif
       
       call shortening_sprPair(sprPair,hm)

       if (decideSpr==1) then
          write(135,*) l0(apclSprNCP21),l0(apclSprNCP22),sprPair
       elseif (decideSpr==2) then
          write(135,*) l0(baslSprNCP21),l0(baslSprNCP22),sprPair
       elseif (decideSpr==3) then
          write(135,*) l0(ltrlSpr2NCP21),l0(ltrlSpr2NCP22),sprPair
       endif
       
       call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame + 1
       call switchto_NI_model_run_and_switchbackto_TN
       
       if (i.ne.Nstp) then
          
          if (decideSpr==1) then
             l0(apclSprNCP11) = Initl0P11
             l0(apclSprNCP12) = Initl0P12
             l0(apclSprNCP21) = Initl0P21
             l0(apclSprNCP22) = Initl0P22  
          elseif (decideSpr==2) then
             l0(baslSprNCP11) = Initl0P11
             l0(baslSprNCP12) = Initl0P12
             l0(baslSprNCP21) = Initl0P21
             l0(baslSprNCP22) = Initl0P22  
          elseif (decideSpr==3) then
             l0(ltrlSpr2NCP11) = Initl0P11
             l0(ltrlSpr2NCP12) = Initl0P12
             l0(ltrlSpr2NCP21) = Initl0P21
             l0(ltrlSpr2NCP22) = Initl0P22  
          endif
          
       endif
       
    enddo
    
    close(135)
    
  end subroutine StepwiseSprShorten_MltplNC_Equilibrate
  
  subroutine Stepwise_ExpandingIC(Exprmnt,Frame,hmuch)
    implicit none
    integer, intent(inout) :: Exprmnt,Frame
    real*8,  intent(in)    :: hmuch
    
    real*8  :: InitA0,FinlA0
    integer :: Nstep
    real*8  :: UnitStep,stepSize
    integer :: i
    
    write(*,*) Exprmnt,Frame,"Exp and Frm in changeA0 in IC"
    
    InitA0 = A0(InitiatorCell)
    FinlA0 = (1.0d0+hmuch)*A0(InitiatorCell)
    write(*,*) InitA0,FinlA0,"InitA0,FinlA0"
    
    Nstep = 5
    Unitstep = (1.0d0)/(Nstep-1)
    write(*,*) Unitstep,"UnitStep"
    
    do i = 1,Nstep
       stepSize = (i-1)*UnitStep
       A0(InitiatorCell) = InitA0 + stepSize*(FinlA0-InitA0)
       write(*,*) A0(InitiatorCell),"A0 of InitiatorCell"
       
       call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
       
       call Equilibrate_system
       call coordntes_to_nodes(coordntes_xy,node_xy)
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
    enddo
    
  end subroutine Stepwise_ExpandingIC

  
  subroutine Stepwise_ExpandingNC(Exprmnt,Frame,NC,hmuch)
    implicit none
    integer, intent(in)    :: Exprmnt
    integer, intent(inout) :: Frame
    integer, intent(in)    :: NC
    real*8,  intent(in)    :: hmuch
    
    integer :: NC1,NC2
    real*8  :: InitA01,InitA02
    real*8  :: FinlA01,FinlA02
    integer :: Nstep
    real*8  :: UnitStep,stepSize
    integer :: i
    real*8  :: TINYval

    open(unit=67,file='ExpandingNC.dat')
    
    write(67,*) Exprmnt,Frame,"Exp and Frm in changeA0 in NC"
    
    TINYval = 1.0d-15
    
    NC1=(ncl)-(NC-1)
    NC2=(ncl+ncr)-(NC-1)
    write(67,*) NC1,NC2,"cell1,cell2"
    
    InitA01 = A0(NC1) ; InitA02 = A0(NC2)
    
    FinlA01 = (1.0d0+hmuch)*A0(NC1)
    FinlA02 = (1.0d0+hmuch)*A0(NC2)
    
    write(67,*) InitA01,FinlA01,"InitA01,FinlA01"
    write(67,*) InitA02,FinlA02,"InitA02,FinlA02"
    
    Nstep = 5
    Unitstep = (1.0d0)/(Nstep-1)
    write(67,*) Unitstep,"UnitStep"
    
    do i = 1,Nstep
       stepSize = (i-1)*UnitStep
       
       A0(NC1) = InitA01 + stepSize*(FinlA01-InitA01)
       A0(NC2) = InitA02 + stepSize*(FinlA02-InitA02)
       
       write(67,*) A0(NC1),"A0 of NC1"
       write(67,*) A0(NC2),"A0 of NC2"
       
       call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
       
       call Equilibrate_system
       call coordntes_to_nodes(coordntes_xy,node_xy)
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
    enddo

    close(67)
    
  end subroutine Stepwise_ExpandingNC
  
  subroutine adjust_BsalLtrlSprIC_and_Equilibrate(Exprmnt,Frame,hmuch)
    implicit none
    integer, intent(in) :: Exprmnt
    integer, intent(inout) :: Frame
    real*8 , intent(in) :: hmuch
    
    real*8  :: hm
    integer :: Nstp
    real*8  :: stepS
    
    integer :: baslSpr,ltrlSpr1,ltrlSpr2
    real*8  :: Initl0bsl,Initl0lt1,Initl0lt2
    integer :: i,imax
    integer :: cnt_ltrl,spr_nm
    
    open(unit=134,file="BsalLtrlSpr_Shortn.dat")
    
    imax = area_spr(InitiatorCell,0)
    cnt_ltrl = 0
    
    do i = 1,imax
       spr_nm = area_spr(InitiatorCell,i)
       
       if (typ_spr(spr_nm)==4) then
          baslSpr = area_spr(InitiatorCell,i)
       elseif (typ_spr(spr_nm)==5) then
          cnt_ltrl=cnt_ltrl+1
          
          if (cnt_ltrl==1) then
             ltrlSpr1 = area_spr(InitiatorCell,i)
          elseif (cnt_ltrl==2) then
             ltrlSpr2 = area_spr(InitiatorCell,i)
          endif
          
       endif
       
    enddo
    
    write(134,fmt=*) baslSpr,ltrlSpr1,ltrlSpr2,"baslSpr,ltrlSpr1and2"

    !k_spr(baslSpr)  =  0.2*k_spr(baslSpr)
    !k_spr(ltrlSpr1) =  8.0*k_spr(ltrlSpr1)
    !k_spr(ltrlSpr2) =  8.0*k_spr(ltrlSpr2)
    
    Nstp  = 5
    stepS = (hmuch)/(Nstp-1)
    
    Initl0bsl = l0(baslSpr)
    Initl0lt1 = l0(ltrlSpr1)
    Initl0lt2 = l0(ltrlSpr2)
    
    
    write(134,fmt=*) Initl0bsl,Initl0lt1,Initl0lt2,"Initl0s"
    
    do i = 1,Nstp
       
       hm = (i-1)*stepS
       call shortening_singleSpr(baslSpr,hm)
       call shortening_singleSpr(ltrlSpr1,hm)
       call shortening_singleSpr(ltrlSpr2,hm)
       
       write(134,*) l0(baslSpr),l0(ltrlSpr1),l0(ltrlSpr2),"l0 aftShrtn insd sb"
       
       call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame + 1
       call switchto_NI_model_run_and_switchbackto_TN
       
       if (i.ne.Nstp) then
          l0(baslSpr)  = Initl0bsl
          l0(ltrlSpr1) = Initl0lt1
          l0(ltrlSpr2) = Initl0lt2
       endif
       
    enddo
    
    write(134,*) Exprmnt,Frame,"Exprmnt,Frame"
    write(134,*) " "
    close(134)
    
  end subroutine adjust_BsalLtrlSprIC_and_Equilibrate
  
  
  subroutine Stepwise_SprShorteningNC_and_Equilibrate_UntilMet(Exprmnt,Frame,NCP,hmMltplL)
    implicit none
    integer, intent(in)    :: Exprmnt
    integer, intent(inout) :: Frame
    integer, intent(in)    :: NCP
    
    real*8, intent(out) :: hmMltplL(1:3)
    
    integer :: Nstp
    real*8  :: stepS
    
    integer :: apclSprNC1,apclSprNC2
    integer :: baslSprNC1,baslSprNC2
    integer :: ltrlSpr1NC1,ltrlSpr1NC2
    
    real*8  :: Initl0,Initl01,Initl02
    integer :: i,j,jmax
    integer :: NC1,NC2
    real*8  :: hmuch,hm
    integer :: spr_nm,sprPair(1:2)
    integer :: cnt_ltrl
    real*8  :: TINYval
    integer :: spr1,spr2,spr3,spr4
    
    logical :: lgcl_Rdfn,lgcl_Rdfn1
    real*8  :: tol_Rdfn
    
    integer :: ProgCycl,strctV1,strctV2,strctV
    integer :: strghtOrNI
    
    TINYval = 1.0d-15
    
    NC1=(ncl)-(NCP-1)
    NC2=(ncl+ncr)-(NCP-1)
    
    write(*,*) NC1,NC2,max_spr_area,"NC1,NC2,max_spr_area"
    write(*,*) area_spr(NC1,0:max_spr_area),"NC1"
    write(*,*) area_spr(NC2,0:max_spr_area),"NC2"
    
    if (area_spr(NC1,0)==4) then
       spr1 = area_spr(NC1,1) ; spr2 = area_spr(NC1,2)
       spr3 = area_spr(NC1,3) ; spr4 = area_spr(NC1,4)  
       write(*,*) typ_spr(spr1),typ_spr(spr2),typ_spr(spr3),typ_spr(spr4),"typ_spr of NC1"
       
    elseif (area_spr(NC1,0)==3) then
       spr1 = area_spr(NC1,1) ; spr2 = area_spr(NC1,2)
       spr3 = area_spr(NC1,3)   
       write(*,*) typ_spr(spr1),typ_spr(spr2),typ_spr(spr3),"typ_spr of NC1"
       
    endif
    
    if (area_spr(NC2,0)==4) then
       spr1 = area_spr(NC2,1) ; spr2 = area_spr(NC2,2)
       spr3 = area_spr(NC2,3) ; spr4 = area_spr(NC2,4)  
       write(*,*) typ_spr(spr1),typ_spr(spr2),typ_spr(spr3),typ_spr(spr4),"typ_spr of NC2"
       
    elseif (area_spr(NC2,0)==3) then
       spr1 = area_spr(NC2,1) ; spr2 = area_spr(NC2,2)
       spr3 = area_spr(NC2,3) 
       write(*,*) typ_spr(spr1),typ_spr(spr2),typ_spr(spr3),"typ_spr of NC2"
       
    endif
    
    !stop
    
    do i = 1,2
       cnt_ltrl = 0
       
       if (i==1) jmax = area_spr(NC1,0)
       if (i==2) jmax = area_spr(NC2,0)
       
       do j = 1,jmax
          
          if (i==1) spr_nm = area_spr(NC1,j) 
          if (i==2) spr_nm = area_spr(NC2,j)
          
          if (typ_spr(spr_nm)==1.or. typ_spr(spr_nm)==3.or.typ_spr(spr_nm)==6) then
             if (i==1) apclSprNC1 = spr_nm
             if (i==2) apclSprNC2 = spr_nm
             
          elseif (typ_spr(spr_nm)==2 .or. typ_spr(spr_nm)==4) then
             if (i==1) baslSprNC1 = spr_nm
             if (i==2) baslSprNC2 = spr_nm
             
          elseif (typ_spr(spr_nm)==5) then
             
             if (i==1) then
                cnt_ltrl = cnt_ltrl+1
                if (cnt_ltrl==1) ltrlSpr1NC1 = spr_nm
                if (cnt_ltrl==2) continue
             elseif (i==2) then
                cnt_ltrl = cnt_ltrl+1
                if (cnt_ltrl==1) ltrlSpr1NC2 = spr_nm
                if (cnt_ltrl==2) continue
             endif
             
          endif
          
       enddo
       
    enddo
    
    write(*,*) apclSprNC1,apclSprNC2,baslSprNC1,baslSprNC2,"apclNC1 & NC2,baslNC1 & NC2"
    write(*,*) ltrlSpr1NC1,ltrlSpr1NC2,"lt1 NC1 & NC2"
    
    !APICAL SPRING SHORTEN(+hm)/RELAX(-hm)
    
    hmuch = hmMltplL(1)
    
    if (abs(hmuch).gt.TINYval) then
       open(unit=134,file="ApclSprNC_Shortn.dat")
       
       Nstp  = 5
       stepS = (hmuch)/(Nstp-1)
       
       Initl01 = l0(apclSprNC1)
       Initl02 = l0(apclSprNC2)
       
       sprPair(1) = apclSprNC1 ; sprPair(2) = apclSprNC2
       
       write(134,fmt=*) Initl01,Initl02,"Initl0s"
       write(134,fmt=*) sprPair(1:2),"sprPair"
       
       do i = 1,Nstp
          
          hm = (i-1)*stepS
          call shortening_sprPair(sprPair,hm)
          write(134,fmt=*) l0(apclSprNC1),l0(apclSprNC2),"l0 aftShrtn insd main sb"
          
          call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
          call Equilibrate_system
          call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
          Frame = Frame + 1
          call switchto_NI_model_run_and_switchbackto_TN
          
          if (i.ne.Nstp) then
             l0(apclSprNC1) = Initl01
             l0(apclSprNC2) = Initl02
          endif
          
          !lgcl_Rdfn = .False.
          !tol_Rdfn = 0.05d0
          !call get_decision_of_redefining(tol_Rdfn,lgcl_Rdfn)
          
          !if (lgcl_Rdfn.eqv..True.) then
           !  write(134,*) (Frame-1),hm,"FrameNo,hm"
            ! exit
          !elseif (lgcl_Rdfn.eqv..False.) then
           !  continue
          !endif
          
          
       enddo
       
       write(134,fmt=*) Exprmnt,Frame,"Exprmnt,Frame"
       write(134,fmt=*) " "
       close(134)
       
    endif
    
    !BASAL SPRING SHORTEN(+hm)/RELAX(-hm)
    
    hmuch = hmMltplL(2)
    
    if (abs(hmuch).gt.TINYval) then
       open(unit=135,file="BaslSprNC_Shortn.dat")
       
       Nstp  = 5
       stepS = (hmuch)/(Nstp-1)
       
       Initl01 = l0(baslSprNC1)
       Initl02 = l0(baslSprNC2)
       
       sprPair(1) = baslSprNC1 ; sprPair(2) = baslSprNC2
       
       write(135,fmt=*) Initl01,Initl02,"Initl0s"
       write(135,fmt=*) sprPair(1:2),"sprPair"
       
       do i = 1,Nstp
          
          hm = (i-1)*stepS
          call shortening_sprPair(sprPair,hm)
          write(135,fmt=*) l0(baslSprNC1),l0(baslSprNC2),"l0 AftShrtn insd main sb"
          
          call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
          call Equilibrate_system
          call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
          Frame = Frame + 1
          call switchto_NI_model_run_and_switchbackto_TN
          
          if (i.ne.Nstp) then
             l0(baslSprNC1) = Initl01
             l0(baslSprNC2) = Initl02
          endif
          
          
          !lgcl_Rdfn = .False.
          !tol_Rdfn = 0.05d0
          !call get_decision_of_redefining(tol_Rdfn,lgcl_Rdfn)
          
          !if (lgcl_Rdfn.eqv..True.) then
           !  write(135,*) (Frame-1),hm,"FrameNo,hm"
            ! exit
          !elseif (lgcl_Rdfn.eqv..False.) then
           !  continue
          !endif
          
       enddo
       
       write(135,fmt=*) Exprmnt,Frame,"Exprmnt,Frame"
       write(135,fmt=*) " "
       close(135)
       
    endif
       
    !LATERAL SPRING1 SHORTEN(+hm)/RELAX(-hm)
    
    hmuch = hmMltplL(3)
    
    if (abs(hmuch).gt.TINYval) then
       
       open(unit=136,file='LtrlSpr1NC_Shortn.dat')
       open(unit=137,file='HmMltpl_Ltrl.dat',position='append')
       
       Nstp  = 15
       stepS = (hmuch)/(Nstp-1)
       
       Initl01 = l0(ltrlSpr1NC1)
       Initl02 = l0(ltrlSpr1NC2)
       
       sprPair(1) = ltrlSpr1NC1 ; sprPair(2) = ltrlSpr1NC2
       
       write(136,fmt=*) Initl01,Initl02,"Initl0s"
       write(136,fmt=*) sprPair(1:2),"sprPair"
       
       do i = 1,Nstp
          
          hm = (i-1)*stepS
          call shortening_sprPair(sprPair,hm)
          write(136,*) l0(ltrlSpr1NC1),l0(ltrlSpr1NC2),"l0 AftShrtn insd main sb"
          
          call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
          call Equilibrate_system
          call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
          Frame = Frame + 1
          !call switchto_NI_model_run_and_switchbackto_TN
          ProgCycl= ProgCyclNo ; strghtOrNI = 2
          strctV1 = strctVno1  ; strctV2    = strctVno2
          call switchto_NI_model_save_and_switchbackto_TN(ProgCycl,strctV1,strctV2,strghtOrNI,Frame)
          
          !if (i.ne.Nstp) then
          !  l0(ltrlSpr1NC1) = Initl01
          !  l0(ltrlSpr1NC2) = Initl02
          !endif
          
          lgcl_Rdfn = .False.
          tol_Rdfn = 0.05d0
          call get_decision_of_redefining(tol_Rdfn,lgcl_Rdfn)
          
          if (lgcl_Rdfn .eqv. .True.) then
             write(137,*) (Frame-1),hm,"FrameNo,hm"
             hmMltplL(3) = hm
             exit
             
          elseif (lgcl_Rdfn .eqv. .False.) then
             l0(ltrlSpr1NC1) = Initl01
             l0(ltrlSpr1NC2) = Initl02
          endif
          
       enddo
       
       
       if (lgcl_Rdfn .eqv. .False.) then !this if block for lateral only
          
          j=1
          
          do
             hm = (Nstp+j)*stepS
             call shortening_sprPair(sprPair,hm)
             write(136,*) l0(ltrlSpr1NC1),l0(ltrlSpr1NC2),"l0 AftShrtn insd main sb infinite do loop"
             
             call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
             call Equilibrate_system
             call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
             Frame = Frame + 1
             
             ProgCycl= ProgCyclNo  ; strghtOrNI=2
             strctV1 = strctVno1   ; strctV2 = strctVno2
             call switchto_NI_model_save_and_switchbackto_TN(ProgCycl,strctV1,strctV2,strghtOrNI,Frame)
             
             
             !l0(ltrlSpr1NC1) = Initl01
             !l0(ltrlSpr1NC2) = Initl02
             
             lgcl_Rdfn1 = .False.
             tol_Rdfn   = 0.05d0
             call get_decision_of_redefining(tol_Rdfn,lgcl_Rdfn1)
             
             if (lgcl_Rdfn1 .eqv. .True.) then
                write(137,*) (Frame-1),hm,"FrameNo,hm"
                hmMltplL(3) = hm
                exit
             elseif (lgcl_Rdfn1 .eqv. .False.) then
                l0(ltrlSpr1NC1) = Initl01
                l0(ltrlSpr1NC2) = Initl02
             endif
             
             j=j+1
             
          enddo
          
          
       elseif (lgcl_Rdfn .eqv. .True.) then
          continue
          
       endif
       
       
       write(136,fmt=*) Exprmnt,Frame,"Exprmnt,Frame"
       write(136,fmt=*) " "
       close(136)
       close(137)
       
    endif
    
    
  end subroutine Stepwise_SprShorteningNC_and_Equilibrate_UntilMet
  
  subroutine get_first_critical_pair_spr(CP1) !ltrl ones closest to Pulley
    implicit none
    integer, intent(out) :: CP1(1:2)
    integer :: i,j
    integer :: spr_num

    write(*,*) ncl,ncr,N_cell,nsecl,nsecr,"nc1"
    
    do i = 1,N_cell
       
       if (i==ncl) then
          do j=1,nsecl
             spr_num = (i-1)*nsecl+j
             
             if (typ_spr(spr_num)==5) then
                CP1(1) = spr_num
             endif
             
          enddo
          
       elseif (i==(ncl+ncr)) then
          do j=1,nsecr
             spr_num = (i-1)*nsecr+j
             
             if (typ_spr(spr_num)==5) then
                CP1(2) = spr_num
             endif
          enddo
       endif
       
    enddo
    
  end subroutine get_first_critical_pair_spr

  subroutine pulling_or_pushing_the_embryo(Exprmnt,Frame,hmuch)
    implicit none
    integer, intent(in)    :: Exprmnt
    integer, intent(inout) :: Frame   
    real*8, intent(in)     :: hmuch
    
    integer :: m,mmax
    real*8  :: TINY=1.0d-15
    
    mmax = int((hmuch+TINY) / ppS)
    write(*,*) hmuch,ppS,mmax,"hmuch,ppS,mmax"
    
    do m = 2,(mmax+1)
       call get_CgXNode_withStep(m)
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame=Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
    enddo
    
    write(*,*) CgXNode(1),"CgXNode"
    write(*,*) Frame-1,"Frame"
    
  end subroutine pulling_or_pushing_the_embryo
  
  
  subroutine taking_lftrghtPulleyNeigh_downwards(Pulley,LftNeigh,RghtNeigh)
    
    implicit none
    integer, intent(in) :: Pulley,LftNeigh,RghtNeigh
    real*8  :: ydis_frmPulley
    
    call get_Pulley_and_Neighbours(Pulley,LftNeigh,RghtNeigh)
    
    ydis_frmPulley = 0.01d0
    
    call coordntes_to_nodes(coordntes_xy,node_xy)
    
    node_xy(LftNeigh,1)  = node_xy(Pulley,1)
    node_xy(LftNeigh,2)  = node_xy(Pulley,2) - ydis_frmPulley
    
    node_xy(RghtNeigh,1) = node_xy(Pulley,1)
    node_xy(RghtNeigh,2) = node_xy(Pulley,2) - ydis_frmPulley
    
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
  end subroutine taking_lftrghtPulleyNeigh_downwards
  
  subroutine taking_lftrghtOriginNeigh_downwards(LftNeigh,RghtNeigh)
    
    implicit none
    integer, intent(inout) :: LftNeigh,RghtNeigh
    real*8                 :: ydis_frmPulley
    
    write(*,*) LftNeigh,RghtNeigh,"LftNeigh,RghtNeigh,bfr"
    call get_LftRghtNeigh_ofOrigin(LftNeigh,RghtNeigh)
    write(*,*) LftNeigh,RghtNeigh,"LftNeigh,RghtNeigh,aft"
    call sleep(1)
    
    ydis_frmPulley = 0.01d0
    
    call coordntes_to_nodes(coordntes_xy,node_xy)
    
    node_xy(LftNeigh,1)  = origin(1)
    node_xy(LftNeigh,2)  = origin(2) - ydis_frmPulley
    
    node_xy(RghtNeigh,1) = origin(1)
    node_xy(RghtNeigh,2) = origin(2) - ydis_frmPulley
    
    write(*,*) node_xy(LftNeigh,1:2),"node_xy LftNeigh"
    write(*,*) node_xy(RghtNeigh,1:2),"node_xy RghtNeigh"
    
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
  end subroutine taking_lftrghtOriginNeigh_downwards
  
  
  subroutine taking_lftrghtOriginNeigh_downwards_diffWay(LftNeigh,RghtNeigh)
    
    implicit none
    integer, intent(inout) :: LftNeigh,RghtNeigh
    real*8                 :: ydis_frmPulley
    
    LftNeigh  = (Hlf_Ncell+1)*2 - 1 - (CellsMeet)*2
    RghtNeigh = LftNeigh + (Hlf_Ncell+1)*2
    write(*,*) LftNeigh,RghtNeigh,"Lft-Rght Neigh"
    
    ydis_frmPulley = 0.01d0
    
    call coordntes_to_nodes(coordntes_xy,node_xy)
    
    node_xy(LftNeigh,1)  = origin(1)
    node_xy(LftNeigh,2)  = origin(2) - ydis_frmPulley
    
    node_xy(RghtNeigh,1) = origin(1)
    node_xy(RghtNeigh,2) = origin(2) - ydis_frmPulley
    
    write(*,*) node_xy(LftNeigh,1:2),"node_xy LftNeigh"
    write(*,*) node_xy(RghtNeigh,1:2),"node_xy RghtNeigh"
    
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
  end subroutine taking_lftrghtOriginNeigh_downwards_diffWay
  
  
  subroutine Extrapolate_and_adjust_strctProp_ifNeeded_in_MAE(Exprmnt,Frame)!MAE=Manipl_&_Equillibrte
    implicit none
    integer, intent(in)    :: Exprmnt
    integer, intent(inout) :: Frame   
    
    real*8  :: tol_Rdfn
    logical :: lgcl_rdfn
    real*8  :: incrSize
    
    open(unit=89,file='Extrapolate_MAE.dat',position='append')
    
    tol_Rdfn = 0.05d0
    lgcl_rdfn = .False.
    
    call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
    
    if (lgcl_rdfn .eqv. .True.) then
       write(89,fmt=*) lgcl_rdfn,"lgcl_rdfn"
       
    elseif (lgcl_rdfn .eqv. .False.) then
       write(89,fmt=*) lgcl_rdfn,"lgcl_rdfn"
       incrSize = 0.025d0
       
       call allocate_and_initialize_StepVars
       call get_strct_step_in_between_strcts(incrSize)
       
       do
          call Extrapolate_A0l0kaks_for_gettingOverPulley
          
          call Equilibrate_system
          call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
          Frame = Frame+1
          call switchto_NI_model_run_and_switchbackto_TN
          
          call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
          write(89,fmt=*) lgcl_rdfn,(Frame-1),"lgcl_rdfn,Frame"
          
          if (lgcl_rdfn.eqv..True.) then
             call adjst_strctF_Props
             !write(*,*) "stop in flnm: Man&Eq,sb: Extrapolate_and_adjust_strctProp_ifNeeded"
             exit
          endif
          
       enddo
       
       call deallocate_StepVars
       
    endif
    
    close(89)
    
  end subroutine Extrapolate_and_adjust_strctProp_ifNeeded_in_MAE
  
  
  subroutine Extrapolate_and_adjust_strctProp_ifNeeded_in_MAE_woSwitchBack(Exprmnt,Frame)!MAE=Manipl_&_Equillibrte
    implicit none
    integer, intent(in)    :: Exprmnt
    integer, intent(inout) :: Frame   
    
    real*8  :: tol_Rdfn
    logical :: lgcl_rdfn
    real*8  :: incrSize
    
    open(unit=89,file='Extrapolate_MAE.dat',position='append')
    
    tol_Rdfn = 0.05d0
    lgcl_rdfn = .False.
    
    call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
    
    if (lgcl_rdfn .eqv. .True.) then
       write(89,fmt=*) lgcl_rdfn,"lgcl_rdfn"
       
    elseif (lgcl_rdfn .eqv. .False.) then
       
       write(89,fmt=*) lgcl_rdfn,"lgcl_rdfn"
       incrSize = 0.025d0
       
       call allocate_and_initialize_StepVars
       call get_strct_step_in_between_strcts(incrSize)
       
       do
          
          call Extrapolate_A0l0_only_for_gettingOverPulley
          
          call Equilibrate_system
          call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
          Frame = Frame+1
          !call switchto_NI_model_run_and_switchbackto_TN
          
          call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
          write(89,fmt=*) lgcl_rdfn,(Frame-1),"lgcl_rdfn,Frame"
          
          if (lgcl_rdfn.eqv..True.) then
             call adjst_strctF_Props
             !write(*,*) "stop in flnm: Man&Eq,sb: Extrapolate_and_adjust_strctProp_ifNeeded"
             exit
          endif
          
       enddo
       
       call deallocate_StepVars
       
    endif
    
    close(89)
    
  end subroutine Extrapolate_and_adjust_strctProp_ifNeeded_in_MAE_woSwitchBack
  
  
  subroutine Extrapolate_Adjust_strctPropWithCond_in_MAE_woSwitchBack(Exprmnt,Frame)
    implicit none
    integer, intent(in)    :: Exprmnt
    integer, intent(inout) :: Frame   
    
    real*8  :: tol_Rdfn
    logical :: lgcl_rdfn
    real*8  :: incrSize
    
    real*8  :: PresVal(1:N_cell),TnsnVal(1:N_spr)
    logical :: PresLgcl(1:N_cell),TnsnLgcl(1:N_spr)
    logical :: PresNegLgcl(1:N_cell),TnsnPosLgcl(1:N_spr)
    real*8  :: ka_P(1:N_cell),A0_P(1:N_cell),ks_P(1:N_spr),l0_P(1:N_spr)
    
    interface
       
       function TnsnComprsn(dum_Nodes,duml0)
         use system_parameters
         use transfrm_info
         implicit none
         real*8 :: TnsnComprsn(1:N_spr)
         real*8 :: dum_Nodes(1:N_node,1:N_dmnsn)
         real*8 :: duml0(1:N_spr)
       end function TnsnComprsn
       
       function Pressure(dum_Nodes,dumA0)
         use system_parameters
         use transfrm_info
         implicit none
         real*8 :: Pressure(1:N_cell)
         real*8 :: dum_Nodes(1:N_node,1:N_dmnsn)
         real*8 :: dumA0(1:N_cell)
       end function Pressure
       
    end interface
    
    
    open(unit=89,file='Extrapolate_MAE.dat',position='append')
    
    tol_Rdfn = 0.05d0
    lgcl_rdfn = .False.
    
    call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
    
    if (lgcl_rdfn .eqv. .True.) then
       write(89,fmt=*) lgcl_rdfn,"lgcl_rdfn"
       
    elseif (lgcl_rdfn .eqv. .False.) then
       
       PresVal(1:N_cell) = Pressure(node_xy,A0)
       TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
       call get_TnsnPresReverseLgcl(PresVal,TnsnVal,PresNegLgcl,TnsnPosLgcl)
       
       write(89,fmt=*) lgcl_rdfn,"lgcl_rdfn"
       incrSize = 0.025d0
       
       call allocate_and_initialize_StepVars
       call get_strct_step_in_between_strcts(incrSize)
       
       do
          
          call store_PreExtraplt_Props(ka_P,A0_P,ks_P,l0_P)
          
          PresVal(1:N_cell) = Pressure(node_xy,A0)
          TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
          call get_TnsnPresLgcls(PresVal,TnsnVal,PresLgcl,TnsnLgcl,Frame)
          
          !call Extrapolate_A0l0_only_for_gettingOverPulley
          call Extrapolate_A0l0_only_withCondtn_for_gettingOverPulley(PresLgcl,TnsnLgcl)
          
          call Equilibrate_system
          call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
          Frame = Frame+1
          !call switchto_NI_model_run_and_switchbackto_TN
          
          call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
          write(89,fmt=*) lgcl_rdfn,(Frame-1),"lgcl_rdfn,Frame"
          
          if (lgcl_rdfn.eqv..True.) then
             call adjst_strctF_Props
             exit
          endif
          
       enddo
       
       call deallocate_StepVars
       
    endif
    
    close(89)
    
  end subroutine Extrapolate_Adjust_strctPropWithCond_in_MAE_woSwitchBack
  
  
  subroutine revExtraplt_Adjust_strctPropWithCond_in_MAE_woSwitchBack(Exprmnt,Frame)
    implicit none
    integer, intent(in)    :: Exprmnt
    integer, intent(inout) :: Frame   
    
    real*8  :: tol_Rdfn
    logical :: lgcl_rdfn
    real*8  :: incrSize
    
    real*8  :: PresVal(1:N_cell),TnsnVal(1:N_spr)
    logical :: PresLgcl(1:N_cell),TnsnLgcl(1:N_spr)
    logical :: PresNegLgcl(1:N_cell),TnsnPosLgcl(1:N_spr)
    real*8  :: ka_P(1:N_cell),A0_P(1:N_cell),ks_P(1:N_spr),l0_P(1:N_spr)
    
    interface
       
       function TnsnComprsn(dum_Nodes,duml0)
         use system_parameters
         use transfrm_info
         implicit none
         real*8 :: TnsnComprsn(1:N_spr)
         real*8 :: dum_Nodes(1:N_node,1:N_dmnsn)
         real*8 :: duml0(1:N_spr)
       end function TnsnComprsn
       
       function Pressure(dum_Nodes,dumA0)
         use system_parameters
         use transfrm_info
         implicit none
         real*8 :: Pressure(1:N_cell)
         real*8 :: dum_Nodes(1:N_node,1:N_dmnsn)
         real*8 :: dumA0(1:N_cell)
       end function Pressure
       
    end interface
    
    
    open(unit=89,file='rev_Extrapolate_MAE.dat',position='append')
    
    tol_Rdfn = 0.05d0
    lgcl_rdfn = .False.
    
    call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
    
    if (lgcl_rdfn.eqv..False.) then
       write(*,*) "lgcl_rdfn cant be False here, stop"
       stop
    endif
    
    PresVal(1:N_cell) = Pressure(node_xy,A0)
    TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
    call get_TnsnPresReverseLgcl(PresVal,TnsnVal,PresNegLgcl,TnsnPosLgcl)
    
    write(89,fmt=*) lgcl_rdfn,"lgcl_rdfn"
    incrSize = 0.025d0
    
    call allocate_and_initialize_StepVars
    call get_strct_step_in_between_strcts(incrSize)
    
    do
       
       call store_PreExtraplt_Props(ka_P,A0_P,ks_P,l0_P)
       
       PresVal(1:N_cell) = Pressure(node_xy,A0)
       TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
       call get_TnsnPresLgcls(PresVal,TnsnVal,PresLgcl,TnsnLgcl,Frame)
       
       call reverse_Extrapolate_A0l0_only_withCondtn_for_gettingOverPulley(PresLgcl,TnsnLgcl)
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       
       call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
       write(89,*) lgcl_rdfn,(Frame-1),"lgcl_rdfn,Frame"
       
       if (lgcl_rdfn.eqv..False.) then
          
          k_area(1:N_cell) = ka_P(1:N_cell)
          A0(1:N_cell)     = A0_P(1:N_cell)
          k_spr(1:N_spr)   = ks_P(1:N_spr)
          l0(1:N_spr)      = l0_P(1:N_spr)
          
          call adjst_strctF_Props
          exit
          
       endif
       
    enddo
    
    call deallocate_StepVars
    
    close(89)
    
  end subroutine revExtraplt_Adjust_strctPropWithCond_in_MAE_woSwitchBack
  
  
  
  subroutine allocate_and_initialize_StepVars
    implicit none
    
    allocate(A0_stepF(1:N_cell),ka_stepF(1:N_cell))
    allocate(l0_stepF(1:N_spr) ,ks_stepF(1:N_spr))
    
    A0_stepF = -1.0d30 ; ka_stepF = -1.0d30
    l0_stepF = -1.0d30 ; ks_stepF = -1.0d30
    
  end subroutine allocate_and_initialize_StepVars
  
  
  subroutine deallocate_StepVars
    implicit none
    
    deallocate(A0_stepF,ka_stepF)
    deallocate(l0_stepF,ks_stepF)
    
  end subroutine deallocate_StepVars
  
  
  subroutine get_strct_step_in_between_strcts(incrSize)
    implicit none
    real*8, intent(in) :: incrSize !can also be done for individ param incrSize
    
    integer :: i,j,jmax
    
    open(unit=193,file='StepValues.dat')
    
    do i = 1,2
       
       if (i==1) jmax = N_cell
       if (i==2) jmax = N_spr
       
       do j = 1,jmax
          
          if (i==1) then
             A0_stepF(j) = incrSize*(A0_F(j)-A0_I(j))  
             ka_stepF(j) = incrSize*(ka_F(j)-ka_I(j))
             write(193,*) A0_stepF(j),ka_stepF(j),j,"A0 and ka stepF,area nm"
             
          elseif (i==2) then
             l0_stepF(j) = incrSize*(l0_F(j)-l0_I(j))  
             ks_stepF(j) = incrSize*(ks_F(j)-ks_I(j))
             write(193,*) l0_stepF(j),ks_stepF(j),j,"l0 and ks stepF,area nm"
             
          endif
          
       enddo
       
    enddo
    
    close(193)
    
  end subroutine get_strct_step_in_between_strcts
  
  subroutine Extrapolate_A0l0kaks_for_gettingOverPulley
    implicit none
    integer :: i,j,jmax
    real*8  :: Tiny
    
    Tiny = -0.0d0
    
    do i = 1,2
       if (i==1) jmax=N_cell
       if (i==2) jmax=N_spr
       
       do j = 1,jmax
          
          if (i==1) then
             
             A0(j)     = A0(j) + A0_stepF(j)
             k_area(j) = k_area(j) + ka_stepF(j)
             
             if ((A0(j).lt.Tiny).or.(k_area(j).lt.Tiny)) then
                write(*,*) "A0/k_area got_negative"
                write(*,*) A0(j),k_area(j),j,"A0/k_area/j"
                stop
             endif
             
          elseif (i==2) then
             
             l0(j)     = l0(j) + l0_stepF(j)
             k_spr(j)  = k_spr(j) + ks_stepF(j)
             
             if ((l0(j).lt.Tiny).or.(k_spr(j).lt.Tiny)) then
                write(*,*) "l0/k_spr got_negative"
                write(*,*) l0(j),k_spr(j),j,"l0/k_spr/j"
                stop
             endif
             
          endif
       enddo
       
    enddo
    
  end subroutine Extrapolate_A0l0kaks_for_gettingOverPulley
  
  
  subroutine Extrapolate_A0l0_only_for_gettingOverPulley
    implicit none
    integer :: i,j,jmax
    real*8  :: Tiny
    real*8  :: A0_stre(1:N_cell),l0_stre(1:N_spr)
    real*8  :: A0_diff,l0_diff
    
    Tiny = -0.0d0
    
    A0_stre(1:N_cell) = A0(1:N_cell) 
    l0_stre(1:N_spr)  = l0(1:N_spr)
    
    do i = 1,2
       
       if (i==1) jmax = N_cell
       if (i==2) jmax = N_spr
       
       do j = 1,jmax
          
          if (i==1) then
             
             A0(j)   = A0(j) + A0_stepF(j)
             
             if ((A0(j).lt.Tiny).or.(k_area(j).lt.Tiny)) then
                write(*,*) "A0/k_area got_negative"
                write(*,*) A0(j),k_area(j),j,"A0/k_area/j"
                
                A0(j) = A0_stre(j)
             endif
             
          elseif (i==2) then
             
             l0(j) = l0(j) + l0_stepF(j)
             
             if ((l0(j).lt.Tiny).or.(k_spr(j).lt.Tiny)) then
                write(*,*) "l0/k_spr got_negative"
                write(*,*) l0(j),k_spr(j),j,"l0/k_spr/j"
                
                l0(j) = l0_stre(j)
             endif
             
          endif
       enddo
       
    enddo
    
  end subroutine Extrapolate_A0l0_only_for_gettingOverPulley
  
  
  subroutine Extrapolate_A0l0_only_withCondtn_for_gettingOverPulley(PresLgcl,TnsnLgcl)
    implicit none
    logical, intent(in) :: PresLgcl(1:N_cell)
    logical, intent(in) :: TnsnLgcl(1:N_spr)

    real*8  :: A0_stre(1:N_cell),l0_stre(1:N_spr)
    integer :: i,j,jmax
    real*8  :: Tiny
    
    Tiny = -0.0d0
    
    A0_stre(1:N_cell) = A0(1:N_cell) 
    l0_stre(1:N_spr)  = l0(1:N_spr)
    
    
    do i = 1,2
       
       if (i==1) jmax=N_cell
       if (i==2) jmax=N_spr
       
       do j = 1,jmax
          
          if (i==1) then
             
             if (PresLgcl(j) .eqv. .True.) A0(j) = A0(j) + A0_stepF(j)
             
             if (A0(j).lt.Tiny) then
                write(*,*) "A0 got_negative"
                write(*,*) A0(j),k_area(j),j,"A0/k_area/j"
                
                A0(j) = A0_stre(j)
             endif
             
          elseif (i==2) then
             
             if (TnsnLgcl(j) .eqv. .True.) l0(j) = l0(j) + l0_stepF(j)
             
             if (l0(j).lt.Tiny) then
                write(*,*) "l0 got_negative"
                write(*,*) l0(j),k_spr(j),j,"l0/k_spr/j"
                
                l0(j) = l0_stre(j)
             endif
             
          endif
       enddo
       
    enddo
    
  end subroutine Extrapolate_A0l0_only_withCondtn_for_gettingOverPulley
  
  
  subroutine reverse_Extrapolate_A0l0_only_withCondtn_for_gettingOverPulley(PresLgcl,TnsnLgcl)
    implicit none
    logical, intent(in) :: PresLgcl(1:N_cell)
    logical, intent(in) :: TnsnLgcl(1:N_spr)
    
    integer :: i,j,jmax
    real*8  :: Tiny
    
    Tiny = -0.0d0
    
    do i = 1,2
       
       if (i==1) jmax=N_cell
       if (i==2) jmax=N_spr
       
       do j = 1,jmax
          
          if (i==1) then
             
             if (PresLgcl(j) .eqv. .True.) A0(j) = A0(j) - A0_stepF(j)
             
             if (A0(j).lt.Tiny) then
                write(*,*) "A0 got_negative"
                write(*,*) A0(j),k_area(j),j,"A0/k_area/j"
                stop
             endif
             
          elseif (i==2) then
             
             if (TnsnLgcl(j) .eqv. .True.) l0(j) = l0(j) - l0_stepF(j)
             
             if (l0(j).lt.Tiny) then
                write(*,*) "l0 got_negative"
                write(*,*) l0(j),k_spr(j),j,"l0/k_spr/j"
                stop
             endif
             
          endif
       enddo
       
    enddo
    
  end subroutine reverse_Extrapolate_A0l0_only_withCondtn_for_gettingOverPulley
  
  
  
  subroutine adjst_strctF_Props
    implicit none
    integer :: i,j,jmax
    
    do i = 1,2
       
       if (i==1) jmax=N_cell
       if (i==2) jmax=N_spr
       
       do j = 1,jmax
          if (i==1) then
             A0_F(j) = A0(j)
             ka_F(j) = k_area(j)
          elseif (i==2) then
             l0_F(j) = l0(j)
             ks_F(j) = k_spr(j)
          endif
       enddo
       
    enddo
    
  end subroutine adjst_strctF_Props
  
  
  
  subroutine pulling_or_pushing_forces_at_CS(strctToRead,cntrlSt,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: strctToRead,cntrlSt,ExpNo
    integer, intent(inout) :: FrmNo
    
    integer :: RghtstrtTN
    real*8  :: E
    integer :: saveWhom
    integer :: ExtrpltUsingWhatProps
    
    real*8, allocatable :: PresVal(:),TnsnVal(:)
    
    interface
       
       function TnsnComprsn(dum_Nodes,duml0)
         use system_parameters
         use transfrm_info
         implicit none
         real*8 :: TnsnComprsn(1:N_spr)
         real*8 :: dum_Nodes(1:N_node,1:N_dmnsn)
         real*8 :: duml0(1:N_spr)
       end function TnsnComprsn
       
       function Pressure(dum_Nodes,dumA0)
         use system_parameters
         use transfrm_info
         implicit none
         real*8 :: Pressure(1:N_cell)
         real*8 :: dum_Nodes(1:N_node,1:N_dmnsn)
         real*8 :: dumA0(1:N_cell)
       end function Pressure
       
    end interface
    
    open(unit=64,file='PullingPushingAtNI.dat')
    
    if (modelID==1) then
       call switch_to_NI_model
       E = Energy(node_xy,l0,A0)
       call get_gradient(node_xy,l0,A0,gradient)   
    elseif (modelID==2) then
       continue
    endif
    
    allocate(PresVal(1:N_cell),TnsnVal(1:N_spr))
    call read_strctProps_withAddedCell(strctToRead)
    
    A0(1:N_cell)     = A0_strctNI(strctToRead,1:N_cell)
    k_area(1:N_cell) = ka_strctNI(strctToRead,1:N_cell)
    l0(1:N_spr)      = l0_strctNI(strctToRead,1:N_spr)
    k_spr(1:N_spr)   = ks_strctNI(strctToRead,1:N_spr)
    CgXNode(1:N_node) = CgX_strctNI(strctToRead,1:N_node)
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    
    RghtstrtTN = (Hlf_Ncell+1)*2 + 1
    CgXNode(1:2) = +0.25d0 ; CgXNode(RghtstrtTN:(RghtstrtTN+1)) = -0.25d0 !Pull
    !CgXNode(1:2) = -0.25d0 ; CgXNode(RghtstrtTN:(RghtstrtTN+1)) = +0.25d0 !Push
    
    write(*,*) node_xy(1,1:2),node_xy(2,1:2),"node_xy 1-2"
    write(*,*) node_xy(RghtstrtTN,1:2),node_xy(RghtstrtTN+1,1:2),"node_xy R 1-2"
    
    if (cntrlSt == 1) then
       
       write(*,*) cntrlSt,(FrmNo-1),"cntrlSt1,FrmNo"
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       
    elseif (cntrlSt == 2) then
       
       write(*,*) cntrlSt,(FrmNo-1),"cntrlSt2,FrmNo"
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       
       allocate(A0_I(1:N_cell),ka_I(1:N_cell),l0_I(1:N_spr),ks_I(1:N_spr))
       allocate(A0_F(1:N_cell),ka_F(1:N_cell),l0_F(1:N_spr),ks_F(1:N_spr))
       
       A0_I(1:N_cell)=1.0d30 ; ka_I(1:N_cell)=1.0d30 ; l0_I(1:N_spr)=1.0d30 ; ks_I(1:N_spr)=1.0d30
       A0_F(1:N_cell)=1.0d30 ; ka_F(1:N_cell)=1.0d30 ; l0_F(1:N_spr)=1.0d30 ; ks_F(1:N_spr)=1.0d30
       
       allocate(A0I_NI(1:N_cell),kaI_NI(1:N_cell),l0I_NI(1:N_spr),ksI_NI(1:N_spr))
       allocate(A0F_NI(1:N_cell),kaF_NI(1:N_cell),l0F_NI(1:N_spr),ksF_NI(1:N_spr))
       
       A0I_NI(1:N_cell)=1.0d30 ; kaI_NI(1:N_cell)=1.0d30 ; l0I_NI(1:N_spr)=1.0d30 ; ksI_NI(1:N_spr)=1.0d30
       A0F_NI(1:N_cell)=1.0d30 ; kaF_NI(1:N_cell)=1.0d30 ; l0F_NI(1:N_spr)=1.0d30 ; ksF_NI(1:N_spr)=1.0d30
       
       A0_F(1:N_cell) = A0_strctNI(strctToRead,1:N_cell) ; A0_I(1:N_cell) = A0_strctNI((strctToRead-1),1:N_cell)
       ka_F(1:N_cell) = ka_strctNI(strctToRead,1:N_cell) ; ka_I(1:N_cell) = ka_strctNI((strctToRead-1),1:N_cell)
       l0_F(1:N_spr)  = l0_strctNI(strctToRead,1:N_spr)  ; l0_I(1:N_spr)  = l0_strctNI((strctToRead-1),1:N_spr)
       ks_F(1:N_spr)  = ks_strctNI(strctToRead,1:N_spr)  ; ks_I(1:N_spr)  = ks_strctNI((strctToRead-1),1:N_spr)
       
       write(*,*) A0_F(5),A0_I(5),ka_F(5),ka_I(5),"A0F+kaF in cell 5"
       write(*,*) l0_F(9),l0_I(9),ks_F(9),ks_I(9),"l0F+ksF in spr  9"
       
       call sleep(5)
       !call Extrapolate_and_adjust_strctProp_ifNeeded_in_MAE_woSwitchBack(ExpNo,FrmNo)
       call Extrapolate_Adjust_strctPropWithCond_in_MAE_woSwitchBack(ExpNo,FrmNo)
       
       A0I_NI(1:N_cell) = A0_I(1:N_cell) ; kaI_NI(1:N_cell) = ka_I(1:N_cell)
       l0I_NI(1:N_spr)  = l0_I(1:N_spr)  ; ksI_NI(1:N_spr)  = ks_I(1:N_spr)
       
       A0F_NI(1:N_cell) = A0_F(1:N_cell) ; kaF_NI(1:N_cell) = ka_F(1:N_cell)
       l0F_NI(1:N_spr)  = l0_F(1:N_spr)  ; ksF_NI(1:N_spr)  = ks_F(1:N_spr)
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       
    elseif (cntrlSt==3) then
       
       write(*,*) cntrlSt,(FrmNo-1),"cntrlSt3,FrmNo"
       
       A0(1:N_cell) = A0F_NI(1:N_cell)   ; k_area(1:N_cell) = kaF_NI(1:N_cell)
       l0(1:N_spr)  = l0F_NI(1:N_spr)    ; k_spr(1:N_spr)   = ksF_NI(1:N_spr)
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       
    elseif (cntrlSt==4) then
       
       write(*,*) cntrlSt,(FrmNo-1),"cntrlSt4,FrmNo"
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       
    endif
    
    !call switchto_NI_model_run_and_switchbackto_TN
    
    !saveWhom = 3 ! not saving
    !call saveA0l0_ofStrct_AddedCell_TNorNI(strctToRead,saveWhom)
    
    if (modelID==1) then
       continue
    elseif (modelID==2) then
       call deallocate_repetitive_arrays
       call switchback_to_TN_model
    endif
    
    close(64)
    
  end subroutine pulling_or_pushing_forces_at_CS
  
  
  subroutine getting_strctPropsBack_pulling_or_pushing_forces_at_CS(strctToRead,cntrlSt,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: strctToRead,cntrlSt,ExpNo
    integer, intent(inout) :: FrmNo
    
    integer :: RghtstrtTN
    real*8  :: E
    integer :: saveWhom
    
    real*8, allocatable :: PresVal(:),TnsnVal(:)
    
    interface
       
       function TnsnComprsn(dum_Nodes,duml0)
         use system_parameters
         use transfrm_info
         implicit none
         real*8 :: TnsnComprsn(1:N_spr)
         real*8 :: dum_Nodes(1:N_node,1:N_dmnsn)
         real*8 :: duml0(1:N_spr)
       end function TnsnComprsn
       
       function Pressure(dum_Nodes,dumA0)
         use system_parameters
         use transfrm_info
         implicit none
         real*8 :: Pressure(1:N_cell)
         real*8 :: dum_Nodes(1:N_node,1:N_dmnsn)
         real*8 :: dumA0(1:N_cell)
       end function Pressure
       
    end interface
    
    open(unit=64,file='PullingPushingAtNI.dat')
    
    if (modelID==1) then
       call switch_to_NI_model
       E = Energy(node_xy,l0,A0)
       call get_gradient(node_xy,l0,A0,gradient)   
    elseif (modelID==2) then
       continue
    endif
    
    allocate(PresVal(1:N_cell),TnsnVal(1:N_spr))
    call read_strctProps_withAddedCell(strctToRead)
    
    A0(1:N_cell)     = A0_strctNI(strctToRead,1:N_cell)
    k_area(1:N_cell) = ka_strctNI(strctToRead,1:N_cell)
    l0(1:N_spr)      = l0_strctNI(strctToRead,1:N_spr)
    k_spr(1:N_spr)   = ks_strctNI(strctToRead,1:N_spr)
    CgXNode(1:N_node) = CgX_strctNI(strctToRead,1:N_node)
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    
    RghtstrtTN = (Hlf_Ncell+1)*2 + 1
    CgXNode(1:2) = +0.00d0 ; CgXNode(RghtstrtTN:(RghtstrtTN+1)) = -0.00d0
    
    if (cntrlSt == 1) then
       
       write(*,*) cntrlSt,(FrmNo-1),"cntrlSt1,FrmNo"
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       
    elseif (cntrlSt == 2) then
       
       write(*,*) cntrlSt,(FrmNo-1),"cntrlSt2,FrmNo"
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       
       allocate(A0_I(1:N_cell),ka_I(1:N_cell),l0_I(1:N_spr),ks_I(1:N_spr))
       allocate(A0_F(1:N_cell),ka_F(1:N_cell),l0_F(1:N_spr),ks_F(1:N_spr))
       
       A0_I(1:N_cell)=1.0d30 ; ka_I(1:N_cell)=1.0d30 ; l0_I(1:N_spr)=1.0d30 ; ks_I(1:N_spr)=1.0d30
       A0_F(1:N_cell)=1.0d30 ; ka_F(1:N_cell)=1.0d30 ; l0_F(1:N_spr)=1.0d30 ; ks_F(1:N_spr)=1.0d30
       
       allocate(A0I_NI(1:N_cell),kaI_NI(1:N_cell),l0I_NI(1:N_spr),ksI_NI(1:N_spr))
       allocate(A0F_NI(1:N_cell),kaF_NI(1:N_cell),l0F_NI(1:N_spr),ksF_NI(1:N_spr))
       
       A0I_NI(1:N_cell)=1.0d30 ; kaI_NI(1:N_cell)=1.0d30 ; l0I_NI(1:N_spr)=1.0d30 ; ksI_NI(1:N_spr)=1.0d30
       A0F_NI(1:N_cell)=1.0d30 ; kaF_NI(1:N_cell)=1.0d30 ; l0F_NI(1:N_spr)=1.0d30 ; ksF_NI(1:N_spr)=1.0d30
       
       A0_F(1:N_cell) = A0_strctNI(strctToRead,1:N_cell) ; A0_I(1:N_cell) = A0_strctNI((strctToRead-1),1:N_cell)
       ka_F(1:N_cell) = ka_strctNI(strctToRead,1:N_cell) ; ka_I(1:N_cell) = ka_strctNI((strctToRead-1),1:N_cell)
       l0_F(1:N_spr)  = l0_strctNI(strctToRead,1:N_spr)  ; l0_I(1:N_spr)  = l0_strctNI((strctToRead-1),1:N_spr)
       ks_F(1:N_spr)  = ks_strctNI(strctToRead,1:N_spr)  ; ks_I(1:N_spr)  = ks_strctNI((strctToRead-1),1:N_spr)
       
       write(*,*) A0_F(5),A0_I(5),ka_F(5),ka_I(5),"A0F+kaF in cell 5"
       write(*,*) l0_F(9),l0_I(9),ks_F(9),ks_I(9),"l0F+ksF in spr  9"
       
       call sleep(5)
       call revExtraplt_Adjust_strctPropWithCond_in_MAE_woSwitchBack(ExpNo,FrmNo)
       
       A0I_NI(1:N_cell) = A0_I(1:N_cell) ; kaI_NI(1:N_cell) = ka_I(1:N_cell)
       l0I_NI(1:N_spr)  = l0_I(1:N_spr)  ; ksI_NI(1:N_spr)  = ks_I(1:N_spr)
       
       A0F_NI(1:N_cell) = A0_F(1:N_cell) ; kaF_NI(1:N_cell) = ka_F(1:N_cell)
       l0F_NI(1:N_spr)  = l0_F(1:N_spr)  ; ksF_NI(1:N_spr)  = ks_F(1:N_spr)
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       
    elseif (cntrlSt==3) then
       
       write(*,*) cntrlSt,(FrmNo-1),"cntrlSt3,FrmNo"
       
       A0(1:N_cell) = A0F_NI(1:N_cell)   ; k_area(1:N_cell) = kaF_NI(1:N_cell)
       l0(1:N_spr)  = l0F_NI(1:N_spr)    ; k_spr(1:N_spr)   = ksF_NI(1:N_spr)
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       
    elseif (cntrlSt==4) then
       
       write(*,*) cntrlSt,(FrmNo-1),"cntrlSt4,FrmNo"
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       
    endif
    
    !call switchto_NI_model_run_and_switchbackto_TN
    
    !saveWhom = 3 ! not saving
    !call saveA0l0_ofStrct_AddedCell_TNorNI(strctToRead,saveWhom)
    
    if (modelID==1) then
       continue
    elseif (modelID==2) then
       call deallocate_repetitive_arrays
       call switchback_to_TN_model
    endif
    
    close(64)
    
  end subroutine getting_strctPropsBack_pulling_or_pushing_forces_at_CS
  
  
   
  subroutine store_PreExtraplt_Props(ka_P,A0_P,ks_P,l0_P)
    implicit none
    real*8, intent(out) :: ka_P(1:N_cell),A0_P(1:N_cell)
    real*8, intent(out) :: ks_P(1:N_spr) ,l0_P(1:N_spr)
    
    ka_P(1:N_cell) = k_area(1:N_cell)
    A0_P(1:N_cell) = A0(1:N_cell)
    ks_P(1:N_spr)  = k_spr(1:N_spr)
    l0_P(1:N_spr)  = l0(1:N_spr)
    
  end subroutine store_PreExtraplt_Props
  
  subroutine get_TnsnPresReverseLgcl(PresVal,TnsnVal,PresNegLgcl,TnsnPosLgcl)
    implicit none
    real*8, intent(in)  :: PresVal(1:N_cell)
    real*8, intent(in)  :: TnsnVal(1:N_spr)
    logical,intent(out) :: PresNegLgcl(1:N_cell)
    logical,intent(out) :: TnsnPosLgcl(1:N_spr)
    
    integer :: i,j,jmax
    real*8  :: ZERO
    integer :: cnt_PresNeg,cnt_TnsnPos
    
    PresNegLgcl(1:N_cell) = .False.
    TnsnPosLgcl(1:N_spr)  = .False.
    
    ZERO = 0.00d0
    cnt_PresNeg=0 ; cnt_TnsnPos=0
    
    do i = 1,2
       
       if (i==1) jmax = N_cell
       if (i==2) jmax = N_spr
       
       do j = 1,jmax
          
          if (i==1) then
             
             if (PresVal(j) .lt. ZERO) then
                PresNegLgcl(j) = .True.
                cnt_PresNeg = cnt_PresNeg+1
             endif
             
          elseif (i==2) then
             
             if (TnsnVal(j) .gt. ZERO) then
                TnsnPosLgcl(j) = .True.
                cnt_TnsnPos = cnt_TnsnPos+1
             endif
             
          endif
             
       enddo
       
    enddo
    
    if ((cnt_PresNeg.ne.0) .or. (cnt_TnsnPos.ne.0)) then
       
       open(unit=76,file='PresTnsnReversibility.dat')
       
       do i = 1,2
          
          if (i==1) jmax=N_cell
          if (i==2) jmax=N_spr
          
          do j = 1,jmax
             if (i==1) write(76,*) PresNegLgcl(j),j
             if (i==2) write(76,*) TnsnPosLgcl(j),j
          enddo
          
          write(76,*) " "
          
       enddo
       
       close(76)
       
       write(*,*) "Pressure Negative or Tension Positive" 
       !stop
       
    endif
    
  end subroutine get_TnsnPresReverseLgcl
  
  subroutine get_TnsnPresLgcls(PresVal,TnsnVal,PresLgcl,TnsnLgcl,Frame)
    implicit none
    real*8,  intent(in)  :: PresVal(1:N_cell)
    real*8,  intent(in)  :: TnsnVal(1:N_spr)
    logical, intent(out) :: PresLgcl(1:N_cell)
    logical, intent(out) :: TnsnLgcl(1:N_spr)
    integer, intent(in)  :: Frame
    
    real*8  :: tol_Pres,tol_Tnsn
    integer :: i,j,jmax,p,pmax
    integer :: nsprsInACell
    integer :: N_itm
    logical :: lgcl_p
    
    integer, allocatable :: Apcls(:),Bsals(:),Ltrls(:)
    integer, allocatable :: diagSprL(:),diagSprR(:)
    
    integer :: spclAttnCond,cellL,cellR
    
    tol_Pres = 0.005d0
    tol_Tnsn = 0.005d0
    
    PresLgcl(1:N_cell) = .True.
    TnsnLgcl(1:N_spr)  = .True.
    
    open(unit=64,file='TnsnPresLgcls.dat',position='append')
    
    do i = 1,2
       
       if (i==1) jmax = N_cell
       if (i==2) jmax = N_spr
       
       do j = 1,jmax
          
          if (i==1) then
             
             if (abs(PresVal(j)) .le. tol_Pres) then
                PresLgcl(j) = .False.
                
             elseif (abs(PresVal(j)) .gt. tol_Pres) then
                PresLgcl(j) = .True.
                
             endif
             
          elseif (i==2) then
             
             if (abs(TnsnVal(j)) .le. tol_Tnsn) then
                TnsnLgcl(j) = .False.
                
             elseif (TnsnVal(j) .gt. tol_Tnsn) then
                TnsnLgcl(j) = .True.
             endif
             
          endif
             
       enddo
       
    enddo
    
    if (modelID==1) continue
    
    if (modelID==2) then

       nsprsInACell = (NAEC_Apcl+1) + (NAEC_Bsal+1) + (NAEC_Ltrl+1)
       write(*,*) "nsprsInACell is =", nsprsInACell
       
       allocate(Apcls(1:(NAEC_Apcl+1)),Bsals(1:(NAEC_Bsal+1)),Ltrls(1:(NAEC_Ltrl+1)))
       N_itm = 3
       
       do i = 1,(N_cell-1)
          
          call get_the_apcls_bsals_and_oneltrls(i,Apcls,Bsals,Ltrls,NAEC_Apcl,NAEC_Bsal,NAEC_Ltrl) 
          
          do j = 1,N_itm
             
             if (j==1) pmax = NAEC_Apcl+1
             if (j==2) pmax = NAEC_Bsal+1
             if (j==3) pmax = NAEC_Ltrl+1

             lgcl_p = .True.
             
             do p = 1,pmax
                
                if (j==1) lgcl_p = TnsnLgcl(Apcls(p))
                if (j==2) lgcl_p = TnsnLgcl(Bsals(p))
                if (j==3) lgcl_p = TnsnLgcl(Ltrls(p))
                
                if (lgcl_p .eqv. .False.) then
                   
                   if (j==1) TnsnLgcl(Apcls(1) : Apcls(NAEC_Apcl+1)) = .False.
                   if (j==2) TnsnLgcl(Bsals(1) : Bsals(NAEC_Bsal+1)) = .False.
                   if (j==3) TnsnLgcl(Ltrls(1) : Ltrls(NAEC_Ltrl+1)) = .False.
                   
                   lgcl_p = .True.
                   exit
                endif
                   
             enddo
             
          enddo
          
       enddo
       
       Ltrls(1) = (N_cell-1)*(nsprsInACell)+1 ; Ltrls(2) = Ltrls(1)+1 ; Ltrls(3) = Ltrls(1)+2
       pmax = NAEC_Ltrl+1
       lgcl_p = .True.
       
       do p = 1,pmax
          
          lgcl_p = TnsnLgcl(Ltrls(p))
          
          if (lgcl_p .eqv. .False.) then
             
             TnsnLgcl(Ltrls(1) : Ltrls(NAEC_Ltrl+1)) = .False.
             lgcl_p = .True.
             exit
          endif
          
       enddo
       
       
       spclAttnCond = 1
       
       if (spclAttnCond == 0) then
          continue 
       elseif (spclAttnCond == 1) then
          
          cellL = invgnting_cells(1) ; cellR = invgnting_cells(2)
          if (cellL .ne. 7) write(*,*) "NOT 7"
          
          allocate(diagSprL(1:(NAEC_Ltrl+1)),diagSprR(1:(NAEC_Ltrl+1)))
          
          do i = 1,(NAEC_Ltrl+1)
             diagSprL(i) = (cellL-1)*(nsprsInACell) + (NAEC_Apcl+1) + (NAEC_Bsal+1) + i
             diagSprR(i) = (cellR-1)*(nsprsInACell) + (NAEC_Apcl+1) + (NAEC_Bsal+1) + i
          enddo
          
          do i = 1,(NAEC_Ltrl+1)
             TnsnLgcl(diagSprL(i)) = .False.
             TnsnLgcl(diagSprR(i)) = .False.
          enddo
          
       elseif (spclAttnCond == 2) then
          continue
       endif
       
       
    endif

    write(64,*) Frame,"Frame"
    
    
    do i = 1,2
       
       if (i==1) jmax=N_cell
       if (i==2) jmax=N_spr
       
       do j = 1,jmax
          if (i==1) write(64,*) PresLgcl(j),j,"Pres"
          if (i==2) write(64,*) TnsnLgcl(j),j,"Tnsn"
       enddo
       
    enddo
    
    write(64,*) " "
    
    close(64)
    
  end subroutine get_TnsnPresLgcls
  
  
  
end module Manipulate_and_equilibrate


