
module spr_PenFunc
  use system_parameters
  use cell_info
  use transfrm_info
  
  implicit none
  
contains
  
  real*8 function spr_PenF(dum_Nodes) !dummy, N_spring = N_spr
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    integer :: i
    character(len=100) :: flnm
    
    spr_PenF = 0.0d0
    
    do i = 1,N_spr
       
       write(flnm,'(a09,i2.2,a10)')"PenF_Spr",i,"_vs_l0.dat" 
       open(unit=41,file=adjustl(trim(flnm)))
       
       spr_PenF = spr_PenF + indvd_spr_PenF(i)
       
       write(unit=41,fmt=*) indvd_spr_PenF(i)
       !write(*,*) spr_PenF,"for spring i =",i
       close(41)
       
    enddo
    !write(*,*) spr_PenF,"sprPenF"

  contains
    
    real*8 function indvd_spr_PenF(spr_nmbr)
      implicit none
      integer :: spr_nmbr
      
      !if (spr_nmbr.eq.18) write(*,*) dflctn(spr_nmbr),"dflctn"
      indvd_spr_PenF = 0.5d0 * alpha(spr_nmbr) * (dflctn(spr_nmbr))**2

      !if (spr_nmbr==1) write(*,*) dflctn(spr_nmbr)
      
    end function indvd_spr_PenF
    
    real*8 function dflctn(spr_nmbr)
      implicit none
      integer :: spr_nmbr
      integer :: Node1,Node2,Node3
      integer :: pulley_typ
      
      pulley_typ = 6
      
      !if (typ_spr(spr_nmbr) .ne. pulley_typ) then
      if (spr_node(spr_nmbr,0) == 2) then
         
         Node1 = spr_node(spr_nmbr,1)
         Node2 = spr_node(spr_nmbr,2)
         !write(*,*) Node1,Node2,"Node1,Node2"
         l(spr_nmbr) = lngth_u_two_nodes(Node1,Node2)
         dflctn      = l(spr_nmbr) - Lt(spr_nmbr)

         !write(*,*) l(spr_nmbr),Lt(spr_nmbr),dflctn,"D1"
         !stop
         !if (spr_nmbr.eq.18) write(*,*) dflctn,l(spr_nmbr),Lt(spr_nmbr),"dflctn,l,Lt for spring i =",spr_nmbr
         
      elseif (spr_node(spr_nmbr,0) == 3) then
         Node1 = spr_node(spr_nmbr,1)
         Node2 = spr_node(spr_nmbr,2)
         Node3 = spr_node(spr_nmbr,3)
         
         l(spr_nmbr) = lngth_u_three_nodes(Node1,Node2,Node3)
         dflctn = l(spr_nmbr) - Lt(spr_nmbr)
         !write(*,*) l(spr_nmbr),Lt(spr_nmbr),dflctn,"D2"
         
      endif
      
    end function dflctn
    
    real*8 function lngth_u_two_nodes(Node1,Node2)
      implicit none
      integer :: Node1,Node2
      real*8  :: N1(1:N_dmnsn),N2(1:N_dmnsn)
      real*8  :: res_sq
      integer :: i

      N1 = 0.0d0 ; N2 = 0.0d0
      
      N1(1:N_dmnsn) = dum_Nodes(Node1,1:N_dmnsn)
      N2(1:N_dmnsn) = dum_Nodes(Node2,1:N_dmnsn)
      !write(*,*) N1(1:2),"N1"
      !write(*,*) N2(1:2),"N2"
      
      lngth_u_two_nodes = 0.0d0
      res_sq = 0.0d0

      !open(unit=137,file='N2_N1.dat',position='append')
      !write(unit=137,fmt=*) N2,"N2"
      !write(unit=137,fmt=*) N1,"N1"
      !close(137)
      
      do i = 1,2
         res_sq = res_sq + (N2(i) - N1(i))**2 
      enddo
      
      lngth_u_two_nodes = sqrt(res_sq) 
      !write(*,*) lngth_u_two_nodes,"lutn"

    end function lngth_u_two_nodes

    real*8 function lngth_u_three_nodes(Node1,Node2,Node3)
      implicit none
      integer :: Node1,Node2,Node3
      real*8  :: N1(1:N_dmnsn),N2(1:N_dmnsn),N3(1:N_dmnsn)
      integer :: prtn,i
      real*8  :: res_sq
      
      N1 = dum_Nodes(Node1,1:N_dmnsn)
      N2 = dum_Nodes(Node2,1:N_dmnsn)
      N3 = dum_Nodes(Node3,1:N_dmnsn)
      
      lngth_u_three_nodes = 0.0d0
      !res_sq = 0.0d0
      
      do prtn = 1,2
         res_sq = 0.0d0
         do i = 1,N_dmnsn
            if(prtn.eq.1) then
               res_sq = res_sq + (N2(i) - N1(i))**2
            else
               res_sq = res_sq + (N3(i) - N2(i))**2
            endif
         enddo
         lngth_u_three_nodes = lngth_u_three_nodes + sqrt(res_sq)
      enddo
      
    end function lngth_u_three_nodes
    
  end function spr_PenF
  
end module spr_PenFunc





module diag_PenFunc
  use system_parameters
  use cell_info
  use transfrm_info
  use triangle_routines
  
contains

  ! real*8 function diag_PenF()
  !   implicit none
  !   real*8  :: dum_Nodes(1:N_node,1:N_dmnsn) 
  !   integer :: i

  !   diag_PenF = 0.0d0

  !   do i = 1,N_cell
  !      if (area_sides(i) .eq. 4) then
  !         diag_PenF = diag_PenF + indvd_areaDiag_PenF(i)
          
  !      endif
  !   enddo


  ! contains

  !   real*8 function indvd_areaDiag_PenF(i)
  !     integer :: j

  !     side = area_sides(i)
      
  !     jmax = (factorial(side)/(factorial(side-2) * factorial(2))) - side
      
  !     do j = 1,jmax
         
  !     enddo
      
      
  !   end function indvd_areaDiag_PenF

    
  ! end function diag_PenF


  real*8 function diag_PenF(dum_Nodes) 
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    integer :: i
    character(len=100) :: flnm
    
    diag_PenF = 0.0d0
    
    do i = 1,N_diag

       write(flnm,'(a10,i2.2,a10)')"PenF_Diag",i,"_vs_l0.dat" 
       open(unit=41,file=adjustl(trim(flnm)),position='append')

       diag_PenF = diag_PenF + indvd_diag_PenF(i)

       write(unit=41,fmt=*) l0(32),indvd_diag_PenF(i)
       !write(*,*) diag_PenF,"for spring i =",i
       close(41)
       
    enddo
    !write(*,*) diag_PenF,"diagPenF"
  contains
    
    real*8 function indvd_diag_PenF(diag_nmbr)
      implicit none
      integer :: diag_nmbr
      
      !if (diag_nmbr.eq.18) write(*,*) dflctn(diag_nmbr),"dflctn"
      indvd_diag_PenF = 0.5d0 * gammaa(diag_nmbr) * (dflctn(diag_nmbr))**2
      
    end function indvd_diag_PenF
    
    real*8 function dflctn(diag_nmbr)
      implicit none
      integer :: diag_nmbr
      integer :: Node1,Node2,Node3
      integer :: pulley_typ
      
      ! pulley_typ = 6
      
      ! if (typ_diag(diag_nmbr) .ne. pulley_typ) then
      !    Node1 = diag_node(diag_nmbr,1)
      !    Node2 = diag_node(diag_nmbr,2)
         
      !    l(diag_nmbr) = lngth_u_two_nodes(Node1,Node2)
      !    dflctn      = l(diag_nmbr) - Lt(diag_nmbr)

      !    !write(*,*) l(diag_nmbr),Lt(diag_nmbr),dflctn,"D1"
      !    !if (diag_nmbr.eq.18) write(*,*) dflctn,l(diag_nmbr),Lt(diag_nmbr),"dflctn,l,Lt for diagonal i =",diag_nmbr
      ! else
      !    Node1 = diag_node(diag_nmbr,1)
      !    Node2 = diag_node(diag_nmbr,2)
      !    Node3 = diag_node(diag_nmbr,3)
         
      !    l(diag_nmbr) = lngth_u_three_nodes(Node1,Node2,Node3)
      !    dflctn = l(diag_nmbr) - Lt(diag_nmbr)
      !    !write(*,*) l(diag_nmbr),Lt(diag_nmbr),dflctn,"D2"
         
      ! endif
      
    end function dflctn
    
    real*8 function lngth_u_two_nodes(Node1,Node2)
      implicit none
      integer :: Node1,Node2
      real*8  :: N1(1:N_dmnsn),N2(1:N_dmnsn)
      real*8  :: res_sq
      integer :: i

      N1 = 0.0d0 ; N2 = 0.0d0
      
      N1(1:N_dmnsn) = dum_Nodes(Node1,1:N_dmnsn)
      N2(1:N_dmnsn) = dum_Nodes(Node2,1:N_dmnsn)
          
      lngth_u_two_nodes = 0.0d0
      res_sq = 0.0d0

      open(unit=138,file='N2_N1_diag.dat',position='append')
      write(unit=138,fmt=*) N2,"N2"
      write(unit=138,fmt=*) N1,"N1"
      close(138)
      
      do i = 1,2
         res_sq = res_sq + (N2(i) - N1(i))**2 
      enddo
      
      lngth_u_two_nodes = sqrt(res_sq) 
      !write(*,*) lngth_u_two_nodes,"lutn"

    end function lngth_u_two_nodes

  end function diag_PenF
  
end module diag_PenFunc





module area_PenFunc
  use system_parameters
  use cell_info
  use transfrm_info
  use triangle_routines
  
contains

  real*8 function area_PenF(dum_Nodes)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    integer :: i
    character(len=100) :: flnm
    
    area_PenF = 0.0d0
    
    do i = 1,N_cell

       write(flnm,'(a10,i2.2,a10)')"PenF_Area",i,"_vs_l0.dat" 
       open(unit=41,file=adjustl(trim(flnm)),position='append')
       
       area_PenF = area_PenF + indvd_area_PenF(i)

       write(unit=41,fmt=*) indvd_area_PenF(i)
       !write(*,*) area_PenF,i,"area_PenF,i"
       close(41)
    enddo
    
    !write(*,*) area_PenF,"area_PenF"

  contains
    
    real*8 function indvd_area_PenF(area_nmbr)
      implicit none
      integer :: area_nmbr
      
      indvd_area_PenF = 0.5 * beta(area_nmbr) * (area_chng(area_nmbr))**2
      
    end function indvd_area_PenF
    
    real*8 function area_chng(area_nmbr)
      implicit none
      integer :: area_nmbr
      
      A(area_nmbr) = A_currnt(area_nmbr)
      area_chng    = A(area_nmbr) - At(area_nmbr)

      !if (area_nmbr .ge. 9)
      !write(*,*) A(area_nmbr),At(area_nmbr),"A_and_At"
      !if (area_nmbr .ge. 9) write(*,*) area_chng,"area_chng"
      
    end function area_chng
    
    real*8 function A_currnt(area_nmbr)
      implicit none
      integer :: area_nmbr
      integer :: Node1,Node2,Node3
      real*8  :: N1(1:3),N2(1:3),N3(1:3)
      integer :: i
      real*8  :: trngl

      N1=0.0d0 ; N2=0.0d0 ; N3=0.0d0
      
      trngl = 0.0d0
      
      A_currnt = 0.0d0
      
      i = 1
      
      open(unit=24,file="check_inside_areaE1.dat")
      
      do
         if ((i+2) .gt. max_node_area) then
            !write(*,*) "exit1"
            exit
         elseif (area_node(area_nmbr,(i+2)) .eq. (-1)) then
            !write(*,*) "exit2"
            exit
         endif
         
         N1 = 0.0d0 ; N2 = 0.0d0 ; N3 =0.0d0
         
         Node1 = area_node(area_nmbr,1)
         Node2 = area_node(area_nmbr,i+1)
         Node3 = area_node(area_nmbr,i+2)
         
         N1(1:2) = dum_Nodes(Node1,1:2)
         N2(1:2) = dum_Nodes(Node2,1:2)
         N3(1:2) = dum_Nodes(Node3,1:2)

         N1(3) = 0.0d0 ; N2(3) = 0.0d0 ; N3(3) = 0.0d0
         !write(unit=24,fmt=*) N1,N2,N3

         !write(*,*) N2(1:3),"N2_inside"
         !write(*,*) N3(1:3),"N3_inside"
         !write(*,*) N1(1:3),"N1_inside"
         
         !write(unit=24,fmt='(12(f8.6x),A)') dum_Nodes(1:6,1:2),"inside area_PenF"

         trngl = calc_triangle(N1,N2,N3)
         !write(*,*) trngl,"trngl"
         write(unit=24,fmt=*) trngl,"calc_triangle"
         
         A_currnt = A_currnt + trngl
         
         i = i+1
         
      enddo
      
      close(24)
      
    end function A_currnt
    
  end function area_PenF
  
end module area_PenFunc



module PenF_module
  use system_parameters
  use spr_PenFunc
  use area_PenFunc
  use conversion_routines

  implicit none
  real*8,allocatable :: grdmv_PenF(:) 
  
contains
  
  real*8 function PenF(dum_Nodes)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    integer :: lam_spr,lam_area
    integer :: check=1
    real*8  :: Ps,Pa
    
    !write(*,*) N_mvCoordnte,"N_mv"
    lam_spr = 0.0d0 ; lam_area = 0.0d0
    
    if (spr_lgcl.eqv..True.)    lam_spr    = 1.0d0
    if (area_lgcl.eqv..True.)   lam_area   = 1.0d0
   
    if (check.eq.1) then
       write(*,*) lam_spr,lam_area,"lams_in_PenF"
       check = check + 1
    endif
    
    !write(*,*) l0(31),l0(32),l0(34),l0(35),"l0 31-32-34-35"
    
    Ps = spr_PenF(dum_Nodes)
    Pa = area_PenF(dum_Nodes)
    
    
    PenF = lam_spr*Ps + lam_area*Pa 

    !write(*,*) Ps,Pa,"Ps,Pa"
 
    open(unit=101,file='All_PenFunc.dat',position='append')
    write(unit=101,fmt=*) Ps,Pa
    
    !write(*,*) Ps,"spr_PenF"
    !write(*,*) Pa,"area_PenF"

    close(101)
    
  end function PenF

end module PenF_module



module Wfunc_and_its_derivative
  use Energy_module
  use PenF_module
  use gradient_module
  use grdPenF_module
  
  implicit none
  real*8  :: alphaW=1.0d0,betaW=1.0d0
 
contains

  real*8 function Wfunc(dum_coordntes)
    implicit none
    real*8  :: dum_coordntes(1:N_mvCoordnte_withl0A0)
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8  :: duml0(1:N_spr),dumA0(1:N_cell)
    
    !real*8  :: Energy_Value,PenF_Value
    
    call coordntes_to_nodesl0A0(dum_coordntes,dum_Nodes,duml0,dumA0)
    
    !Energy_Value = Energy(dum_nodes,duml0,dumA0)
    !PenF_Value   = PenF(dum_Nodes)
    
    Wfunc = alphaW*Energy(dum_nodes,duml0,dumA0) + betaW*PenF(dum_Nodes)
    !Wfunc = alphaW*Energy_Value + betaW*PenF_Value
    
    !write(*,*) Energy_Value,PenF_Value,alphaW,betaW,"E,Penf,alpW,betW"
    
    
  end function Wfunc
  
  
  subroutine grdW_coordntes(dum_coordntesW,dum_grdWC)
    implicit none
    
    real*8,intent(in)  :: dum_coordntesW(1:N_mvCoordnte_withl0A0)
    real*8,intent(out) :: dum_grdWC(1:N_mvCoordnte_withl0A0) !C=coordnte_based

    real*8 :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8 :: duml0(1:N_spr),dumA0(1:N_cell)
    
    real*8 :: dum_grdN(1:N_node,1:N_dmnsn)!N=node_based
    
    real*8 :: coordnteE(1:N_mvCoordnte_withl0A0),coordntePenF(1:N_mvCoordnte)
    real*8 :: grdE_coordnte(1:N_mvCoordnte_withl0A0),grdPenF_coordnte(1:N_mvCoordnte)
    integer :: i
    
    coordnteE(1:N_mvCoordnte_withl0A0) = dum_coordntesW(1:N_mvCoordnte_withl0A0)
    coordntePenF(1:N_mvCoordnte) = dum_coordntesW(1:N_mvCoordnte)

    !do i = 1,N_mvCoordnte_withl0A0
     !  write(*,*) coordnteE(i),i,"coorE" 
    !enddo
    !do i = 1,N_mvCoordnte
     !  write(*,*) coordntePenF(i),i,"coorPenF"
    !enddo
    
    call gradient_coordntes(coordnteE,grdE_coordnte)
    call grdPenF_coordntes(coordntePenF,grdPenF_coordnte)
    call combine_all_variation
    
    
  contains
    
    subroutine combine_all_variation
      implicit none
      integer :: cnt
      integer :: i
     
      !write(*,*) alphaW,betaW
      
      do i = 1,N_mvCoordnte_withl0A0

         if ((l0_variatn_lgcl.eqv..True.).OR.(A0_variatn_lgcl.eqv..True.)) then
            !write(*,*) "Entering 1"
   
            if (i.le.N_mvCoordnte) then
               dum_grdWC(i) = alphaW*grdE_coordnte(i) + betaW*grdPenF_coordnte(i)
               !if (i==3) write(*,*) grdE_coordnte(i),grdPenF_coordnte(i),"i=3"
               
            elseif (i.gt.N_mvCoordnte) then
               dum_grdWC(i) = alphaW*grdE_coordnte(i)
            endif
            
         elseif ((l0_variatn_lgcl.eqv..False.).AND.(A0_variatn_lgcl.eqv..False.)) then
            !write(*,*) "Entering 2"
            dum_grdWC(i) = alphaW*grdE_coordnte(i)
            !write(*,*) dum_grdWC(i),grdE_coordnte(i),"Entering 2"
         endif
         
      enddo
      
    end subroutine combine_all_variation
   
      
  end subroutine grdW_coordntes


  subroutine grdW_Numrcl(dum_coordntesW,dum_grdWC)
    implicit none
    real*8,intent(inout) :: dum_coordntesW(1:N_mvCoordnte_withl0A0)
    real*8,intent(out)   :: dum_grdWC(1:N_mvCoordnte_withl0A0)

    real*8 :: dum_coordntesW_tmp(1:N_mvCoordnte_withl0A0)
    real*8 :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8 :: dum_grdN(1:N_node,1:N_dmnsn)
    real*8 :: duml0(1:N_spr),dumA0(1:N_cell)
    
    real*8  :: str,W_aft,W_bfr
    real*8  :: h=1e-4
    
    integer :: i
    real*8  :: EPS = 1.e-08
    
    dum_grdWC = 0.0d0
    
    dum_coordntesW_tmp = dum_coordntesW
    
    call coordntes_to_nodesl0A0(dum_coordntesW,dum_Nodes,duml0,dumA0)

    !write(*,*) l0(1),l0(2),l0(3),"l0 1-3"
    
    do i = 1,N_mvCoordnte_withl0A0
       
       str = dum_coordntesW(i)
       !write(*,*) str,"str"
       
       dum_coordntesW(i) = dum_coordntesW(i) + h
       write(*,*) dum_coordntesW(i),"dum_coordntesW"
       
       W_aft = Wfunc(dum_coordntesW)
       write(*,*) W_aft,"W_aft"
       
       dum_coordntesW(i) = str
       
       dum_coordntesW(i) = dum_coordntesW(i) - h
       write(*,*) dum_coordntesW(i),"dum_coordntesW"
       
       W_bfr = Wfunc(dum_coordntesW)
       write(*,*) W_bfr,"W_bfr"
       
       dum_coordntesW(i) = str
       
       dum_grdWC(i) = (W_aft-W_bfr)/(2.d0*h)
       
       if (abs(dum_coordntesW(i) - dum_coordntesW_tmp(i)) .gt. EPS) then
          write(*,*) "coordntes chnge in numrcl,flnm:PenF_andWfunc,sbrtn:grdW_Numrcl"
          stop
       endif
       
    enddo
    
    
    call print_numerical_gradients
    
  contains
    
    subroutine print_numerical_gradients
      implicit none
      integer :: i

      do i = 1,N_mvCoordnte_withl0A0
         !write(*,*) dum_grdWC(i),"grdWC_Numerical()"
      enddo
      
    end subroutine print_numerical_gradients
    
  end subroutine grdW_Numrcl


  

  subroutine check_Wfunc
    implicit none
    real*8  :: kspr_tmp(1:N_spr),karea_tmp(1:N_cell)
    real*8  :: l0_tmp(1:N_spr),A0_tmp(1:N_cell)
    real*8  :: Lt_tmp(1:N_spr),At_tmp(1:N_cell)
    real*8  :: l_tmp(1:N_spr), A_tmp(1:N_cell)
    
    real*8  :: W_expct,W_eval,diff
    integer :: N_itm
    integer :: i,j,jmax
    
    kspr_tmp  = k_spr
    karea_tmp = k_area

    l0_tmp = l0
    A0_tmp = A0
    Lt_tmp = Lt
    At_tmp = At
    l_tmp  = l
    A_tmp  = A

    N_itm = 2
    
    call make_all_k_zero

    open(unit=132,file='Wfunc_check.dat')

    do i = 1,N_itm
       if (i==1) jmax=N_spr
       if (i==1) jmax=N_cell
       
       do j = 1,jmax
       
          if (i==1) l0(j) = 1.0d0
          if (i==2) A0(j) = 10.0d0
          
          if (i==1) write(*,*) l(j),l0(j),Lt(j),"l(j),l0(j),Lt(j)"
          if (i==2) write(*,*) A(j),A0(j),At(j),"A(j),A0(j),At(j)"
          
          if (i==1) k_spr(j)  = kspr_tmp(j)
          if (i==2) k_area(j) = karea_tmp(j)
          
          if (i==1) W_expct = 0.5d0 * ( (l(j)-l0(j))**2 + (l(j)-Lt(j))**2 )
          if (i==2) W_expct = 0.5d0 * ( (A(j)-A0(j))**2 + (A(j)-At(j))**2 )
          
          W_eval  = Wfunc(node_xy)
          
          diff = abs(W_expct-W_eval)
          
          write(unit=132,fmt=*) W_expct,W_eval,diff,"W_expct,W_eval,diff"
       
          if(i==1) k_spr(j)  = 0.0d0
          if(i==2) k_area(j) = 0.0d0
          
       enddo
       
    enddo
    
    close(132)
    
  end subroutine check_Wfunc


  
  subroutine check_dWfunc
    implicit none
    
    real*8  :: kspr_tmp(1:N_spr),karea_tmp(1:N_cell)
    real*8  :: l0_tmp(1:N_spr),A0_tmp(1:N_cell)
    real*8  :: Lt_tmp(1:N_spr),At_tmp(1:N_cell)
    real*8  :: l_tmp(1:N_spr), A_tmp(1:N_cell)
    
    real*8  :: dW_analtl(1:N_mvCoordnte_withl0A0)
    real*8  :: dW_numrcl(1:N_mvCoordnte_withl0A0)
    real*8  :: diff(1:N_mvCoordnte_withl0A0)
    real*8  :: coordntesW(1:N_mvCoordnte_withl0A0)

    real*8  :: hm
    integer :: i
    
    kspr_tmp  = k_spr
    karea_tmp = k_area

    l0_tmp = l0
    A0_tmp = A0
    Lt_tmp = Lt
    At_tmp = At
    l_tmp  = l
    A_tmp  = A

    hm = 0.10d0
    
    write(*,*)N_mvCoordnte,N_spr,N_cell,N_mvCoordnte_withl0A0,"N_mvCoord,N_spr,N_cell,N_mv_withl0A0"

    do i = 1,N_node
       !write(*,*) node_xy(i,1:2),"node_xy HAHA"
    enddo
    
    !call make_all_k_zero
    
    !k_spr(1)  = kspr_tmp(1)

    do i = 1,N_spr
       l0(i) = (1.0d0-hm)*l0(i)
    enddo
    
    call nodesl0A0_to_coordntes(node_xy,l0,A0,coordntesW)
    call grdW_coordntes(coordntesW,dW_analtl)
    call grdW_Numrcl(coordntesW,dW_numrcl)

    !k_spr(1)  = 0.0d0
    
    open(unit=133,file='dWfunc_check.dat')
    
    do i = 1,N_mvCoordnte_withl0A0
       !write(*,*) coordntesW(i),"coordntesW"
       
       diff(i) = abs(dW_analtl(i) - dW_numrcl(i))
       write(unit=133,fmt=*) dW_analtl(i),dW_numrcl(i),diff(i),"dW_analtl,dW_numrcl,diff"

    enddo

    close(133)
    
  end subroutine check_dWfunc


  subroutine make_all_k_zero
    implicit none
    
    k_spr(1:N_spr)   = 0.0d0
    k_area(1:N_cell) = 0.0d0
    
  end subroutine make_all_k_zero
  
  
end module Wfunc_and_its_derivative
