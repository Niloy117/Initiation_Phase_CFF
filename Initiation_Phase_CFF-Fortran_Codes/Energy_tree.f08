module triangle_routines
  use system_parameters
  use transfrm_info
  
contains
  
  real*8 function calc_triangle(N1,N2,N3)
    implicit none
    real*8  :: N1(1:3),N2(1:3),N3(1:3)
    real*8  :: vectA(1:3),vectB(1:3),vectC(1:3)
    real*8  :: x1,x2,x3,y1,y2,y3
    integer :: ACW
    
    vectB = 0.0d0 ; vectC = 0.0d0 ; vectA = 0.0d0
    N1(3) = 0.0d0 ; N2(3) = 0.0d0 ; N3(3) = 0.0d0
    
    !write(*,*) N2(1:3),"N2" ;write(*,*) N3(1:3),"N3" ;write(*,*) N1(1:3),"N1"
    
    vectB = N2 - N1  !;write(*,*) vectB,"VECTB"
    vectC = N3 - N1  !;write(*,*) vectC,"VECTC"
    
    vectA         = cross_prdct(vectB,vectC)    !;write(*,*) vectA,"VECTA"
    calc_triangle = 0.5*vector_magnitude(vectA) !;write(*,*) calc_triangle,"calc_tr_inside"
    
    call check_triangle_area_value(calc_triangle)
    
    x1 = N1(1); x2 = N2(1); x3 = N3(1)
    y1 = N1(2); y2 = N2(2); y3 = N3(2)
    
    call check_ACW(x1,x2,x3,y1,y2,y3,ACW)

    if (ACW==1) then
       calc_triangle = calc_triangle
    elseif (ACW==(-1)) then
       calc_triangle = -calc_triangle
    elseif (ACW==0) then
       write(*,*) "ACW can't be zero, Energy_tree"
    endif   
    
  end function calc_triangle
  
  subroutine check_triangle_area_value(calc_triangle)
    implicit none
    real*8 :: calc_triangle
    
    if (calc_triangle .lt. 0.0) then
       !write(*,*) "WARNING:Negative triangle_area, check_again"
       !try to include location of this error, THINK
       !stop
    else
       !write(*,*) "You are Okay"
    endif
    
  end subroutine check_triangle_area_value
  
  function cross_prdct(vectB,vectC)
    
    implicit none
    real*8  :: vectB(1:3), vectC(1:3)
    real*8  :: vectA(1:3)
    real*8  :: cross_prdct(1:3)
    !real*8  :: strB(1:dmnsn), strC(1:dmnsn)
    !real*8  :: vectB_(1:3),vectC_(1:3)
    
    ! if (dmnsn .eq. 1) then
    !    cross_prdct = 0
    ! elseif (dmnsn .eq. 2) then
    !    call TwoD_to_ThreeD_vector(vectB,vectB_)
    !    call TwoD_to_ThreeD_vector(vectC,vectC_)
       
    !    call loop_for_cross_prdct(vectB_,vectC_,vectA)
    !    cross_prdct = vectA
    ! endif
    
    call loop_for_cross_prdct(vectB,vectC,vectA)
    cross_prdct = vectA
    
  end function cross_prdct
  
  subroutine loop_for_cross_prdct(A,B,C)
    implicit none
    real*8              :: LCT(1:3,1:3,1:3)
    real*8, intent(in)  :: A(1:3),B(1:3)
    real*8, intent(out) :: C(1:3)
  
    integer :: j1,j2,j3
    
    C(1:3) = 0.0d0
    
    call levi_civita_tensor(LCT)
    
    do j1 = 1,3
       do j2 = 1,3
          do j3 = 1,3
             C(j1) = C(j1) + LCT(j1,j2,j3) * A(j2) * B(j3) 
          enddo
       enddo
    enddo

  end subroutine loop_for_cross_prdct
  
  subroutine levi_civita_tensor(LCT)

    implicit none
    real*8 :: LCT(1:3,1:3,1:3)
    integer :: j1,j2,j3
    
    LCT(1:3,1:3,1:3) = 0.0d0
    
    do j1 = 1,3
       do j2 = 1,3
          if (j2.eq.j1) cycle
          do j3 = 1,3
             if (j3.eq.j2 .or. j3.eq.j1) cycle
             LCT(j1,j2,j3) = 1.0d0
             !write(*,*) LCT(j1,j2,j3),j1,j2,j3,"LCT_ji_j2_j3"
             
             if(j1.eq.1 .or. j1.eq.3) then
                if (j2 .gt. j3) LCT(j1,j2,j3) = -LCT(j1,j2,j3)
             elseif(j1.eq.2) then
                if (j2 .lt. j3) LCT(j1,j2,j3) = -LCT(j1,j2,j3)
             endif
             !write(*,*) LCT(j1,j2,j3),j1,j2,j3,"LCT_ji_j2_j3"
          enddo
       enddo
    enddo
    
  end subroutine levi_civita_tensor
  
  subroutine TwoD_to_ThreeD_vector(vectA,vectA_)

    implicit none
    real*8, intent(in)  :: vectA(1:2)
    real*8, intent(out) :: vectA_(1:3)
    real*8  :: strA(1:2)
    !integer :: i
    
    strA = vectA
    
    vectA_(1:2) = strA(1:2)
    vectA_(3)   = 0.0d0
    
  end subroutine TwoD_to_ThreeD_vector

  real*8 function vector_magnitude(vectA)
    implicit none
    
    real*8  :: vectA(1:3)
    real*8  :: magnitude_sqr
    integer :: i

  
    !write(*,*) vectA,"vectA"
    vector_magnitude = 0.0d0
    magnitude_sqr    = 0.0d0
    
    do i = 1,3
       magnitude_sqr = magnitude_sqr + vectA(i)*vectA(i)
    enddo

    vector_magnitude = sqrt(magnitude_sqr)
    
  end function vector_magnitude
  
  subroutine check_ACW(x1,x2,x3,y1,y2,y3,ACW)
    implicit none
    real*8, intent(in)  :: x1,x2,x3,y1,y2,y3
    integer,intent(out) :: ACW

    real*8 :: x12,y12,x23,y23
    real*8 :: MatVal
    
    open(unit=919,file='check_ACW.dat')
    
    x12 = x1-x2 ; y12 = y1-y2
    x23 = x2-x3 ; y23 = y2-y3
    
    write(919,*) x1,x2,x3,y1,y2,y3,"x1~"
    write(919,*) x12,y12,x23,y23,"x12~"
    close(919)
    
    MatVal = (x12*y23-y12*x23)
    
    if (MatVal.ge.0.0d0) then
       ACW = 1
    elseif (MatVal.lt.0.0d0) then
       ACW = -1
    endif
    
  end subroutine check_ACW
  
  recursive function factorial(val) result(fact)
    implicit none
    
    integer, intent(in) :: val
    integer :: fact
    
    if (val .eq. 0) then
       Fact = 1
    else
       Fact = val * Factorial(val-1)
    endif
    
  end function factorial
  
  
  subroutine chk_area_sign(node_nm,signV)
    implicit none
    integer, intent(in)  :: node_nm
    integer, intent(out) :: signV
    
    integer :: vrtx1,vrtx2
    real*8  :: x1,y1,x2,y2,xN,yN
    real*8  :: slope,IC
    real*8  :: SL,TOL
    
    interface
       real*8 function strght_line(xd,yd,slope,IC)
         implicit none
         real*8, intent(in) :: xd,yd
         real*8, intent(in) :: slope,IC
       end function strght_line
    end interface
    
    signV = 0
    TOL = 1.d-16
    
    vrtx1 = insrtd_vrtxNode(node_nm,1)
    vrtx2 = insrtd_vrtxNode(node_nm,2)
    
    x1 = node_xy(vrtx1,1) ; y1 = node_xy(vrtx1,2)
    x2 = node_xy(vrtx1,1) ; y2 = node_xy(vrtx1,2)
    
    xN = node_xy(node_nm,1) ; yN = node_xy(node_nm,2)
    
    if (abs(x1-x2).ge.TOL) then
       call get_slope_and_intercept(x1,y1,x2,y2,slope,IC)
       SL =  strght_line(xN,yN,slope,IC)

       if (SL.ge.0.0d0) then
          signV = 1
       elseif (SL.lt.0.0d0) then
          signV = -1
       endif
    else
       if (xN.ge.x1) then
          signV = 1
       elseif (xN.lt.x1) then
          signV = -1
       endif
    endif
    
  end subroutine chk_area_sign
  
end module triangle_routines

module spr_energy
  use system_parameters
  use cell_info
  use transfrm_info
  
  implicit none
  
contains
  
  real*8 function spr_E(dum_Nodes,duml0) !dummy, N_spring = N_spr
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8  :: duml0(1:N_spr)
    integer :: i
    
    spr_E = 0.0d0
    
    do i = 1,N_spr
       spr_E = spr_E + indvd_spr_E(i)
       !write(*,*) spr_E,"for spring i =",i
       
       !if (printIn_sprE==1) then
       !   write(*,*) i,k_spr(i),l(i),l0(i),spr_E,"Esval"
       !endif
       
    enddo
    !write(*,*) spr_E,"sprE"
    
  contains
    
    real*8 function indvd_spr_E(spr_nmbr)
      implicit none
      integer :: spr_nmbr
      
      !if (spr_nmbr.eq.18) write(*,*) dflctn(spr_nmbr),"dflctn"
      indvd_spr_E = 0.5d0 * k_spr(spr_nmbr) * (dflctn(spr_nmbr))**2

      !if (printIn_sprE == 1) then
      !   write(*,*) indvd_spr_E,k_spr(spr_nmbr),dflctn(spr_nmbr),spr_nmbr,"indvd aft print"
      !endif
      
    end function indvd_spr_E
    
    real*8 function dflctn(spr_nmbr)
      implicit none
      integer :: spr_nmbr
      integer :: Node1,Node2,Node3
      integer :: pulley_typ
      
      integer             :: nodeNum,cnt,nodeVal
      real*8, allocatable :: Nodes(:,:)
      
      !open(unit=79,file='len727374.dat',position='append')
      
      pulley_typ = 6
      
      !if (typ_spr(spr_nmbr) .ne. pulley_typ) then
      if (spr_node(spr_nmbr,0)==2) then
         
         Node1 = spr_node(spr_nmbr,1)
         Node2 = spr_node(spr_nmbr,2)
         
         l(spr_nmbr) = lngth_u_two_nodes(Node1,Node2)
         dflctn      = l(spr_nmbr) - duml0(spr_nmbr)
         
         !if ((spr_nmbr==72) .or. (spr_nmbr==73) .or. (spr_nmbr==74)) write(79,*) l(spr_nmbr),Node1,Node2,spr_nmbr,"72/73/74"
         !write(*,*) l(spr_nmbr),duml0(spr_nmbr),dflctn,"D1"
         !if (spr_nmbr.eq.18) write(*,*) dflctn,l(spr_nmbr),duml0(spr_nmbr),"dflctn,l,l0 for spring i =",spr_nmbr
         
      elseif (spr_node(spr_nmbr,0)==3) then
         
         Node1 = spr_node(spr_nmbr,1)
         Node2 = spr_node(spr_nmbr,2)
         Node3 = spr_node(spr_nmbr,3)
         
         l(spr_nmbr) = lngth_u_three_nodes(Node1,Node2,Node3)
         dflctn = l(spr_nmbr) - duml0(spr_nmbr)
         
         !if ((spr_nmbr==72) .or. (spr_nmbr==73) .or. (spr_nmbr==74)) write(79,*) l(spr_nmbr),Node1,Node2,spr_nmbr,"72/73/74"
         !write(*,*) l(spr_nmbr),duml0(spr_nmbr),dflctn,"D2"
         
      elseif (spr_node(spr_nmbr,0).gt.3) then
         
         nodeNum = spr_node(spr_nmbr,0)
         allocate(Nodes(1:nodeNum,1:N_dmnsn))
         
         do cnt = 1,nodeNum
            nodeVal              = spr_node(spr_nmbr,cnt)
            Nodes(cnt,1:N_dmnsn) = dum_Nodes(nodeVal,1:N_dmnsn)
         enddo
         
         l(spr_nmbr) = lngth_u_twoOrMore_nodes(nodeNum,Nodes)
         dflctn = l(spr_nmbr) - duml0(spr_nmbr)
         
         if (abs(dflctn) .le. 1.0d-15) then
            dflctn = 0.0d0
         endif
      endif
      
      !close(79)
      
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
      
      !open(unit=137,file='N2_N1.dat',position='append')
      !write(unit=137,fmt=*) N2,"N2"
      !write(unit=137,fmt=*) N1,"N1"
      !close(137)
      
      do i = 1,2
         res_sq = res_sq + (N2(i) - N1(i))**2 
      enddo

      !open(unit=212,file='res_sq.dat',position='append')
      !write(unit=212,fmt=*) res_sq
      !close(212)
      
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
    
    real*8 function lngth_u_twoOrMore_nodes(nodeNum,Nodes)
      implicit none
      integer, intent(in)  :: nodeNum
      real*8,  intent(in)  :: Nodes(1:nodeNum,1:N_dmnsn) 
      
      integer :: numOfprtns
      integer :: i,j
      real*8  :: res_sq
      real*8  :: N1(1:N_dmnsn),N2(1:N_dmnsn)
      
      numOfprtns              = nodeNum-1
      lngth_u_twoOrMore_nodes = 0.0d0
      
      
      do i = 1,numOfprtns
         
         res_sq        = 0.0d0
         N1(1:N_dmnsn) = Nodes(i,1:N_dmnsn)
         N2(1:N_dmnsn) = Nodes(i+1,1:N_dmnsn)
         
         do j = 1,N_dmnsn
            res_sq = res_sq + (N2(j)-N1(j))**2 
         enddo
         
         lngth_u_twoOrMore_nodes = lngth_u_twoOrMore_nodes + abs(sqrt(res_sq))
         
      enddo
      
    end function lngth_u_twoOrMore_nodes
    
    
  end function spr_E
  
end module spr_energy

module area_energy
  use system_parameters
  use cell_info
  use transfrm_info
  use triangle_routines
  use vitelline_fluid_calctn_mod
  
contains
  
  real*8 function area_E(dum_Nodes,dumA0)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8  :: dumA0(1:N_cell)
    integer :: i
    
    area_E = 0.0d0
    
    do i = 1,N_cell
       area_E = area_E + indvd_area_E(i)
       !write(*,*) area_E,i,"area_E,i"
    enddo
    
    
  contains
    
    real*8 function indvd_area_E(area_nm)
      implicit none
      integer :: area_nm
      
      indvd_area_E = 0.5*k_area(area_nm)*(area_chng(area_nm))**2
      
    end function indvd_area_E
    
    real*8 function area_chng(area_nm)
      implicit none
      integer :: area_nm
      
      if (VF_regionModelled==0) then
         A(area_nm) = A_currnt(area_nm)   
      elseif (VF_regionModelled==1) then
         if (area_nm.ne.N_cell) A(area_nm) = A_currnt(area_nm)
         if (area_nm == N_cell) A(area_nm) = vitelline_fluid_region(area_nm,dum_Nodes)
      endif
      
      area_chng  = A(area_nm) - dumA0(area_nm)
      
    end function area_chng
    
    real*8 function A_currnt(area_nm)
      implicit none
      integer :: area_nm
      integer :: Node1,Node2,Node3
      real*8  :: N1(1:3),N2(1:3),N3(1:3)
      integer :: i
      real*8  :: trngl
      integer :: signV
      
      N1=0.0d0 ; N2=0.0d0 ; N3=0.0d0
      
      trngl = 0.0d0
      signV = 0
      
      A_currnt = 0.0d0
      
      i = 1
      
      open(unit=25,file="check_inside_areaE.dat",position='append')
      
      do
         if ((i+2) .gt. max_node_area) then
            exit
         elseif (area_node(area_nm,(i+2)) .eq. (-1)) then
            exit
         endif
         
         N1 = 0.0d0 ; N2 = 0.0d0 ; N3 =0.0d0
         
         Node1 = area_node(area_nm,1)
         Node2 = area_node(area_nm,i+1)
         Node3 = area_node(area_nm,i+2)
         
         N1(1:2) = dum_Nodes(Node1,1:2)
         N2(1:2) = dum_Nodes(Node2,1:2)
         N3(1:2) = dum_Nodes(Node3,1:2)
         
         N1(3) = 0.0d0 ; N2(3) = 0.0d0 ; N3(3) = 0.0d0
         
         trngl = calc_triangle(N1,N2,N3)
         
         if (modelID==2.and.filechk==1) then
            
            if (area_nm==1) then
               write(25,*) "area_nm=",area_nm
               write(25,*) node1,N1(1:3),"N1"
               write(25,*) node2,N2(1:3),"N2"
               write(25,*) node3,N3(1:3),"N3"
               write(25,*) trngl,i,"calc_triangle,i"
               write(25,*) " "
               
            elseif (area_nm==(Hlf_Ncell+1)) then
               write(25,*) "area_nm=",area_nm
               write(25,*) node1,N1(1:3),"N1"
               write(25,*) node2,N2(1:3),"N2"
               write(25,*) node3,N3(1:3),"N3"
               write(25,*) trngl,i,"calc_triangle,i"
               write(25,*) " "
            endif
            
         endif
         
         A_currnt = A_currnt + trngl
         !write(25,*) A_currnt,trngl,"A_currnt,trngl"
         
         if (modelID==2.and.filechk==1) then
            
            if (area_nm==1) then
               write(25,*) A_currnt,trngl,"A_currnt,trngl"
               write(25,*) " "
               
            elseif (area_nm==(Hlf_Ncell+1)) then
               write(25,*) A_currnt,trngl,"A_currnt,trngl"
               write(25,*) " "
            endif
            
         endif
         
         i = i+1
         
      enddo
      
      if (modelID==2.and.filechk==1) then
         
         if (area_nm==1) then
            write(25,*) A_currnt,"A_currnt"
            write(25,*) " "
            
         elseif (area_nm==(Hlf_Ncell+1)) then
            write(25,*) A_currnt,"A_currnt"
            write(25,*) " "
         endif
         
      endif
      
      close(25)
      
    end function A_currnt
    
  end function area_E
  
  
end module area_energy


module gravitational_energy
  use system_parameters
  use cell_info

  implicit none
  real*8 :: x_grd = 0.0d0
  real*8 :: y_grd = 0.0d0 
  
contains
  
  real*8 function grvtnl_E(dum_Nodes)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    integer :: i
    real*8  :: delta_X=0.0d0, delta_Y=0.0d0
    
    !write(*,*) select_xy,"select_xy"
    
    grvtnl_E = 0.0d0
    
    
    do i = 1,N_node
       !write(*,*) select_xy,"select_xy"
       
       if (select_xy == 1) then !in x direction
          
          delta_X = dum_Nodes(i,1) - x_grd
          
          if (node_cnnctd(i)==0) then
             grvtnl_E = grvtnl_E + CgXNode(i)*delta_X
          elseif (node_cnnctd(i).ne.0) then
             if (count_this_dn(i)==1) then
                grvtnl_E = grvtnl_E + CgXNode(i)*delta_X
             endif
          endif
          
       elseif (select_xy == 2) then !in y direction
          
          delta_Y = dum_Nodes(i,2) - y_grd
          
          if (node_cnnctd(i)==0) then
             grvtnl_E = grvtnl_E + CgYNode(i)*delta_Y
          elseif (node_cnnctd(i).ne.0) then
             if (count_this_dn(i)==1) then
                grvtnl_E = grvtnl_E + CgYNode(i)*delta_Y
             endif
          endif
          
          
       elseif (select_xy == 3) then
          
          delta_X  = dum_Nodes(i,1) - x_grd
          delta_Y  = dum_Nodes(i,2) - y_grd 
          
          if (node_cnnctd(i)==0) then
             
             grvtnl_E = grvtnl_E + CgXNode(i)*delta_X + CgYNode(i)*delta_Y
             
          elseif (node_cnnctd(i).ne.0) then
             
             if (count_this_dn(i)==1) then
                grvtnl_E = grvtnl_E + CgXNode(i)*delta_X + CgYNode(i)*delta_Y
             endif
             
          endif
          
       endif
       
    enddo

    !write(*,*) grvtnl_E,"grvE"
    
  end function grvtnl_E

  
end module gravitational_energy

module SRYP_mod ! SRYP=Short-Range-Yukawa-Potential
  use system_parameters
  use cell_info
  use transfrm_info
  use triangle_routines
  use vitelline_fluid_calctn_mod
  
  implicit none
  
contains
  
  
  real*8 function SRyp_E(dum_Nodes) !SRyp_E=Short Rnge Yukwa Potential
    implicit none
    real*8, intent(in) :: dum_Nodes(1:N_node,1:2)
    real*8             :: yNodes(1:N_node)
    integer            :: i,j
    integer            :: apclNode
    real*8             :: yNodeV,numertr,denumertr
    real*8             :: ZeroV=0.00000000d0
    
    open(unit=872,file='SRyp_Chk.dat')
    
    yNodes(1:N_node) = dum_Nodes(1:N_node,2)
    SRyp_E           = 0.00d0 !; write(*,*) SRyp_E,"SRyp_E begin"
    
    if (SRyp_lgcl.eqv..False.) then
       SRyp_E = 0.0000d0
    elseif (SRyp_lgcl.eqv..True.) then
       
       do i = 1,numApclNodes
          apclNode = apclNodes(i)
          yNodeV   = yNodes(apclNode)
          
          if (activtnFctrApcl(i) == 0) then
             SRyp_E      = SRyp_E + ZeroV 
          elseif (activtnFctrApcl(i) == 1) then
             
             write(872,*)  AmpApclYP(i),AlphaApclYP(i),yNodeV,apclNode,"SRypchange 0"
             
             numertr   = AmpApclYP(i)*exp(AlphaApclYP(i)*(yNodeV))
             denumertr = (-yNodeV+epsApclYP(i))**(powrApclYP(i))
             write(872,*) apclNode,activtnFctrApcl(i),numertr,denumertr,i,"SRypchange1"
             SRyp_E    = SRyp_E + numertr/denumertr !;write(*,*)SRyp_E,nume,denume,i,"SRyp_E chng"
             write(872,*) SRyp_E,i,"SRypchange2"
          endif
          
       enddo
       
    endif
    
    close(872)
    
  end function SRyp_E
  
  subroutine get_yNodesApcl(yNodes,yNodesApcl)
    implicit none
    real*8, intent(in)    :: yNodes(1:N_node)
    real*8, intent(inout) :: yNodesApcl(1:numApclNodes)
    integer               :: i,j
    integer               :: apclNode
    real*8                :: yVal
    
    do i = 1,numApclNodes
       apclNode      = apclNodes(i)     ; write(*,*) apclNode,i,"apclNode,i"
       yNodesApcl(i) = yNodes(apclNode) ; write(*,*) yNodesApcl(i),i,"yNodesApcl,i"
    enddo
    
  end subroutine get_yNodesApcl
  
end module SRYP_mod

module repulsive_energy
  use system_parameters
  use cell_info
  use transfrm_info
  
  implicit none
  real*8 :: delta = 0.1d0


contains
  
  real*8 function rep_E(dum_nodes)
    implicit none
    integer :: p,q,r,s
    integer :: pmax,qmax,rmax,smax
    
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    integer :: node1,node2
    real*8  :: kr !will edit later
    
    real*8  :: vect(1:N_dmnsn),val
    !real*8  :: Heaviside_step
    
    pmax = N_cell        ; qmax = N_cell
    rmax = max_node_area ; smax = max_node_area
    
    rep_E = 0.0d0
    
    do p = 1,pmax

       do q = 1,qmax
          if (q.le.p) cycle
          
          do r = 1,rmax
             
             node1 = area_node(p,r)
             
             do s = 1,smax
                
                node2 = area_node(q,s)
                
                vect(1:2) = dum_nodes(node1,1:2) - dum_nodes(node2,1:2)
                val = vect_mgnitude(vect,N_dmnsn)
                
                rep_E = rep_E + 0.50d0*kr*((delta-val)**2)*Heaviside_step(delta,val)
                
             enddo
             
          enddo
       enddo
    enddo
    
    
  end function rep_E
  
  
  real*8 function vect_mgnitude(vect,dmnsn)
    implicit none
    integer :: dmnsn
    real*8  :: vect(1:dmnsn)
    
    integer :: i
    real*8  :: mgnitude_sq
    
    vect_mgnitude = 0.0d0
    mgnitude_sq = 0.0d0
    
    do i = 1,dmnsn
       mgnitude_sq = mgnitude_sq + vect(i)*vect(i)
    enddo
    
    vect_mgnitude = sqrt(mgnitude_sq)
    
  end function vect_mgnitude
  
  
  real*8 function Heaviside_step(delta,disp)
    implicit none
    real*8, intent(in) :: delta,disp
    real*8 :: EPS 
    
    EPS = 1e-08
    
    if (abs(delta-disp) .lt. EPS) then
       Heaviside_step = 1.0d0
    else
       Heaviside_step = 0.0d0
    endif
    
  end function Heaviside_step
  
  
end module repulsive_energy



module bending_energy
  use system_parameters
  use transfrm_info
  use neighbour_info_and_trnsfrmInfo_dependent_info
  implicit none
  
contains
  
  real*8 function bend_E(dum_Nodes)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    integer :: i
    
    !open(unit=224,file='bendE.dat')
    !open(unit=225,file='area_and_node_in_bendE.dat')
    !open(unit=226,file='cotSq_and_Phi.dat',position='append')
    !open(unit=227,file='bendE_for_onePhi.dat')
    
    bend_E = 0.0d0
      
    do i = 1,N_cell
       !write(226,fmt=*) i,"area_nm"
       bend_E = bend_E + indvd_area_bendE(i)
       !write(224,fmt=*) bend_E,i,"bend_E,i"
    enddo
    
    !close(224)
    !close(225)
    !close(226)
    !close(227)
    
    
  contains
    
    real*8 function indvd_area_bendE(area_nm)
      implicit none
      integer :: area_nm
      integer :: node_nm
      integer :: area_serial
      integer :: j,jmax
      real*8  :: cotSq,tanSq,Phi
      real*8  :: bendE_for_onePhi
      real*8  :: ZERO
      
      ZERO             = 0.0000000000000000d0
      indvd_area_bendE = 0.0d0
      bendE_for_onePhi = 0.0d0
      
      jmax = area_node(area_nm,0)
      !write(225,fmt=*) area_nm,"area_nm"
      
      do j = 1,jmax
         node_nm = area_node(area_nm,j)
         
         call replace_the_dnSecNode_with_fn(node_nm)
         !write(225,fmt=*) node_nm,j,"node_nm,pstn"
         !write(226,fmt=*) node_nm,j,"node_nm,pstn"
         
         call get_area_serial(node_nm,area_nm,area_serial)
         
         
         if (nodePhi_typ(node_nm)==1) then !1=cotSq
            
            if (abs(k_phi(node_nm,area_serial)-ZERO) .le. 1.0d-14) then
               bendE_for_onePhi = 0.0d0
            else
               call get_CotTanSq_and_Phi(node_nm,area_nm,area_serial,dum_Nodes,cotSq,tanSq,Phi)
               bendE_for_onePhi = 0.5d0*k_phi(node_nm,area_serial)*cotSq
            endif
            
         elseif (nodePhi_typ(node_nm)==2) then
            
            if (abs(k_phi(node_nm,area_serial)-ZERO) .le. 1.0d-14) then
               bendE_for_onePhi = 0.0d0
            else
               call get_CotTanSq_and_Phi(node_nm,area_nm,area_serial,dum_Nodes,cotSq,tanSq,Phi)
               bendE_for_onePhi = 0.5d0*k_phi(node_nm,area_serial)*tanSq
            endif
            
         endif
         
         indvd_area_bendE = indvd_area_bendE + bendE_for_onePhi
         
         !write(226,fmt=*) cotSq,tanSq,Phi,"cotSq,Phi"
         !write(227,fmt=*) bendE_for_onePhi,node_nm,area_serial
         
      enddo
      
      !write(225,fmt=*) " "
      !write(226,fmt=*) " "
      
    end function indvd_area_bendE
    
  end function bend_E
  
  
  subroutine get_CotTanSq_and_Phi(node_nm,area_nm,area_serial,dum_Nodes,cotSq,tanSq,Phi)
    implicit none
    integer, intent(in)  :: node_nm,area_nm,area_serial
    real*8,  intent(in)  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8,  intent(out) :: cotSq,tanSq,Phi
    
    real*8  :: CN(1:2),NN1(1:2),NN2(1:2)
    real*8  :: Aa,Bb
    integer :: i
    
    if (modelID==2) then
       !open(unit=276,file='Entered_node.dat')
       
       !do i = 1,N_node
       !   write(276,fmt=*) dum_Nodes(i,1:2),i
       !enddo
       !write(276,fmt=*) " "
       !close(276)
    endif
    
    call get_CN_NN1_NN2_nodes(node_nm,area_serial,dum_Nodes,CN,NN1,NN2)
    
    ! if (CellsMeet==1) then
    !    write(*,*) node_nm,area_serial,"node_nm,area_serial"
       
    !    do i = 1,N_node
    !       write(*,*) dum_Nodes(i,1:2),"dum_Nodes"
    !    enddo
       
    !    write(*,*) CN(1:2),"CN"
    !    write(*,*) NN1(1:2),"NN1"
    !    write(*,*) NN2(1:2),"NN2"
       
    ! endif
    
    call get_A_and_B(CN,NN1,NN2,Aa,Bb)
    
    if (nodePhi_typ(node_nm)==1) then
       
       if (abs(Bb-Aa).le.1.0d-15) then
          !write(*,*) Bb,Aa,abs(Bb-Aa),"Bb,Aa"
          !write(*,*) CN,NN1,NN2
          !stop
       endif
       
       cotSq = Aa / (Bb-Aa)
       tanSq = 10.0d5
    elseif (nodePhi_typ(node_nm)==2) then
       tanSq = (Bb-Aa) / Aa
       cotSq = 10.0d5
    endif
    
    !if(modelID==2) write(*,*) Aa,Bb,"Aa,Bb" 
    call get_Phi(CN,NN1,NN2,Phi)
    !if(modelID==2) write(*,*) Phi,"Phi"
    
  end subroutine get_CotTanSq_and_Phi
  
  subroutine get_CN_NN1_NN2_nodes(node_nm,area_serial,dum_Nodes,CN,NN1,NN2)
    implicit none
    integer, intent(in)  :: node_nm,area_serial
    real*8 , intent(in)  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8 , intent(out) :: CN(1:2),NN1(1:2),NN2(1:2)
    
    integer :: NeighNodes(1:2)
    
    NeighNodes(1:2) = Nlist(node_nm,area_serial,1:2)
    !if (modelID==2) write(*,*) NeighNodes(1:2),node_nm,"node,Neigh"
    
    if (NeighNodes(1).lt.0 .or. NeighNodes(2).lt.0) then
       write(*,*) node_nm,area_serial,NeighNodes(1:2),"node,area_serial,neighs"
       write(*,*)"Neighs can't be (-)ve, sb:get_CotTanSq_and_Phi,fl:Energy_tree"
       stop
    endif
    
    NN1(1:2) = dum_Nodes(NeighNodes(1),1:2)
    NN2(1:2) = dum_Nodes(NeighNodes(2),1:2)
    CN(1:2)  = dum_Nodes(node_nm,1:2)
    
    !if(SystemTyp == 1) then
    !   write(*,*) node_nm,NeighNodes(1),NeighNodes(2),"nodes and neighs"
    !   write(*,*) CN(1:2),"CN"
    !   write(*,*) NN1(1:2),"NN1"
    !   write(*,*) NN2(1:2),"NN2"
    !endif
    
    !if (modelID==2) write(*,*) NN1(1:2),NN2(1:2),CN(1:2),"NN1,NN2,CN"
    
  end subroutine get_CN_NN1_NN2_nodes
  
  
  subroutine get_Phi(CN,NN1,NN2,Phi)
    implicit none
    real*8, intent(in)  :: CN(1:2),NN1(1:2),NN2(1:2)
    real*8, intent(out) :: Phi
    
    real*8 :: Pi,Zero
    real*8 :: NN1x,NN1y,NN2x,NN2y,CNx,CNy
    
    real*8 :: vectA(1:N_dmnsn),vectB(1:N_dmnsn)
    real*8 :: vectA_val,vectB_val
    real*8 :: dotP,mltpAB,dotAB_by_AB
    
    interface
       real*8 function dot_prdct(vectA,vectB,dmnsn)
         implicit none
         integer,intent(in) :: dmnsn
         real*8, intent(in) :: vectA(1:dmnsn),vectB(1:dmnsn)
       end function dot_prdct

       real*8 function Mgnitude(vectA,dmnsn)
         implicit none
         integer,intent(in) :: dmnsn
         real*8 ,intent(in) :: vectA(1:dmnsn)
       end function Mgnitude
    end interface
    
    open(unit=125,file='Phi.dat',position='append')
    
    Pi = 4.0d0*atan(1.0d0)
    Zero = 0.0d0
    
    NN1x = NN1(1) ; NN1y = NN1(2)
    NN2x = NN2(1) ; NN2y = NN2(2)
    CNx  = CN(1)  ; CNy  = CN(2)
    
    vectA(1) = CNx-NN1x ; vectA(2) = CNy-NN1y
    vectB(1) = CNx-NN2x ; vectB(2) = CNy-NN2y
    
    dotP = dot_prdct(vectA,vectB,N_dmnsn)
    vectA_val = Mgnitude(vectA,N_dmnsn)
    vectB_val = Mgnitude(vectB,N_dmnsn)
    
    mltpAB = vectA_val*vectB_val
    dotAB_by_AB = dotP/(mltpAB)

    if ((dotAB_by_AB.gt.1.0d0).AND.(dotAB_by_AB.le.1.0001d0)) then
       dotAB_by_AB = 1.0d0
    elseif ((dotAB_by_AB.lt.(-1.0d0)).AND.(dotAB_by_AB.ge.(-1.0001d0))) then
       dotAB_by_AB = -1.0d0
    endif
    
    write(125,fmt=*) vectA_val,vectB_val,mltpAB,dotAB_by_AB
    
    if (dotAB_by_AB.gt.1.0d0 .or. dotAB_by_AB.lt.(-1.0d0)) then
       write(*,*) CNx,CNy,NN1x,NN1y,NN2x,NN2y,"CNxy,NN1xy,NN2xy"
       write(*,*) vectA(1:2),vectB(1:2),"vectA,vectB"
       write(*,*) vectA_val,vectB_val,"val of vectA,vectB"
       write(*,*) dotP,mltpAB,"dotP,mltAB"
       write(*,*) "out of range for cosine"
       stop
    endif
    Phi = acos(dotAB_by_AB)*(180.0d0/Pi)
    
    write(125,fmt=*) Phi,"Phi"
    close(125)
    
  end subroutine get_Phi
  
  subroutine get_A_and_B(CN,NN1,NN2,Aa,Bb)
    implicit none
    real*8, intent(in)  :: NN1(1:2),NN2(1:2),CN(1:2)
    real*8, intent(out) :: Aa,Bb
    
    real*8  :: NN1x,NN1y,NN2x,NN2y,CNx,CNy
    integer :: i
    
    NN1x = NN1(1) ; NN1y = NN1(2)
    NN2x = NN2(1) ; NN2y = NN2(2)
    CNx  = CN(1)  ; CNy  = CN(2)
    
    Aa = ((CNx-NN1x)*(CNx-NN2x) + (CNy-NN1y)*(CNy-NN2y))**2
    Bb = ((CNx-NN1x)**2+(CNy-NN1y)**2) * ((CNx-NN2x)**2+(CNy-NN2y)**2)
    
    if (Aa.gt.10.d15 .or. Bb.gt.10.d15) then
       write(*,*) CNx,CNy,NN1x,NN1y,NN2x,NN2y,"CurrNodeX-Y,NeighNode1X-Y,NeighNode2X-Y"
       write(*,*) (CNx-NN1x),(CNx-NN2x),(CNy-NN1y),(CNy-NN2y),"diff(1)"
       
       do i = 1,N_node
          write(*,*) node_xy(i,1:2),i,"dum in A_and_B"
       enddo
       
       stop
    endif
       
  end subroutine get_A_and_B
  
  subroutine get_Angle_AtANode(nodeNm,areaNm,Angle)
    implicit none
    integer, intent(in)  :: nodeNm,areaNm
    real*8 , intent(out) :: Angle
    
    integer :: areaSerial
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8  :: CN(1:N_dmnsn),NN1(1:N_dmnsn),NN2(1:N_dmnsn)
    
    dum_Nodes(1:N_node,1:N_dmnsn) = node_xy(1:N_node,1:N_dmnsn)
    
    call get_area_serial(nodeNm,areaNm,areaSerial)
    call get_CN_NN1_NN2_nodes(nodeNm,areaSerial,dum_Nodes,CN,NN1,NN2)
    call get_Phi(CN,NN1,NN2,Angle)
    
  end subroutine get_Angle_AtANode
  
end module bending_energy



module Energy_module
  use system_parameters
  use spr_energy
  use area_energy
  use gravitational_energy
  use bending_energy
  use SRYP_mod
  use conversion_routines
  
contains
  
  real*8 function Energy(dum_Nodes,duml0,dumA0)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8  :: duml0(1:N_spr),dumA0(1:N_cell)
    real*8  :: lam_spr,lam_area,lam_grvtnl,lam_bend,lam_SRyp
    integer :: check=1
    real*8  :: Es,Ea,Eg,Eb,Er
    
    lam_spr=0.0d0 ; lam_area=0.0d0 ; lam_grvtnl=0.0d0 ; lam_bend=0.0d0 ; lam_SRyp=0.0d0
    
    if (spr_lgcl.eqv..True.)    lam_spr    = 1.0d0
    if (area_lgcl.eqv..True.)   lam_area   = 1.0d0
    if (grvtnl_lgcl.eqv..True.) lam_grvtnl = 1.0d0
    if (bend_lgcl.eqv..True.)   lam_bend   = 1.0d0
    if (SRyp_lgcl.eqv..True.)   lam_SRyp   = 1.0d0
    
    !check = 1
    if (check.eq.1) then
       write(*,*) lam_spr,lam_area,lam_grvtnl,lam_bend,"lams_in_Energy"
       check = check + 1
    endif
    
    Es = spr_E(dum_Nodes,duml0)
    Ea = area_E(dum_Nodes,dumA0)
    Eg = grvtnl_E(dum_Nodes)
    Eb = bend_E(dum_Nodes)
    Er = SRyp_E(dum_Nodes)
    
    Energy = lam_spr*Es + lam_area*Ea + lam_grvtnl*Eg + lam_bend*Eb + lam_SRyp*Er
    
    open(unit=101,file='All_Energy.dat',position='append')
    write(101,*) Es,Ea,Eg,Eb,Er
    close(101)
    
  end function Energy
  
  
  real*8 function Energy_mv1(coordnte)
    implicit none
    
    real*8  :: coordnte(1:N_mvCoordnte)
    real*8  :: nodes(1:N_node,1:N_dmnsn)
    
    call coordntes_to_nodes(coordnte,nodes)
    
    Energy_mv1 = Energy(nodes,l0,A0)
    
  end function Energy_mv1
  
  real*8 function Energy_mv(dum_coordntes)
    implicit none
    
    real*8  :: dum_coordntes(1:N_mvCoordnte_withl0A0)
    real*8  :: dum_Nodes(1:N_dmnsn,1:N_dmnsn)
    real*8  :: duml0(1:N_spr)
    real*8  :: dumA0(1:N_cell)
    
    call coordntes_to_nodesl0A0(dum_coordntes,dum_Nodes,duml0,dumA0)
    
    Energy_mv = Energy(dum_Nodes,duml0,dumA0)
    
  end function Energy_mv
  
  
  real*8 function Energy_constrntBasd(dum_Nodes,duml0,dumA0)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8  :: duml0(1:N_spr),dumA0(1:N_cell)
    real*8  :: Es,Ea,Eg,Eb
    
    if (spr_lgcl_cnstr.eqv..True.)     Es = spr_E(dum_Nodes,duml0)
    if (spr_lgcl_cnstr.eqv..False.)    Es = 0.0d0
    
    if (area_lgcl_cnstr.eqv..True.)    Ea = area_E(dum_Nodes,dumA0)
    if (area_lgcl_cnstr.eqv..False.)   Ea = 0.0d0
    
    if (grvtnl_lgcl_cnstr.eqv..True.)  Eg = grvtnl_E(dum_Nodes)
    if (grvtnl_lgcl_cnstr.eqv..False.) Eg = 0.0d0
    
    if (bend_lgcl_cnstr.eqv..True.)    Eb = bend_E(dum_Nodes)
    if (bend_lgcl_cnstr.eqv..False.)   Eb = 0.0d0
    
    Energy_constrntBasd = Es+Ea+Eg+Eb
    
    open(unit=102,file='All_Energy_Cnstrnt.dat',position='append')
    write(unit=102,fmt=*) Es,Ea,Eg,Eb
    close(102)
    
  end function Energy_constrntBasd
  
  
end module Energy_module

