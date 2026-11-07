module curve_varsValue
  use system_parameters
  use transfrm_info
  
  implicit none
  
contains

  
  subroutine get_curveVars(dum_Nodes,Cell_Pressure,Tension)
    implicit none
    real*8, intent(in) :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8, intent(in) :: Cell_Pressure(1:N_cell)
    real*8, intent(in) :: Tension(1:N_spr)
    
    integer :: i,j,jmax
    integer :: spr_nm,typ_nm
    integer :: cnt,cnt_Pspr
    real*8  :: dP_crv,TC_crv,Rdius_crv
    real*8  :: Cx(1:2),Cy(1:2)
    real*8  :: slope,IC
    real*8  :: x1,y1,x2,y2,xc,yc
    real*8  :: Angl1,Angl2
    
    integer :: node1,node2
    integer :: area1,area2
    
    integer :: spr1,spr2
    real*8  :: ks_jnt,l0_jnt
    real*8  :: P1_val,P2_val

    cnt  = 0
    cnt_Pspr = 0
    jmax = 3
    
    do i = 1,N_cell  
       write(*,*) Cell_Pressure(i)
       
       do j = 1,jmax
          spr_nm = (i-1)*jmax + j
          typ_nm = typ_spr(spr_nm)
          
          if (i.le.(ncl+ncr)) then
             
             if (typ_nm==2 .or. typ_nm==4 .or. typ_nm==5) then
                cnt = cnt+1
                write(*,*) cnt,"cnt in 1"
                call get_oneCurveVars
                write(*,*) P1(cnt),P2(cnt),cnt,"P1P2,cnt"
                !stop
             endif
             
          elseif (i.gt.(ncl+ncr) .and. i.le.(ncl+ncr+2)) then
             !if (typ_nm .ne. 6) then
                !cnt = cnt+1
             !endif
             continue
             
          elseif (i.gt.(ncl+ncr+2)) then
             
             if (typ_nm == 3) then
                cnt_Pspr = cnt_Pspr+1
                !write(*,*) cnt_Pspr,"cnt_Pspr"
                if (cnt_Pspr==1) then
                   spr1 = spr_nm
                endif
                
                if (cnt_Pspr==2) then
                   cnt = cnt+1
                   write(*,*) cnt,"cnt in 2"
                   spr2 = spr_nm
                   write(*,*) spr1,spr2,"spr1-2"

                   call joint_spr_Prop(spr1,spr2,ks_jnt,l0_jnt)
                   write(*,*) ks_jnt,l0_jnt,"ks_jnt,l0_jnt"
                   
                   node1 = spr_node(spr1,1)
                   node2 = spr_node(spr1,2)
                   write(*,*) node1,node2,"node1-2"
                   
                   curveX1Y1(cnt,1:2) = dum_Nodes(node1,1:2)
                   curveX2Y2(cnt,1:2) = dum_Nodes(node2,1:2)
                   
                   area1 = spr_area(spr1,1)
                   area2 = spr_area(spr2,1)
                   write(*,*) area1,area2,"area1-2"
                   
                   P1(cnt) = Cell_Pressure(area1)
                   P2(cnt) = Cell_Pressure(area2)
                   
                   P1_val = P1(cnt) ; P2_val = P2(cnt)
                   dP(cnt)= abs(P1_val-P2_val)
                   
                   write(*,*) l(spr1),l(spr2),"l of spr1,spr2"
                   TC(cnt) = ks_jnt*(l0_jnt-l(spr1))
                   
                   dP_crv    = dP(cnt)
                   TC_crv    = TC(cnt)
                   write(*,*) dP_crv,TC_crv,"dP_crv,TC_crv"
                   
                   x1 = curveX1Y1(cnt,1) ; y1 = curveX1Y1(cnt,2)
                   x2 = curveX2Y2(cnt,1) ; y2 = curveX2Y2(cnt,2)
                   
                   !call get_Radius(dP_crv,TC_crv,Rdius_crv)
                   !R(cnt) = Rdius_crv
                   call get_Radius_ArcApproach(P1_val,P2_val,area1,area2,spr_nm,Rdius_crv)
                   R(cnt) = Rdius_crv
                   
                   write(*,*) x1,y1,x2,y2,"x1y1x2y2"
                   
                   call get_circle_center(x1,y1,x2,y2,Rdius_crv,Cx,Cy)
                   write(*,*) Cx(1:2),"Cx"
                   write(*,*) Cy(1:2),"Cy"
                   
                   call get_slope_and_intercept(x1,y1,x2,y2,slope,IC)
                   call get_allowable_center(area1,area2,P1_val,P2_val,Cx,Cy,xc,yc)
                   xc_crv(cnt) = xc ; yc_crv(cnt) = yc
                   
                   call get_CurveStrtEndAngl(x1,y1,x2,y2,xc,yc,Angl1,Angl2)
                   AnglS(cnt) = Angl1 ; AnglE(cnt) = Angl2
                   
                   cnt_Pspr = 0
                   
                endif
                
             elseif (typ_nm==4 .or. typ_nm==5) then
                cnt = cnt+1
                write(*,*) cnt,"cnt in 3"
                call get_oneCurveVars
             endif
             
          endif
          
       enddo
    enddo
    
  contains
  
    subroutine get_oneCurveVars
      implicit none
      
      node1 = spr_node(spr_nm,1)
      node2 = spr_node(spr_nm,2)
      write(*,*) node1,node2,spr_nm,"node1-2,spr"
      write(*,*) dum_Nodes(node1,1:2),dum_Nodes(node2,1:2)
      
      curveX1Y1(cnt,1:2) = dum_Nodes(node1,1:2)
      curveX2Y2(cnt,1:2) = dum_Nodes(node2,1:2)
      
      area1 = spr_area(spr_nm,1)
      area2 = spr_area(spr_nm,2)
      
      P1(cnt) = Cell_Pressure(area1)
      
      if (area2.lt.0) then
         P2(cnt) = 0.0d0
      else
         P2(cnt) = Cell_Pressure(area2)
      endif
      
      P1_val = P1(cnt) ; P2_val = P2(cnt)
      dP(cnt)= abs(P1_val-P2_val)
      
      TC(cnt) = Tension(spr_nm)
      
      dP_crv    = dP(cnt)
      TC_crv    = TC(cnt)
      
      x1 = curveX1Y1(cnt,1) ; y1 = curveX1Y1(cnt,2)
      x2 = curveX2Y2(cnt,1) ; y2 = curveX2Y2(cnt,2)
      
      !call get_Radius(dP_crv,TC_crv,Rdius_crv)
      !R(cnt) = Rdius_crv

      call get_Radius_ArcApproach(P1_val,P2_val,area1,area2,spr_nm,Rdius_crv)
      R(cnt) = Rdius_crv
      
      call get_circle_center(x1,y1,x2,y2,Rdius_crv,Cx,Cy)
      call get_slope_and_intercept(x1,y1,x2,y2,slope,IC)
      
      call get_allowable_center(area1,area2,P1_val,P2_val,Cx,Cy,xc,yc)
      xc_crv(cnt) = xc ; yc_crv(cnt) = yc
      
      call get_CurveStrtEndAngl(x1,y1,x2,y2,xc,yc,Angl1,Angl2)
      AnglS(cnt) = Angl1 ; AnglE(cnt) = Angl2

    end subroutine get_oneCurveVars
    
  end subroutine get_curveVars

  
  subroutine get_Radius(dP_crv,TC_crv,Rdius_crv)
    implicit none
    real*8, intent(in)  :: dP_crv,TC_crv
    real*8, intent(out) :: Rdius_crv

    real*8 :: TOL 

    TOL = 1.d-16
    
    if(dP_crv.lt.TOL) then
       Rdius_crv = 10000.d0
    else
       Rdius_crv = abs(TC_crv)/dP_crv
    endif
    
  end subroutine get_Radius

  subroutine get_Radius_ArcApproach(P1,P2,area1,area2,spr_nm,Rdius_crv)
    implicit none
    real*8,  intent(in)  :: P1,P2
    integer, intent(in)  :: area1,area2
    integer, intent(in)  :: spr_nm
    real*8,  intent(out) :: Rdius_crv

    real*8 :: Neutrl(1:N_dmnsn)
    real*8 :: sprMid(1:N_dmnsn)
    real*8 :: NeutrlAxisDist,Diamtr_crv,dll0 !dll0=diff_in_l_and_l0
    real*8 :: TOL 
    
    interface
       real*8 function distance(Point1,Point2,N_dim)
         implicit none
         integer, intent(in) :: N_dim
         real*8 , intent(in) :: Point1(1:N_dim),Point2(1:N_dim)
       end function distance
    end interface
    
    TOL = 1.d-16

    open(unit=342,file='EqNgtvPressure.dat',position='append')
    
    if (P1.gt.P2) then
       call get_Neutral_Axis(area1,Neutrl)
       
    elseif (P2.gt.P1) then
       
       if (area2.lt.0) then
          write(342,fmt=*) "negativePressure",P1,P2
          write(342,fmt=*) 
       endif
       
       call get_Neutral_Axis(area2,Neutrl)
       
    else
       write(unit=342,fmt=*) "Pressure is equal in both cells",area1,area2
    endif
    
    close(342)
    
    call get_sprMid(spr_nm,sprMid)
    
    NeutrlAxisDist = distance(Neutrl,sprMid,N_dmnsn)

    dll0 = l(spr_nm)-l0(spr_nm)
    
    if(abs(dll0).lt.TOL) then
       Rdius_crv = 10000.d0
    else
       Rdius_crv = abs((l0(spr_nm)/dll0)*NeutrlAxisDist)
    endif
    
    Diamtr_crv = 2.d0*Rdius_crv

    if (Diamtr_crv.lt.l(spr_nm)) then
       write(*,*) "Not structurally possible"
       write(*,*) l(spr_nm),Diamtr_crv,"l,Diamtr"
       stop
    endif
    
    !write(*,*) Rdius_crv,"R"
    
  end subroutine get_Radius_ArcApproach
  
end module curve_varsValue




module neighbour_info_and_trnsfrmInfo_dependent_info
  use system_parameters
  use transfrm_info
  use curve_varsValue
  implicit none
  
contains
  
  subroutine get_Nlist
    implicit none
    integer :: i,j
    integer :: n_areaCnn
    integer :: area_nm,node_nm
    integer :: neighbours(1:2)

    open(unit=121,file='Nlist.dat')
    
    do i = 1,N_node
       n_areaCnn = node_area(i,0)
       write(121,fmt=*) i,"node_nm"
       !write(*,*)       i,"node_nm"
       
       do j = 1,n_areaCnn
          area_nm = node_area(i,j)
          node_nm = i
          write(121,fmt=*) area_nm,j,"area_nm,area serial"
          !write(*,*)       area_nm,j,"area_nm,area serial"
          
          call get_neighbour_nodes(area_nm,node_nm,neighbours)
          call replace_the_dn_with_fn(neighbours)
          
          Nlist(i,j,1:2) = neighbours(1:2)
          write(121,fmt=*) Nlist(i,j,1:2),"neighbours"
          !write(*,*)       Nlist(i,j,1:2),"neighbours"
          
       enddo
       
       write(121,fmt=*) " "
       !write(*,*)       " "
       
    enddo
    
    close(121)
    
  end subroutine get_Nlist
  
  
  subroutine replace_the_dn_with_fn(neighbours)
    implicit none
    integer, intent(inout) :: neighbours(1:2)
    integer :: i
    integer :: node_nm
    logical :: lgcl_dn
    integer :: double_node_num,node_positn
    integer :: other
    
    do i = 1,2
       node_nm = neighbours(i)
       call is_it_a_dn(node_nm,lgcl_dn)

       if (lgcl_dn .eqv. .True.) then
          call first_or_second_node_of_double_node(node_nm,double_node_num,node_positn)
          
          if (node_positn==2) then
             call get_the_other_dn(node_nm,other)
             neighbours(i) = other
          endif
       
       elseif (lgcl_dn .eqv. .False.) then
          continue
       endif
       
    enddo
    
  end subroutine replace_the_dn_with_fn

  subroutine replace_the_dnSecNode_with_fn(node_nm)
    implicit none
    integer, intent(inout) :: node_nm
    integer :: i
    logical :: lgcl_dn
    integer :: double_node_num,node_positn
    integer :: other
    
    
    call is_it_a_dn(node_nm,lgcl_dn)

    if (lgcl_dn .eqv. .True.) then
       call first_or_second_node_of_double_node(node_nm,double_node_num,node_positn)
       
       if (node_positn==2) then
          call get_the_other_dn(node_nm,other)
          node_nm = other
       endif
       
    elseif (lgcl_dn .eqv. .False.) then
       continue
    endif
    
  end subroutine replace_the_dnSecNode_with_fn
  
  subroutine get_area_serial(node_nm,area_nm,area_serial)
    implicit none
    integer, intent(in)  :: node_nm,area_nm
    integer, intent(out) :: area_serial
    
    integer :: i,imax
    integer :: curr_area
    
    open(unit=341,file='area_serial.dat',position='append')

    !write(*,*) node_nm,"node"
    imax = node_area(node_nm,0)
    !write(*,*) node_area(node_nm,0:imax),node_nm,"node_nm"
    
    do i = 1,imax
       curr_area = node_area(node_nm,i)
       
       if (curr_area==area_nm) then
          area_serial = i
          exit
       endif
       
    enddo
    
    write(341,fmt=*) node_nm,area_nm,area_serial,"node,area,area_serial"
    
    close(341)
    
  end subroutine get_area_serial
  
  
  subroutine get_phiInfo
    implicit none
    integer :: i,j,jmax
    integer :: cnt_phi
    integer :: area_nm,node_nm
    
    cnt_phi = 1
    
    do i =1,N_cell
       area_nm = i
       jmax    = area_node(i,0)
       
       do j = 1,jmax
          node_nm = area_node(i,j)
          
          phi_info(cnt_phi,1) = node_nm
          phi_info(cnt_phi,2) = area_nm
          
          cnt_phi = cnt_phi + 1
       enddo
       
    enddo
    
  end subroutine get_phiInfo
  
  
  
end module neighbour_info_and_trnsfrmInfo_dependent_info
