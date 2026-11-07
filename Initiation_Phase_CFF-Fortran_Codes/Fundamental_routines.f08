real*8 function distance(Point1,Point2,N_dim)
  implicit none
  integer, intent(in) :: N_dim
  real*8 , intent(in) :: Point1(1:N_dim),Point2(1:N_dim)
  
  real*8  :: res_sq
  integer :: i
  
  res_sq = 0.0d0
  
  do i = 1,N_dim
     res_sq = res_sq + (Point1(i)-Point2(i))**2
  enddo
  
  distance = sqrt(res_sq)
  distance = abs(distance)
  
end function distance

function TnsnComprsn(dum_Nodes,duml0)
  use system_parameters
  use transfrm_info
  
  implicit none
  real*8  :: TnsnComprsn(1:N_spr)
  real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
  real*8  :: duml0(1:N_spr)
  integer :: i

  if (strctNo==1) then
     open(unit=321,file='TnsnCmpr_in_strct1.dat',position='append')
     open(unit=322,file='ll0_inTnsnCmpr_strct1.dat',position='append')
  
  elseif (strctNo==2) then
     open(unit=321,file='TnsnCmpr_in_strct2.dat',position='append')
     open(unit=322,file='ll0_inTnsnCmpr_strct2.dat',position='append')
  endif
  
  do i = 1,N_spr
     TnsnComprsn(i) = -k_spr(i) * dflctn(i)
     !write(*,*) TnsnComprsn(i),i,"TnsnComprsn"
     !if (strctNo==1.and.i.ge.1) write(321,fmt=*) TnsnComprsn(i),i,"Tnsn,spr_nm"
     if (strctNo.le.2) write(321,fmt=*) TnsnComprsn(i),i,"Tnsn,spr_nm"
     !if (modelID==2) then
        !write(*,*) k_spr(i),dflctn(i),i,"k,dflctn,spr"
     !endif
  enddo

  if (strctNo==1) then
     close(321)
     close(322)
  elseif (strctNo==2) then
     close(321)
     close(322)
  endif
  
contains
  
  real*8 function dflctn(spr_nmbr)
    implicit none
    
    integer :: spr_nmbr
    integer :: Node1,Node2,Node3
    real*8  :: N1(1:N_dmnsn),N2(1:N_dmnsn),N3(1:N_dmnsn)
    integer :: nodeNum,nodeVal,cnt
    
    real*8, allocatable :: Nodes(:,:)
    
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
       
       real*8 function Length_MultiPrtn(nodeNum,Nodes,N_dmnsn)
         implicit none
         integer, intent(in) :: nodeNum
         real*8 , intent(in) :: Nodes(1:nodeNum,1:N_dmnsn)
         integer, intent(in) :: N_dmnsn
       end function Length_MultiPrtn
    end interface
    
    
    if (spr_node(spr_nmbr,0) == 2) then
       !write(*,*) spr_nmbr,"spr_nmbr"
       
       Node1 = spr_node(spr_nmbr,1)
       Node2 = spr_node(spr_nmbr,2)
       
       N1(1:N_dmnsn) = dum_Nodes(Node1,1:N_dmnsn)
       N2(1:N_dmnsn) = dum_Nodes(Node2,1:N_dmnsn)
       
       l(spr_nmbr) = Length(N1,N2,N_dmnsn)
       dflctn      = l(spr_nmbr) - duml0(spr_nmbr)
       
       !if(strctNo==1.and.spr_nmbr.ge.1)write(322,fmt=*) l(spr_nmbr),duml0(spr_nmbr),Lt(spr_nmbr),spr_nmbr,"l,l0,Lt,spr_nm"
       if(strctNo==1.and.spr_nmbr.ge.1) write(322,fmt=*) l(spr_nmbr),duml0(spr_nmbr),spr_nmbr,"l,l0,spr_nm"
       
       if (spr_nmbr==2) then
          !write(*,*) l(spr_nmbr),duml0(spr_nmbr),k_spr(spr_nmbr),spr_nmbr,"ll0"
       endif
       
       !if (spr_nmbr==2) then
        !  write(323,fmt=*) l(spr_nmbr),duml0(spr_nmbr),k_spr(spr_nmbr),spr_nmbr
       !endif
       
    elseif (spr_node(spr_nmbr,0)==3) then
       !write(*,*) spr_nmbr,"spr_nmbr"
       
       Node1 = spr_node(spr_nmbr,1)
       Node2 = spr_node(spr_nmbr,2)
       Node3 = spr_node(spr_nmbr,3)
       
       N1(1:N_dmnsn) = dum_Nodes(Node1,1:N_dmnsn)
       N2(1:N_dmnsn) = dum_Nodes(Node2,1:N_dmnsn)
       N3(1:N_dmnsn) = dum_Nodes(Node3,1:N_dmnsn)
       
       l(spr_nmbr) = Length_2Prtn(N1,N2,N3,N_dmnsn)
       dflctn = l(spr_nmbr) - duml0(spr_nmbr)
       
       if(strctNo==1.and.spr_nmbr.ge.1) write(322,fmt=*) l(spr_nmbr),duml0(spr_nmbr),spr_nmbr,"l,l0,spr_nm"
       
       if (spr_nmbr==11.or.spr_nmbr==23.or.spr_nmbr==27.or.spr_nmbr==30) then
          write(323,fmt=*) l(spr_nmbr),duml0(spr_nmbr),k_spr(spr_nmbr),spr_nmbr
       endif
       
       if (spr_nmbr==21.or.spr_nmbr==24.or.spr_nmbr==26.or.spr_nmbr==29) then
          write(323,fmt=*) l(spr_nmbr),duml0(spr_nmbr),k_spr(spr_nmbr),spr_nmbr
       endif
       
       
    elseif (spr_node(spr_nmbr,0).gt.3) then
       
       nodeNum = spr_node(spr_nmbr,0)
       allocate(Nodes(1:nodeNum,1:N_dmnsn))
       
       do cnt = 1,nodeNum
          nodeVal              = spr_node(spr_nmbr,cnt)
          Nodes(cnt,1:N_dmnsn) = dum_Nodes(nodeVal,1:N_dmnsn)
       enddo
       
       l(spr_nmbr) = Length_MultiPrtn(nodeNum,Nodes,N_dmnsn)
       dflctn = l(spr_nmbr) - duml0(spr_nmbr)
       
    endif
    
    if (abs(dflctn) .le. 1e-16) then
       dflctn = 0.0d0
    endif
       
  end function dflctn
         
end function TnsnComprsn


! function Bending(dum_Nodes,duml0)
!   use system_parameters
!   use transfrm_info
  
!   implicit none
!   real*8  :: Bending(1:N_spr)
!   real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
!   real*8  :: duml0(1:N_spr)
!   integer :: i
  
!   if (strctNo==1) then
!      open(unit=321,file='Bending_in_strct1.dat',position='append')
!      open(unit=322,file='ll0_inBending_strct1.dat',position='append')
!   endif
  
!   do i = 1,N_spr
!      Bending(i) = -k_spr(i) * dflctn(i)
!      !write(*,*) Bending(i),i,"Bending"
!      if (strctNo==1.and.i.ge.1) write(321,fmt=*) Bending(i),i,"Bending,spr_nm"
!   enddo

!   if (strctNo==1) then
!      close(321)
!      close(322)
!   endif
         
! end function Bending
  
real*8 function Length(N1,N2,N_dmnsn)
  
  implicit none
  integer, intent(in) :: N_dmnsn
  real*8 , intent(in) :: N1(1:N_dmnsn),N2(1:N_dmnsn)
  
  real*8  :: res_sq
  integer :: i
  
  Length = 0.0d0
  res_sq = 0.0d0
  
  do i = 1,N_dmnsn
     res_sq = res_sq + (N1(i)-N2(i))**2 
  enddo
  
  Length = sqrt(res_sq) 
  !write(*,*) lngth_u_two_nodes,"lutn"
  
end function Length

real*8 function Length_2Prtn(N1,N2,N3,N_dmnsn)
  implicit none
  real*8 , intent(in) :: N1(1:N_dmnsn),N2(1:N_dmnsn),N3(1:N_dmnsn)
  integer, intent(in) :: N_dmnsn
  real*8  :: LP1,LP2
  
  interface
     real*8 function Length(N1,N2,N_dmnsn)
       implicit none
       integer, intent(in) :: N_dmnsn
       real*8 , intent(in) :: N1(1:N_dmnsn),N2(1:N_dmnsn)
     end function Length
  end interface
  
  Length_2Prtn = 0.0d0
     
  LP1 = Length(N1,N2,N_dmnsn)
  LP2 = Length(N2,N3,N_dmnsn)
  
  Length_2Prtn = LP1 + LP2
   
end function Length_2Prtn


real*8 function Length_MultiPrtn(nodeNum,Nodes,N_dmnsn)
  implicit none
  integer, intent(in)  :: nodeNum
  real*8,  intent(in)  :: Nodes(1:nodeNum,1:N_dmnsn) 
  integer, intent(in)  :: N_dmnsn
  
  integer :: numOfprtns
  integer :: i,j
  real*8  :: res_sq
  real*8  :: N1(1:N_dmnsn),N2(1:N_dmnsn)
  
  numOfprtns       = nodeNum-1
  Length_MultiPrtn = 0.0d0
  
  
  do i = 1,numOfprtns
     
     res_sq        = 0.0d0
     N1(1:N_dmnsn) = Nodes(i,1:N_dmnsn)
     N2(1:N_dmnsn) = Nodes(i+1,1:N_dmnsn)
     
     do j = 1,N_dmnsn
        res_sq = res_sq + (N2(j)-N1(j))**2 
     enddo
     
     Length_MultiPrtn = Length_MultiPrtn + abs(sqrt(res_sq))
         
  enddo
  
end function Length_MultiPrtn


function Pressure(dum_Nodes,dumA0)
  use system_parameters
  use transfrm_info
  
  implicit none
  real*8  :: Pressure(1:N_cell)
  real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
  real*8  :: dumA0(1:N_cell)
  integer :: i
  
  !strctNo=1
  !write(*,*) strctNo,"strctNo"
  
  if (strctNo==1) then
     open(unit=321,file='Pressure_in_strct1.dat',position='append')
     open(unit=322,file='AA0_inPress_strct1.dat',position='append')
     open(unit=323,file='karea_strct1.dat',position='append')
  elseif (strctNo==2) then
     open(unit=421,file='Pressure_in_strct2.dat',position='append')
     open(unit=422,file='AA0_inPress_strct2.dat',position='append')
     open(unit=423,file='karea_strct2.dat',position='append')
  endif
  
  do i = 1,N_cell
     Pressure(i) = -k_area(i) * area_chng(i)
     if (strctNo==1.and.i.ge.1) write(321,fmt=*) Pressure(i),i
     if (strctNo==1.and.i.ge.1) write(323,fmt=*) k_area(i),i
     if (strctNo==2.and.i.ge.1) write(421,fmt=*) Pressure(i),i
     if (strctNo==2.and.i.ge.1) write(423,fmt=*) k_area(i),i
  enddo  

  if (strctNo==1) then
     close(321)
     close(322)
     close(323)
  elseif (strctNo==2) then
     close(421)
     close(422)
     close(423)
  endif
  
contains
  
  real*8 function area_chng(area_nm)
    implicit none
    integer :: area_nm
    integer :: Nnode_poly
    
    real*8,allocatable :: Nodes(:,:)
    
    integer :: j
    integer :: node_nm
    real*8  :: Poly_A
    
    interface
       real*8 function Polygon_Area(Nodes,Nnode_poly)
         use Triangle_routines         
         implicit none
         integer, intent(in) :: Nnode_poly
         real*8 , intent(in) :: Nodes(1:Nnode_poly,1:3)
       end function Polygon_Area
       
       real*8 function vitelline_fluid_region_in_FR(area_nm,dum_Nodes)
         use vitelline_fluid_calctn_mod
         implicit none
         integer, intent(in) :: area_nm
         real*8 , intent(in) :: dum_Nodes(1:N_node,1:N_dmnsn)
       end function vitelline_fluid_region_in_FR
    end interface
    
    Nnode_poly = area_node(area_nm,0)
    allocate(Nodes(1:Nnode_Poly,1:3))
    
    do j = 1,Nnode_poly
       node_nm = area_node(area_nm,j)
       Nodes(j,1:N_dmnsn) = dum_Nodes(node_nm,1:N_dmnsn)
    enddo
    
    if (VF_regionModelled==0) then
       A(area_nm) = Polygon_Area(Nodes,Nnode_poly)
    elseif (VF_regionModelled==1) then
       if (i.ne.N_cell) A(area_nm) = Polygon_Area(Nodes,Nnode_poly)
       if (i == N_cell) A(area_nm) = vitelline_fluid_region_in_FR(area_nm,dum_Nodes)
    endif
       
    area_chng  = A(area_nm) - dumA0(area_nm)
    
    if(strctNo==1.and.area_nm.ge.1) then
       write(322,fmt=*) A(area_nm),dumA0(area_nm),area_nm,"A,A0,area_nm" 
    elseif(strctNo==2.and.area_nm.ge.1) then
       write(422,fmt=*) A(area_nm),dumA0(area_nm),area_nm,"A,A0,area_nm" 
    endif
    
    deallocate(Nodes)
    
  end function area_chng
  
end function Pressure

  
real*8 function Polygon_Area(Nodes,Nnode_poly)
  use Triangle_routines
  
  implicit none
  integer, intent(in) :: Nnode_poly
  real*8 , intent(in) :: Nodes(1:Nnode_poly,1:3)
  
  real*8  :: N1(1:3),N2(1:3),N3(1:3)
  integer :: i
  real*8  :: areachk
  
  Polygon_Area = 0.0d0
  
  open(unit=24,file="check_inside_areaE.dat")
  
  do i = 1,(Nnode_poly-2)
     N1=0.0d0 ; N2=0.0d0 ; N3=0.0d0
     
     N1(1:3) = Nodes(1,1:3)
     N2(1:3) = Nodes(i+1,1:3)
     N3(1:3) = Nodes(i+2,1:3)

     if (N_dmnsn==2) then
        N1(3)=0.0d0 ; N2(3)=0.0d0 ; N3(3)=0.0d0
     endif
     
     !write(unit=24,fmt='(12(f8.6x),A)') dum_Nodes(1:6,1:2),"inside area_E"
     
     areachk = calc_triangle(N1,N2,N3)
     !write(*,*) areachk,"areachk"
     
     Polygon_Area = Polygon_Area + calc_triangle(N1,N2,N3)   
  enddo
  
  close(24)
    
end function Polygon_Area


real*8 function vitelline_fluid_region_in_FR(area_nm,dum_Nodes)
  use vitelline_fluid_calctn_mod
  implicit none
  integer, intent(in) :: area_nm
  real*8 , intent(in) :: dum_Nodes(1:N_node,1:N_dmnsn)
  
  vitelline_fluid_region_in_FR = vitelline_fluid_region(area_nm,dum_Nodes)
  
end function vitelline_fluid_region_in_FR

real*8 function LengthForArbtrySgmnt(Nodes,Nnode_Spr)
  implicit none
  integer, intent(in) :: Nnode_Spr
  real*8 , intent(in) :: Nodes(1:Nnode_Spr,1:3)
  
  real*8  :: N1(1:3),N2(1:3)
  real*8  :: resSq
  integer :: i,j

  LengthForArbtrySgmnt = 0.0d0
  
  do i = 1,(Nnode_Spr-1)
     N1(1:3) = Nodes(i,1:3)
     N2(1:3) = Nodes((i+1),1:3)
     
     resSq = 0.0d0
     
     do j = 1,3
        resSq = resSq + (N2(j) - N1(j))**2
     enddo

     LengthForArbtrySgmnt = LengthForArbtrySgmnt + sqrt(resSq)
  enddo
  
end function LengthForArbtrySgmnt

subroutine get_Neutral_Axis(cell_nm,Neutrl)
  use system_parameters
  use transfrm_info
  
  implicit none
  integer, intent(in)  :: cell_nm
  real*8,  intent(out) :: Neutrl(1:N_dmnsn)
  
  integer :: TopN1,TopN2,BotN1,BotN2
  real*8  :: Tx,Ty,Bx,By
  
  if (typ_area(cell_nm)==1) then

     TopN1 = area_node(cell_nm,1)
     TopN2 = area_node(cell_nm,2)
     BotN1 = area_node(cell_nm,3)
     BotN2 = area_node(cell_nm,4)
     
     Tx = (node_xy(TopN1,1)+node_xy(TopN2,1))/2.d0
     Ty = (node_xy(TopN1,2)+node_xy(TopN2,2))/2.d0
     Bx = (node_xy(BotN1,1)+node_xy(BotN2,1))/2.d0
     By = (node_xy(BotN1,1)+node_xy(BotN2,1))/2.d0
     
     Neutrl(1) = (Tx+Bx)/2.d0
     Neutrl(2) = (Ty+By)/2.d0
     
  elseif (typ_area(cell_nm)==2) then
     continue
  endif
  
end subroutine get_Neutral_Axis

subroutine get_sprMid(spr_nm,sprMid)
  use system_parameters
  use transfrm_info
  
  implicit none
  integer, intent(in)  :: spr_nm
  real*8 , intent(out) :: sprMid(1:N_dmnsn)

  integer :: sprN1,sprN2
  
  sprN1 = spr_node(spr_nm,1)
  sprN2 = spr_node(spr_nm,2)
  
  if (typ_spr(spr_nm).ne.6) then
     sprMid(1) = (node_xy(sprN1,1)+node_xy(sprN2,1))/2.d0
     sprMid(2) = (node_xy(sprN1,2)+node_xy(sprN2,2))/2.d0
  elseif (typ_spr(spr_nm)==6) then
     continue
  endif

end subroutine get_sprMid

subroutine Optimize_system
  use calls_for_tests
  implicit none
  integer :: iter
  real*8  :: fret,fpp
  
  fret = 0.0d0    ; fpp = 0.0d0

  funcChoice = 2 ; N_variabls = N_mvCoordnte_withl0A0 ; ftol=1.d-03
  call frprmn(coordntes,grd_mv,iter,funcChoice,N_variabls,ftol,fret,fpp)
  call coordntes_to_nodesl0A0(coordntes,node_xy,l0,A0)
  coordntes_xy(1:N_mvCoordnte) = coordntes(1:N_mvCoordnte)
  grdmv_xy(1:N_mvCoordnte)     = grd_mv(1:N_mvCoordnte)

  write(*,*) iter,"iter taken"
  
end subroutine Optimize_system


subroutine Equilibrate_system
  use calls_for_tests
  implicit none
  integer :: iter
  real*8  :: fret,fpp
  integer :: i1
  
  
  if (EquilAlgrthm == 1) then
     
     open(unit=74,file='CG_Energy.dat',position='append') ; fret = 0.0d0  ; fpp = 0.0d0
     funcChoice = 1 ; N_variabls = N_mvCoordnte ; ftol=1.d-06    
     call frprmn(coordntes_xy,grdmv_xy,iter,funcChoice,N_variabls,ftol,fret,fpp)
     call coordntes_to_nodes(coordntes_xy,node_xy)
     write(74,fmt=*) fpp,fret,iter
     close(74)
     
     coordntes(1:N_mvCoordnte) = coordntes_xy(1:N_mvCoordnte)
     grd_mv(1:N_mvCoordnte)    = grdmv_xy(1:N_mvCoordnte)
     write(*,*) iter,"iter taken"
     
  elseif (EquilAlgrthm == 2) then
     funcChoice = 1 ; N_variabls = N_mvCoordnte ; ftol=1.0d-06 
     call TimeStepAlg(coordntes_xy,grdmv_xy,funcChoice,N_variabls,ftol)
     
     call coordntes_to_nodes(coordntes_xy,node_xy)
     coordntes(1:N_mvCoordnte) = coordntes_xy(1:N_mvCoordnte)
     grd_mv(1:N_mvCoordnte)    = grdmv_xy(1:N_mvCoordnte)
  endif
     
end subroutine Equilibrate_system


subroutine Equilibrate_and_Optimize_system
  use calls_for_tests
  implicit none
  integer :: iter
  real*8  :: fret,fpp
  
  if (EquilAlgrthm == 1) then
     fret = 0.0d0 ; fpp = 0.0d0
     
     funcChoice = 1 ; N_variabls = N_mvCoordnte ; ftol=1.d-06    
     call frprmn(coordntes_xy,grdmv_xy,iter,funcChoice,N_variabls,ftol,fret,fpp)
     call coordntes_to_nodes(coordntes_xy,node_xy)
     coordntes(1:N_mvCoordnte) = coordntes_xy(1:N_mvCoordnte)
     grd_mv(1:N_mvCoordnte)    = grdmv_xy(1:N_mvCoordnte)
     
     
     funcChoice = 2 ; N_variabls = N_mvCoordnte_withl0A0 ; ftol=1.d-03
     call frprmn(coordntes,grd_mv,iter,funcChoice,N_variabls,ftol,fret,fpp)
     call coordntes_to_nodesl0A0(coordntes,node_xy,l0,A0)
     coordntes_xy(1:N_mvCoordnte) = coordntes(1:N_mvCoordnte)
     grdmv_xy(1:N_mvCoordnte)     = grd_mv(1:N_mvCoordnte)
     
     write(*,*) iter,"iter taken" 
     
     
  elseif (EquilAlgrthm == 2) then
     
     fret = 0.0d0 ; fpp = 0.0d0
     
     funcChoice = 1 ; N_variabls = N_mvCoordnte ; ftol=1.d-06    
     
     call coordntes_to_nodes(coordntes_xy,node_xy)
     coordntes(1:N_mvCoordnte) = coordntes_xy(1:N_mvCoordnte)
     grd_mv(1:N_mvCoordnte)    = grdmv_xy(1:N_mvCoordnte)
     call TimeStepAlg(coordntes_xy,grdmv_xy,funcChoice,N_variabls,ftol)
  
     funcChoice = 2 ; N_variabls = N_mvCoordnte_withl0A0 ; ftol=1.d-03
     call TimeStepAlg(coordntes_xy,grdmv_xy,funcChoice,N_variabls,ftol)
     
     call coordntes_to_nodesl0A0(coordntes,node_xy,l0,A0)
     coordntes_xy(1:N_mvCoordnte) = coordntes(1:N_mvCoordnte)
     grdmv_xy(1:N_mvCoordnte)     = grd_mv(1:N_mvCoordnte)
     
  endif
  
end subroutine Equilibrate_and_Optimize_system



subroutine get_neighbour_nodes(area_nmbr,node_nmbr,neighbour_nodes)
  use system_parameters
  use transfrm_info
  
  implicit none
  integer, intent(in)  :: area_nmbr,node_nmbr
  integer, intent(out) :: neighbour_nodes(1:2)
  
  integer :: i
  logical :: lgcl_nn !nn=neighbouring_node
  integer :: other_dn
  
  lgcl_nn = .False.
  neighbour_nodes(1:2) = -1 !unrealistic value
  !write(*,*) node_nmbr,area_nmbr,"node,area"
  
  do i = 1,max_node_area
     if (node_nmbr .eq. (-1)) then
        write(*,*) "stop reasn: nodes cant be negtve,Flnm:Fundamental_routines"
        stop
     else
        if (node_nmbr .eq. area_node(area_nmbr,i)) then
           call get_first_neighbour(i)
           call get_second_neighbour(i)
           exit
        else
           if (i.eq.max_node_area) then
              lgcl_nn = .True.
              !write(*,*) lgcl_nn,"lgcl_nn"
              call get_the_other_dn(node_nmbr,other_dn)
              !write(*,*) node_nmbr,other_dn,"node,other"
              
           endif
           
        endif
     endif
  enddo
  
  if (lgcl_nn .eqv. .True.) then
     do i = 1,max_node_area
        if (other_dn .eq. area_node(area_nmbr,i)) then
           call get_first_neighbour(i)
           call get_second_neighbour(i)
        endif
     enddo
  endif
  
  
  !write(*,*) node_nmbr,neighbour_nodes(1:2),"node_and_neighbours"
  
contains
  
  subroutine get_first_neighbour(node_pstn_in_array)
    implicit none
    integer :: node_pstn_in_array
    
    if (modelID==1) then
       
       if ((node_pstn_in_array+1) .le. max_node_area) then
          
          if (area_node(area_nmbr,(node_pstn_in_array+1)).ne.(-1)) then
             neighbour_nodes(1) = area_node(area_nmbr,(node_pstn_in_array+1))
          else
             neighbour_nodes(1) = area_node(area_nmbr,1) 
          endif
       else
          neighbour_nodes(1) = area_node(area_nmbr,1)
       endif
       
    elseif (modelID==2) then

       !if (area_nmbr.le.(ncl+ncr+2)) then
          
          if ((node_pstn_in_array+1) .le. max_node_area) then
             
             if (area_node(area_nmbr,(node_pstn_in_array+1)).ne.(-1)) then
                neighbour_nodes(1) = area_node(area_nmbr,(node_pstn_in_array+1))
             else
                neighbour_nodes(1) = area_node(area_nmbr,1) 
             endif
             
          else
             neighbour_nodes(1) = area_node(area_nmbr,1)
          endif
          
       !elseif (area_nmbr.gt.(ncl+ncr+2)) then
          
          !if ((node_pstn_in_array+1) .le. max_node_areaS) then
             
           !  if (area_node(area_nmbr,(node_pstn_in_array+1)).ne.(-1)) then
             !  neighbour_nodes(1) = area_node(area_nmbr,(node_pstn_in_array+1))
            ! else
              !  neighbour_nodes(1) = area_node(area_nmbr,1) 
             !endif
             
          !else
           !  neighbour_nodes(1) = area_node(area_nmbr,1)
          !endif
          
       !endif
       
    endif
    
  end subroutine get_first_neighbour
  
  subroutine get_second_neighbour(node_pstn_in_array)
    implicit none
    integer :: node_pstn_in_array
    integer :: u,cnt_lpInSN,neighNodeVal
    
    if (modelID==1) then
       
       if ((node_pstn_in_array-1) .ne. 0) then
          neighbour_nodes(2) = area_node(area_nmbr,(node_pstn_in_array-1))  
       else
          
          !if(area_node(area_nmbr,max_node_area) .ne. (-1)) then
          !   neighbour_nodes(2) = area_node(area_nmbr,max_node_area)
          !else
          !   if (area_node(area_nmbr,0).ne.3) then
          !      neighbour_nodes(2) = area_node(area_nmbr,(max_node_area-1))
          !   elseif (area_node(area_nmbr,0)==3) then !for triangles
          !      neighbour_nodes(2) = area_node(area_nmbr,(max_node_area-2))
          !   endif
          !endif
          
          !This endless do will help finding out non (-1) node from area_nodes's
          !highest position to reverse direction
          
          cnt_lpInSN   = 0
          do
             neighNodeVal = area_node(area_nmbr,(max_node_area-cnt_lpInSN))
             
             if (neighNodeVal==-1) then
                cnt_lpInSN = cnt_lpInSN+1
             elseif (neighNodeVal.ne.-1) then
                exit
             endif
          enddo
          neighbour_nodes(2) = neighNodeVal
       endif
       
    elseif (modelID==2) then
       
       if ((node_pstn_in_array-1) .ne. 0) then
          neighbour_nodes(2) = area_node(area_nmbr,(node_pstn_in_array-1))  
       else
          
          !if (area_nmbr.le.(ncl+ncr+2)) then
             
             if(area_node(area_nmbr,max_node_area) .ne. (-1)) then
                neighbour_nodes(2) = area_node(area_nmbr,max_node_area)
             else
                
                do u = (max_node_area-1),1,(-1)
                   if (area_node(area_nmbr,u) .ne. (-1)) then
                      neighbour_nodes(2) = area_node(area_nmbr,u)
                      exit
                   endif
                enddo
                
             endif
             
          ! elseif (area_nmbr.gt.(ncl+ncr+2)) then
             
          !    if(area_node(area_nmbr,max_node_areaS) .ne. (-1)) then
          !       neighbour_nodes(2) = area_node(area_nmbr,max_node_areaS)
          !    else
          !       neighbour_nodes(2) = area_node(area_nmbr,(max_node_areaS-1))
          !    endif
             
          ! endif

       endif
       
    endif
    
  end subroutine get_second_neighbour
  
end subroutine get_neighbour_nodes


real*8 function dot_prdct(vectA,vectB,dmnsn)
  implicit none
  integer, intent(in) :: dmnsn
  real*8,  intent(in) :: vectA(1:dmnsn),vectB(1:dmnsn)
  integer :: i
  
  dot_prdct = 0.0d0
  
  do i = 1,dmnsn
     dot_prdct = dot_prdct + vectA(i)*vectB(i)
  enddo
  
end function dot_prdct

real*8 function Mgnitude(vectA,dmnsn)
  implicit none
  integer,intent(in) :: dmnsn
  real*8 ,intent(in) :: vectA(1:dmnsn)
  integer :: i
  real*8  :: MgnitudeSq
  
  Mgnitude = 0.0d0 ; MgnitudeSq = 0.0d0
  
  do i = 1,dmnsn
     MgnitudeSq = MgnitudeSq + vectA(i)*vectA(i)
  enddo
  
  Mgnitude = sqrt(MgnitudeSq)
  
end function Mgnitude


subroutine get_Ranking(inpt_array,N_Elmnt,Ranking)
  implicit none
  integer, intent(in)  :: N_Elmnt
  real*8,  intent(in)  :: inpt_array(1:N_Elmnt)
  integer, intent(out) :: Ranking(1:N_Elmnt)
  
  integer :: i,j
  
  open(unit=164,file='Ranking.dat',position='append')
  
  do i = 1,N_Elmnt
     Ranking(i) = 1
     
     do j = 1,N_Elmnt
        if (j==i) cycle
        
        if(inpt_array(i).gt.inpt_array(j)) then
           continue
        elseif(inpt_array(i).lt.inpt_array(j)) then
           Ranking(i) = Ranking(i) + 1
        else
           write(*,*) inpt_array(i),inpt_array(j),"inpt_array_i,j"
           write(*,*) "cant be equal,Sb:get_rnking,fl:Fundamental_routines"
           stop
        endif
        
     enddo
     
     write(164,fmt=*) Ranking(i),i,"Ranking,ElmntNo" 
     
  enddo
  
  close(164)
  
end subroutine get_Ranking

subroutine get_cubic_fit(x1y1,x2y2,x3y3,N_lon,LoN)
  implicit none
  real*8,  intent(in)  :: x1y1(1:2),x2y2(1:2),x3y3(1:2)
  integer, intent(in)  :: N_lon
  real*8,  intent(out) :: LoN(1:N_lon,1:2) !LoN=List of Nodes
  
  integer :: i
  integer :: caseNo
  real*8  :: x1,y1,x2,y2,x3,y3
  real*8  :: P1,P2,Q1,Q2,R1,R2
  real*8  :: A,B,C
  real*8  :: LoNX(1:N_lon),LoNY(1:N_lon)
  real*8  :: TOL
  
  interface
     real*8 function funcLonX(Y,A,B,C)
       implicit none
       real*8 :: Y,A,B,C
     end function funcLonX
     
     real*8 function funcLonY(X,A,B,C)
       implicit none
       real*8 :: X,A,B,C
     end function funcLonY
  end interface
  
  TOL = 1.0d-15
  
  x1 = x1y1(1) ; y1 = x1y1(2)
  x2 = x2y2(1) ; y2 = x2y2(2)
  x3 = x3y3(1) ; y3 = x3y3(2)

  if (abs(x1-x3).gt.abs(y1-y3)) then
     caseNo = 1
  elseif (abs(x1-x3).le.abs(y1-y3)) then
     caseNo = 2
  endif
  
  call get_PQR(caseNo,x1,y1,x2,y2,x3,y3,P1,P2,Q1,Q2,R1,R2)
  call get_ABval(P1,Q1,R1,P2,Q2,R2,A,B)
  call get_Cval(caseNo,x3,y3,A,B,C)
  
  if (caseNo==1) then
     call get_LoNX(x1,x2,x3,N_lon,LoNX)
     LoN(1:N_lon,1) = LoNX(1:N_lon) 
     
     do i = 1,N_lon
        if (i==1) then
           LoN(i,2) = y1
        elseif (i==((N_lon/2)+1)) then
           LoN(i,2) = y2
        elseif (i==N_lon) then
           LoN(i,2) = y3
        else
           LoN(i,2) = funcLonY(Lon(i,1),A,B,C)
        endif
        !write(*,*) LoN(i,1),LoN(i,2),i,"LoNX,LoNY"
     enddo
     
     
  elseif (caseNo==2) then
     call get_LonY(y1,y2,y3,N_lon,LoNY)
     LoN(1:N_lon,2) = LoNY(1:N_lon)
     
     do i = 1,N_lon
        LoN(i,1) = funcLonX(Lon(i,2),A,B,C) 
     enddo
     
  endif
     
  ! elseif ((abs(x1-x3).le.TOL) .and. (abs(y1-y3).gt.TOL)) then !x1=x3,y1=/y3
     
  !    call get_ABval_SpclCase(x2,y2,A,B)
  !    call get_LoNY(y1,y2,y3,N_lon,LoNY)
  !    call get_LonXspcl(LonY,A,B,LonX)
     
  ! elseif ((abs(x1-x3).gt.TOL) .and. (abs(y1-y3).le.TOL)) then !x1=/x3,y1=y3
     
 
  !    call get_LonX(x1,x2,x3,N_lon,LoNX)
  !    call get_LonYspcl(LonX,A,B,LonY)
     
  !endif
  
end subroutine get_cubic_fit

subroutine get_PQR(caseNo,x1,y1,x2,y2,x3,y3,P1,P2,Q1,Q2,R1,R2)
  implicit none
  integer, intent(in)  :: caseNo
  real*8,  intent(in)  :: x1,y1,x2,y2,x3,y3
  real*8,  intent(out) :: P1,P2,Q1,Q2,R1,R2
  
  if (caseNo==1) then
     P1=(x1**2-x2**2) ; Q1=(x1-x2) ; R1 = -(y1-y2)
     P2=(x2**2-x3**2) ; Q2=(x2-x3) ; R2 = -(y2-y3)
     
  elseif (caseNo==2) then
     P1=(y1**2-y2**2) ; Q1=(y1-y2) ; R1 = -(x1-x2)
     P2=(y2**2-y3**2) ; Q2=(y2-y3) ; R2 = -(x2-x3)
     
  endif
  
end subroutine get_PQR

subroutine get_ABval(P1,Q1,R1,P2,Q2,R2,A,B)
  implicit none
  real*8, intent(in)  :: P1,Q1,R1,P2,Q2,R2
  real*8, intent(out) :: A,B

  A = (Q1*R2-Q2*R1) / (P1*Q2-P2*Q1)
  B = (R1*P2-R2*P1) / (P1*Q2-P2*Q1)
  
end subroutine get_ABval

subroutine get_Cval(caseNo,x3,y3,A,B,C)
  implicit none
  integer, intent(in)  :: caseNo
  real*8 , intent(in)  :: x3,y3
  real*8 , intent(in)  :: A,B
  real*8 , intent(out) :: C

  if (caseNo==1) then
     C = y3-A*(x3)**2-B*x3
  elseif (caseNo==2) then
     C = x3-A*(y3)**2-B*y3
  endif
  
end subroutine get_Cval


subroutine get_LonX(x1,x2,x3,N_lon,LoNX)
  implicit none
  real*8, intent(in)  :: x1,x2,x3
  integer,intent(in)  :: N_lon
  real*8, intent(out) :: LoNX(1:N_lon)

  integer :: intrmd_Nlon
  real*8  :: frst_xrnge,scnd_xrnge
  integer :: i,cnt
  real*8  :: dx
  
  intrmd_Nlon = (N_lon-3)/2
  frst_xrnge  = (x2-x1)
  scnd_xrnge  = (x3-x2)
  
  dx = frst_xrnge/real(intrmd_Nlon+1)
  
  LonX(1) = x1
  cnt = 1
  
  do i = 1,intrmd_Nlon
     cnt = cnt+1
     LonX(cnt) = x1 + real(i)*dx
  enddo
  
  cnt = cnt+1
  LoNX(cnt) = x2
  dx = scnd_xrnge/real(intrmd_Nlon+1)
  
  do i = 1,intrmd_Nlon
     cnt = cnt+1
     LonX(cnt) = x2 + real(i)*dx
  enddo
  
  cnt = cnt+1
  LoNX(cnt) = x3
  
end subroutine get_LonX


subroutine get_LonY(y1,y2,y3,N_lon,LoNY)
  implicit none
  real*8, intent(in)  :: y1,y2,y3
  integer,intent(in)  :: N_lon
  real*8, intent(out) :: LoNY(1:N_lon)

  integer :: intrmd_Nlon
  real*8  :: frst_yrnge,scnd_yrnge
  integer :: i,cnt
  real*8  :: dy
  
  intrmd_Nlon = (N_lon-3)/2
  frst_yrnge  = (y2-y1)
  scnd_yrnge  = (y3-y2)
  
  dy = frst_yrnge/real(intrmd_Nlon+1)

  LonY(1) = y1
  cnt = 1
  
  do i = 1,intrmd_Nlon
     cnt = cnt+1
     LonY(cnt) = y1 + real(i)*dy
  enddo

  cnt = cnt+1
  LoNY(cnt) = y2
  dy = scnd_yrnge/real(intrmd_Nlon+1)
  
  do i = 1,intrmd_Nlon
     cnt = cnt+1
     LonY(cnt) = y2 + real(i)*dy
  enddo
  
  cnt = cnt+1
  LoNY(cnt) = y3
  
end subroutine get_LonY


real*8 function funcLonX(Y,A,B,C)
  implicit none
  real*8 :: Y,A,B,C
  
  funcLoNX = A*(Y**2) + B*Y + C
  
end function funcLonX

real*8 function funcLonY(X,A,B,C)
  implicit none
  real*8 :: X,A,B,C
  
  funcLoNY = A*(X**2) + B*X + C
  
end function funcLonY



function MAT_INV4b4(A1) !! Performs a direct calculation of the inverse of a 4×4 matrix.
  implicit none
  real*8, intent(in)  :: A1(4,4)            !! Matrix
  real*8              :: detinv,detrmnt
  real*8              :: ZeroV = 0.00000000000000001d0
  
  real*8              :: MAT_INV4b4(4,4) 
  
  ! Calculate the inverse determinant of the matrix
  
  detrmnt = (A1(1,1)*(A1(2,2)*(A1(3,3)*A1(4,4)-A1(3,4)*A1(4,3))+A1(2,3)*(A1(3,4)*A1(4,2)-A1(3,2)*A1(4,4))+ &
       A1(2,4)*(A1(3,2)*A1(4,3)-A1(3,3)*A1(4,2))) &
       - A1(1,2)*(A1(2,1)*(A1(3,3)*A1(4,4)-A1(3,4)*A1(4,3))+A1(2,3)*(A1(3,4)*A1(4,1)-A1(3,1)*A1(4,4))+ &
       A1(2,4)*(A1(3,1)*A1(4,3)-A1(3,3)*A1(4,1)))&
       + A1(1,3)*(A1(2,1)*(A1(3,2)*A1(4,4)-A1(3,4)*A1(4,2))+A1(2,2)*(A1(3,4)*A1(4,1)-A1(3,1)*A1(4,4))+ &
       A1(2,4)*(A1(3,1)*A1(4,2)-A1(3,2)*A1(4,1)))&
       - A1(1,4)*(A1(2,1)*(A1(3,2)*A1(4,3)-A1(3,3)*A1(4,2))+A1(2,2)*(A1(3,3)*A1(4,1)-A1(3,1)*A1(4,3))+ &
       A1(2,3)*(A1(3,1)*A1(4,2)-A1(3,2)*A1(4,1))))
  
  write(*,*) detrmnt,"det Val"
  if (abs(detrmnt) .le. ZeroV) stop 'inconsistent determinant' 
  
  detinv  = 1/detrmnt
  
  ! Calculate the inverse of the matrix
  MAT_INV4b4(1,1) = detinv*(A1(2,2)*(A1(3,3)*A1(4,4)-A1(3,4)*A1(4,3))+A1(2,3)*(A1(3,4)*A1(4,2)- &
       A1(3,2)*A1(4,4))+A1(2,4)*(A1(3,2)*A1(4,3)-A1(3,3)*A1(4,2)))
  MAT_INV4b4(2,1) = detinv*(A1(2,1)*(A1(3,4)*A1(4,3)-A1(3,3)*A1(4,4))+A1(2,3)*(A1(3,1)*A1(4,4)- &
       A1(3,4)*A1(4,1))+A1(2,4)*(A1(3,3)*A1(4,1)-A1(3,1)*A1(4,3)))
  MAT_INV4b4(3,1) = detinv*(A1(2,1)*(A1(3,2)*A1(4,4)-A1(3,4)*A1(4,2))+A1(2,2)*(A1(3,4)*A1(4,1)- &
       A1(3,1)*A1(4,4))+A1(2,4)*(A1(3,1)*A1(4,2)-A1(3,2)*A1(4,1)))
  MAT_INV4b4(4,1) = detinv*(A1(2,1)*(A1(3,3)*A1(4,2)-A1(3,2)*A1(4,3))+A1(2,2)*(A1(3,1)*A1(4,3)- &
       A1(3,3)*A1(4,1))+A1(2,3)*(A1(3,2)*A1(4,1)-A1(3,1)*A1(4,2)))
  MAT_INV4b4(1,2) = detinv*(A1(1,2)*(A1(3,4)*A1(4,3)-A1(3,3)*A1(4,4))+A1(1,3)*(A1(3,2)*A1(4,4)- &
       A1(3,4)*A1(4,2))+A1(1,4)*(A1(3,3)*A1(4,2)-A1(3,2)*A1(4,3)))
  MAT_INV4b4(2,2) = detinv*(A1(1,1)*(A1(3,3)*A1(4,4)-A1(3,4)*A1(4,3))+A1(1,3)*(A1(3,4)*A1(4,1)- &
       A1(3,1)*A1(4,4))+A1(1,4)*(A1(3,1)*A1(4,3)-A1(3,3)*A1(4,1)))
  MAT_INV4b4(3,2) = detinv*(A1(1,1)*(A1(3,4)*A1(4,2)-A1(3,2)*A1(4,4))+A1(1,2)*(A1(3,1)*A1(4,4)- &
       A1(3,4)*A1(4,1))+A1(1,4)*(A1(3,2)*A1(4,1)-A1(3,1)*A1(4,2)))
  MAT_INV4b4(4,2) = detinv*(A1(1,1)*(A1(3,2)*A1(4,3)-A1(3,3)*A1(4,2))+A1(1,2)*(A1(3,3)*A1(4,1)- &
       A1(3,1)*A1(4,3))+A1(1,3)*(A1(3,1)*A1(4,2)-A1(3,2)*A1(4,1)))
  MAT_INV4b4(1,3) = detinv*(A1(1,2)*(A1(2,3)*A1(4,4)-A1(2,4)*A1(4,3))+A1(1,3)*(A1(2,4)*A1(4,2)- &
       A1(2,2)*A1(4,4))+A1(1,4)*(A1(2,2)*A1(4,3)-A1(2,3)*A1(4,2)))
  MAT_INV4b4(2,3) = detinv*(A1(1,1)*(A1(2,4)*A1(4,3)-A1(2,3)*A1(4,4))+A1(1,3)*(A1(2,1)*A1(4,4)- &
       A1(2,4)*A1(4,1))+A1(1,4)*(A1(2,3)*A1(4,1)-A1(2,1)*A1(4,3)))
  MAT_INV4b4(3,3) = detinv*(A1(1,1)*(A1(2,2)*A1(4,4)-A1(2,4)*A1(4,2))+A1(1,2)*(A1(2,4)*A1(4,1)- &
       A1(2,1)*A1(4,4))+A1(1,4)*(A1(2,1)*A1(4,2)-A1(2,2)*A1(4,1)))
  MAT_INV4b4(4,3) = detinv*(A1(1,1)*(A1(2,3)*A1(4,2)-A1(2,2)*A1(4,3))+A1(1,2)*(A1(2,1)*A1(4,3)- &
       A1(2,3)*A1(4,1))+A1(1,3)*(A1(2,2)*A1(4,1)-A1(2,1)*A1(4,2)))
  MAT_INV4b4(1,4) = detinv*(A1(1,2)*(A1(2,4)*A1(3,3)-A1(2,3)*A1(3,4))+A1(1,3)*(A1(2,2)*A1(3,4)- &
       A1(2,4)*A1(3,2))+A1(1,4)*(A1(2,3)*A1(3,2)-A1(2,2)*A1(3,3)))
  MAT_INV4b4(2,4) = detinv*(A1(1,1)*(A1(2,3)*A1(3,4)-A1(2,4)*A1(3,3))+A1(1,3)*(A1(2,4)*A1(3,1)- &
       A1(2,1)*A1(3,4))+A1(1,4)*(A1(2,1)*A1(3,3)-A1(2,3)*A1(3,1)))
  MAT_INV4b4(3,4) = detinv*(A1(1,1)*(A1(2,4)*A1(3,2)-A1(2,2)*A1(3,4))+A1(1,2)*(A1(2,1)*A1(3,4)- &
       A1(2,4)*A1(3,1))+A1(1,4)*(A1(2,2)*A1(3,1)-A1(2,1)*A1(3,2)))
  MAT_INV4b4(4,4) = detinv*(A1(1,1)*(A1(2,2)*A1(3,3)-A1(2,3)*A1(3,2))+A1(1,2)*(A1(2,3)*A1(3,1)- &
       A1(2,1)*A1(3,3))+A1(1,3)*(A1(2,1)*A1(3,2)-A1(2,2)*A1(3,1)))
  
end function MAT_INV4b4

subroutine find_the_numOfLines
  implicit none
  
  
end subroutine find_the_numOfLines
