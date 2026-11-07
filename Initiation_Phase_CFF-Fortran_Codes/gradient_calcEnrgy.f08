module gradient_from_total_system  !M1
  use system_parameters
  
  implicit none
  integer :: Analtcl_or_Numrcl
  
  real*8, allocatable :: gradient(:,:),     gradientS(:,:),     gradientS2(:,:),     gradientSNI(:,:)
  real*8, allocatable :: gradient_ana(:,:), gradient_anaS(:,:), gradient_anaS2(:,:), gradient_anaSNI(:,:)
  real*8, allocatable :: gradient_num(:,:), gradient_numS(:,:), gradient_numS2(:,:), gradient_numSNI(:,:)
  real*8, allocatable :: grd_mv(:),         grd_mvS(:),         grd_mvS2(:),         grd_mvSNI(:)
  real*8, allocatable :: grdmv_xy(:),       grdmv_xyS(:),       grdmv_xyS2(:),       grdmv_xySNI(:)
  
  !real*8, allocatable :: l0Vartn(:),A0Vartn(:)
  
contains
  
  subroutine allocate_and_initialize_gradient_system_variables
    implicit none
    
    allocate(gradient(1:N_node,1:N_dmnsn),    gradientS(1:N_node,1:N_dmnsn))    
    allocate(gradient_ana(1:N_node,1:N_dmnsn),gradient_anaS(1:N_node,1:N_dmnsn))
    allocate(gradient_num(1:N_node,1:N_dmnsn),gradient_numS(1:N_node,1:N_dmnsn))
    
    !write(*,*) N_mvCoordnte,"N_mvCoordnte inside alloc"
    
    allocate(grd_mv(1:N_mvCoordnte_withl0A0),grd_mvS(1:N_mvCoordnte_withl0A0))
    allocate(grdmv_xy(1:N_mvCoordnte),grdmv_xyS(1:N_mvCoordnte))
    !allocate(l0Vartn(1:N_spr),A0Vartn(1:N_cell)) !may have to use lt N_spr,N_cell, say Nspr_mv,Ncell_mv
    
    gradient     = 10.d5 ; gradientS     = 10.d5
    gradient_ana = 10.d5 ; gradient_anaS = 10.d5
    gradient_num = 10.d5 ; gradient_numS = 10.d5 
    grd_mv       = 10.d5 ; grd_mvS       = 10.d5  !very high value(unrealistic)
    grdmv_xy     = 10.d5 ; grdmv_xyS     = 10.d5
    
    !l0Vartn  = 10e5
    !A0Vartn  = 10e5
    
  end subroutine allocate_and_initialize_gradient_system_variables
  
  subroutine allocate_and_initialize_gradient_system_variables_wo_StrVars
    implicit none
    
    allocate(gradient(1:N_node,1:N_dmnsn))    
    allocate(gradient_ana(1:N_node,1:N_dmnsn))
    allocate(gradient_num(1:N_node,1:N_dmnsn))
    
    allocate(grd_mv(1:N_mvCoordnte_withl0A0))
    allocate(grdmv_xy(1:N_mvCoordnte))
    
    gradient     = 10.d5
    gradient_ana = 10.d5
    gradient_num = 10.d5
    grd_mv       = 10.d5 !very high value(unrealistic)
    grdmv_xy     = 10.d5
    
  end subroutine allocate_and_initialize_gradient_system_variables_wo_StrVars
  
  subroutine allocate_and_initialize_gradient_system_variables_StrVars
    implicit none
    
    allocate(gradientS(1:N_node,1:N_dmnsn))    
    allocate(gradient_anaS(1:N_node,1:N_dmnsn))
    allocate(gradient_numS(1:N_node,1:N_dmnsn))
    
    allocate(grd_mvS(1:N_mvCoordnte_withl0A0))
    allocate(grdmv_xyS(1:N_mvCoordnte))
    
    gradient     = 10.d5
    gradient_ana = 10.d5
    gradient_num = 10.d5
    grd_mv       = 10.d5 !very high value(unrealistic)
    grdmv_xy     = 10.d5
    
  end subroutine allocate_and_initialize_gradient_system_variables_StrVars
  
  subroutine deallocate_gradient_system_variables
    implicit none
    
    deallocate(gradient,gradientS)    
    deallocate(gradient_ana,gradient_anaS)
    deallocate(gradient_num,gradient_numS)
    deallocate(grd_mv,grd_mvS)
    deallocate(grdmv_xy,grdmv_xyS)
    
  end subroutine deallocate_gradient_system_variables
  
  subroutine deallocate_gradient_system_variables_wo_StrVars
    implicit none
    
    deallocate(gradient)    
    deallocate(gradient_ana)
    deallocate(gradient_num)
    deallocate(grd_mv)
    deallocate(grdmv_xy)
    
  end subroutine deallocate_gradient_system_variables_wo_StrVars
  
  subroutine deallocate_gradient_system_variables_StrVars
    implicit none
    
    deallocate(gradientS)    
    deallocate(gradient_anaS)
    deallocate(gradient_numS)
    deallocate(grd_mvS)
    deallocate(grdmv_xyS)
    
  end subroutine deallocate_gradient_system_variables_StrVars
  
  subroutine store_gradient_system_variables
    implicit none
    
    gradientS     = gradient
    gradient_anaS = gradient_ana
    gradient_numS = gradient_num
    grd_mvS       = grd_mv
    grdmv_xyS     = grdmv_xy
    
  end subroutine store_gradient_system_variables
  
  subroutine restore_gradient_system_variables
    implicit none
    
    gradient     = gradientS
    gradient_ana = gradient_anaS
    gradient_num = gradient_numS
    grd_mv       = grd_mvS
    grdmv_xy     = grdmv_xyS
    
  end subroutine restore_gradient_system_variables
  
  subroutine store_gradient_system_variables_NI
    implicit none
    
    gradientSNI     = gradient
    gradient_anaSNI = gradient_ana
    gradient_numSNI = gradient_num
    grd_mvSNI       = grd_mv
    grdmv_xySNI     = grdmv_xy
    
  end subroutine store_gradient_system_variables_NI
  
  subroutine store_gradient_system_variables_S2
    implicit none
    
    gradientS2     = gradient
    gradient_anaS2 = gradient_ana
    gradient_numS2 = gradient_num
    grd_mvS2       = grd_mv
    grdmv_xyS2     = grdmv_xy
    
  end subroutine store_gradient_system_variables_S2
  
end module gradient_from_total_system




module gradient_contrb_from_spr !M2
  use system_parameters
  use transfrm_info
  use triangle_routines
  implicit none
  
  real*8, allocatable :: grd_spr(:,:),grd_sprS(:,:),grd_sprS2(:,:),grd_sprSNI(:,:)
  real*8, allocatable :: delE_delR(:),delE_delRS(:),delE_delRS2(:),delE_delRSNI(:)
  real*8, allocatable :: delR_delNode(:,:,:),delR_delNodeS(:,:,:),delR_delNodeS2(:,:,:),delR_delNodeSNI(:,:,:)
  
contains
  
  subroutine allocate_and_initialize_gradient_spr_variables
    implicit none
    
    allocate(grd_spr(1:N_node,1:N_dmnsn),grd_sprS(1:N_node,1:N_dmnsn))
    allocate(delE_delR(1:N_spr),delE_delRS(1:N_spr))
    allocate(delR_delNode(1:N_node,1:max_spr_node,1:2),delR_delNodeS(1:N_node,1:max_spr_node,1:2))
    
    grd_spr      = -1.d30 ; grd_sprS      = -1.d30
    delE_delR    = -1.d30 ; delE_delRS    = -1.d30
    delR_delNode = -1.d30 ; delR_delNodeS = -1.d30
    
  end subroutine allocate_and_initialize_gradient_spr_variables
  
  subroutine allocate_and_initialize_gradient_spr_variables_wo_StrVars
    implicit none
    
    allocate(grd_spr(1:N_node,1:N_dmnsn))
    allocate(delE_delR(1:N_spr))
    allocate(delR_delNode(1:N_node,1:max_spr_node,1:2))
    
    grd_spr      = -1.d30
    delE_delR    = -1.d30
    delR_delNode = -1.d30
    
  end subroutine allocate_and_initialize_gradient_spr_variables_wo_StrVars
  
  subroutine allocate_and_initialize_gradient_spr_variables_StrVars
    implicit none
    
    allocate(grd_sprS(1:N_node,1:N_dmnsn))
    allocate(delE_delRS(1:N_spr))
    allocate(delR_delNodeS(1:N_node,1:max_spr_node,1:2))
    
    grd_sprS      = -1.d30
    delE_delRS    = -1.d30
    delR_delNodeS = -1.d30
    
  end subroutine allocate_and_initialize_gradient_spr_variables_StrVars
  
  
  subroutine deallocate_gradient_spr_variables
    implicit none
    
    deallocate(grd_spr,grd_sprS)
    deallocate(delE_delR,delE_delRS)
    deallocate(delR_delNode,delR_delNodeS)
    
  end subroutine deallocate_gradient_spr_variables
  
  
  subroutine deallocate_gradient_spr_variables_wo_StrVars
    implicit none
    
    deallocate(grd_spr)
    deallocate(delE_delR)
    deallocate(delR_delNode)
    
  end subroutine deallocate_gradient_spr_variables_wo_StrVars
  
  subroutine deallocate_gradient_spr_variables_StrVars
    implicit none
    
    deallocate(grd_sprS)
    deallocate(delE_delRS)
    deallocate(delR_delNodeS)
    
  end subroutine deallocate_gradient_spr_variables_StrVars
  
  subroutine store_gradient_spr_variables
    implicit none
    
    grd_sprS      = grd_spr
    delE_delRS    = delE_delR
    delR_delNodeS = delR_delNode
    
  end subroutine store_gradient_spr_variables
  
  subroutine restore_gradient_spr_variables
    implicit none
    
    grd_spr      = grd_sprS
    delE_delR    = delE_delRS
    delR_delNode = delR_delNodeS
    
  end subroutine restore_gradient_spr_variables
  
  subroutine store_gradient_spr_variables_NI
    implicit none
    
    grd_sprSNI      = grd_spr
    delE_delRSNI    = delE_delR
    delR_delNodeSNI = delR_delNode
    
  end subroutine store_gradient_spr_variables_NI
  
  subroutine store_gradient_spr_variables_S2
    implicit none
    
    grd_sprS2      = grd_spr
    delE_delRS2    = delE_delR
    delR_delNodeS2 = delR_delNode
    
  end subroutine store_gradient_spr_variables_S2
  
  
  subroutine gradient_components_frm_spr(dum_Nodes,duml0,dum_grdSpr)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8  :: duml0(1:N_spr)
    real*8  :: dum_grdSpr(1:N_node,1:N_dmnsn)
    integer :: i
    integer :: spr_nmbr
    integer :: totl_spr_cnnctd
    integer :: given_node,other_node
    
    dum_grdSpr   = 0.0d0
    delE_delR    = 0.0d0
    delR_delNode = 0.0d0
    
    given_node = 0 ; other_node = 0
    
    call calc_delE_delR_array(dum_Nodes,duml0)
    call calc_delR_delNode_array(dum_Nodes)
    
    do i = 1,N_node

       if (node_typ(i) .eq. 0) then
          continue
       elseif (node_typ(i) .ne. 0) then
          totl_spr_cnnctd = node_spr(i,0)
          !write(*,*) totl_spr_cnnctd,"tsc
          
          if (node_cnnctd(i) .eq. 0) then
             
             if (node_typ(i) .eq. 1) then !free node
                call grdSpr_lp(i)
             elseif (node_typ(i) .eq. 2) then !x free, y fixed
                call grdSpr_lp(i)
             endif
             
          elseif (node_cnnctd(i) .ne. 0) then
             
             if (count_this_dn(i) .eq. 1) then
                
                if (node_typ(i) .eq. 1) then !free node
                   call grdSpr_lp(i)
                elseif (node_typ(i) .eq. 2) then !x free, y fixed
                   call grdSpr_lp(i)
                endif
                
                
             elseif (count_this_dn(i) .ne. 1) then
                !will have to replace with other double node value
                !has been replaced after calc
                continue
             endif
             
          endif
          
          !write(*,*) dum_grdSpr(i,1:2),"dum_grdSpr and i = ",i
       endif
       
       if (modelID==2.and.filechk==1) then
          write(*,*) dum_grdSpr(i,1:2),"dum_grdSpr and i = ",i
       endif
       !write(*,*) dum_grdSpr(i,1:2),"dum_grdSpr and i = ",i
       
    enddo
    
    
  contains
    
    subroutine grdSpr_lp(node_nm)
      implicit none
      integer :: j
      integer :: node_nm
      
      !if (node_nm .eq. 14) write(*,*) dum_grdSpr(node_nm,1:2),"dum_grdSpr14_bfr"
      
      do j = 1,totl_spr_cnnctd
         
         spr_nmbr = node_spr(node_nm,j)
         
         if (node_typ(node_nm) .eq. 1) then !x,y free
            dum_grdSpr(node_nm,1:2) = dum_grdSpr(node_nm,1:2) &
                 + delE_delR(spr_nmbr) * delR_delNode(node_nm,j,1:2)
            
            ! if (SystemTyp==1) then
            !    if (node_nm==47 .or. node_nm==48) then
            !       write(*,*) "For node_nm =",node_nm
            !       write(*,*) delE_delR(spr_nmbr),spr_nmbr,j,"delE_delR"
            !       write(*,*) delR_delNode(node_nm,j,1:2),"delR_delN"
            !    endif
            ! endif
            
         
         elseif (node_typ(node_nm) .eq. 2) then !x free, y fixed        
            dum_grdSpr(node_nm,1) = dum_grdSpr(node_nm,1) &
                 + delE_delR(spr_nmbr) * delR_delNode(node_nm,j,1)
            
            ! if (SystemTyp==1) then
            !    if (node_nm==49 .or. node_nm==50) then
            !       write(*,*) "For node_nm =",node_nm
            !       write(*,*) delE_delR(spr_nmbr),spr_nmbr,j,"delE_delR"
            !       write(*,*) delR_delNode(node_nm,j,1:2),"delR_delN"
            !    endif
            ! endif
            
         endif
         
      enddo
      
      !if (node_nm .eq. 14) write(*,*) dum_grdSpr(node_nm,1:2),"dum_grdSpr14_aft"
    end subroutine grdSpr_lp
    
  end subroutine gradient_components_frm_spr
  
  
  
  subroutine calc_delE_delR_array(dum_Nodes,duml0)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8  :: duml0(1:N_spr)
    integer :: i
    
    open(unit=290,file="delE_delR.dat")
    
    do i = 1,N_spr
       delE_delR(i) = indvd_delE_delR(i)
       write(290,fmt=*) delE_delR(i),i,k_spr(i),l(i),l0(i),"delE_delR"
    enddo
    
    close(290)
    
  contains
    
    real*8 function indvd_delE_delR(spr_nmbr)
      implicit none
      integer :: spr_nmbr
      
      indvd_delE_delR = k_spr(spr_nmbr) * dflctn(spr_nmbr)
      
      if (modelID==2.and.filechk==1) then
         write(*,*) k_spr(spr_nmbr),l(spr_nmbr),l0(spr_nmbr),dflctn(spr_nmbr),"kll0_dflc"
      endif
      
    end function indvd_delE_delR
    
    
!!!!dflctn and related subroutines are common for both Energy and gradient,find a suitable way to write once(!HERE)
    
    real*8 function dflctn(spr_nmbr)
      implicit none
      integer :: spr_nmbr
      integer :: Node1,Node2,Node3
      integer :: pulley_typ
      
      integer             :: nodeNum,cnt,nodeVal
      real*8, allocatable :: Nodes(:,:)
      
      pulley_typ = 6
      
      if (spr_node(spr_nmbr,0)==2) then
         
         Node1 = spr_node(spr_nmbr,1)
         Node2 = spr_node(spr_nmbr,2)
         
         l(spr_nmbr) = lngth_u_two_nodes(Node1,Node2)
         dflctn      = l(spr_nmbr) - duml0(spr_nmbr)
         
         if (abs(dflctn) .le. 1.0d-16) then
            dflctn = 0.0d0
         endif
         
      elseif (spr_node(spr_nmbr,0)==3) then
         
         Node1 = spr_node(spr_nmbr,1)
         Node2 = spr_node(spr_nmbr,2)
         Node3 = spr_node(spr_nmbr,3)
         
         l(spr_nmbr) = lngth_u_three_nodes(Node1,Node2,Node3)
         dflctn = l(spr_nmbr) - duml0(spr_nmbr)
         
         if (abs(dflctn) .le. 1.0d-15) then
            dflctn = 0.0d0
         endif
         
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
        
    end function dflctn
    
    real*8 function lngth_u_two_nodes(Node1,Node2)
      implicit none
      integer :: Node1,Node2
      real*8  :: N1(1:N_dmnsn),N2(1:N_dmnsn)
      real*8  :: res_sq
      integer :: i
      
      N1 = dum_Nodes(Node1,1:N_dmnsn)
      N2 = dum_Nodes(Node2,1:N_dmnsn)
      
      lngth_u_two_nodes = 0.0d0
      res_sq = 0.0d0
      
      do i = 1,N_dmnsn
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
         
         lngth_u_three_nodes = lngth_u_three_nodes + abs(sqrt(res_sq))
         
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
    
    
    
  !!!! dflctn related sbrtns are finished here  
      
  end subroutine calc_delE_delR_array
  
  
  
  
  subroutine calc_delR_delNode_array(dum_Nodes)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    
    integer :: i,j
    integer :: totl_spr_cnnctd,cnnctd_sprs_node
    integer :: spr_nmbr,prev_spr_nmbr
    integer :: init_node
    integer :: Pulley_typ
    integer :: cnt_NTN1,cnt_NTN2 ! count for NTN=Non-Terminal Node1 and 2
    
    open(unit=327,file='Track_dR_dN.dat')
    
    Pulley_typ    = 6
    prev_spr_nmbr = 0
    cnt_NTN1      = 0 ; cnt_NTN2=0
    
    do i = 1,N_node
       !write(*,*) i,node_spr(i,1:max_spr_node),node_spr(i,0),"node_spr for node i"
       
       init_node = i
       totl_spr_cnnctd = node_spr(i,0)
       
       !if(i.eq.13) write(*,*) totl_spr_cnnctd,"tsc"
       
       do j = 1,totl_spr_cnnctd
          spr_nmbr = node_spr(i,j)
          
          if (spr_nmbr==(-1)) then
             write(*,*) i,j,node_spr(i,j),"i,j,node_spr(i,j)","negative"
             stop
          endif
          
          if (j == 1) prev_spr_nmbr = 0
          if (j.ne.1) prev_spr_nmbr = node_spr(i,j-1)  
          
          
          if (i==103 .or. i==159) write(327,*) spr_nmbr,prev_spr_nmbr,node_spr(i,j-1),"prev"
          
          if (spr_node(spr_nmbr,0)==2) then   
             delR_delNode(i,j,1:2) = indvd_delR_delNode_NonPulley(spr_nmbr,init_node)
             
          elseif (spr_node(spr_nmbr,0)==3) then
             delR_delNode(i,j,1:2) = indvd_delR_delNode_Pulley_3node(spr_nmbr,prev_spr_nmbr,init_node)
             
          elseif (spr_node(spr_nmbr,0)==4) then   
             delR_delNode(i,j,1:2) = indvd_delR_delNode_Pulley_4node(spr_nmbr,cnt_NTN1,cnt_NTN2,init_node)
          endif
          
          if(i==103.or.i==159)write(327,*)delR_delNode(i,j,1:N_dmnsn),i,j,spr_nmbr,typ_spr(spr_nmbr),"delR_delN"
          
       enddo
       
    enddo
    
    close(327)
    
  contains
    
    function indvd_delR_delNode_NonPulley(spr_nmbr,init_node)
      implicit none
      real*8  :: indvd_delR_delNode_NonPulley(1:N_dmnsn)
      integer :: spr_nmbr
      integer :: finl_node,init_node !node_indexes
      integer :: i
      logical :: lgcl_dn,lgcl_sameDn
      
      indvd_delR_delNode_NonPulley = 0.0d0
      lgcl_dn = .False.
      lgcl_sameDn = .False.
      
      !if (typ_spr(spr_nmbr) .ne. pulley_typ) then
      
      do i = 1,max_node_spr
         
         if ((spr_node(spr_nmbr,i).ne.init_node) .and. (spr_node(spr_nmbr,i).ne.(-1))) then
            finl_node = spr_node(spr_nmbr,i)
            
            call are_they_same_dn(lgcl_sameDn,init_node,finl_node)
            !write(*,*) lgcl_sameDn,"lgcl_sameDn"
            
            if (lgcl_sameDn .eqv. .True.) then
               finl_node = 0
               cycle
            endif
            
            indvd_delR_delNode_NonPulley(1:N_dmnsn) = get_unit_vector(finl_node,init_node)
            
         endif
      enddo
      
    end function indvd_delR_delNode_NonPulley
    
    
    function indvd_delR_delNode_Pulley_3node(spr_nmbr,prev_spr_nmbr,init_node)
      implicit none
      
      real*8  :: indvd_delR_delNode_Pulley_3node(1:N_dmnsn)
      integer :: spr_nmbr,prev_spr_nmbr
      integer :: finl_node,init_node !node_indexes
      logical :: lgcl_trminl
      
      indvd_delR_delNode_Pulley_3node = 0.0d0
      
      !if (init_node.eq.17) write(*,*) init_node,spr_nmbr,"init_node,spr_nmbr"
      
      call is_it_a_trminl_node(init_node,spr_nmbr,lgcl_trminl)
      call get_the_final_node_for_Pulley_case_3node(lgcl_trminl,spr_nmbr,prev_spr_nmbr,finl_node)
      
      indvd_delR_delNode_Pulley_3node(1:N_dmnsn) = get_unit_vector(finl_node,init_node)
      
    end function indvd_delR_delNode_Pulley_3node
    
    function indvd_delR_delNode_Pulley_4node(spr_nmbr,cnt_NTN1,cnt_NTN2,init_node)
      implicit none
      
      real*8                 :: indvd_delR_delNode_Pulley_4node(1:N_dmnsn)
      integer, intent(in)    :: spr_nmbr
      integer, intent(inout) :: cnt_NTN1,cnt_NTN2
      
      integer :: finl_node,init_node !node_indexes
      logical :: lgcl_trminl
      
      indvd_delR_delNode_Pulley_4node = 0.0d0
      
      call is_it_a_trminl_node(init_node,spr_nmbr,lgcl_trminl)
      call get_the_final_node_for_Pulley_case_4node(lgcl_trminl,spr_nmbr,cnt_NTN1,cnt_NTN2,finl_node)
      
      indvd_delR_delNode_Pulley_4node(1:N_dmnsn) = get_unit_vector(finl_node,init_node)
      
    end function indvd_delR_delNode_Pulley_4node
    
    
    subroutine get_the_final_node_for_Pulley_case_3node(lgcl_trminl,spr_nmbr,prev_spr_nmbr,finl_node)
      implicit none
      logical, intent(in)  :: lgcl_trminl
      integer, intent(in)  :: spr_nmbr,prev_spr_nmbr
      integer, intent(out) :: finl_node
      
      if (lgcl_trminl .eqv. .True.) then
         finl_node = spr_node(spr_nmbr,2)
         !if (init_node.eq.17) write(*,*) finl_node,spr_node(spr_nmbr,2),"2"
         
      elseif (lgcl_trminl .eqv. .False.) then
         
         if (spr_nmbr .ne. prev_spr_nmbr) then
            finl_node = spr_node(spr_nmbr,1)
            !if (init_node.eq.17) write(*,*) finl_node,spr_node(spr_nmbr,1),"1"
            
         elseif (spr_nmbr .eq. prev_spr_nmbr) then
            finl_node = spr_node(spr_nmbr,3)
            !if (init_node.eq.17) write(*,*) finl_node,spr_node(spr_nmbr,3),"3"
            
         endif
         
      endif
      
      !write(*,*) init_node,finl_node,"init_node,finl_node"
      
    end subroutine get_the_final_node_for_Pulley_case_3node
    
    subroutine get_the_final_node_for_Pulley_case_4node(lgcl_trminl,spr_nmbr,cnt_NTN1,cnt_NTN2,finl_node)
      implicit none
      logical, intent(in)    :: lgcl_trminl
      integer, intent(in)    :: spr_nmbr
      integer, intent(inout) :: cnt_NTN1,cnt_NTN2
      integer, intent(out)   :: finl_node
      
      if (SystemTyp.ne.1) then
         write(*,*) "may or may not work for SytemTyp different than 1"
         write(*,*) "CHECK the Analytical and Numerical gradient MATCHING"
      endif
      
      if (lgcl_trminl .eqv. .True.) then
         
         if (init_node==spr_node(spr_nmbr,1)) then
            finl_node = spr_node(spr_nmbr,2)
            !write(*,*) init_node,finl_node,spr_node(spr_nmbr,2),"lgcl_trm1"
            
         elseif (init_node==spr_node(spr_nmbr,4)) then
            finl_node = spr_node(spr_nmbr,3)
            !write(*,*) init_node,finl_node,spr_node(spr_nmbr,3),"lgcl_trm2"
         endif
         
      elseif (lgcl_trminl .eqv. .False.) then
         
         if (init_node==spr_node(spr_nmbr,2)) then
            
            if (cnt_NTN1==0) then
               finl_node = spr_node(spr_nmbr,1)
               cnt_NTN1  = 1
               !write(*,*) init_node,finl_node,cnt_NTN1,"lgcl_trm3_1"
               
            elseif (cnt_NTN1==1) then
               finl_node = spr_node(spr_nmbr,3)
               cnt_NTN1  = 0
               !write(*,*) init_node,finl_node,cnt_NTN1,"lgcl_trm3_1"
            endif
            
         elseif (init_node==spr_node(spr_nmbr,3)) then
            
            if (cnt_NTN2==0) then
               finl_node = spr_node(spr_nmbr,4)
               cnt_NTN2  = 1
               !write(*,*) init_node,finl_node,cnt_NTN2,"lgcl_trm4_1"
               
            elseif (cnt_NTN2==1) then
               finl_node = spr_node(spr_nmbr,2)
               cnt_NTN2  = 0
               !write(*,*) init_node,finl_node,cnt_NTN2,"lgcl_trm4_2"
            endif
            
         endif
         
      endif
      
    end subroutine get_the_final_node_for_Pulley_case_4node
    
    function get_unit_vector(finl_node,init_node)
      implicit none
      real*8  :: get_unit_vector(1:N_dmnsn)
      integer :: finl_node,init_node
      !real*8  :: vectr(1:2)
      
      real*8  :: r2(1:N_dmnsn),r1(1:N_dmnsn)
      real*8  :: r12(1:N_dmnsn)
      real*8  :: r12_abs
      
      real*8  :: r12_uv(1:N_dmnsn)
      real*16 :: prcsn = 1.0d-30
      
      !write(*,*) init_node,finl_node,"init_finl_node_in_get_unit_vector"
      !write(*,*) dum_Nodes(init_node,1:2),"dum_init"
      !write(*,*) dum_Nodes(finl_node,1:2),"dum_finl"
      
      r1 = dum_Nodes(init_node,1:N_dmnsn)
      r2 = dum_Nodes(finl_node,1:N_dmnsn)
      
      !write(*,*) init_node,finl_node,"init_finl"
      !write(*,*) r1,"r1"
      !write(*,*) r2,"r2"
      
      r12     = r1-r2
      r12_abs = sqrt(r12(1)**2 + r12(2)**2)
      
      if (abs(r12_abs) .lt. prcsn) then
         write(*,*) "check vector_absValue_grdcalcE.dat"
         open(unit=125,file='vector_absValue_grdcalcE.dat')
         
         write(125,*) init_node,finl_node,"init and finl node"
         write(125,*) "less than precision value, flnm:gradient_calcEnrgy,sbrtn: get_unit_vector"
         write(125,*) r1(1:2),r2(1:2),"r1,r2"
         write(125,*) r12(1),r12(2),r12_abs,"r12(1:2)"
         write(125,*) " "
         
         close(125)
         stop
      endif
      !write(*,*) r12,r12_abs,"in get_unit_vector"
      
      r12_uv  = r12/r12_abs
      
      get_unit_vector = r12_uv
      
    end function get_unit_vector
    
  end subroutine calc_delR_delNode_array

  
  subroutine check_Ana_Grd_for_spr
    implicit none
    
    integer :: node
    integer :: i
    integer :: n1,n2,n3,n4
    integer :: s1,s2,s3,s4
    
    real*8  :: dEdR(1:4)
    real*8  :: dRdN(1:4,1:2)
    
    real*8  :: length
    real*8  :: x0,xNn1,xNn2,xNn3,xNn4
    real*8  :: y0,yNn1,yNn2,yNn3,yNn4
    
    real*8  :: grdVx,grdVy
    
    open(unit=61,file='SprGrdChk.dat')
    
    node=29
    
    s1=92 ; s2=96 ; s3=190 ; s4=194 !conn_sprs
    n1=27 ; n2=30 ; n3=57  ; n4=60 !neighNodes
    
    grdVx = 0.0d0 ; grdVy = 0.0d0
    
    xNn1 = node_xy(n1,1) ; yNn1 = node_xy(n1,2)
    xNn2 = node_xy(n2,1) ; yNn2 = node_xy(n2,2)
    xNn3 = node_xy(n3,1) ; yNn3 = node_xy(n3,2)
    xNn4 = node_xy(n4,1) ; yNn4 = node_xy(n4,2)
    x0 = node_xy(node,1) ; y0 = node_xy(node,2)
    
    write(61,*) "Nodes"
    write(61,*) x0,y0,"29"
    write(61,*) xNn1,yNn1,"27"
    write(61,*) xNn2,yNn2,"30"
    write(61,*) xNn3,yNn3,"57"
    write(61,*) xNn4,yNn4,"60"
    
    dEdR(1) = k_spr(s1)*(l(s1)-l0(s1))
    dEdR(2) = k_spr(s2)*(l(s2)-l0(s2))
    dEdR(3) = k_spr(s3)*(l(s3)-l0(s3))
    dEdR(4) = k_spr(s4)*(l(s4)-l0(s4))
    
    write(61,*) dEdR(1:4)
    
    length = sqrt((x0-xNn1)**2 + (y0-yNn1)**2)
    dRdN(1,1) = (x0-xNn1)/length ; dRdN(1,2) = (y0-yNn1)/length
    
    length = sqrt((x0-xNn2)**2 + (y0-yNn2)**2)
    dRdN(2,1) = (x0-xNn2)/length ; dRdN(2,2) = (y0-yNn2)/length
    
    length = sqrt((x0-xNn3)**2 + (y0-yNn3)**2)
    dRdN(3,1) = (x0-xNn3)/length ; dRdN(3,2) = (y0-yNn3)/length
    
    length = sqrt((x0-xNn4)**2 + (y0-yNn4)**2)
    dRdN(4,1) = (x0-xNn4)/length ; dRdN(4,2) = (y0-yNn4)/length
    
    do i = 1,4
       grdVx = grdVx + dEdR(i)*dRdN(i,1)
       grdVy = grdVy + dEdR(i)*dRdN(i,2)
    enddo
    
    write(61,*) grdVx,grdVy,"grdVx,grdVy"

    close(61)
    
  end subroutine check_Ana_Grd_for_spr
  
  
  
  
  subroutine check_Ana_Grd_for_spr_node1
    implicit none
    
    integer :: node
    integer :: i
    integer :: n1,n2,n3,n4
    integer :: s1,s2,s3,s4
    
    real*8  :: dEdR(1:1)
    real*8  :: dRdN(1:1,1:2)
    
    real*8  :: length
    real*8  :: x0,xNn1,xNn2,xNn3,xNn4
    real*8  :: y0,yNn1,yNn2,yNn3,yNn4
    
    real*8  :: grdVx,grdVy
    
    open(unit=61,file='SprGrdChk.dat')
    
    node=1
    
    s1=1    !conn_sprs
    n1=49   !neighNodes
    
    grdVx = 0.0d0 ; grdVy = 0.0d0
    
    xNn1 = node_xy(n1,1)   ; yNn1 = node_xy(n1,2)
    x0   = node_xy(node,1) ; y0   = node_xy(node,2)
    
    write(61,*) "Nodes"
    write(61,*) x0,y0,"1"
    write(61,*) xNn1,yNn1,"49"
    
    dEdR(1) = k_spr(s1)*(l(s1)-l0(s1))
    
    write(61,*) dEdR(1:1)
    
    length = sqrt((x0-xNn1)**2 + (y0-yNn1)**2)
    dRdN(1,1) = (x0-xNn1)/length ; dRdN(1,2) = (y0-yNn1)/length
     
    do i = 1,1
       grdVx = grdVx + dEdR(i)*dRdN(i,1)
       grdVy = grdVy + dEdR(i)*dRdN(i,2)
    enddo
    
    write(61,*) grdVx,grdVy,"grdVx,grdVy"

    close(61)
    
  end subroutine check_Ana_Grd_for_spr_node1
  
  
  
  
  
  
  subroutine check_Ana_Grd_for_spr_47
    implicit none
    
    integer :: node
    integer :: i
    integer :: n1,n2,n3
    integer :: s1,s2,s3
    
    real*8  :: dEdR(1:3)
    real*8  :: dRdN(1:3,1:2)
    
    real*8  :: length
    real*8  :: x0,xNn1,xNn2,xNn3
    real*8  :: y0,yNn1,yNn2,yNn3
    
    real*8  :: grdVx,grdVy
    
    open(unit=61,file='SprGrdChk1.dat')
    
    node=47
    
    s1=64 ; s2=66 ; s3=67 !conn_sprs
    n1=50 ; n2=48 ; n3=50 !neighNodes
    
    grdVx = 0.0d0 ; grdVy = 0.0d0
    
    xNn1 = node_xy(n1,1)   ; yNn1 = node_xy(n1,2)
    xNn2 = node_xy(n2,1)   ; yNn2 = node_xy(n2,2)
    xNn3 = node_xy(n3,1)   ; yNn3 = node_xy(n3,2)
    x0   = node_xy(node,1) ; y0 = node_xy(node,2)
    
    write(61,*) "Nodes"
    write(61,*) x0,y0,"47"
    write(61,*) xNn1,yNn1,"50"
    write(61,*) xNn2,yNn2,"48"
    write(61,*) xNn3,yNn3,"50"
    
    dEdR(1) = k_spr(s1)*(l(s1)-l0(s1))
    dEdR(2) = k_spr(s2)*(l(s2)-l0(s2))
    dEdR(3) = k_spr(s3)*(l(s3)-l0(s3))
    
    write(61,*) dEdR(1:3),"dEdR"
    
    length = sqrt((x0-xNn1)**2 + (y0-yNn1)**2)
    dRdN(1,1) = (x0-xNn1)/length ; dRdN(1,2) = (y0-yNn1)/length
    
    length = sqrt((x0-xNn2)**2 + (y0-yNn2)**2)
    dRdN(2,1) = (x0-xNn2)/length ; dRdN(2,2) = (y0-yNn2)/length
    
    length = sqrt((x0-xNn3)**2 + (y0-yNn3)**2)
    dRdN(3,1) = (x0-xNn3)/length ; dRdN(3,2) = (y0-yNn3)/length
    
    write(61,*) dRdN(1,1:2),"dRdN"
    write(61,*) dRdN(2,1:2),"dRdN"
    write(61,*) dRdN(3,1:2),"dRdN"
    
    do i = 1,3
       grdVx = grdVx + dEdR(i)*dRdN(i,1)
       grdVy = grdVy + dEdR(i)*dRdN(i,2)
    enddo
    
    write(61,*) grdVx,grdVy,"grdVx,grdVy"

    close(61)
    
  end subroutine check_Ana_Grd_for_spr_47

  subroutine check_Ana_Grd_for_spr_48
    implicit none
    
    integer :: node
    integer :: i
    integer :: n1,n2,n3
    integer :: s1,s2,s3
    
    real*8  :: dEdR(1:3)
    real*8  :: dRdN(1:3,1:2)
    
    real*8  :: length
    real*8  :: x0,xNn1,xNn2,xNn3
    real*8  :: y0,yNn1,yNn2,yNn3
    
    real*8  :: grdVx,grdVy
    
    open(unit=61,file='SprGrdChk2.dat')
    
    node=48
    
    s1=65 ; s2=66 ; s3=68 !conn_sprs
    n1=46 ; n2=47 ; n3=24 !neighNodes
    
    grdVx = 0.0d0 ; grdVy = 0.0d0
    
    xNn1 = node_xy(n1,1) ; yNn1 = node_xy(n1,2)
    xNn2 = node_xy(n2,1) ; yNn2 = node_xy(n2,2)
    xNn3 = node_xy(n3,1) ; yNn3 = node_xy(n3,2)
    x0 = node_xy(node,1) ; y0 = node_xy(node,2)
    
    write(61,*) "Nodes"
    write(61,*) x0,y0,"48"
    write(61,*) xNn1,yNn1,"46"
    write(61,*) xNn2,yNn2,"47"
    write(61,*) xNn3,yNn3,"24"
    
    dEdR(1) = k_spr(s1)*(l(s1)-l0(s1))
    dEdR(2) = k_spr(s2)*(l(s2)-l0(s2))
    dEdR(3) = k_spr(s3)*(l(s3)-l0(s3))
    
    write(61,*) dEdR(1:3)
    
    length = sqrt((x0-xNn1)**2 + (y0-yNn1)**2)
    dRdN(1,1) = (x0-xNn1)/length ; dRdN(1,2) = (y0-yNn1)/length
    
    length = sqrt((x0-xNn2)**2 + (y0-yNn2)**2)
    dRdN(2,1) = (x0-xNn2)/length ; dRdN(2,2) = (y0-yNn2)/length
    
    length = sqrt((x0-xNn3)**2 + (y0-yNn3)**2)
    dRdN(3,1) = (x0-xNn3)/length ; dRdN(3,2) = (y0-yNn3)/length

    write(61,*) dRdN(1,1:2),"dRdN"
    write(61,*) dRdN(2,1:2),"dRdN"
    write(61,*) dRdN(3,1:2),"dRdN"
    
    do i = 1,3
       grdVx = grdVx + dEdR(i)*dRdN(i,1)
       grdVy = grdVy + dEdR(i)*dRdN(i,2)
    enddo
    
    write(61,*) grdVx,grdVy,"grdVx,grdVy"

    close(61)
    
  end subroutine check_Ana_Grd_for_spr_48
  
  subroutine check_Ana_Grd_for_spr_49
    implicit none
    
    integer :: node
    integer :: i
    integer :: n1,n2,n3,n4
    integer :: s1,s2,s3,s4
    
    real*8  :: dEdR(1:4)
    real*8  :: dRdN(1:4,1:2)
    
    real*8  :: length
    real*8  :: x0,xNn1,xNn2,xNn3,xNn4
    real*8  :: y0,yNn1,yNn2,yNn3,yNn4
    
    real*8  :: grdVx,grdVy
    
    open(unit=61,file='SprGrdChk3.dat')
    
    node=49
    
    s1=31 ; s2=31 ; s3=67 ; s4=67 !conn_sprs
    n1=21 ; n2=23 ; n3=50 ; n4=23 !neighNodes
    
    grdVx = 0.0d0 ; grdVy = 0.0d0
    
    xNn1 = node_xy(n1,1) ; yNn1 = node_xy(n1,2)
    xNn2 = node_xy(n2,1) ; yNn2 = node_xy(n2,2)
    xNn3 = node_xy(n3,1) ; yNn3 = node_xy(n3,2)
    xNn4 = node_xy(n4,1) ; yNn4 = node_xy(n4,2)
    x0 = node_xy(node,1) ; y0 = node_xy(node,2)
    
    write(61,*) "Nodes"
    write(61,*) x0,y0,"init node"
    write(61,*) xNn1,yNn1,"n1"
    write(61,*) xNn2,yNn2,"n2"
    write(61,*) xNn3,yNn3,"n3"
    write(61,*) xNn4,yNn4,"n4"
    
    dEdR(1) = k_spr(s1)*(l(s1)-l0(s1))
    dEdR(2) = k_spr(s2)*(l(s2)-l0(s2))
    dEdR(3) = k_spr(s3)*(l(s3)-l0(s3))
    dEdR(4) = k_spr(s4)*(l(s4)-l0(s4))
    
    write(61,*) dEdR(1:4)
    
    length = sqrt((x0-xNn1)**2 + (y0-yNn1)**2)
    dRdN(1,1) = (x0-xNn1)/length ; dRdN(1,2) = (y0-yNn1)/length
    
    length = sqrt((x0-xNn2)**2 + (y0-yNn2)**2)
    dRdN(2,1) = (x0-xNn2)/length ; dRdN(2,2) = (y0-yNn2)/length
    
    length = sqrt((x0-xNn3)**2 + (y0-yNn3)**2)
    dRdN(3,1) = (x0-xNn3)/length ; dRdN(3,2) = (y0-yNn3)/length
    
    length = sqrt((x0-xNn4)**2 + (y0-yNn4)**2)
    dRdN(4,1) = (x0-xNn4)/length ; dRdN(4,2) = (y0-yNn4)/length
    
    do i = 1,4
       grdVx = grdVx + dEdR(i)*dRdN(i,1)
       grdVy = grdVy + dEdR(i)*dRdN(i,2)
    enddo
    
    write(61,*) grdVx,grdVy,"grdVx,grdVy"

    close(61)
    
  end subroutine check_Ana_Grd_for_spr_49
  
  subroutine check_Ana_Grd_for_spr_50
    implicit none
    
    integer :: node
    integer :: i
    integer :: n1,n2,n3,n4
    integer :: s1,s2,s3,s4
    
    real*8  :: dEdR(1:4)
    real*8  :: dRdN(1:4,1:2)
    
    real*8  :: length
    real*8  :: x0,xNn1,xNn2,xNn3,xNn4
    real*8  :: y0,yNn1,yNn2,yNn3,yNn4
    
    real*8  :: grdVx,grdVy
    
    open(unit=61,file='SprGrdChk4.dat')
    
    node=50
    
    s1=64 ; s2=64 ; s3=67 ; s4=67 !conn_sprs
    n1=45 ; n2=47 ; n3=49 ; n4=47 !neighNodes
    
    grdVx = 0.0d0 ; grdVy = 0.0d0
    
    xNn1 = node_xy(n1,1) ; yNn1 = node_xy(n1,2)
    xNn2 = node_xy(n2,1) ; yNn2 = node_xy(n2,2)
    xNn3 = node_xy(n3,1) ; yNn3 = node_xy(n3,2)
    xNn4 = node_xy(n4,1) ; yNn4 = node_xy(n4,2)
    x0 = node_xy(node,1) ; y0 = node_xy(node,2)
    
    write(61,*) "Nodes"
    write(61,*) x0,y0,"init node"
    write(61,*) xNn1,yNn1,"n1"
    write(61,*) xNn2,yNn2,"n2"
    write(61,*) xNn3,yNn3,"n3"
    write(61,*) xNn4,yNn4,"n4"
    
    dEdR(1) = k_spr(s1)*(l(s1)-l0(s1))
    dEdR(2) = k_spr(s2)*(l(s2)-l0(s2))
    dEdR(3) = k_spr(s3)*(l(s3)-l0(s3))
    dEdR(4) = k_spr(s4)*(l(s4)-l0(s4))
    
    write(61,*) dEdR(1:4)
    
    length = sqrt((x0-xNn1)**2 + (y0-yNn1)**2)
    dRdN(1,1) = (x0-xNn1)/length ; dRdN(1,2) = (y0-yNn1)/length
    
    length = sqrt((x0-xNn2)**2 + (y0-yNn2)**2)
    dRdN(2,1) = (x0-xNn2)/length ; dRdN(2,2) = (y0-yNn2)/length
    
    length = sqrt((x0-xNn3)**2 + (y0-yNn3)**2)
    dRdN(3,1) = (x0-xNn3)/length ; dRdN(3,2) = (y0-yNn3)/length
    
    length = sqrt((x0-xNn4)**2 + (y0-yNn4)**2)
    dRdN(4,1) = (x0-xNn4)/length ; dRdN(4,2) = (y0-yNn4)/length
    
    do i = 1,4
       grdVx = grdVx + dEdR(i)*dRdN(i,1)
       grdVy = grdVy + dEdR(i)*dRdN(i,2)
    enddo
    
    write(61,*) grdVx,grdVy,"grdVx,grdVy"

    close(61)
    
  end subroutine check_Ana_Grd_for_spr_50
  
  
  subroutine check_Ana_Grd_for_spr_IN(postn)
    implicit none
    integer, intent(in) :: postn
    
    integer :: node
    integer :: i
    integer :: n1,n2,n3,n4
    integer :: s1,s2
    
    real*8  :: dEdR(1:2)
    real*8  :: dRdN(1:2,1:2)
    
    real*8  :: length
    real*8  :: x0,xNn1,xNn2
    real*8  :: y0,yNn1,yNn2
    
    real*8  :: grdVx,grdVy
    integer :: LftOrRght
    
    !open(unit=61,file='SprGrdChk_IN1.dat',position='append') !IN1=InsertedNode1
    
    ! if (postn==1) then
    !    node=78
       
    !    s1=51 ; s2=52  !conn_sprs
    !    n1=16 ; n2=79  !neighNodes
       
    ! elseif (postn==2) then
       
    !    node=102
       
    !    s1=72 ; s2=73  !conn_sprs
    !    n1=22 ; n2=103 !neighNodes
    ! endif
    
    ! grdVx = 0.0d0 ; grdVy = 0.0d0
    
    ! xNn1 = node_xy(n1,1)   ; yNn1 = node_xy(n1,2)
    ! xNn2 = node_xy(n2,1)   ; yNn2 = node_xy(n2,2)
    ! x0   = node_xy(node,1) ; y0   = node_xy(node,2)
    
    ! write(61,*) "Nodes"
    ! write(61,*) x0,y0,"102"
    ! write(61,*) xNn1,yNn1,"22"
    ! write(61,*) xNn2,yNn2,"103"
    
    ! dEdR(1) = k_spr(s1)*(l(s1)-l0(s1))
    ! dEdR(2) = k_spr(s2)*(l(s2)-l0(s2))
    
    ! write(61,*) dEdR(1:2)
    
    ! length = sqrt((x0-xNn1)**2 + (y0-yNn1)**2)
    ! dRdN(1,1) = (x0-xNn1)/length ; dRdN(1,2) = (y0-yNn1)/length
    
    ! length = sqrt((x0-xNn2)**2 + (y0-yNn2)**2)
    ! dRdN(2,1) = (x0-xNn2)/length ; dRdN(2,2) = (y0-yNn2)/length
    
    ! write(61,*)  dRdN(1,1), dRdN(1,2), dRdN(2,1), dRdN(2,2),"dRdN"
    
    ! do i = 1,2
    !    grdVx = grdVx + dEdR(i)*dRdN(i,1)
    !    grdVy = grdVy + dEdR(i)*dRdN(i,2)
    ! enddo
    
    ! write(61,*) grdVx,grdVy,"grdVx,grdVy"
    
    ! close(61)
    
    
    open(unit=62,file='SprGrdChk_IN2.dat')
    
    do LftOrRght = 1,2
       
       if (postn==1) then
          
          if (LftOrRght==1) then
             node=79
             
             s1=52  ; s2=53  !conn_sprs
             n1=78  ; n2=18  !neighNodes
             
          elseif (LftOrRght==2) then
             node = 123
             
             s1=129 ; s2=130
             n1=122 ; n2=42
             
          endif
          
       elseif (postn==2) then
          
          if (LftOrRght==1) then
             write(62,*) Postn,LftOrRght,"Lft/Rght + Postn"
             
             node=103
             
             s1=73  ; s2=74  !conn_sprs
             n1=102 ; n2=24  !neighNodes
             
          elseif (LftOrRght==2) then
             write(62,*) Postn,LftOrRght,"Lft/Rght + Postn"
             
             node=159
             
             s1=171  ; s2=172  !conn_sprs
             n1=158  ; n2=54  !neighNodes
             
          endif
          
       endif
       
       grdVx = 0.0d0 ; grdVy = 0.0d0
    
       xNn1 = node_xy(n1,1)   ; yNn1 = node_xy(n1,2)
       xNn2 = node_xy(n2,1)   ; yNn2 = node_xy(n2,2)
       x0   = node_xy(node,1) ; y0   = node_xy(node,2)
       
       write(62,*) "Nodes"
       write(62,*) x0,y0,"103"
       write(62,*) xNn1,yNn1,"102"
       write(62,*) xNn2,yNn2,"24"
       
       dEdR(1) = k_spr(s1)*(l(s1)-l0(s1))
       dEdR(2) = k_spr(s2)*(l(s2)-l0(s2))
       
       write(62,*) dEdR(1:2)
       
       length = sqrt((x0-xNn1)**2 + (y0-yNn1)**2)
       dRdN(1,1) = (x0-xNn1)/length ; dRdN(1,2) = (y0-yNn1)/length
       
       length = sqrt((x0-xNn2)**2 + (y0-yNn2)**2)
       dRdN(2,1) = (x0-xNn2)/length ; dRdN(2,2) = (y0-yNn2)/length
       
       do i = 1,2
          grdVx = grdVx + dEdR(i)*dRdN(i,1)
          grdVy = grdVy + dEdR(i)*dRdN(i,2)
       enddo
       
       write(62,*) grdVx,grdVy,"grdVx,grdVy"
       
    enddo
       
    close(62)
    
    
  end subroutine check_Ana_Grd_for_spr_IN
  
  subroutine chk_Grd_For_twoEquivNode(postn)
    implicit none
    integer, intent(in) :: postn
    
    integer :: nodeL,nodeR
    integer :: nL1,nL2,nR1,nR2
    integer :: sL1,sL2,sR1,sR2
    
    real*8  :: dEdR_L(1:2),dEdR_R(1:2)
    real*8  :: dRdN_L(1:2,1:2),dRdN_R(1:2,1:2)
    
    real*8  :: length
    real*8  :: x0L,y0L,x0R,y0R
    real*8  :: xNn1L,yNn1L,xNn1R,yNn1R
    real*8  :: xNn2L,yNn2L,xNn2R,yNn2R
    
    real*8  :: grdVxL,grdVyL
    real*8  :: grdVxR,grdVyR

    integer :: i
    
    open(unit=64,file='GrdChk_TwoEquiv.dat')
    
    write(64,*) "START"
    
    if (postn==1) then
       nodeL = 79 ; nodeR = 123
       
       sL1 = 52 ; sL2 = 53  ; sR1 = 129 ; sR2 = 130
       nL1 = 78 ; nL2 = 18  ; nR1 = 122 ; nR2 = 42
       
    elseif (postn==2) then
       nodeL = 103 ; nodeR = 159
       
       sL1 = 73  ; sL2 = 74  ; sR1 = 171 ; sR2 = 172
       nL1 = 102 ; nL2 = 24  ; nR1 = 158 ; nR2 = 54
    endif
       
    grdVxL = 0.0d0 ; grdVyL = 0.0d0
    grdVxR = 0.0d0 ; grdVyR = 0.0d0
    
    xNn1L = node_xy(nL1,1)   ; yNn1L = node_xy(nL1,2)
    xNn2L = node_xy(nL2,1)   ; yNn2L = node_xy(nL2,2)
    x0L   = node_xy(nodeL,1) ; y0L   = node_xy(nodeL,2)
    
    xNn1R = node_xy(nR1,1)   ; yNn1R = node_xy(nR1,2)
    xNn2R = node_xy(nR2,1)   ; yNn2R = node_xy(nR2,2)
    x0R   = node_xy(nodeR,1) ; y0R   = node_xy(nodeR,2)
    
    write(64,*) "Nodes"
    write(64,*) x0L,x0R,"x0 of 102 and 158"
    write(64,*) y0L,y0R,"y0 of 102 and 158"
    
    write(64,*) xNn1L,xNn1R,"xNn1 22/52" 
    write(64,*) xNn2L,xNn2R,"xNn2 103/159"
    
    write(64,*) yNn1L,yNn1R,"yNn1 22/52"
    write(64,*) yNn2L,yNn2R,"yNn2 103/159"

    write(64,*) " "
    
    dEdR_L(1) = k_spr(sL1)*(l(sL1)-l0(sL1)) ; dEdR_R(1) = k_spr(sR1)*(l(sR1)-l0(sR1))
    dEdR_L(2) = k_spr(sL2)*(l(sL2)-l0(sL2)) ; dEdR_R(2) = k_spr(sR2)*(l(sR2)-l0(sR2))

    write(64,*) "spr k"
    write(64,*) k_spr(sL1),k_spr(sR1)
    write(64,*) k_spr(sL2),k_spr(sR2)
    write(64,*) " "
    
    write(64,*) "lengths"
    write(64,*) l(sL1),l(sR1)
    write(64,*) l(sL2),l(sR2)
    write(64,*) " "

    write(64,*) "rest lengths"
    write(64,*) l0(sL1),l0(sR1)
    write(64,*) l0(sL2),l0(sR2)
    write(64,*) " "
    
    write(64,*) "l-l0"
    write(64,*) (l(sL1)-l0(sL1)),(l(sR1)-l0(sR1))
    write(64,*) (l(sL2)-l0(sL2)),(l(sR2)-l0(sR2))
    write(64,*) " "
    
    write(64,*) dEdR_L(1:2),"dEdR_L"
    write(64,*) dEdR_R(1:2),"dEdR_R"


    length      = sqrt((x0L-xNn1L)**2 + (y0L-yNn1L)**2)
    dRdN_L(1,1) = (x0L-xNn1L)/length ; dRdN_L(1,2) = (y0L-yNn1L)/length 
    
    length      = sqrt((x0L-xNn2L)**2 + (y0L-yNn2L)**2)
    dRdN_L(2,1) = (x0L-xNn2L)/length ; dRdN_L(2,2) = (y0L-yNn2L)/length
    
    do i = 1,2
       grdVxL = grdVxL + dEdR_L(i)*dRdN_L(i,1)
       grdVyL = grdVyL + dEdR_L(i)*dRdN_L(i,2)
    enddo

    
    length      = sqrt((x0R-xNn1R)**2 + (y0R-yNn1R)**2)
    dRdN_R(1,1) = (x0R-xNn1R)/length ; dRdN_R(1,2) = (y0R-yNn1R)/length 
    
    length      = sqrt((x0R-xNn2R)**2 + (y0R-yNn2R)**2)
    dRdN_R(2,1) = (x0R-xNn2R)/length ; dRdN_R(2,2) = (y0R-yNn2R)/length
    
    do i = 1,2
       grdVxR = grdVxR + dEdR_R(i)*dRdN_R(i,1)
       grdVyR = grdVyR + dEdR_R(i)*dRdN_R(i,2)
    enddo

    write(64,*) "dRdN"
    write(64,*) dRdN_L(1,1),dRdN_L(1,2),dRdN_L(2,1),dRdN_L(2,2),"dRdN_L"
    write(64,*) dRdN_R(1,1),dRdN_R(1,2),dRdN_R(2,1),dRdN_R(2,2),"dRdN_R"
    write(64,*) " "
    
    write(64,*) "grdVx~grdVy"
    write(64,*) grdVxL,grdVyL,"grdV L"
    write(64,*) grdVxR,grdVyR,"grdV_R"
    
    write(64,*) "END"
    close(64)
    
  end subroutine chk_Grd_For_twoEquivNode
  
  
end module gradient_contrb_from_spr


module gradient_contrb_from_area  !M3
  use system_parameters !cell_info + non_changing_parameters
  use transfrm_info
  use triangle_routines
  use vitelline_fluid_calctn_mod
  
  implicit none
  
  real*8, allocatable :: grd_area(:,:),grd_areaS(:,:),grd_areaS2(:,:),grd_areaSNI(:,:)
  real*8, allocatable :: delE_delA(:) ,delE_delAS(:), delE_delAS2(:) ,delE_delASNI(:)
  real*8, allocatable :: delA_delNode(:,:,:),delA_delNodeS(:,:,:), delA_delNodeS2(:,:,:),delA_delNodeSNI(:,:,:) 
  
contains
  
  subroutine allocate_and_initialize_gradient_area_variables    
    implicit none
    
    allocate(grd_area(1:N_node,1:N_dmnsn),grd_areaS(1:N_node,1:N_dmnsn))
    allocate(delE_delA(1:N_cell),delE_delAS(1:N_cell))
    allocate(delA_delNode(1:N_node,1:max_area_node,1:N_dmnsn),delA_delNodeS(1:N_node,1:max_area_node,1:N_dmnsn))
    
    grd_area     = -1.d30 ; grd_areaS     = -1.d30
    delE_delA    = -1.d30 ; delE_delAS    = -1.d30 !(ridiculous value)
    delA_delNode = -1.d30 ; delA_delNodeS = -1.d30
    
  end subroutine allocate_and_initialize_gradient_area_variables
  
  subroutine allocate_and_initialize_gradient_area_variables_wo_StrVars
    implicit none
    
    allocate(grd_area(1:N_node,1:N_dmnsn))
    allocate(delE_delA(1:N_cell))
    allocate(delA_delNode(1:N_node,1:max_area_node,1:N_dmnsn))
    
    grd_area     = -1.d30
    delE_delA    = -1.d30
    delA_delNode = -1.d30
    
  end subroutine allocate_and_initialize_gradient_area_variables_wo_StrVars
  
  subroutine allocate_and_initialize_gradient_area_variables_StrVars
    implicit none
    
    allocate(grd_areaS(1:N_node,1:N_dmnsn))
    allocate(delE_delAS(1:N_cell))
    allocate(delA_delNodeS(1:N_node,1:max_area_node,1:N_dmnsn))
    
    grd_areaS     = -1.d30
    delE_delAS    = -1.d30
    delA_delNodeS = -1.d30
    
  end subroutine allocate_and_initialize_gradient_area_variables_StrVars
  
  
  subroutine deallocate_gradient_area_variables
    implicit none
    
    deallocate(grd_area,grd_areaS)
    deallocate(delE_delA,delE_delAS)
    deallocate(delA_delNode,delA_delNodeS)
    
  end subroutine deallocate_gradient_area_variables
  
  subroutine deallocate_gradient_area_variables_wo_StrVars
    implicit none
    
    deallocate(grd_area)
    deallocate(delE_delA)
    deallocate(delA_delNode)
    
  end subroutine deallocate_gradient_area_variables_wo_StrVars
  
  subroutine deallocate_gradient_area_variables_StrVars
    implicit none
    
    deallocate(grd_areaS)
    deallocate(delE_delAS)
    deallocate(delA_delNodeS)
    
  end subroutine deallocate_gradient_area_variables_StrVars
  
  
  subroutine store_gradient_area_variables
    implicit none
    
    grd_areaS     = grd_area
    delE_delAS    = delE_delA
    delA_delNodeS = delA_delNode
    
  end subroutine store_gradient_area_variables
  
  subroutine restore_gradient_area_variables
    implicit none
    
    grd_area     = grd_areaS
    delE_delA    = delE_delAS
    delA_delNode = delA_delNodeS
    
  end subroutine restore_gradient_area_variables
  
  
  subroutine store_gradient_area_variables_NI
    implicit none
    
    grd_areaSNI     = grd_area
    delE_delASNI    = delE_delA
    delA_delNodeSNI = delA_delNode
    
  end subroutine store_gradient_area_variables_NI
  
  subroutine store_gradient_area_variables_S2
    implicit none
    
    grd_areaS2     = grd_area
    delE_delAS2    = delE_delA
    delA_delNodeS2 = delA_delNode
    
  end subroutine store_gradient_area_variables_S2
  
  
  subroutine gradient_components_frm_area(dum_Nodes,dumA0,dum_grdArea)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8  :: dumA0(1:N_cell)
    real*8  :: dum_grdArea(1:N_node,1:N_dmnsn)
    integer :: i
    integer :: area_nmbr
    integer :: totl_cell_cnnctd 
    
    dum_grdArea  = 0.0d0
    delE_delA    = 0.0d0
    delA_delNode = 0.0d0
    
    call calc_delE_delA_array(dum_Nodes,dumA0)
    call calc_delA_delNode_array(dum_Nodes)
    
    do i = 1,N_node
       if (node_typ(i) .eq. 0) then
          continue
       elseif (node_typ(i) .ne. 0) then
          totl_cell_cnnctd = node_area(i,0)
          !write(*,*) totl_cell_cnnctd,"tcc"
          
          if (node_cnnctd(i) .eq. 0) then
             
             if (node_typ(i) .eq. 1) then !free node
                call grdArea_lp(i)
             elseif (node_typ(i) .eq. 2) then !x free, y fixed
                call grdArea_lp(i)
             endif
             
          elseif (node_cnnctd(i) .ne. 0) then
             
             if (count_this_dn(i) .eq. 1) then
                
                if (node_typ(i) .eq. 1) then !free node
                   call grdArea_lp(i)
                elseif (node_typ(i) .eq. 2) then !x free, y fixed
                   call grdArea_lp(i)
                endif
                
             elseif (count_this_dn(i) .ne. 1) then
                !will have to replace with other double node value
                
                continue
             endif
             
          endif
          
       endif
       
       !write(*,*) dum_grdArea(i,1:2),"dum_grdArea and i = ",i
       
    enddo
    
    
  contains
    
    subroutine grdArea_lp(node_nm)
      implicit none
      integer :: j
      integer :: node_nm
      
      do j = 1,totl_cell_cnnctd
         area_nmbr = node_area(node_nm,j)
         
         if (node_typ(node_nm) .eq. 1) then  !x,y free
            dum_grdArea(node_nm,1:2) = dum_grdArea(node_nm,1:2) &
                 + delE_delA(area_nmbr) * delA_delNode(node_nm,j,1:2)
            
         elseif (node_typ(node_nm) .eq. 2) then !x free, y fixed        
            dum_grdArea(node_nm,1) = dum_grdArea(node_nm,1) &
                 + delE_delA(area_nmbr) * delA_delNode(node_nm,j,1)
         endif
      enddo
      
    end subroutine grdArea_lp
    
    
  end subroutine gradient_components_frm_area
  
  subroutine calc_delE_delA_array(dum_Nodes,dumA0)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8  :: dumA0(1:N_cell)
    integer :: i

    open(unit=289,file="delE_delA.dat")
    
    do i = 1,N_cell
       delE_delA(i) = indvd_delE_delA(i)
       write(289,fmt=*) delE_delA(i),k_area(i),i
    enddo
    
    close(289)
    
  contains
    
    real*8 function indvd_delE_delA(area_nmbr)
      implicit none
      integer :: area_nmbr
      real*8  :: ka,area_ch
      
      ka              = k_area(area_nmbr)
      area_ch         = area_chng(area_nmbr)
      indvd_delE_delA = ka*area_ch
      
      !HERE A(area_nmbr) and A0(area_nmbr) needs to be accesibile
      
    end function indvd_delE_delA
    
    
    real*8 function area_chng(area_nmbr)
      implicit none
      integer :: area_nmbr
      
      !write(289,*) A_currnt(area_nmbr),dumA0(area_nmbr),"A_curr_A0"
      
      if (VF_regionModelled==0) then
         area_chng  = A_currnt(area_nmbr) - dumA0(area_nmbr)   
      elseif (VF_regionModelled==1) then
         if (area_nmbr.ne.N_cell) area_chng = A_currnt(area_nmbr) - dumA0(area_nmbr)
         if (area_nmbr == N_cell) area_chng = vitelline_fluid_region(area_nmbr,dum_Nodes) - dumA0(area_nmbr)
      endif
      
    end function area_chng
    
    real*8 function A_currnt(area_nmbr)
      implicit none
      integer :: area_nmbr
      integer :: Node1,Node2,Node3
      real*8  :: N1(1:3),N2(1:3),N3(1:3)
      integer :: i
      real*8  :: areachk
      
      A_currnt = 0.0d0
      
      i = 1
      
      open(unit=24,file="check_inside_areaE.dat")
      
      do
         if ((i+2) .gt. max_node_area) then
            exit
         elseif (area_node(area_nmbr,(i+2)) .eq. (-1)) then
            exit
         endif
         
         Node1 = area_node(area_nmbr,1)
         Node2 = area_node(area_nmbr,i+1)
         Node3 = area_node(area_nmbr,i+2)
         
         N1(1:2) = dum_Nodes(Node1,1:2)
         N2(1:2) = dum_Nodes(Node2,1:2)
         N3(1:2) = dum_Nodes(Node3,1:2)
         
         N1(3)=0.0d0 ; N2(3)=0.0d0 ; N3(3)=0.0d0
         
         !write(unit=24,fmt='(12(f8.6x),A)') dum_Nodes(1:6,1:2),"inside area_E"
         areachk = calc_triangle(N1,N2,N3)
         !write(*,*) areachk,"Area_chk"
         A_currnt = A_currnt + calc_triangle(N1,N2,N3)
         !write(*,*) N1,N2,N3,"N1-N3"
         !write(*,*) A_currnt,"A_currnt"
         i = i+1
         
      enddo
      
      close(24)
      
    end function A_currnt
    
    
  end subroutine calc_delE_delA_array
  
  
  subroutine calc_delA_delNode_array(dum_Nodes)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    integer :: i,j
    integer :: totl_area_cnnctd
    integer :: area_nmbr,node_nmbr
    
    open(unit=326,file='Track_dA_dN.dat')
    
    do i = 1,N_node
       node_nmbr        = i
       totl_area_cnnctd = node_area(i,0)
       
       do j = 1,totl_area_cnnctd          
          area_nmbr = node_area(i,j)
          
          if (VF_regionModelled==0) then
             delA_delNode(i,j,1:N_dmnsn) = indvd_delA_delNode(area_nmbr,node_nmbr)
          elseif (VF_regionModelled==1) then
             if (area_nmbr.ne.N_cell) delA_delNode(i,j,1:N_dmnsn) = indvd_delA_delNode(area_nmbr,node_nmbr)
             if (area_nmbr == N_cell) delA_delNode(i,j,1:N_dmnsn) = indvd_delVF_delNode(area_nmbr,node_nmbr)
          endif
          
          if (i==39) write(326,*) delA_delNode(i,j,1:N_dmnsn),i,j,"delA_delN"
       enddo
       
    enddo
    
    close(326)
    
 contains
   
   function indvd_delA_delNode(area_nmbr,node_nmbr)
     implicit none
     real*8  :: indvd_delA_delNode(1:N_dmnsn)
     integer :: area_nmbr,node_nmbr
     integer :: neighbour_nodes(1:2)
     integer :: Neigh1,Neigh2     
     real*8  :: first_xy(1:3),second_xy(1:3),third_xy(1:3)
     real*8  :: triangle_area
     
     indvd_delA_delNode = 0.0d0
     call get_neighbour_nodes(area_nmbr,node_nmbr,neighbour_nodes)
     
     !if(node_nmbr==3) write(*,*) neighbour_nodes(1:2),"neigh_nodes"
     
     Neigh1 = neighbour_nodes(1)
     Neigh2 = neighbour_nodes(2)!;write(*,*)Neigh1,Neigh2,"Neigh1,Neigh2 for (Nn,An)=(",node_nmbr,area_nmbr,")"
     
     first_xy(1:2)  = dum_Nodes(node_nmbr,1:2)
     second_xy(1:2) = dum_Nodes(Neigh1,1:2)
     third_xy(1:2)  = dum_Nodes(Neigh2,1:2)
     
     if (node_nmbr==39) write(326,fmt=*) first_xy(1:2),second_xy(1:2),third_xy(1:2)
     
     first_xy(3)=0.0d0 ; second_xy(3)=0.0d0 ; third_xy(3)=0.0d0
    
     triangle_area = calc_triangle(first_xy,second_xy,third_xy)
     
     call check_triangle_area_value(triangle_area)
     
     indvd_delA_delNode(1) = 0.5 * (dum_Nodes(Neigh1,2) - dum_Nodes(Neigh2,2)) ! (y_Neigh1 - y_Neigh2)
     indvd_delA_delNode(2) = 0.5 * (dum_Nodes(Neigh2,1) - dum_Nodes(Neigh1,1)) ! (x_Neigh2 - x_Neigh1)
     
     ! if (modelID==2) then
     !    if (node_nmbr==3.and.area_nmbr==1) then
     !       indvd_delA_delNode(1) = 0.5 * (dum_Nodes(1,2) - dum_Nodes(4,2))
     !       indvd_delA_delNode(2) = 0.5 * (dum_Nodes(4,1) - dum_Nodes(1,1))
     !    elseif (node_nmbr==3.and.area_nmbr==2) then
     !       indvd_delA_delNode(1) = 0.5 * (dum_Nodes(4,2) - dum_Nodes(5,2))
     !       indvd_delA_delNode(2) = 0.5 * (dum_Nodes(5,1) - dum_Nodes(4,1))
     !    endif
     ! endif
     
     !write(*,*) indvd_delA_delNode(1:2),"indvd_delA_delNode"
     
   end function indvd_delA_delNode
   
   function indvd_delVF_delNode(area_nmbr,node_nmbr)
     implicit none
     real*8  :: indvd_delVF_delNode(1:N_dmnsn)
     integer :: area_nmbr,node_nmbr
     integer :: fn,ln,endPstn
     integer :: iV,chkNode,nodePstn,nc,nb,na
     real*8  :: xc,yc,xb,yb,xa,ya
     
     indvd_delVF_delNode = 0.0d0
     if (area_nmbr.ne.N_cell) stop 'area_nmbr=/N_cell'
     
     endPstn = area_node(N_cell,0) 
     fn = area_node(N_cell,1) ; ln = area_node(N_cell,endPstn) !; write(*,*) endPstn,fn,ln,"endPstn,fn,ln"
     
     do iV = 1,endPstn
        chkNode = area_node(N_cell,iV)
        if (chkNode == node_nmbr) then
           nodePstn = iV
        endif
     enddo
     
     if (nodePstn==1) then
        
        nc=area_node(N_cell,nodePstn) ; na=area_node(N_cell,nodePstn+1) !; write(*,*) nc,na,"nc-na-i" 
        xc= dum_Nodes(nc,1)           ; yc=dum_Nodes(nc,2)              !; write(*,*) xc,yc,"xc-yc"
        xa= dum_Nodes(na,1)           ; ya=dum_Nodes(na,2)              !; write(*,*) xa,ya,"xa-ya"
        
        indvd_delVF_delNode(1) = 0.5000d0 * (yc+ya)
        indvd_delVF_delNode(2) = 0.5000d0 * (xc-xa)
          
     elseif (nodePstn==endPstn) then
        
        nc=area_node(N_cell,nodePstn) ; nb=area_node(N_cell,nodePstn-1) !; write(*,*) nc,nb,"nc-nb-i" 
        xc=dum_Nodes(nc,1)            ; yc=dum_Nodes(nc,2)              !; write(*,*) xc,yc,"xc-yc"
        xb=dum_Nodes(nb,1)            ; yb=dum_Nodes(nb,2)              !; write(*,*) xb,yb,"xb-yb"
        
        indvd_delVF_delNode(1) = 0.5000d0 * (-yc-yb)
        indvd_delVF_delNode(2) = 0.5000d0 * (-xc+xb)
        
     else
        
        nc=area_node(N_cell,nodePstn)   ; nb=area_node(N_cell,nodePstn-1)
        na=area_node(N_cell,nodePstn+1) !; write(*,*) nc,nb,na,"nc-nb-na-i" 
        
        xc=dum_Nodes(nc,1)     ; yc=dum_Nodes(nc,2) !; write(*,*) xc,yc,"xc-yc"
        xb=dum_Nodes(nb,1)     ; yb=dum_Nodes(nb,2) !; write(*,*) xb,yb,"xb-yb"
        xa=dum_Nodes(na,1)     ; ya=dum_Nodes(na,2) !; write(*,*) xa,ya,"xa-ya"
        
        indvd_delVF_delNode(1) = 0.5000d0 * (ya-yb)
        indvd_delVF_delNode(2) = 0.5000d0 * (xb-xa)
        
     endif
     
   end function indvd_delVF_delNode
   
 end subroutine calc_delA_delNode_array
 
 
 
 subroutine get_neighbour_nodes(area_nmbr,node_nmbr,neighbour_nodes)
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
         write(*,*) "nodes cant be negative, flnm: gradient_calc, sbrtn: get_neighbour_nodes"
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
           
         !  if ((node_pstn_in_array+1) .le. max_node_areaS) then
              
          !    if (area_node(area_nmbr,(node_pstn_in_array+1)).ne.(-1)) then
           !      neighbour_nodes(1) = area_node(area_nmbr,(node_pstn_in_array+1))
            !  else
             !    neighbour_nodes(1) = area_node(area_nmbr,1) 
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
           !     neighbour_nodes(2) = area_node(area_nmbr,(max_node_area-1))
           !  elseif (area_node(area_nmbr,0)==3) then !for triangles
           !     neighbour_nodes(2) = area_node(area_nmbr,(max_node_area-2))
           !  endif
           !endif
           
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
              
            !elseif (area_nmbr.gt.(ncl+ncr+2)) then
              
            !  if(area_node(area_nmbr,max_node_areaS) .ne. (-1)) then
             !    neighbour_nodes(2) = area_node(area_nmbr,max_node_areaS)
              !else
               !  neighbour_nodes(2) = area_node(area_nmbr,(max_node_areaS-1))
              !endif
              
           !endif
           
        endif
        
     endif
     
   end subroutine get_second_neighbour
   
 end subroutine get_neighbour_nodes
 
end module gradient_contrb_from_area


module gradient_contrb_frm_gravitation
  use system_parameters
  implicit none
  real*8, allocatable :: grd_grvtnl(:,:),grd_grvtnlS(:,:),grd_grvtnlS2(:,:),grd_grvtnlSNI(:,:)
  
contains
  
  subroutine allocate_and_initialize_gradient_grvtnl  
    implicit none
    
    allocate(grd_grvtnl(1:N_node,1:N_dmnsn),grd_grvtnlS(1:N_node,1:N_dmnsn))
    
    grd_grvtnl = -1.d30 ; grd_grvtnlS = -1.d30
    
  end subroutine allocate_and_initialize_gradient_grvtnl
  
  subroutine allocate_and_initialize_gradient_grvtnl_wo_StrVars  
    implicit none
    
    allocate(grd_grvtnl(1:N_node,1:N_dmnsn))
    
    grd_grvtnl = -1.d30
    
  end subroutine allocate_and_initialize_gradient_grvtnl_wo_StrVars
  
  subroutine allocate_and_initialize_gradient_grvtnl_StrVars  
    implicit none
    
    allocate(grd_grvtnlS(1:N_node,1:N_dmnsn))
    
    grd_grvtnlS = -1.d30
    
  end subroutine allocate_and_initialize_gradient_grvtnl_StrVars
  
  
  subroutine deallocate_gradient_grvtnl
    implicit none
    
    deallocate(grd_grvtnl,grd_grvtnlS)
    
  end subroutine deallocate_gradient_grvtnl
  
  subroutine deallocate_gradient_grvtnl_wo_StrVars
    implicit none
    
    deallocate(grd_grvtnl)
    
  end subroutine deallocate_gradient_grvtnl_wo_StrVars
  
  subroutine deallocate_gradient_grvtnl_StrVars
    implicit none
    
    deallocate(grd_grvtnlS)
    
  end subroutine deallocate_gradient_grvtnl_StrVars
  
  subroutine store_gradient_grvtnl_variables
    implicit none
    
    grd_grvtnlS = grd_grvtnl
    
  end subroutine store_gradient_grvtnl_variables
  
  subroutine restore_gradient_grvtnl_variables
    implicit none
    
    grd_grvtnl = grd_grvtnlS
    
  end subroutine restore_gradient_grvtnl_variables
  
  subroutine store_gradient_grvtnl_variables_NI
    implicit none
    
    grd_grvtnlSNI = grd_grvtnl
    
  end subroutine store_gradient_grvtnl_variables_NI
  
  subroutine store_gradient_grvtnl_variables_S2
    implicit none
    
    grd_grvtnlS2 = grd_grvtnl
    
  end subroutine store_gradient_grvtnl_variables_S2
  
  subroutine gradient_components_frm_grvtnl(dum_Nodes,dum_grdGrv)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8  :: dum_grdGrv(1:N_node,1:N_dmnsn)
    integer :: i
    
    !write(*,*) "Check_delX_delY in gradient_comps_frm_grvtnl"
    
    dum_grdGrv(1:N_node,1:2) = 0.0d0
    
    do i = 1,N_node
       
       if (node_typ(i) .eq. 0) then
          continue
          
       elseif (node_typ(i) .ne. 0) then
          
          if (node_cnnctd(i) .eq. 0) then       
             call grdGrvtnl_lp(i)
             
          elseif (node_cnnctd(i) .ne. 0) then
             
             if (count_this_dn(i) .eq. 1) then  
                call grdGrvtnl_lp(i)
                
             elseif (count_this_dn(i) .ne. 1) then
                !will have to replace with other double node value
                
                continue
             endif
             
          endif
       endif
       
    enddo
    
    
  contains
    
    
    subroutine grdGrvtnl_lp1(node_nm)
      !although no loop here, to maintain similarity, lp in name
      implicit none
      integer :: node_nm
      
      if (node_typ(node_nm) == 1) then !free node
         
         if (select_xy .eq. 1) then
            if (node_nm.eq.(N_node-1) .or. node_nm.eq.(N_node)) then
               dum_grdGrv(node_nm,1) = -CgX*1.0d0
            endif
         elseif (select_xy == 2) then
            dum_grdGrv(node_nm,2) = CgX*1.0d0
            
         elseif (select_xy == 3) then
            if (node_nm.eq.(N_node-1) .or. node_nm.eq.(N_node)) then
               dum_grdGrv(node_nm,1) = -CgX*1.0d0
            endif
            dum_grdGrv(node_nm,2) = CgX*1.0d0
         endif
         
      elseif (node_typ(node_nm) .eq. 2) then ! x free, y fixed
         
         if (select_xy == 1) then
            
            if (node_nm.eq.(N_node-1) .or. node_nm.eq.(N_node)) then
               dum_grdGrv(node_nm,1) = -CgX*1.0d0
            endif
            
         elseif (select_xy == 2) then
            continue
         endif
         
      endif
         
    end subroutine grdGrvtnl_lp1
    
    subroutine grdGrvtnl_lp(node_nm)
      implicit none
      integer :: node_nm
      integer :: j    
      
     !write(*,*) CgXNode(node_nm),CgYNode(node_nm),node_nm,"CgVal"
      
      if (node_typ(node_nm) == 1) then !x,y free [Eg = func(x,y)]
         
         if (select_xy==1) then
            dum_grdGrv(node_nm,1) = CgXNode(node_nm)
            dum_grdGrv(node_nm,2) = 0.00d0
            
         elseif (select_xy==2) then
            dum_grdGrv(node_nm,1) = 0.00d0
            dum_grdGrv(node_nm,2) = CgXNode(node_nm)
            
         elseif (select_xy==3) then
            dum_grdGrv(node_nm,1) = CgXNode(node_nm)
            dum_grdGrv(node_nm,2) = CgYNode(node_nm)
         endif
         
      elseif (node_typ(node_nm) == 2) then !x free, y fixed [Eg = func(x,yf)]       
         
         if (select_xy==1) then ! CgX exists, CgY not exists
            dum_grdGrv(node_nm,1) = CgXNode(node_nm)
            dum_grdGrv(node_nm,2) = 0.0d0
            
         elseif (select_xy==2) then ! CgX not exists, CgY exists
            dum_grdGrv(node_nm,1) = 0.0d0
            dum_grdGrv(node_nm,2) = 0.0d0
            
         elseif (select_xy==3) then ! CgX exists, CgY exists exists
            dum_grdGrv(node_nm,1) = CgXNode(node_nm)
            dum_grdGrv(node_nm,2) = 0.0d0
         endif
         
      endif
      
    end subroutine grdGrvtnl_lp
    
  end subroutine gradient_components_frm_grvtnl
  
  
end module gradient_contrb_frm_gravitation

module gradient_contrb_frm_SRYP
  use SRYP_mod
  implicit none
  real*8, allocatable :: grd_SRyp(:,:),grd_SRypS(:,:)
  
contains
  
  subroutine allocate_and_initialize_gradient_SRyp
    implicit none
    allocate(grd_SRyp(1:N_node,1:N_dmnsn),grd_SRypS(1:N_node,1:N_dmnsn))
    grd_SRyp = -1.d30 ; grd_SRypS = -1.d30
  end subroutine allocate_and_initialize_gradient_SRyp
  
  subroutine allocate_and_initialize_gradient_SRyp_wo_StrVars
    implicit none
    allocate(grd_SRyp(1:N_node,1:N_dmnsn))  
    grd_SRyp = -1.d30
  end subroutine allocate_and_initialize_gradient_SRyp_wo_StrVars
  
  subroutine allocate_and_initialize_gradient_SRyp_StrVars  
    implicit none  
    allocate(grd_SRypS(1:N_node,1:N_dmnsn))
    grd_SRypS = -1.d30
  end subroutine allocate_and_initialize_gradient_SRyp_StrVars
  
  subroutine deallocate_gradient_SRyp
    implicit none  
    deallocate(grd_SRyp,grd_SRypS)
  end subroutine deallocate_gradient_SRyp
  
  subroutine deallocate_gradient_SRyp_wo_StrVars
    implicit none
    deallocate(grd_SRyp)
  end subroutine deallocate_gradient_SRyp_wo_StrVars
  
  subroutine deallocate_gradient_SRyp_StrVars
    implicit none
    deallocate(grd_SRypS)
  end subroutine deallocate_gradient_SRyp_StrVars
  
  !subroutine store_gradient_SRyp_variables
  !  implicit none
  !  grd_SRypS = grd_SRyp
  !end subroutine store_gradient_SRyp_variables
  
  subroutine store_gradient_SRyp_variables
    implicit none
    integer :: cntAllctnM,cntAllctnS,cntAllctnSum
    
    cntAllctnM = 0 ; cntAllctnS = 0 ; cntAllctnSum = 0
    
    if (.not.allocated(grd_SRyp)) then
       continue
    else
       cntAllctnM = cntAllctnM+1
    endif
    if (.not.allocated(grd_SRypS)) then
       continue
    else
       cntAllctnS = cntAllctnS+1
    endif
    
    cntAllctnSum = cntAllctnM+cntAllctnS ; write(*,*) cntAllctnSum,cntAllctnM,cntAllctnS,"cntAllctn"
    
    if (cntAllctnSum==2) then
       grd_SRypS = grd_SRyp
    else
       write(*,*) "Check Val Above line"
    endif
    
  end subroutine store_gradient_SRyp_variables
  
  
  subroutine restore_gradient_SRyp_variables
    implicit none
    grd_SRyp = grd_SRypS
  end subroutine restore_gradient_SRyp_variables
  
  subroutine gradient_components_frm_SRyp(dum_Nodes,dum_grdSRyp)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8  :: dum_grdSRyp(1:N_node,1:N_dmnsn)
    integer :: apclNode,i
    real*8  :: portn1,portn2,yNodeV
    
    dum_grdSRyp(1:N_node,1:2) = 0.0d0
    
    do i = 1,numApclNodes
       apclNode = apclNodes(i)
       yNodeV   = dum_Nodes(apclNode,2)
       
       if (activtnFctrApcl(i) == 0) then
          continue
       elseif (activtnFctrApcl(i) == 1) then
          portn1=(AlphaApclYP(i))*(exp(AlphaApclYP(i)*yNodeV))*((-yNodeV+epsApclYP(i))**(-powrApclYP(i)))
          portn2=(-powrApclYP(i))*(exp(AlphaApclYP(i)*yNodeV))*((-yNodeV+epsApclYP(i))**(-powrApclYP(i)-1))*(-1)
          !write(*,*) AmpApclYP(i),portn1,portn2,(portn1+portn2),i,"prtn"
          dum_grdSRyp(apclNode,2) = AmpApclYP(i) * (portn1+portn2)
       endif
       
    enddo
    
  end subroutine gradient_components_frm_SRyp
  
  
end module gradient_contrb_frm_SRYP

module gradient_contrb_frm_bending
  use system_parameters
  use transfrm_info
  use triangle_routines
  use bending_energy ! get_A_and_B
  
  implicit none
  
  real*8, allocatable :: grd_bend(:,:),grd_bendS(:,:),grd_bendS2(:,:),grd_bendSNI(:,:)
  
  real*8, allocatable :: delB_delCN(:,:,:) ,delB_delCNS(:,:,:) ,delB_delCNS2(:,:,:), delB_delCNSNI(:,:,:)
  real*8, allocatable :: delB_delNN1(:,:,:),delB_delNN1S(:,:,:),delB_delNN1S2(:,:,:),delB_delNN1SNI(:,:,:)
  real*8, allocatable :: delB_delNN2(:,:,:),delB_delNN2S(:,:,:),delB_delNN2S2(:,:,:),delB_delNN2SNI(:,:,:)
  
contains
  
  subroutine allocate_and_initialize_gradient_bend_variables
    implicit none
    
    allocate(grd_bend(1:N_node,1:N_dmnsn),grd_bendS(1:N_node,1:N_dmnsn))
    
    allocate(delB_delCN(1:N_node,1:max_area_node,1:N_dmnsn))
    allocate(delB_delCNS(1:N_node,1:max_area_node,1:N_dmnsn))
    allocate(delB_delNN1(1:N_node,1:max_area_node,1:N_dmnsn))
    allocate(delB_delNN1S(1:N_node,1:max_area_node,1:N_dmnsn))
    allocate(delB_delNN2(1:N_node,1:max_area_node,1:N_dmnsn))
    allocate(delB_delNN2S(1:N_node,1:max_area_node,1:N_dmnsn))

    grd_bend = -1.d30 ; grd_bendS = -1.d30
    
    delB_delCN  = -1.d30 ; delB_delCNS  = -1.d30
    delB_delNN1 = -1.d30 ; delB_delNN1S = -1.d30
    delB_delNN2 = -1.d30 ; delB_delNN2S = -1.d30
    
  end subroutine allocate_and_initialize_gradient_bend_variables
  
  subroutine allocate_and_initialize_gradient_bend_variables_wo_StrVars
    implicit none
    
    allocate(grd_bend(1:N_node,1:N_dmnsn))
    
    allocate(delB_delCN(1:N_node,1:max_area_node,1:N_dmnsn))
    allocate(delB_delNN1(1:N_node,1:max_area_node,1:N_dmnsn))
    allocate(delB_delNN2(1:N_node,1:max_area_node,1:N_dmnsn))
    
    grd_bend = -1.d30
    
    delB_delCN  = -1.d30
    delB_delNN1 = -1.d30
    delB_delNN2 = -1.d30
    
  end subroutine allocate_and_initialize_gradient_bend_variables_wo_StrVars
  
  subroutine allocate_and_initialize_gradient_bend_variables_StrVars
    implicit none
    
    allocate(grd_bendS(1:N_node,1:N_dmnsn))
    
    allocate(delB_delCNS(1:N_node,1:max_area_node,1:N_dmnsn))
    allocate(delB_delNN1S(1:N_node,1:max_area_node,1:N_dmnsn))
    allocate(delB_delNN2S(1:N_node,1:max_area_node,1:N_dmnsn))
    
    grd_bendS    = -1.d30
    
    delB_delCNS  = -1.d30
    delB_delNN1S = -1.d30
    delB_delNN2S = -1.d30
    
  end subroutine allocate_and_initialize_gradient_bend_variables_StrVars
  
  subroutine deallocate_gradient_bend_variables
    implicit none
    
    deallocate(grd_bend,grd_bendS)
    
    deallocate(delB_delCN,delB_delCNS)
    deallocate(delB_delNN1,delB_delNN1S)
    deallocate(delB_delNN2,delB_delNN2S)
    
  end subroutine deallocate_gradient_bend_variables
  
  subroutine deallocate_gradient_bend_variables_wo_StrVars
    implicit none
    
    deallocate(grd_bend)
    
    deallocate(delB_delCN)
    deallocate(delB_delNN1)
    deallocate(delB_delNN2)
    
  end subroutine deallocate_gradient_bend_variables_wo_StrVars
  
  subroutine deallocate_gradient_bend_variables_StrVars
    implicit none
    
    deallocate(grd_bendS)
    
    deallocate(delB_delCNS)
    deallocate(delB_delNN1S)
    deallocate(delB_delNN2S)
    
  end subroutine deallocate_gradient_bend_variables_StrVars
  
  subroutine store_gradient_bend_variables
    implicit none
    
    grd_bendS    = grd_bend
    delB_delCNS  = delB_delCN
    delB_delNN1S = delB_delNN1
    delB_delNN2S = delB_delNN2
    
  end subroutine store_gradient_bend_variables
  
  subroutine restore_gradient_bend_variables
    implicit none
    
    grd_bend    = grd_bendS
    delB_delCN  = delB_delCNS
    delB_delNN1 = delB_delNN1S
    delB_delNN2 = delB_delNN2S
    
  end subroutine restore_gradient_bend_variables
  
  subroutine store_gradient_bend_variables_NI
    implicit none
    
    grd_bendSNI    = grd_bend
    delB_delCNSNI  = delB_delCN
    delB_delNN1SNI = delB_delNN1
    delB_delNN2SNI = delB_delNN2
    
  end subroutine store_gradient_bend_variables_NI
  
  subroutine store_gradient_bend_variables_S2
    implicit none
    
    grd_bendS2    = grd_bend
    delB_delCNS2  = delB_delCN
    delB_delNN1S2 = delB_delNN1
    delB_delNN2S2 = delB_delNN2
    
  end subroutine store_gradient_bend_variables_S2
  
  subroutine gradient_components_frm_bend(dum_Nodes,dum_grdBend)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8  :: dum_grdBend(1:N_node,1:N_dmnsn)
    integer :: i
    integer :: totl_area_cnnctd
    
    dum_grdBend = 0.0d0
    delB_delCN  = 0.0d0
    delB_delNN1 = 0.0d0
    delB_delNN2 = 0.0d0
    
    call calc_delB_delN_array(dum_Nodes)
    
    open(unit=123,file="grdBend.dat")
    
    do i = 1,N_node
       
       if (node_typ(i) == 0) then
          continue
          
       elseif (node_typ(i).ne.0) then
          totl_area_cnnctd = node_area(i,0)
          
          if (node_cnnctd(i) == 0) then
             
             if (node_typ(i) .eq. 1) then !free node
                call grdBend_lp(i)
             elseif (node_typ(i) .eq. 2) then !x free, y fixed
                call grdBend_lp(i)
             endif
             
          elseif (node_cnnctd(i) .ne. 0) then
             
             if (count_this_dn(i) .eq. 1) then
                
                if (node_typ(i) .eq. 1) then
                   call grdBend_lp(i)
                elseif (node_typ(i) .eq. 2) then
                   call grdBend_lp(i)
                endif
                
             elseif (count_this_dn(i) .ne. 1) then
                continue
             endif
             
          endif
          
       endif
       
       write(123,fmt=*) dum_grdBend(i,1:2),i,"grdBend,node_nm"
       
    enddo
    
    close(123)
    
  contains
    
    subroutine grdBend_lp(node_nm)
      implicit none
      integer :: j,p
      integer :: area_nm
      integer :: node_nm,neigh1_nm,neigh2_nm
      integer :: node_elmnt
      real*8  :: store(1:2)
      
      open(unit=223,file='grdBendCheck.dat',position='append')
      open(unit=224,file='store.dat',position='append')
      open(unit=230,file='node_19.dat')
      
      if (node_nm==19) write(230,fmt=*) totl_area_cnnctd,"tAc"
      
      do j = 1,totl_area_cnnctd
         area_nm   = node_area(node_nm,j)
         neigh1_nm = Nlist(node_nm,j,1)
         neigh2_nm = Nlist(node_nm,j,2)
         
         do p = 1,3
            
            if(p==1) then
               node_elmnt = node_nm
               store(1:2) = delB_delCN(node_nm,j,1:2)
            elseif(p==2) then
               node_elmnt = neigh1_nm
               store(1:2) = delB_delNN1(node_nm,j,1:2)
            elseif(p==3) then
               node_elmnt = neigh2_nm
               store(1:2) = delB_delNN2(node_nm,j,1:2)
            endif
            
            write(224,fmt=*) store(1:2),node_nm,j,p,"store,node,j,p"
            
            if (node_typ(node_elmnt) == 1) then!x,y free
               dum_grdBend(node_elmnt,1:2) = dum_grdBend(node_elmnt,1:2) + store(1:2)       
            elseif (node_typ(node_elmnt) == 2) then!xfree,yfixed
               dum_grdBend(node_elmnt,1) = dum_grdBend(node_elmnt,1) + store(1)
            endif
            
            if (node_nm==19 .and. p==1) then
               write(223,fmt=*) store(1:2),delB_delCN(node_nm,j,1:2),node_nm,"node"
            endif
            
            if (neigh1_nm==19 .and. p==2) then
               write(223,fmt=*) store(1:2),delB_delNN1(node_nm,j,1:2),neigh1_nm,node_nm,"n1,node"
            endif
            
            if (neigh2_nm==19 .and. p==3) then
               write(223,fmt=*) store(1:2),delB_delNN2(node_nm,j,1:2),neigh2_nm,node_nm,"n2,node"
            endif
            
         enddo
         
      enddo
      
      close(223)
      close(224)
      close(230)
      
    end subroutine grdBend_lp
    
    
  end subroutine gradient_components_frm_bend
  
  
  subroutine calc_delB_delN_array(dum_Nodes) !B=bending
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    integer :: i,j
    integer :: totl_area_cnnctd
    integer :: area_nm,node_nm,area_serial
    real*8  :: all_delB_delN(1:3,1:2)
    real*8  :: ZERO
    
    ZERO = 0.0000000000000000d0
    
    open(unit=122,file='delB_delN.dat')
    
    do i = 1,N_node
       node_nm = i
       totl_area_cnnctd = node_area(i,0)
       write(122,fmt=*) node_nm,"node_nm"
       
       do j = 1,totl_area_cnnctd
          area_nm = node_area(i,j)
          area_serial = j
          write(122,fmt=*) area_nm,area_serial,"area_nm,area_serial"
          
          if (abs(k_phi(node_nm,area_serial)-ZERO) .le. 1.0d-14) then
             all_delB_delN = 0.000000000d0
          else
             all_delB_delN = indvd_delB_delN(area_nm,area_serial,node_nm)
          endif
          
          delB_delCN(i,j,1:N_dmnsn)  = 0.5d0*k_phi(node_nm,area_serial)*all_delB_delN(1,1:2)
          delB_delNN1(i,j,1:N_dmnsn) = 0.5d0*k_phi(node_nm,area_serial)*all_delB_delN(2,1:2)
          delB_delNN2(i,j,1:N_dmnsn) = 0.5d0*k_phi(node_nm,area_serial)*all_delB_delN(3,1:2)
          
          write(122,fmt=*) delB_delCN(i,j,1:2) ,"delB_delCN"
          write(122,fmt=*) delB_delNN1(i,j,1:2),"delB_delNN1"
          write(122,fmt=*) delB_delNN2(i,j,1:2),"delB_delNN2"
          write(122,fmt=*) " "
       enddo
       
       write(122,fmt=*) " "
       
    enddo
    
    close(122)
    
  contains
    
    function indvd_delB_delN(area_nm,area_serial,node_nm)
      implicit none
      integer :: area_nm,area_serial,node_nm
      real*8  :: indvd_delB_delN(1:3,1:N_dmnsn)
      integer :: NeighNodes(1:2)
      real*8  :: NN1(1:2),NN2(1:2),CN(1:2)
      real*8  :: Aa,Bb
      real*8  :: x12,y12,x13,y13
      real*8  :: dAdCnX,dAdCnY,dBdCnX,dBdCnY
      real*8  :: dAdNn1X,dAdNn1Y,dBdNn1X,dBdNn1Y
      real*8  :: dAdNn2X,dadNn2Y,dBdNn2X,dBdNn2Y
      real*8  :: rtA1,rtA2
      
      indvd_delB_delN = 0.0d0
      
      NeighNodes(1:2) = Nlist(node_nm,area_serial,1:2)
      write(122,fmt=*) NeighNodes(1:2),"Neighbours"
      
      NN1(1:2) = dum_Nodes(NeighNodes(1),1:2)
      NN2(1:2) = dum_Nodes(NeighNodes(2),1:2)
      CN(1:2)  = dum_Nodes(node_nm,1:2)
      
      call get_A_and_B(CN,NN1,NN2,Aa,Bb)
      call get_x12y12_x13y13(CN,NN1,NN2,x12,y12,x13,y13)
      
      call get_dAdCn_and_dBdCn(Aa,Bb,x12,y12,x13,y13,dAdCnX,dAdCnY,dBdCnX,dBdCnY)
      call get_dAdNn1_and_dBdNn1(Aa,Bb,x12,y12,x13,y13,dAdNn1X,dAdNn1Y,dBdNn1X,dBdNn1Y)
      call get_dAdNn2_and_dBdNn2(Aa,Bb,x12,y12,x13,y13,dAdNn2X,dadNn2Y,dBdNn2X,dBdNn2Y)
      
      if (nodePhi_typ(node_nm)==1) then
         
         indvd_delB_delN(1,1) = (Bb*dAdCnX-Aa*dBdCnX)/((Bb-Aa)**2)
         indvd_delB_delN(1,2) = (Bb*dAdCnY-Aa*dBdCnY)/((Bb-Aa)**2)
         
         indvd_delB_delN(2,1) = (Bb*dAdNn1X-Aa*dBdNn1X)/((Bb-Aa)**2)
         indvd_delB_delN(2,2) = (Bb*dAdNn1Y-Aa*dBdNN1Y)/((Bb-Aa)**2)
         
         indvd_delB_delN(3,1) = (Bb*dAdNn2X-Aa*dBdNn2X)/((Bb-Aa)**2)
         indvd_delB_delN(3,2) = (Bb*dAdNn2Y-Aa*dBdNn2Y)/((Bb-Aa)**2)
         
      elseif (nodePhi_typ(node_nm)==2) then
         
         indvd_delB_delN(1,1) = -(Bb*dAdCnX-Aa*dBdCnX)/((Aa)**2)
         indvd_delB_delN(1,2) = -(Bb*dAdCnY-Aa*dBdCnY)/((Aa)**2)
         
         indvd_delB_delN(2,1) = -(Bb*dAdNn1X-Aa*dBdNn1X)/((Aa)**2)
         indvd_delB_delN(2,2) = -(Bb*dAdNn1Y-Aa*dBdNN1Y)/((Aa)**2)
         
         indvd_delB_delN(3,1) = -(Bb*dAdNn2X-Aa*dBdNn2X)/((Aa)**2)
         indvd_delB_delN(3,2) = -(Bb*dAdNn2Y-Aa*dBdNn2Y)/((Aa)**2)
         
      endif
      
    end function indvd_delB_delN
    
  end subroutine calc_delB_delN_array
  
  subroutine get_x12y12_x13y13(CN,NN1,NN2,x12,y12,x13,y13)
    implicit none
    real*8, intent(in) :: CN(1:2)  !x1,y1
    real*8, intent(in) :: NN1(1:2) !x2,y2
    real*8, intent(in) :: NN2(1:2) !x3,y3
    
    real*8, intent(out) :: x12,y12,x13,y13
    
    x12 = CN(1)-NN1(1)
    y12 = CN(2)-NN1(2)
    x13 = CN(1)-NN2(1)
    y13 = CN(2)-NN2(2)
    
  end subroutine get_x12y12_x13y13
  
  subroutine get_dAdCn_and_dBdCn(Aa,Bb,x12,y12,x13,y13,dAdCnX,dAdCnY,dBdCnX,dBdCnY)
    implicit none !A or Aa here is not area, its numerator of bendE 
    real*8, intent(in)  :: Aa,Bb,x12,y12,x13,y13 
    real*8, intent(out) :: dAdCnx,dAdCnY,dBdCnX,dBdCnY
    
    dAdCnX = 2.0d0*(x12*x13+y12*y13)*(x13+x12)
    dAdCnY = 2.0d0*(x12*x13+y12*y13)*(y13+y12)
    dBdCnX = (x12**2+y12**2)*(2.0d0*x13) + (x13**2+y13**2)*(2.0d0*x12)
    dBdCnY = (x12**2+y12**2)*(2.0d0*y13) + (x13**2+y13**2)*(2.0d0*y12)
    
  end subroutine get_dAdCn_and_dBdCn
  
  subroutine get_dAdNn1_and_dBdNn1(Aa,Bb,x12,y12,x13,y13,dAdNn1X,dAdNn1Y,dBdNn1X,dBdNn1Y)
    implicit none
    real*8, intent(in)  :: Aa,Bb,x12,y12,x13,y13
    real*8, intent(out) :: dAdNn1X,dAdNn1Y,dBdNn1X,dBdNn1Y
    
    dAdNn1X = 2.0d0*(x12*x13+y12*y13)*(-x13)
    dAdNn1Y = 2.0d0*(x12*x13+y12*y13)*(-y13)
    dBdNn1X = (x13**2+y13**2)*(-2.0d0*x12)
    dBdNn1Y = (x13**2+y13**2)*(-2.0d0*y12)
    
  end subroutine get_dAdNn1_and_dBdNn1
  
  subroutine get_dAdNn2_and_dBdNn2(Aa,Bb,x12,y12,x13,y13,dAdNn2X,dAdNn2Y,dBdNn2X,dBdNn2Y)
    implicit none
    real*8, intent(in)  :: Aa,Bb,x12,y12,x13,y13
    real*8, intent(out) :: dAdNn2X,dAdNn2Y,dBdNn2X,dBdNn2Y
    
    dAdNn2X = 2.0d0*(x12*x13+y12*y13)*(-x12)
    dAdNn2Y = 2.0d0*(x12*x13+y12*y13)*(-y12)
    dBdNn2X = (x12**2+y12**2)*(-2.0d0*x13)
    dBdNn2Y = (x12**2+y12**2)*(-2.0d0*y13)
    
  end subroutine get_dAdNn2_and_dBdNn2
  
  
end module gradient_contrb_frm_bending



module l0_A0_variatns
  use system_parameters
  use gradient_contrb_from_spr
  use gradient_contrb_from_area
  
  implicit none
  !real*8, allocatable :: grdmv_l0vartn(:) 
  !real*8, allocatable :: grdmv_A0vartn(:)
  
contains

  ! subroutine allocate_and_initialize_grdmv_l0A0vartn
  !   implicit none
  
  !   allocate(grdmv_l0vartn(1:N_spr)) !may have to use lt N_spr, say Nspr_mv
  !   allocate(grdmv_A0vartn(1:N_cell))
  
  !   grdmv_l0vartn = -1.d30
  !   grdmv_A0vartn = -1.d30
  
  ! end subroutine allocate_and_initialize_grdmv_l0A0vartn
  
  subroutine get_grdmv_l0vartn(dum_nodes,duml0,dum_grdmv_l0vartn)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8  :: duml0(1:N_spr)
    real*8  :: dum_grdmv_l0vartn(1:N_spr)
    real*8  :: lam_l0variatn
    integer :: i,cnt
    
    if (l0_variatn_lgcl .eqv. .True.) then
       lam_l0variatn = 1.0d0
    else
       lam_l0variatn = 0.0d0
    endif
    
    delE_delR = -1.d30
    
    call calc_delE_delR_array(dum_Nodes,duml0)
    
    cnt = 0
    
    do i = 1,N_spr
       
       if (optmSpr(i)==0) then
          continue
       elseif (optmSpr(i)==1) then
          cnt = cnt + 1
          dum_grdmv_l0vartn(cnt) = -delE_delR(i) * optmSpr(i)
          
       elseif (optmSpr(i).ne.0 .AND. optmSpr(i).ne.1) then
          write(*,*) "Neither 0 nor 1, sbrtn:get_grdmv_l0vartn"
          stop
       endif
       
    enddo
    
    
    if (cnt.ne.N_optmSpr) then
       write(*,*) cnt,N_optmSpr,"cnt,N_optmSpr"
       write(*,*) "cnt and N_optmSpr must be EQUAL,sbrtn:get_grdmv_l0vartn"
       stop
    endif
    
  end subroutine get_grdmv_l0vartn
  
  
  subroutine get_grdmv_A0vartn(dum_Nodes,dumA0,dum_grdmv_A0vartn)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8  :: dumA0(1:N_cell)
    real*8  :: dum_grdmv_A0vartn(1:N_cell)
    real*8  :: lam_A0variatn
    integer :: i,cnt
    
    if (A0_variatn_lgcl .eqv. .True.) then
       lam_A0variatn = 1.0d0
    else
       lam_A0variatn = 0.0d0
    endif
    
    
    delE_delA = -1.d30
    
    call calc_delE_delA_array(dum_Nodes,dumA0)
    
    cnt = 0
    
    do i = 1,N_cell
       
       if (optmCell(i)==0) then
          continue
       elseif (optmCell(i)==1) then
          cnt = cnt + 1
          dum_grdmv_A0vartn(cnt) = -delE_delA(i) * optmCell(i)
          
       elseif (optmCell(i).ne.0 .AND. optmCell(i).ne.1) then
          write(*,*) "cnt and N_optmCell must be EQUAL,sbrtn:get_grdmv_A0vartn"
          stop
       endif
       
    enddo
    
    if (cnt.ne.N_optmCell) then
       write(*,*) "cnt and N_optmCell must be EQUAL,sbrtn:get_grdmv_A0vartn"
       stop
    endif
    
  end subroutine get_grdmv_A0vartn
  
end module l0_A0_variatns


module gradient_module  !M4
  use system_parameters
  use gradient_from_total_system
  use gradient_contrb_from_spr
  use gradient_contrb_from_area
  use gradient_contrb_frm_gravitation
  use gradient_contrb_frm_bending
  use gradient_contrb_frm_SRYP
  use l0_A0_variatns
  use conversion_routines
  use Energy_module
  
  implicit none  
  
contains
  
  subroutine input_from_bash_about_grdnt_choice
    implicit none
    
    read(*,*) Analtcl_or_Numrcl
    
  end subroutine input_from_bash_about_grdnt_choice
  
  subroutine get_all_gradient_variables
    implicit none
    
    call allocate_and_initialize_gradient_system_variables
    call allocate_and_initialize_gradient_spr_variables
    call allocate_and_initialize_gradient_area_variables
    call allocate_and_initialize_gradient_grvtnl
    call allocate_and_initialize_gradient_bend_variables
    
    if (SRyp_lgcl.eqv..True.) call allocate_and_initialize_gradient_SRyp
    
  end subroutine get_all_gradient_variables
  
  subroutine get_all_gradient_variables_wo_StrVars
    implicit none
    
    call allocate_and_initialize_gradient_system_variables_wo_StrVars
    call allocate_and_initialize_gradient_spr_variables_wo_StrVars
    call allocate_and_initialize_gradient_area_variables_wo_StrVars
    call allocate_and_initialize_gradient_grvtnl_wo_StrVars
    call allocate_and_initialize_gradient_bend_variables_wo_StrVars
    
    if (SRyp_lgcl.eqv..True.) call allocate_and_initialize_gradient_SRyp_wo_StrVars
    
  end subroutine get_all_gradient_variables_wo_StrVars
  
  
  subroutine get_all_gradient_variables_StrVars
    implicit none
    
    call allocate_and_initialize_gradient_system_variables_StrVars
    call allocate_and_initialize_gradient_spr_variables_StrVars
    call allocate_and_initialize_gradient_area_variables_StrVars
    call allocate_and_initialize_gradient_grvtnl_StrVars
    call allocate_and_initialize_gradient_bend_variables_StrVars
    
    if (SRyp_lgcl.eqv..True.) call allocate_and_initialize_gradient_SRyp_StrVars
    
  end subroutine get_all_gradient_variables_StrVars
  
  
  subroutine deallocate_all_gradient_variables
    implicit none
    
    call deallocate_gradient_system_variables
    call deallocate_gradient_spr_variables
    call deallocate_gradient_area_variables
    call deallocate_gradient_grvtnl  
    call deallocate_gradient_bend_variables
    
    if (SRyp_lgcl.eqv..True.) call deallocate_gradient_SRyp
    
  end subroutine deallocate_all_gradient_variables
  
  subroutine deallocate_all_gradient_variables_wo_StrVars
    implicit none
    
    call deallocate_gradient_system_variables_wo_StrVars
    call deallocate_gradient_spr_variables_wo_StrVars
    call deallocate_gradient_area_variables_wo_StrVars
    call deallocate_gradient_grvtnl_wo_StrVars  
    call deallocate_gradient_bend_variables_wo_StrVars
    
    if (SRyp_lgcl.eqv..True.) call deallocate_gradient_SRyp_wo_StrVars
    
  end subroutine deallocate_all_gradient_variables_wo_StrVars
  
  
  subroutine deallocate_all_gradient_variables_StrVars
    implicit none
    
    call deallocate_gradient_system_variables_StrVars
    call deallocate_gradient_spr_variables_StrVars
    call deallocate_gradient_area_variables_StrVars
    call deallocate_gradient_grvtnl_StrVars  
    call deallocate_gradient_bend_variables_StrVars
    
    if (SRyp_lgcl.eqv..True.) call deallocate_gradient_SRyp_StrVars
    
  end subroutine deallocate_all_gradient_variables_StrVars
  
  
  subroutine get_gradient(dum_Nodes,duml0,dumA0,dum_grdN)
    implicit none
    real*8 :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8 :: duml0(1:N_spr)
    real*8 :: dumA0(1:N_cell)
    real*8 :: dum_grdN(1:N_node,1:N_dmnsn)
    
    !write(*,*) dum_Nodes,"node_xy should match"
    !write(*,*) dum_grdN, "gradient"
    
    if (Analtcl_or_Numrcl==1) call gradient_Analtcl(dum_Nodes,duml0,dumA0,dum_grdN)
    if (Analtcl_or_Numrcl==2) call gradient_Numrcl(dum_Nodes,dum_grdN)
    
  end subroutine get_gradient
  
  subroutine get_gradient_with_changing_nodeTyp_test(dum_Nodes,duml0,dumA0,dum_grdN)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8  :: duml0(1:N_spr)
    real*8  :: dumA0(1:N_cell)
    real*8  :: dum_grdN(1:N_node,1:N_dmnsn)
    integer :: node_typValStr(1:N_node)
    
    node_typValStr(1:N_node) = node_typ(1:N_node)
    node_typ(1:N_node)       = 1                   ! making all nodes free
    
    if (Analtcl_or_Numrcl==1) call gradient_Analtcl(dum_Nodes,duml0,dumA0,dum_grdN)
    if (Analtcl_or_Numrcl==2) call gradient_Numrcl(dum_Nodes,dum_grdN)
    
    node_typ(1:N_node)       = node_typValStr(1:N_node) ! changing the node_typ_to_previous value
    
  end subroutine get_gradient_with_changing_nodeTyp_test
  
  subroutine get_gradient_with_changing_nodeTyp_test_withCompnts(dum_Nodes,duml0,dumA0,dum_grdN,&
       dum_grdSpr,dum_grdArea,dum_grdGrvtnl,dum_grdBend,dum_grdSRyp)
    implicit none
    real*8, intent(in)  :: dum_Nodes(1:N_node,1:N_dmnsn),duml0(1:N_spr),dumA0(1:N_cell)
    real*8, intent(out) :: dum_grdN(1:N_node,1:N_dmnsn)
    real*8, intent(out) :: dum_grdSpr(1:N_node,1:N_dmnsn)
    real*8, intent(out) :: dum_grdArea(1:N_node,1:N_dmnsn)
    real*8, intent(out) :: dum_grdGrvtnl(1:N_node,1:N_dmnsn)
    real*8, intent(out) :: dum_grdBend(1:N_node,1:N_dmnsn)
    real*8, intent(out) :: dum_grdSRyp(1:N_node,1:N_dmnsn)
    
    integer             :: node_typValStr(1:N_node)
    
    node_typValStr(1:N_node) = node_typ(1:N_node)
    node_typ(1:N_node)       = 1                   ! making all nodes free
    
    call gradient_Analtcl_withCompnts(dum_Nodes,duml0,dumA0,dum_grdN,dum_grdSpr,&
         dum_grdArea,dum_grdGrvtnl,dum_grdBend,dum_grdSRyp)
    
    node_typ(1:N_node)       = node_typValStr(1:N_node) ! changing the node_typ_to_previous value
    
  end subroutine get_gradient_with_changing_nodeTyp_test_withCompnts
  
  
  
  subroutine complte_the_grdnts(dum_grd)
    implicit none
    real*8, intent(inout)  :: dum_grd(1:N_node,1:N_dmnsn)
    integer :: node_toBe_used, node_toBe_filled
    integer :: i
    
    !write(*,*) N_doubleNode,"N_doubleNode"
    
    do i = 1,N_doubleNode
       node_toBe_used   = double_Node(i,1)
       node_toBe_filled = double_Node(i,2)
       !write(*,*) node_toBe_used,node_toBe_filled,"used-filled"
       
       !write(*,*) dum_grd(node_toBe_used,1:2),"toBe_used gradient"
       !write(*,*) dum_grd(node_toBe_filled,1:2),"toBe_filled_gradient_bfr"
       
       dum_grd(node_toBe_filled,1:2) = dum_grd(node_toBe_used,1:2)
       
       !write(*,*) dum_grd(node_toBe_filled,1:2),"toBe_filled_gradient_aft"
       
    enddo
    
    
  end subroutine complte_the_grdnts
  
  
  subroutine gradient_Analtcl_withCompnts(dum_Nodes,duml0,dumA0,dum_grdN,dum_grdSpr,&
       dum_grdArea,dum_grdGrvtnl,dum_grdBend,dum_grdSRyp)
    implicit none
    real*8, intent(in)  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8, intent(in)  :: duml0(1:N_spr),dumA0(1:N_cell)
    real*8, intent(out) :: dum_grdN(1:N_node,1:N_dmnsn)
    real*8, intent(out) :: dum_grdSpr(1:N_node,1:N_dmnsn)
    real*8, intent(out) :: dum_grdArea(1:N_node,1:N_dmnsn)
    real*8, intent(out) :: dum_grdGrvtnl(1:N_node,1:N_dmnsn)
    real*8, intent(out) :: dum_grdBend(1:N_node,1:N_dmnsn)
    real*8, intent(out) :: dum_grdSRyp(1:N_node,1:N_dmnsn)
    
    integer :: i,j
    real*8  :: lam_spr,lam_area,lam_grvtnl,lam_bend,lam_SRyp
        
    dum_grdN      = 0.0d0
    dum_grdSpr    = 0.0d0
    dum_grdArea   = 0.0d0
    dum_grdGrvtnl = 0.0d0
    dum_grdBend   = 0.0d0
    dum_grdSRyp   = 0.0d0
    
    lam_spr=0.0d0 ; lam_area=0.0d0 ; lam_grvtnl=0.0d0 ; lam_bend=0.0d0 ; lam_SRyp=0.0d0
    
    if (spr_lgcl.eqv..True.) then
       lam_spr    = 1.0d0
       call gradient_components_frm_spr(dum_Nodes,duml0,dum_grdSpr)
       call complte_the_grdnts(dum_grdSpr)
    endif
    
    if (area_lgcl.eqv..True.) then
       lam_area   = 1.0d0
       call gradient_components_frm_area(dum_Nodes,dumA0,dum_grdArea)
       call complte_the_grdnts(dum_grdArea)
    endif
    
    if (grvtnl_lgcl.eqv..True.) then
       lam_grvtnl = 1.0d0
       call gradient_components_frm_grvtnl(dum_Nodes,dum_grdGrvtnl)
       call complte_the_grdnts(dum_grdGrvtnl)
    endif
    
    if (bend_lgcl.eqv..True.) then
       lam_bend = 1.0d0
       call gradient_components_frm_bend(dum_Nodes,dum_grdBend)
       call complte_the_grdnts(dum_grdBend)
    endif
    
    if (SRyp_lgcl.eqv..True.) then
       lam_SRyp = 1.0d0
       call gradient_components_frm_SRyp(dum_Nodes,dum_grdSRyp)
       call complte_the_grdnts(dum_grdSRyp)
    endif
    
    do i = 1,N_node
       do j = 1,N_dmnsn
          dum_grdN(i,j) = lam_spr*dum_grdSpr(i,j) + lam_area*dum_grdArea(i,j) &
               + lam_grvtnl*dum_grdGrvtnl(i,j) + lam_bend*dum_grdBend(i,j) &
               + lam_SRyp*dum_grdSRyp(i,j)
       enddo
    enddo
    
    call get_the_other_dn_grd_values_inCompAnltcl
    
  contains
    
    subroutine get_the_other_dn_grd_values_inCompAnltcl
      implicit none
      integer :: node_nmbr,other_node
      
      do node_nmbr = 1,N_node
         if(node_cnnctd(node_nmbr).ne.0 .and. count_this_dn(node_nmbr).eq.1) then
            other_node                    = node_cnnctd(node_nmbr)
            dum_grdN(other_node,1:2)      = dum_grdN(node_nmbr,1:2)
            dum_grdSpr(other_node,1:2)    = dum_grdSpr(node_nmbr,1:2)
            dum_grdArea(other_node,1:2)   = dum_grdArea(node_nmbr,1:2)
            dum_grdGrvtnl(other_node,1:2) = dum_grdGrvtnl(node_nmbr,1:2)
            dum_grdSRyp(other_node,1:2)   = dum_grdSRyp(node_nmbr,1:2)
         endif
      enddo
      
    end subroutine get_the_other_dn_grd_values_inCompAnltcl
    
  end subroutine gradient_Analtcl_withCompnts
  
  
  subroutine gradient_Analtcl(dum_Nodes,duml0,dumA0,dum_grdN)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8  :: duml0(1:N_spr),dumA0(1:N_cell)
    
    real*8  :: dum_grdN(1:N_node,1:N_dmnsn)
    real*8  :: dum_grdSpr(1:N_node,1:N_dmnsn)
    real*8  :: dum_grdArea(1:N_node,1:N_dmnsn)
    real*8  :: dum_grdGrvtnl(1:N_node,1:N_dmnsn)
    real*8  :: dum_grdBend(1:N_node,1:N_dmnsn)
    real*8  :: dum_grdSRyp(1:N_node,1:N_dmnsn)
    
    integer :: i,j
    real*8  :: lam_spr,lam_area,lam_grvtnl,lam_bend,lam_SRyp
    
    dum_grdN      = 0.0d0
    dum_grdSpr    = 0.0d0
    dum_grdArea   = 0.0d0
    dum_grdGrvtnl = 0.0d0
    dum_grdBend   = 0.0d0
    dum_grdSRyp   = 0.0d0
    
    lam_spr=0.00d0; lam_area=0.00d0; lam_grvtnl=0.00d0; lam_bend=0.00d0; lam_SRyp=0.00d0
    
    if (spr_lgcl.eqv..True.) then
       lam_spr    = 1.0d0
       call gradient_components_frm_spr(dum_Nodes,duml0,dum_grdSpr)
       call complte_the_grdnts(dum_grdSpr)
    endif
    
    if (area_lgcl.eqv..True.) then
       lam_area   = 1.0d0
       call gradient_components_frm_area(dum_Nodes,dumA0,dum_grdArea)
       call complte_the_grdnts(dum_grdArea)
    endif
    
    if (grvtnl_lgcl.eqv..True.) then
       lam_grvtnl = 1.0d0
       call gradient_components_frm_grvtnl(dum_Nodes,dum_grdGrvtnl)
       call complte_the_grdnts(dum_grdGrvtnl)
    endif
    
    if (bend_lgcl.eqv..True.) then
       lam_bend = 1.0d0
       call gradient_components_frm_bend(dum_Nodes,dum_grdBend)
       call complte_the_grdnts(dum_grdBend)
    endif
    
    if (SRyp_lgcl.eqv..True.) then   
       lam_SRyp = 1.0d0
       call gradient_components_frm_SRyp(dum_Nodes,dum_grdSRyp)
       call complte_the_grdnts(dum_grdSRyp)
    endif
    
    do i = 1,N_node
       do j = 1,N_dmnsn
          dum_grdN(i,j) = lam_spr*dum_grdSpr(i,j) + lam_area*dum_grdArea(i,j) &
               + lam_grvtnl*dum_grdGrvtnl(i,j) + lam_bend*dum_grdBend(i,j) &
               + lam_SRyp*dum_grdSRyp(i,j)
       enddo
    enddo
    
    call get_the_other_dn_grd_values ! ;write(*,*) dum_grdN,"dum_grdN_INSIDE"
    call print_analytical_gradients
    
  contains
    
    subroutine get_the_other_dn_grd_values
      implicit none
      integer :: node_nmbr,other_node
      
      do node_nmbr = 1,N_node
         if(node_cnnctd(node_nmbr).ne.0 .and. count_this_dn(node_nmbr).eq.1) then
            other_node = node_cnnctd(node_nmbr)
            dum_grdN(other_node,1:2) = dum_grdN(node_nmbr,1:2)
         endif
      enddo
      
    end subroutine get_the_other_dn_grd_values
    
    subroutine print_analytical_gradients
      implicit none
      integer :: i1
      
      open(unit=25,file='grdSpr_grdArea_grdGrvtnl_grdBend.dat')
      open(unit=26,file='grdAna.dat')
      open(unit=27,file='grdSprAna.dat')
      open(unit=28,file='grdAreaAna.dat')
      open(unit=29,file='grdGrvAna.dat')
      open(unit=30,file='grdBendAna.dat')
      open(unit=31,file='grdSRypAna.dat')
      
      write(25,*) "grdSpr(1:2),grdArea(1:2),grdGrvtnl(1:2),grdBend(1:2)"
      
      do i1 = 1,N_node
         
         write(25,*) dum_grdSpr(i1,1:2),dum_grdArea(i1,1:2),dum_grdGrvtnl(i1,1:2),&
              dum_grdBend(i1,1:2),dum_grdSRyp(i1,1:2)
         
         write(25,*) dum_grdN(i1,1:2),"grdN=grdSpr+grdArea+grdGrvtnl+grdBend+grdSRyp"
         write(26,*) dum_grdN(i1,1:2),i1,"grd_Ana,node"
         write(27,*) dum_grdSpr(i1,1:2),i1,"grd_AnaSpr,node"
         write(28,*) dum_grdArea(i1,1:2),i1,"grd_AnaArea,node"
         write(29,*) dum_grdGrvtnl(i1,1:2),i1,"grd_AnaGrv,node"
         write(30,*) dum_grdBend(i1,1:2),i1,"grd_AnaBend,node"
         write(31,*) dum_grdSRyp(i1,1:2),i1,"grd_AnaSRyp,node"
         
      enddo
      
      close(25)
      close(26)
      close(27)
      close(28)
      close(29)
      close(30)
      close(31)
      
    end subroutine print_analytical_gradients
    
    
  end subroutine gradient_Analtcl
  
  
  
  subroutine gradient_Numrcl(dum_Nodes,dum_grdN)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8  :: dum_grdN(1:N_node,1:N_dmnsn)
    real*8  :: dum_grdSpr(1:N_node,1:N_dmnsn)
    real*8  :: dum_grdArea(1:N_node,1:N_dmnsn)
    real*8  :: dum_grdGrvtnl(1:N_node,1:N_dmnsn)
    real*8  :: dum_grdBend(1:N_node,1:N_dmnsn)
    real*8  :: dum_grdSRyp(1:N_node,1:N_dmnsn)
    
    real*8  :: str,str1,str2
    real*8  :: E_aft,Es_aft,Ea_aft,Eg_aft,Eb_aft,Er_aft
    real*8  :: E_bfr,Es_bfr,Ea_bfr,Eg_bfr,Eb_bfr,Er_bfr
    real*8  :: h=1.0d-4
    integer :: i,j,cnt
    logical :: lgcl_dn
    integer :: node_nm,othr_dn
    
    open(unit=131,file='grdNumCheck.dat')
    
    dum_grdN(1:N_node,1:N_dmnsn)      = 0.0000d0
    dum_grdSpr(1:N_node,1:N_dmnsn)    = 0.0000d0
    dum_grdArea(1:N_node,1:N_dmnsn)   = 0.0000d0
    dum_grdGrvtnl(1:N_node,1:N_dmnsn) = 0.0000d0
    dum_grdBend(1:N_node,1:N_dmnsn)   = 0.0000d0
    dum_grdSRyp(1:N_node,1:N_dmnsn)   = 0.0000d0
    
    do i = 1,N_node
       if (node_typ(i)==0) cycle
        
       lgcl_dn = .False.
       node_nm = i ; call is_it_a_dn(node_nm,lgcl_dn)
       if (lgcl_dn .eqv. .True.) call get_the_other_dn(node_nm,othr_dn)
       
       do j = 1,N_dmnsn
          
          if (node_typ(i)==2 .and. j==N_dmnsn) cycle
          
          if (lgcl_dn.eqv..True.) then
             str1 = dum_Nodes(node_nm,j)
             str2 = dum_Nodes(othr_dn,j)
             
             dum_Nodes(node_nm,j) = dum_Nodes(node_nm,j) + h
             dum_Nodes(othr_dn,j) = dum_Nodes(othr_dn,j) + h
             
          elseif (lgcl_dn.eqv..False.) then
             str = dum_Nodes(i,j)
             dum_Nodes(i,j) = dum_Nodes(i,j) + h
          endif
             
          E_aft  = Energy(dum_Nodes,l0,A0)
          Es_aft = spr_E(dum_Nodes,l0)
          Ea_aft = area_E(dum_Nodes,A0)
          Eg_aft = grvtnl_E(dum_Nodes)
          Eb_aft = bend_E(dum_Nodes)
          Er_aft = SRyp_E(dum_Nodes)
          
          if (lgcl_dn.eqv..True.) then
             dum_Nodes(node_nm,j) = str1
             dum_Nodes(othr_dn,j) = str2
             
             dum_Nodes(node_nm,j) = dum_Nodes(node_nm,j) - h
             dum_Nodes(othr_dn,j) = dum_Nodes(othr_dn,j) - h
             
          elseif (lgcl_dn.eqv..False.) then
             dum_Nodes(i,j) = str
             dum_Nodes(i,j) = dum_Nodes(i,j) - h
          endif
          
          
          E_bfr  = Energy(dum_Nodes,l0,A0)
          Es_bfr = spr_E(dum_Nodes,l0)
          Ea_bfr = area_E(dum_Nodes,A0)
          Eg_bfr = grvtnl_E(dum_Nodes)
          Eb_bfr = bend_E(dum_Nodes)
          Er_bfr = SRyp_E(dum_Nodes)
          
          if (lgcl_dn.eqv..True.) then
             dum_Nodes(node_nm,j) = str1
             dum_Nodes(othr_dn,j) = str2   
          elseif (lgcl_dn .eqv. .False.) then
             dum_Nodes(i,j) = str
          endif
          
          dum_grdN(i,j)      = (E_aft-E_bfr)/(2.d0*h)
          dum_grdSpr(i,j)    = (Es_aft-Es_bfr)/(2.d0*h)
          dum_grdArea(i,j)   = (Ea_aft-Ea_bfr)/(2.d0*h)
          dum_grdGrvtnl(i,j) = (Eg_aft-Eg_bfr)/(2.d0*h)
          dum_grdBend(i,j)   = (Eb_aft-Eb_bfr)/(2.d0*h)
          dum_grdSRyp(i,j)   = (Er_aft-Er_bfr)/(2.d0*h)
          
          if (i==31) then
             write(131,*) E_aft,E_bfr,"Eaft,Ebfr"
             write(131,*) Es_aft,Es_bfr,"Es_aft,Es_bfr"
             write(131,*) Ea_aft,Ea_bfr,"Ea_aft,Ea_bfr"
             write(131,*) Eg_aft,Eg_bfr,"Eg_aft,Eg_bfr"
             write(131,*) Eb_aft,Eb_bfr,"Eb_aft,Eb_bfr"
             write(131,*) Er_aft,Er_bfr,"Er_aft,Er_bfr"
             write(131,*) dum_grdN(i,j),i,j,"grdN,node,j"
          endif
          
       enddo
       
    enddo
    
    close(131)
    
    call print_numerical_gradients
    
  contains
    
    subroutine print_numerical_gradients
      implicit none
      integer :: i
      
      open(unit=126,file='grdNum.dat')
      open(unit=127,file='grdSprNum.dat')
      open(unit=128,file='grdAreaNum.dat')
      open(unit=129,file='grdGrvNum.dat')
      open(unit=130,file='grdBendNum.dat')
      open(unit=119,file='grdSRypNum.dat')
      
      do i = 1,N_node
         write(126,*) dum_grdN(i,1:2),i,"grd_Num,node"
         write(127,*) dum_grdSpr(i,1:2),i,"grd_NumSpr,node"
         write(128,*) dum_grdArea(i,1:2),i,"grd_NumArea,node"
         write(129,*) dum_grdGrvtnl(i,1:2),i,"grd_NumGrv,node"
         write(130,*) dum_grdBend(i,1:2),i,"grd_NumBend,node"
         write(119,*) dum_grdSRyp(i,1:2),i,"grd_NumSRyp,node"
      enddo
      
      close(126)
      close(127)
      close(128)
      close(129)
      close(130)
      close(119)
      
    end subroutine print_numerical_gradients
    
  end subroutine gradient_Numrcl
  
  
  
  subroutine gradient_Numrcl_withCompnts(dum_Nodes,duml0,dumA0,dum_grdN,dum_grdSpr,&
       dum_grdArea,dum_grdGrvtnl,dum_grdBend,dum_grdSRyp)
    implicit none
    real*8, intent(inout) :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8, intent(in)    :: duml0(1:N_spr),dumA0(1:N_cell)
    real*8, intent(out)   :: dum_grdN(1:N_node,1:N_dmnsn)
    real*8, intent(out)   :: dum_grdSpr(1:N_node,1:N_dmnsn)
    real*8, intent(out)   :: dum_grdArea(1:N_node,1:N_dmnsn)
    real*8, intent(out)   :: dum_grdGrvtnl(1:N_node,1:N_dmnsn)
    real*8, intent(out)   :: dum_grdBend(1:N_node,1:N_dmnsn)
    real*8, intent(out)   :: dum_grdSRyp(1:N_node,1:N_dmnsn)
    
    real*8  :: str,str1,str2
    real*8  :: E_aft,Es_aft,Ea_aft,Eg_aft,Eb_aft,Er_aft
    real*8  :: E_bfr,Es_bfr,Ea_bfr,Eg_bfr,Eb_bfr,Er_bfr
    real*8  :: h=1.0d-4
    
    integer :: i,j,cnt
    logical :: lgcl_dn
    integer :: node_nm,othr_dn
    
    dum_grdN(1:N_node,1:N_dmnsn)      = 0.0d0
    dum_grdSpr(1:N_node,1:N_dmnsn)    = 0.0d0
    dum_grdArea(1:N_node,1:N_dmnsn)   = 0.0d0
    dum_grdGrvtnl(1:N_node,1:N_dmnsn) = 0.0d0
    dum_grdBend(1:N_node,1:N_dmnsn)   = 0.0d0
    dum_grdSRyp(1:N_node,1:N_dmnsn)   = 0.0d0
    
    do i = 1,N_node
       if (node_typ(i)==0) cycle
       
       lgcl_dn = .False.
       
       node_nm = i
       call is_it_a_dn(node_nm,lgcl_dn)
       
       if (lgcl_dn .eqv. .True.) then
          call get_the_other_dn(node_nm,othr_dn)
       endif
       
       do j = 1,N_dmnsn
          
          if (node_typ(i)==2 .and. j==N_dmnsn) cycle
          
          if (lgcl_dn.eqv..True.) then
             str1 = dum_Nodes(node_nm,j)
             str2 = dum_Nodes(othr_dn,j)
             
             dum_Nodes(node_nm,j) = dum_Nodes(node_nm,j) + h
             dum_Nodes(othr_dn,j) = dum_Nodes(othr_dn,j) + h
             
          elseif (lgcl_dn.eqv..False.) then
             str = dum_Nodes(i,j)
             dum_Nodes(i,j) = dum_Nodes(i,j) + h
          endif
             
          E_aft  = Energy(dum_Nodes,l0,A0)
          Es_aft = spr_E(dum_Nodes,l0)
          Ea_aft = area_E(dum_Nodes,A0)
          Eg_aft = grvtnl_E(dum_Nodes)
          Eb_aft = bend_E(dum_Nodes)
          Er_aft = SRyp_E(dum_Nodes)
          
          if (lgcl_dn.eqv..True.) then
             dum_Nodes(node_nm,j) = str1
             dum_Nodes(othr_dn,j) = str2
             
             dum_Nodes(node_nm,j) = dum_Nodes(node_nm,j) - h
             dum_Nodes(othr_dn,j) = dum_Nodes(othr_dn,j) - h
             
          elseif (lgcl_dn.eqv..False.) then
             dum_Nodes(i,j) = str
             dum_Nodes(i,j) = dum_Nodes(i,j) - h
          endif
          
          E_bfr  = Energy(dum_Nodes,l0,A0)
          Es_bfr = spr_E(dum_Nodes,l0)
          Ea_bfr = area_E(dum_Nodes,A0)
          Eg_bfr = grvtnl_E(dum_Nodes)
          Eb_bfr = bend_E(dum_Nodes)
          Er_bfr = SRyp_E(dum_Nodes)
          
          if (lgcl_dn.eqv..True.) then
             dum_Nodes(node_nm,j) = str1
             dum_Nodes(othr_dn,j) = str2   
          elseif (lgcl_dn .eqv. .False.) then
             dum_Nodes(i,j) = str
          endif
          
          dum_grdN(i,j)      = (E_aft-E_bfr)/(2.d0*h)
          dum_grdSpr(i,j)    = (Es_aft-Es_bfr)/(2.d0*h)
          dum_grdArea(i,j)   = (Ea_aft-Ea_bfr)/(2.d0*h)
          dum_grdGrvtnl(i,j) = (Eg_aft-Eg_bfr)/(2.d0*h)
          dum_grdBend(i,j)   = (Eb_aft-Eb_bfr)/(2.d0*h)
          dum_grdSRyp(i,j)   = (Er_aft-Er_bfr)/(2.d0*h)
          
       enddo
       
    enddo
    
  end subroutine gradient_Numrcl_withCompnts
  
  
  
  !VECTOR TO COORDINATE TRANSFRM STUFFS FOR GRADIENT ARE BELOW
  
  subroutine gradient_coordntes(dum_coordntes,dum_grdC)
    
    implicit none
    real*8,intent(in)  :: dum_coordntes(1:N_variabls)
    real*8,intent(out) :: dum_grdC(1:N_variabls)!C=coordnte_based
    
    real*8  :: dum_grdC_xy(1:N_mvCoordnte)
    real*8  :: dum_grdC_spr(1:N_spr),dum_grdC_area(1:N_cell)
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8  :: duml0(1:N_spr),dumA0(1:N_cell)
    real*8  :: dum_grdN(1:N_node,1:N_dmnsn)!N=node_based
    integer :: cnt
    
    if (funcChoice==1) then
       call coordntes_to_nodes(dum_coordntes,dum_Nodes)
    elseif (funcChoice==2) then
       call coordntes_to_nodesl0A0(dum_coordntes,dum_Nodes,duml0,dumA0)
       
    elseif (funcChoice.ne.1 .AND. funcChoice.ne.2) then
       write(*,*) "N_variabls not equal to any of vars"
       write(*,*) "flnm:gradient_calcEnrgy.f08,sbrtn:gradient_coordntes"
       stop
    endif
    
    !write(*,*) dum_coordntes,"dum_coordntes inside_gradient_coordntes_sbrtn" !CO1
    !write(*,*) dum_Nodes,"dum_Nodes inside_gradient_coordntes_sbrtn"!CO2
    
    
    if (funcChoice==1) then
       call get_gradient(dum_Nodes,l0,A0,dum_grdN)
    elseif (funcChoice==2) then
       call get_gradient(dum_Nodes,duml0,dumA0,dum_grdN)
    endif
    !write(*,*) dum_grdN,"dum_grdN"!CO3
    
       
    call reduced_to_grd_coordntes_xy(dum_grdN,dum_grdC_xy)
    !write(*,*) dum_grdC,"dum_grdC_out"!CO4
    
    if (funcChoice==2) then
       if (l0_variatn_lgcl.eqv..True.) call get_grdmv_l0vartn(dum_Nodes,duml0,dum_grdC_spr)
       if (A0_variatn_lgcl.eqv..True.) call get_grdmv_A0vartn(dum_Nodes,dumA0,dum_grdC_area)
    endif
    
    call get_complete_grdC
    call print_coordnte_based_gradients(dum_grdC)
    
    !if (CellsMeet==1) then
    !   write(*,*) dum_grdN(1,1:2),dum_grdN(2,1:2),"1-2"
    !   write(*,*) dum_grdN(3,1:2),dum_grdN(4,1:2),"3-4"
    !   write(*,*) dum_grdC(1:6),"aft print_coordnte_based_gradients"
    !endif
    
  contains
    
    subroutine reduced_to_grd_coordntes_xy(dum_grdN,dum_grdC)
      implicit none
      real*8  :: dum_grdN(1:N_node,1:N_dmnsn)
      real*8  :: dum_grdC(1:N_mvCoordnte)
      integer :: i,count_mv
      
      count_mv = 1
      !write(*,*) N_mvCoordnte,"N_mvCoordnte"
      
      dum_grdC(1:N_mvCoordnte) = -1.d30
      
      do i = 1,N_node
         
         if (node_typ(i) .eq. 0) then !fixed node
            continue
         elseif (node_typ(i) .ne. 0) then
            if (node_cnnctd(i) .eq. 0) then
               
               if (node_typ(i) .eq. 1) then !free node
                  dum_grdC(count_mv:(count_mv+1)) = dum_grdN(i,1:N_dmnsn)
                  count_mv = count_mv + 2
                  
               elseif (node_typ(i) .eq. 2) then !x free, y fixed
                  dum_grdC(count_mv) = dum_grdN(i,1)
                  count_mv = count_mv + 1
               endif
               
            elseif (node_cnnctd(i) .ne. 0) then
               if (count_this_dn(i) .eq. 1) then
                  
                  if (node_typ(i) .eq. 1) then !free node
                     dum_grdC(count_mv:(count_mv+1)) = dum_grdN(i,1:N_dmnsn)
                     count_mv = count_mv + 2
                     
                  elseif (node_typ(i) .eq. 2) then !x free, y fixed
                     dum_grdC(count_mv) = dum_grdN(i,1)
                     count_mv = count_mv + 1
                  endif
                  
               elseif (count_this_dn(i) .ne. 1) then
                  continue
               endif
               
            endif
         endif
         
      enddo
      
      
    end subroutine reduced_to_grd_coordntes_xy
    
    subroutine get_complete_grdC
      implicit none
      integer :: cnt
      integer :: lam_l0variatn=0,lam_A0variatn=0
      integer :: i
      
      dum_grdC = -1.d20
      dum_grdC(1:N_mvCoordnte) = dum_grdC_xy(1:N_mvCoordnte)
      cnt = N_mvCoordnte
      
      if (funcChoice==2) then
         
         if (l0_variatn_lgcl .eqv. .True.) then
            
            do i = 1,N_spr
               if (optmSpr(i)==1) then
                  cnt = cnt + 1
                  dum_grdC(cnt) = dum_grdC_spr(i)
                  
               elseif (optmSpr(i)==0) then
                  continue
               elseif (optmSpr(i).ne.1 .AND. optmSpr(i).ne.0) then
                  write(*,*) "Neither 0 nor 1,m1,,sbrtn:get_complete_grdC"
                  stop
               endif
            enddo
            
         endif
         !write(*,*) cnt,"cnt aft N_spr"
         
         if (A0_variatn_lgcl .eqv. .True.) then
            
            do i = 1,N_cell
               
               if (optmCell(i)==1) then
                  cnt = cnt + 1
                  dum_grdC(cnt) = dum_grdC_area(i)
                  
               elseif (optmCell(i)==0) then
                  continue
               elseif(optmCell(i).ne.0 .and. optmCell(i).ne.1) then
                  write(*,*) "Neither 0 nor 1, m2,sbrtn:get_complete_grdC"
                  stop
               endif
               
            enddo
            
         endif
         !write(*,*) cnt,"cnt aft N_cell"
         
         if (cnt.ne.N_mvCoordnte_withl0A0) then
            write(*,*) cnt,N_mvCoordnte_withl0A0,"cnt,N_mvCoordnte_withl0A0"
            write(*,*) "cnt must be equal to N_mvCoordnte_withl0A0"
            write(*,*) "Neither 0 nor 1,m3, sbrtn:get_complete_grdC"
            stop
         endif
         
      endif
      
    end subroutine get_complete_grdC
    
    subroutine print_coordnte_based_gradients(dum_grdC)
      implicit none
      real*8  :: dum_grdC(1:N_variabls)
      integer :: i
      
      open(unit=26,file='co_ordnte_based_gradient.dat')
      
      do i = 1,N_variabls
         write(unit=26,fmt=*) dum_grdC(i), "coordnte_based_gradients"
      enddo
      
      close(26)
      
    end subroutine print_coordnte_based_gradients
    
  end subroutine gradient_coordntes
  
  
  subroutine check_values_of_gradient
    implicit none
    
    integer :: fnode,spr_no
    integer :: i
    real*8  :: dx(1:3),dy(1:3),r(1:3)
    real*8  :: grd_prt_spr(1:3,1:2)
    
    integer :: dnode1=0,dnode2=0,area_nm
    real*8  :: dx_a(1:2),dy_a(1:2)
    real*8  :: grd_prt_area(1:2,1:2)
    
    real*8  :: grd_spr(1:2),grd_area(1:2)
    real*8  :: grd(1:2)
    
    grd_spr(1:2)  = 0.0d0
    grd_area(1:2) = 0.0d0
    grd(1:2)      = 0.0d0
    
    open(unit=22, file='checking_gradient_value.dat')
    write(unit=22,fmt=*) "We will check for Node 3"
    
    write(unit=22, fmt=*) node_spr(3,1:max_spr_node),"node_to_spr"
    write(unit=22, fmt=*) k_spr(1),k_spr(3),k_spr(4),"k_spr"
    write(unit=22, fmt=*) l0(1),l0(3),l0(4),"l0"
    
    do i = 1,3
       if (i.eq.1) then
          fnode = 1; spr_no = 1
       elseif(i.eq.2) then
          fnode = 4; spr_no = 3
       else
          fnode = 5; spr_no = 4
       endif
       
       write(unit=22, fmt=*) node_xy(3,1),node_xy(fnode,1),"x value(3 & fnode)"
       write(unit=22, fmt=*) node_xy(3,2),node_xy(fnode,2),"y value(3 & fnode)"
       
       dx(i) = node_xy(3,1) - node_xy(fnode,1)
       dy(i) = node_xy(3,2) - node_xy(fnode,2)
       r(i)  = sqrt((dx(i))**2 + (dy(i))**2)
       
       write(unit=22, fmt=*) dx(i),dy(i),"dx,dy,for i =",i
       write(unit=22, fmt=*) r(i),"r for i =",i
       
       grd_prt_spr(i,1) = k_spr(spr_no) * (r(i)-l0(spr_no)) * (dx(i) / r(i))
       grd_prt_spr(i,2) = k_spr(spr_no) * (r(i)-l0(spr_no)) * (dy(i) / r(i))
       
       write(unit=22,fmt=*) grd_prt_spr(i,1:2),"grd_prt_spr for i =",i
       write(unit=22,fmt=*) " "
       write(unit=22,fmt=*) " "
       
       grd_spr(1) = grd_spr(1) + grd_prt_spr(i,1)
       grd_spr(2) = grd_spr(2) + grd_prt_spr(i,2)
       
    enddo
    
    do i = 1,2
       if (i.eq.1) then
          dnode1  = 1 ; dnode2 = 4
          area_nm = 1          
       elseif (i.eq.2) then
          dnode1  = 4 ; dnode2 = 5
          area_nm = 2
       endif
       
       dy_a(i) = 0.5 * (node_xy(dnode1,2) - node_xy(dnode2,2))
       dx_a(i) = 0.5 * (node_xy(dnode2,1) - node_xy(dnode1,1))
       
       
       !grd_prt_area(i,1) = k_area(area_nm) * area_chng(area_nm) * (dy_a(i))
       !grd_prt_area(i,2) = k_area(area_nm) * area_chng(area_nm) * (dx_a(i))
       
       
       !area_chng can't be accessed HERE as it is in the contains part of a subroutine
       
       
       !write(unit=22,fmt=*) area_chng(area_nm),"area_change for i=",i
       
       write(unit=22,fmt=*) dy_a(i),dx_a(i),"dy_a,dx_a"
       write(unit=22,fmt=*) grd_prt_area(i,1:2),"grd_prt_area for i =",i
       write(unit=22,fmt=*) " "
       write(unit=22,fmt=*) " "
       
    enddo
    
    grd(1) = grd_spr(1) + grd_area(1)
    grd(2) = grd_spr(2) + grd_area(2)
    
    write(unit=22,fmt=*) grd(1:2),"grd"
    
    close(22)
    
  end subroutine check_values_of_gradient
  
  subroutine turn_on_one_strength_property(spr_or_area,spr_or_area_nmbr) !wont be used now
    implicit none
    integer :: spr_or_area
    integer :: spr_or_area_nmbr
    
    if (spr_or_area .eq. 1) then       
       !k_spr(1:N_spr)   = 0.0d0 !you have to store first
       !k_area(1:N_cell) = 0.0d0
       !k_spr(spr_or_area_nmbr) = 1.0d0
       
    elseif (spr_or_area .eq. 2) then
       !k_spr(1:N_spr)   = 0.0d0
       !k_area(1:N_cell) = 0.0d0
       !k_area(spr_or_area_nmbr) = 1.0d0
       
    endif

  end subroutine turn_on_one_strength_property
  
  subroutine Find_Analytical_and_Numerical_Mismatch
    implicit none
    integer :: i
    real*8  :: grd_mvCaseAna(1:N_mvCoordnte),grd_mvCaseNum(1:N_mvCoordnte)
    real*8  :: grd_mvDiff(1:N_mvCoordnte)
    
    grdMismatch = 1
    
    funcChoice=1 ; N_variabls = N_mvCoordnte
    call gradient_coordntes(coordntes_xy,grdmv_xy)
    grd_mv(1:N_mvCoordnte) = grdmv_xy(1:N_mvCoordnte)
    
    grd_mvCaseAna(1:N_mvCoordnte) = grd_mv(1:N_mvCoordnte) !
    
    Analtcl_or_Numrcl = 2
    call gradient_coordntes(coordntes_xy,grdmv_xy)
    grd_mv(1:N_mvCoordnte) = grdmv_xy(1:N_mvCoordnte)
    
    grd_mvCaseNum(1:N_mvCoordnte) = grd_mv(1:N_mvCoordnte)
    grd_mvDiff                    = grd_mvCaseAna - grd_mvCaseNum 
    
    write(*,*) "cheking grd Mismatch"
    write(*,*) " "
    
    do i = 1,N_mvCoordnte
       write(*,*) grd_mvCaseAna(i),grd_mvCaseNum(i),grd_mvDiff(i),i,"grd_chk"
    enddo
    
    write(*,*) " "
    
    Analtcl_or_Numrcl = 1
    
    grdMismatch = 0
    
    !do i = 1,N_node
     !  write(*,*) node_xy(i,1:2),i,"x,y,i"
    !enddo
    
  end subroutine Find_Analytical_and_Numerical_Mismatch
  
  
  subroutine Find_Analytical_and_Numerical_Mismatch_allInfo(cnt)
    implicit none
    integer, intent(in) :: cnt
    
    real*8  :: grdNV1(1:N_node,1:2),grdSprV1(1:N_node,1:2),grdAreaV1(1:N_node,1:2)
    real*8  :: grdGrvtnlV1(1:N_node,1:2),grdBendV1(1:N_node,1:2),grdSRypV1(1:N_node,1:2)
    real*8  :: grdNV2(1:N_node,1:2),grdSprV2(1:N_node,1:2),grdAreaV2(1:N_node,1:2)
    real*8  :: grdGrvtnlV2(1:N_node,1:2),grdBendV2(1:N_node,1:2),grdSRypV2(1:N_node,1:2)
    
    call gradient_Analtcl_withCompnts(node_xy,l0,A0,grdNV1,grdSprV1,grdAreaV1,&
         grdGrvtnlV1,grdBendV1,grdSRypV1)
    call gradient_Numrcl_withCompnts(node_xy,l0,A0,grdNV2,grdSprV2,grdAreaV2,&
         grdGrvtnlV2,grdBendV2,grdSRypV2)
    
    call cmpreTwoGrdCompnts(grdNV1,grdSprV1,grdAreaV1,grdGrvtnlV1,grdBendV1,grdSRypV1,&
          grdNV2,grdSprV2,grdAreaV2,grdGrvtnlV2,grdBendV2,grdSRypV2,cnt)
    
    
  end subroutine Find_Analytical_and_Numerical_Mismatch_allInfo
  
  
  
  subroutine cmpreTwoGrdCompnts(grdV1,grdSprV1,grdAreaV1,grdGrvtnlV1,grdBendV1,&
       grdSrypV1,grdV2,grdSprV2,grdAreaV2,grdGrvtnlV2,grdBendV2,grdSRypV2,cnt)
    implicit none
    real*8,  intent(in) :: grdV1(1:N_node,1:N_dmnsn),grdSprV1(1:N_node,1:N_dmnsn)
    real*8,  intent(in) :: grdAreaV1(1:N_node,1:N_dmnsn),grdGrvtnlV1(1:N_node,1:N_dmnsn)
    real*8,  intent(in) :: grdBendV1(1:N_node,1:N_dmnsn),grdSRypV1(1:N_node,1:N_dmnsn)
    real*8,  intent(in) :: grdV2(1:N_node,1:N_dmnsn),grdSprV2(1:N_node,1:N_dmnsn)
    real*8,  intent(in) :: grdAreaV2(1:N_node,1:N_dmnsn),grdGrvtnlV2(1:N_node,1:N_dmnsn)
    real*8,  intent(in) :: grdBendV2(1:N_node,1:N_dmnsn),grdSRypV2(1:N_node,1:N_dmnsn)
    integer, intent(in) :: cnt
    
    integer             :: i,j,imax,jmax
    real*8              :: diffInComp1,diffInComp2,compChk
    integer             :: comp1Sign,comp2Sign
    
    if (cnt==1)    open(unit=346,file='cmpre_TwoGrdCompnts_001.dat')
    if (cnt==2)    open(unit=346,file='cmpre_TwoGrdCompnts_002.dat')
    if (cnt==3)    open(unit=346,file='cmpre_TwoGrdCompnts_003.dat')
    if (cnt==4)    open(unit=346,file='cmpre_TwoGrdCompnts_004.dat')
    if (cnt==5)    open(unit=346,file='cmpre_TwoGrdCompnts_005.dat')
    if (cnt==6)    open(unit=346,file='cmpre_TwoGrdCompnts_006.dat')
    if (cnt==7)    open(unit=346,file='cmpre_TwoGrdCompnts_007.dat')
    if (cnt==8)    open(unit=346,file='cmpre_TwoGrdCompnts_008.dat')
    if (cnt==9)    open(unit=346,file='cmpre_TwoGrdCompnts_009.dat')
    if (cnt==10)   open(unit=346,file='cmpre_TwoGrdCompnts_010.dat')
    if (cnt.gt.10) stop 'cntFile not in Rnge in cmpreTwogrd'
    
    imax    = 6
    jmax    = N_node
    compChk = 1.00d-6
    
    do i = 1,imax
       
       if (i==1) write(346,*) "FOR combined grd"
       if (i==2) write(346,*) "FOR Spr    Compnent" 
       if (i==3) write(346,*) "FOR Area   Compnent"
       if (i==4) write(346,*) "FOR Grvtnl Compnent"
       if (i==5) write(346,*) "FOR Bend   Compnent"
       if (i==6) write(346,*) "FOR SRyp   Compnent"
       
       do j = 1,jmax
          
          if (i==1) then
             
             diffInComp1 = (grdV1(j,1)-grdV2(j,1))
             diffInComp2 = (grdV1(j,2)-grdV2(j,2))
             
             if (abs(diffInComp1) .gt. compChk) comp1Sign=1 
             if (abs(diffInComp1) .le. compChk) comp1Sign=0
             if (abs(diffInComp2) .gt. compChk) comp2Sign=1
             if (abs(diffInComp2) .le. compChk) comp2Sign=0
             
             
             write(346,*) j,grdV1(j,1),grdV2(j,1),diffInComp1,comp1Sign,&
                  grdV1(j,2),grdV2(j,2),diffInComp2,comp2Sign
             
          elseif (i==2) then
             
             diffInComp1 = (grdSprV1(j,1)-grdSprV2(j,1))
             diffInComp2 = (grdSprV1(j,2)-grdSprV2(j,2))
             
             if (abs(diffInComp1) .gt. compChk) comp1Sign=1 
             if (abs(diffInComp1) .le. compChk) comp1Sign=0
             if (abs(diffInComp2) .gt. compChk) comp2Sign=1
             if (abs(diffInComp2) .le. compChk) comp2Sign=0
             
             write(346,*) j,grdSprV1(j,1),grdSprV2(j,1),diffInComp1,comp1Sign,&
                  grdSprV1(j,2),grdSprV2(j,2),diffInComp2,comp2Sign
             
          elseif (i==3) then
             
             diffInComp1 = (grdAreaV1(j,1)-grdAreaV2(j,1))
             diffInComp2 = (grdAreaV1(j,2)-grdAreaV2(j,2))
             
             if (abs(diffInComp1) .gt. compChk) comp1Sign=1 
             if (abs(diffInComp1) .le. compChk) comp1Sign=0
             if (abs(diffInComp2) .gt. compChk) comp2Sign=1
             if (abs(diffInComp2) .le. compChk) comp2Sign=0
             
             write(346,*) j,grdAreaV1(j,1),grdAreaV2(j,1),diffInComp1,comp1Sign,&
                  grdAreaV1(j,2),grdAreaV2(j,2),diffInComp2,comp2Sign
             
          elseif (i==4) then
             
             diffInComp1 = (grdGrvtnlV1(j,1)-grdGrvtnlV2(j,1))
             diffInComp2 = (grdGrvtnlV1(j,2)-grdGrvtnlV2(j,2))
             
             if (abs(diffInComp1) .gt. compChk) comp1Sign=1 
             if (abs(diffInComp1) .le. compChk) comp1Sign=0
             if (abs(diffInComp2) .gt. compChk) comp2Sign=1
             if (abs(diffInComp2) .le. compChk) comp2Sign=0
             
             write(346,*) j,grdGrvtnlV1(j,1),grdGrvtnlV2(j,1),diffInComp1,comp1Sign,&
                  grdGrvtnlV1(j,2),grdGrvtnlV2(j,2),diffInComp2,comp2Sign
             
          elseif (i==5) then
             
             diffInComp1 = (grdBendV1(j,1)-grdBendV2(j,1))
             diffInComp2 = (grdBendV1(j,2)-grdBendV2(j,2))
             
             if (abs(diffInComp1) .gt. compChk) comp1Sign=1 
             if (abs(diffInComp1) .le. compChk) comp1Sign=0
             if (abs(diffInComp2) .gt. compChk) comp2Sign=1
             if (abs(diffInComp2) .le. compChk) comp2Sign=0
             
             write(346,*) j,grdBendV1(j,1),grdBendV2(j,1),diffInComp1,comp1Sign,&
                  grdBendV1(j,2),grdBendV2(j,2),diffInComp2,comp2Sign
             
          elseif (i==6) then
                
             diffInComp1 = (grdSRypV1(j,1)-grdSRypV2(j,1))
             diffInComp2 = (grdSRypV1(j,2)-grdSRypV2(j,2))
             
             if (abs(diffInComp1) .gt. compChk) comp1Sign=1 
             if (abs(diffInComp1) .le. compChk) comp1Sign=0
             if (abs(diffInComp2) .gt. compChk) comp2Sign=1
             if (abs(diffInComp2) .le. compChk) comp2Sign=0
             
             write(346,*) j,grdSRypV1(j,1),grdSRypV2(j,1),diffInComp1,comp1Sign,&
                  grdSRypV1(j,2),grdSRypV2(j,2),diffInComp2,comp2Sign
             
          endif
          
       enddo
       
       write(346,*) " "
       
    enddo
    
    close(346)
    
  end subroutine cmpreTwoGrdCompnts
  
end module gradient_module
