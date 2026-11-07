module grdPenF_from_total_system  !M1
  use system_parameters
  use gradient_from_total_system
  
  implicit none
  
  real*8, allocatable :: grdPenF(:,:)
  real*8, allocatable :: grdPenF_ana(:,:),grdPenF_num(:,:)
  real*8, allocatable :: grdPenF_mv(:) !1D array

contains

  subroutine allocate_and_initialize_grdPenF_system_variables
    implicit none
    
    allocate(grdPenF(1:N_node,1:N_dmnsn))    
    allocate(grdPenF_ana(1:N_node,1:N_dmnsn))
    allocate(grdPenF_num(1:N_node,1:N_dmnsn))
    
    !write(*,*) N_mvCoordnte,"N_mvCoordnte inside alloc"
    allocate(grdPenF_mv(1:N_mvCoordnte))
    
    
    grdPenF     = 10e5
    grdPenF_ana = 10e5
    grdPenF_num = 10e5
    grdPenF_mv  = 10e5  !very high value(unrealistic)
  
  end subroutine allocate_and_initialize_grdPenF_system_variables

  subroutine deallocate_grdPenF_system_variables
    implicit none

    deallocate(grdPenF)    
    deallocate(grdPenF_ana)
    deallocate(grdPenF_num)
    deallocate(grdPenF_mv)
    
  end subroutine deallocate_grdPenF_system_variables
  
  
end module grdPenF_from_total_system



module grdPenF_contrb_from_spr !M2
  use system_parameters
  use transfrm_info
  use triangle_routines
  
  implicit none
  
  real*8, allocatable :: grd_sprPf(:,:)
  real*8, allocatable :: delPf_delR(:)
  real*8, allocatable :: delR_delNodePf(:,:,:)
  
contains
  
  subroutine allocate_and_initialize_grdPenF_spr_variables
    implicit none
    
    allocate(grd_sprPf(1:N_node,1:N_dmnsn))
    allocate(delPf_delR(1:N_spr))
    allocate(delR_delNodePf(1:N_node,1:max_spr_node,1:N_dmnsn))
    
    grd_sprPf      = -1e30
    delPf_delR     = -1e30
    delR_delNodePf = -1e30
    
  end subroutine allocate_and_initialize_grdPenF_spr_variables

  subroutine deallocate_grdPenF_spr_variables
    implicit none
    
    deallocate(grd_sprPf)
    deallocate(delPf_delR)
    deallocate(delR_delNodePf)

  end subroutine deallocate_grdPenF_spr_variables

  
  subroutine grdPenF_components_frm_spr(dum_Nodes,dum_grdSprpf)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8  :: dum_grdSprpf(1:N_node,1:N_dmnsn)
    integer :: i
    integer :: spr_nmbr
    integer :: totl_spr_cnnctd

    integer :: given_node,other_node
    
    dum_grdSprpf = 0.0d0
    given_node = 0 ; other_node = 0
    
    call calc_delPf_delR_array(dum_Nodes)
    call calc_delR_delNodePf_array(dum_Nodes)     
    
    do i = 1,N_node

       if (node_typ(i) .eq. 0) then
          continue
       elseif (node_typ(i) .ne. 0) then
          totl_spr_cnnctd = node_spr(i,0)
          !write(*,*) totl_spr_cnnctd,"tsc
          
          if (node_cnnctd(i) .eq. 0) then
             
             if (node_typ(i) .eq. 1) then !free node
                call grdSprpf_lp(i)
             elseif (node_typ(i) .eq. 2) then !x free, y fixed
                call grdSprpf_lp(i)
             endif
             
          elseif (node_cnnctd(i) .ne. 0) then
             
             if (count_this_dn(i) .eq. 1) then
                
                if (node_typ(i) .eq. 1) then !free node
                   call grdSprpf_lp(i)
                elseif (node_typ(i) .eq. 2) then !x free, y fixed
                   call grdSprpf_lp(i)
                endif

                
             elseif (count_this_dn(i) .ne. 1) then
                !will have to replace with other double node value
                !has been replaced after calc
                continue
             endif
             
          endif
          
          !write(*,*) dum_grdSprpf(i,1:2),"dum_grdSprpf and i = ",i
       endif
       
       !write(*,*) dum_grdSprpf(i,1:2),"dum_grdSprpf and i = ",i
       
    enddo
    
    
  contains

    subroutine grdSprpf_lp(node_nm)
      implicit none
      integer :: j
      integer :: node_nm
      
      !if (node_nm .eq. 14) write(*,*) dum_grdSprpf(node_nm,1:2),"dum_grdSprpf14_bfr"

      do j = 1,totl_spr_cnnctd
         spr_nmbr = node_spr(node_nm,j)
         
         if (node_typ(node_nm) .eq. 1) then !x,y free
            dum_grdSprpf(node_nm,1:2) = dum_grdSprpf(node_nm,1:2) &
                 + delPf_delR(spr_nmbr) * delR_delNodePf(node_nm,j,1:2)
            
            !if (node_nm.eq.11.or.node_nm.eq.12)write(*,*)delPf_delR(spr_nmbr),delR_delNodePf(node_nm,j,1:2),"dPf/dR,dR/dNpf(1:2)"
            !if (node_nm .eq. 14) write(*,*) delPf_delR(spr_nmbr),delR_delNodePf(node_nm,j,2),"dPf/dR,dR/dNpf"
            !if (node_nm .eq. 14) write(*,*) dum_grdSprpf(node_nm,1:2),"dum_grdSprpf14_insd"
         
         elseif (node_typ(node_nm) .eq. 2) then !x free, y fixed        
            dum_grdSprpf(node_nm,1) = dum_grdSprpf(node_nm,1) &
                 + delPf_delR(spr_nmbr) * delR_delNodePf(node_nm,j,1)
            !if (node_nm .eq. 14) write(*,*) delPf_delR(spr_nmbr),delR_delNodePf(node_nm,j,1)
            !if (node_nm .eq. 14) write(*,*) dum_grdSprpf(node_nm,1:2),"dum_grdSprpf14_insd"
         endif
      enddo

      !if (node_nm .eq. 14) write(*,*) dum_grdSprpf(node_nm,1:2),"dum_grdSprpf14_aft"
    end subroutine grdSprpf_lp
    
  end subroutine grdPenF_components_frm_spr


  
  subroutine calc_delPf_delR_array(dum_Nodes)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    integer :: i
    
    do i = 1,N_spr
       delPf_delR(i) = indvd_delPf_delR(i)
       !write(*,*) delPf_delR(i),i,"delPf_delR"
    enddo
    
  contains
    
    real*8 function indvd_delPf_delR(spr_nmbr)
      implicit none
      integer :: spr_nmbr
      
      indvd_delPf_delR = alpha(spr_nmbr) * dflctn(spr_nmbr)
      
    end function indvd_delPf_delR
    
    
   !!!!dflctn and related subroutines are common for both Energy and gradient,find a suitable way to write once(!HERE)
    
    real*8 function dflctn(spr_nmbr)
      implicit none
      integer :: spr_nmbr
      integer :: Node1,Node2,Node3
      integer :: pulley_typ
      
      pulley_typ = 6
      
      !if (typ_spr(spr_nmbr) .ne. pulley_typ) then
      if (spr_node(spr_nmbr,0)==2) then

         Node1 = spr_node(spr_nmbr,1)
         Node2 = spr_node(spr_nmbr,2)
         
         l(spr_nmbr) = lngth_u_two_nodes(Node1,Node2)
         dflctn      = l(spr_nmbr) - Lt(spr_nmbr)

         !if (spr_nmbr==1.or.spr_nmbr==3.or.spr_nmbr==4) then
         !   write(*,*) dflctn,"dflctn"
         !endif
         
         if (abs(dflctn) .le. 1e-16) then
            dflctn = 0.0d0
         endif
         
      elseif (spr_node(spr_nmbr,0)==3) then
         
         Node1 = spr_node(spr_nmbr,1)
         Node2 = spr_node(spr_nmbr,2)
         Node3 = spr_node(spr_nmbr,3)
         
         l(spr_nmbr) = lngth_u_three_nodes(Node1,Node2,Node3)
         dflctn = l(spr_nmbr) - Lt(spr_nmbr)
         
         if (abs(dflctn) .le. 1e-15) then
            dflctn = 0.0d0
         endif
         
         !if (spr_nmbr.eq.13) write(*,*) dflctn,"dflctn"
         
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
         
         lngth_u_three_nodes = lngth_u_three_nodes + sqrt(res_sq)
      enddo
      
    end function lngth_u_three_nodes
    
    
!!!! dflctn related sbrtns are finished here  
      
  end subroutine calc_delPf_delR_array
  
  subroutine calc_delR_delNodePf_array(dum_Nodes)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    
    integer :: i,j
    integer :: totl_spr_cnnctd
    integer :: spr_nmbr,prev_spr_nmbr
    integer :: init_node
    integer :: Pulley_typ

    Pulley_typ = 6
    prev_spr_nmbr = 0
    
    do i = 1,N_node
       !write(*,*) i,node_spr(i,1:max_spr_node),node_spr(i,0),"node_spr for node i"
       
       init_node = i
       totl_spr_cnnctd = node_spr(i,0)
       
       !if(i.eq.13) write(*,*) totl_spr_cnnctd,"tsc"
       
       do j = 1,totl_spr_cnnctd
          spr_nmbr = node_spr(i,j)

          if (j.eq.1) then
             prev_spr_nmbr = 0
          elseif (j.ne.1) then
             prev_spr_nmbr = node_spr(i,j-1)  
          endif
          
          !if(i.eq.13) write(*,*) spr_nmbr,"spr_nmbr"
          
          !if (typ_spr(spr_nmbr) .ne. Pulley_typ) then
          if (spr_node(spr_nmbr,0)==2) then
             
             delR_delNodePf(i,j,1:2) = indvd_delR_delNodePf_NonPulley(spr_nmbr,init_node)

          !elseif (typ_spr(spr_nmbr) .eq. Pulley_typ) then
          elseif (spr_node(spr_nmbr,0)==3) then
             delR_delNodePf(i,j,1:2) = indvd_delR_delNodePf_Pulley(spr_nmbr,prev_spr_nmbr,init_node)

          endif
          
          !delR_delNodePf(i,j,1:2) = 0.0d0
       enddo
       
    enddo
    
    
  contains
    
    function indvd_delR_delNodePf_NonPulley(spr_nmbr,init_node)
      implicit none
      real*8  :: indvd_delR_delNodePf_NonPulley(1:N_dmnsn)
      integer :: spr_nmbr
      integer :: finl_node,init_node !node_indexes
      integer :: i
      logical :: lgcl_dn,lgcl_sameDn
      
      indvd_delR_delNodePf_NonPulley = 0.0d0
      lgcl_dn = .False.
      lgcl_sameDn = .False.
      
      !if (typ_spr(spr_nmbr) .ne. pulley_typ) then
      
      do i = 1,max_node_spr
         if ((spr_node(spr_nmbr,i).ne.init_node) .and. (spr_node(spr_nmbr,i).ne.(-1))) then
            finl_node = spr_node(spr_nmbr,i)
            
            !if (init_node.eq.17) then
            !   write(*,*) finl_node,"finl_nodeIN"
            !endif
            
            !call is_it_a_dn(finl_node,lgcl_dn)
            !write(*,*) lgcl_dn,"lgcl_dn"
            call are_they_same_dn(lgcl_sameDn,init_node,finl_node)
            !write(*,*) lgcl_sameDn,"lgcl_sameDn"
            
            if (lgcl_sameDn .eqv. .True.) then
               finl_node = 0
               cycle
            endif
            
            !write(*,*) init_node,finl_node,"init_node,finl_node"
            !write(*,*) lgcl_dn,"lgcl_dn"
            
            indvd_delR_delNodePf_NonPulley(1:N_dmnsn) = get_unit_vector(finl_node,init_node)
            
         endif
      enddo
     
    end function indvd_delR_delNodePf_NonPulley


    function indvd_delR_delNodePf_Pulley(spr_nmbr,prev_spr_nmbr,init_node)
      implicit none
      
      real*8  :: indvd_delR_delNodePf_Pulley(1:N_dmnsn)
      integer :: spr_nmbr,prev_spr_nmbr
      integer :: finl_node,init_node !node_indexes
      logical :: lgcl_trminl
      
      indvd_delR_delNodePf_Pulley = 0.0d0
      
      !if (init_node.eq.17) write(*,*) init_node,spr_nmbr,"init_node,spr_nmbr"
      call is_it_a_trminl_node(init_node,spr_nmbr,lgcl_trminl)
      
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

      indvd_delR_delNodePf_Pulley(1:N_dmnsn) = get_unit_vector(finl_node,init_node)

    end function indvd_delR_delNodePf_Pulley

      
    
    
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
         write(*,*) "Check vector_absValue_grdcalcPenF.dat"
         open(unit=128,file='vector_absValue_grdcalcPenF.dat')
         
         write(128,*) init_node,finl_node,"init and finl node"
         write(128,*) "less than precision value, flnm:gradient_calcEnrgy,sbrtn: get_unit_vector"
         write(128,*) r12(1),r12(2),r12_abs,"r12(1:2)"
         write(128,*) " "
         
         close(128)
         stop
      endif
      !write(*,*) r12,r12_abs,"in get_unit_vector"

      r12_uv  = r12/r12_abs
     
      get_unit_vector = r12_uv
      
    end function get_unit_vector
      
  end subroutine calc_delR_delNodePf_array
       
end module grdPenF_contrb_from_spr
  
  
module grdPenF_contrb_from_area  !M3
  use system_parameters !cell_info + non_changing_parameters
  use transfrm_info
  use triangle_routines
  
  implicit none
  
  real*8, allocatable :: grd_areaPf(:,:)
  real*8, allocatable :: delPf_delA(:)
  real*8, allocatable :: delA_delNodePf(:,:,:) !3D array
  
contains
  
  subroutine allocate_and_initialize_grdPenF_area_variables
    
    implicit none

    allocate(grd_areaPf(1:N_node,1:N_dmnsn))
    allocate(delPf_delA(1:N_cell))
    allocate(delA_delNodePf(1:N_node,1:max_area_node,1:N_dmnsn))

    grd_areaPf     = -1e30
    delPf_delA     = -1e30 
    delA_delNodePf = -1e30 
    
  end subroutine allocate_and_initialize_grdPenF_area_variables

  subroutine deallocate_grdPenF_area_variables
    implicit none

    deallocate(grd_areaPf)
    deallocate(delPf_delA)
    deallocate(delA_delNodePf)

  end subroutine deallocate_grdPenF_area_variables
  
  subroutine grdPenF_components_frm_area(dum_Nodes,dum_grdAreapf)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8  :: dum_grdAreapf(1:N_node,1:N_dmnsn)
    integer :: i
    integer :: area_nmbr
    integer :: totl_cell_cnnctd 

    dum_grdAreapf = 0.0d0
    
    call calc_delPf_delA_array(dum_Nodes)
    call calc_delA_delNodePf_array(dum_Nodes)
    
    do i = 1,N_node
       
       if (node_typ(i) .eq. 0) then
          continue
       elseif (node_typ(i) .ne. 0) then
          totl_cell_cnnctd = node_area(i,0)
          !write(*,*) totl_cell_cnnctd,"tcc"

          if (node_cnnctd(i) .eq. 0) then
             
             if (node_typ(i) .eq. 1) then !free node
                call grdAreapf_lp(i)
             elseif (node_typ(i) .eq. 2) then !x free, y fixed
                call grdAreapf_lp(i)
             endif
             
          elseif (node_cnnctd(i) .ne. 0) then
             
             if (count_this_dn(i) .eq. 1) then
                
                if (node_typ(i) .eq. 1) then !free node
                   call grdAreapf_lp(i)
                elseif (node_typ(i) .eq. 2) then !x free, y fixed
                   call grdAreapf_lp(i)
                endif

             elseif (count_this_dn(i) .ne. 1) then
                !will have to replace with other double node value

                continue
             endif
             
          endif
          
       endif
       
       !write(*,*) dum_grdAreapf(i,1:2),"dum_grdAreapf and i = ",i
       
    enddo

    
    contains

    subroutine grdAreapf_lp(node_nm)
      implicit none
      integer :: j
      integer :: node_nm
      
      do j = 1,totl_cell_cnnctd
         area_nmbr = node_area(node_nm,j)
         
         if (node_typ(node_nm) .eq. 1) then  !x,y free
            dum_grdAreapf(node_nm,1:2) = dum_grdAreapf(node_nm,1:2) &
                 + delPf_delA(area_nmbr) * delA_delNodePf(node_nm,j,1:2)

            !if (node_nm.eq.9) write(*,*) delPf_delA(area_nmbr),delA_delNodePf(node_nm,j,1:2),"dPf/dA,dA/dNpf(1:2)",node_nm
            !if (node_nm.eq.10)write(*,*) delPf_delA(area_nmbr),delA_delNodePf(node_nm,j,1:2),"dPf/dA,dA/dNpf(1:2)",node_nm
            !if (node_nm.eq.11)write(*,*) delPf_delA(area_nmbr),delA_delNodePf(node_nm,j,1:2),"dPf/dA,dA/dNpf(1:2)",node_nm
            !if (node_nm.eq.12)write(*,*) delPf_delA(area_nmbr),delA_delNodePf(node_nm,j,1:2),"dPf/dA,dA/dNpf(1:2)",node_nm
            !if (node_nm.eq.19)write(*,*) delPf_delA(area_nmbr),delA_delNodePf(node_nm,j,1:2),"dPf/dA,dA/dNpf(1:2)",node_nm
            
         elseif (node_typ(node_nm) .eq. 2) then !x free, y fixed        
            dum_grdAreapf(node_nm,1) = dum_grdAreapf(node_nm,1) &
                 + delPf_delA(area_nmbr) * delA_delNodePf(node_nm,j,1)
         endif
      enddo
      
    end subroutine grdAreapf_lp
    
  end subroutine grdPenF_components_frm_area

  
  subroutine calc_delPf_delA_array(dum_Nodes)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    integer :: i
    
    do i = 1,N_cell
       delPf_delA(i) = indvd_delPf_delA(i)
       !write(*,*) delPf_delA(i),i,"DeLE_DeLA_inside"
    enddo


  contains

    real*8 function indvd_delPf_delA(area_nmbr)
      implicit none
      integer :: area_nmbr
      real*8  :: ka,area_ch

      ka = beta(area_nmbr)
      area_ch = area_chng(area_nmbr)

      !write(*,*) ka,area_ch,"Beta and Area Change"
      !indvd_delPf_delA = beta(area_nmbr) * area_chng(area_nmbr)
      indvd_delPf_delA = ka*area_ch
      
      !HERE A(area_nmbr) and At(area_nmbr) needs to be accesibile
      
    end function indvd_delPf_delA


    real*8 function area_chng(area_nmbr)
      implicit none
      integer :: area_nmbr
      
      !A(area_nmbr) = A_currnt(area_nmbr)
      area_chng    = A_currnt(area_nmbr) - At(area_nmbr)
      
      !write(*,*) A_currnt(area_nmbr),At(area_nmbr),"A_curr_At"
      !write(*,*) area_chng,"area_chng"
      !write(*,*) ""
      !if (area_nmbr.eq.4) write(*,*) A_currnt(area_nmbr),At(area_nmbr),beta(area_nmbr),"A_curr_At_beta"
      
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

    
  end subroutine calc_delPf_delA_array

  

  subroutine calc_delA_delNodePf_array(dum_Nodes)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    integer :: i,j
    integer :: totl_area_cnnctd
    integer :: area_nmbr,node_nmbr

    
    do i = 1,N_node
       node_nmbr = i
       totl_area_cnnctd = node_area(i,0)
       
       do j = 1,totl_area_cnnctd          
          area_nmbr = node_area(i,j)
          delA_delNodePf(i,j,1:N_dmnsn) = indvd_delA_delNodePf(area_nmbr,node_nmbr)   
       enddo
    enddo
    
  !end subroutine calc_delA_delNodePf_array

 contains

   function indvd_delA_delNodePf(area_nmbr,node_nmbr)
     implicit none
     real*8  :: indvd_delA_delNodePf(1:N_dmnsn)
     integer :: area_nmbr,node_nmbr
     integer :: neighbour_nodes(1:2)
     !integer :: max_i
     integer :: Neigh1,Neigh2
     
     real*8  :: first_xy(1:3),second_xy(1:3),third_xy(1:3)
     real*8  :: triangle_area
     !max_i = area_node(area_nmbr,0)
     
     indvd_delA_delNodePf = 0.0d0
     
     !write(*,*) dum_Nodes,"Dummy Nodes inside Indvd delA delNodePf"
     call get_neighbour_nodes(area_nmbr,node_nmbr,neighbour_nodes)
     
     !write(*,*) neighbour_nodes(1:2),"neigh_nodes"
     Neigh1 = neighbour_nodes(1)
     Neigh2 = neighbour_nodes(2)
     !write(*,*) Neigh1,Neigh2,"Neigh1,Neigh2 for (Nn,An) = (",node_nmbr,area_nmbr,")"
     
     first_xy(1:2)  = dum_Nodes(node_nmbr,1:2)
     second_xy(1:2) = dum_Nodes(Neigh1,1:2)
     third_xy(1:2)  = dum_Nodes(Neigh2,1:2)

     first_xy(3)=0.0d0 ; second_xy(3)=0.0d0 ; third_xy(3)=0.0d0
    
     triangle_area = calc_triangle(first_xy,second_xy,third_xy)
    
    
     call check_triangle_area_value(triangle_area)
    
     indvd_delA_delNodePf(1) = 0.5 * (dum_Nodes(Neigh1,2) - dum_Nodes(Neigh2,2)) ! (y_Neigh1 - y_Neigh2)
     indvd_delA_delNodePf(2) = 0.5 * (dum_Nodes(Neigh2,1) - dum_Nodes(Neigh1,1)) ! (x_Neigh2 - x_Neigh1)
     
     !write(*,*) indvd_delA_delNodePf(1:2),"indvd_delA_delNodePf"
     
   end function indvd_delA_delNodePf

 end subroutine calc_delA_delNodePf_array


 
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
         write(*,*) "nodes cant be negative, flnm: grdPenF_calc, sbrtn: get_neighbour_nodes"
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
            !   neighbour_nodes(1) = area_node(area_nmbr,(node_pstn_in_array+1))
             !else
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
 
end module grdPenF_contrb_from_area





module grdPenF_module  !M4
  use system_parameters
  use grdPenF_from_total_system
  use grdPenF_contrb_from_spr
  use grdPenF_contrb_from_area

  !use l0_A0_variatns !not needed for grdPenF
  use conversion_routines
  use PenF_module
  
  implicit none  
  
contains

  
  subroutine get_all_grdPenF_variables
    implicit none

    call allocate_and_initialize_grdPenF_system_variables
    call allocate_and_initialize_grdPenF_spr_variables
    call allocate_and_initialize_grdPenF_area_variables
    
  end subroutine get_all_grdPenF_variables

  subroutine deallocate_all_grdPenF_variables
    implicit none

    call deallocate_grdPenF_system_variables
    call deallocate_grdPenF_spr_variables
    call deallocate_grdPenF_area_variables
    
  end subroutine deallocate_all_grdPenF_variables

  
  subroutine get_grdPenF(dum_Nodes,dum_grdN)
    implicit none
    real*8 :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8 :: dum_grdN(1:N_node,1:N_dmnsn)

    !write(*,*) dum_Nodes,"node_xy should match"
    !write(*,*) dum_grdN, "grdPenF"
    
    if (Analtcl_or_Numrcl .eq. 1) then
       call grdPenF_Analtcl(dum_Nodes,dum_grdN)
           
    elseif (Analtcl_or_Numrcl .eq. 2) then
       call grdPenF_Numrcl(dum_Nodes,dum_grdN)
       
    endif
    
  end subroutine get_grdPenF

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

       !write(*,*) dum_grd(node_toBe_used,1:2),"toBe_used grdPenF"
       !write(*,*) dum_grd(node_toBe_filled,1:2),"toBe_filled_grdPenF_bfr"
       
       dum_grd(node_toBe_filled,1:2) = dum_grd(node_toBe_used,1:2)

       !write(*,*) dum_grd(node_toBe_filled,1:2),"toBe_filled_grdPenF_aft"
       
    enddo
    

  end subroutine complte_the_grdnts
  
  subroutine grdPenF_Analtcl(dum_Nodes,dum_grdN)
    implicit none
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8  :: dum_grdN(1:N_node,1:N_dmnsn)
    real*8  :: dum_grdSprpf(1:N_node,1:N_dmnsn)
    real*8  :: dum_grdAreapf(1:N_node,1:N_dmnsn)
   
    integer :: i,j   

    integer :: lam_spr,lam_area
        
    dum_grdN      = 0.0d0
    dum_grdSprpf  = 0.0d0
    dum_grdAreapf = 0.0d0
  
    lam_spr = 0.0d0 ; lam_area = 0.0d0
    
    if (spr_lgcl.eqv..True.) then
       lam_spr    = 1.0d0
       call grdPenF_components_frm_spr(dum_Nodes,dum_grdSprpf)
       call complte_the_grdnts(dum_grdSprpf)
       !write(*,*) dum_grdSprpf,"Dum_GRD_SPRPF"

       do i = 1,N_node
          !write(*,*) dum_grdSprpf(i,1:2),"dum_grdSprpf for Node =",i
       enddo
       
    endif
    
    if (area_lgcl.eqv..True.) then
       lam_area   = 1.0d0
       call grdPenF_components_frm_area(dum_Nodes,dum_grdAreapf)
       call complte_the_grdnts(dum_grdAreapf)
       !write(*,*) dum_grdAreapf,"Dum_GRD_AREA"

       do i = 1,N_node
          !write(*,*) dum_grdAreapf(i,1:2),"dum_grdAreapf for Node =",i
       enddo
       
    endif
    
    
    do i = 1,N_node
       do j = 1,N_dmnsn
          dum_grdN(i,j) = lam_spr*dum_grdSprpf(i,j) + lam_area*dum_grdAreapf(i,j)
       enddo
    enddo
    
    call get_the_other_dn_grd_values
    !write(*,*) dum_grdN,"dum_grdN_INSIDE"
    call print_analytical_grdPenFs
    
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

    subroutine print_analytical_grdPenFs
      implicit none
      integer :: i1
      
      open(unit=25,file="grd_sprPf_grd_areaPf.dat")
      open(unit=26,file='grdAnaPf.dat')
      open(unit=27,file='grdSprAnaPf.dat')
      open(unit=28,file='grdAreaAnaPf.dat')
      
      do i1 = 1,N_node
         write(25,*) dum_grdSprpf(i1,1:2),dum_grdAreaPf(i1,1:2)
         write(25,*) dum_grdN(i1,1:2),"dum_grdN=dum_grdSprpf+dum_grdAreapf"
         write(26,*) dum_grdN(i1,1:2),i1,"grd_AnaPf,node" 
         write(27,*) dum_grdSprpf(i1,1:2),"grd_AnaSprPf,node"
         write(28,*) dum_grdAreaPf(i1,1:2),"grd_AnaAreaPf,node"
      enddo
      
      close(25)
      close(26)
      close(27)
      close(28)
      
    end subroutine print_analytical_grdPenFs
    
  end subroutine grdPenF_Analtcl


  
  subroutine grdPenF_Numrcl(dum_Nodes,dum_grdN)
    implicit none

    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8  :: dum_grdN(1:N_node,1:N_dmnsn)
    real*8  :: str,PenF_aft,PenF_bfr
    real*8  :: h=1e-4
    
    integer :: i,j
    
    dum_grdN(1:N_node,1:N_dmnsn) = 0.d0
    !write(*,*) N_mvCoordnte,"TMP"
    
    do i = 1,N_node
       
       do j = 1,N_dmnsn
          str            = dum_Nodes(i,j)
          !if(i.eq.7) write(*,*) str,"str"
          
          dum_Nodes(i,j) = dum_Nodes(i,j) + h
          !write(*,*) dum_Nodes(i,j)
          
          PenF_aft       = PenF(dum_Nodes)
          dum_Nodes(i,j) = str
          
          dum_Nodes(i,j) = dum_Nodes(i,j) - h
          !write(*,*) dum_Nodes(i,j)
          
          PenF_bfr       = PenF(dum_Nodes)
          
          dum_Nodes(i,j) = str
          
          dum_grdN(i,j) = (PenF_aft-PenF_bfr)/(2.d0*h)
          
          !write(*,*) PenF_aft,PenF_bfr,"PenFaft,PenFbfr"
          !write(*,*) dum_grdN(i,j),"dum_grdN(i,j)"
          !if(i.eq.7) write(*,*) PenF_aft,PenF_bfr,"PenFaft,PenFbfr"
          
       enddo
        
    enddo

    write(*,*) " adjust this as  in gradient_calcEnergy"
    stop
    call print_numerical_grdPenFs
    
  contains
    
    subroutine print_numerical_grdPenFs
      implicit none
      integer :: i
      
      open(unit=132,file='NumgrdPenF.dat')
      
      do i = 1,N_node
         write(132,fmt=*) dum_grdN(i,1:2),"grd_Numerical()"
      enddo
      
      close(132)
      
    end subroutine print_numerical_grdPenFs
    
  end subroutine grdPenF_Numrcl
  


  
  !VECTOR TO COORDINATE TRANSFRM STUFFS FOR GRADIENT ARE BELOW
  
  subroutine grdPenF_coordntes(dum_coordntes,dum_grdC)

    implicit none
    real*8,intent(in) :: dum_coordntes(1:N_mvCoordnte)
    real*8,intent(out) :: dum_grdC(1:N_mvCoordnte)!C=coordnte_based
    real*8 :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8 :: dum_grdN(1:N_node,1:N_dmnsn)!N=node_based
    
    call coordntes_to_nodes(dum_coordntes,dum_Nodes)
    !write(*,*) dum_coordntes,"dum_coordntes inside_grdPenF_coordntes_sbrtn" !CO1
    !write(*,*) dum_Nodes,"dum_Nodes inside_grdPenF_coordntes_sbrtn"!CO2

    call get_grdPenF(dum_Nodes,dum_grdN)
    !write(*,*) dum_grdN,"dum_grdN"!CO3
    
    call reduced_to_grd_coordntes_xy(dum_grdN,dum_grdC)
    !write(*,*) dum_grdC,"dum_grdC_out"!CO4
    
    call print_coordnte_based_grdPenFs(dum_grdC)
    
  contains
    
    subroutine reduced_to_grd_coordntes_xy(dum_grdN,dum_grdC)
      implicit none
      real*8  :: dum_grdN(1:N_node,1:N_dmnsn)
      real*8  :: dum_grdC(1:N_mvCoordnte)
      integer :: i,count_mv
      
      count_mv = 1
      !write(*,*) N_mvCoordnte,"N_mvCoordnte"
      
      dum_grdC(1:N_mvCoordnte) = -1e30
      
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

    subroutine print_coordnte_based_grdPenFs(dum_grdC)
      implicit none
      real*8  :: dum_grdC(1:N_mvCoordnte)
      integer :: i

      open(unit=26,file='co_ordnte_based_grdPenF.dat')
      
      do i = 1,N_mvCoordnte
         write(unit=26,fmt=*) dum_grdC(i), "coordnte_based_grdPenFs"
      enddo

      close(26)
      
    end subroutine print_coordnte_based_grdPenFs
    
  end subroutine grdPenF_coordntes


  subroutine check_values_of_grdPenF !wont be used now
    implicit none
    
    integer :: fnode,spr_no
    integer :: i
    real*8  :: dx(1:3),dy(1:3),r(1:3)
    real*8  :: grd_prt_spr(1:3,1:2)
    
    integer :: dnode1=0,dnode2=0,area_nm
    real*8  :: dx_a(1:2),dy_a(1:2)
    real*8  :: grd_prt_area(1:2,1:2) = 0.0d0
    
    real*8  :: grd_sprPfchk(1:2),grd_areaPfchk(1:2)
    real*8  :: grdchk(1:2)
    
    grd_sprPfchk(1:2)  = 0.0d0
    grd_areaPfchk(1:2) = 0.0d0
    grdchk(1:2)        = 0.0d0
    
    open(unit=22, file='checking_grdPenF_value.dat')
    write(unit=22,fmt=*) "We will check for Node 3"
    
    write(unit=22, fmt=*) node_spr(3,1:max_spr_node),"node_to_spr"
    write(unit=22, fmt=*) alpha(1),alpha(3),alpha(4),"alpha"
    write(unit=22, fmt=*) Lt(1),Lt(3),Lt(4),"Lt"
    
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
       write(unit=22, fmt=*) "dis to fnode  =",r(i)

       grd_prt_spr(i,1) = alpha(spr_no) * (r(i)-Lt(spr_no)) * (dx(i) / r(i))
       grd_prt_spr(i,2) = alpha(spr_no) * (r(i)-Lt(spr_no)) * (dy(i) / r(i))

       write(unit=22,fmt=*) grd_prt_spr(i,1:2),"grd_prt_spr for i =",i
       write(unit=22,fmt=*) " "
       write(unit=22,fmt=*) " "

       grd_sprPfchk(1) = grd_sprPfchk(1) + grd_prt_spr(i,1)
       grd_sprPfchk(2) = grd_sprPfchk(2) + grd_prt_spr(i,2)
       
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

       
       !grd_prt_area(i,1) = beta(area_nm) * area_chng(area_nm) * (dy_a(i))
       !grd_prt_area(i,2) = beta(area_nm) * area_chng(area_nm) * (dx_a(i))

       
       !area_chng can't be accessed HERE as it is in the contains part of a subroutine

       
       !write(unit=22,fmt=*) area_chng(area_nm),"area_change for i=",i

       write(unit=22,fmt=*) dy_a(i),dx_a(i),"dy_a,dx_a"
       write(unit=22,fmt=*) grd_prt_area(i,1:2),"grd_prt_area for i =",i
       write(unit=22,fmt=*) " "
       write(unit=22,fmt=*) " "

    enddo

    grdchk(1) = grd_sprPfchk(1) + grd_areaPfchk(1)
    grdchk(2) = grd_sprPfchk(2) + grd_areaPfchk(2)
    
    write(unit=22,fmt=*) grdchk(1:2),"grdchk"
    
    close(22)
    
  end subroutine check_values_of_grdPenF

  subroutine turn_on_one_strength_property(spr_or_area,spr_or_area_nmbr) !wont be used now
    implicit none
    integer :: spr_or_area
    integer :: spr_or_area_nmbr
    
    if (spr_or_area .eq. 1) then       
       !alpha(1:N_spr)   = 0.0d0 !you have to store first
       !beta(1:N_cell) = 0.0d0
       !alpha(spr_or_area_nmbr) = 1.0d0

    elseif (spr_or_area .eq. 2) then
       !alpha(1:N_spr)   = 0.0d0
       !beta(1:N_cell) = 0.0d0
       !beta(spr_or_area_nmbr) = 1.0d0

    endif

  end subroutine turn_on_one_strength_property
  
end module grdPenF_module


