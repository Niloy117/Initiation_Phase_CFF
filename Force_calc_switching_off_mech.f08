module node_info_tmp
  implicit none
  
  real*8, allocatable :: node_xy_tmp(:,:)
  integer,allocatable :: node_typ_tmp(:),node_cnnctd_tmp(:)
  integer,allocatable :: double_node_tmp(:,:)
  integer,allocatable :: count_this_dn_tmp(:)
  
end module node_info_tmp

module spr_info_tmp
  implicit none
  
  integer,allocatable :: typ_spr_tmp(:)
  real*8, allocatable :: k_spr_tmp(:),l0_tmp(:)
  real*8, allocatable :: l_tmp(:)

end module spr_info_tmp

module area_info_tmp
  implicit none
      
  integer,allocatable :: typ_area_tmp(:)
  real*8, allocatable :: k_area_tmp(:), A0_tmp(:)
  real*8, allocatable :: A_tmp(:)
  
end module area_info_tmp

module grv_info_tmp
  implicit none
  
  real*8, allocatable :: CgXNode_tmp(:)
  
end module grv_info_tmp

module bend_info_tmp
  implicit none

  real*8, allocatable :: k_phi_tmp(:,:)
  
end module bend_info_tmp

module cell_info_tmp
  
  use node_info_tmp
  use spr_info_tmp
  use area_info_tmp
  use grv_info_tmp
  use bend_info_tmp
  use system_parameters
  use neighbour_info_and_trnsfrmInfo_dependent_info
  
contains
  
  subroutine alloc_and_init_tmp_cell_variables
    implicit none
    
    allocate(node_xy_tmp(1:N_node,1:N_dmnsn))
    allocate(node_typ_tmp(1:N_node),node_cnnctd_tmp(1:N_node))
    allocate(double_node_tmp(1:N_node,1:2))
    allocate(count_this_dn_tmp(1:N_node))
    
    node_xy_tmp       = -1e20
    node_typ_tmp      = -1
    node_cnnctd_tmp   =  0
    double_node_tmp   = -1
    count_this_dn_tmp = -1
    
    allocate(typ_spr_tmp(1:N_spr),k_spr_tmp(1:N_spr))
    allocate(l0_tmp(1:N_spr),l_tmp(1:N_spr))
    
    typ_spr_tmp  = -1
    k_spr_tmp    = -1e20
    l0_tmp       = -1e25
    l_tmp        = -1e25
    
    allocate(k_area_tmp(1:N_cell),A0_tmp(1:N_cell))
    allocate(A_tmp(1:N_cell),typ_area_tmp(1:N_cell))
    
    typ_area_tmp  = -1
    k_area_tmp    = -1e20
    A0_tmp        = -1e25
    A_tmp         = -1e25
    
    allocate(CgXNode_tmp(1:N_node))
    
    CgXNode_tmp = -1.d20
    
    allocate(k_phi_tmp(1:N_node,1:max_Phi_node))

    k_phi_tmp = -1.d20
    
  end subroutine alloc_and_init_tmp_cell_variables
  
  subroutine store_coeffs_to_tmp
    implicit none
    
    k_spr_tmp(1:N_spr)    = k_spr(1:N_spr)
    k_area_tmp(1:N_cell)  = k_area(1:N_cell)
    CgXNode_tmp(1:N_node) = CgXNode(1:N_node)
    
    k_phi_tmp(1:N_node,1:max_Phi_node) = k_phi(1:N_node,1:max_Phi_node)
    
  end subroutine store_coeffs_to_tmp
  
  subroutine set_all_coeffs_to_zero_value
    implicit none
    
    k_spr(1:N_spr)    = 0.0d0
    k_area(1:N_cell)  = 0.0d0
    CgXNode(1:N_node) = 0.0d0
    
    k_phi(1:N_node,1:max_Phi_node) = 0.0d0
    
  end subroutine set_all_coeffs_to_zero_value
  
  subroutine set_value_to_one_entity(entity,entity_nmbr)
    implicit none
    integer :: entity
    integer :: entity_nmbr
    
    integer :: node_nm,area_nm,area_serial
    
    if (entity==1) then
       k_spr(entity_nmbr)  = k_spr_tmp(entity_nmbr)
    elseif (entity==2) then
       k_area(entity_nmbr) = k_area_tmp(entity_nmbr)
    elseif (entity==3) then
       CgXNode(entity_nmbr) = CgXNode_tmp(entity_nmbr)
    elseif (entity==4) then
       node_nm = Phi_info(entity_nmbr,1)
       area_nm = Phi_info(entity_nmbr,2)
       call get_area_serial(node_nm,area_nm,area_serial)
       k_phi(node_nm,area_serial) = k_phi_tmp(node_nm,area_serial) 
    else
       write(*,*) "Entitty number must be 1,2,3 or 4 flnm:Force_calc_switching_off_mech, sbrtn_nm:set_value_to_one_entity"
       stop
    endif
    
  end subroutine set_value_to_one_entity
  
  subroutine restore_coeffs_frm_tmp
    implicit none
    
    k_spr(1:N_spr)    = k_spr_tmp(1:N_spr)
    k_area(1:N_cell)  = k_area_tmp(1:N_cell)
    CgXNode(1:N_node) = CgXNode_tmp(1:N_node)
    
    k_phi(1:N_node,1:max_Phi_node) = k_phi_tmp(1:N_node,1:max_Phi_node)
    
  end subroutine restore_coeffs_frm_tmp
  
  subroutine dealloc_tmp_cell_variables
    implicit none
    
    deallocate(node_xy_tmp)
    deallocate(node_typ_tmp,node_cnnctd_tmp)
    deallocate(double_node_tmp)
    deallocate(count_this_dn_tmp)
    
    deallocate(typ_spr_tmp,k_spr_tmp)
    deallocate(l0_tmp,l_tmp)
    
    deallocate(k_area_tmp,A0_tmp)
    deallocate(A_tmp,typ_area_tmp)
    
    deallocate(CgXNode_tmp)
    deallocate(k_phi_tmp)
    
  end subroutine dealloc_tmp_cell_variables
  
end module cell_info_tmp

module force_calc
  use gradient_module
  use cell_info_tmp
  
  implicit none
  real*8, allocatable :: Force_node(:,:)
  real*8, allocatable :: Force_spr(:,:)
  real*8, allocatable :: Force_area(:,:)
  real*8, allocatable :: Force_grvty(:,:)
  real*8, allocatable :: Force_bend(:,:)
  
contains


  subroutine alloc_and_init_force_variables
    implicit none

    allocate(Force_node(1:N_node,1:N_dmnsn))
    allocate(Force_spr(1:N_spr,1:N_dmnsn))
    allocate(Force_area(1:N_cell,1:N_dmnsn))
    allocate(Force_grvty(1:N_node,1:N_dmnsn))
    allocate(Force_bend(1:N_phi,1:N_dmnsn))
    
    Force_node  = -1e30 !(unrealistic_value)
    Force_spr   = -1e30
    Force_area  = -1e30
    Force_grvty = -1e30
    Force_bend  = -1e30
    
  end subroutine alloc_and_init_force_variables
  
  
  subroutine get_Springwise_Areawise_and_Nodewise_forces
    implicit none
    
    call alloc_and_init_force_variables
    
    call get_force_in_all_springs
    call get_force_in_all_areas
    call get_force_in_all_nodes_dueto_grvty
    call get_force_in_all_nodes_dueto_bend
    
    call get_forces_at_nodes
    
  end subroutine get_Springwise_Areawise_and_Nodewise_forces


  
  subroutine get_force_in_all_springs
    implicit none
    integer :: i

    open(unit=307,file='Force_in_spr.dat')
    
    do i = 1,N_spr
       call get_force_in_a_spr(i)
       write(unit=307,fmt=*) Force_spr(i,1:N_dmnsn),i,"spr_nmbr"
    enddo

    close(307)
    
  end subroutine get_force_in_all_springs


  
  subroutine get_force_in_all_areas
    implicit none
    integer :: i
    
    open(unit=307,file='Force_in_area.dat')

    do i = 1,N_cell
       call get_force_in_an_area(i)
       write(unit=307,fmt=*) Force_area(i,1:N_dmnsn),i,"area_nmbr"
    enddo
    
    close(307)
    
  end subroutine get_force_in_all_areas
  
  
  subroutine get_force_in_all_nodes_dueto_grvty
    implicit none
    integer :: i
    
    open(unit=307,file='Force_in_grvty.dat')
    
    do i = 1,N_node
       call get_force_dueto_Cg_atNode(i)
       write(unit=307,fmt=*) Force_grvty(i,1:N_dmnsn),i,"node_nmbr"
    enddo
    
    close(307)
    
  end subroutine get_force_in_all_nodes_dueto_grvty
  
  subroutine get_force_in_all_nodes_dueto_bend
    implicit none
    integer :: i

    open(unit=308,file='Force_in_bend.dat')
    
    do i = 1,N_phi
       call get_force_dueto_bend_atPhi(i)
       write(unit=308,fmt=*) Force_bend(i,1:N_dmnsn),i,"node_nmbr"
    enddo
    
    close(308)
    
  end subroutine get_force_in_all_nodes_dueto_bend
  
  
  
  subroutine get_forces_from_gradient_at_nodes(grd_value)
    
    implicit none
    real*8  :: grd_value(1:N_node,1:N_dmnsn)
    integer :: i

    !open(unit=307,file="Force_in_node.dat")

    do i = 1,N_node
       Force_node(i,1:2) = - grd_value(i,1:2)
       !write(unit=307,fmt=*) Force_node(i,1:N_dmnsn),i,"node_nm"
    enddo

    !close(307)
    
  end subroutine get_forces_from_gradient_at_nodes

  
  subroutine get_forces_at_nodes
    
    implicit none
    integer :: i

    call get_gradient(node_xy,l0,A0,gradient)
    
    open(unit=307,file="Force_in_node.dat")
    
    do i = 1,N_node
       Force_node(i,1:2) = - gradient(i,1:2)
       write(unit=307,fmt=*) Force_node(i,1:N_dmnsn),i,"node_nm"
    enddo
    
    close(307)
    
  end subroutine get_forces_at_nodes

  
  
  subroutine get_force_in_a_spr(spr_nmbr)
    implicit none
    
    integer :: spr_nmbr
    integer :: entity
    integer :: entity_nmbr
    
    entity = 1
    entity_nmbr = spr_nmbr
    
    call switch_off_every_strength_except_one(entity,entity_nmbr)
    call get_gradient(node_xy,l0,A0,gradient)
    
    
    call get_forces_from_gradient_at_nodes(gradient)
    
    Force_spr(spr_nmbr,1:N_dmnsn) =  force_in_a_spr(spr_nmbr)
    
    call restore_coeffs_frm_tmp
    call dealloc_tmp_cell_variables
    
  contains

    function force_in_a_spr(spr_nmbr)
      implicit none
      real*8  :: force_in_a_spr(1:N_dmnsn)
      integer :: spr_nmbr
      
      integer :: i,imax
      integer :: node_nmbr
      
      force_in_a_spr(1:N_dmnsn) = 0.0d0
      
      imax = spr_node(spr_nmbr,0)
      
      do i = 1,imax
         node_nmbr = spr_node(spr_nmbr,i)
         force_in_a_spr(1:2) = force_in_a_spr(1:2) + Force_node(node_nmbr,1:2)
      enddo
      
    end function force_in_a_spr
    
  end subroutine get_force_in_a_spr



  subroutine get_force_in_an_area(area_nmbr)
    implicit none
    integer :: area_nmbr
    integer :: entity
    integer :: entity_nmbr
    
    entity = 2
    entity_nmbr = area_nmbr
    
    call switch_off_every_strength_except_one(entity,entity_nmbr)
    call get_gradient(node_xy,l0,A0,gradient)

    call get_forces_from_gradient_at_nodes(gradient)
    
    Force_area(area_nmbr,1:N_dmnsn) = force_in_an_area(area_nmbr) !force_in_an_area is a function
    
    call restore_coeffs_frm_tmp
    call dealloc_tmp_cell_variables
    
  contains
    
    function force_in_an_area(area_nmbr)
      implicit none
      real*8  :: force_in_an_area(1:N_dmnsn)
      integer :: area_nmbr
      
      integer :: lp_cnt,i
      integer :: node_nmbr
      
      force_in_an_area = 0.0d0
      lp_cnt = area_node(area_nmbr,0)
      
      do i = 1,lp_cnt
         node_nmbr = area_node(area_nmbr,0)
         force_in_an_area(1:2) = force_in_an_area(1:2) + Force_node(node_nmbr,1:2)
      enddo
      
    end function force_in_an_area
    
  end subroutine get_force_in_an_area
  
  
  
  subroutine get_force_dueto_Cg_atNode(node_nmbr)    
    implicit none    
    integer :: node_nmbr
    integer :: entity
    integer :: entity_nmbr
    
    
    entity = 3
    entity_nmbr = node_nmbr
    
    call switch_off_every_strength_except_one(entity,entity_nmbr)
    call get_gradient(node_xy,l0,A0,gradient)
    
    call get_forces_from_gradient_at_nodes(gradient)
    
    Force_grvty(node_nmbr,1:N_dmnsn) = force_at_node_for_Cg(node_nmbr) !force_at_node_for_Cg is a function
    
    call restore_coeffs_frm_tmp
    call dealloc_tmp_cell_variables
    
  contains
    
    function force_at_node_for_Cg(node_nmbr)
      implicit none
      real*8  :: force_at_node_for_Cg(1:N_dmnsn)
      integer :: node_nmbr
      
      force_at_node_for_Cg = 0.0d0
      
      force_at_node_for_Cg(1:2) = Force_node(node_nmbr,1:2)
      
    end function force_at_node_for_Cg
    
  end subroutine get_force_dueto_Cg_atNode
  
  
  subroutine get_force_dueto_bend_atPhi(phi_nm)
    implicit none
    integer :: phi_nm
    integer :: entity
    integer :: entity_nmbr

    entity = 4
    entity_nmbr = phi_nm
    
    call switch_off_every_strength_except_one(entity,entity_nmbr)
    call get_gradient(node_xy,l0,A0,gradient)
    
    call get_forces_from_gradient_at_nodes(gradient)

    
    Force_bend(phi_nm,1:N_dmnsn) = force_dueto_a_Phi(phi_nm)  
    
    call restore_coeffs_frm_tmp
    call dealloc_tmp_cell_variables

  contains

    function force_dueto_a_Phi(phi_nm)
      implicit none
      integer :: phi_nm
      real*8  :: force_dueto_a_Phi(1:N_dmnsn)
      integer :: i,imax
      integer :: node_nm,area_nm,area_serial
      integer :: neigh1,neigh2
      
      force_dueto_a_Phi(1:N_dmnsn) = 0.0d0
      
      node_nm = phi_info(phi_nm,1)
      area_nm = phi_info(phi_nm,2)
      call get_area_serial(node_nm,area_nm,area_serial)
      
      neigh1 = Nlist(node_nm,area_serial,1)
      neigh2 = Nlist(node_nm,area_serial,2)

      imax = 3
      
      do i = 1,imax
         if (i==1) then
           force_dueto_a_Phi(1:2) = force_dueto_a_Phi(1:2) + Force_node(node_nm,1:2)
        elseif (i==2) then
           force_dueto_a_Phi(1:2) = force_dueto_a_Phi(1:2) + Force_node(neigh1,1:2)
        elseif (i==3) then
           force_dueto_a_Phi(1:2) = force_dueto_a_Phi(1:2) + Force_node(neigh2,1:2)
        endif
     enddo
     
    end function force_dueto_a_Phi
    
  end subroutine get_force_dueto_bend_atPhi

  
  subroutine switch_off_every_strength_except_one(entity,entity_nmbr)
    implicit none
    integer :: entity ! entity 1 is spr,2 is area,3 in grvty,4 is bend
    integer :: entity_nmbr

    call alloc_and_init_tmp_cell_variables
    call store_coeffs_to_tmp
    call set_all_coeffs_to_zero_value
    call set_value_to_one_entity(entity,entity_nmbr)
      
  end subroutine switch_off_every_strength_except_one
  
  
end module force_calc









