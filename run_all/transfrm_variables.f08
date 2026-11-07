module transfrm_info
  use sys_building_info
  use system_parameters
  use generating_the_shape!inside that get_cnnctng_node(blck_id,prev_blck_id)
  
  implicit none
  
  integer :: max_node_spr,max_node_sprS,max_node_sprS2,max_node_sprSNI
  integer :: max_area_spr,max_area_sprS,max_area_sprS2,max_area_sprSNI
  integer :: max_node_nospr
  
  integer :: max_spr_node, max_spr_nodeS, max_spr_nodeS2, max_spr_nodeSNI
  integer :: max_area_node,max_area_nodeS,max_area_nodeS2,max_area_nodeSNI
  
  integer :: max_node_area,max_node_areaS,max_node_areaS2,max_node_areaSNI
  integer :: max_spr_area, max_spr_areaS, max_spr_areaS2, max_spr_areaSNI
  
  
  integer, allocatable :: spr_node(:,:),spr_nodeS(:,:),spr_nodeSNI(:,:),spr_nodeS2(:,:)
  integer, allocatable :: spr_area(:,:),spr_areaS(:,:),spr_areaSNI(:,:),spr_areaS2(:,:)
  integer, allocatable :: sprless_node(:,:)
  
  integer, allocatable :: node_spr(:,:), node_sprS(:,:), node_sprSNI(:,:) ,node_sprS2(:,:)
  integer, allocatable :: node_area(:,:),node_areaS(:,:),node_areaSNI(:,:),node_areaS2(:,:)
  
  integer, allocatable :: area_spr(:,:), area_sprS(:,:), area_sprSNI(:,:) , area_sprS2(:,:)
  integer, allocatable :: area_node(:,:),area_nodeS(:,:),area_nodeSNI(:,:), area_nodeS2(:,:)
  
  integer, allocatable :: Nlist(:,:,:),NlistS(:,:,:),NlistSNI(:,:,:),NlistS2(:,:,:)
  integer, allocatable :: phi_info(:,:)
  
  integer, allocatable :: trmnl_intrmdSpr(:,:)
  integer, allocatable :: insrtd_vrtxNode(:,:)
  
contains
  
  subroutine allocate_and_initialize_transfrm_variables
    implicit none
    
    write(*,*) modelID,"modelID in allocate_and_initialize_transfrm_variables"
    if (modelID .ne. 1) stop 'model ID must be 1 in Here in allocate_and_initialize_transfrm_variables'
    
    max_node_spr   = 3
    max_node_nospr = 2
    max_area_spr   = 2
    
    max_spr_node  = 6 !max no of spr  cnnctd to a node
    max_area_node = 4 !max no of area cnnctd to a node
    
    max_node_area = 5
    max_spr_area  = 4 !for everyone 4 except the boundary ones
    
    allocate(spr_node(1:N_spr,0:max_node_spr),spr_nodeS(1:N_spr,0:max_node_spr))
    allocate(sprless_node(1:N_nospr,0:max_node_nospr))
    allocate(spr_area(1:N_spr,0:max_area_spr),spr_areaS(1:N_spr,0:max_area_spr))
    
    allocate(node_spr(1:N_node,0:max_spr_node),  node_sprS(1:N_node,0:max_spr_node))
    allocate(node_area(1:N_node,0:max_area_node),node_areaS(1:N_node,0:max_spr_node))
    
    allocate(area_spr(1:N_cell,0:max_spr_area),  area_sprS(1:N_cell,0:max_spr_area))
    allocate(area_node(1:N_cell,0:max_node_area),area_nodeS(1:N_cell,0:max_node_area))
    
    allocate(Nlist(1:N_node,1:max_area_node,1:2),NlistS(1:N_node,1:max_area_node,1:2))
    
    spr_node  = -1 ; spr_nodeS  = -1
    spr_area  = -1 ; spr_areaS  = -1
    
    node_spr  = -1 ; node_sprS  = -1
    node_area = -1 ; node_areaS = -1 
    
    area_node = -1 ; area_nodeS = -1
    area_spr  = -1 ; area_sprS  = -1
    
    Nlist     = -1 ; NlistS     = -1
    
  end subroutine allocate_and_initialize_transfrm_variables
  
  
  subroutine allocate_and_initialize_phiInfo !trnsfrm_var dependent info 
    implicit none
    
    call get_Nphi
    allocate(phi_info(1:N_phi,1:2))
    
    phi_info = -1
    
  end subroutine allocate_and_initialize_phiInfo


  subroutine deallocate_transfrm_variables
    implicit none
    
    deallocate(spr_node,spr_nodeS)
    deallocate(sprless_node)
    deallocate(spr_area,spr_areaS)
    
    deallocate(node_spr,node_sprS)
    deallocate(node_area,node_areaS)
    
    deallocate(area_spr,area_sprS)
    deallocate(area_node,area_nodeS)
    
    deallocate(Nlist,NlistS)

  end subroutine deallocate_transfrm_variables
  
  
  subroutine deallocate_and_reallocate_transfrmStr_variables
    implicit none
    
    deallocate(spr_nodeS,spr_areaS)
    deallocate(node_sprS,node_areaS)
    deallocate(area_sprS,area_nodeS)
    
    deallocate(NlistS)
    
    write(*,*) max_node_spr,max_area_spr,max_spr_node,max_area_node,max_node_area,&
         max_spr_area,"max values"
    
    allocate(spr_nodeS(1:N_spr,0:max_node_spr))
    allocate(spr_areaS(1:N_spr,0:max_area_spr))
    
    allocate(node_sprS(1:N_node,0:max_spr_node))
    allocate(node_areaS(1:N_node,0:max_spr_node))
    
    allocate(area_sprS(1:N_cell,0:max_spr_area))
    allocate(area_nodeS(1:N_cell,0:max_node_area))
    
    allocate(NlistS(1:N_node,1:max_area_node,1:2))
    
    spr_nodeS  = -1 ; spr_areaS  = -1
    node_sprS  = -1 ; node_areaS = -1 
    area_nodeS = -1 ; area_sprS  = -1
    
    NlistS     = -1
    
  end subroutine deallocate_and_reallocate_transfrmStr_variables
  
  subroutine deallocate_phiInfo !trnsrm_vars_dependent_var
    implicit none

    deallocate(phi_info)

  end subroutine deallocate_phiInfo

  subroutine get_Nphi
    implicit none
    integer :: i,j
    
    N_phi = 0
    
    do i = 1,N_cell
       N_phi = N_phi + area_node(i,0)
    enddo
    
  end subroutine get_Nphi
  
  subroutine is_it_a_trminl_node(node_nm,spr_nm,lgcl_trminl)
    implicit none
    integer :: node_nm,spr_nm
    logical :: lgcl_trminl,lgcl_dn
    
    integer :: i,imax
    integer :: pstn
    integer :: other_node
    
    lgcl_trminl = .False.
    lgcl_dn     = .False.
    
    pstn = 0
    
    imax = spr_node(spr_nm,0)
    
    do i = 1,imax
       if (node_nm .eq. spr_node(spr_nm,i)) then
          pstn = i
          exit
       endif
    enddo
    
    if (spr_node(spr_nm,0) == 3) then
       
       if (pstn==1 .or. pstn==3) then
          lgcl_trminl = .True.
       elseif (pstn==2) then
          lgcl_trminl = .False.
       elseif (pstn==0) then
          call is_it_a_dn(node_nm,lgcl_dn)
       endif
       
    elseif (spr_node(spr_nm,0) == 4) then
       
       if (pstn==1 .or. pstn==4) then
          lgcl_trminl = .True.
       elseif (pstn==2 .or. pstn==3) then
          lgcl_trminl = .False.
       elseif (pstn==0) then
          call is_it_a_dn(node_nm,lgcl_dn)
       endif
       
    endif
       
    
    if (lgcl_dn .eqv. .True.) then
       call get_the_other_dn(node_nm,other_node)
       
       do i = 1,max_node_spr
          if (other_node .eq. spr_node(spr_nm,i)) then
             pstn = i
             exit
          endif
       enddo
       
       if (spr_node(spr_nm,0) == 3) then
          if (pstn==1 .or. pstn==3) then
             lgcl_trminl = .True.
          elseif (pstn==2) then
             lgcl_trminl = .False.
          endif
       elseif (spr_node(spr_nm,0) == 4) then
          if (pstn==1 .or. pstn==4) then
             lgcl_trminl = .True.
          elseif (pstn==2 .or. pstn==3) then
             lgcl_trminl = .False. 
          endif
       endif
       
    endif
    
  end subroutine is_it_a_trminl_node
  
end module transfrm_info
