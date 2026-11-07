module storing_changing_restoring_routines
  use conversion_routines !=cell_info+system_parameters+moving_coordnte_info
  implicit none
  real*8, allocatable  :: node_xy_tmp(:,:)
  real*8, allocatable  :: node_xy_cmpr(:,:)
  real*8, allocatable  :: coordntes_zero_grv(:)
  
  real*8, allocatable  :: A0_tmp(:)
  !real*8, allocatable :: delX(:),delY(:)
  integer, allocatable :: node_typ_tmp(:)
  integer              :: readNodeFilesTyp=0
  
contains
  
  subroutine allocating_variables_scrr
    implicit none
    
    continue
    
  end subroutine allocating_variables_scrr
  
  subroutine store_node_xy
    implicit none
    
    allocate(node_xy_tmp(1:N_node,1:N_dmnsn))
    
    node_xy_tmp = -1.0d30
    node_xy_tmp = node_xy
    
  end subroutine store_node_xy
  
  subroutine restore_node_xy
    implicit none
    
    node_xy = node_xy_tmp
    deallocate(node_xy_tmp)
    
  end subroutine restore_node_xy
  
  subroutine store_and_change_A0_area(prcnt_chng)
    implicit none
    real*8  :: prcnt_chng ! can be (+)ve or (-)ve
    !real*8  :: A0_tmp(1:N_cell)
    
    allocate(A0_tmp(1:N_cell))
    
    A0_tmp = -1e30
    A0_tmp = A0
    A0     = (1.0 + prcnt_chng) * A0
    
  end subroutine store_and_change_A0_area
  
  
  subroutine restore_A0_area
    implicit none
    
    A0 = A0_tmp
    deallocate(A0_tmp)
    
  end subroutine restore_A0_area
  
  
  subroutine store_and_change_node_type(which_type,changed_type)
    implicit none
    integer, intent(in) :: which_type
    integer, intent(in) :: changed_type

    integer :: i
    
    allocate(node_typ_tmp(1:N_node))

    node_typ_tmp = node_typ

    do i = 1,N_node
       if (node_typ(i) .eq. which_type) then
          node_typ(i) = changed_type
       endif  
    enddo

  end subroutine store_and_change_node_type
  
  
  subroutine store_and_change_mltipl_node_types(hmt,which_types,changed_types)
    
    implicit none
    integer, intent(in) :: hmt !hmt=how_many_types_to_be_changed
    integer, intent(in) :: which_types(1:hmt)
    integer, intent(in) :: changed_types(1:hmt)
    
    integer :: i,j
    
    allocate(node_typ_tmp(1:N_node))
    
    node_typ_tmp = node_typ
    
    do i = 1,N_node
       do j = 1,hmt
          if (node_typ(i) .eq. which_types(j)) then
             node_typ(i) = changed_types(j)
          endif
       enddo
    enddo
    
  end subroutine store_and_change_mltipl_node_types

  
  subroutine restore_node_type
    implicit none

    node_typ = node_typ_tmp
    deallocate(node_typ_tmp)

  end subroutine restore_node_type

  
  subroutine alloc_and_init_node_xy_cmpr_coordntes_zero_grv
    implicit none

    allocate(node_xy_cmpr(1:N_node,1:N_dmnsn))
    node_xy_cmpr = -1e30

    allocate(coordntes_zero_grv(1:N_mvCoordnte))
    coordntes_zero_grv = -1e30
    
  end subroutine alloc_and_init_node_xy_cmpr_coordntes_zero_grv

  
  subroutine generate_delX_and_delY(dum_coordntes,test_no,lp_cnt,N_iter)

    implicit none
    integer, intent(in) :: test_no,lp_cnt
    
    real*8  :: dum_coordntes(1:N_mvCoordnte)
    real*8  :: delX(1:N_node),delY(1:N_node)
    real*8  :: nonD_delA0_area
    
    integer, intent(in) :: N_iter
    
    call coordntes_to_nodes(dum_coordntes,node_xy)
   
    !write(*,*) node_xy,"node_xy_in_gen_delx_dely"
    !write(*,*) node_xy_tmp,"node_xy_tmp_in_gen_delx_dely"

    
    nonD_delA0_Area = -(int(N_iter/2))*0.1d0 + (lp_cnt-1)*0.1d0
    
    call get_delX_and_delY(node_xy,node_xy_tmp,delX,delY)
    
    if (lp_cnt .eq. 1) then
       !write(*,*) delX(1:N_node)
       !write(*,*) delY(1:N_node)
    endif
       
    call writing_delX_and_delY(delX,delY,test_no,lp_cnt,nonD_delA0_area)
    
  end subroutine generate_delX_and_delY
  
  
  
  
  
  subroutine save_config_and_generate_data(dum_coordntes,test_no,subtest_no,lp_cnt)
    implicit none
    real*8  :: dum_coordntes(1:N_mvCoordnte)
    integer, intent(in) :: test_no,subtest_no,lp_cnt
    
    !write(*,*) lp_cnt,"LOOP COUNT"

    !if (lp_cnt .eq. 8) then
     !  write(*,*) node_xy(9,1:2),"node9"
      ! write(*,*) node_xy(10,1:2),"node10"
    !endif
    
    call coordntes_to_nodes(dum_coordntes,node_xy)
    
    ! if (lp_cnt .eq. 8) then
    !    write(*,*) node_xy(9,1:2),"node9"
    !    write(*,*) node_xy(10,1:2),"node10"
    !    write(*,*) dum_coordntes,"dum_coor_in_save_config"
    !    write(*,*) node_xy,"node_xy_in_save_config"
    ! endif
    
    call writing_nodes_over_spr(node_xy,test_no,subtest_no,lp_cnt)
    call writing_nodes_with_number_and_type(node_xy,test_no,subtest_no,lp_cnt)
    call writing_nodes_over_area(node_xy,test_no,subtest_no,lp_cnt)
    
  end subroutine save_config_and_generate_data
  
  
  
  subroutine save_config_and_generate_ColorData(dum_coordntes,duml0,dumA0,ExpNo,frmeNo)
    implicit none
    real*8 , intent(in) :: dum_coordntes(1:N_mvCoordnte)
    real*8 , intent(in) :: duml0(1:N_spr),dumA0(1:N_cell)
    integer, intent(in) :: ExpNo,frmeNo
    
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8  :: Cell_Pressure(1:N_cell)
    real*8  :: Tension(1:N_spr)
    integer :: i
    integer :: ExpNoStr,ExpNm
    integer :: N_lon
    
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
    
    call coordntes_to_nodes(dum_coordntes,dum_Nodes)
    
    if (ExpNo==37) then
       write(*,*) dum_coordntes(1),"FrmNo=",frmeNo
       write(*,*) dum_Nodes(1,1:2),"FrmNo=",frmeNo
       write(*,*) node_xy(1,1:2),"FrmNo=",frmeNo
    endif
    
    open(unit=21,file="TnsnCrossCheck.dat",position="append")
    open(unit=22,file="PresCrossCheck.dat",position="append")
    open(unit=323,file='ll0_kspr_Lt21_24_26_29.dat',position='append')
    
    if (ExpNo==16) write(323,fmt=*) ExpNo,frmeNo,"Exp,frme"
    
    Tension(1:N_spr)        = TnsnComprsn(dum_Nodes,duml0) 
    Cell_Pressure(1:N_cell) = Pressure(dum_Nodes,dumA0)  
    
    do i=1,N_spr
       write(21,fmt=*) Tension(i),i,k_spr(i),l0(i),l(i),"Tension,spr_nm"
       !write(*,*) Tension(i),i,k_spr(i),l0(i),l(i),"Tension,spr_nm"
    enddo
    
    do i=1,N_cell
       write(22,fmt=*) Cell_Pressure(i),i,"Pressure,area_nm"
    enddo
    
    close(21)
    close(22)
    close(323)
    
    call writing_nodes_for_Exps(dum_Nodes,ExpNo,frmeNo)
    call writing_nodes_over_spr_single_node(dum_Nodes,Tension,ExpNo,frmeNo)
    call writing_nodes_over_area_single_node(dum_Nodes,Cell_Pressure,ExpNo,frmeNo)
    
    !if (frmeNo==31) then
    !call writing_datas_for_curve(dum_Nodes,Cell_Pressure,Tension,ExpNo,frmeNo)
    !endif
    
    if (modelID==2) then
       !ExpNm = 22
       !N_lon = 9
       !call smooth_image_with_cubic_fit(dum_Nodes,N_lon,Tension,Cell_Pressure,ExpNm,frmeNo)
       continue
    endif
    
  end subroutine save_config_and_generate_ColorData
  
  
  subroutine read_config_and_start_simlnFrm_there(ExpNo,FrmNoToBeRead)
    implicit none
    integer, intent(in) :: ExpNo
    integer, intent(in) :: FrmNoToBeRead
    
    integer :: objTyp,NumChar
    
    character(len=200)              :: flNmN,flNmS,flNmA
    character(len=200), allocatable :: flNmN_arr(:),flNmS_arr(:),flNmA_arr(:)
    
    NumChar = 1
    allocate(flNmN_arr(1:NumChar),flNmS_arr(1:NumChar),flNmA_arr(1:NumChar))
    
    objTyp = 1
    call manage_ExpsName_withlen200(ExpNo,FrmNoToBeRead,objTyp,flnmN)
    flNmN_arr(1)=flNmN
    
    if (readNodeFilesTyp==0) call read_nodeXY_and_CgVals_frmAFile(NumChar,flNmN)
    if (readNodeFilesTyp==1) call read_nodeXY_and_CgVals_frm_PrePulleyToPostPulley(NumChar,flNmN)
    
    objTyp = 2
    call manage_ExpsName_withlen200(ExpNo,FrmNoToBeRead,objTyp,flnmS)
    flNmS_arr(1)=flNmS
    call read_kspr_and_l0vals_frmAFile(NumChar,flNmS)
    
    objTyp = 5
    call manage_ExpsName_withlen200(ExpNo,FrmNoToBeRead,objTyp,flnmA)
    flNmA_arr(1)=flNmA
    call read_karea_and_A0vals_frmAFile(NumChar,flNmA)
    
  end subroutine read_config_and_start_simlnFrm_there
  
  
  subroutine read_config_and_start_simlnFrm_there_flpceLoc(flplce,ExpNo,FrmNoToBeRead)
    implicit none
    character(len=200), intent(in)  :: flplce
    integer,   intent(in)           :: ExpNo
    integer,   intent(in)           :: FrmNoToBeRead
    
    integer                         :: objTyp,NumChar
    character(len=200)              :: flNmN,flNmS,flNmA
    character(len=200), allocatable :: flNmN_arr(:),flNmS_arr(:),flNmA_arr(:)
    
    NumChar = 2
    allocate(flNmN_arr(1:NumChar),flNmS_arr(1:NumChar),flNmA_arr(1:NumChar))
    
    objTyp = 1
    call manage_ExpsName_withlen200(ExpNo,FrmNoToBeRead,objTyp,flNmN)
    flNmN_arr(1)=flplce ; flNmN_arr(2)=flNmN
    
    write(*,*)trim(adjustl(flplce)),"N0"
    write(*,*)trim(adjustl(flNmN_arr(1))),"N1"
    write(*,*)trim(adjustl(flNmN_arr(2))),"N2"
    
    call read_nodeXY_and_CgVals_frmAFile(NumChar,flNmN_arr)
    
    objTyp = 2
    call manage_ExpsName_withlen200(ExpNo,FrmNoToBeRead,objTyp,flnmS)
    flNmS_arr(1)=flplce ; flNmS_arr(2)=flNmS
    
    write(*,*)trim(adjustl(flplce)),"S0"
    write(*,*)trim(adjustl(flNmS_arr(1))),"S1"
    write(*,*)trim(adjustl(flNmS_arr(2))),"S2"
    
    call read_kspr_and_l0vals_frmAFile(NumChar,flNmS_arr)
    
    objTyp = 5
    call manage_ExpsName_withlen200(ExpNo,FrmNoToBeRead,objTyp,flnmA)
    flNmA_arr(1)=flplce ; flNmA_arr(2)=flNmA

    write(*,*)trim(adjustl(flplce)),"A0"
    write(*,*)trim(adjustl(flNmA_arr(1))),"A1"
    write(*,*)trim(adjustl(flNmA_arr(2))),"A2"
    
    call read_karea_and_A0vals_frmAFile(NumChar,flNmA_arr)
    
  end subroutine read_config_and_start_simlnFrm_there_flpceLoc
  
  
  subroutine smooth_image_with_cubic_fit(dum_Nodes,N_lon,Tension,Cell_Pressure,ExpNo,frmeNo)
    use system_parameters
    use transfrm_info
    
    implicit none
    real*8, intent(in) :: dum_Nodes(1:N_node,1:N_dmnsn)
    integer,intent(in) :: N_lon
    real*8, intent(in) :: Tension(1:N_spr)
    real*8, intent(in) :: Cell_Pressure(1:N_cell)
    integer,intent(in) :: ExpNo,frmeNo
    
    character(len=100) :: flnm!1,flnm2
    character(len=100) :: flplce!1,flplce2
    
    integer :: i,j,jmax,m
    integer :: spr1,spr2
    integer :: node1,node2,node3
    real*8  :: N1xy(1:N_dmnsn),N2xy(1:N_dmnsn),N3xy(1:N_dmnsn)
    
    integer :: strtNode,fnshNode
    integer :: spr_nm,node_nm,area_nm
    integer :: objTyp
    
    real*8  :: TnsnWall
    real*8  :: LoN(1:N_lon,1:N_dmnsn)
    
    objTyp = 2
    
    call manage_ExpsName(ExpNo,frmeNo,objTyp,flnm)
    
    write(flplce,*) "/home/rniloy/PROJECT_CELL/2D_codes/unit_blck/restructuring_Data_Structure/"
    
    !open(unit=31,file=trim(adjustl(flplce))//trim(adjustl(flnm)))
    open(unit=31,file=trim(adjustl(directryFlnm))) 
    
    TnsnWall = 0.00d0
    
    do i = 0,(N_sprS+1)

       if (i==0) then
          strtNode = 1 ; fnshNode = 2
          write(unit=31,fmt=*) TnsnWall
          write(unit=31,fmt=*) dum_Nodes(strtNode,1:N_dmnsn)
          write(unit=31,fmt=*) dum_Nodes(fnshNode,1:N_dmnsn)
          write(unit=31,fmt=*) " "
          
       elseif ((i.ne.0) .AND. (i.ne.(N_sprS+1))) then
          
          if (is_it_a_curve(i)==1) then
             
             spr1  = trmnl_intrmdSpr(i,1)
             spr2  = trmnl_intrmdSpr(i,2)
             
             node1 = spr_node(spr1,1) !inside loop if more than 1 node added
             node2 = spr_node(spr1,2)
             node3 = spr_node(spr2,2)
             
             N1xy = dum_Nodes(node1,1:N_dmnsn)
             N2xy = dum_Nodes(node2,1:N_dmnsn)
             N3xy = dum_Nodes(node3,1:N_dmnsn)
             
             call get_cubic_fit(N1xy,N2xy,N3xy,N_lon,LoN)
             
             do j = 1,N_lon
                
                if (j.le.(N_lon/2)) then
                   if (j==1) write(unit=31,fmt=*) Tension(spr1)
                   
                   write(unit=31,fmt=*) LoN(j,1:N_dmnsn)
                   
                   if (j==(N_lon/2)) then
                      write(unit=31,fmt=*) LoN((j+1),1:N_dmnsn)
                      write(unit=31,fmt=*) " "
                   endif
                   
                   
                elseif (j.gt.(N_lon/2)) then
                   if (j==((N_lon/2)+1)) write(unit=31,fmt=*) Tension(spr2)
                   
                   write(unit=31,fmt=*) LoN(j,1:N_dmnsn)
                   
                   if (j==N_lon) write(unit=31,fmt=*) " "
                   
                endif !!!if block will be looped if more than 1 node added
                
             enddo
             
             
          elseif (is_it_a_curve(i)==0) then
             
             spr_nm = trmnl_intrmdSpr(i,1)
             jmax   = spr_node(spr_nm,0)
             
             if (jmax==2) then
                write(unit=31,fmt=*) Tension(spr_nm)
                
                do j = 1,jmax
                   node_nm = spr_node(spr_nm,j)
                   write(unit=31,fmt=*) dum_Nodes(node_nm,1:N_dmnsn)
                enddo
                
                write(unit=31,fmt=*) " "
                
             elseif (jmax==3) then
                
                write(unit=31,fmt=*) Tension(spr_nm)
                
                do j = 1,2
                   node_nm = spr_node(spr_nm,j)
                   write(unit=31,fmt=*) dum_Nodes(node_nm,1:N_dmnsn)
                enddo
                write(unit=31,fmt=*) " "
                
                write(unit=31,fmt=*) Tension(spr_nm)
                do j = 2,3
                   node_nm = spr_node(spr_nm,j)
                   write(unit=31,fmt=*) dum_Nodes(node_nm,1:N_dmnsn)
                enddo
                write(unit=31,fmt=*) " "
                
             endif
          endif
          
       elseif (i==(N_sprS+1)) then
          
          strtNode = (2*nvsl)+1 ; fnshNode = (2*nvsl)+2
          write(unit=31,fmt=*) TnsnWall
          write(unit=31,fmt=*) dum_Nodes(strtNode,1:N_dmnsn)
          write(unit=31,fmt=*) dum_Nodes(fnshNode,1:N_dmnsn)
          write(unit=31,fmt=*) " "
          
       endif
             
       !if (i.ne.0 .AND. i.ne.(N_sprS+1)) then
       !  write(*,*) is_it_a_curve(i),i,"IIAC"
       !endif
       
    enddo
    
    close(31)
    
    objTyp = 3
    
    call manage_ExpsName(ExpNo,frmeNo,objTyp,flnm)
    
    write(flplce,*) "/home/rniloy/PROJECT_CELL/2D_codes/unit_blck/restructuring_Data_Structure/"
    
    !open(unit=31,file=trim(adjustl(flplce))//adjustl(flnm))
    open(unit=31,file=trim(adjustl(directryFlnm)))
    
    do i = 1,N_cell
       write(unit=31,fmt=*) Cell_Pressure(i)
       
       area_nm = i   
       jmax = area_node(i,0)
       
       do j = 1,jmax
          node_nm = area_node(i,j)
          
          if (node_nm .gt. N_nodeS) then
             
             if ((j-1).ne.0) then
                node1 = area_node(i,(j-1))
             elseif ((j-1)==0) then
                node1 = area_node(i,jmax)
             endif
             
             node2 = area_node(i,j)
             
             if ((j+1).le.jmax) then
                node3 = area_node(i,(j+1))
             elseif ((j+1).gt.area_node(i,0)) then
                node3 = area_node(i,1)
             endif
             
             N1xy = dum_Nodes(node1,1:N_dmnsn)
             N2xy = dum_Nodes(node2,1:N_dmnsn)
             N3xy = dum_Nodes(node3,1:N_dmnsn)

             call get_cubic_fit(N1xy,N2xy,N3xy,N_lon,LoN)
             
             do m = 2,(N_lon-1)
                write(unit=31,fmt=*) LoN(m,1:N_dmnsn)
             enddo
             
          elseif (node_nm .le. N_nodeS) then          
             write(unit=31,fmt=*) dum_Nodes(node_nm,1:N_dmnsn)
          endif
          
       enddo
       write(unit=31,fmt=*) " "
    enddo
    
    close(31)
    
    
  end subroutine smooth_image_with_cubic_fit
  
  subroutine save_config_and_generate_data_InsideSubtst(dum_coordntes,test_no,subtest_no,intraSubtst_no,lp_cnt)
    implicit none
    real*8  :: dum_coordntes(1:N_mvCoordnte)
    integer, intent(in) :: test_no,subtest_no,intraSubtst_no,lp_cnt
    
    !write(*,*) lp_cnt,"LOOP COUNT"
    
    call coordntes_to_nodes(dum_coordntes,node_xy)
    
    ! if (lp_cnt .eq. 1) then
    !    write(*,*) dum_coordntes,"dum_coor_in_save_config"
    !    write(*,*) node_xy,"node_xy_in_save_config"
    ! endif
    
    call writing_nodes_over_spr(node_xy,test_no,subtest_no,lp_cnt)
    call writing_nodes_with_number_and_type(node_xy,test_no,subtest_no,lp_cnt)
    call writing_nodes_over_area(node_xy,test_no,subtest_no,lp_cnt)
    
  end subroutine save_config_and_generate_data_InsideSubtst





  
  
  subroutine generate_delX_delY_CgX(dum_coordntes,dum_coordntes_zero_grv,test_no,lp_cnt,N_iter,CgX)
    
    implicit none
    integer, intent(in) :: test_no,lp_cnt
    
    real*8  :: dum_coordntes(1:N_mvCoordnte)
    real*8  :: dum_coordntes_zero_grv(1:N_mvCoordnte)
    real*8  :: delX(1:N_node),delY(1:N_node)
    !real*8  :: nonD_g
    
    real*8,intent(in)   :: CgX
    !real*8,intent(in)   :: CgX_incr
    integer, intent(in) :: N_iter
   
    !write(*,*) " "
    !write(*,*) dum_coordntes,"coor_in delxdely_CgX"
    !write(*,*) " "
    !write(*,*) dum_coordntes_zero_grv,"coor_zero_grv_indelxdely_CgX"
    
    call coordntes_to_nodes(dum_coordntes,node_xy)
    call coordntes_to_nodes(dum_coordntes_zero_grv,node_xy_cmpr)
    
    
    !nonD_g = lp_cnt*(1/N_iter)
    
    !CgX = (lp_cnt-1)*CgX_incr
    
    call get_delX_and_delY(node_xy,node_xy_cmpr,delX,delY)
    
    if (lp_cnt .eq. 1) then
       !write(*,*) delX(1:N_node)
       !write(*,*) delY(1:N_node)
    endif
    
    call writing_delX_delY_and_CgX(delX,delY,test_no,lp_cnt,CgX)
    
  end subroutine generate_delX_delY_CgX
  
  
  subroutine generate_delX_delY_CgX_pairwise(dum_coordntes,dum_coordntes_zero_grv,&
       test_no,lp_cnt,N_iter,CgX)
    
    implicit none
    integer, intent(in) :: test_no,lp_cnt
    
    real*8  :: dum_coordntes(1:N_mvCoordnte)
    real*8  :: dum_coordntes_zero_grv(1:N_mvCoordnte)
    !real*8  :: delX(1:N_node),delY(1:N_node)
    !real*8  :: nonD_g
    
    real*8,intent(in)   :: CgX
    !real*8,intent(in)   :: CgX_incr
    integer, intent(in) :: N_iter
    integer :: hlf_Node
    
    real*8,allocatable  :: nodes_pw(:,:),nodes_cmpr_pw(:,:)
    real*8,allocatable  :: delX_pw(:),delY_pw(:) 
    
    write(*,*) N_iter,"N_iter"
    
    !write(*,*) " "
    !write(*,*) dum_coordntes,"coor_in delxdely_CgX"
    !write(*,*) " "
    write(*,*) dum_coordntes_zero_grv,"coor_zero_grv_indelxdely_CgX"
    
    call coordntes_to_nodes(dum_coordntes,node_xy)
    call coordntes_to_nodes(dum_coordntes_zero_grv,node_xy_cmpr)
    
    if (mod(N_node,2).eq.0) then
       hlf_Node = N_node/2
    elseif (mod(N_node,2).eq.1) then
       hlf_Node = 0
       write(*,*) "Odd nodes are not compatible in sbrtn: generate_delX_delY_CgX_pairwise"
    endif
    
    allocate(nodes_pw(1:hlf_Node,1:N_dmnsn))
    allocate(nodes_cmpr_pw(1:hlf_Node,1:N_dmnsn))
    allocate(delX_pw(1:hlf_Node),delY_pw(1:hlf_Node))
    
    !nonD_g = lp_cnt*(1/N_iter)
    
    !CgX = (lp_cnt-1)*CgX_incr
    call get_pw_nodes!(node_xy,node_xy_cmpr,nodes_pw,nodes_cmpr_pw)
    call get_delX_and_delY_pairwise(nodes_pw,nodes_cmpr_pw,delX_pw,delY_pw,hlf_node)
    
    if (lp_cnt .eq. 1) then
       !write(*,*) delX(1:N_node)
       !write(*,*) delY(1:N_node)
    endif
    
    call writing_delX_delY_and_CgX_pairwise(delX_pw,delY_pw,test_no,lp_cnt,CgX,hlf_Node)
    

  contains
    
    subroutine get_pw_nodes!(node_xy,node_xy_cmpr,nodes_pw,nodes_cmpr_pw)
      
      implicit none
      integer :: i
      integer :: count_node
      
      count_node = 1
      
      do i = 1,hlf_Node         
         nodes_pw(i,1) = (node_xy(count_node,1) + node_xy((count_node+1),1))/2.0d0
         nodes_pw(i,2) = (node_xy(count_node,2) + node_xy((count_node+1),2))/2.0d0
         
         nodes_cmpr_pw(i,1) = (node_xy_cmpr(count_node,1) + node_xy_cmpr((count_node+1),1))/2.0d0
         nodes_cmpr_pw(i,2) = (node_xy_cmpr(count_node,2) + node_xy_cmpr((count_node+1),2))/2.0d0
         
         count_node = count_node + 2
         
         !write(*,*) nodes_pw(i,1:2),"pw"
         !write(*,*) nodes_cmpr_pw(i,1:2),"pw_cmpr"
         
      enddo
      
    end subroutine get_pw_nodes
    
  end subroutine generate_delX_delY_CgX_pairwise
  
  
  subroutine CgX_zero_data(dum_coordntes,dum_coordntes_zero_grv)
    implicit none
    real*8 :: dum_coordntes(1:N_mvCoordnte)
    real*8 :: dum_coordntes_zero_grv(1:N_mvCoordnte)
    
    dum_coordntes_zero_grv = dum_coordntes
    
  end subroutine CgX_zero_data
  
  
  
  
  subroutine print_the_properties
    implicit none
    integer :: i
    
    write(*,*) " "
    do i = 1,N_node
       write(*,*) node_xy(i,1:N_dmnsn),i,"nodes"
    enddo
    write(*,*) " "
    do i = 1,N_spr
       write(*,*) k_spr(i),l(i),l0(i),i,"sprs"
    enddo
    write(*,*) " "
    do i = 1,N_cell
       write(*,*) k_area(i),A(i),A0(i),i,"cells"
    enddo
    write(*,*) " "
    do i = 1,N_mvCoordnte
       write(*,*) coordntes_xy(i),i,"coor"
    enddo
    write(*,*) " "
    do i = 1,N_node
       write(*,*) CgXNode(i),CgYNode(i),i,"Cgs"
    enddo
    write(*,*) ""
    
  end subroutine print_the_properties
  
  
  subroutine get_fileNames_inside_the_modFile(ExpNo,frmeNo,objTyp,flnm)
    implicit none
    integer, intent(in)             :: ExpNo,frmeNo,objTyp
    character(len=100), intent(out) :: flnm
    
    call manage_ExpsName(ExpNo,frmeNo,objTyp,flnm)
    
  end subroutine get_fileNames_inside_the_modFile
  
  
  
  subroutine read_nodeXY_and_CgVals_frmAFile(NumChar,flnm_arr)
    implicit none
    integer, intent(in)            :: NumChar
    character(len=200), intent(in) :: flnm_arr(1:NumChar)
    
    real*8  :: CgXVal,CgXNval(1:N_node)
    real*8  :: CgYVal,CgYNval(1:N_node)
    integer :: nodeTval,nodeT(1:N_node)
    real*8  :: nodeRead(1:N_dmnsn),nodeVal(1:N_node,1:N_dmnsn)
    integer :: cnt,error
    integer :: nodeCnt,i
    integer :: node1,node2
    
    if (NumChar==1) then
       write(*,*)trim(adjustl(flnm_arr(1))),"flnm_arr In read_nodeXY 1"
       open(unit=449,file=trim(adjustl(flnm_arr(1))))
       
    elseif (NumChar==2) then
       write(*,*)trim(adjustl(flnm_arr(1))),"flnm_arr In read_nodeXY 21"
       write(*,*)trim(adjustl(flnm_arr(2))),"flnm_arr In read_nodeXY 22"
       open(unit=449,file=trim(adjustl(flnm_arr(1)))//trim(adjustl(flnm_arr(2))))
    endif
    
    cnt=1 ; nodeCnt=1
    nodeVal = -100.0d0
    
    do   
       read(unit=449,fmt=*,IOSTAT=error) nodeRead(1:N_dmnsn),nodeCnt,nodeTval,CgXval,CgYval
       write(*,*) nodeRead(1:N_dmnsn),nodeCnt,nodeTval,CgXval,CgYval,"p"
       
       if (error.lt.0) exit
       
       nodeVal(nodeCnt,1:N_dmnsn) = nodeRead(1:N_dmnsn)
       nodeT(nodeCnt)             = nodeTval 
       CgXNval(nodeCnt)           = CgXval
       CgYNval(nodeCnt)           = CgYval
       cnt = cnt+1
       
    enddo
    
    close(449)
    
    if (N_doubleNode==0) then
       continue
    elseif (N_doubleNode.gt.0) then
       
       do i = 1,N_doubleNode
          
          node1 = double_node(i,1) ; node2 = double_node(i,2)
          write(*,*) node1,node2,"doubleNode"
          nodeVal(node2,1:N_dmnsn) = nodeVal(node1,1:N_dmnsn)
          nodeT(node2)             = nodeT(node1)
          CgXNval(node2)           = CgXNval(node1)
          CgYNval(node2)           = CgYNval(node1)
          
       enddo
       
    endif
    
    if (nodeCnt.ne.N_node) then
       write(*,*) "nodeCnt MUST BE equal to N_node"
       write(*,*) nodeCnt,N_node,"N_node"
       stop
    endif
    
    node_xy(1:N_node,1:N_dmnsn) = nodeVal(1:N_node,1:N_dmnsn)
    node_typ(1:N_node)          = nodeT(1:N_node)
    CgXNode(1:N_node)           = CgXNval(1:N_node)
    CgYNode(1:N_node)           = CgYNval(1:N_node)
    
    
    do i = 1,N_node
       write(*,*) node_xy(i,1:N_dmnsn),i,node_typ(i),CgXNode(i),CgYNode(i),"aft read"
    enddo
    
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
  end subroutine read_nodeXY_and_CgVals_frmAFile
  
  
  
  subroutine read_nodeXY_and_CgVals_frm_PrePulleyToPostPulley(NumChar,flnm_arr)
    implicit none
    integer,            intent(in) :: NumChar
    character(len=200), intent(in) :: flnm_arr(1:NumChar)
    
    real*8  :: CgXVal,CgXNval(1:N_node)
    real*8  :: CgYVal,CgYNval(1:N_node)
    integer :: nodeTval,nodeT(1:N_node)
    real*8  :: nodeRead(1:N_dmnsn),nodeVal(1:N_node,1:N_dmnsn)
    integer :: cnt,error
    integer :: nodeCnt,i
    integer :: node1,node2
    
    if (NumChar==1) then
       write(*,*)trim(adjustl(flnm_arr(1))),"flnm_arr In read_nodeXY 1 pp"
       open(unit=449,file=trim(adjustl(flnm_arr(1))))
       
    elseif (NumChar==2) then
       write(*,*)trim(adjustl(flnm_arr(1))),"flnm_arr In read_nodeXY 21 pp"
       write(*,*)trim(adjustl(flnm_arr(2))),"flnm_arr In read_nodeXY 22 pp"
       open(unit=449,file=trim(adjustl(flnm_arr(1)))//trim(adjustl(flnm_arr(2))))
    endif
    
    cnt=1 ; nodeCnt=1
    nodeVal = -100.0d0
    
    do   
       read(unit=449,fmt=*,IOSTAT=error) nodeRead(1:N_dmnsn),nodeCnt,nodeTval,CgXval,CgYval
       write(*,*) nodeRead(1:N_dmnsn),nodeCnt,nodeTval,CgXval,CgYval,"p"
       
       if (error.lt.0) exit
       
       nodeVal(nodeCnt,1:N_dmnsn) = nodeRead(1:N_dmnsn)
       nodeT(nodeCnt)             = nodeTval 
       CgXNval(nodeCnt)           = CgXval
       CgYNval(nodeCnt)           = CgYval
       cnt = cnt+1
       
    enddo
    
    close(449)
    
    
    if (N_doubleNode==0) then
       continue
    elseif (N_doubleNode.gt.0) then
       
       do i = 1,N_doubleNode
          
          if (i==1) then
             node1            = double_node(i,1) ; node2            = double_node(i,2)
             nodeT(node1)     = 1                ; nodeT(node2)     = 1
             nodeVal(node1,1) = 0.0000d0         ; nodeVal(node1,2) = -0.0500d0
             nodeVal(node2,1) = nodeVal(node1,1) ; nodeVal(node2,2) = nodeVal(node1,2)
             
             CgXNval(node1)   = CgXNval(node1)   + CgXNval(node2) ! adding as they're connected now
             CgYNval(node1)   = CgYNval(node1)   + CgYNval(node2)
             CgXNval(node2)   = CgXNval(node1)
             CgYNval(node2)   = CgYNval(node1)
          else
             node1 = double_node(i,1) ; node2 = double_node(i,2);write(*,*)node1,node2,"doubleNode"
             nodeVal(node2,1:N_dmnsn) = nodeVal(node1,1:N_dmnsn)
             nodeT(node2)             = nodeT(node1)
             CgXNval(node2)           = CgXNval(node1)
             CgYNval(node2)           = CgYNval(node1)
          endif
          
       enddo
       
    endif
    
    if (nodeCnt.ne.N_node) then
       write(*,*) "nodeCnt MUST BE equal to N_node"
       write(*,*) nodeCnt,N_node,"N_node"
       stop
    endif
    
    node_xy(1:N_node,1:N_dmnsn) = nodeVal(1:N_node,1:N_dmnsn)
    node_typ(1:N_node)          = nodeT(1:N_node)
    CgXNode(1:N_node)           = CgXNval(1:N_node)
    CgYNode(1:N_node)           = CgYNval(1:N_node)
    
    
    do i = 1,N_node
       write(*,*) node_xy(i,1:N_dmnsn),i,node_typ(i),CgXNode(i),CgYNode(i),"aft read post pulley"
    enddo
    
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
  end subroutine read_nodeXY_and_CgVals_frm_PrePulleyToPostPulley
  
  
  subroutine read_kspr_and_l0vals_frmAFile(NumChar,flnm_arr)
    implicit none
    integer, intent(in)            :: NumChar
    character(len=200), intent(in) :: flnm_arr(1:NumChar)
    
    integer :: cntSpr,N_blck,cntLp
    integer :: max_blck
    integer :: i,readTheTnsn
    
    real*8 :: Tnsn(1:N_spr)
    real*8 :: Node1Val(1:N_dmnsn),Node2Val(1:N_dmnsn)
    real*8 :: TnsnV,ksV,lV,l0V
    real*8 :: TINY=1.0d-16
    
    real*8  :: x2,y2
    integer :: readALGO
    
    readALGO = 1 ! 1 with readTheTnsn and 0 without readTheTnsn
    max_blck = 6 ! CHANGE IT IF SYSTEM VARIES (IMPORTANT)
    call find_total_blck(N_blck)
    
    if ((N_blck-N_spr) .gt. max_blck) then
       write(*,*) "Extra Block MUST NOT that high",N_blck,N_spr,"N_blck,N_spr"
    endif
    
    
    if (NumChar==1) then
       write(*,*)trim(adjustl(flnm_arr(1))),"Name In read_kspr 1"
       open(unit=559,file=trim(adjustl(flnm_arr(1))))
       
    elseif (NumChar==2) then
       write(*,*)trim(adjustl(flnm_arr(2))),"Name In read_kspr 2"
       open(unit=559,file=trim(adjustl(flnm_arr(1)))//trim(adjustl(flnm_arr(2))))
    endif
    
    
    cntSpr=0
    
    write(*,*) N_blck,N_spr,"N_lck and N_spr"
    
    do i=1,N_blck
       
       if (readALGO == 0) then
          
          read(559,*) TnsnV,ksV,lV,l0V
          read(559,*) Node1Val(1:N_dmnsn)
          read(559,*) Node2Val(1:N_dmnsn)
          
          if (abs(TnsnV) .le. TINY) then
             continue
          else
             
             if (abs(Node2Val(1)).le.TINY .and. abs(Node2Val(2)).le.TINY) then
                continue
             else
                
                cntSpr=cntSpr+1
                
                Tnsn(cntSpr)   = TnsnV
                k_spr(cntSpr)  = ksV
                l(cntSpr)      = lV
                l0(cntSpr)     = l0V
                
                write(*,*) Tnsn(cntSpr),k_spr(cntSpr),l(cntSpr),l0(cntSpr),cntSpr,"sprPrp aft Read Algo 0"
                
             endif
             
          endif
          
       elseif (readALGO == 1) then
          
          read(559,*) TnsnV,ksV,lV,l0V,readTheTnsn
          read(559,*) Node1Val(1:N_dmnsn)
          read(559,*) Node2Val(1:N_dmnsn)
          
          write(*,*) TnsnV,ksV,lV,l0V,readTheTnsn
          
          if (readTheTnsn==1) then
             
             cntSpr=cntSpr+1
             
             Tnsn(cntSpr)   = TnsnV
             k_spr(cntSpr)  = ksV
             l(cntSpr)      = lV
             l0(cntSpr)     = l0V
             
             write(*,*) Tnsn(cntSpr),k_spr(cntSpr),l(cntSpr),l0(cntSpr),cntSpr,"sprPrp aft read Algo 1"
          
          endif
          
       else
          write(*,*) "ALGO is neither 0 nor 1"
          stop
       endif
       
    enddo
    
    close(559)
    
  end subroutine read_kspr_and_l0vals_frmAFile
  
  
  subroutine find_total_blck(sgmntCnt)
    use transfrm_info
    
    implicit none
    integer :: i,sgmntCnt
    
    sgmntCnt = 0 
    
    do i = 1,N_spr
       sgmntCnt = sgmntCnt + spr_node(i,0) - 1
    enddo
    
    sgmntCnt = sgmntCnt+2
    
  end subroutine find_total_blck
  
  
  subroutine read_karea_and_A0vals_frmAFile(NumChar,flnm_arr)
    implicit none
    integer, intent(in)            :: NumChar
    character(len=200), intent(in) :: flnm_arr(1:NumChar)
    
    real*8  :: PresV,kaV,AV,A0V,x1,y1
    real*8  :: Pres(1:N_cell)
    integer :: endDcsn,i
    
    
    if (NumChar==1) then
       write(*,*)trim(adjustl(flnm_arr(1)))," flnm_arr in read_karea 1"
       open(unit=561,file=trim(adjustl(flnm_arr(1))))
       
    elseif (NumChar==2) then
       write(*,*)trim(adjustl(flnm_arr(2)))," flnm_arr in read_karea 2"
       open(unit=561,file=trim(adjustl(flnm_arr(1)))//trim(adjustl(flnm_arr(2))))
    endif
    
    do i = 1,N_cell
       
       read(561,*) PresV,kaV,AV,A0V
       
       do
          read(561,*) x1,y1,endDcsn
          if (endDcsn==1) exit
       enddo
       
       Pres(i)   = PresV
       k_area(i) = kaV
       A(i)      = AV
       A0(i)     = A0V
       
       write(*,*) Pres(i),k_area(i),A(i),A0(i),i,"areaPrp aft read"
       
    enddo
    
    
    close(561)
    
  end subroutine read_karea_and_A0vals_frmAFile
  
  
end module storing_changing_restoring_routines




!##########################################################################


subroutine writing_nodes_with_number_and_type(dum_nodes,test_no,subtest_no,lp_cnt)
  use system_parameters
  implicit none
  
  real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
  integer :: test_no,subtest_no
  integer :: lp_cnt
  
  integer :: i
  character(len=100) :: flplce
  character(len=100) :: flnm
  character(len=100) :: test_name,subtest_name
  
  
  write(flplce,*) "/home/rniloy/PROJECT_CELL/2D_codes/unit_blck/restructuring_Data_Structure/"
  
  if (test_no .eq. 1) then
     write(test_name,*)"transition_test_NW_" !NW=nodewise
     write(flnm,'(a20,i3.3,a)')test_name,lp_cnt,".dat"
  elseif (test_no .eq. 2) then
     write(test_name,*)"No_partial_fixed_node_NW_"
     write(flnm,'(a26,i3.3,a)')test_name,lp_cnt,".dat"
  elseif (test_no .eq. 3) then
     write(test_name,*)"Effct_of_grvtnl_potntl_NW_"
     write(flnm,'(a27,i3.3,a)')test_name,lp_cnt,".dat"
     
  elseif (test_no .eq. 5) then

     write(test_name,*)"Twelve_cells_NW_"
     
     if (subtest_no .eq. 1) then
        write(subtest_name,*)"St1_"
     elseif (subtest_no .eq. 2) then
        write(subtest_name,*)"St2_"
     endif
     
     write(flnm,'(a21,i2.2,a)')trim(test_name)//adjustl(trim(subtest_name)),lp_cnt,".dat"
     
  elseif (test_no .eq. 6) then
     
     write(test_name,*)"Intrmd_cells_NW_"
     
     if (subtest_no .eq. 1) then
        write(subtest_name,*)"St1_"
     elseif (subtest_no .eq. 2) then
        write(subtest_name,*)"St2_"
     endif
     
     write(flnm,'(a21,i2.2,a)')trim(test_name)//adjustl(trim(subtest_name)),lp_cnt,".dat"   

  elseif (test_no .eq. 7) then
     
     write(test_name,*)"SubIntrmd_cells_NW_"
     
     if (subtest_no .eq. 1) then
        write(subtest_name,*)"St1_"
     elseif (subtest_no .eq. 2) then
        write(subtest_name,*)"St2_"
     endif
     
     write(flnm,'(a24,i2.2,a)')trim(test_name)//adjustl(trim(subtest_name)),lp_cnt,".dat"     
     
  endif

  
  !open(unit=32,file=trim(adjustl(flplce))//adjustl(flnm))
  open(unit=32,file=trim(adjustl(directryFlnm))//adjustl(flnm))
  
  do i = 1,N_node
     
     if (node_cnnctd(i) .eq. 0) then
        write(unit=32,fmt=*) dum_nodes(i,1:N_dmnsn), i, node_typ(i)
     elseif (node_cnnctd(i) .ne. 0) then
        if (count_this_dn(i) .eq. 1) then
           write(unit=32,fmt=*) dum_nodes(i,1:N_dmnsn), i, node_typ(i)
        endif
     endif
     
  enddo
  
  close(32)
  
end subroutine writing_nodes_with_number_and_type

subroutine writing_nodes_for_Exps(dum_Nodes,ExpNo,frmeNo)
  use system_parameters
  implicit none
  
  real*8,  intent(in) :: dum_Nodes(1:N_node,1:N_dmnsn)
  integer, intent(in) :: ExpNo,frmeNo 
  
  character(len=100) :: flplce
  character(len=100) :: flnm
  
  integer :: i
  integer :: objTyp
  
  objTyp = 1 !1 for Nodes
  
  call manage_ExpsName(ExpNo,frmeNo,objTyp,flnm)
  
  write(flplce,*) "/home/rniloy/PROJECT_CELL/2D_codes/unit_blck/restructuring_Data_Structure/"
  
  write(*,*) trim(adjustl(directryFlnm))
  write(*,*) adjustl(flnm),"1"
  
  !open(unit=32,file=trim(adjustl(flplce))//adjustl(flnm))
  open(unit=32,file=trim(adjustl(directryFlnm))//adjustl(flnm))
  
  do i = 1,N_node
     
     if (node_cnnctd(i) .eq. 0) then
        write(32,*) dum_Nodes(i,1:N_dmnsn), i, node_typ(i),CgXNode(i),CgYNode(i)
        
     elseif (node_cnnctd(i) .ne. 0) then
        if (count_this_dn(i) .eq. 1) then
           write(32,*) dum_Nodes(i,1:N_dmnsn), i, node_typ(i),CgXNode(i),CgYNode(i)
        endif
     endif
     
     !if (ExpNo==37) then
     !   if (i==1) then
     !      write(*,*) dum_Nodes(1,1:2),"FrmNo=",frmeNo
     !      write(*,*) node_xy(1,1:2),"FrmNo=",frmeNo
     !   endif
     !endif
     
  enddo
  
  
  close(32)
  
end subroutine writing_nodes_for_Exps


subroutine writing_nodes_over_spr(dum_Nodes,test_no,subtest_no,lp_cnt) !,test_no,N_lp,lp_var)
  
  use system_parameters
  use transfrm_info
  !use loop_info to replace (test_no,N_lp,lp_var)
  
  implicit none
  real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
  integer :: test_no,subtest_no
  integer :: lp_cnt
  !integer :: N_lp
  !integer :: lp_var(1:N_lp)
  
  integer :: i
  integer :: Pulley_type
  
  character(len=100) :: flplce
  character(len=100) :: flnm
  character(len=100) :: test_name,subtest_name

  integer :: node1,node2,node3
  
  !if (lp_cnt .eq. 1) then
     !write(*,*) dum_Nodes,"in writing"
  !endif
  
  Pulley_type = 6
  
  write(flplce,*) "/home/rniloy/PROJECT_CELL/2D_codes/unit_blck/restructuring_Data_Structure/"
  
  if (test_no .eq. 1) then
     write(test_name,*)"transition_test_SW_" !SW=spring wise
     write(flnm,'(a20,i3.3,a)')test_name,lp_cnt,".dat"
  elseif (test_no .eq. 2) then
     write(test_name,*)"No_partial_fixed_node_SW_"
     write(flnm,'(a26,i3.3,a)')test_name,lp_cnt,".dat"
  elseif (test_no .eq. 3) then
     write(test_name,*)"Effct_of_grvtnl_potntl_SW_"
     write(flnm,'(a27,i3.3,a)')test_name,lp_cnt,".dat"
     
  elseif (test_no .eq. 5) then
     write(test_name,*)"Twelve_cells_SW_"

     if (subtest_no .eq. 1) then
        write(subtest_name,*)"St1_"
     elseif (subtest_no .eq. 2) then
        write(subtest_name,*)"St2_"
     endif
     
     write(flnm,'(a21,i2.2,a)')trim(test_name)//adjustl(trim(subtest_name)),lp_cnt,".dat"  

  elseif (test_no .eq. 6) then
     write(test_name,*)"Intrmd_cells_SW_"
     
     if (subtest_no .eq. 1) then
        write(subtest_name,*)"St1_"
     elseif (subtest_no .eq. 2) then
        write(subtest_name,*)"St2_"
     endif
     
     write(flnm,'(a21,i2.2,a)')trim(test_name)//adjustl(trim(subtest_name)),lp_cnt,".dat" 

  elseif (test_no .eq. 7) then
     write(test_name,*)"SubIntrmd_cells_SW_"
     
     if (subtest_no .eq. 1) then
        write(subtest_name,*)"St1_"
     elseif (subtest_no .eq. 2) then
        write(subtest_name,*)"St2_"
     endif
     
     write(flnm,'(a24,i2.2,a)')trim(test_name)//adjustl(trim(subtest_name)),lp_cnt,".dat" 

  endif
  
  
  !open(unit=31,file=trim(adjustl(flplce))//adjustl(flnm))
  open(unit=31,file=trim(adjustl(directryFlnm))//adjustl(flnm))
  
  !write(*,*) N_spr,"N_spr"
  !write(*,*) spr_node(1,1:2),"spr_node"
  
  do i = 0,(N_spr+1)
     if (i.eq.0 .or. i.eq.(N_spr+1)) then
        if (i.eq.0) then
           node1 = 1
           node2 = 2
           !write(unit=31,fmt=*) dum_nodes(node1,1:2),dum_nodes(node2,1:2)
           
        elseif (i.eq.(N_spr+1)) then
           write(*,*) "WARNING :: will be changed with SHAPE,sbrtn: writing_nodes_over_spr"
           node1 = lft_endNode(2) + 1
           node2 = lft_endNode(2) + 2
           !node1 = N_node - 1
           !node2 = N_node
           !write(unit=31,fmt=*) dum_nodes(node1,1:2),dum_nodes(node2,1:2)
        endif
        
        write(unit=31,fmt=*) dum_nodes(node1,1:2),dum_nodes(node2,1:2)
        
     elseif (i.gt.0 .and. i.lt.(N_spr+1)) then
        
        !if (typ_spr(i) .ne. Pulley_type) then
        if (spr_node(i,0)==2) then
           node1 = spr_node(i,1)
           node2 = spr_node(i,2)
           write(unit=31,fmt=*) dum_nodes(node1,1:2),dum_nodes(node2,1:2)
           
        !elseif (typ_spr(i) .eq. Pulley_type) then
        elseif (spr_node(i,0)==3) then
           node1 = spr_node(i,1)
           node2 = spr_node(i,2)
           node3 = spr_node(i,3)
           write(unit=31,fmt=*) dum_nodes(node1,1:2),dum_nodes(node2,1:2)
           write(unit=31,fmt=*) dum_nodes(node2,1:2),dum_nodes(node3,1:2)       
        endif
        
     endif
     
  enddo
  
  close(31)
  
end subroutine writing_nodes_over_spr


subroutine writing_nodes_over_spr_single_node(dum_Nodes,Tension,ExpNo,frmeNo)
  use system_parameters
  use transfrm_info
  use conversion_routines
  use Adding_cells
  
  implicit none
  real*8, intent(in) :: dum_Nodes(1:N_node,1:N_dmnsn)
  real*8, intent(in) :: Tension(1:N_spr)
  integer,intent(in) :: ExpNo,frmeNo
  
  character(len=100) :: flnm
  character(len=100) :: flplce
  
  integer :: i1,j1,j1max
  integer :: prtn,prtnMax
  integer :: strtNode,fnshNode
  integer :: spr_nm,node_nm
  integer :: objTyp
  integer :: j1St,j1En
  
  real*8  :: TnsnWall,k_sprWall,lwall,l0wall
  integer :: readTheTnsn
  
  objTyp = 2 !2 for sprs
  
  call manage_ExpsName(ExpNo,frmeNo,objTyp,flnm)
  
  write(flplce,*) "/home/rniloy/PROJECT_CELL/2D_codes/unit_blck/restructuring_Data_Structure/"
  
  !open(unit=31,file=trim(adjustl(flplce))//trim(adjustl(flnm)))
  open(unit=31,file=trim(adjustl(directryFlnm))//adjustl(flnm))
  
  TnsnWall = 0.00d0 ; k_sprWall = 0.00d0 ; lwall = 0.00d0 ; l0wall = 0.00d0
  
  do i1 = 0,(N_spr+1)
     
     if (i1==0) then
        
        readTheTnsn = 0
        strtNode    = 1 ; fnshNode = 2
        write(31,*) TnsnWall,k_sprWall,lwall,l0wall,readTheTnsn
        write(31,*) dum_Nodes(strtNode,1:N_dmnsn)
        write(31,*) dum_Nodes(fnshNode,1:N_dmnsn)
        write(31,*) " "
        
     elseif (i1.ne.0 .and. i1.ne.(N_spr+1)) then
        
        spr_nm  = i1
        j1max   = spr_node(spr_nm,0)
        
        if (j1max==2) then
           readTheTnsn = 1
           write(31,*) Tension(spr_nm),k_spr(spr_nm),l(spr_nm),l0(spr_nm),readTheTnsn
           
           do j1 = 1,j1max
              node_nm = spr_node(spr_nm,j1)
              write(31,*) dum_Nodes(node_nm,1:N_dmnsn)
           enddo
           
           write(31,*) " "
           
        elseif (j1max.gt.2) then
           
           prtnMax = j1Max-1
           
           do prtn = 1,prtnMax
              
              if (prtn.ne.prtnMax) readTheTnsn=0
              if (prtn == prtnMax) readTheTnsn=1
              
              write(31,*) Tension(spr_nm),k_spr(spr_nm),l(spr_nm),l0(spr_nm),readTheTnsn
              
              if (prtn==1) then
                 j1St = 1
                 j1En = 2
              else
                 j1St = j1St+1
                 j1En = j1En+1
              endif
              
              do j1 = j1St,j1En
                 node_nm = spr_node(spr_nm,j1)
                 write(31,*) dum_Nodes(node_nm,1:N_dmnsn)
              enddo
              
              write(31,*) " "
           enddo
           
        endif
        
     elseif (i1==(N_spr+1)) then
        
        if (AddedCellModelInitiation==0) then
           strtNode = ((2*nvsl)+1)    ; fnshNode = ((2*nvsl)+2)
           !write(*,*) strtNode,fnshNode,nvsl,N_TNsAC,"strtfnshA1"
           
        elseif (AddedCellModelInitiation==1) then
           strtNode = (N_TNsAC/2)+1   ; fnshNode = (N_TNsAC/2)+2
           !write(*,*) strtNode,fnshNode,nvsl,N_TNsAC,"strtfnshA2"
        endif
        
        readTheTnsn = 0
        write(31,*) TnsnWall,k_sprWall,lwall,l0wall,readTheTnsn
        write(31,*) dum_Nodes(strtNode,1:N_dmnsn)
        write(31,*) dum_Nodes(fnshNode,1:N_dmnsn)
        write(31,*) " "
        
        write(*,*) strtNode,fnshNode,nvsl,N_TNsAC,"strtfnshB"
        
     endif
     
  enddo
  
  close(31)
  
end subroutine writing_nodes_over_spr_single_node



subroutine writing_nodes_over_area(dum_Nodes,test_no,subtest_no,lp_cnt)
  use system_parameters
  use transfrm_info

  implicit none
  real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
  integer :: test_no,subtest_no
  integer :: lp_cnt

  integer :: i

  character(len=100) :: flplce
  character(len=100) :: flnm
  character(len=100) :: test_name,subtest_name
  character(len=100) :: fmt
  
  integer :: node(1:max_node_area)
  integer :: node_in_the_area
  
  write(flplce,*) "/home/rniloy/PROJECT_CELL/2D_codes/unit_blck/restructuring_Data_Structure/"
  
  if (test_no .eq. 5) then
     write(test_name,*)"Twelve_cells_AW_" !AW = Area wise
     
     if (subtest_no .eq. 1) then
        write(subtest_name,*)"St1_"
     elseif (subtest_no .eq. 2) then
        write(subtest_name,*)"St2_"
     endif
     
     write(flnm,'(a21,i2.2,a)') trim(test_name)//adjustl(trim(subtest_name)),lp_cnt,".dat"

  elseif (test_no.eq.6) then
     write(test_name,*)"Intrmd_cells_AW_"

     if (subtest_no .eq. 1) then
        write(subtest_name,*)"St1_"
     elseif (subtest_no .eq. 2) then
        write(subtest_name,*)"St2_"
     endif
     
     write(flnm,'(a21,i2.2,a)') trim(test_name)//adjustl(trim(subtest_name)),lp_cnt,".dat"

  elseif (test_no.eq.7) then
     write(test_name,*)"SubIntrmd_cells_AW_"
     
     if (subtest_no .eq. 1) then
        write(subtest_name,*)"St1_"
     elseif (subtest_no .eq. 2) then
        write(subtest_name,*)"St2_"
     endif
     
     write(flnm,'(a24,i2.2,a)') trim(test_name)//adjustl(trim(subtest_name)),lp_cnt,".dat"

  endif
  
  !open(unit=33,file=trim(adjustl(flplce))//adjustl(flnm))
  open(unit=33,file=trim(adjustl(directryFlnm))//adjustl(flnm))
  
  node(1:5) = 0
  
  do i = 1,N_cell
     
     node_in_the_area = area_node(i,0)
     
     node(1) = area_node(i,1)
     node(2) = area_node(i,2)
     node(3) = area_node(i,3)
     node(4) = area_node(i,4)
     
     if (typ_area(i) .eq. 1) then
        continue
     elseif (typ_area(i) .eq. 2) then
        node(5) = area_node(i,5)
     endif

     write(fmt,'(I4,"(f5.2)")') node_in_the_area
     !write(*,*) adjustl(fmt),"FMT"
     !do j = 1,node_in_the_area
     !write(unit=33,'(5(f5.2))') (dum_nodes(node(j)),j=1,node_in_the_area)

     if (typ_area(i) .eq. 1) then
        write(unit=33,fmt='(4(f5.2))') dum_nodes(node(1),1:2),dum_nodes(node(2),1:2),dum_nodes(node(3),1:2),&
             dum_nodes(node(4),1:2)

     elseif (typ_area(i) .eq. 2) then
        write(unit=33,fmt='(5(f5.2))') dum_nodes(node(1),1:2),dum_nodes(node(2),1:2),dum_nodes(node(3),1:2),&
              dum_nodes(node(4),1:2),dum_nodes(node(5),1:2)
     endif
     !write(unit=31,'(5(f5.2))') (node(j),j=1,node_in_the_are
        !if (j.eq.node_in_the_area) write(unit=31,fmt=*) " "
     !enddo
     
  enddo

  close(33)
  
end subroutine writing_nodes_over_area


subroutine writing_nodes_over_area_single_node(dum_Nodes,Cell_Pressure,ExpNo,frmeNo)
  use system_parameters
  use transfrm_info
  use conversion_routines
  
  implicit none
  real*8, intent(in) :: dum_Nodes(1:N_node,1:N_dmnsn)
  real*8, intent(in) :: Cell_Pressure(1:N_cell)
  integer,intent(in) :: ExpNo,frmeNo
  
  character(len=100) :: flplce,flnm1,flnm2
  
  integer :: i1,j1,j1max
  integer :: area_nm,node_nm
  integer :: objTyp
  integer :: zeroI,oneI
  
  zeroI = 0 ; oneI = 1
  
  objTyp = 3
  call manage_ExpsName(ExpNo,frmeNo,objTyp,flnm1)
  objTyp = 5
  call manage_ExpsName(ExpNo,frmeNo,objTyp,flnm2)
  
  write(flplce,*) "/home/rniloy/PROJECT_CELL/2D_codes/unit_blck/restructuring_Data_Structure/"
  
  !open(unit=31,file=trim(adjustl(flplce))//adjustl(flnm))
  open(unit=31,file=trim(adjustl(directryFlnm))//adjustl(flnm1))
  open(unit=41,file=trim(adjustl(directryFlnm))//adjustl(flnm2))
  
  do i1 = 1,N_cell
     write(unit=31,fmt=*) Cell_Pressure(i1),k_area(i1),A(i1),A0(i1)
     write(unit=41,fmt=*) Cell_Pressure(i1),k_area(i1),A(i1),A0(i1)
     
     area_nm = i1   
     j1max = area_node(i1,0)
     
     do j1 = 1,j1max
        node_nm = area_node(i1,j1)
        if (j1.ne.j1max) write(unit=31,fmt=*) dum_Nodes(node_nm,1:N_dmnsn)
        if (j1.ne.j1max) write(unit=41,fmt=*) dum_Nodes(node_nm,1:N_dmnsn),zeroI
        if (j1 == j1max) write(unit=31,fmt=*) dum_Nodes(node_nm,1:N_dmnsn)
        if (j1 == j1max) write(unit=41,fmt=*) dum_Nodes(node_nm,1:N_dmnsn),oneI
     enddo

     write(unit=31,fmt=*) " "
     write(unit=41,fmt=*) " "
  enddo
  
  close(31)
  close(41)
  
end subroutine writing_nodes_over_area_single_node



subroutine writing_datas_for_curve(dum_Nodes,Cell_Pressure,Tension,ExpNo,frmeNo)
  use system_parameters
  use neighbour_info_and_trnsfrmInfo_dependent_info
  use transfrm_info
  use conversion_routines
  
  implicit none
  real*8, intent(in) :: dum_Nodes(1:N_node,1:N_dmnsn)
  real*8, intent(in) :: Cell_Pressure(1:N_cell)
  real*8, intent(in) :: Tension(1:N_spr)
  integer,intent(in) :: ExpNo,frmeNo
  
  character(len=100) :: flplce,flnm
  
  integer :: i,j,jmax
  integer :: area_nm,node_nm
  integer :: objTyp
  
  objTyp = 4
  
  call manage_ExpsName(ExpNo,frmeNo,objTyp,flnm)
  
  write(flplce,*) "/home/rniloy/PROJECT_CELL/2D_codes/unit_blck/restructuring_Data_Structure/"
  
  !open(unit=31,file=trim(adjustl(flplce))//adjustl(flnm))
  open(unit=31,file=trim(adjustl(directryFlnm))//adjustl(flnm))
  
  call get_curveVars(dum_Nodes,Cell_Pressure,Tension)
  
  do i = 1,N_curve
     write(31,fmt=*) P1(i),P2(i),dP(i),TC(i),curveX1Y1(i,1:2),curveX2Y2(i,1:2)
     write(31,fmt=*) xc_crv(i),yc_crv(i),R(i),AnglS(i),AnglE(i)
     write(unit=31,fmt=*) " "
  enddo
  
  
  close(31)
  
end subroutine writing_datas_for_curve



subroutine get_delX_and_delY(dum_Nodes,dum_Nodes_cmpr,delX,delY)
  
  use system_parameters
  
  implicit none
  real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
  real*8  :: dum_Nodes_cmpr(1:N_node,1:N_dmnsn) !cmpr = to be compared
  real*8  :: delX(1:N_node),delY(1:N_node)
  
  integer :: i
  
  delX(1:N_node) = 0.0d0
  delY(1:N_node) = 0.0d0
  
  !write(*,*) node_typ(1:N_node),"node_typ"
  
  do i = 1,N_node
     !write(*,*) i,"Node_number"
     !write(*,*) dum_Nodes(i,1:N_dmnsn),"dum_Nodes"
     !write(*,*) dum_Nodes_cmpr(i,1:N_dmnsn),"dum_Nodes_cmpr"
     
     if (node_typ(i) .eq. 0) then
        continue
     elseif (node_typ(i) .eq. 1) then
        delX(i) = dum_Nodes(i,1) - dum_Nodes_cmpr(i,1)
        delY(i) = dum_Nodes(i,2) - dum_Nodes_cmpr(i,2)
     elseif (node_typ(i) .eq. 2) then
        !write(*,*) dum_Nodes(i,2),dum_Nodes_cmpr(i,2),"dum_Nodes and cmpr"
        delX(i) = dum_Nodes(i,1) - dum_Nodes_cmpr(i,1)
     endif
     !write(*,*) delX(i),delY(i),"delX-delY"
  enddo


end subroutine get_delX_and_delY


subroutine get_delX_and_delY_pairwise(nodes_pw,nodes_cmpr_pw,delX_pw,delY_pw,hlf_node)
  use system_parameters
  
  implicit none
  integer :: hlf_node
  real*8  :: nodes_pw(1:hlf_node,1:N_dmnsn),nodes_cmpr_pw(1:hlf_node,1:N_dmnsn)
  real*8  :: delX_pw(1:hlf_Node),delY_pw(1:hlf_node)
 
  integer :: i

  delX_pw(1:hlf_Node) = 0.0d0
  delY_pw(1:hlf_Node) = 0.0d0
  
  do i = 1,hlf_Node
     delX_pw(i) = nodes_pw(i,1) - nodes_cmpr_pw(i,1)
     delY_pw(i) = nodes_pw(i,2) - nodes_cmpr_pw(i,2)
     write(*,*) nodes_cmpr_pw(i,1:2),"cmpr_pw"
     write(*,*) nodes_pw(i,1:2),"pw"
  enddo

  
end subroutine get_delX_and_delY_pairwise





subroutine writing_delX_and_delY(delX,delY,test_no,lp_cnt,nonD_delA0_area)
  use system_parameters
  
  implicit none
  real*8  :: delX(1:N_node),delY(1:N_node)
  integer :: test_no
  integer :: lp_cnt

  real*8, intent(in) :: nonD_delA0_area
  
  integer :: i

  character(len=100) :: flplce
  character(len=100) :: flnm
  character(len=100) :: test_name
  
  
  write(flplce,*) "/home/rniloy/PROJECT_CELL/2D_codes/unit_blck/restructuring_Data_Structure/"
  
  if (test_no .eq. 1) then
     write(test_name,*)"delX_delY_for_TT_"
     write(flnm,'(a18,i3.3,a)')test_name,lp_cnt,".dat"
     
  elseif (test_no .eq. 2) then
     write(test_name,*)"delX_delY_for_TT_NPFN_"
     write(flnm,'(a23,i3.3,a)')test_name,lp_cnt,".dat"
     
  elseif (test_no .eq. 3) then
     write(test_name,*)"delX_delY_for_grvtnl_potntl_"
     write(flnm,'(a29,i3.3,a)')test_name,lp_cnt,".dat"
  endif
  
  
  !open(unit=32,file=trim(adjustl(flplce))//adjustl(flnm))
  open(unit=32,file=trim(adjustl(directryFlnm))//adjustl(flnm))
  
  do i = 1,N_node
     write(unit=32,fmt=*) delX(i),delY(i),nonD_delA0_area,i
  enddo
  
  close(32)
  
end subroutine writing_delX_and_delY



subroutine create_nonD_delA0_area_vs_coordnte_data(N_iter,test_no)
   !this subroutine will create a nondimensionalized del_A0_area versus del_coordnte data to plot
  use system_parameters
  
  implicit none
  integer :: N_iter
  integer :: test_no
  
  integer :: nlines = 0
  !integer :: node
  real*8  :: dx,dy
  real*8  :: nonD_delA0_Area
  
  character(len=100) :: flnm,created_flnm
  character(len=100) :: flnmbr,created_flnmbr
  character(len=100) :: full_flnm,created_full_flnm
  
  integer :: i,j
  
  created_flnm='nonD_delA0_area_vs_coordnte'
  
  do i = 1,N_node
  
     write(created_flnmbr,'(i3.3,a)') i,'.dat'
     created_full_flnm=trim(adjustl(created_flnm))//trim(adjustl(created_flnmbr))
     
     open(unit=11,file=trim(adjustl(created_full_flnm)))
     
     if (test_no .eq. 1) flnm='delX_delY_for_TT_'
     if (test_no .eq. 2) flnm='delX_delY_for_TT_NPFN_'
     
     do j = 1,N_iter
        
        write(flnmbr,'(i3.3,a)') j,'.dat'
        
        full_flnm=trim(adjustl(flnm))//trim(adjustl(flnmbr))
        
        open(unit=3,file=trim(adjustl(full_flnm)))
        
        do
           read(3,fmt=*) dx,dy,nonD_delA0_Area
           nlines = nlines+1
           
           if (nlines .eq. i) then
              write(unit=11,fmt=*) dx,dy,nonD_delA0_area
              exit
           endif
           
        enddo
        
        close(3)
        
        nlines = 0
        
     enddo
     
     close(11)
  enddo
  
end subroutine create_nonD_delA0_area_vs_coordnte_data




subroutine writing_delX_delY_and_CgX(delX,delY,test_no,lp_cnt,Cg_dum)
  use system_parameters
  
  implicit none
  real*8  :: delX(1:N_node),delY(1:N_node)
  integer :: test_no
  integer :: lp_cnt
  
  !real*8, intent(in) :: nonD_g
  real*8, intent(in) :: Cg_dum
  
  integer :: i
  
  character(len=100) :: flplce
  character(len=100) :: flnm
  character(len=100) :: test_name
  
  !write(*,*) Cg_dum,"Cg_dum"
  
  write(flplce,*) "/home/rniloy/PROJECT_CELL/2D_codes/unit_blck/restructuring_Data_Structure/"
  !work for test_03
  if (test_no .ne. 3) then
     write(*,*) "only for test 3 compatible, sbrtn name: writing_delX_delY_CgX"
     stop
  endif
  
  
  write(test_name,*)"delX_delY_for_grvtnl_potntl_"
  write(flnm,'(a29,i3.3,a)')test_name,lp_cnt,".dat" 

  
  !open(unit=32,file=trim(adjustl(flplce))//adjustl(flnm))
  open(unit=32,file=trim(adjustl(directryFlnm))//adjustl(flnm))
  
  do i = 1,N_node
     !write(*,*) delX(i),delY(i),Cg_dum,i,"delx,dely,Cg_dum,i"
     write(unit=32,fmt=*) delX(i),delY(i),Cg_dum,i
  enddo
  
  close(32)
  
end subroutine writing_delX_delY_and_CgX



subroutine writing_delX_delY_and_CgX_pairwise(delX_pw,delY_pw,test_no,lp_cnt,Cg_dum,hlf_node)
  use system_parameters

  implicit none
  integer :: hlf_node
  real*8  :: delX_pw(1:hlf_node),delY_pw(1:hlf_node)
  integer :: test_no
  integer :: lp_cnt
  
  !real*8, intent(in) :: nonD_g
  real*8, intent(in) :: Cg_dum
  
  integer :: i
  
  character(len=100) :: flplce
  character(len=100) :: flnm
  character(len=100) :: test_name
  
  !write(*,*) Cg_dum,"Cg_dum"
  
  write(flplce,*) "/home/rniloy/PROJECT_CELL/2D_codes/unit_blck/restructuring_Data_Structure/"
  !work for test_03
  if (test_no .ne. 3) then
     write(*,*) "only for test 3 compatible, sbrtn name: writing_delX_delY_CgX"
     stop
  endif
  
  
  write(test_name,*)"delX_delY_for_grvtnl_potntl_"
  write(flnm,'(a29,i3.3,a)')test_name,lp_cnt,".dat" 
  
  
  !open(unit=32,file=trim(adjustl(flplce))//adjustl(flnm))
  open(unit=32,file=trim(adjustl(directryFlnm))//adjustl(flnm))
  
  do i = 1,hlf_node
     !write(*,*) delX(i),delY(i),Cg_dum,i,"delx,dely,Cg_dum,i"
     write(unit=32,fmt=*) delX_pw(i),delY_pw(i),Cg_dum,i
     write(*,*) delX_pw(i),delY_pw(i),Cg_dum,i
  enddo
  
  close(32)
  
end subroutine writing_delX_delY_and_CgX_pairwise





subroutine create_CgX_vs_coordnte_data(N_iter,test_no)
  use system_parameters
  
  implicit none
  integer :: N_iter
  integer :: test_no
  
  integer :: nlines = 0
  !integer :: node
  real*8  :: dx,dy
  !real*8  :: nonD_g
  
  character(len=100) :: flnm,created_flnm
  character(len=100) :: flnmbr,created_flnmbr
  character(len=100) :: full_flnm,created_full_flnm
  integer :: i,j
  integer :: reason
  
  created_flnm='CgX_vs_coordnte'
  
  i=0
  
  do
     i = i+1
     write(created_flnmbr,'(i3.3,a)') i,'.dat'
     created_full_flnm=trim(adjustl(created_flnm))//trim(adjustl(created_flnmbr))
     
     open(unit=11,file=trim(adjustl(created_full_flnm)))
     
     if (test_no .eq. 1) flnm ='delX_delY_for_TT_'
     if (test_no .eq. 2) flnm ='delX_delY_for_TT_NPFN_'
     if (test_no .eq. 3) flnm ='delX_delY_for_grvtnl_potntl_'
     
     do j = 1,(N_iter)
        !write(*,*) j,"j"
        write(flnmbr,'(i3.3,a)') (j-1),'.dat'
        
        full_flnm=trim(adjustl(flnm))//trim(adjustl(flnmbr))
        
        open(unit=3,file=trim(adjustl(full_flnm)))
        
        do
           read(3,fmt=*,IOSTAT=reason) dx,dy,CgX
           !write(*,*) reason,"reason"
           nlines = nlines+1
              
           if (nlines .eq. i) then
              !write(*,*) dx,dy,CgX,"dx-dy-CgX"
              !write(*,*) "E1"
              !write(*,*) reason,"reason"
              if (reason .gt. 0) then
                 write(*,*) "something wrong in delX_delY dat, sbrtn: create_CgX_vs_coordnte_data"
                 stop
              elseif (reason .lt. 0) then
                 !write(*,*) reason,"reason should be -ve"
                 exit
              elseif (reason .eq. 0) then
                 write(unit=11,fmt=*) dx,dy,CgX
                 exit
              endif
              !write(*,*) "After_exit"
           endif
           !write(*,*) "After_exit"
        enddo
        
        
        !write(*,*) "Aft_loop"
        close(3)
        
        nlines = 0
        
     enddo
     
     if (reason .lt. 0) then
        close(11,status='delete')
     else
        close(11)
     endif
     
     if (reason .lt. 0) exit
     
  enddo
  
end subroutine create_CgX_vs_coordnte_data
  

subroutine fret_fpp_data(fret,fpp,iter,test_no,CgX,dA0)
  implicit none
  integer, intent(in) :: test_no
  real*8,  intent(in) :: fret,fpp
  integer, intent(in) :: iter
  real*8 , intent(in) :: CgX
  real*8 , intent(in) :: dA0
  
  character(len=100) :: flplce
  character(len=100) :: flnm,flnmbr,full_flnm
  character(len=100) :: test_name
  
  !write(*,*) Cg_dum,"Cg_dum"
  
  write(flplce,*) "/home/rniloy/PROJECT_CELL/2D_codes/unit_blck/restructuring_Data_Structure/"
  ! mainly work for test_03, but flexibility is maintained to integrate for example, if in case of combining test_02 and test_03
  
  if (test_no .eq. 1) then
     write(test_name,*)"fret_fpp_for_TT_"
     write(flnm,'(a17)')test_name
     
  elseif (test_no .eq. 2) then
     write(test_name,*)"fret_fpp_for_TT_NPFN_"
     write(flnm,'(a22)')test_name
     
  elseif (test_no .eq. 3) then
     write(test_name,*)"fret_fpp_for_grvtnl_potntl_"
     write(flnm,'(a28)')test_name
  endif
  
  write(flnmbr,'(a,f4.2,a)') 'dA0=',dA0,'.dat'
  
  
  
  full_flnm=trim(adjustl(flnm))//trim(adjustl(flnmbr))
  
  open(unit=32,file=trim(adjustl(full_flnm)),position='append')
  
  write(unit=32,fmt=*) CgX,fret,fpp,(fret-fpp),iter
  
  close(32)
  
  
end subroutine fret_fpp_data


subroutine get_analytical_delY_vs_delA0(prcnt_chng)
  use system_parameters
  implicit none
  
  real*8  :: prcnt_chng
  integer :: hlf_Node=1
  integer :: N_stbl !No of Stable config = 4

  !shape1 is corner point up
  !shape2 is corner point down
  !shape3 is striaght up
  !shape4 is straight down
  
  real*8, allocatable :: delY_PwAna(:,:) !Pairwise Analytical
 
  
  if (mod(N_node,2).eq.0) then
     hlf_Node = N_node/2
  elseif (mod(N_node,2).eq.1) then
     write(*,*) "Not compatible with odd N_node, sbrtn: get_analytical_delY_vs_delA0"
     stop
  endif

  N_stbl = 4
  allocate(delY_PwAna(1:N_stbl,1:hlf_Node))
  
  call get_delYpwAna

contains

  subroutine get_delYpwAna
    implicit none
    real*8 :: area_nd
    real*8 :: len
    real*8 :: phi_x !angle with x-axis

    integer :: i
    
    area_nd = 1+prcnt_chng
    len = l0(1) !1/2/3 anyone will work here as of equal length
    
    if (area_nd.gt.1.0d0) then
       write(*,*) "Cant be greater than 1.0; sbrtn : get_for_shape1"
       stop
    endif
    
    phi_x = acos(area_nd)

    open(unit=122,file='delY_PwAna.dat',position='append')
    
    do i = 1,N_stbl
    
       delY_PwAna(i,1) = 0.0d0
       
       if (i.eq.1 .or. i.eq.3) then
          delY_PwAna(i,2) = len*sin(phi_x)
       elseif (i.eq.2 .or. i.eq.4) then
          delY_PwAna(i,2) = -len*sin(phi_x)
       endif
       
       if (i.le.2) then
          delY_PwAna(i,3) = 0.0d0
       elseif (i.eq.3) then
          delY_PwAna(i,3) = 2.0*len*sin(phi_x)
       elseif (i.eq.4) then
          delY_PwAna(i,3) = -2.0*len*sin(phi_x)
       endif
       
       write(unit=122,fmt=*) prcnt_chng,delY_PwAna(i,1:3),i
       
    enddo
    
    close(122)
  end subroutine get_delYpwAna

  
end subroutine get_analytical_delY_vs_delA0


subroutine get_experimntl_delY_vs_delA0(prcnt_chng)
  use system_parameters
  implicit none
  real*8 :: prcnt_chng
  character(len=200) :: flplce,flplce1
  character(len=100) :: flnm,created_flnm
  !character(len=100) :: full_flnm,created_full_flnm
  
  integer :: nlines
  real*8  :: dx,dy,CgValue
  integer :: pw_node
  real*8  :: dy1,dy2,dy3
  
  integer :: i
  integer :: directry_id
  
  !file for writing
  write(flplce,*) "/home/rniloy/PROJECT_CELL/2D_codes/unit_blck/restructuring_Data_Structure/grv_vs_crdnte_diff_dA0/"

  created_flnm='delY_PwExp.dat'
  
  !open(unit=11,file=trim(adjustl(flplce))//adjustl(created_flnm),position='append')
  open(unit=11,file=trim(adjustl(directryFlnm))//adjustl(created_flnm),position='append')
  
  !file for reading
  directry_id = 100 + int(100.00d0*prcnt_chng)
  write(*,*) directry_id

  if(directry_id .lt. 100) then
     !write(flplce1,'(a106,a2,i2.2,a)')flplce,"A_",directry_id,"_prcnt/"
     write(flplce1,'(a106,a2,i2.2,a)')directryFlnm,"A_",directry_id,"_prcnt/"
     
  elseif(directry_id .ge. 100) then
     !write(flplce1,'(a106,a2,i3.3,a)')flplce,"A_",directry_id,"_prcnt/"
     write(flplce1,'(a106,a2,i3.3,a)')directryFlnm,"A_",directry_id,"_prcnt/"
  endif
  
  flnm ='delX_delY_for_grvtnl_potntl_000.dat'      
  
  nlines = 0
  
  do i=1,3
     open(unit=3,file=trim(adjustl(flplce1))//adjustl(flnm))  

     do
        read(3,fmt=*) dx,dy,CgValue,pw_node
        nlines = nlines+1
        
        if (nlines .eq. i) then
     
           if (i.eq.1) dy1=dy
           if (i.eq.2) dy2=dy
           if (i.eq.3) dy3=dy
           exit
        endif
        
     enddo
     close(3)
     nlines = 0
  enddo

  i=3 !stbl_condn 3
  write(unit=11,fmt=*) prcnt_chng,dy1,-dy2,dy3,i
  i=4 !stbl_condn 4
  write(unit=11,fmt=*) prcnt_chng,dy1,dy2,dy3,i
  
  close(11)

  
end subroutine get_experimntl_delY_vs_delA0

subroutine save_l0_A0_nodeXY(updt_cnt)
  use system_parameters
  use conversion_routines
  
  implicit none
  integer,intent(in) :: updt_cnt
  character(len=100) :: flnm
  
  integer :: i,j,jmax
  integer :: N_item

  
  call coordntes_to_nodes(coordntes,node_xy)
  
  write(flnm,'(a14,i2.2,a)')"system_params",updt_cnt,".dat"
  flnm=adjustl(trim(flnm))
  write(*,*) flnm,"filename"
  
  N_item = 3
  jmax = -1
  
  open(unit=114,file=flnm)
  
  do i = 1,N_item
     if (i.eq.1) jmax=N_node
     if (i.eq.2) jmax=N_spr
     if (i.eq.3) jmax=N_cell
     
     do j = 1,jmax
        if (i==1) write(unit=114,fmt=*) node_xy(j,1:2),j
        if (i==2) write(unit=114,fmt=*) l(j),l0(j),j
        if (i==3) write(unit=114,fmt=*) A(j),A0(j),j
        
        if (j.eq.jmax) write(unit=114,fmt=*) " "
     enddo
     
  enddo
  
  close(114)
  
end subroutine save_l0_A0_nodeXY


subroutine use_saved_l0_A0_nodeXY(updt_cnt)
  use system_parameters
  use conversion_routines

  implicit none
  integer, intent(in) :: updt_cnt

  integer :: i,j,jmax
  character(len=100) :: flnm

  integer :: N_item,ilines
  integer :: ioVal,n_skip
  
  write(flnm,'(a14,i2.2,a)')"system_params",updt_cnt,".dat"
  flnm=adjustl(trim(flnm))
  write(*,*) flnm,"filename"
  
 

  N_item = 3
  jmax   = -1
  ioVal  = 0
  n_skip = 1
  
  open(unit=124,file=flnm,action='read')

  do i = 1,N_item
     if (i.eq.1) jmax=N_node
     if (i.eq.2) jmax=N_spr
     if (i.eq.3) jmax=N_cell
     
     do j = 1,jmax
        if (ioVal > 0)  then
           write(*,*) "WARNING in READING"
        elseif (ioval < 0) then
           exit
        else
           continue
        endif
        
        
        if (i==1) read(unit=124,fmt=*,iostat=ioVal) node_xy(j,1:2)
        if (i==2) read(unit=124,fmt=*,iostat=ioVal) l(j),l0(j)
        if (i==3) read(unit=124,fmt=*,iostat=ioVal) A(j),A0(j)
        
        if (j.eq.jmax) then
           do ilines = 1,n_skip
              read(124,'()',advance='yes')
           enddo
        endif
        
     enddo
     
  enddo

  close(124)

  call nodes_to_coordntes(node_xy,coordntes)
  
end subroutine use_saved_l0_A0_nodeXY


subroutine manage_Filename(N_info,test_info,object_info,flnm)
  
  implicit none
  integer, intent(in) :: N_info
  integer, intent(in) :: object_info
  
  integer, intent(in) :: test_info(1:N_info)

  character(len=100), intent(out) :: flnm

  integer :: i
  character(len=100) :: test_name,subtest_name,subsubtest_name
  
  ! do i = 1,N_info
     
  !    if (i==1) then
  !       test_no       = test_info(i)
  !    elseif (i==2) then
  !       subtest_no    = test_info(i)
  !    elseif (i==3) then
  !       subsubtest_no = test_info(i)
  !    endif
     
  ! enddo

  ! if (test_no==5) then
  !    test_name =
  !    if (subtest_no==1) then
  !       subtest_name =
  !    elseif (subtest_no==2) then
  !       subtest_name =
  !    endif
     
  ! elseif (test_no==6) then
  !    test_name =
     
  ! elseif (test_no==7) then
  !    test_name =
  ! endif
  
end subroutine manage_Filename


subroutine manage_ExpsName(ExpNo,frmeNo,objTyp,flnm)
  implicit none
  
  integer,intent(in)             :: ExpNo,objTyp,frmeNo
  character(len=100),intent(out) :: flnm
  character(len=100)             :: mainName,suffxName
  
  open(unit=56,file='NameCheck.dat',position='append')
  
  if (ExpNo==0) then
     write(mainName,*) "Exp0_aftOptmzToIncrseCg"
  elseif (ExpNo==1) then
     write(mainName,*)"Exp1_botCgIncr"
  elseif (ExpNo==2) then
     write(mainName,*)"Exp2_l0ShortningSecndMidSpr"
  elseif (ExpNo==3) then
     write(mainName,*)"Exp3_l0ShorteningPulleySpr"
  elseif (ExpNo==4) then
     write(mainName,*)"Exp4_l0ShorteningApicalSprNxtToPulley"
  elseif (ExpNo==5) then
     write(mainName,*)"Exp5_l0ShorteningDiagonalSpr"
  elseif (ExpNo==6) then
     write(mainName,*)"Exp6_l0ShorteningBasalFurrowSpr"
  elseif (ExpNo==7) then
     write(mainName,*)"Exp7_l0ShorteningThrdMidSpr"
  elseif (ExpNo==8) then
     write(mainName,*)"Exp8_l0ShorteningVertSprBfrPulley"
  elseif (ExpNo==9) then
     write(mainName,*)"Exp9_l0ShorteningThrdSideSpr"
  elseif (ExpNo==10) then
     write(mainName,*)"Exp10_l0ShorteningSecndSideSpr"
  elseif (ExpNo==11) then
     write(mainName,*)"Exp11_l0ShorteningThrdConnSpr"
  elseif (ExpNo==12) then
     write(mainName,*)"Exp12_l0ShorteningBasalSprBfrPulley"
  elseif (ExpNo==13) then
     write(mainName,*)"Exp0_atCg0_A0_incr_thrdPairfrmBottm"
  elseif (ExpNo==14) then
     write(mainName,*)"Exp0_Cg0_bot6_A0incr_strct1_"
  elseif (ExpNo==15) then
     write(mainName,*)"Exp0_Cg0_bot6_A0incr_strct2_"
  elseif (ExpNo==16) then
     write(mainName,*)"Trnsfrmtn"
  elseif (ExpNo==17) then
     write(mainName,*)"Trnsfrmtn_threeStrct"
  elseif (ExpNo==18) then
     write(mainName,*)"tst_Exp0_Cg0_bot6_A0incr_strct1_"
  elseif (ExpNo==19) then
     write(mainName,*)"tst_Exp0_Cg0_bot6_A0incr_strct2_"
  elseif (ExpNo==20) then
     write(mainName,*)"NI_modelPrgrsn"
  elseif (ExpNo==21) then
     write(mainName,*)"stg1E_"
  elseif (ExpNo==22) then
     write(mainName,*)"Smth_model"
  elseif (ExpNo==23) then
     write(mainName,*)"BeginStg_"
  elseif (ExpNo==24) then
     write(mainName,*)"S4T1_strct1_"
  elseif (ExpNo==25) then
     write(mainName,*)"S4T1_strct2_"
  elseif (ExpNo==26) then
     write(mainName,*)"PrgrsnStg_"
  elseif (ExpNo==27) then
     write(mainName,*)"BeginStgAddCell_"
  elseif (ExpNo==28) then
     write(mainName,*)"AddedCellProgrsnStg_"
  elseif (ExpNo==29) then
     write(mainName,*)"ProgrsnStgMovemnt_"
  elseif (ExpNo==30) then
     write(mainName,*)"NI_model"
  elseif (ExpNo==31) then
     write(mainName,*)"PrgrsnStgMovementSeparately_"
  elseif (ExpNo==32) then
     write(mainname,*)"NI_modelAddedCell"
     
  elseif (ExpNo==33) then ! very genralizd naming(33-36) for multpl cases
     write(mainName,*)"Case01"
  elseif (ExpNo==34) then
     write(mainName,*)"Case02"
  elseif (ExpNo==35) then
     write(mainName,*)"Case03"
  elseif (ExpNo==36) then
     write(mainName,*)"Case04"
  elseif (ExpNo==37) then
     write(mainName,*)"NI_modelInitiatn"
  elseif (ExpNo==38) then 
     write(mainName,*)"NI_modelInitiatn_WT_VF"
  elseif (ExpNo==39) then
     !write(mainName,*)""
     continue
  elseif (ExpNo==40) then
     write(mainName,*)"ProductnRunInitiatn"
  endif
  
  !write(*,*) mainName
  
  if (objTyp==1) then
     write(suffxName,*)"NW"
  elseif (objTyp==2) then
     write(suffxName,*)"SW"
  elseif (objTyp==3) then
     write(suffxName,*)"AW"
  elseif (objTyp==4) then
     write(suffxName,*)"CW"
  elseif (objTyp==5) then
     write(suffxName,*)"AP" !AP=Area+Position
  else
     write(*,*) "Not an Object,flnme:writing_sbrtns,sbrtn:manag_ExpName"
     stop
  endif
  
  !write(suffxName,*) trim(adjustl(suffxName))
  !write(*,*) trim(adjustl(suffxName))
  
  !write(flnm,'(a,i3.3,a)')(trim(adjustl(mainName)))//(trim(adjustl(suffxName))),frmeNo,".dat"
  write(flnm,'(a,i4.4,a)')(trim(adjustl(mainName)))//(trim(adjustl(suffxName))),frmeNo,".dat"
  
  write(56,*) trim(adjustl(flnm))
  write(*,*) trim(adjustl(flnm))
  close(56)
  
  !call sleep(0.5)
  
end subroutine manage_ExpsName


subroutine manage_ExpsName_withlen200(ExpNo,frmeNo,objTyp,flnm)
  
  implicit none
  
  integer,intent(in) :: ExpNo,objTyp,frmeNo
  character(len=200),intent(out) :: flnm
  character(len=200) :: mainName,suffxName
  
  open(unit=56,file='NameCheck.dat',position='append')
  
  if (ExpNo==0) then
     write(mainName,*) "Exp0_aftOptmzToIncrseCg"
  elseif (ExpNo==1) then
     write(mainName,*)"Exp1_botCgIncr"
  elseif (ExpNo==2) then
     write(mainName,*)"Exp2_l0ShortningSecndMidSpr"
  elseif (ExpNo==3) then
     write(mainName,*)"Exp3_l0ShorteningPulleySpr"
  elseif (ExpNo==4) then
     write(mainName,*)"Exp4_l0ShorteningApicalSprNxtToPulley"
  elseif (ExpNo==5) then
     write(mainName,*)"Exp5_l0ShorteningDiagonalSpr"
  elseif (ExpNo==6) then
     write(mainName,*)"Exp6_l0ShorteningBasalFurrowSpr"
  elseif (ExpNo==7) then
     write(mainName,*)"Exp7_l0ShorteningThrdMidSpr"
  elseif (ExpNo==8) then
     write(mainName,*)"Exp8_l0ShorteningVertSprBfrPulley"
  elseif (ExpNo==9) then
     write(mainName,*)"Exp9_l0ShorteningThrdSideSpr"
  elseif (ExpNo==10) then
     write(mainName,*)"Exp10_l0ShorteningSecndSideSpr"
  elseif (ExpNo==11) then
     write(mainName,*)"Exp11_l0ShorteningThrdConnSpr"
  elseif (ExpNo==12) then
     write(mainName,*)"Exp12_l0ShorteningBasalSprBfrPulley"
  elseif (ExpNo==13) then
     write(mainName,*)"Exp0_atCg0_A0_incr_thrdPairfrmBottm"
  elseif (ExpNo==14) then
     write(mainName,*)"Exp0_Cg0_bot6_A0incr_strct1_"
  elseif (ExpNo==15) then
     write(mainName,*)"Exp0_Cg0_bot6_A0incr_strct2_"
  elseif (ExpNo==16) then
     write(mainName,*)"Trnsfrmtn"
  elseif (ExpNo==17) then
     write(mainName,*)"Trnsfrmtn_threeStrct"
  elseif (ExpNo==18) then
     write(mainName,*)"tst_Exp0_Cg0_bot6_A0incr_strct1_"
  elseif (ExpNo==19) then
     write(mainName,*)"tst_Exp0_Cg0_bot6_A0incr_strct2_"
  elseif (ExpNo==20) then
     write(mainName,*)"NI_modelPrgrsn"
  elseif (ExpNo==21) then
     write(mainName,*)"stg1E_"
  elseif (ExpNo==22) then
     write(mainName,*)"Smth_model"
  elseif (ExpNo==23) then
     write(mainName,*)"BeginStg_"
  elseif (ExpNo==24) then
     write(mainName,*)"S4T1_strct1_"
  elseif (ExpNo==25) then
     write(mainName,*)"S4T1_strct2_"
  elseif (ExpNo==26) then
     write(mainName,*)"PrgrsnStg_"
  elseif (ExpNo==27) then
     write(mainName,*)"BeginStgAddCell_"
  elseif (ExpNo==28) then
     write(mainName,*)"AddedCellProgrsnStg_"
  elseif (ExpNo==29) then
     write(mainName,*)"ProgrsnStgMovemnt_"
  elseif (ExpNo==30) then
     write(mainName,*)"NI_model"
  elseif (ExpNo==31) then
     write(mainName,*)"PrgrsnStgMovementSeparately_"
  elseif (ExpNo==32) then
     write(mainname,*)"NI_modelAddedCell"
     
  elseif (ExpNo==33) then ! very genralizd naming(33-36) for multpl caswes
     write(mainName,*)"Case01"
  elseif (ExpNo==34) then
     write(mainName,*)"Case02"
  elseif (ExpNo==35) then
     write(mainName,*)"Case03"
  elseif (ExpNo==36) then
     write(mainName,*)"Case04"
  elseif (ExpNo==37) then
     write(mainName,*)"NI_modelInitiatn"
  elseif (ExpNo==38) then 
     write(mainName,*)"NI_modelInitiatn_WT_VF"
  elseif (ExpNo==39) then
     !write(mainName,*)""
     continue
  elseif (ExpNo==40) then
     write(mainName,*)"ProductnRunInitiatn"
  endif
  
  !write(*,*) mainName
  
  if (objTyp==1) then
     write(suffxName,*)"NW"
  elseif (objTyp==2) then
     write(suffxName,*)"SW"
  elseif (objTyp==3) then
     write(suffxName,*)"AW"
  elseif (objTyp==4) then
     write(suffxName,*)"CW"
  elseif (objTyp==5) then
     write(suffxName,*)"AP" !AP=Area+Position
  else
     write(*,*) "Not an Object,flnme:writing_sbrtns,sbrtn:manag_ExpName"
     stop
  endif
  
  !write(suffxName,*) trim(adjustl(suffxName))
  !write(*,*) trim(adjustl(suffxName))
  
  !write(flnm,'(a,i3.3,a)')(trim(adjustl(mainName)))//(trim(adjustl(suffxName))),frmeNo,".dat"
  write(flnm,'(a,i4.4,a)')(trim(adjustl(mainName)))//(trim(adjustl(suffxName))),frmeNo,".dat"
  
  write(56,*) trim(adjustl(flnm))
  write(*,*) trim(adjustl(flnm))
  close(56)
  
  !call sleep(0.5)
  
end subroutine manage_ExpsName_withlen200
