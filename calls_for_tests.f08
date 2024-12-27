module calls_for_tests
  
  use sys_building_info !input_frm_bash_what_to_build.f08
  use cell_info
  use end_info
  use system_parameters 
  
  use generating_the_shape !generating_system.f08
  use perturbing_the_shape !generating_system.f08
  use changing_parameters  !Manipulating_properties.f08
  
  use node_variables
  use spring_variables
  use area_variables
  use grvtnl_variables
  use bend_variables
  use SRyp_variables
  use curve_variables
  
  use spr_to_other_transfrm
  use node_to_other_transfrm
  use area_to_other_transfrm
  use neighbour_info_and_trnsfrmInfo_dependent_info
  use transfrm_info
  
  use conversion_routines
  
  use Energy_module
  use gradient_module
  use PenF_module
  use grdPenF_module
  use PenF_minimizer_mod
  use Wfunc_and_its_derivative
  
  use mltiple_calls_together
  use moving_coordnte_info
  use moving_coordnte_variables
  
  use force_calc  
  use storing_changing_restoring_routines
  use Manipulate_and_equilibrate
  use switch_and_unswitch_models
  
  use redefining_system_module
  use stiff_bottom4cell_holding_it_in_prev_pos
  
  implicit none
  integer :: test_no,subtest_no
  integer :: lp_cnt
  real*8  :: ftol
  
  integer :: Exprmnt,Frame
  logical :: Hold_lgclRdfn

  
contains
  
  subroutine test_01
    !test_01 = In test_01, I will make rest area A0_area from less than La*Lb (l0^2) to greater than and will see if I can get a plot that will create a plot from unstable equilibrium to stable equillibrium through the critical point
    
    implicit none
    real*8  :: prcnt_chng
    integer :: i
    integer :: N_iter
    integer :: lp_strt,lp_end
    integer :: iter
    real*8  :: fret,fpp

    call get_all_moving_coordnte_variables!conversion_for_moving_system.f08  
    call get_all_gradient_variables !gradient_calc.f08
    call get_gradient(node_xy,l0,A0,gradient) !gradient_calc.f08
    
    test_no = 1 ; lp_cnt = 0
    N_iter  = 11
    
    lp_strt = -int(N_iter/2)
    lp_end  = int(N_iter/2)
    
    !write(*,*) lp_strt,lp_end,"lp_strt,lp_end"

    !open(unit=21,file='A0_value_check.dat')
    !open(unit=22,file='see_results.dat')
    
    do i = lp_strt,lp_end
       lp_cnt = lp_cnt + 1
       call store_node_xy
       
       !write(unit=21,fmt=*) i,"lp_no"
       prcnt_chng = (0.1d0*i)
       call store_and_change_A0_area(prcnt_chng)
       !write(unit=21,fmt=*) A0,"A0 aft chnage"
       
       if (i.ne.0) then
          funcChoice = 1 ; ftol=1.d-08 
          call frprmn(coordntes_xy,grdmv_xy,iter,funcChoice,N_mvCoordnte,ftol,fret,fpp)
       endif
       write(*,*) lp_cnt,"lp_cnt_aft_frprmn"
       !write(unit=22,fmt=*) coordntes_xy,"aft_frprmn"
       !write(unit=22,fmt=*) fret,"aft_frprmn"
       
       call generate_delX_and_delY(coordntes_xy,test_no,lp_cnt,N_iter)

!UNCOMMENT       !call save_config_and_generate_data(coordntes_xy,test_no,lp_cnt)
       
       call restore_node_xy
       call restore_A0_area
       !write(unit=21,fmt=*) A0,"A0 aft restoring"  
      
    enddo
    
    call deallocate_moving_coordnte_variables
    call create_nonD_delA0_area_vs_coordnte_data(N_iter,test_no)
    
    !close(21)
    !close(22)
    !call system ('gnuplot transition_test.gnu')

  end subroutine test_01
  
  
  subroutine test_02
    !test_02 = In test_02, I will make rest area A0_area = La*Lb (l0^2) and will see if I can get a stable equilibrium with minimal force at the moving coordinates (critical point)
    
    implicit none
    real*8  :: prcnt_chng
    integer :: i
    integer :: N_iter
    integer :: lp_strt,lp_end
    integer :: iter
    real*8  :: fret,fpp
    
    integer :: which_type,changed_type


    which_type   = 2
    changed_type = 1
    
    call store_and_change_node_type(which_type,changed_type)
    write(*,*) node_typ(3),node_typ(5),"node_typ-3,5"
    
    call get_all_moving_coordnte_variables!conversion_for_moving_system.f08  
    call get_all_gradient_variables !gradient_calc.f08
    call get_gradient(node_xy,l0,A0,gradient) !gradient_calc.f08
    
    
    test_no = 2 ; lp_cnt = 0
    N_iter  = 11
    
    lp_strt = -int(N_iter/2)
    lp_end  = int(N_iter/2)
    
    write(*,*) lp_strt,lp_end,"lp_strt,lp_end"
    
    !open(unit=21,file='A0_value_check.dat') !co1
    !open(unit=22,file='see_results.dat') !co2
    
    !open(unit=41,file='inside_loop.dat') !co3
    
   

    do i = lp_strt,lp_end
       lp_cnt = lp_cnt + 1
       !write(unit=41, fmt=*) lp_cnt,"lp_cnt" !co4
       !write(unit=41, fmt=*) node_xy,"node_xy" !co5
       call store_node_xy
       
       !write(unit=41,fmt=*) i,"lp_no" !co6
       prcnt_chng = (0.1d0*i)
       write(*,*) prcnt_chng,"prcnt_chng"
       call store_and_change_A0_area(prcnt_chng)

       !write(unit=21,fmt=*) A0,"A0 aft chnage" !co7
       
       if (i.ne.0) then
          !call nodes_to_coordntes(node_xy,coordntes_xy) !!!!!!!!!!!DO I NEED THIS?
          funcChoice = 1 ; ftol=1.d-08  
          call frprmn(coordntes_xy,grdmv_xy,iter,funcChoice,N_mvCoordnte,ftol,fret,fpp)
       endif
       
       write(*,*) lp_cnt,"lp_cnt_aft_frprmn"
       
       !write(unit=22,fmt=*) coordntes_xy,"aft_frprmn" !co8
       !write(unit=22,fmt=*) fret,"aft_frprmn" !co9
       
       call generate_delX_and_delY(coordntes_xy,test_no,lp_cnt,N_iter)
!UNCOMMENT       call save_config_and_generate_data(coordntes_xy,test_no,lp_cnt)
       
       call restore_node_xy
       call restore_A0_area
       
       
       !write(unit=21,fmt=*) A0,"A0 aft restoring" !co10  
       
    enddo
    
    call restore_node_type
    call deallocate_moving_coordnte_variables
    call create_nonD_delA0_area_vs_coordnte_data(N_iter,test_no)
    
    !close(21)
    !close(22)
    !call system ('gnuplot transition_test.gnu')
    
    
    
  end subroutine test_02
  
  
  
  subroutine test_03_PosdA0
    !test_03 = In test_03, I will take one specific value of A0(say A0=0.7*La*Lb); then I will check the gravitational potential for different CgXNode value.    
    implicit none
    
    
    real*8  :: prcnt_chng
    integer :: i
    integer :: N_iter
    !integer :: lp_strt,lp_end
    integer :: iter
    real*8  :: fret,fpp
    
    integer :: which_type,changed_type
    !real*8  :: CgX_incr
    
    integer :: node_nmbr
    real*8  :: dx,dy
    
    logical :: spr,area,grv,bend,sryp
    real*8  :: Enrgy_frst
    
    spr=.True. ; area=.True. ; grv=.True. ; bend=.False. ; sryp=.False. 
    !CgX = 0.30d0    
    call get_sys_lgcls(spr,area,grv,bend,sryp)
    
    select_xy = 2 !change here to shift from horizontal to vertical
    write(*,*) select_xy,"WATCH : Check the select xy"
    
    which_type   = 2
    changed_type = 1
    
    call store_and_change_node_type(which_type,changed_type)
    write(*,*) node_typ(3),node_typ(5),"node_typ-3,5"
    !write(*,*) node_typ(3),node_typ(4),"node_typ-3,4"
    
    write(*,*) A0,"A0"
    write(*,*) k_spr,"k_spr"

    !do i = 1,2 !this commntd out blck is to change a node value
       !if (i.eq.1) node_nmbr = 5
       !if (i.eq.2) node_nmbr = 6
       node_nmbr = 5
       dx = 0.00d0
       dy = -0.05d0
       call change_coordnte_of_a_node(node_nmbr,dx,dy)
       write(*,*) node_xy(node_nmbr,1:2) ,"node_xy-5 or 6"
    !enddo
    
    !dx = 0.00d0
    !dy = 0.00d0
    
    !type_of_node = 1

    !call change_coordnte_of_a_type_of_node(type_of_node,dx,dy)
    !write(*,*) node_xy(3,1:2),node_xy(4,1:2),node_xy(5,1:2),node_xy(6,1:2),"node_xy frm 3 to 6"
    
    call get_all_moving_coordnte_variables!conversion_for_moving_system.f08  
    call get_all_gradient_variables !gradient_calc.f08

    call get_gradient(node_xy,l0,A0,gradient) !gradient_calc.f08
    
    !write(*,*) gradient,"grd"
    !write(*,*) " "

    !write(*,*) node_xy,"node_xy"

    Enrgy_frst = Energy(node_xy,l0,A0)
    write(*,*) Enrgy_frst,"Enrgy_frst"

    test_no = 3 ; lp_cnt = -1 
    
    call alloc_and_init_node_xy_cmpr_coordntes_zero_grv
   
    prcnt_chng = 0.10d0 ! 0.50 means A0 = 1.5 LaLb
 



    N_imp_region  = 5
    allocate(imp_region(1:N_imp_region))
    imp_region(1) = 1
    imp_region(2) = 2
    imp_region(3) = 3
    imp_region(4) = 4
    imp_region(5) = 5
    
    !do i = 1,2
    div_imp_region = 20
    
    call rgion_dvsn(N_imp_region,div_imp_region,imp_region)
    
    !enddo
    
    
    !deallocate(imp_region)
    
    
    
    !N_iter = N_CgXEval
    N_iter = 2 !watch out this and prev line
    
   
    ! if (abs(prcnt_chng-1.0d0) .lt. (1e-16)) then
    !    write(*,*) "prcnt_chng can not equal to 1.0d0 in test 3"
    !    stop
    ! endif
    write(*,*) N_iter,"N_iter"
    
    do i = 1,N_iter
       lp_cnt = lp_cnt + 1
       
       CgX = CgX_val(i)
       !write(*,*) CgX,"CgX"
       
       call store_node_xy
       !write(*,*) node_xy(5,1:2),node_xy(6,1:2),"node_xy_5-6"
       !write(*,*) prcnt_chng,"prcnt_chng"! t(-1)
       call store_and_change_A0_area(prcnt_chng) !t0
       
       fret = 0.0d0 ; fpp = 0.0d0
       grdmv_xy = 10e5
       
       !write(*,*) A0(1:2),"insd A0"
       !write(*,*) coordntes_xy,"coordntes_xy_bfr"
       !write(*,*) grdmv_xy,"grdmv_xy_bfr"
       !write(*,*) " "
       
       !if (i.eq.368) then    
       funcChoice = 1 ; ftol=1.d-08  
       call frprmn(coordntes_xy,grdmv_xy,iter,funcChoice,N_mvCoordnte,ftol,fret,fpp)
       !endif
       
       call fret_fpp_data(fret,fpp,iter,test_no,CgX,prcnt_chng)
       !write(*,*) fret,"fret_aft"
       
       !write(*,*) iter,"ITER took to Minimize"
       !write(*,*) coordntes_xy,"coordntes_xy_aft"
       !write(*,*) grdmv_xy,"grdmv_xy_aft"
       
       !write(*,*) lp_cnt,"lp_cnt_aft_frprmn"
       
       !write(unit=22,fmt=*) coordntes_xy,"aft_frprmn" 
       !write(unit=22,fmt=*) fret,"aft_frprmn" 
       
       
       if (i .eq. 1) then
          write(*,*) coordntes_xy,"at CgX=0.00"
          call CgX_zero_data(coordntes_xy,coordntes_zero_grv)
       endif   
       
       !elseif (i .ne. 1) then
          if (i==2) write(*,*) coordntes_xy,"at CgX=0.005"
          
          call generate_delX_delY_CgX(coordntes_xy,coordntes_zero_grv,test_no,lp_cnt,N_iter,CgX) !t1
          
!UNCOMMENT    call save_config_and_generate_data(coordntes_xy,test_no,lp_cnt) !t2
          
          
       !endif
        
       call restore_node_xy !t3
       call nodes_to_coordntes(node_xy,coordntes_xy) !t4
       call restore_A0_area !t5
       
       !write(unit=21,fmt=*) A0,"A0 aft restoring" !co10  
       
    enddo
    
    call restore_node_type
    call deallocate_moving_coordnte_variables
    call create_CgX_vs_coordnte_data(N_iter,test_no)
    
    !!close(21)
    !!close(22)
    !!call system ('gnuplot transition_test.gnu')
    
  end subroutine test_03_PosdA0
  
  
  subroutine test_03_NegdA0
    !test_03 = In test_03, I will take one specific value of A0(say A0=0.7*La*Lb); then I will check the gravitational potential for different CgX value.    
    implicit none
    
    real*8  :: prcnt_chng
    integer :: i
    integer :: N_iter,N_iter_updted
    !integer :: lp_strt,lp_end
    integer :: lp_incr
    integer :: iter
    real*8  :: fret,fpp
    
    integer :: which_type,changed_type
    !real*8  :: CgX_incr

    integer :: node_nmbr
    real*8  :: dx,dy
    
    logical :: spr,area,grv,bend,sryp
    real*8  :: Enrgy_frst
    
    spr=.True. ; area=.True. ; grv=.True. ; bend=.False. ; sryp=.False.
    call get_sys_lgcls(spr,area,grv,bend,sryp)
    
    select_xy = 2 !change here to shift from horizontal to vertical
    write(*,*) select_xy,"WATCH : Check the select xy"
    
    which_type   = 2
    changed_type = 1
    
    call store_and_change_node_type(which_type,changed_type)
    
    
    write(*,*) A0,"A0"
    write(*,*) k_spr,"k_spr"
  
    call get_all_moving_coordnte_variables!conversion_for_moving_system.f08
    call alloc_and_init_node_xy_cmpr_coordntes_zero_grv

    call get_all_gradient_variables !gradient_calc.f08
    call get_gradient(node_xy,l0,A0,gradient) !gradient_calc.f08
    
    !Enrgy_frst = Energy(node_xy,l0,A0)
    write(*,*) Enrgy_frst,"Enrgy_frst"

    
    prcnt_chng = -0.10d0 ! 0.50 means A0 = 1.5 LaLb
    call get_symmetric_config

    do i = 1,2
       if (i.eq.1) node_nmbr = 5
       if (i.eq.2) exit !node_nmbr = 5
       dx = 0.00d0
       dy = -0.05d0
    
       call change_coordnte_of_a_node(node_nmbr,dx,dy)
       write(*,*) node_xy(node_nmbr,1:2) ,"node_xy-5 or 6"
    enddo
    
    call nodes_to_coordntes(node_xy,coordntes_xy)

    test_no = 3 ; lp_cnt = -1 

    N_imp_region  = 3
    allocate(imp_region(1:N_imp_region))
    imp_region(1) = 1
    imp_region(2) = 2
    imp_region(3) = 3
    !imp_region(4) = 4
    !imp_region(5) = 5
    
    div_imp_region = 100
    call rgion_dvsn(N_imp_region,div_imp_region,imp_region)
    
    
    N_iter = N_CgXEval
    !N_iter = 10 !watch out this and prev line
    
    CgX_max = 0.1101d0 !This value will include 0.1100000000000001 also
    
    N_iter_updted = 0
    
    do i = 1,N_iter
       if (CgX_val(i) .gt. CgX_max) then
          N_iter_updted = (i-1)
          exit
       endif
    enddo
    
    N_iter = N_iter_updted
    ! if (abs(prcnt_chng-1.0d0) .lt. (1e-16)) then
    !    write(*,*) "prcnt_chng can not equal to 1.0d0 in test 3"
    !    stop
    ! endif
    !write(*,*) N_iter,"N_iter"
    call store_and_change_A0_area(prcnt_chng)

    lp_cnt  = 0
    lp_incr = -1
    
    do i = N_iter,1,lp_incr
       lp_cnt = i-1
       CgX = CgX_val(i)
       write(*,*) CgX,"CgX"
       
       if (i.eq.1) write(*,*) coordntes_xy,"bfr_at CgX=0.000"
       if (i.eq.2) write(*,*) coordntes_xy,"bfr_at CgX=0.005"
       if (i.eq.3) write(*,*) coordntes_xy,"bfr_at CgX=0.010"
       
       if (i.eq.N_iter) write(*,*) coordntes_xy,"bfr_at CgX(N_iter)"
       if (i.eq.(N_iter-1)) write(*,*) coordntes_xy,"bfr_at CgX(N_iter-1)"
       if (i.eq.(N_iter-2)) write(*,*) coordntes_xy,"bfr_at CgX(N_iter-2)"
       
       !CgX = 0.03
       !if (i.eq.1) CgX = 0.2475d0
       !if (i.eq.2) CgX = 0.2500d0
       !write(*,*) CgX,"CgX"
       
       !call store_node_xy
       !write(*,*) node_xy(5,1:2),node_xy(6,1:2),"node_xy_5-6"
       !write(*,*) prcnt_chng,"prcnt_chng"! t(-1)
       !call store_and_change_A0_area(prcnt_chng) !t0
       
       fret = 0.0d0 ; fpp = 0.0d0
       grdmv_xy = 10e5
       
       write(*,*) A0(1:2),"insd A0"
       !write(*,*) coordntes_xy,"coordntes_xy_bfr"
       !write(*,*) grdmv_xy,"grdmv_xy_bfr"
       !write(*,*) " "
       
       !if (i.eq.368) then    
          funcChoice = 1 ; ftol=1.d-08  
          call frprmn(coordntes_xy,grdmv_xy,iter,funcChoice,N_mvCoordnte,ftol,fret,fpp)
       !endif
          
       call fret_fpp_data(fret,fpp,iter,test_no,CgX,prcnt_chng)
       !write(*,*) fret,"fret_aft"
       
       !write(*,*) iter,"ITER took to Minimize"
       !write(*,*) coordntes_xy,"coordntes_xy_aft"
       !write(*,*) grdmv_xy,"grdmv_xy_aft"
       
       !write(*,*) lp_cnt,"lp_cnt_aft_frprmn"
       
       !write(unit=22,fmt=*) coordntes_xy,"aft_frprmn" 
       !write(unit=22,fmt=*) fret,"aft_frprmn"


       if (i.eq.1) write(*,*) coordntes_xy,"aft_at CgX=0.000"
       if (i.eq.2) write(*,*) coordntes_xy,"aft_at CgX=0.005"
       if (i.eq.3) write(*,*) coordntes_xy,"aft_at CgX=0.010"
       
       
       if (i.eq.(N_iter)) write(*,*) coordntes_xy,"aft_at CgX(N_iter)"
       if (i.eq.(N_iter-1)) write(*,*) coordntes_xy,"aft_at CgX(N_iter-1)"
       if (i.eq.(N_iter-2)) write(*,*) coordntes_xy,"aft_at CgX(N_iter-2)"
       
       
       call generate_delX_delY_CgX_pairwise(coordntes_xy,coordntes_zero_grv,test_no,lp_cnt,N_iter,CgX) !t1
!UNCOMMENT   call save_config_and_generate_data(coordntes_xy,test_no,lp_cnt) !t2

       !call restore_node_xy !t3
       !call nodes_to_coordntes(node_xy,coordntes_xy) !t4
       !call restore_A0_area !t5
       
       
       !write(unit=21,fmt=*) A0,"A0 aft restoring" !co10  
       
    enddo
    
    call restore_A0_area    
    call restore_node_type
    call deallocate_moving_coordnte_variables
    call create_CgX_vs_coordnte_data(N_iter,test_no)
    
    call get_analytical_delY_vs_delA0(prcnt_chng)
    call get_experimntl_delY_vs_delA0(prcnt_chng)
    
    !!close(21)
    !!close(22)
    !!call system ('gnuplot transition_test.gnu')


  contains

    subroutine get_symmetric_config

      implicit none

      !write(*,*) prcnt_chng,"inside get_symmetric"

      !CgX = 0.00d0
      
      !call store_node_xy
      !call store_and_change_A0_area(prcnt_chng) !t0

      !fret = 0.0d0 ; fpp = 0.0d0
      !grdmv_xy = 10e5
      
      !call nodes_to_coordntes(node_xy,coordntes_xy)
      !call frprmn(coordntes_xy,grdmv_xy,iter,fret,fpp,N_mvCoordnte)
      !call CgX_zero_data(coordntes_xy,coordntes_zero_grv)
      
      !write(*,*) coordntes_xy,"coordntes_xy inside symmetric"
      
      !call restore_node_xy
      !call nodes_to_coordntes(node_xy,coordntes_xy)
      !call restore_A0_area
      
      call nodes_to_coordntes(node_xy,coordntes_xy) !unperturbed_sys is taken as symmtric config and be fed in CgX_zero_data 
      call CgX_zero_data(coordntes_xy,coordntes_zero_grv)
      
    end subroutine get_symmetric_config
    
    
  end subroutine test_03_NegdA0
  
  
  
  
  
  subroutine test_04
    
    implicit none
    integer :: iter
    real*8  :: fret,fpp
    
    integer :: which_type,changed_type
    integer :: node_nmbr
    real*8  :: dx,dy
    logical :: spr,area,grv,bend,sryp
    
    real*8  :: Enrgy_frst
    
    integer :: i
    
    spr=.True. ; area=.False. ; grv=.True. ; bend=.False. ; sryp=.False.
    CgX = 0.3d0
    
    call get_sys_lgcls(spr,area,grv,bend,sryp)
    select_xy = 1
    
    which_type   = 2
    changed_type = 1
    
    call store_and_change_node_type(which_type,changed_type)
    
    write(*,*) node_typ(3),node_typ(5),"node_typ-3,5"
    !write(*,*) node_typ(3),node_typ(4),"node_typ-3,4"
    
    write(*,*) A0,"A0"
    write(*,*) k_spr,"k_spr"
    
    do i = 1,2
       if (i.eq.1) node_nmbr = 5
       if (i.eq.2) node_nmbr = 6
       
       dx = 0.05d0
       dy = 0.00d0
       call change_coordnte_of_a_node(node_nmbr,dx,dy)
       write(*,*) node_xy(node_nmbr,1:2) ,"node_xy-5 or 6"
    enddo
    
    call get_all_moving_coordnte_variables!conversion_for_moving_system.f08  
    call get_all_gradient_variables !gradient_calc.f08
    !write(*,*) node_xy
    !write(*,*) grd_spr,"grd_spr"
    !write(*,*) grd_area,"grd_area"
    !write(*,*) grd_grvtnl,"grd_grvtnl"
    
    call get_gradient(node_xy,l0,A0,gradient) !gradient_calc.f08
    write(*,*) gradient,"grd"
    write(*,*) " "
    !write(*,*) grd_spr,"grd_spr_aft"
    !write(*,*) grd_area,"grd_area_aft"
    !write(*,*) grd_grvtnl,"grd_grvtnl_aft"
    
    !write(*,*) gradient(3,1:2),gradient(4,1:2),"grad"
    !write(*,*) gradient(5,1:2),gradient(6,1:2),"grad"
    write(*,*) node_xy,"node_xy"
    Enrgy_frst = Energy(node_xy,l0,A0)
    write(*,*) Enrgy_frst,"Enrgy_frst"
    test_no = 4
    
    call store_node_xy
    
    write(*,*) grdmv_xy,"grdmv_xy_bfr"
    
    funcChoice = 1 ; ftol=1.d-08  
    call frprmn(coordntes_xy,grdmv_xy,iter,funcChoice,N_mvCoordnte,ftol,fret,fpp)
    
    write(*,*) iter,"iter_aft_frprmn"
    write(*,*) coordntes_xy,"aft_frprmn" !co8
    write(*,*) fret,"aft_frprmn" !co9
    
    !call generate_delX_and_delY(coordntes_xy,test_no,lp_cnt,N_iter)
    !call save_config_and_generate_data(coordntes_xy,test_no,lp_cnt)
    
    !call restore_node_xy
    !call restore_A0_area
    
    
    !write(unit=21,fmt=*) A0,"A0 aft restoring" !co10  
    
    
    
    call restore_node_type
    call deallocate_moving_coordnte_variables
    !call create_nonD_delA0_area_vs_coordnte_data(N_iter,test_no)
   
    

  end subroutine test_04
  
  
  
  
  subroutine test_05
    implicit none
    integer :: i,i1
    integer :: N_iter
    
    logical :: spr,area,grv
    real*8  :: Enrgy_frst

    integer :: htu !htu=how to update
    integer :: max_change
    
    integer,allocatable :: type_of_change(:)
    
    real*8  :: Eval,penF_val,Wval
    integer :: updt_cnt
    
    integer :: i2
    real*8  :: Es=0.0d0,Ea=0.0d0,Eg=0.0d0
    real*8  :: Pfs=0.0d0,Pfa=0.0d0
    
    max_change = 10
    allocate(type_of_change(1:max_change))
    type_of_change(1:max_change) = 0
    
    updt_cnt = 0
    
    write(*,*) A0,"A0"
    write(*,*) k_spr,"k_spr"
    
    write(*,*) N_mvCoordnte,"N_mvCooordnte"
    write(*,*) N_mvCoordnte_withl0A0,"N_mvCooordnte_withl0A0"
    
    
    do i = 1,N_node
       write(*,*) node_xy(i,1:2),"node_xy for i =",i
    enddo
    
    Enrgy_frst = Energy(node_xy,l0,A0)
    write(*,*) Enrgy_frst,"Enrgy_frst"
    
    !stop
    
    test_no    = 5
    lp_cnt     = -1
    subtest_no = 2
    N_iter     = 2
    
    grd_mv = 10.0d5 ; grdmv_xy = 10.0d5       
    
    type_of_change(9) = 9    
    
    call nodesl0A0_to_coordntes(node_xy,l0,A0,coordntes)
    coordntes_xy(1:N_mvCoordnte) = coordntes(1:N_mvCoordnte)
    
    !for test this prtn
    funcChoice=1 ; N_variabls = N_mvCoordnte
    call gradient_coordntes(coordntes_xy,grdmv_xy)
    grd_mv(1:N_mvCoordnte) = grdmv_xy(1:N_mvCoordnte)
    !stop
    !for test this prtn
    
    funcChoice=2 ; N_variabls = N_mvCoordnte_withl0A0
    call grdW_coordntes(coordntes,grd_mv)
    grdmv_xy(1:N_mvCoordnte) = grd_mv(1:N_mvCoordnte)
    
    do i = 1,N_mvCoordnte_withl0A0
       write(*,*) coordntes(i),grd_mv(i),"coordntes,grd_mv"
    enddo
    
    write(*,*) N_mvCoordnte_withl0A0,"N_mvCoordnte_withl0A0"

    if (stageNo==4 .and. stageType==1) then
       if (strctNo==1) Exprmnt=24 ; Frame=1
       if (strctNo==2) Exprmnt=25 ; Frame=1
    elseif (stageNo==4 .and. stageType==2) then
       if (strctNo==1) Exprmnt=14 ; Frame=1
       if (strctNo==2) Exprmnt=15 ; Frame=1
    endif
    
    do i = 1,N_iter
       
       lp_cnt = lp_cnt + 1    
       
       call update_and_equilibrate(type_of_change,max_change,updt_cnt)
       updt_cnt = updt_cnt + 1
       
       call write_coordntes_as_nodes(coordntes_xy)   
       
       call coordntes_to_nodesl0A0(coordntes,node_xy,l0,A0)
       call save_l0_A0_nodeXY(updt_cnt)
       
       Wval=Wfunc(coordntes);penF_val=penF(node_xy);Eval=Energy(node_xy,l0,A0) 
       write(*,*) Wval,penF_val,Eval,"Wval,penF_val,Eval"
       
    enddo
    
    
    
    open(unit=133,file="l0A0kskaCoordntes.dat",position='append')
    
    write(133,fmt=*) strctno,"strctNo"
    
    do i = 1,N_spr
       write(133,fmt=*) l0(i),k_spr(i),"l0,kspr for spring nm =",i
    enddo
    write(133,fmt=*) " "
    
    do i = 1,N_cell
       write(133,fmt=*) A0(i),k_area(i),"A0,karea for cell nm =",i
    enddo
    write(133,fmt=*) " "
    
    do i = 1,N_mvCoordnte_withl0A0
       if (i.le.N_mvCoordnte) then
          write(133,fmt=*) coordntes(i),coordntes_xy(i),i,"coordntes,coordnteXY"
       elseif (i.gt.N_mvCoordnte .and. i.le.(N_mvCoordnte+N_spr)) then
          write(133,fmt=*) coordntes(i),l0(i-N_mvCoordnte),i,(i-N_mvCoordnte),"coordnte,l0,i,i-N_mvCoordnte"
       elseif (i.gt.(N_mvCoordnte+N_spr) .and. i.le.(N_mvCoordnte+N_spr+N_cell)) then
          write(133,fmt=*) coordntes(i),A0(i-N_mvCoordnte-N_spr),i,(i-N_mvCoordnte-N_spr),"coordnte,A0,i,i-N_mvCoordnte-N_spr"
       endif
    enddo
    
    write(133,fmt=*) N_mvCoordnte_withl0A0,N_mvCoordnte,N_spr,N_cell,"1=2+3+4"
    write(133,fmt=*) " "
    close(133)
    
  end subroutine test_05
  
  subroutine test_S1T1
    implicit none
    real*8  :: hm,hmSpr,hmArea
    real*8  :: hmMltplL(1:3),hmMltplK(1:3)
    integer :: i,sprNo
    integer :: NCpairs(1:2),NCP,decideSpr
    integer :: NC
    
    real*8  :: Enrgy_frst!,Esf,Eaf,Egf,Ebf,Esindvd
    logical :: lgcl_rdfn
    real*8  :: tol_Rdfn
    integer :: LftNeigh=0,RghtNeigh=0
    
    open(unit=87,file='FrameNoAtDiffStrct.dat',position='append')
    
    if (strctNo==1) then
       write(*,*) N_mvCoordnte,"N_mvCooordnte"
       write(*,*) N_mvCoordnte_withl0A0,"N_mvCooordnte_withl0A0"
       
       Enrgy_frst = Energy(node_xy,l0,A0)
       write(*,*) Enrgy_frst,"Enrgy_frst"
       
       grd_mv = 10.0d5 ; grdmv_xy = 10.0d5       
       
       call nodesl0A0_to_coordntes(node_xy,l0,A0,coordntes)
       coordntes_xy(1:N_mvCoordnte) = coordntes(1:N_mvCoordnte)
       
       Exprmnt=21 ; Frame=1
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
    endif
    
    !No pressure/Tension to Initial Version
    
    if (strctNo==1) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       
       !hm = 0.015d0
       !call pulling_or_pushing_the_embryo(Exprmnt,Frame,hm)
       
       !call successiveChangein_f   
       !call changeA0_inSteps_Equilibrate_and_Optimize(Exprmnt,Frame)
       
       call read_strctProps_strct(strctNo)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
       !Initial version to First Meet (Frankenstein)
       
    elseif (strctNo==2) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       
       !hmMltplL(1)=0.99d0; hmMltplL(2)=-0.40d0 ; hmMltplL(3)=0.00d0
       !hmMltplK(1)=10.0d0; hmMltplK(2)= 0.10d0 ; hmMltplK(3)=0.00d0
       !call Stepwise_SprShorteningIC_and_Equilibrate(Exprmnt,Frame,hmMltplL,hmMltplK)
       !write(*,*) Frame,"Frame aft IC"
       
       !stop (Final)
       
       call read_strctProps_strct(strctNo)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
       
       lgcl_rdfn = .False.
       Hold_lgclRdfn = .False.
       tol_Rdfn  = 0.05d0
       call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
       
       if (lgcl_rdfn.eqv..True.) then
          write(*,*) " "
          write(*,*) "HURRAY, THE CELLS MEET 1st time !!!!!!!"
          write(*,*) " "
          CellsMeet = 1
          Hold_lgclRdfn = lgcl_rdfn
          
          !MeetData=1
          !call saveA0l0_0fS1T1 !will'be shiftd build_structures.f08
          !call sleep(1)
          !call taking_lftrghtOriginNeigh_downwards(LftNeigh,RghtNeigh)
          
       elseif (lgcl_rdfn.eqv..False.) then
          write(*,*) " "
          write(*,*) " :( Cells did not meet yet, stop here to adjust (1st)"
          stop
       endif
       
    elseif (strctNo==3) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       
       lgcl_rdfn = Hold_lgclRdfn
       
       if (lgcl_rdfn.eqv..True.) then
          call taking_lftrghtOriginNeigh_downwards(LftNeigh,RghtNeigh)
       elseif (lgcl_rdfn.eqv..False.) then
          write(*,*) "lgcl is not true,stop"
          stop
       endif
       
       do i=1,N_node
          write(*,*) node_xy(i,1:2),i,"node"
       enddo
       
       !call sleep(1)
       
       sprNo = cntrlCell_Spring(1)
       dmlshDecsn = 1
       call demolish_spring_and_redefine_system(sprNo)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
       
       !*****First Meet(Frankenstein) to Adjusted Version *****
       
    elseif (strctNo==4) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       
       hm=-0.20d0 !(+/-) for (Expn/Sqz)
       call Stepwise_ExpandingIC(Exprmnt,Frame,hm)
       hmSpr = -hm*3.0d0
       call adjust_BsalLtrlSprIC_and_Equilibrate(Exprmnt,Frame,hmSpr)
       
       NCP=1 ; hmMltplL(1)=-1.60d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.00d0
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       NCP=2 ; hmMltplL(1)=-0.60d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.00d0
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       NCP=3 ; hmMltplL(1)=-0.60d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.00d0
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       NCP=4 ; hmMltplL(1)=-0.60d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.00d0
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       NCP=5 ; hmMltplL(1)=-0.60d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.00d0
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       NCP=6 ; hmMltplL(1)=-0.50d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.00d0
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       NCP=7 ; hmMltplL(1)=-0.40d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.00d0
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       NCP=8 ; hmMltplL(1)=-0.30d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.00d0
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       
       NCpairs(1)=2 ; NCpairs(2)=3 ; decideSpr=3 ; hm=0.35d0
       call StepwiseSprShorten_MltplNC_Equilibrate(Exprmnt,Frame,NCpairs,decideSpr,hm)
       
       hmMltplL(1)=0.00d0 ; hmMltplL(2)=-1.20d0 ; hmMltplL(3)=0.00d0
       hmMltplK(1)=0.00d0 ; hmMltplK(2)= 0.00d0 ; hmMltplK(3)=0.00d0
       call Stepwise_SprShorteningIC_and_Equilibrate(Exprmnt,Frame,hmMltplL,hmMltplK)
       
       !stop    
       
    !****** First Meet to Second Meet (Frankenstein)**********
    elseif (strctNo==5) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       write(*,*) Frame,strctNo,"Frme,strct=5"
       call sleep(1)
       
       hmSpr = 0.25d0 !change to 0.99
       call adjust_BsalLtrlSprIC_and_Equilibrate(Exprmnt,Frame,hmSpr) !1

       NCP=1 ; hmArea = 0.25d0 
       call Stepwise_ExpandingNC(Exprmnt,Frame,NCP,hmArea)
       NCP=2 ; hmArea = 0.25d0
       call Stepwise_ExpandingNC(Exprmnt,Frame,NCP,hmArea)
       
       !stop
       
       NCP=1 ; hmMltplL(1)=-0.60d0 ; hmMltplL(2)=0.80d0 ; hmMltplL(3)=0.00d0
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL) !3
       NCP=2 ; hmMltplL(1)=-0.60d0 ; hmMltplL(2)=0.90d0 ; hmMltplL(3)=0.00d0
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL) !5
       
       NCP=1 ; hmMltplL(1)=0.00d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=-0.50d0 
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL) !6
       NCP=2 ; hmMltplL(1)=0.00d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=-0.05d0
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL) !7
       
       NCP=3 ; hmMltplL(1)=-0.35d0 ; hmMltplL(2)=0.35d0 ; hmMltplL(3)=0.00d0 
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL) !9
       NCP=4 ; hmMltplL(1)=-0.35d0 ; hmMltplL(2)=0.35d0 ; hmMltplL(3)=0.00d0
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL) !11
       
       stop !(Final stop)
       
       call read_strctProps_strct(strctNo)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
       
       lgcl_rdfn = .False.
       Hold_lgclRdfn = .False.
       tol_Rdfn  = 0.05d0
       call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
       
       if (lgcl_rdfn.eqv..True.) then
          write(*,*) " "
          write(*,*) "HURRAY, THE CELLS MEET 2nd Time !!!!!!!"
          write(*,*) " "
          CellsMeet = 2
          Hold_lgclRdfn = lgcl_rdfn
          
          !MeetData=2
          !call saveA0l0_0fS1T1 !will'be shiftd build_structures.f08
          !call sleep(1)
          !call taking_lftrghtOriginNeigh_downwards(LftNeigh,RghtNeigh)
          
       elseif (lgcl_rdfn.eqv..False.) then
          write(*,*) " "
          write(*,*) " :( Cells did not meet yet, stop here to adjust (2nd)"
          stop
       endif
       
    elseif (strctNo==6) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       
       lgcl_rdfn = Hold_lgclRdfn
       
       if (lgcl_rdfn.eqv..True.) then
          call taking_lftrghtOriginNeigh_downwards(LftNeigh,RghtNeigh)
       elseif (lgcl_rdfn.eqv..True.) then
          write(*,*) "lgcl is not true,stop"
          stop
       endif
       
       dmlshDecsn = 0
       call redefine_system_wo_demolishSpr
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
       
       ! ****Frankenstein Second Meet to Adjusted Ones
       
    elseif (strctNo==7) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       
       !NCP=1 ; hmMltplL(1)=0.00d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.50d0
       !call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       !NCP=1 ; hmMltplL(1)=0.99d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.00d0
       !call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       
       
       !NCP=3 ; hmMltplL(1)=-0.60d0 ; hmMltplL(2)=0.90d0 ; hmMltplL(3)=0.00d0
       !call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       
       !NCP=2 ; hmMltplL(1)=0.00d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=-0.38d0
       !call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       
       !stop !(Final stop)
       
       call read_strctProps_strct(strctNo)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
       lgcl_rdfn = .False.
       Hold_lgclRdfn = .False.
       tol_Rdfn  = 0.05d0
       call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
       
       if (lgcl_rdfn.eqv..True.) then
          write(*,*) " "
          write(*,*) "HURRAY, THE CELLS MEET 3rd Time !!!!!!!"
          write(*,*) " "
          
          CellsMeet=3
          Hold_lgclRdfn = lgcl_rdfn
          
       elseif (lgcl_rdfn.eqv..False.) then
          write(*,*) " "
          write(*,*) " :( Cells did not meet yet, stop here to adjust (3rd)"
          stop
       endif
       
    elseif (strctNo==8) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"

       lgcl_rdfn = Hold_lgclRdfn
       
       if (lgcl_rdfn.eqv..True.) then
          call taking_lftrghtOriginNeigh_downwards(LftNeigh,RghtNeigh)
       elseif (lgcl_rdfn.eqv..True.) then
          write(*,*) "lgcl is not true,stop"
          stop
       endif
       
       dmlshDecsn = 0
       call redefine_system_wo_demolishSpr
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
    elseif (strctNo==9) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       
       !NCP=2 ; hmMltplL(1)=0.00d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.50d0
       !call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       !NCP=2 ; hmMltplL(1)=0.60d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.00d0
       !call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       
       
       !NCP=4 ; hmMltplL(1)=-0.60d0 ; hmMltplL(2)=0.90d0 ; hmMltplL(3)=0.00d0
       !call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       
       !NCP=3 ; hmMltplL(1)=0.00d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=-0.38d0
       !call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       
       
       !NCP=1 ; hmMltplL(1)=0.00d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.30d0
       !call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       !NCP=2 ; hmMltplL(1)=0.00d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.40d0
       !call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       
       !stop (Final stop)
       
       
       call read_strctProps_strct(strctNo)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
       lgcl_rdfn = .False.
       Hold_lgclRdfn = .False.
       tol_Rdfn  = 0.05d0
       call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
       
       if (lgcl_rdfn.eqv..True.) then
          write(*,*) " "
          write(*,*) "HURRAY, THE CELLS MEET 4th Time !!!!!!!"
          write(*,*) " "
          
          CellsMeet=4
          Hold_lgclRdfn = lgcl_rdfn
          
       elseif (lgcl_rdfn.eqv..False.) then
          write(*,*) " "
          write(*,*) " :( Cells did not meet yet, stop here to adjust (4th)"
          stop
       endif
       
    elseif (strctNo==10) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       
       lgcl_rdfn = Hold_lgclRdfn
       
       if (lgcl_rdfn.eqv..True.) then
          call taking_lftrghtOriginNeigh_downwards(LftNeigh,RghtNeigh)
       elseif (lgcl_rdfn.eqv..True.) then
          write(*,*) "lgcl is not true,stop"
          stop
       endif
       
       dmlshDecsn = 0
       call redefine_system_wo_demolishSpr
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
    elseif (strctNo==11) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       
       !NCP=3 ; hmMltplL(1)=0.00d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.50d0
       !call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       !NCP=3 ; hmMltplL(1)=0.60d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.00d0
       !call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       
       
       !NCP=5 ; hmMltplL(1)=-0.60d0 ; hmMltplL(2)=0.90d0 ; hmMltplL(3)=0.00d0
       !call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       
       !NCP=4 ; hmMltplL(1)=0.00d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=-0.38d0
       !call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       
       
       !NCP=1 ; hmMltplL(1)=0.00d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.30d0
       !call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       !NCP=2 ; hmMltplL(1)=0.00d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.10d0
       !call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       !NCP=3 ; hmMltplL(1)=0.00d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.24d0
       !call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)

       !stop (Final stop)
       
       call read_strctProps_strct(strctNo)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
       lgcl_rdfn = .False.
       Hold_lgclRdfn = .False.
       tol_Rdfn  = 0.05d0
       call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
       
       if (lgcl_rdfn.eqv..True.) then
          write(*,*) " "
          write(*,*) "HURRAY, THE CELLS MEET 5th Time !!!!!!!"
          write(*,*) " "
          
          CellsMeet=5
          Hold_lgclRdfn = lgcl_rdfn
          
       elseif (lgcl_rdfn.eqv..False.) then
          write(*,*) " "
          write(*,*) " :( Cells did not meet yet, stop here to adjust (5th)"
          stop
       endif
       
    elseif (strctNo==12) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"

       lgcl_rdfn = Hold_lgclRdfn
       
       if (lgcl_rdfn.eqv..True.) then
          call taking_lftrghtOriginNeigh_downwards(LftNeigh,RghtNeigh)
       elseif (lgcl_rdfn.eqv..True.) then
          write(*,*) "lgcl is not true,stop"
          stop
       endif
       
       dmlshDecsn = 0
       call redefine_system_wo_demolishSpr
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
    elseif (strctNo==13) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       
       NCP=4 ; hmMltplL(1)=0.00d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.50d0
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       
       NCP=6 ; hmMltplL(1)=-0.80d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.00d0
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       NCP=3 ; hmMltplL(1)=-1.00d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.00d0
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       
       NCP=1 ; hmMltplL(1)=0.00d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.90d0
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       NCP=2 ; hmMltplL(1)=0.00d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.90d0
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       NCP=3 ; hmMltplL(1)=0.00d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.60d0
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       
       NCP=6 ; hmMltplL(1)=0.00d0 ; hmMltplL(2)=0.75d0 ; hmMltplL(3)=0.00d0
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)

       !!comparing frm datas
       
       NCP=3 ; hmMltplL(1)=-1.20d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.00d0
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       NCP=2 ; hmMltplL(1)=-1.20d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.00d0
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)

       NCP=1 ; hmMltplL(1)= 0.00d0 ; hmMltplL(2)=-0.50d0 ; hmMltplL(3)=0.00d0
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       
       k_spr(16) = 0.50*k_spr(16) ; k_spr(40) = 0.50*k_spr(40) !NCP3 apcl
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
       k_spr(20) = 0.50*k_spr(20) ; k_spr(44) = 0.50*k_spr(44) !NCP2 bsal
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
       k_spr(22) = 5.0*k_spr(22) ; k_spr(46) = 5.0*k_spr(46) !NCP1 apcl
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
       k_spr(21) = 3.50*k_spr(21) ; k_spr(45) = 3.50*k_spr(45)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN

       !k_spr(24) = 2.00*k_spr(24) ; k_spr(48) = 2.00*k_spr(48)
       !call Equilibrate_system
       !call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       !Frame = Frame+1
       !call switchto_NI_model_run_and_switchbackto_TN

       A0(8) = 1.25*A0(8) ; A0(16) = 1.25*A0(16)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
    endif
    
    close(87)
    
  end subroutine test_S1T1


  
  subroutine test_S1T1_TargetVersion
    implicit none
    
    real*8  :: hm,hmSpr,hmArea
    real*8  :: hmMltplL(1:3),hmMltplK(1:3)
    integer :: i,sprNo
    integer :: NCpairs(1:2),NCP,decideSpr
    integer :: NC
    
    real*8  :: Enrgy_frst!,Esf,Eaf,Egf,Ebf,Esindvd
    logical :: lgcl_rdfn
    real*8  :: tol_Rdfn
    integer :: LftNeigh=0,RghtNeigh=0
    integer :: conVstrct

    real*8  :: l0_incr,ks_incr,A0_incr,ka_incr
    real*8,allocatable  :: A0F_finl(:),kaF_finl(:),l0F_finl(:),ksF_finl(:)
    real*8,allocatable  :: A0F_init(:),kaF_init(:),l0F_init(:),ksF_init(:)
    
    integer :: SmallV
    
    SmallV=0.0000001d0
    
    open(unit=87,file='FrameNoAtDiffStrct.dat',position='append')
    
    if (strctNo==1) then
       write(*,*) N_mvCoordnte,"N_mvCooordnte"
       write(*,*) N_mvCoordnte_withl0A0,"N_mvCooordnte_withl0A0"
       
       Enrgy_frst = Energy(node_xy,l0,A0)
       write(*,*) Enrgy_frst,"Enrgy_frst"
       
       grd_mv = 10.0d5 ; grdmv_xy = 10.0d5       
       
       call nodesl0A0_to_coordntes(node_xy,l0,A0,coordntes)
       coordntes_xy(1:N_mvCoordnte) = coordntes(1:N_mvCoordnte)
       
       Exprmnt=21 ; Frame=1
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
    endif
    
    !No pressure/Tension to Initial Version
    
    if (strctNo==1) then
       CyclNo = 1
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       
       hm = 0.015d0
       call pulling_or_pushing_the_embryo(Exprmnt,Frame,hm)
       
       call successiveChangein_f   
       call changeA0_inSteps_Equilibrate_and_Optimize(Exprmnt,Frame)
       
       !call read_strctProps_strct(strctNo)
       !call Equilibrate_system
       !call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       !Frame = Frame+1
       !call switchto_NI_model_run_and_switchbackto_TN
       
       !Initial version to First Meet (Frankenstein)
       
    elseif (strctNo==2) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       
       !hmMltplL(1)=0.99d0; hmMltplL(2)=-0.40d0 ; hmMltplL(3)=0.00d0
       !hmMltplK(1)=10.0d0; hmMltplK(2)= 0.10d0 ; hmMltplK(3)=0.00d0
       !call Stepwise_SprShorteningIC_and_Equilibrate(Exprmnt,Frame,hmMltplL,hmMltplK)
       !write(*,*) Frame,"Frame aft IC"
       
       !stop (Final)
       
       call read_strctProps_strct(strctNo)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
       
       lgcl_rdfn = .False.
       Hold_lgclRdfn = .False.
       tol_Rdfn  = 0.05d0
       call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
       
       if (lgcl_rdfn.eqv..True.) then
          write(*,*) " "
          write(*,*) "HURRAY, THE CELLS MEET 1st time !!!!!!!"
          write(*,*) " "
          CellsMeet = 1
          Hold_lgclRdfn = lgcl_rdfn
          
       elseif (lgcl_rdfn.eqv..False.) then
          write(*,*) " "
          write(*,*) " :( Cells did not meet yet, stop here to adjust"
          stop
       endif
       write(*,*) lgcl_rdfn,"lgcl Rdfn2"
       call sleep(1)
       
    elseif (strctNo==3) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       
       lgcl_rdfn = Hold_lgclRdfn
       write(*,*) lgcl_rdfn,"lgcl Rdfn3"
       call sleep(1)
       
       if (lgcl_rdfn.eqv..True.) then
          call taking_lftrghtOriginNeigh_downwards(LftNeigh,RghtNeigh)
       elseif (lgcl_rdfn.eqv..False.) then
          write(*,*) "lgcl is not true,stop"
          stop
       endif
       
       sprNo = cntrlCell_Spring(1)
       dmlshDecsn = 1
       call demolish_spring_and_redefine_system(sprNo)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
       
       !*****First Meet(Frankenstein) to Adjusted Version *****
       
    elseif (strctNo==4) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       call convert_strctData_from_S4_to_S1
       
       CyclNo = 2
       conVstrct = 1
       call read_conVstrctProps_with_adjstmnt(conVstrct)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
       allocate(A0_I(1:N_cell),ka_I(1:N_cell),l0_I(1:N_spr),ks_I(1:N_spr))
       allocate(A0_F(1:N_cell),ka_F(1:N_cell),l0_F(1:N_spr),ks_F(1:N_spr))
       
       A0_I = 1.0d30 ; ka_I = 1.0d30 ; l0_I = 1.0d30 ; ks_I = 1.0d30
       A0_F = 1.0d30 ; ka_F = 1.0d30 ; l0_F = 1.0d30 ; ks_F = 1.0d30
       
       A0_I(1:N_cell) = A0(1:N_cell) ; ka_I(1:N_cell) = k_area(1:N_cell)
       l0_I(1:N_spr)  = l0(1:N_spr)  ; ks_I(1:N_spr)  = k_spr(1:N_spr)
       
       
       !****** First Meet to Second Meet (Frankenstein)**********
       
    elseif (strctNo==5) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       write(*,*) CyclNo,"CyCl must be 2 here"
       
       conVstrct = 3
       call read_conVstrctProps_with_adjstmnt(conVstrct)
       !call read_strctProps_strct(strctNo)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
       write(87,fmt=*) Frame,"Frame at which Extrapolation begins 2"
       
       A0_F(1:N_cell) = A0(1:N_cell) ; ka_F(1:N_cell) = k_area(1:N_cell)
       l0_F(1:N_spr)  = l0(1:N_spr)  ; ks_F(1:N_spr)  = k_spr(1:N_spr)
       
       call Extrapolate_and_adjust_strctProp_ifNeeded_in_MAE(Exprmnt,Frame)
       
       write(87,fmt=*) (Frame-1),"Frame at which Extrapolation ends 2"

       do i=1,N_spr
          write(*,*) l0(i),l0_F(i),k_spr(i),ks_F(i),"l,l0F,ks,ksF,spr"
       enddo
       
       do i=1,N_cell
          write(*,*) A0(i),A0_F(i),k_area(i),ka_F(i),"A,A0F,ka,kaF,area"
       enddo
       
       deallocate(A0_I,ka_I,l0_I,ks_I)
       deallocate(A0_F,ka_F,l0_F,ks_F)
       
       
       lgcl_rdfn = .False.
       Hold_lgclRdfn = .False.
       tol_Rdfn  = 0.051d0
       call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
       
       if (lgcl_rdfn.eqv..True.) then
          write(*,*) " "
          write(*,*) "HURRAY, THE CELLS MEET 2nd Time !!!!!!!"
          write(*,*) " "
          CellsMeet = 2
          Hold_lgclRdfn = lgcl_rdfn
          
       elseif (lgcl_rdfn.eqv..False.) then
          write(*,*) " "
          write(*,*) " :( Cells did not meet yet, stop here to adjust"
          stop
       endif
       
    elseif (strctNo==6) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       
       lgcl_rdfn = Hold_lgclRdfn
       
       if (lgcl_rdfn.eqv..True.) then
          call taking_lftrghtOriginNeigh_downwards(LftNeigh,RghtNeigh)
       elseif (lgcl_rdfn.eqv..False.) then
          write(*,*) "lgcl is not true,stop"
          stop
       endif
       
       dmlshDecsn = 0
       call redefine_system_wo_demolishSpr
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
       
       ! ****Frankenstein Second Meet to Adjusted Ones
       
       
    elseif (strctNo==7) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       
       CyclNo = 3
       conVstrct = 1
       call read_conVstrctProps_with_adjstmnt(conVstrct)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
       allocate(A0_I(1:N_cell),ka_I(1:N_cell),l0_I(1:N_spr),ks_I(1:N_spr))
       allocate(A0_F(1:N_cell),ka_F(1:N_cell),l0_F(1:N_spr),ks_F(1:N_spr))
       
       A0_I = 1.0d30 ; ka_I = 1.0d30 ; l0_I = 1.0d30 ; ks_I = 1.0d30
       A0_F = 1.0d30 ; ka_F = 1.0d30 ; l0_F = 1.0d30 ; ks_F = 1.0d30
       
       A0_I(1:N_cell) = A0(1:N_cell) ; ka_I(1:N_cell) = k_area(1:N_cell)
       l0_I(1:N_spr)  = l0(1:N_spr)  ; ks_I(1:N_spr)  = k_spr(1:N_spr)
       
       
    elseif (strctNo==8) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       write(*,*) CyclNo,"CyCl must be 3 here"
       
       conVstrct = 3
       call read_conVstrctProps_with_adjstmnt(conVstrct)
       !call read_strctProps_strct(strctNo)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN


       write(87,fmt=*) Frame,"Frame at which Extrapolation begins 3"
       
       A0_F(1:N_cell) = A0(1:N_cell) ; ka_F(1:N_cell) = k_area(1:N_cell)
       l0_F(1:N_spr)  = l0(1:N_spr)  ; ks_F(1:N_spr)  = k_spr(1:N_spr)
       
       call Extrapolate_and_adjust_strctProp_ifNeeded_in_MAE(Exprmnt,Frame)
       
       write(87,fmt=*) (Frame-1),"Frame at which Extrapolation ends 3"

       do i=1,N_spr
          write(*,*) l0(i),l0_F(i),k_spr(i),ks_F(i),"l,l0F,ks,ksF,spr"
       enddo
       
       do i=1,N_cell
          write(*,*) A0(i),A0_F(i),k_area(i),ka_F(i),"A,A0F,ka,kaF,area"
       enddo
       
       deallocate(A0_I,ka_I,l0_I,ks_I)
       deallocate(A0_F,ka_F,l0_F,ks_F)
       
       
       lgcl_rdfn = .False.
       Hold_lgclRdfn = .False.
       tol_Rdfn  = 0.05d0
       call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
       
       if (lgcl_rdfn.eqv..True.) then
          write(*,*) " "
          write(*,*) "HURRAY, THE CELLS MEET 3rd Time !!!!!!!"
          write(*,*) " "
          
          CellsMeet=3
          Hold_lgclRdfn = lgcl_rdfn
          
       elseif (lgcl_rdfn.eqv..False.) then
          write(*,*) " "
          write(*,*) " :( Cells did not meet yet, stop here to adjust"
          stop
       endif
       
    elseif (strctNo==9) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       
       lgcl_rdfn = Hold_lgclRdfn
       
       if (lgcl_rdfn.eqv..True.) then
          call taking_lftrghtOriginNeigh_downwards(LftNeigh,RghtNeigh)
       elseif (lgcl_rdfn.eqv..True.) then
          write(*,*) "lgcl is not true,stop"
          stop
       endif
       
       dmlshDecsn = 0
       call redefine_system_wo_demolishSpr
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
       
    elseif (strctNo==10) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       
       CyclNo = 4
       conVstrct = 1
       call read_conVstrctProps_with_adjstmnt(conVstrct)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
       allocate(A0_I(1:N_cell),ka_I(1:N_cell),l0_I(1:N_spr),ks_I(1:N_spr))
       allocate(A0_F(1:N_cell),ka_F(1:N_cell),l0_F(1:N_spr),ks_F(1:N_spr))
       
       A0_I = 1.0d30 ; ka_I = 1.0d30 ; l0_I = 1.0d30 ; ks_I = 1.0d30
       A0_F = 1.0d30 ; ka_F = 1.0d30 ; l0_F = 1.0d30 ; ks_F = 1.0d30
       
       A0_I(1:N_cell) = A0(1:N_cell) ; ka_I(1:N_cell) = k_area(1:N_cell)
       l0_I(1:N_spr)  = l0(1:N_spr)  ; ks_I(1:N_spr)  = k_spr(1:N_spr)
        
    elseif (strctNo==11) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       write(*,*) CyclNo,"CyCl must be 4 here"
       
       conVstrct = 3
       call read_conVstrctProps_with_adjstmnt(conVstrct)
       !call read_strctProps_strct(strctNo)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
       write(87,fmt=*) Frame,"Frame at which Extrapolation begins 4"
       
       A0_F(1:N_cell) = A0(1:N_cell) ; ka_F(1:N_cell) = k_area(1:N_cell)
       l0_F(1:N_spr)  = l0(1:N_spr)  ; ks_F(1:N_spr)  = k_spr(1:N_spr)
       
       call Extrapolate_and_adjust_strctProp_ifNeeded_in_MAE(Exprmnt,Frame)
       
       write(87,fmt=*) (Frame-1),"Frame at which Extrapolation ends 4"
       
       do i=1,N_spr
          write(*,*) l0(i),l0_F(i),k_spr(i),ks_F(i),"l,l0F,ks,ksF,spr"
       enddo
       
       do i=1,N_cell
          write(*,*) A0(i),A0_F(i),k_area(i),ka_F(i),"A,A0F,ka,kaF,area"
       enddo
       
       deallocate(A0_I,ka_I,l0_I,ks_I)
       deallocate(A0_F,ka_F,l0_F,ks_F)
       
       lgcl_rdfn = .False.
       Hold_lgclRdfn = .False.
       tol_Rdfn  = 0.05d0
       call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
       
       if (lgcl_rdfn.eqv..True.) then
          write(*,*) " "
          write(*,*) "HURRAY, THE CELLS MEET 4th Time !!!!!!!"
          write(*,*) " "
          
          CellsMeet=4
          Hold_lgclRdfn = lgcl_rdfn
          
       elseif (lgcl_rdfn.eqv..False.) then
          write(*,*) " "
          write(*,*) " :( Cells did not meet yet, stop here to adjust"
          stop
       endif
       
    elseif (strctNo==12) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       
       lgcl_rdfn = Hold_lgclRdfn
       
       if (lgcl_rdfn.eqv..True.) then
          call taking_lftrghtOriginNeigh_downwards(LftNeigh,RghtNeigh)
       elseif (lgcl_rdfn.eqv..True.) then
          write(*,*) "lgcl is not true,stop"
          stop
       endif
       
       dmlshDecsn = 0
       call redefine_system_wo_demolishSpr
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
       
    elseif (strctNo==13) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       
       CyclNo = 5
       conVstrct = 1
       call read_conVstrctProps_with_adjstmnt(conVstrct)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
    elseif (strctNo==14) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       
       conVstrct = 3
       call read_conVstrctProps_with_adjstmnt(conVstrct)
       !call read_strctProps_strct(strctNo)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN

       NCP=6 ; hmMltplL(1)=-0.70d0 ; hmMltplL(2)=0.70d0 ; hmMltplL(3)=0.00d0
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       
       NCP=4 ; hmMltplL(1)=+0.00d0 ; hmMltplL(2)=0.50d0 ; hmMltplL(3)=0.00d0
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       
       NCP=7 ; hmMltplL(1)=-0.50d0 ; hmMltplL(2)=0.00d0 ; hmMltplL(3)=0.00d0
       call Stepwise_SprShorteningNC_and_Equilibrate(Exprmnt,Frame,NCP,hmMltplL)
       
       
       lgcl_rdfn = .False.
       Hold_lgclRdfn = .False.
       tol_Rdfn  = 0.05d0
       call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
       
       if (lgcl_rdfn.eqv..True.) then
          write(*,*) " "
          write(*,*) "HURRAY, THE CELLS MEET 5th Time !!!!!!!"
          write(*,*) " "
          
          CellsMeet=5
          Hold_lgclRdfn = lgcl_rdfn
          
       elseif (lgcl_rdfn.eqv..False.) then
          write(*,*) " "
          write(*,*) " :( Cells did not meet yet, stop here to adjust"
          stop
       endif
       
    elseif (strctNo==15) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"

       lgcl_rdfn = Hold_lgclRdfn
       
       if (lgcl_rdfn.eqv..True.) then
          call taking_lftrghtOriginNeigh_downwards(LftNeigh,RghtNeigh)
       elseif (lgcl_rdfn.eqv..True.) then
          write(*,*) "lgcl is not true,stop"
          stop
       endif
       
       dmlshDecsn = 0
       call redefine_system_wo_demolishSpr
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
    elseif (strctNo==16) then
       write(87,fmt=*) Frame,strctNo,"Frme,strct"
       
       conVstrct = 2
       call read_conVstrctProps_with_adjstmnt(conVstrct)
       !call read_strctProps_strct(strctNo)
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
    endif
    
    close(87)
    
  end subroutine test_S1T1_TargetVersion
  
  
  subroutine test_S1T2
    implicit none
    real*8 :: hm
    
    Exprmnt=21 ; Frame=1
    
    !hm = 0.20d0
    !call Stepwise_shortening_ApclSprofIC_and_Equilibrate(Exprmnt,Frame,hm)
    
  end subroutine test_S1T2
  
  subroutine saveA0l0_0fS1T1
    implicit none
    integer :: i
    
    if (MeetData==0) then
       write(*,*) "MeetData cant be ZERO"
       stop
    elseif (MeetData==1) then
       open(unit=148,file='IntrmdStrucPropS1T1C_FrstMeet.dat')!C=Checking purpose
       open(unit=149,file='IntrmdStrucPropS1T1R_FrstMeet.dat')!R=Reading  purpose
    elseif (MeetData==2) then
       open(unit=148,file='IntrmdStrucPropS1T1C_ScndMeet.dat')!C=Checking purpose
       open(unit=149,file='IntrmdStrucPropS1T1R_ScndMeet.dat')!R=Reading  purpose
    endif
    
       
    call coordntes_to_nodes(coordntes_xy,node_xy)
    call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
    
    do i=1,N_spr
       write(148,*) k_spr(i),l0(i),i,"k_spr,l0,l,spr_nm"!,l(i) 
       write(149,*) k_spr(i),l0(i)
    enddo
    
    write(148,*) " "
    write(149,*) " "
    
    do i=1,N_cell
       write(148,*) k_area(i),A0(i),i,"k_area,A0,A,area_nm"!,A(i)
       write(149,*) k_area(i),A0(i)
    enddo
    write(148,*) " "
    write(149,*) " "

    do i=1,N_node
       write(148,*) CgXNode(i),i,"CgXNode,node_nm"
       write(149,*) CgXNode(i)
    enddo
    write(148,*) " "
    write(149,*) " "
    
    do i=1,N_mvCoordnte
       write(148,*) coordntes_xy(i),i,"coordntes_xy,i"
       write(149,*) coordntes_xy(i)
    enddo
    write(148,*) " "
    write(149,*) " "

    do i=1,N_node
       write(148,*) node_xy(i,1:N_dmnsn),i,"node_xy,node_nm"
       write(149,*) node_xy(i,1:N_dmnsn)
    enddo
    write(148,*) " "
    write(149,*) " "
    
    close(148)
    close(149)
    
  end subroutine saveA0l0_0fS1T1
  
  subroutine read_A0l0_ofS1T1_MeetData
    implicit none
    integer :: i,j,jmax
    integer :: N_itm
    character(len=100) :: flnm
    
    if (MeetData==1) then
       flnm='IntrmdStrucPropS1T1R_FrstMeet.dat'
       write(*,*) trim(adjustl(flnm)),"flnm"
       call sleep(1)
       
    elseif (MeetData==2) then
       flnm='IntrmdStrucPropS1T1R_ScndMeet.dat'
       write(*,*) trim(adjustl(flnm)),"flnm"
       call sleep(1)
    endif
    
    open(unit=21,file=trim(adjustl(flnm)))
    
    if (MeetData==1) then
       open(unit=22,file='ReadChkS1T1Prop_FrstMeet.dat')
    elseif (MeetData==2) then
       open(unit=22,file='ReadChkS1T1Prop_ScndMeet.dat')
    endif
    
    N_itm = 5
    
    do i=1,N_itm
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       if (i==3) jmax=N_node
       if (i==4) jmax=N_mvCoordnte   
       if (i==5) jmax=N_node
       
       do j = 1,jmax
          if (i==1) read(21,*)  k_spr(j),l0(j)!,l(j)
          if (i==1) write(22,*) k_spr(j),l0(j)!,l(j)
          if (i==2) read(21,*)  k_area(j),A0(j)!,A(j)
          if (i==2) write(22,*) k_area(j),A0(j)!,A(j)
          if (i==3) read(21,*)  CgXNode(j),CgYNode(j)
          if (i==3) write(22,*) CgXNode(j),CgYNode(j)
          if (i==4) read(21,*)  coordntes_xy(j)
          if (i==4) write(22,*) coordntes_xy(j)
          if (i==5) read(21,*)  node_xy(j,1:N_dmnsn)
          if (i==5) write(22,*) node_xy(j,1:N_dmnsn)
       enddo
       
    enddo
    
    close(21)
    close(22)
    
  end subroutine read_A0l0_ofS1T1_MeetData
  
  subroutine read_strctProps_strct(strctNum)
    implicit none
    integer, intent(in) :: strctNum
    integer :: i,j,jmax
    integer :: N_itm
    
    character(len=100) :: flnm
    character(len=100) :: flnmbr
    character(len=100) :: ReadFull_flnm,WriteFull_flnm
    
    flnm='strctProps'
    N_itm = 5
    
    if (strctNum.le.9) then
       write(flnmbr,'(i1.1,a)') strctNum,'S1T1.dat'
    elseif (strctNum.gt.9 .and. strctNum.le.99) then
       write(flnmbr,'(i2.2,a)') strctNum,'S1T1.dat'
    endif
    
    write(*,*) trim(adjustl(flnm)),trim(adjustl(flnmbr))
    ReadFull_flnm =trim(adjustl(flnm))//trim(adjustl(flnmbr))
    
    flnm='writestrctProps'
    WriteFull_flnm=trim(adjustl(flnm))//trim(adjustl(flnmbr))
    
    open(unit=21,file=trim(adjustl(ReadFull_flnm)))
    open(unit=22,file=trim(adjustl(WriteFull_flnm)))
    
    
    do i=1,N_itm
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       if (i==3) jmax=N_node
       if (i==4) jmax=N_mvCoordnte   
       if (i==5) jmax=N_node
       
       do j = 1,jmax
          if (i==1) read(21,*)  k_spr(j),l0(j)!,l(j)
          if (i==1) write(22,*) k_spr(j),l0(j)!,l(j)
          if (i==2) read(21,*)  k_area(j),A0(j)!,A(j)
          if (i==2) write(22,*) k_area(j),A0(j)!,A(j)
          if (i==3) read(21,*)  CgXNode(j),CgYNode(j)
          if (i==3) write(22,*) CgXNode(j),CgYNode(j)
          if (i==4) read(21,*)  coordntes_xy(j)
          if (i==4) write(22,*) coordntes_xy(j)
          if (i==5) read(21,*)  node_xy(j,1:N_dmnsn)
          if (i==5) write(22,*) node_xy(j,1:N_dmnsn)
       enddo
       
    enddo
    
    close(21)
    close(22)
    
  end subroutine read_strctProps_strct
  
  
  subroutine update_and_equilibrate(type_of_change,max_change,updt_cnt)
    !use calls_for_tests
    
    implicit none
    integer :: max_change
    integer :: type_of_change(1:max_change)
    
    integer :: which_typ_spr,typ_cnt
    real*8,  allocatable :: hm_area(:),hm_spr(:)
    integer, allocatable :: which_cells(:),which_sprs(:)
    integer, allocatable :: which(:),what(:)
    
    integer :: cell_cnt,spr_cnt
    real*8  :: hmMltplL(1:3),hmMltplK(1:3)
    real*8  :: max_prcnt
    
    integer :: spr_nmbr
    real*8  :: hm !how much
    integer :: hmn !how many nodes
    integer :: cell_nmbr
    integer :: spr_or_area
    
    integer :: j,m
    integer :: updt_cnt
    
    integer :: serial
  
    open(unit=78,file='FrameCheck.dat',position='append')
    
    max_prcnt = 0.0d0
    
    do j = 1,max_change
       
       
       if (type_of_change(j) .eq. 9) then
          write(*,*) updt_cnt,"updt_cnt"
        
          if (updt_cnt == 0) then
             call Equilibrate_and_Optimize_system
             call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
             Frame = Frame+1
             call switchto_NI_model_run_and_switchbackto_TN
             
          elseif (updt_cnt == 1) then
             write(78,fmt=*) Frame,"Frame bfr CgXNode"
             
             hm = 0.015d0
             call pulling_or_pushing_the_embryo(Exprmnt,Frame,hm)
             
             write(78,fmt=*) Frame,"Frame aft CgXNode"
             write(78,fmt=*) Frame,"Frame bfr A0 change"
             
             call successiveChangein_f   
             call changeA0_inSteps_Equilibrate_and_Optimize(Exprmnt,Frame)
             
             write(78,fmt=*) Frame,"Frame aft A0 change"
             write(78,fmt=*) Frame,"Frame bfr MidSpr"
             
             serial = 2
             hm = 0.15d0
             call stepwise_shortening_MidSpr_frmtop_and_Equilibrate(Exprmnt,Frame,serial,hm)
             write(78,fmt=*) Frame,"Frame aft MidSpr"
             
             if (stageNo==4 .and. stageType==1) then
                hmMltplL(1)=0.00d0; hmMltplL(2)= 0.40d0 ; hmMltplL(3)=0.30d0
                hmMltplK(1)=0.00d0; hmMltplK(2)= 0.00d0 ; hmMltplK(3)=0.00d0
                call Stepwise_SprShorteningIC_and_Equilibrate(Exprmnt,Frame,hmMltplL,hmMltplK)
             endif
             
             call symmetricize_prop
             call Equilibrate_system
             call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
             Frame=Frame+1
             call switchto_NI_model_run_and_switchbackto_TN
             write(78,fmt=*) Frame,"Frame aft Symm"
             
             k_phi(1:N_node,1:max_Phi_node) = 0.00d0
             call Equilibrate_system
             call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
             Frame=Frame+1
             call switchto_NI_model_run_and_switchbackto_TN
             write(78,fmt=*) Frame,"Frame aft k_phi"
             
             k_phi(1:N_node,1:max_Phi_node) = 0.0d0
             
          elseif (updt_cnt == 2) then
             
             !hm = 0.30d0
             !call shortening_the_pulley_sprs(hm)
             !call shortening_the_pulling_sprs(hm)
             !call nodesl0A0_to_coordntes(node_xy,l0,A0,coordntes)
             
             !call switch_off_alpha_beta_excpt_bottom4_cells
             !call minimize_activePairSprs
             continue
             
             
          elseif (updt_cnt == 3) then
             !hm = 0.30d0
             !call shortening_basalFurrow_pair2(hm)
             !call nodesl0A0_to_coordntes(node_xy,l0,A0,coordntes)
             continue
             
          elseif (updt_cnt == 4) then
             
             ! call switch_on_alphabeta_ofFirstPair_ofbottomCell
             
             !call get_trnsfrng_cells(trns_frmCell,trns_toCell)
             !max_prcnt=0.9001d0
             !call TrgtPropChng_and_equilibrate(trns_frmCell,trns_tocell,max_prcnt)
             !call smallChnge_and_equilibrate(trns_frmCell,trns_tocell,max_prcnt)
             !max_prcnt=0.4001d0
             !call fctr_and_l0s_frstPairCntrlBot_cells(max_prcnt)
             
          endif
          
          write(*,*) l0(31),l0(32),l0(34),l0(35),"l0_aft minimizing PenF"
          
          
          
       elseif (type_of_change(j) .eq. 1) then
          cell_cnt = 8
          allocate(which_cells(1:cell_cnt),hm_area(1:cell_cnt))
          which_cells(1) = 5 ; which_cells(2) = 6
          which_cells(3) = 7 ; which_cells(4) = 8
          which_cells(5) = 9 ; which_cells(6) = 10
          which_cells(7) = 11; which_cells(8) = 12
          
          hm_area(1:2) = -0.40d0
          hm_area(3:4) = 0.00d0
          hm_area(5:6) = 3.50d0
          hm_area(7:8) = 3.60d0
          
          call change_A0_of_mltipl_cell(which_cells,cell_cnt,hm_area)
          write(*,*) A0(7),A0(8),A0(9),A0(10),"A0 of 7-10"
          deallocate(which_cells,hm_area)
          
          
       elseif (type_of_change(j) .eq. 2) then
          cell_cnt = 8
          allocate(which_cells(1:cell_cnt),hm_area(1:cell_cnt))
          which_cells(1) = 5 ; which_cells(2) = 6
          which_cells(3) = 7 ; which_cells(4) = 8
          which_cells(5) = 9 ; which_cells(6) = 10
          which_cells(7) = 11; which_cells(8) = 12
          
          hm_area(1:2) = 0.00d0
          hm_area(3:4) = 0.00d0
          hm_area(5:6) = 0.00d0
          hm_area(7:8) = 5.00d0
          
          call change_k_area_of_mltipl_cell(which_cells,cell_cnt,hm_area)
          write(*,*) k_area(7),k_area(8),k_area(9),k_area(10),"k_area of 7-10"
          
          deallocate(which_cells,hm_area)
          
          
       elseif (type_of_change(j) .eq. 3) then
          spr_cnt = nsecclb + nseccrb + 2
          allocate(which_sprs(1:spr_cnt),hm_spr(1:spr_cnt))
          which_sprs(1) = 21 ; which_sprs(2) = 24
          which_sprs(3) = 25 ; which_sprs(4) = 28
          which_sprs(5) = 26 ; which_sprs(6) = 29
          which_sprs(7) = 27 ; which_sprs(8) = 30
          hm_spr = -0.30d0
          call change_l0_of_mltipl_spr(which_sprs,spr_cnt,hm_spr)
          write(*,*) l0(21),l0(30),"l0 21,30"
          
          cell_nmbr = 9
          call get_k_sprs_of_rectangular_shape(cell_nmbr)
          cell_nmbr = 10
          call get_k_sprs_of_rectangular_shape(cell_nmbr)
          
          write(*,*) k_spr(25),k_spr(26),k_spr(21),k_spr(27),"k_spr 25,30"
          deallocate(which_sprs,hm_spr)
          
          
       elseif (type_of_change(j) .eq. 4) then
          
          spr_nmbr = 19 !pulling one
          hm = -0.28d0
          call change_l0_of_one_spr(spr_nmbr,hm)
          write(*,*) l0(spr_nmbr),spr_nmbr,"l0 of spr_nmbr"
          hm = 1.35d0
          call change_k_spr_of_one_spr(spr_nmbr,hm)
          write(*,*) k_spr(spr_nmbr),spr_nmbr,"k_spr of spr_nmbr"
          
          
          spr_nmbr = 22 !pulling one
          hm = -0.28d0
          call change_l0_of_one_spr(spr_nmbr,hm)
          write(*,*) l0(spr_nmbr),spr_nmbr,"l0 of spr_nmbr"
          hm = 1.35d0
          call change_k_spr_of_one_spr(spr_nmbr,hm)
          write(*,*) k_spr(spr_nmbr),spr_nmbr,"k_spr of spr_nmbr"
          
          
       elseif (type_of_change(j) .eq. 5) then
          which_typ_spr = 6 !pulley_spr
          call get_spr_typ_cnt(which_typ_spr,typ_cnt)
          write(*,*) typ_cnt,"typ_cnt6"
          allocate(hm_spr(1:typ_cnt))
          
          hm_spr = 2.00d0 !0.90
          call change_k_spr_of_a_typ_of_spr(which_typ_spr,typ_cnt,hm_spr)
          write(*,*) k_spr(13),k_spr(16),"k_spr of TYPE 6"
          
          hm_spr = -0.90d0    
          call change_l0_of_a_typ(which_typ_spr,typ_cnt,hm_spr)   
          write(*,*) l0(13),l0(16),"l0 of TYPE 6_aft"
          
          deallocate(hm_spr)
          
       elseif (type_of_change(j) .eq. 6) then
          spr_cnt = 10
          allocate(which_sprs(1:spr_cnt))
          which_sprs(1) = 21 ; which_sprs(2) = 24
          which_sprs(3) = 25 ; which_sprs(4) = 28
          which_sprs(5) = 26 ; which_sprs(6) = 29
          which_sprs(7) = 27 ; which_sprs(8) = 30
          which_sprs(9) = 19 ; which_sprs(10) = 22
          !call switch_off_all_sprs_except_some(spr_cnt,which_sprs)
          deallocate(which_sprs)
          write(*,*) k_spr(20),k_spr(21),k_spr(25),k_spr(29),k_spr(23),"k_sprs 20-21-25-29-23"
          
          
       elseif (type_of_change(j) .eq. 7) then
          cell_cnt = 4
          allocate(which_cells(1:cell_cnt))
          which_cells(1) = 7 ; which_cells(2) = 8
          which_cells(3) = 9 ; which_cells(4) = 10
          !call switch_off_all_areas_except_some(cell_cnt,which_cells)
          deallocate(which_cells)
          write(*,*) A0(7),A0(8),A0(9),A0(10),"A0 7-10"
          
       elseif (type_of_change(j) .eq. 8) then
          spr_nmbr = 25 !pulling one
          hm = -0.28d0
          call change_l0_of_one_spr(spr_nmbr,hm)
          write(*,*) l0(spr_nmbr),spr_nmbr,"l0 of spr_nmbr"
          hm = 1.35d0
          call change_k_spr_of_one_spr(spr_nmbr,hm)
          write(*,*) k_spr(spr_nmbr),spr_nmbr,"k_spr of spr_nmbr"
          
          
          spr_nmbr = 28 !pulling one
          hm = -0.28d0
          call change_l0_of_one_spr(spr_nmbr,hm)
          write(*,*) l0(spr_nmbr),spr_nmbr,"l0 of spr_nmbr"
          hm = 1.35d0
          call change_k_spr_of_one_spr(spr_nmbr,hm)
          write(*,*) k_spr(spr_nmbr),spr_nmbr,"k_spr of spr_nmbr"
          
       endif
       
    enddo
    
    write(*,*) type_of_change(8),"INSD the updt_params"
    
    close(78)
    
  end subroutine update_and_equilibrate


end module calls_for_tests

subroutine update_node_types(htu)
  use calls_for_tests
  
  implicit none
  integer :: htu
  
  integer :: hmt
  integer :: which_type,changed_type
  integer,allocatable :: which_types(:),changed_types(:)
 
  
  if (htu .eq. 1) then
     which_type   = 2
     changed_type = 1
     
     call store_and_change_node_type(which_type,changed_type)

  elseif (htu .eq. 2) then
     hmt = 2
     allocate(which_types(1:hmt),changed_types(1:hmt))
     which_types(1)= 0 ; which_types(2)= 2
     changed_types(1:2) = 1
     
     call store_and_change_mltipl_node_types(hmt,which_types,changed_types)
     
     write(*,*) node_typ(1),node_typ(3),node_typ(15),node_typ(7),"node_typ-3,5"
     write(*,*) node_typ(3),node_typ(4),"node_typ-3,4"

  endif
     

end subroutine update_node_types





subroutine testing_gradient
  use calls_for_tests
  
  implicit none
  integer :: spr_or_area

  spr_or_area = 2
  
  if (spr_or_area .eq. 1) then
     call test_for_spr
  elseif (spr_or_area .eq. 2) then    
     call test_for_area
  endif
  

contains


  subroutine test_for_spr
    implicit none

    integer :: i
    real*8,allocatable :: k_spr_tmp(:)

    allocate(k_spr_tmp(1:N_spr))

    !k_spr_tmp = 0.0d0
    k_spr_tmp  = k_spr
    k_spr      = 0.0d0
    k_area_tmp = k_area
    k_area     = 0.0d0
    
    do i = 17,18
       if (i.eq.18) exit
       write(*,*)"SPRING",i,"is ON"
       k_spr(i) = k_spr_tmp(i)
       write(*,*) k_spr(i),l0(i),"k_spr and l0 for i = ",i
       call get_gradient(node_xy,l0,A0,gradient)
       !write(*,*) k_spr,"bfr"
       k_spr(i) = 0.0d0
       !write(*,*) k_spr,"aft"
      
       !write(*,*) gradient(14,1:2)
       !write(*,*) gradient(13,1:2)
    enddo
    
    k_spr  = k_spr_tmp 
    k_area = k_area_tmp
    
    !write(*,*)
  end subroutine test_for_spr


  subroutine test_for_area
    implicit none

    integer :: i
    real*8,allocatable :: k_area_tmp(:)

    allocate(k_area_tmp(1:N_cell))
    
    k_spr_tmp  = k_spr
    k_spr      = 0.0d0
    k_area_tmp = k_area
    k_area     = 0.0d0
    
    do i = 1,N_cell
       write(*,*) i,"Cell_bfr"
       if (i.ne.4) cycle
       write(*,*) i,"Cell_aft"
       
       write(*,*)"AREA",i,"is ON"
       
       k_area(i) = k_area_tmp(i)
       write(*,*) k_area(i),A0(i),"k_area and A0 for i = ",i
       
       call get_gradient(node_xy,l0,A0,gradient)
       !write(*,*) k_area,"bfr"
       k_area(i) = 0.0d0
       !write(*,*) k_area,"aft"
    enddo

    k_spr  = k_spr_tmp 
    k_area = k_area_tmp
    
  end subroutine test_for_area

  
end subroutine testing_gradient


subroutine testing_aftr_A_frame_and_saving_result
  use calls_for_tests
  
  implicit none

  real*8 :: l0_frm(1:N_spr)
  real*8 :: A0_frm(1:N_cell)
  real*8 :: coordntesXY_frm(1:N_mvCoordnte)


  l0_frm  = l0
  A0_frm  = A0
  coordntesXY_frm = coordntes_xy

  !Test what you want to test 


  !call save_config_and_generate_data()
  


end subroutine testing_aftr_A_frame_and_saving_result
