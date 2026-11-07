module area_to_other_transfrm
  use system_parameters
  use transfrm_info
  
  implicit none
  integer :: cnn_lft,cnn_rght
  integer :: nodesRingF=6
  integer :: nodesRingA=7  !nodes in Ring, F=First,A=All other excpt clt/crt
  integer :: nodesRingCT=8 !nodes in Ring for Cntrl Top (lft/rght) 
  integer :: nodesRingER=6 !nodes in Ring for End Rgn
  integer :: nodesRingCR   !nodes in Ring for End Rgn
  
  integer :: nodesRingF_withAS=7 !AS=Apicalside [wont work with pulley concept]
  integer :: nodesRingA_withAS=8
  
contains
  
  subroutine get_area_to_spr_transfrm
    implicit none
    
    if (modelID==1) then
       call get_area_to_spr_transfrm_TN_model
    elseif (modelID==2) then
       call get_area_to_spr_transfrm_NI_model
    endif
    
  end subroutine get_area_to_spr_transfrm
  
  subroutine get_area_to_spr_transfrm_TN_model
    implicit none
    integer :: i1
    integer :: btbb !btbb=blocks_to_be_built
    
    count_area = 0
    
    btbb = wntbb(0)
    write(*,*) btbb,"btbb"
    
    do i1 = 1,max_blck
       if (wntbb(i1) == 1) then
          call lft_area_to_spr
       elseif (wntbb(i1) == 2) then
          call rght_area_to_spr
       elseif (wntbb(i1) == 3) then
          call clft_top_area_to_spr
       elseif (wntbb(i1) == 4) then   
          !call adjstment_for_additnl_spr_in_area_spr
          !not needed, if needed will adjust later
          continue
       elseif (wntbb(i1) == 5) then
          call crght_top_area_to_spr
       elseif (wntbb(i1) == 6) then
          call clft_bot_area_to_spr
       elseif (wntbb(i1) == 7) then
          call crght_bot_area_to_spr
       elseif (wntbb(i1)==8) then
          call cntrlCell_area_to_spr
       elseif (wntbb(i1)==9) then
          call endRgnCell_area_to_spr
       else
          continue
       endif
       
    enddo
    
    call print_area_to_spr_transfrm
    
  contains
    
    subroutine lft_area_to_spr
      implicit none
      integer :: i
      
      do i = 1,ncl !no_of_cells_lft
         count_area = count_area + 1
         
         if (i==1) then
            area_spr(count_area,0) = 3
            area_spr(count_area,1) = 3*count_area-2
            area_spr(count_area,2) = 3*count_area-1
            area_spr(count_area,3) = 3*count_area
            
         else
            area_spr(count_area,0) = 4
            area_spr(count_area,1) = 3*count_area-3
            area_spr(count_area,2) = 3*count_area-2
            area_spr(count_area,3) = 3*count_area-1
            area_spr(count_area,4) = 3*count_area
         endif

      enddo
      
    end subroutine lft_area_to_spr
    
    subroutine rght_area_to_spr
      implicit none
      integer :: i
      
      do i = 1,ncr
         count_area = count_area + 1
         
         if (i == 1) then
            area_spr(count_area,0) = 3
            area_spr(count_area,1) = 3*count_area-2
            area_spr(count_area,2) = 3*count_area-1
            area_spr(count_area,3) = 3*count_area
         else
            area_spr(count_area,0) = 4
            area_spr(count_area,1) = 3*count_area-3
            area_spr(count_area,2) = 3*count_area-2
            area_spr(count_area,3) = 3*count_area-1
            area_spr(count_area,4) = 3*count_area
         endif
      enddo
      
    end subroutine rght_area_to_spr

    subroutine cntrlCell_area_to_spr
      implicit none
      integer :: i
      
      
      do i = 1,ncc
         count_area = count_area+1

         if (i==1) then
            area_spr(count_area,0) = 4
            area_spr(count_area,1) = lft_endSpring(1)
            area_spr(count_area,2) = 3*count_area-2
            area_spr(count_area,3) = 3*count_area-1
            area_spr(count_area,4) = 3*count_area-3
         endif
      
      enddo
      
    end subroutine cntrlCell_area_to_spr
    
    subroutine clft_top_area_to_spr
      implicit none
      
      count_area = count_area + 1
      
      area_spr(count_area,0) = 4
      area_spr(count_area,1) = lft_endSpring(1)
      area_spr(count_area,2) = clft_top_endSpring(1)
      area_spr(count_area,3) = clft_top_endSpring(1) + 1 
      area_spr(count_area,4) = clft_top_endSpring(2)
      
    end subroutine clft_top_area_to_spr
    
    ! subroutine additnl_area_to_spr
    !   implicit none
    !   integer :: i,j
    
    ! end subroutine additnl_area_to_spr

    subroutine adjstment_for_additnl_spr_in_area_spr
      implicit none
      integer :: i
      integer :: cnnctng_area
      
      call get_cnnctng_area(wntbb(i-1),cnnctng_area)

    end subroutine adjstment_for_additnl_spr_in_area_spr
    
    subroutine crght_top_area_to_spr
      implicit none
      
      count_area = count_area + 1
      
      area_spr(count_area,0) = 4
      area_spr(count_area,1) = rght_endSpring(1)
      area_spr(count_area,2) = crght_top_endSpring(1)
      area_spr(count_area,3) = crght_top_endSpring(1) + 1 
      area_spr(count_area,4) = crght_top_endSpring(2)
      
    end subroutine crght_top_area_to_spr

    subroutine clft_bot_area_to_spr
      implicit none
      integer :: i

      count_area = ncl + ncr + ncclt + nccrt
      
      do i = 1,ncclb
         
         if (i==1) then
            count_area = count_area + 1
         else
            count_area = count_area + 2
         endif
         
         area_spr(count_area,0) = 4

         if (i==1) then
            area_spr(count_area,1) = clft_top_endSpring(1)+1
         else
            area_spr(count_area,1) = 3*count_area-6
         endif
         
         area_spr(count_area,2) = 3*count_area-2
         area_spr(count_area,3) = 3*count_area-1
         area_spr(count_area,4) = 3*count_area
        
      enddo
      
    end subroutine clft_bot_area_to_spr
    
    
    subroutine crght_bot_area_to_spr
      implicit none
      integer :: i
      
      count_area = ncl + ncr + ncclt + nccrt + 1
      
      do i = 1,nccrb
         
         if (i==1) then
            count_area = count_area + 1
         else
            count_area = count_area + 2
         endif
         
         area_spr(count_area,0) = 4
         
         if (i==1) then
            area_spr(count_area,1) = crght_top_endSpring(1)+1
         else
            area_spr(count_area,1) = 3*count_area-6
         endif
         
         area_spr(count_area,2) = 3*count_area-2
         area_spr(count_area,3) = 3*count_area-1
         area_spr(count_area,4) = 3*count_area
         
      enddo
      
    end subroutine crght_bot_area_to_spr
    
    
    subroutine endRgnCell_area_to_spr
      implicit none
      integer :: i
      
      write(*,*) clft_bot_endSpring(1),crght_bot_endSpring(1),endCell_Spring(1),"endSpr,clft-rght"
      
      do i = 1,ncer
         count_area = count_area+1

         if (i==1) then
            area_spr(count_area,0) = 3
            area_spr(count_area,1) = clft_bot_endSpring(1) 
            area_spr(count_area,2) = endCell_Spring(1)
            area_spr(count_area,3) = crght_bot_endSpring(1)
         endif
         
      enddo
      
      
    end subroutine endRgnCell_area_to_spr
    
    
  end subroutine get_area_to_spr_transfrm_TN_model !!!!!!
  
  
  
  subroutine get_area_to_spr_transfrm_NI_model
    implicit none
    integer :: i,j,m,jmax,mmax
    integer :: area_nm,spr_nm
    integer :: cnt_curr
    
    
    do i = 1,N_cell
       area_nm = i 
       jmax    = area_sprS(area_nm,0)
       cnt_curr = 0
       
       do j = 1,jmax
          spr_nm = area_sprS(area_nm,j)
          !write(*,*) spr_nm,j,"spr_nm,j"
          
          if (trmnl_intrmdSpr(spr_nm,0)==1) then
             cnt_curr = cnt_curr+1
             area_spr(area_nm,cnt_curr) = trmnl_intrmdSpr(spr_nm,1)
             !write(*,*) area_spr(area_nm,cnt_curr),cnt_curr,"E1"
             
          elseif (trmnl_intrmdSpr(spr_nm,0).ne.1) then
             mmax = trmnl_intrmdSpr(spr_nm,0)
             
             do m = 1,mmax
                cnt_curr = cnt_curr+1
                area_spr(area_nm,cnt_curr) =  trmnl_intrmdSpr(spr_nm,m)
                !write(*,*) area_spr(area_nm,cnt_curr),cnt_curr,"E2"
             enddo
             
          endif
          
       enddo
       
       area_spr(area_nm,0) = cnt_curr
       !write(*,*) area_spr(area_nm,0:max_spr_area),i,"areaspr"
       
    enddo
    
    call print_area_to_spr_transfrm
    
  end subroutine get_area_to_spr_transfrm_NI_model
  
  
  subroutine print_area_to_spr_transfrm
    implicit none
    integer :: i
    
    open(unit=223,file='area_spr.dat')
    open(unit=224,file='area_sprWOtxt.dat')
    
    write(223,*) "area-spr","    ","first one node_area(i,0)"
    
    !write(*,*) N_cell,"area_to_spr"
    
    do i = 1,N_cell
       write(223,*)"area to spr for i =",i,"area are",area_spr(i,0:max_spr_area)
       write(224,*) i,area_spr(i,0:max_spr_area)
    enddo
    
    close(223)
    close(224)
    
  end subroutine print_area_to_spr_transfrm
  
  
  subroutine get_area_to_node_transfrm
    implicit none
    
    if (modelID==1) then
       call get_area_to_node_transfrm_TN_model
    elseif (modelID==2) then
       call get_area_to_node_transfrm_NI_model
    endif
    
  end subroutine get_area_to_node_transfrm

  subroutine get_area_to_node_transfrm_TN_model
    implicit none
    integer :: i,lp_cnt
    integer :: count_area,count_nodes

    count_area  = 0
    count_nodes = 0
    
    lp_cnt = wntbb(0)

    do i = 1,max_blck
       
       if (wntbb(i) == 1) then
          call lft_area_to_node
       elseif (wntbb(i) == 2) then
          call rght_area_to_node
       elseif (wntbb(i) == 3) then  
          call clft_top_area_to_node
       elseif (wntbb(i) == 4) then
          continue
       elseif (wntbb(i) == 5) then
          call crght_top_area_to_node
       elseif (wntbb(i) == 6) then
          call clft_bot_area_to_node
       elseif (wntbb(i) == 7) then
          call crght_bot_area_to_node
       elseif (wntbb(i) == 8) then
          call cntrlCell_area_to_node
       elseif (wntbb(i) == 9) then
          call endRgnCell_area_to_node
       else
          continue
       endif
       
    enddo
    
    call print_area_to_node_trnsfrm
    
    
  contains
    
    subroutine lft_area_to_node
      implicit none
      integer :: i,j
      integer :: nnfbcl(1:nnecl) !nnfbcl=node_nos_frst_blck_cell_lft
      
      !area_node(1:ncl,0) = nnecl !=4
      
      nnfbcl(1) = 1 ; nnfbcl(2) = 2   !ANTI-CLOCKWISE
      nnfbcl(3) = 4 ; nnfbcl(4) = 3
      
      do i = 1,ncl
         count_area = count_area + 1
         area_node(count_area,0) = nnecl !=4
         
         do j = 1,nnecl
            area_node(count_area,j) = nnfbcl(j) + (i-1)*2 ! i here,not count_area in bracket
         enddo
         !write(*,*) area_node(i,1:max_node_area),"Area,i",i
         
      enddo

      count_nodes = count_nodes + nvsl*2
      
    end subroutine lft_area_to_node
    
    subroutine rght_area_to_node
      implicit  none
      integer :: i,j
      integer :: nnfbcr(1:nnecr) !nnfbcr=node_nos_frst_blck_cell_rght
      
      
      nnfbcr(1) = count_nodes + 1 ; nnfbcr(2) = count_nodes + 3  !reverse in lft to maintain anti_clockwise ness
      nnfbcr(3) = count_nodes + 4 ; nnfbcr(4) = count_nodes + 2
      
      !area_node(start:finish,0) = nnecr !=4
      
      do i = 1,ncr
         count_area = count_area + 1
         area_node(count_area,0) = nnecr !=4
         
         do j = 1,nnecr
            area_node(count_area,j) = nnfbcr(j) + (i-1)*2
         enddo
         !write(*,*) area_node(i,1:max_node_area),"Area,i",i
      enddo
      
      count_nodes = count_nodes + nvsr*2
      
    end subroutine rght_area_to_node
    
    subroutine cntrlCell_area_to_node
      implicit none
      integer :: i,j
      
      write(*,*) rght_endNode(0:2),lft_endNode(0:2),nnecc,"rght_endNode,lft_endNode,nnecc"
      
      do i = 1,ncc
         count_area = count_area+1
         
         area_node(count_area,0) = nnecc
            
         do j = 1,nnecc   
            if (j==1) area_node(count_area,j) = rght_endNode(1)
            if (j==2) area_node(count_area,j) = lft_endNode(1)
            if (j==3) area_node(count_area,j) = lft_endNode(2)
            if (j==4) area_node(count_area,j) = rght_endNode(2)
         enddo
         
      enddo
      
      
      
    end subroutine cntrlCell_area_to_node
    
    subroutine clft_top_area_to_node
      implicit none
      !integer :: i
      integer :: j

      !this routine is okay for nccl=1; if nccl > 1, then I have to rewrite it in a loopwise manner.Same thing for the crght_top_area_to_node!!!

      
      count_area = count_area + ncclt
      !write(*,*) nneccl,"nneccl"
      area_node(count_area,0) = nnecclt !=5
      
      do j = 1,nnecclt
         if (j == 1) then
            area_node(count_area,j) = lft_endNode(2)
         elseif (j.gt.1 .AND. j.lt.nnecclt) then
            area_node(count_area,j) = clft_top_endNode(3) - (j-2)!clft_top_end_n-(0/1/2)
         else
            area_node(count_area,j) = lft_endNode(1)
         endif
      enddo
      
    end subroutine clft_top_area_to_node
    
    !NO AREA FOR ADDITNL NODE
    
    subroutine crght_top_area_to_node
      implicit none
      !integer :: i
      integer :: j
      
      count_area = count_area + nccrt      
      area_node(count_area,0) = nneccrt
      
      do j = 1,nnecclt
         if (j == 1) then
            area_node(count_area,j) = rght_endNode(1)
         elseif (j.gt.1 .AND. j.lt.nnecclt) then
           area_node(count_area,j) = crght_top_endNode(1) + (j-2)
         else
            area_node(count_area,j) = rght_endNode(2)
         endif
      enddo
      
    end subroutine crght_top_area_to_node



    subroutine clft_bot_area_to_node
      implicit none
      integer :: i
      
      count_area = ncl + ncr + ncclt + nccrt
      
      do i = 1,ncclb
         
         if (i==1) then
            count_area = count_area + 1
         else
            count_area = count_area + 2
         endif
         
         area_node(count_area,0) = 4
         
         if (i==1) then
            area_node(count_area,1) = clft_top_endNode(2) 
            area_node(count_area,2) = clft_top_endNode(3)
            area_node(count_area,3) = 8 + 2*(count_area-1)
            area_node(count_area,4) = 8 + 2*(count_area-1) - 1
            write(*,*) area_node(count_area,3:4),count_area,"area_node"
            !stop
         else
            area_node(count_area,1) = 3 + 2*(count_area-1)
            area_node(count_area,2) = 3 + 2*(count_area-1) + 1
            area_node(count_area,3) = 7 + 2*(count_area-1) + 1
            area_node(count_area,4) = 7 + 2*(count_area-1)
            
         endif
         
      enddo
      
    end subroutine clft_bot_area_to_node

    

    subroutine crght_bot_area_to_node
      implicit none
      integer :: i
      
      count_area = ncl + ncr + ncclt + nccrt + 1
      
      do i = 1,nccrb
         
         if (i==1) then
            count_area = count_area + 1
         else
            count_area = count_area + 2
         endif
         
         area_node(count_area,0) = 4
         
         if (i==1) then
            area_node(count_area,1) = crght_top_endNode(3) 
            area_node(count_area,2) = crght_top_endNode(2)
            area_node(count_area,3) = 8 + 2*(count_area-1) - 1
            area_node(count_area,4) = 8 + 2*(count_area-1)

         else
            area_node(count_area,1) = 4 + 2*(count_area-1)
            area_node(count_area,2) = 4 + 2*(count_area-1) - 1
            area_node(count_area,3) = 8 + 2*(count_area-1) - 1
            area_node(count_area,4) = 8 + 2*(count_area-1)
            
         endif
         
      enddo
      

    end subroutine crght_bot_area_to_node
    

    subroutine endRgnCell_area_to_node
      implicit none
      integer :: i,j

      write(*,*)clft_bot_endNode(1:2),crght_bot_endNode(1:2),"endNode,clft-rght"
      
      
      do i = 1,ncer
         count_area = count_area+1
         
         area_node(count_area,0) = nnecer
         
         do j = 1,nnecer   
            if (j==1) area_node(count_area,j) = clft_bot_endNode(1)
            if (j==2) area_node(count_area,j) = clft_bot_endNode(2)
            if (j==3) area_node(count_area,j) = crght_bot_endNode(2)
         enddo
         
      enddo
      
    end subroutine endRgnCell_area_to_node
    
    
  end subroutine get_area_to_node_transfrm_TN_model


  subroutine get_area_to_node_transfrm_NI_model
    implicit none
    integer :: i,lp_cnt
    integer :: count_area,count_nodes
    integer :: cnt_addNode,cnt_rglNode
    
    count_area  = 0
    cnt_rglNode = 0
    cnt_addNode = N_nodeS
    
    lp_cnt = wntbb(0)
    
    do i = 1,max_blck
       
       if (wntbb(i) == 1) then
          call lft_area_to_node_NI
       elseif (wntbb(i) == 2) then
          call rght_area_to_node_NI
       elseif (wntbb(i) == 3) then  
          call clft_top_area_to_node_NI
       elseif (wntbb(i) == 4) then
          continue
       elseif (wntbb(i) == 5) then
          call crght_top_area_to_node_NI
       elseif (wntbb(i) == 6) then
          call clft_bot_area_to_node_NI
       elseif (wntbb(i) == 7) then
          call crght_bot_area_to_node_NI
       elseif (wntbb(i)==8) then
          call cntrlRgn_area_to_node_NI
       elseif (wntbb(i)==9) then
          call endRgn_area_to_node_NI
       endif
       
    enddo
    
    if (SystemTyp.ne.0) then
       call adjst_area_to_node_trnsfrm_NI
    endif
    
  contains
    
    subroutine lft_area_to_node_NI
      implicit none
      integer :: i,j
      integer :: nnfbcl(1:max_node_area,1:2) !nnfbcl=node_nos_frst_blck_cell_lft
      integer :: area_indctr
      
      do i = 1,ncl
         call get_nnfbcl(i,nnfbcl)
         count_area = count_area + 1
         
         if (i==1) then
            area_node(count_area,0) = max_node_area-1-1*(NAEC)
         elseif (i.gt.1 .and. i.le.ncl) then
            area_node(count_area,0) = max_node_area-1
         endif
         
         do j = 1,(max_node_area-1)
            if (nnfbcl(j,1).ne.(-1)) then
               
               if (nnfbcl(j,1)==1) then !rgl node
                  area_node(count_area,j) = nnfbcl(j,2) + (i-1)*2
                  
               elseif (nnfbcl(j,1)==2) then !insrtd node
                  area_node(count_area,j) = nnfbcl(j,2) + (i-1)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)
               endif
               
            endif
         enddo
         !write(*,*) area_node(i,1:max_node_area),"Area,i",i

         if (area_nodeS(count_area,0)==5) then
            area_indctr=1
            call adjstmnt_area_to_node_for_S1T1(area_indctr,count_area)
         endif
         
      enddo
      
      cnt_rglNode = cnt_rglNode + nvsl*2
      cnt_addNode = cnt_addNode + ncl*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)
      
      lft_endNode(0)=3 ; lft_endNode(3) = cnt_addNode
      write(*,*) cnt_rglNode,cnt_addNode,"rgl,add"
      
    end subroutine lft_area_to_node_NI
    
    subroutine rght_area_to_node_NI
      implicit  none
      integer :: i,j
      integer :: nnfbcr(1:max_node_area,1:2) !nnfbcr ~ nnfbcl
      integer :: area_indctr
      
      write(*,*) cnt_rglNode,cnt_addNode,"rgl,add"
      
      do i = 1,ncr
         call get_nnfbcr(i,cnt_rglNode,cnt_addNode,nnfbcr)
         
         count_area = count_area + 1
         
         if (i==1) then
            area_node(count_area,0) = max_node_area-1-1*(NAEC)
         elseif (i.gt.1 .and. i.le.ncr) then
            area_node(count_area,0) = max_node_area-1
         endif
         
         do j = 1,(max_node_area-1)
            if (nnfbcr(j,1).ne.(-1)) then
               
               if (nnfbcr(j,1)==1) then !rgl node
                  area_node(count_area,j) = nnfbcr(j,2) + (i-1)*2
                  
               elseif (nnfbcr(j,1)==2) then !insrtd node
                  area_node(count_area,j) = nnfbcr(j,2) + (i-1)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)
               endif
               
            endif
         enddo
         !write(*,*) area_node(i,1:max_node_area),"Area,i",i
         
         if (area_nodeS(count_area,0)==5) then
            area_indctr=2
            call adjstmnt_area_to_node_for_S1T1(area_indctr,count_area)
         endif
         
      enddo
      
      cnt_rglNode = cnt_rglNode + nvsr*2
      cnt_addNode = cnt_addNode + ncr*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)
      
      rght_endNode(0) = 3 ; rght_endNode(3) = cnt_addNode
      
    end subroutine rght_area_to_node_NI
    
    subroutine clft_top_area_to_node_NI
      implicit none
      integer :: j,k
      integer :: grp1(NAEC),grp2(1:NAEC),grp3(1:NAEC)
      integer :: cnt_N
      
      count_area = count_area + ncclt
      area_node(count_area,0) = max_node_area

      cnt_N=1
      
      call clft_top_node_groups(cnt_addNode,grp1,grp2,grp3)
      
      do j = 1,nodesRingCT
         
         if (j==1) then
            area_node(count_area,cnt_N) = lft_endNode(2)
            cnt_N=cnt_N+1
            
         elseif (j==2) then
            
            do k=1,NAEC
               area_node(count_area,cnt_N) = grp1(k)
               cnt_N=cnt_N+1
            enddo
            !area_node(count_area,cnt_N) = cnt_addNode+2
            
         elseif (j==3) then
            area_node(count_area,cnt_N) = clft_top_endNode(3) -(j-3)
            cnt_N=cnt_N+1
            
         elseif (j==4) then
            
            do k=1,NAEC
               area_node(count_area,cnt_N) = grp2(k)
               cnt_N=cnt_N+1
            enddo
            !area_node(count_area,cnt_N) = cnt_addNode+1
            
         elseif (j==5.or.j==6) then
            area_node(count_area,cnt_N) = clft_top_endNode(3) -(j-4)
            cnt_N=cnt_N+1
            
         elseif (j == (nodesRingCT-1)) then
            
            area_node(count_area,cnt_N) = lft_endNode(1)
            cnt_N=cnt_N+1
            
         elseif (j == nodesRingCT) then
            
            do k=1,NAEC
               area_node(count_area,cnt_N) = grp3(k)
               cnt_N=cnt_N+1
            enddo
            !area_node(count_area,cnt_N) = lft_endNode(3)
            
         endif
      enddo
      
      cnn_lft     = cnt_addNode + 1
      cnt_rglNode = cnt_rglNode + 3
      cnt_addNode = cnt_addNode + 2*(NAEC)
      
      lft_endNode(3) = 0
      
    end subroutine clft_top_area_to_node_NI
    
    
    subroutine crght_top_area_to_node_NI
      implicit none
      integer :: j,k
      integer :: grp1(NAEC),grp2(1:NAEC),grp3(1:NAEC)
      integer :: cnt_N
      
      count_area = count_area + nccrt
      area_node(count_area,0) = max_node_area

      cnt_N=1
      
      call crght_top_node_groups(cnt_addNode,grp1,grp2,grp3)
      
      do j = 1,nodesRingCT
         
         if (j==1) then
            area_node(count_area,cnt_N) = rght_endNode(1)
            cnt_N=cnt_N+1
            
         elseif (j==2.or.j==3) then
            area_node(count_area,cnt_N) = crght_top_endNode(1) + (j-2)
            cnt_N=cnt_N+1
            
         elseif (j==4) then
            
            do k=1,NAEC
               area_node(count_area,cnt_N) = grp1(k)
               cnt_N=cnt_N+1
            enddo
            !area_node(count_area,j) = cnt_addNode+1
            
         elseif (j==5) then
            area_node(count_area,cnt_N) = crght_top_endNode(1) + (j-3)
            cnt_N=cnt_N+1
            
         elseif (j==6) then
            
            do k=1,NAEC
               area_node(count_area,cnt_N) = grp2(k)
               cnt_N=cnt_N+1
            enddo
            !area_node(count_area,j) = cnt_addNode+2
            
         elseif (j==(nodesRingCT-1)) then
            area_node(count_area,cnt_N) = rght_endNode(2)
            cnt_N=cnt_N+1
            
         elseif (j==nodesRingCT) then
            
            do k=1,NAEC
               area_node(count_area,cnt_N) = grp3(k)
               cnt_N=cnt_N+1
            enddo
            !area_node(count_area,j) = rght_endNode(3)
            
         endif
      enddo
      
      cnn_rght    = cnt_addNode + 1
      
      cnt_rglNode = cnt_rglNode + 3
      cnt_addNode = cnt_addNode + 2*(NAEC)
      
      rght_endNode(3) = 0
      
    end subroutine crght_top_area_to_node_NI
    
    subroutine clft_bot_area_to_node_NI
      implicit none
      integer :: i,j,k
      integer :: cnt_N
      integer :: grp1(1:NAEC),grp2(1:NAEC),grp3(1:NAEC)
      
      count_area = ncl + ncr + ncclt + nccrt
      
      do i = 1,ncclb
         
         if (i==1) then
            count_area = count_area + 1
            area_node(count_area,0) = max_node_area-1
            cnt_rglNode = clft_top_endNode(3)
            cnt_addNode = cnn_lft-1+ 2*(NAEC)
            
         elseif (i==2) then
            count_area = count_area + 2
            area_node(count_area,0) = max_node_area-1
            cnt_rglNode = cnt_rglNode+5
            cnt_addNode = cnt_addNode+4*(NAEC)
            
         elseif (i.gt.2) then
            count_area = count_area + 2
            area_node(count_area,0) = max_node_area-1
            cnt_rglNode = cnt_rglNode+4
            cnt_addNode = cnt_addNode+4*(NAEC)
         endif
         
         call clft_bot_node_groups(i,cnt_addNode,grp1,grp2,grp3)
               
         cnt_N = 1
         
         do j = 1,nodesRingA
            
            if (j==1) then
               if (i.eq.1) area_node(count_area,cnt_N) = clft_top_endNode(2)
               if (i.gt.1) area_node(count_area,cnt_N) = cnt_rglNode-1
               cnt_N=cnt_N+1
               
            elseif (j==2) then
               
               do k=1,NAEC
                  area_node(count_area,cnt_N) = grp1(k)
                  cnt_N=cnt_N+1
               enddo
               
            elseif (j==3) then
               if (i.eq.1) area_node(count_area,cnt_N) = clft_top_endNode(3)
               if (i.gt.1) area_node(count_area,cnt_N) = cnt_rglNode
               cnt_N=cnt_N+1
               
            elseif (j==4) then  
               
               do k=1,NAEC
                  area_node(count_area,cnt_N) = grp2(k)
                  cnt_N=cnt_N+1
               enddo
               
            elseif (j==5) then
               if (i.eq.1) area_node(count_area,cnt_N) = cnt_rglNode+5
               if (i.gt.1) area_node(count_area,cnt_N) = cnt_rglNode+4
               cnt_N=cnt_N+1
               
            elseif (j==6) then
               
               do k=1,NAEC
                  area_node(count_area,cnt_N) = grp3(k)
                  cnt_N=cnt_N+1
               enddo
               
            elseif (j==nodesRingA) then
               if (i.eq.1) area_node(count_area,cnt_N) = cnt_rglNode+4
               if (i.gt.1) area_node(count_area,cnt_N) = cnt_rglNode+3
               cnt_N = cnt_N+1
            endif
            
         enddo
         
         !if (i==1) then
         !area_node(count_area,1) = clft_top_endNode(2)
         !area_node(count_area,2) = cnn_lft
         !area_node(count_area,3) = clft_top_endNode(3)
         !area_node(count_area,4) = (N_nodeS+2) + 2*(count_area-1) - 1
         !area_node(count_area,5) = 8 + 2*(count_area-1)
         !area_node(count_area,6) = (N_nodeS+2) + 2*(count_area-1)
         !area_node(count_area,7) = 8 + 2*(count_area-1) - 1
         !else
         !area_node(count_area,1) = 3 + 2*(count_area-1)
         !area_node(count_area,2) = (N_nodeS-2) + 2*(count_area-1)
         !area_node(count_area,3) = 3 + 2*(count_area-1) + 1
         !area_node(count_area,4) = (N_nodeS-2) + 2*(count_area-1) + 3
         !area_node(count_area,5) = 7 + 2*(count_area-1) + 1
         !area_node(count_area,6) = (N_nodeS-2) + 2*(count_area-1) + 4
         !area_node(count_area,7) = 7 + 2*(count_area-1)
         !endif
         
      enddo
      
    end subroutine clft_bot_area_to_node_NI
    
    subroutine crght_bot_area_to_node_NI
      implicit none
      integer :: i,j,k
      integer :: cnt_N
      integer :: grp1(1:NAEC),grp2(1:NAEC),grp3(1:NAEC)
      
      count_area = ncl + ncr + ncclt + nccrt + 1
      
      do i = 1,nccrb
         
         if (i==1) then
            count_area = count_area + 1
            area_node(count_area,0) = max_node_area-1
            cnt_rglNode = crght_top_endNode(3)
            cnt_addNode = cnn_rght-1+ 2*(NAEC)
            
         elseif (i==2) then
            count_area = count_area + 2
            area_node(count_area,0) = max_node_area-1
            cnt_rglNode = cnt_rglNode+4
            cnt_addNode = cnt_addNode+4*(NAEC)
            
         elseif (i.gt.2) then
            count_area = count_area + 2
            area_node(count_area,0) = max_node_area-1
            cnt_rglNode = cnt_rglNode+4
            cnt_addNode = cnt_addNode+4*(NAEC)
         endif
         
         call crght_bot_node_groups(i,cnt_addNode,grp1,grp2,grp3)
         
         cnt_N = 1
         
         do j=1,nodesRingA
            
            if (j==1) then
               if (i.eq.1) area_node(count_area,cnt_N) = crght_top_endNode(3)
               if (i.gt.1) area_node(count_area,cnt_N) = cnt_rglNode
               cnt_N=cnt_N+1
               
            elseif (j==2) then
               
               do k=1,NAEC
                  area_node(count_area,cnt_N) = grp1(k)
                  cnt_N=cnt_N+1
               enddo
               
            elseif (j==3) then
               if (i.eq.1) area_node(count_area,cnt_N) = crght_top_endNode(2)
               if (i.gt.1) area_node(count_area,cnt_N) = cnt_rglNode-1
               cnt_N=cnt_N+1
               
            elseif (j==4) then
               area_node(count_area,cnt_N) = cnt_rglNode+3
               cnt_N=cnt_N+1
               
            elseif (j==5) then
               
               do k=1,NAEC
                  area_node(count_area,cnt_N) = grp2(k)
                  cnt_N=cnt_N+1
               enddo
               
            elseif (j==6) then
               area_node(count_area,cnt_N) = cnt_rglNode+4
               cnt_N=cnt_N+1
               
            elseif (j==nodesRingA) then
               
               do k=1,NAEC
                  area_node(count_area,cnt_N) = grp3(k)
                  cnt_N=cnt_N+1
               enddo
               
            endif
            
         enddo
         
         !if (i==1) then
         !area_node(count_area,1) = crght_top_endNode(3)
         !area_node(count_area,2) = cnn_rght
         !area_node(count_area,3) = crght_top_endNode(2)
         !area_node(count_area,4) = 8 + 2*(count_area-1) - 1
         !area_node(count_area,5) = (N_nodeS+2) + 2*(count_area-1)
         !area_node(count_area,6) = 8 + 2*(count_area-1)
         !area_node(count_area,7) = (N_nodeS+2) + 2*(count_area-1) - 1  
         !else
         !area_node(count_area,1) = 4 + 2*(count_area-1)
         !area_node(count_area,2) = (N_nodeS-2) + 2*(count_area-1)
         !area_node(count_area,3) = 4 + 2*(count_area-1) - 1
         !area_node(count_area,4) = 8 + 2*(count_area-1) - 1
         !area_node(count_area,5) = (N_nodeS-2) + 2*(count_area-1) + 4
         !rea_node(count_area,6) = 8 + 2*(count_area-1)
         !area_node(count_area,7) = (N_nodeS-2) + 2*(count_area-1) + 3   
         !endif
         
      enddo
      
    end subroutine crght_bot_area_to_node_NI
    
    subroutine cntrlRgn_area_to_node_NI
      implicit none
      integer :: j,k
      integer :: grp1(1:NAEC_Ltrl),grp2(1:NAEC_Bsal),grp3(1:NAEC_Ltrl),grp4(1:NAEC_Apcl)
      integer :: cnt_N
      
      count_area = count_area + 1
      
      if (CellsMeet==0) then
         area_node(count_area,0) = max_node_area-1
         
         if (NI_incldAS_woPP == 0) nodesRingCR = 7
         if (NI_incldAS_woPP == 1) nodesRingCR = 8
         
      elseif (CellsMeet.gt.0) then
         area_node(count_area,0) = max_node_area-2
         
         if (NI_incldAS_woPP == 0) nodesRingCR = 6
         if (NI_incldAS_woPP == 1) nodesRingCR = 6 !same nodesRingCR 
         
      endif
      
      write(*,*) lft_endNode(1:3),rght_endNode(1:3),cnt_addNode,"lft_endNode,rght_endNode,cnt_addNode"
      
      cnt_N=1
      
      call cntrlRgn_node_groups(cnt_addNode,grp1,grp2,grp3,grp4)
      
      do j = 1,nodesRingCR
         
         if (j==1) then
            area_node(count_area,cnt_N) = lft_endNode(1) + 2*CellsMeet
            cnt_N = cnt_N+1
            
         elseif (j==2) then
            
            do k=1,NAEC
               area_node(count_area,cnt_N) = grp1(k)
               cnt_N = cnt_N+1
            enddo
            
         elseif (j==3) then
            area_node(count_area,cnt_N) = lft_endNode(2) + 2*CellsMeet
            cnt_N = cnt_N+1
            
         elseif (j==4) then
            
            do k=1,NAEC
               area_node(count_area,cnt_N) = grp2(k)
               cnt_N = cnt_N+1
            enddo
            
         elseif (j==5) then
            area_node(count_area,cnt_N) = rght_endNode(2) + 2*CellsMeet
            cnt_N = cnt_N+1
            
         elseif (j==6) then
            
            do k=1,NAEC
               area_node(count_area,cnt_N) = grp3(k)
               cnt_N = cnt_N+1
            enddo
            
         elseif (j==7) then
            area_node(count_area,cnt_N) = rght_endNode(1) + 2*CellsMeet
            cnt_N = cnt_N+1
            
         elseif (j==8) then
            
            do k=1,NAEC
               area_node(count_area,cnt_N) = grp4(k)
               cnt_N = cnt_N+1
            enddo
            
         endif
         
      enddo
      
    end subroutine cntrlRgn_area_to_node_NI
    
    subroutine endRgn_area_to_node_NI
      implicit none
      integer :: j,k
      integer :: grp1(1:NAEC),grp2(1:NAEC),grp3(1:NAEC)
      integer :: cnt_N
      
      count_area = count_area + 1
      area_node(count_area,0) = max_node_area-2
      cnt_addNode = cnt_addNode + 4*(NAEC)
      
      cnt_N=1

      call endRgn_node_groups(cnt_addNode,grp1,grp2,grp3)

      do j = 1,nodesRingER

         if (j==1) then
            area_node(count_area,cnt_N) = clft_bot_endNode(1)
            cnt_N = cnt_N+1
            
         elseif (j==2) then
            
            do k=1,NAEC
               area_node(count_area,cnt_N) = grp1(k)
               cnt_N=cnt_N+1
            enddo
            
         elseif (j==3) then
            area_node(count_area,cnt_N) = clft_bot_endNode(2)
            cnt_N = cnt_N+1
            
         elseif (j==4) then
            
            do k=1,NAEC
               area_node(count_area,cnt_N) = grp2(k)
               cnt_N=cnt_N+1
            enddo
            
         elseif (j==5) then
            area_node(count_area,cnt_N) = crght_bot_endNode(2)
            cnt_N = cnt_N+1
            
         elseif (j==nodesRingER) then
            
            do k=1,NAEC
               area_node(count_area,cnt_N) = grp3(k)
               cnt_N=cnt_N+1
            enddo
            
         endif
         
      enddo
      
      !area_node(count_area,1) = clft_bot_endNode(1)
      !area_node(count_area,2) = (N_nodeS-2) + 2*(count_area-1)
      !area_node(count_area,3) = clft_bot_endNode(2)
      !area_node(count_area,4) = (N_nodeS-2) + 2*(count_area-1) + 3
      !area_node(count_area,5) = crght_bot_endNode(2)
      !area_node(count_area,6) = (N_nodeS-2) + 2*(count_area-1) + 2
      
    end subroutine endRgn_area_to_node_NI
    
    
    subroutine other_area_to_node_NI
      implicit none
      
      area_node((count_area+1):N_cell,0:max_node_areaS) = area_nodeS((count_area+1):N_cell,0:max_node_areaS)

    end subroutine other_area_to_node_NI
    
    subroutine get_a1_of_AP_clft_bot !AP=Arithmatic Progression [an=a1+(n-1)d]
      implicit none
      
      continue
      
    end subroutine get_a1_of_AP_clft_bot
    
  end subroutine get_area_to_node_transfrm_NI_model
  
  subroutine get_nnfbcl(cell_order,nnfbcl)
    implicit none
    integer, intent(in)  :: cell_order
    integer, intent(out) :: nnfbcl(1:max_node_area,1:2)
    
    if (NI_incldAS_woPP == 0) call get_nnfbcl_methd1(cell_order,nnfbcl)
    if (NI_incldAS_woPP == 1) call get_nnfbcl_methd2(cell_order,nnfbcl)
     
  end subroutine get_nnfbcl
  
  subroutine get_nnfbcl_methd1(cell_order,nnfbcl)
    implicit none
    integer, intent(in)  :: cell_order
    integer, intent(out) :: nnfbcl(1:max_node_area,1:2)
    
    integer :: i,j
    integer :: cnt_N
    
    nnfbcl = -1
    
    if (cell_order==1) then
       
       !ANTI-CLOCKWISE
       
       cnt_N=1
       
       do i = 1,nodesRingF
          
          if (i==1) then
             nnfbcl(cnt_N,1) = 1
             nnfbcl(cnt_N,2) = 1
             cnt_N = cnt_N+1
             
          elseif (i==2) then
             nnfbcl(cnt_N,1) = 1
             nnfbcl(cnt_N,2) = 2
             cnt_N = cnt_N+1
             
          elseif (i==3 .or. i==5) then
             
             do j=1,NAEC
                nnfbcl(cnt_N,1) = 2
                
                if (i==3) nnfbcl(cnt_N,2) = N_nodeS+j
                if (i==5) nnfbcl(cnt_N,2) = N_nodeS+(NAEC)+(NAEC-(j-1))
                
                cnt_N=cnt_N+1
             enddo
             
          elseif (i==4) then
             nnfbcl(cnt_N,1) = 1
             nnfbcl(cnt_N,2) = 4
             cnt_N = cnt_N+1
             
          elseif (i==6) then
             nnfbcl(cnt_N,1) = 1
             nnfbcl(cnt_N,2) = 3
             cnt_N = cnt_N+1
          endif
          
       enddo
       
       !write(*,*) nnfbcl(1:max_node_area,1),"nnfbcl typs"
       !write(*,*) nnfbcl(1:max_node_area,2),"nnfbcl vals"
       
    elseif (cell_order.gt.1 .and. cell_order.le.ncl) then
       
       cnt_N=1
       
       do i = 1,nodesRingA
          
          if (i==1) then
             nnfbcl(cnt_N,1) = 1
             nnfbcl(cnt_N,2) = 1
             cnt_N = cnt_N+1
             
          elseif (i==3) then
             nnfbcl(cnt_N,1) = 1
             nnfbcl(cnt_N,2) = 2
             cnt_N = cnt_N+1
             
          elseif (i==5) then
             nnfbcl(cnt_N,1) = 1
             nnfbcl(cnt_N,2) = 4
             cnt_N = cnt_N+1
             
          elseif (i==7) then
             nnfbcl(cnt_N,1) = 1
             nnfbcl(cnt_N,2) = 3
             cnt_N = cnt_N+1
             
          elseif (i==2 .or. i==4 .or. i==6) then
             
             do j=1,NAEC
                nnfbcl(cnt_N,1) = 2
                
                if (i==2) nnfbcl(cnt_N,2) = N_nodeS-(NAEC-j)
                if (i==4) nnfbcl(cnt_N,2) = N_nodeS+(j)
                if (i==6) nnfbcl(cnt_N,2) = N_nodeS+(NAEC)+(NAEC-(j-1))
                
                cnt_N=cnt_N+1
             enddo
             
          endif
       enddo
       
       !write(*,*) nnfbcl(1:max_node_area,1),"nnfbcl typs"
       !write(*,*) nnfbcl(1:max_node_area,2),"nnfbcl vals"
       
    endif
    
  end subroutine get_nnfbcl_methd1
  
   
  subroutine get_nnfbcl_methd2(cell_order,nnfbcl)
    implicit none
    integer, intent(in)  :: cell_order
    integer, intent(out) :: nnfbcl(1:max_node_area,1:2)
    
    integer :: i,j
    integer :: cnt_N
    
    nnfbcl = -1
    
    if (cell_order==1) then
       
       !ANTI-CLOCKWISE
       
       cnt_N=1
       
       do i = 1,nodesRingF_withAS
          
          if (i==1) then
             nnfbcl(cnt_N,1) = 1
             nnfbcl(cnt_N,2) = 1
             cnt_N = cnt_N+1
             
          elseif (i==2) then
             nnfbcl(cnt_N,1) = 1
             nnfbcl(cnt_N,2) = 2
             cnt_N = cnt_N+1
             
          elseif (i==3 .or. i==5 .or. i==7) then
             
             do j=1,NAEC
                nnfbcl(cnt_N,1) = 2
                
                if (i==3) nnfbcl(cnt_N,2) = N_nodeS + (NAEC) + j
                if (i==5) nnfbcl(cnt_N,2) = N_nodeS + (NAEC) + (NAEC) + (NAEC-(j-1))
                if (i==7) nnfbcl(cnt_N,2) = N_nodeS + (NAEC-(j-1)) 
                
                cnt_N=cnt_N+1
             enddo
             
          elseif (i==4) then
             nnfbcl(cnt_N,1) = 1
             nnfbcl(cnt_N,2) = 4
             cnt_N = cnt_N+1
             
          elseif (i==6) then
             nnfbcl(cnt_N,1) = 1
             nnfbcl(cnt_N,2) = 3
             cnt_N = cnt_N+1
          endif
          
       enddo
       
       !write(*,*) nnfbcl(1:max_node_area,1),"nnfbcl typs"
       !write(*,*) nnfbcl(1:max_node_area,2),"nnfbcl vals"
       
    elseif (cell_order.gt.1 .and. cell_order.le.ncl) then
       
       cnt_N=1
       
       do i = 1,nodesRingA_withAS
          
          if (i==1) then
             nnfbcl(cnt_N,1) = 1
             nnfbcl(cnt_N,2) = 1
             cnt_N = cnt_N+1
             
          elseif (i==3) then
             nnfbcl(cnt_N,1) = 1
             nnfbcl(cnt_N,2) = 2
             cnt_N = cnt_N+1
             
          elseif (i==5) then
             nnfbcl(cnt_N,1) = 1
             nnfbcl(cnt_N,2) = 4
             cnt_N = cnt_N+1
             
          elseif (i==7) then
             nnfbcl(cnt_N,1) = 1
             nnfbcl(cnt_N,2) = 3
             cnt_N = cnt_N+1
             
          elseif (i==2 .or. i==4 .or. i==6 .or. i==8) then
             
             do j=1,NAEC
                nnfbcl(cnt_N,1) = 2
                
                if (i==2) nnfbcl(cnt_N,2) = N_nodeS - (NAEC-j)
                if (i==4) nnfbcl(cnt_N,2) = N_nodeS + (NAEC)   + (j)
                if (i==6) nnfbcl(cnt_N,2) = N_nodeS + (NAEC)   + (NAEC) + (NAEC-(j-1))
                if (i==8) nnfbcl(cnt_N,2) = N_nodeS + (NAEC-(j-1))
                
                cnt_N=cnt_N+1
             enddo
             
          endif
       enddo
       
       !write(*,*) nnfbcl(1:max_node_area,1),"nnfbcl typs"
       !write(*,*) nnfbcl(1:max_node_area,2),"nnfbcl vals"
       
    endif
    
  end subroutine get_nnfbcl_methd2
  
  subroutine get_nnfbcr(cell_order,cnt_rglNode,cnt_addNode,nnfbcr)
    implicit none
    integer, intent(in)  :: cell_order
    integer, intent(in)  :: cnt_rglNode,cnt_addNode
    integer, intent(out) :: nnfbcr(1:max_node_area,1:2)
    
    if (NI_incldAS_woPP == 0) call get_nnfbcr_methd1(cell_order,cnt_rglNode,cnt_addNode,nnfbcr)
    if (NI_incldAS_woPP == 1) call get_nnfbcr_methd2(cell_order,cnt_rglNode,cnt_addNode,nnfbcr)
    
  end subroutine get_nnfbcr
  
  subroutine get_nnfbcr_methd1(cell_order,cnt_rglNode,cnt_addNode,nnfbcr)
    implicit none
    integer, intent(in)  :: cell_order
    integer, intent(in)  :: cnt_rglNode,cnt_addNode
    integer, intent(out) :: nnfbcr(1:max_node_area,1:2)
    
    integer :: i,j,cnt_N
    
    nnfbcr = -1
    
    if (cell_order==1) then
       
       !ANTI-CLOCKWISE
       
       cnt_N=1
       
       do i = 1,nodesRingF

          if (i==1) then
             nnfbcr(cnt_N,1) = 1
             nnfbcr(cnt_N,2) = cnt_rglNode+1
             cnt_N = cnt_N+1
             
          elseif (i==2) then
             nnfbcr(cnt_N,1) = 1
             nnfbcr(cnt_N,2) = cnt_rglNode+3
             cnt_N = cnt_N+1
             
          elseif (i==3 .or. i==5) then
             
             do j=1,NAEC
                nnfbcr(cnt_N,1) = 2
                
                if (i==3) nnfbcr(cnt_N,2) = cnt_addNode+(NAEC)+j
                if (i==5) nnfbcr(cnt_N,2) = cnt_addNode+(NAEC-(j-1))
                
                cnt_N=cnt_N+1
             enddo
             
          elseif (i==4) then
             nnfbcr(cnt_N,1) = 1
             nnfbcr(cnt_N,2) = cnt_rglNode+4
             cnt_N = cnt_N+1
             
          elseif (i==6) then
             nnfbcr(cnt_N,1) = 1
             nnfbcr(cnt_N,2) = cnt_rglNode+2
             cnt_N = cnt_N+1
          endif
          
       enddo
       
       !write(*,*) nnfbcr(1:max_node_area,1),"nnfbcr typs"
       !write(*,*) nnfbcr(1:max_node_area,2),"nnfbcr vals"
       
    elseif (cell_order.gt.1.and.cell_order.le.ncr) then
       
       cnt_N = 1
       
       do i = 1,nodesRingA
          
          if (i==1) then
             nnfbcr(cnt_N,1) = 1
             nnfbcr(cnt_N,2) = cnt_rglNode+1
             cnt_N = cnt_N+1
             
          elseif (i==2) then
             nnfbcr(cnt_N,1) = 1
             nnfbcr(cnt_N,2) = cnt_rglNode+3
             cnt_N = cnt_N+1
             
          elseif (i==4) then
             nnfbcr(cnt_N,1) = 1
             nnfbcr(cnt_N,2) = cnt_rglNode+4
             cnt_N = cnt_N+1
             
          elseif (i==6) then
             nnfbcr(cnt_N,1) = 1
             nnfbcr(cnt_N,2) = cnt_rglNode+2
             cnt_N = cnt_N+1
             
          elseif (i==3 .or. i==5 .or. i==7) then
             
             do j=1,NAEC
                nnfbcr(cnt_N,1) = 2
                
                if (i==3) nnfbcr(cnt_N,2) = cnt_addNode+(NAEC)+j
                if (i==5) nnfbcr(cnt_N,2) = cnt_addNode+(NAEC-(j-1))
                if (i==7) nnfbcr(cnt_N,2) = cnt_addNode-(j-1)
                
                cnt_N = cnt_N+1
             enddo
             
          endif
       enddo
       
       !write(*,*) nnfbcr(1:max_node_area,1),"nnfbcr typs"
       !write(*,*) nnfbcr(1:max_node_area,2),"nnfbcr vals"
       
    endif
    
  end subroutine get_nnfbcr_methd1
  
  subroutine get_nnfbcr_methd2(cell_order,cnt_rglNode,cnt_addNode,nnfbcr)
    implicit none
    integer, intent(in)  :: cell_order
    integer, intent(in)  :: cnt_rglNode,cnt_addNode
    integer, intent(out) :: nnfbcr(1:max_node_area,1:2)
    
    integer :: i,j,cnt_N
    
    nnfbcr = -1
    
    if (cell_order==1) then
       
       !ANTI-CLOCKWISE
       
       cnt_N=1
       
       do i = 1,nodesRingF_withAS

          if (i==1) then
             nnfbcr(cnt_N,1) = 1
             nnfbcr(cnt_N,2) = cnt_rglNode+1
             cnt_N = cnt_N+1
             
          elseif (i==3) then
             nnfbcr(cnt_N,1) = 1
             nnfbcr(cnt_N,2) = cnt_rglNode+3
             cnt_N = cnt_N+1
             
          elseif (i==2 .or. i==4 .or. i==6) then
             
             do j=1,NAEC
                nnfbcr(cnt_N,1) = 2
                
                if (i==2) nnfbcr(cnt_N,2) = cnt_addNode +  j
                if (i==4) nnfbcr(cnt_N,2) = cnt_addNode + (NAEC) + (NAEC) + j
                if (i==6) nnfbcr(cnt_N,2) = cnt_addNode + (NAEC) + (NAEC-(j-1))
                
                cnt_N=cnt_N+1
             enddo
             
          elseif (i==5) then
             nnfbcr(cnt_N,1) = 1
             nnfbcr(cnt_N,2) = cnt_rglNode+4
             cnt_N = cnt_N+1
             
          elseif (i==7) then
             nnfbcr(cnt_N,1) = 1
             nnfbcr(cnt_N,2) = cnt_rglNode+2
             cnt_N = cnt_N+1
          endif
          
       enddo
       
       !write(*,*) nnfbcr(1:max_node_area,1),"nnfbcr typs"
       !write(*,*) nnfbcr(1:max_node_area,2),"nnfbcr vals"
       
    elseif (cell_order.gt.1.and.cell_order.le.ncr) then
       
       cnt_N = 1
       
       do i = 1,nodesRingA_withAS
          
          if (i==1) then
             nnfbcr(cnt_N,1) = 1
             nnfbcr(cnt_N,2) = cnt_rglNode+1
             cnt_N = cnt_N+1
             
          elseif (i==3) then
             nnfbcr(cnt_N,1) = 1
             nnfbcr(cnt_N,2) = cnt_rglNode+3
             cnt_N = cnt_N+1
             
          elseif (i==5) then
             nnfbcr(cnt_N,1) = 1
             nnfbcr(cnt_N,2) = cnt_rglNode+4
             cnt_N = cnt_N+1
             
          elseif (i==7) then
             nnfbcr(cnt_N,1) = 1
             nnfbcr(cnt_N,2) = cnt_rglNode+2
             cnt_N = cnt_N+1
             
          elseif (i==2 .or. i==4 .or. i==6 .or. i==8) then
             
             do j=1,NAEC
                nnfbcr(cnt_N,1) = 2

                if (i==2) nnfbcr(cnt_N,2) = cnt_addNode +  j
                if (i==4) nnfbcr(cnt_N,2) = cnt_addNode + (NAEC) + (NAEC) + j
                if (i==6) nnfbcr(cnt_N,2) = cnt_addNode + (NAEC) + (NAEC-(j-1))
                if (i==8) nnfbcr(cnt_N,2) = cnt_addNode-(j-1)
                
                cnt_N = cnt_N+1
             enddo
             
          endif
       enddo
       
       !write(*,*) nnfbcr(1:max_node_area,1),"nnfbcr typs"
       !write(*,*) nnfbcr(1:max_node_area,2),"nnfbcr vals"
       
    endif
    
    
  end subroutine get_nnfbcr_methd2
  
  subroutine clft_top_node_groups(cnt_addNode,grp1,grp2,grp3)
    implicit none
    integer, intent(in)  :: cnt_addNode
    integer, intent(out) :: grp1(1:NAEC),grp2(1:NAEC),grp3(1:NAEC)
    
    integer :: i,j
    
    do i = 1,3 !clft has three curve where nodes to be addded
       
       do j = 1,NAEC
          
          if (i==1) then
             grp1(j) = cnt_addNode + (NAEC) + j
          elseif (i==2) then
             grp2(j) = cnt_addNode + (NAEC-(j-1))
          elseif (i==3) then
             grp3(j) = lft_endNode(3) - (NAEC) + j
          endif
          
       enddo
       
    enddo
    
  end subroutine clft_top_node_groups
  
  
  subroutine crght_top_node_groups(cnt_addNode,grp1,grp2,grp3)
    implicit none
    integer, intent(in)  :: cnt_addNode
    integer, intent(out) :: grp1(1:NAEC),grp2(1:NAEC),grp3(1:NAEC)

    integer :: i,j
    
    do i = 1,3 !crght has three curve where nodes to be addded
       
       do j = 1,NAEC
          
          if (i==1) then
             grp1(j) = cnt_addNode + j 
          elseif (i==2) then
             grp2(j) = cnt_addNode + (NAEC) + (NAEC-(j-1))
          elseif (i==3) then
             grp3(j) = rght_endNode(3) - (j-1)
          endif
          
       enddo
       
    enddo
    
    
  end subroutine crght_top_node_groups
  
  subroutine clft_bot_node_groups(cell_order,cnt_addNode,grp1,grp2,grp3)
    implicit none
    integer, intent(in)  :: cell_order,cnt_addNode
    integer, intent(out) :: grp1(1:NAEC),grp2(1:NAEC),grp3(1:NAEC)
    
    integer :: i,j
    
    do i = 1,3
       
       do j = 1,NAEC
          
          if (i==1) then
             if (cell_order==1) then
                grp1(j) = cnn_lft + (j-1)
             elseif (cell_order.gt.1) then
                grp1(j) = cnt_addNode -(NAEC) + (j)
             endif
          elseif (i==2) then
             grp2(j) = cnt_addNode + (2*NAEC) + j
          elseif (i==3) then
             grp3(j) = cnt_addNode + (3*NAEC) + (NAEC-(j-1))
          endif
          
       enddo
       
    enddo
    
  end subroutine clft_bot_node_groups
  
  subroutine crght_bot_node_groups(cell_order,cnt_addNode,grp1,grp2,grp3)
    implicit none
    integer, intent(in)  :: cell_order,cnt_addNode
    integer, intent(out) :: grp1(1:NAEC),grp2(1:NAEC),grp3(1:NAEC)
    
    integer :: i,j
    
    do i = 1,3
       
       do j = 1,NAEC
          
          if (i==1) then
             if (cell_order==1) then
                grp1(j) = cnn_rght + NAEC - j
             elseif (cell_order.gt.1) then
                grp1(j) = cnt_addNode - (j-1)
             endif
          elseif (i==2) then
             grp2(j) = cnt_addNode + (3*NAEC) + j
          elseif (i==3) then
             grp3(j) = cnt_addNode + (2*NAEC) + (NAEC-(j-1))
          endif
          
       enddo
       
    enddo
    
  end subroutine crght_bot_node_groups
  
  subroutine cntrlRgn_node_groups(cnt_addNode,grp1,grp2,grp3,grp4)
    implicit none
    integer, intent(in)  :: cnt_addNode
    integer, intent(out) :: grp1(1:NAEC_Ltrl),grp2(1:NAEC_Bsal),grp3(1:NAEC_Ltrl),grp4(1:NAEC_Apcl)
    
    integer :: i,imax,j,jmax
    
    
    if (NI_incldAS_woPP == 0) imax = 3
    if (NI_incldAS_woPP == 1) imax = 4
    
    grp1=0 ; grp2=0 ; grp3=0 ; grp4=0
    
    do i = 1,imax
       
       do j = 1,NAEC

          if (i==1) then
             grp1(j) = lft_endNode(3) - (NAEC) + j
          elseif (i==2) then
             grp2(j) = cnt_addNode + (NAEC_Apcl) + j
          elseif (i==3) then
             grp3(j) = cnt_addNode - (j-1)
          elseif (i==4) then
             grp4(j) = cnt_addNode + (NAEC_Ltrl) - (j-1)
          endif
          
       enddo
       
    enddo
    
  end subroutine cntrlRgn_node_groups
  
  subroutine endRgn_node_groups(cnt_addNode,grp1,grp2,grp3)
    implicit none
    integer, intent(in)  :: cnt_addNode
    integer, intent(out) :: grp1(1:NAEC),grp2(1:NAEC),grp3(1:NAEC)
    
    integer :: i,j
    
    do i = 1,3
       
       do j = 1,NAEC
          
          if (i==1) then
             grp1(j) = cnt_addNode - (3*NAEC) + j
          elseif (i==2) then
             grp2(j) = cnt_addNode + j
          elseif(i==3) then
             grp3(j) = cnt_addNode - (j-1)
          endif
          
       enddo
       
    enddo
    
  end subroutine endRgn_node_groups


  subroutine adjstmnt_area_to_node_for_S1T1(area_indctr,area_no)
    implicit none
    integer, intent(in)  :: area_indctr,area_no
    integer, allocatable :: AreaNodeStr(:)
    
    integer :: prvSpce,curSpce
    
    if (NI_incldAS_woPP == 1) then
       write(*,*) "sbrtn: adjstmnt_area_to_node_for_S1T1 -->area_to_other_trnsfrm.f08 needs to be updated";stop
    endif
    
    
    if (area_indctr==1) then ! 1 is for lft
       prvSpce = area_node(area_no,0)
       curSpce = prvSpce+1
       allocate(AreaNodeStr(1:prvSpce))
       
       AreaNodeStr(1:prvSpce) = area_node(area_no,1:prvSpce)
       area_node(area_no,0)   = curSpce
       
       area_node(area_no,1:prvSpce) = AreaNodeStr(1:prvSpce)
       area_node(area_no,curSpce)   = N_nodeS
       deallocate(AreaNodeStr)
       
    elseif (area_indctr==2) then ! 2 is for rght
       prvSpce = area_node(area_no,0)
       curSpce = prvSpce+1
       allocate(AreaNodeStr(1:prvSpce))
       
       AreaNodeStr(1:prvSpce) = area_node(area_no,1:prvSpce) 
       area_node(area_no,0) = curSpce
       area_node(area_no,1) = AreaNodeStr(1)
       area_node(area_no,2) = N_nodeS
       
       area_node(area_no,3:curSpce) = AreaNodeStr(2:prvSpce)
       
       deallocate(AreaNodeStr)
       
    endif
    
  end subroutine adjstmnt_area_to_node_for_S1T1
  
  subroutine adjst_area_to_node_trnsfrm_NI
    implicit none
    integer :: i
    integer :: lim1,lim2,lim3,lim4,lim5

    write(*,*) Hlf_Ncell,"Hlf_Ncell"
    if (Hlf_Ncell .ne. ((N_cell-1)/2)) then
       write(*,*) "stopped due to Hlf_Ncell not equal to what it is"
       stop
    endif
    
    write(*,*) N_nodeS,"N_nodeS" ! 50 (includes two added node)
    
    lim1 = Hlf_Ncell-1 ! = 10
    lim2 = Hlf_Ncell
    lim3 = lim2+(Hlf_Ncell-1)
    lim4 = N_cell-1
    lim5 = N_cell
    
    if (SystemTyp == 1) then
       
       do i = 1,N_cell
          
          if (i.le.lim1) then
             area_node(i,0) = area_node(i,0) - 1
             
          elseif ((i.gt.lim1) .and. (i.le.lim2)) then
             area_node(i,0)               = area_node(i,0) - 1
             area_node(i,max_node_area-1) = area_node(i,max_node_area)-1
             area_node(i,max_node_area)   = -1
             
          elseif ((i.gt.lim2) .and. (i.le.lim3)) then
             area_node(i,0) = area_node(i,0) - 1
             
          elseif ((i.gt.lim3) .and. (i.le.lim4)) then
             area_node(i,0) = area_node(i,0) - 1
             
          elseif ((i.gt.lim4) .and. (i.le.lim5)) then
             area_node(i,0)               = max_node_area
             area_node(i,max_node_area-1) = N_nodeS
             area_node(i,max_node_area)   = N_nodeS-1
          endif
          
       enddo
       
    endif
    
  end subroutine adjst_area_to_node_trnsfrm_NI
  
  subroutine print_area_to_node_trnsfrm
    implicit none
    integer :: i     
    
    open(unit=112,file="area_node.dat")
    open(unit=114,file="area_nodeWOtxt.dat")
    
    write(112,*) "area-node","   ","first one area_node(i,0)"
    !write(*,*) N_cell,"area_to_node"
    
    do i = 1,N_cell
       write(112,*) "area_to_node for area =",i,"are",area_node(i,0:max_node_area)
       write(114,*) i,area_node(i,0:max_node_area)
    enddo
    
    close(112)
    close(114)
    
  end subroutine print_area_to_node_trnsfrm
  
  
end module area_to_other_transfrm












