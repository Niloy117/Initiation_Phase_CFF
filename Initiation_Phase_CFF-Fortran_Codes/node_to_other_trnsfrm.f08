module node_to_other_transfrm
  use system_parameters
  use transfrm_info
  implicit none
  
  integer :: i
  
contains
  
  subroutine get_node_to_spr_transfrm
    implicit none
    integer :: i,lp_cnt
    integer :: count_vs
    
    
    count_nodes = 0
    count_vs    = 0
    
    lp_cnt = wntbb(0)

    if (modelID == 1) then
       
       do i = 1,max_blck
          
          if (wntbb(i) == 1) then
             call lft_node_to_spr
          elseif (wntbb(i) == 2) then
             call rght_node_to_spr
          elseif (wntbb(i) == 3) then
             call clft_top_node_to_spr
          elseif (wntbb(i) == 4) then
             call additnl_node_to_spr
          elseif (wntbb(i) == 5) then
             call crght_top_node_to_spr
          elseif (wntbb(i) == 6) then
             call clft_bot_node_to_spr
          elseif (wntbb(i) == 7) then
             call crght_bot_node_to_spr
          elseif (wntbb(i) == 8) then
             continue
          elseif (wntbb(i) == 9) then
             continue
          else
             continue
          endif
       
       enddo
       
       call trnsfrmtn_adjstmnt_node_spr_for_doubleNode
       
    elseif (modelID == 2) then
       call NI_models_trmnlNode_to_spr
       call Inserted_node_to_spr
       !trnsfrmtn_adjstmnt_node_spr_for_doubleNode not needed
       !in modelID=2
    endif
    
    call print_node_to_spr_transfrm
    
    !if (SystemTyp==1 .and. modelID==2) then
       !write(*,*) "Stoppeded here" 
       !stop
    !endif
    
  contains
    
    subroutine lft_node_to_spr
      implicit none
      integer :: i,j
      
      !write(*,*) nvsl,"nvsl"
      
      do i = 1,nvsl
         count_vs = count_vs + 1
         !write(*,*) count_vs,"count_vs"
         
         if (i==1) then
            do j = (2*count_vs-1),(2*count_vs)
               node_spr(j,0) = 1
               node_spr(j,1) = j
            enddo
            
         elseif ((i.ne.1) .and. (i.ne.nvsl)) then
            
            do j = (2*count_vs-1),(2*count_vs)
               node_spr(j,0) = 3
               
               if (mod(j,2) == 1) then
                  
                  node_spr(j,1) = (-2) + (count_vs-1)*3 ! nth term = a + (n-1)*d
                  node_spr(j,2) = (0)  + (count_vs-1)*3
                  node_spr(j,3) = (1)  + (count_vs-1)*3
                  
               elseif (mod(j,2) == 0) then
                  
                  node_spr(j,1) = (-1) + (count_vs-1)*3
                  node_spr(j,2) = (0) + (count_vs-1)*3
                  node_spr(j,3) = (2) + (count_vs-1)*3
                  
               endif
               
            enddo
            
         elseif (i==nvsl) then
            
            do j = (2*count_vs-1),(2*count_vs)

               if (stageNo==1 .and. stageType==1) then

                  if (cntrl_cell_nmbr .ne. 0) then
                     node_spr(j,0) = 3
                  else
                     node_spr(j,0) = 2
                  endif
                  
                  if (mod(j,2) == 1) then
                     
                     node_spr(j,1) = -2+(count_vs-1)*3  !nth term = a + (n-1)*d 
                     node_spr(j,2) = 0 +(count_vs-1)*3
                     node_spr(j,3) = cntrlCell_Spring(1)
                     !write(*,*) node_spr(j,1:3),"node_spr_lft1"
                     
                  elseif (mod(j,2) == 0) then               
                     node_spr(j,1) = -1+(count_vs-1)*3 ! nth term = a + (n-1)*d
                     node_spr(j,2) = 0 +(count_vs-1)*3
                     node_spr(j,3) = cntrlCell_Spring(2)
                     !write(*,*) node_spr(j,1:3),"node_spr_lft2" 
                  endif
                  
               elseif (stageNo==4) then
                  
                  if (clft_top_cell_nmbr .ne. 0) then
                     node_spr(j,0) = 3
                  else
                     node_spr(j,0) = 2
                  endif
                  
                  if (mod(j,2) == 1) then
                     
                     node_spr(j,1) = -2+(count_vs-1)*3  ! nth term = a + (n-1)*d 
                     node_spr(j,2) = 0 +(count_vs-1)*3
                     node_spr(j,3) =  clft_top_endSpring(1)
                     !write(*,*) node_spr(j,1:3),"node_spr_lft1"
                     
                  elseif (mod(j,2) == 0) then               
                     node_spr(j,1) = -1+(count_vs-1)*3 ! nth term = a + (n-1)*d
                     node_spr(j,2) = 0 +(count_vs-1)*3
                     node_spr(j,3) = clft_top_endSpring(2)
                     !write(*,*) node_spr(j,1:3),"node_spr_lft2"
                     
                  endif
               endif
               
            enddo
         endif
      enddo
      
      count_nodes = 2*count_vs
      
    end subroutine lft_node_to_spr

    subroutine rght_node_to_spr
      implicit none
      integer :: i,j
      
      !write(*,*) nvsr,"nvsr"
      !write(*,*) count_vs,"count_vs"
      
      do i = 1,nvsr
         count_vs = count_vs + 1
         !write(*,*) count_vs,"count_vs inside loop"
         
         if (i == 1) then
            do j = (2*count_vs-1),(2*count_vs)
               node_spr(j,0) = 1
               
               if (mod(j,2)==1) then
                  node_spr(j,1) = (-2) + (count_vs-1)*3
               elseif (mod(j,2)==0) then
                  node_spr(j,1) = (-1) + (count_vs-1)*3
               endif
               
            enddo
            
         elseif ((i.ne.1) .and. (i.ne.nvsr)) then
            
            do j = (2*count_vs-1),(2*count_vs)
               node_spr(j,0) = 3
               
               if (mod(j,2) == 1) then
                  node_spr(j,1) = (-5) + (count_vs-1)*3 ! nth term = a + (n-1)*d
                  node_spr(j,2) = (-3) + (count_vs-1)*3
                  node_spr(j,3) = (-2) + (count_vs-1)*3
                  
               elseif (mod(j,2) == 0) then
                  node_spr(j,1) = (-4) + (count_vs-1)*3
                  node_spr(j,2) = (-3) + (count_vs-1)*3
                  node_spr(j,3) = (-1) + (count_vs-1)*3
               endif
            enddo
            
         elseif (i==nvsr) then
            
            do j = (2*count_vs-1),(2*count_vs)
               
               if (stageNo==1 .and. stageType==1) then
                  
                  if (cntrl_cell_nmbr .ne. 0) then
                     node_spr(j,0) = 3
                  else
                     node_spr(j,0) = 2
                  endif
                  
                  if (mod(j,2) == 1) then
                     
                     node_spr(j,1) = -5+(count_vs-1)*3  !nth term = a + (n-1)*d 
                     node_spr(j,2) = -3+(count_vs-1)*3
                     node_spr(j,3) = cntrlCell_Spring(1)
                     !write(*,*) node_spr(j,1:3),"node_spr_lft1"
                     
                  elseif (mod(j,2) == 0) then               
                     node_spr(j,1) = -4+(count_vs-1)*3 !nth term = a + (n-1)*d
                     node_spr(j,2) = -3+(count_vs-1)*3
                     node_spr(j,3) = cntrlCell_Spring(2)
                     !write(*,*) node_spr(j,1:3),"node_spr_lft2"
                  endif
                  
               elseif (stageNo==4) then
                  
                  if (crght_top_cell_nmbr .ne. 0) then
                     node_spr(j,0) = 3
                  else
                     node_spr(j,0) = 2
                  endif
                  
                  if (mod(j,2) == 1) then
                     node_spr(j,1) = -5+(count_vs-1)*3 ! nth term = a + (n-1)*d
                     node_spr(j,2) = -3+(count_vs-1)*3 
                     node_spr(j,3) =  crght_top_endSpring(1)
                     !write(*,*) node_spr(j,1:3),"node_spr_rght1"
                     
                  elseif (mod(j,2) == 0) then               
                     node_spr(j,1) = -4+(count_vs-1)*3 ! nth term = a + (n-1)*d
                     node_spr(j,2) = -3+(count_vs-1)*3 
                     node_spr(j,3) = crght_top_endSpring(2)
                     !write(*,*) node_spr(j,1:3),"node_spr_rght2"
                     
                  endif
                  
               endif
               
            enddo
         endif
         
      enddo

      count_nodes = 2*count_vs
      
    end subroutine rght_node_to_spr
    
    subroutine clft_top_node_to_spr
      implicit none
      integer :: i
      
      do i = 1,clft_top_node 
         count_nodes = count_nodes + 1
         
         if (i == 1) then
            node_spr(count_nodes,0) = 4  !!!NEEDS TO BE ADDRESSED LATER
            node_spr(count_nodes,1) = clft_top_endSpring(1)
            node_spr(count_nodes,2) = clft_top_endSpring(1)
            node_spr(count_nodes,3) = crght_top_endSpring(1)
            node_spr(count_nodes,4) = crght_top_endSpring(1)
            
         elseif (i == 2) then
           
            if(wntbb(4).ne.0) then
               node_spr(count_nodes,0) = 5
               node_spr(count_nodes,1) = clft_top_endSpring(1)
               node_spr(count_nodes,2) = clft_top_endSpring(1) + 1
               node_spr(count_nodes,3) = additnl_node_endSpring(1)
               node_spr(count_nodes,4) = crght_top_endSpring(1)
               node_spr(count_nodes,5) = crght_top_endSpring(1) + 1
               
            elseif(wntbb(4)==0 .AND. wntbb(6).ne.0) then
               node_spr(count_nodes,0) = 6
               node_spr(count_nodes,1) = clft_top_endSpring(1)
               node_spr(count_nodes,2) = clft_top_endSpring(1) + 1
               node_spr(count_nodes,3) = crght_top_endSpring(1)
               node_spr(count_nodes,4) = crght_top_endSpring(1) + 1
               node_spr(count_nodes,5) = first_spr_of_clft_bot 
               node_spr(count_nodes,6) = first_spr_of_crght_bot
              
            endif
            
         elseif (i == (clft_top_node)) then
            if(wntbb(4).ne.0) then
               node_spr(count_nodes,0) = 2
               node_spr(count_nodes,1) = clft_top_endSpring(1) + 1
               node_spr(count_nodes,2) = clft_top_endSpring(2)
            elseif (wntbb(4)==0 .AND. wntbb(6).ne.0) then
               node_spr(count_nodes,0) = 3
               node_spr(count_nodes,1) = clft_top_endSpring(1) + 1
               node_spr(count_nodes,2) = clft_top_endSpring(2)
               node_spr(count_nodes,3) = first_spr_of_clft_bot + 1
            endif
            
         endif
         
      enddo 
     
    end subroutine clft_top_node_to_spr

    
    subroutine additnl_node_to_spr
      implicit none
      !integer :: i,j

      count_nodes = count_nodes + 1
      
      node_spr(count_nodes,0) = 1
      node_spr(count_nodes,1) = additnl_node_endSpring(1)
      
    end subroutine additnl_node_to_spr

    subroutine crght_top_node_to_spr
      implicit none
      integer :: i!,j

      do i = 1,crght_top_node 
         count_nodes = count_nodes + 1
         
         if (i == 1) then
            node_spr(count_nodes,0) = 4  !!!NEEDS TO BE ADDRESSED LATER
            node_spr(count_nodes,1) = clft_top_endSpring(1)
            node_spr(count_nodes,2) = clft_top_endSpring(1)
            node_spr(count_nodes,3) = crght_top_endSpring(1)
            node_spr(count_nodes,4) = crght_top_endSpring(1)
            
         elseif (i == 2) then
            
            if(wntbb(4).ne.0) then
               node_spr(count_nodes,0) = 5
               node_spr(count_nodes,1) = clft_top_endSpring(1)
               node_spr(count_nodes,2) = clft_top_endSpring(1) + 1
               node_spr(count_nodes,3) = additnl_node_endSpring(1)
               node_spr(count_nodes,4) = crght_top_endSpring(1)
               node_spr(count_nodes,5) = crght_top_endSpring(1) + 1
               
            elseif(wntbb(4)==0 .AND. wntbb(6).ne.0) then
               node_spr(count_nodes,0) = 6
               node_spr(count_nodes,1) = clft_top_endSpring(1)
               node_spr(count_nodes,2) = clft_top_endSpring(1) + 1
               node_spr(count_nodes,3) = crght_top_endSpring(1)
               node_spr(count_nodes,4) = crght_top_endSpring(1) + 1
               node_spr(count_nodes,5) = first_spr_of_clft_bot 
               node_spr(count_nodes,6) = first_spr_of_crght_bot

            endif
            
         elseif (i == (crght_top_node)) then

            if(wntbb(4).ne.0) then
               node_spr(count_nodes,0) = 2
               node_spr(count_nodes,1) = crght_top_endSpring(1) + 1
               node_spr(count_nodes,2) = crght_top_endSpring(2)
            elseif (wntbb(4)==0 .AND. wntbb(6).ne.0) then
               node_spr(count_nodes,0) = 3  
               node_spr(count_nodes,1) = crght_top_endSpring(1) + 1
               node_spr(count_nodes,2) = crght_top_endSpring(2)
               node_spr(count_nodes,3) = first_spr_of_crght_bot + 1
            endif
            
         endif
         
      enddo
      
    end subroutine crght_top_node_to_spr
    
    
    
    subroutine clft_bot_node_to_spr
      implicit none
      integer :: i,j
      integer :: spr1,spr2,spr3
      
      spr1=0 ; spr2=0 ; spr3=0
      write(*,*) nvsclb,"nvsclb"
      
      count_vs = nvsl + nvsr + nvsclt + nvscrt 
      
      do i = 1,nvsclb
         if (i==1) then
            count_vs = count_vs + 1
         elseif (i.ne.1) then
            count_vs = count_vs + 2
         endif
         
         write(*,*) count_vs,"count_vs"
         
         if (i.ne.nvsclb) then
            
            do j = (2*count_vs+1),(2*count_vs+2) !its 2*cnt_vs-1,2*cnt_vs in lft
               node_spr(j,0) = 3
               
               if (mod(j,2) == 1) then
                  write(*,*) count_vs,"count_vs bfr get_first_three"
                  call get_first_three_sprs(count_vs,spr1,spr2,spr3)
                  write(*,*) spr1,spr2,spr3,"spr1-3"
                  
                  node_spr(j,1) = spr1
                  node_spr(j,2) = spr2
                  node_spr(j,3) = spr3
                  
                  write(*,*) node_spr(j,1:3),"node_spr(j,1:3)"
               elseif (mod(j,2) == 0) then
                  
                  node_spr(j,1) = (spr1+1)
                  node_spr(j,2) = spr2 
                  node_spr(j,3) = (spr3+1)
                  write(*,*) node_spr(j,1:3),"node_spr(j,1:3)"
               endif
               
            enddo
            
         elseif (i==nvsclb) then
            
            do j = (2*count_vs+1),(2*count_vs+2)
               
               if (mod(j,2) == 1) then
                  node_spr(j,0) = 2
                  
                  write(*,*) count_vs,"count_vs bfr get_first_three"
                  call get_first_three_sprs(count_vs,spr1,spr2,spr3)
                  node_spr(j,1) = spr1   ! nth term = a + (n-1)*d 
                  node_spr(j,2) = spr2 
                  !write(*,*) node_spr(j,1:3),"node_spr_clb1"
                  
               elseif (mod(j,2) == 0) then
                  
                  if (stageNo==4 .and. stageType==1) then
                     node_spr(j,0) = 3
                     node_spr(j,1) = (spr1+1)  ! nth term = a + (n-1)*d
                     node_spr(j,2) = spr2 
                     node_spr(j,3) = spr3
                     !write(*,*) node_spr(j,1:3),"node_spr_clb2"
                     
                  elseif (stageNo==4 .and. stageType==2) then
                     node_spr(j,0) = 2
                     node_spr(j,1) = (spr1+1)  ! nth term = a + (n-1)*d
                     node_spr(j,2) = spr2
                     
                  endif
                  
               endif
               
            enddo
         endif
      enddo
      
      count_nodes = 2*(count_vs) + 2*(nvsclb-1)
      
            
    end subroutine clft_bot_node_to_spr
    
    
    
    subroutine crght_bot_node_to_spr
      implicit none
      integer :: i,j
      integer :: spr1,spr2,spr3
      
      write(*,*) nvscrb,"nvscrb"
      
      count_vs = nvsl + nvsr + nvsclt + nvscrt + 1 
      
      do i = 1,nvscrb
         if (i==1) then
            count_vs = count_vs + 1
         elseif (i.ne.1) then
            count_vs = count_vs + 2
         endif
         
         write(*,*) count_vs,"count_vs"
         
         if (i.ne.nvscrb) then
            
            do j = (2*count_vs+1),(2*count_vs+2)
               node_spr(j,0) = 3
               
               if (mod(j,2) == 1) then
                  call get_first_three_sprs(count_vs,spr1,spr2,spr3)
                  
                  node_spr(j,1) = spr1
                  node_spr(j,2) = spr2  !not using cnt_vs, using i
                  node_spr(j,3) = spr3 
                  
               elseif (mod(j,2) == 0) then
              
                  node_spr(j,1) = (spr1+1) 
                  node_spr(j,2) = spr2 
                  node_spr(j,3) = (spr3+1) 
                  
               endif
               
            enddo
            
         elseif (i==nvscrb) then
            
            do j = (2*count_vs+1),(2*count_vs+2)
               node_spr(j,0) = 2
               
               if (mod(j,2) == 1) then
                  node_spr(j,0) = 2
                  call get_first_three_sprs(count_vs,spr1,spr2,spr3)
                  node_spr(j,1) = spr1  
                  node_spr(j,2) = spr2 
                  !write(*,*) node_spr(j,1:3),"node_spr_crb1"
                  
               elseif (mod(j,2) == 0) then
                  
                  if (stageNo==4 .and. stageType==1) then
                     node_spr(j,0) = 3
                     node_spr(j,1) = (spr1+1)
                     node_spr(j,2) = spr2
                     node_spr(j,3) = spr3-3
                     
                  elseif (stageNo==4 .and. stageType==2) then
                     node_spr(j,0) = 2
                     node_spr(j,1) = (spr1+1)
                     node_spr(j,2) = spr2
                     !write(*,*) node_spr(j,1:3),"node_spr_crb2" 
                  endif
                  
               endif
               
            enddo
         endif
      enddo

      count_nodes = 2*(count_vs) + 2*(nvscrb-1)
      
    end subroutine crght_bot_node_to_spr
    
    subroutine get_first_three_sprs(vs,frst,scnd,thrd)
      implicit none
      integer :: vs !vertical_side
      integer :: frst,scnd,thrd
      
      frst = (-5) + (vs-1)*3
      scnd = (-3) + (vs-1)*3
      thrd = (1)  + (vs-1)*3
      
    end subroutine get_first_three_sprs

    
    subroutine NI_models_trmnlNode_to_spr
      implicit none
      integer :: i,j,m
      integer :: jmax
      integer :: spr_Prv,cnt,trmnl
      
      integer :: Init_rghtEndNode1,Init_lftEndNode1
      integer :: Init_rghtEndNode2,Init_lftEndNode2
      integer :: lim1,lim2
      integer :: LR !lfr rght node indicator for S1T1
      logical :: lgcl_2nd,lgcl_3rd
      
      !I want the rght_endNode(2) for starting position, this is for getting that rght_endNode(2) at CellsMeet=0, as rght_endNode(2) is changing
      
      trmnl = 1
      
      write(*,*) rght_endNode,"rghtEndNode in fl:node_to_other_trnsfrm ; sb: NI_models_trmnlNode_to_spr"
      write(*,*) lft_endNode, "lftEndNode  in fl:node_to_other_trnsfrm ; sb: NI_models_trmnlNode_to_spr"
      
      if (stageNo==1 .and. stageType==1) then
         
         open(unit=35,file='AlteredNode.dat',position='append')
         
         Init_lftEndNode1  = lft_endNode(1)  + 2*CellsMeet
         Init_rghtEndNode1 = rght_endNode(1) + 2*CellsMeet
         Init_lftEndNode2  = lft_endNode(2)  + 2*CellsMeet
         Init_rghtEndNode2 = rght_endNode(2) + 2*CellsMeet
         
         lim1 = lft_endNode(2) ; lim2 = rght_endNode(2)
         
         write(35,fmt=*) Init_rghtEndNode1,Init_lftEndNode1,"Init1"!47/23
         write(35,fmt=*) Init_rghtEndNode2,Init_lftEndNode2,"Init2"!48/24
         write(35,fmt=*) lft_endNode,rght_endNode,"lft_endNode,rght_endNode"
         close(35)
         
      endif
      
      do i = 1,N_nodeS
         
         jmax          = node_sprS(i,0)
         node_spr(i,0) = node_sprS(i,0)
         
         if (stageNo==1 .and. stageType==1) then
            
            if ((i.le.lim1) .or. (i.gt.Init_lftEndNode2 .and. i.le.lim2)) then
               cnt = 1
               
               do j = 1,jmax
                  spr_Prv = node_sprS(i,j)
                  
                  if (trmnl_intrmdSpr(spr_Prv,0)==1) then
                     node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,1) !; write(*,*) i,j,node_spr(i,cnt),cnt,"i,j,node_spr,cnt"
                     cnt = cnt+1
                     
                  else
                     
                     m = trmnl_intrmdSpr(spr_Prv,0) !; write(*,*) m,i,j,"m,i,j"
                     
                     if (j==1) then
                        
                        if(i==1.or.i==2.or.i==(Init_lftEndNode2+1).or.i==(Init_lftEndNode2+2)) then !%%%%%%%CHANGEd aft adding node in APCL side 
                           node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,1)
                        else
                           node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,m)
                        endif
                        
                        !write(*,*) i,j,node_spr(i,cnt),cnt,"i,j,node_spr,cnt"
                        cnt = cnt+1
                        
                     elseif (j==2) then
                        
                        if (trmnl==1) then
                           node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,1)
                           trmnl=trmnl+1
                        elseif (trmnl==2) then
                           node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,m)
                           trmnl=1
                        endif
                        
                        !write(*,*) i,j,node_spr(i,cnt),cnt,"i,j,node_spr,cnt"
                        cnt = cnt+1
                        
                     elseif (j==3) then
                        
                        if ((i.ne.Init_rghtEndNode1).and.(i.ne.Init_rghtEndNode2)) then !%%%%%%%CHANGEd aft adding node in APCL side
                           node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,1) !;write(*,*) i,j,node_spr(i,cnt),cnt,"i,j,node_spr,cnt"
                           cnt = cnt+1
                           
                        elseif ((i==Init_rghtEndNode1).or.(i==Init_rghtEndNode2)) then !%%%%%%%CHANGEd aft adding node in APCL side
                           node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,m) !; write(*,*) i,j,node_spr(i,cnt),cnt,"i,j,node_spr,cnt"
                           cnt = cnt+1
                        endif
                        
                     endif
                     
                  endif
                  
               enddo
               
            else
               
               if (i.le.Init_lftEndNode2) LR=1
               if (i.gt.Init_lftEndNode2) LR=2
               
               lgcl_2nd=.False.
               lgcl_3rd=.False.
               
               if (LR==1) then
                  
                  if ((i-lim1).le.2) then
                     lgcl_2nd=.True.
                  elseif ((i-lim1).gt.2) then
                     lgcl_3rd=.True.
                  endif
                     
               elseif (LR==2) then
                  
                  if ((i-lim2).le.2) then
                     lgcl_2nd=.True.
                  elseif ((i-lim2).gt.2) then
                     lgcl_3rd=.True.
                  endif
                  
               endif
               
               if (lgcl_2nd.eqv..True.) then
                  cnt = 1
                  
                  do j=1,jmax
                     spr_Prv = node_sprS(i,j)
                     
                     if (trmnl_intrmdSpr(spr_Prv,0)==1) then
                        node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,1)
                        !write(*,*) i,j,node_spr(i,cnt),cnt,"i,j,node_spr,cnt"
                        cnt = cnt+1
                        
                     else
                        
                        m = trmnl_intrmdSpr(spr_Prv,0)
                        
                        if ((i-lim1)==1 .or. (i-lim2)==1) then
                           node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,1)
                           !write(*,*) i,j,node_spr(i,cnt),cnt,"i,j,node_spr,cnt"
                           cnt = cnt+1
                           
                        elseif ((i-lim1)==2 .or. (i-lim2)==2) then
                           
                           if (j==1 .or. j==2) then
                              node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,m)
                              !write(*,*) i,j,node_spr(i,cnt),cnt,"i,j,node_spr,cnt"
                              cnt = cnt+1
                           
                           elseif (j==3) then
                           
                              if (i.ne.Init_rghtEndNode2) then
                                 node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,1)
                                 !write(*,*) i,j,node_spr(i,cnt),cnt,"i,j,node_spr,cnt"
                                 cnt = cnt+1
                                 
                              elseif (i==Init_rghtEndNode2) then
                                 node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,m)
                                 !write(*,*) i,j,node_spr(i,cnt),cnt,"i,j,node_spr,cnt"
                                 cnt = cnt+1
                              endif
                           
                           endif
                           
                        endif
                        
                     endif
                  enddo
                  
               elseif (lgcl_3rd.eqv..True.) then
                  
                  cnt = 1
                  
                  !if (i.ge.17.or.i.ge.35) then
                   !  open(unit=122,file='EntryChk.dat',position='append')
                    ! write(122,*) i,jmax,lim1,lim2,"jmax in 17-35,lims"
                     !close(122)
                  !endif
                  
                  do j = 1,jmax
                     spr_Prv = node_sprS(i,j)
                     
                     if (trmnl_intrmdSpr(spr_Prv,0)==1) then
                        node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,1)
                        !write(*,*) i,j,node_spr(i,cnt),cnt,"i,j,node_spr,cnt"
                        cnt = cnt+1
                        
                     else
                        m = trmnl_intrmdSpr(spr_Prv,0)
                        
                        if (j==1.or.j==4) then
                           
                           if ((i.ne.Init_lftEndNode1).and.(i.ne.Init_rghtEndNode1)) then
                              
                              if (j==1) node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,m)
                              if (j==4) node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,1)
                              
                           elseif ((i==Init_lftEndNode1).or.(i==Init_rghtEndNode1)) then
                              node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,1)
                           endif
                           
                           !write(*,*) i,j,node_spr(i,cnt),cnt,"i,j,node_spr,cnt"
                           cnt = cnt+1
                           
                           
                        elseif (j==2.or.j==5) then
                           !if (i==17) then
                            !  open(unit=122,file='trmnlChk',position='append')
                             ! write(122,*) trmnl,"trmnl in 17"
                              !close(122)
                           !endif
                           
                           if (trmnl==1) then
                              if (mod(i,2)==1) node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,1)
                              if (mod(i,2)==0) node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,m)
                              
                              if ((i.ne.Init_lftEndNode1).and.(i.ne.Init_rghtEndNode1)) then
                                 if(j==5) trmnl=trmnl+1 !othrwise will be incresd twice
                              elseif ((i==Init_lftEndNode1).or.(i==Init_rghtEndNode1)) then
                                 if(j==2) trmnl=trmnl+1 !j wont reach upto 5 here
                              endif
                              
                           elseif (trmnl==2) then
                              node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,m)
                              trmnl=1
                           endif
                           
                           !write(*,*) i,j,node_spr(i,cnt),cnt,"i,j,node_spr,cnt"
                           cnt = cnt+1
                           
                        elseif (j==3.or.j==6) then
                           
                           if ((i.ne.Init_lftEndNode2).and.(i.ne.Init_rghtEndNode2)) then
                              node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,1)
                           elseif (i==Init_lftEndNode2) then
                              node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,1)
                           elseif (i==Init_rghtEndNode2) then
                              node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,m)
                           endif
                           
                           !write(*,*) i,j,node_spr(i,cnt),cnt,"i,j,node_spr,cnt"
                           cnt = cnt+1
                           
                        endif
                        
                     endif
                  
                  enddo
               
               endif
               
            endif
            
         elseif (stageNo==4) then
            
            if (i.le.rght_endNode(2)) then
               cnt = 1
               
               do j = 1,jmax
                  spr_Prv = node_sprS(i,j)
                  
                  if (trmnl_intrmdSpr(spr_Prv,0)==1) then
                     node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,1)
                     !write(*,*) i,j,node_spr(i,cnt),cnt,"i,j,node_spr,cnt"
                     cnt = cnt+1
                     
                  else
                     
                     m = trmnl_intrmdSpr(spr_Prv,0)
                     !write(*,*) m,i,j,"m,i,j"
                     
                     if (j==1) then
                        
                        if(i==2.or.i==(lft_endNode(2)+2)) then
                           node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,1)
                        else
                           node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,m)
                        endif
                        
                        !write(*,*) i,j,node_spr(i,cnt),cnt,"i,j,node_spr,cnt"
                        cnt = cnt+1
                        
                     elseif (j==2) then
                        
                        if (trmnl==1) then
                           node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,1)
                           trmnl=trmnl+1
                        elseif (trmnl==2) then
                           node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,m)
                           trmnl=1
                        endif
                        
                        !write(*,*) i,j,node_spr(i,cnt),cnt,"i,j,node_spr,cnt"
                        cnt = cnt+1
                        
                     elseif (j==3) then
                        node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,1)
                        !write(*,*) i,j,node_spr(i,cnt),cnt,"i,j,node_spr,cnt"
                        cnt = cnt+1
                        
                     endif
                     
                  endif
                  
               enddo
               
            elseif ((i.gt.rght_endNode(2)).and.(i.le.crght_top_endNode(3))) then
               
               cnt = 1
               
               do j = 1,jmax
                  spr_Prv = node_sprS(i,j)
                  
                  if (trmnl_intrmdSpr(spr_Prv,0)==1) then
                     node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,1)
                     !write(*,*) i,j,node_spr(i,cnt),cnt,"i,j,node_spr,cnt"
                     cnt = cnt+1
                     
                  else
                     m = trmnl_intrmdSpr(spr_Prv,0)
                     !write(*,*) m,i,j,"m,i,j"
                     
                     if (i==clft_top_endNode(2).or.i==crght_top_endNode(2)) then
                        node_spr(i,cnt)=trmnl_intrmdSpr(spr_Prv,1)
                        cnt = cnt+1
                        
                     elseif (i==clft_top_endNode(3).or.i==crght_top_endNode(3))then
                        
                        if (j==1.or.j==2) then
                           node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,m)
                           !write(*,*) i,j,node_spr(i,cnt),cnt,"i,j,node_spr,cnt"
                           cnt = cnt+1
                        elseif (j==3) then
                           node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,1)
                           !write(*,*) i,j,node_spr(i,cnt),cnt,"i,j,node_spr,cnt"
                           cnt = cnt+1
                        endif
                        
                     endif
                  endif
                  
               enddo
               
               
            elseif ((i.gt.crght_top_endNode(3)) .and. (i.le.N_nodeS)) then
               
               cnt = 1
               
               do j = 1,jmax
                  spr_Prv = node_sprS(i,j)
                  
                  if (trmnl_intrmdSpr(spr_Prv,0)==1) then
                     node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,1)
                     !write(*,*) i,j,node_spr(i,cnt),cnt,"i,j,node_spr,cnt"
                     cnt = cnt+1
                     
                  else
                     m = trmnl_intrmdSpr(spr_Prv,0)
                     
                     if (j==1.or.j==4) then
                        
                        if ((i.ne.clft_bot_endNode(1)).and.(i.ne.crght_bot_endNode(1))) then
                           node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,m)
                        elseif ((i==clft_bot_endNode(1)).or.(i==crght_bot_endNode(1))) then
                           node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,1)
                        endif
                        !if (i==35) write(*,*) trmnl_intrmdSpr(spr_Prv,0:2)
                        !write(*,*) i,j,node_spr(i,cnt),cnt,"i,j,node_spr,cnt"
                        cnt = cnt+1
                        
                        
                     elseif (j==2.or.j==5) then
                        
                        if (trmnl==1) then
                           node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,1)
                           
                           if ((i.ne.clft_bot_endNode(1)).and.(i.ne.crght_bot_endNode(1))) then
                              if(j==5) trmnl=trmnl+1 !othrwise will be incresd twice
                           elseif ((i==clft_bot_endNode(1)).or.(i==crght_bot_endNode(1))) then
                              if(j==2) trmnl=trmnl+1 !j wont reach upto 5 here
                           endif
                           
                        elseif (trmnl==2) then
                           node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,m)
                           trmnl=1
                        endif
                        
                        write(*,*) i,j,node_spr(i,cnt),cnt,"i,j,node_spr,cnt"
                        cnt = cnt+1
                        
                     elseif (j==3.or.j==6) then
                        
                        if ((i.ne.clft_bot_endNode(2)).and.(i.ne.crght_bot_endNode(2))) then
                           node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,1)
                        elseif (i==clft_bot_endNode(2)) then
                           node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,1)
                        elseif (i==crght_bot_endNode(2)) then
                           node_spr(i,cnt) = trmnl_intrmdSpr(spr_Prv,m)
                        endif
                        
                        !write(*,*) i,j,node_spr(i,cnt),cnt,"i,j,node_spr,cnt"
                        cnt = cnt+1
                        
                     endif
                     
                  endif
                  
               enddo
               
            endif
         endif
         
      enddo
      
    end subroutine NI_models_trmnlNode_to_spr
    
    subroutine Inserted_node_to_spr
      implicit none
      integer :: i,j,rs,re
      integer :: spr_Prv
      integer :: cnt
      integer :: spr_rplcd(1:(NAEC+1))
      
      cnt = N_nodeS
      
      do i = 1,N_curve
         rs = 1 ; re = 2
         spr_Prv = curve_spr(i)
         call get_rplcd_spr(spr_Prv,spr_rplcd)
         
         do j = 1,NAEC
            node_spr((cnt+j),0)   = 2
            node_spr((cnt+j),1:2) = spr_rplcd(rs:re)
            rs = rs+1 ; re = re+1
         enddo
         
         cnt = cnt+NAEC
         
      enddo

    end subroutine Inserted_node_to_spr
    
    subroutine get_rplcd_spr(spr_Prv,spr_rplcd)
      implicit none
      integer, intent(in)  :: spr_Prv
      integer, intent(out) :: spr_rplcd(1:(NAEC+1))

      integer :: i

      do i = 1,(NAEC+1)
         spr_rplcd(i) = trmnl_intrmdSpr(spr_Prv,i)
      enddo
      
    end subroutine get_rplcd_spr
    
    subroutine trnsfrmtn_adjstmnt_node_spr_for_doubleNode
      implicit none
      integer :: i
      integer :: frst_node,scnd_node
      integer :: valued_pos1,valued_pos2
      integer :: N_valued_pos
      
      frst_node = 0 ; scnd_node = 0
      valued_pos1 = 0 ; valued_pos2 = 0
      
      do i = 1,N_doubleNode
         
         if(i.gt.2) then !Will be Edited, can be done for all double_Nodes
            frst_node = double_node(i,1)
            scnd_node = double_node(i,2)

            valued_pos1 = node_spr(frst_node,0)
            valued_pos2 = node_spr(scnd_node,0)
            N_valued_pos =  valued_pos1 + valued_pos2
            write(*,*) N_valued_pos, "N_valued_pos"
            !node_spr(frst_node,((valued_pos1+1):(valued_pos1+valued_pos2))) = node_spr(scnd_node,1:valued_pos2)
            node_spr(frst_node,(valued_pos1+1):N_valued_pos) = node_spr(scnd_node,1:valued_pos2)
            
            !node_spr(scnd_node,1:(valued_pos1+valued_pos2)) = node_spr(frst_node,1:(valued_pos1+valued_pos2))
            node_spr(scnd_node,1:N_valued_pos) = node_spr(frst_node,1:N_valued_pos)
            
            node_spr(frst_node,0) = node_spr(frst_node,0) + node_spr(scnd_node,0)
            node_spr(scnd_node,0) = node_spr(frst_node,0)
         endif
         
      enddo

    end subroutine trnsfrmtn_adjstmnt_node_spr_for_doubleNode

    
    subroutine print_node_to_spr_transfrm
      implicit none
      integer :: i
      
      open(unit=21,file='node_spr.dat')
      open(unit=20,file='node_sprWOtxt.dat')
      
      write(21,*) N_node,"node_to_spr"
      write(21,*) "node-spr,first one node_spr(i,0)"
      
      do i = 1,N_node
         write(21,*) "node to spr for i =",i,"is",node_spr(i,0:max_spr_node)
         write(20,*) i,node_spr(i,0:max_spr_node)
      enddo
      
      close(21)
      close(20)
      
    end subroutine print_node_to_spr_transfrm
    
  end subroutine get_node_to_spr_transfrm !!!!!!

  

  subroutine get_node_to_area_transfrm
    implicit none
    integer :: i,lp_cnt
    integer :: count_vs

    count_nodes = 0
    count_vs    = 0
    count_area  = 0
    
    lp_cnt = wntbb(0)
    
    do i = 1,max_blck
       
       if (wntbb(i) == 1) then
          call lft_node_to_area
       elseif (wntbb(i) == 2) then
          call rght_node_to_area
       elseif (wntbb(i) == 3) then
          call clft_top_node_to_area
       elseif (wntbb(i) == 4) then
          call additnl_node_to_area
       elseif (wntbb(i) == 5) then
          call crght_top_node_to_area
       elseif (wntbb(i) == 6) then
          call clft_bot_node_to_area
       elseif (wntbb(i) == 7) then
          call crght_bot_node_to_area
       elseif (wntbb(i) == 8) then
          continue
       elseif (wntbb(i) == 9) then
          continue
       else
          continue
       endif
       
    enddo
    
    if (modelID==1) then
       call trnsfrmtn_adjstmnt_node_area_for_doubleNode
    elseif (modelID==2) then
       
       call trnsfrmtn_adjstmnt_node_area_for_doubleNode
       
       if (SystemTyp==1) then
          node_area(N_nodeS-1,0:max_area_node) = node_areaS(N_nodeS-1,0:max_area_node)
          node_area(N_nodeS,0:max_area_node)   = node_areaS(N_nodeS,0:max_area_node)
       endif
       
       call Inserted_node_to_area
       
    endif
    
    call print_node_to_area_trnsfrm    
    
  contains
    
    subroutine lft_node_to_area
      implicit none
      integer :: i,j
      
      write(*,*) nvsl,"nvsl"
      
      do i = 1,nvsl
         count_vs   = count_vs + 1
         count_area = count_area + 1
         
         if (i==1) then
            do j = (2*count_vs-1),(2*count_vs)
               node_area(j,0) = 1
               node_area(j,1) = count_area
            enddo
            
         elseif ((i.ne.1) .and. (i.ne.nvsl)) then
            do j = (2*count_vs-1),(2*count_vs)
              node_area(j,0) = 2
              node_area(j,1) = count_area-1
              node_area(j,2) = count_area
           enddo
           
        elseif (i==nvsl) then
           
           if (stageNo==1 .and. stageType==1) then
              
              do j = (2*count_vs-1),(2*count_vs) 
                 node_area(j,0) = 2
                 node_area(j,1) = count_area-1
                 node_area(j,2) = cntrl_cell_nmbr
              enddo
              
              count_area = count_area - 1 !count_area = 2 now
              
           elseif (stageNo==2 .or. stageNo==3) then
              write(*,*) "BRING MY ATTENTION to node_to_other,lft_node_to_area"
           elseif (stageNo==4) then
              
              do j = (2*count_vs-1),(2*count_vs)
                 if (clft_top_cell_nmbr .ne. 0) then   
                    node_area(j,0) = 2
                    node_area(j,1) = count_area-1
                    node_area(j,2) = clft_top_cell_nmbr
                 elseif (clft_top_cell_nmbr == 0) then
                    node_area(j,0) = 1
                    node_area(j,1) = count_area-1
                 endif
              enddo
              count_area = count_area - 1 !count_area = 2 now
              
           endif
        endif
     enddo
     
     count_nodes = count_nodes + 2*nvsl
     !write(*,*) count_nodes,"count_nodes"
     
    end subroutine lft_node_to_area

    subroutine rght_node_to_area
      implicit  none
      integer :: i,j
      
      write(*,*) nvsr,"nvsr"
      
      do i = 1,nvsr
         count_vs   = count_vs + 1
         count_area = count_area + 1
         !write(*,*) count_vs,count_area,"count_vs,count_area"
         
         if (i==1) then
            do j = (2*count_vs-1),(2*count_vs)
               node_area(j,0) = 1
               node_area(j,1) = count_area
            enddo
         elseif ((i.ne.1) .and. (i.ne.nvsr)) then
            do j = (2*count_vs-1),(2*count_vs)
               node_area(j,0) = 2
               node_area(j,1) = count_area-1
               node_area(j,2) = count_area
            enddo
            
         elseif (i==nvsr) then
            write(*,*) (count_area-1),cntrl_cell_nmbr,"count_area,cntrl_cell_nmbr"
            write(*,*) j,stageNo,stageType,"ST"
            
            if (stageNo==1 .and. stageType==1) then
               
               do j = (2*count_vs-1),(2*count_vs)   
                  node_area(j,0) = 2
                  node_area(j,1) = count_area-1
                  node_area(j,2) = cntrl_cell_nmbr
               enddo
               
            elseif (stageNo==2 .or. stageNo==3) then
               write(*,*) "BRING MY ATTENTION to node_to_othr,rght_node_to_area"
            elseif (stageNo==4) then
               
               do j = (2*count_vs-1),(2*count_vs)
                  if (crght_top_cell_nmbr .ne. 0) then   
                     node_area(j,0) = 2
                     node_area(j,1) = count_area-1
                     node_area(j,2) = crght_top_cell_nmbr
                  elseif (crght_top_cell_nmbr == 0) then
                     node_area(j,0) = 1
                     node_area(j,1) = count_area-1
                  endif
               enddo
               
            endif
            
         endif
         
      enddo

      count_nodes = count_nodes + 2*nvsr
      !write(*,*) count_nodes,"count_nodes"
      
    end subroutine rght_node_to_area
    
    subroutine clft_top_node_to_area
      implicit none
      integer :: i

      do i = 1,clft_top_node
         count_nodes = count_nodes + 1
         write(*,*) count_nodes,"count_nodes"
         !write(*,*) wntbb(0:9),"wntbb"
         !stop
         if (wntbb(4).ne.0) then

            if (count_nodes .ne. clft_top_endNode(3))then!clft_top_endNode(3)=15
               node_area(count_nodes,0) = 2 !if possible to be written better
               node_area(count_nodes,1) = clft_top_cell_nmbr
               node_area(count_nodes,2) = crght_top_cell_nmbr
            else
               node_area(count_nodes,0) = 1
               node_area(count_nodes,1) = clft_top_cell_nmbr
            endif
            
         elseif(wntbb(4)==0 .AND. wntbb(6).ne.0) then
            
            if (i==1) then
               node_area(count_nodes,0) = 2 
               node_area(count_nodes,1) = clft_top_cell_nmbr
               node_area(count_nodes,2) = crght_top_cell_nmbr
               
            elseif (i==2) then
               node_area(count_nodes,0) = 4
               node_area(count_nodes,1) = clft_top_cell_nmbr
               node_area(count_nodes,2) = crght_top_cell_nmbr
               node_area(count_nodes,3) = crght_top_cell_nmbr + 1
               node_area(count_nodes,4) = crght_top_cell_nmbr + 2
               
            elseif (i==3) then
               node_area(count_nodes,0) = 2 
               node_area(count_nodes,1) = clft_top_cell_nmbr
               node_area(count_nodes,2) = clft_top_cell_nmbr + 2
               
            endif
            
         endif
         
      enddo
      
    end subroutine clft_top_node_to_area
    
    
    
    subroutine additnl_node_to_area
      implicit none
      
      count_nodes = count_nodes + 1
      
      node_area(count_nodes,0) = 0
      
    end subroutine additnl_node_to_area
    
    
    
    subroutine crght_top_node_to_area
      implicit none
      integer :: i
      
      do i = 1,crght_top_node
         count_nodes = count_nodes + 1
         
         if (wntbb(4).ne.0) then
            if (count_nodes .ne. crght_top_endNode(3)) then !crght_top_endNode(2)=19
               node_area(count_nodes,0) = 2
               node_area(count_nodes,1) = clft_top_cell_nmbr
               node_area(count_nodes,2) = crght_top_cell_nmbr
            else
               node_area(count_nodes,0) = 1
               node_area(count_nodes,1) = crght_top_cell_nmbr
            endif
            
         elseif (wntbb(4)==0 .AND. wntbb(6).ne.0) then
            
            if (i==1) then
               node_area(count_nodes,0) = 2 
               node_area(count_nodes,1) = clft_top_cell_nmbr
               node_area(count_nodes,2) = crght_top_cell_nmbr
               
            elseif (i==2) then
               node_area(count_nodes,0) = 4
               node_area(count_nodes,1) = clft_top_cell_nmbr
               node_area(count_nodes,2) = crght_top_cell_nmbr
               node_area(count_nodes,3) = crght_top_cell_nmbr + 1
               node_area(count_nodes,4) = crght_top_cell_nmbr + 2
               
            elseif (i==3) then
               node_area(count_nodes,0) = 2 
               node_area(count_nodes,1) = crght_top_cell_nmbr
               node_area(count_nodes,2) = crght_top_cell_nmbr + 2
               
            endif
            
         endif
      enddo
      
    end subroutine crght_top_node_to_area
    
    
    subroutine clft_bot_node_to_area
      implicit none
      integer :: i,j
      
      write(*,*) nvsclb,"nvsclb"

      count_vs = nvsl + nvsr + nvsclt + nvscrt
      count_area = ncl + ncr + ncclt + nccrt
      
      do i = 1,nvsclb
         !count_vs   = count_vs + 1
         
         if(i==1) then
            count_vs   = count_vs + 1
            count_area = count_area + 1
         else
            count_vs   = count_vs + 2
            count_area = count_area + 2
         endif
            
         if (i.ne.nvsclb) then
            
            do j = (2*count_vs+1),(2*count_vs+2)
               node_area(j,0) = 2
               node_area(j,1) = count_area !count_area-1 in lft
               node_area(j,2) = count_area + 2
            enddo
            
         elseif (i==nvsclb) then
            
            do j = (2*count_vs+1),(2*count_vs+2)
               
               if (stageNo==4 .and. stageType==1) then
                  node_area(j,0) = 2
                  node_area(j,1) = count_area
                  node_area(j,2) = N_cell
                  
               elseif (stageNo==4 .and. stageType==2) then
                  node_area(j,0) = 1
                  node_area(j,1) = count_area
                  
               endif
               
            enddo
           
            !count_area = count_area - 1 !count_area = 2 now

        endif
        
     enddo
     
     count_nodes = 2*(count_vs) + 2*(nvsclb-1)
     !write(*,*) count_nodes,"count_nodes"
      
      
    end subroutine clft_bot_node_to_area
    
    
    
    subroutine crght_bot_node_to_area
      implicit none
      integer :: i,j
      
      write(*,*) nvscrb,"nvscrb"
      
      count_vs = nvsl + nvsr + nvsclt + nvscrt + 1
      count_area = ncl + ncr + ncclt + nccrt + 1
      
      do i = 1,nvscrb
         !count_vs   = count_vs + 1
         
         if(i==1) then
            count_vs   = count_vs + 1
            count_area = count_area + 1
         else
            count_vs   = count_vs + 2
            count_area = count_area +2
         endif
         
         !write(*,*) count_vs,count_area,"count_vs,count_area"
         
         if (i.ne.nvscrb) then
            
            do j = (2*count_vs+1),(2*count_vs+2)
               node_area(j,0) = 2
               node_area(j,1) = count_area !count_area-1 in lft
               node_area(j,2) = count_area + 2
            enddo
            
         elseif (i==nvsclb) then
            
            do j = (2*count_vs+1),(2*count_vs+2)
               
               if (stageNo==4 .and. stageType==1) then
                  node_area(j,0) = 2
                  node_area(j,1) = count_area
                  node_area(j,2) = N_cell
                  
               elseif (stageNo==4 .and. stageType==2) then
                  node_area(j,0) = 1
                  node_area(j,1) = count_area
                  
               endif
               
            enddo
            
            !count_area = count_area - 1 !count_area = 2 now
            
         endif
         
      enddo
      
      count_nodes =  2*(count_vs) + 2*(nvscrb-1)
      !write(*,*) count_nodes,"count_nodes"
      
      
    end subroutine crght_bot_node_to_area
    
    subroutine inserted_node_to_area
      implicit none
      integer :: cnt_Node,spr_Prv
      integer :: i,j
      
      cnt_Node = N_nodeS

      node_area(N_nodeS,0:max_area_node) = node_areaS(N_nodeS,0:max_area_node)

      open(unit=655,file='inserted_node_to_area.dat')
      
      do i = 1,N_curve
         spr_Prv = curve_spr(i)
         write(655,*) spr_Prv,"spr_Prv"
         
         do j = 1,NAEC
            node_area((cnt_Node+j),0:2) = spr_areaS(spr_Prv,0:2)
            write(655,*) spr_areaS(spr_Prv,0:2),"spr_areaS"
            write(655,*) node_area((cnt_Node+j),0:2),"node_area"
            
         enddo
         cnt_Node = cnt_Node+NAEC
         
      enddo
      
      close(655)
      
    end subroutine inserted_node_to_area
    
    subroutine trnsfrmtn_adjstmnt_node_area_for_doubleNode
      implicit none
      integer :: i
      integer :: frst_node,scnd_node
      integer :: valued_pos1,valued_pos2
      integer :: N_valued_pos
      
      frst_node   = 0 ; scnd_node   = 0
      valued_pos1 = 0 ; valued_pos2 = 0

      if (stageNo==4) then
         
         do i = 1,N_doubleNode
            
            if (i.gt.2 .and. i.ne.N_doubleNode) then
               
               frst_node = double_node(i,1)
               scnd_node = double_node(i,2)
               
               valued_pos1  = node_area(frst_node,0)
               valued_pos2  = node_area(scnd_node,0)
               N_valued_pos = valued_pos1 + valued_pos2
               
               node_area(frst_node,(valued_pos1+1):N_valued_pos) = node_area(scnd_node,1:valued_pos2)
               
               node_area(scnd_node,1:N_valued_pos) = node_area(frst_node,1:N_valued_pos)
               
               node_area(frst_node,0)=node_area(frst_node,0)+node_area(scnd_node,0)
               node_area(scnd_node,0)=node_area(frst_node,0)
               
               
            elseif (i.gt.2 .and. i.eq.N_doubleNode) then
               
               frst_node = double_node(i,1)
               scnd_node = double_node(i,2)
               
               if (stageType==1) then
                  
                  N_valued_pos = 3
                  
                  node_area(frst_node,0) = N_valued_pos
                  node_area(frst_node,2) = node_area(scnd_node,1)
                  node_area(frst_node,3) = node_area(scnd_node,2)
                  
                  node_area(scnd_node,0:N_valued_pos) = node_area(frst_node,0:N_valued_pos)
                  
               elseif (stageType==2) then
                  
                  N_valued_pos = 2
                  
                  node_area(frst_node,0) = N_valued_pos
                  node_area(frst_node,2) = node_area(scnd_node,1)
                  
                  node_area(scnd_node,0:N_valued_pos) = node_area(frst_node,0:N_valued_pos)
                  
               endif
               
            endif
            
         enddo
         
      elseif (stageNo==1 .and. stageType==1) then
         
         do i = 1,N_doubleNode
            
            if (CellsMeet==1) then
               
               if (i.gt.1) then
                  write(*,*) "amount of doubleNode > 1, fl: node_to_other_trnsfrm, sb: trnsfrmtn_adjstmnt_node_area_dn"
                  stop
               endif
               
               frst_node = double_node(i,1)
               scnd_node = double_node(i,2)
               
               N_valued_pos =  3
               
               node_area(frst_node,0) = N_valued_pos
               node_area(frst_node,2) = node_area(scnd_node,1)
               node_area(frst_node,3) = node_area(scnd_node,2)
               
               node_area(scnd_node,0:N_valued_pos) = node_area(frst_node,0:N_valued_pos)
               
            elseif (CellsMeet.gt.1) then
               
               if (i.ne.N_doubleNode) then
                  
                  frst_node = double_node(i,1)
                  scnd_node = double_node(i,2)
                  
                  valued_pos1  = node_area(frst_node,0)
                  valued_pos2  = node_area(scnd_node,0)
                  N_valued_pos = valued_pos1 + valued_pos2
                  
                  node_area(frst_node,(valued_pos1+1):N_valued_pos) = node_area(scnd_node,1:valued_pos2)
                  node_area(scnd_node,1:N_valued_pos)               = node_area(frst_node,1:N_valued_pos)
                  
                  node_area(frst_node,0)=node_area(frst_node,0)+node_area(scnd_node,0)
                  node_area(scnd_node,0)=node_area(frst_node,0)
                  
               elseif (i==N_doubleNode) then
                  
                  frst_node = double_node(i,1)
                  scnd_node = double_node(i,2)
                  
                  N_valued_pos =  3
                  
                  node_area(frst_node,0) = N_valued_pos
                  node_area(frst_node,2) = node_area(scnd_node,1)
                  node_area(frst_node,3) = node_area(scnd_node,2)
               
                  node_area(scnd_node,0:N_valued_pos) = node_area(frst_node,0:N_valued_pos)
                  
               endif
               
            endif
            
         enddo
         
      endif
      
    end subroutine trnsfrmtn_adjstmnt_node_area_for_doubleNode
    
    
    subroutine print_node_to_area_trnsfrm
      implicit none
      integer :: i
      
      open(unit=22,file='node_area.dat')
      open(unit=23,file='node_areaWOtxt.dat')
      
      write(22,*) "node-area","   ","first one node_area(i,0)"
      
      do i = 1,N_node
         write(22,*) "node to area for i =",i,"are",node_area(i,0:max_area_node)
         write(23,*) i,node_area(i,0:max_area_node)
      enddo
      
      close(22)
      close(23)
      
    end subroutine print_node_to_area_trnsfrm
    
  end subroutine get_node_to_area_transfrm
  
end module node_to_other_transfrm





