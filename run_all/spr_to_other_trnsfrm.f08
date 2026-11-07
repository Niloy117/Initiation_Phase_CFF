module spr_to_other_transfrm
  use transfrm_info !transfrm_variables.f08
  
  implicit none
  
contains
  
  subroutine get_spr_to_node_transfrm
    implicit none
    
    if (modelID==1) then
       call get_spr_to_node_transfrm_TN_model
    elseif (modelID==2) then
       call get_spr_to_node_transfrm_NI_model
    endif
    
  end subroutine get_spr_to_node_transfrm
  
  subroutine get_spr_to_node_transfrm_TN_model
    implicit none
    integer :: lp_cnt
    integer :: cnnctng_node
    integer :: CMval,nslStr,nsrStr
    
    count_spr = 0
    
    do lp_cnt = 1,max_blck
       
       if (wntbb(lp_cnt) == 1) then
          call lft_spr_to_node
       elseif (wntbb(lp_cnt) == 2) then
          call rght_spr_to_node
       elseif (wntbb(lp_cnt) == 3) then
          call clft_top_spr_to_node
       elseif (wntbb(lp_cnt) == 4) then
          call get_cnnctng_node(wntbb(lp_cnt-1),cnnctng_node)
          call additnl_spr_to_node(cnnctng_node)
       elseif (wntbb(lp_cnt) == 5) then
          call crght_top_spr_to_node
       elseif (wntbb(lp_cnt) == 6) then
          call clft_bot_spr_to_node
       elseif (wntbb(lp_cnt) == 7) then
          call crght_bot_spr_to_node
       elseif (wntbb(lp_cnt)==8) then
          call cntrlRgn_spr_to_node
       elseif (wntbb(lp_cnt)==9) then
          call endRgn_spr_to_node
       else
          continue
       endif
       
    enddo
       
    call print_spr_to_node_transfrm
       
  contains
    
    subroutine lft_spr_to_node
      implicit none
      integer :: i
      integer :: cell_no
      
      spr_node((count_spr+1):(count_spr+nsl),0) = 2
      
      do i = 1,nsl
         
         count_spr = count_spr + 1
         
         if (mod(i,3) .ne. 0) then
            cell_no = (count_spr/3) + 1
         elseif (mod(i,3) == 0)  then
            cell_no = (count_spr/3)  
         endif
         
         if (mod(i,3) == 1) then
            spr_node(count_spr,1) = 1 + (cell_no-1)*2
            spr_node(count_spr,2) = 3 + (cell_no-1)*2
         elseif (mod(i,3) == 2) then
            spr_node(count_spr,1) = 2 + (cell_no-1)*2        
            spr_node(count_spr,2) = 4 + (cell_no-1)*2
         else
            spr_node(count_spr,1) = 3 + (cell_no-1)*2
            spr_node(count_spr,2) = 4 + (cell_no-1)*2
         endif
         
      enddo
      
    end subroutine lft_spr_to_node
    
    subroutine rght_spr_to_node
      implicit none
      integer :: i
      integer :: cell_no
      
      spr_node((count_spr+1):(count_spr+nsr),0) = 2
      
      do i = 1,nsr
         count_spr = count_spr + 1
         
         if (mod(i,3) .ne. 0) then
            cell_no = (count_spr/3) + 1
         elseif (mod(i,3) == 0)  then
            cell_no = (count_spr/3)  
         endif
         
         if (mod(i,3) == 1) then
            spr_node(count_spr,1) = 1 + cell_no*2
            spr_node(count_spr,2) = 3 + cell_no*2
         elseif (mod(i,3) == 2) then
            spr_node(count_spr,1) = 2 + cell_no*2
            spr_node(count_spr,2) = 4 + cell_no*2
         else
            spr_node(count_spr,1) = 3 + cell_no*2
            spr_node(count_spr,2) = 4 + cell_no*2
         endif
         
      enddo
      
    end subroutine rght_spr_to_node
    
    subroutine clft_top_spr_to_node
      implicit none
      integer :: i
      
      spr_node((count_spr+1),0)                          = 3
      spr_node((count_spr+2):(count_spr+clft_top_spr),0) = 2
      
      do i = 1,clft_top_spr !clft_top_spr = nsclt
         count_spr = count_spr + 1
         
         if (i == 1) then
            spr_node(count_spr,1) = lft_endNode(1)
            spr_node(count_spr,2) = clft_top_endNode(1)
            spr_node(count_spr,3) = clft_top_endNode(2)
         elseif (i == 2) then
            spr_node(count_spr,1) = clft_top_endNode(2)
            spr_node(count_spr,2) = clft_top_endNode(3)
         elseif (i == 3) then
            spr_node(count_spr,1) = clft_top_endNode(3)
            spr_node(count_spr,2) = lft_endNode(2)
         endif
         
      enddo
      
    end subroutine clft_top_spr_to_node

    subroutine additnl_spr_to_node(cnnctng_node)
      implicit none
      !integer :: i
      integer :: cnnctng_node

      !write(*,*) clft_top_end_spr,"clft_top_end_spr"
      count_spr = count_spr + 1
      
      spr_node(count_spr,0) = 2
      spr_node(count_spr,1) = cnnctng_node !blck_node + 2
      spr_node(count_spr,2) = additnl_node_endNode(1)
      
    end subroutine additnl_spr_to_node

    subroutine crght_top_spr_to_node
      implicit none
      integer :: i

      spr_node((count_spr+1),0)                           = 3
      spr_node((count_spr+2):(count_spr+crght_top_spr),0) = 2

      !write(*,*) rght_end_node,"rght_end_node"
      
      do i = 1,crght_top_spr !crght_top_spr = nscrt
         count_spr = count_spr + 1
         
         if (i == 1) then
            spr_node(count_spr,1) = rght_endNode(1)
            spr_node(count_spr,2) = crght_top_endNode(1)
            spr_node(count_spr,3) = crght_top_endNode(2)
         elseif (i == 2) then
            spr_node(count_spr,1) = crght_top_endNode(2)
            spr_node(count_spr,2) = crght_top_endNode(3)
         elseif (i == 3) then
            spr_node(count_spr,1) = crght_top_endNode(3)
            spr_node(count_spr,2) = rght_endNode(2)
         endif
         
      enddo
      
    end subroutine crght_top_spr_to_node



    subroutine clft_bot_spr_to_node
      implicit none
      integer :: i,j
      integer :: cell_no
      integer :: prv_cell_cnnctng_node(1:2)
      integer :: prv_blckCell,cell_diff

      cell_no = -1
      prv_blckCell = ncl + ncr + ncclt + nccrt
      cell_diff = -1
      
      !write(*,*) count_spr,"clft_bot_spr_to_node"
      
      prv_cell_cnnctng_node(1:2) = 0


      do i = 1,ncclb

         do j = 1,nsecclb
            count_spr = count_spr + 1
            !write(*,*) j,count_spr,"count_spr"
            spr_node(count_spr,0) = 2
            
            if (mod(count_spr,3) .ne. 0) then
               cell_no = (count_spr/3) + 1
            elseif (mod(count_spr,3) == 0)  then
               cell_no = (count_spr/3)  
            endif
            
            cell_diff = cell_no - prv_blckCell
            
            if (cell_diff .lt. 0) then
               write(*,*) "Warning : Cell_diff cant be negative, sbrtn: clft_bot_spr_to_node"
               stop
            endif
            
            if (mod(count_spr,3) == 1) then
               
               if (cell_diff==1) then
                  spr_node(count_spr,1) = clft_top_endNode(2)
               else
                  spr_node(count_spr,1) = prv_cell_cnnctng_node(1)
               endif
               
               spr_node(count_spr,2)    = 3 + (cell_no+1)*2!lft_blck its cell_no-1 
               prv_cell_cnnctng_node(1) = spr_node(count_spr,2)

               !write(*,*) prv_cell_cnnctng_node(1),"prv_cell_cnnctng(1)"
               
               !write(*,*) count_spr,"count_spr"
               !write(*,*) spr_node(count_spr,1:2),"spr_node"
               
            elseif (mod(count_spr,3) == 2) then
               
               if (cell_diff==1) then
                  spr_node(count_spr,1) = clft_top_endNode(3)
                  !write(*,*) clft_top_endNode(3),"clft_top_endNode(3)"
               else
                  spr_node(count_spr,1) = prv_cell_cnnctng_node(2)
                  !write(*,*) spr_node(count_spr,2),count_spr,"count_spr=26"
               endif
               
               spr_node(count_spr,2)    = 4 + (cell_no+1)*2
               prv_cell_cnnctng_node(2) = spr_node(count_spr,2)

               !write(*,*) prv_cell_cnnctng_node(2),"prv_cell_cnnctng(2)"
               
            elseif (mod(count_spr,3) == 0) then
               spr_node(count_spr,1) = 3 + (cell_no+1)*2
               spr_node(count_spr,2) = 4 + (cell_no+1)*2
            endif
            
         enddo
         count_spr = count_spr + 3
      enddo
      
    end subroutine clft_bot_spr_to_node
    
    
    
    subroutine crght_bot_spr_to_node
      implicit none
      integer :: i,j
      integer :: cell_no
      integer :: prv_cell_cnnctng_node(1:2)
      integer :: prv_blckCell,cell_diff

      count_spr = nsl + nsr + nsclt + nscrt + nsecclb
      
      cell_no = -1
      prv_blckCell = ncl + ncr + ncclt + nccrt
      cell_diff = -1
      
      write(*,*) count_spr,"crght_bot_spr_to_node"
      
      prv_cell_cnnctng_node(1:2) = 0


      do i = 1,nccrb
      
         do j = 1,nseccrb
            count_spr = count_spr + 1
            !write(*,*) i,count_spr,"count_spr"
            
            spr_node(count_spr,0) = 2
            
            if (mod(count_spr,3) .ne. 0) then
               cell_no = (count_spr/3) + 1
            elseif (mod(count_spr,3) == 0)  then
               cell_no = (count_spr/3)  
            endif
            
            cell_diff = cell_no - prv_blckCell
            
            if (cell_diff .lt. 0) then
               write(*,*) "Warning : Cell_diff cant be negative, sbrtn: crght_bot_spr_to_node"
               stop
            endif
            
            
            if (mod(count_spr,3) == 1) then
               
               if (cell_diff==2) then
                  spr_node(count_spr,1) = crght_top_endNode(2)
               else
                  spr_node(count_spr,1) = prv_cell_cnnctng_node(1)
               endif
               
               spr_node(count_spr,2)    = 3 + (cell_no+1)*2!lft_blck its cell_no-1 
               prv_cell_cnnctng_node(1) = spr_node(count_spr,2)
               
               !write(*,*) count_spr,"count_spr"
               !write(*,*) spr_node(count_spr,1:2),"spr_node"
               
            elseif (mod(count_spr,3) == 2) then
               
               if (cell_diff==2) then
                  spr_node(count_spr,1) = crght_top_endNode(3)
                  !write(*,*) crght_top_endNode(3),"crght_top_endNode(3)"
               else
                  spr_node(count_spr,1) = prv_cell_cnnctng_node(2)
               endif
               
               spr_node(count_spr,2)   = 4 + (cell_no+1)*2
               prv_cell_cnnctng_node(2) = spr_node(count_spr,2)
               
            elseif (mod(count_spr,3) == 0) then
               spr_node(count_spr,1) = 3 + (cell_no+1)*2
               spr_node(count_spr,2) = 4 + (cell_no+1)*2
            endif
            
         enddo
         count_spr = count_spr + 3
      enddo
      
    end subroutine crght_bot_spr_to_node

    subroutine cntrlRgn_spr_to_node
      implicit none
      integer :: i,j

      if (stageNo==1 .and. stageType==1) then
         continue
      else
         stop
         write(*,*) "diff stg-diff typ,fl:spr_to_other,sb:cntrlRgn_spr_to_node "
      endif
      
      count_spr = nsl + nsr
      
      write(*,*) lft_endNode,rght_endNode,"lft-rght end node"
      
      do i = 1,ncc
         
         do j = 1,nsecc
            count_spr = count_spr+1
            spr_node(count_spr,0) = 2 
            
            if (j==1) then
               spr_node(count_spr,1) = lft_endNode(1)
               spr_node(count_spr,2) = rght_endNode(1)
            elseif (j==2) then
               spr_node(count_spr,1) = lft_endNode(2)
               spr_node(count_spr,2) = rght_endNode(2)
            endif
            
         enddo
         
      enddo
      
    end subroutine cntrlRgn_spr_to_node
    
    subroutine endRgn_spr_to_node
      implicit none
      integer :: i
      
      count_spr = nsl + nsr + nsclt + nscrt + nsclb + nscrb
      
      do i = 1,nser
         count_spr = count_spr+1
         
         if (i==1) then
            spr_node(count_spr,0) = 2
            spr_node(count_spr,1) = N_node-2
            spr_node(count_spr,2) = N_node
         endif
         !write(*,*)spr_node(count_spr,0:max_node_spr),count_spr,"endRgn spr_nde"
         
      enddo
      
    end subroutine endRgn_spr_to_node
    
    subroutine springless_boundary_to_node
      implicit none
      integer :: i
      
      sprless_node(1:N_nospr,0) = 2
      
      do i = 1,N_nospr
         if (i==1) then
            sprless_node(i,1) = 1
            sprless_node(i,2) = 2
         elseif (i==N_nospr) then

            if (wntbb(0) == 1) then !only lft_blck
               sprless_node(i,1) = lft_endNode(2)
               sprless_node(i,2) = lft_endNode(2)
            else
               sprless_node(i,1) = lft_endNode(2) + 1
               sprless_node(i,2) = lft_endNode(2) + 2
            endif
            
         endif
      enddo
      
    end subroutine springless_boundary_to_node
    
  end subroutine get_spr_to_node_transfrm_TN_model !!!!!!
  
  
  subroutine get_spr_to_node_transfrm_NI_model
    implicit none
    
    if (NI_incldAS_woPP == 0) call get_spr_to_node_transfrm_NI_model_methd1
    if (NI_incldAS_woPP == 1) call get_spr_to_node_transfrm_NI_model_methd2

       
         
  end subroutine get_spr_to_node_transfrm_NI_model
  
  subroutine get_spr_to_node_transfrm_NI_model_methd1
    implicit none
    integer :: i,j,jmax,jEnd,k
    integer :: spr_nm,typ_nm
    integer :: cnt,cnt_Pspr !Pspr=P
    integer :: spr_Prv,spr_Curr
    integer :: cnt_addedNode
    
    jmax = 3
    spr_node = -1
    spr_Prv = 0 ; spr_Curr = 0
    cnt_addedNode = 1
    
    do i = 1,N_cell
       
       if (i==N_cell) then
          
          if (stageNo==1 .and. stageType==1) then
             if (CellsMeet==0)   jEnd = 2
             if (CellsMeet.gt.0) jEnd = 1
          elseif (stageNo==4 .and. stageType==1) then
             jEnd = 1
          elseif (stageNo==4 .and. stageType==2) then
             jEnd = 3
          endif
          
       elseif (i.ne.N_cell) then
          jEnd = jmax
       endif
       
       do j = 1,jEnd
          spr_Prv = (i-1)*jmax + j
          typ_nm  = typ_sprStr(spr_Prv)
          
          if (i.le.(ncl+ncr)) then
             
             if (typ_nm==2 .or. typ_nm==4 .or. typ_nm==5) then
                
                do k = 1,(NAEC+1)
                   
                   if (k==1) then
                      spr_node(spr_Curr+k,0) = 2
                      spr_node(spr_Curr+k,1) = spr_nodeS(spr_Prv,1)
                      spr_node(spr_Curr+k,2) = N_nodeS+(cnt_addedNode)
                      
                      cnt_addedNode = cnt_addedNode+1

                   elseif (k.gt.1 .and. k.ne.(NAEC+1)) then
                      spr_node(spr_Curr+k,0) = 2
                      spr_node(spr_Curr+k,1) = spr_node(spr_Curr+(k-1),2) !not spr_nodeS
                      spr_node(spr_Curr+k,2) = N_nodeS+(cnt_addedNode)
                      
                      cnt_addedNode = cnt_addedNode+1
                      
                   elseif (k==(NAEC+1)) then
                      spr_node(spr_Curr+k,0) = 2
                      spr_node(spr_Curr+k,1) = spr_node(spr_Curr+(k-1),2)
                      spr_node(spr_Curr+k,2) = spr_nodeS(spr_Prv,2)

                      spr_Curr = spr_Curr+(NAEC+1)
                   endif
                   
                enddo
                
             else
                !spr_node(spr_curr+1,0) = 2 
                !spr_node(spr_curr+1,1:2) = spr_nodeS(spr_Prv,1:2)
                !spr_Curr = spr_Curr+1
                spr_node(spr_Curr+1,0:max_node_spr) = spr_nodeS(spr_Prv,0:max_node_spr)
                spr_Curr = spr_Curr+1
             endif
             
          elseif (i.gt.(ncl+ncr) .and. i.le.(ncl+ncr+2)) then
             
             if (stageNo==1 .and. stageType==1) then

                if (typ_nm==2 .or. typ_nm==4 .or. typ_nm==5) then
                   
                   do k = 1,(NAEC+1)
                      
                      if (k==1) then
                         spr_node(spr_Curr+k,0) = 2
                         spr_node(spr_Curr+k,1) = spr_nodeS(spr_Prv,1)
                         spr_node(spr_Curr+k,2) = N_nodeS+(cnt_addedNode)
                         
                         cnt_addedNode = cnt_addedNode+1
                         
                      elseif (k.gt.1 .and. k.ne.(NAEC+1)) then
                         spr_node(spr_Curr+k,0) = 2
                         spr_node(spr_Curr+k,1) = spr_node(spr_Curr+(k-1),2) !not spr_nodeS
                         spr_node(spr_Curr+k,2) = N_nodeS+(cnt_addedNode)
                         
                         cnt_addedNode = cnt_addedNode+1
                         
                      elseif (k==(NAEC+1)) then
                         spr_node(spr_Curr+k,0) = 2
                         spr_node(spr_Curr+k,1) = spr_node(spr_Curr+(k-1),2)
                         spr_node(spr_Curr+k,2) = spr_nodeS(spr_Prv,2)
                         
                         spr_Curr = spr_Curr+(NAEC+1)
                      endif
                      
                   enddo
                
                else
                   !spr_node(spr_curr+1,0) = 2 
                   !spr_node(spr_curr+1,1:2) = spr_nodeS(spr_Prv,1:2)
                   
                   spr_node(spr_Curr+1,0:max_node_spr) = spr_nodeS(spr_Prv,0:max_node_spr)
                   spr_Curr = spr_Curr+1
                endif
                
                
             elseif (stageNo==4) then
                
                if (typ_nm==7) then
                   
                   do k = 1,(NAEC+1)
                      
                      if (k==1) then
                         spr_node(spr_Curr+k,0) = 2
                         spr_node(spr_Curr+k,1) = spr_nodeS(spr_Prv,1)
                         spr_node(spr_Curr+k,2) = N_nodeS+(cnt_addedNode)
                         
                         cnt_addedNode = cnt_addedNode+1
                         
                      elseif (k.gt.1 .and. k.ne.(NAEC+1)) then
                         spr_node(spr_Curr+k,0) = 2
                         spr_node(spr_Curr+k,1) = spr_node(spr_Curr+(k-1),2) !not spr_nodeS
                         spr_node(spr_Curr+k,2) = N_nodeS+(cnt_addedNode)
                         
                         cnt_addedNode = cnt_addedNode+1
                         
                      elseif (k==(NAEC+1)) then
                         spr_node(spr_Curr+k,0) = 2
                         spr_node(spr_Curr+k,1) = spr_node(spr_Curr+(k-1),2)
                         spr_node(spr_Curr+k,2) = spr_nodeS(spr_Prv,2)
                         
                         spr_Curr = spr_Curr+(NAEC+1)
                      endif
                   
                   enddo
                   
                   !spr_node(spr_Curr+1,0) = 2
                   !spr_node(spr_Curr+1,1) = spr_nodeS(spr_Prv,1)
                   !spr_node(spr_Curr+1,2) = N_nodeS+(cnt_addedNode)
                   
                   !cnt_addedNode = cnt_addedNode+1
                   
                   !spr_node(spr_Curr+2,0) = 2
                   !spr_node(spr_Curr+2,1) = spr_node(spr_Curr+1,2) !not spr_nodeS
                   !spr_node(spr_Curr+2,2) = spr_nodeS(spr_Prv,2)
                   
                   !spr_Curr = spr_Curr+2
                   
                elseif (typ_nm==8) then
                   
                   do k = 1,(NAEC+1)
                      
                      if (k==1) then
                         spr_node(spr_Curr+k,0) = 2
                         spr_node(spr_Curr+k,1) = spr_nodeS(spr_Prv,2)
                         spr_node(spr_Curr+k,2) = N_nodeS+(cnt_addedNode)
                         
                         cnt_addedNode = cnt_addedNode+1
                         
                      elseif (k.gt.1 .and. k.ne.(NAEC+1)) then
                         spr_node(spr_Curr+k,0) = 2
                         spr_node(spr_Curr+k,1) = spr_node(spr_Curr+(k-1),2) !not spr_nodeS
                         spr_node(spr_Curr+k,2) = N_nodeS+(cnt_addedNode)
                         
                         cnt_addedNode = cnt_addedNode+1
                         
                      elseif (k==(NAEC+1)) then
                         spr_node(spr_Curr+k,0) = 2
                         spr_node(spr_Curr+k,1) = spr_node(spr_Curr+(k-1),2)
                         spr_node(spr_Curr+k,2) = spr_nodeS(spr_Prv,1)
                         
                         spr_Curr = spr_Curr+(NAEC+1)
                      endif
                      
                   enddo
                   
                else
                   spr_node(spr_Curr+1,0:max_node_spr) = spr_nodeS(spr_Prv,0:max_node_spr)
                   spr_Curr = spr_Curr+1
                endif
                
             endif
             
          elseif (i.gt.(ncl+ncr+2) .and. i.le.N_cell) then
             
             if (typ_nm==4 .or. typ_nm==5) then

                do k = 1,(NAEC+1)
                   
                   if (k==1) then
                      spr_node(spr_Curr+k,0) = 2
                      spr_node(spr_Curr+k,1) = spr_nodeS(spr_Prv,1)
                      spr_node(spr_Curr+k,2) = N_nodeS+(cnt_addedNode)
                      
                      cnt_addedNode = cnt_addedNode+1
                      
                   elseif (k.gt.1 .and. k.ne.(NAEC+1)) then
                      spr_node(spr_Curr+k,0) = 2
                      spr_node(spr_Curr+k,1) = spr_node(spr_Curr+(k-1),2) !not spr_nodeS
                      spr_node(spr_Curr+k,2) = N_nodeS+(cnt_addedNode)
                      
                      cnt_addedNode = cnt_addedNode+1
                      
                   elseif (k==(NAEC+1)) then
                      spr_node(spr_Curr+k,0) = 2
                      spr_node(spr_Curr+k,1) = spr_node(spr_Curr+(k-1),2)
                      spr_node(spr_Curr+k,2) = spr_nodeS(spr_Prv,2)
                      
                      spr_Curr = spr_Curr+(NAEC+1)
                   endif
                   
                enddo
                
                !spr_node(spr_Curr+1,0) = 2
                !spr_node(spr_Curr+1,1) = spr_nodeS(spr_Prv,1)
                !spr_node(spr_Curr+1,2) = N_nodeS+(cnt_addedNode)
                
                !cnt_addedNode = cnt_addedNode+1
                
                !spr_node(spr_Curr+2,0) = 2
                !spr_node(spr_Curr+2,1) = spr_node(spr_Curr+1,2) !not spr_nodeS
                !spr_node(spr_Curr+2,2) = spr_nodeS(spr_Prv,2)
                
                !spr_Curr = spr_Curr+2
                
             else
                spr_node(spr_curr+1,0) = 2 
                spr_node(spr_curr+1,1:2) = spr_nodeS(spr_Prv,1:2)
                
                spr_Curr = spr_Curr+1
             endif
             
          endif
          !write(*,*) spr_Curr,"spr_Curr"
       enddo
    enddo
    
    call print_spr_to_node_transfrm
    
  end subroutine get_spr_to_node_transfrm_NI_model_methd1
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  subroutine get_spr_to_node_transfrm_NI_model_methd2
    implicit none
    integer :: i,j,jmax,jEnd,k
    integer :: spr_nm,typ_nm
    integer :: cnt,cnt_Pspr      ! Pspr=P
    integer :: spr_Prv,spr_Curr
    integer :: cnt_addedNode,chk_val
    
    jmax          =  3
    spr_node      = -1
    spr_Prv       =  0  ; spr_Curr = 0
    cnt_addedNode =  1
    
    do i = 1,N_cell
       
       if (i==N_cell) then
          if (CellsMeet==0)   jEnd = 2
          if (CellsMeet.gt.0) jEnd = 1
       elseif (i.ne.N_cell) then
          jEnd = jmax
       endif
       
       do j = 1,jEnd
          spr_Prv = (i-1)*jmax + j
          typ_nm  = typ_sprStr(spr_Prv)
          
          chk_val=0 ; call typ_spr_chk(typ_nm,chk_val)
          
          if (chk_val == 1) then
             
             do k = 1,(NAEC+1)
                
                if (k==1) then
                   spr_node(spr_Curr+k,0) = 2
                   spr_node(spr_Curr+k,1) = spr_nodeS(spr_Prv,1)
                   spr_node(spr_Curr+k,2) = N_nodeS+(cnt_addedNode)
                   
                   cnt_addedNode = cnt_addedNode+1
                   
                elseif (k.gt.1 .and. k.ne.(NAEC+1)) then
                   spr_node(spr_Curr+k,0) = 2
                   spr_node(spr_Curr+k,1) = spr_node(spr_Curr+(k-1),2) !not spr_nodeS
                   spr_node(spr_Curr+k,2) = N_nodeS+(cnt_addedNode)
                   
                   cnt_addedNode = cnt_addedNode+1
                   
                elseif (k==(NAEC+1)) then
                   spr_node(spr_Curr+k,0) = 2
                   spr_node(spr_Curr+k,1) = spr_node(spr_Curr+(k-1),2)
                   spr_node(spr_Curr+k,2) = spr_nodeS(spr_Prv,2)
                   
                   spr_Curr = spr_Curr+(NAEC+1)
                endif
             enddo
             
          elseif (chk_val==0) then
             spr_node(spr_Curr+1,0:max_node_spr) = spr_nodeS(spr_Prv,0:max_node_spr)
             spr_Curr = spr_Curr+1
          endif
          
       enddo
    enddo
    
    
  end subroutine get_spr_to_node_transfrm_NI_model_methd2  
  
  
  subroutine typ_spr_chk(typ_sprVal,chk_val)
    implicit none
    integer, intent(in)  :: typ_sprVal
    integer, intent(out) :: chk_val
    
    integer :: cnt_apcl,cnt_bsal,cnt_ltrl
    integer :: chk_valBfr,i
    
    cnt_apcl=0 ; cnt_bsal=0 ; cnt_ltrl=0
    
    if (NAEC_Apcl.gt.0) cnt_apcl=1
    if (NAEC_Bsal.gt.0) cnt_bsal=1
    if (NAEC_Ltrl.gt.0) cnt_ltrl=1
    
    chk_val = 0 ; chk_valBfr = chk_val
    
    do i = 1,3
       
       if (i==1) then
          if (cnt_apcl==1) then
             if (typ_sprVal==1 .or. typ_sprVal==3) then
                chk_val=1
                exit
             endif
          endif
       elseif (i==2) then
          if (cnt_bsal==1) then
             if (typ_sprVal==2 .or. typ_sprVal==4) then
                chk_val=1
                exit
             endif
          endif
       elseif (i==3) then
          if (cnt_ltrl==1) then
             if (typ_sprVal==5) then
                chk_val=1
                exit
             endif
          endif
       endif
       
    enddo
    
    write(*,*) typ_sprVal,chk_val,chk_valBfr,"typ"
    
  end subroutine typ_spr_chk
  
  
  
  subroutine print_spr_to_node_transfrm
    implicit none
    integer :: i
      
    open(unit=22,file='spr_node.dat')
    open(unit=23,file='spr_nodeWOtxt.dat')
    
    write(22,*) "spr-node","   ","first one spr_node(i,0)"
    !write(*,*) N_spr,"spr_to_node"
    
    do i = 1,N_spr
       write(22,*) "spr to node for i =",i,"are",spr_node(i,0:max_node_spr)
       write(23,*) i,spr_node(i,0:max_node_spr)
    enddo
    
    close(22)
    close(23)
    
  end subroutine print_spr_to_node_transfrm
  
  subroutine get_spr_to_area_transfrm
    implicit none
    
    if (modelID==1) then
       call get_spr_to_area_transfrm_TN_model
    elseif (modelID==2) then
       call get_spr_to_area_transfrm_NI_model
    endif
    
  end subroutine get_spr_to_area_transfrm
  
  subroutine get_spr_to_area_transfrm_TN_model
    implicit none
    integer :: lp_cnt
    
    count_spr = 0

    !write(*,*) count_area,"inside_area_transform"

    !lp_cnt = wntbb(0)
    
    do lp_cnt = 1,max_blck
    
       if (wntbb(lp_cnt) == 1) then
          call lft_spr_to_area
       elseif (wntbb(lp_cnt) == 2) then
          call rght_spr_to_area
       elseif (wntbb(lp_cnt) == 3) then
          call clft_top_spr_to_area
       elseif (wntbb(lp_cnt) == 4) then
          !No_area
          call additnl_spr_to_area
       elseif (wntbb(lp_cnt) == 5) then
          call crght_top_spr_to_area
       elseif (wntbb(lp_cnt) == 6) then
          call clft_bot_spr_to_area
       elseif (wntbb(lp_cnt) == 7) then
          call crght_bot_spr_to_area
       elseif (wntbb(lp_cnt) == 8) then
          call cntrlRgn_spr_to_area
       elseif (wntbb(lp_cnt) == 9) then
          call endRgn_spr_to_area
       else
          continue
       endif
       
    enddo
    
    call print_spr_to_area_trnsfrm
    
  contains
    
    subroutine lft_spr_to_area
      implicit none
      integer :: i
      !integer :: cell_no
      
      !write(*,*) clft_top_cell_nmbr,"clft_top_cell_nmbr"
      
      do i = 1,nsl
         count_spr = count_spr + 1
         
         if (mod(i,3) .ne. 0) then
            spr_area(count_spr,0) = 1
            spr_area(count_spr,1) = (count_spr/3) + 1
            
         elseif (mod(i,3) == 0) then

            if ((i/3) .ne. ncl) then
               spr_area(count_spr,0) = 2
               spr_area(count_spr,1) = (count_spr/3) 
               spr_area(count_spr,2) = (count_spr/3) + 1
               
            elseif ((i/3) == ncl) then
               
               if (stageNo==1 .and. stageType==1) then
                  
                  if (cntrl_cell_nmbr .le. 0) then
                     write(*,*)"fl:spr_to_othr,sb:lft_spr_to_area"
                     stop
                  endif
                  
                  spr_area(count_spr,0) = 2
                  spr_area(count_spr,1) = (count_spr/3)
                  spr_area(count_spr,2) = cntrl_cell_nmbr
                  
                  
               elseif (stageNo==4) then
                  
                  if (clft_top_cell_nmbr .le. 0) then
                     write(*,*)"fl:spr_to_othr,sb:lft_spr_to_area"
                     stop
                  endif
                  
                  spr_area(count_spr,0) = 2
                  spr_area(count_spr,1) = (count_spr/3)
                  spr_area(count_spr,2) = clft_top_cell_nmbr
                  
               else
                  
                  spr_area(count_spr,0) = 1
                  spr_area(count_spr,1) = (count_spr/3)
                  
               endif
               
            endif
         endif   
      enddo
      
    end subroutine lft_spr_to_area

    subroutine rght_spr_to_area
      implicit  none
      integer :: i
      
      do i = 1,nsr
         count_spr = count_spr + 1
         
         if (mod(i,3) .ne. 0) then
            spr_area(count_spr,0) = 1
            spr_area(count_spr,1) = (count_spr/3) + 1
            
         elseif (mod(count_spr,3) == 0) then
            
            if ((i/3) .ne. ncr) then
               spr_area(count_spr,0) = 2
               spr_area(count_spr,1) = (count_spr/3) 
               spr_area(count_spr,2) = (count_spr/3) + 1
               
            elseif ((i/3) == ncr) then

               if (stageNo==1 .and. stageType==1) then
                  
                  if (cntrl_cell_nmbr .le. 0) then
                     write(*,*)"fl:spr_to_othr,sb:rght_spr_to_area"
                     stop
                  endif
                  
                  spr_area(count_spr,0) = 2
                  spr_area(count_spr,1) = (count_spr/3)
                  spr_area(count_spr,2) = cntrl_cell_nmbr
                  
                  
               elseif (stageNo==4) then
                  
                  if (clft_top_cell_nmbr .le. 0) then
                     write(*,*)"fl:spr_to_othr,sb:rght_spr_to_area"
                     stop
                  endif
                  
                  spr_area(count_spr,0) = 2
                  spr_area(count_spr,1) = (count_spr/3)
                  spr_area(count_spr,2) = crght_top_cell_nmbr
                  
               else
                  
                  spr_area(count_spr,0) = 1
                  spr_area(count_spr,1) = (count_spr/3)
                  
               endif
               
            endif
            
         endif
      enddo
      
    end subroutine rght_spr_to_area
    
    subroutine clft_top_spr_to_area
      implicit none
      integer :: i
      
      !count_spr = count_spr + 1
      !write(*,*) count_spr,clft_top_spr,"count_spr,clft_top_spr"
      
      !spr_area(count_spr:(count_spr + (clft_top_spr-1)), 0) = 1
      
      !write(*,*) spr_area(13:15,0),"spr_area(13:15,0)"
      !write(*,*) clft_top_cell_nmbr,"clft_top_cell_nmbr"
      
      !spr_area(count_spr:(count_spr + (clft_top_spr-1)), 1) =clft_top_cell_nmbr
      
      !count_spr = count_spr + (nscl-1)

      do i = 1,nsecclt
         count_spr = count_spr + 1
         
         if (mod(i,3)==1 .or. mod(i,3)==0) then
            spr_area(count_spr,0) = 1

            if (mod(i,3)==1) then
               spr_area(count_spr,1) = (count_spr/3) + 1
            elseif (mod(i,3)==0) then
               spr_area(count_spr,1) = (count_spr/3)
            endif
            
         elseif (mod(i,3)==2) then

            if (wntbb(6) .ne. 0) then
               spr_area(count_spr,0) = 2
               spr_area(count_spr,1) = clft_top_cell_nmbr
               spr_area(count_spr,2) = clft_top_cell_nmbr + 2
               
            elseif (wntbb(6) == 0) then
               spr_area(count_spr,0) = 1
               spr_area(count_spr,1) = clft_top_cell_nmbr
            endif
            
         endif
         
      enddo
      
      
    end subroutine clft_top_spr_to_area

    
    subroutine additnl_spr_to_area
      implicit none

      count_spr = count_spr + 1
      spr_area(count_spr, 0) = 0
      
    end subroutine additnl_spr_to_area

    subroutine crght_top_spr_to_area
      implicit none
      integer :: i
      
      !count_spr = count_spr + 1
      
      !spr_area(count_spr:(count_spr + (crght_top_spr-1)), 0) = 1
      !spr_area(count_spr:(count_spr + (crght_topspr-1)), 1)=crght_top_cell_nmbr

      !count_spr = count_spr + (nscr-1)

      do i = 1,nseccrt
         count_spr = count_spr + 1
         
         if (mod(i,3)==1 .or. mod(i,3)==0) then
            spr_area(count_spr,0) = 1
            
            if (mod(i,3)==1) then
               spr_area(count_spr,1) = (count_spr/3) + 1
            elseif (mod(i,3)==0) then
               spr_area(count_spr,1) = (count_spr/3)
            endif
            
         elseif (mod(i,3)==2) then

            if (wntbb(7) .ne. 0) then
               spr_area(count_spr,0) = 2
               spr_area(count_spr,1) = crght_top_cell_nmbr
               spr_area(count_spr,2) = crght_top_cell_nmbr + 2
               
            elseif (wntbb(6) == 0) then
               spr_area(count_spr,0) = 1
               spr_area(count_spr,1) = crght_top_cell_nmbr
               
            endif
            
         endif
         
      enddo
      
    end subroutine crght_top_spr_to_area
    
    
    subroutine clft_bot_spr_to_area
      implicit none
      integer :: i
      
      !write(*,*)  clft_bot_cell_nmbr,"clft_bot_cell"
      write(*,*) clft_bot_cell_nmbr,"clbcn"
      
      do i = 1,nsclb
         if (i.ne.1 .and. mod(i,nsecclb)==1) then !nsecclb=3
            count_spr = count_spr + 4
         else
            count_spr = count_spr + 1 
         endif
         
         if (mod(i,3) .ne. 0) then
            spr_area(count_spr,0) = 1
            spr_area(count_spr,1) = (count_spr/3) + 1
            
         elseif (mod(i,3) == 0) then
            !write(*,*) count_spr,"count_spr_in clft_bot_spr_to_area"
            
            if ((count_spr/3) .ne. clft_bot_cell_nmbr ) then
               spr_area(count_spr,0) = 2
               spr_area(count_spr,1) = (count_spr/3) 
               spr_area(count_spr,2) = (count_spr/3) + 2 !its + 1 in lft_blcks
               
            elseif ((count_spr/3) == clft_bot_cell_nmbr) then
               write(*,*) clft_bot_cell_nmbr,"clbcn"

               if (stageNo==4 .and. stageType==1) then
                  spr_area(count_spr,0) = 2
                  spr_area(count_spr,1) = (count_spr/3)
                  spr_area(count_spr,2) = N_cell
                  
               elseif (stageNo==4 .and. stageType==2) then
                  spr_area(count_spr,0) = 1
                  spr_area(count_spr,1) = (count_spr/3)
               endif
               
            endif
            
         endif
         
      enddo
      
      
    end subroutine clft_bot_spr_to_area
    
    
    subroutine crght_bot_spr_to_area
      implicit none
      integer :: i

      count_spr = nsl + nsr + nsclt + nscrt + nsecclb
      
      do i = 1,nscrb
        
         if (i.ne.1 .and. mod(i,nseccrb)==1) then !nseccrb=3
            count_spr = count_spr + 4
         else
            count_spr = count_spr + 1
         endif
         
         if (mod(i,3) .ne. 0) then
            spr_area(count_spr,0) = 1
            spr_area(count_spr,1) = (count_spr/3) + 1
            
         elseif (mod(i,3) == 0) then

            if ((count_spr/3) .ne. crght_bot_cell_nmbr) then
               spr_area(count_spr,0) = 2
               spr_area(count_spr,1) = (count_spr/3) 
               spr_area(count_spr,2) = (count_spr/3) + 2 !its + 1 in lft_blcks
               
            elseif ((count_spr/3) == crght_bot_cell_nmbr) then
               write(*,*) crght_bot_cell_nmbr,"crbcn"

               if (stageNo==4 .and. stageType==1) then
                  spr_area(count_spr,0) = 2
                  spr_area(count_spr,1) = (count_spr/3)
                  spr_area(count_spr,2) = N_cell
                  
               elseif (stageNo==4 .and. stageType==2) then  
                  spr_area(count_spr,0) = 1
                  spr_area(count_spr,1) = (count_spr/3)
               endif
               
            endif
            
         endif
         
      enddo
      
    end subroutine crght_bot_spr_to_area
    
    
    subroutine cntrlRgn_spr_to_area
      implicit none
      integer :: i
      
      !write(*,*) nsc,"nsc"
      !call sleep(2)

      count_spr = nsl+nsr
      
      do i = 1,nsc
         count_spr = count_spr+1
         
         spr_area(count_spr,0) = 1
         spr_area(count_spr,1) = cntrl_cell_nmbr
      enddo
      
    end subroutine cntrlRgn_spr_to_area
    
    
    
    subroutine endRgn_spr_to_area
      implicit none
      integer :: i
      
      count_spr = nsl + nsr + nsclt + nscrt + nsclb + nscrb
      
      do i = 1,nser
         count_spr = count_spr+1
         
         if (i==1) then
            spr_area(count_spr,0) = 1
            spr_area(count_spr,1) = N_cell
         endif
         
      enddo
      
    end subroutine endRgn_spr_to_area
    
  end subroutine get_spr_to_area_transfrm_TN_model
  
  subroutine get_spr_to_area_transfrm_NI_model
    implicit none
    integer :: i,j,jmax
    integer :: spr_Prv,spr_Curr
    
    do i = 1,N_sprS
       spr_Prv = i
       jmax = trmnl_intrmdSpr(i,0)
       
       do j = 1,jmax
          spr_Curr = trmnl_intrmdSpr(i,j)
          spr_area(spr_Curr,0:max_area_spr) = spr_areaS(spr_Prv,0:max_area_spr)
       enddo
       
    enddo
    
    call print_spr_to_area_trnsfrm
    
  end subroutine get_spr_to_area_transfrm_NI_model
  
  
  subroutine print_spr_to_area_trnsfrm
    implicit none
    integer :: i
    
    open(unit=22,file='spr_area.dat')
    open(unit=23,file='spr_areaWOtxt.dat')
    
    write(22,*) "spr-area","   ","first one spr_area(i,0)"
    write(22,*) N_spr,"spr_to_area"
    
    do i = 1,N_spr
       write(22,*)"sprToArea for i =",i,"and they're",spr_area(i,0:max_area_spr)
       write(23,*) i,spr_area(i,0:max_area_spr)
    enddo
    
    close(22)
    close(23)
    
  end subroutine print_spr_to_area_trnsfrm
  
end module spr_to_other_transfrm

