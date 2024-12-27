module readjust_struct
  use build_struct
  use nodeInsrtd_cntrlStates
  use Adding_cells
  
contains
  
  subroutine readjust_strctProp
    implicit none
    integer :: i

    if (modelID==2) then
       write(*,*) "routine: readjust_strctProp isn't for modelID=2"
       stop
    endif
    
    open(unit=811,file='ks_strct.dat')
    open(unit=812,file='l0_strct.dat')
    open(unit=813,file='ka_strct.dat')
    open(unit=814,file='A0_strct.dat')
    
    do i = 1,N_spr
       write(811,fmt=*) ks_strctTN(1,i),ks_strctTN(2,i),ks_strctTN(3,i),i
       write(812,fmt=*) l0_strctTN(1,i),l0_strctTN(2,i),l0_strctTN(3,i),i
    enddo
    
    do i = 1,N_cell
       write(813,fmt=*) ka_strctTN(1,i),ka_strctTN(2,i),ka_strctTN(3,i),i
       write(814,fmt=*) A0_strctTN(1,i),A0_strctTN(2,i),A0_strctTN(3,i),i
    enddo
    
    write(811,fmt=*)
    write(812,fmt=*)
    write(813,fmt=*)
    write(814,fmt=*)
    
    call adjstd_strct1
    call adjstd_strct2
    call adjstd_strct3
    
    do i = 1,N_spr
       write(811,fmt=*) ks_strctTN(1,i),ks_strctTN(2,i),ks_strctTN(3,i),i
       write(812,fmt=*) l0_strctTN(1,i),l0_strctTN(2,i),l0_strctTN(3,i),i
    enddo
    
    do i = 1,N_cell
       write(813,fmt=*) ka_strctTN(1,i),ka_strctTN(2,i),ka_strctTN(3,i),i
       write(814,fmt=*) A0_strctTN(1,i),A0_strctTN(2,i),A0_strctTN(3,i),i
    enddo
    
    close(811)
    close(812)
    close(813)
    close(814)
    
  end subroutine readjust_strctProp
  
  subroutine adjstd_strct1
    implicit none
    integer :: strctNum
    
    if (modelID==2) then
       write(*,*) "routine: adjstd_strct1 isn't for modelID=2"
       stop
    endif
    
    strctNum = 1
    
    l0_strctTN(1,1:N_spr)  = l0_strctTN(2,1:N_spr)
    ks_strctTN(1,1:N_spr)  = ks_strctTN(2,1:N_spr)
    A0_strctTN(1,1:N_cell) = A0_strctTN(2,1:N_cell)
    ka_strctTN(1,1:N_cell) = ka_strctTN(2,1:N_cell)
    CgX_strctTN(1,1:N_node) = CgX_strctTN(2,1:N_node)
    
    write(*,*) N_mvCoordnte1,N_mvCoordnte2,"Nmv1,Nmv2"
    write(*,*) N_mvCoordnte,"Nmv"
    
    call coordntes_to_nodes(coordntes_xy,node_xy)
    N_mvCoordnte1 = N_mvCoordnte2
    deallocate(coordntesXY_strct1)
    allocate(coordntesXY_Strct1(1:N_mvCoordnte1))
    
    coordntesXY_strct1(1:N_mvCoordnte1)  = coordntes_xy(1:N_mvCoordnte2)
    nodeXY_strctTN(1,1:N_node,1:N_dmnsn) = node_xy(1:N_node,1:N_dmnsn)
    
  end subroutine adjstd_strct1
  
  subroutine adjstd_strct2
    implicit none
    integer :: areas_Readjst(1:N_cell)
    integer :: sprs_Readjst(1:N_spr)
    
    real*8  :: A0_strct2(1:N_cell),ka_strct2(1:N_cell)
    real*8  :: l0_strct2(1:N_spr),ks_strct2(1:N_spr)
    
    integer :: strctNum
    integer :: i
    integer :: trgt_spr,trgt_cell
    integer :: n1,n2
    real*8  :: CgXV
    
    if (modelID==2) then
       write(*,*) "routine: adjstd_strct2 isn't for modelID=2"
       stop
    endif
    
    strctNum = 2
    
    A0_strct2(1:N_cell) = A0_strctTN(2,1:N_cell)
    ka_strct2(1:N_cell) = ka_strctTN(2,1:N_cell)
    l0_strct2(1:N_spr)  = l0_strctTN(2,1:N_spr)
    ks_strct2(1:N_spr)  = ks_strctTN(2,1:N_spr)
    
    call get_readjstAreas(areas_Readjst,strctNum)!cyclNo be inclded later
    call get_readjstSprs(sprs_Readjst,strctNum)
    
    do i = 1,N_cell
       trgt_cell = areas_Readjst(i)
       A0_strctTN(2,i) = A0_strct2(trgt_cell)
       ka_strctTN(2,i) = ka_strct2(trgt_cell)
       !ka_strctTN(2,i) = ka_strctTN(1,i)
       write(*,*) ka_strctTN(2,i),A0_strctTN(2,i),i,"area"
    enddo
    
    do i = 1,N_spr
       trgt_spr = sprs_Readjst(i)
       l0_strctTN(2,i) = l0_strct2(trgt_spr)
       ks_strctTN(2,i) = ks_strct2(trgt_spr)
       !ks_strctTN(2,i) = ks_strctTN(1,i)
       write(*,*) ks_strctTN(2,i),ks_strct2(trgt_spr),l0_strctTN(2,i),l0_strct2(trgt_spr),i,trgt_spr,"spr"
    enddo
    
    open(unit=308,file='readjst.dat')
    
    CgX_strctTN(2,1:N_node) = 0.0d0
    CgXV = CgXNode(1)
    CgX_strctTN(2,1:2) = CgXV
    n1 = 2*nvsl-1 ; n2 = 2*nvsl
    CgX_strctTN(2,n1:n2) = CgXV 
    write(308,fmt=*) n1,n2,"n1,n2"
    write(308,fmt=*) CgXV,"CgXV"
    
    close(308)
    
  end subroutine adjstd_strct2
  
  subroutine adjstd_strct3
    implicit none
    integer :: areas_Readjst(1:N_cell)
    integer :: sprs_Readjst(1:N_spr)
    
    real*8  :: A0_strct3(1:N_cell),ka_strct3(1:N_cell)
    real*8  :: l0_strct3(1:N_spr),ks_strct3(1:N_spr)
    
    integer :: strctNum
    integer :: i
    integer :: trgt_spr,trgt_cell
    
    if (modelID==2) then
       write(*,*) "routine: adjstd_strct3 isn't for modelID=2"
       stop
    endif
    
    strctNum = 3
    
    A0_strct3(1:N_cell) = A0_strctTN(3,1:N_cell)
    ka_strct3(1:N_cell) = ka_strctTN(3,1:N_cell)
    l0_strct3(1:N_spr)  = l0_strctTN(3,1:N_spr)
    ks_strct3(1:N_spr)  = ks_strctTN(3,1:N_spr)
    
    call get_readjstAreas(areas_Readjst,strctNum)!cyclNo be included latr
    call get_readjstSprs(sprs_Readjst,strctNum)
    
    do i = 1,N_cell
       trgt_cell = areas_Readjst(i)
       A0_strctTN(3,i) = A0_strct3(trgt_cell)
       !ka_strctTN(3,i) = ka_strct3(trgt_cell)
       ka_strctTN(3,i) = ka_strctTN(1,i)
    enddo
    
    do i = 1,N_spr
       trgt_spr = sprs_Readjst(i)
       l0_strctTN(3,i) = l0_strct3(trgt_spr)
       !ks_strctTN(3,i) = ks_strct3(trgt_spr)
       ks_strctTN(3,i) = ks_strctTN(1,i)
    enddo
    
  end subroutine adjstd_strct3
  
  subroutine get_readjstAreas(areas_Readjst,strctNum)
    implicit none
    integer, intent(inout) :: areas_Readjst(1:N_cell)
    integer, intent(in)    :: strctNum
    
    integer :: i
    integer :: PrevNcl,PrevNcr,FutrNcl,FutrNcr
    integer :: CurrNcl,CurrNcr
    integer :: cclr,Futr_cclr
    
    CurrNcl = ncl ; PrevNcl = ncl+1 ; FutrNcl = ncl-1
    CurrNcr = ncr ; PrevNcr = ncr+1 ; FutrNcr = ncr-1
    cclr = CurrNcl + CurrNcr
    Futr_cclr = FutrNcl+FutrNcr 
    
    if (strctNum==2) then
       open(unit=305,file='areasReadjst_strct2.dat')
       
       do i = 1,N_cell
             
          if (i.le.FutrNcl) then
                areas_Readjst(i) = i+1
          elseif (i.gt.FutrNcl.and.i.le.Futr_cclr) then
                areas_Readjst(i) = i+2
          elseif (i.gt.Futr_cclr.and.i.le.(N_cell-2)) then
             areas_Readjst(i) = i+2
          elseif ((i.gt.(N_cell-2)).and.i.le.N_cell) then
             areas_Readjst(i) = i
          endif
          
          write(305,fmt=*)"Radjst area is =",areas_Readjst(i), "for area old =",i
          
       enddo
       
       close(305)
       
    elseif (strctNum==3) then
       open(unit=306,file='areasReadjst_strct3.dat')
       
       do i = 1,N_cell
          
          if (i.le.cclr) then
             
             if ((i.ne.CurrNcl) .and. (i.ne.cclr)) then
                areas_Readjst(i) = i+1
             elseif (i==CurrNcl) then
                areas_Readjst(i) = (ncl+ncr)+1
             elseif (i==cclr) then
                areas_Readjst(i) = (ncl+ncr)+2
             endif
             
          elseif ((i.gt.cclr).and.(i.le.(N_cell-2))) then
             areas_Readjst(i) = i+2
          else
             areas_Readjst(i) = i
          endif
          
          
          write(306,fmt=*)"Radjst area is =",areas_Readjst(i), "for area old =",i
       enddo
       
       close(306)
       
    endif

  end subroutine get_readjstAreas
  
  subroutine get_readjstSprs(sprs_Readjst,strctNum)
    implicit none
    integer, intent(inout) :: sprs_Readjst(1:N_spr)
    integer, intent(in)    :: strctNum
    
    integer :: i
    integer :: PrevNsl,PrevNsr,FutrNsl
    integer :: CurrNsl,CurrNsr,FutrNsr
    integer :: PrevNsCntrlTop
    integer :: cslr,Futr_cslr
    
    CurrNsl = nsl ; PrevNsl = nsl+nsecl*1 ; FutrNsl = nsl-nsecl*1
    CurrNsr = nsr ; PrevNsr = nsr+nsecr*1 ; FutrNsr = nsl-nsecr*1
    
    PrevNsCntrlTop = PrevNsl + PrevNsr + nsclt + nscrt
    cslr = CurrNsl+CurrNsr
    Futr_cslr = FutrNsl+FutrNsr
    
    if (strctNum==2) then
       open(unit=305,file='sprsReadjst_strct2.dat')
       
       do i = 1,N_spr
          
          if (i.le.FutrNsl) then
             sprs_Readjst(i) = i+nsecl
             
          elseif (i.gt.FutrNsl.and.i.le.Futr_cslr) then
             sprs_Readjst(i) = i+nsecl+nsecr
             
          elseif ((i.gt.Futr_cslr).and.(i.le.(N_spr-6))) then
             sprs_Readjst(i) = i+nsecl+nsecr
             
          elseif (i.gt.(N_spr-6).and.i.le.N_spr) then
             sprs_Readjst(i) = i
          endif
          
          write(305,fmt=*) "Radjst spr is =",sprs_Readjst(i), "for area old =",i
          
       enddo
       
       close(305)
       
    elseif (strctNum==3) then
       open(unit=306,file='sprsReadjst_strct3.dat')
       
       do i = 1,N_spr
          
          if (i.le.(CurrNsl-3)) then
             sprs_Readjst(i) = i+nsecl
             
          elseif ((i.gt.(CurrNsl-3)).and.(i.le.CurrNsl)) then
             if (i==(CurrNsl-2)) sprs_Readjst(i) = cslr+1
             if (i==(CurrNsl-1)) sprs_Readjst(i) = cslr+3
             if (i==(CurrNsl-0)) sprs_Readjst(i) = cslr+2
             
          elseif ((i.gt.CurrNsl).and.(i.le.(cslr-3))) then
             sprs_Readjst(i) = i+nsecr
             
          elseif ((i.gt.(cslr-3)).and.i.le.cslr) then
             if (i==(cslr-2)) sprs_Readjst(i) = cslr+nsecclt+1
             if (i==(cslr-1)) sprs_Readjst(i) = cslr+nsecclt+3
             if (i==(cslr-0)) sprs_Readjst(i) = cslr+nsecclt+2
             
          elseif ((i.gt.cslr).and.(i.le.(cslr+3))) then
             if (i==(cslr+1)) sprs_Readjst(i) = i+(2*nsecclt)
             if (i==(cslr+2)) sprs_Readjst(i) = i+(2*nsecclt)+1
             if (i==(cslr+3)) sprs_Readjst(i) = i+(2*nsecclt)-1
             
          elseif ((i.gt.cslr+3).and.(i.le.(cslr+6))) then
             if (i==(cslr+4)) sprs_Readjst(i) = i+(2*nsecclt)
             if (i==(cslr+5)) sprs_Readjst(i) = i+(2*nsecclt)+1
             if (i==(cslr+6)) sprs_Readjst(i) = i+(2*nsecclt)-1
             
          elseif ((i.gt.cslr+6).and.i.le.(N_spr-6)) then
             sprs_Readjst(i) = i+(2*nsecclt)
          elseif (i.gt.(N_spr-6).and.i.le.N_spr) then
             sprs_Readjst(i) = i
          endif
          
          write(306,fmt=*) "Radjst spr is =",sprs_Readjst(i), "for area old =",i
          
       enddo

       close(306)
    endif
    
    
  end subroutine get_readjstSprs
  
  subroutine get_adjstNodes(nodes_Readjst)
    implicit none
    integer, intent(inout) :: nodes_Readjst(1:N_node)
    
    integer :: i
    integer :: CurrNvsl,PrevNvsl
    integer :: CurrNvsr,PrevNvsr
    integer :: CurrNnl,PrevNnl
    integer :: CurrNnr,PrevNnr
    integer :: cnlr,PrevCnlr
    integer :: PrevNnCntrlTop
    
    open(unit=305,file='nodesReadjst.dat')
    
    CurrNvsl = nvsl ; PrevNvsl = nvsl+1
    CurrNvsr = nvsr ; PrevNvsr = nvsr+1

    CurrNnl = 2*CurrNvsl ; PrevNnl = 2*PrevNvsl 
    CurrNnr = 2*CurrNvsr ; PrevNnr = 2*PrevNvsr

    cnlr = CurrNnl+CurrNnr
    PrevCnlr = PrevNnl + PrevNnr
    
    PrevNnCntrlTop = 2*PrevNvsl + 2*PrevNvsr + clft_top_node + crght_top_node

    do i = 1,N_node
       
       if (i.le.CurrNnl) then
          nodes_Readjst(i) = i+2
       endif
       
    enddo
    
    close(305)
    
  end subroutine get_adjstNodes
  
  
end module readjust_struct

! do i = 1,N_spr
   
!    if (i.le.(CurrNsl-3)) then
!       sprs_Readjst(i) = i+nsecl
      
!    elseif ((i.gt.(CurrNsl-3)).and.(i.le.CurrNsl)) then
!       if (i==(CurrNsl-2)) sprs_Readjst(i) = cslr+1
!       if (i==(CurrNsl-1)) sprs_Readjst(i) = cslr+3
!       if (i==(CurrNsl-0)) sprs_Readjst(i) = cslr+2
      
!    elseif ((i.gt.CurrNsl).and.(i.le.(cslr-3))) then
!       sprs_Readjst(i) = i+nsecr
      
!    elseif ((i.gt.(cslr-3)).and.i.le.cslr) then
!       if (i==(cslr-2)) sprs_Readjst(i) = cslr+nsecclt+1
!       if (i==(cslr-1)) sprs_Readjst(i) = cslr+nsecclt+3
!       if (i==(cslr-0)) sprs_Readjst(i) = cslr+nsecclt+2
      
!    elseif ((i.gt.cslr).and.(i.le.(cslr+3))) then
!       if (i==(cslr+1)) sprs_Readjst(i) = i+(2*nsecclt)
!       if (i==(cslr+2)) sprs_Readjst(i) = i+(2*nsecclt)+1
!       if (i==(cslr+3)) sprs_Readjst(i) = i+(2*nsecclt)-1
             
!    elseif ((i.gt.cslr+3).and.(i.le.(cslr+6))) then
!       if (i==(cslr+4)) sprs_Readjst(i) = i+(2*nsecclt)
!       if (i==(cslr+5)) sprs_Readjst(i) = i+(2*nsecclt)+1
!       if (i==(cslr+6)) sprs_Readjst(i) = i+(2*nsecclt)-1
      
!    elseif ((i.gt.cslr+6).and.i.le.(N_spr-6)) then
!       sprs_Readjst(i) = i+(2*nsecclt)
!    elseif (i.gt.(N_spr-6).and.i.le.N_spr) then
!       sprs_Readjst(i) = i
!    endif
   
!    write(305,fmt=*) "Radjst spr is =",sprs_Readjst(i), "for area old =",i
   
! enddo

 
          !if (i.le.CurrNcl) then
           !  areas_Readjst(i) = i+1 
          !elseif (i.gt.CurrNcl.and.i.le.cclr) then
           !  areas_Readjst(i) = i+2
          !elseif ((i.gt.cclr).and.(i.le.(N_cell-2))) then
           !  areas_Readjst(i) = i+2
          !elseif ((i.gt.(N_cell-2)).and.(i.le.N_cell)) then
           !  areas_Readjst(i) = i
!endif

! if (i.le.CurrNsl) then
          !    sprs_Readjst(i) = i+nsecl
          ! elseif (i.gt.CurrNsl .and. i.le.cslr) then
          !    sprs_Readjst(i) = i+nsecl+nsecr
          ! elseif ((i.gt.cslr).and.(i.le.(N_spr-(nsecl+nsecr)))) then
          !    sprs_Readjst(i) = i+nsecl+nsecr
          ! elseif (i.gt.(N_spr-(nsecl+nsecr))) then
          !    sprs_Readjst(i) = i
          ! endif
