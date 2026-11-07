
  
subroutine frprmn(Mc,grd_M,iter,funcChoice,N_variabls,ftol,fret,fpp)
  implicit none
  
  integer,intent(in) :: N_variabls,funcChoice
  real*8, intent(in) :: ftol
  
  real*8  :: Mc(1:N_variabls),grd_M(1:N_variabls)
  real*8  :: dgg,fpp,gam,gg,g(1:(N_variabls)),h(1:(N_variabls))
  integer :: iter
  
  integer :: NMAX,ITMAX
  
  real*8  :: fret
  real*8  :: EPS
  integer :: its,j
  integer :: iVar
  
  !integer :: gam_flg
  real*8  :: max_TOL=1.d20
  
  PARAMETER(ITMAX=500,NMAX=50)
  
  interface

     real*8 function func_Minimizer(mv_pnts,N_variabls,funcChoice)
       implicit none
       integer :: N_variabls,funcChoice
       real*8  :: mv_pnts(1:N_variabls)
       
     end function func_Minimizer
     
  end interface
  
  
  open(unit=71,file='inside_minimizer.dat',position='append')
  write(unit=71,fmt=*) Mc,"what Mc minimizer has taken"
  
  EPS  = 1.d-08
  
  write(unit=71,fmt=*) Mc,"bfr fpp"
  !write(*,*) Mc(1:6),"Mc P1"
  !write(*,*) grd_M(1:6),"grd_M P1"
  fpp = func_Minimizer(Mc,N_variabls,funcChoice)
  !write(*,*) "AFT FPP"
  !write(*,*) Mc(1:6),"Mc P2"
  
  !write(*,*) Mc,"aft fpp"
  
  write(unit=71,fmt=*) fpp,"Energy_calc_with_initial_Mc"
  
  call dfunc(Mc,grd_M,N_variabls,funcChoice)

  !write(*,*) Mc(1:6),"Mc P3"
  !write(*,*) grd_M(1:6),"grd_M P3"
  
  !do j = 1,N_variabls
   !  write(*,*) grd_M(j),"grd_M"
  !enddo
  
  !if (N_variabls==96) stop
  
  write(unit=71,fmt=*) grd_M,"gradient_using_initial_Mc"

  do j = 1,(N_variabls)
     g(j)      = -grd_M(j)
     h(j)      = g(j)
     grd_M(j)  = h(j)
  enddo
  
  !open(unit=72,file='E_vs_Niter.dat')
  open(unit=73,file='fret_fpp.dat',position='append')
  
  do its = 1,ITMAX
     !write(*,*) its,"its"
     iter = its
     
     !if(its.le.2) write(unit=71,fmt=*) iter, "bfr"
     !if(its.le.2) write(unit=71,fmt='(8(f4.2x))') 
     
     write(unit=71,fmt=*) Mc,"Mc bfr_lnmn"
     write(unit=71,fmt=*) grd_M,"grd bfr_lnmn"

     !write(*,*) Mc(1:6),"Mc P4"
     call linmin(Mc,grd_M,fret,N_variabls,funcChoice)
     !write(*,*) Mc(1:6),"Mc P5"
     
     write(unit=71,fmt=*) Mc,"aft_linmin"
     write(unit=71,fmt=*) grd_M,"grd_aft_lnmn"
     
     !if(its.le.2)write(*,'(14(f4.2x))') Mc
     
     write(unit=73,fmt=*) its,fret,fpp,"fret_fpp"
     write(*,*) its,fret,fpp,"fret_fpp"

     !if (abs(fret) .gt. 30.00) then
      !  do iVar = 1,N_variabls
       !    write(*,*) Mc(iVar),grd_M(iVar),iVar,"unusual values"
       ! enddo
     !endif
     
     if (fret.gt.fpp) then
        write(*,*) "fret is getting big,not small"
        !exit
        !fret = fpp

        if (fpp==0.0d0) exit
        if ((fret/fpp) .gt. max_TOL) then
           write(*,*) "(fret/fpp)>max_TOL"
           stop
        endif
     endif
     
     !endif
     
     !if (fret.lt.fpp) then
        if ((2.d0*abs(fret-fpp)) .le. (ftol*(abs(fret)+abs(fpp)+EPS))) return
     !endif
     
     ! if ((2.d0*abs(fret-fpp)) .le. (ftol*(abs(fret)+abs(fpp)+EPS))) then
     !    write(*,*) "exit condn satisfied"
     !    exit
     ! endif
     
     fpp=fret
     
     call dfunc(Mc,grd_M,N_variabls,funcChoice)

     !do j = 1,N_variabls
      !  write(*,*) grd_M(j),"grd_M"
     !enddo
     
     
     gg  = 0.d0
     dgg = 0.d0
     
     !write(*,*) g,"g"
     !write(*,*) grd_M,"grd_M"
     
     do j = 1,(N_variabls)
        !write(*,*) g(j),grdnt(j)
        gg  = gg + g(j)**2
        dgg = dgg + (grd_M(j)+g(j)) * grd_M(j)   !!POLAK-RIBIERE
     enddo
     
     !write(*,*) gg,dgg,"gg,dgg"
     
     
     !if(gg .eq. 0.d0) exit
     if(gg .eq. 0.d0) return

     gam = dgg/gg
     !gam_flg = 0
     
     if (gam .lt. 0.0d0) then
        gam = 0.0d0
     endif
     
     
     open(unit=111,file='gam.dat',position='append')
     write(unit=111,fmt=*) gam,its
     close(111)
     
     !if (gam_flg .eq. 0) then
     
     do j = 1,(N_variabls)
        g(j)      = -grd_M(j)
        h(j)      = g(j) + gam*h(j)
        grd_M(j)  = h(j)
     enddo
     
     !elseif (gam_flg .eq. 1) then
     !grd_M = -1.10d0 + its*0.01
     !endif
     
        
  enddo
  
  if (iter .eq. ITMAX) then
     write(*,*) "WARNING: MAXIMUM ITERATION REACHED"
  endif
  
  write(unit=73,fmt=*) " "
  !write(unit=72,fmt=*) iter,fret

  close(71)
  close(73)
  
  
end subroutine frprmn



subroutine linmin(Mc,grd_M,fret,N_variabls,funcChoice)  !!!Line MINIMIZATION routine   
  implicit none
  integer :: N_variabls,funcChoice
  real*8  :: fret,Mc(1:N_variabls),grd_M(1:N_variabls)!,f1dim,df1dim
  integer :: j
  real*8  :: ax,bx,fa,fb,fx,xmin,xx,tol
  
  ax = 0.01d0
  xx = 1.01d0
  !write(*,*) Mc(1:6),"Mc in linmin1"
  !write(*,*) grd_M(1:6),"grd_M in linmin1"
  
  call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
  !write(*,*) ax,xx,bx,"abx"
  
  !write(*,*) Mc(1:6),"Mc in linmin2"
  !write(*,*) grd_M(1:6),"grd_M in linmin2"
  
  !write(*,*) fret,"bfr_dbrent"
  fret = dbrent(ax,xx,bx,f1dim,df1dim,tol,xmin)
  !write(*,*) fret,"aft_dbrent"
  
  !write(*,*) fret,xmin,"fr_xmn_lnmn"
  !write(*,*) xmin,"xmin"
  
  do j = 1,(N_variabls)
     grd_M(j) = xmin*grd_M(j)
     !write(*,*) grdnt(j)
     Mc(j) = Mc(j) + grd_M(j)
  enddo
  !write(*,*) Mc(1:6),"Mc in linmin3"
  !write(*,*) grd_M(1:6),"grd_M in linmin3"
  
contains
  
  subroutine mnbrak(ax,bx,cx,fa,fb,fc,func)
    implicit none
    real*8 :: ax,bx,cx,fa,fb,fc,func
    real*8 :: temp,denom,num,fu,q,r,u,ulim
    
    real*8,parameter :: GOLD=1.618034, GLIMIT=100., TINY=1.e-16

    !write(*,*) ax,bx,"ax bx"
    fa=func(ax)
    fb=func(bx)
    
    !write(*,*) ax,bx,fa,fb,"fa_fb"
    
    if (abs(fb-fa) .lt. TINY) then
     write(*,*) "fb and fa can not be equal,flnm:Minimizer_restr, sbrtn:mnbrak"
     !stop
     return
    endif
    
    if(fb.gt.fa)then
       temp=ax
       ax=bx
       bx=temp
       temp=fb
       fb=fa
       fa=temp
    endif
    
    cx=bx+GOLD*(bx-ax)
    fc=func(cx)
    
    !write(*,*) ax,bx,fa,fb,"ax_bx_fa_fb"

    do
       if(fb.ge.fc) then
          r=(bx-ax)*(fb-fc)
          q=(bx-cx)*(fb-fa)
          denom=2.0d0*(sign(max(abs(q-r),TINY),q-r))
          num  = (bx-cx)*q-(bx-ax)*r
          u    = bx-(num/denom)
          ulim = bx+GLIMIT*(cx-bx)
          
          open(unit=120,file='u_ulim.dat')
          write(unit=120,fmt=*) u,ulim,"u_ulim"
          close(120)
          
          if(((bx-u)*(u-cx)).gt.0.) then
             fu = func(u)
             if(fu.lt.fc) then
                ax = bx
                fa = fb
                bx = u
                fb = fu
                return
             elseif(fu.gt.fb) then
                cx = u
                fc = fu
                return
             endif
             u  = cx+GOLD*(cx-bx) !!Parabolic fit is of no use,so default magnify
             fu = func(u)
          elseif(((u-ulim)*(ulim-cx)).ge.0.) then
             fu=func(u)
             if(fu.lt.fc) then
                bx = cx
                cx = u
                u  = cx + GOLD*(cx-bx)
                fb = fc
                fc = fu
                fu = func(u)
             endif
          elseif(((u-ulim)*(ulim-cx)).ge.0.) then
             u  = ulim
             fu = func(u)
          else
             u = cx + GOLD*(cx-bx)
             fu = func(u)
          endif
          ax = bx
          bx = cx
          cx = u
          fa = fb
          fb = fc
          fc = fu
       else
          !write(*,*) "BRACKETED"
          exit
       endif
    enddo
    
    !write(*,*) "aft_loop_of_mnbrk"
    
    return
  end subroutine mnbrak
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!MNBRAK ENDS!!!!!!!!!DBRENT STARTS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  function dbrent(ax,bx,cx,f,df,tol,xmin)
    implicit none
    integer,parameter :: ITMAX=100
    real*8,parameter  :: ZEPS=1.e-10
    
    real*8  :: dbrent,ax,bx,cx,tol,xmin
    integer :: iter
    real*8  :: a=0.d0,b=0.d0,d=0.d0,d1=0.d0,d2=0.d0,du=0.d0,dv=0.d0,dw=0.d0,dx=0.d0,e
    real*8  :: fu=0.d0,fv=0.d0,fw=0.d0,fx=0.d0,olde=0.d0,tol1,tol2
    real*8  :: f,df
    real*8  :: u=0.d0,u1=0.d0,u2=0.d0,v=0.d0,w=0.d0,x=0.d0,xm=0.d0
    logical :: ok1,ok2
    integer :: flag1=0,flag2=0,flag3=0
    
    tol=1.e-08
    a = min(ax,cx)
    b = max(ax,cx)    
    v = bx
    w = v
    x = v
    e = 0.0d0
    
    fx = f(x)
    !write(*,*) fx,"fx"
    fv = fx
    fw = fx
    dx = df(x)
    dv = dx
    dw = dx
      
    do iter = 1,100
       xm = 0.5d0*(a+b)
       !write(*,*) xm,a,b,"xm_a_b" !2
       tol1 = tol*abs(x)+ZEPS
       tol2 = 2.d0*tol1
       !write(*,*) iter
       !write(*,*) (x-xm),(tol2-0.5d0*(b-a))
       if(abs(x-xm).le.(tol2-0.5d0*(b-a))) then !1
          !write(*,*) "CONDN1"
          exit
       endif !1
       
       if(abs(e).gt.tol1) then  !2
          !write(*,*) iter,"CONDN2" !3  !1 e enter kore na
          d1=2.d0*(b-a)
          d2=d1
          !write(*,*) d1,d2,"d1,d2" !4
          !write(*,*) dw,dv,dx,"dw,dv,dx" !5
          if(dw.ne.dx) d1=(w-x)*dx/(dx-dw)
          if(dv.ne.dx) d2=(v-x)*dx/(dx-dv)
          !write(*,*) dw,dv,"dw_dv" !6
          u1=x+d1
          u2=x+d2
          !write(*,*) u1,u2,"u1u2" !7
          ok1=((a-u1)*(u1-b).gt.0.d0) .AND. (dx*d1.le.0.d0)
          ok2=((a-u2)*(u2-b).gt.0.d0) .AND. (dx*d2.le.0.d0)
          !write(*,*) ok1,ok2,"ok1ok2" !8
          olde=e
          !write(*,*) olde,e,"olde_e"
          e=d
          !write(*,*) e,d,"e_d"
          if(.not.(ok1 .OR. ok2)) then
             flag1 = 1
             !write(*,*) "THE FLAG1 IS ON"
          elseif (ok1.AND.ok2) then
             if(abs(d1).lt.abs(d2)) then
                d = d1
             else
                d = d2
             endif
          elseif(ok1) then
             d = d1
          else
             d = d2
          endif
          !write(*,*) d !12
          if(flag1.eq.0) then 
             if(abs(d).gt.abs(0.5d0*olde)) then
                flag2 = 1
                !write(*,*) "THE FLAG2 IS ON" !13
             endif
             
             if(flag2.eq.0) then
                u = x + d
                if((u-a).lt.tol2 .OR. (b-u).lt.tol2) d=sign(tol1,xm-x)
                flag3 = 1
                !write(*,*) "THE FLAG3 IS ON"
                !write(*,*) u,"flag2_0_hole_u " !14
             endif
          endif
          
       endif     !2
       
       if(flag2.eq.1..OR.flag2.eq.0.OR.flag1.eq.0.OR.flag1.eq.1.OR.flag3.eq.0) then  !3
          !write(*,*) "CONDN3"
          !write(*,*) iter,"iter"
          !write(*,*) dx,"dx" !17
          if(dx.ge.0.d0) then
             e = a-x
             
             !write(*,*) iter,"+"
             !write(*,*) e,a,x
          else
             e = b-x
             !write(*,*) iter,"-"
             !write(*,*) e,b,x
          endif
          d = 0.5d0*e
          !write(*,*) d,"d" !22
       endif    !3
       !if(flag3.eq.1) then    !4
       !write(*,*) "CRAPP__CONDN4"
       !write(*,*) d,tol1,"d_tol1" !23
       if((abs(d)) .ge. tol1) then
          !write(*,*) "No1"
          u  = x + d
          fu = f(u)
       else
          u = x + sign(tol1,d)
          !write(*,*) "No2" !25
          !write(*,*) u
          fu = f(u)
          !write(*,*) u,x,fu,fx,"u_x_fu_fx" !26
          if(fu.gt.fx) exit
       endif
       !write(*,*) u,fu,"u_fu" !27
       !endif                  !4
       !write(*,*) u,d
       !write(*,*) flag1,flag2,flag3
       du = df(u)
       !write(*,*) u,du,"u_du" !28
       !write(*,*) fu,fx,"fu_fx"
       if(fu.le.fx) then     !!!5
            !write(*,*) "CONDN5",u,x
          if(u.ge.x) then
             a = x
             !write(*,*) a,x,"a_x"
          else
             b = x
             !write(*,*) b,x,"b_x"
          endif
          !write(*,*) w,fw,dw,"w_fw_dw" !33
          v  = w
          fv = fw
          dv = dw
          !write(*,*) v,fv,dv,"v_fv_dv"
          !write(*,*) x,fx,dx,"x_fx_dx"
          w  = x
          fw = fx
          dw = dx
          !write(*,*) w,fw,dw,"w_fw_dw"
          !write(*,*) u,fu,du,"u,fu,du"
          x  = u
          fx = fu
          dx = du
          !write(*,*) x,fx,dx,"x_fx_dx"
       else
          !write(*,*) "CONDN5_else"
          !write(*,*) u,x,"u_x" !40
          if(u.lt.x) then
             !write(*,*) "CONDN5_else11"
             a = u
             !write(*,*) a,"a" !41
          else
             b = u
             !write(*,*) b,"b"
          endif
          
          if((fu.le.fw) .OR. (w.eq.x)) then
             !write(*,*) w,fw,dw,"w_fw_dw"
             v  = w
             fv = fw
             dv = dw
             !write(*,*) v,fv,dv,"v_fv_dv"
             !write(*,*) u,fu,du,"u_fu_du"
             w  = u
             fw = fu
             dw = du
             !write(*,*) w,fw,dw,"w_fw_dw"
          elseif(fu.le.fv .OR. v.eq.x .OR. v.eq.w) then
             !write(*,*) u,fu,du,"u_fu_du"
             v  = u
             fv = fu
             dv = du
             !write(*,*) v,fv,dv,"v_fv_dv" !48
          endif
       endif   !!!5
       !write(*,*) e
       !endif !6
    enddo
    xmin   = x
    dbrent = fx
    !write(*,*) xmin,dbrent,"xmin_dbrent"
    return
    
  end function dbrent
  
 

  real*8 function f1dim(lam)
    implicit none
    real*8 :: lam
    integer :: j
    real*8 :: xt(1:(N_variabls))
    
    interface
       real*8 function func_Minimizer(mv_pnts,N_variabls,funcChoice)
         implicit none
         integer :: N_variabls,funcChoice
         real*8  :: mv_pnts(1:N_variabls)
       end function func_Minimizer    
    end interface
    
    !write(*,*) grd_M,"grd_M_in_f1dim"
    
    do j = 1,N_variabls
       xt(j) = Mc(j) + lam*grd_M(j)
    enddo
    
    !write(*,*) xt(1:6),"xt"
    !write(*,*) Mc(1:6),"Mc in f1dim"
    !write(*,*) lam,"lam"
    !write(*,*) grd_M(1:6),"grd_M"
    
    f1dim = func_Minimizer(xt,N_variabls,funcChoice)   
    !write(*,*) f1dim,"f1dim"
    return
    
  end function f1dim
  
  real*8 function df1dim(lam)
    implicit none
    real*8::lam
    integer :: j
    real*8 :: df(1:N_variabls),xt(1:N_variabls)
    
    
    do j = 1,(N_variabls)
       xt(j) = Mc(j) + lam*grd_M(j)
       !write(*,*) xt(j),Mc(j),grd_M(j),"xt,Mc,grd_M for lam=",lam
    enddo
    
    call dfunc(xt,df,N_variabls,funcChoice)  
    !write(*,*) df,"df"

    df1dim=0.d0
    
    do j=1,N_variabls
       df1dim = df1dim + df(j)*grd_M(j)
    enddo
    !write(*,*) df1dim,"df1dim"
    return
  end function df1dim
  
end subroutine linmin



real*8 function func_Minimizer(mv_pnts,N_variabls,funcChoice)
  implicit none
  integer :: N_variabls,funcChoice
  real*8  :: mv_pnts(1:N_variabls)
  
  if (funcChoice.eq.1) then
     func_Minimizer = Energy_Minimizer(mv_pnts,N_variabls)
  elseif (funcChoice.eq.2) then
     func_Minimizer = Wfunc_Minimizer(mv_pnts,N_variabls)
  endif
  

contains

  real*8 function Energy_Minimizer(dum_coordntes,dum_Nvars)
    use Energy_module
    implicit none
    integer :: dum_Nvars
    real*8  :: dum_coordntes(1:dum_Nvars)
    real*8  :: dum_nodes(1:N_node,1:N_dmnsn)
    integer :: iVal
    real*8  :: Es,Ea
    real*8  :: Es_prv,Ea_prv
    
    if (dum_Nvars.ne.N_mvCoordnte) then
       write(*,*) "N_variabls =/ N_mvCoordnte"
       stop
    endif
    
    !do iVal = 1,10
    !   write(*,*) dum_coordntes(iVal),"coor PB"
    !enddo
    
    call coordntes_to_nodes(dum_coordntes,dum_nodes)
    
    !do iVal = 1,10
    !   write(*,*) dum_nodes(iVal,1:N_dmnsn),ival,"nodeP"
    !enddo
    !do iVal = 1,10
    !   write(*,*) dum_coordntes(iVal),"coor P"
    !enddo
    
    Energy_Minimizer = Energy(dum_nodes,l0,A0)
    
    ! if (abs(Energy_Minimizer).gt.30.00) then
    !    open(unit=45,file='UnusualEval.dat')

    !    Es=0.00d0 ; Ea=0.00d0
    !    Es_prv=0.0d0 ; Ea_prv=0.0d0
       
    !    do iVal = 1,N_node
    !       write(45,*) dum_nodes(iVal,1:N_dmnsn),ival,"nodeV"
    !    enddo
       
    !    do iVal = 1,N_spr
    !       write(45,*) k_spr(iVal),l(iVal),l0(iVal),iVal,"sprV"
    !       Es_prv = Es
    !       Es = Es + 0.5d0*k_spr(iVal)*(l(iVal)-l0(iVal))**2
    !       write(*,*) Es,Es_prv,(Es-Es_prv),iVal,"sprE"
    !    enddo
       
    !    do iVal = 1,N_cell
    !       write(45,*) k_area(iVal),A(iVal),A0(iVal),iVal,"areaV"
    !       Ea_prv = Ea
    !       Ea = Ea + 0.5d0*k_area(iVal)*(A(iVal)-A0(iVal))**2
    !       write(*,*) Ea,Ea_prv,(Ea-Ea_prv),"areaE"
    !    enddo

    !    do iVal = 1,N_node
    !       write(45,*) CgXNode(iVal),iVal,"CgXNode"
    !    enddo
       
    !    do iVal = 1,N_node
    !       write(45,*) k_phi(iVal,1:max_Phi_node),iVal,"k_phi"
    !    enddo
       
    !    close(45)
    ! endif
    
  end function Energy_Minimizer
  
  
  real*8 function Wfunc_Minimizer(dum_coordntes,dum_Nvars)
    use Wfunc_and_its_derivative
    implicit none
    integer :: dum_Nvars
    real*8  :: dum_coordntes(1:dum_Nvars)
    
    if (dum_Nvars.ne.N_mvCoordnte_withl0A0) then
       write(*,*) "N_variabls =/ N_mvCoordnte_withl0A0"
       stop
    endif
    
    Wfunc_Minimizer = Wfunc(dum_coordntes)
    
  end function Wfunc_Minimizer
  
end function func_Minimizer


subroutine dfunc(dum_coordntes,dum_grdC,dum_Nvars,dum_funcC)
  
  implicit none
  integer :: dum_funcC
  integer :: dum_Nvars
  real*8  :: dum_coordntes(1:dum_Nvars)
  real*8  :: dum_grdC(1:dum_Nvars)
  
  if (dum_funcC.eq.1) then
     call grdE_coordntes_Minimizer(dum_coordntes,dum_grdC,dum_Nvars)
  elseif (dum_funcC.eq.2) then
     call grdW_coordntes_Minimizer(dum_coordntes,dum_grdC,dum_Nvars)
  endif
  
end subroutine dfunc


subroutine grdE_coordntes_Minimizer(dum_coordntes,dum_grdC,dum_Nvars)
  use gradient_module
    
  implicit none
  integer :: dum_Nvars
  real*8  :: dum_coordntes(1:dum_Nvars)
  real*8  :: dum_grdC(1:dum_Nvars)
  
  dum_grdC = 0.0d0
  
  !write(*,*) dum_coordntes,"dum_c"
  !write(*,*) dum_grdC,"dum_grdC"
  call gradient_coordntes(dum_coordntes,dum_grdC)
  
  !write(*,*) dum_grdC,"dum_grdC_aft"
  
end subroutine grdE_coordntes_Minimizer



subroutine grdW_coordntes_Minimizer(dum_coordntes,dum_grdC,dum_Nvars)
  use Wfunc_and_its_derivative
  
  implicit none
  integer :: dum_Nvars
  real*8  :: dum_coordntes(1:dum_Nvars)
  real*8  :: dum_grdC(1:dum_Nvars)
  
  dum_grdC = 0.0d0
  
  
  !write(*,*) dum_coordntes,"dum_c"
  !write(*,*) dum_grdC,"dum_grdC"
  call grdW_coordntes(dum_coordntes,dum_grdC)
  
  !write(*,*) dum_grdC,"dum_grdC_aft"
    
end subroutine grdW_coordntes_Minimizer

