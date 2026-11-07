module PenF_minimizer_mod

  implicit none
  
contains

  
  subroutine minimize_PenF_wrt_l0s(hm,nmbrs,spr_or_area)
    use PenF_module
    use system_parameters
    use gradient_module
    !use storing_changing_restoring_routines
    
    implicit none
    integer :: hm
    integer :: nmbrs(1:hm)
    integer :: spr_or_area
    real*8  :: ax,bx,fa,fb,fx,xmin,xx
    real*8  :: tol
    
    real*8  :: PenF_minimized
    integer :: cnt
    integer :: spr_nmbr!,area_nmbr
    real*8  :: l0_val,fctr_val!not same as areas
    
    cnt = 0
    tol = 1e-05

    spr_nmbr = nmbrs(1)
    
    l0_val   = l0(spr_nmbr)
    fctr_val = 1.02d0
    
    ax = l0_val
    xx = l0_val/fctr_val
    !bx = l0_val/(fctr_val)**2  !needed in brak_with_itertn
    
    call mnbrak(ax,xx,bx,fa,fx,fb,PenF_Minmztn)
    write(*,*) ax,xx,bx,"ax_xx_bx_"
    
    ! call brak_with_itertn(ax,xx,bx,fa,fx,fb,PenF_Minmztn,fctr_val)
    ! write(*,*) ax,xx,bx,"ax_xx,bx"
    
    PenF_minimized = brent(ax,xx,bx,PenF_Minmztn,tol,xmin)
    write(*,*) PenF_minimized,"PenF_Minimized"
    
    !call PenF_vs_l0
    
  contains
    
    subroutine mnbrak(ax,bx,cx,fa,fb,fc,func)
      implicit none
      real*8 :: ax,bx,cx,fa,fb,fc,func
      real*8 :: temp,denom,num,fu,q,r,u,ulim
      
      real*8,parameter :: GOLD=1.618034, GLIMIT=100., TINY=1.e-16
      
      
      fa=func(ax)
      fb=func(bx)

      !write(*,*) ax,bx,fa,fb,"ax_bx_fa_fb"
      
      !if (abs(fb-fa) .lt. TINY) then
      !write(*,*) "fb and fa can not be equal,flnm:Minimizer_restr, sbrtn:mnbrak"
      !stop
      !endif
      
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
      !write(*,*) cx,fc,"cx_fc"
      
      !write(*,*) "bfr_loop_in_mnbrk"    
      !write(*,*) ax,bx,fa,fb,"ax_bx_fa_fb"
      
      do
         if(fb.ge.fc) then
            r=(bx-ax)*(fb-fc)
            q=(bx-cx)*(fb-fa)
            denom=2.0d0*(sign(max(abs(q-r),TINY),q-r))
            num  = (bx-cx)*q-(bx-ax)*r
            u    = bx-(num/denom)
            ulim = bx+GLIMIT*(cx-bx)
            
            !open(unit=120,file='u_ulim.dat')
            !write(unit=120,fmt=*) u,ulim,"u_ulim"
            !close(120)
            !write(*,*) u,ulim,"u,u_lim"
            
            if(((bx-u)*(u-cx)).gt.0.) then
               fu = func(u)
               !write(*,*) fu,"fu1 in mnbrak"
               
               if(fu.lt.fc) then
                  ax = bx
                  fa = fb
                  bx = u
                  fb = fu
              
               elseif(fu.gt.fb) then
                  cx = u
                  fc = fu
                  
               endif
               u  = cx+GOLD*(cx-bx) !!Parabolic fit is of no use,so default magnify
               fu = func(u)
               !write(*,*) fu,"fu2 in mnbrak"
               
            elseif(((u-ulim)*(ulim-cx)).ge.0.) then
               fu=func(u)
               !write(*,*) fu,"fu3 in mnbrak"
               
               if(fu.lt.fc) then
                  bx = cx
                  cx = u
                  u  = cx + GOLD*(cx-bx)
                  fb = fc
                  fc = fu
                  fu = func(u)
                  !write(*,*) fu,"fu4 in mnbrak"
                  
               endif
            elseif(((u-ulim)*(ulim-cx)).ge.0.) then
               u  = ulim
               fu = func(u)
               !write(*,*) fu,"fu5 in mnbrak"
               
            else
               u = cx + GOLD*(cx-bx)
               fu = func(u)
               !write(*,*) fu,"fu6 in mnbrak"
               
            endif
            ax = bx
            bx = cx
            cx = u
            fa = fb
            fb = fc
            fc = fu
         else
            !write(*,*) ax,bx,cx,"ax_bx_cx BRACKETED"
            exit
         endif
         
         !write(*,*) ax,bx,cx,"l0s in every step"
      enddo
      
      
      return
    end subroutine mnbrak



    
    real*8 function brent(ax,bx,cx,f,tol,xmin)
      implicit none
      integer :: ITMAX
      real*8  :: ax,bx,cx
      real*8  :: tol,xmin
      real*8  :: CGOLD,ZEPS
      real*8  :: f
      
      PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.0e-10)
      
      integer :: iter
      real*8  :: a,b,d,e,etemp,p,q,r
      real*8  :: u,v,w,x,xm
      real*8  :: fu,fv,fw,fx
      real*8  :: tol1,tol2


      a = min(ax,cx)
      b = max(ax,cx)
      
      v=bx
      w=v
      x=v
      e=0.0d0
      
      fx=f(x)
      !write(*,*) fx,"fx"
      
      fv=fx
      fw=fx
      
      do iter=1,ITMAX
         
         xm   = 0.5*(a+b)
         tol1 = tol*abs(x)+ZEPS
         tol2 = 2.0*tol1
         
         if(abs(x-xm).le.(tol2-.5*(b-a))) exit
       
         if(abs(e).gt.tol1) then
            r=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=2.*(q-r)
            if(q.gt.0.) p=-p
            q=abs(q)
            etemp=e
            e=d
          
            if(abs(p).ge.abs(.5*q*etemp) .or. p.le.q*(a-x) .or. p.ge.q*(b-x)) goto 1
            d=p/q
            u=x+d
            !write(*,*) u,"u1"
            if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
            goto 2
         endif
         
1        if(x.ge.xm) then
            e=a-x
         else
            e=b-x
         endif
         d=CGOLD*e
         
         !write(*,*) abs(d),"abs d"
         
2        if(abs(d).ge.tol1) then
            u=x+d
            !write(*,*) u,"u2"
         else
            u=x+sign(tol1,d)
            !write(*,*) u,"u3"
         endif
         
         !write(*,*) u,"brent"
         fu=f(u)
         !write(*,*) fu,"fu"
         
         if(fu.le.fx) then
            if(u.ge.x) then
               a=x
            else
               b=x
            endif
            v=w
            fv=fw
            w=x
            fw=fx
            x=u
            fx=fu
            
         else
            if(u.lt.x) then
               a=u
            else
               b=u
            endif
            
            if(fu.le.fw .or. w.eq.x) then
               v=w
               fv=fw
               w=u
               fw=fu
            elseif(fu.le.fv .or. v.eq.x .or. v.eq.w) then
               v=u
               fv=fu
            endif
         endif
         
      enddo
      
      if (iter .eq. ITMAX) then
         write(*,*) 'brent exceed maximum iterations'
         stop
      endif
      
      xmin=x
      brent=fx
      write(*,*) xmin,brent,"xmin_brent"
      return
      
    end function brent
  
    real*8 function brent1(ax,bx,cx,f,tol,xmin)
      
      implicit none
      integer :: ITMAX
      real*8  :: ax,bx,cx
      real*8  :: tol,xmin
      real*8  :: CGOLD,ZEPS
      real*8  :: f
      
      PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.0e-10)
    
      integer :: iter
      real*8  :: a,b,d,e,etemp,p,q,r
      real*8  :: u,v,w,x,xm
      real*8  :: fu,fv,fw,fx
      real*8  :: tol1,tol2
      
      integer :: choice_flg

      
      a = min(ax,cx)
      b = max(ax,cx)
      
      v=bx
      w=v
      x=v
      e=0.d0
      
      fx=f(x)
      fv=fx
      fw=fx
      
      choice_flg = -1
      
      do iter=1,ITMAX
         
         xm   = 0.5*(a+b)
         tol1 = tol*abs(x)+ZEPS
         tol2 = 2.0*tol1
         
         if(abs(x-xm).le.(tol2-.5*(b-a))) exit
         
         if(abs(e) .gt. tol1) then
            
            r = (x-w)*(fx-fv)
            q = (x-v)*(fx-fw)
            p = (x-v)*q-(x-w)*r
            q = 2.0d0*(q-r)
            
            if (q .gt. 0.0d0) p = -p 
            
            q = abs(q)
            etemp=e
            e = d
            
            if(abs(p).ge.abs(.5*q*etemp) .or. p.le.q*(a-x) .or. p.ge.q*(b-x)) then
               choice_flg = 1
               if(x.ge.xm) then
                  e = a-x
               else
                  e = b-x
               endif
               
               d = CGOLD*e
            else
               choice_flg = 0
               d=p/q
               u=x+d
               !write(*,*) u,"u1"
               if((u-a).lt.tol2 .or. (b-u).lt.tol2) d=sign(tol1,xm-x)
               
            endif
            
         else
            if(x.ge.xm) then
               e=a-x
            else
               e=b-x
            endif
            d=CGOLD*e
            
         endif
         
         
         !write(*,*) abs(d),"abs d"   
         
         if(abs(d).ge.tol1) then
            u=x+d
            !write(*,*) u,"u2"
         else
            u=x+sign(tol1,d)
            !write(*,*) u,"u3"
         endif
         
         !write(*,*) u,"brent1"
         fu=f(u)
         
         if(fu.le.fx) then
          
            if(u.ge.x) then  
               a=x
            else
               b=x
            endif
            
            v=w
            fv=fw
            w=x
            fw=fx
            x=u
            fx=fu
          
         else
            if(u.lt.x) then
               a=u
            else
               b=u
            endif
            
            if(fu.le.fw .or. w.eq.x) then
               v=w
               fv=fw
               w=u
               fw=fu
            elseif (fu.le.fv .or. v.eq.x .or. v.eq.w) then
               v=u
               fv=fu
            endif
            
         endif
         
      enddo
      
      if (iter .eq. ITMAX) then
         write(*,*) 'brent exceed maximum iterations'
         stop
      endif
      
      xmin=x 
      brent1=fx
      write(*,*) xmin,brent1,"xmin_brent1"
      return
      
    end function brent1
    
    
    real*8 function PenF_Minmztn(Val)
      
      implicit none
      real*8  :: Val
      integer :: iter
      real*8  :: fret,fpp
      !real*8  :: PenF

      integer :: i
 
      
      fret = 0.0d0; fpp = 0.0d0
      
      if (spr_or_area .eq. 1) then
         do i = 1,hm
            l0(nmbrs(i)) = Val
         enddo
         !write(*,*) l0(nmbrs(1)),l0(nmbrs(2)),"los"
         
      elseif (spr_or_area .eq. 2) then
         do i = 1,hm
            A0(nmbrs(i)) = Val
         enddo
         
      endif
      
      call frprmn(coordntes,grd_mv,iter,fret,fpp,N_mvCoordnte)
      
      cnt = cnt + 1
      !write(*,*) cnt,"cnt in PenF"
      
      call coordntes_to_nodes(coordntes,node_xy)
      !write(*,*) l(31),l(32),l(34),l(35),"l 31-32-34-35"
      
      PenF_Minmztn = PenF(node_xy)
      !write(*,*) PenF_Minmztn,"PenF_Minmztn"
      
      
    end function PenF_Minmztn



    subroutine brak_with_itertn(ax,bx,cx,fa,fb,fc,func,fctr_val)
      implicit none
      real*8 :: ax,bx,cx,fa,fb,fc,func
      real*8 :: fctr_val
      

      do
         fa = func(ax)
         fb = func(bx)
         fc = func(cx)
         
         if (fb.lt.fa .and. fb.lt.fc) then
            exit
            
         elseif (fc.lt.fa .and. fc.lt.fb) then
            ax = cx
            bx = cx/fctr_val
            cx = cx/(fctr_val)**2
            
         elseif (fa.lt.fb .and. fa.lt.fc) then
            ax = ax
            bx = ax*fctr_val
            cx = ax*(fctr_val)**2
         endif
         
      enddo

      
    end subroutine brak_with_itertn


    subroutine PenF_vs_l0
      implicit none
      integer :: i
      real*8  :: l0_val,PenF_val
      integer :: lp_cnt
      
      lp_cnt = 0
      open(unit=15,file='PenF_vs_l0.dat')
      
      do i = -5,+5
         l0_val = 1.0d0 + (i*0.1d0)
         PenF_val = PenF_Minmztn(l0_val)
       
         write(unit=15,fmt=*) l0_val,PenF_val
         
         if (i.eq.(-5) .or. i.eq.(+5)) then
            !call save_config_and_generate_data(coordntes,5,2,lp_cnt)
            lp_cnt = lp_cnt + 1
         endif
      enddo
      
      close(15)
      
    end subroutine PenF_vs_l0

    

    
  end subroutine minimize_PenF_wrt_l0s


 
end module PenF_minimizer_mod
