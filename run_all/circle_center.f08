
subroutine get_circle_center(x1,y1,x2,y2,R,Cx,Cy)
  implicit none
  real*8, intent(in)  :: x1,y1,x2,y2
  real*8, intent(in)  :: R
  real*8, intent(out) :: Cx(1:2),Cy(1:2)

  real*8 :: k1,k2,l1,l2
  real*8 :: a,b,c
  real*8 :: f1,f2,g1,g2,c1,c2
  real*8 :: TOL
  
  TOL = 1.d-15
  
  if (abs(x1-x2).gt.TOL .and. abs(y1-y2).gt.TOL) then
     call get_k1_k2(x1,y1,x2,y2,k1,k2)
     write(*,*) k1,k2,"k1k2"
     call get_l1_l2(x1,y1,k1,k2,l1,l2)
     write(*,*) l1,l2,"l1l2"
     call get_abc(k1,k2,l1,l2,R,a,b,c)
     write(*,*) a,b,c,"abc"
     
     call get_root(a,b,c,f1,f2)
     !write(*,*) f1,f2,"f1f2"
     call get_g(f1,f2,k1,k2,g1,g2)
     call get_c(f1,f2,l1,l2,c1,c2)
     
     Cx(1) = -g1 ; Cy(1) = -f1
     Cx(2) = -g2 ; Cy(2) = -f2
     
  elseif (abs(x1-x2).lt.TOL .and. abs(y1-y2).gt.TOL) then
     call get_f_equalX(y1,y2,f1,f2)
     call get_l1_l2_equalX(x1,y1,x2,y2,l1,l2)
     call get_abc_equalX(l1,l2,y1,y2,R,a,b,c)
     
     call get_root(a,b,c,g1,g2)
     call get_c_equalX(g1,g2,l1,l2,c1,c2)

     Cx(1) = -g1 ; Cy(1) = -f1
     Cx(2) = -g2 ; Cy(2) = -f2
     
  elseif (abs(y1-y2).lt.TOL .and. abs(x1-x2).gt.TOL) then
     call get_g_equalY(x1,x2,g1,g2)
     call get_l1_l2_equalY(x1,y1,x2,y2,l1,l2)
     call get_abc_equalY(l1,l2,x1,x2,R,a,b,c)
     
     call get_root(a,b,c,f1,f2)
     call get_c_equalY(f1,f2,l1,l2,c1,c2)

     Cx(1) = -g1 ; Cy(1) = -f1
     Cx(2) = -g2 ; Cy(2) = -f2
  else
     write(*,*) "both node are the same,can't happen"
     write(*,*) x1,y1,x2,y2,"x1y1x2y2"
     stop
  endif
  
end subroutine get_circle_center

subroutine get_k1_k2(x1,y1,x2,y2,k1,k2)
  implicit none
  real*8, intent(in)  :: x1,y1,x2,y2
  real*8, intent(out) :: k1,k2
  
  real*8 :: k2Num,k2dNum
  
  k1 = -(y1-y2)/(x1-x2)
  k2Num  = ((x2)**2-(x1)**2) + ((y2)**2-(y1)**2)
  k2dNum = 2.0*(x1-x2)
 
  k2 = k2Num/k2dNum
  
end subroutine get_k1_k2

subroutine get_l1_l2(x1,y1,k1,k2,l1,l2)
  implicit none
  real*8, intent(in)  :: x1,y1,k1,k2
  real*8, intent(out) :: l1,l2

  l1 = -(2.0d0*k1*x1+2.0d0*y1)
  l2 = -(x1)**2-(y1)**2-(2.0d0*k2*x1)
  
end subroutine get_l1_l2

subroutine get_l1_l2_equalX(x1,y1,x2,y2,l1,l2)
  implicit none
  real*8, intent(in)  :: x1,y1,x2,y2
  real*8, intent(out) :: l1,l2

  l1 = -(2.0d0*x1)
  l2 = -x1**2-y1**2+y1*(y1+y2)
  
end subroutine get_l1_l2_equalX

subroutine get_l1_l2_equalY(x1,y1,x2,y2,l1,l2)
  implicit none
  real*8, intent(in)  :: x1,y1,x2,y2
  real*8, intent(out) :: l1,l2

  l1 = -(2.0d0*y1)
  l2 =  -x1**2-y1**2+x1*(x1+x2)
  
end subroutine get_l1_l2_equalY

subroutine get_abc(k1,k2,l1,l2,R,a,b,c)
  implicit none
  real*8, intent(in)  :: k1,k2,l1,l2
  real*8, intent(in)  :: R
  real*8, intent(out) :: a,b,c
  
  a = 1.0d0+(k1)**2
  b = 2.0d0*k1*k2-l1
  c = (k2)**2-(R)**2-l2
  
end subroutine get_abc

subroutine get_abc_equalX(l1,l2,y1,y2,R,a,b,c)
  implicit none
  real*8, intent(in)  :: l1,l2,y1,y2,R
  real*8, intent(out) :: a,b,c
  
  a = 1.0d0
  b = -l1
  c = ((y1+y2)/2.0d0)**2-l2-R**2
  
end subroutine get_abc_equalX

subroutine get_abc_equalY(l1,l2,x1,x2,R,a,b,c)
  implicit none
  real*8, intent(in)  :: l1,l2,x1,x2,R
  real*8, intent(out) :: a,b,c
  
  a = 1.0d0
  b = -l1
  c = ((x1+x2)/2.0d0)**2-l2-R**2
  
end subroutine get_abc_equalY

subroutine get_root(a,b,c,rt1,rt2)
  implicit none
  real*8, intent(in)  :: a,b,c
  real*8, intent(out) :: rt1,rt2
  real*8 :: det,det_sq

  det_sq = b**2-4.0d0*a*c
  write(*,*) a,b,c,det_sq,"det_sq"
  det = sqrt(b**2-4.0d0*a*c) 
  rt1  = (-b+det)/(2.0d0*a)
  rt2  = (-b-det)/(2.0d0*a)
  
end subroutine get_root

subroutine get_f_equalX(y1,y2,f1,f2)
  implicit none
  real*8, intent(in)  :: y1,y2 
  real*8, intent(out) :: f1,f2
  
  f1 = -(y1+y2)/2.0d0
  f2 = f1
  
end subroutine get_f_equalX
  
subroutine get_g(f1,f2,k1,k2,g1,g2)
  implicit none
  real*8, intent(in)  :: f1,f2,k1,k2
  real*8, intent(out) :: g1,g2
  
  g1 = k1*f1+k2
  g2 = k1*f2+k2
  
end subroutine get_g

subroutine get_g_equalY(x1,x2,g1,g2)
  implicit none
  real*8, intent(in)  :: x1,x2
  real*8, intent(out) :: g1,g2

  g1 = -(x1+x2)/2.0d0
  g2 = g1
  
end subroutine get_g_equalY

subroutine get_c(f1,f2,l1,l2,c1,c2)
  implicit none
  real*8, intent(in)  :: f1,f2,l1,l2 
  real*8, intent(out) :: c1,c2

  c1 = l1*f1+l2
  c2 = l1*f2+l2

end subroutine get_c

subroutine get_c_equalX(g1,g2,l1,l2,c1,c2)
  implicit none
  real*8, intent(in)  :: g1,g2,l1,l2
  real*8, intent(out) :: c1,c2
  
  c1 = l1*g1 + l2
  c2 = l1*g2 + l2
  
end subroutine get_c_equalX

subroutine get_c_equalY(f1,f2,l1,l2,c1,c2)
  implicit none
  real*8, intent(in)  :: f1,f2,l1,l2
  real*8, intent(out) :: c1,c2

  c1 = l1*f1 + l2
  c2 = l1*f2 + l2

end subroutine get_c_equalY

subroutine get_allowable_center(cellNo1,cellNo2,P1,P2,Cx,Cy,xc,yc)
  implicit none
  integer,intent(in)  :: cellNo1,cellNo2
  real*8 ,intent(in)  :: P1,P2
  real*8, intent(in)  :: Cx(1:2),Cy(1:2)
  real*8, intent(out) :: xc,yc

  real*8 :: xdyd(1:2)
  real*8 :: d1,d2
  real*8 :: Center1(1:2),Center2(1:2)  

  interface
     function diag_meeting_point(cell)
       implicit none
       real*8  :: diag_meeting_point(1:2)
       integer :: cell
     end function diag_meeting_point
     
     real*8 function distance(A,B,N_dim)
       implicit none
       real*8  :: A(1:2),B(1:2)
       integer :: N_dim
     end function distance
  end interface
  
  if (P1.ge.P2) then
     xdyd(1:2) = diag_meeting_point(cellNo1)
  elseif (P1.lt.P2) then
     xdyd(1:2) = diag_meeting_point(cellNo2)
  else
     write(*,*) "P1=P2"
  endif
  
  Center1(1) = Cx(1) ; Center1(2) = Cy(1)
  Center2(1) = Cx(2) ; Center2(2) = Cy(2)
  
  d1 = distance(xdyd,Center1,2) !2 for N_dmnsn
  d2 = distance(xdyd,Center2,2)
  
  if(d1.lt.d2) then
     xc = Cx(1) ; yc = Cy(1)
  elseif (d2.lt.d1) then
     xc = Cx(2) ; yc = Cy(2)
  else
     write(*,*) "d1=d2"
  endif
  
end subroutine get_allowable_center

subroutine get_slope_and_intercept(x1,y1,x2,y2,slope,IC)
  implicit none
  real*8, intent(in)  :: x1,y1,x2,y2
  real*8, intent(out) :: slope,IC
  real*8 :: TOL

  TOL = 1.d-16
  
  if(abs(x1-x2).lt.TOL) then
     slope = 10000.d0
     IC    = -10000.d0
  else
     slope = (y2-y1)/(x2-x1)
     IC    = (x1*y2-x2*y1)/(x1-x2)
  endif
  
end subroutine get_slope_and_intercept

real*8 function strght_line(xd,yd,slope,IC)
  implicit none
  real*8, intent(in) :: xd,yd
  real*8, intent(in) :: slope,IC

  strght_line = yd-(slope*xd)-IC

end function strght_line

function diag_meeting_point(cellNo)
  use system_parameters
  use transfrm_info
  
  implicit none
  real*8  :: diag_meeting_point(1:2)
  integer :: cellNo
  integer :: totlNodeCnnctd
  
  integer :: D1,D2,D3,D4
  real*8  :: c1,c2,m1,m2
  real*8  :: x1,y1,x2,y2
  
  totlNodeCnnctd = area_node(cellNo,0)
  
  if (totlNodeCnnctd==4) then
     D1 = area_node(cellNo,1)
     D2 = area_node(cellNo,3)
     D3 = area_node(cellNo,2)
     D4 = area_node(cellNo,4)
     
     x1=node_xy(D1,1) ; y1=node_xy(D1,2)
     x2=node_xy(D2,1) ; y2=node_xy(D2,2)
     call get_slope_and_intercept(x1,y1,x2,y2,m1,c1)
     
     x1=node_xy(D3,1) ; y1=node_xy(D3,2)
     x2=node_xy(D4,1) ; y2=node_xy(D4,2)
     call get_slope_and_intercept(x1,y1,x2,y2,m2,c2)
     
     diag_meeting_point(1) = (c2-c1)/(m1-m2)
     diag_meeting_point(2) = ((m1*(c2-c1))/(m1-m2)) + c1
  else
     continue
  endif
  
end function diag_meeting_point

subroutine get_CurveStrtEndAngl(x1,y1,x2,y2,xc,yc,Angl1,Angl2)
  implicit none
  real*8, intent(in)  :: x1,y1,x2,y2,xc,yc
  real*8, intent(out) :: Angl1,Angl2

  real*8 :: xtrns,xs,ys
  real*8 :: diffangl
  real*8 :: dotP,dotab_by_ab
  real*8 :: mltpAB,vectA_val,vectB_val

  real*8  :: vectA(1:2),vectB(1:2)
  integer :: quardrnt
  real*8  :: Pi
  real*8  :: AnglStr
  
  interface
     real*8 function dot_prdct(vectA,vectB,dmnsn)
       implicit none
       integer, intent(in) :: dmnsn
       real*8,  intent(in) :: vectA(1:dmnsn),vectB(1:dmnsn)
     end function dot_prdct
     real*8 function Mgnitude(vectA,dmnsn)
       implicit none
       integer,intent(in) :: dmnsn
       real*8 ,intent(in) :: vectA(1:dmnsn)
     end function Mgnitude     
  end interface
  
  open(unit=198,file='CurveStrtEndAngl.dat',position='append')
  
  Pi = atan(1.0d0)*4.0d0
  AnglStr = 0.0d0
  
  xtrns = 0.1d0
  xs = xc + xtrns
  ys = yc

  vectA(1) = (xc-xs) ; vectA(2) = (yc-ys)
  vectB(1) = (xc-x1) ; vectB(2) = (yc-y1)
  
  write(198,fmt=*) vectA(1:2),vectB(1:2),"vectAB"
  write(198,fmt=*) xc,yc,"xcyc"
  write(198,fmt=*) x1,y1,x2,y2,"xy12"
  
  dotP      = dot_prdct(vectA,vectB,2)
  vectA_val = Mgnitude(vectA,2)
  vectB_val = Mgnitude(vectB,2)
  mltpAB    = vectA_val*vectB_val
  
  dotAB_by_AB = dotP/(mltpAB)
  Angl1 = acos(dotAB_by_AB)*(180.0d0/Pi)
  
  call get_quardrnt(xc,yc,x1,y1,quardrnt)

  if (quardrnt==1.or.quardrnt==2) then
     continue
  elseif (quardrnt==3) then
     diffAngl = 180.0d0-Angl1
     Angl1 = 180.0d0+diffAngl
  elseif (quardrnt==4) then
     diffAngl = Angl1
     Angl1 = 360.0d0-diffAngl
  endif

  vectB(1)  = (xc-x2) ; vectB(2) = (yc-y2)
  dotP      = dot_prdct(vectA,vectB,2)
  vectB_val = Mgnitude(vectB,2)
  mltpAB    = vectA_val*vectB_val
  
  dotAB_by_AB = dotP/(mltpAB)
  
  Angl2 = acos(dotAB_by_AB)*(180.0d0/Pi)
  call get_quardrnt(xc,yc,x2,y2,quardrnt)

  if (quardrnt==1.or.quardrnt==2) then
     continue
  elseif (quardrnt==3) then
     diffAngl = 180.0d0-Angl2
     Angl2 = 180.0d0+diffAngl
  elseif (quardrnt==4) then
     diffAngl = Angl2
     Angl2 = 360.0d0-diffAngl
  endif

  if (Angl2.gt.Angl1) then
     continue
  elseif (Angl2.le.Angl1) then
     AnglStr = Angl2
     Angl2   = Angl1
     Angl1   = Anglstr
  endif

  if (abs(Angl1-Angl2).gt.180.0d0) then
     AnglStr = Angl2
     Angl2   = Angl1
     Angl1   = Anglstr
  endif
  
  close(198)
  
end subroutine get_CurveStrtEndAngl

subroutine get_quardrnt(xc,yc,xa,ya,quardrnt)
  implicit none
  real*8, intent(in)  :: xc,yc,xa,ya
  integer,intent(out) :: quardrnt

  !real*8 :: orginX,orginY
  real*8 :: xaU,yaU
  
  !orginX = 0.0d0 ; orginY = 0.0d0
  xaU = xa-xc    ; yaU = ya-yc

  if (xaU.ge.0.0d0 .and. yaU.ge.0) then
     quardrnt = 1
  elseif (xaU.lt.0.0d0 .and. yaU.ge.0) then
     quardrnt = 2
  elseif (xaU.lt.0.0d0 .and. yaU.lt.0) then
     quardrnt = 3
  elseif (xaU.ge.0.0d0 .and. yaU.lt.0) then
     quardrnt = 4
  endif
  
end subroutine get_quardrnt



