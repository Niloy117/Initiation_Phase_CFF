
subroutine TimeStepAlg(Mc,grd_M,funcChoice,N_variabls,ftol)
    
  implicit none
  integer,intent(in) :: funcChoice,N_variabls
  real*8, intent(in) :: ftol
  
  real*8, intent(inout) :: Mc(1:N_variabls),grd_M(1:N_variabls)
  
  real*8  :: v0(1:N_variabls),vf(1:N_variabls)
  
  real*8  :: order_E,order_EPC,order_RK4
  real*8  :: tmn,tmx
  
  tmn=0.0d0
  tmx=15.0d0
  
  v0(1:N_variabls) = Mc(1:N_variabls)
  call dfunc(v0,grd_M,N_variabls,funcChoice)
    
  call Euler(v0,vf,grd_M,tmx,tmn,N_variabls,funcChoice,order_E)
  !call Euler_PC(v0,vf,grd_M,tmx,tmn,N_variabls,funcChoice,order_EPC)
  
  Mc(1:N_variabls) = vf(1:N_variabls)
  
end subroutine TimeStepAlg


subroutine Euler(v0,vf,grd_M,tmx,tmn,N_variabls,funcChoice,order)
  !use Energy_module
  implicit none
  
  real*8, intent(in)  :: v0(1:N_variabls)
  real*8, intent(out) :: vf(1:N_variabls)
  real*8, intent(in)  :: grd_M(1:N_variabls)
  real*8, intent(in)  :: tmx,tmn
  integer,intent(in)  :: N_variabls
  integer,intent(in)  :: funcChoice
  real*8, intent(out) :: order
    
  integer :: N,i,j
  integer :: stpIncr
  
  real*8  :: dt
  real*8  :: vi(1:N_variabls)
  real*8  :: ti,tf
  real*8  :: Res(1:3)
  real*8  :: ResV
  
  real*8  :: slope(1:N_variabls)
  real*8  :: fprv,fcur
  logical :: lgcl_cmpr
  real*8  :: time
  
  open(unit=16,file='TimeStepEnergy.dat',position='append')
  open(unit=17,file='TimeVsEnergy.dat',position='append')
  
  slope = grd_M
  
  
  stpIncr=1
    
  N=5*(2**(stpIncr-1))*1000
  dt=(tmx-tmn)/N
    
  vi=v0
    
  !do i = 1,N
  i=1
  do
     do j = 1,N_variabls
        vf(j) = vi(j) - slope(j)*dt
     enddo
     
     call dfunc(vf,slope,N_variabls,funcChoice)
     call Energy_comp(vi,vf,N_variabls,fprv,fcur,lgcl_cmpr)

     time = tmn+i*dt
     write(16,fmt=*) time,fcur,lgcl_cmpr
     
     if (lgcl_cmpr.eqv..True.) then
        write(*,*) i,"i"
        write(17,fmt=*) time,fcur
        exit
     endif
     
     tf=tmn+i*dt
       
     vi=vf
     i=i+1
  enddo
  
  ResV = vi(1)
  !Res(stpIncr)=vi(1)
  
  write(16,fmt=*) " "
  close(16)
  close(17)
  
end subroutine Euler
    
  
subroutine Euler_PC(v0,vf,grd_M,tmx,tmn,N_variabls,funcChoice,order)
  
  implicit none
  
  real*8, intent(in)  :: v0(1:N_variabls)
  real*8, intent(out) :: vf(1:N_variabls)
  real*8, intent(in)  :: grd_M(1:N_variabls)
  real*8, intent(in)  :: tmx,tmn
  integer,intent(in)  :: N_variabls
  integer,intent(in)  :: funcChoice
  real*8, intent(out) :: order
  
  integer :: N,i,j
  integer :: stpIncr
  
  real*8 :: dt
  real*8 :: vi(1:N_variabls),vp(1:N_variabls),vc(1:N_variabls)
  
  real*8 :: ti
  real*8 :: k1,k2,l1,l2
  real*8 :: Res(1:3),ResV
  
  real*8 :: slope1(1:N_variabls),slope2(1:N_variabls)
  
  slope1 = grd_M
  
  !do stepIncr=1,3
  stpIncr=1
  
  N=(2**(stpIncr-1))*1000
  dt=(tmx-tmn)/N
  vi=v0
  
  do i=1,N
     
     ti=tmn+i*dt
     
     do j = 1,N_variabls
        vp(j)=vi(j)-slope1(j)*dt
     enddo
     
     call dfunc(vp,slope2,N_variabls,funcChoice)
     
     do j = 1,N_variabls
        vc(j)=vi(j)-0.5d0*(slope1(j)+slope2(j))*dt 
     enddo
     
     vi=vc
     
  enddo
  
  ResV = vi(1)
  !Res(j)=si
  
  vf = vc
  
  
end subroutine Euler_PC


subroutine Energy_comp(coorPrv,coorCur,N_variabls1,fprv,fcur,lgcl_cmpr)
  use Energy_module
  use calls_for_tests
  
  implicit none
  real*8,  intent(in)  :: coorPrv(1:N_variabls1),coorCur(1:N_variabls1)
  integer, intent(in)  :: N_variabls1
  real*8,  intent(out) :: fprv,fcur
  logical, intent(out) :: lgcl_cmpr
  
  real*8 :: nodesPrv(1:N_node,1:N_dmnsn),nodesCur(1:N_node,1:N_dmnsn)
  real*8 :: EPS
  
  if (N_variabls1==N_mvCoordnte) then
     call coordntes_to_nodes(coorPrv,nodesPrv)
     call coordntes_to_nodes(coorCur,nodesCur)
     
  elseif (N_variabls1.ne.N_mvCoordnte) then
     write(*,*) "N_variabls not equal N_mvCoordnte"
     write(*,*) "flnm:Time_Stepping, sbrtn:Energy_comp"
     stop
  endif
     
  fprv = Energy(nodesPrv,l0,A0)
  fcur = Energy(nodesCur,l0,A0)
  
  EPS  = 1.0d-08
  
  lgcl_cmpr = .False.

  !write(*,*) ftol,"ftol"
  !stop
  
  if ((2.d0*abs(fcur-fprv)) .le. (ftol*(abs(fcur)+abs(fprv)+EPS))) then
     lgcl_cmpr=.True.
  endif
  
end subroutine Energy_comp

  
 
 
