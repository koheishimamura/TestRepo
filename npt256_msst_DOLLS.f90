      module common_variables
!---npart=128→250---------------
!npart=128:rcoff=1.d-9,side=4.d0*a,hmax=11,bcc=3,kkk=200,Pe=1.d0
!npart=250:rcoff=1.3d-9,side=5.d0*a,hmax=12,bcc=4,kkk=250,Pe=0.6d0,1000step40min
!--------------------------------
integer, parameter :: npart=256          !Argon粒子数 --
integer, parameter :: nt=3*npart
integer, parameter :: nnos=1
integer, parameter :: kkk=256
integer,parameter::Termo=10   !--熱浴の数設定--
character(5) :: CHECK
        real*8 a0,e0,m0,t0,boltz,pi,e
        real*8 a,tref,Pe,dt,rcoff,M,Q,mar,rar
        real*8 side,sideh,rcoffs,tscale,V,pv,ps,s,P
        real*8 s7ii,s7ia,s7aa,alpha,w,qe2,eps,side2
        real*8 ekin,epot,vir,rds,r148,rd,hx,hy,hz
        real*8 x(npart),y(npart),z(npart),x0(npart),y0(npart),z0(npart)
        real*8 vx(npart),vy(npart),vz(npart)
        real*8 fx(npart),fy(npart),fz(npart),qe(nt)
        real*8 vn,vn1,vn2,vn3
real*8 cstress(3,3),cstress0(3,3),cstress1(3,3),cstress2(3,3),cstress3(3,3),alp
real*8 qi(3,npart),qvi(3,npart),boxinv(3,3),abcreal3,creal3(3,3)
real*8 ek,etot2,etott,temp,pres,hen,etot,abcreal,virp
real*8 veigv(3,3),veig(3),EPSI,W1R(3),W1I(3),work(3)
real*8 nf,gboxg(3,3),vboxg(3,3),box(3,3),ekinb,ekins,hen1,hen2,hen3
real*8 glogs(nnos),vlogs(nnos),xlogs(nnos),Qq(nnos)
real*8 vgtvg(3,3),pext,pint(3,3),B,Qq2
real*8 dr,gagag,gagi,gii,ra,VV,absxii(nt),absyii(nt),abszii(nt)
real*8 nagag(kkk),nagi(kkk),nii(kkk),qi0(3,npart),qia(3,nt),sp2(3),trefc
integer NN,MTRX,IFVEC,IFSRT

        integer arg,vmd,avg,gn,conti
integer timemx,clock,hmax,chi,timea,clocks
integer chtemp,clockss

!--11/12日追加フーバー形式M≧2--
integer hoover
real*8 Qw(Termo),psw(Termo),sw(Termo),G

real*8 V0,tmass,Vs

real*8 Vfin,teps
integer Maxtime

      end module

!---------------------------------------------------------------------
      PROGRAM Argon
!---------------------------------------------------------------------

      use common_variables

      !-- Hartree単位換算の基準値（現実系） --

      a0=5.29177d-11     !-- 長さ(m) --
      m0=9.10939d-31     !-- 質量(kg) --
      e0=4.35975d-18     !-- エネルギー(J) --
      t0=2.41888d-17     !-- 時間(s) --
      boltz=1.38065d-23  !-- ボルツマン定数(J/K) --

      pi=acos(-1.d0)     !-- 円周率 --
      e=1.60218d-19      !-- 電気素量 --

      !-- 初期値の設定（現実系） --

      a=5.739925d-10        !-- Argonfccの格子定数(資料より) --
      tref=50.d0        !-- 温度(K) --
trefc=tref
      Pe=0.d0            !-- 圧力(GPa) --
      dt= 2.0d-15        !-- タイムステップ --
      timemx=10000         !-- ステップ数 --
      rcoff=1.0d-9       !-- ポテンシャルカット（セルの一辺の長さの半分より少し小さい:適当) --
      mar=6.69d-26       !--Argon質量(資料より)--
      rar=1.91d-10       !--Argon半径(使ってない)--
      Vs = 000.d0      !shock speed m/s

      Maxtime = 10000 
      Vfin    = 0.d0

      !-- アルゴリズムの決定 --

      arg=10  !-- NVE=0,NVT=1,NPH=2,NPT=3,MSST=9,DOLLS=10 --
      vmd=1  !-- VMD用ファイルを作成するか(yes=1) --
      avg=0  !-- 3D AVS Player用ファイルを作成するか(yes=1) --
      gn=0   !-- Gnuplot用ファイルを作成するか(yes=1) --


      !-- 単位の換算（現実系→仮想系） --
      a=a/a0
      side=4.d0*a        !-- 基本セルの一辺の長さ １セル中の粒子数を256個にするために必要--
      sideh=side/2.d0
      side2=side**2.d0
      rcoff=rcoff/a0
      rcoffs=rcoff*rcoff
      dt=dt/t0
      mar=mar/m0
      rar=rar/a0
      tscale=1.d0/(3.d0*npart)
      tmass = mar*dble(npart)
      M=1.d9
      Q=300000.d0
      Qw=(/300000.0,300000.0,300000.0,300000.0,300000.0,300000.0,300000.0,300000.0,300000.0,300000.0/)  !--熱浴の質量の値--      
      tref=tref*boltz/e0
      Pe=Pe*(1.d9*a0**3)/e0
      V=side*side*side   !-- 基本セルの体積 --
      V0 = V
      Vs = Vs*t0/a0
      Vfin = Vfin/(a0**3)
      Vfin = V0*2.d0
      teps = dlog(Vfin/V0)/dble(Maxtime)/dt
      write(*,*)"teps=",teps
      !stop
      pv=0.d0
      hoover=2  !--フーバー形式でM=1(=1)かM≧2(=2)を使うかを決める--
      s=1.d0  !--s=0.0からはじめるべきでは。フーバー形式M=1。熱浴の座標初期値。--
      sw=1.0  !--フーバー形式M≧2。各熱浴の座標初期値--
      ps=0.d0 !--フーバー形式M=1。熱浴の運動量初期値--
      psw=0.0  !--フーバー形式M≧2。各熱浴の運動量初期値--
      clock=0            !-- ステップ数のカウンター --
clocks=0
clockss=0
chi=4.d0
hmax=12
alp=chi/rcoff
nf=3.d0*nt
cstress=0.0
cstress1=0.0
cstress2=0.0
cstress3=0.0

!---SUBROUTINE RsEIGQR---
NN=3
MTRX=3
EPSI=1.d-30
IFVEC=1
IFSRT=2
CHECK='check'

gboxg=0.d0
vboxg=0.d0
vboxg(1,1) = -0.0000
ekinb=0.d0
glogs=0.d0
vlogs=0.d0
xlogs=0.d0

!--セルの形を決定今回は立方体なので以下のような設定--
box=0.0
box(1,1)=side
box(2,2)=side
box(3,3)=side

qia=0.d0

!---subroutine cf---
nagag=0
nagi=0
nii=0
gagag=0.d0
gagi=0.d0
gii=0.d0
ra=0.d0
VV=0.d0
dr=0.05d-10
dr=dr/a0
timea=1000


      !-- 力の計算の初期値(Hartree) --

      s7ii=(ri+ri)**7
      s7ia=(rag+ri)**9
      s7aa=(rag+rag)**11
      eps=5.42186d-3
      alpha=44.d0
      w=167.044d0
      do i=1,nt
        if (i>npart) then
          qe(i)=-0.6d0
        else
          qe(i)=0.6d0
        end if
      end do
      qe2=0.6d0*0.6d0

write(*,*)'npart=',npart
write(*,*)'rcoff=',rcoff
write(*,*)'side=',side
write(*,*)'hmax=',hmax
write(*,*)'Pe=',Pe

!write(*,*)'continue? 0=no 1=yes'
!read(*,*)conti
!write(*,*)'timemx='
!read(*,*)timemx
conti=0
!timemx=200

      !-- 演算開始 --
open(52,file='bep.dat',status='unknown')  !データを書き込むファイル
open(90,file='cf.dat',status='unknown')
open(91,file='msd.dat',status='unknown')
open(92,file='analyze_box.dat',status='unknown')
open(50,file='parameter.dat',status='unknown') !データを書き込むファイル
!open(150,file='absp.xyz',status='unknown')
if(vmd == 1) then
 open(51,file='Argon.xyz',status='unknown')
end if

pext=Pe
Qq(1)=300000.0
!B=100.d0*(1.d9*a0**3)/e0/1.4d0
!Qq2=tref/4.d0/pi/pi*(100.d0*dt)**2
!Qq=Qq2
!Qq(1)=Qq2*3.d0*nt
!M=3.d0*B*V*(100.d0*dt)**2/4.d0/pi/pi
!Qq(2)=10.d30
!Qq(1)=10.d30
!M=10.d50
!--q=L^-1rの計算をするために必要な逆行列L^-1を作るサブルチ。--
call inversebox

if(conti==0)then
 call fcc                         !初期配置の設定(今回は面心立方格子の配列を初期配置とする)--
 call mxwell                      !初期速度の設定(マクスウェル分布での初期速度の設定)--
 call force                       !初期力の設定(レナードジョーンズポテンシャルでの力の計算)--
end if

if(conti==1)then
 call continue
 clock=clocks
end if

if(arg==1) then
   write(6,'(1a6,5a15)')'step','conserve','temp','pres','etot'
end if
call result                      !初期状況の書き込み

open(80,file='con.dat',status='unknown')

!--以下より数値積分の開始--
 do clock=clocks+1,clocks+timemx
  if(clock<0 .and. mod(clock,10)==0)then
   call scaling     
  end if
  if(arg.eq.0) call NVE          !選んだアルゴリズムを開始する
  if(arg.eq.1) call NVT
  if(arg.eq.2) call NPH
  if(arg.eq.3) call NPT
  if(arg.eq.9) call MSST
  if(arg.eq.10) call DOLLS 
  call result                    !結果の書き込み
end do
! call cfb
! call absposition

write(80,*)x,y,z,x0,y0,z0,vx,vy,vz,box,vboxg,xlogs,vlogs,sp2,clock-1,qia
      stop
      end
!---------------------------------------------------------------------
subroutine scaling

use common_variables

write(6,*) 'velocity adjustment'
ts=tscale*ekin
write(6,*) 'temperature before scaling is ',ts
sc=tref/ts
sc=sqrt(sc)
write(6,*) 'scale factor is ',sc
do i=1,npart
 vx(i)=vx(i)*sc
 vy(i)=vy(i)*sc
 vz(i)=vz(i)*sc
end do

call rouekin

vboxg=0.d0
ekinb=0.d0
vlogs=0.d0
xlogs=0.d0

return
end
!---------------------------------------------------------------------
subroutine continue

use common_variables

open(100,file='con0.dat')
read(100,*)x,y,z,x0,y0,z0,vx,vy,vz,box,vboxg,xlogs,vlogs,sp2,clocks,qia
close(100)

call inversebox

do i=1,npart
 qi0(1,i)=boxinv(1,1)*x(i)+boxinv(1,2)*y(i)+boxinv(1,3)*z(i)
 qi0(2,i)=boxinv(2,1)*x(i)+boxinv(2,2)*y(i)+boxinv(2,3)*z(i)
 qi0(3,i)=boxinv(3,1)*x(i)+boxinv(3,2)*y(i)+boxinv(3,3)*z(i)
end do

qi=qi0

call force

call rouekin

!---get the box kinetic energy---
do i=1,3
 do j=1,3
  vgtvg(i,j)=vboxg(1,i)*vboxg(1,j)+vboxg(2,i)*vboxg(2,j)+vboxg(3,i)*vboxg(3,j)
 end do
end do
ekinb=(vgtvg(1,1)+vgtvg(2,2)+vgtvg(3,3))*M

return
end
!---------------------------------------------------------------------
subroutine inversebox

use common_variables

boxinv(1,1)=box(2,2)*box(3,3)-box(3,2)*box(2,3)
boxinv(1,2)=box(2,1)*box(3,3)-box(2,3)*box(3,1)
boxinv(1,3)=box(2,1)*box(3,2)-box(2,2)*box(3,1)
boxinv(2,1)=box(1,2)*box(3,3)-box(1,3)*box(3,2)
boxinv(2,2)=box(1,1)*box(3,3)-box(1,3)*box(3,1)
boxinv(2,3)=box(1,1)*box(3,2)-box(1,2)*box(3,1)
boxinv(3,1)=box(1,2)*box(2,3)-box(1,3)*box(2,2)
boxinv(3,2)=box(1,1)*box(2,3)-box(3,1)*box(2,1)
boxinv(3,3)=box(1,1)*box(2,2)-box(1,2)*box(2,1)
V=boxinv(1,1)*box(1,1)-boxinv(1,2)*box(1,2)+boxinv(1,3)*box(1,3) !--Vはセルの体積でもあるが、行列Lの行列式でもある--
boxinv=boxinv/V                                                  !--V=a・(b×c)=det|a,b,c|より--

return
end

!---------------------------------------------------------------------
subroutine rouekin !--ビリアル（≡運動エネルギー×２）の計算（正しい計算はtextのように運動エネルギーで行なうが、1/2を入れるのが面倒なので最初から省いておく）--
use common_variables

ekin=0.d0
do i=1,npart
 ekin=ekin+mar*(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
end do


call cvn0

cstress=cstress0+cstress1+cstress2+cstress3
pint=cstress

return
end

!--------------------------------------------------------------------

      subroutine bcc

      use common_variables

      integer n,ph

      n=0

      do ph=0,1
      do i=0,4
      do j=0,4
      do k=0,4
        n=n+1
        x(n)=a*i+0.5d0*a*ph
        y(n)=a*j+0.5d0*a*ph
        z(n)=a*k+0.5d0*a*ph
      end do
      end do
      end do
      end do


      do ph=0,1
      do i=0,4
      do j=0,4
      do k=0,4
        n=n+1
        x(n)=0.5d0*a+a*i-0.5d0*a*ph
        y(n)=0.25d0*a+a*j+0.5d0*a*ph
        z(n)=a*k+0.5d0*a*ph
      end do
      end do
      end do
      end do

x0=x
y0=y
z0=z

do i=1,npart
 qi0(1,i)=boxinv(1,1)*x0(i)+boxinv(1,2)*y0(i)+boxinv(1,3)*z0(i)
 qi0(2,i)=boxinv(2,1)*x0(i)+boxinv(2,2)*y0(i)+boxinv(2,3)*z0(i)
 qi0(3,i)=boxinv(3,1)*x0(i)+boxinv(3,2)*y0(i)+boxinv(3,3)*z0(i)
end do

qi=qi0

sp2=0.d0
do i=1,npart
 sp2(1)=sp2(1)+mag*x(i)+mi*x(i+npart)
 sp2(2)=sp2(2)+mag*y(i)+mi*y(i+npart)
 sp2(3)=sp2(3)+mag*z(i)+mi*z(i+npart)
end do
sp2(1)=sp2(1)/(npart*(mag+mi))
sp2(2)=sp2(2)/(npart*(mag+mi))
sp2(3)=sp2(3)/(npart*(mag+mi))

      return
      end
!---------------------------------------------------------------------
subroutine fcc   !--二つに分ける。最初128個を配置する。次に最初に配置した粒子の間に128個配置する。--

use common_variables

integer n,lg
n=0

do lg=0,1
do i=0,3
do j=0,3
do k=0,3

   n=n+1
   x(n)=i*a+lg*a*0.5
   y(n)=j*a+lg*a*0.5
   z(n)=k*a
end do
end do
end do
end do

do lg=1,2
do i=0,3
do j=0,3
do k=0,3
  n=n+1
  x(n)=i*a+(2-lg)*a*0.5
  y(n)=j*a+(lg-1)*a*0.5
  z(n)=k*a+a*0.5
end do
end do
end do
end do


x0=x
y0=y
z0=z
!--先ほど計算した逆行列を使ってq=L^-1rの計算を行なっている。立方体の場合対角項のみの計算で、残りの非対角項は０になる。--
do i=1,npart
 qi0(1,i)=boxinv(1,1)*x0(i)+boxinv(1,2)*y0(i)+boxinv(1,3)*z0(i)
 qi0(2,i)=boxinv(2,1)*x0(i)+boxinv(2,2)*y0(i)+boxinv(2,3)*z0(i)
 qi0(3,i)=boxinv(3,1)*x0(i)+boxinv(3,2)*y0(i)+boxinv(3,3)*z0(i)
end do

qi=qi0

return
end
!---------------------------------------------------------------------
      subroutine mxwell !--資料参考--

      use common_variables

      real u1,u2,v1,v2,sa,r,sc
      real*8 vh(1:nt),arand
      integer ir

      ir=1

      do i=1,nt,2
 1      u1=arand(ir)
        u2=arand(ir)
        v1=2.d0*u1-1.d0
        v2=2.d0*u2-1.d0
        sa=v1*v1+v2*v2
        if (sa.ge.1.0) goto 1
        r=-2.d0*alog(sa)/sa
        vh(i)=v1*sqrt(r)
        vh(i+1)=v2*sqrt(r)
      end do

      do i=1,npart
         vx(i)=vh(i)
         vy(i)=vh(i+npart)
         vz(i)=vh(i+npart*2)
      end do

      ekin=0.d0

      sp=0.d0
      do i=1,npart
        sp=sp+mar*vx(i)
      end do
      write(6,*) 'total linear momentum in x direction is ',sp
      sp=sp/(npart*mar)
      do i=1,npart
        vx(i)=vx(i)-sp
        ekin=ekin+mar*vx(i)*vx(i)
      end do

      sp=0.d0
      do i=1,npart
        sp=sp+mar*vy(i)
      end do
      write(6,*) 'total linear momentum in y direction is ',sp
      sp=sp/(npart*mar)
      do i=1,npart
        vy(i)=vy(i)-sp
        ekin=ekin+mar*vy(i)*vy(i)
      end do

      sp=0.d0
      do i=1,npart
        sp=sp+mar*vz(i)
      end do
      write(6,*) 'total linear momentum in z direction is ',sp
      sp=sp/(npart*mar)
      do i=1,npart
        vz(i)=vz(i)-sp
        ekin=ekin+mar*vz(i)*vz(i)
      end do

      write(6,*) 'velocity adjustment'
      ts=tscale*ekin
      write(6,*) 'temperature before scaling is ',ts
      sc=tref/ts
      sc=sqrt(sc)
      write(6,*) 'scale factor is ',sc
      do i=1,npart
        vx(i)=vx(i)*sc
        vy(i)=vy(i)*sc
        vz(i)=vz(i)*sc
      end do

call rouekin

return
      end

!=====================================================================
real*8 function arand(ir)!--資料参考--

      implicit integer*4 (i-n)
real*8 anorm
      parameter (ml=1664501)
      parameter (ia=1229,ic=351750)
      parameter (anorm=1.d0/ml)

      ir=mod(ia*ir+ic,ml)
      arand=ir*anorm

      return
      end

!--------------------------------------------------------------------
      subroutine force

      use common_variables

      real*8 kx,ky,kz
real*8 xi,yi,zi,xx,yy,zz,xxx,yyy,zzz
real*8 sig,ip,frc,Crc,rcdpdr

sig=3.4d-10 !--レナードポテンシャルの定数について。Argonの場合のものを資料より--
ip=1.65d-21

sig=sig/a0 !--ハートレー換算で無次元化--
ip=ip/e0

vn=0.d0
call cvn0

do i=1,npart
 fx(i)=0.d0
 fy(i)=0.d0
 fz(i)=0.d0
end do
vir=0.d0
epot=0.d0
rcdpdr=-48.0*ip*(sig**12.0*rcoffs**(-7.0)-0.5*sig**6.0*rcoffs**(-4.0))*rcoff
Crc=-(4.0*ip*(sig**12.0*rcoffs**(-6.0)-sig**6.0*rcoffs**(-3.0)))

cstress3=0.d0
do i=1,npart
 xi=qi(1,i)
 yi=qi(2,i)
 zi=qi(3,i)
do j=i+1,npart
 xx=xi-qi(1,j)
 yy=yi-qi(2,j)
 zz=zi-qi(3,j)

        if (xx.lt.-0.5d0) xx=xx+1.d0
        if (xx.gt. 0.5d0) xx=xx-1.d0

        if (yy.lt.-0.5d0) yy=yy+1.d0
        if (yy.gt. 0.5d0) yy=yy-1.d0

        if (zz.lt.-0.5d0) zz=zz+1.d0
        if (zz.gt. 0.5d0) zz=zz-1.d0

xxx=box(1,1)*xx+box(1,2)*yy+box(1,3)*zz
yyy=box(2,1)*xx+box(2,2)*yy+box(2,3)*zz
zzz=box(3,1)*xx+box(3,2)*yy+box(3,3)*zz
xx=xxx
yy=yyy
zz=zzz
        rd=xx*xx+yy*yy+zz*zz
       
        if (rd.gt.rcoffs) goto 270 !--カットオフ以上はスルー--
!--ポテンシャルの計算（レナードジョーンズの式より）--
rds=sqrt(rd)
frc=-(rds-rcoff)*rcdpdr
epot=epot+4.0*ip*(sig**12.0*rd**(-6.0)-sig**6.0*rd**(-3.0))+frc+Crc
!--レナードポテンシャルの微分式より力を算出--
r148=48.0*ip*(sig**12.0*rd**(-7.0)-0.5*sig**6.0*rd**(-4.0))+rcdpdr/rds
vir=vir-rd*r148

kx=xx*r148 !--実際はkx,ky,kzの計算で力が算出できる。ただし一回目の呼び出しで算出しているのはt-⊿tの時の力であることに注意。二回目の呼び出しで⊿t時の計算。--
fx(i)=fx(i)+kx
fx(j)=fx(j)-kx
ky=yy*r148
fy(i)=fy(i)+ky
fy(j)=fy(j)-ky
kz=zz*r148
fz(i)=fz(i)+kz
fz(j)=fz(j)-kz

!--パリネロ・ラーマンの方法時のビリアル定理第二項（ポテンシャル依存項）の算出。ただし、二体ポテンシャル形式である。時間は一回目でt-⊿t、二回目でt=tである。--
!--式としてはtextにも載っていて(8,8)式である。--
cstress3(1,1)=cstress3(1,1)+r148*xx*xx
cstress3(1,2)=cstress3(1,2)+r148*xx*yy
cstress3(1,3)=cstress3(1,3)+r148*xx*zz
cstress3(2,2)=cstress3(2,2)+r148*yy*yy
cstress3(2,3)=cstress3(2,3)+r148*yy*zz
cstress3(3,3)=cstress3(3,3)+r148*zz*zz


  270 end do
      end do
      epot=epot+vn
!--上のvnはエワルドより。下では圧力テンソルは対称行列であることから同じ値を代入している。（一般に角運動量保存則に従っているならば対称性がいえる。参考wikipedia「応力]）
cstress3(2,1)=cstress3(1,2)
cstress3(3,1)=cstress3(1,3)
cstress3(3,2)=cstress3(2,3)
cstress3=cstress3/V

!--cstress1,2に=0.0を代入しておく必要はあるか?(上で代入した)--
cstress=cstress0+cstress1+cstress2+cstress3
pint=cstress

!do i=1,npart
!   write(*,*)fx(i),fy(i),fz(i)
!end do

      return
      end
!---------------------------------------------------------------------
      subroutine ewald

      use common_variables
!      vir=0.d0
        call cvn1
        call cvn2
        call cvn3

      vn=vn1+vn2+vn3
      return
      end
!---------------------------------------------------------------------
subroutine cvn0

use common_variables

!--パリネロ・ラーマン法でのビリアル定理の第一項（ビリアル依存項）の算出。--
cstress=0.d0

do i=1,npart
 cstress0(1,1)=cstress0(1,1)+mar*vx(i)*vx(i)
 cstress0(1,2)=cstress0(1,2)+mar*vx(i)*vy(i)
 cstress0(1,3)=cstress0(1,3)+mar*vx(i)*vz(i)
 cstress0(2,2)=cstress0(2,2)+mar*vy(i)*vy(i)
 cstress0(2,3)=cstress0(2,3)+mar*vy(i)*vz(i)
 cstress0(3,3)=cstress0(3,3)+mar*vz(i)*vz(i)
end do

!--第二項同様に圧力テンソルは対称行列であるので--
cstress0(2,1)=cstress0(1,2)
cstress0(3,1)=cstress0(1,3)
cstress0(3,2)=cstress0(2,3)
cstress0=cstress0/V
return
end

!---------------------------------------------------------------------
      subroutine cvn1


      use common_variables
      real*8 xa,erfc,xx,yy,zz,xi,yi,zi,xxx,yyy,zzz
      real*8 kx,ky,kz,r149

cstress1=0.d0
vn1=0.d0

do i=1,npart
 xi=qi(1,i)
 yi=qi(2,i)
 zi=qi(3,i)
do j=i+1,npart
 xx=xi-qi(1,j)
 yy=yi-qi(2,j)
 zz=zi-qi(3,j)

        if (xx.lt.-0.5d0) xx=xx+1.d0
        if (xx.gt. 0.5d0) xx=xx-1.d0

        if (yy.lt.-0.5d0) yy=yy+1.d0
        if (yy.gt. 0.5d0) yy=yy-1.d0

        if (zz.lt.-0.5d0) zz=zz+1.d0
        if (zz.gt. 0.5d0) zz=zz-1.d0

xxx=box(1,1)*xx+box(1,2)*yy+box(1,3)*zz
yyy=box(2,1)*xx+box(2,2)*yy+box(2,3)*zz
zzz=box(3,1)*xx+box(3,2)*yy+box(3,3)*zz
xx=xxx
yy=yyy
zz=zzz
        rd=xx*xx+yy*yy+zz*zz
        rds=sqrt(rd)

        if(rds<=sideh)then
        xa=alp*rds
        erfc=1.d0-erf(xa)

vn1=vn1+qe(i)*qe(j)*erfc/rds
        r148=-qe(i)*qe(j)*(erfc+2.d0/sqrt(pi)*xa*exp(-xa*xa))/rd/rds

r149=-r148
cstress1(1,1)=cstress1(1,1)+r149*xx*xx
cstress1(1,2)=cstress1(1,2)+r149*xx*yy
cstress1(1,3)=cstress1(1,3)+r149*xx*zz
cstress1(2,2)=cstress1(2,2)+r149*yy*yy
cstress1(2,3)=cstress1(2,3)+r149*yy*zz
cstress1(3,3)=cstress1(3,3)+r149*zz*zz

        kx=xx*r148
        fx(i)=fx(i)-kx
        fx(j)=fx(j)+kx
        ky=yy*r148
        fy(i)=fy(i)-ky
        fy(j)=fy(j)+ky
        kz=zz*r148
        fz(i)=fz(i)-kz
        fz(j)=fz(j)+kz
        end if

      end do
      end do
cstress1(2,1)=cstress1(1,2)
cstress1(3,1)=cstress1(1,3)
cstress1(3,2)=cstress1(2,3)
cstress1=cstress1/V
      return
      end
!---------------------------------------------------------------------
      subroutine cvn2


      use common_variables

real*8 xc,h2,LL,MM,exc,LM,LLI,MMI,LL2,MM2,boxreg(3),boxreg2(3,3),xi,yi,zi
real*8 sini(3,-20:20,nt),cosi(3,-20:20,nt),cosxb(nt),sinxb(nt),cos23,sin23
integer hx2,hy2,hz2
cstress2=0.d0
vn2=0.d0
sini=0.d0
cosi=0.d0
do i=1,npart
 xi=qi(1,i)
 yi=qi(2,i)
 zi=qi(3,i)

 sini(1,0,i)=0.d0
 cosi(1,0,i)=1.d0
 sini(2,0,i)=0.d0
 cosi(2,0,i)=1.d0
 sini(3,0,i)=0.d0
 cosi(3,0,i)=1.d0

 sini(1,1,i)=sin(2.d0*pi*xi)
 cosi(1,1,i)=cos(2.d0*pi*xi)
 sini(2,1,i)=sin(2.d0*pi*yi)
 cosi(2,1,i)=cos(2.d0*pi*yi)
 sini(3,1,i)=sin(2.d0*pi*zi)
 cosi(3,1,i)=cos(2.d0*pi*zi)
 do j=1,3
  do k=2,hmax
   sini(j,k,i)=sini(j,k-1,i)*cosi(j,1,i)+cosi(j,k-1,i)*sini(j,1,i)
   cosi(j,k,i)=cosi(j,k-1,i)*cosi(j,1,i)-sini(j,k-1,i)*sini(j,1,i)
  end do
  do k=-hmax,-1
   sini(j,k,i)=-sini(j,-k,i)
   cosi(j,k,i)=cosi(j,-k,i)
  end do
 end do
end do
      do hx2=-hmax,hmax
      do hy2=-hmax,hmax
      do hz2=-hmax,hmax

boxreg=0.d0

      LL=0.d0
      MM=0.d0

      if(hx2/=0.or.hy2/=0.or.hz2/=0)then

boxreg(1)=boxinv(1,1)*hx2+boxinv(2,1)*hy2+boxinv(3,1)*hz2
boxreg(2)=boxinv(1,2)*hx2+boxinv(2,2)*hy2+boxinv(3,2)*hz2
boxreg(3)=boxinv(1,3)*hx2+boxinv(2,3)*hy2+boxinv(3,3)*hz2
boxreg=2.d0*pi*boxreg

      do i=1,npart
cos23=cosi(2,hy2,i)*cosi(3,hz2,i)-sini(2,hy2,i)*sini(3,hz2,i)
sin23=sini(2,hy2,i)*cosi(3,hz2,i)+cosi(2,hy2,i)*sini(3,hz2,i)
cosxb(i)=cosi(1,hx2,i)*cos23-sini(1,hx2,i)*sin23
sinxb(i)=sini(1,hx2,i)*cos23+cosi(1,hx2,i)*sin23
        LL=LL+qe(i)*cosxb(i)
        MM=MM+qe(i)*sinxb(i)

      end do

      LL2=LL*LL
      MM2=MM*MM

      LM=LL2+MM2

h2=boxreg(1)*boxreg(1)+boxreg(2)*boxreg(2)+boxreg(3)*boxreg(3)
xc=-1.d0*h2/(alp*alp*4.d0)
      exc=exp(xc)

vn2=vn2+2.d0*pi/V/h2*exc*LM

r148=exc/h2*LM

boxreg2=0.d0
boxreg2(1,1)=boxreg(1)*boxreg(1)
boxreg2(1,2)=boxreg(1)*boxreg(2)
boxreg2(1,3)=boxreg(1)*boxreg(3)
boxreg2(2,2)=boxreg(2)*boxreg(2)
boxreg2(2,3)=boxreg(2)*boxreg(3)
boxreg2(3,3)=boxreg(3)*boxreg(3)
boxreg2(2,1)=boxreg2(1,2)
boxreg2(3,1)=boxreg2(1,3)
boxreg2(3,2)=boxreg2(2,3)

cstress2(1,1)=cstress2(1,1)+(1-0.5d0*boxreg2(1,1)/(alp*alp)-2.d0*boxreg2(1,1)/h2)*r148
cstress2(1,2)=cstress2(1,2)+(-0.5d0*boxreg2(1,2)/(alp*alp)-2.d0*boxreg2(1,2)/h2)*r148
cstress2(1,3)=cstress2(1,3)+(-0.5d0*boxreg2(1,3)/(alp*alp)-2.d0*boxreg2(1,3)/h2)*r148
cstress2(2,2)=cstress2(2,2)+(1-0.5d0*boxreg2(2,2)/(alp*alp)-2.d0*boxreg2(2,2)/h2)*r148
cstress2(2,3)=cstress2(2,3)+(-0.5d0*boxreg2(2,3)/(alp*alp)-2.d0*boxreg2(2,3)/h2)*r148
cstress2(3,3)=cstress2(3,3)+(1-0.5d0*boxreg2(3,3)/(alp*alp)-2.d0*boxreg2(3,3)/h2)*r148

      do i=1,npart

        LLI=qe(i)*cosxb(i)
        MMI=qe(i)*sinxb(i)

r148=4.d0*pi*exc/h2/V*(MMI*LL-LLI*MM)

fx(i)=fx(i)+boxreg(1)*r148
fy(i)=fy(i)+boxreg(2)*r148
fz(i)=fz(i)+boxreg(3)*r148
        
      end do

      end if

      end do
      end do
      end do
cstress2(2,1)=cstress2(1,2)
cstress2(3,1)=cstress2(1,3)
cstress2(3,2)=cstress2(2,3)
cstress2=cstress2*2.d0*pi/V/V

      return
      end
!---------------------------------------------------------------------
      subroutine cvn3

      use common_variables

vn3=0.d0

      do i=1,npart

      vn3=vn3-1.d0*alp/sqrt(pi)*qe2

      end do
      return
      end

!--------------------------------------------------------------------
      subroutine NPT !--吉村さんの卒論中の式を用いて説明する--
      use common_variables

call npt1
call npt2
call npt3
call npt2
call npt1

      return
      end
!--------------------------------------------------------------------
subroutine npt1

use common_variables

integer, parameter :: nyosh=1
integer iyosh,ne,nresn
real*8 wdti2(nyosh),wdti4(nyosh),wdti8(nyosh),ww(nyosh),trvg,aa
real*8 sc1,sc2,sc3,uv1,uv2,uv3

ww(1)=1.d0
!ww(1)=1.d0/(2.d0-2.d0**(1.d0/3.d0))
!ww(2)=1.d0-2.d0*ww(1)
!ww(3)=ww(1)
!ww(1)=1.d0/(4.d0-4.d0**(1.d0/3.d0))
!ww(2)=ww(1)
!ww(3)=1.d0-4.d0*ww(1)
!ww(4)=ww(1)
!ww(5)=ww(1)
nresn=1
sc1=0.d0
sc2=0.d0
sc3=0.d0


!---start the multiple time step procedure---
do ne=1,nresn
do iyosh=1,nyosh

wdti2(iyosh)=ww(iyosh)*dt/ne/2.d0 !--=⊿t/2--
wdti4(iyosh)=wdti2(iyosh)/2.d0    !--=⊿t/4--
wdti8(iyosh)=wdti4(iyosh)/2.d0    !--=⊿t/8--

!--(2.32)式。Qq()を用いている。フーバー形式M≧2に発展させるためにはif文で分けた方がいいかも--
!---update the thermostat forces---
glogs(1)=(ekin+ekinb-(nt+9.d0)*tref)/Qq(1)
!do i=2,nnos
! glogs(i)=(Qq(i-1)*vlogs(i-1)*vlogs(i-1)-tref)/Qq(i)
!end do

!--(2.31)式--
!---update the box forces---
gboxg(1,1)=(ekin/nt+(pint(1,1)-pext)*V)/M
gboxg(1,2)=pint(1,2)*V/M
gboxg(1,3)=pint(1,3)*V/M
gboxg(2,2)=(ekin/nt+(pint(2,2)-pext)*V)/M
gboxg(2,3)=pint(2,3)*V/M
gboxg(3,3)=(ekin/nt+(pint(3,3)-pext)*V)/M
gboxg(2,1)=gboxg(1,2)    !--実対称行列であることより--
gboxg(3,1)=gboxg(1,3)
gboxg(3,2)=gboxg(2,3)

!--P6の下から３行目の真ん中の計算。psの⊿t/4の時間発展--
!---update the thermostat velocities---
vlogs(nnos)=vlogs(nnos)+glogs(nnos)*wdti4(iyosh)
!do i=1,nnos-1
! aa=exp(-wdti8(iyosh)*vlogs(nnos+1-i))
! vlogs(nnos-i)=vlogs(nnos-i)*aa*aa+wdti4(iyosh)*glogs(nnos-i)*aa
!end do

!--P6の下から２行目の計算。三つまとめて計算している。--
!---update the box velocities---
aa=exp(-wdti8(iyosh)*vlogs(1))
do i=1,3
 do j=1,3
  vboxg(i,j)=(vboxg(i,j)*aa+wdti4(iyosh)*gboxg(i,j))*aa
 end do
end do

!--P7の上から１行目の計算。熱浴座標の時間発展--
!---update the thermostat positions---
do i=1,nnos
 xlogs(i)=xlogs(i)+wdti2(iyosh)*vlogs(i)
end do

!--速度の圧力による補正。P6の最後の行--
!---update the particle velocities---
!--viの係数となる行列の成分計算--
trvg=(vboxg(1,1)+vboxg(2,2)+vboxg(3,3))/nt   !--圧力の速度の対角和を自由度で割っている。P6の下から１行目のexp中第三項の係数計算--
veigv(1,1)=vboxg(1,1)+trvg+vlogs(1)
veigv(1,2)=vboxg(1,2)
veigv(1,3)=vboxg(1,3)
veigv(2,2)=vboxg(2,2)+trvg+vlogs(1)
veigv(2,3)=vboxg(2,3)
veigv(3,3)=vboxg(3,3)+trvg+vlogs(1)
veigv(2,1)=veigv(1,2)
veigv(3,1)=veigv(1,3)
veigv(3,2)=veigv(2,3)

!--対角化実行。対角化の理由は計算時間削減の為。対角化したい行列はveigvで、対角化されて帰ってくる行列の固有値がveigだと思われる--
veig=0.d0
call RsEIGQR( veigv, veig, NN, MTRX, EPSI, IFVEC, IFSRT,  &
    &                   CHECK, W1R, W1I, work )

sc1=exp(-veig(1)*wdti2(iyosh)) 
sc2=exp(-veig(2)*wdti2(iyosh))
sc3=exp(-veig(3)*wdti2(iyosh))
do i=1,npart    !--後ろから順次計算P7(2.40)式--
 uv1=vx(i)*veigv(1,1)+vy(i)*veigv(2,1)+vz(i)*veigv(3,1) !--直交行列の逆行列（対称行列)×ベクトルviの計算--
 uv2=vx(i)*veigv(1,2)+vy(i)*veigv(2,2)+vz(i)*veigv(3,2)
 uv3=vx(i)*veigv(1,3)+vy(i)*veigv(2,3)+vz(i)*veigv(3,3)
 uv1=uv1*sc1
 uv2=uv2*sc2
 uv3=uv3*sc3
 vx(i)=uv1*veigv(1,1)+uv2*veigv(1,2)+uv3*veigv(1,3)
 vy(i)=uv1*veigv(2,1)+uv2*veigv(2,2)+uv3*veigv(2,3)
 vz(i)=uv1*veigv(3,1)+uv2*veigv(3,2)+uv3*veigv(3,3)
end do

!--ビリアル更新のため再計算--
!---get the particle kinetic energy---
call rouekin

!--ビリアル更新のためにこれも再計算--
!---update the box forces---
gboxg(1,1)=(ekin/nt+(pint(1,1)-pext)*V)/M
gboxg(1,2)=pint(1,2)*V/M
gboxg(1,3)=pint(1,3)*V/M
gboxg(2,2)=(ekin/nt+(pint(2,2)-pext)*V)/M
gboxg(2,3)=pint(2,3)*V/M
gboxg(3,3)=(ekin/nt+(pint(3,3)-pext)*V)/M
gboxg(2,1)=gboxg(1,2)
gboxg(3,1)=gboxg(1,3)
gboxg(3,2)=gboxg(2,3)

!--先ほどの逆の手順で計算--
!---update the box velocities---
aa=exp(-wdti8(iyosh)*vlogs(1))
do i=1,3
 do j=1,3
  vboxg(i,j)=(vboxg(i,j)*aa+wdti4(iyosh)*gboxg(i,j))*aa
 end do
end do
!---get the box kinetic energy---
do i=1,3
 do j=1,3
  vgtvg(i,j)=vboxg(1,i)*vboxg(1,j)+vboxg(2,i)*vboxg(2,j)+vboxg(3,i)*vboxg(3,j)
 end do
end do

!--圧力の速度が時間発展したために圧力ビリアルの再計算--
ekinb=(vgtvg(1,1)+vgtvg(2,2)+vgtvg(3,3))*M

!---update the thermostat forces---
glogs(1)=(ekin+ekinb-(nt+9.d0)*tref)/Qq(1)

!---update the thermostat velocities---
!do i=1,nnos-1
! aa=exp(-wdti8(iyosh)*vlogs(i+1))
! vlogs(i)=vlogs(i)*aa*aa+wdti4(iyosh)*glogs(i)*aa
!---update the thermostat forces---
! glogs(i+1)=(Qq(i)*vlogs(i)*vlogs(i)-tref)/Qq(i+1)
!!----------------------------------
!end do
vlogs(nnos)=vlogs(nnos)+glogs(nnos)*wdti4(iyosh)
!--------------------------------------

end do
end do
return
end
!--------------------------------------------------------------------
subroutine npt2 !--速度⊿t/2時間発展--

use common_variables

!---update the particle velocities---
do i=1,npart
 vx(i)=vx(i)+fx(i)*(dt/2.d0/mar)
 vy(i)=vy(i)+fy(i)*(dt/2.d0/mar)
 vz(i)=vz(i)+fz(i)*(dt/2.d0/mar)
end do

call rouekin

return
end
!--------------------------------------------------------------------
subroutine npt3

use common_variables

real*8 e2,e4,e6,e8
real*8 aa,aa2(3),arg2,poly,bb(3),ubox(3,3),u1,u2,u3,uv1,uv2,uv3

e2=1.d0/6.d0 !--３の階乗分の１--
e4=e2/20.d0 !--５の階乗分の１--
e6=e4/42.d0 !--７の階乗分の１--
e8=e6/72.d0 !--９の階乗分の１--

!---update the particle positions---
veigv=vboxg
veig=0.d0
aa=0.d0
call RsEIGQR( veigv, veig, NN, MTRX, EPSI, IFVEC, IFSRT,  &
    &                   CHECK, W1R, W1I, work )
do i=1,3
 aa=exp(veig(i)*dt/2.d0)
 aa2(i)=aa*aa !--P8の上から５行目の計算--
 arg2=(veig(i)*dt/2.d0)*(veig(i)*dt/2.d0) !--ここから双曲線関数の計算--
 poly=(((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.d0 !--双曲線関数マクローリン展開５項目までの計算--
 bb(i)=aa*poly*dt !--P8の上から６行目の算出--
end do
do i=1,npart 
 u1=x(i)*veigv(1,1)+y(i)*veigv(2,1)+z(i)*veigv(3,1)  !--P8(2.42)式の第一項の後ろから順次計算--
 u2=x(i)*veigv(1,2)+y(i)*veigv(2,2)+z(i)*veigv(3,2)
 u3=x(i)*veigv(1,3)+y(i)*veigv(2,3)+z(i)*veigv(3,3)
 uv1=vx(i)*veigv(1,1)+vy(i)*veigv(2,1)+vz(i)*veigv(3,1) !--P8(2.42)式の第二項の後ろから順次計算--
 uv2=vx(i)*veigv(1,2)+vy(i)*veigv(2,2)+vz(i)*veigv(3,2)
 uv3=vx(i)*veigv(1,3)+vy(i)*veigv(2,3)+vz(i)*veigv(3,3)
 u1=u1*aa2(1)+uv1*bb(1) !--この計算で(2.42)式の直交行列計算前までいく--
 u2=u2*aa2(2)+uv2*bb(2)
 u3=u3*aa2(3)+uv3*bb(3)
 x(i)=u1*veigv(1,1)+u2*veigv(1,2)+u3*veigv(1,3)  !--直交行列との積の計算で最後--
 y(i)=u1*veigv(2,1)+u2*veigv(2,2)+u3*veigv(2,3)
 z(i)=u1*veigv(3,1)+u2*veigv(3,2)+u3*veigv(3,3)
end do
!---update the box---
!--P8(2.43)式の後ろから計算。逆行列をかけていてその逆行列は直交行列のものであるから対称行列である。故に下のように計算するので注意--
do i=1,3
 do j=1,3
  ubox(i,j)=veigv(1,i)*box(1,j)+veigv(2,i)*box(2,j)+veigv(3,i)*box(3,j) 
 end do
end do

!--Ieをかける--
do i=1,3
 do j=1,3
  ubox(i,j)=ubox(i,j)*aa2(i)
 end do
end do

!--直交行列をかける。ここは普通の行列積の計算--
do i=1,3
 do j=1,3
  box(i,j)=veigv(i,1)*ubox(1,j)+veigv(i,2)*ubox(2,j)+veigv(i,3)*ubox(3,j)
 end do
end do
call inversebox !--逆行列の計算--
!--周期境界条件の適用。セル内に全ての粒子を戻す。--
do i=1,npart
 qi(1,i)=boxinv(1,1)*x(i)+boxinv(1,2)*y(i)+boxinv(1,3)*z(i)
 qi(2,i)=boxinv(2,1)*x(i)+boxinv(2,2)*y(i)+boxinv(2,3)*z(i)
 qi(3,i)=boxinv(3,1)*x(i)+boxinv(3,2)*y(i)+boxinv(3,3)*z(i)
 if (qi(1,i).lt.0) qi(1,i)=qi(1,i)+1.d0
 if (qi(1,i).gt.1.d0) qi(1,i)=qi(1,i)-1.d0
 if (qi(2,i).lt.0) qi(2,i)=qi(2,i)+1.d0
 if (qi(2,i).gt.1.d0) qi(2,i)=qi(2,i)-1.d0
 if (qi(3,i).lt.0) qi(3,i)=qi(3,i)+1.d0
 if (qi(3,i).gt.1.d0) qi(3,i)=qi(3,i)-1.d0
 x(i)=box(1,1)*qi(1,i)+box(1,2)*qi(2,i)+box(1,3)*qi(3,i)  !--サブルチforceに渡すのはqiである。xは粒子位置として力計算には用いない。--
 y(i)=box(2,1)*qi(1,i)+box(2,2)*qi(2,i)+box(2,3)*qi(3,i)
 z(i)=box(3,1)*qi(1,i)+box(3,2)*qi(2,i)+box(3,3)*qi(3,i)
end do
call force
return
end
!--------------------------------------------------------------------
subroutine MSST 
use common_variables

call msst1
call npt2
call msst2
call npt2
call msst1

return
end




!--------------------------------------------------------------------
subroutine msst1
use common_variables

integer, parameter :: nyosh=1
integer iyosh,ne,nresn
real*8 wdti2(nyosh),wdti4(nyosh),wdti8(nyosh),ww(nyosh),trvg,aa
real*8 sc1,sc2,sc3,uv1,uv2,uv3
real*8 apint

ww(1)=1.d0
nresn=1
sc1=0.d0
sc2=0.d0
sc3=0.d0

!---start the multiple time step procedure---
do ne=1,nresn
do iyosh=1,nyosh

wdti2(iyosh)=ww(iyosh)*dt/ne/2.d0 !--=⊿t/2--
wdti4(iyosh)=wdti2(iyosh)/2.d0    !--=⊿t/4--
wdti8(iyosh)=wdti4(iyosh)/2.d0    !--=⊿t/8--

!---update the box forces---
gboxg = 0.d0
!gboxg(1,1)=(ekin/dble(nt)+(pint(1,1)-pext)*V+tmass*Vs**2*(V0-V)*V/V0**2)/M
!apint = (pint(1,1)+pint(2,2)+pint(3,3))/3.d0
!gboxg(1,1)=((apint-pext)*V+tmass*Vs**2*(V0-V)*V/V0**2)/M
gboxg(1,1)=((pint(1,1)-pext)*V-tmass*Vs**2*(V0-V)*V/V0**2)/M
!write(*,*)"b",pint(1,1)

!---update the box velocities---
i=1
j=1
vboxg(i,j)=vboxg(i,j)+wdti4(iyosh)*gboxg(i,j)

!---update the particle velocities---
!trvg=vboxg(1,1)/dble(nt) 
!veigv(1,1)=trvg
veigv(1,1)=vboxg(1,1)!+trvg

!write(*,*)"b",trvg
!do i=1,3
!   do j=1,3
!      write(2222,*)"b",i,j,veigv(i,j)
!   end do
!end do

do i=1,npart    
   vx(i)=vx(i)*exp(-veigv(1,1)*wdti2(iyosh)) 
!   write(2222,*)"b",i,vx(i)
end do

!---get the particle kinetic energy---
call rouekin

!---update the box forces---
gboxg = 0.d0
!gboxg(1,1)=(ekin/dble(nt)+(pint(1,1)-pext)*V+tmass*Vs**2*(V0-V)*V/V0**2)/M
!apint = (pint(1,1)+pint(2,2)+pint(3,3))/3.d0
!gboxg(1,1)=((apint-pext)*V+tmass*Vs**2*(V0-V)*V/V0**2)/M
gboxg(1,1)=((pint(1,1)-pext)*V-tmass*Vs**2*(V0-V)*V/V0**2)/M

!---update the box velocities---
i=1
j=1
vboxg(i,j)=vboxg(i,j)+wdti4(iyosh)*gboxg(i,j)

!---get the box kinetic energy---
vgtvg(1,1)=vboxg(1,1)*vboxg(1,1)
ekinb=(vgtvg(1,1))*M

end do
end do

return
end subroutine




!--------------------------------------------------------------------
subroutine msst2
use common_variables

real*8 aa

!---update the particle positions---
do i=1,npart 
 x(i)=x(i) + vx(i)*dt/4.d0
 y(i)=y(i) + vy(i)*dt/4.d0
 z(i)=z(i) + vz(i)*dt/4.d0

 aa = vboxg(1,1)*dt/2.d0 
 !x(i)=x(i)*exp(aa) 

 x(i)=x(i) + vx(i)*dt/4.d0
 y(i)=y(i) + vy(i)*dt/4.d0
 z(i)=z(i) + vz(i)*dt/4.d0

 qi(1,i)=boxinv(1,1)*x(i)+boxinv(1,2)*y(i)+boxinv(1,3)*z(i)
 qi(2,i)=boxinv(2,1)*x(i)+boxinv(2,2)*y(i)+boxinv(2,3)*z(i)
 qi(3,i)=boxinv(3,1)*x(i)+boxinv(3,2)*y(i)+boxinv(3,3)*z(i)
 if (qi(1,i).lt.0) qi(1,i)=qi(1,i)+1.d0
 if (qi(1,i).gt.1.d0) qi(1,i)=qi(1,i)-1.d0
 if (qi(2,i).lt.0) qi(2,i)=qi(2,i)+1.d0
 if (qi(2,i).gt.1.d0) qi(2,i)=qi(2,i)-1.d0
 if (qi(3,i).lt.0) qi(3,i)=qi(3,i)+1.d0
 if (qi(3,i).gt.1.d0) qi(3,i)=qi(3,i)-1.d0
end do

!---update the box---
i=1
j=1
aa = vboxg(1,1)*dt
box(i,j)=box(i,j)*exp(aa)
!write(1111,*)clock,exp(aa),exp((pint(1,1)-pext)*V/M*dt/2.d0)
call inversebox 

do i=1,npart
 x(i)=box(1,1)*qi(1,i)+box(1,2)*qi(2,i)+box(1,3)*qi(3,i)  
 y(i)=box(2,1)*qi(1,i)+box(2,2)*qi(2,i)+box(2,3)*qi(3,i)
 z(i)=box(3,1)*qi(1,i)+box(3,2)*qi(2,i)+box(3,3)*qi(3,i)

!---update the particle positions---
 x(i)=x(i) + vx(i)*dt/4.d0
 y(i)=y(i) + vy(i)*dt/4.d0
 z(i)=z(i) + vz(i)*dt/4.d0

 aa = vboxg(1,1)*dt/2.d0
 !x(i)=x(i)*exp(aa)

 x(i)=x(i) + vx(i)*dt/4.d0
 y(i)=y(i) + vy(i)*dt/4.d0
 z(i)=z(i) + vz(i)*dt/4.d0

 qi(1,i)=boxinv(1,1)*x(i)+boxinv(1,2)*y(i)+boxinv(1,3)*z(i)
 qi(2,i)=boxinv(2,1)*x(i)+boxinv(2,2)*y(i)+boxinv(2,3)*z(i)
 qi(3,i)=boxinv(3,1)*x(i)+boxinv(3,2)*y(i)+boxinv(3,3)*z(i)
 if (qi(1,i).lt.0) qi(1,i)=qi(1,i)+1.d0
 if (qi(1,i).gt.1.d0) qi(1,i)=qi(1,i)-1.d0
 if (qi(2,i).lt.0) qi(2,i)=qi(2,i)+1.d0
 if (qi(2,i).gt.1.d0) qi(2,i)=qi(2,i)-1.d0
 if (qi(3,i).lt.0) qi(3,i)=qi(3,i)+1.d0
 if (qi(3,i).gt.1.d0) qi(3,i)=qi(3,i)-1.d0

 x(i)=box(1,1)*qi(1,i)+box(1,2)*qi(2,i)+box(1,3)*qi(3,i)  
 y(i)=box(2,1)*qi(1,i)+box(2,2)*qi(2,i)+box(2,3)*qi(3,i)
 z(i)=box(3,1)*qi(1,i)+box(3,2)*qi(2,i)+box(3,3)*qi(3,i)

end do

call force

return
end




!--------------------------------------------------------------------
subroutine msst3
use common_variables

real*8 e2,e4,e6,e8
real*8 aa,aah

e2=1.d0/6.d0 
e4=e2/20.d0 
e6=e4/42.d0 
e8=e6/72.d0 

!---update the particle positions---
do i=1,npart 
!    write(1111,*)"b",i,x(i),y(i)
 aa = vboxg(1,1)*dt 
 aah = aa/2.d0

! write(3333,*)"b",clock,i,x(i)

 !x(i)=x(i)*exp(aa)+vx(i)*dt*exp(aah)*sinh(aah)/aah

 !write(3333,*)"a",clock,i,x(i)
 x(i)=x(i) + vx(i)*dt
 y(i)=y(i) + vy(i)*dt
 z(i)=z(i) + vz(i)*dt
 qi(1,i)=boxinv(1,1)*x(i)+boxinv(1,2)*y(i)+boxinv(1,3)*z(i)
 qi(2,i)=boxinv(2,1)*x(i)+boxinv(2,2)*y(i)+boxinv(2,3)*z(i)
 qi(3,i)=boxinv(3,1)*x(i)+boxinv(3,2)*y(i)+boxinv(3,3)*z(i)
 if (qi(1,i).lt.0) qi(1,i)=qi(1,i)+1.d0
 if (qi(1,i).gt.1.d0) qi(1,i)=qi(1,i)-1.d0
 if (qi(2,i).lt.0) qi(2,i)=qi(2,i)+1.d0
 if (qi(2,i).gt.1.d0) qi(2,i)=qi(2,i)-1.d0
 if (qi(3,i).lt.0) qi(3,i)=qi(3,i)+1.d0
 if (qi(3,i).gt.1.d0) qi(3,i)=qi(3,i)-1.d0
! qvi(1,i)=boxinv(1,1)*vx(i)+boxinv(1,2)*vy(i)+boxinv(1,3)*vz(i)
! qvi(2,i)=boxinv(2,1)*vx(i)+boxinv(2,2)*vy(i)+boxinv(2,3)*vz(i)
! qvi(3,i)=boxinv(3,1)*vx(i)+boxinv(3,2)*vy(i)+boxinv(3,3)*vz(i)
!    write(1111,*)"a",i,x(i),y(i)
end do

 write(3333,*)clock,sinh(aah)/aah,exp(aah)

!---update the box---
i=1
j=1
aa = vboxg(1,1)*dt
box(i,j)=box(i,j)*exp(aa)
write(1111,*)clock,exp(aa),exp((pint(1,1)-pext)*V/M*dt/2.d0)

call inversebox 
do i=1,npart
 x(i)=box(1,1)*qi(1,i)+box(1,2)*qi(2,i)+box(1,3)*qi(3,i)  
 y(i)=box(2,1)*qi(1,i)+box(2,2)*qi(2,i)+box(2,3)*qi(3,i)
 z(i)=box(3,1)*qi(1,i)+box(3,2)*qi(2,i)+box(3,3)*qi(3,i)
! vx(i)=box(1,1)*qvi(1,i)+box(1,2)*qvi(2,i)+box(1,3)*qvi(3,i)  
! vy(i)=box(2,1)*qvi(1,i)+box(2,2)*qvi(2,i)+box(2,3)*qvi(3,i)
! vz(i)=box(3,1)*qvi(1,i)+box(3,2)*qvi(2,i)+box(3,3)*qvi(3,i)
end do

!write(*,*)"ok"

call force
return
end




!--------------------------------------------------------------------
subroutine DOLLS
use common_variables

call dolls1
call npt2
call dolls2
call npt2
call dolls1

return
end




!--------------------------------------------------------------------
subroutine dolls1
use common_variables

integer, parameter :: nyosh=1
integer iyosh,ne,nresn
real*8 wdti2(nyosh),wdti4(nyosh),wdti8(nyosh),ww(nyosh),trvg,aa
real*8 sc1,sc2,sc3,uv1,uv2,uv3
real*8 apint

ww(1)=1.d0
nresn=1
sc1=0.d0
sc2=0.d0
sc3=0.d0

!---start the multiple time step procedure---
do ne=1,nresn
do iyosh=1,nyosh

wdti2(iyosh)=ww(iyosh)*dt/ne/2.d0 !--=<87><99>t/2--
wdti4(iyosh)=wdti2(iyosh)/2.d0    !--=<87><99>t/4--
wdti8(iyosh)=wdti4(iyosh)/2.d0    !--=<87><99>t/8--

!---update the box velocities---
i=1
j=1
vboxg(i,j) = teps
veigv(i,j) = vboxg(i,j)

!---update the atomic velocities---
do i=1,npart
   vx(i)=vx(i)*exp(-veigv(1,1)*wdti2(iyosh))
end do

!---get the particle kinetic energy---
call rouekin

!---update the box velocities---
i=1
j=1
vboxg(i,j)= teps

!---get the box kinetic energy---
vgtvg(1,1)=vboxg(1,1)*vboxg(1,1)
ekinb=(vgtvg(1,1))*M

end do
end do

return
end subroutine




!--------------------------------------------------------------------
subroutine dolls2
use common_variables

real*8 aa

!---update the atomic positions---
do i=1,npart
 x(i)=x(i) + vx(i)*dt/4.d0
 y(i)=y(i) + vy(i)*dt/4.d0
 z(i)=z(i) + vz(i)*dt/4.d0

 aa = vboxg(1,1)*dt/2.d0
 x(i)=x(i)*exp(aa)

 x(i)=x(i) + vx(i)*dt/4.d0
 y(i)=y(i) + vy(i)*dt/4.d0
 z(i)=z(i) + vz(i)*dt/4.d0

 qi(1,i)=boxinv(1,1)*x(i)+boxinv(1,2)*y(i)+boxinv(1,3)*z(i)
 qi(2,i)=boxinv(2,1)*x(i)+boxinv(2,2)*y(i)+boxinv(2,3)*z(i)
 qi(3,i)=boxinv(3,1)*x(i)+boxinv(3,2)*y(i)+boxinv(3,3)*z(i)
 if (qi(1,i).lt.0) qi(1,i)=qi(1,i)+1.d0
 if (qi(1,i).gt.1.d0) qi(1,i)=qi(1,i)-1.d0
 if (qi(2,i).lt.0) qi(2,i)=qi(2,i)+1.d0
 if (qi(2,i).gt.1.d0) qi(2,i)=qi(2,i)-1.d0
 if (qi(3,i).lt.0) qi(3,i)=qi(3,i)+1.d0
 if (qi(3,i).gt.1.d0) qi(3,i)=qi(3,i)-1.d0
end do

!---update the box---
i=1
j=1
aa = vboxg(1,1)*dt
box(i,j)=box(i,j)*exp(aa)
!write(1111,*)clock,exp(aa),exp((pint(1,1)-pext)*V/M*dt/2.d0)
call inversebox

do i=1,npart
 x(i)=box(1,1)*qi(1,i)+box(1,2)*qi(2,i)+box(1,3)*qi(3,i)
 y(i)=box(2,1)*qi(1,i)+box(2,2)*qi(2,i)+box(2,3)*qi(3,i)
 z(i)=box(3,1)*qi(1,i)+box(3,2)*qi(2,i)+box(3,3)*qi(3,i)

!---update the atomic positions---
 x(i)=x(i) + vx(i)*dt/4.d0
 y(i)=y(i) + vy(i)*dt/4.d0
 z(i)=z(i) + vz(i)*dt/4.d0

 aa = vboxg(1,1)*dt/2.d0
 x(i)=x(i)*exp(aa)

 x(i)=x(i) + vx(i)*dt/4.d0
 y(i)=y(i) + vy(i)*dt/4.d0
 z(i)=z(i) + vz(i)*dt/4.d0

 qi(1,i)=boxinv(1,1)*x(i)+boxinv(1,2)*y(i)+boxinv(1,3)*z(i)
 qi(2,i)=boxinv(2,1)*x(i)+boxinv(2,2)*y(i)+boxinv(2,3)*z(i)
 qi(3,i)=boxinv(3,1)*x(i)+boxinv(3,2)*y(i)+boxinv(3,3)*z(i)
 if (qi(1,i).lt.0) qi(1,i)=qi(1,i)+1.d0
 if (qi(1,i).gt.1.d0) qi(1,i)=qi(1,i)-1.d0
 if (qi(2,i).lt.0) qi(2,i)=qi(2,i)+1.d0
 if (qi(2,i).gt.1.d0) qi(2,i)=qi(2,i)-1.d0
 if (qi(3,i).lt.0) qi(3,i)=qi(3,i)+1.d0
 if (qi(3,i).gt.1.d0) qi(3,i)=qi(3,i)-1.d0

 x(i)=box(1,1)*qi(1,i)+box(1,2)*qi(2,i)+box(1,3)*qi(3,i)
 y(i)=box(2,1)*qi(1,i)+box(2,2)*qi(2,i)+box(2,3)*qi(3,i)
 z(i)=box(3,1)*qi(1,i)+box(3,2)*qi(2,i)+box(3,3)*qi(3,i)

end do

call force

return
end




!---------------------------------------------------------------------
      subroutine result

      use common_variables

real*8 creal(3,3),pots,realV
!--11/12追加NVTフーバー形式--
real*8 sumsw,sumpsw

      if ((clock/1)*1.eq.clock) then
        ek=0.5d0*ekin
        epot=epot
        etot=ek+epot
        etot2=(pv*pv)/(2.d0*M)+Pe*V+ps*ps/(2.d0*Q) &
     &        +((3.d0*npart)*tref*s)
        etott=etot+etot2
        temp=tscale*ekin*e0/boltz
creal=cstress*e0/(1.d9*a0**3)
abcreal=(creal(1,1)+creal(2,2)+creal(3,3))/3.d0
pres=abcreal

sumsw=0.0
sumpsw=0.0
do i=1,Termo
  if(i==1)then
       sumsw=sumsw+((3.d0*npart)*tref*sw(i))
  else 
       sumsw=sumsw+tref*sw(i)
  end if
  sumpsw=sumpsw+(psw(i)**2.0)/(2.0*Qw(i))
end do
        if(arg.eq.0) then
          if (clock.eq.0) &
     &    write(6,'(1a6,5a15)') 'step','ek','epot','etot', &
     &                          'temp[K]','pres[GPa]'
          write(6,'(1i6,6f15.6)')clock,ek,epot,etot,temp,pres
!write(6,'(1a6,5a15)')'x(1)','vx(1)','xi(256)','vx(256)'
!write(6,'(1i6,6f15.6)')x(2),vx(1),x(254),vx(256)
!write(6,'(1a6,5a15)')'qi(1,1)','qi(1,256)','qi(2,256)','qi(3,256)'
!write(6,'(1i6,6f15.6)')qi(1,1),qi(1,256),qi(2,256),qi(3,256)
!write(6,'(1a6,5a15)')'fx(1)','fy(1)','fx(256)','fz(256)'
!write(6,'(1i6,6f15.6)')fx(1),fy(1),fx(256),fz(256)
        end if

if(arg.eq.1 .and. hoover==1) then
  hen=etot+ps*ps/(2.d0*Q)+((3.d0*npart)*tref*s)
  write(6,'(1i6,6f15.6)')clock,hen,temp,pres,etot
end if
if(arg.eq.1 .and. hoover==2)then
  hen=etot+sumpsw+sumsw
  write(6,'(1i6,6f15.6)')clock,hen,temp,pres,etot
end if
        
        if(arg.eq.2) then
          hen=etot+(pv*pv)/(2.d0*M)+Pe*V
          write(6,'(1i6,3f15.6)')clock,hen,temp,pres
        end if
        if(arg.eq.3) then
realV=0.d0
ekins=0.d0
pots=0.d0
hen1=0.d0
hen2=0.d0
hen3=0.d0
do i=1,nnos !--熱浴全ビリアル算出--
 ekins=ekins+Qq(i)*vlogs(i)*vlogs(i)
end do

!--圧力運動エネルギー算出--
ekinb=ekinb*0.5d0
!--熱浴全運動エネルギー算出--
ekins=ekins*0.5d0
!do i=2,nnos
! pots=pots+xlogs(i)
!end do
!--P5(2.25)式の５項目、圧力ポテンシャルの計算--
hen1=Pe*V
!--６項目、熱浴１の座標エネルギーの算出--
hen2=(nt+9.d0)*tref*xlogs(1)
!--７項目、熱浴が２個以上あった場合の座標エネルギーのM=1を省いた総和--
hen3=tref*pots
!--NPT,Nose-Hoover方式の保存量--
hen=etot+ekinb+ekins+hen1+hen2+hen3
!hen=etot+ekins+hen2+hen3
realV=V*a0*a0*a0*1.d30
write(*,6000) clock,ek,epot,etot,hen,ekinb,ekins,hen1,hen2,hen3,pres,temp,realV
        end if
write(50,6000) clock,ek,epot,etot,hen,ekinb,ekins,hen1,hen2,hen3,pres,temp,realV
 6000   format(1i6,30f15.8,f15.6)
!write(*,'(3f15.6)') ((creal(i,j),j=1,3),i=1,3)
      end if

if(arg.eq.9)then
ekinb=ekinb*0.5d0
hen1=Pe*(V-V0)
hen=-tmass*Vs**2*(V0-V)**2/2.d0/V0**2
hen=etot+ekinb+hen1+hen
write(*,6000) clock,ek,epot,etot,hen,ekinb,-tmass*Vs**2*(V0-V)**2/2.d0/V0**2,hen1,box(1,1),pint(1,1),temp,V/V0
end if

if(arg.eq.10)then
ekinb=ekinb*0.5d0
hen1=0.d0
hen=0.d0
hen=etot+ekinb+hen1+hen
write(*,6000) clock,ek,epot,etot,hen,ekinb,-tmass*Vs**2*(V0-V)**2/2.d0/V0**2,hen1,box(1,1),pint(1,1),temp,V/V0
end if

call MSD
call analyze_box
if(clock .gt. clocks+timemx-timea .and. mod(clock,10)==1)then
 call cfa
end if
if(vmd.eq.1 .and. clock == 0) call vmdfile
if(vmd.eq.1 .and. clock>clocks+timemx-2000 .and. mod(clock,10)==1) call vmdfile

      return
      end
!--------------------------------------------------------------------
      subroutine vmdfile
      use common_variables
      real arx,ary,arz
real*8 sigma

      sigma=a0*1.d10

 write(51,*) npart
 write(51,*) nt

      do i=1,npart
        arx=x(i)*sigma
        ary=y(i)*sigma
        arz=z(i)*sigma

  write(51,'(a,3f7.3)') 'Ar',arx,ary,arz

      end do

      return
      end

!--------------------------------------------------------------------
subroutine cfa

use common_variables

real*8 xij,yij,zij,cfr,cfrs
real*8 q1,q2,q3

xij=0.d0
yij=0.d0
zij=0.d0
cfr=0.d0
cfrs=0.d0

do i=1,npart
 do j=i+1,npart
  q1=qi(1,i)-qi(1,j)
  q2=qi(2,i)-qi(2,j)
  q3=qi(3,i)-qi(3,j)
  if(abs(q1)>0.5d0)q1=q1-dsign(1.d0,q1)
  if(abs(q2)>0.5d0)q2=q2-dsign(1.d0,q2)
  if(abs(q3)>0.5d0)q3=q3-dsign(1.d0,q3)
  xij=box(1,1)*q1+box(1,2)*q2+box(1,3)*q3
  yij=box(2,1)*q1+box(2,2)*q2+box(2,3)*q3
  zij=box(3,1)*q1+box(3,2)*q2+box(3,3)*q3
  cfr=xij*xij+yij*yij+zij*zij
  cfrs=sqrt(cfr)
  k=(cfrs+0.5d0*dr)/dr
  if(k .le. kkk)then
   if((i .le. npart) .and. (j .le. npart))then
    nagag(k)=nagag(k)+1
   else if((i .le. npart) .and. (j .gt. npart))then
    nagi(k)=nagi(k)+1
   else if((i .gt. npart) .and. (j .gt. npart))then
    nii(k)=nii(k)+1
   end if
  end if
 end do
end do
VV=VV+V
return
end
!--------------------------------------------------------------------
subroutine cfb

use common_variables

real*8 dnagag(kkk),dnagi(kkk),dnii(kkk)

VV=VV/(timea/10)
do i=1,kkk
 ra=dr*i
 dnagag(i)=2.d0*nagag(i)/(timea/10)
 dnagi(i)=nagi(i)/(timea/10)
 dnii(i)=2.d0*nii(i)/(timea/10)
 gagag=VV/(4.d0*pi*ra*ra*dr)/npart/(npart-1)*dnagag(i)
 gagi=VV/(4.d0*pi*ra*ra*dr)/npart/npart*dnagi(i)
 gii=VV/(4.d0*pi*ra*ra*dr)/npart/(npart-1)*dnii(i)
 ra=ra*a0*1.d10

  write(90,*)ra,gagag,gagi,gii

end do

return
end
!--------------------------------------------------------------------
subroutine MSD

use common_variables

real*8 xii,yii,zii,rrr,msdag,msdi,xl(nt),yl(nt),zl(nt),sp3(3)

msdag=0.d0
msdi=0.d0
xii=0.d0
yii=0.d0
zii=0.d0
rrr=0.d0

do i=1,npart
 xii=qi(1,i)-qi0(1,i)
 yii=qi(2,i)-qi0(2,i)
 zii=qi(3,i)-qi0(3,i)

 if(abs(xii).gt.0.5d0) qia(1,i)=qia(1,i)-dsign(1.d0,xii)
 if(abs(yii).gt.0.5d0) qia(2,i)=qia(2,i)-dsign(1.d0,yii)
 if(abs(zii).gt.0.5d0) qia(3,i)=qia(3,i)-dsign(1.d0,zii)

 if(clock>clocks+timemx-1000)then
  if(abs(xii).gt.0.5d0)then
   xii=qi(1,i)-dsign(1.d0,xii)
  else
   xii=qi(1,i)
  end if
  if(abs(yii).gt.0.5d0)then
   yii=qi(2,i)-dsign(1.d0,yii)
  else
   yii=qi(2,i)
  end if
  if(abs(zii).gt.0.5d0)then
   zii=qi(3,i)-dsign(1.d0,zii)
  else
   zii=qi(3,i)
  end if
  absxii(i)=absxii(i)+xii
  absyii(i)=absyii(i)+yii
  abszii(i)=abszii(i)+zii
 end if

 xii=qi(1,i)+qia(1,i)
 yii=qi(2,i)+qia(2,i)
 zii=qi(3,i)+qia(3,i)

 xl(i)=box(1,1)*xii+box(1,2)*yii+box(1,3)*zii
 yl(i)=box(2,1)*xii+box(2,2)*yii+box(2,3)*zii
 zl(i)=box(3,1)*xii+box(3,2)*yii+box(3,3)*zii

end do

sp3=0.d0
do i=1,npart
 sp3(1)=sp3(1)+mar*xl(i)
 sp3(2)=sp3(2)+mar*yl(i)
 sp3(3)=sp3(3)+mar*zl(i)
end do
sp3(1)=sp3(1)/(npart*mar)
sp3(2)=sp3(2)/(npart*mar)
sp3(3)=sp3(3)/(npart*mar)

do i=1,npart
 xii=xl(i)-x0(i)-sp3(1)+sp2(1)
 yii=yl(i)-y0(i)-sp3(2)+sp2(2)
 zii=zl(i)-z0(i)-sp3(3)+sp2(3)

 rrr=xii*xii+yii*yii+zii*zii
 if(i .le. npart)then
  msdag=msdag+rrr
 else
  msdi=msdi+rrr
 end if
 qi0(1,i)=qi(1,i)
 qi0(2,i)=qi(2,i)
 qi0(3,i)=qi(3,i)
end do
msdag=msdag/npart
msdi=msdi/npart
write(91,*)clock,msdag,msdi
return
end

!---------------------------------------------------------------------
subroutine analyze_box

use common_variables

real*8 absbox(3),boxtheta(3),inprobox(3)

absbox=0.d0
boxtheta=0.d0

do i=1,3
 do j=1,3
  absbox(i)=absbox(i)+box(j,i)*box(j,i)
 end do
end do

absbox=sqrt(absbox)

inprobox(1)=box(1,1)*box(1,2)+box(2,1)*box(2,2)+box(3,1)*box(3,2)
inprobox(2)=box(1,2)*box(1,3)+box(2,2)*box(2,3)+box(3,2)*box(3,3)
inprobox(3)=box(1,3)*box(1,1)+box(2,3)*box(2,1)+box(3,3)*box(3,1)

boxtheta(1)=acos(inprobox(1)/(absbox(1)*absbox(2)))
boxtheta(2)=acos(inprobox(2)/(absbox(2)*absbox(3)))
boxtheta(3)=acos(inprobox(3)/(absbox(3)*absbox(1)))

boxtheta=boxtheta*180.d0/pi
absbox=absbox*a0*1.d10

write(92,*)clock,absbox,boxtheta

return
end
!---------------------------------------------------------------------
subroutine absposition
use common_variables

real*8 absxiir(nt),absyiir(nt),absziir(nt)

absxiir=0.d0
absyiir=0.d0
absziir=0.d0

absxii=absxii/1000.d0
absyii=absyii/1000.d0
abszii=abszii/1000.d0

do i=1,npart
 absxiir(i)=box(1,1)*absxii(i)+box(1,2)*absyii(i)+box(1,3)*abszii(i)
 absyiir(i)=box(2,1)*absxii(i)+box(2,2)*absyii(i)+box(2,3)*abszii(i)
 absziir(i)=box(3,1)*absxii(i)+box(3,2)*absyii(i)+box(3,3)*abszii(i)
end do

absxiir=absxiir*a0*1.d10
absyiir=absyiir*a0*1.d10
absziir=absziir*a0*1.d10


 write(150,*) npart
 write(150,*) npart



 do i=1,npart
   write(150,'(a,3f7.3)') 'Ar ',absxiir(i),absyiir(i),absziir(i)
 end do


return
end
!---------------------------------------------------------------------

      subroutine NVE
      use common_variables

      call ver

      return
      end
!---------------------------------------------------------------------
      subroutine ver  !--marをfに含める予定--
      use common_variables

      do i=1,npart !--速度ベルレ法による数値積分１回目(t=t+⊿t時の計算)--
        vx(i)=vx(i)+fx(i)*(dt/2.d0/mar)
        vy(i)=vy(i)+fy(i)*(dt/2.d0/mar)
        vz(i)=vz(i)+fz(i)*(dt/2.d0/mar)
      end do

      do i=1,npart !--速度ベルレ法の位置算出には第三項目がある(P75(5.8)式)が⊿tの２乗のオーダーなので無視している。そのため式自体は⊿tの精度しかもっていない--
        x(i)=x(i)+vx(i)*dt
        y(i)=y(i)+vy(i)*dt
        z(i)=z(i)+vz(i)*dt
      end do

      do i=1,npart !--周期条件--
        if (x(i).lt.0) x(i)=x(i)+side
        if (x(i).gt.side) x(i)=x(i)-side
        if (y(i).lt.0) y(i)=y(i)+side
        if (y(i).gt.side) y(i)=y(i)-side
        if (z(i).lt.0) z(i)=z(i)+side
        if (z(i).gt.side) z(i)=z(i)-side
      end do

do i=1,npart      !--サブルチforceで使うので位置rをqに再び変換している。--
 qi0(1,i)=boxinv(1,1)*x(i)+boxinv(1,2)*y(i)+boxinv(1,3)*z(i)
 qi0(2,i)=boxinv(2,1)*x(i)+boxinv(2,2)*y(i)+boxinv(2,3)*z(i)
 qi0(3,i)=boxinv(3,1)*x(i)+boxinv(3,2)*y(i)+boxinv(3,3)*z(i)
end do

qi=qi0

      call force !--t=tにおける力を算出する--

      do i=1,npart !--速度ベルレ法での数値積分２回目(t=t+⊿t時)--
        vx(i)=vx(i)+fx(i)*(dt/2.d0/mar)
        vy(i)=vy(i)+fy(i)*(dt/2.d0/mar)
        vz(i)=vz(i)+fz(i)*(dt/2.d0/mar)
      end do
      
call rouekin


      return
      end
!--------------------------------------------------------------------
      subroutine NVT
      use common_variables
   if(hoover==1)then
      call nvt1
      call ver
      call nvt1
   end if
   
   if(hoover==2)then
      call nvt2
      call ver
      call nvt2
   end if

      return
      end
!---------------------------------------------------------------------
      subroutine nvt1  !--フーバー形式（M=1の場合が下記。ps,sではなくpη,ηと表記したらtext通り)--
      use common_variables
!--①psをt=t+⊿t/4に時間発展--⑤二回目ではt=t+3⊿t/4に時間発展--
      ps=ps+(ekin-(3.d0*npart)*tref)*(dt/4.d0)
!--②sをt=t+⊿t/2に時間発展--⑥二回目ではt=t+⊿tに時間発展
      s=s+ps/Q*(dt/2.d0)
!--③（および二回目⑦）vx,vy,vzを補正。！--       !--②、③（または⑥、⑦）はexp同士が可換である事から、どちらからでも計算していい。--
      do i=1,npart
        vx(i)=vx(i)*exp(-ps/Q*(dt/2.d0))
        vy(i)=vy(i)*exp(-ps/Q*(dt/2.d0))
        vz(i)=vz(i)*exp(-ps/Q*(dt/2.d0))
      end do
!--vx,vy,vzが時間発展したのでt=t+⊿t/2時のビリアル(ekin)を計算--
call rouekin
!--④計算したビリアルを使ってpsをt=t+⊿t/2に時間発展--⑧二回目ではpsをt=t+⊿tに時間発展--
      ps=ps+(ekin-(3.d0*npart)*tref)*(dt/4.d0)


      return
      end
!--------------------------------------------------------------------
      subroutine nvt2 !--うまく動いたらM=1の場合も使えるはず--
      use common_variables
!--フーバー形式M≧2によるNVT--
!--①一番外側の熱浴の運動量の時間発展。ここだけ最終項が付かない。一項だけの計算になる。--
G=(((psw(Termo-1)**2.0)/Qw(Termo-1))-tref)
psw(Termo)=psw(Termo)+G*dt/4.0
!--②外側の熱浴から段階を踏んで熱浴の運動量を時間発展。最終項があるので各iについて項が２項あるのでRESPA法の時間発展に従って各iで３回ずつの計算を行なう。i同士も隣接するものが可換でないのでこの順番で計算を行なう。
do i=Termo-1,1,-1
 psw(i)=psw(i)*exp(-psw(i+1)/Qw(i+1)*(dt/8.0))
   if(i==1) then
        G=(ekin-(3.d0*npart)*tref)
   else 
        G=(((psw(i-1)**2.0)/Qw(i-1))-tref)
   end if
 psw(i)=psw(i)+G*(dt/4.0)
 psw(i)=psw(i)*exp(-psw(i+1)/Qw(i+1)*(dt/8.0))
end do

!--③熱浴の座標を時間発展させる。③、④は可換なのでどちらから始めても良い。--
do i=1,Termo
 sw(i)=sw(i)+psw(i)/Qw(i)*(dt/2.0)
end do

!--④速度補正--
do i=1,npart
        vx(i)=vx(i)*exp(-psw(1)/Qw(1)*(dt/2.d0))
        vy(i)=vy(i)*exp(-psw(1)/Qw(1)*(dt/2.d0))
        vz(i)=vz(i)*exp(-psw(1)/Qw(1)*(dt/2.d0))
end do

!--速度補正をしたため、ビリアル計算--
call rouekin

!--⑤内側から熱浴の運動量を時間発展--
do i=1,Termo-1
 psw(i)=psw(i)*exp(-psw(i+1)/Qw(i+1)*(dt/8.0))
   if(i==1) then
         G=(ekin-(3.d0*npart)*tref)
   else 
         G=(((psw(i-1)**2.0)/Qw(i-1))-tref)
   end if
 psw(i)=psw(i)+G*(dt/4.0)
 psw(i)=psw(i)*exp(-psw(i+1)/Qw(i+1)*(dt/8.0))
end do

!--⑥RESPA法に従って最後に一番外側の熱浴の運動量を時間発展--
G=(((psw(Termo-1)**2.0)/Qw(Termo-1))-tref)
psw(Termo)=psw(Termo)+G*dt/4.0
      return
      end
!--------------------------------------------------------------------
      subroutine NPH
      use common_variables

      call nph1
      call npt3
      call nph1

      return
      end
!--------------------------------------------------------------------
      subroutine nph1
      use common_variables

      pv=pv+(3.d0*V*(((1.d0/(3.d0*V))*(ekin-vir))-Pe))*(dt/4.d0)

      V=V*exp(3.d0*pv/M*(dt/2.d0))

      do i=1,nt
        vx(i)=vx(i)*exp((-(1+1.0/nt)*pv/M)*(dt/2.d0))
        vy(i)=vy(i)*exp((-(1+1.0/nt)*pv/M)*(dt/2.d0))
        vz(i)=vz(i)*exp((-(1+1.0/nt)*pv/M)*(dt/2.d0))
      end do

call rouekin
      pv=pv+(3.d0*V*(((1.d0/(3.d0*V))*(ekin-vir))-Pe))*(dt/4.d0)

      P=(ekin-vir)/(3.d0*V)

      return
      end



!***********************************************************************
      SUBROUTINE RsEIGQR( AR, EVR, N, MTRX, EPSI, IFVEC, IFSRT,  &
     &                   CHECK, W1R, W1I, work )
!-----------------------------------------------------------------------
!    Eigenvalues and eigenvectors of real symmetric  matrix   1992/10/16
!
!     Householders method for reducing the matrix to tridiagonal form
!     QR method for obtaining the eigenvalues
!
!      up-date  1992/11/12
!-----------------------------------------------------------------------
!  ( input )
!     AR    ...... Real parts of the Hermitian matrix elements
!       ( AR_ij ) for i>=j must be set correctly.
!     N     ...... Order of the Hermitian matrix
!     IFVEC ......  = 1 : obtain eigenvectors
!                   = 0 : do not obtain eigenvectors
!     IFSRT ......  = 2 : sorting eigenvalues and eigenvectors
!                   = 1 : sorting eigenvalues only ( W1I: sorting list )
!                   = 0 : do not sorting
!     EPSI  ...... Accuracy
!  ( output )
!     AR    ...... Real parts of the eigenvectors
!      ( AR_ij ) for i=1...N is eigenvector of eigenvalue, EVR_j.
!     EVR   ...... The eigenvalues
!  *****  W1R, W1I ...... work area  *****
!-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-R,T-Z)
      CHARACTER CHECK*5
      LOGICAL   VEC
      DIMENSION      AR(MTRX,*),EVR(*)
      DIMENSION      W1R(*),W1I(*),work(*)
      DATA  ZERO / 1.0D-30 /
      DATA  ITMIN, KMAX / 2, 50 /

!      WRITE(*,*) 'ITMIN, KMAX'
!      READ(*,*) ITMIN,KMAX
      VEC = IFVEC .EQ. 1


!-----------------------------------------------------------------------
!----    Householders method for reducing the matrix to tridiagonal form


      DO 1000 K = 1, N-2
         K1 = K + 1

         T1 = 0.0
         DO 1100 I = K1, N
            T1 = T1 + AR(I,K)*AR(I,K)
 1100    CONTINUE
         IF( ABS(T1).LT.ZERO ) THEN
             DO 1110 J = K1, N
                AR(K,J) = 0.0
 1110        CONTINUE
             GO TO 1000
         ENDIF

         TS = SQRT(T1) * SIGN( 1.0D0, AR(K1,K) )

         A1 = 1.0D0/( T1 + AR(K1,K)*TS )
         A1 = SQRT(A1)

         W1R(K) = 0.0
         W1R(K1) = (AR(K1,K) + TS)*A1
         DO 1200 I = K+2, N
            W1R(I) = AR(I,K)*A1
 1200    CONTINUE
         DO 1210 I = K1, N
            AR(K,I) = W1R(I)
 1210    CONTINUE

!         DO 1300 I = K, N
!            W1I(I) = 0.0
!            DO 1310 J = K1, I
!               W1I(I) = W1I(I) + AR(I,J)*W1R(J)
! 1310       CONTINUE
!            DO 1320 J = I+1, N
!               W1I(I) = W1I(I) + AR(J,I)*W1R(J)
! 1320       CONTINUE
! 1300    CONTINUE
         DO I = K, N
            W1I(I) = 0.d0
         end do
         DO J = K1, N
            w1rj = W1R(J)
            DO I = J, N
               W1I(I) = W1I(I) + AR(I,J)*w1rj
            end do
         end do
         DO I = K, N
            tsum = W1I(I)
            DO J = I+1, N
               tsum = tsum + AR(J,I)*W1R(J)
            end do
            W1I(I) = tsum
         end do

         U1R = 0.0
         DO 1400 I = K1, N
            U1R = U1R + W1I(I)*W1R(I)
 1400    CONTINUE
         U1R = -0.5D0*U1R

         DO 1500 I = K, N
            W1I(I) = W1I(I) + U1R*W1R(I)
 1500    CONTINUE

         DO 1600 I = K, K+1
            AR(I,K) = AR(I,K) - W1R(I)*W1I(K) - W1I(I)*W1R(K)
 1600    CONTINUE
         DO 1610 I = K+2, N
            AR(I,K) = 0.0
 1610    CONTINUE
         DO 1620 J = K1, N
            wij = W1I(J)
            wrj = W1R(J)
            DO 1620 I = J, N
               !AR(I,J) = AR(I,J) - W1R(I)*W1I(J) - W1I(I)*W1R(J)
               AR(I,J) = AR(I,J) - W1R(I)*wij - W1I(I)*wrj
 1620    CONTINUE


!check
!         WRITE(2,*) 'K=',K
!         DO 9700 L1 = 1, N
!            WRITE(2,9000) (AR(L1,L2),AI(L1,L2),L2=1,N)
! 9700    CONTINUE
! 9000    FORMAT(100('(',2E12.4') '))

 1000 CONTINUE




!-----------------------------------------------------------------------
!----   Store matrix elements
!----      EVR        :  diagonal elements ( real )
!----      W1R        :  non-diagonal (i,i+1)-elements


      DO 2000 I = 1, N
         EVR(I) = AR(I,I)
 2000 CONTINUE
      DO 2010 I = 1, N-1
         W1R(I) = AR(I+1,I)
 2010 CONTINUE


      IF( .NOT.VEC ) GO TO 10
!-----------------------------------------------------------------------
!----  Set the Unitary matrix elements


      DO 2100 I = 1, N
         AR(I,I) = 1.0D0
 2100 CONTINUE
!      DO 2110 I = 2, N
!      DO 2110 J = 1, I - 1
!         AR(I,J) = 0.0
! 2110 CONTINUE
      DO J = 1, N - 1
      DO I = J + 1, N
         AR(I,J) = 0.d0
      end do
      end do


      DO 2200 K = 1, N-2
         KN  = N - K
         KN1 = KN - 1

         DO 2205 J = KN+1, N
            AR(KN,J) = 0.0
 2205    CONTINUE
         DO J = KN, N
            w1i(j) = AR(KN1,J)
         enddo
         DO 2220 I = KN, N
            !work(i) = 0.0
            tsum = 0.d0
            DO J = KN, N
               !work(i) = work(i) + w1i(j)*AR(J,I)
               tsum = tsum + w1i(j)*AR(J,I)
            end do
            work(i) = tsum
 2220    CONTINUE
         DO 2210 I = KN, N
            VVR = -work(i)
            DO 2210 J = KN, N
               AR(J,I) = AR(J,I) + w1i(j)*VVR
 2210    CONTINUE
 2200 CONTINUE

      DO 2300 K = 2, N
         AR(1,K) = 0.0
 2300 CONTINUE


!check
!         WRITE(2,*) 'U'
!         DO 710 L1 = 1, N
!            WRITE(2,9000) (AR(L1,L2),L2=1,N)
!  710    CONTINUE





   10 CONTINUE
!-----------------------------------------------------------------------
!----               Givens' QR method


      DMX = 0.0
      WMX = 0.0
      DO 3000 I = 1, N
         DMX = MAX( DMX, ABS(EVR(I)) )
 3000 CONTINUE
      DO 3010 I = 1, N-1
         WMX = MAX( WMX, ABS(W1R(I)) )
 3010 CONTINUE

      ESHIFT = MAX( WMX, DMX ) * 1.1D0
      DO 3020 I = 1, N
         EVR(I) = EVR(I) + ESHIFT
 3020 CONTINUE

!check
!         WRITE(*,*) 'A_initial'
!         WRITE(*,9002) (EVR(L1),L1=1,N)
!         WRITE(*,9002) (W1R(L1),L1=1,N-1)
!         WRITE(*,9002) (W1I(L1),L1=1,N-1)
! 9002    FORMAT(6E12.4)


      KAI   = 0
      KNTSF = 0
      KNM   = N

    1 CONTINUE
      EVRSFT = 0.0
      GO TO 3
    2 CONTINUE
!----   eigenvalue shift
        EVRSFT = EVR(KNM)
        DO 3030 I = 1, KNM
           EVR(I) = EVR(I) - EVRSFT
 3030   CONTINUE


    3 CONTINUE
      KAI   = KAI   + 1
      KNTSF = KNTSF + 1

      COLD = 1.0D0
      DO 3100 K = 2, KNM
         K1 = K - 1
         AA   = EVR(K1)*EVR(K1) + W1R(K1)*W1R(K1)
         AA   = SQRT(AA)
!           IF( ABS(AA).LT.ZERO ) THEN
!               WRITE(*,*) AA
!               CHECK = 'ZERO '
!               RETURN
!           ENDIF
         DCOSR = EVR(K1)/AA
         DSINR = W1R(K1)/AA

         EK1R = W1R(K1)*COLD

         ARK  = EVR(K)
         W1NR = EK1R*DCOSR + ARK*DSINR
         EVR(K) = -EK1R*DSINR + ARK*DCOSR

         BK1 = COLD*DCOSR

         EVR(K1) = AA*BK1 + W1NR*DSINR
         IF( K.GE.3 ) THEN
             W1R(K-2) = AA*DSOLR
         ENDIF

         COLD  = DCOSR
         DSOLR = DSINR
         IF( VEC ) THEN
             DO 3150 I = 1, N
                ARK1 = AR(I,K-1)
                ARK  = AR(I,K)
                AR(I,K-1) = ARK1*DCOSR + ARK*DSINR
                AR(I,K)   = -ARK1*DSINR + ARK*DCOSR
 3150        CONTINUE
         ENDIF
 3100 CONTINUE
      ARK1 = EVR(KNM)
      EVR(KNM)   = ARK1*COLD
      W1R(KNM-1) = ARK1*DSOLR



!----   diagonal part

      DO 3200 I = 1, KNM
         EVR(I) = EVR(I) + EVRSFT
 3200 CONTINUE

!      WRITE(*,*) 'kaisuu',KAI
!      WRITE(*,9002) (EVR(L1),L1=1,N)
!      WRITE(*,9002) (EVR(L1)-ESHIFT,L1=1,N)
!C         WRITE(*,9002) (W1R(L1),L1=1,N-1)
!C         WRITE(*,9002) (W1I(L1),L1=1,N-1)


!----   hantei

      IF( KNTSF.GT.KMAX ) THEN
          CHECK = 'KAISU'
          RETURN
      ENDIF
      DMX = 0.0
      DO 3300 I = 1, KNM
         DMX = MAX( DMX, ABS(EVR(I)) )
 3300 CONTINUE
      DMX = DMX*EPSI

      IF( ABS(W1R(KNM-1)).LT.DMX ) THEN
          KNM = KNM - 1
          KNTSF = 0
      ENDIF

      DO 3310 I = 1, KNM - 1
         IF( ABS(W1R(I)).GT.DMX ) then
         IF( KAI.lE.ITMIN ) then
             go to 1
          else
             go to 2
         endif
         endif
 3310 CONTINUE


!      WRITE(*,*) '*** iteration :',KAI
!-----------------------------------------------------------------------

      DO 3400 I = 1, N
         EVR(I) = EVR(I) - ESHIFT
 3400 CONTINUE


      IF( IFSRT.LE.0 ) RETURN
!-----------------------------------------------------------------------
!----    sorting eigenvalues


      W1I(1)= 1

      DO 3500 I = 2, N
         DO 3510 J = 1, I-1
         IF( EVR(I).GE.EVR(J) ) GO TO 3510

            EVRI = EVR(I)
            DO 3520 L = I, J+1, -1
               EVR(L) = EVR(L-1)
               W1I(L) = W1I(L-1)
 3520       CONTINUE
            EVR(J) = EVRI
            W1I(J) = I
            GO TO 3500

 3510    CONTINUE
         W1I(I) = I
 3500 CONTINUE



      IF( .NOT.VEC .OR. IFSRT.LE.1 ) RETURN

!----    sorting eigenvectors


      DO 3600 I = 1, N
         W1R(I) = 0
 3600 CONTINUE
      DO 3610 I = 1, N
         IF( NINT(W1R(I)).NE.0 ) GO TO 3610
         W1R(I) = 10.0
         IT = NINT(W1I(I))
         IF( I.EQ.IT ) GO TO 3610
             IOR = I
             IO  = I

 3611        CONTINUE
             DO 3620 L=1,N
                ARI = AR(L,IO)
                AR(L,IO)  = AR(L,IT)
                AR(L,IT) = ARI
 3620        CONTINUE
             IO = IT
             W1R(IO) = 10.0
             IT = NINT(W1I(IO))
             IF( IT.EQ.IOR ) GO TO 3610
             GO TO 3611

 3610 CONTINUE


      RETURN
      END
