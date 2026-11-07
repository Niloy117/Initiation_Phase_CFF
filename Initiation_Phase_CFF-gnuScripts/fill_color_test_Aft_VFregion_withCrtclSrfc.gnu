
max(x,y)=(x<y)?y:x
ftr(x,alpha)=abs(x)**(1.+alpha)/x
alpha=1.0

#set term x11 enhanced font "arial,15" 

Lx=640
xmin=-10.50 #-42.00
xmax=10.50 #+42.00
ymin=-7.00 #-44.50
ymax=+2.00
Ly=floor(100*(ymax-ymin)/(xmax-xmin))
print Lx,Ly

set terminal postscript eps enhanced color dashed lw 1 "Times-New-Roman" 9 #size 1800, 1200
set size Lx,Ly

set size ratio -1
set angles degrees

unset key
#unset border
#set notics
#set noxtics
#set noytics

set xrange [-22.00:22.00]
set yrange [-20.50:2.00]


#set ytics -20.50,22.50,2.00 offset 5

#set xtics nomirror
#set xtics -22.00,44.00,22.00 scale 0.1 font "Times-New-Roman,1" #THIS LINE VERY IMPORTANT FOR canvas SIZE
set xtics

# define the location of your plot:
bm = 0.15
lm = 0.12
rm = 0.75
tm = 0.90

# letter spacing - play with this as needed:
STRDIST = 0.022

# set up the plot window:
set lmargin at screen lm
set rmargin at screen rm
set bmargin at screen bm
set tmargin at screen tm


set palette model RGB defined (-1 "red",0 "white", 1 "blue")
cbr=1.
set cbrange [-cbr:cbr]
set cbtics ('Compression' 1.0, '0' 0, 'Tension' -1.0)
set colorbox vertical user origin rm+0.01,0.41 size .015,0.22*(tm-bm)


# your label
LABEL1 = "Tension    "
LABEL2 = "Compression"

N_iter = 145
nstrt  = 1
nend   = N_iter

Hlf_Ncell       = 11
ncells          = (2*Hlf_Ncell)+1+1
NCP_CrtclApSrfc = 1

NAEC_apcl = 2 ; NAEC_bsal=2 ; NAEC_ltrl=2 ; NAEC_apclCrtclApSrfc=5

nsprsInACellwithoutNode           = 3
nsprsInACellwiththeNode           = (NAEC_apcl+1)+(NAEC_bsal+1)+(NAEC_ltrl+1) 
nsprsInACellwithtNodeAndCrtclNode = (NAEC_apclCrtclApSrfc+1)+(NAEC_bsal+1)+(NAEC_ltrl+1)

curvesInACell  = 3
nsprsInACell   = nsprsInACellwiththeNode
nsprsInACellT2 = nsprsInACellwithtNodeAndCrtclNode

#nsprsBfrMet = (ncells-2)*(nsprsInACell) + 2 + 1*(NAEC_apcl+1) + 1*(NAEC_bsal+1)
nsprsBfrMet = (ncells-2-2*NCP_CrtclApSrfc)*(nsprsInACell) + (2*NCP_CrtclApSrfc)*(nsprsInACellT2) + 2 + 1*(NAEC_apclCrtclApSrfc+1) + 1*(NAEC_bsal+1) 
nsprsAftMet = (ncells-2)*(nsprsInACell) + 3 + 1*(NAEC_bsal+1) + 2 -1 + (1)

nsprs       = nsprsBfrMet
ncurve      = (ncells-2)*(curvesInACell) + 1
print nsprsBfrMet,nsprsAftMet," nsprsBfr and nsprsAft"

numRglrNode           = (2*(Hlf_Ncell+1))*2
numNodeBFRcrtclApSrfc = 2*(Hlf_Ncell-NCP_CrtclApSrfc)*(nsprsInACellwiththeNode-3)
numNodeAFTcrtclApSrfc = (2*NCP_CrtclApSrfc*(nsprsInACellwithtNodeAndCrtclNode-3)) + (NAEC_apclCrtclApSrfc+NAEC_bsal)
Nnode                 = numRglrNode + numNodeBFRcrtclApSrfc + numNodeAFTcrtclApSrfc
print numRglrNode,numNodeBFRcrtclApSrfc,numNodeAFTcrtclApSrfc,Nnode," node stuff"

array fcolorArea[ncells]

ItrnNoBfrMeet=200

stressScaleA = 0.38 #0.23 #0.24 #0.30 #0.38 #0.37 #0.17 #0.37 #0.71 #0.21
stressScaleS = 1.58 #1.45 #1.47 #1.52 #1.61 #1.10 #0.64 #1.18 #2.40 #0.87

stressMaxA=0.
stressMaxS=0.

lpMax=int((N_iter-nstrt-1)/20)
print lpMax

lpEnd=(lpMax*20)+(nstrt-1)
print lpEnd,"lpEnd"

diffLpEnd=N_iter-lpEnd
print N_iter,lpEnd,diffLpEnd,"chk"

if (diffLpEnd==0) {gnuStepMax=lpMax}
if (diffLpEnd>0)  {gnuStepMax=lpMax+1}

print diffLpEnd,gnuStepMax," ll"

do for [gnuStep=1:gnuStepMax]{
   
   nstrtLp = nstrt  + (20*(gnuStep-1))
   if (gnuStep!=gnuStepMax) {nend = (nstrt-1+20) + (20*(gnuStep-1))}
   if (gnuStep==gnuStepMax) {
      if(diffLpEnd!=0) {nend = (nstrt-1) + diffLpEnd + (20*(gnuStep-1))}
      if(diffLpEnd==0) {nend = (nstrt-1+20) + (20*(gnuStep-1))}
      }
   print nstrtLp,diffLpEnd,nstrt,gnuStep,nend,"val"
   
   do for [n=nstrtLp:nend]{
      key_val=0+(n-1)*2
	
      if (n<=ItrnNoBfrMeet) {nsprs=nsprsBfrMet}
      else {
      nsprs = nsprsAftMet}
      	  
      print n,nsprs
      array fcolorSpr[nsprs]
   	
      set term x11
      	
      datadir='./'
	
      datafileArea = datadir.sprintf('NI_modelInitiatn_WT_VFAW%04d.dat',n)
      datafileSpr  = datadir.sprintf('NI_modelInitiatn_WT_VFSW%04d.dat',n)   
      datafileNode = datadir.sprintf('NI_modelInitiatn_WT_VFNW%04d.dat',n)
      
      do for [i=1:ncells]{
      	 plot datafileArea every 1:1:0:i-1:0:i-1 u (stressA=$1)
      	 plot datafileArea every 1:1:0:i-1:0:i-1 u (fcolorArea[i]=(1.+ftr($1/stressScaleA,alpha))/2.)
      	 stressMaxA=max(stressMaxA,abs(stressA))
	 }
      
      do for [i=1:nsprs]{
      	 plot datafileSpr every 1:1:0:i-1:0:i-1 u (stressS=$1)
      	 plot datafileSpr every 1:1:0:i-1:0:i-1 u (fcolorSpr[i]=(1.+ftr($1/stressScaleS,alpha))/2.)
      	 stressMaxS=max(stressMaxS,abs(stressS))
      	 }
     	 	 	 	 	 
   	print n,stressMaxA,stressMaxS
   	
   	#set terminal postscript landscape enhanced color dashed lw 1 "Times-New-Roman" 9
   	set terminal postscript eps enhanced color dashed lw 1 "Times-New-Roman" 9 #size 1800, 1200
  	 
	outputfile=sprintf('Exp0_%04d.eps',n)
   	set output outputfile
   
	set multiplot
   	#set title sprintf("frame%02d",(n+0))
   		
   	plot for [i=1:ncells] datafileArea every 1:1:1:i-1::i-1 w filledcurves lc palette frac fcolorArea[i] 
   
	plot \
   	for [i=1:nsprs] datafileSpr every 1:1:1:i-1::i-1 w l lt 1 lc rgb "grey" lw 4,\
   	for [i=1:nsprs] datafileSpr every 1:1:1:i-1::i-1 w filledcurves lc palette frac fcolorSpr[i] lw 2

	#plot datafileNode every ::1::48 w circle 
   
	unset multiplot
   	
   }
   pause 1
}

pause -1
