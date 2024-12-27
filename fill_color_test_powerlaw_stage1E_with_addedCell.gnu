
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


set xrange [-62.00:62.00]
set yrange [-20.50:2.00]

#set xrange [-11.00:11.00]
#set yrange [-10.25:2.00]


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


N_iter = 2
nstrt  = 1
nend   = N_iter

numOfAddedCell = 20
Hlf_Ncells     = 11 + numOfAddedCell
ncells         = (2*Hlf_Ncells)+1
print Hlf_Ncells,ncells,numOfAddedCell," hlf-Ncells"

nodeInsrtdOrNot = 1

if (nodeInsrtdOrNot == 0) {NAEC_apcl=0 ; NAEC_bsal=0 ; NAEC_ltrl=0}
if (nodeInsrtdOrNot == 1) {NAEC_apcl=0 ; NAEC_bsal=2 ; NAEC_ltrl=2}

nsprsInACellwithoutNode = 3
nsprsInACellwiththeNode = (NAEC_apcl+1)+(NAEC_bsal+1)+(NAEC_ltrl+1) 

curvesInACell = 2 #1 for basal and 1 for lateral

if (nodeInsrtdOrNot==0) {nsprsInACell=nsprsInACellwithoutNode}
if (nodeInsrtdOrNot==1) {nsprsInACell=nsprsInACellwiththeNode}

nsprsBfrMet = (ncells-1)*(nsprsInACell) + 3 + 1*(NAEC_bsal+1)
nsprsAftMet = (ncells-1)*(nsprsInACell) + 3 + 1*(NAEC_bsal+1) + 2 -1
nsprs       = nsprsBfrMet
ncurve      = (ncells-1)*(curvesInACell) + 1

print nsprsBfrMet,nsprsAftMet," nsprsBfr and nsprsAft"

array fcolorArea[ncells]

ItrnNoBfrMeet=11

stressScaleA = 0.36 #1 #0.28 #0.33 #0.37 #0.37
stressScaleS = 0.97 #1.05 #1.06 #1.01 #1.51 #0.76 #1.20 #1.12

stressMaxA=0.
stressMaxS=0.

do for [n=1:N_iter]{
   key_val=0+(n-1)*2
   
   if (n<=ItrnNoBfrMeet) {
      if (n<=1) {nsprs=nsprsBfrMet
      } else{
      nsprs=nsprsBfrMet+3
      }
   } else {
       nsprs = nsprsAftMet
       }
   
   print n,nsprs
   array fcolorSpr[nsprs]
   
   set term x11 #size 1800, 1200
   
   datadir='./'
   
   datafileArea = datadir.sprintf('NI_modelAddedCellAW%03d.dat',n)
   datafileSpr  = datadir.sprintf('NI_modelAddedCellSW%03d.dat',n)
   
   
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
   #set term gif
   
   outputfile=sprintf('Exp0_%03d.eps',n)  # eps, not gif 
   set output outputfile
   
   set multiplot
   #set title sprintf("frame%02d",(n+0))
   		
   plot for [i=1:ncells] datafileArea every 1:1:1:i-1::i-1 w filledcurves lc palette frac fcolorArea[i] 
   
   #lp pt 6
   
   plot \
   for [i=1:nsprs] datafileSpr every 1:1:1:i-1::i-1 w l lt 1 lc rgb "grey" lw 4,\
   for [i=1:nsprs] datafileSpr every 1:1:1:i-1::i-1 w filledcurves lc palette frac fcolorSpr[i] lw 2
   
   unset multiplot
   
}

pause -1

#set label i LABEL1[i:i] at screen 0.7+((i-1.)*STRDIST),bm-0.01 textcolor palette frac 0.0\