module param
implicit none
integer, parameter :: dp = selected_real_kind(10,100)
integer :: i, j, k, l, m, n, totatom, tbond, tmolec, last_pos=0, n_neighbor=0, a_neighbor=0, e_mapper=0, a_mapper=0, slt=1
integer :: ipsbond=0, ipsang=0, ipsdih=0, ipsmolecule=0, tstep=0, axsfstep=0, totsolute=0, atomorig=1, e_charge=0, a_charge=0
integer :: ntbody=0, cartcharge=0
real(kind=dp) :: last, cellparam=0.0d0, rcut=0.0d0, orig_x=0.0d0, orig_y=0.0d0, orig_z=0.0d0, cellfactor=0.0d0, tmp_mapper=0.0d0
real(kind=dp) :: cartmass=0.0d0
real(kind=dp), dimension(:), allocatable :: rxcoord, rycoord, rzcoord, rmap
real(kind=dp), dimension(:), allocatable :: xcoord, ycoord, zcoord, amap, samap
real(kind=dp), dimension(:), allocatable :: cmassx, cmassy, cmassz
real(kind=dp), dimension(:), allocatable :: molecmass
integer, dimension(:), allocatable :: molec, imolec, iamap, moleccharge, iacharge, sacharge
integer, dimension(:), allocatable :: acharge
integer, dimension(:,:), allocatable :: sbond, sang, sdih, itbody
real(kind=dp), dimension(:,:), allocatable :: bond, molec_bond
integer, dimension(:,:), allocatable :: molec_map
logical, dimension(:), allocatable :: ltbody
logical, dimension(:,:), allocatable :: lbond, molec_mask
logical :: full_cell=.False., ptrajec=.False., pcell=.False., pnns=.False., in_cell=.False., debug=.False., within=.False.
logical :: cbonds=.False.,  pbonds=.False., pangs=.False., pdih=.False., molecule=.False., in_cellcmass=.False., linp=.False.
logical :: pmolecule=.False., neighb=.False., aneighb=.False., pcart=.False., psolute=.False., paxsf=.False., topos=.False.
logical :: psbond=.False., psang=.False., psdih=.False., psmolecule=.False., pstrajec=.False., cmass=.False., psaxsf=.False.
logical :: lcharge=.False., tbody=.False.
character (len=20) :: keywords
character (len=3) :: ci
character (len=2), dimension(:), allocatable :: rsymb
character (len=2), dimension(:), allocatable :: symb
character (len=80) :: arg, pinp2
!
integer, parameter :: totalelem=118
real(kind=dp), parameter :: pi=3.14159265358979
character (len=2), dimension(totalelem) :: elem
real(kind=dp), dimension(totalelem) :: radiuscov
real(kind=dp), dimension(totalelem) :: mass
real(kind=dp), dimension(totalelem) :: mapper=1.100D0
integer, dimension(totalelem) :: charge=0

! SYMBOLS OF PERIODIC TABLE OF ELEMENTS  Z=1-118
data (elem(l), l=1,totalelem)/&
"H ","He",&
"Li","Be","B ","C ","N ","O ","F ","Ne",&
"Na","Mg","Al","Si","P ","S ","Cl","Ar",&
"K ","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn",&
"Ga","Ge","As","Se", "Br","Kr",&
"Rb","Sr","Y ","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",&
"In","Sn","Sb","Te","I ","Xe",&
"Cs","Ba",&
"La","Ce","Pr","Nd","Pm","Sm","Eu",&
"Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",&
"Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",&
"Tl","Pb","Bi","Po","At","Rn",&
"Fr","Ra",&
"Ac","Th","Pa","U ","Np","Pu","Am",&
"Cm","Bk","Cf","Es","Fm","Md","No","Lr",&
"Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn",&
"Nh","Fl","Mc","Lv","Ts","Og"/
! COVALENT RADII OF SINGLE BONDS IN ANGSTRONS
data (radiuscov(l), l=1,totalelem)/&
0.310,0.280,&
1.280,0.960,0.840,0.760,0.710,0.660,0.570,0.580,&
1.660,1.410,1.210,1.110,1.070,1.050,1.020,1.060,&
2.030,1.760,1.700,1.600,1.530,1.390,1.500,1.420,1.380,1.240,1.320,1.220,&
1.220,1.200,1.190,1.200,1.200,1.160,&
2.200,1.950,1.900,1.750,1.640,1.540,1.470,1.460,1.420,1.390,1.450,1.440,&
1.420,1.390,1.390,1.380,1.390,1.400,&
2.440,2.150,&
2.070,2.040,2.030,2.010,1.990,1.980,1.980,&
1.960,1.940,1.920,1.920,1.890,1.900,1.870,1.870,&
1.750,1.700,1.620,1.510,1.440,1.410,1.360,1.360,1.320,&
1.450,1.460,1.480,1.400,1.500,1.500,&
2.600,2.210,&
2.150,2.060,2.000,1.960,1.900,1.870,1.800,&
1.690,1.690,1.700,1.660,1.680,1.730,1.760,1.610,&
1.570,1.490,1.430,1.410,1.340,1.290,1.280,1.210,1.220,&
1.360,1.430,1.620,1.750,1.650,1.570/
!
! ATOMIC MASS
data (mass(l), l=1,totalelem)/&
1.00797,4.00260,&                                                                   
6.941,9.01218,10.81,12.011,14.0067,15.9994,18.998403,20.179,&                       
22.98977,24.305,26.98154,28.0855,30.97376,32.06,35.453,39.948,&                     
39.0983,40.08,44.9559,47.90,50.9415,51.996,54.9380,55.847,58.9332,58.70,63.546,65.38,&
69.72,72.59,74.9216,78.96,79.904,83.80,&                                            
85.4678,87.62,88.9059,91.22,92.9064,95.94,98.00,101.07,102.9055,106.4,107.868,112.41,&
114.82,118.69,121.75,127.60,126.9045,131.30,&                                       
132.9054,137.33,&                                                                   
138.9055,140.12,140.9077,144.24,145.00,150.4,151.96,&                               
157.25,158.9254,162.50,164.9304,167.26,168.9342,173.04,174.967,&                    
178.49,180.9479,183.85,186.207,190.2,192.22,195.09,196.9665,200.59,&                
204.37,207.2,208.9804,209.00,210.00,222.00,&                                        
223.00,226.0254,227.0278,232.0381,231.0359,238.029,237.0482,242.00,243.00,&         
247.07,247.07,251.07,252.08,257.09,258.09,259.10,266.12,&
267.12,268.12,269.12,270.13,269.13,277.15,282.16,282.16,286.17,&
286.18,290.19,290.19,293.20,294.21,295.21/
!
contains
!
!
!
! Calculates the distance between to atoms
real(kind=dp) function distan(na,nb)
integer :: na, nb
distan=sqrt((xcoord(na)-xcoord(nb))**2+(ycoord(na)-ycoord(nb))**2+(zcoord(na)-zcoord(nb))**2)
end function distan

! Calculates angle between na and nc: na -- nb -- nc 
real(kind=dp) function ang(na,nb,nc)
! na # of Atom a
! nb # of Atom b
! nc # of Atom c
integer :: na, nb, nc
! Angle A--B--C
! Vector #1 A-B  B->A
! Vector #2 C-B  B->C
ang=(((xcoord(na)-xcoord(nb))*(xcoord(nc)-xcoord(nb)))+&
    ((ycoord(na)-ycoord(nb))*(ycoord(nc)-ycoord(nb)))+&
    (((zcoord(na)-zcoord(nb))*(zcoord(nc)-zcoord(nb)))))
ang=ang/(distan(na,nb)*distan(nc,nb))
ang=acos(ang)*180/pi
end function ang

! Calculates dihedral angle between na and nd: na -- nb -- nc -- nd 
real(kind=dp) function dih(na,nb,nc,nd)
! na # of Atom a
! nb # of Atom b
! nc # of Atom c
! nd # of Atom d
integer :: na, nb, nc, nd
!real(kind=dp) :: a, b, c, d, e, f
real(kind=dp) :: abx, aby, abz, cbx, cby, cbz, bcx, bcy, bcz, dcx, dcy, dcz
real(kind=dp) ::  nx, ny, nz, mx, my, mz, ox, oy, oz 
! dihedral A--B--C--D
! Vector a: B->A
abx=xcoord(na)-xcoord(nb)
aby=ycoord(na)-ycoord(nb)
abz=zcoord(na)-zcoord(nb)
! Vector b: B->C
cbx=xcoord(nc)-xcoord(nb)
cby=ycoord(nc)-ycoord(nb)
cbz=zcoord(nc)-zcoord(nb)
! Vector b: C->B
bcx=-cbx
bcy=-cby
bcz=-cbz
! Vector c: C->D
dcx=xcoord(nd)-xcoord(nc)
dcy=ycoord(nd)-ycoord(nc)
dcz=zcoord(nd)-zcoord(nc)
! CODE OF NORMAL VECTOR N = B->A X B->C
nx=(aby*cbz)-(abz*cby)
ny=(abz*cbx)-(abx*cbz)
nz=(abx*cby)-(aby*cbx)
! CODE OF NORMAL VECTOR M = C->B X C->D
mx=(bcy*dcz)-(bcz*dcy)
my=(bcz*dcx)-(bcx*dcz)
mz=(bcx*dcy)-(bcy*dcx)
! CODE OF NORMAL VECTOR O = N X M
ox=(ny*mz)-(nz*my)
oy=(nz*mx)-(nx*mz)
oz=(nx*my)-(ny*mx)
! COS THETA = n.m / |n|.|m|
! n.m / |n|.|m|
dih=((nx*mx)+(ny*my)+(nz*mz))/((sqrt((nx*nx)+(ny*ny)+(nz*nz)))*(sqrt((mx*mx)+(my*my)+(mz*mz))))
! FROM RAD TO GRAUS
dih=acos(dih)*180/pi
! func sign(A,B) gives the A value with the sign of B
dih=dsign(dih,((ox*cbx)+(oy*cby)+(oz*cbz)))
end function dih

! gives the covalent radius 
real (kind=dp) function radius(name)                                     
character (len=2) :: name
integer :: i
do i=1,totalelem
  if(elem(i) == name) radius=radiuscov(i)
end do
end function radius

! gives the atomic number
integer function atonum(name)                                     
character (len=2) :: name
integer :: i
do i=1,totalelem
  if(elem(i) == name) atonum=i
end do
end function atonum

! gives the atomic mass
real(kind=dp) function atomass(name)                                     
character (len=2) :: name
integer :: i
do i=1,totalelem
  if(elem(i) == name) atomass=mass(i)
end do
end function atomass

!gives the mapper value for an atom 
real(kind=dp) function map(name)                                     
character (len=2) :: name
integer :: i
do i=1,totalelem
  if(elem(i) == name) map=mapper(i)
end do
end function map

! subrotina read_inp
subroutine read_inp()
character (len=2) :: name
real(kind=dp) :: val
integer :: i, j, m, c
!
! Setting False for input loop
debug=.False.; pcart=.False.; cbonds=.False.; cmass=.False.; molecule=.False.
psbond=.False.; pbonds=.False.; psang=.False.; pangs=.False.; psdih=.False.
pdih=.False.; full_cell=.False.; pcell=.False.; in_cell=.False.; within=.False.
psmolecule=.False.; in_cellcmass=.False.; pmolecule=.False.; neighb=.False.; aneighb=.False.
psolute=.False.; pnns=.False.; pstrajec=.False.; ptrajec=.False.; psaxsf=.False.
paxsf=.False.; linp=.False.; lcharge=.False.
! Setting null for input loop
n_neighbor=0; e_mapper=0; a_mapper=0; slt=1; totsolute=0; atomorig=1
ipsbond=0; ipsang=0; ipsdih=0; ipsmolecule=0; tstep=0; axsfstep=0; a_neighbor=0
!  Reading options
print '(" Reading options ... ")'
do
  read(1,FMT='(*(A))',END=101) keywords
  n=len_trim(keywords)
!!
  if(keywords(:n) == 'DEBUG')then
    print '("  Printing debug information in the output files")'
    debug=.True.
    pcart=.True.
!!
  else if(keywords(:n) == 'TWO BODY')then
    read(1,FMT='(i2)',ERR=129,IOSTAT=m,END=101) ntbody
    129 if(ntbody == 0)then
      print '("  Error reading the number of two body property pairs")'
    else
      allocate(itbody(ntbody,2),stat=m)
      allocate(ltbody(ntbody),stat=m)
      itbody(:,:)=0
      ltbody(:)=.False.
      print '("  Number of two body property pairs: ", i2)', ntbody
      do i=1,ntbody
        read(1,FMT=*,ERR=130,IOSTAT=m,END=101) itbody(i,1), itbody(i,2), ltbody(i)
        print '("   Two body property between atoms ",i4,",",i4," with ", l," value")',&
        & itbody(i,1), itbody(i,2), ltbody(i)
        130 if(itbody(i,1) == 0 .or. itbody(i,2) == 0)then
          print '("  Error Reading Atom indexes: ")'
          STOP
        end if
      end do
    end if
    tbody=.True.
    cbonds=.True.
!!
  else if(keywords(:n) == 'CENTER OF MASS')then
    print '("  Calculating molecular center of mass")'
    cmass=.True.
    molecule=.True.
    cbonds=.True.
    pcart=.True.
!!
  else if(keywords(:n) == 'PRINT BOND')then
    read(1,FMT='(i2)',ERR=112,IOSTAT=m,END=101) ipsbond
    112 if(ipsbond == 0)then
      print '("  All bond lengths will be printed")'
    else
      print '("  Selected", i4," Bond length(s) will be printed")', ipsbond
      allocate(sbond(ipsbond,2),stat=m)
      sbond(:,:)=0d0
      psbond=.True.
      do i=1,ipsbond
        read(1,FMT=*,ERR=113,IOSTAT=m,END=101)  sbond(i,1), sbond(i,2)
        print '("   Printing Bond length between atoms ", i4,",",i4)',  sbond(i,1), sbond(i,2)
        113 if(sbond(i,1) == 0 .or. sbond(i,2) == 0)then
          print '("  Error Reading Atom indexes: ")'
          STOP
        end if
      end do
    end if
    pcart=.True.
    cbonds=.True.
    pbonds=.True.
!! 
  else if(keywords(:n) == 'MAPPER')then
    read(1,FMT='(f5.3)',ERR=131,IOSTAT=m,END=101) tmp_mapper
    131 if(tmp_mapper == 0.00)then
      print '("  Error Reading the new mapper value: ")'
      mapper=1.100d0
      print '("   It continues with default value of ", f5.3)', 1.100
    else
      print '("  New mapper value found: ", f5.3)', tmp_mapper
      do i=1, totalelem
        if(mapper(i) /= 1.100d0) cycle
        mapper(i)=tmp_mapper
      end do
    end if
!!
  else if(keywords(:n) == 'ELEMENT MAPPER')then
    read(1,FMT='(i3)',ERR=105,IOSTAT=m,END=101) e_mapper
    105 if(e_mapper == 0)then
      print '("  Error Reading Number of mapper(s): ")'
      STOP
    end if
    print '("  Number of element mapper(s): ",i2)',  e_mapper
    do i=1,e_mapper
      read(1,FMT='(a2,f6.3)',ERR=106,IOSTAT=m,END=101) name, val
      106 if(val == 0.0d0 .or. name == '')then
         print '("  Error Reading Name and mapper value(s): ")'
         STOP
      end if
      print '("    New element mapper for ",a2," = ",f6.3)',  name, val
      do j=1,totalelem
        if(elem(j) == name) mapper(j)=val
      end do
    end do
!!
  else if(keywords(:n) == 'ATOM MAPPER')then
    read(1,FMT='(i3)',ERR=122,IOSTAT=m,END=101) a_mapper
    122 if(a_mapper == 0)then
      print '("  Error Reading Number of mapper(s): ")'
      STOP
    end if
    print '("  Number of atom mapper(s): ",i2)',  a_mapper
    allocate(iamap(a_mapper),stat=m)
    allocate(samap(a_mapper),stat=m)
        iamap(:)=0
        samap(:)=0.0d0
    do i=1,a_mapper
      read(1,FMT=*,ERR=123,IOSTAT=m,END=101) iamap(i), samap(i)
      123 if(iamap(i) == 0 .or. samap(i) == 0.0d0)then
        if(iamap(i) == 0) print '("  Error Reading Atom indexes: ")'
        if(samap(i) == 0.0d0) print '("  Error Reading Atom mappers: ")'
        STOP
      end if
      print '("    New mapper for atom ",i3," = ",f6.3)',  iamap(i), samap(i)
    end do
!!
  else if(keywords(:n) == 'ELEMENT CHARGE')then
    lcharge=.True.
    read(1,FMT='(i3)',ERR=120,IOSTAT=m,END=101) e_charge
    120 if(e_charge == 0)then
      print '("  Error Reading element charge(s): ")'
      STOP
    end if
    print '("  Number of element Charge(s): ",i2)',  e_charge
    do i=1,e_charge
      read(1,FMT='(a2,i3)',ERR=121,IOSTAT=m,END=101) name, c
      121 if(m /= 0 .or. name == '')then
         print '("  Error Reading Name and Charge value(s): ")'
         STOP
      end if
      print '("    Charge for element ",a2," = ",i3)',  name, c
      do j=1,totalelem
        if(elem(j) == name) charge(j)=c
      end do
    end do
!!
  else if(keywords(:n) == 'ATOM CHARGE')then
    lcharge=.True.
    read(1,FMT='(i3)',ERR=127,IOSTAT=m,END=101) a_charge
    127 if(a_charge == 0)then
      print '("  Error Reading atom charge(s): ")'
      STOP
    end if
    print '("  Number of atom charge(s): ",i2)',  a_charge
    allocate(iacharge(a_charge),stat=m)
    allocate(sacharge(a_charge),stat=m)
        iacharge(:)=0
        sacharge(:)=0
    do i=1,a_charge
      read(1,FMT=*,ERR=128,IOSTAT=m,END=101) iacharge(i), sacharge(i)
      128 print '("    New Charge for atom ",i3," = ",i3)',  iacharge(i), sacharge(i)
    end do
!!
  else if(keywords(:n) == 'PRINT ANGLE')then
    read(1,FMT='(i2)',ERR=114,IOSTAT=m,END=101) ipsang
    114 if(ipsang == 0)then
            print '("  All bond angles will be printed")'
    else
            print '( "  Selected", i4," Bond angles will be printed")', ipsang
            allocate(sang(ipsang,3),stat=m)
            sang(:,:)=0d0
            psang=.True.
            do i=1,ipsang
            read(1,FMT=*,ERR=115,IOSTAT=m,END=101)  sang(i,1), sang(i,2), sang(i,3)
            print '("   Printing Bond angle of", i4,",",i4,",",i4)',  sang(i,1), sang(i,2), sang(i,3)
          115 if(sang(i,1) == 0 .or. sang(i,2) == 0 .or. sang(i,3) == 0)then
              print '("  Error Reading Atom indexes: ")'
              STOP
            end if
        end do
    end if
    pcart=.True.
    cbonds=.True.
    pangs=.True.
!!
  else if(keywords(:n) == 'PRINT DIHEDRAL')then
    read(1,FMT='(i2)',ERR=116,IOSTAT=m,END=101) ipsdih
    116 if(ipsdih == 0)then
            print '("  All Dihedral angles will be printed")'
    else
            print '( "  Selected", i4," Dihedral angles will be printed")', ipsdih
            allocate(sdih(ipsdih,4),stat=m)
            sdih(:,:)=0d0
            psdih=.True.
            do i=1,ipsdih
            read(1,FMT=*,ERR=117,IOSTAT=m,END=101)  sdih(i,1), sdih(i,2), sdih(i,3), sdih(i,4)
            print '("   Printing Dihedral angle of", i4,",",i4,",",i4, ",",i4)',  sdih(i,1), sdih(i,2), sdih(i,3), sdih(i,4)
          117 if(sdih(i,1) == 0 .or. sdih(i,2) == 0 .or. sdih(i,3) == 0 .or. sdih(i,4) == 0)then
              print '("  Error Reading Atom indexes: ")'
              STOP
            end if
        end do
    end if
    pcart=.True.
    cbonds=.True.
    pdih=.True.
!!
  else if(keywords(:n) == 'PRINT CARTESIAN')then
    print '("  Cartesian coordinates will be printed")'
    pcart=.True.
!!
  else if(keywords(:n) == 'ORIGIN ATOM')then
    read(1,FMT='(i3)',ERR=110,IOSTAT=m,END=101) atomorig
    print '("  Atom will be located at center",i3)', atomorig
      110 if(atomorig == 0)then
           print '("  Error Reading origin atom: ")'
           STOP
    end if
!!
  else if(keywords(:n) == 'ORIGIN COORDINATES')then
    read(1,FMT=*,ERR=125,IOSTAT=m,END=101) orig_x, orig_y, orig_z
    print '("  Atom will be located at X, Y, Z coords",f12.8,2X,f12.8,2X,f12.8)', orig_x, orig_y, orig_z
    125 if(m /= 0)then
        print '("  Error Reading coordinates of origin atom: ")'
        STOP
    end if
!!
  else if(keywords(:n) == 'CELL')then
    print '("  Super cell will be generated")'
    read(1,FMT='(f6.3)',ERR=102,IOSTAT=m,END=101) cellparam
    full_cell=.True.
    print '("   Cell Parameter value: ",f6.3)',  cellparam
    102 if(full_cell.eqv..False..or.cellparam == 0.0d0)then
    print '("  Error Reading Cell Parameter value: ")'
    STOP
  end if
!!
  else if(keywords(:n) == 'PRINT SCELL')then
    if(full_cell)then
      print '("  Super Cell file will be printed")'
      pcell=.True.
  else
    print '("  Cannot print Super Cell file without generate it")'
    print '("  Please specify CELL parameter in the input file")'
    print '("  or specify SCELL after CELL parameter in the input file")'
    STOP
  end if
!!
  else if(keywords(:n) == 'INCELL')then
    print '("  Putting all atoms inside the cell")'
    read(1,FMT='(f6.3)',ERR=107,IOSTAT=m,END=101) cellparam
    in_cell=.True.
    pcart=.True.
    print '("  Found Cell Parameter value: ",f6.3)',  cellparam
    107 if(in_cell.eqv..False..or.cellparam == 0.0d0)then
        print '("  Error Reading Cell Parameter value: ")'
        STOP
    end if
!!
  else if(keywords(:n) == 'WITHIN')then
    print '("  Printing solvent within a radii cutoff from solute")'
    cbonds=.True.
    within=.True.
    read(1,FMT='(i3)',ERR=108,IOSTAT=m,END=101) totsolute
    print '("    Number of solute atoms",i3)', totsolute
      108 if(totsolute == 0)then
           print '("  Error Reading Number of total of solute atoms: ")'
           STOP
    end if
    read(1,FMT='(f6.3)',ERR=109,IOSTAT=m,END=101) rcut
      print '("    Printing solvent within a radii cutoff of",f6.3," Angstron(s)")', rcut
      109 if(rcut == 0.0d0)then
           print '("  Error Reading radius cut: ")'
           STOP
    end if
!!
  else if(keywords(:n) == 'MOLECULE')then
    print '("  Computing molecules")'
    cbonds=.True.
    molecule=.True.
    cmass=.True.
!!
  else if(keywords(:n) == 'PRINT MOLECULE')then
    print '("  Printing molecules")'
    read(1,FMT='(i4)',ERR=111,IOSTAT=m,END=101) ipsmolecule
    111 if(ipsmolecule == 0)then
            print '("  All molecules will be printed")'
          else
            print '("  Molecule ",i4," will be printed")', ipsmolecule
            psmolecule=.True.
    end if
    pmolecule=.True.
    cbonds=.True.
    molecule=.True.
    cmass=.True.
    pcart=.True.
!!
  else if(keywords(:n) == 'NEIGHBOR')then
    print '("  Neighbors will be computed ")'
    neighb=.True.
    cbonds=.True.
    molecule=.True.
    cmass=.True.
    read(1,FMT='(i3)',ERR=104,IOSTAT=m,END=101) n_neighbor
    print '("   Number of nearest neighbor(s) requested: ",i2)',  n_neighbor
    104 if(n_neighbor == 0)then
      print '("  All nearest neighbor(s) will be printed: ")'
    end if
!!
  else if(keywords(:n) == 'ATOM NEIGHBOR')then
    print '("  Atom neighbors will be computed ")'
    aneighb=.True.
    neighb=.True.
    cbonds=.True.
    molecule=.True.
    cmass=.True.
    read(1,FMT='(i3)',ERR=124,IOSTAT=m,END=101) a_neighbor
    print '("   Requested atom number for nearest neighbor(s): ",i2)',  a_neighbor
    124 if(a_neighbor == 0)then
      aneighb=.False.
      print '("  Warning: Impossible to read selected atom neighbor! ")'
      print '("  All nearest neighbor(s) for solute will be printed: ")'
    end if
!!
   else if(keywords(:n) == 'PRINT SOLUTE')then
     print '("  Solute geometry will be printed ")'
     psolute=.True.
     cbonds=.True.
     molecule=.True.
     pmolecule=.True.
     cmass=.True.
     pcart=.True.
!!
  else if(keywords(:n) == 'SOLUTE')then
    print '("  Setting up new solute: ")'
    read(1,FMT=*,ERR=103,IOSTAT=m,END=101) slt
    cbonds=.True.
    molecule=.True.
    cmass=.True.
    103 if(slt == 0)then
      print '("  Error Reading Solute Number: ")'
      STOP
    end if
    print '("    Number of solute: ",i5)',  slt
!!
  else if(keywords(:n) == 'PRINT NEIGHBOR')then
    print '("  Nearest Neighbor(s) file will be printed")'
    cbonds=.True.
    molecule=.True.
    cmass=.True.
    neighb=.True.
    pnns=.True.
!!
  else if(keywords(:n) == 'PRINT TRAJECTORY')then
    read(1,FMT='(i2)',ERR=118,IOSTAT=m,END=101) tstep
    118 if(tstep == 0)then
      print '("  Trajectory file will be printed")'
    else
      print '("  Trajectory file will be printed")'
      print '("    Trajectory coords will be printed in steps of ",i3)', tstep
      pstrajec=.True.
    end if
    cbonds=.True.
    molecule=.True.
    cmass=.True.
    neighb=.True.
    ptrajec=.True.
!!
  else if(keywords(:n) == 'PRINT AXSF')then
    read(1,FMT='(i2)',ERR=119,IOSTAT=m,END=101) axsfstep
    119 if(axsfstep == 0)then
      print '("  AXSF file will be printed")'
    else
      print '("  AXFS file will be printed")'
      print '("    AXFS coords will be printed in steps of ",i3)', axsfstep
      psaxsf=.True.
    end if
    if(cellparam == 0.0d0)then
        print '("  Cannot print AXSF file without cell parameter")'
        print '("  Please specify CELL parameter in the input file")'
        print '("  or specify PRINT AXSF after CELL parameter in the input file")'
        STOP
    end if
    cbonds=.True.
    molecule=.True.
    cmass=.True.
    neighb=.True.
    paxsf=.True.
!!
  else if(keywords(:n) == 'H TOPOLOGY')then
    print '("  Considering only Hydrogen with the shortest contact")'
    topos=.True.
!!
  else if(keywords(:n) == 'INCELL CMASS')then
    print '("  Putting all molecules inside the cell based on center of mass")'
    read(1,FMT='(f6.3)',ERR=126,IOSTAT=m,END=101) cellparam
    pcart=.True.
    cbonds=.True.
    molecule=.True.
    cmass=.True.
    in_cellcmass=.True.
    print '("  Found Cell Parameter value: ",f6.3)',  cellparam
    126 if(in_cellcmass.eqv..False..or.cellparam == 0.0d0)then
        print '("  Error Reading Cell Parameter value: ")'
        STOP
    end if
  end if
end do
! read(1,FMT='(A)',END=101)
101 rewind(1)
print '(" End of reading options ... ")'
end subroutine read_inp

! subrotina read_cart
subroutine read_cart()
integer :: j, m
! 
print '(" Calling read_cart subroutine")'
! Reading total number of atoms and saving in totatom variable
read(1,*) totatom
read(1,*) ! reading comment line
!
! Allocating matrices 
allocate(rsymb(totatom),stat=m)
if(m /= 0) print *, "allocating rsymb failed"
allocate(rxcoord(totatom),stat=m)
if(m /= 0) print *, "allocating rxcoord failed"
allocate(rycoord(totatom),stat=m)
if(m /= 0) print *, "allocating rycoord failed"
allocate(rzcoord(totatom),stat=m)
if(m /= 0) print *, "allocating rzcoord failed"
allocate(rmap(totatom),stat=m)
if(m /= 0) print *, "allocating rmap failed"
!
!
if(full_cell)then
  ! Allocating matrices with totatom*27 size duo to full_cell request
  allocate(symb(totatom*27),stat=m)
  if(m /= 0) print *, "allocating symb failed"
  allocate(xcoord(totatom*27),stat=m)
  if(m /= 0) print *, "allocating xcoord failed"
  allocate(ycoord(totatom*27),stat=m)
  if(m /= 0) print *, "allocating ycoord failed"
  allocate(zcoord(totatom*27),stat=m)
  if(m /= 0) print *, "allocating zcoord failed"
  allocate(amap(totatom*27),stat=m)
  if(m /= 0) print *, "allocating amap failed"
  if(lcharge)then
    allocate(acharge(totatom*27),stat=m)
    acharge=0
    if(m /= 0) print *, "allocating acharge failed"
  end if
else
  ! Allocating matrices with totatom size
  allocate(symb(totatom),stat=m)
  if(m /= 0) print *, "allocating symb failed"
  allocate(xcoord(totatom),stat=m)
  if(m /= 0) print *, "allocating xcoord failed"
  allocate(ycoord(totatom),stat=m)
  if(m /= 0) print *, "allocating ycoord failed"
  allocate(zcoord(totatom),stat=m)
  if(m /= 0) print *, "allocating zcoord failed"
  allocate(amap(totatom),stat=m)
  if(m /= 0) print *, "allocating amap failed"
  if(lcharge)then
    allocate(acharge(totatom),stat=m)
    if(m /= 0) print *, "allocating acharge failed"
    acharge=0
  end if
end if
! 
print '("  Number of Atoms found", i6)', totatom
!
! Setting rmpa to zero before copying sampa into it
rmap(:)=0.0d0
!
! reading symbols and cartersian coordinates from input file
do j=1, totatom
  read(1,*) rsymb(j), rxcoord(j), rycoord(j), rzcoord(j)
end do
!
! Copying read mappers to rmap
if(a_mapper /= 0)then
  do i=1,a_mapper
    if(iamap(i) <= totatom)then
      rmap(iamap(i))=samap(i)
    else
      print '("Requested index atom mapper > total # atoms")'
    end if
  end do
end if
! Filling the remaining rmap and filling the acharge
do j=1, totatom
  if(rmap(j) == 0.0d0) rmap(j)=map(rsymb(j))
  if(lcharge) acharge(j)=charge(atonum(rsymb(j)))
end do
! Copying read atom charge to acharge
if(a_charge /= 0 .and. lcharge)then
  do i=1,a_charge
    if(iacharge(i) <= totatom)then
      acharge(iacharge(i))=sacharge(i)  
    else
      print '("Requested index atom charge > total # atoms")'
    end if
  end do
end if
!
! reading and putting atom atomorig (default #1) in the center of the supercell (cellfactor=1.5d0) or not using cellfactor=0.5d0
if(full_cell) cellfactor=1.5d0
!
if(in_cell .or. in_cellcmass) cellfactor=0.5d0
!
!
do j=1, totatom
symb(j)=rsymb(j) ! Delta_i+Shift+DeltaCell
  xcoord(j)=rxcoord(j)-rxcoord(atomorig)+orig_x+(cellparam*cellfactor)
  ycoord(j)=rycoord(j)-rycoord(atomorig)+orig_y+(cellparam*cellfactor)
  zcoord(j)=rzcoord(j)-rzcoord(atomorig)+orig_z+(cellparam*cellfactor)
  amap(j)=rmap(j)
end do
!
if(debug)then
  write (4,*) 'DEBUG INFO FROM READ_CART SUBROUTINE'
  write (4,*) cellfactor, 'cellfactor from read_cart'
  do j=1, totatom
    write (4,'(A2,2X,f12.8,2X,f12.8,2X,f12.8,2X,f6.4,2X,f8.4,2X,f8.4,2X,f8.4,A8)') symb(j), xcoord(j), ycoord(j), zcoord(j), &
    radius(symb(j)), map(symb(j)), rmap(j), amap(j), ' Leitura'
    if(.not.lcharge) write (4,'(A2,2X,f12.8,2X,f12.8,2X,f12.8,2X,f6.4,2X,f6.4,2X,f6.4)') symb(j), xcoord(j), ycoord(j), zcoord(j), &
    radius(symb(j)), amap(j), map(symb(j))
    if(lcharge) write (4,'(A2,2X,f12.8,2X,f12.8,2X,f12.8,2X,f6.4,2X,f6.4,2X,f6.4,2X,I2)') symb(j), xcoord(j), ycoord(j), &
    zcoord(j), radius(symb(j)), amap(j), map(symb(j)), acharge(j)
  end do
  write (4,*) 'END DEBUG INFO FROM READ_CART SUBROUTINE'
end if
!
end subroutine read_cart

! subrotina print_cart
subroutine print_cart()
integer :: j
do j=1, totatom
 if(lcharge) cartcharge=cartcharge+acharge(j)
 cartmass=cartmass+atomass(symb(j))
end do
print '(" Calling print_cart subroutine")'
write (2,'(i5)') totatom
!write (2,*)
if(lcharge) then
  write(2,'(" Mass ", f8.3, " Charge", i4)') cartmass, cartcharge
else
  write(2,'("mass ", f8.3)') cartmass
end if

!
do j=1, totatom
  write (2,'(A2,2X,f12.8,2X,f12.8,2X,f12.8)') symb(j), xcoord(j), ycoord(j), zcoord(j)
end do
!
end subroutine print_cart

! subroutine incell
subroutine incell()
integer :: i
print '(" Calling incell subroutine")'
!
do i=1, totatom
! Coords larger than the cell parameter
  if(xcoord(i) > cellparam) xcoord(i)=xcoord(i)-(int(xcoord(i)/cellparam)*cellparam)
  if(ycoord(i) > cellparam) ycoord(i)=ycoord(i)-(int(ycoord(i)/cellparam)*cellparam)
  if(zcoord(i) > cellparam) zcoord(i)=zcoord(i)-(int(zcoord(i)/cellparam)*cellparam)
! Coords smaller than zero
  if(xcoord(i) < 0.D0) xcoord(i)=cellparam-(abs(xcoord(i))-int(abs(xcoord(i))/cellparam)*cellparam)
  if(ycoord(i) < 0.D0) ycoord(i)=cellparam-(abs(ycoord(i))-int(abs(ycoord(i))/cellparam)*cellparam)
  if(zcoord(i) < 0.D0) zcoord(i)=cellparam-(abs(zcoord(i))-int(abs(zcoord(i))/cellparam)*cellparam)
  if(full_cell)then
    xcoord(i)=xcoord(i)+(cellparam)
    ycoord(i)=ycoord(i)+(cellparam)
    zcoord(i)=zcoord(i)+(cellparam)
  end if
end do
!
end subroutine incell

! subroutine scell
subroutine scell()
integer i, l, m, n, j
cellfactor=1.5d0
print '(" Calling scell subroutine")'
if(debug)then
  write(4,*) 'DEBUG INFO FROM SCELL SUBROUTINE'
  write(4,*) 'From scell: '
  write(4,'("   Cell parameter is ",f6.3)')  cellparam
  write(4,*) totatom, 'totatom before changing  scell'
  write(4,*) size(symb), size(xcoord)
  write(4,*) cellparam, 'cellparam from scell'
end if
! using the original value of totatom read from input file. in this way no big changes are needed
! in read_cart and scell subroutines 
j=totatom+1
do l=-1,1
  do m=-1,1
    do n=-1,1
      if(l /= 0.or.m /= 0.or.n /= 0)then
          if(debug) write(4,*) l, m, n, 'Indexes from super cell'
          do i=1, totatom
            symb(j)=symb(i)
            amap(j)=amap(i)
            if(lcharge) acharge(j)=acharge(i)
            xcoord(j)=xcoord(i)-xcoord(atomorig)+((cellfactor+l)*cellparam)+orig_x
            ycoord(j)=ycoord(i)-ycoord(atomorig)+((cellfactor+m)*cellparam)+orig_y
            zcoord(j)=zcoord(i)-zcoord(atomorig)+((cellfactor+n)*cellparam)+orig_z
            !if(debug) write (4,'(A2,2X,f12.8,2X,f12.8,2X,f12.8,2X,f12.8,2X,f12.8,2X,f8.4)') symb(j),xcoord(atomorig),xcoord(atomorig), &
            !& xcoord(j),ycoord(j),zcoord(j), amap(j)
            j=j+1
          end do
      end if
    end do  
  end do  
end do
! setting totatom 27 x the original value for the next subroutines 
totatom=totatom*27
if(debug)then
  write(4,*) totatom, 'totatom from scell 27 times more'
  write (4,*) 'END DEBUG INFO FROM SCELL SUBROUTINE'
end if
end subroutine scell

! subroutine bonds
subroutine bonds()
print '(" Calling bonds subroutine")'
!
allocate(bond(totatom,totatom),stat=m)
if(m /= 0) print *, "allocating bond failed"
allocate(lbond(totatom,totatom),stat=m)
if(m /= 0) print *, "allocating lbond failed"
!
! Setting null value for bond
bond(:,:)=0.0d0
! Setting FALSE value for lbond
lbond(:,:)=.FALSE.
! Setting null value for tbond
tbond=0
! Looking for atom connections
!$OMP PARALLEL SHARED(bond,lbond,amap,totatom) PRIVATE(j,k)
!$OMP DO
do j=1, totatom-1
  do k=j+1, totatom
    bond(j,k)=distan(j,k)
    bond(k,j)=distan(j,k)
    if(((radius(symb(j))*amap(j))+(radius(symb(k))*amap(k))) >= distan(j,k))then
      lbond(j,k)=.TRUE.
      lbond(k,j)=.TRUE.
    end if
  end do
end do
!$OMP END DO
!$OMP END PARALLEL
!
tbond=count(lbond)/2
!
if(tbody)then
    print '("  Changing The Number of Bond(s) Found by Two Body Property")'
  do j=1,ntbody
    if(itbody(j,1).gt.totatom.or.itbody(j,2).gt.totatom)then
      print '("   Atom Index is larger than the number of atoms")'
      STOP
    end if  
    if(lbond(itbody(j,1),itbody(j,2)).eqv.ltbody(j))then
      print '("   Warning: Atoms", i4, " and", i4, " already is equal to ", L)',&
      & itbody(j,1), itbody(j,2), lbond(itbody(j,1),itbody(j,2))
      cycle
    else
      !if(ltbody(j)) tbond=tbond+1
      !if(.not.ltbody(j)) tbond=tbond-1
      lbond(itbody(j,1),itbody(j,2))=ltbody(j)
      lbond(itbody(j,2),itbody(j,1))=ltbody(j)
    end if
  end do
  tbond=count(lbond)/2
  print '("   Number of Bond(s) Found Changed by Two Body Property: ",i4)', tbond
else
  print '("  Number of Bond(s) Found: ",i4)', tbond
end if
!
if(debug)then
  write(4,*) 'DEBUG INFO FROM BOND SUBROUTINE'
  write(4,*) 'TOTAL NUMBER OF BONDS', tbond
  write(4,*) 'LOOP INDEXES OF LBOND AND BOND - SQUARE SIMMETRIC MATRICES'
  do j=1, totatom-1
    do k=j+1, totatom
      write(4,*) j, k, lbond(j,k), lbond(k,j), bond(j,k), bond(k,j)
    end do
  end do
  write(4,*) 'END DEBUG INFO FROM BOND SUBROUTINE'
end if
!
end subroutine bonds

! subroutine topology
subroutine topology()
print '(" Calling topology subroutine")'
print '("  Fixing H atoms with two or more bonds")'
!
! Looking for H connections
!$OMP PARALLEL SHARED(bond,lbond,symb,tbond) PRIVATE(j)
!$OMP DO 
do j=1, totatom
  if(symb(j).ne.'H ') cycle
    if((minloc(bond(:,j),DIM=1,mask=lbond(:,j))).ne.(maxloc(bond(:,j),DIM=1,mask=lbond(:,j))))then
      lbond(j,(maxloc(bond(:,j),DIM=1,mask=lbond(:,j))))=.False.
      lbond((maxloc(bond(:,j),DIM=1,mask=lbond(:,j))),j)=.False.
      !tbond=tbond-1
      if(debug)then
        write(4,*) 'DEBUG INFO FROM TOPOLOGY SUBROUTINE'
        write(4,*) 'lbond', lbond(j,(maxloc(bond(:,j),DIM=1,mask=lbond(:,j)))), lbond((maxloc(bond(:,j),DIM=1,mask=lbond(:,j))),j)
        tbond=count(lbond)/2
        write(4,*) 'tbond', tbond 
        write(4,*) 'Miniloc', minloc(bond(:,j),DIM=1,mask=lbond(:,j))
        write(4,*) 'Minival', minval(bond(:,j),DIM=1,mask=lbond(:,j))
        write(4,*) 'Maxloc', maxloc(bond(:,j),DIM=1,mask=lbond(:,j))
        write(4,*) 'Maxval', maxval(bond(:,j),DIM=1,mask=lbond(:,j))
        write(4,*) 'END DEBUG INFO FROM TOPOLOGY SUBROUTINE'
      end if
    end if
end do
!$OMP END DO
!$OMP END PARALLEL
!
tbond=count(lbond)/2
print '("  New Number of Bond(s) Found: ",i4)', tbond
!
end subroutine topology

! subroutine print_bonds
subroutine print_bonds()
print '(" Calling print_bonds subroutine")'
!
if(psbond.eqv..False.)then
!
  write(2,*)
  write(2,'(" Chemical Bond(s):", i4)') tbond
  do j=1, totatom-1
    do k=j+1, totatom
      if(lbond(k,j)) write(2,'(A,i4,2X,A,i4,2X," = ", f9.5)') &
      &symb(j), j, symb(k), k, bond(k,j)
    end do
  end do
end if
!
if(psbond)then
  if(ipsbond == 0)then
    print *, ' ipsbond cannot be zero '
    STOP
  end if
  write(2,*)
  write(2,'(" Selected Chemical Bond(s):", i4)') ipsbond
  do k=1, ipsbond
      write(2,'(A,i4,2X,A,i4,2X," = ", f9.5,2X,L)') &
      &symb(sbond(k,1)), sbond(k,1), symb(sbond(k,2)), sbond(k,2), bond(sbond(k,1),sbond(k,2)), lbond(sbond(k,1),sbond(k,2))
  end do
end if
end subroutine print_bonds

! subroutine print_angle
subroutine print_angle()
integer :: j, k, l
print '(" Calling print_angle subroutine")'
!
if(psang)then
  if(ipsang == 0)then
    print *, ' ipsang cannot be zero '
    STOP
  end if
  write(2,'(" Selected Bond Angles: ",i4)') ipsang
  do j=1,ipsang
    write (2,'(A,i4,2X,A,i4,2X,A,i4,2X," = ", f7.3)') symb(sang(j,1)), sang(j,1), &
    &symb(sang(j,2)), sang(j,2), symb(sang(j,3)), sang(j,3), &
    &ang(sang(j,1),sang(j,2),sang(j,3))
  end do
! 
else
!
write(2,'(" Angles:")')
  if(full_cell)then
  !$OMP PARALLEL SHARED(lbond) PRIVATE(j,k,l)
    do j=1, totatom/27
      !$OMP DO
      do k=1, totatom/27
        do l=1, totatom/27
          if(lbond(k,j).and.lbond(l,k).and.j /= l) write (2,'(A,i4,2X,A,i4,2X,A,i4,2X," = ", f7.3)') &
          &symb(j), j, symb(k), k, symb(l), l, ang(j,k,l)
        end do
      end do
      !$OMP END DO
    end do
  !$OMP END PARALLEL
  else
  !$OMP PARALLEL SHARED(lbond) PRIVATE(j,k,l)
    do j=1, totatom
      !$OMP DO
      do k=1, totatom
        do l=1, totatom
          if(lbond(k,j).and.lbond(l,k).and.j /= l) write (2,'(A,i4,2X,A,i4,2X,A,i4,2X," = ", f7.3)') &
          &symb(j), j, symb(k), k, symb(l), l, ang(j,k,l)
        end do
      end do
      !$OMP END DO
    end do
  !$OMP END PARALLEL
  end if
end if
end subroutine print_angle

! subroutine print_dihedral
subroutine print_dihedral()
integer :: j, k, l, m
print '(" Calling print_dihedral subroutine")'

if(psdih)then
  write(2,'(" Selected Dihedral Angles: ",i4)') ipsdih
  do j=1,ipsdih
    write (2,'(A,i4,2X,A,i4,2X,A,i4,2X,A,i4,2X," = ", f8.3)') &
    &symb(sdih(j,1)), sdih(j,1), symb(sdih(j,2)), sdih(j,2), &
    &symb(sdih(j,3)), sdih(j,3), symb(sdih(j,4)), sdih(j,4), &
    &dih(sdih(j,1), sdih(j,2), sdih(j,3), sdih(j,4))
  end do
else
  write(2,'(" Dihedral Angles:")')
!
  if(full_cell)then
  !$OMP PARALLEL SHARED(lbond) PRIVATE(j,k,l,m)
    do j=1, totatom/27
      do k=1, totatom/27
        !$OMP DO
        do l=1, totatom/27
          do m=1, totatom/27
            if(lbond(k,j).and.lbond(l,k).and.lbond(m,l).and.j /= l.and.k /= m.and.j /= m) &
            & write (2,'(A,i4,2X,A,i4,2X,A,i4,2X,A,i4,2X," = ", f8.3)') &
            &symb(j), j, symb(k), k, symb(l), l, symb(m), m, dih(j,k,l,m)
          end do
        end do
        !$OMP END DO
      end do
    end do
  !$OMP END PARALLEL
  else
  !$OMP PARALLEL SHARED(lbond) PRIVATE(j,k,l,m)
    do j=1, totatom
      do k=1, totatom
        !$OMP DO
        do l=1, totatom
          do m=1, totatom
            if(lbond(k,j).and.lbond(l,k).and.lbond(m,l).and.j /= l.and.k /= m.and.j /= m) &
            & write (2,'(A,i4,2X,A,i4,2X,A,i4,2X,A,i4,2X," = ", f8.3)') &
            &symb(j), j, symb(k), k, symb(l), l, symb(m), m, dih(j,k,l,m)
          end do
        end do
        !$OMP END DO
      end do
    end do
  !$OMP END PARALLEL
  end if
end if
end subroutine print_dihedral


! subroutine molecules
subroutine molecules()
integer :: i
print '(" Calling molecules subroutine")'
allocate(molec(totatom),stat=m)
if(m /= 0) print *, "allocating molec failed"
! Setting null value for counter
i=0
! Setting null value for molec
molec=0
! Setting null value for tmolec
tmolec=0
!
if(debug)then
  write(4,*) 'DEBUG INFO FROM MOLECULE SUBROUTINE'
end if  
!
! Calculating the total number of molecules
!$OMP PARALLEL SHARED(lbond,molec,tmolec) PRIVATE(j,k)
!$OMP DO ORDERED
do j=1, totatom-1
        !$OMP ORDERED
  do k=j+1, totatom
        if(.not.lbond(k,j)) cycle
      if(debug)then
        write(4,*) j, k, lbond(k,j), molec(j), molec(k)
        !write(4,*) ' FUNC COUNTS T', count(lbond(:,j)), count(lbond(j,:)), count(lbond)
      end if
    !if(lbond(k,j))then
      if(molec(j) == 0.and.molec(k) == 0)then
        if(debug) write(4,*) 'Zero and Zero'
        tmolec=tmolec+1
        molec(j)=tmolec
        molec(k)=tmolec
      else if(molec(j) /= 0.and.molec(k) == 0)then
        if(debug) write(4,*) 'Non Zero and Zero'
        molec(k)=molec(j)
      else if(molec(j) == 0.and.molec(k) /= 0)then
        if(debug) write(4,*) 'Zero and  Non Zero'
        molec(j)=molec(k)
      !!!!!!
      else if(molec(j) /= 0.and.molec(k) /= 0.and.molec(j) /= molec(k))then
        if(debug) write(4,*) 'NON Zero and NON Zero', molec(j), molec(k), tmolec, 'before'
        if(molec(j) < molec(k))then
          if(debug) write(4,*) 'molec(j) < molec(k)', molec(j), molec(k)
          do l=1, totatom
            if(l == k) cycle
            !print *, l , 'L'
            if(molec(l) == molec(k))then
              if(debug) write(4,*)  'molec(l) == molec(k)'
              if(debug) write(4,*) l, 'L', molec(j), 'molecJ', molec(k), 'molecK', molec(l), 'molecL'
              molec(l)=molec(j)
              if(debug) write(4,*) l, 'L', molec(j), 'molecJ', molec(k), 'molecK', molec(l), 'molecL*'
            else if(molec(l) > molec(k))then
              if(debug) write(4,*) 'molec(l) > molec(k)'
              if(debug) write(4,*)  l, 'L', molec(j), 'molecJ', molec(k), 'molecK', molec(l), 'molecL'
              molec(l)=molec(l)-1
              if(debug) write(4,*) l, 'L', molec(j), 'molecJ', molec(k), 'molecK', molec(l), 'molecL*'
            end if
          end do
          molec(k)=molec(j)!; print *, 'mudei', molec(k), 'molecK'
      !!!!!!
        else if(molec(j) > molec(k))then
          if(debug) write(4,*) 'molec(j) > molec(k)', molec(j), molec(k)
          do l=1, totatom
            if(l == j) cycle
            !print *, l , 'L', molec(j), 'molecJ'
            if(molec(l) == molec(j))then
              if(debug) write(4,*)  'molec(l) == molec(j)'
              if(debug) write(4,*)  l, 'L', molec(j), 'molecJ', molec(k), 'molecK', molec(l), 'molecL'
              molec(l)=molec(k)
              if(debug) write(4,*)  l, 'L', molec(j), 'molecJ', molec(k), 'molecK', molec(l), 'molecL*'
            else if(molec(l) > molec(j))then
              if(debug) write(4,*) 'molec(l) > molec(j)'
              if(debug) write(4,*) l, 'L', molec(j), 'molecJ', molec(k), 'molecK', molec(l), 'molecL'
              molec(l)=molec(l)-1
              if(debug) write(4,*)  l, 'L', molec(j), 'molecJ', molec(k), 'molecK', molec(l), 'molecL*'
            end if
          end do
          molec(j)=molec(k)!; print *, 'mudei', molec(j), 'molecJ'
        end if
        tmolec=tmolec-1
        if(debug) write(4,*)  molec(j), molec(k), tmolec, 'depois'
      end if
  end do
  ! Verifying if there is no molec(j)=0 left from 1 to totatom-1 range
  if(molec(j) == 0)then
      tmolec=tmolec+1
      molec(j)=tmolec
      i=i+1
      if(debug) write(4,*) molec(j), j, i, tmolec, 'from molecules'
  end if
       !$OMP END ORDERED
end do
!$OMP END DO  
!$OMP END PARALLEL
! Veryfing the last one atom is in the molec list
if(molec(totatom) == 0)then
  tmolec=tmolec+1
  molec(totatom)=tmolec
  i=i+1
end if
!
print '("  Number of Molecule(s): ",i4)', tmolec
!
if(i == 0)then
  print '("  No Single species found")'
else if(i == 1)then
  print '("  One Single specie found")'
else if(i > 1)then
  print '("  Number of Single species found: ",i4)', i
end if
!
allocate(imolec(tmolec),stat=m)
if(m /= 0) print *, "allocating imolec failed"
imolec=0
!
allocate(molecmass(tmolec),stat=m)
if(m /= 0) print *, "allocating molecmass failed"
molecmass=0.0d0

if(lcharge)then
  allocate(moleccharge(tmolec),stat=m)
  if(m /= 0) print *, "allocating moleccharge failed"
  moleccharge=0
end if

! Calculating the number of atoms in each molecule, its mass and charge
do j=1, tmolec
  do k=1, totatom
    if(j == molec(k))then
      if(debug) write (4,'(A2,2X,f12.8,2X,f12.8,2X,f12.8,2X,I5)') &
      &symb(k),xcoord(k),ycoord(k),zcoord(k),molec(k)
      molecmass(j)=molecmass(j)+atomass(symb(k))
      imolec(j)=imolec(j)+1
      if(lcharge) moleccharge(j)=moleccharge(j)+acharge(k)
    end if
  end do
  if(debug.and.lcharge) write (4,*) '# atoms ', imolec(j), 'MASS ', molecmass(j), 'CHARGE ', moleccharge(j)
end do  
!
if(debug)then
    write(4,*) 'END DEBUG INFO FROM MOLECULE SUBROUTINE'
end if  
!
end subroutine molecules

! subroutine centermass
subroutine centermass()
print '(" Calling center of mass subroutine")'
!
allocate(cmassx(tmolec),stat=m)
if(m /= 0) print *, "allocating cmassx failed"
cmassx=0.0d0
allocate(cmassy(tmolec),stat=m)
if(m /= 0) print *, "allocating cmassy failed"
cmassy=0.0d0
allocate(cmassz(tmolec),stat=m)
if(m /= 0) print *, "allocating cmassz failed"
cmassz=0.0d0

! Calculating center of mass for each molecule
do i=1, tmolec
  do k=1, totatom
    if(i == molec(k))then
      cmassx(i)=cmassx(i)+(atomass(symb(k))*xcoord(k))
      cmassy(i)=cmassy(i)+(atomass(symb(k))*ycoord(k))
      cmassz(i)=cmassz(i)+(atomass(symb(k))*zcoord(k))
    end if
  end do
  cmassx(i)=cmassx(i)/molecmass(i)
  cmassy(i)=cmassy(i)/molecmass(i)
  cmassz(i)=cmassz(i)/molecmass(i)
end do  
!
end subroutine centermass

! subroutine incelcmass
subroutine incellcmass()
integer :: i, j
print '(" Calling incellcmass subroutine")'
!
! multiplying cellparam by 3 due to supercell request 
if(full_cell) cellparam=cellparam*3
!
do i=1, tmolec
! CM Coords larger than the cell parameter
  if(cmassx(i) > cellparam)then
    do j=1, totatom
      if(molec(j) == i)then
        xcoord(j)=xcoord(j)-cmassx(i)
        xcoord(j)=(cmassx(i)-(int(cmassx(i)/cellparam)*cellparam))+xcoord(j)
      end if
    end do
    cmassx(i)=cmassx(i)-(int(cmassx(i)/cellparam)*cellparam)
  end if
  if(cmassy(i) > cellparam)then
    do j=1, totatom
      if(molec(j) == i)then
        ycoord(j)=ycoord(j)-cmassy(i)
        ycoord(j)=(cmassy(i)-(int(cmassy(i)/cellparam)*cellparam))+ycoord(j)
      end if
    end do
    cmassy(i)=cmassy(i)-(int(cmassy(i)/cellparam)*cellparam)
  end if
  if(cmassz(i) > cellparam)then
    do j=1, totatom
      if(molec(j) == i)then
        zcoord(j)=zcoord(j)-cmassz(i)
        zcoord(j)=(cmassz(i)-(int(cmassz(i)/cellparam)*cellparam))+zcoord(j)
      end if
    end do  
    cmassz(i)=cmassz(i)-(int(cmassz(i)/cellparam)*cellparam)
  end if
! CM Coords smaller than zero
  if(cmassx(i) < 0.D0)then
    do j=1, totatom
      if(molec(j) == i)then
        xcoord(j)=xcoord(j)-cmassx(i)
        xcoord(j)=(cellparam-(abs(cmassx(i))-int(abs(cmassx(i))/cellparam)*cellparam))+xcoord(j)
      end if
    end do
    cmassx(i)=cellparam-(abs(cmassx(i))-int(abs(cmassx(i))/cellparam)*cellparam)
  end if
  if(cmassy(i) < 0.D0)then
    do j=1, totatom
      if(molec(j) == i)then
        ycoord(j)=ycoord(j)-cmassy(i)
        ycoord(j)=(cellparam-(abs(cmassy(i))-int(abs(cmassy(i))/cellparam)*cellparam))+ycoord(j)
      end if
    end do
    cmassy(i)=cellparam-(abs(cmassy(i))-int(abs(cmassy(i))/cellparam)*cellparam)
  end if
  if(cmassz(i) < 0.D0)then
    do j=1, totatom
      if(molec(j) == i)then
        zcoord(j)=zcoord(j)-cmassz(i)
        zcoord(j)=(cellparam-(abs(cmassz(i))-int(abs(cmassz(i))/cellparam)*cellparam))+zcoord(j)
      end if
    end do  
    cmassz(i)=cellparam-(abs(cmassz(i))-int(abs(cmassz(i))/cellparam)*cellparam)
  end if
end do
!
!  backing the cellparam value to the original for the remaining routines
if(full_cell) cellparam=cellparam/3
end subroutine incellcmass



! subroutine print_molecules
subroutine print_molecules()
print '(" Calling print_molecules subroutine")'
! Printing molecular coordinates
  write(2,*)
  write(2,'(" ********** SUBROUTINE MOLECULES ********** ")')
  write(2,*)
!
if(debug)then
  write(4,*) 'DEBUG INFO FROM PRINT MOLECULES SUBROUTINE'
end if  
!
! Priting selected molecules
if(psmolecule)then
  if(ipsmolecule > tmolec)then
    print '("  Molecule number is larger than the total molecules found ")'
    print '("  Printing molecule number 1 ")'
    ipsmolecule=1
  end if
  if(lcharge)then
    write(2,'(" Printing selected molecule #",i5, " with ", i4, " atom(s), mass ", f8.3, " and charge ", i4)') ipsmolecule,&
    imolec(ipsmolecule), molecmass(ipsmolecule), moleccharge(ipsmolecule)
  else
    write(2,'(" Printing selected molecule #",i5, " with ", i4, " atom(s), mass ", f8.3)') ipsmolecule,&
    imolec(ipsmolecule), molecmass(ipsmolecule)
  end if  
  write(2,'("   Center of mass: ",f12.8,1x,f12.8,1x,f12.8)') cmassx(ipsmolecule), cmassy(ipsmolecule), cmassz(ipsmolecule)  
  write(2,*)
else
  write(2,'(" Printing",i5," molecule(s)")') tmolec
  write(2,*)
end if
!
! Printing all molecules in unit 2
!
do j=1, tmolec
  if(psmolecule.eqv..False.)then
    if(lcharge)then
      write(2,'(" Printing molecule #",i5, " with ", i4, " atom(s), mass ", f8.3, " and charge", i4)')&
      j, imolec(j), molecmass(j), moleccharge(j)
      write(2,'("   Center of mass: ",f12.8,1x,f12.8,1x,f12.8)') cmassx(j), cmassy(j), cmassz(j)
    else
      write(2,'(" Printing molecule #",i5, " with ", i4, " atom(s), mass ", f8.3)') j, imolec(j), molecmass(j)
      write(2,'("   Center of mass: ",f12.8,1x,f12.8,1x,f12.8)') cmassx(j), cmassy(j), cmassz(j)
    end if
  end if
!
  do k=1, totatom
    if(j == molec(k))then
      if(psmolecule.eqv..False.) write (2,'(A2,2X,f12.8,2X,f12.8,2X,f12.8,2X)') &
      &symb(k),xcoord(k),ycoord(k),zcoord(k)
      ! Priting just the requested molecule
      if(psmolecule.and.molec(k) == ipsmolecule) write (2,'(A2,2X,f12.8,2X,f12.8,2X,f12.8,2X)') &
      &symb(k),xcoord(k),ycoord(k),zcoord(k)
      ! Debug
      if(debug) write (4,'(A2,2X,f12.8,2X,f12.8,2X,f12.8,2X,I5)') &
      &symb(k),xcoord(k),ycoord(k),zcoord(k),molec(k)
    end if
  end do
end do  
!
if(debug)then
    write(4,*) 'END DEBUG INFO FROM PRINT MOLECULES SUBROUTINE'
end if  
end subroutine print_molecules

! subroutine print_solute
subroutine print_solute()
print '(" Calling print_solute subroutine")'
!
! Priting solute molecule in inpu_solute.xyz
pinp2=arg(1:index(arg,'.')-1)//'_solute.xyz'
print *,'  New info file... ', pinp2
open(unit=3,file=pinp2,status='replace')
write(3,'(i5)') imolec(slt)
write(3,'(a15)') 'Solute Molecule' 
do k=1, totatom
  if(slt == molec(k))then
    write (3,'(A2,2X,f12.8,2X,f12.8,2X,f12.8,2X)')&
    & symb(k),xcoord(k),ycoord(k),zcoord(k)
  end if
end do
!
close(3,status='keep')
end subroutine print_solute

! subroutine neighbor
subroutine neighbor()
print '(" Calling neighbor subroutine")'
! Allocating Matrices with tmolec size
! Allocating molec_bond matrix
allocate(molec_bond(tmolec,tmolec),stat=m)
if(m /= 0) print *, "allocating molec_bond failed"
! Allocating molec_map matrix
allocate(molec_map(tmolec,tmolec),stat=m)
if(m /= 0) print *, "allocating molec_map failed"
! Allocating molec_mask matrix
allocate(molec_mask(tmolec,tmolec),stat=m)
if(m /= 0) print *, "allocating molec_map failed"
!
! Setting null value for molec_bond matrix elements
molec_bond=0.d0
molec_map=0
! J is setting up for molecule one for setting the n-1 neighbors
j=1
!
if(debug)then
  write(4,*) 'DEBUG INFO FROM NEIGHBOR SUBROUTINE'
  if(aneighb) write(4,*) slt, 'slt from neighbor', molec(a_neighbor), 'molec(a_neighbor)', a_neighbor, 'a_neighbor'
end if  
!
! Verifying if atom neighbor is in solute, otherwise will continue the search based on all atoms of solute
if(aneighb)then
  if (molec(a_neighbor) /= slt)then
    print '("  Warning: Selected atom does not belong to solute molecules! ")' 
    print '("           We will continue with NNS of solute atoms!  ")'
    aneighb=.False.
  end if
end if
! Verifying the smaller contact between molecule one and the n-1 remaining
!$OMP PARALLEL SHARED(molec,molec_bond,bond,lbond) PRIVATE(k,m,n)
!$OMP DO
!
! Code to find the smallest contact of n-1 neighbors of molecule slt
do k=1, tmolec
  if(k == slt) cycle
  do m=1, totatom-1
    do n=m+1, totatom
      if(aneighb .and. (m /= a_neighbor .and. n /= a_neighbor)) cycle
      if(lbond(n,m)) cycle
      if(molec(m) == slt.and.molec(n) == k.or.molec(m) == k.and.molec(n) == slt)then
        if(molec_bond(k,slt) == 0.d0)then
          molec_bond(slt,k)=bond(n,m)
          molec_bond(k,slt)=bond(n,m)
          !if(debug) write(4,*) k, m, n, molec_bond(slt,k), molec_bond(k, slt), 'from neighbor', '*' 
        else if(bond(n,m) < molec_bond(k,slt))then
          molec_bond(k,slt)=bond(n,m)
          molec_bond(slt,k)=bond(n,m)
          !if(debug) write(4,*) k, m, n, molec_bond(slt,k), molec_bond(k, slt), 'from neighbor', '***'
        end if
      end if
    end do
  end do
end do
!$OMP END DO
!$OMP END PARALLEL
!
if(debug)then
  do k=1, tmolec
    if(molec_bond(k,slt) == 0.0d0) write(4,*) k, molec_bond(k,slt)
  end do
end if
!
molec_mask=.TRUE.
!
if(debug) write(4,*) 'molec_mask shows all T elements: ', molec_mask(:,slt)
!
do k=1, tmolec
  if(debug) write(4,*) minloc(molec_bond(:,slt),DIM=1,mask=molec_mask(:,slt))
  if(debug) write(4,*) minval(molec_bond(:,slt),DIM=1,mask=molec_mask(:,slt))
  last_pos=minloc(molec_bond(:,slt),DIM=1,mask=molec_mask(:,slt))
  molec_mask(last_pos,slt)=.FALSE.
  molec_mask(slt,last_pos)=.FALSE.
  molec_map(k,slt)=last_pos
  molec_map(slt,k)=last_pos
  if(debug) write(4,*) k, slt, last_pos, molec_bond(last_pos,slt),&
  molec_map(k,slt), molec_map(slt,k),molec_mask(last_pos,slt), molec_mask(slt,last_pos)
end do
!
if(debug) write(4,*) 'END DEBUG INFO FROM NEIGHBOR SUBROUTINE'
end subroutine neighbor

! subroutine print_trajec
subroutine print_trajectory()
integer :: i, j, k, l, m, n
print '(" Calling print_trajectory subroutine")'
! Code to build XYZ trajectory file like
! opening output file 3
pinp2=arg(1:index(arg,'.')-1)//'_trajectory.xyz'
print *,'  New info file... ', pinp2
open(unit=3,file=pinp2,status='replace')
!
m=0
! if(debug)  write(4,*) ' printing ', n_neighbor, ' molecule(s) - n_neighbor'
!
! Needs to set up tmolec-1 due to the following loop
if(n_neighbor == 0) n_neighbor=tmolec-1
if(n_neighbor > tmolec)then
  print '("  Number of requested is larger than total molecules found")'
  print '("  Printing Total Number of Molecules ",i3)', tmolec
  n_neighbor=tmolec-1
end if
!
!if(debug) write(4,*) slt, 'slt from print_trajec'
do l=1, n_neighbor+1
  m=m+imolec(molec_map(l,slt))
  write(3,'(i5)') m
  if(lcharge)then
    write(3,'(a15,1x,i5,1x,a4,i3,1x,a10,f9.5,1x,a10,1x,i3)') 'Adding Molecule', l, 'with', imolec(molec_map(l,slt)),&
    'atom(s) R=', molec_bond(molec_map(l,slt),slt), 'and charge', moleccharge(molec_map(l,slt))
  else
    write(3,'(a15,1x,i5,1x,a4,i3,1x,a10,f9.5)') 'Adding Molecule', l, 'with', imolec(molec_map(l,slt)),&
    'atom(s) R=', molec_bond(molec_map(l,slt),slt)
  end if
  do k=1,l
    do j=1, totatom
      if(molec(j) == molec_map(k,slt).and.(full_cell)) write (3,'(A2,2X,f12.8,2X,f12.8,2X,f12.8)')&
      symb(j),xcoord(j)-cellparam,ycoord(j)-cellparam,zcoord(j)-cellparam
      if(molec(j) == molec_map(k,slt).and.(full_cell.eqv..False.)) write (3,'(A2,2X,f12.8,2X,f12.8,2X,f12.8)')&
      symb(j),xcoord(j),ycoord(j),zcoord(j)
    end do
  end do
end do
!
close(3,status='keep')

if(pstrajec)then
  do n=1, n_neighbor+1,tstep
    write(ci,FMT='(i3.3)') n-1
    pinp2=arg(1:index(arg,'.')-1)//'_'//ci//'_trajectory.xyz'
    print *,'  New info file... ', pinp2
    open(unit=3,file=pinp2,status='replace')
    m=0
    i=0
    do l=1, n_neighbor+1
      m=m+imolec(molec_map(l,slt))
      if(lcharge) i=i+moleccharge(molec_map(l,slt))
      if(l == n)then
        write(3,'(i5)') m
        if(lcharge)then
          write(3,'(a10,1x,i5,1x,a31,1x,i3)') 'Solute and', l-1, 'nearest neighbor(s) with charge', i   
        else
          write(3,'(a10,1x,i5,1x,a19)') 'Solute and', l-1, 'nearest neighbor(s)'
        end if
        do k=1,l
          do j=1, totatom
            if(molec(j) == molec_map(k,slt).and.(full_cell)) write (3,'(A2,2X,f12.8,2X,f12.8,2X,f12.8)')&
            symb(j),xcoord(j)-cellparam,ycoord(j)-cellparam,zcoord(j)-cellparam
            if(molec(j) == molec_map(k,slt).and.(full_cell.eqv..False.)) write (3,'(A2,2X,f12.8,2X,f12.8,2X,f12.8)')&
            symb(j),xcoord(j),ycoord(j),zcoord(j)
          end do
        end do
      end if
    end do
    close(3,status='keep')
  end do
end if
!
end subroutine print_trajectory

! subroutine print_axsf
subroutine print_axsf()
integer :: j, k, l, m, n
print '(" Calling print_axsf subroutine")'
! Code to build Trajectory AXSF
! opening output file 3
pinp2=arg(1:index(arg,'.')-1)//'_trajectory.axsf'
print *,'  New info file... ', pinp2
open(unit=3,file=pinp2,status='replace')
m=0
! if(debug)  write(4,*) ' printing ', n_neighbor, ' molecule(s) - n_neighbor'
! Needs to set up tmolec-1 due to the following loop
if(n_neighbor == 0) n_neighbor=tmolec-1
if(n_neighbor > tmolec)then
  print '("  Number of requested is larger than total molecules found")'
  print '("  Printing Total Number of Molecules ",i3)', tmolec
  n_neighbor=tmolec-1
end if
!

write(3,'(1x,a9,5x,i7)') 'ANIMSTEPS', n_neighbor+1
write(3,'(1x,a7)') 'CRYSTAL'
do l=1, n_neighbor+1
  write(3,'(1x,a7,5x,i7)') 'PRIMVEC', l
  write(3,'(3f15.10)') cellparam, 0.0d0, 0.0d0
  write(3,'(3f15.10)') 0.0d0, cellparam, 0.0d0
  write(3,'(3f15.10)') 0.0d0, 0.0d0, cellparam
  write(3,'(1x,a7,5x,i7)') 'CONVVEC', l
  write(3,'(3f15.10)') cellparam, 0.0d0, 0.0d0
  write(3,'(3f15.10)') 0.0d0, cellparam, 0.0d0
  write(3,'(3f15.10)') 0.0d0, 0.0d0, cellparam
  write(3,'(1x,a9,5x,i7)') 'PRIMCOORD', l
  m=m+imolec(molec_map(l,slt))
  write(3,'(5x,i7,5x,i7)') m, 1   
  do k=1,l
    do j=1, totatom
      if(molec(j) == molec_map(k,slt).and.(full_cell)) write (3,'(i2,2X,f12.8,2X,f12.8,2X,f12.8)')&
      atonum(symb(j)),xcoord(j)-cellparam,ycoord(j)-cellparam,zcoord(j)-cellparam
      if(molec(j) == molec_map(k,slt).and.(full_cell.eqv..False.)) write (3,'(i2,2X,f12.8,2X,f12.8,2X,f12.8)')&
      atonum(symb(j)),xcoord(j),ycoord(j),zcoord(j)
    end do
  end do
end do
!
close(3,status='keep')

if(psaxsf)then
  do n=1, n_neighbor+1,axsfstep
    write(ci,FMT='(i3.3)') n-1
    pinp2=arg(1:index(arg,'.')-1)//'_'//ci//'_trajectory.axsf'
    print *,'  New info file... ', pinp2
    open(unit=3,file=pinp2,status='replace')
    write(3,'(1x,a9,5x,i7)') 'ANIMSTEPS', 1
    write(3,'(1x,a7)') 'CRYSTAL'
    m=0d0
    do l=1, n_neighbor+1
      m=m+imolec(molec_map(l,slt))
      if(l == n)then
        write(3,'(1x,a7,5x,i7)') 'PRIMVEC', 1
        write(3,'(3f15.10)') cellparam, 0.0d0, 0.0d0
        write(3,'(3f15.10)') 0.0d0, cellparam, 0.0d0
        write(3,'(3f15.10)') 0.0d0, 0.0d0, cellparam
        write(3,'(1x,a7,5x,i7)') 'CONVVEC', 1
        write(3,'(3f15.10)') cellparam, 0.0d0, 0.0d0
        write(3,'(3f15.10)') 0.0d0, cellparam, 0.0d0
        write(3,'(3f15.10)') 0.0d0, 0.0d0, cellparam
        write(3,'(1x,a9,5x,i7)') 'PRIMCOORD', 1
        write(3,'(5x,i7,5x,i7)') m, 1
        do k=1,l
          do j=1, totatom
            if(molec(j) == molec_map(k,slt).and.(full_cell)) write (3,'(i2,2X,f12.8,2X,f12.8,2X,f12.8)')&
            atonum(symb(j)),xcoord(j)-cellparam,ycoord(j)-cellparam,zcoord(j)-cellparam
            if(molec(j) == molec_map(k,slt).and.(full_cell.eqv..False.)) write (3,'(i2,2X,f12.8,2X,f12.8,2X,f12.8)')&
            atonum(symb(j)),xcoord(j),ycoord(j),zcoord(j)
          end do
        end do
      end if
    end do
    close(3,status='keep')
  end do
end if
!
end subroutine print_axsf


! subroutine print_within
subroutine print_within()
print '(" Calling print_within subroutine")'
print '("  Total number solute atoms: ",i4)', totsolute
print '("  Cutoff distance in Angstrons: ",f6.2)', rcut
! Code to build XYZ of Solute + atoms in giving cutoff distance
! opening output file 3
pinp2=arg(1:index(arg,'.')-1)//'_within.xyz'
print *,'  New info file... ', pinp2
open(unit=3,file=pinp2,status='replace')
!
m=0
do l=1, totsolute
  do k=totsolute+1,totatom
    if(bond(l,k) <= rcut)then
      m=m+1
    end if
  end do
end do
!
write (3,'(i5)') m+totsolute
write (3,'(a31,1x,f6.2,1x,a11)') 'Solute and solvent atoms within',rcut,'Angstron(s)'
! printing solute atoms
do j=1, totsolute
  write (3,'(A2,2X,f12.8,2X,f12.8,2X,f12.8)') symb(j),xcoord(j),ycoord(j),zcoord(j)
end do
! printing solvent atoms
do l=1, totsolute
  do k=totsolute+1,totatom
    if(bond(l,k) <= rcut)then
      write (3,'(A2,2X,f12.8,2X,f12.8,2X,f12.8)') symb(k),xcoord(k),ycoord(k),zcoord(k)
    end if
  end do
end do
!
close(3,status='keep')
end subroutine print_within

! subroutine print_scell
subroutine print_scell()
print '(" Calling print_scell subroutine")'

! printing geometry in cartesian format
!opening output file
pinp2=arg(1:index(arg,'.')-1)//'_supercell.xyz'
print *,'  New info file... ', pinp2
open(unit=3,file=pinp2,status='replace')
!
write(3,'(i5)') totatom
write(3,'(" Super Cell Cartesian Coordiantes")')
do j=1, totatom
  write (3,'(A2,2X,f12.8,2X,f12.8,2X,f12.8)') symb(j),xcoord(j),ycoord(j),zcoord(j)
end do
close(3,status='keep')
!
! printing geometry in axsf format
!opening output file
pinp2=arg(1:index(arg,'.')-1)//'_supercell.axsf'
print *,'  New info file... ', pinp2
open(unit=3,file=pinp2,status='replace')
!
write(3,'(1x,a9,5x,i7)') 'ANIMSTEPS', 1
write(3,'(1x,a7)') 'CRYSTAL'
write(3,'(1x,a7,5x,i7)') 'PRIMVEC', 1
write(3,'(3f15.10)') cellparam*3, 0.0d0, 0.0d0
write(3,'(3f15.10)') 0.0d0, cellparam*3, 0.0d0
write(3,'(3f15.10)') 0.0d0, 0.0d0, cellparam*3
write(3,'(1x,a7,5x,i7)') 'CONVVEC', 1
write(3,'(3f15.10)') cellparam*3, 0.0d0, 0.0d0
write(3,'(3f15.10)') 0.0d0, cellparam*3, 0.0d0
write(3,'(3f15.10)') 0.0d0, 0.0d0, cellparam*3
write(3,'(1x,a9,5x,i7)') 'PRIMCOORD', 1
write(3,'(5x,i7,5x,i7)') totatom, 1
do j=1, totatom
  write (3,'(i2,2X,f12.8,2X,f12.8,2X,f12.8)') atonum(symb(j)),xcoord(j),ycoord(j),zcoord(j)
end do
close(3,status='keep')
end subroutine print_scell

! subroutine print_nns
subroutine print_neighbor()
print '(" Calling print_neighbor subroutine")'
! opening output file 3
pinp2=arg(1:index(arg,'.')-1)//'_neighbor.xyz'
print *,'  New info file... ', pinp2
open(unit=3,file=pinp2,status='replace')
!
i=0
m=0
! if(debug) write(4,*) ' printing ', n_neighbor, ' molecule(s) - n_neighbor'
!
! Needs to set up tmolec-1 due to the following loop
!if(debug) write(4,*) slt, 'slt from print_nns'
if(n_neighbor == 0) n_neighbor=tmolec-1
if(n_neighbor >= tmolec)then
  print '("  Requested number is larger than total molecules found")'
  print '("  Printing Total Number of Molecules ",i3)', tmolec
  n_neighbor=tmolec-1
end if
! if(debug)  write(4,*) ' printing ', tmolec, ' molecule(s) - tmolec' 
!
! Calculating the number of atoms of n_neighbor+1
do l=1, n_neighbor+1
  m=m+imolec(molec_map(l,slt))
  if(lcharge) i=i+moleccharge(molec_map(l,slt))
end do
write(3,'(i5)') m
if(lcharge)then
  write(3,'(a20,1x,i5,1x,a31,i3)') 'Solute Molecule with', n_neighbor, 'Nearest Neighbor(s) with charge', i
else
  write(3,'(a20,1x,i5,1x,a19)') 'Solute Molecule with', n_neighbor, 'Nearest Neighbor(s)'
end if  
do k=1,n_neighbor+1
  do j=1, totatom
    if(molec(j) == molec_map(k,slt).and.(full_cell)) write (3,'(A2,2X,f12.8,2X,f12.8,2X,f12.8)')&
    symb(j),xcoord(j)-cellparam,ycoord(j)-cellparam,zcoord(j)-cellparam
    if(molec(j) == molec_map(k,slt).and.(full_cell.eqv..False.)) write (3,'(A2,2X,f12.8,2X,f12.8,2X,f12.8)')&
    symb(j),xcoord(j),ycoord(j),zcoord(j)
  end do
end do
close(3,status='keep')
end subroutine print_neighbor
!
end module param
