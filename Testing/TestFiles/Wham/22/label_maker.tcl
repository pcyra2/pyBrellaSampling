mol new ../complex.parm7 waitfor -1
mol addfile prod_1.22.dcd waitfor -1
label add Bonds 0/13730 0/13717
label graph Bonds 0 C4-C8.dat
label delete Bonds 0

label add Dihedrals 0/13741 0/13738 0/13735 0/13730
label graph Dihedrals 0 SobVert.dat
label delete Dihedrals 0

quit
