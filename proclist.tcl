# Basic definitions
set bohrtoang 0.52917721092
set angtobohr [expr 1/0.52917721092]
set ehtokcal 627.509469
set ehtoev 27.21138505
set evtokcal 23.06054888
set degtorad [expr 3.1415926535897931/180]
set radtodeg [expr 180/3.1415926535897931]
set tcl_precision 14

#########################
# QM/MM-MD-only additions
#########################

# Finding (and writing) phenol H-bond donating/accepting information for steps in MD
## This can be rewritten by slurping contents into list and doing math/rewriting directly
## See how extrapcut.tcl works, should be able to do something very similar for some speedup
proc write_hbonds {coords} {
    # Setting some basic definitions
    global istep ohfile
    set numatoms [get_number_of_atoms coords=$coords]
    set ophe 12
    set hphe 13
    set odist 100
    set hdist 100
    set ohdist [interatomic_distance coords=$coords i=$hphe j=$ophe unit=angstrom]
    # Loop to find closest O to phenol O-H hydrogen
    push_banner_flag 0
    for {set n [expr $hphe+1]} {$n <= $numatoms} {incr n} {
        set type [lindex [get_atom_entry atom_number=$n coords=$coords] 0]
        if {"$type" == "O"} {
            set distance [interatomic_distance coords=$coords i=$hphe j=$n unit=angstrom]
            if {$distance < $odist} {
                set odist $distance
                set onum $n
            }
        }
    }
    pop_banner_flag
    # Loop to find closest H to phenol O
    push_banner_flag 0
    for {set m [expr $hphe+1]} {$m <= $numatoms} {incr m} {
        set type [lindex [get_atom_entry atom_number=$m coords=$coords] 0]
        if {"$type" == "H"} {
            set distance [interatomic_distance coords=$coords i=$ophe j=$m unit=angstrom]
            if {$distance < $hdist} {
                set hdist $distance
                set hnum $m
            }
        }
    }
    pop_banner_flag
    # Writing data to hbonds file with flags for unusual numbers
    set oflag ""
    set hflag ""
    if {$ohdist > $odist} {
        set oflag "O*"
    }
    if {$ohdist > $hdist} {
        set hflag "H*"
    }
    puts $ohfile "$istep $onum $odist $hnum $hdist $oflag $hflag"
}

# SBC restraint on individual atom. Not documented I think but in code:
# /users/home/ragnarbj/chemshell-360/chemsh/src/chemsh/restraint.c
proc set_sbc_indiv_restraint {qmatoms centre sbc_indiv_radius sbc_indiv_restraint_k} {
    # sbc_individual atom_number x y z r0 K
    set centre_x [lindex $centre 0]; set centre_y [lindex $centre 1]; set centre_z [lindex $centre 2]
    set restraintsettings {}
    # Setting Cartesian restraint for each QM atom and adding to restraintsettings
    for {set atom 1} {$atom <= [llength $qmatoms]} {incr atom 1} {
        set newrest "sbc_individual $atom $centre_x $centre_y $centre_z $sbc_indiv_radius $sbc_indiv_restraint_k"
        lappend restraintsettings $newrest
    }
    return $restraintsettings
}

# Proper sequence from begin to end: seq 5 .. 10 creates: 5 6 7 8 9 10
proc seq {start to end} {
    set result []
    for {set i $start} {$i <= $end} {incr i} {
        lappend result $i
    }
    return $result
}

# Define frozen region as distance from centre of cluster (in angstrom)
proc def_frozen_region {qmregion cluster radbohr} {
    global bohrtoang
    set numatoms [get_number_of_atoms coords=$cluster]
    set centre [get_molecule_centre coords=$cluster]
    puts "Setting frozen region as distance $radbohr bohr from centre $centre (bohr)"
    set centx [lindex $centre 0]
    set centy [lindex $centre 1]
    set centz [lindex $centre 2]
    puts "centx is $centx; centy is $centy; centz is $centz"
    set act { }
    set frozen { }
    for { set a 1 } { $a <= [get_number_of_atoms coords=$cluster] } { incr a } {
        # Making sure "a" (atom number) isn't in qmatoms list
        set runflag "on"
        if { [lsearch $qmregion $a] != "-1" } {
            set runflag "off"
        }
        # Guts of the program: finding distance from centre and adding to frozen if needed
        if { $runflag == "on" } {
            push_banner_flag 0
            set atomx [lindex [get_atom_entry coords=$cluster atom_number=$a] 1]
            set atomy [lindex [get_atom_entry coords=$cluster atom_number=$a] 2]
            set atomz [lindex [get_atom_entry coords=$cluster atom_number=$a] 3]
            set distance [expr sqrt(pow(($centx-$atomx),2)+pow(($centy-$atomy),2)+pow(($centz-$atomz),2))]
            if { $distance > $radbohr } {
                set molecule [get_molecule_members coords=$cluster atom_number=$a]
                for { set c 0 } { $c < [llength $molecule] } {incr c } {
                    lappend frozen [lindex $molecule $c]
                }
            }
            pop_banner_flag
        }
    }
    # Cleaning up frozen list and setting act list as everything not frozen
    set frozen [lsort -unique -integer $frozen]
    set all [seq 1 to $numatoms]
    set act [listcomp $all $frozen]

    return [list $frozen $act]
}

# Huge loop: sets lists for frozen/active regions, sets constraints for active, types, groups
# Assumes all solvent molecules are the same as first in list
# Should give proper lists for con/groups irrespective of Chemshell connectivity scale/tolerance
proc set_system_lists {qmatoms solutetypes system x solvmemb qmsbcrad} {

    global bohrtoang frozrad frozsbcrad
    set centre [get_molecule_centre coords=$system]
    set sysatoms [get_number_of_atoms coords=$system]
    set solvstart [expr [llength $qmatoms]+1]

    ## Quick and dirty fix: if doing periodic MD, set frozrad to be huge to it does all constraints
    ## This needs to be fixed later to be cleaner/more general when making fused input file
    global solvgeom
    if {$solvgeom == "box"} {
        set frozrad 10000
        puts "Updating frozrad for periodic"
    }

    # Setting frozen and active regions:
    set frozact [def_frozen_region $qmatoms $system $frozrad]
    set froz [lindex $frozact 0]
    set act [lindex $frozact 1]
    puts "frozen region for solute $x ([llength $froz]) is $froz"
    puts "active region for solute $x ([llength $act]) is $act"
    # Setting solvent-only active region:
    set actnoqm [lreplace $act 0 [expr [llength $qmatoms]-1]]
    puts "active solvents for solute $x ([llength $actnoqm]) is $actnoqm"

    # Writing SBC restraints:
    set froz_sbc "sbc $centre $frozsbcrad 3"

    set qm_sbc [set_sbc_indiv_restraint $qmatoms $centre $qmsbcrad 3]

    set restraints [lappend qm_sbc $froz_sbc]
    puts "restraints for solute $x are $restraints"

    # Making types list with all atoms (froz/act)
    # QM atoms are lowercase, solvent are xT for clarity
    set types $solutetypes
    for {set i $solvstart} {$i <= $sysatoms} {incr i} {
        set atom [lindex [get_atom_entry coords=$system atom_number=$i] 0]
        if {"$atom" == "Bq"} {set atom "E"}
        append atom "T"
        lappend types $atom
    }
    puts "types for solute $x ([llength $types]) are $types"

    # Making constraints/distance list (conX) for active-region solvent molecules
    # Everything is defined based on the first solvent molecule in the merged system
    # Solvent molecule number data (solvmemb) pulled from insertsolute.tcl
    puts "solvmemb is $solvmemb"
    set solvnum [llength $solvmemb]
    puts "solvnum is $solvnum"
    # Ensuring a whole number of solvent molecules in active region
    set actsolvnum [expr [llength $act]-[llength $qmatoms]]
    puts "actsolvnum is $actsolvnum atoms"
    if {$actsolvnum%$solvnum} {
        puts "Error!  Non-integer number of solvent molecules in active region."
        exit
    }
    set actsolvmol [expr $actsolvnum/$solvnum]
    puts "actsolvmol is $actsolvmol molecules"
    # Getting atom combinations and distance values for first solvent in active region:
    foreach i $solvmemb {
        lappend actsolvmemb [lindex $actnoqm [expr $i-1]]
    }
    puts "actsolvmemb is $actsolvmemb"
    set actsolvcomb [combinations $actsolvmemb 2]
    puts "actsolvcomb ([llength $actsolvcomb]) is $actsolvcomb"
    foreach i $actsolvcomb {
        set num1 [lindex $i 0]
        set num2 [lindex $i 1]
        set dist [string trim [interatomic_distance coords=$system i=$num1 j=$num2]]
        lappend firstcon [list $num1 $num2 $dist]
    }
    puts "firstcon is $firstcon"
    # Applying these combinations to all remaining atoms in order
    # Convoluted but think it's the most effective/general way
    set con $firstcon
    for {set i [expr $solvnum+1]} {$i <= [llength $actnoqm]} {incr i $solvnum} {
    # Finding atom numbers and combinations of next solvent molecule in active list
        set range [lrange $actnoqm $i-1 [expr $i+$solvnum-2]]
        set rangecomb [combinations $range 2]
    # Writing these combinations over the initial ones to copy initial distances
        foreach j $firstcon {
            set index [lsearch $firstcon $j]
            set newatoms [lrange [lindex $rangecomb $index] 0 1]
            set j [lreplace $j 0 1 [lindex $newatoms 0] [lindex $newatoms 1]]
            lappend con $j
        }
    }
    puts "con for solute $x ([llength $con]) is $con"

    # Redoing connections block in .c file to ensure there are no errors
    # First have to wipe the existing block
    connect coords=$system conn=$system scale=0 toler=0
    # Adding solute connections: maybe unnecessary for QM/MM?  Good to have for periodic MM anyway
    set qmcomb [combinations $qmatoms 2]
    foreach i $qmcomb {
        set num1 [lindex $i 0]
        set num2 [lindex $i 1]
        add_connection coords=$system i=$num1 j=$num2
    }
    # Looping through all solvent molecules, making combinations and adding
    for {set i $solvstart} {$i <= $sysatoms} {incr i $solvnum} {
        set connrange [seq $i to [expr $i+$solvnum-1]]
        set connrangecomb [combinations $connrange 2]
        foreach j $connrangecomb {
            set num1 [lindex $j 0]
            set num2 [lindex $j 1]
            add_connection coords=$system i=$num1 j=$num2
        }
    }

    # Making groups list for all atoms, size based on first solvent molecule size
    set groups [list $qmatoms]
    for {set i $solvstart} {$i <= $sysatoms} {incr i $solvnum} {
        set grp [seq $i to [expr $i+$solvnum-1]]
        lappend groups $grp
    }
    puts "groups for solute $x ([llength $groups]) is $groups"

pop_banner_flag

set file3Id [open constraints$x "w"]
puts $file3Id "Constraints for state $x ([llength $con]) applied are: $con"
close $file3Id
set file3Id [open groups$x "w"]
puts $file3Id "Groups for state $x ([llength $groups]) are: $groups"
close $file3Id
set file3Id [open types$x "w"]
puts $file3Id "Types for state $x ([llength $types]) are: $types"
close $file3Id
return [list $froz $act $restraints $types $con $groups]
}

####################################
# Ending QM/MM-MD-specific additions
####################################

# Simple lmap. From http://wiki.tcl.tk/13920
proc lmap {_var list body} {
     upvar 1 $_var var
     set res {}
     foreach var $list {lappend res [uplevel 1 $body]}
     set res
}

# Some vector arithmetic procs below. From http://wiki.tcl.tk/14022
# We need basic scalar operators from [expr] factored out:
foreach op {+ - * / % ==} {proc $op {a b} "expr {\$a $op \$b}"}
proc vec {op a b} {
    if {[llength $a] == 1 && [llength $b] == 1} {
        $op $a $b
    } elseif {[llength $a]==1} {
        lmap i $b {vec $op $a $i}
    } elseif {[llength $b]==1} {
        lmap i $a {vec $op $i $b}
    } elseif {[llength $a] == [llength $b]} {
        set res {}
        foreach i $a j $b {lappend res [vec $op $i $j]}
        set res
    } else {error "length mismatch [llength $a] != [llength $b]"}
}

# Proc to calculate quadratic mean (RMS). Used to calculate RMSD.
proc qmean list {
    set sum 0.0
    foreach value $list { set sum [expr {$sum + $value**2}] }
    return [expr { sqrt($sum / [llength $list]) }]
}

proc CM5convert {frag Hcharges} {
push_banner_flag 0
set numatoms [get_number_of_atoms coords=$frag]
set qmatoms [iota 1 $numatoms]
set bohrtoang 0.52917721092
# CM5 parameters from paper: http://pubs.acs.org/doi/pdfplus/10.1021/ct200866d
set alpha 2.474
set C 0.705
# Full DZlist in SI: http://pubs.acs.org/doi/suppl/10.1021/ct200866d/suppl_file/ct200866d_si_001.pdf
set DZlist {0.0056 -0.1543 0.000 0.0333 -0.1030 -0.0446 -0.1072 -0.0802 -0.0629 -0.1088 0.0184 0.000 -0.0726 -0.0790 -0.0756 -0.0565 -0.0444 -0.0767 0.0130 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 -0.0512 -0.0557 -0.0533 -0.0399 -0.0313 -0.0541 0.0092 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 -0.0361 -0.0393 -0.0376 -0.0281 -0.0220 -0.0381 0.0065 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 -0.0255 -0.0277 -0.0265 -0.0198 -0.0155 -0.0269 0.0046 0.000 0.000 0.000 0.000 0.000}
# Table in paper
set DZZ_16 0.0502
set DZZ_17 0.1747
set DZZ_18 0.1671
set DZZ_67 0.0556
set DZZ_68 0.0234
set DZZ_78 -0.0346
set DZZ_61 [expr -1 * $DZZ_16]
set DZZ_71 [expr -1 * $DZZ_17]
set DZZ_81 [expr -1 * $DZZ_18]
set DZZ_76 [expr -1 * $DZZ_67]
set DZZ_86 [expr -1 * $DZZ_68]
set DZZ_87 [expr -1 * $DZZ_78]

# Covalent radii (not vdW), tabulated by Cramer and Truhlar in CRC book: https://books.google.is/books?id=MmDOBQAAQBAJ&pg=SA9-PA49&lpg=SA9-PA49&dq=van+der+Waals+Radii+truhlar&source=bl&ots=64ynIx3euo&sig=qEk-0HTSoYIeaGbBUGTDmioeoJE&hl=en&sa=X&ei=3-VVVYqYBYHRsgGcwYH4CA&redir_esc=y#v=onepage&q=van%20der%20Waals%20Radii%20truhlar&f=false
# Periodic order from hydrogen to uranium(92)
set RZklist {0.32 0.37 1.30 0.99 0.84 0.75 0.71 0.64 0.60 0.62 1.60 1.40 1.24 1.14 1.09 1.04 1.00 1.01 2.00 1.74 1.59 1.48 1.44 1.30 1.29 1.24 1.18 1.17 1.22 1.20 1.23 1.20 1.20 1.18 1.17 1.16 2.15 1.90 1.76 1.64 1.56 1.46 1.38 1.36 1.34 1.30 1.36 1.40 1.42 1.40 1.40 1.37 1.36 1.36 2.38 2.06 1.94 1.84 1.90 1.88 1.86 1.85 1.83 1.82 1.81 1.80 1.79 1.77 1.77 1.78 1.74 1.64 1.58 1.50 1.41 1.36 1.32 1.30 1.30 1.32 1.44 1.45 1.50 1.42 1.48 1.46 2.42 2.11 2.01 1.90 1.84 1.83}

set CM5charges {}
for { set i 0 } { $i < [llength $qmatoms] } { incr i 1 } {
	set sumtb 0
	set bkk_ 0
	set tkk_ 0
	set d [expr $i + 1]
	set Z [get_atom_znum coords=$frag atom_number=[expr $i + 1]]
	set elem [lindex [get_atom_entry coords=$frag atom_number=[expr $i + 1]] 0]
	set qhpa [lindex $Hcharges $i]
	set RZk [lindex $RZklist [expr $Z - 1]]
		for { set f 1 } { $f <= [llength $qmatoms] } { incr f 1 } {
			if {$f == $d} {
			} else {
				set Z_ [get_atom_znum coords=$frag atom_number=$f]
				set RZk_ [lindex $RZklist [expr $Z_ - 1]]
				set rkk_ [ expr [lindex [ interatomic_distance coords=$frag i=$d j=$f ] 0] * $bohrtoang ]
				set bkk_ [expr  exp(-1 * $alpha * ($rkk_ - $RZk - $RZk_ ))]
				if {[lsearch {1 6 7 8} $Z] >= 0 && [lsearch {1 6 7 8} $Z_] >= 0 }  {
					if  {$Z == $Z_} {
						set DZkZk_ 0
						set tkk_ $DZkZk_
					} else {
						set DZkZk_ [eval set z $[set ly "DZZ_$Z$Z_"]]
						set tkk_ $DZkZk_
					}
				} else {
					set tkk_ [expr [lindex $DZlist [expr $Z - 1]]  - [lindex $DZlist [expr $Z_ - 1]]]
				}
				set sumtb [expr $sumtb + [expr $bkk_ * $tkk_ ]]
			}
		}
	set qcm5 [expr $qhpa + $sumtb]
	lappend CM5charges $qcm5
}
pop_banner_flag
return $CM5charges
}

proc grabHirs {qmatoms} {
set fpout [ open orca1.out r ]
set natoms [llength $qmatoms]
set moltype {}
set code [ gets $fpout line ]
    while { $line != "HIRSHFELD ANALYSIS" && $code != -1 } {
		set code [ gets $fpout line ]
    }
if { $code == -1 } then { chemerr "Couldn't find HIRSHFELD output" }
    # skip lines
    set code [ gets $fpout line ]
    set code [ gets $fpout line ]
    set code [ gets $fpout line ]
    set code [ gets $fpout line ]
    set code [ gets $fpout line ]
    set code [ gets $fpout line ]
    set code [ gets $fpout line ]
    for {set i 0} {$i < $natoms} { incr i } {
		lappend moltype [lindex $line 2]
        set code [ gets $fpout line ]
    }
    close $fpout
        return $moltype
}

# Sums everything in list
proc lsum L {
	expr [join $L +]+0
}

# Average everything in list
proc average list {
	expr ([join $list +])/[llength $list].
}

# Mean squared
proc mean2 list {
	set sum 0
	foreach i $list {set sum [expr {$sum+$i*$i}]}
	expr {double($sum)/[llength $list]}
}

# Unbiased standard deviation
proc stdev list {
    set x [average $list]
    set sum 0
    foreach i $list {set sum [expr {$sum+($i-$x)**2}]}
    expr {sqrt($sum/([llength $list]-1))}
}

# Set list of snapshots to grab
proc snaplist { numsnapshots rangebegin rangeend } {
    # Create list of numbers
	set delsnap [expr ($rangeend - $rangebegin) / $numsnapshots]
    if {! "$delsnap" >= 1 && $numsnapshots != 1} {
        puts "Wanting more snapshots than possible in range, exiting."
        exit
    }
	lappend snaps $rangebegin
	for { set i 1 } { $i < $numsnapshots } { incr i 1 } {
		lappend snaps [expr $rangebegin + $delsnap*$i]
	}
    return $snaps
}

# Changes masses of all H in QM region (solute) to deuterium
proc prepmass {sys qmregion} {
    for {set i 1} {$i <= [llength $qmregion]} {incr i} {
        set sym [lindex [get_atom_entry coords=$sys atom_number=$i] 0]
        if {"$sym" == "H" || "$sym" == "h"} {
            lappend masslist [list $i 2.014]
        }
    }
    return $masslist
}

# Returns all elements in a that are not in b
proc listcomp {a b} {
	set diff {}
	foreach i $a {
		if {[lsearch -exact $b $i]==-1} {
			lappend diff $i
		}
	}
return $diff
}

proc iota {base n} {
	set res {}
    for {set i $base} {$i<$n+$base} {incr i} {
		lappend res $i
	}
    set res
}

# Proc to grab NPA charges from ORCA output
proc grabNPA {qmatoms moltype} {
set fpout [ open orca1.out r ]
set natoms [llength $qmatoms]
set $moltype {}
set code [ gets $fpout line ]
while { $line != " Summary of Natural Population Analysis:" && $code != -1 } {
	set code [ gets $fpout line ]
}
if { $code == -1 } then { chemerr "Couldn't find NPA output" }
# Skip lines
	set code [ gets $fpout line ]
	set code [ gets $fpout line ]
	set code [ gets $fpout line ]
	set code [ gets $fpout line ]
	set code [ gets $fpout line ]
	set code [ gets $fpout line ]
for {set i 0} {$i < $natoms} { incr i } {
	lappend $moltype [lindex $line 2]
    set code [ gets $fpout line ]
}
close $fpout
return [subst $$moltype]
}

# Finds dot product of two vectors
proc dotprod {vec1 vec2} {
    set vec1x [lindex $vec1 0]
    set vec1y [lindex $vec1 1]
    set vec1z [lindex $vec1 2]

    set vec2x [lindex $vec2 0]
    set vec2y [lindex $vec2 1]
    set vec2z [lindex $vec2 2]

    return [expr ($vec1x*$vec2x)+($vec1y*$vec2y)+($vec1z*$vec2z)]
}

# Finds cross product of two vectors
proc xprod {vec1 vec2} {
    set vec1x [lindex $vec1 0]
    set vec1y [lindex $vec1 1]
    set vec1z [lindex $vec1 2]

    set vec2x [lindex $vec2 0]
    set vec2y [lindex $vec2 1]
    set vec2z [lindex $vec2 2]

# Finding x, y, z (i, j, k) components
    lappend prodlist [expr ($vec1y*$vec2z)-($vec2y*$vec1z)]
    lappend prodlist [expr ($vec2x*$vec1z)-($vec1x*$vec2z)]
    lappend prodlist [expr ($vec1x*$vec2y)-($vec2x*$vec1y)]

    return $prodlist
}

# Finds length of a vector (or coordinates)
proc lvec {vec} {
    set vecx [lindex $vec 0]
    set vecy [lindex $vec 1]
    set vecz [lindex $vec 2]

    return [expr sqrt((pow($vecx,2))+(pow($vecy,2))+(pow($vecz,2)))]
}

# Finds distance between two 3D points
proc distance {pt1 pt2} {
    set vecx [expr [lindex $pt2 0]-[lindex $pt1 0]]
    set vecy [expr [lindex $pt2 1]-[lindex $pt1 1]]
    set vecz [expr [lindex $pt2 2]-[lindex $pt1 2]]

    return [expr sqrt((pow($vecx,2))+(pow($vecy,2))+(pow($vecz,2)))]
}

# Adds some number to each value in a list
proc ladd {skra num} {
    set count 0
    foreach i $skra {
        set skra [lreplace $skra $count $count [expr $i+$num]]
        incr count
    }
    return $skra
}

# Multiplies each value in a list by a constant
proc listmult {skra num} {
    set count 0
    foreach i $skra {
        set skra [lreplace $skra $count $count [expr $i*$num]]
        incr count
    }
    return $skra
}

# Removes multiples from list while preserving order
proc listuniq {l} {
    foreach i $l {
        dict set tmp $i 1
    }
    set unique [dict keys $tmp]
    return $unique
}

# Finds all combinations of size from a list
proc combinations { list size } {
    if { $size == 0 } {
        return [list [list]]
    }
    set retval {}
    for { set i 0 } { ($i + $size) <= [llength $list] } { incr i } {
        set firstElement [lindex $list $i]
        set remainingElements [lrange $list [expr { $i + 1 }] end]
        foreach subset [combinations $remainingElements [expr { $size - 1 }]] {
            lappend retval [linsert $subset 0 $firstElement]
        }
    }
    return $retval
}

# Finds all combinations of 2 including matching with self (for ffgen.tcl)
proc ffcomb { list } {
    set comblist {}
    for {set i 0} {$i <= [llength $list]} {incr i} {
        set firstElement [lindex $list $i]
        set remainingElements [lrange $list $i end]
        foreach subset [combinations $remainingElements 1] {
            lappend retval [linsert $subset 0 $firstElement]
        }
    }
    return $retval
}

# Finds Cartesian product of two lists (https://en.wikipedia.org/wiki/Cartesian_product)
# Matches each element of 1 with each of 2: {1 2} {a b} => {1 a} {1 b} {2 a} {2 b}
proc cartprod {list1 list2} {
    foreach i $list1 {
        foreach j $list2 {
            lappend combs [list $i $j]
        }
    }
    return $combs
}

# Strips data from files matching string (can be used on updir) and concatenates into list
proc datastrip {files} {
    foreach i [exec ls {*}[glob $files]] {
        set datfile [open $i r]
        set datlist [read $datfile]
        foreach i $datlist {
            lappend finallist $i
        }
        close $datfile
    }
    return $finallist
}

# Turns a dictionary into a list in order of increasing key values
proc dictolist {dict} {
    for {set i 1} {$i <= [dict size $dict]} {incr i} {
        set val [dict get $dict $i]
        lappend finlist $val
    }
    return $finlist
}

# Sorts dictionaries so values (not keys) increase smallest to largest
# Taken from: https://stackoverflow.com/questions/5726938/sort-tcl-dict-by-value
proc dictvalsort {dict args} {
    set lst {}
    dict for {k v} $dict {lappend lst [list $k $v]}
    return [concat {*}[lsort -index 1 {*}$args $lst]]
}

# Proc lifted from: http://www.flooxs.tec.ufl.edu/FLOOXS%20Manual/Post/Fit.html
# Usage: l is a list of x,y pairs: l = {x1 y1 x2 y2 x3 y3 ...}
# Returns: {slope yint r}

proc bestfitline {l} {
    set sumx 0.0
    set sumy 0.0
    set sum2x 0.0
    set sum2y 0.0
    set sumxy 0.0
    set num 0

    for {set i 0} {$i < [llength $l]/2} {incr i} {
        set x [lindex $l [expr 2*$i]]
        set y [lindex $l [expr 2*$i+1]]
        set sumx [expr $x + $sumx]
        set sum2x [expr $x * $x + $sum2x]
        set sumy [expr $y + $sumy]
        set sum2y [expr $y * $y + $sum2y]
        set sumxy [expr $x * $y + $sumxy]
        incr num
    }

    set denx [expr $num * $sum2x - $sumx * $sumx]
    set deny [expr $num * $sum2y - $sumy * $sumy]
    set top  [expr $num * $sumxy - $sumx * $sumy]

    set corr [expr $top / sqrt($denx * $deny)]
    set slp [expr $top / $denx]
    set inter [expr ($sumy - $sumx * $slp) / $num]

    return "$slp $inter $corr"
}
