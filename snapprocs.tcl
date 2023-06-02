##############################
# GAS PHASE ENERGY CALCULATION
##############################

proc gassnap {x cm type} {
    # Killing if "type" is "none"; allows skipping of high-level gas calculations
    if {"$type" == "none"} {
        return
    }

    global masterdir charge$cm mult$cm ${type}theory$cm

    set charge [set charge$cm]
    set mult [set mult$cm]
    set theory [set ${type}theory$cm]

    puts "Running gas calculation on gas-mol$x.c with charge: $charge and mult: $mult"

    energy list_option=none energy=e coords=$masterdir/snaps/gas-mol$x.c $theory
    exec cp orca1.out $masterdir/calc-snaps/gasmol$x-$cm-$type.out
    # Getting final single point energy (FPSE) from output matrix
    set fspe [get_matrix_element matrix=e indices= {0 0}]
    
    return $fspe
}

###########################
# FULL CLUSTER CALCULATIONS
###########################

proc calcsnapsenergy {x snapdict snapfolder qmtype usepol mmodel pmodel qmsolvradius polradius args} {
global numatomsolute$x solutetypes$x solvwitht solvtypes solvtypesp drudeq drudepol wid rOMa rOMap rOLap aLLdp solvtype masterdir

puts "starting calcsnapsenergy"

# Setting a few variables
set solutetypes [set solutetypes$x]
set numatomsolute [set numatomsolute$x]

# This ensures args doesn't break stuff if it's blank or have any spaces
set args [join $args ""]

# Opening files to print energies
set xafile [open $masterdir/${args}${x}Avie$wid.en w]
set xbfile [open $masterdir/${args}${x}Bvie$wid.en w]
set ipfile [open $masterdir/${args}${x}ip$wid.en w]

# Going through snapdict to choose correct snapshot to use
dict for {key val} $snapdict {
set snapclust $val
# Changing because files are stored in snapfolder now
set snapcluster $snapfolder/$snapclust
fragment $snapcluster old

# Setting pol radius to 0 if pol isn't on, makes doing everything generally much easier
if {$usepol=="no"} {set polradius 0}

# Setting number of atoms in snapshot
set numatomcluster [get_number_of_atoms coords=$snapcluster]
puts "At start there are $numatomcluster atoms ([expr $numatomcluster-$numatomsolute] solv)"

# Setting QM, pol, and MM shells
# Making QM region shell
set solutelist [seq 1 to $numatomsolute]
set qmatoms [makeshell $snapcluster 0 $qmsolvradius $solutelist]

# Making polarizable shell
set polatoms [makeshell $snapcluster 0 [expr $polradius+$qmsolvradius] $solutelist $qmatoms]

# Making list of MM atoms
set fullist [seq 1 to $numatomcluster]
set used [concat $qmatoms $polatoms]
set mmatoms [listcomp $fullist $used]

puts "starting qmatoms ([llength $qmatoms]) are $qmatoms"
puts "starting polatoms ([llength $polatoms]) are $polatoms"
puts "starting mmatoms ([llength $mmatoms]) are $mmatoms"
puts "usepol is $usepol, pmodel is $pmodel, solvtype is $solvtype"

# Only do this stuff if polarization is on
if {$usepol=="yes"} {
    # Adding/moving polarizable M charges if needed
    if {$pmodel=="swm4dp" || $pmodel=="swm4ndp" || $pmodel=="swm6"} {
        if {$solvtype=="tip4p"} {
            snapmovem $snapcluster $rOMap $polatoms
        } else {
            snapaddm $snapcluster $rOMap polatoms qmatoms mmatoms
        }
    } elseif {$pmodel=="pspc" && $solvtype=="tip4p"} {
        snapdelm $snapcluster $polatoms
    }

    # Adding L charges if needed
    if {$pmodel=="swm6"} {
        snapaddl $snapcluster $rOLap $aLLdp polatoms qmatoms mmatoms
    }
}

# Adding/removing MM M charges if needed
if {$solvtype!="tip4p" && $mmodel=="tip4p"} {
    snapaddm $snapcluster $rOMa mmatoms qmatoms polatoms
} elseif {$solvtype=="tip4p" && $mmodel!="tip4p"} {
    snapdelm $snapcluster $mmatoms
}

puts "qmatoms ([llength $qmatoms]) are $qmatoms"
puts "polatoms ([llength $polatoms]) are $polatoms"
puts "mmatoms ([llength $mmatoms]) are $mmatoms"

# Setting Drude list from polarized oxygens
if {$usepol=="yes"} {
if {$pmodel=="swm4dp" || $pmodel=="swm4ndp" || $pmodel=="swm6"} {
    set cospara {}
    push_banner_flag 0
    foreach n $polatoms {
        set elem [lindex [get_atom_entry coords=$snapcluster atom_number=$n] 0]
        if {$elem=="O"} {
            lappend cospara [list $n $drudepol $drudeq]
        }
    }
    pop_banner_flag
    puts "cospara is $cospara"
}
}

# Setting charge list for EFPs
## Really ugly and hardset for TIP3P, needs fixing eventually (reading forcefield.ff file?)
if {$pmodel=="efp"} {
    set solvchrg {-.834 .417 .417}
    set pccharges {}
    set snum 0
    foreach v $mmatoms {
        lappend pccharges [lindex $solvchrg $snum]
        if {$snum<2} {
            incr snum
        } else {
            set snum 0
        }
    }
    puts "pccharges is $pccharges"
}

# Recalculating number of atoms in snapshot after additions/deletions
set numatomcluster [get_number_of_atoms coords=$snapcluster]
puts "After changes there are $numatomcluster atoms ([expr $numatomcluster-$numatomsolute] solv)"

# Making types list, has to figure out which shell each atom/molecule is in and add correct types
# Starting with solute types
set types $solutetypes

for {set i [expr $numatomsolute+1]} {$i<=$numatomcluster} {} {
    if {[lsearch $qmatoms $i]!="-1"} {
        set types [concat $types $solvwitht]
    }
    if {[lsearch $polatoms $i]!="-1"} {
        if {$pmodel=="efp"} {
            set types [concat $types $solvwitht]
        } else {
            set types [concat $types $solvtypesp]
        }
    }
    if {[lsearch $mmatoms $i]!="-1"} {
        set types [concat $types $solvtypes]
    }
    push_banner_flag 0
    set numconn [llength [get_molecule_members coords=$snapcluster atom_number=$i]]
    pop_banner_flag
    incr i $numconn
}
puts "types ([llength $types]) are $types"

if {[llength $types]!=$numatomcluster} {
    puts "types length [llength $types] != total atoms $numatomcluster.  Exiting..."
    exit
}

# If pol list is empty do non-pol energy calculations
if {[llength $polatoms]==0} {
    puts "Running non-polarized calculation for snap $snapclust"
    clustersnap $x $key $snapfolder $snapclust $qmtype $qmatoms $mmatoms $types $xafile $xbfile $ipfile $args
# Otherwise do polarized: Drude or EFP
} elseif {$pmodel=="swm4dp" || $pmodel=="swm4ndp" || $pmodel=="swm6"} {
    puts "Running Drude polarized calculation for snap $snapclust"
    clustersnappol $x $key $snapfolder $snapclust $qmtype $qmatoms $polatoms $cospara $mmatoms $types $xafile $xbfile $ipfile $args
} elseif {$pmodel=="efp"} {
    puts "Running EFP calculation for snap $snapclust"
    clustersnapefp $x $key $snapfolder $snapclust $qmatoms $polatoms $mmatoms $pccharges $types $xafile $xbfile $ipfile $args
}
delete_object $snapcluster
}

# Closing files with lists and moving them up a directory (from workgroup to main)
close $xafile
close $xbfile
close $ipfile
#set objlist [list "${args}${x}Avie$wid.en" "${args}${x}Bvie$wid.en" "${args}${x}ip$wid.en"]
#foreach k $objlist {exec cp $k ..}
#exec cp -r ${args}calc-snaps ..
}

#################################################
# FULL CLUSTER DELTAE CALCULATION (NON-POLARIZED)
#################################################

proc clustersnap {x key snapfolder snapclust qmtype qmatoms mmatoms types xafile xbfile ipfile args} {
global cutoff no_elec coupling chargeA multA chargeB multB ehtokcal masterdir
global ${qmtype}theoryA ${qmtype}theoryB

# Setting qmtheory variable
set qmtheoryA [set ${qmtype}theoryA]
set qmtheoryB [set ${qmtype}theoryB]

set args [join $args ""]

# Setting snapclust to folder
set snapcluster $snapfolder/$snapclust

# Setting mxlists based on number of atoms
set mxlist [expr [get_number_of_atoms coords=$snapcluster]/2]
set mxexcl [llength $qmatoms]

# Setting basic snaptheory stuff
set mmargs [ list \
    list_option=none \
    cutoff=$cutoff \
    no_elec=$no_elec \
    use_pairlist=yes \
    atom_types= $types \
    mxlist=$mxlist \
    mxexcl=$mxexcl \
    mm_defs=forcefield.ff ]

puts "Running calculation on $snapcluster with charge $chargeA for state A trajectory"

set ahybargs [ list \
    theory=hybrid : [ list \
        coupling=$coupling \
        list_option=none \
        qm_region= $qmatoms \
        qm_$qmtheoryA \
        mm_theory=dl_poly : $mmargs ]]

energy list_option=none energy=e coords=$snapcluster $ahybargs
if {[string match *orca* $qmtype]} {
    set outfile orca1.out
} elseif {[string match *xtb* $qmtype]} {
    set outfile xtb1.out
} elseif {[string match *m3* $qmtype]} {
    set outfile mndo.out
}
exec cp $outfile $masterdir/${args}calc-snaps/$snapclust-$qmtype-$x-chargeA.out

if { $coupling == "shift" && [string match *orca* $qmtype] } {
    exec cp pointcharges.xyz $masterdir/${args}calc-snaps/$snapclust-PC-$x-chargeA.xyz
}
set en1 [ get_matrix_element matrix=e indices= {0 0} ]
set en1 [expr $en1*$ehtokcal]

puts "Running calculation on $snapcluster with charge $chargeB for state B trajectory"
set bhybargs [ list \
    theory=hybrid : [ list \
        coupling=$coupling \
        list_option=none \
        qm_region= $qmatoms \
        qm_$qmtheoryB \
        mm_theory=dl_poly : $mmargs ]]

energy list_option=none energy=e coords=$snapcluster $bhybargs
if {[string match *orca* $qmtype]} {
    set outfile orca1.out
} elseif {[string match *xtb* $qmtype]} {
    set outfile xtb1.out
} elseif {[string match *m3* $qmtype]} {
    set outfile mndo.out
}
exec cp $outfile $masterdir/${args}calc-snaps/$snapclust-$qmtype-$x-chargeB.out

if { $coupling == "shift" && $qmtype == "orca" } {
    exec cp pointcharges.xyz $masterdir/${args}calc-snaps/$snapclust-PC-$x-chargeB.xyz
}
set en2 [ get_matrix_element matrix=e indices= {0 0} ]
set en2 [expr $en2*$ehtokcal]

set IP [expr ($en2 - $en1)]
puts "deltaE for snap $snapclust is $IP kcal/mol"

# Writing data to appropriate files
puts $xafile "$key $en1"
puts $xbfile "$key $en2"
puts $ipfile "$key $IP"

delete_object $snapcluster
}

#############################################
# FULL CLUSTER DELTAE CALCULATION (POLARIZED)
#############################################

proc clustersnappol {x key snapfolder snapclust qmtype qmatoms polatoms cospara mmatoms types xafile xbfile ipfile args} {
global cutoff no_elec coupling chargeA multA chargeB multB ehtokcal masterdir
global ${qmtype}theoryA ${qmtype}theoryB

# Setting snapclust to folder
set snapcluster $snapfolder/$snapclust

# Setting qmtheory variable
set qmtheoryA [set ${qmtype}theoryA]
set qmtheoryB [set ${qmtype}theoryB]

set args [join $args ""]

# Setting mxlists based on number of atoms
set mxlist [expr [get_number_of_atoms coords=$snapcluster]/2]
set mxexcl [llength $qmatoms]

set topfile newtop.top
set prmfile spc_water_param.inp

set mmargs [ list \
    list_option=none \
    cutoff=$cutoff \
    no_elec=$no_elec \
    atom_types= $types \
    scale14= {1.0 0.0} \
    mxexcl=$mxexcl \
    mxlist=$mxlist \
    exact_srf= yes \
    use_charmm_psf=no \
    mm_defs=forcefield.ff \
    charmm_parameter_file=$prmfile \
    charmm_mass_file=$topfile ]

set polarg [ list \
    polcos_scale14=1.0 \
    polcos_maxcycle=50 \
    polcos_toler_energy=5e-4 \
    polcos_maxdx=2.0e-4 \
    polcos_rmsdx=1.0e-4 \
    polcos_efield=drude \
    polcos_inmaxcyc=1000 \
    polcos_atom_polcosq= $cospara ]

puts "Running calculation on $snapcluster with charge $chargeA for state A trajectory"

set ahybargs [ list \
    theory=hybrid : [ list \
        coupling=$coupling \
        qm_region= $qmatoms \
        list_option=none \
        mm_polcos=yes : [ list \
            $polarg ] \
        qm_$qmtheoryA \
        mm_theory=dl_poly : $mmargs ]]

energy list_option=full energy=e coords= $snapcluster $ahybargs
if {[string match *orca* $qmtype]} {
    set outfile orca1.out
} elseif {[string match *xtb* $qmtype]} {
    set outfile xtb1.out
} elseif {[string match *m3* $qmtype]} {
    set outfile mndo.out
}
exec cp $outfile $masterdir/${args}calc-snaps/$snapclust-$qmtype-$x-chargeA.out

if { $coupling == "shift" && [string match *orca* $qmtype] } {
    exec cp pointcharges.xyz $masterdir/${args}calc-snaps/$snapclust-PC-$x-chargeA.xyz
}
set en1 [ get_matrix_element matrix=e indices= {0 0} ]
set en1 [expr $en1*$ehtokcal]

puts "Running calculation on $snapcluster with charge $chargeB for state B trajectory"

set bhybargs [ list \
    theory=hybrid : [ list \
        coupling=$coupling \
        qm_region= $qmatoms \
        list_option=none \
        mm_polcos=yes : [ list \
            $polarg ] \
        qm_$qmtheoryB \
        mm_theory=dl_poly : $mmargs ]]

energy list_option=full energy=e coords= $snapcluster $bhybargs
if {[string match *orca* $qmtype]} {
    set outfile orca1.out
} elseif {[string match *xtb* $qmtype]} {
    set outfile xtb1.out
} elseif {[string match *m3* $qmtype]} {
    set outfile mndo.out
}
exec cp $outfile $masterdir/${args}calc-snaps/$snapclust-$qmtype-$x-chargeB.out

if { $coupling == "shift" && [string match *orca* $qmtype] } {
    exec cp pointcharges.xyz $masterdir/${args}calc-snaps/$snapclust-PC-$x-chargeB.xyz
}
set en2 [ get_matrix_element matrix=e indices= {0 0} ]
set en2 [expr $en2*$ehtokcal]

set IP [expr ($en2 - $en1)]
puts "deltaE for snap $snapclust is $IP kcal/mol"

# Writing data to appropriate files
puts $xafile "$key $en1"
puts $xbfile "$key $en2"
puts $ipfile "$key $IP"

delete_object $snapcluster
}

#######################################
# FULL CLUSTER DELTAE CALCULATION (EFP)
#######################################

proc clustersnapefp {x key snapfolder snapclust qmatoms polatoms mmatoms pccharges types xafile xbfile ipfile args} {
global cutoff no_elec coupling chargeA multA chargeB multB ehtokcal psi4path molsolventdir psi4cores masterdir

# Manually sourcing Psi4 Chemshell interface
source $molsolventdir/chemshint-psi4.tcl

set args [join $args ""]

# Setting snapclust to folder
set snapcluster $snapfolder/$snapclust

# Setting mxlists based on number of atoms
set mxlist [expr [get_number_of_atoms coords=$snapcluster]/2]
set mxexcl [llength $qmatoms]

# Setting dl-poly MM arguments
set mmargs [ list \
    list_option=none \
    cutoff=$cutoff \
    no_elec=$no_elec \
    atom_types= $types \
    scale14= {1.0 0.0} \
    mxexcl=$mxexcl \
    mxlist=$mxlist \
    exact_srf= yes \
    use_charmm_psf=no \
    mm_defs=forcefield.ff ]

# Setting variables for EFP calculations
set psi4method wb97x
set psi4input "
set basis def2-SVP
set dft_spherical_points 434
set dft_radial_points 100
set diis_max_vecs 20
set maxiter 300
"

puts "Running calculation on $snapcluster with charge $chargeA for state A trajectory"

set hybridargs [ list \
    coupling=$coupling \
    qm_region= $qmatoms \
    list_option=none \
    qm_theory=psi4 : [ list \
        executable=$psi4path/psi4 \
        method=$psi4method \
        psi4input= $psi4input \
        psi4cores= $psi4cores \
        guess=auto \
        charge=$chargeA \
        mult=$multA \
        efp=yes \
        efpatoms= $polatoms \
        pcatoms= $mmatoms \
        pccharges= $pccharges \
        fullcoords= $snapcluster ] \
    mm_theory=dl_poly : $mmargs ]

energy list_option=full energy=e coords= $snapcluster theory= hybrid : [ list $hybridargs ]
exec cp psi41.out ./calc-snaps/$snapclust-psi4-$x-chargeA.out

set en1 [ get_matrix_element matrix=e indices= {0 0} ]
set en1 [expr $en1*$ehtokcal]

puts "Running calculation on $snapcluster with charge $chargeB for state B trajectory"

set hybridargs [ list \
    coupling=$coupling \
    qm_region= $qmatoms \
    list_option=none \
    qm_theory=psi4 : [ list \
        executable=$psi4path/psi4 \
        method=$psi4method \
        psi4input= $psi4input \
        guess=sad \
        charge=$chargeB \
        mult=$multB \
        efp=yes \
        efpatoms= $polatoms \
        pcatoms= $mmatoms \
        pccharges= $pccharges \
        fullcoords= $snapcluster ] \
    mm_theory=dl_poly : $mmargs ]

energy list_option=full energy=e coords= $snapcluster theory= hybrid : [ list $hybridargs ]
exec cp psi41.out $masterdir/${args}calc-snaps/$snapclust-psi4-$x-chargeB.out

set en2 [ get_matrix_element matrix=e indices= {0 0} ]
set en2 [expr $en2*$ehtokcal]

set IP [expr ($en2 - $en1)]
puts "deltaE for snap $snapclust is $IP kcal/mol"

# Writing data to appropriate files
puts $xafile "$key $en1"
puts $xbfile "$key $en2"
puts $ipfile "$key $IP"

delete_object $snapcluster
}

##############################################
# Other procs used in snap energy calculations
##############################################

proc centerox {currO currH1 currH2} {
# Setting the initial coordinates
    set Oxinit [lindex $currO 1]
    set Oyinit [lindex $currO 2]
    set Ozinit [lindex $currO 3]
    set H1xinit [lindex $currH1 1]
    set H1yinit [lindex $currH1 2]
    set H1zinit [lindex $currH1 3]
    set H2xinit [lindex $currH2 1]
    set H2yinit [lindex $currH2 2]
    set H2zinit [lindex $currH2 3]
# Setting all the coords with O at origin
    set Ox [expr $Oxinit-$Oxinit]
    set Oy [expr $Oyinit-$Oyinit]
    set Oz [expr $Ozinit-$Ozinit]
    set H1x [expr $H1xinit-$Oxinit]
    set H1y [expr $H1yinit-$Oyinit]
    set H1z [expr $H1zinit-$Ozinit]
    set H2x [expr $H2xinit-$Oxinit]
    set H2y [expr $H2yinit-$Oyinit]
    set H2z [expr $H2zinit-$Ozinit]
# Returning coords
    set Ofix [list $Ox $Oy $Ox]
    set H1fix [list $H1x $H1y $H1z]
    set H2fix [list $H2x $H2y $H2z]
    return [list $Ofix $H1fix $H2fix]
}

proc mcoords {currO currH1 currH2 rOM} {
# Setting the initial coordinates
    set Oxinit [lindex $currO 1]
    set Oyinit [lindex $currO 2]
    set Ozinit [lindex $currO 3]
    set H1xinit [lindex $currH1 1]
    set H1yinit [lindex $currH1 2]
    set H1zinit [lindex $currH1 3]
    set H2xinit [lindex $currH2 1]
    set H2yinit [lindex $currH2 2]
    set H2zinit [lindex $currH2 3]
# Setting all the coords with O at origin
    set Ox [expr $Oxinit-$Oxinit]
    set Oy [expr $Oyinit-$Oyinit]
    set Oz [expr $Ozinit-$Ozinit]
    set H1x [expr $H1xinit-$Oxinit]
    set H1y [expr $H1yinit-$Oyinit]
    set H1z [expr $H1zinit-$Ozinit]
    set H2x [expr $H2xinit-$Oxinit]
    set H2y [expr $H2yinit-$Oyinit]
    set H2z [expr $H2zinit-$Ozinit]
# Finding the midpoint Hm of H1 and H2
    set Hmx [expr ($H1x+$H2x)/2]
    set Hmy [expr ($H1y+$H2y)/2]
    set Hmz [expr ($H1z+$H2z)/2]
# Finding length from O to Hm
    set rOHm [expr sqrt($Hmx**2+$Hmy**2+$Hmz**2)]
    if {$rOM==0} {set rOM $rOHm}
# Finding the recalibrated M coords, then setting it back to the original coords
    set currM {Bq}
    lappend currM [expr ($Hmx*($rOM/$rOHm))+$Oxinit]
    lappend currM [expr ($Hmy*($rOM/$rOHm))+$Oyinit]
    lappend currM [expr ($Hmz*($rOM/$rOHm))+$Ozinit]
# Returning M coordinates
    return $currM
}

proc snapaddm {snap rOMa srange args} {
global angtobohr
upvar 1 $srange range

# Copying snapshot for backup file of original geometry
exec cp $snap ${snap}.mbk
set frag $snap
fragment $frag old persistent

# Setting rOMa to bohr
set rOM [expr $rOMa*$angtobohr]

puts "Now adding M sites for snapshot $snap"

for {set n 0} {$n<[llength $range]} {incr n} {
    push_banner_flag 0
    set i [lindex $range $n]
    set elem [lindex [get_atom_entry coords=$frag atom_number=$i] 0]
    if {$elem == "O"} {
        set currO [get_atom_entry coords=$frag atom_number=$i]
        set currH1 [get_atom_entry coords=$frag atom_number=[expr $i+1]]
        set currH2 [get_atom_entry coords=$frag atom_number=[expr $i+2]]
        set currM [mcoords $currO $currH1 $currH2 $rOM]
# Adding the M particle to the water molecule and updating range list
        if {[expr $i+3]<[get_number_of_atoms coords=$frag]} {
            add_atom_entry coords=$frag atom_entry= $currM atom_number=[expr $i+3]
        } elseif {[expr $i+3]>=[get_number_of_atoms coords=$frag]} {
            add_atom_entry coords=$frag atom_entry= $currM
        }
# Updating range of shell and increasing values of other shell atoms
        set range2 {}
        foreach l $range {
            if {$l>[expr $i+2]} {
                lappend range2 [expr $l+1]
            } else {
                lappend range2 $l
            }
        }
        set range $range2
        set range [lsort -integer [lappend range [expr $i+3]]]
        foreach j $args {
            set jlist2 {}
            upvar $j jlist$j
            foreach k [set jlist$j] {
                if {$k>$i} {
                    lappend jlist2 [expr $k+1]
                } else {
                    lappend jlist2 $k
                }
            }
            set jlist$j $jlist2
        }
    }
    pop_banner_flag
}
}

proc snapmovem {snap rOMa range} {
global angtobohr

# Copying snapshot for backup file of original geometry
exec cp $snap ${snap}.mmbk
set frag $snap
fragment $frag old persistent
source ~/protocols-qm.mm/proclist.tcl

# Setting rOMa to bohr
set rOM [expr $rOMa*$angtobohr]

# Cycling through oxygens, finding attached M site, moving to correct position
foreach i $range {
    set elem [lindex [get_atom_entry coords=$frag atom_number=$i] 0]
    if {$elem == "O"} {
        set currO [get_atom_entry coords=$frag atom_number=$i]
        set currH1 [get_atom_entry coords=$frag atom_number=[expr $i+1]]
        set currH2 [get_atom_entry coords=$frag atom_number=[expr $i+2]]
        set currM [mcoords $currO $currH1 $currH2 $rOM]
# Moving M particle position
        if {[get_atom_entry coords=$frag atom_number=[expr $i+3]] != "Bq"} {
            puts "Oxygen no. $i does not have attached Bq atom!"
        } else {
            replace_bq_entry coords=$frag bq_entry= $currM bq_number=[expr $i+3]
        }
    }
}
}

proc snapdelm {snap range args} {
foreach i $range {
    if {[lindex [get_atom_entry coords=$snap atom_number=$i] 0]=="Bq"} {
        delete_atom_entry coords=$snap atom_number=$i
    }
}
}

proc snapaddl {snap rOLa aLLd srange args} {
global angtobohr degtorad
upvar 1 $srange range

# Copying snapshot for backup file of original geometry
exec cp $snap ${snap}.lbk
set frag $snap
fragment $snap old persistent

# Setting rOLa to bohr and aLLd to rad
set rOL [expr $rOLa*$angtobohr]
set aLL [expr $aLLd*$degtorad/2]
set rLip [expr $rOL*cos($aLL)]
set rLoop [expr $rOL*sin($aLL)]

puts "Now adding L sites for snapshot $snap"

for {set n 0} {$n<[llength $range]} {incr n} {
    push_banner_flag 0
    set i [lindex $range $n]
    set elem [lindex [get_atom_entry coords=$frag atom_number=$i] 0]
    if {$elem == "O"} {
# Calculating coords of L particles
# Redefining coords with O at origin
        set currO [get_atom_entry coords=$frag atom_number=$i]
        set currH1 [get_atom_entry coords=$frag atom_number=[expr $i+1]]
        set currH2 [get_atom_entry coords=$frag atom_number=[expr $i+2]]
        set newcoords [centerox $currO $currH1 $currH2]
        set Ocent [lindex $newcoords 0]
        set H1cent [lindex $newcoords 1]
        set H2cent [lindex $newcoords 2]
        set Ocentm [linsert [lindex $newcoords 0] 0 O]
        set H1centm [linsert [lindex $newcoords 1] 0 H]
        set H2centm [linsert [lindex $newcoords 2] 0 H]
# Finding L particles in plane coords
        set Hmid [lreplace [mcoords $Ocentm $H1centm $H2centm 0] 0 0]
        set Hmidl [lvec $Hmid]
        set Lip [listmult $Hmid [expr $rLip/$Hmidl]]
        set Lip [listmult $Lip -1]
# Finding L particles out of plane coords
        set Hcross [xprod $H1cent $H2cent]
        set Hcrossl [lvec $Hcross]
        set Loop1 [listmult $Hcross [expr $rLoop/$Hcrossl]]
        set Loop2 [listmult $Loop1 -1]
# Combining coords into final L coords
        set currL1 Bq
        set currL2 Bq
        for {set p 0} {$p<3} {incr p} {
            lappend currL1 [expr [lindex $Loop1 $p]+[lindex $Lip $p]+[lindex $currO $p+1]]
        }
        for {set q 0} {$q<3} {incr q} {
            lappend currL2 [expr [lindex $Loop2 $q]+[lindex $Lip $q]+[lindex $currO $q+1]]
        }
# Adding the L particles to the water molecule and updating range list
        if {[expr $i+4]<[get_number_of_atoms coords=$frag]} {
            add_atom_entry coords=$frag atom_entry= $currL1 atom_number=[expr $i+4]
        } elseif {[expr $i+4]>=[get_number_of_atoms coords=$frag]} {
            add_atom_entry coords=$frag atom_entry= $currL1
        }
        if {[expr $i+5]<[get_number_of_atoms coords=$frag]} {
            add_atom_entry coords=$frag atom_entry= $currL2 atom_number=[expr $i+5]
        } elseif {[expr $i+5]>=[get_number_of_atoms coords=$frag]} {
            add_atom_entry coords=$frag atom_entry= $currL2
        }
# Updating range of shell and increasing values of other shell atoms
        set range2 {}
        foreach l $range {
            if {$l>[expr $i+2]} {
                lappend range2 [expr $l+2]
            } else {
                lappend range2 $l
            }
        }
        set range $range2
        set range [lsort -integer [lappend range [expr $i+3] [expr $i+4]]]
        foreach j $args {
            set jlist2 {}
            upvar $j jlist$j
            foreach k [set jlist$j] {
                if {$k>$i} {
                    lappend jlist2 [expr $k+2]
                } else {
                    lappend jlist2 $k
                }
            }
            set jlist$j $jlist2
        }
    }
    pop_banner_flag
}
}

# Makes shells defined by distance from solute (in anstroms)
proc makeshell {snap start end solutelist {remove {}}} {
fragment $snap old persistent

# Finding number of atoms in the snapshot
set numatoms [get_number_of_atoms coords=$snap]
# Finding start of solvent atoms
set solvstart [expr [llength $solutelist]+1]

# If starting at 0 angstrom, include solute atoms in solvshell list
if {$start=="0"} {
    set solvshell $solutelist
} else {
    set solvshell {}
}

# Cycling through each solute atom
foreach i $solutelist {
    for {set j $solvstart} {$j<=$numatoms} {incr j} {
        push_banner_flag 0
        set dist [interatomic_distance coords=$snap i=$i j=$j unit=angstrom]
        if {$dist>$start && $dist<=$end} {
            set solvx [get_molecule_members coords=$snap atom_number=$j]
            set solvshell [concat $solvshell $solvx]
        }
        pop_banner_flag
    }
}

# Sorting list
set solvshell [lsort -unique -integer $solvshell]

# Removing already assigned atoms if needed
foreach k $remove {
    set index [lsearch $solvshell $k]
    if {$index!="-1"} {
        set solvshell [lreplace $solvshell $index $index]
    }
}

concat $remove $solvshell
return $solvshell
}
