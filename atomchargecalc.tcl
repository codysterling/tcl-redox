proc xtbchrg {x} {
    push_banner_flag 0

    global molsys$x solvsys$x xtbtheory$x chrgscale

    set molsys [set molsys$x]
    set solvsys [set solvsys$x]
    set xtbtheory [set xtbtheory$x]

    puts "Doing xTB charges for geometry $x"

    # Running xTB calculation to generate xtb1.out
    energy coords=$molsys $xtbtheory list_option=full

    # Getting solute CM5 charges from output file and scaling
    set xtboutfile [open xtb1.out r]
    set natoms [get_number_of_atoms coords=$molsys]
    set grab undefined
    set count undefined
    global charges
    set charges {}
    while {[gets $xtboutfile line]>=0} {
        if {$grab == yes && $count <= $natoms} {
            set linelist [lreplace $line 0 -1]
            lappend charges [lindex $linelist 3]
            set count [expr $count+1]
        }
        if {[lsearch $line Mulliken/CM5]>=0} {
            set grab yes
            set count 1
        }
    }
    set charges [listmult $charges $chrgscale]

    # Have to add "unset" charges to make the list the length of the whole system
    # This will use .ff file values for solvent charges like OT and HT
    for {set i [llength $charges]} {$i<[get_number_of_atoms coords=$solvsys]} {incr i} {
        lappend charges "unset"
    }
    puts "xTB CM5 charges for $x ([llength $charges]) are $charges"

    # Setting charges variable for use in main file
    variable charges$x $charges

    pop_banner_flag
}

proc xtbmulchrg {x} {
    push_banner_flag 0

    global molsys$x solvsys$x xtbtheory$x chrgscale

    set molsys [set molsys$x]
    set solvsys [set solvsys$x]
    set xtbtheory [set xtbtheory$x]

    puts "Doing xTB charges for geometry $x"

    # Running xTB calculation to generate xtb1.out
    energy coords=$molsys $xtbtheory list_option=full

    # Getting solute CM5 charges from output file and scaling
    set xtboutfile [open xtb1.out r]
    set natoms [get_number_of_atoms coords=$molsys]
    set grab undefined
    set count undefined
    global charges
    set charges {}
    while {[gets $xtboutfile line]>=0} {
        if {$grab == yes && $count <= $natoms} {
            set linelist [lreplace $line 0 -1]
            lappend charges [lindex $linelist 2]
            set count [expr $count+1]
        }
        if {[lsearch $line Mulliken/CM5]>=0} {
            set grab yes
            set count 1
        }
    }
    set charges [listmult $charges $chrgscale]

    # Have to add "unset" charges to make the list the length of the whole system
    # This will use .ff file values for solvent charges like OT and HT
    for {set i [llength $charges]} {$i<[get_number_of_atoms coords=$solvsys]} {incr i} {
        lappend charges "unset"
    }
    puts "xTB Mulliken charges for $x ([llength $charges]) are $charges"

    # Setting charges variable for use in main file
    variable charges$x $charges

    pop_banner_flag
}

proc orcachrg {x} {
    push_banner_flag 0

    global molsolventdir molsys$x solvsys$x orcatheory$x solutelist$x chrgscale

    source $molsolventdir/proclist.tcl

    set molsys [set molsys$x]
    set solvsys [set solvsys$x]
    set orcatheory [set orcatheory$x]
    set solutelist [set solutelist$x]

    puts "Doing ORCA CM5 charges for geometry $x"

    # Running energy to generate orca1.out
    energy coords=$molsys $orcatheory

    # Finding Hirshfeld charges and converting to CM5 and scaling
    set charges [grabHirs $solutelist]
    set charges [CM5convert $molsys $charges]
    set charges [listmult $charges $chrgscale]

    # Have to add "unset" charges to make the list the length of the whole system
    # This will use .ff file values for solvent charges like OT and HT
    for {set i [llength $charges]} {$i<[get_number_of_atoms coords=$solvsys]} {incr i} {
        lappend charges "unset"
    }
    puts "ORCA CM5 charges for $x ([llength $charges]) are $charges"

    # Setting charges variable for use in main file
    variable charges$x $charges

    pop_banner_flag
}

proc blankcharge {x} {
    push_banner_flag 0

    global solvsys$x

    set solvsys [set solvsys$x]

    # Have to add "unset" charges for length of whole system
    for {set i 0} {$i<[get_number_of_atoms coords=$solvsys]} {incr i} {
        lappend charges "unset"
    }
    puts "blanked charges for $x ([llength $charges]) are $charges"

    # Setting charges variable for use in main file
    variable charges$x $charges

    pop_banner_flag
}
