push_banner_flag 0
# Takes L-J parameters for the QM molecule and combine with L-J parameters for the solvent.
# Uses m_n_vdw (from http://www.cse.scitech.ac.uk/ccg/software/chemshell/manual/dl_poly.html)
# Syntax: m_n_vdw <type1><type2><m><n><r0><epsilon> where r0 = sigma*[2^(1/6)]

# L-J parameters from UFF forcefield. Taken from Chemshell UFF file
set elem {h he li be b c n o f ne na mg al si p s cl ar k ca sc ti v cr mn fe co ni cu zn ga ge as se br kr rb sr y zr nb mo tc ru rh pd ag cd in sn sb te i xe cs ba la ce pr nd pm sm eu gd tb dy ho er tm yb lu hf ta w re os ir pt au hg tl pb bi po at rn fr ra ac th pa u np pu am cm bk cf es fm md no lr}
set vdwr0 {2.886 2.362 2.451 2.745 4.083 3.851 3.660 3.500 3.364 3.243 2.983 3.021 4.499 4.295 4.147 4.035 3.947 3.868 3.812 3.399 3.295 3.175 3.144 3.023 2.961 2.912 2.872 2.834 3.495 2.763 4.383 4.280 4.230 4.205 4.189 4.141 4.114 3.641 3.345 3.124 3.165 3.052 2.998 2.963 2.929 2.899 3.148 2.848 4.463 4.392 4.420 4.470 4.50  4.404 4.517 3.703 3.522 3.556 3.606 3.575 3.547 3.520 3.493 3.368 3.451 3.428 3.409 3.391 3.374 3.355 3.640 4.141 3.170 3.069 2.954 3.120 2.840 2.754 3.293 2.705 4.347 4.297 4.370 4.709 4.750 4.765 4.90  3.677 3.478 3.396 3.424 3.395 3.424 3.424 3.381 3.326 3.339 3.313 3.299 3.286 3.274 3.248 3.236}
set vdweps {0.044 0.056 0.025 0.085 0.180 0.105 0.069 0.060 0.050 0.042 0.030 0.111 0.505 0.402 0.305 0.274 0.227 0.185 0.035 0.238 0.019 0.017 0.016 0.015 0.013 0.013 0.014 0.015 0.005 0.124 0.415 0.379 0.309 0.291 0.251 0.220 0.04 0.235 0.072 0.069 0.059 0.056 0.048 0.056 0.053 0.048 0.036 0.228 0.599 0.567 0.449 0.398 0.339 0.332 0.045 0.364 0.017 0.013 0.010 0.010 0.009 0.008 0.008 0.009 0.007 0.007 0.007 0.007 0.006 0.228 0.041 0.072 0.081 0.067 0.066 0.037 0.073 0.080 0.039 0.385 0.680 0.663 0.518 0.325 0.284 0.248 0.050 0.404 0.033 0.026 0.022 0.022 0.019 0.016 0.014 0.013 0.013 0.013 0.012 0.012 0.011 0.011 0.011}

set polarlist {O N}

#####################
# Solvent definitions
#####################

# Defining point charges
if {"$solvtype"=="spc" || "$solvtype"=="tip3p" || "$solvtype"=="tips3p"} {
set solvparcharges "
charge OT $OTcharge
charge HT $HTcharge
"
} elseif {"$solvtype"=="tip4p"} {
set solvparcharges "
charge OT $OTcharge
charge HT $HTcharge
charge ET $ETcharge
"
}

# Fixing the water model sigma value for Chemshell use
set conv [expr pow(2,1/6.)]
foreach i $solvsig {
    set type [lindex $i 0]
    set isig [lindex $i 1]
    set ir0 [expr $isig*$conv]
    lappend fftypes $type
    lappend solvr0 [list $type $ir0]
}

puts "fftypes is $fftypes"
set fftypescomb [ffcomb $fftypes]
puts "fftypescomb is $fftypescomb"
puts "solvr0 is $solvr0"

# Defining solvent-solvent LJ repulsion
# Looping through each combination of solvent LJ types:
foreach i $fftypescomb {
    set type1 [lindex $i 0]
    set type2 [lindex $i 1]
    set type1r0 [lindex [lindex $solvr0 [lsearch $fftypes $type1]] 1]
    set type2r0 [lindex [lindex $solvr0 [lsearch $fftypes $type2]] 1]
    set type1eps [lindex [lindex $solveps [lsearch $fftypes $type1]] 1]
    set type2eps [lindex [lindex $solveps [lsearch $fftypes $type2]] 1]
    # Unsure how TIP3P/TIP4P calculate repulsion, so just defaulting to CHARMM way:
    # r0 = (r01+r02)/2, eps = sqrt(eps1*eps2)
    set ljr0 [expr ($type1r0+$type2r0)/2]
    set ljeps [expr sqrt($type1eps*$type2eps)]
    append solvparLJ "\nm_n_vdw $type1 $type2 12 6 $ljr0 $ljeps"
}
puts "solvparLJ is $solvparLJ"

##########################
# Making solute-solvent LJ
##########################

proc soluteffgen {x} {
    global solutetype solutetypes$x elem vdwr0 vdweps charmmelem charmmvdwr0 charmmvdweps solvparcharges fftypes solvr0 solveps solvparLJ molsolventdir polarscale polarlist
    source $molsolventdir/proclist.tcl
    set solutetypes [set solutetypes$x]

    set soluteuniq [listuniq $solutetypes]
    puts "soluteuniq is $soluteuniq"

    set solutesolvcombs [cartprod $soluteuniq $fftypes]
    puts "solutesolvcombs is $solutesolvcombs"

    # Declaring all solute atom types to make sure they can be read by dl_poly
    foreach i [listuniq $solutetypes] {
        set solvparcharges [append solvparcharges "declare $i\n"]
    }

    # Cycling through all the solute-solv combinations, finding LJ parameters
    foreach i $solutesolvcombs {
        set soluteatom [lindex $i 0]
        set solventatom [lindex $i 1]
        puts "soluteatom is $soluteatom and solvent atom is $solventatom"
        # Finding solvent parameters, since these are the same regardless of solute
        set solventr0 [lindex [lindex $solvr0 [lsearch $fftypes $solventatom]] 1]
        puts "solventr0 is $solventr0"
        set solventeps [lindex [lindex $solveps [lsearch $fftypes $solventatom]] 1]
        puts "solventeps is $solventeps"
        # Now finds solute parameters and combines them (depends on UFF/CHARMM/OPLS)
        if {"$solutetype"=="uff"} {
            set soluter0 [lindex $vdwr0 [lsearch $elem $soluteatom]]
            puts "soluter0 is $soluter0"
            set soluteeps [lindex $vdweps [lsearch $elem $soluteatom]]
            puts "soluteeps is $soluteeps"
            set r0comb [expr sqrt($soluter0*$solventr0)]
            puts "r0comb is $r0comb"
            set epscomb [expr sqrt($soluteeps*$solventeps)]
            puts "epscomb is $epscomb"
        } elseif {"$solutetype"=="charmm"} {
            # First have to open the cgenff file and find the line that matches soluteatom
            set fp [open $molsolventdir/cgenff-all36-nbtrim.prm]
            set found 0
            foreach line [split [read $fp] "\n"] {
                set data [regexp -all -inline {\S+} $line]
                if {"$soluteatom"=="[lindex $data 0]"} {
                    puts "$soluteatom found on line: $data"
                    set soluter0 [expr [lindex $data 3]*2]
                    puts "soluter0 is $soluter0"
                    set soluteeps [expr [lindex $data 2]*-1]
                    puts "soluteeps is $soluteeps"
                    set found 1
                }
            }
            # Checking that new atom was found
            if {! $found} {
                puts "Solute atom: $soluteatom was not found in list: $molsolventdir/cgenff-all36-nbtrim.prm.  Exiting."
                exit
            }
            # Setting r0comb and epscomb values
            set r0comb [expr ($soluter0+$solventr0)/2]
            puts "r0comb is $r0comb"
            set epscomb [expr sqrt($soluteeps*$solventeps)]
            puts "epscomb is $epscomb"
        } elseif {"$solutetype"=="opls"} {
            # Do something similar to CHARMM: read .prm file and match solutetypes to parameters
            global oplsfile$x
            set fp [open [set oplsfile$x]]
            set write 0
            set found 0
            foreach line [split [read -nonewline $fp] "\n"] {
                # This ensures we don't start taking data until after "NONBONDED" line of file
                if {[lindex $line 0]=="NONBONDED"} {
                    set write 1
                }
                # This ensures we only write from atom information block
                if {$write && [llength $line] == 7} {
                    # Now we match OPLS type to its epsilon and r/0 parameters
                    if {"$soluteatom"=="[lindex $line 0]"} {
                        puts "$soluteatom found on line: $line"
                        set soluter0 [expr [lindex $line 3]*2]
                        puts "soluter0 is $soluter0"
                        set soluteeps [expr [lindex $line 2]*-1]
                        puts "soluteeps is $soluteeps"
                        set found 1
                    }
                }
            }
            # Checking that new atom was found
            if {! $found} {
                puts "Solute atom: $soluteatom was not found in list: [set oplsfile$x].  Exiting."
                exit
            }
            # Setting r0comb and epscomb values according to how OPLS combines
            set r0comb [expr sqrt($soluter0*$solventr0)]
            puts "r0comb is $r0comb"
            set epscomb [expr sqrt($soluteeps*$solventeps)]
            puts "epscomb is $epscomb"
        } else {
            puts "Invalid solutetype: $solutetype.  Exiting..."
            exit
        }
        # Now we scale r0comb if solute atom is polar list
        ## Can be improved using regex and cleaner list searching
        foreach atom $polarlist {
            # Testing if any atom in list matches beginning of soluteatom name
            if {"[string range $soluteatom 0 [expr [string length $atom]-1]]" == "$atom"} {
                puts "Solute atom $soluteatom is in polarlist {$polarlist}, scaling r0comb by $polarscale"
                set r0comb [expr sqrt($soluter0*$solventr0)*$polarscale]
                puts "New r0comb is $r0comb"
                break
            }
        }
        append solutesolvlj "\nm_n_vdw $soluteatom $solventatom 12 6 $r0comb $epscomb"
    }
    puts "solutesolvlj is $solutesolvlj"
    # Now that all the proper LJ interaction lists are done, writing to forcefield file
    set fp [open md-FF$x.ff "w"]
    puts $fp "# Solvent charge and solute type parameters:"
    puts $fp $solvparcharges
    puts $fp "# Solvent-solvent m_n_vdw parameters:"
    puts $fp "# r0 and epsilon are geometric means of type1 and type2 values"
    puts $fp "# m_n_vdw <type1> <type2> <m> <n> <r0> <epsilon>"
    puts $fp $solvparLJ
    puts $fp "\n# Solute-solvent m_n_vdw parameters:"
    puts $fp $solutesolvlj
    close $fp
}

pop_banner_flag
