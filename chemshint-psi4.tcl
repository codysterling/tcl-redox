#
# Psi4 1.0 interface. Based on modified ORCA interface
#
# CHanging all psi4 to psi4

proc new_psi4 { name args } {
    global psi4
    # Option defaults
    set executable  psi4
    set jobname     $name
    set restart     no
    set moinp       save.gbw
    set charge      0
    set mult        1
    set functional  undefined
    set list_option medium
    set unique_listing no
    set dispersion_correction undefined
    set gridsize    undefined
    set finalgridsize undefined
    set use_finalgrid yes
    set convergence NormalSCF
    set guess undefined
    set maxcyc undefined
    set optstr undefined
    set excited no
    set estate 1
    set eroots 1
    set tda yes
    set image undefined
    set efp undefined
    set efpatoms undefined
    set reference undefined
    set fullcoords undefined
    set nproc 1
    set efpregion undefined
    set pcatoms undefined
    set pccharges undefined
    set energy undefined
    set gradient undefined
    set hessian undefined
    
    # Psi4 input options
    set psi4(in_args) { \
                coords energy gradient hessian \
                executable jobname restart moinp \
                charge mult functional \
                gridsize finalgridsize use_finalgrid
        list_option unique_listing \
                dispersion_correction \
                gridsize finalgridsize use_finalgrid \
                convergence guess nproc \
                maxcyc optstr psi4input method\
                excited estate eroots tda image efp reference efpatoms fullcoords pcatoms pccharges \
            }
    
    # All variables including internal ones
    set psi4(variables) [ concat $psi4(in_args) \
            verbose debug quiet restart_psi4  listing_index ]
    
    # Assign name if not given
    if { "$name" == "*" } {
        if { [ info exists psi4(index)  ] } {
            set index $psi4(index)
            incr psi4(index)
        } else {
            set index 1
            set psi4(index) 2
        }
        set name psi4$index
    }
    
    if { [ parsearg psi4.init $psi4(in_args) $args ] != 0 } {
        chemerror "psi4 failed due to bad argument"
    }
    
    set restart_psi4 [ check_boolean restart ]
    set unique_listing [ check_boolean unique_listing ]
    set excited [ check_boolean excited ]
    set tda [ check_boolean tda ]
    
    
    set verbose 0
    set debug 0
    set quiet 0
    switch $list_option {
        medium {
            set verbose 0
        }
        full {
            set verbose 1
        }
        debug {
            set verbose 1
            set debug 1
        }
        none {
            set verbose 0
            set quiet 1
        }
        default {
            chemerr "bad setting for list_option argument to psi4"
        }
    }
    
    # for unique_listing
    set listing_index 1
    
    
    # process dispersion correction
    if { "$dispersion_correction" != "undefined" } {
        # functional first verbose coords energy gradient
        dispersion_correction_calc $dispersion_correction 1 $verbose $coords $energy $gradient
    }
    
    # All variables including internal ones
    set psi4(variables) [ concat $psi4(in_args) \
            verbose debug quiet restart_psi4  listing_index ]
    
    # Store data for use at each point
    foreach var $psi4(variables) {
        set psi4($name,$var) [ set $var ]
    }
    
    # Create the new command
    proc $name { args } [ concat { return [ eval psi4_objcmd } $name {$args ] } ]
    
    end_module
    return $name
    
}

proc psi4_objcmd { name mode args } {
    
    global psi4
    
    # Recover variables for this object
    foreach var $psi4(variables) {
        set $var $psi4($name,$var)
    }
    
    # image>0 is used to keep track of DL-FIND NEB images,
    # so that separate wavefunction guesses etc. can be used for each image.
    set image 0
    
    # Parse any args provided
    # Note this will parse the image argument for DL-FIND NEB calcs
    if { [ parsearg psi4 $psi4(in_args) $args ] != 0 } then {
        chemerror "psi4 failed due to bad argument"
    }
    
    set gradflag 0
    set enerflag 0
    
    switch $mode {
        {energy}   {
            case $energy in {
                {undefined} {
                    chemerr " you must provide an energy= argument to psi4 "
                }
                {default} {
                    set enerflag 1
                }
            }
        }
        eandg -
        gradient {
            case $energy in {
                {undefined} {
                    chemerr " you must provide an energy= argument to psi4 "
                }
                {default} {
                    set enerflag 1
                }
            }
            case $gradient in {
                {undefined} {
                    chemerr " you must provide a gradient= argument to psi4 "
                }
                {default} {
                    set gradflag 1
                }
            }
        }
    }
    
    #  ================= Generate psi4 input file =================
    #
    
    push_banner_flag 0
    set fp [ open $jobname.inp w ]
    
    set natoms [ get_number_of_atoms coords=$coords ]
    #puts "natoms is $natoms"
    set nbq [ get_number_of_bqs coords=$coords ]
    set ncent  [ expr $natoms + $nbq]
    
    puts $fp "# psi4 input file generated by ChemShell"
    
    # Geometry
    puts $fp "molecule chemshellmol {"
        puts $fp "$charge $mult"
        for {set i 1} {$i <= $natoms} { incr i } {
            puts $fp [ get_atom_entry coords=$coords atom_number= $i ]
        }
        #Units, nosymmetry and no reorientation
        puts $fp "units au"
        puts $fp "symmetry c1"
        puts $fp "no_reorient"
        puts $fp "no_com"
        if { $efp == "yes" } {
            
            if { $efpatoms == "undefined" } {
                # No EFP region defined. Hence all MM atoms will be EFPs
                set efpregion undefined
                puts "unset"
                set countbq 0
                puts $fp "  --"
                puts $fp "  efp h2o"
                for { set w 1 } { $w <= $nbq } { incr w } {
                    if {$countbq == 3 } {
                        puts $fp "  --"
                        puts $fp "  efp h2o"
                        set countbq 0
                    }
                    set entry [ get_bq_entry coords=$coords bq_number=$w unit=au ]
                    puts $fp "[ lindex $entry 1 ] [ lindex $entry 2 ] [ lindex $entry 3 ]"
                    set countbq [expr $countbq + 1]
                }
            } else {
                # Here there is an EFP region. Different kind of loop
                set efpregion "yes"
                set countbq 0
                puts $fp "  --"
                puts $fp "  efp h2o"
                for { set w 0 } { $w < [llength $efpatoms] } { incr w } {
                    set cv [lindex $efpatoms $w]
                    if {$countbq == 3 } {
                        puts $fp "  --"
                        puts $fp "  efp h2o"
                        set countbq 0
                    }
                    set entry [ get_atom_entry coords=$fullcoords atom_number=$cv unit=au ]
                    puts $fp "[ lindex $entry 1 ] [ lindex $entry 2 ] [ lindex $entry 3 ]"
                    set countbq [expr $countbq + 1]
                }
            }
        }
        
    puts $fp "}"
    
    
    
    # Adding extra blocks from Chemshell inputfile
    puts $fp "$psi4input"
    
    #Setting
    if { $mult == 1} {
        puts $fp "set reference rks"
    }
    if { $mult > 1} {
        puts $fp "set reference uks"
    }
    
    #Preserve Checkpoint file
    puts $fp "psi4_io.set_specific_retention(PSIF_SCF_MOS,True)"
    #psi4_io.set_specific_path(PSIF_CHKPT, './')
    #psi4_io.set_specific_retention(PSIF_CHKPT, True)
    
    
    # If Point charges and efpregion is still undefined (meaning whole MM region should be EFP and hence no PC).
    if { $nbq > 0 && $efp != "yes"  } then {
        puts "Regular Elstat"
        puts $fp "Chrgfield = QMMM()"
        for { set i 1 } { $i <= $nbq } { incr i } {
            
            # Check if unit for pointcharge coordinates here is correct
            set entry [ get_bq_entry coords=$coords bq_number=$i unit=au ]
            puts $fp "Chrgfield.extern.addCharge( [ lindex $entry 4 ], [ lindex $entry 1 ], [ lindex $entry 2 ], [ lindex $entry 3 ] )"
        }
        puts $fp "psi4.set_global_option_python('EXTERN',Chrgfield.extern)"
    }
    
    # If both EFPs and Point charges. efpregion has to be set to yes
    if { $efp == "yes" } {
        set numqmatoms [get_number_of_atoms coords=$coords]
        if { $nbq > 0 && $efpregion == "yes" } then {
            puts "EFP + Elstat"
            puts $fp "Chrgfield = QMMM()"
            #puts "pcatoms is $pcatoms"
            for { set i 0 } { $i < [llength $pcatoms] } { incr i } {
                set pc [lindex $pcatoms $i]
                set pccharge [lindex $pccharges $i]
                set bqnum [expr $pc - $numqmatoms ]
                # puts "i is $i and pc is $pc and bqnum is $bqnum"
                #Taking PC atom coords using get_atom_entry of fullcoords
                set entry [ get_atom_entry coords=$fullcoords atom_number=$pc unit=au ]
                puts $fp "Chrgfield.extern.addCharge( $pccharge, [ lindex $entry 1 ], [ lindex $entry 2 ], [ lindex $entry 3 ] )"
            }
            puts $fp "psi4.set_global_option_python('EXTERN',Chrgfield.extern)"
        }
    }
    
    
    # "Simple input"
    if { $gradflag } then {
        if { $restart_psi4 } {
            puts $fp "set guess $read"
            puts $fp
            #Weirdly gradient restart does not work.
            #Have to do energy step first with restart, then maxiter 0
            puts $fp "energy('$method',restart_file='psiorbs.180')"
            puts $fp "set maxiter 0"
            puts $fp "set FAIL_ON_MAXITER false"
            puts $fp "gradient('$method')"
        } else {
            puts $fp "set guess $guess"
        puts $fp "gradient('$method')" }
    } else {
        
        if { $restart_psi4 } {
            puts $fp "set guess $guess"
            puts $fp "energy('$method',restart_file='psiorbs.180')"
        } else {
            puts $fp "set guess $guess"
        puts $fp "energy('$method')" }
        
    }
    
    
    # User-defined options
    if { $optstr != "undefined" } {
        puts $fp $optstr
    }
    
    puts $fp "copy_file_from_scratch('psiorbs.180', 'psi', 'chemshellmol', 180)"
    
    
    close $fp
    #
    # =========================== Run psi4 ===========================
    #
    #puts "Starting psi4 calculation"
    # Delete old output file if it exists from a previous run
    file delete $jobname.out
    if { $restart_psi4 == 0 } {
        # Discard restart information
        file delete $moinp
    }
    
    if {$image > 0 } {
        # Specify a file suffix to separate image files
        set imageid ".img$image"
        # Copy input file for actual image if present.
        if { [ file exists "$moinp$imageid" ] } {
            catch { file copy -force -- "$moinp$imageid" $moinp }
            # Check whether file has really been copied.
            set status [catch {exec diff "$moinp$imageid" $moinp} result]
            if {$status == 1} {
                puts "** WARNING, copy process failed !!! psi4 will try to restart from whatever is left in $moinp **"
            }
        }
    }
    
    set code 0
    set code [ exec_module psi4 "$executable $jobname.inp -a -o $jobname.out " NULL ]
    
    #  set orbfile psiorbs
    #  exec cp /tmp/psi.*.chemshellmol.180 .
    
    
    if { "$code" != "0" } then {
        puts stdout "error detected, code $code"
        error "psi4 failed"
    }
    
    #  file copy -force /tmp/psi.*.chemshellmol.180 .
    
    # if { $image > 0 } {
    #      catch { file copy -force -- $jobname.gbw "$moinp$imageid" }
    #      # Check whether file has really been copied back.
    #      set status [catch {exec diff $jobname.gbw "$moinp$imageid"} result]
    #           if {$status == 1} {
    #           puts "** WARNING, process to copy back the actual wave function failed !!! "
    #           puts "This might lead to psi4 failing in the next cycle for this image **"
    #	   }
    #  }
    
    # restart on future calcs
    set psi4($name,restart_psi4) 1
    
    #
    # ======================= Read psi4 output =======================
    #
    if { $gradflag } then {
        puts "Still have not coded the Psi4 gradient. Also pc gradient not possible. RB"
        exit
        set fpout [ open $jobname.out r ]
        if { $enerflag } then {
            # search for "Total Energy =" in $jobname.engrad
            set code [ gets $fpout line ]
            while { [string first "Total Energy =" $line] == -1 && $code != -1 } {
                #puts "line is $line"
                set code [ gets $fpout line ]
            }
            if { $code == -1 } then { chemerr "Couldn't find energy in psi4 output file" }
            set val [ lindex $line 3 ]
            set_matrix_size datatype=real dimensions= "1 1" matrix=$energy name="energy"
            set_matrix_element indices= "0 0" value = $val matrix=$energy
        }
        # search for "Gradient info"
        #puts "We are searching $jobname.out for gradient"
        set code [ gets $fpout line ]
        while { $line != "  -Total Gradient:" && $code != -1 } {
            set code [ gets $fpout line ]
            #puts "line is $line"
        }
        if { $code == -1 } then { chemerr "Couldn't find gradient in psi4 .engrad file" }
        # skip a line
        set code [ gets $fpout line ]
        set code [ gets $fpout line ]
        set code [ gets $fpout line ]
        # read in data
        set_matrix_size datatype=real dimensions= "3 $ncent" matrix=$gradient name= "gradient"
        for {set i 0} {$i < $natoms} { incr i } {
            #puts "line is $line"
            set val [lindex $line 1]
            set_matrix_element indices= "0 $i" value = $val matrix=$gradient
            set val [lindex $line 2]
            set_matrix_element indices= "1 $i" value = $val matrix=$gradient
            set val [lindex $line 3]
            set_matrix_element indices= "2 $i" value = $val matrix=$gradient
            set code [ gets $fpout line ]
        }
        #puts "Gradient found is:"
        #print_matrix matrix=$gradient
        close $fpout
        # get point charge gradient if there is one
        if { $nbq != 0 } then {
            set fpout [ open $jobname.pcgrad r ]
            gets $fpout line
            set val [lindex $line 0]
            if { $val != $nbq } { chemerr "Inconsistent nbq in psi4 .pcgrad file" }
            set ncent [ expr $natoms + $nbq ]
            for { set i $natoms } { $i < $ncent } { incr i } {
                gets $fpout line
                set val [lindex $line 0]
                set_matrix_element indices= "0 $i" value = $val matrix=$gradient
                set val [lindex $line 1]
                set_matrix_element indices= "1 $i" value = $val matrix=$gradient
                set val [lindex $line 2]
                set_matrix_element indices= "2 $i" value = $val matrix=$gradient
            }
            close $fpout
        }
    } elseif { $enerflag } then {
        set fpout [ open $jobname.out r ]
        # search for "FINAL SINGLE POINT ENERGY" in $jobname.out
        # (.engrad has more sig figs but is only written out for gradient calcs)
        set code [ gets $fpout line ]
        while { [ string first "Total Energy =" $line ] == -1  && $code != -1 } {
            set code [ gets $fpout line ]
        }
        if { $code == -1 } then { chemerr "Couldn't find energy in psi4 output file" }
        set val [ lindex $line 3 ]
        set_matrix_size datatype=real dimensions= "1 1" matrix=$energy name="energy"
        set_matrix_element indices= "0 0" value = $val matrix=$energy
        close $fpout
    }
    
    # save i/o files if requested
    if { $unique_listing } {
        file rename -force $jobname.inp $jobname.inp.$listing_index
        file rename -force $jobname.out $jobname.out.$listing_index
        incr psi4($name,listing_index)
    }
    
    pop_banner_flag
    
    # process dispersion_correction
    if { "$dispersion_correction" != "undefined" } {
        # functional first verbose coords energy gradient
        dispersion_correction_calc $dispersion_correction 0 $verbose $coords $energy $gradient
    }
    
    end_module
    
    return 0
    
}

#
# Procedural interfaces
#
proc psi4.init args {
    global psi4
    new_psi4 psi41 $args
    return 0
}
proc psi4.energy { args } {
    return [ eval psi41 energy $args ]
}
proc psi4.gradient { args } {
    return [ eval psi41 gradient $args ]
}
proc psi4.eandg { args } {
    return [ eval psi41 eandg $args ]
}
proc psi4.update { args } {
    return 0
}
proc psi4.kill args {
    return 0
}

proc psi42.init args {
    global psi4
    new_psi4 psi42 $args
    return 0
}
proc psi42.energy { args } {
    return [ eval psi42 energy $args ]
}
proc psi42.gradient { args } {
    return [ eval psi42 gradient $args ]
}
proc psi42.eandg { args } {
    return [ eval psi42 eandg $args ]
}
proc psi42.update { args } {
    return 0
}
proc psi42.kill args {
    return 0
}

proc psi43.init args {
    global psi4
    new_psi4 psi43 $args
    return 0
}
proc psi43.energy { args } {
    return [ eval psi43 energy $args ]
}
proc psi43.gradient { args } {
    return [ eval psi43 gradient $args ]
}
proc psi43.eandg { args } {
    return [ eval psi43 eandg $args ]
}
proc psi43.update { args } {
    return 0
}
proc psi43.kill args {
    return 0
}

