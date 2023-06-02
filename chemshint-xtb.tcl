#
# xtb  interface. Based on modified ORCA interface
#

proc new_xtb { name args } {
    global xtb
    # Option defaults
    set executable  xtb
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
    set xtbupdatecharges undefined
    # xtb input options
    set xtb(in_args) { \
                coords energy gradient hessian \
                executable jobname restart moinp \
                charge mult functional \
                gridsize finalgridsize use_finalgrid
        list_option unique_listing \
                dispersion_correction \
                gridsize finalgridsize use_finalgrid \
                convergence guess nproc \
                maxcyc optstr method\
                excited estate eroots tda image  \
            }
    # All variables including internal ones
    set xtb(variables) [ concat $xtb(in_args) \
            verbose debug quiet restart_xtb  listing_index ]
    
    # Assign name if not given
    if { "$name" == "*" } {
        if { [ info exists xtb(index)  ] } {
            set index $xtb(index)
            incr xtb(index)
        } else {
            set index 1
            set xtb(index) 2
        }
        set name xtb$index
    }
    
    if { [ parsearg xtb.init $xtb(in_args) $args ] != 0 } {
        chemerror "xtb failed due to bad argument"
    }
    
    set restart_xtb [ check_boolean restart ]
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
            chemerr "bad setting for list_option argument to xtb"
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
    set xtb(variables) [ concat $xtb(in_args) \
            verbose debug quiet restart_xtb  listing_index ]
    
    # Store data for use at each point
    foreach var $xtb(variables) {
        set xtb($name,$var) [ set $var ]
    }
    
    # Create the new command
    proc $name { args } [ concat { return [ eval xtb_objcmd } $name {$args ] } ]
    
    end_module
    return $name
    
}

proc xtb_objcmd { name mode args } {
    
    global xtb
    
    # Recover variables for this object
    foreach var $xtb(variables) {
        set $var $xtb($name,$var)
    }
    
    # image>0 is used to keep track of DL-FIND NEB images,
    # so that separate wavefunction guesses etc. can be used for each image.
    set image 0
    
    # Parse any args provided
    # Note this will parse the image argument for DL-FIND NEB calcs
    if { [ parsearg xtb $xtb(in_args) $args ] != 0 } then {
        chemerror "xtb failed due to bad argument"
    }
    
    set gradflag 0
    set enerflag 0
    set hessflag 0
    switch $mode {
        {energy}   {
            case $energy in {
                {undefined} {
                    chemerr " you must provide an energy= argument to xtb "
                }
                {default} {
                    set enerflag 1
                }
            }
        }
        {hessian}   {
            case $hessian in {
                {undefined} {
                    chemerr " you must provide a hessian= argument to xtb "
                }
                {default} {
                    set hessflag 1
                }
            }
        }
        eandg -
        gradient {
            case $energy in {
                {undefined} {
                    chemerr " you must provide an energy= argument to xtb "
                }
                {default} {
                    set enerflag 1
                }
            }
            case $gradient in {
                {undefined} {
                    chemerr " you must provide a gradient= argument to xtb "
                }
                {default} {
                    set gradflag 1
                }
            }
        }
    }
    
    #  ================= Generate xtb input file =================
    # Very simple because there is no inputfile. Just need to create xyz coordinate file
    #Point charges are automatically read if pcharge file is present in dir.
    push_banner_flag 0
    #set fp [ open $jobname.xyz w ]
    write_xyz coords=$coords file=$jobname.xyz
    set natoms [ get_number_of_atoms coords=$coords ]
    #set bohrang 1.88972598857892
    set nbq [ get_number_of_bqs coords=$coords ]
    set ncent  [ expr $natoms + $nbq]
    
    # Point charge file. Called pcharge in xtb, same format as ORCA
    global bohrtoang
    file delete pcharge
    
    if { $nbq > 0 } then {
        format_xtb_bq_list file pcharge
    }
    
    # =========================== Run xtb ===========================

    # Delete old output file if it exists from a previous run
    file delete gradient
    file delete energy
    if { $restart_xtb == 0 } {
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
                puts "** WARNING, copy process failed !!! xtb will try to restart from whatever is left in $moinp **"
            }
        }
    }
    
    set code 0
    
    set unpel [expr $mult - 1]
    #energy+gradient
    
    # Running various types of xTB based on inputs
    if {$hessflag == 1} {
        set code {exec xtb $jobname.xyz -gfn -chrg $charge -uhf $unpel -hess > $jobname.out}
    } elseif { $gradflag } {
        set code {exec xtb $jobname.xyz -gfn -chrg $charge -uhf $unpel -grad > $jobname.out}
        set code3 {exec xtb $jobname.xyz -gfn -chrg $charge -uhf $unpel -grad -etemp 20000 > $jobname.out}
    } else {
        set code {exec xtb $jobname.xyz -gfn -chrg $charge -uhf $unpel > $jobname.out}
    }

    # This runs the above functions and catches errors if they exist
#    global istep
#    if {$istep == 1} {
#        puts "Doing special 20000K for step 1"
#        catch $code3 err3
#        file delete energy
#        file delete gradient
#        puts "Now doing step 1 at 300K"
#        catch $code err
#    } else {
        catch $code err
#    }

    # Here we check if xtb throws "SCC not converged" error and do high-temp/low-temp calculation to fix
    while { "$err" == "SCC not converged" } {
        puts "SCC not converged, redoing $code with high etemp (4000K)"
        file delete xtbrestart
        if {$hessflag == 1} {
            set code2 {exec xtb $jobname.xyz -gfn -chrg $charge -uhf $unpel -hess -etemp 4000 > $jobname.out}
        } elseif { $gradflag } {
            set code2 {exec xtb $jobname.xyz -gfn -chrg $charge -uhf $unpel -grad -etemp 4000 > $jobname.out}
        } else {
            set code2 {exec xtb $jobname.xyz -gfn -chrg $charge -uhf $unpel -etemp 4000 > $jobname.out}
        }
        catch $code2 err2
        file delete energy
        file delete gradient
        puts "Now doing again with default etemp (300K)"
        catch $code err
    }

#    if { $gradflag } then {
#        if {$unpel > 0 } {
#            set code [ exec_module xtb "$executable $jobname.xyz -gfn -chrg $charge -uhf $unpel -grad > $jobname.out " NULL ]
#            puts "**now finished xtb 0, code is $code"
#        } else {
#            set code [ exec_module xtb "$executable $jobname.xyz -gfn -chrg $charge -grad > $jobname.out " NULL ]
#            puts "**now finished xtb 1, code is $code"
#        }
#    } else {
#        if {$unpel > 0 } {  
#            set code [ exec_module xtb "$executable $jobname.xyz -gfn -chrg $charge -uhf $unpel > $jobname.out " NULL ]
#        } else {    
#            set code [ exec_module xtb "$executable $jobname.xyz -gfn -chrg $charge > $jobname.out " NULL ]
#        }
#    }
#    
#    if { $hessflag == 1 } then {
#        if {$unpel > 0 } {
#            set code [ exec_module xtb "$executable $jobname.xyz -gfn -chrg $charge -uhf $unpel -hess > $jobname.out " NULL ]
#        } else {
#            set code [ exec_module xtb "$executable $jobname.xyz -gfn -chrg $charge -hess > $jobname.out " NULL ]
#        }
#    }
    
    if { [string length $err] } then {
        puts stdout "error detected, error: $err"
        error "xtb failed: $err"
    }

    # Deleting xtbrestart afterwards
#   puts "now deleting xtbrestart"
#   file delete xtbrestart
    
    # restart on future calcs
    set xtb($name,restart_xtb) 1
    
    # ======================= Read xtb output =======================
    #
    if { $hessflag == "1"} {
        puts "Now finding xTB Hessian in output"
        set hesssize [expr $natoms *3]
        puts "hesssize is $hesssize"
        set_matrix_size datatype=real dimensions= "$hesssize $hesssize" matrix=$hessian name= "hessian"
        
        set hfpout [ open hessian r ]
        set file_data [read $hfpout]
        close $hfpout
        set data [split $file_data "\n"]
        set data [lreplace $data 0 0]
        set counthess 0
        set counter 0
        set i 0
        foreach line $data {
            if {[llength $line] == 0} {
                break
            }
            set line [lindex $data $counter]
            set linelength [llength $line]
            for {set k 0 } {$k < $linelength} {incr k 1} {
                set val [lindex $line $k]
                #puts "counthess is $counthess and i is $i and val is $val"
                set_matrix_element indices= "$counthess $i" value = $val matrix=$hessian
                set counthess [expr $counthess + 1]
            }
            set counter [expr $counter + 1]
            if {$counthess == [expr $hesssize - 0] } {
                #puts "reset"
                set counthess 0
                set i [expr $i + 1]
            }
            
        }
        puts "Printing xTB Hessian:"
        print_matrix matrix=$hessian
    }
    
    if { $gradflag } then {
        #xtb_read arguments: natoms enerflag energymatrix gradientflag gradientmatrix nbq
        xtb_read $natoms 1 $energy 1 $gradient $nbq
    } else {
        # Otherwise we read from energy file directly for single point calculations
        set_matrix_size datatype=real dimensions= "1 1" matrix=$energy name= "energy"
        set enfile [open energy r]
        set filedata [split [read -nonewline $enfile] "\n"]
        set curr_energy [lindex [lindex $filedata 1] 1]
        set_matrix_element indices= "0 0" value= $curr_energy matrix=$energy
    }

#    puts "energy matrix here is:"
#    print_matrix matrix=$energy
   #puts "Gradient after pcgrad"
    #print_matrix matrix=$gradient format=%12.9f
    
    # If QM/MM mechanical embedding take the xtb atom charges and update.
    # Only if xtbupdatecharges is active
    #if {$xtbupdatecharges == yes } {
    #set xtboutfile [ open xtb1.out r ]
    #set grab undefined
    #set count undefined
    #global cm5charges
    #set cm5charges {}
    #while {[gets $xtboutfile line]>=0} {
    #    if {$grab == yes && $count <= $natoms } {
    #    set linelist [lreplace $line 0 -1]
    #    lappend cm5charges [lindex $linelist 3]
    #    set count [expr $count + 1]
    #    }
    #    if {[lsearch $line Mulliken/CM5]>=0 } {
    #    set grab yes
    #    set count 1
    #    }
    #
    #}
    #
    #}
    
    # save i/o files if requested
    # if { $unique_listing } {
    #   file rename -force $jobname.inp $jobname.inp.$listing_index
    #   file rename -force $jobname.out $jobname.out.$listing_index
    #   incr xtb($name,listing_index)
    # }
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

proc xtb.hess args {
    return [ eval xtb1 hessian $args ]
}

proc xtb.init args {
    global xtb
    new_xtb xtb1 $args
    return 0
}
proc xtb.energy { args } {
    return [ eval xtb1 energy $args ]
}
proc xtb.gradient { args } {
    return [ eval xtb1 gradient $args ]
}
proc xtb.eandg { args } {
    return [ eval xtb1 eandg $args ]
}
proc xtb.update { args } {
    return 0
}
proc xtb.kill args {
    return 0
}

proc xtb2.init args {
    global xtb
    new_xtb xtb2 $args
    return 0
}
proc xtb2.energy { args } {
    return [ eval xtb2 energy $args ]
}
proc xtb2.gradient { args } {
    return [ eval xtb2 gradient $args ]
}
proc xtb2.eandg { args } {
    return [ eval xtb2 eandg $args ]
}
proc xtb2.update { args } {
    return 0
}
proc xtb2.kill args {
    return 0
}

proc xtb3.init args {
    global xtb
    new_xtb xtb3 $args
    return 0
}
proc xtb3.energy { args } {
    return [ eval xtb3 energy $args ]
}
proc xtb3.gradient { args } {
    return [ eval xtb3 gradient $args ]
}
proc xtb3.eandg { args } {
    return [ eval xtb3 eandg $args ]
}
proc xtb3.update { args } {
    return 0
}
proc xtb3.kill args {
    return 0
}

