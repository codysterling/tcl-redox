proc repsnaplist {type num snapslist endict} {
    global molsolventdir 

    source $molsolventdir/proclist.tcl

    # Printing stuff
    puts "type is $type"
    puts "num is $num"
    puts "snapslist is $snapslist"
    puts "endict is $endict"

    # Calculate average of full dictionary
    set avgip [average [dict values $endict]]
    puts "avgip is $avgip"

    # Make lower-than- and greater-than-average dictionaries
    # Key = snapshot num (taken from snapslist), value = difference from avg (always positive)
    foreach key [dict keys $endict] {
        set snapshot [lindex $snapslist [expr $key - 1]]
        set deltae [expr [dict get $endict $key] - $avgip]
        if { $deltae < 0 } {
            dict append ltdict $snapshot [expr -$deltae]
        } elseif { $deltae > 0} {
            dict append gtdict $snapshot $deltae
        }
    }
    # Now we append and make full list
    set fulldict [dict merge $ltdict $gtdict]
    # Now we order the dictionaries so they are increasing difference from average
    set ltdictsort [dictvalsort $ltdict]
    set gtdictsort [dictvalsort $gtdict]
    set fulldictsort [dictvalsort $fulldict]

    # Here we use the $type and $num arguments to find the correct representative snapshots
    if {$type == "num"} {
        # Just need to get first $num snapshots from fulldict
        set replist [lrange [dict keys $fulldictsort] 0 [expr $num - 1]]
    } elseif {$type == "numalt"} {
        # Here need to figure out which list to start with based on first element
        # We are only interested in keys (snapshot numbers)
        if {[lindex $ltdictsort 1] < [lindex $gtdictsort 1]} {
            set list1 [dict keys $ltdictsort]
            set list2 [dict keys $gtdictsort]
        } elseif {[lindex $gtdictsort 1] < [lindex $ltdictsort 1]} {
            set list1 [dict keys $gtdictsort]
            set list2 [dict keys $ltdictsort]
        } else {
            puts "Something is really wrong with numalt sorting, exiting."
            exit
        }
        # Now we loop through these dictionaries to find representative snapshots
        set iter 1
        set list1index 0
        set list2index 0
        while {$iter <= $num} {
            # If iteration number is odd, use list1
            if {[expr $iter%2]==1} {
                lappend replist [lindex $list1 $list1index]
                incr list1index
            } else {
            # Otherwise use list2
                lappend replist [lindex $list2 $list2index]
                incr list2index
            }
            incr iter
        }
    } elseif {$type == "cut"} {
        # Here we go through fulldictsort and add all keys (snapshots) where the value is below cutoff
        dict for {key val} $fulldictsort {
            if {$val <= $num} {
                lappend replist $key
            } else {
                break
            }
        }
        # Here we make sure that there are actually snapshots in this cutoff or exit
        if {! [info exists replist]} {
            puts "There are no values in representative list, exiting."
            exit
        }
    } else {
        puts "Type "$type" is an invalid keyword for finding representative snapshots, exiting."
        exit
    }

    # Returning the replist
    return $replist
}
