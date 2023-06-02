# Cuts all $clusters from $startloc to radius $size angstroms and puts them into $endloc
proc clustercut {clusters startloc endloc size} {
    global molsolventdir
    source $molsolventdir/proclist.tcl
    set sizebohr [expr $size*$angtobohr]

    foreach i $clusters {
        puts "Now cutting cluster $i to $size angstrom radius"
        push_banner_flag 0
        set addcoords {}
        # Getting info from original cluster
        set frag $startloc/$i
        fragment $frag old persistent
        set numatoms [get_number_of_atoms coords=$frag]
        set centre [get_molecule_centre coords=$frag]
        # Get all the coordinates data and turn into list, cut to only include coords
        set coords [lrange [split [string trim [c_prepare_input coords=$frag]] "\n"] 3 end]
        # Go through list and find coordinates that are within range
        for {set m 1} {$m <= $numatoms} {incr m [llength $mollist]} {
            # Get list of atoms in molecule (assumes proper order but should be safe)
            set write 0
            set mollist [get_molecule_members coords=$frag atom_number=$m]
            # Go through each atom in mollist and find distance to centre, set write flag
            foreach j $mollist {
                set atomcoords [lrange [lindex $coords [expr $j-1]] 1 end]
                if {[distance $atomcoords $centre] <= $sizebohr} {
                    set write 1
                    break
                }
            }
            # Add atoms to endcoords list if write flag is on
            if {$write} {
                foreach k $mollist {
                    lappend addcoords [lindex $coords [expr $k-1]]
                }
            }
        }
        set writecoords [join $addcoords "\n"]
        # Making blank file to hold new coordinates and getting info
        c_create coords=temp {}
        set fp [open temp r]
        set origdata [read -nonewline $fp]
        close $fp
        # Adding original block stuff and coords
        set fp [open $endloc/$i w]
        puts $fp "$origdata"
        puts $fp "block = coordinates records = [llength $addcoords]"
        puts $fp "$writecoords"
        close $fp
        # Making connectivity block (maybe not stritcly necessary)
        connect coords=$endloc/$i conn=$endloc/$i
        pop_banner_flag
   }
}
