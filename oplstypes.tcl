proc setoplstypes {x oplsfile} {
    global molsolventdir

    upvar 1 solutetypes$x solutetypes

    # First have to wipe solutetypes so we can write the appropriate OPLS ones
    # This is done where each atom has its own type even though there will be redundancy here
    # However it should automatically preserve order for types list and stuff
    set solutetypes {}

    # Opening .prm file and looking through each line for data
    set fp [open $oplsfile r]
    set write 0
    foreach line [split [read -nonewline $fp] "\n"] {
        # This ensures we don't start taking data until after "NONBONDED" line of file
        if {[lindex $line 0]=="NONBONDED"} {
            set write 1
        }
        # This ensures we only write from the atom information block
        if {$write && [llength $line] == 7} {
            lappend solutetypes [lindex $line 0]
        }
    }

    puts "After reading OPLS file $oplsfile for molecule $x, new solutetypes ([llength $solutetypes]) are: $solutetypes"
}
