# This procedure splits a list evenly over the workgroups
proc gensnaplist {snapshots nwg wid} {
    set numsnaps [llength $snapshots]
    set snaptasks [expr $numsnaps/$nwg]
    set remainder [expr $numsnaps%$nwg]
    set start [expr $wid*$snaptasks+min($wid,$remainder)]
    set stop [expr $start+$snaptasks-1]
    if {$wid<$remainder} {
        set stop [expr $stop+1]
    }
    set range [seq $start to $stop]
    foreach i $range {
        lappend snaplist [lindex $snapshots $i]
    }
    foreach i $snaplist {
        dict lappend snapdict [expr [lsearch $snapshots $i]+1] $i
    }
    puts "snapdict is $snapdict"
    return $snapdict
}
