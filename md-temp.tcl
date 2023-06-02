# Sourcing vars from input file
global title nstep mdrestart mdrandomize pairupd writetra writerstrt istep fprot ohfile

# Molecular dynamics procs
proc mdstep { } {
    global istep fprot writetra writerstrt pairupd stepb stepc stepd stepe tempb tempc tempd tempe snapslist
    if { [ expr $istep % $pairupd ] == 0 } { dyn1 update }
    if { ($istep == 1) } { dyn1 output }
    dyn1 force
    dyn1 step
    dyn1 printe
    #if { [expr $istep%10] == 0 } { write_hbonds hybrid.dl_poly.coords }

    if { ! ($istep % $writerstrt) } { dyn1 dump }
    if { ! ($istep % $writetra) } { writetra full }
    if { ($istep == $stepb) } { dyn1 configure temperature = $tempb ensemble=NVT }
    if { ($istep == $stepc) } { dyn1 configure temperature = $tempc ensemble=NVT }
    if { ($istep == $stepd) } { dyn1 configure temperature = $tempd ensemble=NVT }
    if { ($istep == $stepe) } { dyn1 configure temperature = $tempe ensemble=NVT }
    if { 0 <= [lsearch $snapslist $istep] } {grabsnap}
    puts $fprot [ format "%8d %14.4f %12.6f %12.6f %12.6f %14.2f %13.5e %12.3e " \
                      $istep \
                      [ dyn1 get time ] \
                      [ dyn1 get potential_energy ] \
                      [ dyn1 get kinetic_energy ] \
                      [ dyn1 get total_energy ] \
                      [ dyn1 get temperature ] \
                      [ dyn1 get constraint_force ] \
                      [ dyn1 get friction ] ]
    flush $fprot
    return 0
}

proc grabsnap { } {
	global istep snapshotsA snapshotsB state
	puts "Grabbing snapshot"
	# dyn1.tempc is the temporary coordinate file
	flush_object dyn1.tempc
	exec cp dyn1.tempc snap$state-$istep.c
	lappend snapshots$state snap$state-$istep.c
}

proc writetra { {type full} } {
    switch $type {
        full {
            # The standard ChemShell trajectory
            dyn1 trajectory
            puts "--- ChemShell trajectory written."
        }
        none { }
        default {
            puts "*** Invalid trajectory type: $type."
            puts "*** No trajectory written."
            return 1
        }
    }
    return 0
}

# MD Loop
# Organize restart, protocol, etc.
if { $mdrestart == "T" } {
    dyn1 load
    # Get the last step no. from protocol (if it exists)
    if { [ file exists $title.prot ] && [ file writable $title.prot ] } {
		set fprot [ open "|tail -1 $title.prot" r ]
		gets $fprot line
		close $fprot
		scan $line %d last
		set first [ expr $last + 1 ]
		set fmode a
    } else {
        set first 1
        set fmode w
    }
} else {
    # No restart
    set first 1
    set fmode w
}

set fprot [ open $title-$state.prot $fmode ]
puts $fprot \
{  #Step        time[ps]    E_pot[au]    E_kin[au]    E_tot[au]           T[K]     F_c[au]       Fric     }

# Remove existing EXIT file
if [ file exists EXIT ] { exec rm EXIT }

# Initial velocity distribution
if { $mdrandomize == "T" } { dyn1 initvel }

#################
# SIMULATION LOOP
#################

for { set istep $first } { $istep < $first + $nstep } { incr istep } {
    mdstep
    # Allow graceful termination via EXIT file (1 if file exists)
    set exitflag [ file exists EXIT ]
    if { $exitflag } {
        puts "*** EXIT file detected. Terminating..."
        exit
    }
}
