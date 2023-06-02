proc qmmm-md {x} {
    global molsolventdir molsys$x solvsys$x con$x groups$x froz$x masses$x restraints$x solutelist$x types$x radius$x tempa dyncode charges$x coupling
    global ${dyncode}theory$x

    variable state $x
    set molsys [set molsys$x]
    set solvsys [set solvsys$x]
    set con [set con$x]
    set groups [set groups$x]
    set froz [set froz$x]
    set masses [set masses$x]
    set restraints [set restraints$x]
    set solutelist [set solutelist$x]
    set types [set types$x]
    set radius [set radius$x]
    set theory [set ${dyncode}theory$x]
    set charges [set charges$x]

    puts "Running $dyncode QM/MM MD trajectory for state $x"

    dynamics dyn1 coords=$solvsys \
        ensemble = NVT \
        temperature = $tempa \
        nosehoover=4 \
        taut=.02 \
        energy_unit = "kcal mol-1" \
        timestep = .001 \
        constraints = $con \
        groups = $groups \
        frozen = $froz \
        masses = $masses \
        restraints = $restraints \
        theory=hybrid : [ list \
            cutoff=1000 \
            groups = $groups \
            qm_region= $solutelist \
            atom_charges = $charges \
            coupling=$coupling \
            list_option=none \
            qm_$theory \
            mm_theory=dl_poly : [ list \
                list_option=none \
                exact_srf=on \
                cutoff=1000 \
                no_elec=no \
                use_pairlist=no \
                atom_types= $types \
                mxlist=15000 \
                mxexcl=15000 \
                mm_defs=md-FF$x.ff ]]

    global ohfile
    set ohfile [open hbonds-$x "w"]
    puts $ohfile "# O#     Odist      H#     Hdist       O* H*"

    dyn1 trajectory
    source $molsolventdir/md-hbonds-temp.tcl
    exec cp dynamics.trj trajectory$x.trj
    dyn1 destroy
    close $ohfile

    # Moving snapshot .c files to snaps folder
    puts "State $x MD done; moving snaphot files to snaps folder"
    exec mkdir -p snaps
    exec mv {*}[glob snap$x*.c] snaps

    # Moving pre-MD solute geometries to snaps folder and deleting objects
    exec mv $molsys snaps/gas-mol$x.c
#    delete_object $molsys
    delete_object $solvsys
}

proc permm-md {x} {
    global molsolventdir molsys$x solvsys$x con$x solutelist$x types$x tempa charges$x masses$x cutradius

    variable state $x
    set molsys [set molsys$x]
    set solvsys [set solvsys$x]
    set con [set con$x]
    set solutelist [set solutelist$x]
    set types [set types$x]
    set charges [set charges$x]
    set masses [set masses$x]

    puts "Running MM MD trajectory for state $x"

    dynamics dyn1 \
        coords=$solvsys \
        temperature = $tempa \
        ensemble = NVT \
        energy_unit= "kcal mol-1" \
        list_option=medium \
        timestep = 0.001 \
        frozen = $solutelist \
        constraints = $con \
        masses= $masses \
        theory=dl_poly : [ list \
            list_option=full \
            spme=yes \
            use_pairlist=yes \
            conn=$solvsys \
            atom_types= $types \
            atom_charges= $charges \
            mm_defs=md-FF$x.ff ]

    global ohfile
    set ohfile [open hbonds-$x "w"]
    puts $ohfile "# O#     Odist      H#     Hdist       O* H*"

    dyn1 trajectory
    source $molsolventdir/md-temp.tcl
    exec cp dynamics.trj trajectory$x.trj
    dyn1 destroy
    close $ohfile

    # Cutting periodic system into spheres
    push_banner_flag 0
    exec mkdir -p snaps

    global snapshots$x soluteatomcentre$x
    set snapshots [set snapshots$x]
    set soluteatomcentre [set soluteatomcentre$x]

    puts "Now cutting snapshots into spheres with radius $cutradius bohr"

    foreach i $snapshots {
        set snap [string trimright $i .c]
        cluster_cut coords=$i radius_cluster=$cutradius cluster=$snap-cluster.c crystal_type=molecular origin_atom= $soluteatomcentre
        unset_charges coords=$snap-cluster.c
        exec mv $snap-cluster.c snaps/$snap.c
    }

    puts "Spherical clusters are cut and moved to snaps folder"
    pop_banner_flag

    # Moving pre-MD solute geometries to snaps folder and deleting solvsys objects
    exec mv $molsys snaps/gas-mol$x.c
    delete_object $molsys
    delete_object $solvsys
}
