proc setqmtheory {x} {
    global charge$x mult$x orcapath orcainput orcablocks xtbpath

    set charge [set charge$x]
    set mult [set mult$x]

    # These theories can be used in MD input file
    set orcatheory [ list \
        theory=orca : [ list \
            executable=$orcapath/orca \
            charge=$charge \
            mult=$mult \
            orcasimpleinput= $orcainput \
            orcablocks= $orcablocks \
            optstr= "%output Print\[ P_Hirshfeld \] 1 end" ]]
    set om3theory [ list \
        theory=mndo : [ list \
            hamiltonian=om3 \
            maxcyc=1000 \
            charge=$charge \
            mult=$mult ]]
    set pm3theory [ list \
        theory=mndo : [ list \
            hamiltonian=pm3 \
            maxcyc=1000 \
            charge=$charge \
            mult=$mult ]]
    set xtbtheory [ list \
        theory=xtb : [ list \
            executable=$xtbpath/xtb \
            method=gfn \
            charge=$charge \
            mult=$mult ]]

    variable orcatheory$x $orcatheory
    variable om3theory$x $om3theory
    variable pm3theory$x $pm3theory
    variable xtbtheory$x $xtbtheory
}
