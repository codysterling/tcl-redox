proc setqmtheory {x} {
    global charge$x mult$x orcapath orcainput orcablocks xtbpath gasorcainput gasorcablocks estimorcainput estimorcablocks polmmorcainput polmmorcablocks restA restB

    set charge [set charge$x]
    set mult [set mult$x]
    set rest [set rest$x]

    set orcatheory [ list \
        theory=orca : [ list \
            executable=$orcapath/orca \
            charge=$charge \
            mult=$mult \
            orcasimpleinput= $orcainput \
            orcablocks= $orcablocks ]]
    set gasorcatheory [ list \
        theory=orca : [ list \
            executable=$orcapath/orca \
            charge=$charge \
            mult=$mult \
            orcasimpleinput= $gasorcainput \
            orcablocks= $gasorcablocks ]]
    set estimorcatheory [ list \
        theory=orca : [ list \
            executable=$orcapath/orca \
            charge=$charge \
            mult=$mult \
            orcasimpleinput= $estimorcainput \
            orcablocks= $estimorcablocks ]]
    set polmmorcatheory [ list \
        theory=orca : [ list \
            executable=$orcapath/orca \
            charge=$charge \
            mult=$mult \
            restart=$rest \
            orcasimpleinput= $polmmorcainput \
            orcablocks= $polmmorcablocks ]]
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
    variable gasorcatheory$x $gasorcatheory
    variable estimorcatheory$x $estimorcatheory
    variable polmmorcatheory$x $polmmorcatheory
    variable om3theory$x $om3theory
    variable pm3theory$x $pm3theory
    variable xtbtheory$x $xtbtheory
}
