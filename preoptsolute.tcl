# Pre-optimizes solute geometries using whatever theory from mdcode (defined in md.chmb)

proc preopt {x} {
    fragment optmol$x.c persistent

    global molsys$x optcode
    global ${optcode}theory$x

    set qmtheory [set ${optcode}theory$x]

    dl-find coords=mol$x.c $qmtheory result=optmol$x.c coordinates=dlc maxcycle=700 maxene=300 tolerance=0.00045 \

    set molsys$x optmol$x.c
}
