source $molsolventdir/proclist.tcl

if {"$solvtype"=="spc"} {
    # Setting point charges for water model
	set OTcharge -.82
	set HTcharge .41

    # Setting LJ parameters for oxygen
    # Values correspond to function for 4e(()^12-()^6)
    # These are fixed for Chemshell use later in ffgen.tcl
    # Units: sigma = A	epsilon = kcal/mol
    set solvsig [list {OT 3.1655}]
    set solveps [list {OT .1554}]

    # Sourcing the solvent system and definitions
    if {"$solvgeom"=="box"} {
        set solvent $molsolventdir/solventlibrary/spcbox.c
    } elseif {"$solvgeom"=="sphere"} {
    	if {"$solvsize"=="9k"} {
            set solvent $molsolventdir/solventlibrary/spc9kdrop.c
        } elseif {"$solvsize"=="big"} {
	    	set solvent $molsolventdir/solventlibrary/spcbigdrop.c
    	} elseif {"$solvsize"=="small"} {
    		set solvent $molsolventdir/solventlibrary/spcsmalldrop.c
    	} elseif {"$solvsize"=="tiny"} {
            set solvent $molsolventdir/solventlibrary/spctinydrop.c
        }
    }
}

if {"$solvtype"=="tip3p" || "$solvtype"=="tips3p"} {
    # Setting point charges for water model
	set OTcharge -.834
	set HTcharge .417

    # Setting LJ parameters for oxygen
    # Values correspond to function for 4e(()^12-()^6)
    # These are fixed for Chemshell use later in ffgen.tcl
    # Units: sigma = A	epsilon = kcal/mol
    set solvsig [list {OT 3.15066}]
    set solveps [list {OT .1521}]
    if {"$solvtype"=="tips3p"} {
        lappend solvsig {HT .40001}
        lappend solveps {HT .0460}
    }

    # Sourcing the solvent system and definitions
    if {"$solvgeom"=="box"} {
        set solvent $molsolventdir/solventlibrary/tip3pbox.c
    } elseif {"$solvgeom"=="sphere"} {
        if {"$solvsize"=="9k"} {
            set solvent $molsolventdir/solventlibrary/tip3p9kdrop.c
        } elseif {"$solvsize"=="big"} {
		    set solvent $molsolventdir/solventlibrary/tip3pbigdrop.c
	    } elseif {"$solvsize"=="small"} {
		    set solvent $molsolventdir/solventlibrary/tip3psmalldrop.c
	    } elseif {"$solvsize"=="tiny"} {
            set solvent $molsolventdir/solventlibrary/tip3ptinydrop.c
        }
    }
}

if {"$solvtype"=="tip4p"} {
    # Setting point charges for water model
	set OTcharge 0.0000008
	set HTcharge .52
	set ETcharge -1.04

    # Setting LJ parameters for oxygen
    # Values correspond to function for 4e(()^12-()^6)
    # These are fixed for Chemshell use later in ffgen.tcl
    # Units: sigma = A	epsilon = kcal/mol
    set solvsig [list {OT 3.154}]
    set solveps [list {OT .155}]

    # Sourcing the solvent system and definitions
    if {"$solvgeom"=="box"} {
        set solvent $molsolventdir/solventlibrary/tip3pbox.c
    } elseif {"$solvgeom"=="sphere"} {
        if {"$solvsize"=="9k"} {
            set solvent $molsolventdir/solventlibrary/tip4p9kdrop.c
    	} elseif {"$solvsize"=="big"} {
    		set solvent $molsolventdir/solventlibrary/tip4pbigdrop.c
    	} elseif {"$solvsize"=="small"} {
	    	set solvent $molsolventdir/solventlibrary/tip4psmalldrop.c
        } elseif {"$solvsize"=="tiny"} {
            set solvent $molsolventdir/solventlibrary/tip4ptinydrop.c
    	}
    }
}

# Setting a bunch of variables that will be used with insertsolute.tcl and ffgen.tcl
# Centre of solvent system
set centre [get_molecule_centre coords=$solvent]
puts "centre is $centre"
# Setting frozen radius and sbc radius (in bohr), leaves 3 angstrom frozen with 1 angstrom before SBC
set frozrad [expr [get_molecule_radius coords=$solvent] - 3*$angtobohr]
puts "frozrad is $frozrad"
set frozsbcrad [expr $frozrad + 1*$angtobohr]
puts "frozsbcrad is $frozsbcrad"
# Number of atoms in each solvent molecule (from connections to atom 1)
set solvmemb [lsort [get_molecule_members coords=$solvent atom_number=1]]
set solvatoms [llength $solvmemb]
puts "solvatoms is $solvatoms"
puts "solvmemb is $solvmemb"

foreach i $solvmemb {
    set atom [lindex [get_atom_entry coords=$solvent atom_number=$i] 0]
    append atom "T"
    lappend solvwitht $atom
}

puts "solvwitht is $solvwitht"
