# Opening forcefield.ff file
set ff [open forcefield.ff "w"]

# Writing charge information
if {"$mmodel"=="spc"} {
    set solvtypes {OT HT HT}
    puts $ff "charge OT -.82"
    puts $ff "charge HT .41"
} elseif {"$mmodel"=="tip3p"} {
    set solvtypes {OT HT HT}
    puts $ff "charge OT -.834"
    puts $ff "charge HT .417"
} elseif {"$mmodel"=="tip4p"} {
    set solvtypes {OT HT HT ET}
    puts $ff "charge OT 0"
    puts $ff "charge HT .52"
    puts $ff "charge ET -1.04"
    set rOMa .15
} elseif {"$mmodel"=="blank"} {
    set solvtypes {OT HT HT}
    puts $ff "charge OT 0"
    puts $ff "charge HT 0"
    puts $ff "charge ET 0"
}

if {"$pmodel"=="pspc"} {
    set solvtypesp {OTp HTp HTp}
    puts $ff "charge OTp -.669"
    puts $ff "charge HTp .3345"
    set drudeq 2.08241
    set drudepol 9.717560
} elseif {"$pmodel"=="swm4dp"} {
    set solvtypesp {OTp HTp HTp ETp}
    puts $ff "charge OTp -.0000008"
    puts $ff "charge HTp .5537"
    puts $ff "charge ETp -1.1074"
    set rOMap .23808
    set drudeq 1.77185
    set drudepol 7.035272
} elseif {"$pmodel"=="swm4ndp"} {
    set solvtypesp {OTp HTp HTp ETp}
    puts $ff "charge OTp .0000008"
    puts $ff "charge HTp .55733"
    puts $ff "charge ETp -1.11466"
    set rOMap .24034
    set drudeq -1.71636
    set drudepol 6.601557
} elseif {"$pmodel"=="swm6"} {
    set solvtypesp {OTp HTp HTp ETp LTp LTp}
    puts $ff "charge OTp .288"
    puts $ff "charge HTp .5307"
    puts $ff "charge ETp -1.1334"
    puts $ff "charge LTp -.108"
    set rOMap .247
    set rOLap .315
    set aLLdp 101.098
    set drudeq -1.62789
    set drudepol 6.680850
}

# Closing forcefield.ff file to write data
close $ff
