proc insertsolute {x} {

puts "Starting insert solute for $x"

# Defining needed variables from MD shell
global centre solvent solvsys$x solvmemb solvatoms angtobohr molsys$x soluteatomcentre$x solutetypes$x

set sys [set molsys$x]
puts "sys is $sys"
set mergesys [set solvsys$x]
puts "mergesys is $mergesys"
set solutetypes [set solutetypes$x]
puts "solutetypes here is $solutetypes"

push_banner_flag 0
set numatomsolute [get_number_of_atoms coords=$sys]
puts "initial numatomsolute for $x is $numatomsolute"

set molcentre [get_molecule_centre coords=$sys]
puts "initial molcentre is $molcentre"
puts "initial solv centre is $centre"

set abssmall 100
# Loop to orient mol-atoms relative to $centre
for { set i 1 } { $i <= $numatomsolute } { incr i 1 } {
    set atom [get_atom_entry coords=$sys atom_number=$i]
    set elem [lindex $atom 0]
    # Setting relative to molcentre
    set relx [expr [lindex $atom 1] - [lindex $molcentre 0] ]
    set rely [expr [lindex $atom 2] - [lindex $molcentre 1] ]
    set relz [expr [lindex $atom 3] - [lindex $molcentre 2] ]
    # This finds the atom that is closest to $molcentre
    set absdist [expr pow((pow($relx,2)+pow($rely,2)+pow($relz,2)),.5)]
    if { $absdist < $abssmall } {
        set abssmall $absdist
        set soluteatomcentre$x $i
    }

    # Sets the solute coordinates relative to solvent centre (so solute is in middle of solvent)
    set coordx [expr [lindex $centre 0] + $relx]
    set coordy [expr [lindex $centre 1] + $rely]
    set coordz [expr [lindex $centre 2] + $relz]
    replace_atom_entry coords=$sys atom_number=$i atom_entry= "$elem $coordx $coordy $coordz"
}

# Now we have mol fragment that has the correct center coordinates
set molcentre [get_molecule_centre coords=$sys]
puts "edited molcentre is $molcentre"
puts "soluteatomcentre$x is [set soluteatomcentre$x]"

# Setting atom numbers of the solute (should be 1 ... # solute atoms)
set solutelist [lsort -integer [get_molecule_members coords=$sys atom_number=1]]
puts "solutelist is $solutelist"

# Finding radius for solute SBC, based on radius of solute molecule (in bohr)
set qmsbcrad [expr [get_molecule_radius coords=$sys] + 3*$angtobohr]

# Fragments are merged and new connecitivity is automatically created
# One can modify the scale/toler connectivity before so that more or fewer atoms become connected
# This would increase/decrease the clashing solvent molecules to be deleted
# R_ij < scale*(rad_i + rad_j) + toler.
set chemsh_default_connectivity_scale 1.0
set chemsh_default_connectivity_toler 1.5

# Merging fragments (solute first in .c)
c_merge_fragments base=$sys add=$solvent result=$mergesys
# Saving a copy of the initial merged geomertry
write_xyz coords=$mergesys file=system${x}merge.xyz
# Adding periodic information back in if it exists in solvent file
if {[string bytelength [get_cell coords=$solvent]]} {
    puts "cell exists with [get_cell coords=$solvent]"
    set_cell coords=$mergesys cell={[get_cell coords=$solvent]}
}

# Defining number of atoms in initial merged system
set numatomsys [get_number_of_atoms coords=$mergesys]
puts "Numatoms after merge is is $numatomsys"

# Delete clashing water molecules within 2 angstroms of solute atoms
# Here we have to be careful as the numbering changes with each deletion
for {set i 1} {$i <= $numatomsolute} {incr i} {
    set soluteentry [get_atom_entry coords=$mergesys atom_number=$i]
    set solutecoords [lreplace $soluteentry 0 0]
    for {set j [expr $numatomsolute+1]} {$j <= $numatomsys} {incr j} {
        set jentry [get_atom_entry coords=$mergesys atom_number=$j]
        set jcoords [lreplace $jentry 0 0]
        set dist [distance $solutecoords $jcoords]
        if {$dist <= [expr 2*$angtobohr]} {
            set solvmollist [get_molecule_members coords=$mergesys atom_number=$j]
            foreach k $solvmollist {
                lappend clashes $k
            }
        }
    }
    puts "clashes is $clashes"
    set molclash [listuniq $clashes]
    puts "molclash is $molclash"
}
set clash [listcomp $molclash $solutelist]
set clash_sorted [lsort -integer -decreasing $clash]

puts "clash_sorted is $clash_sorted"
set delmolnum [expr [llength $clash_sorted] / $solvatoms]
puts "delmolnum is $delmolnum"
for { set i 0 } { $i < [llength $clash_sorted] } { incr i 1 } {
set b [lindex $clash_sorted $i]
set c [expr $b - 1]
puts "Deleting atom $b"
delete_atom_entry coords=$mergesys atom_number=$b
}

# Saving copy of geometry after cleanup
write_xyz coords=$mergesys file=sys${x}afterdel.xyz

# Defining number of atoms in final merged system
set numatomafterdel [get_number_of_atoms coords=$mergesys]
puts "numatomafterdel is $numatomafterdel"
# Finding centre of final merged system
set centre [get_molecule_centre coords=$mergesys]
puts "final centre is $centre"
# Finding radius of merged system
set radius [get_molecule_radius coords=$mergesys]
puts "merged radius for $x is $radius"

# Calculating froz/act lists, types, constraints, groups
set lists [set_system_lists $solutelist $solutetypes $mergesys $x $solvmemb $qmsbcrad]
set froz [lindex $lists 0]
set act [lindex $lists 1]
set restraints [lindex $lists 2]
set types [lindex $lists 3]
set con [lindex $lists 4]
set groups [lindex $lists 5]

# Setting QM region hydrogen masses to deuterium
set masses [prepmass $mergesys $solutelist]

#######################
#pop_banner_flag

# Exporting variables to be used later in input file
variable numatomsolute$x $numatomsolute
puts "numatomsolute$x is $numatomsolute"
variable solutelist$x $solutelist
puts "solutelist$x is $solutelist"
variable solutetypes$x $solutetypes
puts "solutetypes$x is $solutetypes"
variable froz$x $froz
puts "froz$x is $froz"
variable restraints$x $restraints
puts "restraints$x is $restraints"
variable types$x $types
puts "types$x is $types"
variable con$x $con
puts "con$x is $con"
variable groups$x $groups
puts "groups$x is $groups"
variable masses$x $masses
puts "masses$x is $masses"
variable radius$x $radius
puts "radius$x is $radius"
}
