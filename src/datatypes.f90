module datatypes

!This module defines datatypes and constant variables that should NEVER change between runs or different PepBD programs

implicit none
!groupdetails contains the following properties for a molecule: 
    !number of atoms
    !atom numes
    !residue name
    !atom coordinates
    !atom masses
!break residue into 3 groups based on amino acid structure. 
    !Group 1 is the N and alpha Carbon, plus their associated hydrogens
    !Group 2 is the amino acid sidechain
    !Group 3 is the C terminus of the residue
type groupdetails
    integer            :: cnum1, cnum2, cnum3 !number of atoms in groups 1 through 3
    character*4        :: atype1(20), atype2(60), atype3(20) !atom types in groups 1 through 3
    character*5        :: gtype !name of residue
    real               :: coo1(3,20), coo2(3,60), coo3(3,20) !coordinates of atoms in groups 1 through 3
    !real               :: bead_loc(3) !coordinate of coarse-grained bead at side-chain position
    real               :: masses1(20), masses2(60), masses3(20)
end type

type energyparameters
!energyparameters stores the following parameters for a residue:
    !atomnumber in system (helpful for checking which atoms are bonded to each other)
    !partial charges
    !Lennard Jones epsilon values
    !atomic radii
    !generalized born radii
    !fs - used to calculate radius in GB method
    !dieletric constants 
    !index into atom types for LCPO method
!as with "groupdetails" data type, the parameters are split into 3 groups
    integer            :: atomid1(20), atomid2(60), atomid3(20) !
    real               :: charge1(20), epsion1(20), r1(20), rborn1(20), fs1(20), dielecons1(20) !partial charge, LJ depth, LJ radius, GB radius, fs (not used), and dielectric constants for group 1 atoms
    real               :: charge2(60), epsion2(60), r2(60), rborn2(60), fs2(60), dielecons2(60)
    real               :: charge3(20), epsion3(20), r3(20), rborn3(20), fs3(20), dielecons3(20)
    integer            :: LCPO1(20), LCPO2(60), LCPO3(3)
end type

type databackup
    !this module is used to store the backup coordinates of alpha carbons needed to insert or remove proline properly
    real               :: coo(3) 
end type

type lib4aminoacid !contains rotamer information for all amino acids
    integer            :: cnum1, cnum2, cnum3 !number of atoms in groups 1 through 3
    integer            :: dihedralangle(6,35) !set of all possible dihedrals for sidechain
    character*4        :: atype1(20), atype2(60), atype3(20) !atom types for groups 1 through 3
    character*4        :: gtype !residue name
    real               :: coo1(3,20), coo2(3,60), coo3(3,20) !coordinates of atoms (defaults to values in rotamer 1)
    integer            :: num_sc_angles, rotanum !number of variable angles; number of rotamers
end type

type index4sidechain
    integer            :: class_No, member_No
end type

type conformer4sidechain
    real               :: member(3,15)
end type

type dihedralparameters
    integer            :: iph(36), jph(36), kph(36), lph(36), multiply(36)
    real               :: pk(4,36), pn(4,36), phase(4,36)
end type

type atomlink            ! The Data type "atomlink" is used to store the neighboring atoms for each atom.
    integer            :: linknum
    integer            :: linkindex(2,4)
end type

type LCPO_params
    character*4        :: atype
    integer            :: n
    real               :: r, p1, p2, p3, p4, p1_NLR, p2_NLR, p3_NLR, p4_NLR
endtype

!different "kind" parameters for defining variable precision

!> Single precision real numbers, 6 digits, range 10⁻³⁷ to 10³⁷-1; 32 bits
integer, parameter :: sp = selected_real_kind(6, 37)
!> Double precision real numbers, 15 digits, range 10⁻³⁰⁷ to 10³⁰⁷-1; 64 bits
integer, parameter :: dp = selected_real_kind(15, 307)
!> Quadruple precision real numbers, 33 digits, range 10⁻⁴⁹³¹ to 10⁴⁹³¹-1; 128 bits
integer, parameter :: qp = selected_real_kind(33, 4931)

!> Char length for integers, range -2⁷ to 2⁷-1; 8 bits
integer, parameter :: i1 = selected_int_kind(2)
!> Short length for integers, range -2¹⁵ to 2¹⁵-1; 16 bits
integer, parameter :: i2 = selected_int_kind(4)
!> Length of default integers, range -2³¹ to 2³¹-1; 32 bits
integer, parameter :: i4 = selected_int_kind(9)
!> Long length for integers, range -2⁶³ to 2⁶³-1; 64 bits
integer, parameter :: i8 = selected_int_kind(18)

!how to use the above "kind" parameters

end module datatypes