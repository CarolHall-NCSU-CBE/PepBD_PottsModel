module database

use datatypes
use sys_vars
use utilities
use math

contains

subroutine PH_checking(group_pep, group_rec)
implicit none
integer                     :: i, ic, count
type(groupdetails)          :: group_pep(pep_res), group_rec(rec_res)

count=0
if(ph_value.le.3.9) then
    if(receptor_name=="rna") then
        count=pep_res
    elseif(receptor_name=="peptide") then
        count=gnum
    elseif(receptor_name=="molecule") then
        count=pep_res
    endif
    do i=1, count
        ic=i
        if(group_pep(ic)%gtype=="GLY".or.group_pep(ic)%gtype=="NGLY".or.group_pep(ic)%gtype=="CGLY".or.group_pep(ic)%gtype=="LEU".or.group_pep(ic)%gtype=="NLEU".or.group_pep(ic)%gtype=="CLEU" &
            .or.group_pep(ic)%gtype=="VAL".or.group_pep(ic)%gtype=="NVAL".or.group_pep(ic)%gtype=="CVAL".or.group_pep(ic)%gtype=="ILE".or.group_pep(ic)%gtype=="NILE".or.group_pep(ic)%gtype=="CILE" &
            .or.group_pep(ic)%gtype=="MET".or.group_pep(ic)%gtype=="NMET".or.group_pep(ic)%gtype=="CMET".or.group_pep(ic)%gtype=="PHE".or.group_pep(ic)%gtype=="NPHE".or.group_pep(ic)%gtype=="CPHE" &
            .or.group_pep(ic)%gtype=="TYR".or.group_pep(ic)%gtype=="NTYR".or.group_pep(ic)%gtype=="CTYR".or.group_pep(ic)%gtype=="TRP".or.group_pep(ic)%gtype=="NTRP".or.group_pep(ic)%gtype=="CTRP" &
            .or.group_pep(ic)%gtype=="ARG".or.group_pep(ic)%gtype=="NARG".or.group_pep(ic)%gtype=="CARG".or.group_pep(ic)%gtype=="LYS".or.group_pep(ic)%gtype=="NLYS".or.group_pep(ic)%gtype=="CLYS" &
            .or.group_pep(ic)%gtype=="SER".or.group_pep(ic)%gtype=="NSER".or.group_pep(ic)%gtype=="CSER".or.group_pep(ic)%gtype=="THR".or.group_pep(ic)%gtype=="NTHR".or.group_pep(ic)%gtype=="CTHR" &
            .or.group_pep(ic)%gtype=="ASN".or.group_pep(ic)%gtype=="NASN".or.group_pep(ic)%gtype=="CASN".or.group_pep(ic)%gtype=="GLN".or.group_pep(ic)%gtype=="NGLN".or.group_pep(ic)%gtype=="CGLN" &
            .or.group_pep(ic)%gtype=="HIP".or.group_pep(ic)%gtype=="NHIP".or.group_pep(ic)%gtype=="CHIP".or.group_pep(ic)%gtype=="PRO".or.group_pep(ic)%gtype=="NPRO".or.group_pep(ic)%gtype=="CPRO" &
            .or.group_pep(ic)%gtype=="CYS".or.group_pep(ic)%gtype=="NCYS".or.group_pep(ic)%gtype=="CCYS".or.group_pep(ic)%gtype=="ALA".or.group_pep(ic)%gtype=="NALA".or.group_pep(ic)%gtype=="CALA" &
            .or.group_pep(ic)%gtype=="GLH".or.group_pep(ic)%gtype=="NGLH".or.group_pep(ic)%gtype=="CGLH".or.group_pep(ic)%gtype=="ASH".or.group_pep(ic)%gtype=="NASH".or.group_pep(ic)%gtype=="CASH") then
        else
            write(*,*) group_pep(ic)%gtype, "of the", ic, "residue could not be found in the PH of", ph_value
            write(*,*) "Please check whether the group_pep name is right or not!"
            stop
        endif
    enddo

elseif(ph_value.gt.3.9.and.ph_value.le.4.3) then
    if(receptor_name=="rna") then
        count=pep_res
    elseif(receptor_name=="peptide") then
        count=gnum
    elseif(receptor_name=="molecule") then
        count=pep_res
    endif
    do i=1, count
        ic=i
        if(group_pep(ic)%gtype=="GLY".or.group_pep(ic)%gtype=="NGLY".or.group_pep(ic)%gtype=="CGLY".or.group_pep(ic)%gtype=="LEU".or.group_pep(ic)%gtype=="NLEU".or.group_pep(ic)%gtype=="CLEU" &
            .or.group_pep(ic)%gtype=="VAL".or.group_pep(ic)%gtype=="NVAL".or.group_pep(ic)%gtype=="CVAL".or.group_pep(ic)%gtype=="ILE".or.group_pep(ic)%gtype=="NILE".or.group_pep(ic)%gtype=="CILE" &
            .or.group_pep(ic)%gtype=="MET".or.group_pep(ic)%gtype=="NMET".or.group_pep(ic)%gtype=="CMET".or.group_pep(ic)%gtype=="PHE".or.group_pep(ic)%gtype=="NPHE".or.group_pep(ic)%gtype=="CPHE" &
            .or.group_pep(ic)%gtype=="TYR".or.group_pep(ic)%gtype=="NTYR".or.group_pep(ic)%gtype=="CTYR".or.group_pep(ic)%gtype=="TRP".or.group_pep(ic)%gtype=="NTRP".or.group_pep(ic)%gtype=="CTRP" &
            .or.group_pep(ic)%gtype=="ARG".or.group_pep(ic)%gtype=="NARG".or.group_pep(ic)%gtype=="CARG".or.group_pep(ic)%gtype=="LYS".or.group_pep(ic)%gtype=="NLYS".or.group_pep(ic)%gtype=="CLYS" &
            .or.group_pep(ic)%gtype=="SER".or.group_pep(ic)%gtype=="NSER".or.group_pep(ic)%gtype=="CSER".or.group_pep(ic)%gtype=="THR".or.group_pep(ic)%gtype=="NTHR".or.group_pep(ic)%gtype=="CTHR" &
            .or.group_pep(ic)%gtype=="ASN".or.group_pep(ic)%gtype=="NASN".or.group_pep(ic)%gtype=="CASN".or.group_pep(ic)%gtype=="GLN".or.group_pep(ic)%gtype=="NGLN".or.group_pep(ic)%gtype=="CGLN" &
            .or.group_pep(ic)%gtype=="HIP".or.group_pep(ic)%gtype=="NHIP".or.group_pep(ic)%gtype=="CHIP".or.group_pep(ic)%gtype=="PRO".or.group_pep(ic)%gtype=="NPRO".or.group_pep(ic)%gtype=="CPRO" &
            .or.group_pep(ic)%gtype=="CYS".or.group_pep(ic)%gtype=="NCYS".or.group_pep(ic)%gtype=="CCYS".or.group_pep(ic)%gtype=="ALA".or.group_pep(ic)%gtype=="NALA".or.group_pep(ic)%gtype=="CALA" &
            .or.group_pep(ic)%gtype=="GLH".or.group_pep(ic)%gtype=="NGLH".or.group_pep(ic)%gtype=="CGLH".or.group_pep(ic)%gtype=="ASP".or.group_pep(ic)%gtype=="NASP".or.group_pep(ic)%gtype=="CASP") then
        else
            write(*,*) group_pep(ic)%gtype, "of the", ic, "residue could not be found in the PH of", ph_value
            write(*,*) "Please check whether the group_pep name is right or not!"
            stop
        endif
    enddo

elseif(ph_value.gt.4.3.and.ph_value.le.6.0) then
    if(receptor_name=="rna") then
        count=pep_res
    elseif(receptor_name=="peptide") then
        count=gnum
    elseif(receptor_name=="molecule") then
        count=pep_res
    endif
    do i=1, count
        ic=i
        if(group_pep(ic)%gtype=="GLY".or.group_pep(ic)%gtype=="NGLY".or.group_pep(ic)%gtype=="CGLY".or.group_pep(ic)%gtype=="LEU".or.group_pep(ic)%gtype=="NLEU".or.group_pep(ic)%gtype=="CLEU" &
            .or.group_pep(ic)%gtype=="VAL".or.group_pep(ic)%gtype=="NVAL".or.group_pep(ic)%gtype=="CVAL".or.group_pep(ic)%gtype=="ILE".or.group_pep(ic)%gtype=="NILE".or.group_pep(ic)%gtype=="CILE" &
            .or.group_pep(ic)%gtype=="MET".or.group_pep(ic)%gtype=="NMET".or.group_pep(ic)%gtype=="CMET".or.group_pep(ic)%gtype=="PHE".or.group_pep(ic)%gtype=="NPHE".or.group_pep(ic)%gtype=="CPHE" &
            .or.group_pep(ic)%gtype=="TYR".or.group_pep(ic)%gtype=="NTYR".or.group_pep(ic)%gtype=="CTYR".or.group_pep(ic)%gtype=="TRP".or.group_pep(ic)%gtype=="NTRP".or.group_pep(ic)%gtype=="CTRP" &
            .or.group_pep(ic)%gtype=="ARG".or.group_pep(ic)%gtype=="NARG".or.group_pep(ic)%gtype=="CARG".or.group_pep(ic)%gtype=="LYS".or.group_pep(ic)%gtype=="NLYS".or.group_pep(ic)%gtype=="CLYS" &
            .or.group_pep(ic)%gtype=="SER".or.group_pep(ic)%gtype=="NSER".or.group_pep(ic)%gtype=="CSER".or.group_pep(ic)%gtype=="THR".or.group_pep(ic)%gtype=="NTHR".or.group_pep(ic)%gtype=="CTHR" &
            .or.group_pep(ic)%gtype=="ASN".or.group_pep(ic)%gtype=="NASN".or.group_pep(ic)%gtype=="CASN".or.group_pep(ic)%gtype=="GLN".or.group_pep(ic)%gtype=="NGLN".or.group_pep(ic)%gtype=="CGLN" &
            .or.group_pep(ic)%gtype=="HIP".or.group_pep(ic)%gtype=="NHIP".or.group_pep(ic)%gtype=="CHIP".or.group_pep(ic)%gtype=="PRO".or.group_pep(ic)%gtype=="NPRO".or.group_pep(ic)%gtype=="CPRO" &
            .or.group_pep(ic)%gtype=="CYS".or.group_pep(ic)%gtype=="NCYS".or.group_pep(ic)%gtype=="CCYS".or.group_pep(ic)%gtype=="ALA".or.group_pep(ic)%gtype=="NALA".or.group_pep(ic)%gtype=="CALA" &
            .or.group_pep(ic)%gtype=="GLU".or.group_pep(ic)%gtype=="NGLU".or.group_pep(ic)%gtype=="CGLU".or.group_pep(ic)%gtype=="ASP".or.group_pep(ic)%gtype=="NASP".or.group_pep(ic)%gtype=="CASP") then
        else
            !open(20, file="error.txt", access="append")
                write(*,*) group_pep(ic)%gtype, "of the", ic, "residue could not be found in the PH of", ph_value
                write(*,*) "Please check whether the group_pep name is right or not!"
            !close(20)
            stop
        endif
    enddo

elseif(ph_value.gt.6.0.and.ph_value.lt.8.3) then
    do i=1, pep_res
        if(group_pep(i)%gtype=="GLY".or.group_pep(i)%gtype=="NGLY".or.group_pep(i)%gtype=="CGLY".or.group_pep(i)%gtype=="LEU".or.group_pep(i)%gtype=="NLEU".or.group_pep(i)%gtype=="CLEU" &
            .or.group_pep(i)%gtype=="VAL".or.group_pep(i)%gtype=="NVAL".or.group_pep(i)%gtype=="CVAL".or.group_pep(i)%gtype=="ILE".or.group_pep(i)%gtype=="NILE".or.group_pep(i)%gtype=="CILE" &
            .or.group_pep(i)%gtype=="MET".or.group_pep(i)%gtype=="NMET".or.group_pep(i)%gtype=="CMET".or.group_pep(i)%gtype=="PHE".or.group_pep(i)%gtype=="NPHE".or.group_pep(i)%gtype=="CPHE" &
            .or.group_pep(i)%gtype=="TYR".or.group_pep(i)%gtype=="NTYR".or.group_pep(i)%gtype=="CTYR".or.group_pep(i)%gtype=="TRP".or.group_pep(i)%gtype=="NTRP".or.group_pep(i)%gtype=="CTRP" &
            .or.group_pep(i)%gtype=="ARG".or.group_pep(i)%gtype=="NARG".or.group_pep(i)%gtype=="CARG".or.group_pep(i)%gtype=="LYS".or.group_pep(i)%gtype=="NLYS".or.group_pep(i)%gtype=="CLYS" &
            .or.group_pep(i)%gtype=="SER".or.group_pep(i)%gtype=="NSER".or.group_pep(i)%gtype=="CSER".or.group_pep(i)%gtype=="THR".or.group_pep(i)%gtype=="NTHR".or.group_pep(i)%gtype=="CTHR" &
            .or.group_pep(i)%gtype=="ASN".or.group_pep(i)%gtype=="NASN".or.group_pep(i)%gtype=="CASN".or.group_pep(i)%gtype=="GLN".or.group_pep(i)%gtype=="NGLN".or.group_pep(i)%gtype=="CGLN" &
            .or.group_pep(i)%gtype=="HIE".or.group_pep(i)%gtype=="NHIE".or.group_pep(i)%gtype=="CHIE".or.group_pep(i)%gtype=="PRO".or.group_pep(i)%gtype=="NPRO".or.group_pep(i)%gtype=="CPRO" &
            .or.group_pep(i)%gtype=="CYS".or.group_pep(i)%gtype=="NCYS".or.group_pep(i)%gtype=="CCYS".or.group_pep(i)%gtype=="ALA".or.group_pep(i)%gtype=="NALA".or.group_pep(i)%gtype=="CALA" &
            .or.group_pep(i)%gtype=="GLU".or.group_pep(i)%gtype=="NGLU".or.group_pep(i)%gtype=="CGLU".or.group_pep(i)%gtype=="ASP".or.group_pep(i)%gtype=="NASP".or.group_pep(i)%gtype=="CASP") then
        else
            !open(20, file="error.txt", access="append")
            write(*, "(A,A,A,I0,A,F6.2)") "Residue type ", group_pep(i)%gtype, " at receptor position ", i, " cannot be found at a pH of ", ph_value
            write(*,*) "Please check the pdb file and input file!"
            !close(20)
            stop
        endif
    enddo

    !check the receptor residues, if it is a protein
    if(receptor_name=="peptide") then
        do i=1, rec_res
            if(group_rec(i)%gtype=="GLY".or.group_rec(i)%gtype=="NGLY".or.group_rec(i)%gtype=="CGLY".or.group_rec(i)%gtype=="LEU".or.group_rec(i)%gtype=="NLEU".or.group_rec(i)%gtype=="CLEU" &
            .or.group_rec(i)%gtype=="VAL".or.group_rec(i)%gtype=="NVAL".or.group_rec(i)%gtype=="CVAL".or.group_rec(i)%gtype=="ILE".or.group_rec(i)%gtype=="NILE".or.group_rec(i)%gtype=="CILE" &
            .or.group_rec(i)%gtype=="MET".or.group_rec(i)%gtype=="NMET".or.group_rec(i)%gtype=="CMET".or.group_rec(i)%gtype=="PHE".or.group_rec(i)%gtype=="NPHE".or.group_rec(i)%gtype=="CPHE" &
            .or.group_rec(i)%gtype=="TYR".or.group_rec(i)%gtype=="NTYR".or.group_rec(i)%gtype=="CTYR".or.group_rec(i)%gtype=="TRP".or.group_rec(i)%gtype=="NTRP".or.group_rec(i)%gtype=="CTRP" &
            .or.group_rec(i)%gtype=="ARG".or.group_rec(i)%gtype=="NARG".or.group_rec(i)%gtype=="CARG".or.group_rec(i)%gtype=="LYS".or.group_rec(i)%gtype=="NLYS".or.group_rec(i)%gtype=="CLYS" &
            .or.group_rec(i)%gtype=="SER".or.group_rec(i)%gtype=="NSER".or.group_rec(i)%gtype=="CSER".or.group_rec(i)%gtype=="THR".or.group_rec(i)%gtype=="NTHR".or.group_rec(i)%gtype=="CTHR" &
            .or.group_rec(i)%gtype=="ASN".or.group_rec(i)%gtype=="NASN".or.group_rec(i)%gtype=="CASN".or.group_rec(i)%gtype=="GLN".or.group_rec(i)%gtype=="NGLN".or.group_rec(i)%gtype=="CGLN" &
            .or.group_rec(i)%gtype=="HIE".or.group_rec(i)%gtype=="NHIE".or.group_rec(i)%gtype=="CHIE".or.group_rec(i)%gtype=="PRO".or.group_rec(i)%gtype=="NPRO".or.group_rec(i)%gtype=="CPRO" &
            .or.group_rec(i)%gtype=="CYS".or.group_rec(i)%gtype=="NCYS".or.group_rec(i)%gtype=="CCYS".or.group_rec(i)%gtype=="ALA".or.group_rec(i)%gtype=="NALA".or.group_rec(i)%gtype=="CALA" &
            .or.group_rec(i)%gtype=="GLU".or.group_rec(i)%gtype=="NGLU".or.group_rec(i)%gtype=="CGLU".or.group_rec(i)%gtype=="ASP".or.group_rec(i)%gtype=="NASP".or.group_rec(i)%gtype=="CASP") then
                continue
            else
                write(*, "(A,A,A,I0,A,F6.2)") "Residue type ", group_rec(i)%gtype, " at receptor position ", i, " cannot be found at a pH of ", ph_value
                write(*,*) "Please check the pdb file and input file!"
                stop
            endif
        enddo
    endif

elseif(ph_value.ge.8.3.and.ph_value.lt.10.1) then
    if(receptor_name=="rna") then
        count=pep_res
    elseif(receptor_name=="peptide") then
        count=gnum
    elseif(receptor_name=="molecule") then
        count=pep_res
    endif
    do i=1, count
        ic=i
        if(group_pep(ic)%gtype=="GLY".or.group_pep(ic)%gtype=="NGLY".or.group_pep(ic)%gtype=="CGLY".or.group_pep(ic)%gtype=="LEU".or.group_pep(ic)%gtype=="NLEU".or.group_pep(ic)%gtype=="CLEU" &
            .or.group_pep(ic)%gtype=="VAL".or.group_pep(ic)%gtype=="NVAL".or.group_pep(ic)%gtype=="CVAL".or.group_pep(ic)%gtype=="ILE".or.group_pep(ic)%gtype=="NILE".or.group_pep(ic)%gtype=="CILE" &
            .or.group_pep(ic)%gtype=="MET".or.group_pep(ic)%gtype=="NMET".or.group_pep(ic)%gtype=="CMET".or.group_pep(ic)%gtype=="PHE".or.group_pep(ic)%gtype=="NPHE".or.group_pep(ic)%gtype=="CPHE" &
            .or.group_pep(ic)%gtype=="TYR".or.group_pep(ic)%gtype=="NTYR".or.group_pep(ic)%gtype=="CTYR".or.group_pep(ic)%gtype=="TRP".or.group_pep(ic)%gtype=="NTRP".or.group_pep(ic)%gtype=="CTRP" &
            .or.group_pep(ic)%gtype=="ARG".or.group_pep(ic)%gtype=="NARG".or.group_pep(ic)%gtype=="CARG".or.group_pep(ic)%gtype=="LYS".or.group_pep(ic)%gtype=="NLYS".or.group_pep(ic)%gtype=="CLYS" &
            .or.group_pep(ic)%gtype=="SER".or.group_pep(ic)%gtype=="NSER".or.group_pep(ic)%gtype=="CSER".or.group_pep(ic)%gtype=="THR".or.group_pep(ic)%gtype=="NTHR".or.group_pep(ic)%gtype=="CTHR" &
            .or.group_pep(ic)%gtype=="ASN".or.group_pep(ic)%gtype=="NASN".or.group_pep(ic)%gtype=="CASN".or.group_pep(ic)%gtype=="GLN".or.group_pep(ic)%gtype=="NGLN".or.group_pep(ic)%gtype=="CGLN" &
            .or.group_pep(ic)%gtype=="HIE".or.group_pep(ic)%gtype=="NHIE".or.group_pep(ic)%gtype=="CHIE".or.group_pep(ic)%gtype=="PRO".or.group_pep(ic)%gtype=="NPRO".or.group_pep(ic)%gtype=="CPRO" &
            .or.group_pep(ic)%gtype=="CYT".or.group_pep(ic)%gtype=="NCYT".or.group_pep(ic)%gtype=="CCYT".or.group_pep(ic)%gtype=="ALA".or.group_pep(ic)%gtype=="NALA".or.group_pep(ic)%gtype=="CALA" &
            .or.group_pep(ic)%gtype=="GLU".or.group_pep(ic)%gtype=="NGLU".or.group_pep(ic)%gtype=="CGLU".or.group_pep(ic)%gtype=="ASP".or.group_pep(ic)%gtype=="NASP".or.group_pep(ic)%gtype=="CASP") then
        else
            !open(20, file="error.txt", access="append")
                write(*,*) group_pep(ic)%gtype, "of the", ic, "residue could not be found in the PH of", ph_value
                write(*,*) "Please check whether the group_pep name is right or not!"
            !close(20)
            stop
        endif
    enddo

elseif(ph_value.ge.10.1.and.ph_value.lt.10.5) then
    if(receptor_name=="rna") then
        count=pep_res
    elseif(receptor_name=="peptide") then
        count=gnum
    elseif(receptor_name=="molecule") then
        count=pep_res
    endif
    do i=1, count
        ic=i
        if(group_pep(ic)%gtype=="GLY".or.group_pep(ic)%gtype=="NGLY".or.group_pep(ic)%gtype=="CGLY".or.group_pep(ic)%gtype=="LEU".or.group_pep(ic)%gtype=="NLEU".or.group_pep(ic)%gtype=="CLEU" &
            .or.group_pep(ic)%gtype=="VAL".or.group_pep(ic)%gtype=="NVAL".or.group_pep(ic)%gtype=="CVAL".or.group_pep(ic)%gtype=="ILE".or.group_pep(ic)%gtype=="NILE".or.group_pep(ic)%gtype=="CILE" &
            .or.group_pep(ic)%gtype=="MET".or.group_pep(ic)%gtype=="NMET".or.group_pep(ic)%gtype=="CMET".or.group_pep(ic)%gtype=="PHE".or.group_pep(ic)%gtype=="NPHE".or.group_pep(ic)%gtype=="CPHE" &
            .or.group_pep(ic)%gtype=="TYX".or.group_pep(ic)%gtype=="NTYX".or.group_pep(ic)%gtype=="CTYX".or.group_pep(ic)%gtype=="TRP".or.group_pep(ic)%gtype=="NTRP".or.group_pep(ic)%gtype=="CTRP" &
            .or.group_pep(ic)%gtype=="ARG".or.group_pep(ic)%gtype=="NARG".or.group_pep(ic)%gtype=="CARG".or.group_pep(ic)%gtype=="LYS".or.group_pep(ic)%gtype=="NLYS".or.group_pep(ic)%gtype=="CLYS" &
            .or.group_pep(ic)%gtype=="SER".or.group_pep(ic)%gtype=="NSER".or.group_pep(ic)%gtype=="CSER".or.group_pep(ic)%gtype=="THR".or.group_pep(ic)%gtype=="NTHR".or.group_pep(ic)%gtype=="CTHR" &
            .or.group_pep(ic)%gtype=="ASN".or.group_pep(ic)%gtype=="NASN".or.group_pep(ic)%gtype=="CASN".or.group_pep(ic)%gtype=="GLN".or.group_pep(ic)%gtype=="NGLN".or.group_pep(ic)%gtype=="CGLN" &
            .or.group_pep(ic)%gtype=="HIE".or.group_pep(ic)%gtype=="NHIE".or.group_pep(ic)%gtype=="CHIE".or.group_pep(ic)%gtype=="PRO".or.group_pep(ic)%gtype=="NPRO".or.group_pep(ic)%gtype=="CPRO" &
            .or.group_pep(ic)%gtype=="CYT".or.group_pep(ic)%gtype=="NCYT".or.group_pep(ic)%gtype=="CCYT".or.group_pep(ic)%gtype=="ALA".or.group_pep(ic)%gtype=="NALA".or.group_pep(ic)%gtype=="CALA" &
            .or.group_pep(ic)%gtype=="GLU".or.group_pep(ic)%gtype=="NGLU".or.group_pep(ic)%gtype=="CGLU".or.group_pep(ic)%gtype=="ASP".or.group_pep(ic)%gtype=="NASP".or.group_pep(ic)%gtype=="CASP") then
        else
            !open(20, file="error.txt", access="append")
                write(*,*) group_pep(ic)%gtype, "of the", ic, "residue could not be found in the PH of", ph_value
                write(*,*) "Please check whether the group_pep name is right or not!"
            !close(20)
            stop
        endif
    enddo

elseif(ph_value.ge.10.5.and.ph_value.lt.12.5) then
    if(receptor_name=="rna") then
        count=pep_res
    elseif(receptor_name=="peptide") then
        count=gnum
    elseif(receptor_name=="molecule") then
        count=pep_res
    endif
    do i=1, count
        ic=i
        if(group_pep(ic)%gtype=="GLY".or.group_pep(ic)%gtype=="NGLY".or.group_pep(ic)%gtype=="CGLY".or.group_pep(ic)%gtype=="LEU".or.group_pep(ic)%gtype=="NLEU".or.group_pep(ic)%gtype=="CLEU" &
            .or.group_pep(ic)%gtype=="VAL".or.group_pep(ic)%gtype=="NVAL".or.group_pep(ic)%gtype=="CVAL".or.group_pep(ic)%gtype=="ILE".or.group_pep(ic)%gtype=="NILE".or.group_pep(ic)%gtype=="CILE" &
            .or.group_pep(ic)%gtype=="MET".or.group_pep(ic)%gtype=="NMET".or.group_pep(ic)%gtype=="CMET".or.group_pep(ic)%gtype=="PHE".or.group_pep(ic)%gtype=="NPHE".or.group_pep(ic)%gtype=="CPHE" &
            .or.group_pep(ic)%gtype=="TYX".or.group_pep(ic)%gtype=="NTYX".or.group_pep(ic)%gtype=="CTYX".or.group_pep(ic)%gtype=="TRP".or.group_pep(ic)%gtype=="NTRP".or.group_pep(ic)%gtype=="CTRP" &
            .or.group_pep(ic)%gtype=="ARG".or.group_pep(ic)%gtype=="NARG".or.group_pep(ic)%gtype=="CARG".or.group_pep(ic)%gtype=="LYN".or.group_pep(ic)%gtype=="NLYN".or.group_pep(ic)%gtype=="CLYN" &
            .or.group_pep(ic)%gtype=="SER".or.group_pep(ic)%gtype=="NSER".or.group_pep(ic)%gtype=="CSER".or.group_pep(ic)%gtype=="THR".or.group_pep(ic)%gtype=="NTHR".or.group_pep(ic)%gtype=="CTHR" &
            .or.group_pep(ic)%gtype=="ASN".or.group_pep(ic)%gtype=="NASN".or.group_pep(ic)%gtype=="CASN".or.group_pep(ic)%gtype=="GLN".or.group_pep(ic)%gtype=="NGLN".or.group_pep(ic)%gtype=="CGLN" &
            .or.group_pep(ic)%gtype=="HIE".or.group_pep(ic)%gtype=="NHIE".or.group_pep(ic)%gtype=="CHIE".or.group_pep(ic)%gtype=="PRO".or.group_pep(ic)%gtype=="NPRO".or.group_pep(ic)%gtype=="CPRO" &
            .or.group_pep(ic)%gtype=="CYT".or.group_pep(ic)%gtype=="NCYT".or.group_pep(ic)%gtype=="CCYT".or.group_pep(ic)%gtype=="ALA".or.group_pep(ic)%gtype=="NALA".or.group_pep(ic)%gtype=="CALA" &
            .or.group_pep(ic)%gtype=="GLU".or.group_pep(ic)%gtype=="NGLU".or.group_pep(ic)%gtype=="CGLU".or.group_pep(ic)%gtype=="ASP".or.group_pep(ic)%gtype=="NASP".or.group_pep(ic)%gtype=="CASP") then
        else
            !open(20, file="error.txt", access="append")
                write(*,*) group_pep(ic)%gtype, "of the", ic, "residue could not be found in the PH of", ph_value
                write(*,*) "Please check whether the group_pep name is right or not!"
            !close(20)
            stop
        endif
    enddo

elseif(ph_value.ge.12.5) then
    if(receptor_name=="rna") then
        count=pep_res
    elseif(receptor_name=="peptide") then
        count=gnum
    elseif(receptor_name=="molecule") then
        count=pep_res
    endif
    do i=1, count
        ic=i
        if(group_pep(ic)%gtype=="GLY".or.group_pep(ic)%gtype=="NGLY".or.group_pep(ic)%gtype=="CGLY".or.group_pep(ic)%gtype=="LEU".or.group_pep(ic)%gtype=="NLEU".or.group_pep(ic)%gtype=="CLEU" &
            .or.group_pep(ic)%gtype=="VAL".or.group_pep(ic)%gtype=="NVAL".or.group_pep(ic)%gtype=="CVAL".or.group_pep(ic)%gtype=="ILE".or.group_pep(ic)%gtype=="NILE".or.group_pep(ic)%gtype=="CILE" &
            .or.group_pep(ic)%gtype=="MET".or.group_pep(ic)%gtype=="NMET".or.group_pep(ic)%gtype=="CMET".or.group_pep(ic)%gtype=="PHE".or.group_pep(ic)%gtype=="NPHE".or.group_pep(ic)%gtype=="CPHE" &
            .or.group_pep(ic)%gtype=="TYX".or.group_pep(ic)%gtype=="NTYX".or.group_pep(ic)%gtype=="CTYX".or.group_pep(ic)%gtype=="TRP".or.group_pep(ic)%gtype=="NTRP".or.group_pep(ic)%gtype=="CTRP" &
            .or.group_pep(ic)%gtype=="ARN".or.group_pep(ic)%gtype=="NARN".or.group_pep(ic)%gtype=="CARN".or.group_pep(ic)%gtype=="LYN".or.group_pep(ic)%gtype=="NLYN".or.group_pep(ic)%gtype=="CLYN" &
            .or.group_pep(ic)%gtype=="SER".or.group_pep(ic)%gtype=="NSER".or.group_pep(ic)%gtype=="CSER".or.group_pep(ic)%gtype=="THR".or.group_pep(ic)%gtype=="NTHR".or.group_pep(ic)%gtype=="CTHR" &
            .or.group_pep(ic)%gtype=="ASN".or.group_pep(ic)%gtype=="NASN".or.group_pep(ic)%gtype=="CASN".or.group_pep(ic)%gtype=="GLN".or.group_pep(ic)%gtype=="NGLN".or.group_pep(ic)%gtype=="CGLN" &
            .or.group_pep(ic)%gtype=="HIE".or.group_pep(ic)%gtype=="NHIE".or.group_pep(ic)%gtype=="CHIE".or.group_pep(ic)%gtype=="PRO".or.group_pep(ic)%gtype=="NPRO".or.group_pep(ic)%gtype=="CPRO" &
            .or.group_pep(ic)%gtype=="CYT".or.group_pep(ic)%gtype=="NCYT".or.group_pep(ic)%gtype=="CCYT".or.group_pep(ic)%gtype=="ALA".or.group_pep(ic)%gtype=="NALA".or.group_pep(ic)%gtype=="CALA" &
            .or.group_pep(ic)%gtype=="GLU".or.group_pep(ic)%gtype=="NGLU".or.group_pep(ic)%gtype=="CGLU".or.group_pep(ic)%gtype=="ASP".or.group_pep(ic)%gtype=="NASP".or.group_pep(ic)%gtype=="CASP") then
        else
            !open(20, file="error.txt", access="append")
                write(*,*) group_pep(ic)%gtype, "of the", ic, "residue could not be found in the PH of", ph_value
                write(*,*) "Please check whether the group_pep name is right or not!"
            !close(20)
            stop
        endif
    enddo

endif

return
end subroutine PH_checking

subroutine rotamerlib
    !loads rotamers for all amino acids 
    !note that N- and C-terminus amino acids are not loaded separately
implicit none
integer                     :: status, num_sc_angles, rotanum, anum, num, i, j
real                        :: x, y, z
integer                     :: cur_dihedrals(5)
character*4                 :: char, atype, name

!for the moment, I am not bothering to fix this function to avoid loading in rotamers for excluded amino acids, since this is only done once during the program and isn't time consuming.
if(ph_value.le.3.9) then
   open(10, file=trim(lib_path)//"/rotamer", status="old")
      do while(.true.)
         read(10, "(a,i3,i3)", iostat=status) name, num_sc_angles, rotanum
         if(status.ne.0) exit
         i=0
         if(name=="GLY") then
            i=1
         elseif(name=="LEU") then
            i=2
         elseif(name=="VAL") then
            i=3
         elseif(name=="ILE") then
            i=4
         elseif(name=="MET") then
            i=5
         elseif(name=="PHE") then
            i=6
         elseif(name=="TYR") then
            i=7
         elseif(name=="TRP") then
            i=8
         elseif(name=="ARG") then
            i=9
         elseif(name=="LYS") then
            i=10
         elseif(name=="SER") then
            i=11
         elseif(name=="THR") then
            i=12
         elseif(name=="ASN") then
            i=13
         elseif(name=="GLN") then
            i=14
         elseif(name=="HIP") then
            i=15
         elseif(name=="PRO") then
            i=16
         elseif(name=="CYS") then
            i=17
         elseif(name=="ALA") then
            i=18
         elseif(name=="GLH") then
            i=19
         elseif(name=="ASH") then
            i=20
         endif

         if(i.ne.0) then
            aa_lib(i)%gtype=name; aa_lib(i)%num_sc_angles=num_sc_angles; aa_lib(i)%rotanum=rotanum
            if(num_sc_angles.ne.0) then
               do j=1, rotanum
                  read(10,*) cur_dihedrals(1:num_sc_angles)
                  aa_lib(i)%dihedralangle(1:num_sc_angles, j) = cur_dihedrals(1:num_sc_angles)
               enddo
            endif
         endif
      enddo
   close(10)
elseif(ph_value.gt.3.9.and.ph_value.le.4.3) then
   open(10, file=trim(lib_path)//"/rotamer", status="old")
      do while(.true.)
         read(10, "(a,i3,i3)", iostat=status) name, num_sc_angles, rotanum
         if(status.ne.0) exit
         i=0
         if(name=="GLY") then
            i=1
         elseif(name=="LEU") then
            i=2
         elseif(name=="VAL") then
            i=3
         elseif(name=="ILE") then
            i=4
         elseif(name=="MET") then
            i=5
         elseif(name=="PHE") then
            i=6
         elseif(name=="TYR") then
            i=7
         elseif(name=="TRP") then
            i=8
         elseif(name=="ARG") then
            i=9
         elseif(name=="LYS") then
            i=10
         elseif(name=="SER") then
            i=11
         elseif(name=="THR") then
            i=12
         elseif(name=="ASN") then
            i=13
         elseif(name=="GLN") then
            i=14
         elseif(name=="HIP") then
            i=15
         elseif(name=="PRO") then
            i=16
         elseif(name=="CYS") then
            i=17
         elseif(name=="ALA") then
            i=18
         elseif(name=="GLH") then
            i=19
         elseif(name=="ASP") then
            i=20
         endif

         if(i.ne.0) then
            aa_lib(i)%gtype=name; aa_lib(i)%num_sc_angles=num_sc_angles; aa_lib(i)%rotanum=rotanum
            if(num_sc_angles.ne.0) then
               do j=1, rotanum
                  read(10,*) cur_dihedrals(1:num_sc_angles)
                  aa_lib(i)%dihedralangle(1:num_sc_angles, j) = cur_dihedrals(1:num_sc_angles)
               enddo
            endif
         endif
      enddo
   close(10)
elseif(ph_value.gt.4.3.and.ph_value.le.6.0) then
   open(10, file=trim(lib_path)//"/rotamer", status="old")
      do while(.true.)
         read(10, "(a,i3,i3)", iostat=status) name, num_sc_angles, rotanum
         if(status.ne.0) exit
         i=0
         if(name=="GLY") then
            i=1
         elseif(name=="LEU") then
            i=2
         elseif(name=="VAL") then
            i=3
         elseif(name=="ILE") then
            i=4
         elseif(name=="MET") then
            i=5
         elseif(name=="PHE") then
            i=6
         elseif(name=="TYR") then
            i=7
         elseif(name=="TRP") then
            i=8
         elseif(name=="ARG") then
            i=9
         elseif(name=="LYS") then
            i=10
         elseif(name=="SER") then
            i=11
         elseif(name=="THR") then
            i=12
         elseif(name=="ASN") then
            i=13
         elseif(name=="GLN") then
            i=14
         elseif(name=="HIP") then
            i=15
         elseif(name=="PRO") then
            i=16
         elseif(name=="CYS") then
            i=17
         elseif(name=="ALA") then
            i=18
         elseif(name=="GLU") then
            i=19
         elseif(name=="ASP") then
            i=20
         endif

         if(i.ne.0) then
            aa_lib(i)%gtype=name; aa_lib(i)%num_sc_angles=num_sc_angles; aa_lib(i)%rotanum=rotanum
            if(num_sc_angles.ne.0) then
               do j=1, rotanum
                  read(10,*) cur_dihedrals(1:num_sc_angles)
                  aa_lib(i)%dihedralangle(1:num_sc_angles, j) = cur_dihedrals(1:num_sc_angles)
               enddo
            endif
         endif
      enddo
   close(10)
elseif(ph_value.gt.6.0.and.ph_value.lt.8.3) then
   open(10, file=trim(lib_path)//"/rotamer", status="old")
      do while(.true.)
         read(10, "(a,i3,i3)", iostat=status) name, num_sc_angles, rotanum
         if(status.ne.0) exit
         i=0
         if(name=="GLY") then
            i=1
         elseif(name=="LEU") then
            i=2
         elseif(name=="VAL") then
            i=3
         elseif(name=="ILE") then
            i=4
         elseif(name=="MET") then
            i=5
         elseif(name=="PHE") then
            i=6
         elseif(name=="TYR") then
            i=7
         elseif(name=="TRP") then
            i=8
         elseif(name=="ARG") then
            i=9
         elseif(name=="LYS") then
            i=10
         elseif(name=="SER") then
            i=11
         elseif(name=="THR") then
            i=12
         elseif(name=="ASN") then
            i=13
         elseif(name=="GLN") then
            i=14
         elseif(name=="HIE") then
            i=15
         elseif(name=="PRO") then
            i=16
         elseif(name=="CYS") then
            i=17
         elseif(name=="ALA") then
            i=18
         elseif(name=="GLU") then
            i=19
         elseif(name=="ASP") then
            i=20
         endif

         if(i.ne.0) then
            aa_lib(i)%gtype=name; aa_lib(i)%num_sc_angles=num_sc_angles; aa_lib(i)%rotanum=rotanum
            if(num_sc_angles.ne.0) then
               do j=1, rotanum
                  read(10,*) cur_dihedrals(1:num_sc_angles)
                  aa_lib(i)%dihedralangle(1:num_sc_angles, j) = cur_dihedrals(1:num_sc_angles)
               enddo
            endif
         endif
      enddo
   close(10)
elseif(ph_value.ge.8.3.and.ph_value.lt.10.1) then
   open(10, file=trim(lib_path)//"/rotamer", status="old")
      do while(.true.)
         read(10, "(a,i3,i3)", iostat=status) name, num_sc_angles, rotanum
         if(status.ne.0) exit
         i=0
         if(name=="GLY") then
            i=1
         elseif(name=="LEU") then
            i=2
         elseif(name=="VAL") then
            i=3
         elseif(name=="ILE") then
            i=4
         elseif(name=="MET") then
            i=5
         elseif(name=="PHE") then
            i=6
         elseif(name=="TYR") then
            i=7
         elseif(name=="TRP") then
            i=8
         elseif(name=="ARG") then
            i=9
         elseif(name=="LYS") then
            i=10
         elseif(name=="SER") then
            i=11
         elseif(name=="THR") then
            i=12
         elseif(name=="ASN") then
            i=13
         elseif(name=="GLN") then
            i=14
         elseif(name=="HIE") then
            i=15
         elseif(name=="PRO") then
            i=16
         elseif(name=="CYT") then
            i=17
         elseif(name=="ALA") then
            i=18
         elseif(name=="GLU") then
            i=19
         elseif(name=="ASP") then
            i=20
         endif

         if(i.ne.0) then
            aa_lib(i)%gtype=name; aa_lib(i)%num_sc_angles=num_sc_angles; aa_lib(i)%rotanum=rotanum
            if(num_sc_angles.ne.0) then
               do j=1, rotanum
                  read(10,*) cur_dihedrals(1:num_sc_angles)
                  aa_lib(i)%dihedralangle(1:num_sc_angles, j) = cur_dihedrals(1:num_sc_angles)
               enddo
            endif
         endif
      enddo
   close(10)
elseif(ph_value.ge.10.1.and.ph_value.lt.10.5) then
   open(10, file=trim(lib_path)//"/rotamer", status="old")
      do while(.true.)
         read(10, "(a,i3,i3)", iostat=status) name, num_sc_angles, rotanum
         if(status.ne.0) exit
         i=0
         if(name=="GLY") then
            i=1
         elseif(name=="LEU") then
            i=2
         elseif(name=="VAL") then
            i=3
         elseif(name=="ILE") then
            i=4
         elseif(name=="MET") then
            i=5
         elseif(name=="PHE") then
            i=6
         elseif(name=="TYX") then
            i=7
         elseif(name=="TRP") then
            i=8
         elseif(name=="ARG") then
            i=9
         elseif(name=="LYS") then
            i=10
         elseif(name=="SER") then
            i=11
         elseif(name=="THR") then
            i=12
         elseif(name=="ASN") then
            i=13
         elseif(name=="GLN") then
            i=14
         elseif(name=="HIE") then
            i=15
         elseif(name=="PRO") then
            i=16
         elseif(name=="CYT") then
            i=17
         elseif(name=="ALA") then
            i=18
         elseif(name=="GLU") then
            i=19
         elseif(name=="ASP") then
            i=20
         endif

         if(i.ne.0) then
            aa_lib(i)%gtype=name; aa_lib(i)%num_sc_angles=num_sc_angles; aa_lib(i)%rotanum=rotanum
            if(num_sc_angles.ne.0) then
               do j=1, rotanum
                  read(10,*) cur_dihedrals(1:num_sc_angles)
                  aa_lib(i)%dihedralangle(1:num_sc_angles, j) = cur_dihedrals(1:num_sc_angles)
               enddo
            endif
         endif
      enddo
   close(10)
elseif(ph_value.ge.10.5.and.ph_value.lt.12.5) then
   open(10, file=trim(lib_path)//"/rotamer", status="old")
      do while(.true.)
         read(10, "(a,i3,i3)", iostat=status) name, num_sc_angles, rotanum
         if(status.ne.0) exit
         i=0
         if(name=="GLY") then
            i=1
         elseif(name=="LEU") then
            i=2
         elseif(name=="VAL") then
            i=3
         elseif(name=="ILE") then
            i=4
         elseif(name=="MET") then
            i=5
         elseif(name=="PHE") then
            i=6
         elseif(name=="TYX") then
            i=7
         elseif(name=="TRP") then
            i=8
         elseif(name=="ARG") then
            i=9
         elseif(name=="LYN") then
            i=10
         elseif(name=="SER") then
            i=11
         elseif(name=="THR") then
            i=12
         elseif(name=="ASN") then
            i=13
         elseif(name=="GLN") then
            i=14
         elseif(name=="HIE") then
            i=15
         elseif(name=="PRO") then
            i=16
         elseif(name=="CYT") then
            i=17
         elseif(name=="ALA") then
            i=18
         elseif(name=="GLU") then
            i=19
         elseif(name=="ASP") then
            i=20
         endif

         if(i.ne.0) then
            aa_lib(i)%gtype=name; aa_lib(i)%num_sc_angles=num_sc_angles; aa_lib(i)%rotanum=rotanum
            if(num_sc_angles.ne.0) then
               do j=1, rotanum
                  read(10,*) cur_dihedrals(1:num_sc_angles)
                  aa_lib(i)%dihedralangle(1:num_sc_angles, j) = cur_dihedrals(1:num_sc_angles)
               enddo
            endif
         endif
      enddo
   close(10)
elseif(ph_value.ge.12.5) then
   open(10, file=trim(lib_path)//"/rotamer", status="old")
      do while(.true.)
         read(10, "(a,i3,i3)", iostat=status) name, num_sc_angles, rotanum
         if(status.ne.0) exit
         i=0
         if(name=="GLY") then
            i=1
         elseif(name=="LEU") then
            i=2
         elseif(name=="VAL") then
            i=3
         elseif(name=="ILE") then
            i=4
         elseif(name=="MET") then
            i=5
         elseif(name=="PHE") then
            i=6
         elseif(name=="TYX") then
            i=7
         elseif(name=="TRP") then
            i=8
         elseif(name=="ARN") then
            i=9
         elseif(name=="LYN") then
            i=10
         elseif(name=="SER") then
            i=11
         elseif(name=="THR") then
            i=12
         elseif(name=="ASN") then
            i=13
         elseif(name=="GLN") then
            i=14
         elseif(name=="HIE") then
            i=15
         elseif(name=="PRO") then
            i=16
         elseif(name=="CYT") then
            i=17
         elseif(name=="ALA") then
            i=18
         elseif(name=="GLU") then
            i=19
         elseif(name=="ASP") then
            i=20
         endif

         if(i.ne.0) then
            aa_lib(i)%gtype=name; aa_lib(i)%num_sc_angles=num_sc_angles; aa_lib(i)%rotanum=rotanum
            if(num_sc_angles.ne.0) then
               do j=1, rotanum
                  read(10,*) cur_dihedrals(1:num_sc_angles)
                  aa_lib(i)%dihedralangle(1:num_sc_angles, j) = cur_dihedrals(1:num_sc_angles)
               enddo
            endif
         endif
      enddo
   close(10)
endif

aa_lib%cnum1=0
aa_lib%cnum2=0
aa_lib%cnum3=0

do i=1, 20
   open (10, file=trim(lib_path)//'/RotamerLibrary/'//trim(aa_lib(i)%gtype), status="old")
   do while(.true.)
      read(10, *, iostat=status) char, anum, atype, name, char, num, x, y, z
      if(status.ne.0) exit
      if(atype=="N".or.atype=="H".or.atype=="H1".or.atype=="H2".or.atype=="H3".or.atype=="CA" &
         .or.atype=="HA".or.atype=="HA2".or.atype=="HA3") then
         aa_lib(i)%cnum1=aa_lib(i)%cnum1+1
         aa_lib(i)%atype1(aa_lib(i)%cnum1)=atype
         aa_lib(i)%coo1(1,aa_lib(i)%cnum1)=x
         aa_lib(i)%coo1(2,aa_lib(i)%cnum1)=y
         aa_lib(i)%coo1(3,aa_lib(i)%cnum1)=z
      elseif(atype=="C".or.atype=="O".or.atype=="OXT") then
         aa_lib(i)%cnum3=aa_lib(i)%cnum3+1
         aa_lib(i)%atype3(aa_lib(i)%cnum3)=atype
         aa_lib(i)%coo3(1,aa_lib(i)%cnum3)=x
         aa_lib(i)%coo3(2,aa_lib(i)%cnum3)=y
         aa_lib(i)%coo3(3,aa_lib(i)%cnum3)=z
      else
         aa_lib(i)%cnum2=aa_lib(i)%cnum2+1
         aa_lib(i)%atype2(aa_lib(i)%cnum2)=atype
         aa_lib(i)%coo2(1,aa_lib(i)%cnum2)=x
         aa_lib(i)%coo2(2,aa_lib(i)%cnum2)=y
         aa_lib(i)%coo2(3,aa_lib(i)%cnum2)=z
      endif
   end do
   close(10)
enddo

return
end subroutine rotamerlib

subroutine findrotamer(res_i, group_pep, name_original, rotanum, aa_group)
!this subroutine gets the coordinates for all rotamers for amino acid and translates/rotates them to fit current position of peptide
implicit none
integer                         :: res_i, rotanum, i, ip
integer                         :: num_sc_angles, sc_groups(6), monitor(6)
real                            :: nr(3), car(3), cr(3), r_norpep(3), dihedrals(6)
real                            :: aa_nr(3), aa_car(3), aa_cr(3), r_norrot(3)
real                            :: r_nca(3), aa_r_nca(3), r_trans(3)
real                            :: CA(3), m(3,3)
character*4                     :: name_original
type(groupdetails)              :: group_pep(pep_res), aa_group(40)
type(index4sidechain)           :: index(60)
type(conformer4sidechain)       :: Iclass(6), Tclass(6)

!determine what amino acid we have
if(name_original=="GLY".or.name_original=="NGLY".or.name_original=="CGLY") then
    ip=1
elseif(name_original=="LEU".or.name_original=="NLEU".or.name_original=="CLEU") then
    ip=2
elseif(name_original=="VAL".or.name_original=="NVAL".or.name_original=="CVAL") then
    ip=3
elseif(name_original=="ILE".or.name_original=="NILE".or.name_original=="CILE") then
    ip=4
elseif(name_original=="MET".or.name_original=="NMET".or.name_original=="CMET") then
    ip=5
elseif(name_original=="PHE".or.name_original=="NPHE".or.name_original=="CPHE") then
    ip=6
elseif(name_original=="TYR".or.name_original=="NTYR".or.name_original=="CTYR".or. &
      name_original=="TYX".or.name_original=="NTYX".or.name_original=="CTYX") then
    ip=7
elseif(name_original=="TRP".or.name_original=="NTRP".or.name_original=="CTRP") then
    ip=8
elseif(name_original=="ARG".or.name_original=="NARG".or.name_original=="CARG".or. &
      name_original=="ARN".or.name_original=="NARN".or.name_original=="CARN") then
    ip=9
elseif(name_original=="LYS".or.name_original=="NLYS".or.name_original=="CLYS".or. &
      name_original=="LYN".or.name_original=="NLYN".or.name_original=="CLYN") then
    ip=10
elseif(name_original=="SER".or.name_original=="NSER".or.name_original=="CSER") then
    ip=11
elseif(name_original=="THR".or.name_original=="NTHR".or.name_original=="CTHR") then
    ip=12
elseif(name_original=="ASN".or.name_original=="NASN".or.name_original=="CASN") then
    ip=13
elseif(name_original=="GLN".or.name_original=="NGLN".or.name_original=="CGLN") then
    ip=14
elseif(name_original=="HIE".or.name_original=="NHIE".or.name_original=="CHIE".or. &
      name_original=="HIP".or.name_original=="NHIP".or.name_original=="CHIP") then
    ip=15
elseif(name_original=="PRO".or.name_original=="NPRO".or.name_original=="CPRO") then
    ip=16
elseif(name_original=="CYS".or.name_original=="NCYS".or.name_original=="CCYS".or. &
      name_original=="CYT".or.name_original=="NCYT".or.name_original=="CCYT") then
    ip=17
elseif(name_original=="ALA".or.name_original=="NALA".or.name_original=="CALA") then
    ip=18
elseif(name_original=="GLU".or.name_original=="NGLU".or.name_original=="CGLU".or. &
      name_original=="GLH".or.name_original=="NGLH".or.name_original=="CGLH") then
    ip=19
elseif(name_original=="ASP".or.name_original=="NASP".or.name_original=="CASP".or. &
      name_original=="ASH".or.name_original=="NASH".or.name_original=="CASH") then
    ip=20
endif

rotanum=aa_lib(ip)%rotanum !get number of rotamers for this amino acid

!get information for each rotamer: number of atoms in each group, the amino acid name, etc. Note that cnum1, cnum3, atype1, and atype3 not accurate for N- and C-terminal residues (this doesn't affect program, these values not ever put into group_pep variable)
do i=1, rotanum
    aa_group(i)%cnum1=aa_lib(ip)%cnum1; aa_group(i)%cnum2=aa_lib(ip)%cnum2; aa_group(i)%cnum3=aa_lib(ip)%cnum3  
    aa_group(i)%gtype=name_original
    aa_group(i)%atype1=aa_lib(ip)%atype1; aa_group(i)%atype2=aa_lib(ip)%atype2; aa_group(i)%atype3=aa_lib(ip)%atype3
enddo
aa_group(1)%coo1=aa_lib(ip)%coo1; aa_group(1)%coo2=aa_lib(ip)%coo2; aa_group(1)%coo3=aa_lib(ip)%coo3

!get position of N, alpha carbon, and carboxyl carbon to correct position rotamers in peptide
do i=1, group_pep(res_i)%cnum1
    if(group_pep(res_i)%atype1(i)=="N") then
        nr(1:3)=group_pep(res_i)%coo1(1:3,i)
    elseif(group_pep(res_i)%atype1(i)=="CA") then
        car(1:3)=group_pep(res_i)%coo1(1:3,i)
    endif
enddo
do i=1, group_pep(res_i)%cnum3
    if(group_pep(res_i)%atype3(i)=="C") then
        cr(1:3)=group_pep(res_i)%coo3(1:3,i)
        exit
    endif
enddo

!get vector normal to plane formed by N and two carbons
call normalvector(nr, car, cr, r_norpep)

!get vector pointing from alpha carbon to nitrogen
r_nca(1:3)=nr(1:3)-car(1:3)

!!STEP ONE: get rotamer into same plane as peptide
!for first rotamer, get the N, CA, and C coordinates from first rotamer in rotamer library file
do i=1, aa_group(1)%cnum1
    if(aa_group(1)%atype1(i)=="N") then
        aa_nr(:)=aa_group(1)%coo1(:,i)
    elseif(aa_group(1)%atype1(i)=="CA") then
        aa_car(:)=aa_group(1)%coo1(:,i)
    endif
enddo
do i=1, aa_group(1)%cnum3
    if(aa_group(1)%atype3(i)=="C") then
        aa_cr(:)=aa_group(1)%coo3(:,i)
    endif
enddo

!get vector normal to N, CA, and C coordinates of rotamer
call normalvector(aa_nr, aa_car, aa_cr, r_norrot)

!get matrix to rotate rotamer vector into reference frame of peptide
call vectorrotation(r_norrot, r_norpep, m)

!rotate rotamer coordinates into reference frame of peptide
aa_group(1)%coo1=matmul(m,aa_group(1)%coo1)
aa_group(1)%coo2=matmul(m,aa_group(1)%coo2)
aa_group(1)%coo3=matmul(m,aa_group(1)%coo3)

!!STEP TWO: get the N-CA bond of rotamer to be parallel with that of peptide
!get coordinates of N and CA of rotamer after the rotation
do i=1, aa_group(1)%cnum1
    if(aa_group(1)%atype1(i)=="N") then
        aa_nr(:)=aa_group(1)%coo1(:,i)
    elseif(aa_group(1)%atype1(i)=="CA") then
        aa_car(:)=aa_group(1)%coo1(:,i)
    endif
enddo

!get vector pointing from alpha carbon to nitrogen in the rotamer
aa_r_nca(1:3)=aa_nr(1:3)-aa_car(1:3)

!get rotation matrix to rotate rotamer vector to be parallel to peptide vector
call vectorrotation(aa_r_nca, r_nca, m)

!perform the rotation
aa_group(1)%coo1=matmul(m, aa_group(1)%coo1)
aa_group(1)%coo2=matmul(m, aa_group(1)%coo2)
aa_group(1)%coo3=matmul(m, aa_group(1)%coo3)

!!STEP THREE: translate the rotamer so it's CA is coincident with CA in peptide.
!!After this, the rotamer has been successfully moved into the reference frame of peptide
!find alpha carbon of rotamer after rotation
do i=1, aa_group(1)%cnum1
    if(aa_group(1)%atype1(i)=="CA") then
        aa_car(1:3)=aa_group(1)%coo1(1:3,i)
    endif
enddo

!get vector to translate alpha carbon of rotamer to alpha carbon of peptide
r_trans(1:3)=car(1:3)-aa_car(1:3)

!translate rotamer atoms by translation vector
do i=1, aa_group(1)%cnum1
    aa_group(1)%coo1(1:3,i)=aa_group(1)%coo1(1:3,i)+r_trans(1:3)
enddo
do i=1, aa_group(1)%cnum2
     aa_group(1)%coo2(1:3,i)=aa_group(1)%coo2(1:3,i)+r_trans(1:3)
enddo
do i=1, aa_group(1)%cnum3
    aa_group(1)%coo3(1:3,i)=aa_group(1)%coo3(1:3,i)+r_trans(1:3)
enddo

if(rotanum.le.1) goto 10 !nothing more to do here, we're done

!Next, partition the sidechain atoms into groups and find the atoms that define the sidechain dihedrals.
!This is specific to each amino acid, which is why this code block is ENORMOUS
!notes for here: num_sc_angles is number of dihedrals in sidechain, sc_groups counts atoms in each section separated by dihedrals, Iclass stores coordinates for sidechain, and monitor stores index to atoms that define dihedrals
sc_groups=0
if(aa_group(1)%gtype=="VAL".or.aa_group(1)%gtype=="NVAL".or.aa_group(1)%gtype=="CVAL") then
    num_sc_angles=1
    do i=1, aa_group(1)%cnum2
        if(aa_group(1)%atype2(i)=="CB") then
            sc_groups(1)=sc_groups(1)+1
            Iclass(1)%member(1:3,sc_groups(1))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=1; index(i)%member_No=sc_groups(1)
            monitor(1)=sc_groups(1)
        else
            sc_groups(2)=sc_groups(2)+1
            Iclass(2)%member(1:3,sc_groups(2))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=2; index(i)%member_No=sc_groups(2)
        endif
    enddo
elseif(aa_group(1)%gtype=="LEU".or.aa_group(1)%gtype=="NLEU".or.aa_group(1)%gtype=="CLEU") then
    num_sc_angles=2
    do i=1, aa_group(1)%cnum2
        if(aa_group(1)%atype2(i)=="CB") then
            sc_groups(1)=sc_groups(1)+1
            Iclass(1)%member(1:3,sc_groups(1))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=1; index(i)%member_No=sc_groups(1)
            monitor(1)=sc_groups(1)
        elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
            sc_groups(2)=sc_groups(2)+1
            Iclass(2)%member(1:3,sc_groups(2))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=2; index(i)%member_No=sc_groups(2)
            if(aa_group(1)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        else
            sc_groups(3)=sc_groups(3)+1
            Iclass(3)%member(1:3,sc_groups(3))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        endif
    enddo
elseif(aa_group(1)%gtype=="ILE".or.aa_group(1)%gtype=="NILE".or.aa_group(1)%gtype=="CILE") then
    num_sc_angles=2
    do i=1, aa_group(1)%cnum2
        if(aa_group(1)%atype2(i)=="CB") then
            sc_groups(1)=sc_groups(1)+1
            Iclass(1)%member(1:3,sc_groups(1))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=1; index(i)%member_No=sc_groups(1)
            monitor(1)=sc_groups(1)
        elseif(aa_group(1)%atype2(i)=="HB".or.aa_group(1)%atype2(i)=="CG2".or.aa_group(1)%atype2(i)=="HG21".or. &
            aa_group(1)%atype2(i)=="HG22".or.aa_group(1)%atype2(i)=="HG23".or.aa_group(1)%atype2(i)=="CG1") then
            sc_groups(2)=sc_groups(2)+1
            Iclass(2)%member(1:3,sc_groups(2))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=2; index(i)%member_No=sc_groups(2)
            if(aa_group(1)%atype2(i)=="CG1") monitor(2)=sc_groups(2)
        else
            sc_groups(3)=sc_groups(3)+1
            Iclass(3)%member(1:3,sc_groups(3))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        endif
    enddo
elseif(aa_group(1)%gtype=="PHE".or.aa_group(1)%gtype=="NPHE".or.aa_group(1)%gtype=="CPHE") then
    num_sc_angles=2
    do i=1, aa_group(1)%cnum2
        if(aa_group(1)%atype2(i)=="CB") then
            sc_groups(1)=sc_groups(1)+1
            Iclass(1)%member(1:3,sc_groups(1))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=1; index(i)%member_No=sc_groups(1)
            monitor(1)=sc_groups(1)
        elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
            sc_groups(2)=sc_groups(2)+1
            Iclass(2)%member(1:3,sc_groups(2))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=2; index(i)%member_No=sc_groups(2)
            if(aa_group(1)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        else
            sc_groups(3)=sc_groups(3)+1
            Iclass(3)%member(1:3,sc_groups(3))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        endif
    enddo
elseif(aa_group(1)%gtype=="TRP".or.aa_group(1)%gtype=="NTRP".or.aa_group(1)%gtype=="CTRP") then
    num_sc_angles=2
    do i=1, aa_group(1)%cnum2
        if(aa_group(1)%atype2(i)=="CB") then
            sc_groups(1)=sc_groups(1)+1
            Iclass(1)%member(1:3,sc_groups(1))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=1; index(i)%member_No=sc_groups(1)
            monitor(1)=sc_groups(1)
        elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
            sc_groups(2)=sc_groups(2)+1
            Iclass(2)%member(1:3,sc_groups(2))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=2; index(i)%member_No=sc_groups(2)
            if(aa_group(1)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        else
            sc_groups(3)=sc_groups(3)+1
            Iclass(3)%member(1:3,sc_groups(3))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        endif
    enddo
elseif(aa_group(1)%gtype=="TYR".or.aa_group(1)%gtype=="NTYR".or.aa_group(1)%gtype=="CTYR".or. &
      aa_group(1)%gtype=="TYX".or.aa_group(1)%gtype=="NTYX".or.aa_group(1)%gtype=="CTYX") then
    num_sc_angles=3
    do i=1, aa_group(1)%cnum2
        if(aa_group(1)%atype2(i)=="CB") then
            sc_groups(1)=sc_groups(1)+1
            Iclass(1)%member(1:3,sc_groups(1))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=1; index(i)%member_No=sc_groups(1)
            monitor(1)=sc_groups(1)
        elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
            sc_groups(2)=sc_groups(2)+1
            Iclass(2)%member(1:3,sc_groups(2))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=2; index(i)%member_No=sc_groups(2)
            if(aa_group(1)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        elseif(aa_group(1)%atype2(i)=="HH") then
            sc_groups(4)=sc_groups(4)+1
            Iclass(4)%member(1:3,sc_groups(4))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=4; index(i)%member_No=sc_groups(4)
        else
            sc_groups(3)=sc_groups(3)+1
            Iclass(3)%member(1:3,sc_groups(3))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=3; index(i)%member_No=sc_groups(3)
            if(aa_group(1)%atype2(i)=="OH") monitor(3)=sc_groups(3)
        endif
    enddo
elseif(aa_group(1)%gtype=="SER".or.aa_group(1)%gtype=="NSER".or.aa_group(1)%gtype=="CSER") then
    num_sc_angles=2
    do i=1, aa_group(1)%cnum2
        if(aa_group(1)%atype2(i)=="CB") then
            sc_groups(1)=sc_groups(1)+1
            Iclass(1)%member(1:3,sc_groups(1))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=1; index(i)%member_No=sc_groups(1)
            monitor(1)=sc_groups(1)
        elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="OG") then
            sc_groups(2)=sc_groups(2)+1
            Iclass(2)%member(1:3,sc_groups(2))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=2; index(i)%member_No=sc_groups(2)
            if(aa_group(1)%atype2(i)=="OG") monitor(2)=sc_groups(2)
        else
            sc_groups(3)=sc_groups(3)+1
            Iclass(3)%member(1:3,sc_groups(3))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        endif
    enddo
elseif(aa_group(1)%gtype=="THR".or.aa_group(1)%gtype=="NTHR".or.aa_group(1)%gtype=="CTHR") then
    num_sc_angles=2
    do i=1, aa_group(1)%cnum2
        if(aa_group(1)%atype2(i)=="CB") then
            sc_groups(1)=sc_groups(1)+1
            Iclass(1)%member(1:3,sc_groups(1))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=1; index(i)%member_No=sc_groups(1)
            monitor(1)=sc_groups(1)
        elseif(aa_group(1)%atype2(i)=="HG1") then
            sc_groups(3)=sc_groups(3)+1
            Iclass(3)%member(1:3,sc_groups(3))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        else
            sc_groups(2)=sc_groups(2)+1
            Iclass(2)%member(1:3,sc_groups(2))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=2; index(i)%member_No=sc_groups(2)
            if(aa_group(1)%atype2(i)=="OG1") monitor(2)=sc_groups(2)
        endif
    enddo
elseif(aa_group(1)%gtype=="CYS".or.aa_group(1)%gtype=="NCYS".or.aa_group(1)%gtype=="CCYS".or. &
      aa_group(1)%gtype=="CYT".or.aa_group(1)%gtype=="NCYT".or.aa_group(1)%gtype=="CCYT") then
    num_sc_angles=2
    do i=1, aa_group(1)%cnum2
        if(aa_group(1)%atype2(i)=="CB") then
            sc_groups(1)=sc_groups(1)+1
            Iclass(1)%member(1:3,sc_groups(1))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=1; index(i)%member_No=sc_groups(1)
            monitor(1)=sc_groups(1)
        elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="SG") then
            sc_groups(2)=sc_groups(2)+1
            Iclass(2)%member(1:3,sc_groups(2))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=2; index(i)%member_No=sc_groups(2)
            if(aa_group(1)%atype2(i)=="SG") monitor(2)=sc_groups(2)
        else
            sc_groups(3)=sc_groups(3)+1
            Iclass(3)%member(1:3,sc_groups(3))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        endif
    enddo
elseif(aa_group(1)%gtype=="MET".or.aa_group(1)%gtype=="NMET".or.aa_group(1)%gtype=="CMET") then
    num_sc_angles=4
    do i=1, aa_group(1)%cnum2
        if(aa_group(1)%atype2(i)=="CB") then
            sc_groups(1)=sc_groups(1)+1
            Iclass(1)%member(1:3,sc_groups(1))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=1; index(i)%member_No=sc_groups(1)
            monitor(1)=sc_groups(1)
        elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
            sc_groups(2)=sc_groups(2)+1
            Iclass(2)%member(1:3,sc_groups(2))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=2; index(i)%member_No=sc_groups(2)
            if(aa_group(1)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        elseif(aa_group(1)%atype2(i)=="HG2".or.aa_group(1)%atype2(i)=="HG3".or.aa_group(1)%atype2(i)=="SD") then
            sc_groups(3)=sc_groups(3)+1
            Iclass(3)%member(1:3,sc_groups(3))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=3; index(i)%member_No=sc_groups(3)
            if(aa_group(1)%atype2(i)=="SD") monitor(3)=sc_groups(3)
        elseif(aa_group(1)%atype2(i)=="CE") then
            sc_groups(4)=sc_groups(4)+1
            Iclass(4)%member(1:3,sc_groups(4))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=4; index(i)%member_No=sc_groups(4)
            monitor(4)=sc_groups(4)
        else
            sc_groups(5)=sc_groups(5)+1
            Iclass(5)%member(1:3,sc_groups(5))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=5; index(i)%member_No=sc_groups(5)
        endif
    enddo
elseif(aa_group(1)%gtype=="ASN".or.aa_group(1)%gtype=="NASN".or.aa_group(1)%gtype=="CASN") then
    num_sc_angles=2
    do i=1, aa_group(1)%cnum2
        if(aa_group(1)%atype2(i)=="CB") then
            sc_groups(1)=sc_groups(1)+1
            Iclass(1)%member(1:3,sc_groups(1))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=1; index(i)%member_No=sc_groups(1)
            monitor(1)=sc_groups(1)
        elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
            sc_groups(2)=sc_groups(2)+1
            Iclass(2)%member(1:3,sc_groups(2))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=2; index(i)%member_No=sc_groups(2)
            if(aa_group(1)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        else
            sc_groups(3)=sc_groups(3)+1
            Iclass(3)%member(1:3,sc_groups(3))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        endif
    enddo
elseif(aa_group(1)%gtype=="GLN".or.aa_group(1)%gtype=="NGLN".or.aa_group(1)%gtype=="CGLN") then
    num_sc_angles=3
    do i=1, aa_group(1)%cnum2
        if(aa_group(1)%atype2(i)=="CB") then
            sc_groups(1)=sc_groups(1)+1
            Iclass(1)%member(1:3,sc_groups(1))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=1; index(i)%member_No=sc_groups(1)
            monitor(1)=sc_groups(1)
        elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
            sc_groups(2)=sc_groups(2)+1
            Iclass(2)%member(1:3,sc_groups(2))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=2; index(i)%member_No=sc_groups(2)
            if(aa_group(1)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        elseif(aa_group(1)%atype2(i)=="HG2".or.aa_group(1)%atype2(i)=="HG3".or.aa_group(1)%atype2(i)=="CD") then
            sc_groups(3)=sc_groups(3)+1
            Iclass(3)%member(1:3,sc_groups(3))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=3; index(i)%member_No=sc_groups(3)
            if(aa_group(1)%atype2(i)=="CD") monitor(3)=sc_groups(3)
        else
            sc_groups(4)=sc_groups(4)+1
            Iclass(4)%member(1:3,sc_groups(4))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=4; index(i)%member_No=sc_groups(4)
        endif
    enddo
elseif(aa_group(1)%gtype=="ASP".or.aa_group(1)%gtype=="NASP".or.aa_group(1)%gtype=="CASP".or. &
      aa_group(1)%gtype=="ASH".or.aa_group(1)%gtype=="NASH".or.aa_group(1)%gtype=="CASH") then
    num_sc_angles=2
    do i=1, aa_group(1)%cnum2
        if(aa_group(1)%atype2(i)=="CB") then
            sc_groups(1)=sc_groups(1)+1
            Iclass(1)%member(1:3,sc_groups(1))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=1; index(i)%member_No=sc_groups(1)
            monitor(1)=sc_groups(1)
        elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
            sc_groups(2)=sc_groups(2)+1
            Iclass(2)%member(1:3,sc_groups(2))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=2; index(i)%member_No=sc_groups(2)
            if(aa_group(1)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        else
            sc_groups(3)=sc_groups(3)+1
            Iclass(3)%member(1:3,sc_groups(3))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        endif
    enddo
elseif(aa_group(1)%gtype=="GLU".or.aa_group(1)%gtype=="NGLU".or.aa_group(1)%gtype=="CGLU".or. &
      aa_group(1)%gtype=="GLH".or.aa_group(1)%gtype=="NGLH".or.aa_group(1)%gtype=="CGLH") then
    num_sc_angles=3
    do i=1, aa_group(1)%cnum2
        if(aa_group(1)%atype2(i)=="CB") then
            sc_groups(1)=sc_groups(1)+1
            Iclass(1)%member(1:3,sc_groups(1))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=1; index(i)%member_No=sc_groups(1)
            monitor(1)=sc_groups(1)
        elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
            sc_groups(2)=sc_groups(2)+1
            Iclass(2)%member(1:3,sc_groups(2))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=2; index(i)%member_No=sc_groups(2)
            if(aa_group(1)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        elseif(aa_group(1)%atype2(i)=="HG2".or.aa_group(1)%atype2(i)=="HG3".or.aa_group(1)%atype2(i)=="CD") then
            sc_groups(3)=sc_groups(3)+1
            Iclass(3)%member(1:3,sc_groups(3))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=3; index(i)%member_No=sc_groups(3)
            if(aa_group(1)%atype2(i)=="CD") monitor(3)=sc_groups(3)
        else
            sc_groups(4)=sc_groups(4)+1
            Iclass(4)%member(1:3,sc_groups(4))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=4; index(i)%member_No=sc_groups(4)
        endif
    enddo
elseif(aa_group(1)%gtype=="HIE".or.aa_group(1)%gtype=="NHIE".or.aa_group(1)%gtype=="CHIE".or. &
      aa_group(1)%gtype=="HIP".or.aa_group(1)%gtype=="NHIP".or.aa_group(1)%gtype=="CHIP") then
    num_sc_angles=2
    do i=1, aa_group(1)%cnum2
        if(aa_group(1)%atype2(i)=="CB") then
            sc_groups(1)=sc_groups(1)+1
            Iclass(1)%member(1:3,sc_groups(1))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=1; index(i)%member_No=sc_groups(1)
            monitor(1)=sc_groups(1)
        elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
            sc_groups(2)=sc_groups(2)+1
            Iclass(2)%member(1:3,sc_groups(2))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=2; index(i)%member_No=sc_groups(2)
            if(aa_group(1)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        else
            sc_groups(3)=sc_groups(3)+1
            Iclass(3)%member(1:3,sc_groups(3))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        endif
    enddo
elseif(aa_group(1)%gtype=="LYS".or.aa_group(1)%gtype=="NLYS".or.aa_group(1)%gtype=="CLYS".or. &
      aa_group(1)%gtype=="LYN".or.aa_group(1)%gtype=="NLYN".or.aa_group(1)%gtype=="CLYN") then
    num_sc_angles=5
    do i=1, aa_group(1)%cnum2
        if(aa_group(1)%atype2(i)=="CB") then
            sc_groups(1)=sc_groups(1)+1
            Iclass(1)%member(1:3,sc_groups(1))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=1; index(i)%member_No=sc_groups(1)
            monitor(1)=sc_groups(1)
        elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
            sc_groups(2)=sc_groups(2)+1
            Iclass(2)%member(1:3,sc_groups(2))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=2; index(i)%member_No=sc_groups(2)
            if(aa_group(1)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        elseif(aa_group(1)%atype2(i)=="HG2".or.aa_group(1)%atype2(i)=="HG3".or.aa_group(1)%atype2(i)=="CD") then
            sc_groups(3)=sc_groups(3)+1
            Iclass(3)%member(1:3,sc_groups(3))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=3; index(i)%member_No=sc_groups(3)
            if(aa_group(1)%atype2(i)=="CD") monitor(3)=sc_groups(3)
        elseif(aa_group(1)%atype2(i)=="HD2".or.aa_group(1)%atype2(i)=="HD3".or.aa_group(1)%atype2(i)=="CE") then
            sc_groups(4)=sc_groups(4)+1
            Iclass(4)%member(1:3,sc_groups(4))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=4; index(i)%member_No=sc_groups(4)
            if(aa_group(1)%atype2(i)=="CE") monitor(4)=sc_groups(4)
        elseif(aa_group(1)%atype2(i)=="HE2".or.aa_group(1)%atype2(i)=="HE3".or.aa_group(1)%atype2(i)=="NZ") then
            sc_groups(5)=sc_groups(5)+1
            Iclass(5)%member(1:3,sc_groups(5))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=5; index(i)%member_No=sc_groups(5)
            if(aa_group(1)%atype2(i)=="NZ") monitor(5)=sc_groups(5)
        else
            sc_groups(6)=sc_groups(6)+1
            Iclass(6)%member(1:3,sc_groups(6))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=6; index(i)%member_No=sc_groups(6)
        endif
    enddo
elseif(aa_group(1)%gtype=="ARG".or.aa_group(1)%gtype=="NARG".or.aa_group(1)%gtype=="CARG".or. &
      aa_group(1)%gtype=="ARN".or.aa_group(1)%gtype=="NARN".or.aa_group(1)%gtype=="CARN") then
    num_sc_angles=4
    do i=1, aa_group(1)%cnum2
        if(aa_group(1)%atype2(i)=="CB") then
            sc_groups(1)=sc_groups(1)+1
            Iclass(1)%member(1:3,sc_groups(1))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=1; index(i)%member_No=sc_groups(1)
            monitor(1)=sc_groups(1)
        elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
            sc_groups(2)=sc_groups(2)+1
            Iclass(2)%member(1:3,sc_groups(2))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=2; index(i)%member_No=sc_groups(2)
            if(aa_group(1)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        elseif(aa_group(1)%atype2(i)=="HG2".or.aa_group(1)%atype2(i)=="HG3".or.aa_group(1)%atype2(i)=="CD") then
            sc_groups(3)=sc_groups(3)+1
            Iclass(3)%member(1:3,sc_groups(3))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=3; index(i)%member_No=sc_groups(3)
            if(aa_group(1)%atype2(i)=="CD") monitor(3)=sc_groups(3)
        elseif(aa_group(1)%atype2(i)=="HD2".or.aa_group(1)%atype2(i)=="HD3".or.aa_group(1)%atype2(i)=="NE") then
            sc_groups(4)=sc_groups(4)+1
            Iclass(4)%member(1:3,sc_groups(4))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=4; index(i)%member_No=sc_groups(4)
            if(aa_group(1)%atype2(i)=="NE") monitor(4)=sc_groups(4)
        else
            sc_groups(5)=sc_groups(5)+1
            Iclass(5)%member(1:3,sc_groups(5))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=5; index(i)%member_No=sc_groups(5)
        endif
    enddo
elseif(aa_group(1)%gtype=="PRO".or.aa_group(1)%gtype=="NPRO".or.aa_group(1)%gtype=="CPRO") then
    num_sc_angles=3
    do i=1, aa_group(1)%cnum2
        if(aa_group(1)%atype2(i)=="CB") then
            sc_groups(1)=sc_groups(1)+1
            Iclass(1)%member(1:3,sc_groups(1))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=1; index(i)%member_No=sc_groups(1)
            monitor(1)=sc_groups(1)
        elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
            sc_groups(2)=sc_groups(2)+1
            Iclass(2)%member(1:3,sc_groups(2))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=2; index(i)%member_No=sc_groups(2)
            if(aa_group(1)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        elseif(aa_group(1)%atype2(i)=="HG2".or.aa_group(1)%atype2(i)=="HG3".or.aa_group(1)%atype2(i)=="CD") then
            sc_groups(3)=sc_groups(3)+1
            Iclass(3)%member(1:3,sc_groups(3))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=3; index(i)%member_No=sc_groups(3)
            if(aa_group(1)%atype2(i)=="CD") monitor(3)=sc_groups(3)
        else
            sc_groups(4)=sc_groups(4)+1
            Iclass(4)%member(1:3,sc_groups(4))=aa_group(1)%coo2(1:3,i)
            index(i)%class_No=4; index(i)%member_No=sc_groups(4)
        endif
    enddo
endif

!make sure we found all the rotamers
if(num_sc_angles.ne.aa_lib(ip)%num_sc_angles) then
    !open(10, file="error.txt", access="append")
        write(*,*) aa_lib(ip)%gtype
        write(*,*) "num_sc_angles=", num_sc_angles
        write(*,*) "aa_lib(",ip,")%num_sc_angles=", aa_lib(ip)%num_sc_angles
        write(*,*) "They are not equal with each other!"
    !close(10)
    stop
endif

!get alpha carbon position in peptide system
do i=1, aa_group(1)%cnum1
    if(aa_group(1)%atype1(i)=="CA") then
        CA(1:3)=aa_group(1)%coo1(1:3,i)
    endif
enddo

!calculate sidechain atom positions for all rotamers. Note that rotamer one is inherently taken care of by how rotamer file set up, so we start at rotamer two
!NOTE: The below rotation method called multiple places in PepBD, perhaps should split off as a subroutine and put in "utilities"
do i=2, rotanum
    !store basic info into aa_group. Most important here is storing the coordinate info
    aa_group(i)%cnum2=aa_group(1)%cnum2
    aa_group(i)%atype2=aa_group(1)%atype2
    aa_group(i)%coo2=aa_group(1)%coo2
    aa_group(i)%gtype=name_original

    !don't actually need the below values since we only want the sidechain coordinates/atom types/atom numbers
    !aa_group(i)%cnum1=aa_group(1)%cnum1; aa_group(i)%cnum3=aa_group(1)%cnum3
    !aa_group(i)%atype1=aa_group(1)%atype1;  aa_group(i)%atype3=aa_group(1)%atype3
    !aa_group(i)%coo1=aa_group(1)%coo1; aa_group(i)%coo3=aa_group(1)%coo3

   Tclass=Iclass !store sidechain positions from rotamer 1 in Tclass

   dihedrals=REAL(aa_lib(ip)%dihedralangle(:,i)-aa_lib(ip)%dihedralangle(:,1)) !get rotation angle from rotamer 1
   call sidechain_rotation(num_sc_angles, sc_groups, dihedrals, Tclass, CA, aa_group, monitor, index, i)  !store different sidechain conformations

enddo
10   continue

return
end subroutine findrotamer

subroutine ramachandranmap
implicit none
integer                         :: i, j

open(10, file=trim(lib_path)//"/rama_map", status="old")
   do j=-179, 180
      read(10,"(360i2)") (rama_map(i,j), i=-179, 180)
   end do
close(10)

flavoredregion_number=0

do i=-179, 180
   do j=-179, 180
      if(rama_map(i,j)==1.or.rama_map(i,j)==2) then
         flavoredregion_number=flavoredregion_number+1
         flavored_region(flavoredregion_number,1)=REAL(i)
         flavored_region(flavoredregion_number,2)=REAL(j)
      endif
   enddo
enddo

return
end subroutine ramachandranmap

subroutine load_LCPO_params()
implicit none
character*10       :: atype
integer            :: i, status, neighbors
real               :: r, p1, p2, p3, p4, p1_NLR, p2_NLR, p3_NLR, p4_NLR

open(10, file=trim(lib_path)//"/LCPO_Parameters", status="old")
read(10, *) !skip the header info in the file
do i=1, 20
   read(10, 20, iostat=status) atype, neighbors, r, p1, p2, p3, p4, p1_NLR, p2_NLR, p3_NLR, p4_NLR
   LCPO(i)%atype = trim(atype)
   LCPO(i)%n = neighbors
   LCPO(i)%r = r + probe_r
   LCPO(i)%p1 = p1
   LCPO(i)%p2 = p2
   LCPO(i)%p3 = p3
   LCPO(i)%p4 = p4
   LCPO(i)%p1_NLR = p1_NLR
   LCPO(i)%p2_NLR = p2_NLR
   LCPO(i)%p3_NLR = p3_NLR
   LCPO(i)%p4_NLR = p4_NLR
enddo

close(10)

20   format(a10, i10, e10.9, 9e12.10)

end subroutine load_LCPO_params

subroutine energy_parameter(group, group_para, group_rec, rec_para)
!This subroutine reads in the force field parameters for all residues in the system
implicit none
integer                     :: i
type(groupdetails)          :: group(pep_res), group_rec(rec_res)
type(energyparameters)      :: group_para(pep_res), rec_para(rec_res)

do i=1, pep_res
   call energy_parameter_single(group, group_para, i)
enddo

do i=1, rec_res
    call energy_parameter_single(group_rec, rec_para, i)
 enddo

return
end subroutine energy_parameter

subroutine energy_parameter_single(group, group_para, i)
!This subroutine reads in the force field parameters a single amino acid at position "i"
!subroutine also fixes the ordering of atoms in "group" to match the order in the parameter files
implicit none
integer                     :: i, j, status, atomid, index
real                        :: charge, epsion, r, rborn, fs, dielecons
character*4                 :: lbres, igraph
integer                     :: LCPO_index
type(groupdetails)          :: group(:), res_reordered
type(energyparameters)      :: group_para(:)

open(10, file=trim(lib_path)//'/ForceField/'//trim(group(i)%gtype), status="old")
read(10, *)
do while(.true.)
    read(10, 20, iostat=status) lbres, igraph, charge, epsion, r, rborn, fs, dielecons, atomid, LCPO_index
    if(status.lt.0) goto 30
    do j=1, group(i)%cnum1
        if(group(i)%atype1(j)==igraph) then
            !the below if accounts for atom numbering of proline being different from other residues
            if (group(i)%gtype.eq.'NPRO' .or. group(i)%gtype.eq.'PRO' .or. group(i)%gtype.eq.'CPRO') then
                if (group(i)%atype1(j).eq.'CA' .or. group(i)%atype1(j).eq.'HA') then
                    index = atomid-group(i)%cnum2
                else
                    index = atomid
                endif
            else
                index = atomid
            endif
            group_para(i)%charge1(index)=charge
            group_para(i)%epsion1(index)=epsion
            group_para(i)%r1(index)=r
            group_para(i)%rborn1(index)=rborn
            group_para(i)%fs1(index)=fs
            group_para(i)%dielecons1(index)=dielecons
            group_para(i)%atomid1(index)=atomid
            group_para(i)%LCPO1(index) = LCPO_index

            res_reordered%coo1(:,index) = group(i)%coo1(:,j)
            res_reordered%atype1(index) = group(i)%atype1(j)
            goto 40
        endif
    enddo
    do j=1, group(i)%cnum2
        if(group(i)%atype2(j)==igraph) then
            !the below if accounts for atom numbering of proline being different from other residues
            if (group(i)%gtype.eq.'NPRO' .or. group(i)%gtype.eq.'PRO' .or. group(i)%gtype.eq.'CPRO') then
                index = atomid-group(i)%cnum1+2
            else
                index = atomid-group(i)%cnum1
            endif
            group_para(i)%charge2(index)=charge
            group_para(i)%epsion2(index)=epsion
            group_para(i)%r2(index)=r
            group_para(i)%rborn2(index)=rborn
            group_para(i)%fs2(index)=fs
            group_para(i)%dielecons2(index)=dielecons
            group_para(i)%atomid2(index)=atomid
            group_para(i)%LCPO2(index) = LCPO_index

            res_reordered%coo2(:,index) = group(i)%coo2(:,j)
            res_reordered%atype2(index) = group(i)%atype2(j)
            goto 40
        endif
    enddo
    do j=1, group(i)%cnum3
        if(group(i)%atype3(j)==igraph) then
            index = atomid-group(i)%cnum1-group(i)%cnum2
            group_para(i)%charge3(index)=charge
            group_para(i)%epsion3(index)=epsion
            group_para(i)%r3(index)=r
            group_para(i)%rborn3(index)=rborn
            group_para(i)%fs3(index)=fs
            group_para(i)%dielecons3(index)=dielecons
            group_para(i)%atomid3(index)=atomid
            group_para(i)%LCPO3(index) = LCPO_index

            res_reordered%coo3(:,index) = group(i)%coo3(:,j)
            res_reordered%atype3(index) = group(i)%atype3(j)
            goto 40
        endif
    enddo

    write(*,"(A,A,A,A,A,I0)") "During parameter loading, unable to find atom ", igraph, " in ", group(i)%gtype, " of residue", i
    write(*,"(A)") "Issue is likely due to atom missing in pdbfile or not having a name that matches the library file"
    write(*,"(A)") "Terminating program!"
    stop
    40   continue
enddo

30   continue
close(10)

group(i)%coo1(:,1:group(i)%cnum1) = res_reordered%coo1(:,1:group(i)%cnum1)
group(i)%coo2(:,1:group(i)%cnum2) = res_reordered%coo2(:,1:group(i)%cnum2)
group(i)%coo3(:,1:group(i)%cnum3) = res_reordered%coo3(:,1:group(i)%cnum3)
group(i)%atype1(1:group(i)%cnum1) = res_reordered%atype1(1:group(i)%cnum1)
group(i)%atype2(1:group(i)%cnum2) = res_reordered%atype2(1:group(i)%cnum2)
group(i)%atype3(1:group(i)%cnum3) = res_reordered%atype3(1:group(i)%cnum3)

20   format(2a4, 6e16.8, i8, i8)
  
do j=1, group(i)%cnum1
   if(group_para(i)%dielecons1(j)<=0.1) then
      write(*,"(A,A,A,A,A,I0,A,F6.3)") "Atom ", group(i)%atype1(j), " in ", group(i)%gtype, " of residue", i, " has bad force field parameter", group_para(i)%dielecons1(j)
      write(*,"(A)") "Please check whether the atom type of PDB file matches the atom type of Force Field LIB or not!"
      stop
   endif
enddo
do j=1, group(i)%cnum2
   if(group_para(i)%dielecons2(j)<=0.1) then
      write(*,"(A,A,A,A,A,I0,A,F6.3)") "Atom ", group(i)%atype2(j), "in ", group(i)%gtype, " of residue", i, " has bad force field parameter", group_para(i)%dielecons2(j)
      write(*,"(A)") "Please check whether the atom type of PDB file matches the atom type of Force Field LIB or not!"
      stop
   endif
enddo
do j=1, group(i)%cnum3
   if(group_para(i)%dielecons3(j)<=0.1) then
      write(*,"(A,A,A,A,A,I0,A,F6.3)") "Atom ", group(i)%atype2(j), "in ", group(i)%gtype, " of residue", i, " has bad force field parameter", group_para(i)%dielecons3(j)
      write(*,"(A)") "Please check whether the atom type of PDB file matches the atom type of Force Field LIB or not!"
      stop
   endif
enddo  
return

end subroutine energy_parameter_single

subroutine find_hydration_category(current_AA, hyd_cat)
implicit none
character*4         :: current_AA
integer             :: hyd_cat, i

do hyd_cat =1,6
    do i=1, num_AAs_available(hyd_cat)
        if (current_AA .eq. trim(AAs_available(i,hyd_cat))) return
    enddo
enddo

end subroutine find_hydration_category

subroutine groupinfo(name, group_name, flag)
implicit none
integer                     :: i, flag
character*4                 :: name, group_name(3)

if(name=="GLY".or.name=="NGLY".or.name=="CGLY") then
    group_name(1)="GLY"
    group_name(2)="NGLY"
    group_name(3)="CGLY"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="LEU".or.name=="NLEU".or.name=="CLEU") then
    group_name(1)="LEU"
    group_name(2)="NLEU"
    group_name(3)="CLEU"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="VAL".or.name=="NVAL".or.name=="CVAL") then
    group_name(1)="VAL"
    group_name(2)="NVAL"
    group_name(3)="CVAL"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="ILE".or.name=="NILE".or.name=="CILE") then
    group_name(1)="ILE"
    group_name(2)="NILE"
    group_name(3)="CILE"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="MET".or.name=="NMET".or.name=="CMET") then
    group_name(1)="MET"
    group_name(2)="NMET"
    group_name(3)="CMET"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="PHE".or.name=="NPHE".or.name=="CPHE") then
    group_name(1)="PHE"
    group_name(2)="NPHE"
    group_name(3)="CPHE"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="TYR".or.name=="NTYR".or.name=="CTYR") then
    group_name(1)="TYR"
    group_name(2)="NTYR"
    group_name(3)="CTYR"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="TYX".or.name=="NTYX".or.name=="CTYX") then
    group_name(1)="TYX"
    group_name(2)="NTYX"
    group_name(3)="CTYX"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="TRP".or.name=="NTRP".or.name=="CTRP") then
    group_name(1)="TRP"
    group_name(2)="NTRP"
    group_name(3)="CTRP"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="ARG".or.name=="NARG".or.name=="CARG") then
    group_name(1)="ARG"
    group_name(2)="NARG"
    group_name(3)="CARG"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="ARN".or.name=="NARN".or.name=="CARN") then
    group_name(1)="ARN"
    group_name(2)="NARN"
    group_name(3)="CARN"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="LYN".or.name=="NLYN".or.name=="CLYN") then
    group_name(1)="LYN"
    group_name(2)="NLYN"
    group_name(3)="CLYN"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="LYS".or.name=="NLYS".or.name=="CLYS") then
    group_name(1)="LYS"
    group_name(2)="NLYS"
    group_name(3)="CLYS"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="SER".or.name=="NSER".or.name=="CSER") then
    group_name(1)="SER"
    group_name(2)="NSER"
    group_name(3)="CSER"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="THR".or.name=="NTHR".or.name=="CTHR") then
    group_name(1)="THR"
    group_name(2)="NTHR"
    group_name(3)="CTHR"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="ASN".or.name=="NASN".or.name=="CASN") then
    group_name(1)="ASN"
    group_name(2)="NASN"
    group_name(3)="CASN"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="GLN".or.name=="NGLN".or.name=="CGLN") then
    group_name(1)="GLN"
    group_name(2)="NGLN"
    group_name(3)="CGLN"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="HIE".or.name=="NHIE".or.name=="CHIE") then
    group_name(1)="HIE"
    group_name(2)="NHIE"
    group_name(3)="CHIE"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="HIP".or.name=="NHIP".or.name=="CHIP") then
    group_name(1)="HIP"
    group_name(2)="NHIP"
    group_name(3)="CHIP"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="PRO".or.name=="NPRO".or.name=="CPRO") then
    group_name(1)="PRO"
    group_name(2)="NPRO"
    group_name(3)="CPRO"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="CYS".or.name=="NCYS".or.name=="CCYS") then
    group_name(1)="CYS"
    group_name(2)="NCYS"
    group_name(3)="CCYS"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="CYT".or.name=="NCYT".or.name=="CCYT") then
    group_name(1)="CYT"
    group_name(2)="NCYT"
    group_name(3)="CCYT"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="ALA".or.name=="NALA".or.name=="CALA") then
    group_name(1)="ALA"
    group_name(2)="NALA"
    group_name(3)="CALA"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="GLH".or.name=="NGLH".or.name=="CGLH") then
    group_name(1)="GLH"
    group_name(2)="NGLH"
    group_name(3)="CGLH"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="GLU".or.name=="NGLU".or.name=="CGLU") then
    group_name(1)="GLU"
    group_name(2)="NGLU"
    group_name(3)="CGLU"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="ASH".or.name=="NASH".or.name=="CASH") then
    group_name(1)="ASH"
    group_name(2)="NASH"
    group_name(3)="CASH"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
elseif(name=="ASP".or.name=="NASP".or.name=="CASP") then
    group_name(1)="ASP"
    group_name(2)="NASP"
    group_name(3)="CASP"
    do i=1, 3
        if(name==group_name(i)) then
        flag=i
        goto 5
        endif
    enddo
endif
5   continue

return
end subroutine groupinfo

subroutine scmf_choose_aminoacid(ip, aminoacid_number, aminoacid_name)
implicit none
integer                :: aminoacid_number
integer                :: ip, ip1, ip2, i
real                   :: ran2
character*4            :: char, aminoacid_name(10)

if(ph_value.le.3.9) then
    if(ip==1) then
        aminoacid_number=1
        aminoacid_name(1)="GLY"
    elseif(ip==2) then
        aminoacid_number=7
        aminoacid_name(1)="LEU"
        aminoacid_name(2)="VAL"
        aminoacid_name(3)="ILE"
        aminoacid_name(4)="MET"
        aminoacid_name(5)="PHE"
        aminoacid_name(6)="TYR"
        aminoacid_name(7)="TRP"
    elseif(ip==3) then
        aminoacid_number=3
        aminoacid_name(1)="ARG"
        aminoacid_name(2)="LYS"
        aminoacid_name(3)="HIP"
    elseif(ip==4) then
        aminoacid_number=6
        aminoacid_name(1)="SER"
        aminoacid_name(2)="THR"
        aminoacid_name(3)="ASN"
        aminoacid_name(4)="GLN"
        aminoacid_name(5)="GLH"
        aminoacid_name(6)="ASH"
    elseif(ip==5) then
        aminoacid_number=3
        aminoacid_name(1)="PRO"
        aminoacid_name(2)="ALA"
        aminoacid_name(3)="CYS"
    endif

elseif(ph_value.gt.3.9.and.ph_value.le.4.3) then
    if(ip==1) then
        aminoacid_number=1
        aminoacid_name(1)="GLY"
    elseif(ip==2) then
        aminoacid_number=7
        aminoacid_name(1)="LEU"
        aminoacid_name(2)="VAL"
        aminoacid_name(3)="ILE"
        aminoacid_name(4)="MET"
        aminoacid_name(5)="PHE"
        aminoacid_name(6)="TYR"
        aminoacid_name(7)="TRP"
    elseif(ip==3) then
        aminoacid_number=3
        aminoacid_name(1)="ARG"
        aminoacid_name(2)="LYS"
        aminoacid_name(3)="HIP"
    elseif(ip==4) then
        aminoacid_number=5
        aminoacid_name(1)="SER"
        aminoacid_name(2)="THR"
        aminoacid_name(3)="ASN"
        aminoacid_name(4)="GLN"
        aminoacid_name(5)="GLH"
    elseif(ip==5) then
        aminoacid_number=3
        aminoacid_name(1)="PRO"
        aminoacid_name(2)="ALA"
        aminoacid_name(3)="CYS"
    elseif(ip==6) then
        aminoacid_number=1
        aminoacid_name(1)="ASP"
    endif

elseif(ph_value.gt.4.3.and.ph_value.le.6.0) then
    if(ip==1) then
        aminoacid_number=1
        aminoacid_name(1)="GLY"
    elseif(ip==2) then
        aminoacid_number=7
        aminoacid_name(1)="LEU"
        aminoacid_name(2)="VAL"
        aminoacid_name(3)="ILE"
        aminoacid_name(4)="MET"
        aminoacid_name(5)="PHE"
        aminoacid_name(6)="TYR"
        aminoacid_name(7)="TRP"
    elseif(ip==3) then
        aminoacid_number=3
        aminoacid_name(1)="ARG"
        aminoacid_name(2)="LYS"
        aminoacid_name(3)="HIP"
    elseif(ip==4) then
        aminoacid_number=4
        aminoacid_name(1)="SER"
        aminoacid_name(2)="THR"
        aminoacid_name(3)="ASN"
        aminoacid_name(4)="GLN"
    elseif(ip==5) then
!         aminoacid_number=3
        aminoacid_number=2
        aminoacid_name(1)="PRO"
        aminoacid_name(2)="ALA"
!         aminoacid_name(3)="CYS"
    elseif(ip==6) then
        aminoacid_number=2
        aminoacid_name(1)="GLU"
        aminoacid_name(2)="ASP"
    endif

elseif(ph_value.gt.6.0.and.ph_value.lt.8.3) then
    aminoacid_number = num_AAs_available(ip)
    aminoacid_name = AAs_available(:,ip)
    
elseif(ph_value.ge.8.3.and.ph_value.lt.10.1) then
    if(ip==1) then
        aminoacid_number=1
        aminoacid_name(1)="GLY"
    elseif(ip==2) then
        aminoacid_number=7
        aminoacid_name(1)="LEU"
        aminoacid_name(2)="VAL"
        aminoacid_name(3)="ILE"
        aminoacid_name(4)="MET"
        aminoacid_name(5)="PHE"
        aminoacid_name(6)="TYR"
        aminoacid_name(7)="TRP"
    elseif(ip==3) then
        aminoacid_number=2
        aminoacid_name(1)="ARG"
        aminoacid_name(2)="LYS"
    elseif(ip==4) then
        aminoacid_number=5
        aminoacid_name(1)="SER"
        aminoacid_name(2)="THR"
        aminoacid_name(3)="ASN"
        aminoacid_name(4)="GLN"
        aminoacid_name(5)="HIE"
    elseif(ip==5) then
        aminoacid_number=2
        aminoacid_name(1)="PRO"
        aminoacid_name(2)="ALA"
    elseif(ip==6) then
        aminoacid_number=3
        aminoacid_name(1)="GLU"
        aminoacid_name(2)="ASP"
        aminoacid_name(3)="CYT"
    endif

elseif(ph_value.ge.10.1.and.ph_value.lt.10.5) then
    if(ip==1) then
        aminoacid_number=1
        aminoacid_name(1)="GLY"
    elseif(ip==2) then
        aminoacid_number=6
        aminoacid_name(1)="LEU"
        aminoacid_name(2)="VAL"
        aminoacid_name(3)="ILE"
        aminoacid_name(4)="MET"
        aminoacid_name(5)="PHE"
        aminoacid_name(6)="TRP"
    elseif(ip==3) then
        aminoacid_number=2
        aminoacid_name(1)="ARG"
        aminoacid_name(2)="LYS"
    elseif(ip==4) then
        aminoacid_number=5
        aminoacid_name(1)="SER"
        aminoacid_name(2)="THR"
        aminoacid_name(3)="ASN"
        aminoacid_name(4)="GLN"
        aminoacid_name(5)="HIE"
    elseif(ip==5) then
        aminoacid_number=2
        aminoacid_name(1)="PRO"
        aminoacid_name(2)="ALA"
    elseif(ip==6) then
        aminoacid_number=4
        aminoacid_name(1)="GLU"
        aminoacid_name(2)="ASP"
        aminoacid_name(3)="TYX"
        aminoacid_name(4)="CYT"
    endif

elseif(ph_value.ge.10.5.and.ph_value.lt.12.5) then
    if(ip==1) then
        aminoacid_number=1
        aminoacid_name(1)="GLY"
    elseif(ip==2) then
        aminoacid_number=6
        aminoacid_name(1)="LEU"
        aminoacid_name(2)="VAL"
        aminoacid_name(3)="ILE"
        aminoacid_name(4)="MET"
        aminoacid_name(5)="PHE"
        aminoacid_name(6)="TRP"
    elseif(ip==3) then
        aminoacid_number=1
        aminoacid_name(1)="ARG"
    elseif(ip==4) then
        aminoacid_number=6
        aminoacid_name(1)="SER"
        aminoacid_name(2)="THR"
        aminoacid_name(3)="ASN"
        aminoacid_name(4)="GLN"
        aminoacid_name(5)="HIE"
        aminoacid_name(6)="LYN"
    elseif(ip==5) then
        aminoacid_number=2
        aminoacid_name(1)="PRO"
        aminoacid_name(2)="ALA"
    elseif(ip==6) then
        aminoacid_number=4
        aminoacid_name(1)="GLU"
        aminoacid_name(2)="ASP"
        aminoacid_name(3)="TYX"
        aminoacid_name(4)="CYT"
    endif

elseif(ph_value.ge.12.5) then
    if(ip==1) then
        aminoacid_number=1
        aminoacid_name(1)="GLY"
    elseif(ip==2) then
        aminoacid_number=6
        aminoacid_name(1)="LEU"
        aminoacid_name(2)="VAL"
        aminoacid_name(3)="ILE"
        aminoacid_name(4)="MET"
        aminoacid_name(5)="PHE"
        aminoacid_name(6)="TRP"
    elseif(ip==4) then
        aminoacid_number=7
        aminoacid_name(1)="SER"
        aminoacid_name(2)="THR"
        aminoacid_name(3)="ASN"
        aminoacid_name(4)="GLN"
        aminoacid_name(5)="HIE"
        aminoacid_name(6)="LYN"
        aminoacid_name(7)="ARN"
    elseif(ip==5) then
        aminoacid_number=2
        aminoacid_name(1)="PRO"
        aminoacid_name(2)="ALA"
    elseif(ip==6) then
        aminoacid_number=4
        aminoacid_name(1)="GLU"
        aminoacid_name(2)="ASP"
        aminoacid_name(3)="TYX"
        aminoacid_name(4)="CYT"
    endif

endif

do i=1, (aminoacid_number-1)
    call RANDOM_NUMBER(ran2)
    ip1=int(ran2*aminoacid_number-1.0e-3)+1
    if(ip1.gt.aminoacid_number) ip1=aminoacid_number
    call RANDOM_NUMBER(ran2)
    ip2=int(ran2*aminoacid_number-1.0e-3)+1
    if(ip2.gt.aminoacid_number) ip2=aminoacid_number

    char=aminoacid_name(ip2)
    aminoacid_name(ip2)=aminoacid_name(ip1)
    aminoacid_name(ip1)=char
enddo

return
end subroutine scmf_choose_aminoacid

subroutine sidechain_category(ic, group, Iclass, num_sc_angles, sc_groups, index, monitor)
implicit none
integer                          :: num_sc_angles, sc_groups(6), monitor(6), i, ic
type(groupdetails)               :: group(gnum)
type(index4sidechain)            :: index(60)
type(conformer4sidechain)        :: Iclass(6)

!this subroutine gets the sc_groups (i.e. number of sidechain dihedral angles) and atom locations for sidechain
!each amino acid type has to be treated uniquely
sc_groups=0
if(group(ic)%gtype=="VAL".or.group(ic)%gtype=="NVAL".or.group(ic)%gtype=="CVAL") then
    num_sc_angles=1
    do i=1, group(ic)%cnum2
        if(group(ic)%atype2(i)=="CB") then
        sc_groups(1)=sc_groups(1)+1
        Iclass(1)%member(1:3,sc_groups(1))=group(ic)%coo2(1:3,i)
        index(i)%class_No=1; index(i)%member_No=sc_groups(1)
        monitor(1)=sc_groups(1)
        else
        sc_groups(2)=sc_groups(2)+1
        Iclass(2)%member(1:3,sc_groups(2))=group(ic)%coo2(1:3,i)
        index(i)%class_No=2; index(i)%member_No=sc_groups(2)
        endif
    enddo
elseif(group(ic)%gtype=="LEU".or.group(ic)%gtype=="NLEU".or.group(ic)%gtype=="CLEU") then
    num_sc_angles=2
    do i=1, group(ic)%cnum2
        if(group(ic)%atype2(i)=="CB") then
        sc_groups(1)=sc_groups(1)+1
        Iclass(1)%member(1:3,sc_groups(1))=group(ic)%coo2(1:3,i)
        index(i)%class_No=1; index(i)%member_No=sc_groups(1)
        monitor(1)=sc_groups(1)
        elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
            sc_groups(2)=sc_groups(2)+1
            Iclass(2)%member(1:3,sc_groups(2))=group(ic)%coo2(1:3,i)
            index(i)%class_No=2; index(i)%member_No=sc_groups(2)
            if(group(ic)%atype2(i)=="CG") monitor(2)=sc_groups(2)
            else
            sc_groups(3)=sc_groups(3)+1
            Iclass(3)%member(1:3,sc_groups(3))=group(ic)%coo2(1:3,i)
            index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        endif
    enddo
elseif(group(ic)%gtype=="ILE".or.group(ic)%gtype=="NILE".or.group(ic)%gtype=="CILE") then
    num_sc_angles=2
    do i=1, group(ic)%cnum2
        if(group(ic)%atype2(i)=="CB") then
        sc_groups(1)=sc_groups(1)+1
        Iclass(1)%member(1:3,sc_groups(1))=group(ic)%coo2(1:3,i)
        index(i)%class_No=1; index(i)%member_No=sc_groups(1)
        monitor(1)=sc_groups(1)
        elseif(group(ic)%atype2(i)=="HB".or.group(ic)%atype2(i)=="CG2".or.group(ic)%atype2(i)=="HG21" &
        .or.group(ic)%atype2(i)=="HG22".or.group(ic)%atype2(i)=="HG23".or.group(ic)%atype2(i)=="CG1") then
        sc_groups(2)=sc_groups(2)+1
        Iclass(2)%member(1:3,sc_groups(2))=group(ic)%coo2(1:3,i)
        index(i)%class_No=2; index(i)%member_No=sc_groups(2)
        if(group(ic)%atype2(i)=="CG1") monitor(2)=sc_groups(2)
        else
        sc_groups(3)=sc_groups(3)+1
        Iclass(3)%member(1:3,sc_groups(3))=group(ic)%coo2(1:3,i)
        index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        endif
    enddo
elseif(group(ic)%gtype=="PHE".or.group(ic)%gtype=="NPHE".or.group(ic)%gtype=="CPHE") then
    num_sc_angles=2
    do i=1, group(ic)%cnum2
        if(group(ic)%atype2(i)=="CB") then
        sc_groups(1)=sc_groups(1)+1
        Iclass(1)%member(1:3,sc_groups(1))=group(ic)%coo2(1:3,i)
        index(i)%class_No=1; index(i)%member_No=sc_groups(1)
        monitor(1)=sc_groups(1)
        elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
        sc_groups(2)=sc_groups(2)+1
        Iclass(2)%member(1:3,sc_groups(2))=group(ic)%coo2(1:3,i)
        index(i)%class_No=2; index(i)%member_No=sc_groups(2)
        if(group(ic)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        else
        sc_groups(3)=sc_groups(3)+1
        Iclass(3)%member(1:3,sc_groups(3))=group(ic)%coo2(1:3,i)
        index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        endif
    enddo
elseif(group(ic)%gtype=="TRP".or.group(ic)%gtype=="NTRP".or.group(ic)%gtype=="CTRP") then
    num_sc_angles=2
    do i=1, group(ic)%cnum2
        if(group(ic)%atype2(i)=="CB") then
        sc_groups(1)=sc_groups(1)+1
        Iclass(1)%member(1:3,sc_groups(1))=group(ic)%coo2(1:3,i)
        index(i)%class_No=1; index(i)%member_No=sc_groups(1)
        monitor(1)=sc_groups(1)
        elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
        sc_groups(2)=sc_groups(2)+1
        Iclass(2)%member(1:3,sc_groups(2))=group(ic)%coo2(1:3,i)
        index(i)%class_No=2; index(i)%member_No=sc_groups(2)
        if(group(ic)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        else
        sc_groups(3)=sc_groups(3)+1
        Iclass(3)%member(1:3,sc_groups(3))=group(ic)%coo2(1:3,i)
        index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        endif
    enddo
elseif(group(ic)%gtype=="TYR".or.group(ic)%gtype=="NTYR".or.group(ic)%gtype=="CTYR") then
    num_sc_angles=3
    do i=1, group(ic)%cnum2
        if(group(ic)%atype2(i)=="CB") then
        sc_groups(1)=sc_groups(1)+1
        Iclass(1)%member(1:3,sc_groups(1))=group(ic)%coo2(1:3,i)
        index(i)%class_No=1; index(i)%member_No=sc_groups(1)
        monitor(1)=sc_groups(1)
        elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
        sc_groups(2)=sc_groups(2)+1
        Iclass(2)%member(1:3,sc_groups(2))=group(ic)%coo2(1:3,i)
        index(i)%class_No=2; index(i)%member_No=sc_groups(2)
        if(group(ic)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        elseif(group(ic)%atype2(i)=="HH") then
        sc_groups(4)=sc_groups(4)+1
        Iclass(4)%member(1:3,sc_groups(4))=group(ic)%coo2(1:3,i)
        index(i)%class_No=4; index(i)%member_No=sc_groups(4)
        else
        sc_groups(3)=sc_groups(3)+1
        Iclass(3)%member(1:3,sc_groups(3))=group(ic)%coo2(1:3,i)
        index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        if(group(ic)%atype2(i)=="OH") monitor(3)=sc_groups(3)
        endif
    enddo
elseif(group(ic)%gtype=="TYX".or.group(ic)%gtype=="NTYX".or.group(ic)%gtype=="CTYX") then
    num_sc_angles=2
    do i=1, group(ic)%cnum2
        if(group(ic)%atype2(i)=="CB") then
        sc_groups(1)=sc_groups(1)+1
        Iclass(1)%member(1:3,sc_groups(1))=group(ic)%coo2(1:3,i)
        index(i)%class_No=1; index(i)%member_No=sc_groups(1)
        monitor(1)=sc_groups(1)
        elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
        sc_groups(2)=sc_groups(2)+1
        Iclass(2)%member(1:3,sc_groups(2))=group(ic)%coo2(1:3,i)
        index(i)%class_No=2; index(i)%member_No=sc_groups(2)
        if(group(ic)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        else
        sc_groups(3)=sc_groups(3)+1
        Iclass(3)%member(1:3,sc_groups(3))=group(ic)%coo2(1:3,i)
        index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        if(group(ic)%atype2(i)=="OH") monitor(3)=sc_groups(3)
        endif
    enddo
elseif(group(ic)%gtype=="SER".or.group(ic)%gtype=="NSER".or.group(ic)%gtype=="CSER") then
    num_sc_angles=2
    do i=1, group(ic)%cnum2
        if(group(ic)%atype2(i)=="CB") then
        sc_groups(1)=sc_groups(1)+1
        Iclass(1)%member(1:3,sc_groups(1))=group(ic)%coo2(1:3,i)
        index(i)%class_No=1; index(i)%member_No=sc_groups(1)
        monitor(1)=sc_groups(1)
        elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="OG") then
        sc_groups(2)=sc_groups(2)+1
        Iclass(2)%member(1:3,sc_groups(2))=group(ic)%coo2(1:3,i)
        index(i)%class_No=2; index(i)%member_No=sc_groups(2)
        if(group(ic)%atype2(i)=="OG") monitor(2)=sc_groups(2)
        else
        sc_groups(3)=sc_groups(3)+1
        Iclass(3)%member(1:3,sc_groups(3))=group(ic)%coo2(1:3,i)
        index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        endif
    enddo
elseif(group(ic)%gtype=="THR".or.group(ic)%gtype=="NTHR".or.group(ic)%gtype=="CTHR") then
    num_sc_angles=2
    do i=1, group(ic)%cnum2
        if(group(ic)%atype2(i)=="CB") then
        sc_groups(1)=sc_groups(1)+1
        Iclass(1)%member(1:3,sc_groups(1))=group(ic)%coo2(1:3,i)
        index(i)%class_No=1; index(i)%member_No=sc_groups(1)
        monitor(1)=sc_groups(1)
        elseif(group(ic)%atype2(i)=="HG1") then
        sc_groups(3)=sc_groups(3)+1
        Iclass(3)%member(1:3,sc_groups(3))=group(ic)%coo2(1:3,i)
        index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        else
        sc_groups(2)=sc_groups(2)+1
        Iclass(2)%member(1:3,sc_groups(2))=group(ic)%coo2(1:3,i)
        index(i)%class_No=2; index(i)%member_No=sc_groups(2)
        if(group(ic)%atype2(i)=="OG1") monitor(2)=sc_groups(2)
        endif
    enddo
elseif(group(ic)%gtype=="CYS".or.group(ic)%gtype=="NCYS".or.group(ic)%gtype=="CCYS") then
    num_sc_angles=2
    do i=1, group(ic)%cnum2
        if(group(ic)%atype2(i)=="CB") then
        sc_groups(1)=sc_groups(1)+1
        Iclass(1)%member(1:3,sc_groups(1))=group(ic)%coo2(1:3,i)
        index(i)%class_No=1; index(i)%member_No=sc_groups(1)
        monitor(1)=sc_groups(1)
        elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="SG") then
        sc_groups(2)=sc_groups(2)+1
        Iclass(2)%member(1:3,sc_groups(2))=group(ic)%coo2(1:3,i)
        index(i)%class_No=2; index(i)%member_No=sc_groups(2)
        if(group(ic)%atype2(i)=="SG") monitor(2)=sc_groups(2)
        else
        sc_groups(3)=sc_groups(3)+1
        Iclass(3)%member(1:3,sc_groups(3))=group(ic)%coo2(1:3,i)
        index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        endif
    enddo
elseif(group(ic)%gtype=="CYT".or.group(ic)%gtype=="NCYT".or.group(ic)%gtype=="CCYT") then
    num_sc_angles=1
    do i=1, group(ic)%cnum2
        if(group(ic)%atype2(i)=="CB") then
        sc_groups(1)=sc_groups(1)+1
        Iclass(1)%member(1:3,sc_groups(1))=group(ic)%coo2(1:3,i)
        index(i)%class_No=1; index(i)%member_No=sc_groups(1)
        monitor(1)=sc_groups(1)
        elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="SG") then
        sc_groups(2)=sc_groups(2)+1
        Iclass(2)%member(1:3,sc_groups(2))=group(ic)%coo2(1:3,i)
        index(i)%class_No=2; index(i)%member_No=sc_groups(2)
        if(group(ic)%atype2(i)=="SG") monitor(2)=sc_groups(2)
        endif
    enddo
elseif(group(ic)%gtype=="MET".or.group(ic)%gtype=="NMET".or.group(ic)%gtype=="CMET") then
    num_sc_angles=3
    do i=1, group(ic)%cnum2
        if(group(ic)%atype2(i)=="CB") then
        sc_groups(1)=sc_groups(1)+1
        Iclass(1)%member(1:3,sc_groups(1))=group(ic)%coo2(1:3,i)
        index(i)%class_No=1; index(i)%member_No=sc_groups(1)
        monitor(1)=sc_groups(1)
        elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
        sc_groups(2)=sc_groups(2)+1
        Iclass(2)%member(1:3,sc_groups(2))=group(ic)%coo2(1:3,i)
        index(i)%class_No=2; index(i)%member_No=sc_groups(2)
        if(group(ic)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        elseif(group(ic)%atype2(i)=="HG2".or.group(ic)%atype2(i)=="HG3".or.group(ic)%atype2(i)=="SD") then
        sc_groups(3)=sc_groups(3)+1
        Iclass(3)%member(1:3,sc_groups(3))=group(ic)%coo2(1:3,i)
        index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        if(group(ic)%atype2(i)=="SD") monitor(3)=sc_groups(3)
        else
        sc_groups(4)=sc_groups(4)+1
        Iclass(4)%member(1:3,sc_groups(4))=group(ic)%coo2(1:3,i)
        index(i)%class_No=4; index(i)%member_No=sc_groups(4)
        endif
    enddo
elseif(group(ic)%gtype=="ASN".or.group(ic)%gtype=="NASN".or.group(ic)%gtype=="CASN") then
    num_sc_angles=2
    do i=1, group(ic)%cnum2
        if(group(ic)%atype2(i)=="CB") then
        sc_groups(1)=sc_groups(1)+1
        Iclass(1)%member(1:3,sc_groups(1))=group(ic)%coo2(1:3,i)
        index(i)%class_No=1; index(i)%member_No=sc_groups(1)
        monitor(1)=sc_groups(1)
        elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
        sc_groups(2)=sc_groups(2)+1
        Iclass(2)%member(1:3,sc_groups(2))=group(ic)%coo2(1:3,i)
        index(i)%class_No=2; index(i)%member_No=sc_groups(2)
        if(group(ic)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        else
        sc_groups(3)=sc_groups(3)+1
        Iclass(3)%member(1:3,sc_groups(3))=group(ic)%coo2(1:3,i)
        index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        endif
    enddo
elseif(group(ic)%gtype=="GLN".or.group(ic)%gtype=="NGLN".or.group(ic)%gtype=="CGLN") then
    num_sc_angles=3
    do i=1, group(ic)%cnum2
        if(group(ic)%atype2(i)=="CB") then
        sc_groups(1)=sc_groups(1)+1
        Iclass(1)%member(1:3,sc_groups(1))=group(ic)%coo2(1:3,i)
        index(i)%class_No=1; index(i)%member_No=sc_groups(1)
        monitor(1)=sc_groups(1)
        elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
        sc_groups(2)=sc_groups(2)+1
        Iclass(2)%member(1:3,sc_groups(2))=group(ic)%coo2(1:3,i)
        index(i)%class_No=2; index(i)%member_No=sc_groups(2)
        if(group(ic)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        elseif(group(ic)%atype2(i)=="HG2".or.group(ic)%atype2(i)=="HG3".or.group(ic)%atype2(i)=="CD") then
        sc_groups(3)=sc_groups(3)+1
        Iclass(3)%member(1:3,sc_groups(3))=group(ic)%coo2(1:3,i)
        index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        if(group(ic)%atype2(i)=="CD") monitor(3)=sc_groups(3)
        else
        sc_groups(4)=sc_groups(4)+1
        Iclass(4)%member(1:3,sc_groups(4))=group(ic)%coo2(1:3,i)
        index(i)%class_No=4; index(i)%member_No=sc_groups(4)
        endif
    enddo
elseif(group(ic)%gtype=="ASP".or.group(ic)%gtype=="NASP".or.group(ic)%gtype=="CASP") then
    num_sc_angles=2
    do i=1, group(ic)%cnum2
        if(group(ic)%atype2(i)=="CB") then
        sc_groups(1)=sc_groups(1)+1
        Iclass(1)%member(1:3,sc_groups(1))=group(ic)%coo2(1:3,i)
        index(i)%class_No=1; index(i)%member_No=sc_groups(1)
        monitor(1)=sc_groups(1)
        elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
        sc_groups(2)=sc_groups(2)+1
        Iclass(2)%member(1:3,sc_groups(2))=group(ic)%coo2(1:3,i)
        index(i)%class_No=2; index(i)%member_No=sc_groups(2)
        if(group(ic)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        else
        sc_groups(3)=sc_groups(3)+1
        Iclass(3)%member(1:3,sc_groups(3))=group(ic)%coo2(1:3,i)
        index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        endif
    enddo
elseif(group(ic)%gtype=="ASH".or.group(ic)%gtype=="NASH".or.group(ic)%gtype=="CASH") then
    num_sc_angles=2
    do i=1, group(ic)%cnum2
        if(group(ic)%atype2(i)=="CB") then
        sc_groups(1)=sc_groups(1)+1
        Iclass(1)%member(1:3,sc_groups(1))=group(ic)%coo2(1:3,i)
        index(i)%class_No=1; index(i)%member_No=sc_groups(1)
        monitor(1)=sc_groups(1)
        elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
        sc_groups(2)=sc_groups(2)+1
        Iclass(2)%member(1:3,sc_groups(2))=group(ic)%coo2(1:3,i)
        index(i)%class_No=2; index(i)%member_No=sc_groups(2)
        if(group(ic)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        else
        sc_groups(3)=sc_groups(3)+1
        Iclass(3)%member(1:3,sc_groups(3))=group(ic)%coo2(1:3,i)
        index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        endif
    enddo
elseif(group(ic)%gtype=="GLU".or.group(ic)%gtype=="NGLU".or.group(ic)%gtype=="CGLU") then
    num_sc_angles=3
    do i=1, group(ic)%cnum2
        if(group(ic)%atype2(i)=="CB") then
        sc_groups(1)=sc_groups(1)+1
        Iclass(1)%member(1:3,sc_groups(1))=group(ic)%coo2(1:3,i)
        index(i)%class_No=1; index(i)%member_No=sc_groups(1)
        monitor(1)=sc_groups(1)
        elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
        sc_groups(2)=sc_groups(2)+1
        Iclass(2)%member(1:3,sc_groups(2))=group(ic)%coo2(1:3,i)
        index(i)%class_No=2; index(i)%member_No=sc_groups(2)
        if(group(ic)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        elseif(group(ic)%atype2(i)=="HG2".or.group(ic)%atype2(i)=="HG3".or.group(ic)%atype2(i)=="CD") then
        sc_groups(3)=sc_groups(3)+1
        Iclass(3)%member(1:3,sc_groups(3))=group(ic)%coo2(1:3,i)
        index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        if(group(ic)%atype2(i)=="CD") monitor(3)=sc_groups(3)
        else
        sc_groups(4)=sc_groups(4)+1
        Iclass(4)%member(1:3,sc_groups(4))=group(ic)%coo2(1:3,i)
        index(i)%class_No=4; index(i)%member_No=sc_groups(4)
        endif
    enddo
elseif(group(ic)%gtype=="GLH".or.group(ic)%gtype=="NGLH".or.group(ic)%gtype=="CGLH") then
    num_sc_angles=3
    do i=1, group(ic)%cnum2
        if(group(ic)%atype2(i)=="CB") then
        sc_groups(1)=sc_groups(1)+1
        Iclass(1)%member(1:3,sc_groups(1))=group(ic)%coo2(1:3,i)
        index(i)%class_No=1; index(i)%member_No=sc_groups(1)
        monitor(1)=sc_groups(1)
        elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
        sc_groups(2)=sc_groups(2)+1
        Iclass(2)%member(1:3,sc_groups(2))=group(ic)%coo2(1:3,i)
        index(i)%class_No=2; index(i)%member_No=sc_groups(2)
        if(group(ic)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        elseif(group(ic)%atype2(i)=="HG2".or.group(ic)%atype2(i)=="HG3".or.group(ic)%atype2(i)=="CD") then
        sc_groups(3)=sc_groups(3)+1
        Iclass(3)%member(1:3,sc_groups(3))=group(ic)%coo2(1:3,i)
        index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        if(group(ic)%atype2(i)=="CD") monitor(3)=sc_groups(3)
        else
        sc_groups(4)=sc_groups(4)+1
        Iclass(4)%member(1:3,sc_groups(4))=group(ic)%coo2(1:3,i)
        index(i)%class_No=4; index(i)%member_No=sc_groups(4)
        endif
    enddo
elseif(group(ic)%gtype=="HIE".or.group(ic)%gtype=="NHIE".or.group(ic)%gtype=="CHIE") then
    num_sc_angles=2
    do i=1, group(ic)%cnum2
        if(group(ic)%atype2(i)=="CB") then
        sc_groups(1)=sc_groups(1)+1
        Iclass(1)%member(1:3,sc_groups(1))=group(ic)%coo2(1:3,i)
        index(i)%class_No=1; index(i)%member_No=sc_groups(1)
        monitor(1)=sc_groups(1)
        elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
        sc_groups(2)=sc_groups(2)+1
        Iclass(2)%member(1:3,sc_groups(2))=group(ic)%coo2(1:3,i)
        index(i)%class_No=2; index(i)%member_No=sc_groups(2)
        if(group(ic)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        else
        sc_groups(3)=sc_groups(3)+1
        Iclass(3)%member(1:3,sc_groups(3))=group(ic)%coo2(1:3,i)
        index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        endif
    enddo
elseif(group(ic)%gtype=="HIP".or.group(ic)%gtype=="NHIP".or.group(ic)%gtype=="CHIP") then
    num_sc_angles=2
    do i=1, group(ic)%cnum2
        if(group(ic)%atype2(i)=="CB") then
        sc_groups(1)=sc_groups(1)+1
        Iclass(1)%member(1:3,sc_groups(1))=group(ic)%coo2(1:3,i)
        index(i)%class_No=1; index(i)%member_No=sc_groups(1)
        monitor(1)=sc_groups(1)
        elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
        sc_groups(2)=sc_groups(2)+1
        Iclass(2)%member(1:3,sc_groups(2))=group(ic)%coo2(1:3,i)
        index(i)%class_No=2; index(i)%member_No=sc_groups(2)
        if(group(ic)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        else
        sc_groups(3)=sc_groups(3)+1
        Iclass(3)%member(1:3,sc_groups(3))=group(ic)%coo2(1:3,i)
        index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        endif
    enddo
elseif(group(ic)%gtype=="LYS".or.group(ic)%gtype=="NLYS".or.group(ic)%gtype=="CLYS") then
    num_sc_angles=4
    do i=1, group(ic)%cnum2
        if(group(ic)%atype2(i)=="CB") then
        sc_groups(1)=sc_groups(1)+1
        Iclass(1)%member(1:3,sc_groups(1))=group(ic)%coo2(1:3,i)
        index(i)%class_No=1; index(i)%member_No=sc_groups(1)
        monitor(1)=sc_groups(1)
        elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
        sc_groups(2)=sc_groups(2)+1
        Iclass(2)%member(1:3,sc_groups(2))=group(ic)%coo2(1:3,i)
        index(i)%class_No=2; index(i)%member_No=sc_groups(2)
        if(group(ic)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        elseif(group(ic)%atype2(i)=="HG2".or.group(ic)%atype2(i)=="HG3".or.group(ic)%atype2(i)=="CD") then
        sc_groups(3)=sc_groups(3)+1
        Iclass(3)%member(1:3,sc_groups(3))=group(ic)%coo2(1:3,i)
        index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        if(group(ic)%atype2(i)=="CD") monitor(3)=sc_groups(3)
        elseif(group(ic)%atype2(i)=="HD2".or.group(ic)%atype2(i)=="HD3".or.group(ic)%atype2(i)=="CE") then
        sc_groups(4)=sc_groups(4)+1
        Iclass(4)%member(1:3,sc_groups(4))=group(ic)%coo2(1:3,i)
        index(i)%class_No=4; index(i)%member_No=sc_groups(4)
        if(group(ic)%atype2(i)=="CE") monitor(4)=sc_groups(4)
        else
        sc_groups(5)=sc_groups(5)+1
        Iclass(5)%member(1:3,sc_groups(5))=group(ic)%coo2(1:3,i)
        index(i)%class_No=5; index(i)%member_No=sc_groups(5)
        endif
    enddo
elseif(group(ic)%gtype=="LYN".or.group(ic)%gtype=="NLYN".or.group(ic)%gtype=="CLYN") then
    num_sc_angles=4
    do i=1, group(ic)%cnum2
        if(group(ic)%atype2(i)=="CB") then
        sc_groups(1)=sc_groups(1)+1
        Iclass(1)%member(1:3,sc_groups(1))=group(ic)%coo2(1:3,i)
        index(i)%class_No=1; index(i)%member_No=sc_groups(1)
        monitor(1)=sc_groups(1)
        elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
        sc_groups(2)=sc_groups(2)+1
        Iclass(2)%member(1:3,sc_groups(2))=group(ic)%coo2(1:3,i)
        index(i)%class_No=2; index(i)%member_No=sc_groups(2)
        if(group(ic)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        elseif(group(ic)%atype2(i)=="HG2".or.group(ic)%atype2(i)=="HG3".or.group(ic)%atype2(i)=="CD") then
        sc_groups(3)=sc_groups(3)+1
        Iclass(3)%member(1:3,sc_groups(3))=group(ic)%coo2(1:3,i)
        index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        if(group(ic)%atype2(i)=="CD") monitor(3)=sc_groups(3)
        elseif(group(ic)%atype2(i)=="HD2".or.group(ic)%atype2(i)=="HD3".or.group(ic)%atype2(i)=="CE") then
        sc_groups(4)=sc_groups(4)+1
        Iclass(4)%member(1:3,sc_groups(4))=group(ic)%coo2(1:3,i)
        index(i)%class_No=4; index(i)%member_No=sc_groups(4)
        if(group(ic)%atype2(i)=="CE") monitor(4)=sc_groups(4)
        else
        sc_groups(5)=sc_groups(5)+1
        Iclass(5)%member(1:3,sc_groups(5))=group(ic)%coo2(1:3,i)
        index(i)%class_No=5; index(i)%member_No=sc_groups(5)
        endif
    enddo
elseif(group(ic)%gtype=="ARG".or.group(ic)%gtype=="NARG".or.group(ic)%gtype=="CARG") then
    num_sc_angles=4
    do i=1, group(ic)%cnum2
        if(group(ic)%atype2(i)=="CB") then
        sc_groups(1)=sc_groups(1)+1
        Iclass(1)%member(1:3,sc_groups(1))=group(ic)%coo2(1:3,i)
        index(i)%class_No=1; index(i)%member_No=sc_groups(1)
        monitor(1)=sc_groups(1)
        elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
        sc_groups(2)=sc_groups(2)+1
        Iclass(2)%member(1:3,sc_groups(2))=group(ic)%coo2(1:3,i)
        index(i)%class_No=2; index(i)%member_No=sc_groups(2)
        if(group(ic)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        elseif(group(ic)%atype2(i)=="HG2".or.group(ic)%atype2(i)=="HG3".or.group(ic)%atype2(i)=="CD") then
        sc_groups(3)=sc_groups(3)+1
        Iclass(3)%member(1:3,sc_groups(3))=group(ic)%coo2(1:3,i)
        index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        if(group(ic)%atype2(i)=="CD") monitor(3)=sc_groups(3)
        elseif(group(ic)%atype2(i)=="HD2".or.group(ic)%atype2(i)=="HD3".or.group(ic)%atype2(i)=="NE") then
        sc_groups(4)=sc_groups(4)+1
        Iclass(4)%member(1:3,sc_groups(4))=group(ic)%coo2(1:3,i)
        index(i)%class_No=4; index(i)%member_No=sc_groups(4)
        if(group(ic)%atype2(i)=="NE") monitor(4)=sc_groups(4)
        else
        sc_groups(5)=sc_groups(5)+1
        Iclass(5)%member(1:3,sc_groups(5))=group(ic)%coo2(1:3,i)
        index(i)%class_No=5; index(i)%member_No=sc_groups(5)
        endif
    enddo
elseif(group(ic)%gtype=="ARN".or.group(ic)%gtype=="NARN".or.group(ic)%gtype=="CARN") then
    num_sc_angles=4
    do i=1, group(ic)%cnum2
        if(group(ic)%atype2(i)=="CB") then
        sc_groups(1)=sc_groups(1)+1
        Iclass(1)%member(1:3,sc_groups(1))=group(ic)%coo2(1:3,i)
        index(i)%class_No=1; index(i)%member_No=sc_groups(1)
        monitor(1)=sc_groups(1)
        elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
        sc_groups(2)=sc_groups(2)+1
        Iclass(2)%member(1:3,sc_groups(2))=group(ic)%coo2(1:3,i)
        index(i)%class_No=2; index(i)%member_No=sc_groups(2)
        if(group(ic)%atype2(i)=="CG") monitor(2)=sc_groups(2)
        elseif(group(ic)%atype2(i)=="HG2".or.group(ic)%atype2(i)=="HG3".or.group(ic)%atype2(i)=="CD") then
        sc_groups(3)=sc_groups(3)+1
        Iclass(3)%member(1:3,sc_groups(3))=group(ic)%coo2(1:3,i)
        index(i)%class_No=3; index(i)%member_No=sc_groups(3)
        if(group(ic)%atype2(i)=="CD") monitor(3)=sc_groups(3)
        elseif(group(ic)%atype2(i)=="HD2".or.group(ic)%atype2(i)=="HD3".or.group(ic)%atype2(i)=="NE") then
        sc_groups(4)=sc_groups(4)+1
        Iclass(4)%member(1:3,sc_groups(4))=group(ic)%coo2(1:3,i)
        index(i)%class_No=4; index(i)%member_No=sc_groups(4)
        if(group(ic)%atype2(i)=="NE") monitor(4)=sc_groups(4)
        else
        sc_groups(5)=sc_groups(5)+1
        Iclass(5)%member(1:3,sc_groups(5))=group(ic)%coo2(1:3,i)
        index(i)%class_No=5; index(i)%member_No=sc_groups(5)
        endif
    enddo
endif

return
end subroutine sidechain_category

subroutine dihedralangle_reading(gtype, dihedral_num_options, dihedral)
!this subroutine gets the atom indices defining each dihedral, and associate dihedral properties (perioidicity, spring constant, and equlibrium value)
implicit none
integer                         :: dihedral_num_options, i, j
character*4                     :: gtype
type(dihedralparameters)        :: dihedral

!this funciton
open(10, file=trim(lib_path)//'/DihedralAngle/'//trim(gtype), status="old")
    read(10, "(i8)") dihedral_num_options
    do i=1, dihedral_num_options
        read(10,"(5i8)") dihedral%iph(i), dihedral%jph(i), dihedral%kph(i), dihedral%lph(i), dihedral%multiply(i)
        do j=1, dihedral%multiply(i)
            read(10,"(3e16.8)") dihedral%pk(j,i), dihedral%pn(j,i), dihedral%phase(j,i)
        enddo
    enddo
close(10)

return
end subroutine dihedralangle_reading

end module database