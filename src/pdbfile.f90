module pdbfile

use datatypes
use sys_vars
contains

subroutine calc_system_size()
implicit none
integer                      :: res_num, status, iunit, anum, flag, pep_atom, atom_tracker
real                         :: x, y, z
character*4                  :: char, atom_name, name

status = 0
res_num = 0
anum = 0
flag = 0
open(newunit = iunit, file=trim(system_pdb))
do while(status.ge.0)
    if (status.eq.0) then
        gnum = res_num; atom_tracker=anum
    endif
    if (flag.eq.0 .and. gnum.gt.pep_res) then
        flag = 1
        pep_atom = atom_tracker-1
    endif
    read(iunit, *, iostat=status) char, anum, atom_name, name, res_num, x, y, z
enddo

atom_num_pep = pep_res * 30
rec_res = gnum - pep_res
atom_num_rec = atom_tracker - pep_atom
atom_num = atom_num_pep + atom_num_rec
rewind(iunit)
close(iunit)
return
end subroutine calc_system_size

subroutine calc_pep_rec_size()
implicit none
integer                      :: anum, status, iunit, res_num, atom_tracker
real                         :: x, y, z
character*4                  :: char, atom_name, name

!first get the system size for the peptide
status = 0
anum = 0
res_num = 0
open(newunit = iunit, file=trim(pep_pdb))
do while(status.eq.0)
    gnum = res_num
    read(iunit, *, iostat=status) char, anum, atom_name, name, res_num, x, y, z
enddo
pep_res = gnum
atom_num_pep = pep_res * 30
rewind(iunit) !be sure to rewind the file, else the pdb file won't be read later to get the system coordinates
close(iunit)

!now get the receptor size
status = 0
res_num = 0
open(newunit = iunit, file=trim(rec_pdb))
do while(status.eq.0)
    gnum = res_num; ; atom_tracker = anum
    read(iunit, *, iostat=status) char, anum, atom_name, name, res_num, x, y, z
enddo
rec_res = gnum
atom_num_rec = atom_tracker
rewind(iunit) !be sure to rewind the file, else the pdb file won't be read later to get the system coordinates
close(iunit)

!specify total system size
gnum = rec_res + pep_res
atom_num = atom_num_pep + atom_num_rec
return
end subroutine calc_pep_rec_size

subroutine load_peptide(group_pep)
implicit none
integer                         :: anum, status=1, flag1, res_num
real                            :: x, y, z, net_mass_sidechain
character(4)                    :: char, atom_name, name
type(groupdetails)              :: group_pep(pep_res)

group_pep%cnum1=0
group_pep%cnum2=0
group_pep%cnum3=0
flag1=0
num_atoms = 0 
net_mass_backbone = 0 
net_mass_sidechain = 0

!tell the program that all coordinates and parameters need to be reloaded into the array format when doing energy calculations
res_params_flags = 1
res_coords_flags = 1
atom_links_flag = 1

open(10, file=trim(pep_pdb))
do while(.true.)
    read(10, *, iostat=status) char, anum, atom_name, name, res_num, x, y, z
    if(status.lt.0) exit
    if(name=="ALA".or.name=="ARG".or.name=="ASN".or.name=="ASP".or.name=="CYS".or.name=="GLN" &
    .or.name=="GLU".or.name=="GLY".or.name=="HIE".or.name=="ILE".or.name=="LEU".or.name=="LYS" &
    .or.name=="MET".or.name=="PHE".or.name=="PRO".or.name=="SER".or.name=="THR".or.name=="TRP" &
    .or.name=="TYR".or.name=="VAL") then
        if (res_num.eq.pep_res) then
            name = "C" // name
        elseif (res_num.eq.1) then
            name = "N" // name
        endif
        group_pep(res_num)%gtype=name
        if(atom_name=="N".or.atom_name=="H".or.atom_name=="H1".or.atom_name=="H2".or.atom_name=="H3".or.atom_name=="CA" &
            .or.atom_name=="HA".or.atom_name=="HA2".or.atom_name=="HA3") then
            group_pep(res_num)%cnum1=group_pep(res_num)%cnum1+1
            group_pep(res_num)%atype1(group_pep(res_num)%cnum1)=atom_name
            group_pep(res_num)%coo1(1,group_pep(res_num)%cnum1)=x
            group_pep(res_num)%coo1(2,group_pep(res_num)%cnum1)=y
            group_pep(res_num)%coo1(3,group_pep(res_num)%cnum1)=z
            num_atoms = num_atoms + 1
            if (atom_name(1:1).eq."H") then
                group_pep(res_num)%masses1(group_pep(res_num)%cnum1) = 1.0
                net_mass_backbone = net_mass_backbone + 1.0
            elseif (atom_name(1:1).eq."N") then
                group_pep(res_num)%masses1(group_pep(res_num)%cnum1) = 14.0
                net_mass_backbone = net_mass_backbone + 14.0
            elseif(atom_name(1:1).eq."C") then
                group_pep(res_num)%masses1(group_pep(res_num)%cnum1) = 12.0
                net_mass_backbone = net_mass_backbone + 12.0
            endif
        elseif(atom_name=="C".or.atom_name=="O".or.atom_name=="OXT") then
            group_pep(res_num)%cnum3=group_pep(res_num)%cnum3+1
            group_pep(res_num)%atype3(group_pep(res_num)%cnum3)=atom_name
            group_pep(res_num)%coo3(1,group_pep(res_num)%cnum3)=x
            group_pep(res_num)%coo3(2,group_pep(res_num)%cnum3)=y
            group_pep(res_num)%coo3(3,group_pep(res_num)%cnum3)=z
            num_atoms = num_atoms + 1
            if (atom_name(1:1).eq."O") then
                group_pep(res_num)%masses3(group_pep(res_num)%cnum3) = 16.0
                net_mass_backbone = net_mass_backbone + 16.0
            elseif(atom_name(1:1).eq."C") then
                group_pep(res_num)%masses3(group_pep(res_num)%cnum3) = 12.0
                net_mass_backbone = net_mass_backbone + 12.0
            endif
        else
            group_pep(res_num)%cnum2=group_pep(res_num)%cnum2+1
            group_pep(res_num)%atype2(group_pep(res_num)%cnum2)=atom_name
            group_pep(res_num)%coo2(1,group_pep(res_num)%cnum2)=x
            group_pep(res_num)%coo2(2,group_pep(res_num)%cnum2)=y
            group_pep(res_num)%coo2(3,group_pep(res_num)%cnum2)=z
            num_atoms = num_atoms + 1
            if (atom_name(1:1).eq."H") then
                group_pep(res_num)%masses2(group_pep(res_num)%cnum2) = 1.0
                net_mass_sidechain = net_mass_sidechain + 1.0
            elseif (atom_name(1:1).eq."N") then
                group_pep(res_num)%masses2(group_pep(res_num)%cnum2) = 14.0
                net_mass_sidechain = net_mass_sidechain + 14.0
            elseif(atom_name(1:1).eq."C") then
                group_pep(res_num)%masses2(group_pep(res_num)%cnum2) = 12.0
                net_mass_sidechain = net_mass_sidechain + 12.0
            elseif(atom_name(1:1).eq."O") then
                group_pep(res_num)%masses2(group_pep(res_num)%cnum2) = 16.0
                net_mass_sidechain = net_mass_sidechain + 16.0
            elseif(atom_name(1:1).eq."S") then
                group_pep(res_num)%masses2(group_pep(res_num)%cnum2) = 32.0
                net_mass_sidechain = net_mass_sidechain + 32.0
            endif
        endif
    endif
enddo
net_mass_pep = net_mass_backbone + net_mass_sidechain

endsubroutine load_peptide

subroutine load_peptide_from_ensemble(group_pep, file)
implicit none
integer                         :: anum, status=1, flag1, res_num
real                            :: x, y, z, net_mass_sidechain
character(4)                    :: char, atom_name, name
character*50                    :: file
character*100                   :: full_file
type(groupdetails)              :: group_pep(pep_res)

group_pep%cnum1=0
group_pep%cnum2=0
group_pep%cnum3=0
flag1=0
num_atoms = 0 
net_mass_backbone = 0 
net_mass_sidechain = 0

full_file = trim(backbones_path)//trim(file)
!tell the program that all coordinates and parameters need to be reloaded into the array format when doing energy calculations
res_params_flags = 1
res_coords_flags = 1
atom_links_flag = 1

open(10, file=trim(full_file))
do while(.true.)
    read(10, *, iostat=status) char, anum, atom_name, name, res_num, x, y, z
    if(status.lt.0) exit
    if(name=="ALA".or.name=="ARG".or.name=="ASN".or.name=="ASP".or.name=="CYS".or.name=="GLN" &
    .or.name=="GLU".or.name=="GLY".or.name=="HIE".or.name=="ILE".or.name=="LEU".or.name=="LYS" &
    .or.name=="MET".or.name=="PHE".or.name=="PRO".or.name=="SER".or.name=="THR".or.name=="TRP" &
    .or.name=="TYR".or.name=="VAL") then
        if (res_num.eq.pep_res) then
            name = "C" // name
        elseif (res_num.eq.1) then
            name = "N" // name
        endif
        group_pep(res_num)%gtype=name
        if(atom_name=="N".or.atom_name=="H".or.atom_name=="H1".or.atom_name=="H2".or.atom_name=="H3".or.atom_name=="CA" &
            .or.atom_name=="HA".or.atom_name=="HA2".or.atom_name=="HA3") then
            group_pep(res_num)%cnum1=group_pep(res_num)%cnum1+1
            group_pep(res_num)%atype1(group_pep(res_num)%cnum1)=atom_name
            group_pep(res_num)%coo1(1,group_pep(res_num)%cnum1)=x
            group_pep(res_num)%coo1(2,group_pep(res_num)%cnum1)=y
            group_pep(res_num)%coo1(3,group_pep(res_num)%cnum1)=z
            num_atoms = num_atoms + 1
            if (atom_name(1:1).eq."H") then
                group_pep(res_num)%masses1(group_pep(res_num)%cnum1) = 1.0
                net_mass_backbone = net_mass_backbone + 1.0
            elseif (atom_name(1:1).eq."N") then
                group_pep(res_num)%masses1(group_pep(res_num)%cnum1) = 14.0
                net_mass_backbone = net_mass_backbone + 14.0
            elseif(atom_name(1:1).eq."C") then
                group_pep(res_num)%masses1(group_pep(res_num)%cnum1) = 12.0
                net_mass_backbone = net_mass_backbone + 12.0
            endif
        elseif(atom_name=="C".or.atom_name=="O".or.atom_name=="OXT") then
            group_pep(res_num)%cnum3=group_pep(res_num)%cnum3+1
            group_pep(res_num)%atype3(group_pep(res_num)%cnum3)=atom_name
            group_pep(res_num)%coo3(1,group_pep(res_num)%cnum3)=x
            group_pep(res_num)%coo3(2,group_pep(res_num)%cnum3)=y
            group_pep(res_num)%coo3(3,group_pep(res_num)%cnum3)=z
            num_atoms = num_atoms + 1
            if (atom_name(1:1).eq."O") then
                group_pep(res_num)%masses3(group_pep(res_num)%cnum3) = 16.0
                net_mass_backbone = net_mass_backbone + 16.0
            elseif(atom_name(1:1).eq."C") then
                group_pep(res_num)%masses3(group_pep(res_num)%cnum3) = 12.0
                net_mass_backbone = net_mass_backbone + 12.0
            endif
        else
            group_pep(res_num)%cnum2=group_pep(res_num)%cnum2+1
            group_pep(res_num)%atype2(group_pep(res_num)%cnum2)=atom_name
            group_pep(res_num)%coo2(1,group_pep(res_num)%cnum2)=x
            group_pep(res_num)%coo2(2,group_pep(res_num)%cnum2)=y
            group_pep(res_num)%coo2(3,group_pep(res_num)%cnum2)=z
            num_atoms = num_atoms + 1
            if (atom_name(1:1).eq."H") then
                group_pep(res_num)%masses2(group_pep(res_num)%cnum2) = 1.0
                net_mass_sidechain = net_mass_sidechain + 1.0
            elseif (atom_name(1:1).eq."N") then
                group_pep(res_num)%masses2(group_pep(res_num)%cnum2) = 14.0
                net_mass_sidechain = net_mass_sidechain + 14.0
            elseif(atom_name(1:1).eq."C") then
                group_pep(res_num)%masses2(group_pep(res_num)%cnum2) = 12.0
                net_mass_sidechain = net_mass_sidechain + 12.0
            elseif(atom_name(1:1).eq."O") then
                group_pep(res_num)%masses2(group_pep(res_num)%cnum2) = 16.0
                net_mass_sidechain = net_mass_sidechain + 16.0
            elseif(atom_name(1:1).eq."S") then
                group_pep(res_num)%masses2(group_pep(res_num)%cnum2) = 32.0
                net_mass_sidechain = net_mass_sidechain + 32.0
            endif
        endif
    endif
enddo
net_mass_pep = net_mass_backbone + net_mass_sidechain

endsubroutine load_peptide_from_ensemble

subroutine load_receptor(group_rec)

implicit none
integer                                 :: anum, status=1, flag1, res_num
real                                    :: x, y, z
character(4)                            :: char, atom_name, res_name
type(groupdetails)                      :: group_rec(:) 

group_rec%cnum1=0
group_rec%cnum2=0
group_rec%cnum3=0
flag1=0
num_atoms = 0 
net_mass_rec = 0 
open(10, file=trim(rec_pdb))
do while(.true.)
    read(10, *, iostat=status) char, anum, atom_name, res_name, res_num, x, y, z
    if(status.lt.0) exit
    group_rec(res_num)%gtype=res_name
    group_rec(res_num)%cnum2=group_rec(res_num)%cnum2+1
    group_rec(res_num)%atype2(group_rec(res_num)%cnum2)=atom_name
    group_rec(res_num)%coo2(1,group_rec(res_num)%cnum2)=x
    group_rec(res_num)%coo2(2,group_rec(res_num)%cnum2)=y
    group_rec(res_num)%coo2(3,group_rec(res_num)%cnum2)=z
    if (atom_name(1:1).eq."H") then
        group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 1.0
        net_mass_rec = net_mass_rec + 1.0
    elseif (atom_name(1:1).eq."N") then
        group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 14.0
        net_mass_rec = net_mass_rec + 14.0
    elseif(atom_name(1:1).eq."C") then
        group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 12.0
        net_mass_rec = net_mass_rec + 12.0
    elseif(atom_name(1:1).eq."O") then
        group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 16.0
        net_mass_rec = net_mass_rec + 16.0
    endif
enddo

end subroutine load_receptor

subroutine load_system(group_pep, group_rec, original_group, groupdata_backup)
! Read the pdbfile
implicit none
integer                      :: anum, status=1, flag, flag1, i, res_num
real                         :: x, y, z, net_mass_sidechain
character*4                  :: char, atom_name, name
type(groupdetails)           :: group_pep(pep_res), group_rec(rec_res), original_group(pep_res)
type(databackup)             :: groupdata_backup(pep_res)

group_pep%cnum1=0
group_pep%cnum2=0
group_pep%cnum3=0
group_rec%cnum1=0
group_rec%cnum2=0
group_rec%cnum3=0
fragmentnum=0
flag=0
flag1=0
net_mass_sidechain=0
net_mass_backbone=0

!tell the program that all coordinates and parameters need to be reloaded into the array format when doing energy calculations
res_params_flags = 1
res_coords_flags = 1
atom_links_flag = 1

flag1=0
open(10, file=trim(system_pdb))
do while(.true.)
    read(10, *, iostat=status) char, anum, atom_name, name, res_num, x, y, z
    if(status.lt.0) exit
    if (status.eq.0) then ! skip over lines that do not fit the desired format (this will be an issue if the pdb file has a chain or extra stuff after the coordinates, should fix eventually)
        if(res_num.le.pep_res) then
            if(name=="ALA".or.name=="ARG".or.name=="ASN".or.name=="ASP".or.name=="CYS".or.name=="GLN" &
                .or.name=="GLU".or.name=="GLY".or.name=="HIE".or.name=="ILE".or.name=="LEU".or.name=="LYS" &
                .or.name=="MET".or.name=="PHE".or.name=="PRO".or.name=="SER".or.name=="THR".or.name=="TRP" &
                .or.name=="TYR".or.name=="VAL") then
                if (res_num.eq.pep_res) then
                    name = "C" // name
                elseif (res_num.eq.1) then
                    name = "N" // name
                endif
                group_pep(res_num)%gtype=name
                if(atom_name=="N".or.atom_name=="H".or.atom_name=="H1".or.atom_name=="H2".or.atom_name=="H3".or.atom_name=="CA" &
                    .or.atom_name=="HA".or.atom_name=="HA2".or.atom_name=="HA3") then
                    group_pep(res_num)%cnum1=group_pep(res_num)%cnum1+1
                    group_pep(res_num)%atype1(group_pep(res_num)%cnum1)=atom_name
                    group_pep(res_num)%coo1(1,group_pep(res_num)%cnum1)=x
                    group_pep(res_num)%coo1(2,group_pep(res_num)%cnum1)=y
                    group_pep(res_num)%coo1(3,group_pep(res_num)%cnum1)=z
                    if(name=="PRO".or.name=="NPRO".or.name=="CPRO") then
                        if(restart_flag==0) flag1=1
                    else
                        if(atom_name=="H".or.atom_name=="H1") then
                            groupdata_backup(res_num)%coo(1)=x
                            groupdata_backup(res_num)%coo(2)=y
                            groupdata_backup(res_num)%coo(3)=z
                        endif
                    endif
                    if (atom_name(1:1).eq."H") then
                        group_pep(res_num)%masses1(group_pep(res_num)%cnum1) = 1.0
                        net_mass_backbone = net_mass_backbone + 1.0
                    elseif (atom_name(1:1).eq."N") then
                        group_pep(res_num)%masses1(group_pep(res_num)%cnum1) = 14.0
                        net_mass_backbone = net_mass_backbone + 14.0
                    elseif(atom_name(1:1).eq."C") then
                        group_pep(res_num)%masses1(group_pep(res_num)%cnum1) = 12.0
                        net_mass_backbone = net_mass_backbone + 12.0
                    endif
                elseif(atom_name=="C".or.atom_name=="O".or.atom_name=="OXT") then
                    group_pep(res_num)%cnum3=group_pep(res_num)%cnum3+1
                    group_pep(res_num)%atype3(group_pep(res_num)%cnum3)=atom_name
                    group_pep(res_num)%coo3(1,group_pep(res_num)%cnum3)=x
                    group_pep(res_num)%coo3(2,group_pep(res_num)%cnum3)=y
                    group_pep(res_num)%coo3(3,group_pep(res_num)%cnum3)=z
                    if (atom_name(1:1).eq."O") then
                        group_pep(res_num)%masses3(group_pep(res_num)%cnum3) = 16.0
                        net_mass_backbone = net_mass_backbone + 16.0
                    elseif(atom_name(1:1).eq."C") then
                        group_pep(res_num)%masses3(group_pep(res_num)%cnum3) = 12.0
                        net_mass_backbone = net_mass_backbone + 12.0
                    endif
                else
                    group_pep(res_num)%cnum2=group_pep(res_num)%cnum2+1
                    group_pep(res_num)%atype2(group_pep(res_num)%cnum2)=atom_name
                    group_pep(res_num)%coo2(1,group_pep(res_num)%cnum2)=x
                    group_pep(res_num)%coo2(2,group_pep(res_num)%cnum2)=y
                    group_pep(res_num)%coo2(3,group_pep(res_num)%cnum2)=z
                
                    if (atom_name(1:1).eq."H") then
                        group_pep(res_num)%masses2(group_pep(res_num)%cnum2) = 1.0
                        net_mass_sidechain = net_mass_sidechain + 1.0
                    elseif (atom_name(1:1).eq."N") then
                        group_pep(res_num)%masses2(group_pep(res_num)%cnum2) = 14.0
                        net_mass_sidechain = net_mass_sidechain + 14.0
                    elseif(atom_name(1:1).eq."C") then
                        group_pep(res_num)%masses2(group_pep(res_num)%cnum2) = 12.0
                        net_mass_sidechain = net_mass_sidechain + 12.0
                    elseif(atom_name(1:1).eq."O") then
                        group_pep(res_num)%masses2(group_pep(res_num)%cnum2) = 16.0
                        net_mass_sidechain = net_mass_sidechain + 16.0
                    elseif(atom_name(1:1).eq."S") then
                        group_pep(res_num)%masses2(group_pep(res_num)%cnum2) = 32.0
                        net_mass_sidechain = net_mass_sidechain + 32.0
                    endif
                endif

                if(res_num.ne.flag) then
                    flag=res_num
                    do i=1, helix_num
                        if(res_num.ge.helix_start(i).and.res_num.le.helix_end(i)) goto 10
                    enddo
                    fragmentnum=fragmentnum+1
                    residuesite(fragmentnum)=res_num
    10                    continue
                endif
            else
                write(*,"(A)") "Unknown atom in the input pdb file!"
                write(*,"(A,A,A,A,A,I0)") "Issue is atom name: ", atom_name, " in residue: ", name, ", atom number: ", anum
                write(*,"(A)") "Please double check the pdb file atom names"
                stop
            endif

        else
            res_num=res_num-pep_res !shift residue number so it starts at 1 for beginning of receptor
            if(receptor_name=="rna") then
                if(name=="RC".or.name=="RA".or.name=="RU".or.name=="RG".or.name=="STA".or.name=="SMU".or.name=="1MA".or.name=="PSU"  &
                        .or.name=="5MU".or.name=="RC5".or.name=="RA5".or.name=="RU5".or.name=="RG5".or.name=="RC3".or.name=="RA3" &
                    .or.name=="RU3".or.name=="RG3".or.name=="2SU".or.name=="6TA") then
                    group_rec(res_num)%gtype=name
                    if(atom_name=="HO5'".or.atom_name=="H5T".or.atom_name=="P".or.atom_name=="O1P".or.atom_name=="OP1".or.atom_name=="O2P".or.atom_name=="OP2" &
                    .or.atom_name=="O5'".or.atom_name=="C5'".or.atom_name=="1H5'".or.atom_name=="H5'".or.atom_name=="2H5'".or.atom_name=="H5''".or. &
                    atom_name=="C4'".or.atom_name=="H4'".or.atom_name=="O4'".or.atom_name=="C1'".or.atom_name=="H1'".or.atom_name=="H5'1".or.atom_name=="H5'2") then
                        if(atom_name=="H5'1") then
                            group_rec(res_num)%cnum1=group_rec(res_num)%cnum1+1
                            group_rec(res_num)%atype1(group_rec(res_num)%cnum1)="1H5'"
                            group_rec(res_num)%coo1(1,group_rec(res_num)%cnum1)=x
                            group_rec(res_num)%coo1(2,group_rec(res_num)%cnum1)=y
                            group_rec(res_num)%coo1(3,group_rec(res_num)%cnum1)=z
                        elseif(atom_name=="H5'2") then
                            group_rec(res_num)%cnum1=group_rec(res_num)%cnum1+1
                            group_rec(res_num)%atype1(group_rec(res_num)%cnum1)="2H5'"
                            group_rec(res_num)%coo1(1,group_rec(res_num)%cnum1)=x
                            group_rec(res_num)%coo1(2,group_rec(res_num)%cnum1)=y
                            group_rec(res_num)%coo1(3,group_rec(res_num)%cnum1)=z
                        else
                            group_rec(res_num)%cnum1=group_rec(res_num)%cnum1+1
                            group_rec(res_num)%atype1(group_rec(res_num)%cnum1)=atom_name
                            group_rec(res_num)%coo1(1,group_rec(res_num)%cnum1)=x
                            group_rec(res_num)%coo1(2,group_rec(res_num)%cnum1)=y
                            group_rec(res_num)%coo1(3,group_rec(res_num)%cnum1)=z
                        endif
                        if (atom_name(1:1).eq."H") then
                            group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 1.0
                            net_mass_rec = net_mass_rec + 1.0
                        elseif (atom_name(1:1).eq."N") then
                            group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 14.0
                            net_mass_rec = net_mass_rec + 14.0
                        elseif(atom_name(1:1).eq."C") then
                            group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 12.0
                            net_mass_rec = net_mass_rec + 12.0
                        elseif(atom_name(1:1).eq."O") then
                            group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 16.0
                            net_mass_rec = net_mass_rec + 16.0
                        elseif(atom_name(1:1).eq."S") then
                            group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 32.0
                            net_mass_rec = net_mass_rec + 32.0
                        endif
                    elseif(atom_name=="C3'".or.atom_name=="H3'".or.atom_name=="C2'".or.atom_name=="1H2'".or.atom_name=="H2'".or.atom_name=="O2'" &
                        .or.atom_name=="2HO'".or.atom_name=="HO2'".or.atom_name=="H3T".or.atom_name=="HO3'".or.atom_name=="O3'".or.atom_name=="H2'1".or.atom_name=="HO'2") then
                        if(atom_name=="H2'1") then
                            group_rec(res_num)%cnum3=group_rec(res_num)%cnum3+1
                            group_rec(res_num)%atype3(group_rec(res_num)%cnum3)="1H2'"
                            group_rec(res_num)%coo3(1,group_rec(res_num)%cnum3)=x
                            group_rec(res_num)%coo3(2,group_rec(res_num)%cnum3)=y
                            group_rec(res_num)%coo3(3,group_rec(res_num)%cnum3)=z
                        elseif(atom_name=="HO'2") then
                            group_rec(res_num)%cnum3=group_rec(res_num)%cnum3+1
                            group_rec(res_num)%atype3(group_rec(res_num)%cnum3)="2HO'"
                            group_rec(res_num)%coo3(1,group_rec(res_num)%cnum3)=x
                            group_rec(res_num)%coo3(2,group_rec(res_num)%cnum3)=y
                            group_rec(res_num)%coo3(3,group_rec(res_num)%cnum3)=z
                        else
                            group_rec(res_num)%cnum3=group_rec(res_num)%cnum3+1
                            group_rec(res_num)%atype3(group_rec(res_num)%cnum3)=atom_name
                            group_rec(res_num)%coo3(1,group_rec(res_num)%cnum3)=x
                            group_rec(res_num)%coo3(2,group_rec(res_num)%cnum3)=y
                            group_rec(res_num)%coo3(3,group_rec(res_num)%cnum3)=z
                            if (atom_name(1:1).eq."H") then
                                group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 1.0
                                net_mass_rec = net_mass_rec + 1.0
                            elseif (atom_name(1:1).eq."N") then
                                group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 14.0
                                net_mass_rec = net_mass_rec + 14.0
                            elseif(atom_name(1:1).eq."C") then
                                group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 12.0
                                net_mass_rec = net_mass_rec + 12.0
                            elseif(atom_name(1:1).eq."O") then
                                group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 16.0
                                net_mass_rec = net_mass_rec + 16.0
                            elseif(atom_name(1:1).eq."S") then
                                group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 32.0
                                net_mass_rec = net_mass_rec + 32.0
                            endif
                        endif
                        
                    else
                        group_rec(res_num)%cnum2=group_rec(res_num)%cnum2+1
                        group_rec(res_num)%atype2(group_rec(res_num)%cnum2)=atom_name
                        group_rec(res_num)%coo2(1,group_rec(res_num)%cnum2)=x
                        group_rec(res_num)%coo2(2,group_rec(res_num)%cnum2)=y
                        group_rec(res_num)%coo2(3,group_rec(res_num)%cnum2)=z

                        if (atom_name(1:1).eq."H") then
                            group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 1.0
                            net_mass_rec = net_mass_rec + 1.0
                        elseif (atom_name(1:1).eq."N") then
                            group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 14.0
                            net_mass_rec = net_mass_rec + 14.0
                        elseif(atom_name(1:1).eq."C") then
                            group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 12.0
                            net_mass_rec = net_mass_rec + 12.0
                        elseif(atom_name(1:1).eq."O") then
                            group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 16.0
                            net_mass_rec = net_mass_rec + 16.0
                        elseif(atom_name(1:1).eq."S") then
                            group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 32.0
                            net_mass_rec = net_mass_rec + 32.0
                        endif
                    endif
                else
                    write(*,"(A)") "Unknown atom in the input pdb file!"
                    write(*,"(A,A,A,A,A,I0)") "Issue is atom name: ", atom_name, " in residue: ", name, ", atom number: ", anum
                    write(*,"(A)") "Please double check the pdb file atom names"
                    stop
                endif

            elseif(receptor_name=="protein") then
                if(name=="ALA".or.name=="ARG".or.name=="ASN".or.name=="ASP".or.name=="CYS".or.name=="GLN" &
                    .or.name=="GLU".or.name=="GLY".or.name=="HIE".or.name=="ILE".or.name=="LEU".or.name=="LYS" &
                    .or.name=="MET".or.name=="PHE".or.name=="PRO".or.name=="SER".or.name=="THR".or.name=="TRP" &
                    .or.name=="TYR".or.name=="VAL".or.name=="NALA".or.name=="NARG".or.name=="NASN".or.name=="NASP" &
                    .or.name=="NCYS".or.name=="NGLN".or.name=="NGLU".or.name=="NGLY".or.name=="NHIE".or.name=="NILE" &
                    .or.name=="NLEU".or.name=="NLYS".or.name=="NMET".or.name=="NPHE".or.name=="NPRO".or.name=="NSER" &
                    .or.name=="NTHR".or.name=="NTRP".or.name=="NTYR".or.name=="NVAL".or.name=="CALA".or.name=="CARG" &
                    .or.name=="CASN".or.name=="CASP".or.name=="CCYS".or.name=="CGLN".or.name=="CGLU".or.name=="CGLY" &
                    .or.name=="CHIE".or.name=="CILE".or.name=="CLEU".or.name=="CLYS".or.name=="CMET".or.name=="CPHE" &
                    .or.name=="CPRO".or.name=="CSER".or.name=="CTHR".or.name=="CTRP".or.name=="CTYR".or.name=="CVAL" &
                    .or.name=="TYX".or.name=="ARN".or.name=="HIP".or.name=="CYT".or.name=="LYN".or.name=="GLH".or.name=="ASH" &
                    .or.name=="NTYX".or.name=="NARN".or.name=="NHIP".or.name=="NCYT".or.name=="NLYN".or.name=="NGLH".or.name=="NASH" &
                    .or.name=="CTYX".or.name=="CARN".or.name=="CHIP".or.name=="CCYT".or.name=="CLYN".or.name=="CGLH".or.name=="CASH" &
                    .or.name=="CYX") then
                    group_rec(res_num)%gtype=name
                    if(atom_name=="N".or.atom_name=="H".or.atom_name=="H1".or.atom_name=="H2".or.atom_name=="H3".or.atom_name=="CA" &
                        .or.atom_name=="HA".or.atom_name=="HA2".or.atom_name=="HA3") then
                        group_rec(res_num)%cnum1=group_rec(res_num)%cnum1+1
                        group_rec(res_num)%atype1(group_rec(res_num)%cnum1)=atom_name
                        group_rec(res_num)%coo1(1,group_rec(res_num)%cnum1)=x
                        group_rec(res_num)%coo1(2,group_rec(res_num)%cnum1)=y
                        group_rec(res_num)%coo1(3,group_rec(res_num)%cnum1)=z
                        if (atom_name(1:1).eq."H") then
                            group_rec(res_num)%masses1(group_rec(res_num)%cnum1) = 1.0
                            net_mass_rec = net_mass_rec + 1.0
                        elseif (atom_name(1:1).eq."N") then
                            group_rec(res_num)%masses1(group_rec(res_num)%cnum1) = 14.0
                            net_mass_rec = net_mass_rec + 14.0
                        elseif(atom_name(1:1).eq."C") then
                            group_rec(res_num)%masses1(group_rec(res_num)%cnum1) = 12.0
                            net_mass_rec = net_mass_rec + 12.0
                        endif
                    elseif(atom_name=="C".or.atom_name=="O".or.atom_name=="OXT") then
                        group_rec(res_num)%cnum3=group_rec(res_num)%cnum3+1
                        group_rec(res_num)%atype3(group_rec(res_num)%cnum3)=atom_name
                        group_rec(res_num)%coo3(1,group_rec(res_num)%cnum3)=x
                        group_rec(res_num)%coo3(2,group_rec(res_num)%cnum3)=y
                        group_rec(res_num)%coo3(3,group_rec(res_num)%cnum3)=z
                        if(atom_name(1:1).eq."C") then
                            group_rec(res_num)%masses3(group_rec(res_num)%cnum3) = 12.0
                            net_mass_rec = net_mass_rec + 12.0
                        elseif(atom_name(1:1).eq."O") then
                            group_rec(res_num)%masses3(group_rec(res_num)%cnum3) = 16.0
                            net_mass_rec = net_mass_rec + 16.0
                        endif
                    else
                        group_rec(res_num)%cnum2=group_rec(res_num)%cnum2+1
                        group_rec(res_num)%atype2(group_rec(res_num)%cnum2)=atom_name
                        group_rec(res_num)%coo2(1,group_rec(res_num)%cnum2)=x
                        group_rec(res_num)%coo2(2,group_rec(res_num)%cnum2)=y
                        group_rec(res_num)%coo2(3,group_rec(res_num)%cnum2)=z
                        if (atom_name(1:1).eq."H") then
                            group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 1.0
                            net_mass_rec = net_mass_rec + 1.0
                        elseif (atom_name(1:1).eq."N") then
                            group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 14.0
                            net_mass_rec = net_mass_rec + 14.0
                        elseif(atom_name(1:1).eq."C") then
                            group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 12.0
                            net_mass_rec = net_mass_rec + 12.0
                        elseif(atom_name(1:1).eq."O") then
                            group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 16.0
                            net_mass_rec = net_mass_rec + 16.0
                        elseif(atom_name(1:1).eq."S") then
                            group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 32.0
                            net_mass_rec = net_mass_rec + 32.0
                        endif
                    endif
                else
                    write(*,"(A)") "Unknown atom in the input pdb file!"
                    write(*,"(A,A,A,A,A,I0)") "Issue is atom name: ", atom_name, " in residue: ", name, ", atom number: ", anum
                    write(*,"(A)") "Please double check the pdb file atom names"
                    stop
                endif

            elseif(receptor_name=="molecule") then
                group_rec(res_num)%gtype=name
                group_rec(res_num)%cnum2=group_rec(res_num)%cnum2+1
                group_rec(res_num)%atype2(group_rec(res_num)%cnum2)=atom_name
                group_rec(res_num)%coo2(1,group_rec(res_num)%cnum2)=x
                group_rec(res_num)%coo2(2,group_rec(res_num)%cnum2)=y
                group_rec(res_num)%coo2(3,group_rec(res_num)%cnum2)=z
                if (atom_name(1:1).eq."H") then
                    group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 1.0
                    net_mass_rec = net_mass_rec + 1.0
                elseif (atom_name(1:1).eq."N") then
                    group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 14.0
                    net_mass_rec = net_mass_rec + 14.0
                elseif(atom_name(1:1).eq."C") then
                    group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 12.0
                    net_mass_rec = net_mass_rec + 12.0
                elseif(atom_name(1:1).eq."O") then
                    group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 16.0
                    net_mass_rec = net_mass_rec + 16.0
                elseif(atom_name(1:1).eq."S") then
                    group_rec(res_num)%masses2(group_rec(res_num)%cnum2) = 32.0
                    net_mass_rec = net_mass_rec + 32.0
                endif
            endif
        endif
    endif
enddo
close(10)

!if(fragmentnum.lt.3) then
!    open(20, file="error.txt", access="append")
!        write(*,*) "The code works only for the situation there are at leat 3 amino &
!        acids as the pivots to move the backbone conformations!"
!    close(20)
!    stop
!endif

if(restart_flag==0) then 
    original_group=group_pep  ! original group used for calculating RMSD
    if(flag1==1) then
        write(*,"(A)") "WARNING: There are Pro residues on the original peptide chain!"
        write(*,"(A)") "WARNING: A supplemental _backup4backbone.txt_ file is required for the program to execute properly"
        open(20, file="backup4backbone.txt", status="old")
            do while(.true.)
                read(20, *, iostat=status) char, anum, atom_name, name, res_num, x, y, z
                if(status.ne.0) exit
                groupdata_backup(res_num)%coo(1)=x
                groupdata_backup(res_num)%coo(2)=y
                groupdata_backup(res_num)%coo(3)=z
            enddo
        close(20)
    endif
else !We are restarting design in middle, so we need to also load the original peptide's structure to do RMSD calculations
    original_group%cnum1=0
    original_group%cnum2=0
    original_group%cnum3=0

    open(10, file=restart_pdb)
    do while(.true.)
        read(10, *, iostat=status) char, anum, atom_name, name, res_num, x, y, z
        if(status.ne.0) exit
        if(res_num.le.pep_res) then
            if(name=="ALA".or.name=="ARG".or.name=="ASN".or.name=="ASP".or.name=="CYS".or.name=="GLN" &
                .or.name=="GLU".or.name=="GLY".or.name=="HIE".or.name=="ILE".or.name=="LEU".or.name=="LYS" &
                    .or.name=="MET".or.name=="PHE".or.name=="PRO".or.name=="SER".or.name=="THR".or.name=="TRP" &
                .or.name=="TYR".or.name=="VAL") then
                if (res_num.eq.1) then
                    name = "N" // name
                elseif (res_num.eq.pep_res) then
                    name = "C" // name
                endif
                original_group(res_num)%gtype=name
                if(atom_name=="N".or.atom_name=="H".or.atom_name=="H1".or.atom_name=="H2".or.atom_name=="H3".or.atom_name=="CA" &
                    .or.atom_name=="HA".or.atom_name=="HA2".or.atom_name=="HA3") then
                    original_group(res_num)%cnum1=original_group(res_num)%cnum1+1
                    original_group(res_num)%atype1(original_group(res_num)%cnum1)=atom_name
                    original_group(res_num)%coo1(1,original_group(res_num)%cnum1)=x
                    original_group(res_num)%coo1(2,original_group(res_num)%cnum1)=y
                    original_group(res_num)%coo1(3,original_group(res_num)%cnum1)=z
                elseif(atom_name=="C".or.atom_name=="O".or.atom_name=="OXT") then
                    original_group(res_num)%cnum3=original_group(res_num)%cnum3+1
                    original_group(res_num)%atype3(original_group(res_num)%cnum3)=atom_name
                    original_group(res_num)%coo3(1,original_group(res_num)%cnum3)=x
                    original_group(res_num)%coo3(2,original_group(res_num)%cnum3)=y
                    original_group(res_num)%coo3(3,original_group(res_num)%cnum3)=z
                else
                    original_group(res_num)%cnum2=original_group(res_num)%cnum2+1
                    original_group(res_num)%atype2(original_group(res_num)%cnum2)=atom_name
                    original_group(res_num)%coo2(1,original_group(res_num)%cnum2)=x
                    original_group(res_num)%coo2(2,original_group(res_num)%cnum2)=y
                    original_group(res_num)%coo2(3,original_group(res_num)%cnum2)=z
                endif
            else
                write(*,"(A)") "Unknown atom in the input pdb file!"
                write(*,"(A,A,A,A,A,I0)") "Issue is atom name: ", atom_name, " in residue: ", name, ", atom number: ", anum
                write(*,"(A)") "Please double check the pdb file atom names"
                stop
            endif
        endif
    enddo
    close(10)
    open(5, file="backup4backbone.txt", status="old")
        do i=1, pep_res
            read(5,*) groupdata_backup(i)%coo(1), groupdata_backup(i)%coo(2), groupdata_backup(i)%coo(3)
        enddo
        read(5,*) flag4conformer
        read(5,*) energy_min
    close(5)
endif

net_mass_pep = net_mass_backbone + net_mass_sidechain

return
end subroutine load_system

subroutine gen_pdb_name(cycle, step, attempt, name)
implicit none
integer             :: cycle, step, attempt
character*5         :: stepchar, attemptchar, cyclechar
character*50        :: name

write(cyclechar,"(i3)") cycle
write(stepchar,"(i5)") step
write(attemptchar,"(i4)") attempt
name = trim(adjustl(cyclechar))//'_'//trim(adjustl(stepchar))//'_'//trim(adjustl(attemptchar))
return

end subroutine gen_pdb_name

subroutine writepdb_pep(group_pep, name)
implicit none
integer                 :: i, j, atomnum
character*3             :: res_type
character*50            :: name
type(groupdetails)      :: group_pep(pep_res)

open(10,file=trim(pdb_path)//trim(name)//'.pdb') 

atomnum=1
do i=1, pep_res

    !trim residue name for N- and C-termini
    if (i.eq.1.or.i.eq.pep_res) then
        res_type = group_pep(i)%gtype(2:)
    else
        res_type = group_pep(i)%gtype
    endif

    do j=1, group_pep(i)%cnum1
        write(10,1) "ATOM  ", atomnum, " ", group_pep(i)%atype1(j), " ", res_type, " A", i, "    ",&
        group_pep(i)%coo1(1,j), group_pep(i)%coo1(2,j), group_pep(i)%coo1(3,j)
        atomnum=atomnum+1
    enddo
    do j=1, group_pep(i)%cnum2
        write(10,1) "ATOM  ", atomnum, " ", group_pep(i)%atype2(j), " ", res_type, " A", i, "    ",&
        group_pep(i)%coo2(1,j), group_pep(i)%coo2(2,j), group_pep(i)%coo2(3,j)
        atomnum=atomnum+1
    enddo
    do j=1, group_pep(i)%cnum3
        write(10,1) "ATOM  ", atomnum, " ", group_pep(i)%atype3(j), " ", res_type, " A", i, "    ",&
        group_pep(i)%coo3(1,j), group_pep(i)%coo3(2,j), group_pep(i)%coo3(3,j)
        atomnum=atomnum+1
    enddo

enddo
close(10)
1 format(a6, i5, a1, a4, a1, a3, a2, i4, a4, 3f8.3)

end subroutine writepdb_pep

subroutine writepdb_sys(group_pep, group_rec, name)
implicit none
integer                 :: i, j, atomnum
character*3             :: res_type
character*50            :: name
type(groupdetails)      :: group_pep(pep_res), group_rec(rec_res)

open(10,file=trim(pdb_path)//trim(name)//'.pdb') 

atomnum=1
do i=1, pep_res
    !trim residue name for N- and C-termini
    if (i.eq.1.or.i.eq.pep_res) then
        res_type = group_pep(i)%gtype(2:)
    else
        res_type = group_pep(i)%gtype
    endif

    do j=1, group_pep(i)%cnum1
        write(10,1) "ATOM  ", atomnum, " ", group_pep(i)%atype1(j), " ", res_type, " A", i, "    ",&
        group_pep(i)%coo1(1,j), group_pep(i)%coo1(2,j), group_pep(i)%coo1(3,j)
        atomnum=atomnum+1
    enddo
    do j=1, group_pep(i)%cnum2
        write(10,1) "ATOM  ", atomnum, " ", group_pep(i)%atype2(j), " ", res_type, " A", i, "    ",&
        group_pep(i)%coo2(1,j), group_pep(i)%coo2(2,j), group_pep(i)%coo2(3,j)
        atomnum=atomnum+1
    enddo
    do j=1, group_pep(i)%cnum3
        write(10,1) "ATOM  ", atomnum, " ", group_pep(i)%atype3(j), " ", res_type, " A", i, "    ",&
        group_pep(i)%coo3(1,j), group_pep(i)%coo3(2,j), group_pep(i)%coo3(3,j)
        atomnum=atomnum+1
    enddo
enddo

!separate peptide from receptor
write(10,"(a3)") "TER"

!now write receptor coordinates
do i=1, rec_res
    res_type = group_rec(i)%gtype
    do j=1, group_rec(i)%cnum1
        write(10,1) "ATOM  ", atomnum, " ", group_rec(i)%atype1(j), " ", res_type, " B", i, "    ",&
        group_rec(i)%coo1(1,j), group_rec(i)%coo1(2,j), group_rec(i)%coo1(3,j)
        atomnum=atomnum+1
    enddo
    do j=1, group_rec(i)%cnum2
        write(10,1) "ATOM  ", atomnum, " ", group_rec(i)%atype2(j), " ", res_type, " B", i, "    ",&
        group_rec(i)%coo2(1,j), group_rec(i)%coo2(2,j), group_rec(i)%coo2(3,j)
        atomnum=atomnum+1
    enddo
    do j=1, group_rec(i)%cnum3
        write(10,1) "ATOM  ", atomnum, " ", group_rec(i)%atype3(j), " ", res_type, " B", i, "    ",&
        group_rec(i)%coo3(1,j), group_rec(i)%coo3(2,j), group_rec(i)%coo3(3,j)
        atomnum=atomnum+1
    enddo
enddo
close(10)

1 format(a6, i5, a1, a4, a1, a3, a2, i4, a4, 3f8.3)

return
end subroutine writepdb_sys

end module pdbfile
