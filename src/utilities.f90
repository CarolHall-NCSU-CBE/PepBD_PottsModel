module utilities
use datatypes
use sys_vars
use math
use omp_lib

contains

subroutine translate_group(group, vector, num_res)
type(groupdetails)            :: group(:)
real                          :: vector(3)
integer                       :: i,j, num_res

if (num_res.eq.pep_res) then
    do i=1,pep_res
        do j=1, group(i)%cnum1
            group(i)%coo1(:,j) = group(i)%coo1(:,j) + vector
        enddo
    
        do j=1, group(i)%cnum2
            group(i)%coo2(:,j) = group(i)%coo2(:,j) + vector
        enddo
    
        do j=1, group(i)%cnum3
            group(i)%coo3(:,j) = group(i)%coo3(:,j) + vector
        enddo
    enddo
else
    do i=1,num_res
        do j=1, group(i)%cnum1
            group(i)%coo1(:,j) = group(i)%coo1(:,j) + vector
        enddo

        do j=1, group(i)%cnum2
            group(i)%coo2(:,j) = group(i)%coo2(:,j) + vector
        enddo

        do j=1, group(i)%cnum3
            group(i)%coo3(:,j) = group(i)%coo3(:,j) + vector
        enddo
    enddo
endif

end subroutine translate_group

subroutine rotate_peptide(group_pep, m, com)
implicit none
integer                     :: i, j
real                        :: m(3,3), com(3)
type(groupdetails)          :: group_pep(pep_res)

!rotate the peptide
do i=1, pep_res
    do j=1,group_pep(i)%cnum1
        group_pep(i)%coo1(:,j) = matmul(m,group_pep(i)%coo1(:,j)-com)+com
    enddo
    do j=1,group_pep(i)%cnum2
        group_pep(i)%coo2(:,j) = matmul(m,group_pep(i)%coo2(:,j)-com)+com
    enddo
    do j=1,group_pep(i)%cnum3
        group_pep(i)%coo3(:,j) = matmul(m,group_pep(i)%coo3(:,j)-com)+com
    enddo
enddo

res_coords_flags = 1
end subroutine rotate_peptide

subroutine rotate_translate_peptide(group_pep, m, vec, com)
implicit none
real                        :: m(3,3), vec(3), com(3)
type(groupdetails)          :: group_pep(:)

call rotate_peptide(group_pep, m, com)
call translate_group(group_pep, vec, pep_res)

end subroutine rotate_translate_peptide

subroutine random_peptide_translation(group_pep)
implicit none
type(groupdetails)      :: group_pep(pep_res)
real                    :: translate_vec(3), val1, val2, val3
REAL                    :: lims(3,2), ran, mag, com(3)
integer                 :: i, sel_index, scalar

!*******NOTE**********!
!in future, may want to change this function to translation occurs only in one direction at at time, similar to rotation function

!generate random numbers between -1 and 1
call RANDOM_NUMBER(val1)
val1 = val1*2 - 1
call RANDOM_NUMBER(val2)
val2 = val2*2 - 1
call RANDOM_NUMBER(val3)
val3 = val3*2 - 1

!determine current center of mass of peptide (current value may not be accurate since center of mass not calculated after backbone or sequence changes)
call calc_com(group_pep, com, pep_res)
!determine max translation in each direction allowed while still being within cutoff of starting center of mass position
!z-direction is special for surface design, so we have extra control over translation distance in z direction

lims(1,:)= (/com(1) - start_pep_com(1) + max_dev_from_center, start_pep_com(1) - com(1) + max_dev_from_center/)
lims(2,:)= (/com(2) - start_pep_com(2) + max_dev_from_center, start_pep_com(2) - com(2) + max_dev_from_center/)
!if (translate_z_flag.eq.1) then
!    lims(3,:)= (/(com(3)-surface_z) - min_zdist_from_surface, max_zdist_from_surface - (com(3)-surface_z)/)
!else 
lims(3,:)= (/com(3) - start_pep_com(3) + max_dev_from_center, start_pep_com(3) - com(3) + max_dev_from_center/)
!endif

translate_vec = 0; scalar=0
do i =1,3
    sel_index = 0
    !pick directions of translations, only moving in sladirections where minimum translation distance is permissible
    if (lims(i,1).gt.translate_min .and. lims(1,2).gt.translate_min) then
        !either direction works, pick one randomly
        call RANDOM_NUMBER(ran)
        if (ran.lt.0.5) then !go in negative direction
            sel_index = 1
            scalar = -1
        else !go in positive direction
            sel_index = 2
            scalar = 1
        endif
    elseif(lims(i,1).gt.translate_min) then
        sel_index = 1
        scalar = -1
    elseif(lims(i,2).gt.translate_min) then
        sel_index = 2
        scalar = 1
    endif

    !now calculate translation vector (if possible to move in this direction).
    !This translation is done such that a random distance is picked in negative/postive direction (selected above) between
    !translate_min and the minimum of translate_max and limit of peptide displacement from starting center of mass
    if (sel_index.ne.0) then
        mag = min(translate_max, lims(i,sel_index))
        call RANDOM_NUMBER(ran) 
        translate_vec(i) = translate_min + scalar*ran*(mag-translate_min)
    endif
enddo

!translate peptide using translation vector calculated in above loop
call translate_group(group_pep, translate_vec, pep_res)

end subroutine random_peptide_translation

subroutine random_peptide_rotation(group_pep)
implicit none
type(groupdetails)          :: group_pep(pep_res)
real                        :: ran, theta, m(3,3)
integer                     :: dir

!pick an axis to rotate about: x, y, or z
call RANDOM_NUMBER(ran)
if (ran.lt.0.3333) then
    dir = 0 ! x-direction
elseif(ran.lt.0.66667) then
    dir = 1 !y-direction
else
    dir = 2 !z-direction
endif

!get rotation matrix
call RANDOM_NUMBER(theta)
theta = rotate_max*(theta*2 - 1)
call rotation_matrix(dir, theta, m)

!get peptide center of mass (axis of rotation passes through this point)
call calc_com(group_pep, pep_com, pep_res)

!perform rotation
call rotate_peptide(group_pep, m, pep_com)

end subroutine random_peptide_rotation

subroutine shift_to_origin(group,num_res)
implicit none

real                        :: com(3)
integer                     :: num_res
type(groupdetails)          :: group(:)

call calc_com(group, com, num_res)
com = -1*com
call translate_group(group, com, num_res)

end subroutine shift_to_origin

subroutine min_peptide_z(group_pep)
implicit none
integer                 :: res, atom
type(groupdetails)      :: group_pep(pep_res)
 
pep_min_z = 1000000
do res=1,pep_res
    do atom=1,group_pep(res)%cnum1
        pep_min_z = min(group_pep(res)%coo1(3,atom), pep_min_z)
    enddo

    do atom=1,group_pep(res)%cnum2
        pep_min_z = min(group_pep(res)%coo2(3,atom), pep_min_z)
    enddo

    do atom=1,group_pep(res)%cnum3
        pep_min_z = min(group_pep(res)%coo3(3,atom), pep_min_z)
    enddo
enddo 

end subroutine min_peptide_z

subroutine min_peptide_z_backbone(group_pep)
implicit none
integer                 :: res, atom
type(groupdetails)      :: group_pep(pep_res)
 
pep_min_z = 1000000
do res=1,pep_res
    do atom=1,group_pep(res)%cnum1
        pep_min_z = min(group_pep(res)%coo1(3,atom), pep_min_z)
    enddo

    do atom=1,group_pep(res)%cnum3
        pep_min_z = min(group_pep(res)%coo3(3,atom), pep_min_z)
    enddo
enddo 

end subroutine min_peptide_z_backbone

subroutine max_receptor_z(group_rec)
implicit none
integer                 :: res, atom
type(groupdetails)      :: group_rec(rec_res)

rec_max_z = -100000
do res=1,rec_res
    do atom=1,group_rec(res)%cnum1
        rec_max_z = max(group_rec(res)%coo1(3,atom), rec_max_z)
    enddo

    do atom=1,group_rec(res)%cnum2
        rec_max_z = max(group_rec(res)%coo2(3,atom), rec_max_z)
    enddo

    do atom=1,group_rec(res)%cnum3
        rec_max_z = max(group_rec(res)%coo3(3,atom), rec_max_z)
    enddo
enddo 

end subroutine max_receptor_z

subroutine update_backup_coords(group_pep, backbone_backup)

implicit none   
type(groupdetails)      :: group_pep(:)
type(databackup)        :: backbone_backup(:)
integer                 :: i,j 

do i=1,pep_res
    do j=1,group_pep(i)%cnum1
        if(group_pep(i)%atype1(j)=="H".or.group_pep(i)%atype1(j)=="H1") then
        backbone_backup(i)%coo(1:3)=group_pep(i)%coo1(1:3,j)
        endif
    enddo
enddo

end subroutine update_backup_coords

subroutine pick_dihedral_pair(phi, psi)
implicit none
integer                     :: index
real                        :: ran, phi, psi

call RANDOM_NUMBER(ran)
index = 1+INT(ran*flavoredregion_number)
phi = REAL(flavored_region(index,1))
psi = REAL(flavored_region(index,2))

return
end subroutine pick_dihedral_pair

subroutine residue_replace(res_i, group, backbone_backup, ip, aa_group, temp_group)
implicit none
integer                         :: res_i, ip, i, j, flag
type(groupdetails)              :: group(:), temp_group(:), aa_group(40)
type(databackup)                :: backbone_backup(:)

temp_group(res_i)=group(res_i)
res_coords_flags(res_i) = 1
GB_area_calc_flags(res_i) = 1
if (temp_group(res_i)%gtype.ne.aa_group(ip)%gtype) then
    res_params_flags(res_i) = 1
    atom_links_flag = 1
endif
if(temp_group(res_i)%gtype=="PRO".or.temp_group(res_i)%gtype=="NPRO".or.temp_group(res_i)%gtype=="CPRO") then
    if(aa_group(ip)%gtype=="PRO".or.aa_group(ip)%gtype=="NPRO".or.aa_group(ip)%gtype=="CPRO") then
        temp_group(res_i)%gtype=aa_group(ip)%gtype
        temp_group(res_i)%cnum2=aa_group(ip)%cnum2
        do i=1, aa_group(ip)%cnum2
            temp_group(res_i)%atype2(i)=aa_group(ip)%atype2(i)
            temp_group(res_i)%coo2(1:3,i)=aa_group(ip)%coo2(1:3,i)
        enddo
    elseif(aa_group(ip)%gtype=="GLY".or.aa_group(ip)%gtype=="NGLY".or.aa_group(ip)%gtype=="CGLY") then
        temp_group(res_i)%gtype=aa_group(ip)%gtype
        flag=0
        do i=1, temp_group(res_i)%cnum1
            if(temp_group(res_i)%atype1(i)=="H2".or.temp_group(res_i)%atype1(i)=="H3") then
                flag=1
            endif
        enddo

        do i=1, aa_group(ip)%cnum1
            if(aa_group(ip)%atype1(i)=="H") then
                temp_group(res_i)%cnum1=temp_group(res_i)%cnum1+1
                do j=(temp_group(res_i)%cnum1-1), 1, -1
                    if(temp_group(res_i)%atype1(j)=="N") then
                        if(flag==1) then
                            temp_group(res_i)%atype1(j+1)="H1"
                        else
                            temp_group(res_i)%atype1(j+1)=aa_group(ip)%atype1(i)
                        endif
                        temp_group(res_i)%coo1(1:3,(j+1))=backbone_backup(res_i)%coo(1:3)
                        goto 10
                    else
                        temp_group(res_i)%atype1(j+1)=temp_group(res_i)%atype1(j)
                        temp_group(res_i)%coo1(1:3,(j+1))=temp_group(res_i)%coo1(1:3,j)
                    endif
                enddo
            elseif(aa_group(ip)%atype1(i)=="HA2") then
                temp_group(res_i)%cnum1=temp_group(res_i)%cnum1+1
                do j=(temp_group(res_i)%cnum1-1), 1, -1
                    if(temp_group(res_i)%atype1(j)=="CA") then
                        temp_group(res_i)%atype1(j+1)=aa_group(ip)%atype1(i)
                        temp_group(res_i)%coo1(1:3,(j+1))=aa_group(ip)%coo1(1:3,i)
                        goto 10
                    else
                        temp_group(res_i)%atype1(j+1)=temp_group(res_i)%atype1(j)
                        temp_group(res_i)%coo1(1:3,(j+1))=temp_group(res_i)%coo1(1:3,j)
                    endif
                enddo
            elseif(aa_group(ip)%atype1(i)=="HA3") then
                temp_group(res_i)%cnum1=temp_group(res_i)%cnum1+1
                do j=(temp_group(res_i)%cnum1-1), 1, -1
                    if(temp_group(res_i)%atype1(j)=="HA2") then
                        temp_group(res_i)%atype1(j+1)=aa_group(ip)%atype1(i)
                        temp_group(res_i)%coo1(1:3,(j+1))=aa_group(ip)%coo1(1:3,i)
                        goto 10
                    else
                        temp_group(res_i)%atype1(j+1)=temp_group(res_i)%atype1(j)
                        temp_group(res_i)%coo1(1:3,(j+1))=temp_group(res_i)%coo1(1:3,j)
                    endif
                enddo
            endif
10            continue
        enddo
        temp_group(res_i)%cnum1=temp_group(res_i)%cnum1-1
        temp_group(res_i)%cnum2=aa_group(ip)%cnum2
    else
        temp_group(res_i)%gtype=aa_group(ip)%gtype
        flag=0
        do i=1, temp_group(res_i)%cnum1
            if(temp_group(res_i)%atype1(i)=="H2".or.temp_group(res_i)%atype1(i)=="H3") then
                flag=1
            endif
        enddo

        do i=1, aa_group(ip)%cnum1
            if(aa_group(ip)%atype1(i)=="H") then
                temp_group(res_i)%cnum1=temp_group(res_i)%cnum1+1
                do j=(temp_group(res_i)%cnum1-1), 1, -1
                    if(temp_group(res_i)%atype1(j)=="N") then
                        if(flag==1) then
                            temp_group(res_i)%atype1(j+1)="H1"
                        else
                            temp_group(res_i)%atype1(j+1)=aa_group(ip)%atype1(i)
                        endif
                        temp_group(res_i)%coo1(1:3,(j+1))=backbone_backup(res_i)%coo(1:3)
                        goto 20
                    else
                        temp_group(res_i)%atype1(j+1)=temp_group(res_i)%atype1(j)
                        temp_group(res_i)%coo1(1:3,(j+1))=temp_group(res_i)%coo1(1:3,j)
                    endif
                enddo
            endif
        enddo
20         continue
        temp_group(res_i)%cnum2=aa_group(ip)%cnum2
        do i=1, aa_group(ip)%cnum2
            temp_group(res_i)%atype2(i)=aa_group(ip)%atype2(i)
            temp_group(res_i)%coo2(1:3,i)=aa_group(ip)%coo2(1:3,i)
        enddo
    endif
elseif(temp_group(res_i)%gtype=="GLY".or.temp_group(res_i)%gtype=="NGLY".or.temp_group(res_i)%gtype=="CGLY") then
    if(aa_group(ip)%gtype=="PRO".or.aa_group(ip)%gtype=="NPRO".or.aa_group(ip)%gtype=="CPRO") then
        temp_group(res_i)%gtype=aa_group(ip)%gtype
        flag=0
        do i=1, (temp_group(res_i)%cnum1-1)
            if(temp_group(res_i)%atype1(i)=="H1".or.temp_group(res_i)%atype1(i)=="H".or.flag==1) then
                temp_group(res_i)%atype1(i)=temp_group(res_i)%atype1(i+1)
                temp_group(res_i)%coo1(1:3,i)=temp_group(res_i)%coo1(1:3,(i+1))
                flag=1
            endif
        enddo
        temp_group(res_i)%cnum1=temp_group(res_i)%cnum1-1

        flag=0
        do i=1, (temp_group(res_i)%cnum1-1)
            if(temp_group(res_i)%atype1(i)=="HA2") then
                temp_group(res_i)%atype1(i)="HA"
            elseif(temp_group(res_i)%atype1(i)=="HA3".or.flag==1) then
                temp_group(res_i)%atype1(i)=temp_group(res_i)%atype1(i+1)
                temp_group(res_i)%coo1(1:3,i)=temp_group(res_i)%coo1(1:3,(i+1))
                flag=1
            endif
        enddo
        temp_group(res_i)%cnum1=temp_group(res_i)%cnum1-1
        temp_group(res_i)%cnum2=aa_group(ip)%cnum2
        do i=1, aa_group(ip)%cnum2
            temp_group(res_i)%atype2(i)=aa_group(ip)%atype2(i)
            temp_group(res_i)%coo2(1:3,i)=aa_group(ip)%coo2(1:3,i)
        enddo
    elseif(aa_group(ip)%gtype=="GLY".or.aa_group(ip)%gtype=="NGLY".or.aa_group(ip)%gtype=="CGLY") then
        temp_group(res_i)%gtype=aa_group(ip)%gtype
        temp_group(res_i)%cnum2=aa_group(ip)%cnum2
        do i=1, aa_group(ip)%cnum2
            temp_group(res_i)%atype2(i)=aa_group(ip)%atype2(i)
            temp_group(res_i)%coo2(1:3,i)=aa_group(ip)%coo2(1:3,i)
        enddo
    else
        temp_group(res_i)%gtype=aa_group(ip)%gtype
        flag=0
        do i=1, (temp_group(res_i)%cnum1-1)
            if(temp_group(res_i)%atype1(i)=="HA2") then
                temp_group(res_i)%atype1(i)="HA"
            elseif(temp_group(res_i)%atype1(i)=="HA3".or.flag==1) then
                temp_group(res_i)%atype1(i)=temp_group(res_i)%atype1(i+1)
                temp_group(res_i)%coo1(1:3,i)=temp_group(res_i)%coo1(1:3,(i+1))
                flag=1
            endif
        enddo
        temp_group(res_i)%cnum1=temp_group(res_i)%cnum1-1
        temp_group(res_i)%cnum2=aa_group(ip)%cnum2
        do i=1, aa_group(ip)%cnum2
            temp_group(res_i)%atype2(i)=aa_group(ip)%atype2(i)
            temp_group(res_i)%coo2(1:3,i)=aa_group(ip)%coo2(1:3,i)
        enddo
    endif
else
    if(aa_group(ip)%gtype=="PRO".or.aa_group(ip)%gtype=="NPRO".or.aa_group(ip)%gtype=="CPRO") then
        temp_group(res_i)%gtype=aa_group(ip)%gtype
        flag=0
        do i=1, (temp_group(res_i)%cnum1-1)
            if(temp_group(res_i)%atype1(i)=="H1".or.temp_group(res_i)%atype1(i)=="H".or.flag==1) then
                temp_group(res_i)%atype1(i)=temp_group(res_i)%atype1(i+1)
                temp_group(res_i)%coo1(1:3,i)=temp_group(res_i)%coo1(1:3,(i+1))
                flag=1
            endif
        enddo
        temp_group(res_i)%cnum1=temp_group(res_i)%cnum1-1
        temp_group(res_i)%cnum2=aa_group(ip)%cnum2
        do i=1, aa_group(ip)%cnum2
            temp_group(res_i)%atype2(i)=aa_group(ip)%atype2(i)
            temp_group(res_i)%coo2(1:3,i)=aa_group(ip)%coo2(1:3,i)
        enddo
    elseif(aa_group(ip)%gtype=="GLY".or.aa_group(ip)%gtype=="NGLY".or.aa_group(ip)%gtype=="CGLY") then
        temp_group(res_i)%gtype=aa_group(ip)%gtype
        do i=1, aa_group(ip)%cnum1
            if(aa_group(ip)%atype1(i)=="HA2") then
                temp_group(res_i)%cnum1=temp_group(res_i)%cnum1+1
                do j=(temp_group(res_i)%cnum1-1), 1, -1
                    if(temp_group(res_i)%atype1(j)=="CA") then
                        temp_group(res_i)%atype1(j+1)=aa_group(ip)%atype1(i)
                        temp_group(res_i)%coo1(1:3,(j+1))=aa_group(ip)%coo1(1:3,i)
                        goto 30
                    else
                        temp_group(res_i)%atype1(j+1)=temp_group(res_i)%atype1(j)
                        temp_group(res_i)%coo1(1:3,(j+1))=temp_group(res_i)%coo1(1:3,j)
                    endif
                enddo
            elseif(aa_group(ip)%atype1(i)=="HA3") then
                temp_group(res_i)%cnum1=temp_group(res_i)%cnum1+1
                do j=(temp_group(res_i)%cnum1-1), 1, -1
                    if(temp_group(res_i)%atype1(j)=="HA2") then
                        temp_group(res_i)%atype1(j+1)=aa_group(ip)%atype1(i)
                        temp_group(res_i)%coo1(1:3,(j+1))=aa_group(ip)%coo1(1:3,i)
                        goto 30
                    else
                        temp_group(res_i)%atype1(j+1)=temp_group(res_i)%atype1(j)
                        temp_group(res_i)%coo1(1:3,(j+1))=temp_group(res_i)%coo1(1:3,j)
                    endif
                enddo
            endif
30      continue
        enddo
        temp_group(res_i)%cnum1=temp_group(res_i)%cnum1-1
        temp_group(res_i)%cnum2=aa_group(ip)%cnum2
    else
        temp_group(res_i)%gtype=aa_group(ip)%gtype
        temp_group(res_i)%cnum2=aa_group(ip)%cnum2
        do i=1, aa_group(ip)%cnum2
            temp_group(res_i)%atype2(i)=aa_group(ip)%atype2(i)
            temp_group(res_i)%coo2(1:3,i)=aa_group(ip)%coo2(1:3,i)
        enddo
    endif
endif

return
end subroutine residue_replace

subroutine residue_replace_inplace(res_i, group_pep, backbone_backup, ip, aa_group)
!!!
! this subroutine inserts rotamer "ip" of residue in aa_group into residue res_i of group, storing result in group (i.e. group IS modified)
!!!
implicit none
integer                         :: res_i, ip, i, j, flag
type(groupdetails)              :: group_pep(pep_res), temp_aa, aa_group(40)
type(databackup)                :: backbone_backup(pep_res)

temp_aa = group_pep(res_i)
res_coords_flags(res_i) = 1
if (temp_aa%gtype.ne.aa_group(ip)%gtype) then
    res_params_flags(res_i) = 1
    atom_links_flag = 1
endif
if(temp_aa%gtype=="PRO".or.temp_aa%gtype=="NPRO".or.temp_aa%gtype=="CPRO") then
    if(aa_group(ip)%gtype=="PRO".or.aa_group(ip)%gtype=="NPRO".or.aa_group(ip)%gtype=="CPRO") then
        temp_aa%cnum2=aa_group(ip)%cnum2
        do i=1, aa_group(ip)%cnum2
            temp_aa%atype2(i)=aa_group(ip)%atype2(i)
            temp_aa%coo2(1:3,i)=aa_group(ip)%coo2(1:3,i)
        enddo
    elseif(aa_group(ip)%gtype=="GLY".or.aa_group(ip)%gtype=="NGLY".or.aa_group(ip)%gtype=="CGLY") then
        flag=0
        do i=1, temp_aa%cnum1
            if(temp_aa%atype1(i)=="H2".or.temp_aa%atype1(i)=="H3") then
                flag=1
            endif
        enddo
        do i=1, aa_group(ip)%cnum1
            if(aa_group(ip)%atype1(i)=="H") then
                temp_aa%cnum1=temp_aa%cnum1+1
                do j=(temp_aa%cnum1-1), 1, -1
                    if(temp_aa%atype1(j)=="N") then
                        if(flag==1) then
                            temp_aa%atype1(j+1)="H1"
                        else
                            temp_aa%atype1(j+1)=aa_group(ip)%atype1(i)
                        endif
                        temp_aa%coo1(1:3,(j+1))=backbone_backup(res_i)%coo(1:3)
                        goto 10
                    else
                        temp_aa%atype1(j+1)=temp_aa%atype1(j)
                        temp_aa%coo1(1:3,(j+1))=temp_aa%coo1(1:3,j)
                    endif
                enddo
            elseif(aa_group(ip)%atype1(i)=="HA2") then
                temp_aa%cnum1=temp_aa%cnum1+1
                do j=(temp_aa%cnum1-1), 1, -1
                    if(temp_aa%atype1(j)=="CA") then
                        temp_aa%atype1(j+1)=aa_group(ip)%atype1(i)
                        temp_aa%coo1(1:3,(j+1))=aa_group(ip)%coo1(1:3,i)
                        goto 10
                    else
                        temp_aa%atype1(j+1)=temp_aa%atype1(j)
                        temp_aa%coo1(1:3,(j+1))=temp_aa%coo1(1:3,j)
                    endif
                enddo
            elseif(aa_group(ip)%atype1(i)=="HA3") then
                temp_aa%cnum1=temp_aa%cnum1+1
                do j=(temp_aa%cnum1-1), 1, -1
                    if(temp_aa%atype1(j)=="HA2") then
                        temp_aa%atype1(j+1)=aa_group(ip)%atype1(i)
                        temp_aa%coo1(1:3,(j+1))=aa_group(ip)%coo1(1:3,i)
                        goto 10
                    else
                        temp_aa%atype1(j+1)=temp_aa%atype1(j)
                        temp_aa%coo1(1:3,(j+1))=temp_aa%coo1(1:3,j)
                    endif
                enddo
            endif
10            continue
        enddo
        temp_aa%cnum1=temp_aa%cnum1-1
        temp_aa%cnum2=aa_group(ip)%cnum2
    else
        flag=0
        do i=1, temp_aa%cnum1
            if(temp_aa%atype1(i)=="H2".or.temp_aa%atype1(i)=="H3") then
                flag=1
            endif
        enddo

        do i=1, aa_group(ip)%cnum1
            if(aa_group(ip)%atype1(i)=="H") then
                temp_aa%cnum1=temp_aa%cnum1+1
                do j=(temp_aa%cnum1-1), 1, -1
                    if(temp_aa%atype1(j)=="N") then
                        if(flag==1) then
                            temp_aa%atype1(j+1)="H1"
                        else
                            temp_aa%atype1(j+1)=aa_group(ip)%atype1(i)
                        endif
                        temp_aa%coo1(1:3,(j+1))=backbone_backup(res_i)%coo(1:3)
                        goto 20
                    else
                        temp_aa%atype1(j+1)=temp_aa%atype1(j)
                        temp_aa%coo1(1:3,(j+1))=temp_aa%coo1(1:3,j)
                    endif
                enddo
            endif
        enddo
20         continue
        temp_aa%cnum2=aa_group(ip)%cnum2
        do i=1, aa_group(ip)%cnum2
            temp_aa%atype2(i)=aa_group(ip)%atype2(i)
            temp_aa%coo2(1:3,i)=aa_group(ip)%coo2(1:3,i)
        enddo
    endif
elseif(temp_aa%gtype=="GLY".or.temp_aa%gtype=="NGLY".or.temp_aa%gtype=="CGLY") then
    if(aa_group(ip)%gtype=="PRO".or.aa_group(ip)%gtype=="NPRO".or.aa_group(ip)%gtype=="CPRO") then
        flag=0
        do i=1, (temp_aa%cnum1-1)
            if(temp_aa%atype1(i)=="H1".or.temp_aa%atype1(i)=="H".or.flag==1) then
                temp_aa%atype1(i)=temp_aa%atype1(i+1)
                temp_aa%coo1(1:3,i)=temp_aa%coo1(1:3,(i+1))
                flag=1
            endif
        enddo
        temp_aa%cnum1=temp_aa%cnum1-1
        flag=0
        do i=1, (temp_aa%cnum1-1)
            if(temp_aa%atype1(i)=="HA2") then
                temp_aa%atype1(i)="HA"
            elseif(temp_aa%atype1(i)=="HA3".or.flag==1) then
                temp_aa%atype1(i)=temp_aa%atype1(i+1)
                temp_aa%coo1(1:3,i)=temp_aa%coo1(1:3,(i+1))
                flag=1
            endif
        enddo
        temp_aa%cnum1=temp_aa%cnum1-1
        temp_aa%cnum2=aa_group(ip)%cnum2
        do i=1, aa_group(ip)%cnum2
            temp_aa%atype2(i)=aa_group(ip)%atype2(i)
            temp_aa%coo2(1:3,i)=aa_group(ip)%coo2(1:3,i)
        enddo
    elseif(aa_group(ip)%gtype=="GLY".or.aa_group(ip)%gtype=="NGLY".or.aa_group(ip)%gtype=="CGLY") then
        temp_aa%cnum2=aa_group(ip)%cnum2
        do i=1, aa_group(ip)%cnum2
            temp_aa%atype2(i)=aa_group(ip)%atype2(i)
            temp_aa%coo2(1:3,i)=aa_group(ip)%coo2(1:3,i)
        enddo
    else
        flag=0
        do i=1, (temp_aa%cnum1-1)
            if(temp_aa%atype1(i)=="HA2") then
                temp_aa%atype1(i)="HA"
            elseif(temp_aa%atype1(i)=="HA3".or.flag==1) then
                temp_aa%atype1(i)=temp_aa%atype1(i+1)
                temp_aa%coo1(1:3,i)=temp_aa%coo1(1:3,(i+1))
                flag=1
            endif
        enddo
        temp_aa%cnum1=temp_aa%cnum1-1
        temp_aa%cnum2=aa_group(ip)%cnum2
        do i=1, aa_group(ip)%cnum2
            temp_aa%atype2(i)=aa_group(ip)%atype2(i)
            temp_aa%coo2(1:3,i)=aa_group(ip)%coo2(1:3,i)
        enddo
    endif
else
    if(aa_group(ip)%gtype=="PRO".or.aa_group(ip)%gtype=="NPRO".or.aa_group(ip)%gtype=="CPRO") then
        flag=0
        do i=1, (temp_aa%cnum1-1)
            if(temp_aa%atype1(i)=="H1".or.temp_aa%atype1(i)=="H".or.flag==1) then
                temp_aa%atype1(i)=temp_aa%atype1(i+1)
                temp_aa%coo1(1:3,i)=temp_aa%coo1(1:3,(i+1))
                flag=1
            endif
        enddo
        temp_aa%cnum1=temp_aa%cnum1-1
        temp_aa%cnum2=aa_group(ip)%cnum2
        do i=1, aa_group(ip)%cnum2
            temp_aa%atype2(i)=aa_group(ip)%atype2(i)
            temp_aa%coo2(1:3,i)=aa_group(ip)%coo2(1:3,i)
        enddo
    elseif(aa_group(ip)%gtype=="GLY".or.aa_group(ip)%gtype=="NGLY".or.aa_group(ip)%gtype=="CGLY") then
        do i=1, aa_group(ip)%cnum1
            if(aa_group(ip)%atype1(i)=="HA2") then
                temp_aa%cnum1=temp_aa%cnum1+1
                do j=(temp_aa%cnum1-1), 1, -1
                    if(temp_aa%atype1(j)=="CA") then
                        temp_aa%atype1(j+1)=aa_group(ip)%atype1(i)
                        temp_aa%coo1(1:3,(j+1))=aa_group(ip)%coo1(1:3,i)
                        goto 30
                    else
                        temp_aa%atype1(j+1)=temp_aa%atype1(j)
                        temp_aa%coo1(1:3,(j+1))=temp_aa%coo1(1:3,j)
                    endif
                enddo
            elseif(aa_group(ip)%atype1(i)=="HA3") then
                temp_aa%cnum1=temp_aa%cnum1+1
                do j=(temp_aa%cnum1-1), 1, -1
                    if(temp_aa%atype1(j)=="HA2") then
                        temp_aa%atype1(j+1)=aa_group(ip)%atype1(i)
                        temp_aa%coo1(1:3,(j+1))=aa_group(ip)%coo1(1:3,i)
                        goto 30
                    else
                        temp_aa%atype1(j+1)=temp_aa%atype1(j)
                        temp_aa%coo1(1:3,(j+1))=temp_aa%coo1(1:3,j)
                    endif
                enddo
            endif
30      continue
        enddo
        temp_aa%cnum1=temp_aa%cnum1-1
        temp_aa%cnum2=aa_group(ip)%cnum2
    else
        temp_aa%cnum2=aa_group(ip)%cnum2
        do i=1, aa_group(ip)%cnum2
            temp_aa%atype2(i)=aa_group(ip)%atype2(i)
            temp_aa%coo2(1:3,i)=aa_group(ip)%coo2(1:3,i)
        enddo
    endif
endif

group_pep(res_i) = temp_aa
group_pep(res_i)%gtype=aa_group(ip)%gtype
return
end subroutine residue_replace_inplace

subroutine check_overlap(flag, res1, res2, group_pep, para_pep, feedback)
!flags and their functions:
    !flag=0: check sidechain of res num1 with rest of system
    !flag=1: check sidechain of res num1 with rest of system, while ignoring sidechain-sidechain interactions between residues num1 and num2
    !flag=2: check for overlaps between all residue pairs
    !flag=3: check for overlaps between peptide and receptor (peptide-peptide overlaps NOT checked)
implicit none
integer                         :: atom_i, atom_j, i, res1, res2, res_i, res_j, flag, feedback
real                            :: rij, cur_r(3)
real                            :: r_i, r_j, r_ij2
type(groupdetails)              :: group_pep(:)
type(energyparameters)          :: para_pep(:)

feedback=1
if(flag.eq.0 .or. flag.eq.1) then
    do atom_i=1, group_pep(res1)%cnum2
        cur_r = group_pep(res1)%coo2(1:3,atom_i)
        r_i = para_pep(res1)%r2(atom_i)
        do res_j=1, pep_res
            if (res_j.ne.res1) then
                do atom_j=1, group_pep(res_j)%cnum1
                    if(res_j==(res1+1).and.group_pep(res1)%atype2(atom_i)=="CB".and.group_pep(res_j)%atype1(atom_j)=="N") goto 53
                    if(group_pep(res1)%gtype=="PRO".or.group_pep(res1)%gtype=="NPRO".or.group_pep(res1)%gtype=="CPRO") then
                        if(res_j==(res1-1).and.group_pep(res1)%atype2(atom_i)=="CD".and.group_pep(res_j)%atype1(atom_j)=="CA") then
                            goto 53
                        endif
                    endif
                    call DIST_SQUARE(cur_r, group_pep(res_j)%coo1(1:3,atom_j), rij)
                    r_j = para_pep(res_j)%r1(atom_j)
                    if(rij<2.4025) then !2.4025 = 1.55^2
                        feedback=0
                        return
                    endif  
                53   continue
                enddo

                do atom_j=1, group_pep(res_j)%cnum3
                    if(res_j==(res1-1).and.group_pep(res1)%atype2(atom_i)=="CB".and.group_pep(res_j)%atype3(atom_j)=="C") goto 57
                    if((group_pep(res1)%gtype=="PRO".or.group_pep(res1)%gtype=="NPRO".or.group_pep(res1)%gtype=="CPRO").and.res_j==(res1-1)) then
                        if(group_pep(res1)%atype2(atom_i)=="CD") then
                            goto 57
                        elseif(group_pep(res1)%atype2(atom_i)=="HD2".or.group_pep(res1)%atype2(atom_i)=="HD3") then
                            if(group_pep(res_j)%atype3(atom_j)=="C") goto 57
                        endif
                    endif
                    call DIST_SQUARE(cur_r, group_pep(res_j)%coo3(:,atom_j), rij)
                    r_j = para_pep(res_j)%r3(atom_j)
                    if(rij<2.4025) then !2.4025 = 1.55^2
                        feedback=0
                        return
                    endif  
                57 continue
                enddo

                if (flag.eq.0 .or. res_j.ne.res2) then
                    do atom_j=1, group_pep(res_j)%cnum2
                        call DIST_SQUARE(cur_r, group_pep(res_j)%coo2(:,atom_j), rij)
                        r_j = para_pep(res_j)%r2(atom_j)
                        if(rij<2.4025) then !2.4025 = 1.55^2
                            feedback=0
                            return
                        endif  
                    enddo
                endif
            endif
        enddo

        !now check for overlaps with receptor 
        do atom_j=1, atom_num_rec
            r_j = r_rec(atom_j)
            r_ij2 = (r_j + r_i)**2
            call DIST_SQUARE(cur_r, mdcrd_rec(:,atom_j), rij)
            if(rij<2.4025) then !2.4025 = 1.55^2
                feedback=0
                return
            endif 
        enddo
    enddo

else
    do res_i=1, pep_res
        !go through group 1 atoms of residue
        do atom_i=1, group_pep(res_i)%cnum1
            cur_r = group_pep(res_i)%coo1(:,atom_i)
            r_i = para_pep(res_i)%r1(atom_i)
            !first check for overlaps with the peptide (only if flag is 2)

            if (flag.eq.2) then
                do res_j=res_i+1, pep_res 
                    do atom_j=1, group_pep(res_j)%cnum2
                        call DIST_SQUARE(cur_r, group_pep(res_j)%coo2(:,atom_j), rij)
                        r_j = para_pep(res_j)%r2(atom_j)
                        if(rij<2.4025) then !2.4025 = 1.55^2
                            feedback=0
                            return
                        endif  
                    enddo

                    do i=1, helix_num
                        if(res_i.ge.helix_start(i).and.res_i.le.helix_end(i).and.res_j.ge.helix_start(i).and.res_j.le.helix_end(i)) goto 97
                    enddo

                    do atom_j=1, group_pep(res_j)%cnum1
                        if(res_j==(res_i+1)) goto 95
                        call DIST_SQUARE(cur_r, group_pep(res_j)%coo1(:,atom_j), rij)
                        r_j = para_pep(res_j)%r1(atom_j)
                        if(rij<2.4025) then !2.4025 = 1.55^2
                            feedback=0
                            return
                        endif  
                    95  continue
                    enddo

                    do atom_j=1, group_pep(res_j)%cnum3
                        call DIST_SQUARE(cur_r, group_pep(res_j)%coo3(:,atom_j), rij)
                        r_j = para_pep(res_j)%r3(atom_j)
                        if(rij<2.4025) then !2.4025 = 1.55^2
                            feedback=0
                            return
                        endif  
                    enddo
                97 continue
                enddo
            endif

            !now check for overlaps with receptor 
            do atom_j=1, atom_num_rec
                r_j = r_rec(atom_j)
                r_ij2 = (r_j + r_i)**2
                call DIST_SQUARE(cur_r, mdcrd_rec(:,atom_j), rij)
                if(rij<2.4025) then !2.4025 = 1.55^2
                    feedback=0
                    return
                endif 
            enddo
        enddo

        !now move onto sidechain of residue
        do atom_i=1, group_pep(res_i)%cnum2
            r_i = para_pep(res_i)%r2(atom_i)
            cur_r = group_pep(res_i)%coo2(:,atom_i)

            !first check for overlaps with peptide (only if flag is 2)
            if (flag.eq.2) then
                do res_j=res_i+1, pep_res
                    do atom_j=1, group_pep(res_j)%cnum1
                        if(res_j==(res_i+1).and.group_pep(res_i)%atype2(atom_i)=="CB".and.group_pep(res_j)%atype1(atom_j)=="N") goto 100
                        call DIST_SQUARE(cur_r, group_pep(res_j)%coo1(:,atom_j), rij)
                        r_j = para_pep(res_j)%r1(atom_j)
                        if(rij<2.4025) then !2.4025 = 1.55^2
                            feedback=0
                            return
                        endif  
                    100 continue
                    enddo

                    do atom_j=1, group_pep(res_j)%cnum2
                        call DIST_SQUARE(cur_r, group_pep(res_j)%coo2(:,atom_j), rij)
                        r_j = para_pep(res_j)%r2(atom_j)
                        if(rij<2.4025) then !2.4025 = 1.55^2
                            feedback=0
                            return
                        endif  
                    enddo

                    do atom_j=1, group_pep(res_j)%cnum3
                        call DIST_SQUARE(cur_r, group_pep(res_j)%coo3(:,atom_j), rij)
                        r_j = para_pep(res_j)%r3(atom_j)
                        if(rij<2.4025) then !2.4025 = 1.55^2
                            feedback=0
                            return
                        endif  
                    enddo
                enddo
            endif

            !now check for overlaps with receptor 
            do atom_j=1, atom_num_rec
                r_j = r_rec(atom_j)
                r_ij2 = (r_j + r_i)**2
                call DIST_SQUARE(cur_r, mdcrd_rec(:,atom_j), rij)
                if(rij<2.4025) then !2.4025 = 1.55^2
                    feedback=0
                    return
                endif 
            enddo
        enddo

        do atom_i=1, group_pep(res_i)%cnum3
            r_i = para_pep(res_i)%r3(atom_i)
            cur_r = group_pep(res_i)%coo3(:,atom_i)

            !first check for overlaps with peptide (only if flag is 2)
            if (flag.eq.2) then
                do res_j=res_i+1, pep_res
                    do atom_j=1, group_pep(res_j)%cnum2
                        if(res_j==(res_i+1)) goto 110
                        call DIST_SQUARE(cur_r, group_pep(res_j)%coo2(:,atom_j), rij)
                        r_j = para_pep(res_j)%r2(atom_j)
                        if(rij<2.4025) then !2.4025 = 1.55^2
                            feedback=0
                            return
                        endif  
                    110 continue
                    enddo

                    do i=1, helix_num
                        if(res_i.ge.helix_start(i).and.res_i.le.helix_end(i).and.res_j.ge.helix_start(i).and.res_j.le.helix_end(i)) goto 120
                    enddo

                    do atom_j=1, group_pep(res_j)%cnum1
                        if(res_j==(res_i+1)) goto 105
                        call DIST_SQUARE(cur_r, group_pep(res_j)%coo1(:,atom_j), rij)
                        r_j = para_pep(res_j)%r1(atom_j)
                        if(rij<2.4025) then !2.4025 = 1.55^2
                            feedback=0
                            return
                        endif  
                    105  continue
                    enddo
                    
                    do atom_j=1, group_pep(res_j)%cnum3
                        if(res_j==(res_i+1)) goto 115
                        call DIST_SQUARE(cur_r, group_pep(res_j)%coo3(:,atom_j), rij)
                        r_j = para_pep(res_j)%r3(atom_j)
                        if(rij<2.4025) then !2.4025 = 1.55^2
                            feedback=0
                            return
                        endif  
                    115  continue
                    enddo
                120 continue
                enddo
            endif

            !now check for overlaps with receptor 
            do atom_j=1, atom_num_rec
                r_j = r_rec(atom_j)
                r_ij2 = (r_j + r_i)**2
                call DIST_SQUARE(cur_r, mdcrd_rec(:,atom_j), rij)
                if(rij<2.4025) then !2.4025 = 1.55^2
                    feedback=0
                    return
                endif 
            enddo
        enddo
    enddo
endif

return
end subroutine check_overlap

subroutine check_backbone_overlap(group_pep, para_pep, feedback)
!this subroutine looks for steric clashes between peptide backbone with either itself or the receptor. Feedback is 0 if no clashes, 1 if there is a clash
implicit none
integer                 :: feedback
integer                 :: i, j, res_i, res_j
real                    :: rij, r_i, r_j, r_ij2, cur_r(3)
type(groupdetails)      :: group_pep(pep_res)
type(energyparameters)  :: para_pep(pep_res)

feedback=1
do res_i=1, pep_res
    do i=1, group_pep(res_i)%cnum1
        r_i = para_pep(res_i)%r1(i)
        cur_r = group_pep(res_i)%coo1(1:3,i)
        
        !first check clashes with other peptide residues
        do res_j=res_i+2, pep_res !skip two residues ahead when looking for overlaps (backbone atoms for adjacent residues will always be within cutoff distance)
            do j=1, helix_num
                if(res_i.ge.helix_start(j).and.res_i.le.helix_end(j).and.res_j.ge.helix_start(j).and.res_j.le.helix_end(j)) goto 20
            enddo

            do j=1, group_pep(res_j)%cnum1
                r_j = para_pep(res_j)%r1(j)
                r_ij2 = (r_j + r_i)**2
                call DIST_SQUARE(cur_r, group_pep(res_j)%coo1(:,j), rij)
                if (rij*min_rij2_ratio < r_ij2) then
                    feedback=0
                    return
                endif
            enddo

            do j=1, group_pep(res_j)%cnum3
                r_j = para_pep(res_j)%r3(j)
                r_ij2 = (r_j + r_i)**2
                call DIST_SQUARE(cur_r, group_pep(res_j)%coo3(:,j), rij)
                if (rij*min_rij2_ratio < r_ij2) then
                    feedback=0
                    return
                endif
            enddo

20      continue
        enddo

        !now check for clashes with the receptor
        do j=1, atom_num_rec
            r_j = r_rec(j)
            r_ij2 = (r_j + r_i)**2
            call DIST_SQUARE(cur_r, mdcrd_rec(:,j), rij)
            if (rij*min_rij2_ratio < r_ij2) then
                feedback=0
                return
            endif
        enddo
    enddo

    do i=1, group_pep(res_i)%cnum3
        r_i = para_pep(res_i)%r3(i)
        cur_r = group_pep(res_i)%coo3(:,i)

        !first check clashes with other peptide residues
        do res_j=res_i+2, pep_res
            do j=1, helix_num
                if(res_i.ge.helix_start(j).and.res_i.le.helix_end(j).and.res_j.ge.helix_start(j).and.res_j.le.helix_end(j)) goto 30
            enddo

            do j=1, group_pep(res_j)%cnum1
                r_j = para_pep(res_j)%r1(j)
                r_ij2 = (r_j + r_i)**2
                call DIST_SQUARE(cur_r, group_pep(res_j)%coo1(:,j), rij)
                if (rij*min_rij2_ratio < r_ij2) then
                    feedback=0
                    return
                endif
            enddo

            do j=1, group_pep(res_j)%cnum3
                r_j = para_pep(res_j)%r3(j)
                r_ij2 = (r_j + r_i)**2
                call DIST_SQUARE(cur_r, group_pep(res_j)%coo3(:,j), rij)
                if (rij*min_rij2_ratio < r_ij2) then
                    feedback=0
                    return
                endif
            enddo

30      continue
        enddo

        !now check for clashes with the receptor
        do j=1, atom_num_rec
            r_j = r_rec(j)
            r_ij2 = (r_j + r_i)**2
            call DIST_SQUARE(cur_r, mdcrd_rec(:,j), rij)
            if (rij*min_rij2_ratio < r_ij2) then
                feedback=0
                return
            endif
        enddo
    enddo
enddo


return
end subroutine check_backbone_overlap

subroutine read_sequences(sequences, count)
implicit none
integer                     :: count, i, ios
character*4,allocatable     :: sequences(:,:)

count = 0
!count number of entries in file
open(101, file=sequences_file)

do 
read (101, '(A)', iostat=ios)
if (ios .ne. 0) exit
count = count + 1
enddo

rewind(unit=101)
allocate(sequences(pep_res,count))

do i=1, count
    read(101,*) sequences(:,i)
enddo

close(101)

end subroutine read_sequences

subroutine read_backbone_files(backbones, count)
implicit none
integer                        :: count, i, ios
character*50, allocatable      :: backbones(:)
character*100                  :: file

count = 0
file = trim(backbones_path)//trim(backbones_file)

!count number of entries in file
open(101, file=file)

do 
read (101, '(A)', iostat=ios)
if (ios .ne. 0) exit
count = count + 1
enddo

rewind(unit=101)
allocate(backbones(count))

do i=1, count
    read(101,"(A)") backbones(i)
enddo

close(101)

end subroutine read_backbone_files

subroutine calc_sidechain_dihedrals(group, res_i, dihedrals, num_dihedrals)
!this subroutine calculates the dihedrals of residue "res_i" in peptide "group" 
!NOTE: The final dihedral angle for proline are NOT calculated correctly at moment, I'm not sure what atoms used to define this dihedral
implicit none
type(groupdetails)              :: group(:)
real                            :: atom_locs(3,9), angle
integer                         :: dihedrals(:,:), res_i, num_dihedrals, atom_i, i
character*4                     :: res_name

!first get the Calpha and N positions (always used for sidechain dihedrals)
do atom_i=1, group(res_i)%cnum1
    if (group(res_i)%atype1(atom_i) .eq. "N") then
        atom_locs(1:3,1) = group(res_i)%coo1(1:3,atom_i) !nitrogen goes into position 1 of atom locs
    elseif(group(res_i)%atype1(atom_i) .eq. "CA") then
        atom_locs(1:3,2) = group(res_i)%coo1(1:3,atom_i) !alpha carbon goes into position 2 of atom locs
    endif
enddo

res_name = trim(group(res_i)%gtype)
if (res_i.eq.1 .or. res_i .eq. pep_res) res_name = res_name(2:4) !trim residue names if at N- or C-terminus

!now get dihedral angles (this is amino acid specific)
if (res_name.eq.'ALA' .or. res_name.eq.'GLY') then
    num_dihedrals = 0  !these residues have no sidechains

else if (res_name.eq.'ARG') then !Confirmed to calculate angles correctly!
    num_dihedrals = 4 !total number of sidechain dihedral angles
    !get position of relevant sidechain atoms that define dihedral angles
    do atom_i=1, group(res_i)%cnum2
        if (group(res_i)%atype2(atom_i) .eq. "CB") then
            atom_locs(:,3) = group(res_i)%coo2(:,atom_i) 
        elseif (group(res_i)%atype2(atom_i) .eq. "CG") then
            atom_locs(:,4) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "CD") then
            atom_locs(:,5) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "NE") then
            atom_locs(:,6) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "CZ") then
            atom_locs(:,7) = group(res_i)%coo2(:,atom_i)
        endif
    enddo
    
else if (res_name.eq.'HIE') then !Confirmed to calculate angles correctly!
    num_dihedrals = 2!total number of sidechain dihedral angles
    !get position of relevant sidechain atoms that define dihedral angles
    do atom_i=1, group(res_i)%cnum2
        if (group(res_i)%atype2(atom_i) .eq. "CB") then
            atom_locs(:,3) = group(res_i)%coo2(:,atom_i) 
        elseif (group(res_i)%atype2(atom_i) .eq. "CG") then
            atom_locs(:,4) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "ND1") then
            atom_locs(:,5) = group(res_i)%coo2(:,atom_i)
        endif
    enddo

else if (res_name.eq.'LYS') then
    num_dihedrals = 5!total number of sidechain dihedral angles
    !get position of relevant sidechain atoms that define dihedral angles
    do atom_i=1, group(res_i)%cnum2
        if (group(res_i)%atype2(atom_i) .eq. "CB") then
            atom_locs(:,3) = group(res_i)%coo2(:,atom_i) 
        elseif (group(res_i)%atype2(atom_i) .eq. "CG") then
            atom_locs(:,4) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "CD") then
            atom_locs(:,5) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "CE") then
            atom_locs(:,6) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "NZ") then
            atom_locs(:,7) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "HZ1") then !need to check if this is the "correct" hydrogen used to define final dihedral
            atom_locs(:,8) = group(res_i)%coo2(:,atom_i)
        endif
    enddo

else if (res_name.eq.'ASP') then
    num_dihedrals = 2 !total number of sidechain dihedral angles
    !get position of relevant sidechain atoms that define dihedral angles
    do atom_i=1, group(res_i)%cnum2
        if (group(res_i)%atype2(atom_i) .eq. "CB") then
            atom_locs(:,3) = group(res_i)%coo2(:,atom_i) 
        elseif (group(res_i)%atype2(atom_i) .eq. "CG") then
            atom_locs(:,4) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "OD1") then !need to check if this is the "correct" oxygen used to define dihedral
            atom_locs(:,5) = group(res_i)%coo2(:,atom_i)
        endif
    enddo

else if (res_name.eq.'GLU') then
    num_dihedrals = 3 !total number of sidechain dihedral angles
    !get position of relevant sidechain atoms that define dihedral angles
    do atom_i=1, group(res_i)%cnum2
        if (group(res_i)%atype2(atom_i) .eq. "CB") then
            atom_locs(:,3) = group(res_i)%coo2(:,atom_i) 
        elseif (group(res_i)%atype2(atom_i) .eq. "CG") then
            atom_locs(:,4) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "CD") then
            atom_locs(:,5) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "OE1") then !need to check if this is the "correct" oxygen used to define dihedral
            atom_locs(:,6) = group(res_i)%coo2(:,atom_i)
        endif
    enddo

else if (res_name.eq.'SER') then
    num_dihedrals = 2 !total number of sidechain dihedral angles
    !get position of relevant sidechain atoms that define dihedral angles
    do atom_i=1, group(res_i)%cnum2
        if (group(res_i)%atype2(atom_i) .eq. "CB") then
            atom_locs(:,3) = group(res_i)%coo2(:,atom_i) 
        elseif (group(res_i)%atype2(atom_i) .eq. "OG") then
            atom_locs(:,4) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "HG") then
            atom_locs(:,5) = group(res_i)%coo2(:,atom_i)
        endif
    enddo

else if (res_name.eq.'THR') then
    num_dihedrals = 2 !total number of sidechain dihedral angles
    !get position of relevant sidechain atoms that define dihedral angles
    do atom_i=1, group(res_i)%cnum2
        if (group(res_i)%atype2(atom_i) .eq. "CB") then
            atom_locs(:,3) = group(res_i)%coo2(:,atom_i) 
        elseif (group(res_i)%atype2(atom_i) .eq. "OG1") then
            atom_locs(:,4) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "HG1") then
            atom_locs(:,5) = group(res_i)%coo2(:,atom_i)
        endif
    enddo

else if (res_name.eq.'ASN') then
    num_dihedrals = 2 !total number of sidechain dihedral angles
    !get position of relevant sidechain atoms that define dihedral angles
    do atom_i=1, group(res_i)%cnum2
        if (group(res_i)%atype2(atom_i) .eq. "CB") then
            atom_locs(:,3) = group(res_i)%coo2(:,atom_i) 
        elseif (group(res_i)%atype2(atom_i) .eq. "CG") then
            atom_locs(:,4) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "OD1") then
            atom_locs(:,5) = group(res_i)%coo2(:,atom_i)
        endif
    enddo

else if (res_name.eq.'GLN') then
    num_dihedrals = 3 !total number of sidechain dihedral angles
    !get position of relevant sidechain atoms that define dihedral angles
    do atom_i=1, group(res_i)%cnum2
        if (group(res_i)%atype2(atom_i) .eq. "CB") then
            atom_locs(:,3) = group(res_i)%coo2(:,atom_i) 
        elseif (group(res_i)%atype2(atom_i) .eq. "CG") then
            atom_locs(:,4) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "CD") then
            atom_locs(:,5) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "OE1") then
            atom_locs(:,6) = group(res_i)%coo2(:,atom_i)
        endif
    enddo

else if (res_name.eq.'CYS') then
    num_dihedrals = 2 !total number of sidechain dihedral angles
    !get position of relevant sidechain atoms that define dihedral angles
    do atom_i=1, group(res_i)%cnum2
        if (group(res_i)%atype2(atom_i) .eq. "CB") then
            atom_locs(:,3) = group(res_i)%coo2(:,atom_i) 
        elseif (group(res_i)%atype2(atom_i) .eq. "SG") then
            atom_locs(:,4) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "HG") then
            atom_locs(:,5) = group(res_i)%coo2(:,atom_i)
        endif
    enddo

else if (res_name.eq.'PRO') then
    num_dihedrals = 3 !total number of sidechain dihedral angles
    !get position of relevant sidechain atoms that define dihedral angles
    atom_locs(:,6) = atom_locs(:,1) !last dihedral wraps back around to nitrogen
    do atom_i=1, group(res_i)%cnum2
        if (group(res_i)%atype2(atom_i) .eq. "CB") then
            atom_locs(:,3) = group(res_i)%coo2(:,atom_i) 
        elseif (group(res_i)%atype2(atom_i) .eq. "CG") then
            atom_locs(:,4) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "CD") then
            atom_locs(:,5) = group(res_i)%coo2(:,atom_i)
        endif
    enddo

else if (res_name.eq.'ILE') then
    num_dihedrals = 2 !total number of sidechain dihedral angles
    !get position of relevant sidechain atoms that define dihedral angles
    do atom_i=1, group(res_i)%cnum2
        if (group(res_i)%atype2(atom_i) .eq. "CB") then
            atom_locs(:,3) = group(res_i)%coo2(:,atom_i) 
        elseif (group(res_i)%atype2(atom_i) .eq. "CG1") then
            atom_locs(:,4) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "CD1") then
            atom_locs(:,5) = group(res_i)%coo2(:,atom_i)
        endif
    enddo

else if (res_name.eq.'LEU') then
    num_dihedrals = 2 !total number of sidechain dihedral angles
    !get position of relevant sidechain atoms that define dihedral angles
    do atom_i=1, group(res_i)%cnum2
        if (group(res_i)%atype2(atom_i) .eq. "CB") then
            atom_locs(:,3) = group(res_i)%coo2(:,atom_i) 
        elseif (group(res_i)%atype2(atom_i) .eq. "CG") then
            atom_locs(:,4) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "CD1") then
            atom_locs(:,5) = group(res_i)%coo2(:,atom_i)
        endif
    enddo

else if (res_name.eq.'MET') then
    num_dihedrals = 4 !total number of sidechain dihedral angles
    !get position of relevant sidechain atoms that define dihedral angles
    do atom_i=1, group(res_i)%cnum2
        if (group(res_i)%atype2(atom_i) .eq. "CB") then
            atom_locs(:,3) = group(res_i)%coo2(:,atom_i) 
        elseif (group(res_i)%atype2(atom_i) .eq. "CG") then
            atom_locs(:,4) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "SD") then
            atom_locs(:,5) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "CE") then
            atom_locs(:,6) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "HE1") then
            atom_locs(:,7) = group(res_i)%coo2(:,atom_i)
        endif
    enddo

else if (res_name.eq.'PHE') then
    num_dihedrals = 2 !total number of sidechain dihedral angles
    !get position of relevant sidechain atoms that define dihedral angles
    do atom_i=1, group(res_i)%cnum2
        if (group(res_i)%atype2(atom_i) .eq. "CB") then
            atom_locs(:,3) = group(res_i)%coo2(:,atom_i) 
        elseif (group(res_i)%atype2(atom_i) .eq. "CG") then
            atom_locs(:,4) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "CD1") then
            atom_locs(:,5) = group(res_i)%coo2(:,atom_i)
        endif
    enddo

else if (res_name.eq.'TRP') then 
    num_dihedrals = 2 !total number of sidechain dihedral angles
    !get position of relevant sidechain atoms that define dihedral angles
    do atom_i=1, group(res_i)%cnum2
        if (group(res_i)%atype2(atom_i) .eq. "CB") then
            atom_locs(:,3) = group(res_i)%coo2(:,atom_i) 
        elseif (group(res_i)%atype2(atom_i) .eq. "CG") then
            atom_locs(:,4) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "CD1") then
            atom_locs(:,5) = group(res_i)%coo2(:,atom_i)
        endif
    enddo

else if (res_name.eq.'TYR') then !NEED TO DOUBLE CHECK IF BELOW IS CORRECT!
    num_dihedrals = 3 !total number of sidechain dihedral angles
    !get position of relevant sidechain atoms that define dihedral angles
    do atom_i=1, group(res_i)%cnum2
        if (group(res_i)%atype2(atom_i) .eq. "CB") then
            atom_locs(:,3) = group(res_i)%coo2(:,atom_i) 
        elseif (group(res_i)%atype2(atom_i) .eq. "CG") then
            atom_locs(:,4) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "CD1") then
            atom_locs(:,5) = group(res_i)%coo2(:,atom_i)    
        elseif (group(res_i)%atype2(atom_i) .eq. "CE1") then
            atom_locs(:,6) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "CZ") then
            atom_locs(:,7) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "OH") then
            atom_locs(:,8) = group(res_i)%coo2(:,atom_i)
        elseif (group(res_i)%atype2(atom_i) .eq. "HH") then
            atom_locs(:,9) = group(res_i)%coo2(:,atom_i)
        endif
    enddo

else if (res_name.eq.'VAL') then
    num_dihedrals = 1 !total number of sidechain dihedral angles
    !get position of relevant sidechain atoms that define dihedral angles
    do atom_i=1, group(res_i)%cnum2
        if (group(res_i)%atype2(atom_i) .eq. "CB") then
            atom_locs(:,3) = group(res_i)%coo2(:,atom_i) 
        elseif (group(res_i)%atype2(atom_i) .eq. "CG1") then
            atom_locs(:,4) = group(res_i)%coo2(:,atom_i) 
        endif
    enddo
    
else
    write(*,*) "Did not recognize residue when trying to calculate sidechain dihedrals!"
    write(*,*) "This error should not normally happen, there is likely a problem with the code."
    write(*,*) "Terminating program."
endif

!now calculate the dihedral angles
if (res_name .eq. 'ARG') then !have a special case for arginine as order for final dihedral atoms is reversed to align with values in rotamer file
    do i=1, num_dihedrals-1
        call calc_dihedral(atom_locs(:,i), atom_locs(:,i+1), atom_locs(:,i+2), atom_locs(:,i+3), angle)
        dihedrals(i, res_i) = -1*NINT(angle) !flip sign of angle to align with tabulated values
    enddo

    call calc_dihedral(atom_locs(:,7), atom_locs(:,6), atom_locs(:,5), atom_locs(:,4), angle)
    dihedrals(num_dihedrals, res_i) = -1*NINT(angle) !flip sign of angle to align with tabulated values

elseif(res_name .eq. 'PRO') then !have special case for tyrosine as last dihedral defined by 4 atoms unique from all other dihedrals

    write(*,*) "NOTE: the final dihedral angle for proline is NOT correct!"
    write(*,*) "I can't find out what combination of atoms is used to define this angle..."
    do i=1, num_dihedrals
        call calc_dihedral(atom_locs(:,i), atom_locs(:,i+1), atom_locs(:,i+2), atom_locs(:,i+3), angle)
        dihedrals(i, res_i) = -1*NINT(angle) !flip sign of angle to align with tabulated values
    enddo

    call calc_dihedral(atom_locs(:,6), atom_locs(:,5), atom_locs(:,4), atom_locs(:,3), angle)
    dihedrals(num_dihedrals, res_i) = -1*NINT(angle) !flip sign of angle to align with tabulated values

elseif(res_name .eq. 'TRP') then !have special case for tyrosine as last dihedral defined by 4 atoms unique from all other dihedrals

    do i=1, num_dihedrals
        call calc_dihedral(atom_locs(:,i), atom_locs(:,i+1), atom_locs(:,i+2), atom_locs(:,i+3), angle)
        dihedrals(i, res_i) = -1*NINT(angle) !flip sign of angle to align with tabulated values
    enddo

elseif(res_name .eq. 'TYR') then !have special case for tyrosine as last dihedral defined by 4 atoms unique from all other dihedrals

    do i=1, num_dihedrals-1
        call calc_dihedral(atom_locs(:,i), atom_locs(:,i+1), atom_locs(:,i+2), atom_locs(:,i+3), angle)
        dihedrals(i, res_i) = -1*NINT(angle) !flip sign of angle to align with tabulated values
    enddo

    call calc_dihedral(atom_locs(:,9), atom_locs(:,8), atom_locs(:,7), atom_locs(:,6), angle)
    dihedrals(num_dihedrals, res_i) = -1*NINT(angle) !flip sign of angle to align with tabulated values
else

    do i=1, num_dihedrals
    call calc_dihedral(atom_locs(:,i), atom_locs(:,i+1), atom_locs(:,i+2), atom_locs(:,i+3), angle)
    dihedrals(i, res_i) = -1*NINT(angle) !flip sign of angle to align with tabulated values

enddo

endif

end subroutine calc_sidechain_dihedrals

subroutine calc_backbone_dihedrals(group, phis, psis)
implicit none
integer                     :: res_i, k
real                        :: phi, psi, phis(pep_res), psis(pep_res), C_prev(3), N(3), CA(3), C(3), N_next(3)
type(groupdetails)          :: group(pep_res)

do res_i=1, pep_res
    !!get coordinates of C_prev, if not residue 1 
    if (res_i.eq.1) then
        phi = 100000 !1st residue does not have this atom, so does not have a phi angle
    else 
        do k=1, group(res_i-1)%cnum3
            if(group(res_i-1)%atype3(k)=="C") then
                C_prev=group(res_i-1)%coo3(:,k)
            endif
        enddo
    endif

    !get coordinate of N and alpha carbon
    do k=1, group(res_i)%cnum1
        if(group(res_i)%atype1(k)=="N") then
            N=group(res_i)%coo1(:,k)
        elseif(group(res_i)%atype1(k)=="CA") then
            CA=group(res_i)%coo1(:,k)
        endif
    enddo

    !Get coordinate of current C atom
    do k=1, group(res_i)%cnum3
        if(group(res_i)%atype3(k)=="C") then
            C=group(res_i)%coo3(:,k)
        endif
    enddo

    !get coordinate of next N atom if not the final residue
    if (res_i.eq.pep_res) then
        psi = 100000
    else
        do k=1, group(res_i+1)%cnum1
            if(group(res_i+1)%atype1(k)=="N") then
                N_next=group(res_i+1)%coo1(:,k)
            endif
        enddo
    endif

    !calculate phi and psi
    if (res_i.ne.1) call calc_dihedral(C_prev, N, CA, C, phi)
    if (res_i.ne.pep_res) call calc_dihedral(N, CA, C, N_next, psi)

    !store the dihedral angles
    phis(res_i) = phi
    psis(res_i) = psi
enddo
return

end subroutine calc_backbone_dihedrals

subroutine load_params_single(group, group_para, flag, natom, res_i)
implicit none
integer                         :: natom, atom_offset
integer                         :: res_i, j, flag
type(groupdetails)              :: group(:)
type(energyparameters)          :: group_para(:)

if (flag.eq.1) then !loading parameters for peptide 
    natom = 0
    if (res_i.eq.1) then
        atom_offset = 0
    else
        atom_offset = res_ends(res_i-1)
    endif   
    res_params_flags(res_i) = 0
    do j=1, group(res_i)%cnum1
        natom=natom+1
        epsion_pep(natom, res_i)=group_para(res_i)%epsion1(j)
        r_pep(natom, res_i)=group_para(res_i)%r1(j)
        charge_pep(natom, res_i)=group_para(res_i)%charge1(j)
        rborn_pep(natom, res_i)=group_para(res_i)%rborn1(j)
        fs_pep(natom, res_i)=group_para(res_i)%fs1(j)
        dielecons_pep(natom, res_i)=group_para(res_i)%dielecons1(j)
        atomid_pep(natom, res_i)=atom_offset+group_para(res_i)%atomid1(j)
    enddo

    do j=1, group(res_i)%cnum2
        natom=natom+1
        epsion_pep(natom, res_i)=group_para(res_i)%epsion2(j)
        r_pep(natom, res_i)=group_para(res_i)%r2(j)
        charge_pep(natom, res_i)=group_para(res_i)%charge2(j)
        rborn_pep(natom, res_i)=group_para(res_i)%rborn2(j)
        fs_pep(natom, res_i)=group_para(res_i)%fs2(j)
        dielecons_pep(natom, res_i)=group_para(res_i)%dielecons2(j)
        atomid_pep(natom, res_i)=atom_offset+group_para(res_i)%atomid2(j)
    enddo

    do j=1, group(res_i)%cnum3
        natom=natom+1
        epsion_pep(natom, res_i)=group_para(res_i)%epsion3(j)
        r_pep(natom, res_i)=group_para(res_i)%r3(j)
        charge_pep(natom, res_i)=group_para(res_i)%charge3(j)
        rborn_pep(natom, res_i)=group_para(res_i)%rborn3(j)
        fs_pep(natom, res_i)=group_para(res_i)%fs3(j)
        dielecons_pep(natom, res_i)=group_para(res_i)%dielecons3(j)
        atomid_pep(natom, res_i)=atom_offset+group_para(res_i)%atomid3(j)
    enddo
    
else !loading parameters receptor
    atom_offset = natom
    do j=1, group(res_i)%cnum1 
        natom=natom+1
        epsion_rec(natom)=group_para(res_i)%epsion1(j)
        r_rec(natom)=group_para(res_i)%r1(j)
        charge_rec(natom)=group_para(res_i)%charge1(j)
        rborn_rec(natom)=group_para(res_i)%rborn1(j)
        fs_rec(natom)=group_para(res_i)%fs1(j)
        dielecons_rec(natom)=group_para(res_i)%dielecons1(j)
        !atomid not needed for receptor, don't need to worry about atom bonding
    enddo

    do j=1, group(res_i)%cnum2
        natom=natom+1
        epsion_rec(natom)=group_para(res_i)%epsion2(j)
        r_rec(natom)=group_para(res_i)%r2(j)
        charge_rec(natom)=group_para(res_i)%charge2(j)
        rborn_rec(natom)=group_para(res_i)%rborn2(j)
        fs_rec(natom)=group_para(res_i)%fs2(j)
        dielecons_rec(natom)=group_para(res_i)%dielecons2(j)
        !atomid not needed for receptor, don't need to worry about atom bonding
    enddo

    do j=1, group(res_i)%cnum3
        natom=natom+1
        epsion_rec(natom)=group_para(res_i)%epsion3(j)
        r_rec(natom)=group_para(res_i)%r3(j)
        charge_rec(natom)=group_para(res_i)%charge3(j)
        rborn_rec(natom)=group_para(res_i)%rborn3(j)
        fs_rec(natom)=group_para(res_i)%fs3(j)
        dielecons_rec(natom)=group_para(res_i)%dielecons3(j)
        !atomid not needed for receptor, don't need to worry about atom bonding
    enddo
endif

end subroutine load_params_single

subroutine load_params(group, group_para, flag)
implicit none
integer                         :: natom
integer                         :: res_i, flag
type(groupdetails)              :: group(:)
type(energyparameters)          :: group_para(:)

natom=0
if (flag.eq.1) then !loading parameters for peptide
    do res_i=1,pep_res
        call load_params_single(group, group_para, flag, natom, res_i)
    enddo
    atom_num_pep = natom
else !loading parameters for receptor
    do res_i=1,rec_res
        call load_params_single(group, group_para, flag, natom, res_i)
    enddo
endif

end subroutine load_params

subroutine load_coords_single(group, flag, natom, res_i)
implicit none
integer                         :: res_i, j, flag, natom
type(groupdetails)              :: group(:)

if (flag.eq.1) then !store coordinates for peptide
    natom = 0
    res_coords_flags(res_i) = 0
    do j=1, group(res_i)%cnum1  
        natom=natom+1
        mdcrd_pep(:,natom, res_i)=group(res_i)%coo1(:,j)
    enddo
    
    do j=1, group(res_i)%cnum2
        natom=natom+1
        mdcrd_pep(:,natom, res_i)=group(res_i)%coo2(:,j)
    enddo
    
    do j=1, group(res_i)%cnum3
        natom=natom+1
        mdcrd_pep(:,natom, res_i)=group(res_i)%coo3(:,j)
    enddo
else 
    do j=1, group(res_i)%cnum1  
        natom=natom+1
        mdcrd_rec(:,natom)=group(res_i)%coo1(:,j)
    enddo
    
    do j=1, group(res_i)%cnum2
        natom=natom+1
        mdcrd_rec(:,natom)=group(res_i)%coo2(:,j)
    enddo
    
    do j=1, group(res_i)%cnum3
        natom=natom+1
        mdcrd_rec(:,natom)=group(res_i)%coo3(:,j)
    enddo
endif

end subroutine load_coords_single

subroutine load_coords(group, flag)
implicit none
integer                         :: natom
integer                         :: res_i, flag, cg_flags(pep_res)
type(groupdetails)              :: group(:)

if (flag.eq.1) then !loading coordinates for peptide
    cg_flags = 0
    call load_res_cutoffs(group, cg_flags)
    atom_num_pep = res_ends(pep_res)
    do res_i=1, pep_res
        natom=0
        call load_coords_single(group,flag, natom, res_i)
    enddo
    GB_area_calc_flags = 1
else !loading coordinates for receptor
    natom=0
    do res_i=1, rec_res
        call load_coords_single(group,flag, natom, res_i)
    enddo
endif

end subroutine load_coords

subroutine load_peptide_info(group_pep, para_pep)
implicit none
integer                         :: res_i, cg_flags(pep_res), natom
type(groupdetails)              :: group_pep(pep_res)
type(energyparameters)          :: para_pep(pep_res)

!get starts/ends for whole peptide and sidechain of peptide residues
cg_flags = 0
call load_res_cutoffs(group_pep, cg_flags)
atom_num_pep = res_ends(pep_res)
do res_i=1, pep_res
    if (res_coords_flags(res_i).eq.1) then 
        !load coordinates
        call load_coords_single(group_pep,1, natom, res_i)
        res_coords_flags(res_i) = 1
    endif
    if (res_params_flags(res_i).eq.1) then
        call load_params_single(group_pep,para_pep, 1, natom, res_i)
    endif
enddo

do res_i=1, pep_res
    if (res_coords_flags(res_i).eq.1) then 
        call calc_pairwise_distances(res_i)
        res_coords_flags(res_i) = 0
    endif
enddo

!form atom links for peptide
if (atom_links_flag.eq.1) then
    call atom_links(group_pep, cg_flags)
    atom_links_flag = 0
endif

end subroutine load_peptide_info

subroutine load_params_single_cg(group, group_para, atom_offset, res_i, cg_flag)
!load parameters for peptide, using coarse-grained bead at beta-carbon if cg_flag is 0
implicit none
integer                         :: natom, atom_offset
integer                         :: res_i, j, cg_flag
type(groupdetails)              :: group(:)
type(energyparameters)          :: group_para(:)

res_params_flags(res_i) = 0
natom = 0
do j=1, group(res_i)%cnum1 
    natom=natom+1
    epsion_pep(natom, res_i)=group_para(res_i)%epsion1(j)
    r_pep(natom, res_i)=group_para(res_i)%r1(j)
    charge_pep(natom, res_i)=group_para(res_i)%charge1(j)
    rborn_pep(natom, res_i)=group_para(res_i)%rborn1(j)
    fs_pep(natom, res_i)=group_para(res_i)%fs1(j)
    dielecons_pep(natom, res_i)=group_para(res_i)%dielecons1(j)
    atomid_pep(natom, res_i)=atom_offset+group_para(res_i)%atomid1(j)
enddo

if (cg_flag.eq.0) then !we are not coarse-graining this residue's sidechain
    do j=1, group(res_i)%cnum2
        natom=natom+1
        epsion_pep(natom, res_i)=group_para(res_i)%epsion2(j)
        r_pep(natom, res_i)=group_para(res_i)%r2(j)
        charge_pep(natom, res_i)=group_para(res_i)%charge2(j)
        rborn_pep(natom, res_i)=group_para(res_i)%rborn2(j)
        fs_pep(natom, res_i)=group_para(res_i)%fs2(j)
        dielecons_pep(natom, res_i)=group_para(res_i)%dielecons2(j)
        atomid_pep(natom, res_i)=atom_offset+group_para(res_i)%atomid2(j)
    enddo
else !we are coarse graining the sidechain, so load parameters for the coarse-grained bead
    natom=natom+1
    epsion_pep(natom, res_i)=cg_epsilon 
    r_pep(natom, res_i)=cg_r 
    charge_pep(natom, res_i)= cg_q 
    rborn_pep(natom, res_i)=cg_rborn 
    fs_pep(natom, res_i)=cg_f 
    dielecons_pep(natom, res_i)=cg_dielec
    atomid_pep(natom, res_i)=atom_offset+group_para(res_i)%atomid2(1)
endif

do j=1, group(res_i)%cnum3
    natom=natom+1
    epsion_pep(natom, res_i)=group_para(res_i)%epsion3(j)
    r_pep(natom, res_i)=group_para(res_i)%r3(j)
    charge_pep(natom, res_i)=group_para(res_i)%charge3(j)
    rborn_pep(natom, res_i)=group_para(res_i)%rborn3(j)
    fs_pep(natom, res_i)=group_para(res_i)%fs3(j)
    dielecons_pep(natom, res_i)=group_para(res_i)%dielecons3(j)
    atomid_pep(natom, res_i)=atom_offset+group_para(res_i)%atomid3(j)
enddo

!update atom offset
atom_offset = atom_offset + natom

end subroutine load_params_single_cg

subroutine load_coords_single_cg(group, res_i, cg_flag)
implicit none
integer                         :: natom
integer                         :: res_i, j, cg_flag
type(groupdetails)              :: group(:)

natom = 0
do j=1, group(res_i)%cnum1  
    natom=natom+1
    mdcrd_pep(:, natom, res_i)=group(res_i)%coo1(:,j)
enddo

sc_starts(res_i) = natom+1
if (cg_flag.eq.0) then !store coordinates for all atoms in the sidechain
    do j=1, group(res_i)%cnum2
        natom=natom+1
        mdcrd_pep(:,natom, res_i)=group(res_i)%coo2(:,j)
    enddo
else !we are coarse-graining the sidechain, so only store coordinates of the beta-carbon
    do j=1, group(res_i)%cnum2
        if (group(res_i)%atype2(j).eq."CB") then
            natom=natom+1
            mdcrd_pep(:,natom, res_i)=group(res_i)%coo2(:,j)
            exit
        endif
    enddo
endif
sc_ends(res_i) = natom

do j=1, group(res_i)%cnum3
    natom=natom+1
    mdcrd_pep(:,natom, res_i)=group(res_i)%coo3(:,j)
enddo

res_ends(res_i)=natom

end subroutine load_coords_single_cg

subroutine load_peptide_info_cg(group_pep, para_pep, cg_flags)
implicit none
integer                         :: res_i, cg_flags(pep_res), atom_offset
type(groupdetails)              :: group_pep(pep_res)
type(energyparameters)          :: para_pep(pep_res)

call load_res_cutoffs(group_pep, cg_flags)
atom_offset = 0
do res_i=1, pep_res
    if (res_coords_flags(res_i).eq.1) then 
        call load_coords_single_cg(group_pep, res_i, cg_flags(res_i))
        res_coords_flags(res_i) = 1 !if we are reloading coordinates, then also need to recalculate pairwise distances. Do this in a separate loop to ensure the correct coordinates are used
    endif
    if (res_params_flags(res_i).eq.1) then
        call load_params_single_cg(group_pep,para_pep, atom_offset, res_i, cg_flags(res_i))
    endif
enddo
atom_num_pep = res_ends(pep_res)

do res_i=1, pep_res
    if (res_coords_flags(res_i).eq.1) then 
        call calc_pairwise_distances(res_i)
        res_coords_flags(res_i) = 0
    endif
enddo

if (atom_links_flag.eq.1) then
    call atom_links(group_pep, cg_flags)
    atom_links_flag = 0
endif

end subroutine load_peptide_info_cg

subroutine load_res_cutoffs(group_pep, cg_flags)
implicit none
integer                         :: natom, res_i, cg_flags(:)
type(groupdetails)              :: group_pep(:)

natom = 0
do res_i=1, pep_res
    res_starts(res_i) = natom+1
    sc_starts(res_i) = res_starts(res_i) + group_pep(res_i)%cnum1
    if (cg_flags(res_i).eq.0) then !residue is not being coarse-grained, so include all atoms in sidechain
        sc_ends(res_i) = sc_starts(res_i) + group_pep(res_i)%cnum2-1
    else !we are coarse-graining sidechain, so only put one atom in sidechain (the coarse-grained bead)
        sc_ends(res_i) = sc_starts(res_i) 
    endif
    res_ends(res_i) = sc_ends(res_i) + group_pep(res_i)%cnum3
    atoms_per_res(res_i) = res_ends(res_i) - res_starts(res_i) + 1
    natom = res_ends(res_i)
enddo


end subroutine load_res_cutoffs

subroutine calc_pairwise_distances(res_i)
!!NOTE: option to optimize the distance calculation is to use a matrix indicating which pairs of residues require distance calculations. Current method can sometimes count distances twice
implicit none
integer         :: res_i, atom_i, atom_j, res_j, i_index, j_index
real            :: rij, i_loc(3)

do atom_i=1, atoms_per_res(res_i)
    i_index = atom_i
    i_loc = mdcrd_pep(:,atom_i, res_i)

    do res_j = 1, pep_res
        do atom_j=1, atoms_per_res(res_j)
            j_index = atom_j

            !get intermolecular distance
            call DIST_SQUARE(i_loc, mdcrd_pep(:,atom_j, res_j), rij)
            rij_matrix_pep(j_index, i_index, res_j, res_i) = rij 
            rij_matrix_pep(i_index, j_index, res_i, res_j) = rij
        enddo
    enddo

    !get peptide-receptor distances
    do atom_j=1, atom_num_rec
        call DIST_SQUARE(i_loc, mdcrd_rec(:,atom_j), rij)
        rij_matrix_sys(atom_j, i_index, res_i) = rij 
    enddo
enddo
end subroutine calc_pairwise_distances

subroutine calc_pairwise_distances_iso_rec
implicit none
integer         :: i_index, j_index
real            :: i_loc(3), rij

do i_index=1, atom_num_rec
    i_loc = mdcrd_rec(:,i_index)
    do j_index=i_index+1, atom_num_rec
        !get intermolecular distances
        call DIST_SQUARE(i_loc, mdcrd_rec(:,j_index), rij)
        rij_matrix_rec(j_index, i_index) = rij 
        rij_matrix_rec(i_index, j_index) = rij
    enddo
enddo
end subroutine calc_pairwise_distances_iso_rec

subroutine atom_links(group_pep, cg_flags)
implicit none
integer              :: status, linkindex(4), anum, resnum, z
integer              :: linknum, cg_flags(:)
integer              :: res_i, res_j, atom_i, atom_j, link_i, link_j, link_i_2, atom_i_ex, atom_i_ex2
type(groupdetails)   :: group_pep(:)

!first read in the atom_links info for each residue in the peptide
do res_i=1, pep_res !only need to create atom_links for peptide (atom_links control whether van der Waals or Coloumbic interactions measured, and these are never calculated for receptor-receptor atom pairs)
    if (cg_flags(res_i).eq.0) then
        open(10, file=trim(lib_path)//'/Atomlink/'//trim(group_pep(res_i)%gtype), status="old")
    else
        open(10, file=trim(lib_path)//'/Atomlink_CG/'//trim(group_pep(res_i)%gtype), status="old")
    endif
    do while(.true.)
        read(10, 20, iostat=status) atom_i, linknum, (linkindex(z), z=1, linknum)
        if(status.ne.0) goto 30
        atom(atom_i, res_i)%linknum=linknum
        do link_i=1, linknum
            anum=linkindex(link_i) 
            if (anum.lt.0) then
                resnum = res_i-1
                anum = atoms_per_res(res_i-1)
            elseif (anum.gt.atoms_per_res(res_i)) then
                resnum = res_i+1
                anum = 1
            else
                resnum = res_i
            endif
            atom(atom_i, res_i)%linkindex(:,link_i)=(/anum, resnum/)
        enddo
    enddo
30      continue
    close(10)
enddo
20  format(i6, i7, 4i3)

do res_i=1, pep_res !go through peptide residues
    do atom_i=1, atoms_per_res(res_i) !go through all atoms in the current peptide

        !find atoms that are only one bond away from current atom, and to numex (excluded list for VDW and ELE interactions)
        numex_pep(atom_i, res_i)=0
        do link_i=1, atom(atom_i, res_i)%linknum
            numex_pep(atom_i, res_i)=numex_pep(atom_i, res_i)+1
            inb_pep(:, numex_pep(atom_i, res_i), atom_i, res_i)=atom(atom_i, res_i)%linkindex(:, link_i)
        enddo

        !find atoms that are two bonds away from current atom, and to numex (excluded list for VDW and ELE interactions)
        do link_i=1, atom(atom_i, res_i)%linknum
            atom_j=atom(atom_i, res_i)%linkindex(1, link_i)
            res_j=atom(atom_i, res_i)%linkindex(2, link_i)
            do link_j=1, atom(atom_j, res_j)%linknum
                if(atom(atom_j, res_j)%linkindex(1, link_j).eq.atom_i .and. atom(atom_j, res_j)%linkindex(2, link_j).eq.res_i) goto 40 !we went back to atom_i, so don't bother evaluating its neighbors
                do link_i_2=1, atom(atom_i, res_i)%linknum !make sure this isn't an atom that both atom_i and atom_j are directly bonded to
                    if(atom(atom_j, res_j)%linkindex(1,link_j).eq.atom(atom_i, res_j)%linkindex(1,link_i_2) .and. &
                    atom(atom_j, res_j)%linkindex(2,link_j).eq.atom(atom_i, res_j)%linkindex(2,link_i_2)) goto 40
                enddo
                numex_pep(atom_i, res_i)=numex_pep(atom_i, res_i)+1
                inb_pep(:, numex_pep(atom_i, res_i), atom_i, res_i)=atom(atom_j, res_j)%linkindex(:, link_j)
                40 continue
            enddo
        enddo
    enddo
enddo

!find atoms that are three bondsaway from current atom, and to numex4 (1-4 list for VDW and ELE interactions)
do res_i=1, pep_res !go through peptide residues
    do atom_i=1, atoms_per_res(res_i) !go through all atoms in the current peptide
        numex4_pep(atom_i, res_i)=0
        do atom_i_ex=1, numex_pep(atom_i, res_i)
            atom_j=inb_pep(1, atom_i_ex, atom_i, res_i)
            res_j=inb_pep(2, atom_i_ex, atom_i, res_i)
            do link_j=1, atom(atom_j, res_j)%linknum
                if(atom(atom_j, res_j)%linkindex(1, link_j).eq.atom_i .and. &
                atom(atom_j, res_j)%linkindex(2, link_j).eq.res_i) goto 50
                do atom_i_ex2=1, numex_pep(atom_i, res_i)
                    if(atom(atom_j, res_j)%linkindex(1,link_j).eq.inb_pep(1,atom_i_ex2, atom_i,res_i) .and. &
                    atom(atom_j, res_j)%linkindex(2,link_j).eq.inb_pep(2,atom_i_ex2, atom_i,res_i)) goto 50
                enddo
                numex4_pep(atom_i, res_i)=numex4_pep(atom_i, res_i)+1
                inb4_pep(:,numex4_pep(atom_i, res_i), atom_i, res_i)=atom(atom_j, res_j)%linkindex(:,link_j)
            50  continue
            enddo
        enddo
    enddo
enddo

return
end subroutine atom_links

!subroutine atom_links4sidechain(res_i, group_pep, numex_pep, inb_pep, numex4_pep, inb4_pep)
!implicit none
!integer                    :: j, status, i1, j1, j2, k, atomid, natom, res_i
!integer                    :: linkindex(4)
!integer                    :: id, linknum, i2
!integer                    :: inb_pep(20,60), inb4_pep(60, 60), numex_pep(60), numex4_pep(60)
!type(groupdetails)         :: group_pep(:)
!type(atomlink)             :: atom(60)
!
!!first read sidechain atom connectivity from file
!open(10, file=trim(lib_path)//'/Atomlink/'//trim(group_pep(res_i)%gtype), status="old")
!do while(.true.)
!    read(10, 20, iostat=status) id, linknum, (linkindex(j), j=1, linknum)
!    if(status.ne.0) goto 30
!    atomid=id
!    atom(atomid)%linknum=linknum
!    do j=1, linknum
!        atom(atomid)%linkindex(j)=linkindex(j)
!    enddo
!enddo
!30   continue
!close(10)
!natom=atomid
!20  format(i6, i7, 4i3)
!
!!Next, add atoms one bond away from each other to exclusion list (don't calculate VDW or ELE energy for these pairs)
!do i1=1, natom
!    numex_pep(i1)=0
!    do j1=1, atom(i1)%linknum
!        numex_pep(i1)=numex_pep(i1)+1
!        inb_pep(numex_pep(i1), i1)=atom(i1)%linkindex(j1)
!    enddo
!
!    !Next, add atoms two bonds from each to exclusion list (don't calculate VDW or ELE energy for these pairs)
!    do j1=1, atom(i1)%linknum
!        i2=atom(i1)%linkindex(j1)
!        if((i2.ge.1).and.(i2.le.natom)) then
!            do j2=1, atom(i2)%linknum
!                if(atom(i2)%linkindex(j2).eq.i1) goto 40
!                do k=1, atom(i1)%linknum
!                    if(atom(i2)%linkindex(j2).eq.atom(i1)%linkindex(k)) goto 40
!                enddo
!                numex_pep(i1)=numex_pep(i1)+1
!                inb_pep(numex_pep(i1),i1)=atom(i2)%linkindex(j2)
!                40  continue
!            enddo
!        endif
!    enddo
!enddo
!
!!Next, add atoms three bonds from each to 1-4 list (rescale the VDW or ELE energy for these pairs)
!do i1=1, natom
!    numex4_pep(i1)=0
!    do j1=1, numex_pep(i1)
!        i2=inb_pep(j1,i1)
!        if((i2.ge.1).and.(i2.le.natom)) then
!            do j2=1, atom(i2)%linknum
!                if(atom(i2)%linkindex(j2).eq.i1) goto 50
!                do k=1, numex_pep(i1)
!                    if(atom(i2)%linkindex(j2).eq.inb_pep(k,i1)) goto 50
!                enddo
!                numex4_pep(i1)=numex4_pep(i1)+1
!                inb4_pep(numex4_pep(i1), i1)=atom(i2)%linkindex(j2)
!                50 continue
!            enddo
!        endif
!    enddo
!enddo

!return
!end subroutine atom_links4sidechain

subroutine res_atom_indices(group, res_i, res_i_start, res_i_end, res_i_sc_start, res_i_sc_end)
implicit none
integer                    :: res_i_start, res_i_end, res_i_sc_start, res_i_sc_end
integer                    :: natom, i, res_i
type(groupdetails)         :: group(:)

natom=0
do i=1, res_i-1
    natom = natom + group(i)%cnum1+group(i)%cnum2+group(i)%cnum3
enddo
res_i_start = natom+1
res_i_sc_start = res_i_start + group(res_i)%cnum1
res_i_sc_end = res_i_sc_start + group(res_i)%cnum2-1
res_i_end = res_i_sc_end + 1 + group(res_i)%cnum3

end subroutine res_atom_indices

subroutine update_rotamer_coords(aa_group, group, res_i)
implicit none
integer                           :: res_i, i
type(groupdetails)                :: group(:), aa_group

do i=1, group(res_i)%cnum2
    aa_group%coo2(:,i) = group(res_i)%coo2(:,i)
enddo

end subroutine update_rotamer_coords

subroutine sidechain_rotation_rodrig(num_sc_angles, sc_groups, dihedrals, sc_coords, CA, group_pep, monitor, index, loc)
!this subroutine uses the Rodrigues rotation formula to rotate sidechain angles by values in "dihedrals". Formula/variables shown here: https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
integer                         :: num_sc_angles, sc_groups(6), monitor(6), i, j, l, loc
real                            :: dihedrals(:), sin, cos, CA(3)
real                            :: k(3), k_cross(3), k_dot, v(3), origin(3)
type(conformer4sidechain)       :: sc_coords(6)
type(groupdetails)              :: group_pep(:)
type(index4sidechain)           :: index(60)

do i=1, num_sc_angles
    cos=cosd(dihedrals(i)); sin=sind(dihedrals(i))
    if(i==1) then
        origin = CA
    else
        origin = sc_coords(i-1)%member(:,monitor(i-1))
    endif
    k = sc_coords(i)%member(:,monitor(i)) - origin
    call normalize_vector(k)

    do j=(i+1), (num_sc_angles+1)
        do l=1, sc_groups(j)
            v = sc_coords(j)%member(:,l) - origin
            call cross_product(k, v, k_cross)
            k_dot =  DOT_PRODUCT(k,v)
            sc_coords(j)%member(:,l) = cos*v + sin*k_cross +(1-cos)*k_dot*k
        enddo
    enddo
enddo

do j=1, group_pep(loc)%cnum2
    group_pep(loc)%coo2(:,j)=sc_coords(index(j)%class_No)%member(:,index(j)%member_No)
enddo

end subroutine sidechain_rotation_rodrig

subroutine sidechain_rotation_4_gradient_rodrig(num_sc_angles, sc_angle, sc_groups, theta, sc_coords, CA, group_pep, monitor, index, loc)
!this subroutine uses the Rodrigues rotation formula to rotate sidechain angles by values in "dihedrals". Formula/variables shown here: https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
!This subroutine differs from "sidechain_rotation" subroutine in that it only performs 1 sidechain rotation, whereas "sidechain_rotation" rotates all rotable bonds in sidechain
implicit none       
integer                         :: num_sc_angles, sc_angle, sc_groups(6), monitor(6), i, j, loc
real                            :: CA(3), sin, cos, theta
real                            :: k(3), k_cross(3), k_dot, v(3), origin(3)
type(conformer4sidechain)       :: sc_coords(6)
type(index4sidechain)           :: index(60)
type(groupdetails)              :: group_pep(:)


cos = cosd(theta)
sin = sind(theta)

!first get rotation axis
if(sc_angle==1) then !special case since one reference point is alpha carbon
    origin = CA
else !more general case where both atoms are stored in "monitor"
    origin = sc_coords(sc_angle-1)%member(:,monitor(sc_angle-1))
endif
k=sc_coords(sc_angle)%member(:,monitor(sc_angle)) - origin
call normalize_vector(k)

!rotate all atoms affected by the rotation (i.e. all atoms farther away from peptide backbone)
do i=(sc_angle+1), (num_sc_angles+1)
    do j=1, sc_groups(i)
        v = sc_coords(i)%member(:,j) - origin
        call cross_product(k, v, k_cross)
        k_dot =  DOT_PRODUCT(k,v)
        sc_coords(i)%member(:,j) = cos*v + sin*k_cross + (1-cos)*k_dot*k
    enddo
enddo

do i=1, group_pep(loc)%cnum2
    group_pep(loc)%coo2(:,i)=sc_coords(index(i)%class_No)%member(:,index(i)%member_No) !put atom positions into group_pep
enddo

end subroutine sidechain_rotation_4_gradient_rodrig

subroutine sidechain_rotation(num_sc_angles, sc_groups, dihedrals, sc_coords, CA, group_pep, monitor, index, loc)
implicit none       

integer                         :: num_sc_angles, sc_groups(6), monitor(6), i, j, k, loc
real                            :: dihedrals(:), sin_angle, cos_angle, CA(3), m(3,3), rotaxis(3)
type(conformer4sidechain)       :: sc_coords(6)
type(groupdetails)              :: group_pep(:)
type(index4sidechain)           :: index(60)

do i=1, num_sc_angles
    cos_angle=cosd(dihedrals(i)); sin_angle=sind(dihedrals(i))
    if(i==1) then
        rotaxis=sc_coords(i)%member(:,monitor(i))-CA
    else
        rotaxis=sc_coords(i)%member(:,monitor(i))-sc_coords(i-1)%member(:,monitor(i-1))
    endif

    !normalize the rotation axis
    call normalize_vector(rotaxis)

    call axisrotation(rotaxis, cos_angle, sin_angle, m)

    do j=(i+1), (num_sc_angles+1)
        do k=1, sc_groups(j)
            sc_coords(j)%member(:,k)=matmul(m,sc_coords(j)%member(:,k)-sc_coords(i)%member(:,monitor(i)))+sc_coords(i)%member(:,monitor(i))
        enddo
    enddo
enddo

do j=1, group_pep(loc)%cnum2
    group_pep(loc)%coo2(:,j)=sc_coords(index(j)%class_No)%member(:,index(j)%member_No)
enddo

end subroutine sidechain_rotation

subroutine sidechain_rotation_4_gradient(num_sc_angles, sc_angle, sc_groups, theta, sc_coords, CA, group_pep, monitor, index, loc)
implicit none       
integer                         :: num_sc_angles, sc_angle, sc_groups(6), monitor(6), j, k, loc
real                            :: theta, cos_angle, sin_angle, CA(3), rotaxis(3), m(3,3)
type(conformer4sidechain)       :: sc_coords(6)
type(index4sidechain)           :: index(60)
type(groupdetails)              :: group_pep(:)

!first get rotation axis
if(sc_angle==1) then !special case since one reference point is alpha carbon
    rotaxis(:)=sc_coords(sc_angle)%member(:,monitor(sc_angle))-CA
else !more general case where both atoms are stored in "monitor"
    rotaxis(:)=sc_coords(sc_angle)%member(:,monitor(sc_angle))-sc_coords(sc_angle-1)%member(:,monitor(sc_angle-1))
endif

!normalize the rotation axis
call normalize_vector(rotaxis)

cos_angle=cosd(theta); sin_angle=sind(theta)
call axisrotation(rotaxis, cos_angle, sin_angle, m) !get matrix that corresponds to rotation about rotation axis by delta_chi

!rotate all atoms affected by the rotation (i.e. all atoms farther away from peptide backbone)
do j=(sc_angle+1), (num_sc_angles+1)
    do k=1, sc_groups(j)
        sc_coords(j)%member(:,k)=matmul(m, sc_coords(j)%member(:,k)-sc_coords(sc_angle)%member(:,monitor(sc_angle)))+sc_coords(sc_angle)%member(:,monitor(sc_angle)) 
    enddo
enddo

do j=1, group_pep(loc)%cnum2
    group_pep(loc)%coo2(:,j)=sc_coords(index(j)%class_No)%member(:,index(j)%member_No) !put atom positions into group_pep
enddo

end subroutine sidechain_rotation_4_gradient

subroutine random_rotation_matrix(dir, rot_m)
!!!
! this subroutine performs a random rotatation for peptide in either x, y, or z direction. Receptor does not move.
!!!
implicit none
integer                         :: dir
real                            :: theta, s, c, rot_m(3,3), m(3,3)

!get random rotation angle
call RANDOM_NUMBER(theta)
theta = PI*(theta*2 - 1)

s = sin(theta); c = cos(theta)

!get rotation matrix
if (dir.eq.1) then
m = reshape((/ 1.,   0.,     0., &
            0.,     c,          s, &
            0.,    -s,            c/), shape(m))
elseif(dir.eq.2) then
m = reshape((/ s,      0.,   0.,      c, &
            0.,  1.,   0., &
            c,        0.,     -s/), shape(m))
else
m = reshape((/ c    ,   s,      0., &
               -s,        c,      0., &
              0.,  0.,  1./), shape(m))
endif

!multiply current rotation matrix by additional matrix
rot_m = MATMUL(m, rot_m)
end subroutine random_rotation_matrix

subroutine position_peptide(group_pep, coords, moments)
implicit none

real                        :: coords(6), delta_z
real, parameter             :: PI = 3.141592653589793238462643383, abserr = 0.00001
real                        :: rot_m(3,3), obj_axes(3,3), shift(3), r, moments(3)
integer                     :: i
type(groupdetails)          :: group_pep(pep_res)

res_coords_flags = 1

!!! STEP ONE: get random orientation and location of peptide over polymer surface

!calculate peptide's center of mass
call calc_com(group_pep, pep_com, pep_res)

!get matrix corresponding to three random rotations about x, y, then z axis
rot_m= reshape((/1,0,0, &
        0,1,0, &
        0,0,1/), shape(rot_m))
do i=1,3
    call random_rotation_matrix(i, rot_m)
enddo

!rotate the peptide
call rotate_peptide(group_pep, rot_m, pep_com)

!calculate vector to shift peptide to system origin
shift(1:2) = -pep_com(1:2)
shift(3) = rec_max_z - pep_com(3) !note this places the peptide directly on the surface. The extra shift above the polymer surface to z_min is done in the do loop below

!add a random additional translation to the "shift" vector
do i=1,3
    call RANDOM_NUMBER(r)
    coords(i) = r
    shift(i) = shift(i) + limits(1,i) + r*(limits(2,i)-limits(1,i))
enddo

!determine if additional offset in z-direction needed to avoid steric clash
call min_peptide_z(group_pep)
delta_z = (pep_min_z + shift(3)) - (rec_max_z + z_min)  !if delta_z is less than or equal to 0, then the distance between peptide and receptor is below the allowed limit and an extra shift is needed!
if (delta_z .lt. 0) then
    shift(3) = shift(3) - delta_z
    coords(3) = coords(3) - delta_z/(limits(2,3)-limits(1,3))
endif 

!translate peptide by the shift vector
call translate_group(group_pep, shift, pep_res)
pep_com = pep_com+shift

!!! STEP TWO: Get coordinates of peptide in terms of (x,y,z) location and Euler angles
!!! NOTE: I use radius of gyration axis and x,y,z axis of system to define the Euler angles
call calc_Euler_angles(group_pep, coords(4:6), obj_axes, moments)

return
end subroutine position_peptide

subroutine reposition_peptide(group_pep, phis, psis, loc, ori, moments)
implicit none
integer                         :: res_i, cutoff
real                            :: loc(3), ori(3), moments(3)
real                            :: phis(pep_res), psis(pep_res), phis_init(pep_res), psis_init(pep_res)
real                            :: zshift1, zshift2
real                            :: shift(3), obj_axes(3,3)
type(groupdetails)              :: group_pep(:)

!     .. Intrinsic Functions ..
INTRINSIC                       :: INT, ABS, SQRT

!debug variables
real                            :: angle_cur(3)

res_coords_flags = 1

!!STEP ONE: put the desired dihedrals into the peptide
!calculate the current peptide dihedrals
call calc_backbone_dihedrals(group_pep, phis_init, psis_init)

!insert the target phi and psis into the peptide
cutoff = pep_res/2
do res_i=1, pep_res
    call insert_dihedral(group_pep, res_i, phis(res_i), psis(res_i), phis_init(res_i), psis_init(res_i), cutoff)
enddo

!confirm the dihedrals are correct, if debugging
if (debug_flag.eq.1) call calc_backbone_dihedrals(group_pep, phis_init, psis_init)

!!STEP TWO: put the peptide in the desired orientation 
!now rotate peptide so Euler angles match those in input (following Euler angle definitions on wikipedia page here: https://en.wikipedia.org/wiki/Euler_angles)
call align_euler_angles(group_pep, obj_axes, ori)

!if debugging, then verify the peptide's Euler angles now match the target values
if (debug_flag.eq.1) then
    write (*,*) "***Checking Euler angle rotation correct***"
    call calc_Euler_angles(group_pep, angle_cur, obj_axes, moments)
    if (abs(angle_cur(1)-ori(1)) > 0.1 .or. abs(angle_cur(2)-ori(2)) > 0.1 .or. abs(angle_cur(3)-ori(3)) > 0.1) then
        write (*,*) "Rotation to target Euler angles not done properly!"
    else
        write (*,*) "Rotation to target Euler angles successful!"
    endif
endif

!!STEP THREE: translate peptide to the desired position
shift(1) = -pep_com(1) + limits(1,1) + ranges(1)*loc(1)
shift(2) = -pep_com(2) + limits(1,2) + ranges(2)*loc(2)

call min_peptide_z_backbone(group_pep)
zshift1 = (rec_max_z - pep_com(3)) + limits(1,3) + ranges(3)*loc(3)
zshift2 = rec_max_z + limits(1,3) - (pep_min_z + zshift1)
if (zshift2.gt.0) zshift1 = zshift1 + zshift2
shift(3) = zshift1
call translate_group(group_pep, shift, pep_res)

return
end subroutine reposition_peptide

subroutine reposition_peptide_gridsearch(group_pep, group_rec, loc, ori)   !, phis, psis, ori, moments)
implicit none
integer                         :: cutoff, res_i 
real                            :: phis_init(pep_res), psis_init(pep_res)
real                            :: shift(3), obj_axes(3,3), loc(3), ori(3)
type(groupdetails)              :: group_pep(:), group_rec(:)

!     .. Intrinsic Functions ..
INTRINSIC                       :: INT, ABS, SQRT

!debug variables
real                            :: angle_cur(3), moments(3)

res_coords_flags = 1

!!STEP ONE: put the desired dihedrals into the peptide, 
!calculate the current peptide dihedrals
call calc_backbone_dihedrals(group_pep, phis_init, psis_init)

!insert the target phi and psis into the peptide
cutoff = pep_res/2
do res_i=1, pep_res
    if (phis(res_i).gt.-200 .and. psis(res_i).gt.-200) then
        call insert_dihedral(group_pep, res_i, phis(res_i), psis(res_i), phis_init(res_i), psis_init(res_i), cutoff)
    endif
enddo

!confirm the dihedrals are correct, if debugging
if (debug_flag.eq.1) call calc_backbone_dihedrals(group_pep, phis_init, psis_init)

!!STEP TWO: put the peptide in the desired orientation, if an orientation is specified
! Euler angle definitions follow those on wikipedia page here: https://en.wikipedia.org/wiki/Euler_angles)
if (alpha.ge.0 .and. beta.ge.0 .and. gamma.ge.0) then
    call align_euler_angles(group_pep, obj_axes, ori)

    !if debugging, then verify the peptide's Euler angles now match the target values
    if (debug_flag.eq.1) then
        write (*,*) "***Checking Euler angle rotation correct***"
        call calc_Euler_angles(group_pep, angle_cur, obj_axes, moments)
        if (abs(angle_cur(1)-ori(1)) > 0.1 .or. abs(angle_cur(2)-ori(2)) > 0.1 .or. abs(angle_cur(3)-ori(3)) > 0.1) then
            write (*,*) "Rotation to target Euler angles not done properly!"
        else
            write (*,*) "Rotation to target Euler angles successful!"
        endif
    endif
endif

!!STEP THREE: translate peptide to the desired position
call max_receptor_z(group_rec)
call calc_com(group_pep, pep_com, pep_res)
shift(1) = loc(1)
shift(2) = loc(2)
shift(3) = -pep_com(3) + rec_max_z  + loc(3)
call translate_group(group_pep, shift, pep_res)

if (debug_flag.eq.1) call calc_com(group_pep, pep_com, pep_res)

return
end subroutine reposition_peptide_gridsearch

subroutine align_euler_angles(group_pep, obj_axes, ori)
implicit none
real                            :: rot_m(3,3), angles_cur(3), ori(3), obj_axes(3,3)
real                            :: moments(3)
type(groupdetails)              :: group_pep(:)

INTRINSIC       :: TRANSPOSE, MATMUL

!debug variables
!real                            :: z_axis(3) = (/0.,0.,1./), x_axis(3) = (/1.,0.,0./),y_axis(3) = (/0.,1.,0./), alt_m(3,3)

!get decomposition of moment of inertia tensor
call decompose_mom_tensor(group_pep, obj_axes, moments)

!if debugging, do an initial rotation to align peptide axes with system axes
if (debug_flag .eq.1) then
    rot_m = TRANSPOSE(obj_axes)
    call rotate_peptide(group_pep, rot_m, pep_com)
    call calc_Euler_angles(group_pep, angles_cur, obj_axes, moments)
endif

!calculate rotation matrix to align Euler angles
call Euler_Rotation_Matrix(ori, rot_m)

!do the rotation
if (debug_flag.eq.1) then
    call rotate_peptide(group_pep, rot_m, pep_com)
    !call calc_Euler_angles(group_pep, angles_cur, obj_axes, moments)
else
    rot_m = MATMUL(rot_m, TRANSPOSE(obj_axes)) 
    call rotate_peptide(group_pep, rot_m, pep_com)
endif

end subroutine align_euler_angles

subroutine insert_dihedral(group_pep, res, phi, psi, phi_init, psi_init, cutoff)
implicit none
integer                     :: res, cutoff
real                        :: phi, psi, phi_init, psi_init, delta_phi, delta_psi
type(groupdetails)          :: group_pep(pep_res)

!debug variables
real                        :: phis_new(pep_res), psis_new(pep_res)

if (debug_flag.eq.1) then
    call calc_backbone_dihedrals(group_pep, phis_new, psis_new)
endif

delta_phi = phi - phi_init
delta_psi = psi - psi_init
if (res .lt. cutoff) then
    call backbone_rotation_Nterm(group_pep, res, delta_phi, delta_psi)
else
    call backbone_rotation_Cterm(group_pep, res, delta_phi, delta_psi)
endif

if (debug_flag.eq.1) then
    call calc_backbone_dihedrals(group_pep, phis_new, psis_new)
endif

end subroutine insert_dihedral

subroutine backbone_rotation_Cterm(group_pep, res, delta_phi, delta_psi)
implicit none
integer                             :: res, j, k, l, temp_num
real                                :: N(3), CA(3), C(3)
real                                :: rotaxis(3), m(3,3)
real                                :: delta_phi, delta_psi
real                                :: cos_angle, sin_angle
type(groupdetails)                  :: group_pep(:)

res_coords_flags(res:pep_res) = 1 !tell system to reload peptide coordinates into array format prior to calculating energies
!first rotate by phi
cos_angle=COS(deg_2_rad*delta_phi)
sin_angle=SIN(deg_2_rad*delta_phi)

!get N and CA positions
do k=1, group_pep(res)%cnum1
    if(group_pep(res)%atype1(k)=="N") then
        N=group_pep(res)%coo1(:,k)
    elseif(group_pep(res)%atype1(k)=="CA") then
        CA=group_pep(res)%coo1(:,k)
    endif
enddo

!get rotation matrix about the CA-N bond (the rotation axis)
rotaxis=N-CA    
call normalize_vector(rotaxis)
call axisrotation(rotaxis, cos_angle, sin_angle, m)

!move all atoms to reference frame where N is at the origin
temp_num=0
do j=1, group_pep(res)%cnum1
    if(group_pep(res)%atype1(j)=="HA".or.group_pep(res)%atype1(j)=="HA2" &
    .or.group_pep(res)%atype1(j)=="HA3") then
        group_pep(res)%coo1(:,j)=matmul(m,group_pep(res)%coo2(:,j) - N) + N
    endif
enddo

do j=1, group_pep(res)%cnum2
    group_pep(res)%coo2(:,j)=matmul(m,group_pep(res)%coo2(:,j) - N) + N
enddo
do j=1, group_pep(res)%cnum3
    group_pep(res)%coo3(:,j)=matmul(m,group_pep(res)%coo3(:,j) - N) + N
enddo

!now move all atoms for residues up to the C-terminus
do j=res+1, pep_res
    !move atoms to reference frame where N is at origin
    do l=1, group_pep(j)%cnum1
        group_pep(j)%coo1(:,l)=matmul(m,group_pep(j)%coo1(:,l)-N)+N
    enddo
    do l=1, group_pep(j)%cnum2
        group_pep(j)%coo2(:,l)=matmul(m,group_pep(j)%coo2(:,l)-N)+N
    enddo
    do l=1, group_pep(j)%cnum3
        group_pep(j)%coo3(:,l)=matmul(m,group_pep(j)%coo3(:,l)-N)+N
    enddo
enddo

!now repeat for the psi angle
!exit if we are at the C-terminus (no psi defined here)
if (res .eq. pep_res) return

!get sine and cosine values of delta_psi
cos_angle=COS(deg_2_rad*delta_psi)
sin_angle=SIN(deg_2_rad*delta_psi)

!get location of C coordinate
do k=1, group_pep(res)%cnum3
    if(group_pep(res)%atype3(k)=="C") then
        C=group_pep(res)%coo3(:,k)
    endif
enddo

!get rotation matrix about C-CA bond (the axis of rotation)
rotaxis=CA-C
call normalize_vector(rotaxis)
call axisrotation(rotaxis, cos_angle, sin_angle, m)

!perform rotation for group3 atoms in current residue
do j=1, group_pep(res)%cnum3
    group_pep(res)%coo3(:,j)=matmul(m,group_pep(res)%coo3(:,j)-CA)+CA
enddo

!do rotation for all residues, moving towards C-terminus
do j=res+1, pep_res
    !move atoms to reference frame where CA is at origin
    do l=1, group_pep(j)%cnum1
        group_pep(j)%coo1(:,l)=matmul(m,group_pep(j)%coo1(:,l)-CA)+CA
    enddo
    do l=1, group_pep(j)%cnum2
        group_pep(j)%coo2(:,l)=matmul(m,group_pep(j)%coo2(:,l)-CA)+CA
    enddo
    do l=1, group_pep(j)%cnum3
        group_pep(j)%coo3(:,l)=matmul(m,group_pep(j)%coo3(:,l)-CA)+CA
    enddo
enddo

return
end subroutine backbone_rotation_Cterm

subroutine backbone_rotation_Nterm(group_pep, res, delta_phi, delta_psi)
implicit none
integer                 :: res, j, k, l, temp_num
real                    :: N(3), CA(3), C(3)
real                    :: rotaxis(3), m(3,3)
real                    :: delta_phi, delta_psi
real                    :: cos_angle, sin_angle
type(groupdetails)      :: group_pep(pep_res)

res_coords_flags(1:res) = 1 !tell system to reload peptide coordinates into array format prior to calculating energies

!get cosine and sine values for delta_psi
cos_angle=COS(deg_2_rad*delta_psi) 
sin_angle=SIN(deg_2_rad*delta_psi)

!get C and CA positions, which is axis for doing psi rotation
do k=1, group_pep(res)%cnum3
    if(group_pep(res)%atype3(k)=="C") then
        C=group_pep(res)%coo3(:,k)
        exit
    endif
enddo
do k=1, group_pep(res)%cnum1
    if(group_pep(res)%atype1(k)=="CA") then
        CA=group_pep(res)%coo1(:,k)
        exit
    endif
enddo

!get rotation rotation axis about CA-C bond
rotaxis = C - CA
call normalize_vector(rotaxis)
call axisrotation(rotaxis, cos_angle, sin_angle, m)

!rotate atoms in current residue
!Don't need to do this for group 3 atoms, since they aren't being moved by this rotation
do j=1, group_pep(res)%cnum2
    group_pep(res)%coo2(:,j)=matmul(m,group_pep(res)%coo2(:,j)-C)+C
enddo
do j=1, group_pep(res)%cnum1
    group_pep(res)%coo1(:,j)=matmul(m,group_pep(res)%coo1(:,j)-C)+C
enddo

!repeat the rotation procedure for all residues going to N-terminus. Note that we need to do this for group 3 atoms now, since they WILL be rotated
do j=res-1, 1, -1
    do l=1, group_pep(j)%cnum3
        group_pep(j)%coo3(:,l)=matmul(m,group_pep(j)%coo3(:,l)-C)+C
    enddo
    do l=1, group_pep(j)%cnum2
        group_pep(j)%coo2(:,l)=matmul(m,group_pep(j)%coo2(:,l)-C)+C
    enddo
    do l=1, group_pep(j)%cnum1
        group_pep(j)%coo1(:,l)=matmul(m,group_pep(j)%coo1(:,l)-C)+C
    enddo
enddo

!now do the rotation for the phi change, if not at the N-terminus
if (res.eq.1) return 
cos_angle=COS(deg_2_rad*delta_phi)
sin_angle=SIN(deg_2_rad*delta_phi)

!get the location of the N atom
do k=1, group_pep(res)%cnum1
    if(group_pep(res)%atype1(k)=="N") then
        N = group_pep(res)%coo1(:,k)
        exit
    endif
enddo

!get the rotation axis about the N-CA bond
rotaxis = CA - N 
call normalize_vector(rotaxis)
call axisrotation(rotaxis, cos_angle, sin_angle, m)

!move group1 atoms on current residue to frame where CA carbon is at origin
temp_num=0
do j=1, group_pep(res)%cnum1
    if(group_pep(res)%atype1(j)=="N".or.group_pep(res)%atype1(j)=="H".or.group_pep(res)%atype1(j)=="H1".or. &
    group_pep(res)%atype1(j)=="H2".or.group_pep(res)%atype1(j)=="H3") then
        group_pep(res)%coo1(:,j) = matmul(m,group_pep(res)%coo1(:,j)-CA)+CA
    endif
enddo

!now repeat rotation for all atoms until N-terminus is reached
do j=res-1, 1, -1
    do l=1, group_pep(j)%cnum3
        group_pep(j)%coo3(:,l)=matmul(m,group_pep(j)%coo3(:,l)-CA)+CA
    enddo
    do l=1, group_pep(j)%cnum2
        group_pep(j)%coo2(:,l)=matmul(m,group_pep(j)%coo2(:,l)-CA)+CA
    enddo
    do l=1, group_pep(j)%cnum1
        group_pep(j)%coo1(:,l)=matmul(m,group_pep(j)%coo1(:,l)-CA)+CA
    enddo
enddo

return
end subroutine backbone_rotation_Nterm

subroutine backbone_rotation_Nterm_multiangle(group_pep, backbone_backup, ran_resi, delta_psi3, delta_phi3, &
    delta_psi2, delta_phi2, Tgroup, Tbackbone_backup)
implicit none
integer                             :: ran_resi(3), i, j, k, l, temp_num, ic
real                                :: ran2
real                                :: CA1(3), C1(3), N(3), CA(3), C(3)
real                                :: rotaxis(3), m(3,3)
real                                :: delta_psi3, delta_phi3, delta_psi2, delta_phi2, delta_psi1
real                                :: delta_psi(3), delta_phi(3), cos_angle, sin_angle
REAL                                :: coo(3,20)
character*4                         :: atype(20)
type(groupdetails)                  :: group_pep(pep_res), Tgroup(pep_res)
type(databackup)                    :: backbone_backup(pep_res), Tbackbone_backup(pep_res)


res_coords_flags(1:ran_resi(3)) = 1 !tell system to reload peptide coordinates into array format prior to calculating energies
delta_psi(3)=-delta_psi3; delta_psi(2)=-delta_psi2
delta_phi(3)=-delta_phi3; delta_phi(2)=-delta_phi2

!transfer data to temporary group
Tgroup=group_pep
Tbackbone_backup=backbone_backup

!do the following for the first two psi/phi pairs
do i=3, 2, -1
    !### first do the rotation for the psi angle ###
    cos_angle=cosd(delta_psi(i))
    sin_angle=sind(delta_psi(i))
    !get CA-C axis
    do k=1, Tgroup(ran_resi(i))%cnum3
        if(Tgroup(ran_resi(i))%atype3(k)=="C") then
            C=Tgroup(ran_resi(i))%coo3(:,k)
        endif
    enddo
    do k=1, Tgroup(ran_resi(i))%cnum1
        if(Tgroup(ran_resi(i))%atype1(k)=="CA") then
            CA=Tgroup(ran_resi(i))%coo1(:,k)
            exit
        endif
    enddo
    rotaxis=CA-C
    call normalize_vector(rotaxis)
    call axisrotation(rotaxis, cos_angle, sin_angle, m)

    !move current residue into reference frame where C is at the origin, rotate peptide, then go back to original reference frame
    do j=1, Tgroup(ran_resi(i))%cnum2
        Tgroup(ran_resi(i))%coo2(:,j)=matmul(m,Tgroup(ran_resi(i))%coo2(:,j)-C)+C
    enddo
    do j=1, Tgroup(ran_resi(i))%cnum1
        Tgroup(ran_resi(i))%coo1(:,j)=matmul(m,Tgroup(ran_resi(i))%coo1(:,j)-C)+C
    enddo
    Tbackbone_backup(ran_resi(i))%coo=matmul(m,Tbackbone_backup(ran_resi(i))%coo-C)+C

    !now do the rotation for all residues, moving towards first pivot in CONROT (same procedure as above)
    do j=ran_resi(i)-1, ran_resi(1), -1
        do l=1, Tgroup(j)%cnum3   
            Tgroup(j)%coo3(:,l)=matmul(m,Tgroup(j)%coo3(:,l)-C)+C
        enddo
        do l=1, Tgroup(j)%cnum2
            Tgroup(j)%coo2(:,l)=matmul(m,Tgroup(j)%coo2(:,l)-C)+C
        enddo
        do l=1, Tgroup(j)%cnum1
            Tgroup(j)%coo1(:,l)=matmul(m,Tgroup(j)%coo1(:,l)-C)+C
        enddo
        Tbackbone_backup(j)%coo=matmul(m,Tbackbone_backup(j)%coo-C)+C
    enddo

    !### now do the rotation for the phi angle ###
    cos_angle=cosd(delta_phi(i))
    sin_angle=sind(delta_phi(i))

    do k=1, Tgroup(ran_resi(i))%cnum1
        if(Tgroup(ran_resi(i))%atype1(k)=="N") then
            N=Tgroup(ran_resi(i))%coo1(:,k)
            exit
        endif
    enddo

    rotaxis=N-CA
    call normalize_vector(rotaxis)
    call axisrotation(rotaxis, cos_angle, sin_angle, m)

    temp_num=0
    do j=1, Tgroup(ran_resi(i))%cnum1
        if(Tgroup(ran_resi(i))%atype1(j)=="N".or.Tgroup(ran_resi(i))%atype1(j)=="H".or.& 
        Tgroup(ran_resi(i))%atype1(j)=="H1".or.Tgroup(ran_resi(i))%atype1(j)=="H2".or.&
        Tgroup(ran_resi(i))%atype1(j)=="H3") then
            temp_num=temp_num+1
            atype(temp_num)=Tgroup(ran_resi(i))%atype1(j)
            coo(:,temp_num)=Tgroup(ran_resi(i))%coo1(:,j)-CA
        endif
    enddo
    Tbackbone_backup(ran_resi(i))%coo=matmul(m,Tbackbone_backup(ran_resi(i))%coo-CA)+CA
    coo=matmul(m,coo) !don't add +CA here since we don't need all the terms in coo usually

    do l=1, temp_num
        do j=1, Tgroup(ran_resi(i))%cnum1
            if(atype(l)==Tgroup(ran_resi(i))%atype1(j)) then
                Tgroup(ran_resi(i))%coo1(:,j)=coo(:,l)+CA
            endif
        enddo
    enddo

    do j=ran_resi(i)-1, ran_resi(1), -1
        do l=1, Tgroup(j)%cnum3
            Tgroup(j)%coo3(:,l)=matmul(m,Tgroup(j)%coo3(:,l)-CA)+CA
        enddo
        do l=1, Tgroup(j)%cnum2
            Tgroup(j)%coo2(:,l)=matmul(m,Tgroup(j)%coo2(:,l)-CA)+CA
        enddo
        do l=1, Tgroup(j)%cnum1
            Tgroup(j)%coo1(:,l)=matmul(m,Tgroup(j)%coo1(:,l)-CA)+CA
        enddo
        Tbackbone_backup(j)%coo=matmul(m,Tbackbone_backup(j)%coo-CA)+CA
    enddo
enddo

do i=1, Tgroup(ran_resi(1))%cnum3
    if(Tgroup(ran_resi(1))%atype3(i)=="C") then
        C1=Tgroup(ran_resi(1))%coo3(:,i)
        exit
    endif
enddo

do i=1, Tgroup(ran_resi(1))%cnum1
    if(Tgroup(ran_resi(1))%atype1(i)=="CA") then
        CA1=Tgroup(ran_resi(1))%coo1(:,i)
        exit
    endif
enddo

call RANDOM_NUMBER(ran2)
ic=anint(ran2*180.)
delta_psi1=-(ic-90)
cos_angle=cosd(delta_psi1)
sin_angle=sind(delta_psi1)

rotaxis=CA1-C1
call normalize_vector(rotaxis)
call axisrotation(rotaxis, cos_angle, sin_angle, m)

do j=1, Tgroup(ran_resi(1))%cnum2
    Tgroup(ran_resi(1))%coo2(:,j)=matmul(m,Tgroup(ran_resi(1))%coo2(:,j)-C1)+C1
enddo
do j=1, Tgroup(ran_resi(1))%cnum1
    Tgroup(ran_resi(1))%coo1(:,j)=matmul(m,Tgroup(ran_resi(1))%coo1(:,j)-C1)+C1
enddo
Tbackbone_backup(ran_resi(1))%coo=matmul(m,Tbackbone_backup(ran_resi(1))%coo-C1)+C1

return
end subroutine backbone_rotation_Nterm_multiangle

subroutine backbone_rotation_Cterm_multiangle(group, backbone_backup, ran_resi, delta_phi1, delta_psi1,  &
delta_phi2, delta_psi2, Tgroup, Tbackbone_backup)
implicit none
integer                             :: ran_resi(3), i, j, k, l, temp_num, ic
real                                :: ran2
real                                :: N3(3), CA3(3), C3(3), N(3), CA(3), C(3)
real                                :: rotaxis(3), m(3,3)
real                                :: delta_phi1, delta_psi1, delta_phi2, delta_psi2, delta_phi3, delta_psi3
real                                :: delta_phi(2), delta_psi(2), cos_angle, sin_angle
real                                :: coo(3,20)
character*4                         :: atype(20)
type(groupdetails)                  :: group(pep_res), Tgroup(pep_res)
type(databackup)                    :: backbone_backup(pep_res), Tbackbone_backup(pep_res)


res_coords_flags(ran_resi(1):pep_res) = 1 !tell system to reload peptide coordinates into array format prior to calculating energies
delta_phi(1)=-delta_phi1; delta_phi(2)=-delta_phi2
delta_psi(1)=-delta_psi1; delta_psi(2)=-delta_psi2

Tgroup=group
Tbackbone_backup=backbone_backup

do i=1, 2
    cos_angle=cosd(delta_phi(i))
    sin_angle=sind(delta_phi(i))

    do k=1, Tgroup(ran_resi(i))%cnum1
        if(Tgroup(ran_resi(i))%atype1(k)=="N") then
            N=Tgroup(ran_resi(i))%coo1(:,k)
        elseif(Tgroup(ran_resi(i))%atype1(k)=="CA") then
            CA=Tgroup(ran_resi(i))%coo1(:,k)
        endif
    enddo

    rotaxis=CA-N
    call normalize_vector(rotaxis)
    call axisrotation(rotaxis, cos_angle, sin_angle, m)

    temp_num=0
    do j=1, Tgroup(ran_resi(i))%cnum1
        if(Tgroup(ran_resi(i))%atype1(j)=="HA".or.Tgroup(ran_resi(i))%atype1(j)=="HA2" &
        .or.Tgroup(ran_resi(i))%atype1(j)=="HA3") then
            temp_num=temp_num+1
            atype(temp_num)=Tgroup(ran_resi(i))%atype1(j)
            coo(:,temp_num)=Tgroup(ran_resi(i))%coo1(:,j)-N
        endif
    enddo

    coo=matmul(m,coo)
    do l=1, temp_num
        do j=1, Tgroup(ran_resi(i))%cnum1
            if(atype(l)==Tgroup(ran_resi(i))%atype1(j)) then
                Tgroup(ran_resi(i))%coo1(:,j)=coo(:,l)+N
            endif
        enddo
    enddo

    do j=1, Tgroup(ran_resi(i))%cnum2
        Tgroup(ran_resi(i))%coo2(:,j)=matmul(m,Tgroup(ran_resi(i))%coo2(:,j)-N)+N
    enddo
    do j=1, Tgroup(ran_resi(i))%cnum3
        Tgroup(ran_resi(i))%coo3(:,j)=matmul(m,Tgroup(ran_resi(i))%coo3(:,j)-N)+N
    enddo

    do j=ran_resi(i)+1, ran_resi(3)
        do l=1, Tgroup(j)%cnum1
            Tgroup(j)%coo1(:,l)=matmul(m,Tgroup(j)%coo1(:,l)-N)+N
        enddo
        do l=1, Tgroup(j)%cnum2
            Tgroup(j)%coo2(:,l)=matmul(m,Tgroup(j)%coo2(:,l)-N)+N
        enddo
        do l=1, Tgroup(j)%cnum3
            Tgroup(j)%coo3(:,l)=matmul(m,Tgroup(j)%coo3(:,l)-N)+N
        enddo
        Tbackbone_backup(j)%coo=matmul(m,Tbackbone_backup(j)%coo-N)+N
    enddo

    cos_angle=cosd(delta_psi(i))
    sin_angle=sind(delta_psi(i))

    do k=1, Tgroup(ran_resi(i))%cnum3
        if(Tgroup(ran_resi(i))%atype3(k)=="C") then
            C=Tgroup(ran_resi(i))%coo3(:,k)
        endif
    enddo

    rotaxis=C-CA
    call normalize_vector(rotaxis)
    call axisrotation(rotaxis, cos_angle, sin_angle, m)

    do j=1, Tgroup(ran_resi(i))%cnum3
        Tgroup(ran_resi(i))%coo3(:,J)=matmul(m,Tgroup(ran_resi(i))%coo3(:,j)-CA)+CA
    enddo

    do j=ran_resi(i)+1, ran_resi(3)
        do l=1, Tgroup(j)%cnum1
            Tgroup(j)%coo1(:,l)=matmul(m,Tgroup(j)%coo1(:,l)-CA)+CA
        enddo
        do l=1, Tgroup(j)%cnum2
            Tgroup(j)%coo2(:,l)=matmul(m,Tgroup(j)%coo2(:,l)-CA)+CA
        enddo
        do l=1, Tgroup(j)%cnum3
            Tgroup(j)%coo3(:,l)=matmul(m,Tgroup(j)%coo3(:,l)-CA)+CA
        enddo
        Tbackbone_backup(j)%coo=matmul(m,Tbackbone_backup(j)%coo-CA)+CA
    enddo
enddo

do i=1, Tgroup(ran_resi(3))%cnum1
    if(Tgroup(ran_resi(3))%atype1(i)=="N") then
        N3=Tgroup(ran_resi(3))%coo1(:,i)
    elseif(Tgroup(ran_resi(3))%atype1(i)=="CA") then
        CA3=Tgroup(ran_resi(3))%coo1(:,i)
    endif
enddo

if(Tgroup(ran_resi(3))%gtype.ne."CPRO") then
    call RANDOM_NUMBER(ran2)
    ic=anint(ran2*180)
    delta_phi3=-(ic-90)
    cos_angle=cosd(delta_phi3)
    sin_angle=sind(delta_phi3)

    rotaxis=CA3-N3
    call normalize_vector(rotaxis)
    call axisrotation(rotaxis, cos_angle, sin_angle, m)

    temp_num=0
    do j=1, Tgroup(ran_resi(3))%cnum1
        if(Tgroup(ran_resi(3))%atype1(j)=="HA".or.Tgroup(ran_resi(3))%atype1(j)=="HA2" & 
        .or.Tgroup(ran_resi(3))%atype1(j)=="HA3") then
            temp_num=temp_num+1
            atype(temp_num)=Tgroup(ran_resi(3))%atype1(j)
            coo(:,temp_num)=Tgroup(ran_resi(3))%coo1(:,j)-N3
        endif
    enddo

    coo=matmul(m,coo)

    do l=1, temp_num
        do j=1, Tgroup(ran_resi(3))%cnum1
            if(atype(l)==Tgroup(ran_resi(3))%atype1(j)) then
                Tgroup(ran_resi(3))%coo1(:,j)=coo(:,l)+N3
            endif
        enddo
    enddo

    do j=1, Tgroup(ran_resi(3))%cnum2
        Tgroup(ran_resi(3))%coo2(:,j)=matmul(m,Tgroup(ran_resi(3))%coo2(:,j)-N3)+N3
    enddo
    do j=1, Tgroup(ran_resi(3))%cnum3
        Tgroup(ran_resi(3))%coo3(:,j)=matmul(m,Tgroup(ran_resi(3))%coo3(:,j)-N3)+N3
    enddo
endif

do i=1, Tgroup(ran_resi(3))%cnum3
    if(Tgroup(ran_resi(3))%atype3(i)=="C") then
        C3=Tgroup(ran_resi(3))%coo3(:,i)
        exit
    endif
enddo

call RANDOM_NUMBER(ran2)
ic=anint(ran2*180)
delta_psi3=-(ic-90)
cos_angle=cosd(delta_psi3)
sin_angle=sind(delta_psi3)

rotaxis=C3-CA3
call normalize_vector(rotaxis)
call axisrotation(rotaxis, cos_angle, sin_angle, m)

do j=1, Tgroup(ran_resi(3))%cnum3
    Tgroup(ran_resi(3))%coo3(:,j)=matmul(m,Tgroup(ran_resi(3))%coo3(:,j)-CA3)+CA3
enddo

return
end subroutine backbone_rotation_Cterm_multiangle

subroutine backbone_rotation_center_multiangle(group_pep, backbone_backup, ran_resi, phi1, psi1, &
    phi2, psi2, Tgroup_pep, Tbackbone_backup)
implicit none
integer                 :: ran_resi(3), i, j, k, l, temp_num
real                    :: C0(3), N1(3), CA1(3), C1(3), N2(3), CA2(3), C2(3), N3(3), N(3), CA(3), C(3)
REAL                    :: rotaxis(3), m(3,3)
real                    :: phi_old(2), psi_old(2), phi_new(2), psi_new(2), phi1, psi1, phi2, psi2
real                    :: delta_phi, delta_psi, cos_angle, sin_angle
real                    :: coo(3,20)
character*4             :: atype(20)
type(groupdetails)      :: group_pep(pep_res), Tgroup_pep(pep_res)
type(databackup)        :: backbone_backup(pep_res), Tbackbone_backup(pep_res)


res_coords_flags(ran_resi(1):ran_resi(3)) = 1 !tell system to reload peptide coordinates into array format prior to calculating energies
phi_new(1)=phi1; phi_new(2)=phi2
psi_new(1)=psi1; psi_new(2)=psi2

do i=1, group_pep(ran_resi(1)-1)%cnum3
    if(group_pep(ran_resi(1)-1)%atype3(i)=="C") then
        C0=group_pep(ran_resi(1)-1)%coo3(:,i)
        exit
    endif
enddo

do i=1, group_pep(ran_resi(1))%cnum1
    if(group_pep(ran_resi(1))%atype1(i)=="N") then
        N1=group_pep(ran_resi(1))%coo1(:,i)
    elseif(group_pep(ran_resi(1))%atype1(i)=="CA") then
        CA1=group_pep(ran_resi(1))%coo1(:,i)
    endif
enddo

do i=1, group_pep(ran_resi(1))%cnum3
    if(group_pep(ran_resi(1))%atype3(i)=="C") then
        C1=group_pep(ran_resi(1))%coo3(:,i)
        exit
    endif
enddo

do i=1, group_pep(ran_resi(2))%cnum1
    if(group_pep(ran_resi(2))%atype1(i)=="N") then
        N2=group_pep(ran_resi(2))%coo1(:,i)
    elseif(group_pep(ran_resi(2))%atype1(i)=="CA") then
        CA2=group_pep(ran_resi(2))%coo1(:,i)
    endif
enddo

do i=1, group_pep(ran_resi(2))%cnum3
    if(group_pep(ran_resi(2))%atype3(i)=="C") then
        C2=group_pep(ran_resi(2))%coo3(:,i)
        exit
    endif
enddo

do i=1, group_pep(ran_resi(3))%cnum1
    if(group_pep(ran_resi(3))%atype1(i)=="N") then
        N3=group_pep(ran_resi(3))%coo1(:,i)
        exit
    endif
enddo

call calc_dihedral(C0, N1, CA1, C1, phi_old(1))
call calc_dihedral(N1, CA1, C1, N2, psi_old(1))
call calc_dihedral(C1, N2, CA2, C2, phi_old(2))
call calc_dihedral(N2, CA2, C2, N3, psi_old(2))

Tgroup_pep=group_pep
Tbackbone_backup=backbone_backup

do i=1, 2
    delta_phi=-(phi_new(i)-phi_old(i))

    cos_angle=cosd(delta_phi)
    sin_angle=sind(delta_phi)

    do k=1, Tgroup_pep(ran_resi(i))%cnum1
        if(Tgroup_pep(ran_resi(i))%atype1(k)=="N") then
            N=Tgroup_pep(ran_resi(i))%coo1(:,k)
        elseif(Tgroup_pep(ran_resi(i))%atype1(k)=="CA") then
            CA=Tgroup_pep(ran_resi(i))%coo1(:,k)
        endif
    enddo

    rotaxis=CA-N
    call normalize_vector(rotaxis)
    call axisrotation(rotaxis, cos_angle, sin_angle, m)

    temp_num=0
        do j=1, Tgroup_pep(ran_resi(i))%cnum1
        if(Tgroup_pep(ran_resi(i))%atype1(j)=="HA".or.Tgroup_pep(ran_resi(i))%atype1(j)=="HA2".or.Tgroup_pep(ran_resi(i))%atype1(j)=="HA3") then
            temp_num=temp_num+1
            atype(temp_num)=Tgroup_pep(ran_resi(i))%atype1(j)
            coo(:,temp_num)=Tgroup_pep(ran_resi(i))%coo1(:,j)-N
        endif
    enddo

    coo=matmul(m,coo)
    do l=1, temp_num
        do j=1, Tgroup_pep(ran_resi(i))%cnum1
            if(atype(l)==Tgroup_pep(ran_resi(i))%atype1(j)) then
                Tgroup_pep(ran_resi(i))%coo1(:,j)=coo(:,l)+N
            endif
        enddo
    enddo

    do j=1, Tgroup_pep(ran_resi(i))%cnum2
        Tgroup_pep(ran_resi(i))%coo2(:,j)=matmul(m,Tgroup_pep(ran_resi(i))%coo2(:,j)-N)+N
    enddo
    do j=1, Tgroup_pep(ran_resi(i))%cnum3
        Tgroup_pep(ran_resi(i))%coo3(:,j)=matmul(m,Tgroup_pep(ran_resi(i))%coo3(:,j)-N)+N
    enddo

    do j=ran_resi(i)+1, ran_resi(3)
        do l=1, Tgroup_pep(j)%cnum1
            Tgroup_pep(j)%coo1(:,j)=matmul(m,Tgroup_pep(j)%coo1(:,j)-N)+N
        enddo
        do l=1, Tgroup_pep(j)%cnum2
            Tgroup_pep(j)%coo2(:,j)=matmul(m,Tgroup_pep(j)%coo2(:,j)-N)+N
        enddo
        if(j.ne.ran_resi(3)) then
            do l=1, Tgroup_pep(j)%cnum3
                Tgroup_pep(j)%coo3(:,j)=matmul(m,Tgroup_pep(j)%coo3(:,j)-N)+N
            enddo
        endif
        Tbackbone_backup(j)%coo=matmul(m,Tbackbone_backup(j)%coo-N)+N
    enddo

    delta_psi=-(psi_new(i)-psi_old(i))

    cos_angle=cosd(delta_psi)
    sin_angle=sind(delta_psi)

    do k=1, Tgroup_pep(ran_resi(i))%cnum3
        if(Tgroup_pep(ran_resi(i))%atype3(k)=="C") then
            C=Tgroup_pep(ran_resi(i))%coo3(:,k)
            exit
        endif
    enddo

    rotaxis=C-CA
    call normalize_vector(rotaxis)
    call axisrotation(rotaxis, cos_angle, sin_angle, m)

    do j=1, Tgroup_pep(ran_resi(i))%cnum3
        Tgroup_pep(ran_resi(i))%coo3(:,j)=matmul(m,Tgroup_pep(ran_resi(i))%coo3(:,j)-CA)+CA
    enddo

    do j=ran_resi(i)+1, ran_resi(3)
        do l=1, Tgroup_pep(j)%cnum1
            Tgroup_pep(j)%coo1(:,j)=matmul(m,Tgroup_pep(j)%coo1(:,j)-CA)+CA
        enddo
        do l=1, Tgroup_pep(j)%cnum2
            Tgroup_pep(j)%coo2(:,j)=matmul(m,Tgroup_pep(j)%coo2(:,j)-CA)+CA
        enddo
        if(j.ne.ran_resi(3)) then
            do l=1, Tgroup_pep(j)%cnum3
                Tgroup_pep(j)%coo3(:,j)=matmul(m,Tgroup_pep(j)%coo3(:,j)-CA)+CA
            enddo
        endif
        Tbackbone_backup(j)%coo=matmul(m,Tbackbone_backup(j)%coo-CA)+CA
    enddo
enddo

return
end subroutine backbone_rotation_center_multiangle

subroutine backbonemove_criterion(rmsd, feedback)
implicit none
integer                     :: feedback
real                        :: rmsd

feedback=0
if(rmsd.le.rmsd_max.and.rmsd.gt.rmsd_min) feedback=1

return
end subroutine backbonemove_criterion

subroutine mc_choose_aminoacid(ic, group_pep, aminoacid_name, res_type_old, res_type_new, trp_delta)
implicit none
integer                     :: ic, ip, ip1, i, trp_delta
integer                     :: res_type_old, res_type_new, res_avail_count, cutoffs(6)
real                        :: ran2
character*4                 :: aminoacid_name
character*3                 :: res_avail(20)
type(groupdetails)          :: group_pep(gnum)

if(ph_value.le.3.9) then
    if(group_pep(ic)%gtype=="GLY") then
        ip=10
    elseif(group_pep(ic)%gtype=="LEU".or.group_pep(ic)%gtype=="VAL".or.group_pep(ic)%gtype=="ILE".or.group_pep(ic)%gtype=="MET".or. &
        group_pep(ic)%gtype=="PHE".or.group_pep(ic)%gtype=="TYR".or.group_pep(ic)%gtype=="TRP") then
        ip=20
    elseif(group_pep(ic)%gtype=="ARG".or.group_pep(ic)%gtype=="LYS".or.group_pep(ic)%gtype=="HIP") then
        ip=30
    elseif(group_pep(ic)%gtype=="SER".or.group_pep(ic)%gtype=="THR".or.group_pep(ic)%gtype=="ASN".or.group_pep(ic)%gtype=="GLN".or. &
            group_pep(ic)%gtype=="GLH".or.group_pep(ic)%gtype=="ASH") then
        ip=40
    elseif(group_pep(ic)%gtype=="PRO".or.group_pep(ic)%gtype=="CYS".or.group_pep(ic)%gtype=="ALA") then
        ip=50
    elseif(group_pep(ic)%gtype=="NGLY") then
        ip=11
    elseif(group_pep(ic)%gtype=="NLEU".or.group_pep(ic)%gtype=="NVAL".or.group_pep(ic)%gtype=="NILE".or.group_pep(ic)%gtype=="NMET".or. &
        group_pep(ic)%gtype=="NPHE".or.group_pep(ic)%gtype=="NTYR".or.group_pep(ic)%gtype=="NTRP") then
        ip=21
    elseif(group_pep(ic)%gtype=="NARG".or.group_pep(ic)%gtype=="NLYS".or.group_pep(ic)%gtype=="NHIP") then
        ip=31
    elseif(group_pep(ic)%gtype=="NSER".or.group_pep(ic)%gtype=="NTHR".or.group_pep(ic)%gtype=="NASN".or.group_pep(ic)%gtype=="NGLN".or. &
            group_pep(ic)%gtype=="NGLH".or.group_pep(ic)%gtype=="NASH") then
        ip=41
    elseif(group_pep(ic)%gtype=="NPRO".or.group_pep(ic)%gtype=="NCYS".or.group_pep(ic)%gtype=="NALA") then
        ip=51
    elseif(group_pep(ic)%gtype=="CGLY") then
        ip=12
    elseif(group_pep(ic)%gtype=="CLEU".or.group_pep(ic)%gtype=="CVAL".or.group_pep(ic)%gtype=="CILE".or.group_pep(ic)%gtype=="CMET".or. &
        group_pep(ic)%gtype=="CPHE".or.group_pep(ic)%gtype=="CTYR".or.group_pep(ic)%gtype=="CTRP") then
        ip=22
    elseif(group_pep(ic)%gtype=="CARG".or.group_pep(ic)%gtype=="CLYS".or.group_pep(ic)%gtype=="CHIP") then
        ip=32
    elseif(group_pep(ic)%gtype=="CSER".or.group_pep(ic)%gtype=="CTHR".or.group_pep(ic)%gtype=="CASN".or.group_pep(ic)%gtype=="CGLN".or. &
            group_pep(ic)%gtype=="CGLH".or.group_pep(ic)%gtype=="CASH") then
        ip=42
    elseif(group_pep(ic)%gtype=="CPRO".or.group_pep(ic)%gtype=="CCYS".or.group_pep(ic)%gtype=="CALA") then
        ip=52
    endif

10 continue
    if(ip.eq.10) then
        aminoacid_name="GLY"
        goto 100
    elseif(ip.eq.20) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*7-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="LEU"
        elseif(ip1.eq.2) then
            aminoacid_name="VAL"
        elseif(ip1.eq.3) then
            aminoacid_name="ILE"
        elseif(ip1.eq.4) then
            aminoacid_name="MET"
        elseif(ip1.eq.5) then
            aminoacid_name="PHE"
        elseif(ip1.eq.6) then
            aminoacid_name="TYR"
        elseif(ip1.eq.7) then
            aminoacid_name="TRP"
        endif
    elseif(ip.eq.30) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*3-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="ARG"
        elseif(ip1.eq.2) then
            aminoacid_name="LYS"
        elseif(ip1.eq.3) then
            aminoacid_name="HIP"
        endif
    elseif(ip.eq.40) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*6-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="SER"
        elseif(ip1.eq.2) then
            aminoacid_name="THR"
        elseif(ip1.eq.3) then
            aminoacid_name="ASN"
        elseif(ip1.eq.4) then
            aminoacid_name="GLN"
        elseif(ip1.eq.5) then
            aminoacid_name="GLH"
        elseif(ip1.eq.6) then
            aminoacid_name="ASH"
        endif
    elseif(ip.eq.50) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*3-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="PRO"
        elseif(ip1.eq.2) then
            aminoacid_name="ALA"
        elseif(ip1.eq.3) then
            aminoacid_name="CYS"
        endif
    elseif(ip.eq.11) then
        aminoacid_name="NGLY"
        goto 100
    elseif(ip.eq.21) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*7-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NLEU"
        elseif(ip1.eq.2) then
            aminoacid_name="NVAL"
        elseif(ip1.eq.3) then
            aminoacid_name="NILE"
        elseif(ip1.eq.4) then
            aminoacid_name="NMET"
        elseif(ip1.eq.5) then
            aminoacid_name="NPHE"
        elseif(ip1.eq.6) then
            aminoacid_name="NTYR"
        elseif(ip1.eq.7) then
            aminoacid_name="NTRP"
        endif
    elseif(ip.eq.31) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*3-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NARG"
        elseif(ip1.eq.2) then
            aminoacid_name="NLYS"
        elseif(ip1.eq.3) then
            aminoacid_name="NHIP"
        endif
    elseif(ip.eq.41) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*6-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NSER"
        elseif(ip1.eq.2) then
            aminoacid_name="NTHR"
        elseif(ip1.eq.3) then
            aminoacid_name="NASN"
        elseif(ip1.eq.4) then
            aminoacid_name="NGLN"
        elseif(ip1.eq.5) then
            aminoacid_name="NGLH"
        elseif(ip1.eq.6) then
            aminoacid_name="NASH"
        endif
    elseif(ip.eq.51) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*3-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NPRO"
        elseif(ip1.eq.2) then
            aminoacid_name="NALA"
        elseif(ip1.eq.3) then
            aminoacid_name="NCYS"
        endif
    elseif(ip.eq.12) then
        aminoacid_name="CGLY"
        goto 100
    elseif(ip.eq.22) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*7-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CLEU"
        elseif(ip1.eq.2) then
            aminoacid_name="CVAL"
        elseif(ip1.eq.3) then
            aminoacid_name="CILE"
        elseif(ip1.eq.4) then
            aminoacid_name="CMET"
        elseif(ip1.eq.5) then
            aminoacid_name="CPHE"
        elseif(ip1.eq.6) then
            aminoacid_name="CTYR"
        elseif(ip1.eq.7) then
            aminoacid_name="CTRP"
        endif
    elseif(ip.eq.32) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*3-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CARG"
        elseif(ip1.eq.2) then
            aminoacid_name="CLYS"
        elseif(ip1.eq.3) then
            aminoacid_name="CHIP"
        endif
    elseif(ip.eq.42) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*6-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CSER"
        elseif(ip1.eq.2) then
            aminoacid_name="CTHR"
        elseif(ip1.eq.3) then
            aminoacid_name="CASN"
        elseif(ip1.eq.4) then
            aminoacid_name="CGLN"
        elseif(ip1.eq.5) then
            aminoacid_name="CGLH"
        elseif(ip1.eq.6) then
            aminoacid_name="CASH"
        endif
    elseif(ip.eq.52) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*3-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CPRO"
        elseif(ip1.eq.2) then
            aminoacid_name="CALA"
        elseif(ip1.eq.3) then
            aminoacid_name="CCYS"
        endif
    endif
    if(aminoacid_name==group_pep(ic)%gtype) goto 10

elseif(ph_value.gt.3.9.and.ph_value.le.4.3) then
    if(group_pep(ic)%gtype=="GLY") then
        ip=10
    elseif(group_pep(ic)%gtype=="LEU".or.group_pep(ic)%gtype=="VAL".or.group_pep(ic)%gtype=="ILE".or.group_pep(ic)%gtype=="MET".or. &
        group_pep(ic)%gtype=="PHE".or.group_pep(ic)%gtype=="TYR".or.group_pep(ic)%gtype=="TRP") then
        ip=20
    elseif(group_pep(ic)%gtype=="ARG".or.group_pep(ic)%gtype=="LYS".or.group_pep(ic)%gtype=="HIP") then
        ip=30
    elseif(group_pep(ic)%gtype=="SER".or.group_pep(ic)%gtype=="THR".or.group_pep(ic)%gtype=="ASN".or.group_pep(ic)%gtype=="GLN".or. &
            group_pep(ic)%gtype=="GLH") then
        ip=40
    elseif(group_pep(ic)%gtype=="PRO".or.group_pep(ic)%gtype=="CYS".or.group_pep(ic)%gtype=="ALA") then
        ip=50
    elseif(group_pep(ic)%gtype=="ASP") then
        ip=60
    elseif(group_pep(ic)%gtype=="NGLY") then
        ip=11
    elseif(group_pep(ic)%gtype=="NLEU".or.group_pep(ic)%gtype=="NVAL".or.group_pep(ic)%gtype=="NILE".or.group_pep(ic)%gtype=="NMET".or. &
        group_pep(ic)%gtype=="NPHE".or.group_pep(ic)%gtype=="NTYR".or.group_pep(ic)%gtype=="NTRP") then
        ip=21
    elseif(group_pep(ic)%gtype=="NARG".or.group_pep(ic)%gtype=="NLYS".or.group_pep(ic)%gtype=="NHIP") then
        ip=31
    elseif(group_pep(ic)%gtype=="NSER".or.group_pep(ic)%gtype=="NTHR".or.group_pep(ic)%gtype=="NASN".or.group_pep(ic)%gtype=="NGLN".or. &
            group_pep(ic)%gtype=="NGLH") then
        ip=41
    elseif(group_pep(ic)%gtype=="NPRO".or.group_pep(ic)%gtype=="NCYS".or.group_pep(ic)%gtype=="NALA") then
        ip=51
    elseif(group_pep(ic)%gtype=="NASP") then
        ip=61
    elseif(group_pep(ic)%gtype=="CGLY") then
        ip=12
    elseif(group_pep(ic)%gtype=="CLEU".or.group_pep(ic)%gtype=="CVAL".or.group_pep(ic)%gtype=="CILE".or.group_pep(ic)%gtype=="CMET".or. &
        group_pep(ic)%gtype=="CPHE".or.group_pep(ic)%gtype=="CTYR".or.group_pep(ic)%gtype=="CTRP") then
        ip=22
    elseif(group_pep(ic)%gtype=="CARG".or.group_pep(ic)%gtype=="CLYS".or.group_pep(ic)%gtype=="CHIP") then
        ip=32
    elseif(group_pep(ic)%gtype=="CSER".or.group_pep(ic)%gtype=="CTHR".or.group_pep(ic)%gtype=="CASN".or.group_pep(ic)%gtype=="CGLN".or. &
            group_pep(ic)%gtype=="CGLH") then
        ip=42
    elseif(group_pep(ic)%gtype=="CPRO".or.group_pep(ic)%gtype=="CCYS".or.group_pep(ic)%gtype=="CALA") then
        ip=52
    elseif(group_pep(ic)%gtype=="CASP") then
        ip=62
    endif

20 continue
    if(ip.eq.10) then
        aminoacid_name="GLY"
        goto 100
    elseif(ip.eq.20) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*7-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="LEU"
        elseif(ip1.eq.2) then
            aminoacid_name="VAL"
        elseif(ip1.eq.3) then
            aminoacid_name="ILE"
        elseif(ip1.eq.4) then
            aminoacid_name="MET"
        elseif(ip1.eq.5) then
            aminoacid_name="PHE"
        elseif(ip1.eq.6) then
            aminoacid_name="TYR"
        elseif(ip1.eq.7) then
            aminoacid_name="TRP"
        endif
    elseif(ip.eq.30) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*3-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="ARG"
        elseif(ip1.eq.2) then
            aminoacid_name="LYS"
        elseif(ip1.eq.3) then
            aminoacid_name="HIP"
        endif
    elseif(ip.eq.40) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*5-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="SER"
        elseif(ip1.eq.2) then
            aminoacid_name="THR"
        elseif(ip1.eq.3) then
            aminoacid_name="ASN"
        elseif(ip1.eq.4) then
            aminoacid_name="GLN"
        elseif(ip1.eq.5) then
            aminoacid_name="GLH"
        endif
    elseif(ip.eq.50) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*3-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="PRO"
        elseif(ip1.eq.2) then
            aminoacid_name="ALA"
        elseif(ip1.eq.3) then
            aminoacid_name="CYS"
        endif
    elseif(ip.eq.60) then
        aminoacid_name="ASP"
        goto 100
    elseif(ip.eq.11) then
        aminoacid_name="NGLY"
        goto 100
    elseif(ip.eq.21) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*7-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NLEU"
        elseif(ip1.eq.2) then
            aminoacid_name="NVAL"
        elseif(ip1.eq.3) then
            aminoacid_name="NILE"
        elseif(ip1.eq.4) then
            aminoacid_name="NMET"
        elseif(ip1.eq.5) then
            aminoacid_name="NPHE"
        elseif(ip1.eq.6) then
            aminoacid_name="NTYR"
        elseif(ip1.eq.7) then
            aminoacid_name="NTRP"
        endif
    elseif(ip.eq.31) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*3-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NARG"
        elseif(ip1.eq.2) then
            aminoacid_name="NLYS"
        elseif(ip1.eq.3) then
            aminoacid_name="NHIP"
        endif
    elseif(ip.eq.41) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*5-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NSER"
        elseif(ip1.eq.2) then
            aminoacid_name="NTHR"
        elseif(ip1.eq.3) then
            aminoacid_name="NASN"
        elseif(ip1.eq.4) then
            aminoacid_name="NGLN"
        elseif(ip1.eq.5) then
            aminoacid_name="NGLH"
        endif
    elseif(ip.eq.51) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*3-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NPRO"
        elseif(ip1.eq.2) then
            aminoacid_name="NALA"
        elseif(ip1.eq.3) then
            aminoacid_name="NCYS"
        endif
    elseif(ip.eq.61) then
        aminoacid_name="NASP"
        goto 100
    elseif(ip.eq.12) then
        aminoacid_name="CGLY"
        goto 100
    elseif(ip.eq.22) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*7-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CLEU"
        elseif(ip1.eq.2) then
            aminoacid_name="CVAL"
        elseif(ip1.eq.3) then
            aminoacid_name="CILE"
        elseif(ip1.eq.4) then
            aminoacid_name="CMET"
        elseif(ip1.eq.5) then
            aminoacid_name="CPHE"
        elseif(ip1.eq.6) then
            aminoacid_name="CTYR"
        elseif(ip1.eq.7) then
            aminoacid_name="CTRP"
        endif
    elseif(ip.eq.32) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*3-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CARG"
        elseif(ip1.eq.2) then
            aminoacid_name="CLYS"
        elseif(ip1.eq.3) then
            aminoacid_name="CHIP"
        endif
    elseif(ip.eq.42) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*5-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CSER"
        elseif(ip1.eq.2) then
            aminoacid_name="CTHR"
        elseif(ip1.eq.3) then
            aminoacid_name="CASN"
        elseif(ip1.eq.4) then
            aminoacid_name="CGLN"
        elseif(ip1.eq.5) then
            aminoacid_name="CGLH"
        endif
    elseif(ip.eq.52) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*3-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CPRO"
        elseif(ip1.eq.2) then
            aminoacid_name="CALA"
        elseif(ip1.eq.3) then
            aminoacid_name="CCYS"
        endif
    elseif(ip.eq.62) then
        aminoacid_name="CASP"
        goto 100
    endif
    if(aminoacid_name==group_pep(ic)%gtype) goto 20

elseif(ph_value.gt.4.3.and.ph_value.le.6.0) then
    if(group_pep(ic)%gtype=="GLY") then
        ip=10
    elseif(group_pep(ic)%gtype=="LEU".or.group_pep(ic)%gtype=="VAL".or.group_pep(ic)%gtype=="ILE".or.group_pep(ic)%gtype=="MET".or. &
        group_pep(ic)%gtype=="PHE".or.group_pep(ic)%gtype=="TYR".or.group_pep(ic)%gtype=="TRP") then
        ip=20
    elseif(group_pep(ic)%gtype=="ARG".or.group_pep(ic)%gtype=="LYS".or.group_pep(ic)%gtype=="HIP") then
        ip=30
    elseif(group_pep(ic)%gtype=="SER".or.group_pep(ic)%gtype=="THR".or.group_pep(ic)%gtype=="ASN".or.group_pep(ic)%gtype=="GLN") then
        ip=40
    elseif(group_pep(ic)%gtype=="PRO".or.group_pep(ic)%gtype=="CYS".or.group_pep(ic)%gtype=="ALA") then
        ip=50
    elseif(group_pep(ic)%gtype=="GLU".or.group_pep(ic)%gtype=="ASP") then
        ip=60
    elseif(group_pep(ic)%gtype=="NGLY") then
        ip=11
    elseif(group_pep(ic)%gtype=="NLEU".or.group_pep(ic)%gtype=="NVAL".or.group_pep(ic)%gtype=="NILE".or.group_pep(ic)%gtype=="NMET".or. &
        group_pep(ic)%gtype=="NPHE".or.group_pep(ic)%gtype=="NTYR".or.group_pep(ic)%gtype=="NTRP") then
        ip=21
    elseif(group_pep(ic)%gtype=="NARG".or.group_pep(ic)%gtype=="NLYS".or.group_pep(ic)%gtype=="NHIP") then
        ip=31
    elseif(group_pep(ic)%gtype=="NSER".or.group_pep(ic)%gtype=="NTHR".or.group_pep(ic)%gtype=="NASN".or.group_pep(ic)%gtype=="NGLN") then
        ip=41
    elseif(group_pep(ic)%gtype=="NPRO".or.group_pep(ic)%gtype=="NCYS".or.group_pep(ic)%gtype=="NALA") then
        ip=51
    elseif(group_pep(ic)%gtype=="NGLU".or.group_pep(ic)%gtype=="NASP") then
        ip=61
    elseif(group_pep(ic)%gtype=="CGLY") then
        ip=12
    elseif(group_pep(ic)%gtype=="CLEU".or.group_pep(ic)%gtype=="CVAL".or.group_pep(ic)%gtype=="CILE".or.group_pep(ic)%gtype=="CMET".or. &
        group_pep(ic)%gtype=="CPHE".or.group_pep(ic)%gtype=="CTYR".or.group_pep(ic)%gtype=="CTRP") then
        ip=22
    elseif(group_pep(ic)%gtype=="CARG".or.group_pep(ic)%gtype=="CLYS".or.group_pep(ic)%gtype=="CHIP") then
        ip=32
    elseif(group_pep(ic)%gtype=="CSER".or.group_pep(ic)%gtype=="CTHR".or.group_pep(ic)%gtype=="CASN".or.group_pep(ic)%gtype=="CGLN") then
        ip=42
    elseif(group_pep(ic)%gtype=="CPRO".or.group_pep(ic)%gtype=="CCYS".or.group_pep(ic)%gtype=="CALA") then
        ip=52
    elseif(group_pep(ic)%gtype=="CGLU".or.group_pep(ic)%gtype=="CASP") then
        ip=62
    endif

30 continue
    if(ip.eq.10) then
        aminoacid_name="GLY"
        goto 100
    elseif(ip.eq.20) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*7-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="LEU"
        elseif(ip1.eq.2) then
            aminoacid_name="VAL"
        elseif(ip1.eq.3) then
            aminoacid_name="ILE"
        elseif(ip1.eq.4) then
            aminoacid_name="MET"
        elseif(ip1.eq.5) then
            aminoacid_name="PHE"
        elseif(ip1.eq.6) then
            aminoacid_name="TYR"
        elseif(ip1.eq.7) then
            aminoacid_name="TRP"
        endif
    elseif(ip.eq.30) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*3-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="ARG"
        elseif(ip1.eq.2) then
            aminoacid_name="LYS"
        elseif(ip1.eq.3) then
            aminoacid_name="HIP"
        endif
    elseif(ip.eq.40) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*4-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="SER"
        elseif(ip1.eq.2) then
            aminoacid_name="THR"
        elseif(ip1.eq.3) then
            aminoacid_name="ASN"
        elseif(ip1.eq.4) then
            aminoacid_name="GLN"
        endif
    elseif(ip.eq.50) then
        call RANDOM_NUMBER(ran2)
!			ip1=int(ran2*3-1.0e-3)+1
        ip1=int(ran2*2-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="PRO"
        elseif(ip1.eq.2) then
            aminoacid_name="ALA"
!			elseif(ip1.eq.3) then
!				aminoacid_name="CYS"
        endif
    elseif(ip.eq.60) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*2-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="GLU"
        elseif(ip1.eq.2) then
            aminoacid_name="ASP"
        endif
    elseif(ip.eq.11) then
        aminoacid_name="NGLY"
        goto 100
    elseif(ip.eq.21) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*7-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NLEU"
        elseif(ip1.eq.2) then
            aminoacid_name="NVAL"
        elseif(ip1.eq.3) then
            aminoacid_name="NILE"
        elseif(ip1.eq.4) then
            aminoacid_name="NMET"
        elseif(ip1.eq.5) then
            aminoacid_name="NPHE"
        elseif(ip1.eq.6) then
            aminoacid_name="NTYR"
        elseif(ip1.eq.7) then
            aminoacid_name="NTRP"
        endif
    elseif(ip.eq.31) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*3-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NARG"
        elseif(ip1.eq.2) then
            aminoacid_name="NLYS"
        elseif(ip1.eq.3) then
            aminoacid_name="NHIP"
        endif
    elseif(ip.eq.41) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*4-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NSER"
        elseif(ip1.eq.2) then
            aminoacid_name="NTHR"
        elseif(ip1.eq.3) then
            aminoacid_name="NASN"
        elseif(ip1.eq.4) then
            aminoacid_name="NGLN"
        endif
    elseif(ip.eq.51) then
        call RANDOM_NUMBER(ran2)
!			ip1=int(ran2*3-1.0e-3)+1
        ip1=int(ran2*2-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NPRO"
        elseif(ip1.eq.2) then
            aminoacid_name="NALA"
!			elseif(ip1.eq.3) then
!				aminoacid_name="NCYS"
        endif
    elseif(ip.eq.61) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*2-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NGLU"
        elseif(ip1.eq.2) then
            aminoacid_name="NASP"
        endif
    elseif(ip.eq.12) then
        aminoacid_name="CGLY"
        goto 100
    elseif(ip.eq.22) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*7-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CLEU"
        elseif(ip1.eq.2) then
            aminoacid_name="CVAL"
        elseif(ip1.eq.3) then
            aminoacid_name="CILE"
        elseif(ip1.eq.4) then
            aminoacid_name="CMET"
        elseif(ip1.eq.5) then
            aminoacid_name="CPHE"
        elseif(ip1.eq.6) then
            aminoacid_name="CTYR"
        elseif(ip1.eq.7) then
            aminoacid_name="CTRP"
        endif
    elseif(ip.eq.32) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*3-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CARG"
        elseif(ip1.eq.2) then
            aminoacid_name="CLYS"
        elseif(ip1.eq.3) then
            aminoacid_name="CHIP"
        endif
    elseif(ip.eq.42) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*4-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CSER"
        elseif(ip1.eq.2) then
            aminoacid_name="CTHR"
        elseif(ip1.eq.3) then
            aminoacid_name="CASN"
        elseif(ip1.eq.4) then
            aminoacid_name="CGLN"
        endif
    elseif(ip.eq.52) then
        call RANDOM_NUMBER(ran2)
!			ip1=int(ran2*3-1.0e-3)+1
        ip1=int(ran2*2-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CPRO"
        elseif(ip1.eq.2) then
            aminoacid_name="CALA"
!			elseif(ip1.eq.3) then
!				aminoacid_name="CCYS"
        endif
    elseif(ip.eq.62) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*2-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CGLU"
        elseif(ip1.eq.2) then
            aminoacid_name="CASP"
        endif
    endif
    if(aminoacid_name==group_pep(ic)%gtype) goto 30

elseif(ph_value.gt.6.0.and.ph_value.lt.8.3) then

    !determine what type of amino acid is being replaced
    if (group_pep(ic)%gtype=="GLY".or.group_pep(ic)%gtype=="NGLY".or.group_pep(ic)%gtype=="CGLY") then
        res_type_old = 1

    elseif(group_pep(ic)%gtype=="LEU".or.group_pep(ic)%gtype=="NLEU".or.group_pep(ic)%gtype=="CLEU".or. &
    group_pep(ic)%gtype=="VAL".or.group_pep(ic)%gtype=="NVAL".or.group_pep(ic)%gtype=="CVAL".or. &
    group_pep(ic)%gtype=="ILE".or.group_pep(ic)%gtype=="NILE".or.group_pep(ic)%gtype=="CILE".or. &
    group_pep(ic)%gtype=="MET".or.group_pep(ic)%gtype=="NMET".or.group_pep(ic)%gtype=="CMET".or. &
    group_pep(ic)%gtype=="PHE".or.group_pep(ic)%gtype=="NPHE".or.group_pep(ic)%gtype=="CPHE".or. &
    group_pep(ic)%gtype=="TYR".or.group_pep(ic)%gtype=="NTYR".or.group_pep(ic)%gtype=="CTYR".or. &
    group_pep(ic)%gtype=="TRP".or.group_pep(ic)%gtype=="NTRP".or.group_pep(ic)%gtype=="CTRP") then
        res_type_old = 2

    elseif (group_pep(ic)%gtype=="ARG".or.group_pep(ic)%gtype=="NARG".or.group_pep(ic)%gtype=="CARG".or. &
    group_pep(ic)%gtype=="LYS".or.group_pep(ic)%gtype=="NLYS".or.group_pep(ic)%gtype=="CLYS") then
        res_type_old = 3

    elseif(group_pep(ic)%gtype=="SER".or.group_pep(ic)%gtype=="NSER".or.group_pep(ic)%gtype=="CSER".or. &
    group_pep(ic)%gtype=="THR".or.group_pep(ic)%gtype=="NTHR".or.group_pep(ic)%gtype=="CTHR".or. &
    group_pep(ic)%gtype=="ASN".or.group_pep(ic)%gtype=="NASN".or.group_pep(ic)%gtype=="CASN".or. &
    group_pep(ic)%gtype=="GLN".or.group_pep(ic)%gtype=="NGLN".or.group_pep(ic)%gtype=="CGLN".or. &
    group_pep(ic)%gtype=="HIE".or.group_pep(ic)%gtype=="NHIE".or.group_pep(ic)%gtype=="CHIE") then
        res_type_old = 4

    elseif(group_pep(ic)%gtype=="PRO".or.group_pep(ic)%gtype=="NPRO".or.group_pep(ic)%gtype=="CPRO".or. &
    group_pep(ic)%gtype=="ALA".or.group_pep(ic)%gtype=="NALA".or.group_pep(ic)%gtype=="CALA".or. &
    group_pep(ic)%gtype=="CYS".or.group_pep(ic)%gtype=="NCYS".or.group_pep(ic)%gtype=="CCYS") then
        res_type_old = 5

    elseif(group_pep(ic)%gtype=="GLU".or.group_pep(ic)%gtype=="NGLU".or.group_pep(ic)%gtype=="CGLU".or. &
    group_pep(ic)%gtype=="ASP".or.group_pep(ic)%gtype=="NASP".or.group_pep(ic)%gtype=="CASP") then
        res_type_old = 6

    else
        !open(20, file="error.txt", access="append")
            write(20,"(A)") "Residue name not recognized during residue replacement, terminating program"
            close(20)
        !stop
    endif
    
!determine which amino acid types can replace the selected amino acid
    res_type_new = 0
    if (below_hyd(res_type_old).eq.0) then
        !Have to replace the amino acid with one of the same type to stay within hydration range
        if (res_type_old.eq.1) then
            res_avail(1) = "GLY"
            res_avail_count = 1
            res_type_new = 1

        elseif (res_type_old.eq.2) then
            !need an extra if statement to not go over the tryptophan limit
            if (trp_count.ge.trp_limit.and.(group_pep(ic)%gtype.ne."NTRP".or.group_pep(ic)%gtype.ne."TRP".or.group_pep(ic)%gtype.ne."CTRP")) then
                res_avail(1:6) = (/"LEU", "ILE", "MET", "VAL", "TYR", "PHE"/)
                res_avail_count = 6
                res_type_new = 2
            else
                res_avail(1:7) = (/"LEU", "ILE", "MET", "VAL", "TYR", "PHE", "TRP"/)
                res_avail_count = 7
                res_type_new = 2
            endif

        elseif (res_type_old.eq.3) then
            res_avail(1:2) = (/"ARG", "LYS"/)
            res_avail_count = 2
            res_type_new = 3

        elseif (res_type_old.eq.4) then
            res_avail(1:5) = (/"SER", "THR", "ASN", "GLN", "HIE"/)
            res_avail_count = 5
            res_type_new = 4

        elseif (res_type_old.eq.5) then
            res_avail(1) = "ALA"
            res_avail_count = 1
            res_type_new = 5

        else 
            res_avail(1:2) = (/"ASP", "GLU"/)
            res_avail_count = 2
            res_type_new = 6
        endif
    
    else
    !see which amino acid types we can increase while staying in the hydration range
        res_avail_count = 0  
        cutoffs = 0
        if (above_hyd(1).lt.0) then
            res_avail(res_avail_count+1) = ("GLY")
            res_avail_count = res_avail_count + 1
            cutoffs(1) = res_avail_count
        endif

        if (above_hyd(2).lt.0) then
            if (trp_count.ge.trp_limit.and.(group_pep(ic)%gtype.ne."NTRP".or.group_pep(ic)%gtype.ne."TRP".or.group_pep(ic)%gtype.ne."CTRP")) then
                res_avail(res_avail_count+1:res_avail_count+6) = (/"LEU", "ILE", "MET", "VAL", "TYR", "PHE"/)
                res_avail_count = res_avail_count + 6
            else
                res_avail(res_avail_count+1:res_avail_count+7) = (/"LEU", "ILE", "MET", "VAL", "TYR", "PHE", "TRP"/)
                res_avail_count = res_avail_count + 7
            endif
            cutoffs(2) = res_avail_count
        endif

        if (above_hyd(3).lt.0) then
            res_avail(res_avail_count+1:res_avail_count+2) = (/"ARG", "LYS"/)
            res_avail_count = res_avail_count + 2
            cutoffs(3) = res_avail_count
        endif

        if (above_hyd(4).lt.0) then
            res_avail(res_avail_count+1:res_avail_count+5) = (/"SER", "THR", "ASN", "GLN", "HIE"/)
            res_avail_count = res_avail_count + 5
            cutoffs(4) = res_avail_count
        endif

        if (above_hyd(5).lt.0) then
            res_avail(res_avail_count+1) = "ALA"
            res_avail_count = res_avail_count + 1
            cutoffs(5) = res_avail_count
        endif

        if (above_hyd(6).lt.0) then
            res_avail(res_avail_count+1:res_avail_count+2) = (/"ASP", "GLU"/)
            res_avail_count = res_avail_count + 2
            cutoffs(6) = res_avail_count
        endif

    endif

    !now pick one of the possible amino acids randomly
    call RANDOM_NUMBER(ran2)
    ip1=min(int(ran2*res_avail_count-1.0e-3)+1, res_avail_count)
    aminoacid_name = res_avail(ip1)
    
    !Check if tryptophan count is changing
    trp_delta = 0 
    if (group_pep(ic)%gtype.eq."TRP".or.group_pep(ic)%gtype.eq."NTRP".or.group_pep(ic)%gtype.eq."CTRP") then
        if (aminoacid_name.ne."TRP") trp_delta = -1
    else 
        if (aminoacid_name.eq."TRP") trp_delta = 1
    endif

    !if we are at one of the peptide termini, then add an N or C to amino acid name
    if (ic.eq.1) then
        aminoacid_name = "N" // aminoacid_name
    elseif (ic.eq.pep_res) then
        aminoacid_name = "C" // aminoacid_name
    endif

    !figure out the type of the replacement residue, if not already known
    if (res_type_new.eq.0) then
        i = 1
        do while (res_type_new.eq.0)
            if (cutoffs(i).ge.ip1) then
                res_type_new = i
            endif
            i = i + 1 
        enddo
    endif 

elseif(ph_value.ge.8.3.and.ph_value.lt.10.1) then
    if(group_pep(ic)%gtype=="GLY") then
        ip=10
    elseif(group_pep(ic)%gtype=="LEU".or.group_pep(ic)%gtype=="VAL".or.group_pep(ic)%gtype=="ILE".or.group_pep(ic)%gtype=="MET".or. &
        group_pep(ic)%gtype=="PHE".or.group_pep(ic)%gtype=="TYR".or.group_pep(ic)%gtype=="TRP") then
        ip=20
    elseif(group_pep(ic)%gtype=="ARG".or.group_pep(ic)%gtype=="LYS") then
        ip=30
    elseif(group_pep(ic)%gtype=="SER".or.group_pep(ic)%gtype=="THR".or.group_pep(ic)%gtype=="ASN".or.group_pep(ic)%gtype=="GLN".or. &
        group_pep(ic)%gtype=="HIE") then
        ip=40
    elseif(group_pep(ic)%gtype=="PRO".or.group_pep(ic)%gtype=="ALA") then
        ip=50
    elseif(group_pep(ic)%gtype=="GLU".or.group_pep(ic)%gtype=="ASP".or.group_pep(ic)%gtype=="CYT") then
        ip=60
    elseif(group_pep(ic)%gtype=="NGLY") then
        ip=11
    elseif(group_pep(ic)%gtype=="NLEU".or.group_pep(ic)%gtype=="NVAL".or.group_pep(ic)%gtype=="NILE".or.group_pep(ic)%gtype=="NMET".or. &
        group_pep(ic)%gtype=="NPHE".or.group_pep(ic)%gtype=="NTYR".or.group_pep(ic)%gtype=="NTRP") then
        ip=21
    elseif(group_pep(ic)%gtype=="NARG".or.group_pep(ic)%gtype=="NLYS") then
        ip=31
    elseif(group_pep(ic)%gtype=="NSER".or.group_pep(ic)%gtype=="NTHR".or.group_pep(ic)%gtype=="NASN".or.group_pep(ic)%gtype=="NGLN".or. &
        group_pep(ic)%gtype=="NHIE") then
        ip=41
    elseif(group_pep(ic)%gtype=="NPRO".or.group_pep(ic)%gtype=="NALA") then
        ip=51
    elseif(group_pep(ic)%gtype=="NGLU".or.group_pep(ic)%gtype=="NASP".or.group_pep(ic)%gtype=="NCYT") then
        ip=61
    elseif(group_pep(ic)%gtype=="CGLY") then
        ip=12
    elseif(group_pep(ic)%gtype=="CLEU".or.group_pep(ic)%gtype=="CVAL".or.group_pep(ic)%gtype=="CILE".or.group_pep(ic)%gtype=="CMET".or. &
        group_pep(ic)%gtype=="CPHE".or.group_pep(ic)%gtype=="CTYR".or.group_pep(ic)%gtype=="CTRP") then
        ip=22
    elseif(group_pep(ic)%gtype=="CARG".or.group_pep(ic)%gtype=="CLYS") then
        ip=32
    elseif(group_pep(ic)%gtype=="CSER".or.group_pep(ic)%gtype=="CTHR".or.group_pep(ic)%gtype=="CASN".or.group_pep(ic)%gtype=="CGLN".or. &
        group_pep(ic)%gtype=="CHIE") then
        ip=42
    elseif(group_pep(ic)%gtype=="CPRO".or.group_pep(ic)%gtype=="CALA") then
        ip=52
    elseif(group_pep(ic)%gtype=="CGLU".or.group_pep(ic)%gtype=="CASP".or.group_pep(ic)%gtype=="CCYT") then
        ip=62
    endif

50 continue
    if(ip.eq.10) then
        aminoacid_name="GLY"
        goto 100
    elseif(ip.eq.20) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*7-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="LEU"
        elseif(ip1.eq.2) then
            aminoacid_name="VAL"
        elseif(ip1.eq.3) then
            aminoacid_name="ILE"
        elseif(ip1.eq.4) then
            aminoacid_name="MET"
        elseif(ip1.eq.5) then
            aminoacid_name="PHE"
        elseif(ip1.eq.6) then
            aminoacid_name="TYR"
        elseif(ip1.eq.7) then
            aminoacid_name="TRP"
        endif
    elseif(ip.eq.30) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*2-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="ARG"
        elseif(ip1.eq.2) then
            aminoacid_name="LYS"
        endif
    elseif(ip.eq.40) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*5-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="SER"
        elseif(ip1.eq.2) then
            aminoacid_name="THR"
        elseif(ip1.eq.3) then
            aminoacid_name="ASN"
        elseif(ip1.eq.4) then
            aminoacid_name="GLN"
        elseif(ip1.eq.5) then
            aminoacid_name="HIE"
        endif
    elseif(ip.eq.50) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*2-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="PRO"
        elseif(ip1.eq.2) then
            aminoacid_name="ALA"
        endif
    elseif(ip.eq.60) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*3-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="GLU"
        elseif(ip1.eq.2) then
            aminoacid_name="ASP"
        elseif(ip1.eq.3) then
            aminoacid_name="CYT"
        endif
    elseif(ip.eq.11) then
        aminoacid_name="NGLY"
        goto 100
    elseif(ip.eq.21) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*7-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NLEU"
        elseif(ip1.eq.2) then
            aminoacid_name="NVAL"
        elseif(ip1.eq.3) then
            aminoacid_name="NILE"
        elseif(ip1.eq.4) then
            aminoacid_name="NMET"
        elseif(ip1.eq.5) then
            aminoacid_name="NPHE"
        elseif(ip1.eq.6) then
            aminoacid_name="NTYR"
        elseif(ip1.eq.7) then
            aminoacid_name="NTRP"
        endif
    elseif(ip.eq.31) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*2-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NARG"
        elseif(ip1.eq.2) then
            aminoacid_name="NLYS"
        endif
    elseif(ip.eq.41) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*5-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NSER"
        elseif(ip1.eq.2) then
            aminoacid_name="NTHR"
        elseif(ip1.eq.3) then
            aminoacid_name="NASN"
        elseif(ip1.eq.4) then
            aminoacid_name="NGLN"
        elseif(ip1.eq.5) then
            aminoacid_name="NHIE"
        endif
    elseif(ip.eq.51) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*2-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NPRO"
        elseif(ip1.eq.2) then
            aminoacid_name="NALA"
        endif
    elseif(ip.eq.61) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*3-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NGLU"
        elseif(ip1.eq.2) then
            aminoacid_name="NASP"
        elseif(ip1.eq.3) then
            aminoacid_name="NCYT"
        endif
    elseif(ip.eq.12) then
        aminoacid_name="CGLY"
        goto 100
    elseif(ip.eq.22) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*7-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CLEU"
        elseif(ip1.eq.2) then
            aminoacid_name="CVAL"
        elseif(ip1.eq.3) then
            aminoacid_name="CILE"
        elseif(ip1.eq.4) then
            aminoacid_name="CMET"
        elseif(ip1.eq.5) then
            aminoacid_name="CPHE"
        elseif(ip1.eq.6) then
            aminoacid_name="CTYR"
        elseif(ip1.eq.7) then
            aminoacid_name="CTRP"
        endif
    elseif(ip.eq.32) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*2-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CARG"
        elseif(ip1.eq.2) then
            aminoacid_name="CLYS"
        endif
    elseif(ip.eq.42) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*5-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CSER"
        elseif(ip1.eq.2) then
            aminoacid_name="CTHR"
        elseif(ip1.eq.3) then
            aminoacid_name="CASN"
        elseif(ip1.eq.4) then
            aminoacid_name="CGLN"
        elseif(ip1.eq.5) then
            aminoacid_name="CHIE"
        endif
    elseif(ip.eq.52) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*2-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CPRO"
        elseif(ip1.eq.2) then
            aminoacid_name="CALA"
        endif
    elseif(ip.eq.62) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*3-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CGLU"
        elseif(ip1.eq.2) then
            aminoacid_name="CASP"
        elseif(ip1.eq.3) then
            aminoacid_name="CCYT"
        endif
    endif
    if(aminoacid_name==group_pep(ic)%gtype) goto 50

elseif(ph_value.ge.10.1.and.ph_value.lt.10.5) then
    if(group_pep(ic)%gtype=="GLY") then
        ip=10
    elseif(group_pep(ic)%gtype=="LEU".or.group_pep(ic)%gtype=="VAL".or.group_pep(ic)%gtype=="ILE".or.group_pep(ic)%gtype=="MET".or. &
        group_pep(ic)%gtype=="PHE".or.group_pep(ic)%gtype=="TRP") then
        ip=20
    elseif(group_pep(ic)%gtype=="ARG".or.group_pep(ic)%gtype=="LYS") then
        ip=30
    elseif(group_pep(ic)%gtype=="SER".or.group_pep(ic)%gtype=="THR".or.group_pep(ic)%gtype=="ASN".or.group_pep(ic)%gtype=="GLN".or. &
        group_pep(ic)%gtype=="HIE") then
        ip=40
    elseif(group_pep(ic)%gtype=="PRO".or.group_pep(ic)%gtype=="ALA") then
        ip=50
    elseif(group_pep(ic)%gtype=="GLU".or.group_pep(ic)%gtype=="ASP".or.group_pep(ic)%gtype=="CYT".or.group_pep(ic)%gtype=="TYX") then
        ip=60
    elseif(group_pep(ic)%gtype=="NGLY") then
        ip=11
    elseif(group_pep(ic)%gtype=="NLEU".or.group_pep(ic)%gtype=="NVAL".or.group_pep(ic)%gtype=="NILE".or.group_pep(ic)%gtype=="NMET".or. &
        group_pep(ic)%gtype=="NPHE".or.group_pep(ic)%gtype=="NTRP") then
        ip=21
    elseif(group_pep(ic)%gtype=="NARG".or.group_pep(ic)%gtype=="NLYS") then
        ip=31
    elseif(group_pep(ic)%gtype=="NSER".or.group_pep(ic)%gtype=="NTHR".or.group_pep(ic)%gtype=="NASN".or.group_pep(ic)%gtype=="NGLN".or. &
        group_pep(ic)%gtype=="NHIE") then
        ip=41
    elseif(group_pep(ic)%gtype=="NPRO".or.group_pep(ic)%gtype=="NALA") then
        ip=51
    elseif(group_pep(ic)%gtype=="NGLU".or.group_pep(ic)%gtype=="NASP".or.group_pep(ic)%gtype=="NCYT".or.group_pep(ic)%gtype=="NTYX") then
        ip=61
    elseif(group_pep(ic)%gtype=="CGLY") then
        ip=12
    elseif(group_pep(ic)%gtype=="CLEU".or.group_pep(ic)%gtype=="CVAL".or.group_pep(ic)%gtype=="CILE".or.group_pep(ic)%gtype=="CMET".or. &
        group_pep(ic)%gtype=="CPHE".or.group_pep(ic)%gtype=="CTRP") then
        ip=22
    elseif(group_pep(ic)%gtype=="CARG".or.group_pep(ic)%gtype=="CLYS") then
        ip=32
    elseif(group_pep(ic)%gtype=="CSER".or.group_pep(ic)%gtype=="CTHR".or.group_pep(ic)%gtype=="CASN".or.group_pep(ic)%gtype=="CGLN".or. &
        group_pep(ic)%gtype=="CHIE") then
        ip=42
    elseif(group_pep(ic)%gtype=="CPRO".or.group_pep(ic)%gtype=="CALA") then
        ip=52
    elseif(group_pep(ic)%gtype=="CGLU".or.group_pep(ic)%gtype=="CASP".or.group_pep(ic)%gtype=="CCYT".or.group_pep(ic)%gtype=="CTYX") then
        ip=62
    endif

60 continue
    if(ip.eq.10) then
        aminoacid_name="GLY"
        goto 100
    elseif(ip.eq.20) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*6-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="LEU"
        elseif(ip1.eq.2) then
            aminoacid_name="VAL"
        elseif(ip1.eq.3) then
            aminoacid_name="ILE"
        elseif(ip1.eq.4) then
            aminoacid_name="MET"
        elseif(ip1.eq.5) then
            aminoacid_name="PHE"
        elseif(ip1.eq.6) then
            aminoacid_name="TRP"
        endif
    elseif(ip.eq.30) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*2-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="ARG"
        elseif(ip1.eq.2) then
            aminoacid_name="LYS"
        endif
    elseif(ip.eq.40) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*5-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="SER"
        elseif(ip1.eq.2) then
            aminoacid_name="THR"
        elseif(ip1.eq.3) then
            aminoacid_name="ASN"
        elseif(ip1.eq.4) then
            aminoacid_name="GLN"
        elseif(ip1.eq.5) then
            aminoacid_name="HIE"
        endif
    elseif(ip.eq.50) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*2-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="PRO"
        elseif(ip1.eq.2) then
            aminoacid_name="ALA"
        endif
    elseif(ip.eq.60) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*4-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="GLU"
        elseif(ip1.eq.2) then
            aminoacid_name="ASP"
        elseif(ip1.eq.3) then
            aminoacid_name="TYX"
        elseif(ip1.eq.4) then
            aminoacid_name="CYT"
        endif
    elseif(ip.eq.11) then
        aminoacid_name="NGLY"
        goto 100
    elseif(ip.eq.21) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*6-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NLEU"
        elseif(ip1.eq.2) then
            aminoacid_name="NVAL"
        elseif(ip1.eq.3) then
            aminoacid_name="NILE"
        elseif(ip1.eq.4) then
            aminoacid_name="NMET"
        elseif(ip1.eq.5) then
            aminoacid_name="NPHE"
        elseif(ip1.eq.6) then
            aminoacid_name="NTRP"
        endif
    elseif(ip.eq.31) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*2-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NARG"
        elseif(ip1.eq.2) then
            aminoacid_name="NLYS"
        endif
    elseif(ip.eq.41) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*5-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NSER"
        elseif(ip1.eq.2) then
            aminoacid_name="NTHR"
        elseif(ip1.eq.3) then
            aminoacid_name="NASN"
        elseif(ip1.eq.4) then
            aminoacid_name="NGLN"
        elseif(ip1.eq.5) then
            aminoacid_name="NHIE"
        endif
    elseif(ip.eq.51) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*2-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NPRO"
        elseif(ip1.eq.2) then
            aminoacid_name="NALA"
        endif
    elseif(ip.eq.61) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*4-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NGLU"
        elseif(ip1.eq.2) then
            aminoacid_name="NASP"
        elseif(ip1.eq.3) then
            aminoacid_name="NTYX"
        elseif(ip1.eq.4) then
            aminoacid_name="NCYT"
        endif
    elseif(ip.eq.12) then
        aminoacid_name="CGLY"
        goto 100
    elseif(ip.eq.22) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*6-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CLEU"
        elseif(ip1.eq.2) then
            aminoacid_name="CVAL"
        elseif(ip1.eq.3) then
            aminoacid_name="CILE"
        elseif(ip1.eq.4) then
            aminoacid_name="CMET"
        elseif(ip1.eq.5) then
            aminoacid_name="CPHE"
        elseif(ip1.eq.6) then
            aminoacid_name="CTRP"
        endif
    elseif(ip.eq.32) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*2-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CARG"
        elseif(ip1.eq.2) then
            aminoacid_name="CLYS"
        endif
    elseif(ip.eq.42) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*5-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CSER"
        elseif(ip1.eq.2) then
            aminoacid_name="CTHR"
        elseif(ip1.eq.3) then
            aminoacid_name="CASN"
        elseif(ip1.eq.4) then
            aminoacid_name="CGLN"
        elseif(ip1.eq.5) then
            aminoacid_name="CHIE"
        endif
    elseif(ip.eq.52) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*2-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CPRO"
        elseif(ip.eq.2) then
            aminoacid_name="CALA"
        endif
    elseif(ip.eq.62) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*4-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CGLU"
        elseif(ip1.eq.2) then
            aminoacid_name="CASP"
        elseif(ip1.eq.3) then
            aminoacid_name="CTYX"
        elseif(ip1.eq.4) then
            aminoacid_name="CCYT"
        endif
    endif
    if(aminoacid_name==group_pep(ic)%gtype) goto 60

elseif(ph_value.ge.10.5.and.ph_value.lt.12.5) then
    if(group_pep(ic)%gtype=="GLY") then
        ip=10
    elseif(group_pep(ic)%gtype=="LEU".or.group_pep(ic)%gtype=="VAL".or.group_pep(ic)%gtype=="ILE".or.group_pep(ic)%gtype=="MET".or. &
        group_pep(ic)%gtype=="PHE".or.group_pep(ic)%gtype=="TRP") then
        ip=20
    elseif(group_pep(ic)%gtype=="ARG") then
        ip=30
    elseif(group_pep(ic)%gtype=="SER".or.group_pep(ic)%gtype=="THR".or.group_pep(ic)%gtype=="ASN".or.group_pep(ic)%gtype=="GLN".or. &
        group_pep(ic)%gtype=="HIE".or.group_pep(ic)%gtype=="LYN") then
        ip=40
    elseif(group_pep(ic)%gtype=="PRO".or.group_pep(ic)%gtype=="ALA") then
        ip=50
    elseif(group_pep(ic)%gtype=="GLU".or.group_pep(ic)%gtype=="ASP".or.group_pep(ic)%gtype=="CYT".or.group_pep(ic)%gtype=="TYX") then
        ip=60
    elseif(group_pep(ic)%gtype=="NGLY") then
        ip=11
    elseif(group_pep(ic)%gtype=="NLEU".or.group_pep(ic)%gtype=="NVAL".or.group_pep(ic)%gtype=="NILE".or.group_pep(ic)%gtype=="NMET".or. &
        group_pep(ic)%gtype=="NPHE".or.group_pep(ic)%gtype=="NTRP") then
        ip=21
    elseif(group_pep(ic)%gtype=="NARG") then
        ip=31
    elseif(group_pep(ic)%gtype=="NSER".or.group_pep(ic)%gtype=="NTHR".or.group_pep(ic)%gtype=="NASN".or.group_pep(ic)%gtype=="NGLN".or. &
        group_pep(ic)%gtype=="NHIE".or.group_pep(ic)%gtype=="NLYN") then
        ip=41
    elseif(group_pep(ic)%gtype=="NPRO".or.group_pep(ic)%gtype=="NALA") then
        ip=51
    elseif(group_pep(ic)%gtype=="NGLU".or.group_pep(ic)%gtype=="NASP".or.group_pep(ic)%gtype=="NCYT".or.group_pep(ic)%gtype=="NTYX") then
        ip=61
    elseif(group_pep(ic)%gtype=="CGLY") then
        ip=12
    elseif(group_pep(ic)%gtype=="CLEU".or.group_pep(ic)%gtype=="CVAL".or.group_pep(ic)%gtype=="CILE".or.group_pep(ic)%gtype=="CMET".or. &
        group_pep(ic)%gtype=="CPHE".or.group_pep(ic)%gtype=="CTRP") then
        ip=22
    elseif(group_pep(ic)%gtype=="CARG") then
        ip=32
    elseif(group_pep(ic)%gtype=="CSER".or.group_pep(ic)%gtype=="CTHR".or.group_pep(ic)%gtype=="CASN".or.group_pep(ic)%gtype=="CGLN".or. &
        group_pep(ic)%gtype=="CHIE".or.group_pep(ic)%gtype=="CLYN") then
        ip=42
    elseif(group_pep(ic)%gtype=="CPRO".or.group_pep(ic)%gtype=="CALA") then
        ip=52
    elseif(group_pep(ic)%gtype=="CGLU".or.group_pep(ic)%gtype=="CASP".or.group_pep(ic)%gtype=="CCYT".or.group_pep(ic)%gtype=="CTYX") then
        ip=62
    endif

70 continue
    if(ip.eq.10) then
        aminoacid_name="GLY"
        goto 100
    elseif(ip.eq.20) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*6-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="LEU"
        elseif(ip1.eq.2) then
            aminoacid_name="VAL"
        elseif(ip1.eq.3) then
            aminoacid_name="ILE"
        elseif(ip1.eq.4) then
            aminoacid_name="MET"
        elseif(ip1.eq.5) then
            aminoacid_name="PHE"
        elseif(ip1.eq.6) then
            aminoacid_name="TRP"
        endif
    elseif(ip.eq.30) then
        aminoacid_name="ARG"
        goto 100
    elseif(ip.eq.40) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*6-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="SER"
        elseif(ip1.eq.2) then
            aminoacid_name="THR"
        elseif(ip1.eq.3) then
            aminoacid_name="ASN"
        elseif(ip1.eq.4) then
            aminoacid_name="GLN"
        elseif(ip1.eq.5) then
            aminoacid_name="HIE"
        elseif(ip1.eq.6) then
            aminoacid_name="LYN"
        endif
    elseif(ip.eq.50) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*2-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="PRO"
        elseif(ip1.eq.2) then
            aminoacid_name="ALA"
        endif
    elseif(ip.eq.60) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*4-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="GLU"
        elseif(ip1.eq.2) then
            aminoacid_name="ASP"
        elseif(ip1.eq.3) then
            aminoacid_name="TYX"
        elseif(ip1.eq.4) then
            aminoacid_name="CYT"
        endif
    elseif(ip.eq.11) then
        aminoacid_name="NGLY"
        goto 100
    elseif(ip.eq.21) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*6-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NLEU"
        elseif(ip1.eq.2) then
            aminoacid_name="NVAL"
        elseif(ip1.eq.3) then
            aminoacid_name="NILE"
        elseif(ip1.eq.4) then
            aminoacid_name="NMET"
        elseif(ip1.eq.5) then
            aminoacid_name="NPHE"
        elseif(ip1.eq.6) then
            aminoacid_name="NTRP"
        endif
    elseif(ip.eq.31) then
        aminoacid_name="NARG"
        goto 100
    elseif(ip.eq.41) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*6-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NSER"
        elseif(ip1.eq.2) then
            aminoacid_name="NTHR"
        elseif(ip1.eq.3) then
            aminoacid_name="NASN"
        elseif(ip1.eq.4) then
            aminoacid_name="NGLN"
        elseif(ip1.eq.5) then
            aminoacid_name="NHIE"
        elseif(ip1.eq.6) then
            aminoacid_name="NLYN"
        endif
    elseif(ip.eq.51) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*2-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NPRO"
        elseif(ip1.eq.2) then
            aminoacid_name="NALA"
        endif
    elseif(ip.eq.61) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*4-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NGLU"
        elseif(ip1.eq.2) then
            aminoacid_name="NASP"
        elseif(ip1.eq.3) then
            aminoacid_name="NTYX"
        elseif(ip1.eq.4) then
            aminoacid_name="NCYT"
        endif
    elseif(ip.eq.12) then
        aminoacid_name="CGLY"
        goto 100
    elseif(ip.eq.22) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*6-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CLEU"
        elseif(ip1.eq.2) then
            aminoacid_name="CVAL"
        elseif(ip1.eq.3) then
            aminoacid_name="CILE"
        elseif(ip1.eq.4) then
            aminoacid_name="CMET"
        elseif(ip1.eq.5) then
            aminoacid_name="CPHE"
        elseif(ip1.eq.6) then
            aminoacid_name="CTRP"
        endif
    elseif(ip.eq.32) then
        aminoacid_name="CARG"
        goto 100
    elseif(ip.eq.42) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*6-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CSER"
        elseif(ip1.eq.2) then
            aminoacid_name="CTHR"
        elseif(ip1.eq.3) then
            aminoacid_name="CASN"
        elseif(ip1.eq.4) then
            aminoacid_name="CGLN"
        elseif(ip1.eq.5) then
            aminoacid_name="CHIE"
        elseif(ip1.eq.6) then
            aminoacid_name="CLYN"
        endif
    elseif(ip.eq.52) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*2-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CPRO"
        elseif(ip.eq.2) then
            aminoacid_name="CALA"
        endif
    elseif(ip.eq.62) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*4-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CGLU"
        elseif(ip1.eq.2) then
            aminoacid_name="CASP"
        elseif(ip1.eq.3) then
            aminoacid_name="CTYX"
        elseif(ip1.eq.4) then
            aminoacid_name="CCYT"
        endif
    endif
    if(aminoacid_name==group_pep(ic)%gtype) goto 70

elseif(ph_value.ge.12.5) then
    if(group_pep(ic)%gtype=="GLY") then
        ip=10
    elseif(group_pep(ic)%gtype=="LEU".or.group_pep(ic)%gtype=="VAL".or.group_pep(ic)%gtype=="ILE".or.group_pep(ic)%gtype=="MET".or. &
        group_pep(ic)%gtype=="PHE".or.group_pep(ic)%gtype=="TRP") then
        ip=20
    elseif(group_pep(ic)%gtype=="SER".or.group_pep(ic)%gtype=="THR".or.group_pep(ic)%gtype=="ASN".or.group_pep(ic)%gtype=="GLN".or. &
        group_pep(ic)%gtype=="HIE".or.group_pep(ic)%gtype=="LYN".or.group_pep(ic)%gtype=="ARN") then
        ip=40
    elseif(group_pep(ic)%gtype=="PRO".or.group_pep(ic)%gtype=="ALA") then
        ip=50
    elseif(group_pep(ic)%gtype=="GLU".or.group_pep(ic)%gtype=="ASP".or.group_pep(ic)%gtype=="CYT".or.group_pep(ic)%gtype=="TYX") then
        ip=60
    elseif(group_pep(ic)%gtype=="NGLY") then
        ip=11
    elseif(group_pep(ic)%gtype=="NLEU".or.group_pep(ic)%gtype=="NVAL".or.group_pep(ic)%gtype=="NILE".or.group_pep(ic)%gtype=="NMET".or. &
        group_pep(ic)%gtype=="NPHE".or.group_pep(ic)%gtype=="NTRP") then
        ip=21
    elseif(group_pep(ic)%gtype=="NSER".or.group_pep(ic)%gtype=="NTHR".or.group_pep(ic)%gtype=="NASN".or.group_pep(ic)%gtype=="NGLN".or. &
        group_pep(ic)%gtype=="NHIE".or.group_pep(ic)%gtype=="NLYN".or.group_pep(ic)%gtype=="NARN") then
        ip=41
    elseif(group_pep(ic)%gtype=="NPRO".or.group_pep(ic)%gtype=="NALA") then
        ip=51
    elseif(group_pep(ic)%gtype=="NGLU".or.group_pep(ic)%gtype=="NASP".or.group_pep(ic)%gtype=="NCYT".or.group_pep(ic)%gtype=="NTYX") then
        ip=61
    elseif(group_pep(ic)%gtype=="CGLY") then
        ip=12
    elseif(group_pep(ic)%gtype=="CLEU".or.group_pep(ic)%gtype=="CVAL".or.group_pep(ic)%gtype=="CILE".or.group_pep(ic)%gtype=="CMET".or. &
        group_pep(ic)%gtype=="CPHE".or.group_pep(ic)%gtype=="CTRP") then
        ip=22
    elseif(group_pep(ic)%gtype=="CSER".or.group_pep(ic)%gtype=="CTHR".or.group_pep(ic)%gtype=="CASN".or.group_pep(ic)%gtype=="CGLN".or. &
        group_pep(ic)%gtype=="CHIE".or.group_pep(ic)%gtype=="CLYN".or.group_pep(ic)%gtype=="CARN") then
        ip=42
    elseif(group_pep(ic)%gtype=="CPRO".or.group_pep(ic)%gtype=="CALA") then
        ip=52
    elseif(group_pep(ic)%gtype=="CGLU".or.group_pep(ic)%gtype=="CASP".or.group_pep(ic)%gtype=="CCYT".or.group_pep(ic)%gtype=="CTYX") then
        ip=62
    endif

80 continue
    if(ip.eq.10) then
        aminoacid_name="GLY"
        goto 100
    elseif(ip.eq.20) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*6-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="LEU"
        elseif(ip1.eq.2) then
            aminoacid_name="VAL"
        elseif(ip1.eq.3) then
            aminoacid_name="ILE"
        elseif(ip1.eq.4) then
            aminoacid_name="MET"
        elseif(ip1.eq.5) then
            aminoacid_name="PHE"
        elseif(ip1.eq.6) then
            aminoacid_name="TRP"
        endif
    elseif(ip.eq.40) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*7-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="SER"
        elseif(ip1.eq.2) then
            aminoacid_name="THR"
        elseif(ip1.eq.3) then
            aminoacid_name="ASN"
        elseif(ip1.eq.4) then
            aminoacid_name="GLN"
        elseif(ip1.eq.5) then
            aminoacid_name="HIE"
        elseif(ip1.eq.6) then
            aminoacid_name="LYN"
        elseif(ip1.eq.7) then
            aminoacid_name="ARN"
        endif
    elseif(ip.eq.50) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*2-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="PRO"
        elseif(ip1.eq.2) then
            aminoacid_name="ALA"
        endif
    elseif(ip.eq.60) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*4-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="GLU"
        elseif(ip1.eq.2) then
            aminoacid_name="ASP"
        elseif(ip1.eq.3) then
            aminoacid_name="TYX"
        elseif(ip1.eq.4) then
            aminoacid_name="CYT"
        endif
    elseif(ip.eq.11) then
        aminoacid_name="NGLY"
        goto 100
    elseif(ip.eq.21) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*6-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NLEU"
        elseif(ip1.eq.2) then
            aminoacid_name="NVAL"
        elseif(ip1.eq.3) then
            aminoacid_name="NILE"
        elseif(ip1.eq.4) then
            aminoacid_name="NMET"
        elseif(ip1.eq.5) then
            aminoacid_name="NPHE"
        elseif(ip1.eq.6) then
            aminoacid_name="NTRP"
        endif
    elseif(ip.eq.41) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*7-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NSER"
        elseif(ip1.eq.2) then
            aminoacid_name="NTHR"
        elseif(ip1.eq.3) then
            aminoacid_name="NASN"
        elseif(ip1.eq.4) then
            aminoacid_name="NGLN"
        elseif(ip1.eq.5) then
            aminoacid_name="NHIE"
        elseif(ip1.eq.6) then
            aminoacid_name="NLYN"
        elseif(ip1.eq.7) then
            aminoacid_name="NARN"
        endif
    elseif(ip.eq.51) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*2-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NPRO"
        elseif(ip1.eq.2) then
            aminoacid_name="NALA"
        endif
    elseif(ip.eq.61) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*4-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="NGLU"
        elseif(ip1.eq.2) then
            aminoacid_name="NASP"
        elseif(ip1.eq.3) then
            aminoacid_name="NTYX"
        elseif(ip1.eq.4) then
            aminoacid_name="NCYT"
        endif
    elseif(ip.eq.12) then
        aminoacid_name="CGLY"
        goto 100
    elseif(ip.eq.22) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*6-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CLEU"
        elseif(ip1.eq.2) then
            aminoacid_name="CVAL"
        elseif(ip1.eq.3) then
            aminoacid_name="CILE"
        elseif(ip1.eq.4) then
            aminoacid_name="CMET"
        elseif(ip1.eq.5) then
            aminoacid_name="CPHE"
        elseif(ip1.eq.6) then
            aminoacid_name="CTRP"
        endif
    elseif(ip.eq.42) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*7-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CSER"
        elseif(ip1.eq.2) then
            aminoacid_name="CTHR"
        elseif(ip1.eq.3) then
            aminoacid_name="CASN"
        elseif(ip1.eq.4) then
            aminoacid_name="CGLN"
        elseif(ip1.eq.5) then
            aminoacid_name="CHIE"
        elseif(ip1.eq.6) then
            aminoacid_name="CLYN"
        elseif(ip1.eq.7) then
            aminoacid_name="CARN"
        endif
    elseif(ip.eq.52) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*2-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CPRO"
        elseif(ip.eq.2) then
            aminoacid_name="CALA"
        endif
    elseif(ip.eq.62) then
        call RANDOM_NUMBER(ran2)
        ip1=int(ran2*4-1.0e-3)+1
        if(ip1.eq.1) then
            aminoacid_name="CGLU"
        elseif(ip1.eq.2) then
            aminoacid_name="CASP"
        elseif(ip1.eq.3) then
            aminoacid_name="CTYX"
        elseif(ip1.eq.4) then
            aminoacid_name="CCYT"
        endif
    endif
    if(aminoacid_name==group_pep(ic)%gtype) goto 80

endif

100 continue

return
end subroutine mc_choose_aminoacid

end module utilities