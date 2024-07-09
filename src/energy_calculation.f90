module energy_calculation

use datatypes
use sys_vars
use math
use utilities
use database
use surface_area

contains

subroutine self_interaction_energy(res_i, GB_flag, pep_vdw, pep_ele, pep_gb_iso, pep_gb_sys)
implicit none
integer                         :: k, res_i, flag, GB_flag, atom_i, atom_j
double precision                :: pep_vdw, pep_ele, pep_gb_iso, pep_gb_sys, gb, vdw, ele
real                            :: rij

pep_vdw=0; pep_ele=0; pep_gb_iso=0; pep_gb_sys=0 !initialize energies to 0
do atom_i=1,atoms_per_res(res_i)

    !get self-interaction energy for the atom
    if (GB_flag.eq.1) then
        call gbcontri(0., alpha_pep_iso(atom_i, res_i), alpha_pep_iso(atom_i, res_i), charge_pep(atom_i, res_i), charge_pep(atom_i, res_i),&
        dielecons_pep(atom_i, res_i), dielecons_pep(atom_i, res_i), gb)
        pep_gb_iso = pep_gb_iso + gb/2.0

        call gbcontri(0., alpha_pep_sys(atom_i, res_i), alpha_pep_sys(atom_i, res_i), charge_pep(atom_i, res_i), charge_pep(atom_i, res_i),&
        dielecons_pep(atom_i, res_i), dielecons_pep(atom_i, res_i), gb)
        pep_gb_sys = pep_gb_sys + gb/2.0
    endif

    !get interaction of atom with all atoms in residue
    do atom_j=atom_i+1, atoms_per_res(res_i)
        rij = rij_matrix_pep(atom_j, atom_i, res_i, res_i)

        !calculate polar solvation energy
        if (GB_flag.eq.1) then
            !isolated peptide energy
            call gbcontri(rij,alpha_pep_iso(atom_i, res_i),alpha_pep_iso(atom_j, res_i),charge_pep(atom_i, res_i), charge_pep(atom_j, res_i),&
            dielecons_pep(atom_i, res_i), dielecons_pep(atom_j, res_i),gb)
            pep_gb_iso = pep_gb_iso + gb

            !peptide energy when bound to receptor
            call gbcontri(rij,alpha_pep_sys(atom_i, res_i),alpha_pep_sys(atom_j, res_i),charge_pep(atom_i, res_i), charge_pep(atom_j, res_i),&
            dielecons_pep(atom_i, res_i), dielecons_pep(atom_j, res_i),gb)
            pep_gb_sys = pep_gb_sys + gb
        endif
            
        !now move on to van der Waals and electrostatics calculations 

        !check if atom pair are in excluded list
        do k=1,numex_pep(atom_i, res_i)
            if(atom_j.eq.inb_pep(1,k,atom_i, res_i) .and. res_i.eq.inb_pep(2,k,atom_i, res_i)) goto 40
        enddo

        !check if atom in 1-4 connection
        flag=0
        do k=1,numex4_pep(atom_i, res_i)
            if(atom_j.eq.inb4_pep(1,k,atom_i, res_i) .and. res_i.eq.inb4_pep(1,k,atom_i, res_i)) then
                flag=1
                exit
            endif
        enddo

        !calculate VDW and ELE energy
        call vdwcontri(rij, epsion_pep(atom_i, res_i), epsion_pep(atom_j, res_i), r_pep(atom_i, res_i), r_pep(atom_j, res_i), vdw)
        call elecontri(rij, charge_pep(atom_i, res_i), charge_pep(atom_j, res_i), dielecons_pep(atom_i, res_i), dielecons_pep(atom_j, res_i), ele)

        !scale if 1-4 connection
        if(flag==1) then
            vdw=vdw/vdw14_coeff
            ele=ele/ele14_coeff
        endif
        
        !store energy
        pep_vdw = pep_vdw + vdw
        pep_ele = pep_ele + ele

        40 continue
    enddo
enddo

end subroutine self_interaction_energy

subroutine pairwise_res_energy(res_i, res_j, GB_flag, vdw_bind, ele_bind, gb_iso_bind, gb_sys_bind)
implicit none
integer                         :: k, res_i, res_j, atom_i, atom_j, flag, GB_flag
double precision                :: vdw, ele, gb, vdw_bind, ele_bind, gb_iso_bind, gb_sys_bind
real                            :: rij

gb_iso_bind=0.0; gb_sys_bind = 0.0; vdw_bind=0.0; ele_bind=0.0
do atom_i=1,atoms_per_res(res_i)
    do atom_j=1, atoms_per_res(res_j)
        rij = rij_matrix_pep(atom_j, atom_i, res_j, res_i)

        !calculate polar solvation energy
        if (GB_flag.eq.1) then
            !isolated peptide energy
            call gbcontri(rij,alpha_pep_iso(atom_i, res_i),alpha_pep_iso(atom_j, res_j),charge_pep(atom_i, res_i), &
            charge_pep(atom_j, res_j),dielecons_pep(atom_i, res_i), dielecons_pep(atom_j, res_j),gb)
            gb_iso_bind = gb_iso_bind + gb

            !peptide energy when bound to receptor
            call gbcontri(rij,alpha_pep_sys(atom_i, res_i),alpha_pep_sys(atom_j, res_j),charge_pep(atom_i, res_i), &
            charge_pep(atom_j, res_j),dielecons_pep(atom_i, res_i), dielecons_pep(atom_j, res_j),gb)
            gb_sys_bind = gb_sys_bind + gb
        endif

        !calculate van der Waals and electrostatic interactions
        !check if atom pair are in excluded list
        do k=1,numex_pep(atom_i, res_i)
            if(atom_j.eq.inb_pep(1, k, atom_i, res_i) .and. res_j.eq.inb_pep(2, k, atom_i, res_i)) goto 40
        enddo
        flag=0

        !check if atom in 1-4 connection
        do k=1,numex4_pep(atom_i, res_i)
            if(atom_j.eq.inb4_pep(1, k, atom_i, res_i) .and. res_j.eq.inb4_pep(2, k, atom_i, res_i)) then
                flag=1
                exit
            endif
        enddo

        !calculate VDW and ELE energy
        call vdwcontri(rij, epsion_pep(atom_i, res_i), epsion_pep(atom_j, res_j), r_pep(atom_i, res_i), r_pep(atom_j, res_j), vdw)
        call elecontri(rij, charge_pep(atom_i, res_i), charge_pep(atom_j, res_j), dielecons_pep(atom_i, res_i), dielecons_pep(atom_j, res_j), ele)

        !scale if 1-4 connection
        if(flag==1) then
            vdw=vdw/vdw14_coeff
            ele=ele/ele14_coeff
        endif
        
        !store energy
        vdw_bind = vdw_bind + vdw
        ele_bind = ele_bind + ele
        40 continue
    enddo
enddo

end subroutine pairwise_res_energy

subroutine res_receptor_energy(res_i, GB_flag, vdw_bind, ele_bind, gb_bind)
!GB energy calculation for this seems good
implicit none
integer                         :: atom_i, res_i, atom_j, GB_flag
double precision                :: vdw, ele, gb, gb_bind, vdw_bind, ele_bind
real                            :: rij

vdw_bind = 0; ele_bind = 0; gb_bind = 0 !initialize energies to 0 

do atom_i=1, atoms_per_res(res_i)
    do atom_j=1, atom_num_rec
        rij = rij_matrix_sys(atom_j, atom_i, res_i)

        !calculate van der Waals and electrostatics energy
        call vdwcontri(rij, epsion_pep(atom_i, res_i), epsion_rec(atom_j), r_pep(atom_i, res_i), r_rec(atom_j), vdw)
        call elecontri(rij, charge_pep(atom_i, res_i), charge_rec(atom_j), dielecons_pep(atom_i, res_i), dielecons_rec(atom_j), ele)

        !store energy of peptide-peptide or peptide-receptor interactions in bound state
        vdw_bind = vdw_bind + vdw
        ele_bind = ele_bind + ele

        !calculate peptide-receptor polar solvation energy
        if (GB_flag.eq.1) then
            call gbcontri(rij, alpha_pep_sys(atom_i, res_i), alpha_rec_sys(atom_j), charge_pep(atom_i, res_i), charge_rec(atom_j), dielecons_pep(atom_i, res_i), dielecons_rec(atom_j), gb)
            gb_bind = gb_bind + gb
        endif
    enddo
enddo
end subroutine res_receptor_energy

subroutine sidechain_dihedral_energy(res_i, dihedral, num_sc_angles, dihedral_energy)
implicit none
type(dihedralparameters)   :: dihedral
integer                    :: i, j, num_sc_angles, ip, jp, kp, lp, res_i
real                       :: p1(3), p2(3), p3(3), p4(3), angle
double precision           :: potential, dihedral_energy

dihedral_energy=0
do i=1, num_sc_angles
    ip=dihedral%iph(i)+res_starts(res_i)-1; jp=dihedral%jph(i)+res_starts(res_i)-1
    kp=dihedral%kph(i)+res_starts(res_i)-1; lp=dihedral%lph(i)+res_starts(res_i)-1
    p1 = mdcrd_pep(:,ip, res_i)
    p2 = mdcrd_pep(:,jp, res_i)
    p3 = mdcrd_pep(:,kp, res_i)
    p4 = mdcrd_pep(:,lp, res_i)
    call calc_dihedral(p1,p2,p3,p4,angle)
    do j=1, dihedral%multiply(i)
        potential=dihedral%pk(j,i)*(1+cosd(dihedral%pn(j,i)*angle-dihedral%phase(j,i)))
        dihedral_energy=dihedral_energy+potential
    enddo
enddo
dihedral_energy = dihedral_energy*dihedral_weighting_factor
return

end subroutine sidechain_dihedral_energy

subroutine calc_rec_iso_gb_energy
implicit none
integer             :: atom_i, atom_j
double precision    :: GB
real                :: rij

!load generalized born parameters
call effgbradi_iso_rec

!go through all atom pairs and calculating GB energy
rec_gb = 0
do atom_i=1, atom_num_rec
    !first get self-interaction energy
    call gbcontri(0.,alpha_rec_iso(atom_i),alpha_rec_iso(atom_i),charge_rec(atom_i), charge_rec(atom_i),dielecons_rec(atom_i), dielecons_rec(atom_i), GB)
    rec_gb = rec_gb + GB/2.0

    !now calculate interaction with all downstream atoms in receptor
    do atom_j=atom_i+1, atom_num_rec
        rij = rij_matrix_rec(atom_j,atom_i)
        call gbcontri(rij, alpha_rec_iso(atom_i), alpha_rec_iso(atom_j), charge_rec(atom_i), charge_rec(atom_j), dielecons_rec(atom_i), dielecons_rec(atom_j), GB)
        rec_gb = rec_gb + GB
    enddo
enddo

end subroutine calc_rec_iso_gb_energy

subroutine calc_rec_sys_gb_energy(GB_complex)
!GB calculation for this seems to be good!
implicit none
integer             :: atom_i, atom_j
double precision    :: GB, GB_complex
real                :: rij

!go through all atom pairs and calculating GB energy
GB_complex = 0
do atom_i=1, atom_num_rec
    !first get self-interaction energy
    call gbcontri(0.,alpha_rec_sys(atom_i),alpha_rec_sys(atom_i),charge_rec(atom_i), charge_rec(atom_i),dielecons_rec(atom_i), dielecons_rec(atom_i),GB)
    GB_complex = GB_complex + GB/2.0
    !now calculate interaction with all downstream atoms in receptor
    do atom_j=atom_i+1, atom_num_rec
        rij = rij_matrix_rec(atom_j,atom_i)
        call gbcontri(rij,alpha_rec_sys(atom_i),alpha_rec_sys(atom_j),charge_rec(atom_i), charge_rec(atom_j), dielecons_rec(atom_i), dielecons_rec(atom_j),GB)
        GB_complex = GB_complex + GB
    enddo
enddo

end subroutine calc_rec_sys_gb_energy

subroutine backbone_energy()
!Energy terms: VDW and ELE
!Group 1 in energy calculations: all atoms in residue i backbone (i ranging from 1 to pep_res)
!Group 2 in energy calculations: all atoms in residue j backbone (j ranging from 1 to pep_res, j not equal to i)
implicit none
!integer                       :: flag1, natom, atom_offset, res_i, res_j, i_index, j_index, cg_flags(pep_res)
!integer                       :: i, j, k
!double precision              :: energy, vdw, ele, totdevdw, totdeele
!real                          :: rij !, ecutoff=1000
!type(groupdetails)            :: group_pep(pep_res)
!type(energyparameters)        :: para_pep(pep_res)


write(*,*) "This subroutine not currently working, please do not use! Terminating program"
stop

!go through all atom pairs in peptide and calculate VDW and ELE energy
!do res_i=1, pep_res
!    
!    do j=i+1, atom_num_pep
!        !book keeping on indices into distance matrix
!        if (j.gt.res_ends(res_j)) res_j = res_j + 1
!        j_index = j - res_starts(res_j) + 1
!
!        !check if atom pair are in excluded list
!        do k=1,numex_pep(atomid_pep(i, res_i))
!            if(atomid_pep(j, res_j).eq.inb_pep(k,atomid_pep(i, res_i))) goto 50
!        enddo
!        flag1=0
!
!        !check if atom in 1-4 connection
!        do k=1,numex4_pep(atomid_pep(i, res_i))
!            if(atomid_pep(j, res_j).eq.inb4_pep(k,atomid_pep(i, res_i))) then
!                flag1=1
!                goto 55
!            endif
!        enddo
!
!    55  continue !go here when atom pair i-j are in numex4_pep list 
!        !call DIST_SQUARE(mdcrd_pep(:,i), mdcrd_pep(:,j), rij)
!        rij = rij_matrix_pep(j_index, i_index, res_j, res_i)
!
!        !calculate VDW and ELE energy
!        call vdwcontri(rij,epsion_pep(i, res_i), epsion_pep(j, res_j),r_pep(i, res_i),r_pep(j, res_j), vdw)
!        call elecontri(rij,charge_pep(i, res_i),charge_pep(j, res_j),dielecons_pep(i, res_i),dielecons_pep(j, res_j),ele)!
!
!        !scale if 1-4 connection
!        if(flag1==1) then
!            vdw=vdw/vdw14_coeff
!            ele=ele/ele14_coeff
!        endif
!
!        !store energy of peptide-peptide or peptide-receptor interactions in bound state
!        totdevdw = totdevdw + vdw
!        totdeele = totdeele + ele
!    50 continue !go here when atom pair i-j are in exclusion list (i.e. they are 3 or fewer atoms away)
!    enddo
!enddo

!energy = totdeele + totdevdw

return
end subroutine backbone_energy

subroutine sidechain_energy(res_i, res_j, group_pep, para_pep, calc_type, cg_flags, ele_flag, energy)
!Energy terms: VDW
!Group 1 in energy calculations: all atoms in sidechain of residue i
!Group 2 in energy calculations: all atoms in receptor (calc_type=0); entire system (calc_type=1); or all atoms in residue j (calc_type=2)
implicit none
integer                    :: i, j, k, calc_type, ele_flag, i_index, j_index, cg_flags(pep_res)
integer                    :: res_i, res_j
double precision           :: energy,  vdw, ele
real                       :: rij, epsion, r0, acoeff, bcoeff
type(groupdetails)         :: group_pep(pep_res)
type(energyparameters)     :: para_pep(pep_res)

!load system parameters and coordinates
cg_flags = 0
call load_peptide_info(group_pep, para_pep)

energy=0.0 !initalize energy to 0

if(calc_type==0) then !get interactions of sidechain with all receptor atoms
    do i=1, atoms_per_res(res_i)
        do j=1, atom_num_rec
            rij = rij_matrix_sys(j, i, res_i)
    
            !calculate the VDW and ELE energy
            call vdwcontri(rij,epsion_pep(i, res_i), epsion_rec(j),r_pep(i, res_i),r_rec(j), vdw)
            energy=energy+vdw
            if (ele_flag.eq.1) then
                call elecontri(rij,charge_pep(i, res_i),charge_rec(j),dielecons_pep(i, res_i),dielecons_rec(j),ele)
                energy=energy+ele
            endif
        enddo
    enddo

elseif(calc_type==1) then !get interactions of sidechain with all other residues in system
    do i=1, group_pep(res_i)%cnum2
        i_index = i + group_pep(res_i)%cnum1
        !get interactions with peptide
        do k=1, pep_res
            if (k.eq.res_i.or.k.eq.res_j) goto 40
            do j=1, group_pep(k)%cnum1
                if(k==(res_i+1)) goto 50
                j_index = j
                if(group_pep(res_i)%gtype=="PRO".or.group_pep(res_i)%gtype=="NPRO".or.group_pep(res_i)%gtype=="CPRO") then
                    if(k==(res_i-1).and.group_pep(res_i)%atype2(i)=="CD".and.group_pep(k)%atype1(j)=="CA") then
                        goto 50
                    endif
                endif
                !call DIST_SQUARE(group_pep(res_i)%coo2(:,i), group_pep(k)%coo1(:,j), rij)
                rij = rij_matrix_pep(j_index, i_index, res_j, res_i)
                epsion=sqrt(para_pep(res_i)%epsion2(i)*para_pep(k)%epsion1(j))
                r0=para_pep(res_i)%r2(i)+para_pep(k)%r1(j)
                acoeff=epsion*(r0**12)
                bcoeff=epsion*2*(r0**6)
                if(group_pep(res_i)%gtype=="PRO".or.group_pep(res_i)%gtype=="NPRO".or.group_pep(res_i)%gtype=="CPRO".or. &
                group_pep(res_i+1)%gtype=="PRO".or.group_pep(res_i+1)%gtype=="NPRO".or.group_pep(res_i+1)%gtype=="CPRO") then
                    if(rij<10.0) rij=12.25
                endif
                vdw=acoeff/(rij**6)-bcoeff/(rij**3)
                energy=energy + weighting_factor*vdw

                if (ele_flag.eq.1) then
                    call elecontri(rij,para_pep(res_i)%charge2(i),para_pep(res_j)%charge1(j), &
                                        para_pep(res_i)%dielecons2(i),para_pep(res_j)%dielecons1(j),ele)
                    energy=energy + weighting_factor*ele
                endif
                50 continue
            enddo

            do j=1, group_pep(k)%cnum2
                j_index = j+group_pep(k)%cnum1
                !call DIST_SQUARE(group_pep(res_i)%coo2(:,i), group_pep(k)%coo2(:,j), rij)
                rij = rij_matrix_pep(j_index, i_index, res_j, res_i)
                epsion=sqrt(para_pep(res_i)%epsion2(i)*para_pep(k)%epsion2(j))
                r0=para_pep(res_i)%r2(i)+para_pep(k)%r2(j)
                acoeff=epsion*(r0**12)
                bcoeff=epsion*2*(r0**6)
                if(group_pep(res_i)%gtype=="PRO".or.group_pep(res_i)%gtype=="NPRO".or.group_pep(res_i)%gtype=="CPRO".or. &
                group_pep(res_i+1)%gtype=="PRO".or.group_pep(res_i+1)%gtype=="NPRO".or.group_pep(res_i+1)%gtype=="CPRO") then
                    if(rij<10.0) rij=12.25
                endif
                vdw=acoeff/(rij**6)-bcoeff/(rij**3)
                energy=energy + weighting_factor*vdw
                if (ele_flag.eq.1) then
                    call elecontri(rij,para_pep(res_i)%charge2(i),para_pep(res_j)%charge2(j), &
                                        para_pep(res_i)%dielecons2(i),para_pep(res_j)%dielecons2(j),ele)
                    energy=energy + weighting_factor*ele
                endif
            enddo

            do j=1, group_pep(k)%cnum3
                j_index = j+group_pep(k)%cnum1+group_pep(k)%cnum2
                if(k==(res_i-1)) goto 60
                if((group_pep(res_i)%gtype=="PRO".or.group_pep(res_i)%gtype=="NPRO".or.group_pep(res_i)%gtype=="CPRO").and.k==(res_i-1)) then
                    if(group_pep(res_i)%atype2(i)=="CD") then
                        goto 60
                    elseif(group_pep(res_i)%atype2(i)=="HD2".or.group_pep(res_i)%atype2(i)=="HD3") then
                        if(group_pep(k)%atype3(j)=="C") goto 60
                    endif
                endif
                !call DIST_SQUARE(group_pep(res_i)%coo2(:,i), group_pep(k)%coo3(:,j), rij)
                rij = rij_matrix_pep(j_index, i_index, res_j, res_i)
                epsion=sqrt(para_pep(res_i)%epsion2(i)*para_pep(k)%epsion3(j))
                r0=para_pep(res_i)%r2(i)+para_pep(k)%r3(j)
                acoeff=epsion*(r0**12)
                bcoeff=epsion*2*(r0**6)
                if(group_pep(res_i)%gtype=="PRO".or.group_pep(res_i)%gtype=="NPRO".or.group_pep(res_i)%gtype=="CPRO".or. &
                group_pep(res_i+1)%gtype=="PRO".or.group_pep(res_i+1)%gtype=="NPRO".or.group_pep(res_i+1)%gtype=="CPRO") then
                    if(rij<10.0) rij=12.25
                endif
                vdw=acoeff/(rij**6)-bcoeff/(rij**3)
                energy=energy + weighting_factor*vdw
                if (ele_flag.eq.1) then
                    call elecontri(rij,para_pep(res_i)%charge2(i),para_pep(res_j)%charge3(j), &
                                        para_pep(res_i)%dielecons2(i),para_pep(res_j)%dielecons3(j),ele)
                    energy=energy +  weighting_factor*ele
                endif
            60 continue
            enddo
40      continue
        enddo
        
        !now get interactions with receptor
        do j=1, atom_num_rec
            rij = rij_matrix_sys(j, i_index, res_i)
    
            !calculate the VDW and ELE energy
            call vdwcontri(rij,epsion_pep(i, res_i), epsion_rec(j),r_pep(i, res_i),r_rec(j), vdw)
            energy=energy+vdw

            if (ele_flag.eq.1) then
                call elecontri(rij,charge_pep(i, res_i),charge_rec(j),dielecons_pep(i, res_i),dielecons_rec(j),ele)
                energy=energy+ele
            endif
        enddo
    enddo

elseif(calc_type==2) then !get interactions of sidechain with all atoms in residue j
    do i=1, group_pep(res_i)%cnum2
        i_index = i + group_pep(res_i)%cnum1
        k=res_j
        do j=1, group_pep(k)%cnum1
            j_index = j
            if(k==(res_i+1)) goto 70
            if(group_pep(res_i)%gtype=="PRO".or.group_pep(res_i)%gtype=="NPRO".or.group_pep(res_i)%gtype=="CPRO") then
                if(k==(res_i-1).and.group_pep(res_i)%atype2(i)=="CD".and.group_pep(k)%atype1(j)=="CA") then
                    goto 70
                endif
            endif
            !call DIST_SQUARE(group_pep(res_i)%coo2(:,i), group_pep(k)%coo1(:,j), rij)
            rij = rij_matrix_pep(j_index, i_index, res_j, res_i)
            epsion=sqrt(para_pep(res_i)%epsion2(i)*para_pep(k)%epsion1(j))
            r0=para_pep(res_i)%r2(i)+para_pep(k)%r1(j)
            acoeff=epsion*(r0**12)
            bcoeff=epsion*2*(r0**6)
            if(group_pep(res_i)%gtype=="PRO".or.group_pep(res_i)%gtype=="NPRO".or.group_pep(res_i)%gtype=="CPRO".or. &
            group_pep(res_i+1)%gtype=="PRO".or.group_pep(res_i+1)%gtype=="NPRO".or.group_pep(res_i+1)%gtype=="CPRO") then
            if(rij<10.0) rij=12.25
            endif
            vdw=acoeff/(rij**6)-bcoeff/(rij**3)
            energy=energy+vdw
            if (ele_flag.eq.1) then
                call elecontri(rij,para_pep(res_i)%charge2(i),para_pep(res_j)%charge1(j), &
                                    para_pep(res_i)%dielecons2(i),para_pep(res_j)%dielecons1(j),ele)
                energy=energy+ele
            endif
            70 continue
        enddo

        do j=1, group_pep(k)%cnum2
            j_index = j + group_pep(res_j)%cnum1
            !call DIST_SQUARE(group_pep(res_i)%coo2(:,i), group_pep(k)%coo2(:,j), rij)
            rij = rij_matrix_pep(j_index, i_index, res_j, res_i)
            epsion=sqrt(para_pep(res_i)%epsion2(i)*para_pep(k)%epsion2(j))
            r0=para_pep(res_i)%r2(i)+para_pep(k)%r2(j)
            acoeff=epsion*(r0**12)
            bcoeff=epsion*2*(r0**6)
            if(group_pep(res_i)%gtype=="PRO".or.group_pep(res_i)%gtype=="NPRO".or.group_pep(res_i)%gtype=="CPRO".or. &
            group_pep(res_i+1)%gtype=="PRO".or.group_pep(res_i+1)%gtype=="NPRO".or.group_pep(res_i+1)%gtype=="CPRO") then
                if(rij<10.0) rij=12.25
            endif
            vdw=acoeff/(rij**6)-bcoeff/(rij**3)
            energy=energy+vdw
            if (ele_flag.eq.1) then
                call elecontri(rij,para_pep(res_i)%charge2(i),para_pep(res_j)%charge2(j), &
                                    para_pep(res_i)%dielecons2(i),para_pep(res_j)%dielecons2(j),ele)
                energy=energy+ele
            endif
        enddo

        do j=1, group_pep(k)%cnum3
            j_index = j + group_pep(res_j)%cnum1 + group_pep(res_j)%cnum2
            if(k==(res_i-1)) goto 80
            if((group_pep(res_i)%gtype=="PRO".or.group_pep(res_i)%gtype=="NPRO".or.group_pep(res_i)%gtype=="CPRO").and.k==(res_i-1)) then
                if(group_pep(res_i)%atype2(i)=="CD") then
                    goto 80
                elseif(group_pep(res_i)%atype2(i)=="HD2".or.group_pep(res_i)%atype2(i)=="HD3") then
                    if(group_pep(k)%atype3(j)=="C") goto 80
                endif
            endif
            !call DIST_SQUARE(group_pep(res_i)%coo2(:,i), group_pep(k)%coo3(:,j), rij)
            rij = rij_matrix_pep(j_index, i_index, res_j, res_i)
            epsion=sqrt(para_pep(res_i)%epsion2(i)*para_pep(k)%epsion3(j))
            r0=para_pep(res_i)%r2(i)+para_pep(k)%r3(j)
            acoeff=epsion*(r0**12)
            bcoeff=epsion*2*(r0**6)
            if(group_pep(res_i)%gtype=="PRO".or.group_pep(res_i)%gtype=="NPRO".or.group_pep(res_i)%gtype=="CPRO".or. &
            group_pep(res_i+1)%gtype=="PRO".or.group_pep(res_i+1)%gtype=="NPRO".or.group_pep(res_i+1)%gtype=="CPRO") then
                if(rij<10.0) rij=12.25
            endif
            vdw=acoeff/(rij**6)-bcoeff/(rij**3)
            energy=energy+vdw
            if (ele_flag.eq.1) then
                call elecontri(rij,para_pep(res_i)%charge2(i),para_pep(res_j)%charge3(j), &
                                    para_pep(res_i)%dielecons2(i),para_pep(res_j)%dielecons3(j),ele)
                energy=energy+ele
            endif
            80 continue
        enddo
    enddo
endif

return
end subroutine sidechain_energy
    
subroutine bindingenergy_sidechainoptimization(group_pep, para_pep, res_i, dihedral, num_sc_angles, GB_flag, CG_flag, CG_flags, binding_energy)
implicit none
integer                         :: res_i, res_j, CG_flag, CG_flags(pep_res)
integer                         :: num_sc_angles
integer                         :: GB_flag
double precision                :: vdw, ele, dihedral_energy
double precision                :: binding_vdw, binding_ele, binding_gb, binding_energy, pep_vdw, pep_ele, pep_gb, sys_gb, pep_energy, gb_iso, gb_sys
type(groupdetails)              :: group_pep(pep_res)
type(dihedralparameters)        :: dihedral
type(energyparameters)          :: para_pep(pep_res)

!load system parameters and coordinates
if (CG_flag.eq.0) then
    call load_peptide_info(group_pep, para_pep)
else
    call load_peptide_info_cg(group_pep, para_pep, cg_flags)
endif

!load generalized born parameters, if including polar solvation energy in calculation
if (GB_flag.eq.1) then
    call effgbradi_complex
endif

!!$OMP PARALLEL PRIVATE(i, j, k, vdw, ele, gb, rij, flag1)
!!$OMP DO REDUCTION (+ : totdevdw, totdeele, totdegb, energy_per_res, rec_gb)

!get interactions of residue with peptide, going towards N-terminus
pep_ele=0.0; pep_gb=0.0; pep_vdw=0.0; pep_gb = 0; sys_gb = 0
binding_ele=0.0; binding_gb=0; binding_vdw=0.0
do res_j=1, res_i-1
    call pairwise_res_energy(res_i, res_j, GB_flag, vdw, ele, gb_iso, gb_sys)
    pep_vdw = pep_vdw + vdw
    pep_ele = pep_ele + ele
    pep_gb = pep_gb + gb_iso
    sys_gb = sys_gb + gb_sys
enddo

!get interaction of residue with itself
call self_interaction_energy(res_i, GB_flag, vdw, ele, gb_iso, gb_sys)
pep_vdw = pep_vdw + vdw
pep_ele = pep_ele + ele
pep_gb = pep_gb + gb_iso
sys_gb = sys_gb + gb_sys

!get interactions of residue with peptide, going towards C-terminus
do res_j=res_i+1, pep_res
    call pairwise_res_energy(res_i, res_j, GB_flag, vdw, ele, gb_iso, gb_sys)
    pep_vdw = pep_vdw + vdw
    pep_ele = pep_ele + ele
    pep_gb = pep_gb + gb_iso
    sys_gb = sys_gb + gb_sys
enddo

pep_energy = weighting_factor*(pep_vdw + pep_ele + pep_gb)

!get interaction of residue with receptor
call res_receptor_energy(res_i, GB_flag, binding_vdw, binding_ele, gb_sys)
sys_gb = sys_gb + gb_sys

!get dihedral energy (deactivating for now since its energy contribution is typically negligible)
!call sidechain_dihedral_energy(res_i, dihedral, num_sc_angles, dihedral_energy)
!dihedral_energy = dihedral_weighting_factor*dihedral_energy
dihedral_energy = 0

!get net binding energy
binding_gb = sys_gb - pep_gb
binding_energy = binding_vdw + binding_ele + binding_gb + pep_energy + dihedral_energy

return
end subroutine bindingenergy_sidechainoptimization

subroutine bindingenergy_full(group_pep, para_pep, score, delta_vdw, delta_ele, delta_gb, delta_np, pep_energy)
!calculates full energy of peptide (VDW and ELE interactions)
implicit none
integer                    :: res_i, res_j
double precision           :: vdw, ele, gb_iso, gb_sys
double precision           :: delta_vdw, delta_ele, delta_gb, delta_np
double precision           :: pep_vdw, pep_ele, pep_gb, pep_energy, comp_gb
double precision           :: score
type(groupdetails)         :: group_pep(:)
type(energyparameters)     :: para_pep(:)

!load peptide parameters and coordinates
call load_peptide_info(group_pep, para_pep)

!load generalized born parameters, if including polar solvation energy in calculation
if (solvation_flag.eq.1) call effgbradi_complex

!set all energy values to 0 initially
delta_vdw=0.0; delta_ele=0.0; delta_np = 0.0; !delta_gb doesn't need to be set to 0 here
pep_vdw = 0.0; pep_ele = 0.0; pep_gb = 0.0; pep_energy = 0; comp_gb=0.0; 

!calculate peptide-peptide interactions
do res_i=1, pep_res
    call self_interaction_energy(res_i, solvation_flag, vdw, ele, gb_iso, gb_sys)
    pep_vdw = pep_vdw + vdw
    pep_ele = pep_ele + ele
    pep_gb = pep_gb + gb_iso
    comp_gb = comp_gb + gb_sys
    do res_j=res_i+1, pep_res
        call pairwise_res_energy(res_i, res_j, solvation_flag, vdw, ele, gb_iso, gb_sys)
        pep_vdw = pep_vdw + vdw
        pep_ele = pep_ele + ele
        pep_gb = pep_gb + gb_iso
        comp_gb = comp_gb + gb_sys
    enddo
enddo
pep_energy = pep_vdw + pep_ele + pep_gb

!get interaction of peptide with receptor
do res_i=1, pep_res
    call res_receptor_energy(res_i, solvation_flag, vdw, ele, gb_sys)
    delta_vdw = delta_vdw + vdw
    delta_ele = delta_ele + ele
    comp_gb = comp_gb + gb_sys
enddo

!get receptor solvation energy
if (solvation_flag.eq.1) then
    call calc_rec_sys_gb_energy(gb_sys)
    comp_gb = comp_gb + gb_sys
endif

!get net binding energy
delta_gb = comp_gb - rec_gb - pep_gb
score = delta_vdw + delta_ele + delta_gb + delta_np + weighting_factor*pep_energy

return
end subroutine bindingenergy_full

subroutine bindingenergy_single(group_pep, para_pep, res_i, GB_flag, CG_flag, energy)
!get the binding energy of a single residue (res_i) with the receptor. No interactions between peptide different peptide residues calculated
!!NOTE: Need to fix this subroutine so it properly calculates energies for glycine, since interactions of the C-alpha hydrogen not calculated (not too urgent, I suspect energy of a single hydrogen won't be significant)
implicit none
integer                    :: GB_flag, CG_flag
integer                    :: res_i, cg_flags(pep_res)
double precision           :: pep_vdw, pep_ele, pep_gb_sys, pep_gb_iso !volume, surfarea
double precision           :: energy, peptide_energy
double precision           :: totdevdw, totdeele, totdegb
type(groupdetails)         :: group_pep(pep_res)
type(energyparameters)     :: para_pep(pep_res)

!load system parameters and coordinates
if (CG_flag.eq.0) then
    call load_peptide_info(group_pep, para_pep)
else
    cg_flags = 1
    cg_flags(res_i) = 0
    call load_peptide_info_cg(group_pep, para_pep, cg_flags)
endif

!load generalized born parameters, if including polar solvation energy in calculation
if (GB_flag.eq.1) call effgbradi_complex

call self_interaction_energy(res_i, GB_flag, pep_vdw, pep_ele, pep_gb_iso, pep_gb_sys)
peptide_energy = pep_vdw + pep_ele + pep_gb_iso

!get interaction of residue with receptor
call res_receptor_energy(res_i, GB_flag, totdevdw, totdeele, totdegb)
totdegb = totdegb + pep_gb_sys - pep_gb_iso

!call sidechain_dihedral_energy(res_i, dihedral, num_sc_angles, dihedral_energy)

energy = totdevdw + totdeele + totdegb + weighting_factor*(peptide_energy)

return
end subroutine bindingenergy_single

subroutine bindingenergy_pair(group_pep, para_pep, res_i, res_j, GB_flag, CG_flag, binding_energy, vdw)
implicit none
integer                    :: res_i, res_j, GB_flag, CG_flag, cg_flags(pep_res)
double precision           :: vdw, ele !, gb !volume, surfarea
double precision           :: binding_energy
double precision           :: gb_iso, gb_sys
type(groupdetails)         :: group_pep(:)
type(energyparameters)     :: para_pep(:)

!load peptide parameters and coordinates
if (CG_flag.eq.0) then
    call load_peptide_info(group_pep, para_pep)
else !coarse-grain all sidechains except for res_i and res_j
    cg_flags = 1
    cg_flags(res_i) = 0; cg_flags(res_j) = 0
    call load_peptide_info_cg(group_pep, para_pep, cg_flags)
endif

!load generalized born parameters, if including polar solvation energy in calculation
if (GB_flag.eq.1) call effgbradi_complex

call pairwise_res_energy(res_i, res_j, solvation_flag, vdw, ele, gb_iso, gb_sys)
binding_energy = gb_sys - gb_iso + weighting_factor*(vdw + ele + gb_iso) !first two terms contribute to total binding energy, while term in parentheses is peptide intramolecular energy in isolated state

return
end subroutine bindingenergy_pair

subroutine vdwcontri(rxy, epsion_x, epsion_y, r_x, r_y,vdw)
implicit none
double precision               :: vdw
real                           :: epsion_xy, r_xy, epsion_x, r_x, epsion_y, r_y, rxy
real                           :: acoeff, bcoeff

epsion_xy=sqrt(epsion_x*epsion_y)
r_xy=r_x+r_y
acoeff=epsion_xy*(r_xy**12)
bcoeff=epsion_xy*2*(r_xy**6)
vdw=acoeff/(rxy**6)-bcoeff/(rxy**3)

return
end subroutine vdwcontri

subroutine elecontri(rxy,qx, qy, dielecons_x, dielecons_y ,ele)
implicit none
double precision               :: ele
real                           :: qx, qy, dielecons4solute, rxy
real                           :: dielecons_x, dielecons_y 

dielecons4solute = MAX(dielecons_x, dielecons_y)
ele=(qx*qy)/(dielecons4solute*sqrt(rxy))

return
end subroutine elecontri

subroutine   gbcontri(rxy,alphax,alphay,qx, qy,dielecons_x, dielecons_y,gb)
implicit none
double precision         :: gb
real                     :: rxy
real                     :: fgb, alphax, alphay
real                     :: gb1, gb2
real                     :: qx, qy, dielecons4solute
real                     :: dielecons_x, dielecons_y


fgb=sqrt(rxy+alphax*alphay*exp(-rxy/(4*alphax*alphay)))
dielecons4solute = MAX(dielecons_x, dielecons_y)
gb1=(1.0/dielecons4solute)-(1.0/dielecons4water)
gb2=qx*qy/fgb
gb=-gb1*gb2

return
end subroutine gbcontri

subroutine effgbradi_complex
implicit none
integer                  :: res_i, atom_i
real                     :: integrals_rec_sys(atom_num_rec), integrals_pep_iso(atom_per_pep_res, pep_res), integrals_pep_sys(atom_per_pep_res, pep_res)
real                     :: alpha=0.8, gamma=2.91
real                     :: psi, i_rborn, i_redborn


!calculate the overlaps for isolated peptide and complex
call areafract_complex(integrals_rec_sys, integrals_pep_iso, integrals_pep_sys)

!calculate alphas for peptide
do res_i=1, pep_res
    do atom_i=1, atoms_per_res(res_i)

        !get parameters for atom
        i_rborn = rborn_pep(atom_i, res_i)
        i_redborn = i_rborn-0.09

        !get alpha for complex
        psi=integrals_pep_sys(atom_i, res_i)*i_redborn
        alpha_pep_sys(atom_i, res_i) = 1.0/(1.0/i_redborn-tanh(alpha*psi + gamma*psi*psi*psi)/i_rborn)

        !get alpha for isolated peptide
        psi=integrals_pep_iso(atom_i, res_i)*i_redborn
        alpha_pep_iso(atom_i, res_i) = 1.0/(1.0/i_redborn-tanh(alpha*psi + gamma*psi*psi*psi)/i_rborn)
    enddo
enddo

!calculate alphas for receptor
do atom_i=1, atom_num_rec
    !get parameters for atom
    i_rborn = rborn_rec(atom_i)
    i_redborn = i_rborn-0.09

    !get alpha for complex
    psi=integrals_rec_sys(atom_i)*i_redborn
    alpha_rec_sys(atom_i) = 1.0/(1.0/i_redborn-tanh(alpha*psi + gamma*psi*psi*psi)/i_rborn)

    !note: don't need to get alphas for isolated receptor, that is already done at beginning of program
enddo

return
end subroutine effgbradi_complex

subroutine areafract_complex(integrals_rec_comp, integrals_pep_iso, integrals_pep_comp)
!this subroutine calculates net atomic shielding for all atoms based on interatomic distances and sizes of all atom pairs
!NOTE: perhaps could speed up this function by skipping sum evaluation if rij is large? need to check how value scales as rij increases (not worrying about for now)
implicit none
integer                  :: atom_i, atom_j, res_i, res_j
real                     :: integrals_rec_comp(atom_num_rec), integrals_pep_iso(atom_per_pep_res, pep_res), integrals_pep_comp(atom_per_pep_res, pep_res)
real                     :: rij, sum_i, sum_j
real                     :: i_rborn, i_redborn, j_rborn, j_redborn

!initialize integrals to 0 for peptide, and to isolated values for receptor (note that values for isolated receptor are never changed during PepBD)
integrals_rec_comp=integrals_rec_iso; integrals_pep_iso=0.0

!get shielding from all peptide-peptide atom pairs
do res_i=1, pep_res
    do atom_i=1, atoms_per_res(res_i)
        !get parameters for atom
        i_rborn = rborn_pep(atom_i, res_i)
        i_redborn = fs_pep(atom_i, res_i)*(i_rborn-0.09)

        !first get shielding from atoms in same residue
        res_j = res_i
        do atom_j = atom_i+1, atoms_per_res(res_j)
            !get parameters for atom
            j_redborn=fs_pep(atom_j, res_j)*(rborn_pep(atom_j, res_j)-0.09)
            j_rborn = rborn_pep(atom_j, res_j)
            
            !look up interatomic distance
            rij = sqrt(rij_matrix_pep(atom_j, atom_i, res_j, res_i))
            
            !update area fraction overlap
            call areafract_helper(rij, i_rborn, i_redborn, j_rborn, j_redborn, sum_i, sum_j)
            integrals_pep_iso(atom_i, res_i) = integrals_pep_iso(atom_i, res_i) + sum_i
            integrals_pep_iso(atom_j, res_j) = integrals_pep_iso(atom_j, res_j) + sum_j
        enddo

        !next get shielding from atoms in other peptide residues
        do res_j = res_i+1, pep_res
            do atom_j=1, atoms_per_res(res_j)
                !get parameters for atom
                j_redborn=fs_pep(atom_j, res_j)*(rborn_pep(atom_j, res_j)-0.09)
                j_rborn = rborn_pep(atom_j, res_j)
                
                !look up interatomic distance
                rij = sqrt(rij_matrix_pep(atom_j, atom_i, res_j, res_i))
                
                !update area fraction overlap
                call areafract_helper(rij, i_rborn, i_redborn, j_rborn, j_redborn, sum_i, sum_j)
                integrals_pep_iso(atom_i, res_i) = integrals_pep_iso(atom_i, res_i) + sum_i
                integrals_pep_iso(atom_j, res_j) = integrals_pep_iso(atom_j, res_j) + sum_j
            enddo
        enddo
    enddo
enddo

integrals_pep_comp = integrals_pep_iso

!update shielding from atom-receptor atom pairs
do res_i=1, pep_res
    if (GB_area_calc_flags(res_i).eq.1) then !if 1, then we need to recalculate the areas for/due to this residue
        integrals_comp_pep_from_rec(:, res_i) = 0 !reset integral to 0 for all receptor atoms due to interactions with this peptide residue 
        integrals_comp_rec_from_pep(:, res_i) = 0 !reset integral to 0 for all atoms in this peptide residue due to interactions with receptor 
        do atom_i=1, atoms_per_res(res_i)
            !get parameters for atom
            i_rborn = rborn_pep(atom_i, res_i)
            i_redborn = fs_pep(atom_i, res_i)*(i_rborn-0.09)

            !go through all atoms in receptor
            do atom_j=1, atom_num_rec
                !get parameters for atom
                j_redborn=fs_rec(atom_j)*(rborn_rec(atom_j)-0.09)
                j_rborn = rborn_rec(atom_j)
                
                !look up interatomic distance
                rij = sqrt(rij_matrix_sys(atom_j, atom_i, res_i))
                
                !calculate area fraction overlap
                call areafract_helper(rij, i_rborn, i_redborn, j_rborn, j_redborn, sum_i, sum_j)
                integrals_comp_pep_from_rec(atom_i, res_i) = integrals_comp_pep_from_rec(atom_i, res_i) + sum_i
                integrals_comp_rec_from_pep(atom_j, res_i) = integrals_comp_rec_from_pep(atom_j, res_i) + sum_j
            enddo
        enddo
        GB_area_calc_flags(res_i) = 0
    endif

    !accumulate values for peptide
    do atom_i=1, atoms_per_res(res_i)
        integrals_pep_comp(atom_i, res_i) = integrals_pep_comp(atom_i, res_i) + integrals_comp_pep_from_rec(atom_i, res_i)
    enddo
enddo

!accumulate values for receptor
do atom_i=1, atom_num_rec
    integrals_rec_comp(atom_i) = integrals_rec_comp(atom_i) + sum(integrals_comp_rec_from_pep(atom_i, :))
enddo

!rescale values
integrals_pep_comp=integrals_pep_comp/2.0
integrals_pep_iso=integrals_pep_iso/2.0
integrals_rec_comp=integrals_rec_comp/2.0

return

end subroutine areafract_complex

subroutine areafract_helper(rij, i_rborn, i_redborn, j_rborn, j_redborn, sum_i, sum_j)
implicit none
real                     :: rij, sum_i, sum_j
real                     :: i_rborn, i_redborn, j_rborn, j_redborn
real                     :: lij, uij

!get shielding of atom_i by atom_j
call evalualij(rij,j_redborn,i_rborn,lij)
call evaluauij(rij,j_redborn,i_rborn,uij)
sum_i=(1.0/lij)-(1.0/uij)+(1.0/(uij*uij)-1.0/(lij*lij))*rij/4.0+log(lij/uij)/(2.0*rij)+ &
    (1.0/(lij*lij)-1.0/(uij*uij))*j_redborn*j_redborn/(4*rij)

!get shielding of atom_j by atom_i
call evalualij(rij,i_redborn,j_rborn,lij)
call evaluauij(rij,i_redborn,j_rborn,uij)
sum_j=(1.0/lij)-(1.0/uij)+(1.0/(uij*uij)-1.0/(lij*lij))*rij/4.0+log(lij/uij)/(2.0*rij)+ &
    (1.0/(lij*lij)-1.0/(uij*uij))*i_redborn*i_redborn/(4*rij)

end subroutine areafract_helper

subroutine effgbradi_iso_rec
implicit none
integer                  :: i
real                     :: alpha=0.8, gamma=2.91
real                     :: psi, i_rborn, i_redborn

!get pairwise distances for all atoms in receptor
call calc_pairwise_distances_iso_rec

!calculate the overlaps for isolated receptor
call areafract_iso_rec

!calculate alphas for isolated receptor
do i=1, atom_num_rec
    !get parameters for atom
    i_rborn = rborn_rec(i)
    i_redborn = i_rborn-0.09

    !get alpha for complex
    psi=integrals_rec_iso(i)*i_redborn
    alpha_rec_iso(i) = 1.0/(1.0/i_redborn-tanh(alpha*psi + gamma*psi*psi*psi)/i_rborn)
enddo

return
end subroutine effgbradi_iso_rec

subroutine areafract_iso_rec
implicit none
integer                  :: i, j
real                     :: rij, sum
real                     :: i_rborn, i_redborn, j_rborn, j_redborn
real                     :: lij, uij

integrals_rec_iso = 0

do i=1,atom_num_rec

    !get parameters for atom
    i_rborn = rborn_rec(i)
    i_redborn = fs_rec(i)*(i_rborn-0.09)

    do j=i+1,atom_num_rec
        !get parameters for atom
        j_rborn = rborn_rec(j)
        j_redborn = fs_rec(j)*(j_rborn-0.09)

        !lookup interatomic distances
        rij = sqrt(rij_matrix_rec(i,j))

        !get shielding of atom_i by atom_j
        call evalualij(rij,j_redborn,i_rborn,lij)
        call evaluauij(rij,j_redborn,i_rborn,uij)
        sum=(1.0/lij)-(1.0/uij)+(1.0/(uij*uij)-1.0/(lij*lij))*rij/4.0+log(lij/uij)/(2.0*rij)+ &
            (1.0/(lij*lij)-1.0/(uij*uij))*j_redborn*j_redborn/(4*rij)
        integrals_rec_iso(i)=integrals_rec_iso(i) + sum

        !get shielding of atom_j by atom_i
        call evalualij(rij,i_redborn,j_rborn,lij)
        call evaluauij(rij,i_redborn,j_rborn,uij)
        sum=(1.0/lij)-(1.0/uij)+(1.0/(uij*uij)-1.0/(lij*lij))*rij/4.0+log(lij/uij)/(2.0*rij)+ &
            (1.0/(lij*lij)-1.0/(uij*uij))*i_redborn*i_redborn/(4*rij)
        integrals_rec_iso(j)=integrals_rec_iso(j) + sum
    enddo
enddo

integrals_rec_iso = integrals_rec_iso/2.0

end subroutine areafract_iso_rec

subroutine evalualij(rij,redborn,rborn_x,lij)
implicit none
real                     :: rij, redborn, rborn_x
real                     :: lij

if(rborn_x.le.(rij-redborn)) then
   lij=rij-redborn
elseif(rborn_x.le.(rij+redborn)) then
   lij=rborn_x-0.09
else
   lij=1.0
endif

return
end subroutine evalualij

subroutine evaluauij(rij,redborn,rborn_x,uij)
implicit none
real                           :: rij, redborn
real                           :: uij
real                           :: rborn_x

if(rborn_x.lt.(rij+redborn)) then
   uij=rij+redborn
else
   uij=1.0
endif

return
end subroutine evaluauij

end module