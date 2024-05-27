module surface_area

use datatypes
use sys_vars
use utilities
use math

contains 

!subroutine isNeighbor(atom_i, atom_k, ns, num_ns, flag)
!implicit none
!integer         :: atom_i, atom_k, ns(max_neighbors,atom_num), num_ns(atom_num), i, flag
!
!if (atom_i.eq.atom_k) then 
!    flag=0
!    return
!else
!    do i=1, num_ns(atom_i)
!        if (ns(i,atom_i).eq.atom_k) then
!            flag=1
!            return
!        endif
!    enddo
!    flag=0
!    return
!endif
!
!end subroutine isNeighbor

subroutine A_IJ_eval(ri, rj, dij, aij)
implicit none
real                :: ri, rj, dij, aij

aij = TWO_PI*ri*(ri - 0.5*dij - (ri**2 - rj**2)/(2*dij))
return

end subroutine A_IJ_eval

subroutine eval_sasa (natom, recepsta, SASA_complex, SASA_pep, SASA_rec) 
implicit none
integer                      :: atom_i, atom_j, atom_k, natom, nj, nk, LCPO_i, LCPO_j, recepsta
real                         :: A_ij_sum, A_jk_sum, A_ij_jk_sum, A_jk_j
real                         :: A_ij_sum_lone, A_jk_sum_lone, A_ij_jk_sum_lone, A_jk_j_lone
real                         :: ri(3), rj(3), rk(3), ri_VDW, rj_VDW, rk_VDW, rij2, rjk2, rij2_VDW, rjk2_VDW, A_ij, A_jk_temp
real                         :: SASA_complex, SASA_pep, SASA_rec !, delta_SASA
integer                      :: ns(max_neighbors), num_ns !list of atom neighbors, and list of number of neighbors for each atom
real                         :: ns_dist(max_neighbors)
integer                      :: i_flag, j_flag, k_flag !if flag 0, then atom belongs to peptide; if flag 1, then atom belongs to receptor

!!Below structure of this subroutine is essentially taken from the cpptraj method
SASA_complex = 0; SASA_pep=0; SASA_rec=0 !set all surface areas equal to 0 initially 

!!for all atoms in the system...
!do atom_i=1, natom
!    num_ns = 0 !set number of neighbors equal for atom i equal to 0
!    LCPO_i = LCPO_indices(atom_i) !get the index of atom i into the LCPO data structure
!    if (LCPO_i.ne.-1) then !skip over atom i if it is a hydrogen
!        ri = mdcrd(1:3,atom_i) !get coordinates of atom i
!        ri_VDW = LCPO(LCPO_i)%r !get vdw radius of atom i
!        do atom_j=1, natom !for all other atoms in the system...
!            LCPO_j = LCPO_indices(atom_j) !get LCPO index of atom j into the LCPO data structure
!            if (LCPO_j.ne.-1.and.atom_i.ne.atom_j) then !skip over atom j if it is a hydrogen
!                rj = mdcrd(:,atom_j) !get coordinates of atom j
!                call DIST_SQUARE(ri, rj, rij2) !get distance between atom i and atom j
!                rij2_VDW =(ri_VDW + LCPO(LCPO_j)%r)**2 !get sum of VDW radii squared
!                if (rij2 < rij2_VDW) then !if true, then atom i and atom j are overlapping!
!                    num_ns = num_ns + 1 !increase neighbor counter by 1
!                    
!                    !make sure we aren't exceeding maximum storage for neighbor list
!                    if (num_ns.gt.max_neighbors) then 
!                        open(5, file="error.txt", access="append")
!                        write(5,*) "More than maximum neighbors found for an atom! This exceeds storage allocated for neighbor list."
!                        write(5,*) "Please increase size of neighbor list and restart design."
!                        write(5,*) "Terminating program."
!                        close(5)
!                        stop
!                    endif
!                    ns(num_ns) = atom_j !add atom j to the neighbor list 
!                    ns_dist(num_ns) = sqrt(rij2) !add distance between atom i and j to distance list
!                endif
!            endif
!        enddo 
!
!        !now calculate the SASA for the current atom
!        A_ij_sum=0; A_jk_sum = 0; A_ij_jk_sum = 0 !accumulators for SA in peptide + receptor complex
!        A_ij_sum_lone=0; A_jk_sum_lone = 0; A_ij_jk_sum_lone = 0 !accumulators for lone peptide or lone polymer (add to correct value depending on atom i)
!        if (atom_i.lt.recepsta) then
!            i_flag = 0 !atom i is in the peptide, so "lone" corresponds to peptide
!        else
!            i_flag = 1 !atom i is in the receptor, so "lone" corresponds to receptor
!        endif
!
!        do nj=1, num_ns !go through all neighbors of atom i and...
!            atom_j = ns(nj) !get atom number for the neighbor 
!            if (atom_j.lt.recepsta) then !figure out if atom j belongs to peptide or receptor
!                j_flag = 0 
!            else
!                j_flag = 1
!            endif
!
!            rj_VDW = LCPO(LCPO_indices(atom_j))%r  !get vdw radius of atom j
!            rj = mdcrd(:,atom_j)                   !get coordinates of atom j
!            call A_IJ_EVAL(ri_VDW, rj_VDW, ns_dist(nj), A_ij) !calculate overlap between atom i and atom j
!            !iterate through all neighbors of atom_i and see if they are also overlap with atom j
!            A_jk_j = 0; A_jk_j_lone = 0 !accumulate sum of overlaps of atom j with atoms that overlap with both atom i and and itself (in complex and lone system)
!            do nk=1, num_ns
!                if (nk.ne.nj) then !skip over atom j (don't count overlaps of atom with itself)
!                    atom_k = ns(nk)  !get atom k from neighbor list
!                    if (atom_k.lt.recepsta) then !determine if atom k belongs to peptide or receptor
!                        k_flag = 0
!                    else
!                        k_flag = 1
!                    endif
!
!                    rk = mdcrd(:,atom_k) !get coordinates of atom k
!                    rk_VDW = LCPO(LCPO_indices(atom_k))%r  !get vdw radius of atom k
!                    call DIST_SQUARE(rj, rk, rjk2) !get square distance between atom j and atom k
!                    rjk2 = sqrt(rjk2) !get distance between atom j and k
!                    rjk2_VDW = rj_VDW + rk_VDW !get sum of vdw radii
!                    if (rjk2 < rjk2_VDW) then !if true, then atom j and atom k are overlapping!
!                        call  A_IJ_EVAL(rj_VDW, rk_VDW, rjk2, A_jk_temp) !calculate overlap between the spheres
!                        A_jk_j = A_jk_j + A_jk_temp !accumulate A_jk value
!                        if (i_flag.eq.k_flag .and. j_flag.eq.k_flag) A_jk_j_lone = A_jk_j_lone + A_jk_temp !if all atoms belong to either the peptide or the receptor, then include overlap for lone system
!                    endif
!                endif
!            enddo
!
!            !accumulate values to sums
!            A_ij_sum = A_ij_sum + A_ij
!            A_jk_sum = A_jk_sum + A_jk_j
!            A_ij_jk_sum = A_ij_jk_sum + A_ij*A_jk_j
!            if (i_flag.eq.j_flag) then
!                A_ij_sum_lone = A_ij_sum_lone + A_ij
!                A_jk_sum_lone = A_jk_sum_lone + A_jk_j_lone
!                A_ij_jk_sum_lone = A_ij_jk_sum_lone + A_ij*A_jk_j_lone
!            endif
!        enddo
!
!        !accumulate to SASA of full system
!        SASA_complex = SASA_complex + LCPO(LCPO_i)%p1*(FOUR_PI * ri_VDW**2) &  
!                                    + LCPO(LCPO_i)%p2*A_ij_sum &
!                                    + LCPO(LCPO_i)%p3*A_jk_sum &
!                                    + LCPO(LCPO_i)%p4*A_ij_jk_sum 
!        if (i_flag.eq.0) then 
!            !accumulate to SASA of peptide
!            SASA_pep = SASA_pep + LCPO(LCPO_i)%p1*(FOUR_PI * ri_VDW**2) &  
!                                + LCPO(LCPO_i)%p2*A_ij_sum_lone &
!                                + LCPO(LCPO_i)%p3*A_jk_sum_lone &
!                                + LCPO(LCPO_i)%p4*A_ij_jk_sum_lone 
!        else
!            !accumulate to SASA of receptor
!            SASA_rec = SASA_rec + LCPO(LCPO_i)%p1*(FOUR_PI * ri_VDW**2) &  
!                                + LCPO(LCPO_i)%p2*A_ij_sum_lone &
!                                + LCPO(LCPO_i)%p3*A_jk_sum_lone &
!                                + LCPO(LCPO_i)%p4*A_ij_jk_sum_lone 
!        endif
!    endif
!enddo

!delta_SASA = SASA_complex - SASA_pep - SASA_rec

!! BELOW IS FOR DEBUGGIN' PURPOSES !!
!write (*,*) "Total surface area:", SASA_complex
!write (*,*) "Peptide surface area:", SASA_pep
!write (*,*) "Plastic surface area:", SASA_rec
!write (*,*) "Change in solvent accessible surface area:", delta_SASA
!stop

end subroutine eval_sasa

subroutine eval_SASA_energy(natom, recepsta, SASA_energy)
implicit none
integer                 :: natom, recepsta
real                    :: SASA_energy, SASA_complex, SASA_pep, SASA_rec

call eval_SASA(natom, recepsta, SASA_complex, SASA_pep, SASA_rec)

SASA_energy = surftens_rec*(SASA_complex - SASA_rec) - surftens_solv*SASA_pep

end subroutine eval_SASA_energy

end module surface_area