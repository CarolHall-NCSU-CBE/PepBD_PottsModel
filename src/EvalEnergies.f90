Program ProteinDesign

use database
use sys_vars
use utilities
use input
use pdbfile
use database
use energy_calculation
use advanced_function
use optimization_techniques

implicit none
integer                                         :: res_i, res_j, rotanum_i, rotanum_j, rot_i, feedback!, num_dihedrals
integer                                         :: AA_i, AA_j
integer, allocatable                            :: dihedral_net_nums(:)
integer                                         :: best_rot_i
double precision                                :: binding_energy, best_energy_rot_i, vdw_energy
character*20                                    :: outfile_single = 'SingleEnergy.txt', outfile_pairs = 'PairwiseEnergy.txt'
integer  , allocatable                          :: best_rot_array(:, :)
integer  , allocatable                          :: dihedrals(:,:)
type(groupdetails), allocatable                 :: aa_groups_opt(:,:)
type(databackup), allocatable                   :: groupdata_backup(:) 
type(groupdetails), allocatable                 :: group_pep(:), group_rec(:), original_pep(:)
type(energyparameters), allocatable             :: para_pep(:), para_rec(:)

!set up list of AA
AA_names(1:num_AAs) = (/'GLY', &
                'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'TYR', 'VAL', &
                'ARG', 'LYS', &
                'SER', 'THR', 'ASN', 'GLN', 'HIE', &
                'ALA', 'CYS', &
                'ASP', 'GLU'/)

call inputfile !read in the input file
call calc_system_size
call rotamerlib !load rotamers for all amino acids

!do all memory allocations (there are a lot!)

!system parameters
allocate(res_coords_flags(pep_res))
allocate(res_params_flags(pep_res))
allocate(phis(pep_res))
allocate(psis(pep_res))
allocate(best_rot_array(num_AAs,pep_res))
allocate(energy_per_res(pep_res, 4))
allocate(dihedral_net_nums(pep_res))
allocate(aa_groups_opt(num_AAs, pep_res))

!Peptide parameters
allocate(res_starts(pep_res))
allocate(res_ends(pep_res))
allocate(sc_starts(pep_res))
allocate(sc_ends(pep_res))
allocate(atoms_per_res(pep_res))
allocate(groupdata_backup(pep_res))
allocate(para_pep(pep_res))
allocate(group_pep(pep_res))
allocate(original_pep(pep_res))
allocate(Tgroup_pep(pep_res))
allocate(residuesite(pep_res))
allocate(atom(atom_per_pep_res, pep_res))
allocate(numex_pep(atom_per_pep_res, pep_res))
allocate(numex4_pep(atom_per_pep_res, pep_res))
allocate(inb_pep(2,20,atom_per_pep_res, pep_res))
allocate(inb4_pep(2,60,atom_per_pep_res, pep_res))
allocate(atomid_pep(atom_per_pep_res, pep_res))
allocate(mdcrd_pep(3,atom_per_pep_res, pep_res))
allocate(charge_pep(atom_per_pep_res, pep_res)) 
allocate(epsion_pep(atom_per_pep_res, pep_res))
allocate(r_pep(atom_per_pep_res, pep_res))
allocate(rborn_pep(atom_per_pep_res, pep_res))
allocate(fs_pep(atom_per_pep_res, pep_res))
allocate(dielecons_pep(atom_per_pep_res, pep_res))
allocate(dihedrals(9,pep_res))
allocate(LCPO_indices(atom_num))
if (solvation_flag.eq.1) then
    allocate(alpha_pep_sys(atom_per_pep_res, pep_res))
    allocate(alpha_pep_iso(atom_per_pep_res, pep_res))
endif

!receptor parameters
allocate(para_rec(rec_res))
allocate(group_rec(rec_res))
allocate(mdcrd_rec(3,atom_num_rec))
allocate(charge_rec(atom_num_rec)) 
allocate(epsion_rec(atom_num_rec))
allocate(r_rec(atom_num_rec))
allocate(rborn_rec(atom_num_rec))
allocate(fs_rec(atom_num_rec))
allocate(dielecons_rec(atom_num_rec))
if (solvation_flag.eq.1) then
    allocate(integrals_rec_iso(atom_num_rec))
    allocate(alpha_rec_sys(atom_num_rec))
    allocate(alpha_rec_iso(atom_num_rec))
    allocate(integrals_comp_rec_from_pep(atom_num_rec, pep_res)) 
    allocate(integrals_comp_pep_from_rec(atom_per_pep_res, pep_res))
    allocate(GB_area_calc_flags(pep_res)) 
endif

!pairwise distances
allocate(rij_matrix_pep(atom_per_pep_res, atom_per_pep_res, pep_res, pep_res))
allocate(rij_matrix_sys(atom_num_rec, atom_per_pep_res, pep_res))
allocate(rij_matrix_rec(atom_num_rec, atom_num_rec))

call load_system(group_pep, group_rec, original_pep, groupdata_backup) !load pdb file for system

!transfer coordinates and parameters for receptor into array format for energy calculations
call energy_parameter(group_pep, para_pep, group_rec, para_rec)
call load_params(group_rec, para_rec, 2)
call load_coords(group_rec, 2)

!get polar solvation energy parameters and energies for isolated receptor 
if (solvation_flag.eq.1) call calc_rec_iso_gb_energy

!check if the peptide backbone does not have an overlap
call check_backbone_overlap(group_pep, para_pep, feedback)
if (feedback.eq.0) then
    write (*,"(A)") "Starting peptide has overlap in backbone! Please provide a different backbone"
    write (*,"(A)") "Terminating program"
    stop
endif

!**************************************************************************
!      MAIN LOOP FOR CALCULATING ALL ONE-BODY AND TWO-BODY ENERGIES
!**************************************************************************

!calculate single body energies
open(301, file=outfile_single)
do res_i = 1, pep_res !for all residue positions, do...
    write (*, "(A,I3)") "Calculating Single Body Energies for Residue ", res_i
    write(301, 11) "#########    Residue ", res_i, " data    #########"
    do AA_i = 1, num_AAs !for all allowed residue types, do...

        !get the name of the residue
        AA_name_i = trim(AA_names(AA_i))
        res_name = AA_name_i
        if (res_i .eq. 1) then
            AA_name_i = "N" // AA_name_i
        elseif (res_i.eq.pep_res) then
            AA_name_i = "C" // AA_name_i
        endif

        call findrotamer(res_i, group_pep, AA_name_i, rotanum_i, aa_group_i)  !find rotamers for the residue

        !go through all rotamers at residue i and find the one with the lowest interaction energy between itself and the receptor
        best_rot_i = -1
        best_energy_rot_i = e_keep_sc_cutoff
        do rot_i=1, max(1,rotanum_i)
            call residue_replace_inplace(res_i, group_pep, groupdata_backup, rot_i, aa_group_i) !insert rotamer i into residue i

            !load parameters and form atom links, if first time seeing this residue
            if (rot_i.eq.1) then 
                call energy_parameter_single(group_pep, para_pep, res_i)
            endif

            !optimize sidechain (if applicable) and calculate binding energy
            if (res_name.eq.'ALA' .or. res_name.eq.'GLY') then
                call bindingenergy_single(group_pep, para_pep, res_i, solvation_flag, 1, binding_energy)
            else
                call sidechain_optimization(group_pep, para_pep, res_i, 2, solvation_flag, 1, binding_energy) !optimize sidechain for current residue
            endif
            
            if (binding_energy .lt. best_energy_rot_i) then
                best_rot_i = rot_i
                best_energy_rot_i = binding_energy

                !store updated coordinates of rotamer after sidechain optimization
                aa_groups_opt(AA_i, res_i) = group_pep(res_i)
            endif
        enddo

        !store index of best rotamer
        best_rot_array(AA_i, res_i) = best_rot_i

        !write to output file the single body energy
        write(301, 12) "Res", res_i, " Type ", group_pep(res_i)%gtype, " ", best_energy_rot_i

        !if coarse graining sidechains, then tell system to reload peptide data for current residue to end of peptide
        res_coords_flags(res_i) = 1
        res_params_flags(res_i) = 1
        atom_links_flag = 1
    enddo
enddo
close(301)

!now calculate pair-wise energies between all of the residues
if (eval_twobody_flag.eq.1 .and. weighting_factor.gt.0) then
    open(302, file=outfile_pairs)
    do res_i=1,pep_res-1 !the first residue in the pair
        write (*, "(A,I3)")   "Calculating Two Body Energies for Residue ", res_i
        write(302, 11) "#########    Residue ", res_i, " data    #########"
        do res_j=res_i+1, pep_res !the second residue in the pair
            do AA_i = 1,num_AAs !amino acid at first residue
                !get name of the amino acid
                AA_name_i = trim(AA_names(AA_i))
                if (res_i .eq. 1) then
                    AA_name_i = "N" // AA_name_i
                elseif (res_i.eq.pep_res) then
                    AA_name_i = "C" // AA_name_i
                endif
                
                if (best_rot_array(AA_i, res_i) .ne. -1) then !if we found a good rotamer for this residue at position i, then will calculate pairwise energies
                    !insert best rotamer for this amino acid at this residue (could probably store best rotamer in above loop so we don't need to call "findrotamer")
                    call findrotamer(res_i, group_pep, AA_name_i, rotanum_i, aa_group_i)
                    aa_group_i(1) = aa_groups_opt(AA_i,res_i)
                    call residue_replace_inplace(res_i, group_pep, groupdata_backup,  one, aa_group_i) 

                    do AA_j = 1, num_AAs !amino acid at the second residue
                        AA_name_j = trim(AA_names(AA_j))
                        if (res_j .eq. 1) then
                            AA_name_j = "N" // AA_name_j
                        elseif (res_i.eq.pep_res) then
                            AA_name_j = "C" // AA_name_j
                        endif

                        if (best_rot_array(AA_j, res_j) .ne. -1) then !only calculate pairwise energy if we found a good rotamer for this residue at position j
                            !insert best rotamer for this amino acid at this residue (could probably store best rotamer in above loop so we don't need to call "findrotamer")
                            call findrotamer(res_j, group_pep, AA_name_j, rotanum_j, aa_group_j)  !find rotamers for the residue
                            aa_group_j(1) = aa_groups_opt(AA_j,res_j)
                            call residue_replace_inplace(res_j, group_pep, groupdata_backup, one, aa_group_j) !insert best rotamer for AA at residue i

                            call bindingenergy_pair(group_pep, para_pep, res_i, res_j, solvation_flag, 1, binding_energy, vdw_energy) !calculate binding energy between res_i and res_j
                            binding_energy = min(binding_energy, e_keep_sc_cutoff)
                        else
                            binding_energy = e_keep_sc_cutoff
                        endif

                        !write pairwise energy to output file
                        write(302,13) "Res", res_i, " Type ", AA_name_i, " Res ", res_j," Type ", AA_name_j, " ", binding_energy  
                    enddo

                else !did not find a good rotamer for this amino acid at this residue position, so don't bother evaluating the pairwise energies
                    do AA_j = 1, num_AAs !amino acid at the second residue
                        AA_name_j = trim(AA_names(AA_j)) !get amino acid name
                        if (res_j .eq. 1) then
                            AA_name_j = "N" // AA_name_j
                        elseif (res_i.eq.pep_res) then
                            AA_name_j = "C" // AA_name_j
                        endif

                        !write pairwise energy to output file
                        write(302,13) "Res", res_i, " Type ", AA_name_i, " Res ", res_j, " Type ", AA_name_j, " ", e_keep_sc_cutoff  
                    enddo
                endif
            enddo
            !if coarse graining sidechains, then tell system to reload peptide data for current residue to end of peptide
            res_coords_flags(res_j:pep_res) = 1
            res_params_flags(res_j:pep_res) = 1
            atom_links_flag = 1
        enddo
        !if coarse graining sidechains, then tell system to reload peptide data for current residue to end of peptide
        res_coords_flags(res_i:pep_res) = 1
        res_params_flags(res_i:pep_res) = 1
        atom_links_flag = 1
    enddo
    close(302)
endif

11 format(a21, i2, a18)
12 format(a7, i4, a6, a4, a1, f8.3)
13 format(a4, i4, a6, a5, a5, i4, a6, a5, a1, f8.3)

deallocate(best_rot_array)
deallocate(dihedral_net_nums)
deallocate(energy_per_res)
!deallocate(rot_indices)
deallocate(dihedrals)
deallocate(para_rec)
deallocate(group_rec)
deallocate(para_pep)
deallocate(group_pep)
deallocate(groupdata_backup)
deallocate(original_pep)
deallocate(phis)
deallocate(psis)
deallocate(residuesite)
deallocate(alpha_pep_sys)
deallocate(alpha_pep_iso)
deallocate(alpha_rec_sys)
deallocate(alpha_rec_iso)
deallocate(numex_pep)
deallocate(numex4_pep)
deallocate(inb_pep)
deallocate(inb4_pep)
deallocate(atomid_pep)
deallocate(mdcrd_pep)
deallocate(charge_pep) 
deallocate(epsion_pep)
deallocate(r_pep)
deallocate(rborn_pep)
deallocate(fs_pep)
deallocate(dielecons_pep)
deallocate(mdcrd_rec)
deallocate(charge_rec) 
deallocate(epsion_rec)
deallocate(r_rec)
deallocate(rborn_rec)
deallocate(fs_rec)
deallocate(dielecons_rec)
deallocate(res_coords_flags)
deallocate(res_params_flags)
deallocate(res_starts)
deallocate(res_ends)
deallocate(sc_starts)
deallocate(sc_ends)
if (solvation_flag.eq.1) then
    deallocate(integrals_rec_iso)
    deallocate(integrals_comp_rec_from_pep)
    deallocate(integrals_comp_pep_from_rec)
    deallocate(GB_area_calc_flags)
endif
deallocate(rij_matrix_pep)
deallocate(rij_matrix_sys)
deallocate(rij_matrix_rec)
deallocate(atoms_per_res)

end program