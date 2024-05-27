module optimization_techniques

use datatypes
use sys_vars
use math
use utilities
use pdbfile
use database
use energy_calculation
use advanced_function
use omp_lib

contains

subroutine Eval_1body_energies(group_pep, para_pep, backbone_backup, aa_groups, rot_per_res, GB_flag, CG_flag, net_score, E_i)
implicit none
double precision        :: net_score, best_res_score, best_rot_score, E_i(40,pep_res), energy
integer                 :: res_i, rot_per_res(:), rot_i, GB_flag, CG_flag, num_good_rotamers
type(groupdetails)      :: group_pep(pep_res)
type(groupdetails)      :: aa_groups(:,:)
type(energyparameters)  :: para_pep(pep_res)
type(databackup)        :: backbone_backup(pep_res)

E_i = 0
if (CG_flag.eq.1) then !tell program to reload coordinates/parameters for residue 3 onward, since they are about to be coarse-grained
    res_coords_flags(3:) = 1
    res_params_flags(3:) = 1
    GB_area_calc_flags(3:) = 1
endif

net_score = 0
do res_i = 1, pep_res !for all residue positions, do...
    if (verbose.eq.1) write(*,"(A,I0)") "Calculating one-body energies for residue ", res_i
    best_res_score = e_keep_sc_cutoff

    !get the name of the residue
    AA_name_i = group_pep(res_i)%gtype
    if (res_i.eq.1 .or. res_i.eq.pep_res) AA_name_i = AA_name_i(2:4)
    
    !go through all rotamers and find the one with the lowest interaction energy between itself and the receptor
    best_rot_score=e_keep_sc_cutoff
    num_good_rotamers = 0
    do rot_i=1, rot_per_res(res_i)
        !insert rotamer i into residue i
        call residue_replace_inplace(res_i, group_pep, backbone_backup, rot_i, aa_groups(:,res_i))

        !optimize sidechain (if applicable) and calculate binding energy
        if (AA_name_i.eq.'ALA' .or. AA_name_i.eq.'GLY' .or. insert_pep_rigor.eq.0) then
            call bindingenergy_single(group_pep, para_pep, res_i, GB_flag, CG_flag, energy)
        else
            call sidechain_optimization(group_pep, para_pep, res_i, 2, GB_flag, CG_flag, energy) !optimize sidechain for current residue
            aa_groups(rot_i,res_i)%coo2 = group_pep(res_i)%coo2 !udpate sidechain coordinates in aa_groups to optimized values
        endif

        if (energy.lt.e_keep_sc_cutoff) then !keep rotamer as option if energy is below a threshold value or if rotamers are not being removed
        
            !!update rotamer info
            num_good_rotamers = num_good_rotamers + 1
        
            !store updated coordinates of rotamer
            call update_rotamer_coords(aa_groups(num_good_rotamers,res_i), group_pep, res_i)
        
            !store single body energy
            E_i(num_good_rotamers, res_i) = energy
        endif
    enddo

    !check if at least one good rotamer found for residue
    if (num_good_rotamers.eq.0) then
        if (verbose.eq.1) write(*,"(A, I2, A)") "!!!!!! All rotamers at residues at ", res_i, " have a steric overlap with the receptor !!!!!!"
        net_score = net_score + 100000
    else
        !add best rotamer energy to the net score for the peptide
        net_score = net_score + MINVAL(E_i(1:num_good_rotamers, res_i))
    endif


    !Update rotamer count to number that can be inserted without steric overlap with receptor
    rot_per_res(res_i) = num_good_rotamers
    if (CG_flag.eq.1) then
        !tell system to reload coords/parameters 
        res_coords_flags(res_i) = 1
        res_params_flags(res_i) = 1
        GB_area_calc_flags(res_i) = 1
        atom_links_flag = 1
    endif
enddo

end subroutine Eval_1body_energies

subroutine Eval_2body_energies(group_pep, para_pep, backbone_backup, aa_groups, rot_per_res, GB_flag, CG_flag, E_ij)
implicit none
double precision        :: E_ij(40,40,pep_res,pep_res), energy, vdw_energy
integer                 :: res_i, res_j, rot_i, rot_j, rot_per_res(:), GB_flag, CG_flag, CG_flags(pep_res), overlap_flag
integer                 :: good_rot_check(pep_res)
type(groupdetails)      :: group_pep(pep_res), aa_groups(:,:)
type(energyparameters)  :: para_pep(pep_res)
type(databackup)        :: backbone_backup(pep_res)

E_ij = 0
if (CG_flag.eq.1) then
    res_coords_flags(2:) = 1
    res_params_flags(2:) = 1
    GB_area_calc_flags(2:) = 1
endif
do res_i = 1, pep_res !for all residues in peptide
    if (verbose.eq.1) write(*,"(A,I0)") "Calculating two-body energies for residue ", res_i
    res_name = trim(group_pep(res_i)%gtype)
    if (res_i.eq.1 .or. res_i.eq.pep_res) res_name = res_name(2:4)
    overlap_flag = 1
    do rot_i=1, max(1,rot_per_res(res_i))  !for all rotamers at this residue (use the max statement since ALA and GLY listed as having 0 rotamers)
        good_rot_check = 0
        if (rot_per_res(res_i).gt.1) call residue_replace_inplace(res_i, group_pep, backbone_backup, rot_i, aa_groups(:,res_i)) !insert rotamer i into residue i
        do res_j=res_i+1, pep_res !for all other residues in the peptide upstream of current residue 
            do rot_j=1, max(1,rot_per_res(res_j)) !For all possible rotamers at this position, do... NOTE: use max statement to ensure pairwise energy is still calculated for residues with one rotamer option

                !insert rotamer j into residue j
                if (rot_per_res(res_j).gt.1) call residue_replace_inplace(res_j, group_pep, backbone_backup, rot_j, aa_groups(:,res_j)) 

                !calculate interaction energy between the two residues
                call bindingenergy_pair(group_pep, para_pep, res_i, res_j, GB_flag, CG_flag, energy, vdw_energy)

                !check if there is a steric overlap; if so, try to optimize the sidechain to remove the overlap
                if (vdw_energy.gt.1000) then 
                    if (res_name .ne. 'ALA' .and. res_name .ne. 'GLY') then
                        call sidechain_optimization4binding(group_pep, para_pep, res_i, GB_flag, CG_flag, CG_flags, energy)
                        aa_groups(rot_i,res_i)%coo2 = group_pep(res_i)%coo2 !udpate sidechain coordinates in aa_groups to optimized values
                        call bindingenergy_pair(group_pep, para_pep, res_i, res_j, GB_flag, CG_flag, energy, vdw_energy)
                    endif
                endif

                !store energy in matrix
                E_ij(rot_i, rot_j, res_i, res_j) =  energy 
                E_ij(rot_j, rot_i, res_j, res_i)  = energy

                !if energy below threshold, then indicate we found a pair of rotamers between residues with an overlap
                if (energy.lt.e_keep_sc_cutoff) good_rot_check(res_j) = 1
            enddo
            !tell system to reload coordinates and parameters (needed for accurate indexing into arrays since res_j will be coarse-grained)
            if (CG_flag.eq.1) then
                atom_links_flag = 1
                res_coords_flags(res_j) = 1
                res_params_flags(res_j) = 1
                GB_area_calc_flags(res_j) = 1
                CG_flags(res_j) = 1 !note that this doesn't do anything if CG_flag is set to 0
            endif
        enddo
        !check if current rotamer could be inserted without any steric overlaps
        if (sum(good_rot_check(res_i+1:pep_res)).eq.pep_res-res_i) overlap_flag = 0
        
    enddo
    !check if ANY rotamer could be inserted without any steric overlaps
    if (overlap_flag.eq.1 .and. verbose.eq.1) write(*,"(A, I2, A, I2, A)") "!!!!!! All rotamer pairs between residues ", res_i, " and ", res_j, " have a steric overlap !!!!!!"

    !tell system to reload coordinates and parameters (needed for accurate indexing into arrays since res_i will be coarse-grained)
    if (CG_flag.eq.1) then
        res_coords_flags(res_i) = 1
        res_params_flags(res_i) = 1
        GB_area_calc_flags(res_i) = 1
        CG_flags(res_i) = 1 !note that this doesn't do anything if CG_flag is set to 0
        atom_links_flag = 1
    endif
enddo

end subroutine Eval_2body_energies

subroutine MC_sequence(energy_new, energy_old, vdw_new, vdw_old, ele_new, ele_old, gb_new, gb_old, &
                       np_new, np_old, pep_energy_new, pep_energy_old, feedback)
implicit none
integer                     :: feedback
real                        :: ran2
double precision            :: energy_new, energy_old, energy_change
double precision            :: vdw_old, ele_old, gb_old, np_old, pep_energy_old
double precision            :: vdw_new, ele_new, gb_new, np_new, pep_energy_new

feedback=0
energy_change=energy_new-energy_old
if(energy_change.le.0) then
    energy_old = energy_new
    vdw_old    = vdw_new
    ele_old    = ele_new
    gb_old    = gb_new
    np_old   = np_new
    pep_energy_old = pep_energy_new
    feedback   = 1
else
    call RANDOM_NUMBER(ran2)
    if(ran2.le.exp(-energy_change/ekt_seq)) then
        energy_old = energy_new
        vdw_old    = vdw_new
        ele_old    = ele_new
        gb_old    = gb_new
        np_old   = np_new
        pep_energy_old =  pep_energy_new
        feedback   = 1
    endif
endif

return
end subroutine MC_sequence

subroutine MC_backbone(energy_new, energy_old, vdw_new, vdw_old, ele_new, ele_old, gb_new, gb_old, &
                        np_new, np_old, pep_energy_new, pep_energy_old, feedback)
implicit none
integer                       :: feedback
real                          :: ran2
double precision              :: energy_new, energy_old, energy_change
double precision              :: vdw_old, ele_old, gb_old, np_old, pep_energy_old
double precision              :: vdw_new, ele_new, gb_new, np_new, pep_energy_new

feedback=0
energy_change=energy_new-energy_old
if(energy_change.le.0) then
    energy_old=energy_new
    vdw_old=vdw_new
    ele_old=ele_new
    gb_old=gb_new
    np_old=np_new
    pep_energy_old = pep_energy_new
    feedback=1
else
    call RANDOM_NUMBER(ran2)
    if(ran2.le.exp(-energy_change/ekt_backbone)) then
        energy_old=energy_new
        vdw_old=vdw_new
        ele_old=ele_new
        gb_old=gb_new
        np_old=np_new
        pep_energy_old = pep_energy_new
        feedback=1
    endif
endif

return
end subroutine MC_backbone

subroutine SCMF_technique_pairwise(rotanum_1, rotanum_2, u, uij, i1_max, j2_max)
implicit none
integer                     :: i, j, rotanum_1, rotanum_2
integer                     :: i1_max, j2_max
double precision            :: u(2,40), uij(40,40), pmatrix_max
double precision            :: pmatrix_old(2,40), pmatrix_new(2,40), effenergy(2,40), effenergy_total, error

i1_max=0
j2_max=0
do i=1, rotanum_1
    pmatrix_old(1,i)=1.0/REAL(rotanum_1)
enddo
do j=1, rotanum_2
    pmatrix_old(2,j)=1.0/REAL(rotanum_2)
enddo

10   continue
pmatrix_new=pmatrix_old
effenergy=0.0
do i=1, rotanum_1
    effenergy(1,i)=u(1,i)
    do j=1, rotanum_2
        effenergy(1,i)=effenergy(1,i)+pmatrix_new(2,j)*uij(i,j)
    enddo
enddo
effenergy_total=0.0000001 !making slightly greater than 0 to prevent divide by 0 errors
do i=1, rotanum_1
    effenergy_total=effenergy_total+exp(-effenergy(1,i)/ekt_scmf)
enddo
do i=1, rotanum_1
    pmatrix_new(1,i)=exp(-effenergy(1,i)/ekt_scmf)/effenergy_total
enddo

do j=1, rotanum_2
    effenergy(2,j)=u(2,j)
    do i=1, rotanum_1
        effenergy(2,j)=effenergy(2,j)+pmatrix_new(1,i)*uij(i,j)
    enddo
enddo
effenergy_total=0.0000001 !making slightly greater than 0 to prevent divide by 0 errors
do j=1, rotanum_2
    effenergy_total=effenergy_total+exp(-effenergy(2,j)/ekt_scmf)
enddo
do j=1, rotanum_2
    pmatrix_new(2,j)=exp(-effenergy(2,j)/ekt_scmf)/effenergy_total
enddo

error=0.0
do i=1, rotanum_1
    error=error+abs(pmatrix_new(1,i)-pmatrix_old(1,i))
enddo
do j=1, rotanum_2
    error=error+abs(pmatrix_new(2,j)-pmatrix_old(2,j))
enddo

if(error.gt.0.01) then
    pmatrix_old=lampda*pmatrix_new+(1-lampda)*pmatrix_old
    goto 10
endif

pmatrix_max=-1.0
do i=1, rotanum_1
    if(pmatrix_max.le.pmatrix_new(1,i)) then
        pmatrix_max=pmatrix_new(1,i)
        i1_max=i
    endif
enddo

pmatrix_max=-1.0
do j=1, rotanum_2
    if(pmatrix_max.le.pmatrix_new(2,j)) then
        pmatrix_max=pmatrix_new(2,j)
        j2_max=j
    endif
enddo

return
end subroutine SCMF_technique_pairwise

subroutine SCMF_optimization_full(rot_per_res, u, uij, best_rotamers)
implicit none
integer                             :: res_i, res_j, rot_i, rot_j, rot_per_res(pep_res), num_iters, best_rotamers(pep_res)
double precision                    :: u(40,pep_res), uij(40,40,pep_res,pep_res), effenergy_total, error, temp
double precision                    :: pmatrix_old(40,pep_res), pmatrix_new(40,pep_res), effenergy(40,pep_res)
real            , parameter         :: threshold=0.01
double precision                    :: shift

!initially set probability matrix to uniform probability for all residue rotamers
do res_i=1, pep_res
    if (rot_per_res(res_i).gt.1) then
        pmatrix_old(1:rot_per_res(res_i),res_i) =1.0/REAL(rot_per_res(res_i))
    else
        pmatrix_old(:, res_i) = 1 !if only one rotamer for residue, then set p matrix to 1 since this rotamer HAS to exist
    endif
enddo

error=1.0
num_iters = 0
do while (error.gt.threshold) !main SCMF loop
    pmatrix_new=pmatrix_old !set probability matrix to last iteration values
    effenergy=0.0 !set energy matrix to 0

    !get the single body energy for each residue
    do res_i=1, pep_res
        if (rot_per_res(res_i).gt.1) then
            effenergy(1:rot_per_res(res_i), res_i) = u(1:rot_per_res(res_i),res_i)
        endif
    enddo

    !get two-body energy for each residue pairing. Note that we don't skip residues with one rotamer here, as doing so would mess up calculation for other residue
    do res_i=1, pep_res-1
        do res_j=res_i+1, pep_res
            do rot_i=1, max(1,rot_per_res(res_i))
                do rot_j=1, max(1,rot_per_res(res_j))
                    effenergy(rot_i,res_i) = effenergy(rot_i,res_i) + pmatrix_new(rot_j,res_j)*uij(rot_j,rot_i, res_j, res_i)
                    effenergy(rot_j,res_j) = effenergy(rot_j,res_j) + pmatrix_new(rot_i,res_i)*uij(rot_i,rot_j, res_i, res_j) 
                enddo
            enddo
        enddo
    enddo

    !use energies to get probabilities for each rotamer
    do res_i=1, pep_res
        if (rot_per_res(res_i).gt.1) then !skip residues with only one rotamer option
            effenergy_total=0.0000001 !making slightly greater than 0 to prevent divide by 0 errors
            shift = MINVAL(effenergy(1:rot_per_res(res_i),res_i))
            do rot_i=1, rot_per_res(res_i)
                temp = effenergy(rot_i,res_i)-shift

                !have a check below to prevent overflow errors
                if (temp.lt.50) then
                    temp = exp(-(effenergy(rot_i,res_i)-shift)/ekt_scmf)
                    pmatrix_new(rot_i, res_i) = temp
                    effenergy_total=effenergy_total+temp
                else
                    pmatrix_new(rot_i, res_i) = 0 !if energy if 50kcal/mol greater, then essentially 0 probability of this state existing
                endif
            enddo
            pmatrix_new(1:rot_per_res(res_i), res_i) = pmatrix_new(1:rot_per_res(res_i), res_i)/effenergy_total
        endif
    enddo

    !check if SCMF has converged
    error=0.0
    do res_i=1, pep_res
        if (rot_per_res(res_i).gt.1) then !skip residues with only one rotamer option
            do rot_i=1, rot_per_res(res_i)
                error=error+abs(pmatrix_new(rot_i,res_i)-pmatrix_old(rot_i,res_i))
            enddo
        endif
    enddo

    num_iters = num_iters + 1
    if(error.gt.0.01) then
        pmatrix_old=lampda*pmatrix_new+(1-lampda)*pmatrix_old
    endif
enddo

!now that SCMF has converged, find most likely rotamer for each residue
do res_i=1, pep_res
    if (rot_per_res(res_i).gt.1) then
        best_rotamers(res_i) = MAXLOC(pmatrix_new(1:rot_per_res(res_i), res_i), 1)
    else
        best_rotamers(res_i) = -1
    endif
enddo

return
end subroutine SCMF_optimization_full

subroutine sidechain_optimization(group_pep, para_pep, res_i, stage, GB_flag, CG_flag, score)
!Use BFGS algorithm to optimize sidechain of res_i in peptide
implicit none
integer                                             :: num_sc_angles, sc_groups(6), monitor(6)
integer                                             :: i, j, res_i, account_num, flag, trial_count, dihedral_num_options
integer                                             :: stage, GB_flag, CG_flag
real                                                :: delta_chi
double precision                                    :: energy_min, energy, score, error
real                                                :: CA(3)
double precision                                    :: h2_denominator, h3_denominator
double precision, dimension(:), allocatable         :: energy_forward, energy_backward
double precision, dimension(:,:), allocatable       :: gradient_old, gradient_new
double precision  , dimension(:,:), allocatable     :: Hessian_old, Hessian_new
double precision  , dimension(:,:), allocatable     :: H2, H3, H31
double precision, dimension(:,:), allocatable       :: d, y
real, dimension(:,:), allocatable                   :: s, Tchi
type(groupdetails)                                  :: group_pep(pep_res), Tgroup_pep(pep_res)
type(index4sidechain)                               :: index(60)
type(conformer4sidechain)                           :: sc_coords_storage(6), sc_coords_trial(6), sc_coords_min(6)
type(dihedralparameters)                            :: dihedral
type(energyparameters)                              :: para_pep(:)

!get sidechain info
call sidechain_category(res_i, group_pep, sc_coords_storage, num_sc_angles, sc_groups, index, monitor)
call dihedralangle_reading(group_pep(res_i)%gtype, dihedral_num_options, dihedral)

!allocate memory
allocate(energy_forward(num_sc_angles)); allocate(energy_backward(num_sc_angles))
allocate(gradient_old(num_sc_angles,1)); allocate(gradient_new(num_sc_angles,1))
allocate(Hessian_old(num_sc_angles,num_sc_angles)); allocate(Hessian_new(num_sc_angles,num_sc_angles))
allocate(H2(num_sc_angles,num_sc_angles)); allocate(H3(num_sc_angles,num_sc_angles)); allocate(H31(num_sc_angles,1))
allocate(d(num_sc_angles,1)); allocate(y(num_sc_angles,1)); allocate(s(num_sc_angles,1))
allocate(Tchi(num_sc_angles,1)) 

Tgroup_pep=group_pep !store group data in temporary structure

!search for alpha carbon position in residue
do i=1, Tgroup_pep(res_i)%cnum1
    if(Tgroup_pep(res_i)%atype1(i)=="CA") then
        CA(:)=Tgroup_pep(res_i)%coo1(1:3,i)
    endif
enddo

30   continue

!determine magnitude of change to sidechain angle. big changes used for coarse optimization, small changes used for fine tuning
if(stage.eq.0 .or. stage.eq.2) then
    delta_chi=5
elseif(stage.eq.1) then
    delta_chi=1
endif

!get energy of starting sidechain
call bindingenergy_single(Tgroup_pep, para_pep, res_i, GB_flag, CG_flag, energy_min)
sc_coords_min = sc_coords_storage

!get numerical derivative of binding energy for each dihedral angle
call dihedral_gradient_sidechain(Tgroup_pep, para_pep, res_i, sc_coords_storage, &
num_sc_angles, sc_groups, monitor, index, delta_chi, GB_flag, CG_flag, CA, energy_forward, energy_backward, gradient_new)

!initialize Hessian matrix as identity matrix
do i=1, num_sc_angles
    do j=1, num_sc_angles
        if(i==j) then
            Hessian_new(i,j)=1
        else
            Hessian_new(i,j)=0
        endif
    enddo
enddo

account_num=1 !tracks how often angles exceed the max dihedral change limits 
s=0.0         !the net "s"tep change for each dihedral angle, accumulated over the entire sidechain optimization process

!the main loop of BFGS algorithm 
do while(.true.)
    !set old values to the values generated in previous iteration (used to assess convergence)
    Hessian_old=Hessian_new !store Hessian in previous iteration
    gradient_old=gradient_new !store gradient in previous iteration
    d=-matmul(Hessian_old, gradient_old) !get direction to perform line search

    trial_count=0 !number of steps attempted
    Tchi=d*delta_chi !current search step size

    !check if changes to the dihedral are within the specified limits defined by current stage
    flag=0
    do i=1, num_sc_angles
        if(Tchi(i,1).gt.delta_chi) then
            Tchi(i,1)=delta_chi
            flag=1
        elseif(Tchi(i,1).lt.-delta_chi) then
            Tchi(i,1)=-delta_chi
            flag=1
        endif
    enddo

    do while(.true.) !the inner loop of BFGS algorithm. Performs linear search to find optimal rotation angle given current search direction
        sc_coords_trial=sc_coords_storage !restore to original sidechain values  (each step in loop takes step size with reference to initial sidechain coordinates)

        if(flag==1) account_num=account_num+1 !increase counter if dihedral found outside limits

        s=s+Tchi !the new change to dihedrals

        !rotate the entire sidechain by the new dihedral change
        call sidechain_rotation(num_sc_angles, sc_groups, s(:,1), sc_coords_trial, CA, Tgroup_pep, monitor, index, res_i)

        !tell program to reload coords for current residue into single array format
        res_coords_flags(res_i) = 1; GB_area_calc_flags(res_i) = 1

        !calculate energy at new position
        call bindingenergy_single(Tgroup_pep, para_pep, res_i, GB_flag, CG_flag, energy)

        !check if change to dihedrals led to a better binding energy
        if(energy - energy_min < -0.01) then !change to dihedral was an improvement!
            energy_min=energy
            sc_coords_min=sc_coords_trial
            trial_count=trial_count+1
        else !score did not improve, so backtrack to previous coordinates
            s=s-Tchi 
            exit
        endif
    enddo

    !if we tried to exceed the bounds too many times, then exit (the sidechain is deviating siginificantly (~100 degrees) from starting rotamer, so likely not good)
    if(stage==0) then
        if(account_num.gt.20) exit
    elseif(stage==1) then
        if(account_num.gt.40) exit
    endif

    if(trial_count==0) exit  !if we never found a better score during last iteration, then exit because we converged

    !check if we have converged, and exit if we have (convergence happens when magnitude of d (dot product of gradient and Hessian) is very small)
    error=0.0   
    do i=1, num_sc_angles
        error=error+abs(d(i,1))
    enddo
    if(error.lt.0.01) exit

    !get gradient at new position
    call dihedral_gradient_sidechain(Tgroup_pep, para_pep, res_i, sc_coords_min, &
    num_sc_angles, sc_groups, monitor, index, delta_chi, GB_flag, CG_flag, CA, energy_forward, energy_backward, gradient_new)

    !get change in gradient from previous iteration
    y=gradient_new-gradient_old

    !update Hessian per the BFGS algorithm
    H2=matmul(s, transpose(s))
    h2_denominator=0.0
    do i=1, num_sc_angles
        h2_denominator=h2_denominator+s(i,1)*y(i,1)
    enddo

    H3=matmul(Hessian_old, matmul(y, matmul(transpose(y), Hessian_old)))
    H31=matmul(Hessian_old, y)
    h3_denominator=0.0
    do i=1, num_sc_angles
        h3_denominator=h3_denominator+y(i,1)*H31(i,1)
    enddo

    Hessian_new=Hessian_old+H2/h2_denominator-H3/h3_denominator
enddo

if(stage==0.and.energy_min.lt.100.0) then
    !if refining sidechains with high precision, then continue optimization with smaller sidechain change
    stage=1
    goto 30
endif

!store energy, if a good change to dihedral angles was found
if(stage.eq.0) then
    score=1000.0
elseif(stage.eq.1 .or. stage.eq.2) then
    score=energy_min
endif

!deallocate memory
deallocate(energy_forward); deallocate(energy_backward)
deallocate(gradient_old); deallocate(gradient_new)
deallocate(Hessian_old); deallocate(Hessian_new)
deallocate(H2); deallocate(H3); deallocate(H31)
deallocate(d); deallocate(y); deallocate(s)
deallocate(Tchi)

!store the final coordinates into group
res_coords_flags(res_i) = 1; GB_area_calc_flags(res_i) = 1
do i=1, group_pep(res_i)%cnum2
    group_pep(res_i)%coo2(1:3,i)=sc_coords_min(index(i)%class_No)%member(1:3,index(i)%member_No)
enddo

return
end subroutine sidechain_optimization

subroutine sidechain_optimization4binding(group_pep, para_pep, res_i, GB_flag, CG_flag, CG_flags, energy_min)
implicit none
integer                                         :: num_sc_angles, sc_groups(6), monitor(6), dihedral_num_options, GB_flag, CG_flag, CG_flags(pep_res)
integer                                         :: i, j, res_i, account_num, flag, trial_count
real                                            :: delta_chi, error, t
double precision                                :: energy, energy_min
real                                            :: CA(3)
double precision                                :: h2_denominator, h3_denominator
double precision, dimension(:), allocatable     :: energy_forward, energy_backward
double precision, dimension(:,:), allocatable   :: gradient_old, gradient_new, Hessian_old, Hessian_new
double precision, dimension(:,:), allocatable   :: H2, H3, H31
double precision, dimension(:,:), allocatable   :: d, y
real, dimension(:,:), allocatable               :: s, Tchi
type(groupdetails)                              :: group_pep(:), Tgroup_pep(pep_res)
type(energyparameters)                          :: para_pep(:)
type(index4sidechain)                           :: index(60)
type(conformer4sidechain)                       :: sc_coords_storage(6), sc_coords_trial(6), sc_coords_min(6)
type(dihedralparameters)                        :: dihedral

!get sidechain info
call sidechain_category(res_i, group_pep, sc_coords_storage, num_sc_angles, sc_groups, index, monitor) 
call dihedralangle_reading(group_pep(res_i)%gtype, dihedral_num_options, dihedral)

!allocate memory for doing BFGS optimization
allocate(energy_forward(num_sc_angles)); allocate(energy_backward(num_sc_angles))
allocate(gradient_old(num_sc_angles,1)); allocate(gradient_new(num_sc_angles,1))
allocate(Hessian_old(num_sc_angles,num_sc_angles)); allocate(Hessian_new(num_sc_angles,num_sc_angles))
allocate(H2(num_sc_angles,num_sc_angles)); allocate(H3(num_sc_angles,num_sc_angles)); allocate(H31(num_sc_angles,1))
allocate(d(num_sc_angles,1)); allocate(y(num_sc_angles,1)); allocate(s(num_sc_angles,1))
allocate(Tchi(num_sc_angles,1))

Tgroup_pep=group_pep !store coordinate group_pep info into a temporary group_pep variable

!find the alpha carbon coordinates 
do i=1, Tgroup_pep(res_i)%cnum1
    if(Tgroup_pep(res_i)%atype1(i)=="CA") then
        CA=Tgroup_pep(res_i)%coo1(:,i)
    endif
enddo

!set step size for dihedral angle changes
delta_chi=5

!get energy of starting sidechain
call bindingenergy_sidechainoptimization(Tgroup_pep, para_pep, res_i, dihedral, dihedral_num_options, GB_flag, CG_flag, CG_flags, energy_min)
sc_coords_min = sc_coords_storage

!calculate numerical gradient for all rotamers
call dihedral_gradient_whole(Tgroup_pep, para_pep, res_i, sc_coords_storage, num_sc_angles, sc_groups, &
     dihedral, monitor, index, delta_chi, GB_flag, CG_flag, CG_flags, CA, energy_forward, energy_backward, gradient_new)

!initialize Hessian as identity matrix
do i=1, num_sc_angles
    do j=1, num_sc_angles
        if(i==j) then
            Hessian_new(i,j)=1
        else
            Hessian_new(i,j)=0
        endif
    enddo
enddo

account_num=0 !tracks number of times the dihedral angle change goes outside allowed limits
t=delta_chi        !stores initial magnitude of dihedral angle change for current search direction (NOTE: original program has this outside the DO loop, but should it be inside to get reset to 0 each time a new search direction is tried?)
s=0.0         !holds the step to make for each dihedral angle

do while(.true.) !outer loop of BFGS algorithm! This iterates until either energy no longer improved, or too many dihedral angle changes want to go outside limits
    30 continue
    Hessian_old=Hessian_new !store Hessian from previous iteration
    gradient_old=gradient_new !store gradient from previous iteration
    d=-matmul(Hessian_old, gradient_old) !get the direction to take a new step
    trial_count=0
    Tchi=t*d !the step size for current search direction

    !check if any of the dihedral angle changes are outside the allowed limits
    flag=0
    do i=1, num_sc_angles 
        if(Tchi(i,1).gt.delta_chi) then
            Tchi(i,1)=delta_chi
            flag=1
        elseif(Tchi(i,1).lt.(-delta_chi)) then
            Tchi(i,1)=-delta_chi
            flag=1
        endif
    enddo

    do while(.true.) !inner loop of the BFGS algorithm! This iterates until moving in current search direction no longer improves binding energy
        sc_coords_trial=sc_coords_storage !restore to original sidechain values  (each step in loop takes step size with reference to initial sidechain coordinates)

        if(flag==1) account_num=account_num+1

        s=s+Tchi
        !rotate sidechain by current dihedral angles
        call sidechain_rotation(num_sc_angles, sc_groups, s(:,1), sc_coords_trial, CA, Tgroup_pep, monitor, index, res_i)

        res_coords_flags(res_i) = 1; GB_area_calc_flags(res_i) = 1

        !calculate binding energy with new dihedrals
        call bindingenergy_sidechainoptimization(Tgroup_pep, para_pep, res_i, dihedral, dihedral_num_options, GB_flag, CG_flag, CG_flags, energy)

        !compare binding energy to previous iterations
        if(energy - energy_min < -0.01) then !taking another step in this direction improved binding energy, so store value and keep going
            energy_min=energy
            sc_coords_min=sc_coords_trial
            trial_count=trial_count+1
        else  !score did not improve, so backtrack to previous coordinates and exit inner loop
            s=s-Tchi
            exit
        endif
    enddo

    if(account_num.gt.40) exit !if too many dihedral angle changes outside limit, then exit (the sidechain is deviating siginificantly from starting rotamer, so likely not good)

    if(trial_count==0) then
        if (delta_chi.eq.5) then
            t = 1; delta_chi = 1 !same as delta_chi = 1
            goto 30
        else
            exit !if no steps in current direction improved energy, then appear to be a minimum so exit
        endif
    endif

    error=0.0
    do i=1, num_sc_angles
        error=error+abs(d(i,1))
    enddo
    if(error.lt.0.01) exit !we have converged, so exit

    !recalculate numerical gradient at new dihedral values
    call dihedral_gradient_whole(Tgroup_pep, para_pep, res_i, sc_coords_storage, num_sc_angles, sc_groups, &
    dihedral, monitor, index, delta_chi, GB_flag, CG_flag, CG_flags, CA, energy_forward, energy_backward, gradient_new)

    !calculate change in gradient
    y=gradient_new-gradient_old

    !update the hessian per the BFGS algorithm
    H2=matmul(s, transpose(s))
    h2_denominator=0.0
    do i=1, num_sc_angles
        h2_denominator=h2_denominator+s(i,1)*y(i,1)
    enddo

    H3=matmul(Hessian_old, matmul(y, matmul(transpose(y), Hessian_old)))
    H31=matmul(Hessian_old, y)
    h3_denominator=0.0
    do i=1, num_sc_angles
        h3_denominator=h3_denominator+y(i,1)*H31(i,1)
    enddo

    if (h3_denominator.eq.0 .or. h2_denominator.eq.0) exit 

    Hessian_new=Hessian_old+H2/h2_denominator-H3/h3_denominator
enddo

!deallocate memory
deallocate(energy_forward); deallocate(energy_backward)
deallocate(gradient_old); deallocate(gradient_new)
deallocate(Hessian_old); deallocate(Hessian_new)
deallocate(H2); deallocate(H3); deallocate(H31)
deallocate(d); deallocate(y); deallocate(s)
deallocate(Tchi)

!put optimized sidechain coordinates back into "group_pep" variable
res_coords_flags(res_i) = 1; GB_area_calc_flags(res_i) = 1
do i=1, group_pep(res_i)%cnum2
    group_pep(res_i)%coo2(1:3,i)=sc_coords_min(index(i)%class_No)%member(1:3,index(i)%member_No)
enddo

return
end subroutine sidechain_optimization4binding

subroutine sidechain_repacking_vdwele_dihedral(group_pep, para_pep, backbone_backup, Tresista, Tresiend, GB_flag, CG_flag)
implicit none
integer                         :: attempt, i, rotanum_1, rotanum_2, feedback_1, GB_flag, CG_flag
integer                         :: res_i, res_j, ip, flag, flag1, stage
integer                         :: Tresista, Tresiend, Tresinum, rounds
real                            :: ran2
double precision                :: energy, energy_min, score, Tscore
type(groupdetails)              :: group_pep(pep_res)
type(databackup)                :: backbone_backup(pep_res)
type(groupdetails)              :: res_i_backup, res_j_backup !stores original coordinates of residue i and j before changing sidechain
type(groupdetails)              :: res_i_best, res_j_best !stores best rotamer coordinates at residue i and j 
type(energyparameters)          :: para_pep(pep_res)

Tresinum=Tresiend-Tresista+1
if((Tresinum*3).lt.(pep_res*2)) then
    rounds=Tresinum*3
else
    rounds=pep_res*2
endif
do attempt=1, rounds  !!pep_res is equal to the number amino acids in the peptide that you wanna study
    call RANDOM_NUMBER(ran2)
    res_i=int(ran2*Tresinum)+1  !! res_i is randomly picked first
    if(res_i.gt.Tresinum) res_i=Tresinum
    flag = 0
    do while(flag.eq.0)
        call RANDOM_NUMBER(ran2)
        res_j=int(ran2*Tresinum)+1 !!res_j is the second amino acid selected for exchange, randomly chosen
        if(res_j.gt.Tresinum) res_j=Tresinum
        if(res_i.ne.res_j) then !!this makes sure that you don't pick the same residue twice
        flag = 1
        endif
    enddo
    res_i=res_i+Tresista-1; res_j=res_j+Tresista-1

    AA_name_i=group_pep(res_i)%gtype
    AA_name_j=group_pep(res_j)%gtype

    call findrotamer(res_i, group_pep, AA_name_i, rotanum_1, aa_group_i) !load rotamers for residue i
    call findrotamer(res_j, group_pep, AA_name_j, rotanum_2, aa_group_j) !load rotamers for residue j

    res_i_backup = group_pep(res_i) !backup original coordinates for residue i
    res_j_backup = group_pep(res_j) !backup original coordinates for residue j

    flag=0
    if(AA_name_i=="GLY".or.AA_name_i=="ALA".or.AA_name_i=="PRO".or. & 
    AA_name_i=="NGLY".or.AA_name_i=="NALA".or.AA_name_i=="NPRO".or. & 
    AA_name_i=="CGLY".or.AA_name_i=="CALA".or.AA_name_i=="CPRO") then
        energy_min=200.0   !initialize
        ip=0 !initialize the rotamer index number
        do i=1, rotanum_1
            call residue_replace_inplace(res_i, group_pep, backbone_backup, i, aa_group_i) !replaces the residue in the temporary group
            call check_overlap(0, res_i, zero, group_pep, para_pep, feedback_1) !this checks for atomic overlaps and gives answer in feedback_1
            if(feedback_1==1) then !this decision will calculate the VdW energy if there are no atomic overlaps
                call bindingenergy_single( group_pep, para_pep, res_i, GB_flag, CG_flag, energy)
                if(energy.lt.energy_min) then !if the calculated energy is lower than the previous minimum, this stores the new minumum and also stores the index to the chosen rotamer
                    energy_min=energy
                    ip=i
                endif
            endif
        enddo

        if(ip>0) then ! we found a good rotamer! so put into peptide and move onto residue j
            call residue_replace_inplace(res_i, group_pep, backbone_backup, ip, aa_group_i) !put best rotamer back into group
            flag=1
        else
            group_pep(res_i) = res_i_backup !restore original coordinates as a good rotamer not found 
        endif

    else !residue is not alanine, glycine, or proline
        flag1=0
        Tscore=1000
        do i=1, rotanum_1
            call residue_replace_inplace(res_i, group_pep, backbone_backup, i, aa_group_i) !replaces the residue
            stage=0
            call sidechain_optimization(group_pep, para_pep, res_i, stage, GB_flag, CG_flag, score)
            if(stage==1) then
                if(score.le.Tscore) then
                    Tscore=score
                    res_i_best = group_pep(res_i)
                endif
                flag1=1
            endif
        enddo

        if(flag1==1) then
            group_pep(res_i) = res_i_best
            flag=1
        else
            group_pep(res_i) = res_i_backup !restore original coordinates as a good rotamer not found
        endif
    endif

    !if we found a good replacement for residue i, then move on to residue j
    if(flag==1) then
        if(AA_name_j=="GLY".or.AA_name_j=="ALA".or.AA_name_j=="PRO".or. & 
        AA_name_j=="NGLY".or.AA_name_j=="NALA".or.AA_name_j=="NPRO".or. & 
        AA_name_j=="CGLY".or.AA_name_j=="CALA".or.AA_name_j=="CPRO") then
            energy_min=200.0   !initialize
            ip=0 !initialize the rotamer index number
            do i=1, rotanum_2
                call residue_replace_inplace(res_j, group_pep, backbone_backup, i, aa_group_j) !replaces the residue
                call check_overlap(0, res_j, zero, group_pep, para_pep, feedback_1) !this checks for atomic overlaps and gives answer in feedback_1
                if(feedback_1==1) then !this decision will calculate the VdW energy if there are no atomic overlaps
                    call bindingenergy_single(group_pep, para_pep, res_i, GB_flag, CG_flag, energy)
                    if(energy.lt.energy_min) then !if the calculated energy is lower than the previous minimum, this stores the new minumum and also stores the index to the chosen rotamer
                        energy_min=energy
                        ip=i
                    endif
                endif
            enddo

            if(ip>0) then !!this decision checks for an accepted rotamer and effects the actual transplant onto the molecule. then calculates new binding energy
                call residue_replace_inplace(res_j, group_pep, backbone_backup, ip, aa_group_j) !put in the best rotamer into the group
            else
                group_pep(res_j) = res_j_backup !restore original coordinates as a good rotamer not found
            endif

        else !the residue is not alanine, glycine, or proline
            flag1=0
            Tscore=2000
            do i=1, rotanum_2
                call residue_replace_inplace(res_j, group_pep, backbone_backup, i, aa_group_j) !replaces the residue
                stage=0
                call sidechain_optimization(group_pep, para_pep, res_j, stage, GB_flag, CG_flag, score)
                if(stage==1) then
                    if(score.le.Tscore) then
                        Tscore=score
                        res_j_best = group_pep(res_j)
                    endif
                    flag1=1
                endif
            enddo

            if(flag1==1) then
                group_pep(res_j) = res_j_best !put best rotamer into group_pep, if a good rotamer found
            else
                group_pep(res_j) = res_j_backup !restore original coordinates as a good rotamer not found 
            endif
        endif
    endif
enddo
return
end subroutine sidechain_repacking_vdwele_dihedral

subroutine sidechain_repacking_old(group_pep, para_pep, Tresista, Tresiend, GB_flag, binding_energy_min, vdw_min, ele_min, gb_min, np_min, pep_energy_min)
implicit none
integer                       :: attempt, feedback_1, MC_feedback, GB_flag, CG_flag, CG_flags(pep_res)
integer                       :: res_i, res_j
integer                       :: Tresista, Tresiend, Tresinum, rounds
character*4                   :: AA_name
real                          :: ran2
double precision              :: binding_energy_old, vdw_old, ele_old, gb_old, np_old, pep_energy_old
double precision              :: binding_energy_new, vdw_new, ele_new, gb_new, np_new, pep_energy_new
double precision              :: binding_energy_min, vdw_min, ele_min, gb_min, np_min, pep_energy_min
double precision              :: energy_min
type(groupdetails)            :: group_pep(pep_res)
type(groupdetails)            :: res_i_backup, res_j_backup
type(energyparameters)        :: para_pep(pep_res)

call bindingenergy_full(group_pep, para_pep, binding_energy_old, vdw_old, ele_old, gb_old, np_old, pep_energy_old)
CG_flag = 0
CG_flags = 0
binding_energy_min=binding_energy_old
vdw_min=vdw_old
ele_min=ele_old
gb_min=gb_old
np_min = np_old
pep_energy_min =  pep_energy_old

!select number of rounds for rotamer optimization
Tresinum=Tresiend-Tresista+1
if((Tresinum*5).lt.(pep_res*4)) then
    rounds=Tresinum*5
else
    rounds=pep_res*4
endif

do attempt=1, rounds 
    !choose two residues for optimizing sidechains
    call RANDOM_NUMBER(ran2)
    res_i = CEILING(ran2*Tresinum)

    res_j=res_i
    do while(res_i.eq.res_j)
        call RANDOM_NUMBER(ran2)
        res_j = CEILING(ran2*Tresinum)
    enddo
    res_i=res_i+Tresista-1; res_j=res_j+Tresista-1

    binding_energy_new=binding_energy_old+100.0

    feedback_1=0 !tracks if sidechain for res_i doesn't have overlaps
    AA_name = group_pep(res_i)%gtype
    if (res_i.eq.1 .or. res_i.eq.pep_res) AA_name(1:3) = AA_name(2:4)

    res_i_backup = group_pep(res_i)
    if(AA_name=="GLY" .or. AA_name=="ALA" .or. AA_name=="PRO") then !no sidechain to optimize for this residue
        feedback_1=1
    else
        call sidechain_optimization4binding(group_pep, para_pep, res_i, GB_flag, CG_flag, CG_flags, energy_min) !optimize side chain for this residue
        !call check_overlap(0, res_i, zero, group_pep, para_pep, feedback_1) !check if there is steric overlap in system; NOTE: perhaps could use energy from sidechain_optimization4binding to determine if there is an overlap
    endif

    if(energy_min.lt.100) then
        res_j_backup = group_pep(res_j)
        AA_name = group_pep(res_j)%gtype
        if (res_j.eq.1 .or. res_j.eq.pep_res) AA_name(1:3) = AA_name(2:4)
        if(AA_name=="GLY" .or. AA_name=="ALA" .or. AA_name=="PRO") then !no sidechain to optimize for this residue
            feedback_1=1
        else
            call sidechain_optimization4binding(group_pep, para_pep, res_j, 0, CG_flag, CG_flags, energy_min)
            !call check_overlap(0, res_j, zero, group_pep, para_pep, feedback_1)
        endif

        if (energy_min.lt.100) then
            !NOTE: we probably don't need to use Metropolis criterion here to accept/reject change, only essential to do that after backbone fully optimized. Leaving this untouched for now, though
            call bindingenergy_full(group_pep, para_pep, binding_energy_new, vdw_new, ele_new, gb_new, np_new,  pep_energy_new)
            call MC_sequence(binding_energy_new, binding_energy_old, vdw_new, vdw_old, &
                        ele_new, ele_old, gb_new, gb_old, np_new, np_old, pep_energy_new, pep_energy_old, MC_feedback)
            if(MC_feedback==1) then
                if(binding_energy_old.lt.binding_energy_min) then
                    binding_energy_min=binding_energy_old
                    vdw_min=vdw_old
                    ele_min=ele_old
                    gb_min=gb_old
                    pep_energy_min= pep_energy_old
                    np_min=np_old
                endif
            else !failed optimization of Metropolis check, so go back to original rotamers
                group_pep(res_i) = res_i_backup
                group_pep(res_j) = res_j_backup 
            endif
        else !failed at insertion of residue j rotamer, so go back to original rotamers
            group_pep(res_i) = res_i_backup
            group_pep(res_j) = res_j_backup 
        endif
    else !failed at insertion of rotamer of residue i, so go back to original rotamers
        group_pep(res_i) = res_i_backup 
    endif
enddo

return
end subroutine sidechain_repacking_old

subroutine sidechain_repacking(group_pep, para_pep, backbone_backup, GB_flag, CG_flag, CG_flags, binding_energy_min, vdw_min, ele_min, gb_min, np_min, pep_energy_min)
implicit none
integer                         :: res_i, rot_per_res(pep_res), best_rotamers(pep_res), iter, GB_flag, CG_flag, CG_flags(pep_res)
double precision                :: binding_energy_min, vdw_min, ele_min, gb_min, np_min, pep_energy_min
double precision                :: net_1body_score, E_i(40,pep_res), E_ij(40,40,pep_res,pep_res), E_pre
type(groupdetails)              :: group_pep(pep_res)
type(energyparameters)          :: para_pep(pep_res)
type(databackup)                :: backbone_backup(pep_res)

!load the rotamers
do res_i=1, pep_res
    AA_name_i = group_pep(res_i)%gtype
    call findrotamer(res_i, group_pep, AA_name_i, rot_per_res(res_i), aa_groups(:,res_i)) 
enddo

!get one- and two-body energies 
call Eval_1body_energies(group_pep, para_pep, backbone_backup, aa_groups, rot_per_res, GB_flag, CG_flag, net_1body_score, E_i)
call Eval_2body_energies(group_pep, para_pep, backbone_backup, aa_groups, rot_per_res, GB_flag, CG_flag, E_ij)

!use SCMF to find optimal sidechain set
call SCMF_optimization_full(rot_per_res, E_i, E_ij, best_rotamers) !find the optimal rotamers using SCMF

!insert the best rotamers into the peptide
do res_i=1, pep_res
    if (best_rotamers(res_i).ne.-1) then
        call residue_replace_inplace(res_i, group_pep, backbone_backup, best_rotamers(res_i), aa_groups(:, res_i))
    endif
enddo

!get full binding energy
call bindingenergy_full(group_pep, para_pep, binding_energy_min, vdw_min, ele_min, gb_min, np_min, pep_energy_min) 

!now gradient descent sidechain positions (not currently used)
!E_pre = binding_energy_min + 1
!iter=0
!do while(E_pre - binding_energy_min .gt. 0.5)
!    E_pre = binding_energy_min
!    !optimize all sidechains
!    do res_i=1, pep_res
!        res_name = trim(group_pep(res_i)%gtype)
!        if (res_i.eq.1 .or. res_i.eq.pep_res) res_name = res_name(2:4)
!        if (res_name .ne. 'ALA' .and. res_name .ne. 'GLY') then
!            call sidechain_optimization4binding(group_pep, para_pep, res_i, 0, CG_flag, CG_flags, binding_energy_min)
!        endif
!    enddo
!    !recalculate binding energy
!    call bindingenergy_full(group_pep, para_pep, binding_energy_min, vdw_min, ele_min, gb_min, np_min, pep_energy_min) 
!    iter = iter+1
!enddo
end subroutine sidechain_repacking

subroutine sequence_randomization(group_pep, para_pep, backbone_backup)
implicit none
integer                                          :: i, rotanum_1, feedback_1, trp_delta, success, attempts
integer                                          :: res_i, ip, res_type_old, res_type_new, CG_flag, GB_flag
double precision                                 :: energy, energy_min=10.0
type(groupdetails)                               :: group_pep(pep_res), res_backup
type(databackup)                                 :: backbone_backup(pep_res)
type(energyparameters)                           :: para_pep(pep_res), para_backup

CG_flag = 0; GB_flag = 0
do res_i=1, pep_res  !go through each residue in peptide
    !insert a random amino acid (as long as hydration constraints satisfied), keep trying until it succeeds
    success = 0; attempts=0
    res_backup = group_pep(res_i)
    para_backup = para_pep(res_i)
    do while (success.eq.0)   
        call mc_choose_aminoacid(res_i, group_pep, AA_name_i, res_type_old, res_type_new, trp_delta) !this picks the residue which will replace the test residue
        call findrotamer(res_i, group_pep, AA_name_i, rotanum_1, aa_group_i)  !find rotamers for the selected amino acid
        ip=0 !stores index of rotamer successfully inserted into peptide
        res_backup = group_pep(res_i)
        do i=1, rotanum_1
            call residue_replace_inplace(res_i, group_pep, backbone_backup, i, aa_group_i) !do a trial of replacing the residue
            if(i==1) then
                para_backup = para_pep(res_i)
                call energy_parameter_single(group_pep, para_pep, res_i) !loads the force field paramters for the amino acid (only needs to happen once)
            endif
            call check_overlap(0, res_i, -1, group_pep, para_pep, feedback_1) !this checks for atomic overlaps and gives answer in feedback_1
            if(feedback_1==1) then !this decision will calculate the VdW energy if there are no atomic overlaps
                call bindingenergy_single(group_pep, para_pep, res_i, GB_flag, CG_flag, energy)
                if(energy.lt.energy_min) then !if the calculated energy is lower than the previous minimum, this stores the new minumum and also stores the index to the chosen rotamer
                    ip=i
        
                    cur_hyd(res_type_old) = cur_hyd(res_type_old) - 1
                    above_hyd(res_type_old) = above_hyd(res_type_old) - 1
                    below_hyd(res_type_old) = below_hyd(res_type_old) + 1

                    cur_hyd(res_type_new) = cur_hyd(res_type_new) + 1
                    above_hyd(res_type_new) = above_hyd(res_type_new) + 1
                    below_hyd(res_type_new) = below_hyd(res_type_new) - 1

                    trp_count = trp_count + trp_delta
                    success = 1
                    exit
                endif
            endif
        enddo

        if (ip.eq.0) then
            attempts = attempts + 1
            if (attempts.gt.20) then
                write(*,"(A,I0,A,I0,A)") "Unable to radomize amino acid at position at residue ", res_i, " after ", attempts, " attempts!"
                write(*,"(A)") "Will not randomize this residue"
                exit
            endif
        endif
    enddo
enddo

return
end subroutine sequence_randomization

subroutine insert_random_dihedrals(group_pep, para_pep, temp_group, aa_groups, aa_groups_temp, cutoff, rot_per_res, CG_flag, CG_flags)
implicit none
integer                         :: res_i, res_j, feedback, flag, start, end, CG_flag, CG_flags(pep_res)
real                            :: phi, psi
double precision                :: binding_energy_min, vdw_min, ele_min, gb_min, np_min, pep_energy_min
integer                         :: cutoff
integer                         :: rot_per_res(pep_res)
type(groupdetails)              :: group_pep(pep_res), temp_group(pep_res), aa_groups(40, pep_res), aa_groups_temp(40,pep_res)
type(energyparameters)          :: para_pep(pep_res)
type(databackup)                :: dummy_backbone_backup(pep_res)

!debug variables
!real                            :: phis_post(pep_res), psis_post(pep_res)
!real                            :: phis_pre(pep_res), psis_pre(pep_res)

do res_i=1, pep_res
    flag = 0
    do while (flag.eq.0) !will keep trying dihedrals at current residue until we find one that works

        call pick_dihedral_pair(phi,psi) !pick a random pair of dihedral angles 
        
        !modify backbone angles per selected dihedral angles
        call insert_dihedral(temp_group, res_i, phi, psi, phis(res_i), psis(res_i), cutoff)

        !identify residues affected by rotation, so the rotamers for the residues can be updated
        if (res_i.le.cutoff) then 
            start = 1
            end = res_i
        else
            start = res_i
            end = pep_res
        endif

        !check for backbone overlaps
        call check_backbone_overlap(temp_group, para_pep, feedback)
        if (feedback.eq.0) then !backbone has overlap :'( undo change
            temp_group = group_pep
        else
            !adjust rotamers for amino acids that were affected by backbone rotation
            do res_j=start, end
                AA_name_j = group_pep(res_j)%gtype
                call findrotamer(res_j, temp_group, AA_name_j, rot_per_res(res_j), aa_groups_temp(:,res_j))
            enddo

            !repack sidechains
            call sidechain_repacking(group_pep, para_pep, dummy_backbone_backup, 0, CG_flag, CG_flags, binding_energy_min, vdw_min, ele_min, gb_min, np_min, pep_energy_min)

            !determine if sidechains could be repacked
            if (binding_energy_min .lt. 10) then !we succeeded!!
                group_pep = temp_group
                aa_groups = aa_groups_temp
                flag = 1
                phis(res_i) = phi
                psis(res_i) = psi
            
            else !we failed to repack the sidechains :'( Undo change 
                temp_group = group_pep
                aa_groups_temp = aa_groups
            endif
        endif
    enddo
enddo

end subroutine insert_random_dihedrals

subroutine sequence_optimization(group_pep, group_rec, para_pep, SA_best_pep, SA_best_para, backbone_backup, step, cycle, binding_energy_old, vdw_old, &
                        ele_old, gb_old, np_old, pep_energy_old)
implicit none
integer                    :: step, attempt, i, j, rotanum_1, rotanum_2, feedback_1, MC_feedback, cycle !, search_flag !basic variables and iterators
integer                    :: res_i, res_j, ip, flag, flag1, flag2, stage, res_type_old, res_type_new, trp_delta !basic variables
real                       :: ran2
double precision           :: binding_energy_old, vdw_old, ele_old, gb_old, np_old, pep_energy_old !stores energy values
double precision           :: binding_energy_new, vdw_new, ele_new, gb_new, np_new, pep_energy_new !stores energy values
double precision           :: energy, energy_min, score, Tscore !stores energy values
character*4                :: group_name_1(3), group_name_2(3) !names of amino acids
character*50               :: pdb_name
!real                       :: res_i_prob, p_shift 
integer                    :: pdb_flag
integer                    :: CG_flag, CG_flags(pep_res)
!integer                    :: TABU_swap, TABU_pick1, TABU_pick2 !variables for TABU search (not used at the moment)
type(groupdetails)         :: group_pep(pep_res), group_rec(rec_res), SA_best_pep(pep_res) !stores system coordinates 
type(databackup)           :: backbone_backup(pep_res)    !backup needed for proline residues 
type(groupdetails)         :: res_backup_i, res_backup_j !these variables backup original coordinates of residue i and residue j before attempting mutation
type(groupdetails)         :: res_i_best, res_j_best !these variables store best rotamers of residue i and residue j during search of all rotamers
type(energyparameters)     :: para_pep(pep_res), SA_best_para(pep_res), para_backup_i, para_backup_j !these variables backup original parameters of residue i and residue j before attempting mutation

!at present, disabling sidechain coarse-graining and polfor this step
CG_flag = 0

print ("(A,I0,A)"), "Step ", step, ": Sequence optimization"
cg_flags=0
do attempt=1, pep_res  !each sequence optimization step attempts as many modifications as there are amino acids in the peptide
    print ("(A,I0)"), "Sequence optimization attempt ", attempt
    pdb_flag = 0

    !tell system that atom links need to be reformed (true regardless of what type of mutation is made)
    atom_links_flag = 1

    !for current iteration, determine if we are doing a point mutation or swapping two amino acids
    call RANDOM_NUMBER(ran2)
    if(ran2.le.res_mutate_prob) then !we will try to do a point mutation

        !if (decomp_mutate_flag.eq.1) then !select residue to mutate based on contribution of residue to peptide score
        !    call RANDOM_NUMBER(ran2)
        !    search_flag = 0; res_i = 1
        !    do while(search_flag.eq.0)
        !        if (ran2.gt.res_select_p_vector(res_i)) then
        !            res_i = res_i+1   
        !        else
        !            search_flag = 1   
        !        endif
        !    enddo

        !elseif (TABU_flag.eq.1) then !select all residues with equal likelihood, but prevent picking residues that were changed recently
        !    call RANDOM_NUMBER(ran2)
        !    res_i=CEILING(ran2*TABU_allowed_length)
        !    res_i = TABU_allowed_residues(TABU_pick1) !get the actual residue number this corresponds to   
        !else !mutate any residue with equal likelihood
            !call RANDOM_NUMBER(ran2)
            !res_i=CEILING(ran2*pep_res)
        !endif

        call RANDOM_NUMBER(ran2)
        res_i=CEILING(ran2*pep_res)

        !update tally of how often we tried to mutate certain amino acids
        SEQ_POINT_attempts = SEQ_POINT_attempts + 1
        SEQ_PerRes_attempts(res_i) = SEQ_PerRes_attempts(res_i) + 1

        !pick a replacement residue and get its rotamers
        call mc_choose_aminoacid(res_i, group_pep, AA_name_i, res_type_old, res_type_new, trp_delta)
        call findrotamer(res_i, group_pep, AA_name_i, rotanum_1, aa_group_i)

        !initialize vdw energy to some arbitrary large number, and index to best rotamer for replacement amino acid to 0
        energy_min=200.0
        ip=0       

        !backup coordinates of residue at res_i
        res_backup_i = group_pep(res_i) 

        if(AA_name_i=="ALA".or.AA_name_i=="PRO".or.AA_name_i=="NALA".or.AA_name_i=="NPRO".or. &
        AA_name_i=="CALA".or.AA_name_i=="CPRO".or.AA_name_i=="GLY".or.AA_name_i=="NGLY".or. &
        AA_name_i=="CGLY") then

            !have special instructions if replacement amino acid is alanine, proline, or glycine, since these only have one rotamer
            do i=1, rotanum_1 !go through all rotamers of new amino acid and ...
                call residue_replace_inplace(res_i, group_pep, backbone_backup, i, aa_group_i) !insert the residue

                if(i==1) then
                    para_backup_i = para_pep(res_i) !store backup of original group's energy parameters
                    call energy_parameter_single(group_pep, para_pep, res_i) !load energy parameters for replacement amino acid
                endif

                call check_overlap(zero, res_i, zero, group_pep, para_pep, feedback_1) !check if inserted amino acid has steric clashes
                if(feedback_1==1) then 
                    !no steric clashes found, woohoo!
                    call bindingenergy_single(group_pep, para_pep, res_i, 0, 0, energy)
                    if(energy.lt.energy_min) then !store vdw energy if it is the lowest of all rotamers tested so far
                        energy_min=energy
                        ip=i !this is index to best rotamer
                    endif
                endif
            enddo

            if(ip>0) then !if we found a rotamer that worked, then...
                call residue_replace_inplace(res_i, group_pep, backbone_backup, ip, aa_group_i) !put that best rotamer back into the peptide 
                call bindingenergy_full(group_pep, para_pep, binding_energy_new, vdw_new, ele_new, gb_new, np_new, pep_energy_new) !calculate binding energy of the new peptide binding_energy_new, vdw_new, ele_new, gb_new, np_new, pep_energy_new)
                call MC_sequence(binding_energy_new, binding_energy_old, vdw_new, vdw_old, & !use Monte Carlo procedure to determine if we accept or reject the new sequence 
                                ele_new, ele_old, gb_new, gb_old, np_new, np_old,  pep_energy_new, pep_energy_old, MC_feedback)

                if (MC_feedback.eq.1) then
                    !we actually changed the peptide!!!

                    !update the hydration properties of the peptide
                    cur_hyd(res_type_old) = cur_hyd(res_type_old) - 1
                    above_hyd(res_type_old) = above_hyd(res_type_old) - 1
                    below_hyd(res_type_old) = below_hyd(res_type_old) + 1

                    cur_hyd(res_type_new) = cur_hyd(res_type_new) + 1
                    above_hyd(res_type_new) = above_hyd(res_type_new) + 1
                    below_hyd(res_type_new) = below_hyd(res_type_new) - 1

                    !update the tryptophan count
                    trp_count = trp_count + trp_delta

                    !update mutation attempt probabilities
                    SEQ_POINT_success = SEQ_POINT_success + 1
                    SEQ_PerRes_success(res_i) = SEQ_PerRes_success(res_i) + 1
                    
                    !update TABU search info, if using TABU search
                    !if (TABU_flag.eq.1) then
                    !    TABU_swap = TABU_allowed_residues(TABU_pick1)
                    !    TABU_allowed_residues(TABU_pick1) = TABU_allowed_residues(TABU_allowed_length)
                    !    TABU_allowed_residues(TABU_allowed_length) = TABU_swap
                    !    TABU_ban_lengths(TABU_pick1) = TABU_ban_lengths(TABU_allowed_length)
                    !    TABU_ban_lengths(TABU_allowed_length) = TABU_ban_length+1 !add 1 extra, since all values in array decreased by 1 before going back through loop
                    !    TABU_allowed_length = TABU_allowed_length-1               
                    !endif

                    !if this is best peptide found so far, then store the information
                    if (energy_min - binding_energy_old.lt.-.1) then
                        energy_min = binding_energy_old
                        if (SA_flag.eq.1) then
                            !store group and hydration properties of current peptide
                            SA_best_pep = group_pep
                            SA_best_para = para_pep
                            SA_best_hyd = cur_hyd
                            SA_above_hyd = above_hyd
                            SA_below_hyd = below_hyd
                            SA_trp_count = trp_count
                        endif
                        open(2, file="minimum_energy.txt", access="append")
                        write(2,40) "cycle", cycle, "step", step, "attempt", attempt, "EBind", binding_energy_old
                        close(2)
                        if (pdb_flag.eq.0) then
                            write(pdb_name, "(I0,A1,I0,A1,I0)")  cycle, '_', step, '_', attempt
                            call writepdb_sys(group_pep, group_rec, pdb_name)
                        endif
                    endif

                else !replacement failed at monte carlo step                    
                    !replace original amino acid parameters and coordinates
                    para_pep(res_i) = para_backup_i 
                    group_pep(res_i) = res_backup_i

                    !tell program that peptide parameters and coordinates need to be reloaded
                    res_coords_flags(res_i) = 1
                    res_params_flags(res_i) = 1
                    GB_area_calc_flags(res_i) = 1
                    atom_links_flag = 1
                endif
            else !replacement failed at check_overlap step
                !replace original amino acid parameters and coordinates
                para_pep(res_i) = para_backup_i 
                group_pep(res_i) = res_backup_i

                !tell program that peptide parameters and coordinates need to be reloaded
                res_coords_flags(res_i) = 1
                res_params_flags(res_i) = 1
                GB_area_calc_flags(res_i) = 1
                atom_links_flag = 1
            endif

        else !the replacement amino acid is neither alanine, glycine, nor proline, so use the general amino acid replacement method
            flag1=0
            do i=1, rotanum_1 !iterate through all rotamers of the replacement amino acid and see if it can fit 
                call residue_replace_inplace(res_i, group_pep, backbone_backup, i, aa_group_i) !insert rotamer of new amino acid into peptide
                if(i==1) then
                    para_backup_i = para_pep(res_i)                                                  !store backup of original group's energy parameters
                    call energy_parameter_single(group_pep, para_pep, res_i)                     !load the force field paramters for new amino acid
                endif

                stage=0 !stage determines how rigorous the sidechain_optimization is. Stage=0 is the faster method used at beginning, and optimization which switch to stage 1 if good energy found
                call sidechain_optimization(group_pep, para_pep, res_i, stage, 0, 0, score) !optimze sidechain, polar solvation energies turned off
                if(stage==1) then !sidechain optimization successful
                    !if best energy seen so far, then store the energy and corresponding group corrdinates
                    if(flag1==0) then !this is the first rotamer to work
                        Tscore=score
                        res_i_best=group_pep(res_i)
                        flag1=1
                    else !we have already found another rotamer that worked, so compare results between two rotamers
                        if(score.le.Tscore) then
                            Tscore=score
                            res_i_best=group_pep(res_i)
                        endif
                    endif
                endif
            enddo

            if(flag1==1) then !we found one rotamer of replacement amino acid that worked, so keep going!
                group_pep(res_i) = res_i_best
                call sidechain_optimization4binding(group_pep, para_pep, res_i, 0, CG_flag, CG_flags, energy_min) !optimize sidechain of best rotamer, still no GB energies
                call check_overlap(0, res_i, zero, group_pep, para_pep, feedback_1)                         !check sidechain has no steric overlaps 
                if(feedback_1==1) then 
                    !no steric overlaps found, so keep going!
                    call bindingenergy_full(group_pep, para_pep, binding_energy_new, vdw_new, ele_new, gb_new, np_new, pep_energy_new) !calculate energy of peptide 
                    call MC_sequence(binding_energy_new, binding_energy_old, vdw_new, vdw_old, & !use Monte Carlo to accept/reject change
                    ele_new, ele_old, gb_new, gb_old, np_new, np_old, pep_energy_new,  pep_energy_old, MC_feedback)

                    if (MC_feedback.eq.1) then !we actually changed the peptide!!!
                        
                        !update the hydration properties
                        cur_hyd(res_type_old) = cur_hyd(res_type_old) - 1
                        above_hyd(res_type_old) = above_hyd(res_type_old) - 1
                        below_hyd(res_type_old) = below_hyd(res_type_old) + 1

                        cur_hyd(res_type_new) = cur_hyd(res_type_new) + 1
                        above_hyd(res_type_new) = above_hyd(res_type_new) + 1
                        below_hyd(res_type_new) = below_hyd(res_type_new) - 1

                        !update tryptophan count
                        trp_count = trp_count + trp_delta

                        !update mutation attempt probabilities
                        SEQ_POINT_success = SEQ_POINT_success + 1
                        SEQ_PerRes_success(res_i) = SEQ_PerRes_success(res_i) + 1

                        !update TABU search info
                        !if (TABU_flag.eq.1) then
                        !    TABU_swap = TABU_allowed_residues(TABU_pick1)
                        !    TABU_allowed_residues(TABU_pick1) = TABU_allowed_residues(TABU_allowed_length)
                        !    TABU_allowed_residues(TABU_allowed_length) = TABU_swap
                        !    TABU_ban_lengths(TABU_pick1) = TABU_ban_lengths(TABU_allowed_length)
                        !    TABU_ban_lengths(TABU_allowed_length) = TABU_ban_length+1 !add 1 extra, since all values in array decreased by 1 before going back through loop
                        !    TABU_allowed_length = TABU_allowed_length-1               
                        !endif

                        !check if this the lowest energy we have found, and update simulated annealing best group (if applicable)
                        if (energy_min - binding_energy_old.lt.-.1) then
                            energy_min = binding_energy_old
                            if (SA_flag.eq.1) then
                                SA_best_pep = group_pep
                                SA_best_para = para_pep
                                SA_best_hyd = cur_hyd
                                SA_above_hyd = above_hyd
                                SA_below_hyd = below_hyd
                                SA_trp_count = trp_count
                            endif
                            open(2, file="minimum_energy.txt", access="append")
                                write(2,40) "cycle", cycle, "step", step, "attempt", attempt, "EBind", binding_energy_old
                            close(2)
                            if (pdb_flag.eq.0) then
                                write(pdb_name, "(I0,A1,I0,A1,I0)")  cycle, '_', step, '_', attempt
                                call writepdb_sys(group_pep, group_rec, pdb_name) 
                            endif
                        endif
                    else !replacement failed at monte carlo step
                        
                        !replace original amino acid parameters and coordinates
                        para_pep(res_i) = para_backup_i 
                        group_pep(res_i) = res_backup_i

                        !tell program that peptide parameters and coordinates need to be reloaded
                        res_coords_flags(res_i) = 1
                        res_params_flags(res_i) = 1
                        GB_area_calc_flags(res_i) = 1
                        atom_links_flag = 1
                    endif
                else !replacement failed at check_overlap step

                    !replace original amino acid parameters and coordinates
                    para_pep(res_i) = para_backup_i
                    group_pep(res_i) = res_backup_i

                    !tell program that peptide parameters and coordinates need to be reloaded
                    res_coords_flags(res_i) = 1
                    res_params_flags(res_i) = 1
                    GB_area_calc_flags(res_i) = 1
                    atom_links_flag = 1
                endif
            else !replacement failed at sidechain_optimization step 

                !replace original amino acid parameters and coordinates
                para_pep(res_i) = para_backup_i 
                group_pep(res_i) = res_backup_i

                !tell program that peptide parameters and coordinates need to be reloaded
                res_coords_flags(res_i) = 1
                res_params_flags(res_i) = 1
                GB_area_calc_flags(res_i) = 1
                atom_links_flag = 1
            endif
        endif

        open(3, file="energydetails.txt",  access="append")
        write(3,4) cycle, step, attempt, binding_energy_old, (binding_energy_old-pep_energy_old), "Mutation" !!Mutation indicates that it was a one amino change instead of swapping amino acids
        !write(3,"(<pep_res>(a4))") (group_pep(j)%gtype, j=1, pep_res)
        write(3,*) (group_pep(j)%gtype, j=1, pep_res)
        write(3,*) "*******************************"
        close(3)

    else !we are going to swap two amino acids already in the peptide
        !if (decomp_mutate_flag.eq.1) then
        !    !select the first residue per probability vector
        !    call RANDOM_NUMBER(ran2)
        !    search_flag = 0; res_i = 1
        !    do while(search_flag.eq.0)
        !        if (ran2.gt.res_select_p_vector(res_i)) then
        !            res_i = res_i+1   
        !        else
        !            search_flag = 1   
        !        endif
        !    enddo
        !
        !    !select second residue per probability vector (ensuring second residue is different from first residue)
        !    !do this in a way such that residue 1 is not reselected and are guaranteed to select a new residue 
        !    call RANDOM_NUMBER(ran2)
        !    search_flag = 0; res_j = 1; 
        !    if (res_i.eq.1) then
        !        res_i_prob = res_select_p_vector(res_i)
        !    else
        !        res_i_prob = res_select_p_vector(res_i) - res_select_p_vector(res_i-1)
        !    endif
        !    ran2 = ran2*(1 - res_i_prob) !this shifts cap on probability so we always are able to select a second residue
        !    p_shift = 0.0
        !    do while(search_flag.eq.0)
        !        if (ran2.gt.(res_select_p_vector(res_j)-p_shift).and.res_j.ne.res_i) then
        !            res_j = res_j + 1
        !            if (res_j.eq.res_i) then
        !                res_j = res_j + 1   
        !                p_shift = res_i_prob
        !            endif
        !        else
        !            search_flag = 1   
        !        endif
        !    enddo
        !
        !!elseif (TABU_flag.eq.1) then         !select res_i and res_j based on TABU search 
        !    call RANDOM_NUMBER(ran2)
        !    TABU_pick1 = CEILING(ran2*TABU_allowed_length) !index into TABU array to get res_i
        !    do while(.true.)
        !        call RANDOM_NUMBER(ran2)
        !        TABU_pick2 = CEILING(ran2*TABU_allowed_length)
        !        if(TABU_pick1.ne.TABU_pick2) then !check if res_i and res_j are actually different residues
        !            exit
        !        endif
        !    enddo
        !    res_i = TABU_allowed_residues(TABU_pick1)
        !    res_j = TABU_allowed_residues(TABU_pick2)

        !else  !randomly select the two amino acids
        !    call RANDOM_NUMBER(ran2)
        !    res_i = CEILING(ran2*pep_res) !res_i selection
        !    do while(.true.)
        !        call RANDOM_NUMBER(ran2)
        !        res_j = CEILING(ran2*pep_res) !res_j selection
        !        if(res_i.ne.res_j) then !check if res_i and res_j are actually different residues
        !            exit
        !        endif
        !    enddo
        !endif

        call RANDOM_NUMBER(ran2)
        res_i = CEILING(ran2*pep_res) !res_i selection
        do while(.true.)
            call RANDOM_NUMBER(ran2)
            res_j = CEILING(ran2*pep_res) !res_j selection
            if(res_i.ne.res_j) exit
        enddo

        !update move attempts for res_i and res_j
        SEQ_SWAP_attempts = SEQ_SWAP_attempts + 1
        SEQ_PerRes_attempts(res_i) = SEQ_PerRes_attempts(res_i) + 1
        SEQ_PerRes_attempts(res_j) = SEQ_PerRes_attempts(res_j) + 1

        !store backups for coordinates and parameters at two selected residues
        res_backup_i = group_pep(res_i)
        res_backup_j = group_pep(res_j)
        para_backup_i = para_pep(res_i) !store energy parameters for residue 1
        para_backup_j = para_pep(res_j) !store energy parameters for residue 

        !get different names for selected current amino acid
        call groupinfo(group_pep(res_i)%gtype, group_name_1, flag1)
        call groupinfo(group_pep(res_j)%gtype, group_name_2, flag2)

        AA_name_i=group_name_2(flag1) !load the name of the second residue into the slot for the first residue, and vice versa
        AA_name_j=group_name_1(flag2)
        call findrotamer(res_i, group_pep, AA_name_i, rotanum_1, aa_group_i) !get rotamers for residue 1
        call findrotamer(res_j, group_pep, AA_name_j, rotanum_2, aa_group_j) !get rotamers for residue 2

        flag=0
        energy_min=200.0 !initialize minimum energy to large number
        ip=0 !initialize best rotamer index

        !step 1: try to insert a rotamer for res_j into position of res_i
        if(AA_name_i=="GLY".or.AA_name_i=="ALA".or.AA_name_i=="PRO".or. &
        AA_name_i=="NGLY".or.AA_name_i=="NALA".or.AA_name_i=="NPRO".or.&
        AA_name_i=="CGLY".or.AA_name_i=="CALA".or.AA_name_i=="CPRO") then !Special instructions for glycine, alanine, and proline
            do i=1, rotanum_1
                call residue_replace_inplace(res_i, group_pep, backbone_backup, i, aa_group_i) !insert rotamer into residue
                if (i.eq.1) then
                    call energy_parameter_single(group_pep, para_pep, res_i) !load the force field paramters for new amino acid
                endif
                call check_overlap(0, res_i, zero, group_pep, para_pep, feedback_1) !this checks for atomic overlaps and gives answer in feedback_1
                if(feedback_1.eq.1) then !no overlaps found! now calculate vdw energy
                    call bindingenergy_single(group_pep, para_pep, res_i, 0, 0, energy)
                    if(energy.lt.energy_min) then !if calculated vdw energy less than previous minimum, stores the new minumum and index to the chosen rotamer
                        energy_min=energy
                        ip=i
                    endif
                endif
            enddo

            if(ip.gt.0) then ! found a good rotamer, so insert into peptide
                call residue_replace_inplace(res_i, group_pep, backbone_backup, ip, aa_group_i) !put in the rotamer
                res_i_best = group_pep(res_i) !store coordinates after inserting new rotamer 
                flag=1 !let program now that we found a good rotamer
            endif

            else !res_i is not alanine, glycine, or proline, so do the general swap method
            flag1=0
            do i=1, rotanum_1
                call residue_replace_inplace(res_i, group_pep, backbone_backup, i, aa_group_i) !insert rotamer i into the peptide
                if(i.eq.1) then
                    call energy_parameter_single(group_pep, para_pep, res_i)                       !load the force field paramters for new amino acid
                endif

                stage=0
                call sidechain_optimization(group_pep, para_pep, res_i, stage, 0, 0, score)    !optimize sidechain rotamer, no polar solvation energies
                if(stage==1) then !rotamer successfully inserted into peptide!
                    if(flag1==0) then !this is first good rotamer we've found
                        Tscore=score
                        res_i_best=group_pep(res_i)
                        flag1=1
                    else        !already found a good rotamer, so compare to previous best rotamer
                        if(score.le.Tscore) then
                        Tscore=score
                        res_i_best=group_pep(res_i)
                        endif
                    endif
                endif
            enddo

            if(flag1.eq.1) then !we found a good rotamer for res_i, so optimize that rotamer a bit more and make sure there aren't overlaps 
                group_pep(res_i) = res_i_best
                call sidechain_optimization4binding(group_pep, para_pep, res_i, 0, CG_flag, CG_flags, energy_min)
                call check_overlap(0, res_i, zero, group_pep, para_pep, feedback_1) !NOTE: I don't think this step necessary - can determine overlap based on energy from sidechain_optimization4binding
                if(feedback_1==1) then !no overlaps, so we have succesfully put a rotamer for res_j in new position!
                    res_i_best=group_pep(res_i) !store backup of new peptide 
                    flag=1 !let program know that we were able to insert a rotamer
                endif
            endif
        endif

        if(flag.eq.1) then !found a good rotamer of residue 2 at position 1, yay! Now look for good rotamer for residue 1 at position 2
            energy_min=200.0 !initialize VDW energy to high number 
            ip=0 !initialize the rotamer index number
            if(AA_name_j=="GLY".or.AA_name_j=="ALA".or.AA_name_j=="PRO" &
            .or.AA_name_j=="NGLY".or.AA_name_j=="NALA".or.AA_name_j=="NPRO".or. & 
            AA_name_j=="CGLY".or.AA_name_j=="CALA".or.AA_name_j=="CPRO") then
                do i=1, rotanum_2
                    call residue_replace_inplace(res_j, group_pep, backbone_backup, i, aa_group_j) !replaces the residue
                    if (i.eq.1) then
                        call energy_parameter_single(group_pep, para_pep,res_j) !load the force field paramters for new amino acid
                    endif
                    call check_overlap(0, res_j, zero, group_pep, para_pep, feedback_1) !checks for atomic overlaps
                    if(feedback_1==1) then !no overlaps found
                        call bindingenergy_single(group_pep, para_pep, res_i, 0, 0, energy)
                        if(energy.lt.energy_min) then !if the calculated energy is lower than the previous minimum vdw energy, stores the new minumum and index to the chosen rotamer
                            energy_min=energy 
                            ip=i
                        endif
                    endif
                enddo

                if(ip>0) then !we found a good rotamer for the residue, so keep going!
                    call residue_replace_inplace(res_j, group_pep, backbone_backup, ip, aa_group_j) !put the best rotamer back in the peptide
                    call bindingenergy_full(group_pep, para_pep, binding_energy_new, vdw_new, ele_new, gb_new, np_new, pep_energy_new)   !calculate binding energy of system
                    call MC_sequence(binding_energy_new, binding_energy_old, vdw_new, vdw_old, &       !reject or accept amino acid swap using Metropolis Monte Carlo scheme
                                ele_new, ele_old, gb_new, gb_old, np_new, np_old, pep_energy_new, pep_energy_old, MC_feedback)

                    if (MC_feedback.eq.1) then !the amino acid swap was accepted!!
                        !update mutation attempt probabilities
                        SEQ_SWAP_success = SEQ_SWAP_success + 1
                        SEQ_PerRes_success(res_i) = SEQ_PerRes_success(res_i) + 1
                        SEQ_PerRes_success(res_j) = SEQ_PerRes_success(res_j) + 1
                        
                        !update TABU search info
                        !if (TABU_flag.eq.1) then
                        !    TABU_swap = TABU_allowed_residues(TABU_pick1)
                        !    TABU_allowed_residues(TABU_pick1) = TABU_allowed_residues(TABU_allowed_length)
                        !    TABU_allowed_residues(TABU_allowed_length) = TABU_swap
                        !    TABU_swap = TABU_allowed_residues(TABU_pick2)
                        !    TABU_allowed_residues(TABU_pick2) = TABU_allowed_residues(TABU_allowed_length-1)
                        !    TABU_allowed_residues(TABU_allowed_length-1) = TABU_swap
                        !    TABU_ban_lengths(TABU_allowed_length) = TABU_ban_length+1 !add 1 extra, since all values in array decreased by 1 before going back through loop
                        !    TABU_ban_lengths(TABU_allowed_length-1) = TABU_ban_length+1 !add 1 extra, since all values in array decreased by 1 before going back through loop
                        !    TABU_allowed_length = TABU_allowed_length-2             
                        !endif

                        !check if best peptide score should be updated
                        if (energy_min - binding_energy_old.lt.-.1) then
                            energy_min = binding_energy_old
                            !update simulated annealing best group info
                            if (SA_flag.eq.1) then
                                SA_best_pep = group_pep
                                SA_best_para = para_pep
                                SA_best_hyd = cur_hyd
                                SA_above_hyd = above_hyd
                                SA_below_hyd = below_hyd
                                SA_trp_count = trp_count
                            endif
                            open(2, file="minimum_energy.txt", access="append")
                                write(2,40) "cycle", cycle, "step", step, "attempt", attempt, "EBind", binding_energy_old
                            close(2)
                            if (pdb_flag.eq.0) then
                                write(pdb_name, "(I0,A1,I0,A1,I0)")  cycle, '_', step, '_', attempt
                                call writepdb_sys(group_pep, group_rec, pdb_name)
                            endif
                        endif
                    else !residue swap failed at Monte Carlo step
                        !restore energy parameters and coordinates
                        para_pep(res_i) = para_backup_i
                        para_pep(res_j) = para_backup_j
                        group_pep(res_i) = res_backup_i
                        group_pep(res_j) = res_backup_j

                        !tell program that peptide parameters and coordinates need to be reloaded
                        res_coords_flags(res_i) = 1
                        res_params_flags(res_i) = 1
                        GB_area_calc_flags(res_i) = 1
                        res_coords_flags(res_j) = 1
                        res_params_flags(res_j) = 1
                        GB_area_calc_flags(res_j) = 1
                        atom_links_flag = 1
                    endif

                else
                    !residue swap failed at finding allowable rotamer for residue 2 at residue 1's position
                    
                    !restore energy parameters and coordinates
                    para_pep(res_i) = para_backup_i
                    para_pep(res_j) = para_backup_j
                    group_pep(res_i) = res_backup_i
                    group_pep(res_j) = res_backup_j

                    !tell program that peptide parameters and coordinates need to be reloaded
                    res_coords_flags(res_i) = 1
                    res_params_flags(res_i) = 1
                    GB_area_calc_flags(res_i) = 1
                    res_coords_flags(res_j) = 1
                    res_params_flags(res_j) = 1
                    GB_area_calc_flags(res_j) = 1
                    atom_links_flag = 1
                endif

            else !replacement residue is not alanine, glycine, or proline, so use more general replacement method
                flag1=0
                do i=1, rotanum_2
                    call residue_replace_inplace(res_j, group_pep, backbone_backup, i, aa_group_j) !insert rotamer into the peptide
                    if(i==1) then
                        call energy_parameter_single(group_pep, para_pep, res_j)                       !load the force field paramters for new amino acid
                    endif

                    stage=0
                    call sidechain_optimization(group_pep, para_pep, res_j, stage, 0, 0, score)    !optimize sidechain, no polar_solvation energies
                    if(stage==1) then !sidechain was able to inserted into peptide!
                        if(flag1==0) then !this is the first rotamer that was successfully inserted into peptide
                            Tscore=score
                            res_j_best = group_pep(res_j)
                            flag1=1
                        else !another rotamer was successfully inserted, so determine which rotamer is better
                            if(score.le.Tscore) then
                                Tscore=score
                                res_j_best = group_pep(res_j)
                            endif
                        endif
                    endif
                enddo

                if(flag1==1) then !a rotamer successfully inserted! so keep going
                    group_pep(res_j) = res_j_best
                    call sidechain_optimization4binding(group_pep, para_pep, res_j, 0, CG_flag, CG_flags, energy_min)  !optimize sidechain of best rotamer a little more
                    call check_overlap(0, res_j, zero, group_pep, para_pep, feedback_1)                          !check for overlaps in system
                    if(feedback_1==1) then !no overlaps found! so keep going
                        call bindingenergy_full(group_pep, para_pep, binding_energy_new, vdw_new, ele_new, gb_new, np_new, pep_energy_new)                    !get binding energy after making the swap
                        call MC_sequence(binding_energy_new, binding_energy_old, vdw_new, vdw_old, &                       !accept or reject change using Metropolis Monte Carlo method
                        ele_new, ele_old, gb_new, gb_old, np_new, np_old, pep_energy_new, pep_energy_old,MC_feedback)

                        if (MC_feedback.eq.1) then !The amino acid swap was successful!

                            !update mutation attempt probabilities
                            SEQ_SWAP_success = SEQ_SWAP_success + 1
                            SEQ_PerRes_success(res_i) = SEQ_PerRes_success(res_i) + 1
                            SEQ_PerRes_success(res_j) = SEQ_PerRes_success(res_j) + 1
                            
                            !update TABU search info
                            !if (TABU_flag.eq.1) then
                            !    TABU_swap = TABU_allowed_residues(TABU_pick1)
                            !    TABU_allowed_residues(TABU_pick1) = TABU_allowed_residues(TABU_allowed_length)
                            !    TABU_allowed_residues(TABU_allowed_length) = TABU_swap
                            !    TABU_swap = TABU_allowed_residues(TABU_pick2)
                            !    TABU_allowed_residues(TABU_pick2) = TABU_allowed_residues(TABU_allowed_length-1)
                            !    TABU_allowed_residues(TABU_allowed_length-1) = TABU_swap
                            !    TABU_ban_lengths(TABU_allowed_length) = TABU_ban_length+1 !add 1 extra, since all values in array decreased by 1 before going back through loop
                            !    TABU_ban_lengths(TABU_allowed_length-1) = TABU_ban_length+1 !add 1 extra, since all values in array decreased by 1 before going back through loop
                            !    TABU_allowed_length = TABU_allowed_length-2   
                            !endif

                            if (energy_min - binding_energy_old.lt.-.1) then
                                energy_min = binding_energy_old
                                if (SA_flag.eq.1) then
                                    SA_best_pep = group_pep
                                    SA_best_para = para_pep
                                    SA_best_hyd = cur_hyd
                                    SA_above_hyd = above_hyd
                                    SA_below_hyd = below_hyd
                                    SA_trp_count = trp_count
                                endif
                                open(2, file="minimum_energy.txt", access="append")
                                    write(2,40) "cycle", cycle, "step", step, "attempt", attempt, "EBind", binding_energy_old
                                close(2)
                                if (pdb_flag.eq.0) then
                                    write(pdb_name, "(I0,A1,I0,A1,I0)")  cycle, '_', step, '_', attempt
                                    call writepdb_sys(group_pep, group_rec, pdb_name)
                                endif
                            endif
                        else !residue swap failed at Monte Carlo step
                            
                            !restore energy parameters and coordinates
                            para_pep(res_i) = para_backup_i
                            para_pep(res_j) = para_backup_j
                            group_pep(res_i) = res_backup_i
                            group_pep(res_j) = res_backup_j

                            !tell program that peptide parameters and coordinates need to be reloaded
                            res_coords_flags(res_i) = 1
                            res_params_flags(res_i) = 1
                            GB_area_calc_flags(res_i) = 1
                            res_coords_flags(res_j) = 1
                            res_params_flags(res_j) = 1
                            GB_area_calc_flags(res_j) = 1
                            atom_links_flag = 1
                        endif
                    else !residue swap failed at finding rotamer for residue 2 at residue 1's position without steric clash
                        !restore energy parameters and coordinates
                        para_pep(res_i) = para_backup_i
                        para_pep(res_j) = para_backup_j
                        group_pep(res_i) = res_backup_i
                        group_pep(res_j) = res_backup_j

                        !tell program that peptide parameters and coordinates need to be reloaded
                        res_coords_flags(res_i) = 1
                        res_params_flags(res_i) = 1
                        GB_area_calc_flags(res_i) = 1
                        res_coords_flags(res_j) = 1
                        res_params_flags(res_j) = 1
                        GB_area_calc_flags(res_j) = 1
                        atom_links_flag = 1
                    endif
                else  !residue swap failed at finding allowable rotamer for residue 2 at residue 1's position
                    !restore energy parameters and coordinates
                    para_pep(res_i) = para_backup_i
                    para_pep(res_j) = para_backup_j
                    group_pep(res_i) = res_backup_i
                    group_pep(res_j) = res_backup_j

                    !tell program that peptide parameters and coordinates need to be reloaded
                    res_coords_flags(res_i) = 1
                    res_params_flags(res_i) = 1
                    GB_area_calc_flags(res_i) = 1
                    res_coords_flags(res_j) = 1
                    res_params_flags(res_j) = 1
                    GB_area_calc_flags(res_j) = 1
                    atom_links_flag = 1
                endif
            endif
        else !residue swap failed at finding allowable rotamer for residue 2 at residue 1's position
            
            !restore energy parameters and coordinates
            para_pep(res_i) = para_backup_i
            para_pep(res_j) = para_backup_j
            group_pep(res_i) = res_backup_i
            group_pep(res_j) = res_backup_j

            !tell program that peptide parameters and coordinates need to be reloaded
            res_coords_flags(res_i) = 1
            res_params_flags(res_i) = 1
            GB_area_calc_flags(res_i) = 1
            res_coords_flags(res_j) = 1
            res_params_flags(res_j) = 1
            GB_area_calc_flags(res_j) = 1
            atom_links_flag = 1
        endif

        !output result of sequence optimization attempt
        open(3, file="energydetails.txt",  access="append")
        write(3,4) cycle,  step, attempt, binding_energy_old, (binding_energy_old-pep_energy_old), "Swap"
        !write(3,"(<pep_res>(a4))") (group_pep(j)%gtype, j=1, pep_res)
        write(3,*) (group_pep(j)%gtype, j=1, pep_res)
        write(3,*) "*******************************"
        close(3)
    endif
    !if (TABU_flag.eq.1.and.MC_feedback.eq.1) call update_TABU_ban_list() !only update TABU list whenever a sequence change is accepted
enddo
4   format(i4, i7, i7, 2f20.4, a15)
40  format(a8, i5, a8, i5, a10, i3, a8, f9.3)

return
end subroutine sequence_optimization

subroutine backbone_optimization(group_pep, original_pep, group_rec, para_pep, SA_best_pep, backbone_backup, step, cycle, binding_energy_old, &
    vdw_old, ele_old, gb_old, np_old, pep_energy_old)
implicit none
integer                    :: step, attempt, i, j, k, phipsi_num, cycle
integer                    :: feedback_1, MC_feedback, flag, flag_Nterm, flag_Cterm
integer                    :: ic(3), ic_min, ic_flag, ran_resi(3)
integer                    :: new_min_flag
integer                    :: GB_flag, CG_flag, CG_flags(pep_res)
character*10               :: move_type
character*50               :: pdb_name
double precision           :: binding_energy_old, vdw_old, ele_old, gb_old, np_old, pep_energy_old
double precision           :: binding_energy_min, vdw_min, ele_min, gb_min, np_min, pep_energy_min
double precision           :: Tbinding_energy_min, Tvdw_min, Tele_min, Tgb_min, Tnp_min, Tpep_energy_min
real                       :: ran2
type(groupdetails)         :: group_pep(pep_res), group_rec(rec_res), SA_best_pep(pep_res)
type(databackup)           :: backbone_backup(pep_res) 
type(databackup)           :: backbone_backup_candidates(pep_res, 20)
type(groupdetails)         :: group_candidates(pep_res, 20)
type(groupdetails)         :: pep_backup(pep_res), pep_best(pep_res), original_pep(pep_res)
type(energyparameters)     :: para_pep(pep_res)
type(databackup)           :: backbone_backup_backup(pep_res), backbone_backup_best(pep_res)

print ("(A,I0,A)"), "Step ", step, ": Backbone optimization"

new_min_flag = 0
attempt=1

!at present, disable sidechain coarse-graining and polar solvation energy
GB_flag = 0; CG_flag = 0; CG_flags = 0

!decide what kind of CONROT move to do (all internal, one termini, or both termini)
call RANDOM_NUMBER(ran2)
if(ran2.le.CONROT_two_term) then !two of the pivots for CONROT are the peptide termini
    flag_Nterm=1
    flag_Cterm=1
    ic(1)=1
    ic(2)=fragmentnum
    ic(3)=int(ran2*(fragmentnum-2))+2 !this expression picks any residue not at N- or C-terminus

elseif(ran2.le.CONROT_one_term)   then  !one of the pivots will be a terminus
    call RANDOM_NUMBER(ran2)
    if(ran2.le.0.5) then  !will use the N-terminus for the terminal pivot
        flag_Nterm=1
        flag_Cterm=0
        ic(1)=1
    else !will use the C-terminus for the terminal pivot
        flag_Nterm=0
        flag_Cterm=1
        ic(1)=fragmentnum
    endif

    call RANDOM_NUMBER(ran2)
    ic(2)=int(ran2*(fragmentnum-2))+2 !this expression picks any residue not at N- or C-terminus
    flag = 0
    do while(flag.eq.0)
        call RANDOM_NUMBER(ran2)
        ic(3)=int(ran2*(fragmentnum-2))+2
        if (ic(3).ne.ic(2)) then
            flag = 1 
        endif
    enddo

else  !pick three internal peptide residues randomly
    flag_Nterm=0
    flag_Cterm=0
    call RANDOM_NUMBER(ran2)
    ic(1)=int(ran2*(fragmentnum-2))+2 !this expression picks any residue not at N- or C-terminus
    flag = 0 
    do while(flag.eq.0)
        call RANDOM_NUMBER(ran2)
        ic(2)=int(ran2*(fragmentnum-2))+2 !this expression picks any residue not at N- or C-terminus
        if (ic(2).ne.ic(1)) then
            flag = 1
        endif
    enddo

    flag = 0
    do while(flag.eq.0)
        call RANDOM_NUMBER(ran2)
        ic(3)=int(ran2*(fragmentnum-2))+2 !this expression picks any residue not at N- or C-terminus
        if(ic(3).ne.ic(1).and.ic(3).ne.ic(2)) then
            flag=1
        endif
    enddo
endif

!order the three residues in res_i from smallest to largest
do i=1, 2
    ic_min=ic(i)
    ic_flag=i
    do j=i+1, 3
        if(ic_min.gt.ic(j)) then
        ic_min=ic(j)
        ic_flag=j
        endif
    enddo
    ic(ic_flag)=ic(i)
    ic(i)=ic_min
enddo

do i=1, 3
    ran_resi(i)=residuesite(ic(i))
enddo

!backup current peptide coordinates
pep_backup(ran_resi(1) : ran_resi(3)) = group_pep(ran_resi(1) : ran_resi(3))

if(flag_Nterm.ne.1.and.flag_Cterm.ne.1) then !doing an internal CONROT move
    CONF_Internal_attempts = CONF_Internal_attempts+1
    move_type = "CRA_mid"
    call concertedrotation_center(group_pep, original_pep, backbone_backup, ran_resi, phipsi_num, group_candidates, backbone_backup_candidates)
elseif(flag_Nterm.eq.1.and.flag_Cterm.eq.1) then !doing CONROT move where both peptide termini are pivots
    CONF_N2C_attempts = CONF_N2C_attempts + 1
    move_type = "CRA_whole"
    call concertedrotation_whole(group_pep,original_pep, backbone_backup, ran_resi, phipsi_num, group_candidates, backbone_backup_candidates)
elseif(flag_Nterm.eq.1.and.flag_Cterm.ne.1) then !N-terminus will be one of the three pivot points for CONROT
    CONF_N_attempts = CONF_N_attempts + 1
    move_type = "CRA_Nterm"
    call concertedrotation_Nterm(group_pep, original_pep, backbone_backup, ran_resi, phipsi_num, group_candidates, backbone_backup_candidates)
elseif(flag_Nterm.ne.1.and.flag_Cterm.eq.1) then
    CONF_C_attempts = CONF_C_attempts + 1
    move_type = "CRA_Cterm"
    call concertedrotation_Cterm(group_pep, original_pep, backbone_backup, ran_resi, phipsi_num, group_candidates, backbone_backup_candidates)
endif

if(phipsi_num.gt.0) then !found a possible CONROT move or moves, so keep going
    binding_energy_min=10000.0
    backbone_backup_backup = backbone_backup
    flag = 0
    do k=1, phipsi_num
        group_pep=group_candidates(:,k)
        res_coords_flags(ran_resi(1):ran_resi(3)) = 1; GB_area_calc_flags(ran_resi(1):ran_resi(3)) = 1
        backbone_backup=backbone_backup_candidates(:,k)
        call check_backbone_overlap(group_pep, para_pep, feedback_1) !make sure backbone doesn't have overlaps
        !call calc_com(group_pep, pep_com, pep_res) !check backbone change doesn't make peptide bend too far from surface
        if(feedback_1.ne.0 ) then !.and. (pep_com(3)-start_pep_com(3)).lt.max_zdist_from_surface) then !no overlaps and within allowable distance from receptor, so keep going
            !call sidechain_repacking_vdwele_dihedral(group_pep, para_pep, backbone_backup, ran_resi(1), ran_resi(3), 0) !do initial cheaper/faster repacking of sidechains, no polar solvation energies included
            !call check_overlap(2, zero, zero, group_pep, para_pep, feedback_1)
            !if(feedback_1.ne.0) then
                call sidechain_repacking(group_pep, para_pep, backbone_backup, GB_flag, CG_flag, CG_flags, Tbinding_energy_min, Tvdw_min, Tele_min, Tgb_min, Tnp_min, Tpep_energy_min)
                if(Tbinding_energy_min.lt.binding_energy_min) then
                    flag = 1
                    pep_best = group_pep
                    backbone_backup_best = backbone_backup
                    binding_energy_min=Tbinding_energy_min
                    vdw_min=Tvdw_min
                    ele_min=Tele_min
                    gb_min=Tgb_min
                    np_min=Tnp_min
                    pep_energy_min= Tpep_energy_min
                endif
            !endif
        endif
    enddo

    if (flag.eq.1) then !a good option found in the CONROT moves, so now move on to Metropolis check
        call MC_backbone(binding_energy_min, binding_energy_old, vdw_min, vdw_old, ele_min, ele_old, &
        gb_min, gb_old, np_min, np_old, pep_energy_min, pep_energy_old, MC_feedback)
        if(MC_feedback==1) then
            !put in best peptide and backup coordinates
            group_pep = pep_best
            backbone_backup = backbone_backup_best

            flag4conformer=1 
            rmsd_flag = 1 !indicate the backbone has changed so the rmsd should be calculated again
            res_coords_flags(ran_resi(1):ran_resi(3)) = 1; GB_area_calc_flags(ran_resi(1):ran_resi(3)) = 1 !tell program to reload peptide coordinates before calculating an energy

            !always generate pdb whenever a backbone change is made
            write(pdb_name, "(I0,A1,I0,A1,I0)")  cycle, '_', step, '_', attempt
            call writepdb_sys(group_pep, group_rec, pdb_name)

            !update probabilities of making successful CONROT move
            if(flag_Nterm.ne.1.and.flag_Cterm.ne.1) then 
                CONF_Internal_success = CONF_Internal_success+1
            elseif(flag_Nterm.eq.1.and.flag_Cterm.eq.1) then
                CONF_N2C_success = CONF_N2C_success + 1
            elseif(flag_Nterm.eq.1.and.flag_Cterm.ne.1) then 
                CONF_N_success = CONF_N_success + 1
            elseif(flag_Nterm.ne.1.and.flag_Cterm.eq.1) then
                CONF_C_success = CONF_C_success + 1
            endif

            !check if this is the best score found so far
            if (binding_energy_old.lt.energy_min-.01) then
                energy_min = binding_energy_old
                !if doing simulated annealing, store current group in SA_best
                if (SA_flag.eq.1) then
                    SA_best_pep = group_pep
                    SA_best_hyd = cur_hyd
                    SA_above_hyd = above_hyd
                    SA_below_hyd = below_hyd
                    SA_trp_count = trp_count
                endif
                
                !write data to minimum energy file
                open(2, file="minimum_energy.txt", access="append")
                    write(2,41) "cycle", cycle, "step", step, "attempt", attempt, "EBind", binding_energy_old
                close(2)
            endif
        else !Metropolis criterion not satisfied
            
            !restore original peptide coordinates and backup coordinates
            group_pep(ran_resi(1) : ran_resi(3)) = pep_backup(ran_resi(1) : ran_resi(3))
            backbone_backup= backbone_backup_backup

            !tell the program that the peptide coordinates need to be reloaded
            res_coords_flags(ran_resi(1):ran_resi(3)) = 1; GB_area_calc_flags(ran_resi(1):ran_resi(3)) = 1
        endif
    else !none of the CONROT options were good
        
        !restore original peptide coordinates
        group_pep(ran_resi(1) : ran_resi(3)) = pep_backup(ran_resi(1) : ran_resi(3))
        backbone_backup = backbone_backup_backup

        !tell the program that the peptide coordinates need to be reloaded
        res_coords_flags(ran_resi(1):ran_resi(3)) = 1; GB_area_calc_flags(ran_resi(1):ran_resi(3)) = 1
    endif
else
    !no CONROT options found, so restore original coordinates
    group_pep(ran_resi(1) : ran_resi(3)) = pep_backup(ran_resi(1) : ran_resi(3))  
    
    !tell the program that the peptide coordinates need to be reloaded
    res_coords_flags(ran_resi(1):ran_resi(3)) = 1; GB_area_calc_flags(ran_resi(1):ran_resi(3)) = 1
endif

!write result of CONROT move (do this even when no options were found)
open(3, file="energydetails.txt",  access="append")
write(3,4) cycle, step, attempt, binding_energy_old, (binding_energy_old-pep_energy_old), move_type
write(3,*) (group_pep(j)%gtype, j=1, pep_res)
write(3,*) "*******************************"
close(3)

4   format(i5,i7,i7,2f20.13,a15)
41   format(a8, i5, a8, i5, a10, i3, a8, f9.3)
return
end subroutine backbone_optimization

subroutine rigid_body_optimization(group_pep, group_rec, para_pep, SA_best_pep, backbone_backup, step, cycle, binding_energy_old, vdw_old, ele_old, gb_old, np_old, pep_energy_old)
implicit none

integer                      :: step, attempt, j, cycle, cg_flags(pep_res)
integer                      :: feedback
integer                      :: move_type
character*50                 :: pdb_name
double precision             :: binding_energy_old, vdw_old, ele_old, gb_old, np_old, pep_energy_old
double precision             :: binding_energy_new, vdw_new, ele_new, gb_new, np_new, pep_energy_new
real                         :: ran2
type(groupdetails)           :: group_pep(pep_res), group_rec(rec_res), pep_backup(pep_res), SA_best_pep(pep_res)
type(databackup)             :: backbone_backup(pep_res) 
type(energyparameters)       :: para_pep(pep_res)

pep_backup = group_pep(1:pep_res)

!pick either a translation or rotation rigid body move
call RANDOM_NUMBER(ran2)
if (ran2.lt.rotation_switch) then
    call random_peptide_translation(group_pep)
    move_type=1
else
    call random_peptide_rotation(group_pep)
    move_type=2
endif

!create atom links for proper evaluation of binding energy
cg_flags = 0
atom_links_flag = 1

!calculate new binding energy and determine if rigid body move should be accepted
call bindingenergy_full(group_pep, para_pep, binding_energy_new, vdw_new, ele_new, gb_new, np_new, pep_energy_new)
call MC_backbone(binding_energy_new, binding_energy_old, vdw_new, vdw_old, &
    ele_new, ele_old, gb_new, gb_old, np_new, np_old, pep_energy_new, pep_energy_old, feedback)

if (feedback.eq.1) then !change accepted
    !update peptide center of mass position (only necessary for translation, since rotations are done about peptide center of mass)
    !cur_pep_com = cur_pep_com + translate_vec !currently not using this option, since center of mass not updated after backbone or sequence changes

    !update backup coordinates (necessary for mutations involving proline)
    call update_backup_coords(group_pep, backbone_backup)

    !write a pdb file
    write(pdb_name, "(I0,A1,I0,A1,I0)")  cycle, '_', step, '_', 1
    call writepdb_sys(group_pep, group_rec, pdb_name)

    !check if we found the minimum energy during design
    if (binding_energy_old.lt.energy_min-.01) then
        energy_min = binding_energy_old
        if (SA_flag.eq.1) then
            SA_best_pep = group_pep
            SA_best_hyd = cur_hyd
            SA_above_hyd = above_hyd
            SA_below_hyd = below_hyd
            SA_trp_count = trp_count
        endif

        open(2, file="minimum_energy.txt", access="append")
            write(2,40) "cycle", cycle, "step", step, "attempt", attempt, "EBind", binding_energy_old
        close(2)
    endif
else
    group_pep = pep_backup
endif   
40  format(a8, i5, a8, i5, a10, i3, a8, f9.3)

open(3, file="energydetails.txt",  access="append")
    if(move_type.eq.1) then
        write(3,4) cycle, step, binding_energy_old, (binding_energy_old-pep_energy_old), "        TRANSLATE" 
    else
        write(3,4) cycle, step, binding_energy_old, (binding_energy_old-pep_energy_old), "        ROTATE" 
    endif
    write(3,*) (group_pep(j)%gtype, j=1, pep_res)
    write(3,*) "*******************************"
close(3)
4   format(i4, i7, 2f20.13, a15)

end subroutine rigid_body_optimization

end module optimization_techniques