module sys_vars

use datatypes

implicit none

!File names and paths
character*100                    :: pep_pdb='pep.pdb' !path to pdb file of isolated peptide, if used by program
character*100                    :: rec_pdb = 'rec.pdb' !path to pdb file of isolated receptor, if used by program
character*100                    :: system_pdb = 'sys.pdb' !path to pdb file of peptide-receptor complex, if used by program
character*100                    :: restart_pdb='restart.pdb'  !path to pdb file of peptide-receptor complex for restarting design, if used by program
character*100                    :: lib_path='~/PepBD/lib/' !path to PepBD library files. Defaults to ~/PepBD/lib
character*100                    :: pdb_path = './pdbfiles/' !path to directory where generated pdb files will be placed
character*100                    :: backbones_file = 'backbones.txt' !file containing paths to different peptide pdb files to be evaluated in quantum annealing programs
character*100                    :: backbones_path = 'backbones.txt' !file containing paths to different peptide pdb files to be evaluated in quantum annealing programs
character*100                    :: sequences_file = 'sequences.txt' !file containing different peptide sequences to be evaluated in quantum annealing programs
logical                          :: system_found, pep_found, rec_found !flags that indicate if a pdb file has been found

!system variables
integer                          :: gnum  ! Total number of monomers in the system (peptide + receptor)
integer                          :: atom_num, atom_num_rec, atom_num_pep !stores info on all atoms in the system, peptide, and receptor
integer,parameter                :: atom_per_pep_res=30 !max number of atoms in a peptide residue (based on tryptohan at N-terminus, which is actually only 26 atoms but 30 is a nice round number)
integer                          :: pep_res, rec_res !number of residues in the peptide and receptor
integer  ,allocatable            :: residuesite(:) !stores the location of helices in the peptide (I think?)
integer, allocatable             :: res_starts(:), res_ends(:), sc_starts(:), sc_ends(:) !stores start and end indices of each residue in peptide
integer                          :: fragmentnum !number of freely rotatable residues in the peptide
integer                          :: helix_num !number of helices in the peptdie
integer  , parameter             :: maximum_helix_num=5 !max number of helices in the peptide

!peptide variables
integer                          :: helix_start(maximum_helix_num), helix_end(maximum_helix_num) !residues numbers for start/end of all helices
real            , allocatable    :: alpha_pep_sys(:,:), alpha_pep_iso(:,:) !stores parameters for system's born radii. _sys is for complex, _iso is for isolated
integer  , allocatable           :: numex_pep(:,:), numex4_pep(:,:), inb_pep(:,:,:,:), inb4_pep(:,:,:,:), atomid_pep(:,:)
real            , allocatable    :: mdcrd_pep(:,:,:), charge_pep(:,:), epsion_pep(:,:)
real            , allocatable    :: r_pep(:,:), rborn_pep(:,:), fs_pep(:,:), dielecons_pep(:,:)
integer, allocatable             :: atoms_per_res(:)

!receptor variables
real            , allocatable    :: mdcrd_rec(:,:), charge_rec(:), epsion_rec(:)
real            , allocatable    :: alpha_rec_sys(:), alpha_rec_iso(:) !stores parameters for system's born radii. _sys is for complex, _iso is for isolated
real            , allocatable    :: r_rec(:), rborn_rec(:), fs_rec(:), dielecons_rec(:)
integer  , allocatable           :: LCPO_indices(:)
real                             :: pep_com(3), pol_com(3) !center of mass of polymer and peptide
double precision                 :: rec_gb=0.0, rec_np=0.0 !stores energy info on the receptor 
real                             :: rec_max_z=-10000, pep_min_z=1000000 !stores z coordinate of surface and minimum z coordinate of peptide
real, allocatable, dimension(:)  :: rec_top_z_arr
integer                          :: num_atoms
real                             :: net_mass_backbone, net_mass_rec, net_mass_pep
real, allocatable                :: phis(:), psis(:) !stores the phi and psi dihedral angles for all peptide residues
character*4,allocatable          :: lbres(:) !stores residue types
real                             :: r_gyr   !radius of gyration for the peptide
real                             :: start_pep_com(3) !stores starting center of mass coordinates of the peptide
integer, allocatable             :: res_coords_flags(:), res_params_flags(:) !tells program if coords/parameters need to be reloaded for system. 1 = need to reload; 0 = no reload needed
real, allocatable                :: rij_matrix_pep(:,:,:,:) !stores pairwise distances between atoms in system. Entry (atom_j, atom_i, res_j, res_i) corresponds to rij distance between atom_j in residue res_j of peptide and atom_i of residue res_i in peptide
real, allocatable                :: rij_matrix_sys(:,:,:) !stores pairwise distances between atoms in system. Entry (atom_j, atom_i, res_i) corresponds to rij distance between atom_j in receptor and atom_i of residue res_i in peptide
real, allocatable                :: rij_matrix_rec(:,:) !stores pairwise distances between atoms in system. Entry (atom_j, atom_i) corresponds to rij distance between atom_j and atom_i in receptor
integer                          :: atom_links_flag !tells program if peptide atom bonds need to be created before calculating energies
type(atomlink), allocatable      :: atom(:,:) !this stores information used to form the exclusions lists 
type(groupdetails), allocatable  :: Tgroup_pep(:) !backup for main group variable
type(lib4aminoacid)              :: aa_lib(20) !this stores the default amino acid rotamer for all amino acid types
type(LCPO_params)                :: LCPO(21) !parameters for calculating non-polar solvation energy from surface accessible area via LCPO method
type(groupdetails)               :: aa_group_i(40), aa_group_j(40) !stores rotamer info for individual amino acids 
type(groupdetails), allocatable  :: aa_groups(:,:) !this stores rotamer options for all residues in the peptide

!generalized born (GB) solvation energy calculation parameters
real, allocatable                :: integrals_rec_iso(:) !stores atomic shielding for isolated receptor, used when calculating generalized born radii
real, allocatable                :: integrals_pep_iso(:,:) !stores atomic shielding for isolated receptor, used when calculating generalized born radii
real, allocatable                :: integrals_comp_rec_from_pep(:,:) !stores atomic shielding of receptor due to peptide; ; index into array is (receptor atom, peptide res)
real, allocatable                :: integrals_comp_pep_from_rec(:,:) !stores atomic shielding of each peptide atom due to receptor; index into array is (atom,res)
integer, allocatable             :: GB_area_calc_flags(:) !indicates if a residue's coordinates have changed so the GB integrals need to be recalculated 

!ramachandran map variables
integer                          :: flavoredregion_number, rama_map(-179:180, -179:180) ! The defined varivables are for the ramachandran map.
real                             :: flavored_region(60000,2) ! tuples of allowed ramachandran values

!General constants
real            , parameter      :: PI=3.14159265358979323846264, TWO_PI = 6.28318530718, FOUR_PI=12.5663706144 !constant pi
integer                          :: one=1, zero=0 !used in a couple locations to avoid integer size errors
real                             :: one_real=1, zero_real=0
real                             :: deg_2_rad = 0.01745329251
real                             :: dielecons4water=80.0 !dielectric constant for water, used in GB solvation energy calculation

!design flags
integer                          :: initial_randomize=1
integer                          :: debug_flag=0
integer                          :: solvation_flag = 0 !if 1, then will calculate solvation energies
integer                          :: end_point_flag = 1 !if 1, then peptides scored based on energy difference in bound vs. unbound state; if 0, then peptides scored on interaction energy
integer                          :: move_pep_over_surface_flag = 0 !if 1, then translate peptide in z-direction over surface before starting design. if 0, then do nothing
integer                          :: restart_flag = 0 !if set to 1, then we are restarting design from some part in middle
integer                          :: NP_flag = 0 !determines if surface area energy calculated in design. 0: do not calculate (default); 1: do calculate 
integer                          :: print_flag = 0
integer                          :: rmsd_flag = 0 !tells program if rmsd should be calculated during design. 1 = calculate, 0 = don't calculate 
integer                          :: SA_flag = 0 !default to no simulated annealing
integer                          :: decomp_mutate_flag = 0 !if this flag is 1, then select residues to mutate with likelihood based on their contribution to binding energy
integer                          :: neighbor_flag = 0 !if 1, then use neighbor cutoff. If 0, then don't use neighbor cutoff
integer                          :: TABU_flag = 0 !if 1, then tabu search used for sequence changes
integer                          :: smart_restart_flag = 0 !controls whether program actively tracks if a new cycle should be started; Not used at present 
integer                          :: eval_twobody_flag=1 !if 1, then two-body (i.e. peptide-peptide) energies will be calculated (not used for "EvalPep" program)
integer                          :: eval_onebody_flag=1 !if 1, then one-body (i.e. peptide-receptor) energies will be calculated (not used for "EvalPep" program)
integer                          :: check_init_structure = 1 !if 1, then check initial structure for steric overlaps with receptor. if 0, then do not check
integer                          :: overlap_type = 1 !determines method of checking for overlaps. Not used at present

!design variables
real                             :: target_pep_distance = 4.5 !if translating peptide over surface, this value is target distance between peptide center of mass and top of surface 
integer                          :: num_ex_AAs = 0 !counts number of excluded amino acids
character*4                      :: AAs_available(10,6) = '' !contains amino acids that can be used during design, where each column corresponds to a certain amino acid category. 
                                                    !As with hydration properties, order is: glycine, hydrophobic, positive, hydrophilic, other, and negative
integer                          :: num_AAs_available(6) = 0 !contains length of each column for the AAs_available array
character*4                      :: excluded_AAs(20) !residues that will not be allowed during design
real                             :: ph_value=7        ! The setting of PH
integer                          :: max_repack_iters = 5
integer                          :: num_backbones
integer  , parameter             :: number4peptide_candidates=6 !max number of conformation changes found by CONROT before it evaluates which of the 6 gives the lowest energy
integer                          :: steps_before_conf_change = 0 !number of steps taken in each cycle before conformation changes are allowed
double precision                 :: energy_min=10000 !stores minimum energy of peptide found during design (set to an arbitrary large positive number initially)
real                             :: rmsd_max    !The maximum restraints of rmsd
real                             :: rmsd_min    !The minimum restraints of rmsd
integer                          :: PDB_gen_freq = -1 !PDB files generated after this many steps during design. Default value means only structures of peptides better than all previous peptides will be stored
real                             :: delta_phi_min = -6.0, delta_phi_max = 6.0, delta_phi_step = 1.0
integer                          :: ran_seed(1) = -1
double precision, parameter      :: e_keep_sc_cutoff = 50 !only keep sidechain rotamer during rotamer search if its interaction with receptor is less than this value
double precision, parameter      :: vdw_cutoff = 20 !if vdw energy between two atoms exceeds this value, then assume there is an overlap/steric clash
real, parameter                  :: min_rij2_ratio = 2.8 !if square distance between two atoms times this value is less than square of van der Waals radii, then we have a steric overlap
double precision, parameter      :: e_sidechain_opt_cut = 0.1 !when optimizing sidechain dihedrals, optimization stops once change in energy falls below this value
double precision, parameter      :: pep_sc_repack_cutoff_energy = 10 !when repacking peptide sidechains, not considering interactions with the receptor
integer                          :: insert_pep_rigor = 1 !determines how hard we try to calculate energy of given peptide sequence by finding best sidechain packing over different backbones
character*30                     :: receptor_name
real, parameter                  :: overlap_dist = 2.4025 ! 2.4025 = 1.55^2
real, parameter                  :: overlap_score = 102.1
integer                          :: max_bb_failed_attempts = 40
real                             :: max_QUBO_score = 50
integer                          :: nstep_start, nstep_terminal, idum,  rec_calculation_flag, flag4conformer
real                             :: ekt_seq = 1, ekt_backbone_ratio = 1,  ekt_backbone,  ekt_rigid_ratio =1,  ekt_rigid, ekt_scmf=0.6 !reference temperatures for sequence and conformation changes, and temperature for SCMF method. These can be changed in input file
real            , parameter      :: lampda=0.5  ! The convergence coefficient
real                             :: weighting_factor=0.01  ! weighting factor is used to adjust the importance to the ligand within the score function
real            , parameter      :: dihedral_weighting_factor=0.50 !weights dihedral energy in score function
real            , parameter      :: vdw14_coeff=2.0  !1-4 VDW scale factor
real            , parameter      :: ele14_coeff=1.2  !1-4 ELE scale factor
real                             :: translate_max = 2.0, translate_min = 0.5 !min/max translation amount in single step
real                             :: rotate_max = pi !max rotation amount
real                             :: max_dev_from_center = 5 !number of Angstroms peptide's center of mass can deviate from original peptide center of mass during design
real                             :: min_zdist_from_surface = 0, max_zdist_from_surface = 1000 !minimum and maximum allowed distance between peptide center of mass and surface (reqires surface normal to be in positive z-direction)
integer                          :: verbose = 0 !determines how much info will be printed out to terminal during program execution              

!variables for non-polar solvation energy
real                             :: probe_r = 1.4 !radius of probe water used to determine solvent accesible surface area (not currently used)
real                             :: surftens_solv = 0.06  !The surface tension in the bulk solvent for purely hydrophobic cavity, in kcal/mol/Angstrom-squared
real                             :: surftens_rec  = 0.03  !The surface tension of the receptor in the solvent for a purely hydrophobic cavity, in kcal/mol/Angstrom-squared
real            , parameter      :: SASA_offset=0.0 !constant added to surface tension energy. Since we are interested in relative binding energy, can set to 0
integer                          :: max_neighbors = 60 !storage allocated for neighbors of each heavy atom

!other parameters
integer  ,parameter              :: num_AAs = 19
character*4                      :: AA_names(num_AAs), AA_name_i, AA_name_j, res_name
integer                          :: backbone_changes_per_iter=10

!Sidechain coarse-grained bead parameters (used for rapid optimization of sidechains when solvation energies are included)
integer                             :: cg_flag=0 !if 1, then will use coarse-graining in program for peptide sidechains
real                                :: cg_rborn = 3
real                                :: cg_f = 0.8
real                                :: cg_dielec=0.1
real                                :: cg_epsilon = 0 !don't include van der Waals interactions to the sidechain
real                                :: cg_q = 0 !don't include electrostatic interactions to the sidechain
real                                :: cg_r = 0 !don't include van der Waals interactions to the sidechain

!parameters for positioning peptide on polymer surface
real                             :: limits(2,3)
real                             :: ranges(3)
real                             :: x_buffer = 4, y_buffer=4, z_min=2, z_max=15
real                             :: x_target, y_target, z_target
real                             :: alpha, beta, gamma !values for the peptide Euler angles

!Parameters for fixing the number of tryptophans
integer                  :: trp_limit = 3, trp_count !default to only allowing 3 tryptophans in a peptide

!Parameters for variable hydration
integer                  :: min_hyd(6), max_hyd(6), cur_hyd(6), above_hyd(6), below_hyd(6), temp_hyd(6) !order is glycine, hydrophobic, positive, hydrophilic, other, and negative

!probabilities governing peptide changes
real                     :: backbone_prob=0.5, backbone_switch  !Probability of doing a backbone change, and corresponding number for 
real                     :: seq_prob = 0.5, seq_switch !Probability of doing a sequence change.
real                     :: res_mutate_prob=0.8  ! probabililty of doing a residue mutation. (1-scmfswitch) is probability of exchanging two residues
real                     :: CONROT_two_term = 0.1, CONROT_one_term = 0.4
real                     :: rotation_switch = 0.5 !probability of doing rigid body rotation, vs. a rigid body translation

!Variables for tracking mutation probability success. Purposes should be apparent from the variable names
integer                  :: SEQ_SWAP_attempts=0, SEQ_POINT_attempts=0, CONF_N2C_attempts=0, CONF_N_attempts=0, CONF_C_attempts=0
integer                  :: CONF_Internal_attempts=0, rotate_attempts=0, translate_attempts=0
integer  , allocatable   :: SEQ_PerRes_attempts(:), SEQ_PerRes_success(:)
integer                  :: SEQ_SWAP_success=0, SEQ_POINT_success=0, CONF_N2C_success=0, CONF_N_success=0, CONF_C_success=0
integer                  :: CONF_Internal_success=0, rotate_success=0, translate_success=0

!Parameters for restart method; NOTE: this method not recommended, doesn't make design better
integer                   :: steps_before_randomize=50 !how many steps will be taken with initial sequence before randomizing peptide, if initial_randomize == 0
integer                   :: randomize_length=10 !number of chnages to make when randomizing peptide 
!integer                  :: num_restarts, steps_per_cycle, randomize_length=10 !number of restarts in 
!integer                  :: initial_randomize = 0 !determines whether peptide sequence is randomized before starting design; 1 = randomize, 0 = don't randomize
!integer                  :: smart_restart_len = 100 !number of steps to average the scores over
!real                     :: smart_restart_offset=5.0 !max offset between current score and best score to keep searching
!integer                  :: smart_restart_wait_len = 5 !wait this many steps after randomizing sequence to start accumulating score (scores at beginning will be very high)
!integer                  :: smart_restart_score_method = 1 !method for determining when to restart. 1 = compare average scores of intervals, 2 = compare best scores of intervals, 3 = compare average score of current interval to best score

!parameters for simulated annealing --> currently using geometric series for temeprature schedule (T_{t+1}  = T_{t} * r, 0 < r < 1)
integer                  :: steps_per_cycle
integer                  :: SA_interval = 100 !after SA_interval steps, reduce the temperature
real                     :: SA_fraction = 0.9 !value of r
real                     :: SA_T_start = 10, SA_T_end = 0.5 !default starting and end temperature
integer                  :: SA_num_cycles !number of intervals needed to go from SA_T_start to SA_T_end (will overshoot if necessary)
integer                  :: SA_start_cycle = 1 !set to 1 if starting design from beginning, otherwise set to cycle at which design should start
integer                  :: SA_best_hyd(6), SA_above_hyd(6), SA_below_hyd(6), SA_trp_count !stores hydration properties of best peptide in current cycle

!variables for per-residue binding energy
double precision, allocatable      :: energy_per_res(:, :) !stores binding energy per residue
!real            , allocatable      :: res_select_p_vector(:) !stores probability of selecting a residue to be mutated
!character*4                        :: decomp_score_type = "INV" 

!variables for neighbor list (currently not used)
real                     :: neighbor_cutoff2 = 9999 !beyond this distance, two atoms not considered neighbors and will not calculate pair-wise energy. By default, this value is large so all atoms are neighbors
!integer                 :: neighbor_length = 200 !should be more than enough to store all neighbors for each residue in the peptide
!integer  , allocatable  :: neighbor_list(:,:)   !will store index to neighbor residues for each residue in the ligand (NOTE: doing this on a per-residue basis) 
!integer  , allocatable  :: neighbor_list_lens(:) !will store how many neighbors each ligand residue has
!real *8, allocatable    :: neighbor_rec_res_box_edges(:,:) !will store edges of each receptor's box edges to determine distance from each ligand residue' box edges

!Tabu search parameters; NOTE: Not recommended to use this, doesn't make design more effective
!integer, allocatable    :: TABU_allowed_residues(:) !list of amino acid positions that CAN be mutated
!integer                 :: TABU_ban_length = 0 !number of mutations that mutated residue cannot be changed. NOT EQUAL TO STEPS!! There are n sequence mutations per step, where n is the number of amino acids in the peptide
!integer, allocatable    :: TABU_ban_lengths(:) !stores list of how long tabu residues have been on the tabu list. list decremented, so residue can be mutated again once the length = 0
!integer                 :: TABU_allowed_length !stores number of residues currently allowed by TABU search

!genetic algorithm docking parameters 
real                            :: mutation_probability = 0.25
real                            :: fraction_retain = 0.4
integer                         :: num_retain 
integer                         :: dim
real                            :: min_improvement = 2.0
integer                         :: patience = 10 !max number of evolution rounds with no improvement before docking ends
real                            :: z_penalty = 0 !multiply z distance between peptide and top of surface by this value and add to score


end module sys_vars