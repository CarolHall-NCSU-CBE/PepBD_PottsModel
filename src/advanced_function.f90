module advanced_function

use datatypes
use sys_vars
use math
use utilities
use input
use database
use energy_calculation

contains

subroutine replace_extra_tryptophans(group_pep, para_pep, backbone_backup)
implicit none
integer                     :: non_trp_count=0, trp_locs(40), non_trp_locs(40) 
integer                     :: i, j, k, num_over, index_trp=0, flag, flag1, flag2, rotanum, feedback
type(groupdetails)          :: group_pep(pep_res)
real                        :: ran
character*4                 :: group_name_1(3), group_name_2(3), aminoacid_name_1
type(databackup)            :: backbone_backup(pep_res)
type(energyparameters)      :: para_pep(pep_res)

!calculate number of tryptophans
trp_count = 0
do i=1, pep_res
   if (group_pep(i)%gtype == "NTRP".or.group_pep(i)%gtype == "TRP".or.group_pep(i)%gtype == "CTRP") then
      trp_locs(trp_count+1) = i
      trp_count = trp_count + 1
   else
      non_trp_locs(non_trp_count+1) = i
      non_trp_count = non_trp_count + 1
   endif
enddo

num_over = trp_count - trp_limit
Tgroup_pep = group_pep
do i=1, num_over
   !pick a tryptophan randomly and replace with a valine (arbitrarily pick valine since it is small and shouldn't have overlaps and is of same type)
   call RANDOM_NUMBER(ran)
   index_trp=int(ran*trp_count-1.0e-3)+1
   if (index_trp.gt.trp_count) index_trp=trp_count

   flag=0

   call groupinfo(group_pep(trp_locs(index_trp))%gtype, group_name_1, flag1)
   call groupinfo("NVAL", group_name_2, flag2)
   aminoacid_name_1=group_name_2(flag1)
   call findrotamer(trp_locs(index_trp), Tgroup_pep, aminoacid_name_1, rotanum, aa_group_i)

   do j=1, rotanum
      call residue_replace(trp_locs(index_trp), group_pep, backbone_backup, j, aa_group_i, Tgroup_pep)
      if (j.eq.1) then
         !NOTE: not necessary to reload all energy parameters for entire system, in future should fix so only parameters for given AA are loaded
         call energy_parameter_single(Tgroup_pep, para_pep, trp_locs(index_trp)) !loads the force field paramters, the force field parameter doesnt change between rotamers so this only happens once
      endif
      call check_overlap(0, i, 0, Tgroup_pep, para_pep, feedback)

      if(feedback==1) then
         !once we find a rotamer without steric clashes, accept the transplant and move on after updating tryptophan info
         group_pep=Tgroup_pep
         non_trp_locs(non_trp_count+1) = trp_locs(index_trp)
         do k=index_trp, trp_count-1
            trp_locs(k) = trp_locs(k+1)
         enddo
         trp_count = trp_count - 1
         non_trp_count = non_trp_count + 1
         goto 10
      endif
   enddo

   !if we reach this part, then we couldn't replace the tryptophan with valine
   !for now, will have program exit since I don't expect this will happen often.
   !if this is a problem, then will add a solution 
   !open(20, file="error.txt", access="append")
      write(20,"(A)") "Unable to replace tryptophan with valine! Terminating program."
      close(20)
   !stop

   10 continue

enddo

return
end subroutine replace_extra_tryptophans

subroutine replace_excluded_AAs(group_pep, para_pep, backbone_backup)
implicit none 
integer                 :: i, j, k, l, flag1, flag2, rotanum, feedback
character*4             :: group_name_1(3), group_name_2(3), current_AA, new_AA, new_AA_options(20)
integer                 :: hyd_cat, new_hyd, num_replacements, replacement_types(20)
type(energyparameters)  :: para_pep(pep_res)
type(databackup)        :: backbone_backup(pep_res)
type(groupdetails)      :: group_pep(pep_res)

Tgroup_pep = group_pep
do i=1,num_ex_AAs
   do j=1,pep_res
      if (j.eq.1 .or. j.eq.pep_res) current_AA = group_pep(i)%gtype(2:4)

      if (trim(current_AA).eq.trim(excluded_AAs(i))) then
        !find type of amino acid that is in exclusion list
        call find_hydration_category(current_AA, hyd_cat)
        call groupinfo(group_pep(j)%gtype, group_name_1, flag1)

        !find possible replacements for amino acids
        replacement_types = 0
        !first type of replacement is amino acids of same type
        num_replacements = num_AAs_available(hyd_cat)
        new_AA_options(1:num_replacements) = AAs_available(1:num_replacements,hyd_cat)
        replacement_types(1:num_replacements) = hyd_cat
        !now go through all other hydration categories and append to list, if possible
        do k=1,6
            if (k.ne.hyd_cat.and.above_hyd(k).lt.0) then
                new_AA_options(num_replacements+1: num_replacements+num_AAs_available(k)) = AAs_available(1:num_AAs_available(k),k)
                replacement_types(num_replacements+1: num_replacements+num_AAs_available(k)) = k
                num_replacements = num_replacements + num_AAs_available(k)
            endif
        enddo

        !now go through replacement AA options one at a time and try replace amino acid. Stop once a good replacement is found
        do k=1, num_replacements
            new_AA = new_AA_options(k)
            call groupinfo(new_AA, group_name_2, flag2)
            new_AA = group_name_2(flag1)
            call findrotamer(j, Tgroup_pep, new_AA, rotanum, aa_group_i)

            do l=1, rotanum
                call residue_replace(j, group_pep, backbone_backup, l, aa_group_i, Tgroup_pep)
                if (l.eq.1) then
                    call energy_parameter_single(Tgroup_pep, para_pep, j) !loads the force field parameters for new residue
                endif
                call check_overlap(0, j, 0, Tgroup_pep, para_pep, feedback)
                if(feedback==1) then
                    !once we find a rotamer without steric clashes, accept the transplant and move on after updating tryptophan info
                    group_pep=Tgroup_pep
                    cur_hyd(hyd_cat) = cur_hyd(hyd_cat) - 1
                    above_hyd(hyd_cat) = above_hyd(hyd_cat) - 1
                    below_hyd(hyd_cat) = below_hyd(hyd_cat) + 1

                    new_hyd = replacement_types(k)
                    cur_hyd(new_hyd) = cur_hyd(new_hyd) + 1
                    above_hyd(new_hyd) = above_hyd(new_hyd) + 1
                    below_hyd(new_hyd) = below_hyd(new_hyd) - 1
                    goto 11
                endif
            enddo
        enddo
        !If we get here, then none of the replacements seem to work. Dang.
        !open(20, file="error.txt", access="append")
        write(20,"(A)") "Unable to replace tryptophan with valine! Terminating program."
        !close(20)
        stop
    endif

11  continue      

   enddo
enddo

return
end subroutine replace_excluded_AAs

subroutine fit_to_hydration_range(group_pep, para_pep, backbone_backup, sub_cycle)
implicit none
integer                 :: sub_cycle, i, j, k, l, flag, flag1, flag2 
integer                 :: positive_number, negative_number, aminoacid_number, rotanum, feedback
integer                 :: icpointer(pep_res, 6), pos_type(pep_res), neg_type(pep_res)
integer                 :: hyd_incr_types(20), hyd_decr_types(20), hyd_incr, hyd_decr
character*4             :: aminoacid_name(20), aminoacid_name_1, group_name_1(3), group_name_2(3)
type(groupdetails)      :: group_pep(pep_res), res_backup
type(databackup)        :: backbone_backup(pep_res)
integer                 :: pos_neg_flag !1 if pos > neg, 0 if pos < neg
integer                 :: replace_incr_count, replace_decr_count, replace_locs(pep_res)
type(energyparameters)  :: para_pep(pep_res), para_backup

cur_hyd=0
icpointer=0
!determine number of amino acids in the six categories used by PepBD. 
!Have different checks for different pH ranges
if(ph_value.le.3.9) then
   do i=1, pep_res
      if(group_pep(i)%gtype=="GLY".or.group_pep(i)%gtype=="NGLY".or.group_pep(i)%gtype=="CGLY") then
         cur_hyd(1)=cur_hyd(1)+1
         icpointer(cur_hyd(1),1)=i
      elseif(group_pep(i)%gtype=="LEU".or.group_pep(i)%gtype=="VAL".or.group_pep(i)%gtype=="ILE".or. &
         group_pep(i)%gtype=="MET".or.group_pep(i)%gtype=="PHE".or.group_pep(i)%gtype=="TYR".or. &
         group_pep(i)%gtype=="TRP".or.group_pep(i)%gtype=="NLEU".or.group_pep(i)%gtype=="NVAL".or. &
         group_pep(i)%gtype=="NILE".or.group_pep(i)%gtype=="NMET".or.group_pep(i)%gtype=="NPHE".or. &
         group_pep(i)%gtype=="NTYR".or.group_pep(i)%gtype=="NTRP".or.group_pep(i)%gtype=="CLEU" &
         .or.group_pep(i)%gtype=="CVAL".or.group_pep(i)%gtype=="CILE".or.group_pep(i)%gtype=="CMET" &
         .or.group_pep(i)%gtype=="CPHE".or.group_pep(i)%gtype=="CTYR".or.group_pep(i)%gtype=="CTRP") then
         cur_hyd(2)=cur_hyd(2)+1
         icpointer(cur_hyd(2),2)=i
      elseif(group_pep(i)%gtype=="ARG".or.group_pep(i)%gtype=="LYS".or.group_pep(i)%gtype=="HIP".or. &
         group_pep(i)%gtype=="NARG".or.group_pep(i)%gtype=="NLYS".or.group_pep(i)%gtype=="NHIP".or. & 
         group_pep(i)%gtype=="CARG".or.group_pep(i)%gtype=="CLYS".or.group_pep(i)%gtype=="CHIP") then
         cur_hyd(3)=cur_hyd(3)+1
         icpointer(cur_hyd(3),3)=i
      elseif(group_pep(i)%gtype=="SER".or.group_pep(i)%gtype=="THR".or.group_pep(i)%gtype=="ASN".or. & 
         group_pep(i)%gtype=="GLN".or.group_pep(i)%gtype=="GLH".or.group_pep(i)%gtype=="ASH".or. &
         group_pep(i)%gtype=="NSER".or.group_pep(i)%gtype=="NTHR".or.group_pep(i)%gtype=="NASN".or. &
         group_pep(i)%gtype=="NGLN".or.group_pep(i)%gtype=="NGLH".or.group_pep(i)%gtype=="NASH".or. &
         group_pep(i)%gtype=="CSER".or.group_pep(i)%gtype=="CTHR".or.group_pep(i)%gtype=="CASN" &
         .or.group_pep(i)%gtype=="CGLN".or.group_pep(i)%gtype=="CGLH".or.group_pep(i)%gtype=="CASH") then
         cur_hyd(4)=cur_hyd(4)+1
         icpointer(cur_hyd(4),4)=i
      elseif(group_pep(i)%gtype=="PRO".or.group_pep(i)%gtype=="CYS".or.group_pep(i)%gtype=="ALA".or. &
         group_pep(i)%gtype=="NPRO".or.group_pep(i)%gtype=="NCYS".or.group_pep(i)%gtype=="NALA" &
         .or.group_pep(i)%gtype=="CPRO".or.group_pep(i)%gtype=="CCYS".or.group_pep(i)%gtype=="CALA") then
         cur_hyd(5)=cur_hyd(5)+1
         icpointer(cur_hyd(5),5)=i
      endif
   enddo

elseif(ph_value.gt.3.9.and.ph_value.le.4.3) then
   do i=1, pep_res
      if(group_pep(i)%gtype=="GLY".or.group_pep(i)%gtype=="NGLY".or.group_pep(i)%gtype=="CGLY") then
         cur_hyd(1)=cur_hyd(1)+1
         icpointer(cur_hyd(1),1)=i
      elseif(group_pep(i)%gtype=="LEU".or.group_pep(i)%gtype=="VAL".or.group_pep(i)%gtype=="ILE".or. & 
         group_pep(i)%gtype=="MET".or.group_pep(i)%gtype=="PHE".or.group_pep(i)%gtype=="TYR".or. & 
         group_pep(i)%gtype=="TRP".or.group_pep(i)%gtype=="NLEU".or.group_pep(i)%gtype=="NVAL".or. &
         group_pep(i)%gtype=="NILE".or.group_pep(i)%gtype=="NMET".or.group_pep(i)%gtype=="NPHE".or. & 
         group_pep(i)%gtype=="NTYR".or.group_pep(i)%gtype=="NTRP".or.group_pep(i)%gtype=="CLEU" &
         .or.group_pep(i)%gtype=="CVAL".or.group_pep(i)%gtype=="CILE".or.group_pep(i)%gtype=="CMET" &
         .or.group_pep(i)%gtype=="CPHE".or.group_pep(i)%gtype=="CTYR".or.group_pep(i)%gtype=="CTRP") then
         cur_hyd(2)=cur_hyd(2)+1
         icpointer(cur_hyd(2),2)=i
      elseif(group_pep(i)%gtype=="ARG".or.group_pep(i)%gtype=="LYS".or.group_pep(i)%gtype=="HIP".or. &
         group_pep(i)%gtype=="NARG".or.group_pep(i)%gtype=="NLYS".or.group_pep(i)%gtype=="NHIP".or. & 
         group_pep(i)%gtype=="CARG".or.group_pep(i)%gtype=="CLYS".or.group_pep(i)%gtype=="CHIP") then
         cur_hyd(3)=cur_hyd(3)+1
         icpointer(cur_hyd(3),3)=i
      elseif(group_pep(i)%gtype=="SER".or.group_pep(i)%gtype=="THR".or.group_pep(i)%gtype=="ASN".or. &
         group_pep(i)%gtype=="GLN".or.group_pep(i)%gtype=="GLH".or.group_pep(i)%gtype=="NSER" &
         .or.group_pep(i)%gtype=="NTHR".or.group_pep(i)%gtype=="NASN".or.group_pep(i)%gtype=="NGLN" & 
         .or.group_pep(i)%gtype=="NGLH".or.group_pep(i)%gtype=="CSER".or.group_pep(i)%gtype=="CTHR".or. &
         group_pep(i)%gtype=="CASN".or.group_pep(i)%gtype=="CGLN".or.group_pep(i)%gtype=="CGLH") then
         cur_hyd(4)=cur_hyd(4)+1
         icpointer(cur_hyd(4),4)=i
      elseif(group_pep(i)%gtype=="PRO".or.group_pep(i)%gtype=="CYS".or.group_pep(i)%gtype=="ALA".or. &
         group_pep(i)%gtype=="NPRO".or.group_pep(i)%gtype=="NCYS".or.group_pep(i)%gtype=="NALA" &
         .or.group_pep(i)%gtype=="CPRO".or.group_pep(i)%gtype=="CCYS".or.group_pep(i)%gtype=="CALA") then
         cur_hyd(5)=cur_hyd(5)+1
         icpointer(cur_hyd(5),5)=i
      elseif(group_pep(i)%gtype=="ASP".or.group_pep(i)%gtype=="NASP".or.group_pep(i)%gtype=="CASP") then
         cur_hyd(6)=cur_hyd(6)+1
         icpointer(cur_hyd(6),6)=i
      endif
   enddo

elseif(ph_value.gt.4.3.and.ph_value.le.6.0) then
   do i=1, pep_res
      if(group_pep(i)%gtype=="GLY".or.group_pep(i)%gtype=="NGLY".or.group_pep(i)%gtype=="CGLY") then
         cur_hyd(1)=cur_hyd(1)+1
         icpointer(cur_hyd(1),1)=i
      elseif(group_pep(i)%gtype=="LEU".or.group_pep(i)%gtype=="VAL".or.group_pep(i)%gtype=="ILE".or. &
         group_pep(i)%gtype=="MET".or.group_pep(i)%gtype=="PHE".or.group_pep(i)%gtype=="TYR".or. &
         group_pep(i)%gtype=="TRP".or.group_pep(i)%gtype=="NLEU".or.group_pep(i)%gtype=="NVAL".or. &
         group_pep(i)%gtype=="NILE".or.group_pep(i)%gtype=="NMET".or.group_pep(i)%gtype=="NPHE".or. &
         group_pep(i)%gtype=="NTYR".or.group_pep(i)%gtype=="NTRP".or.group_pep(i)%gtype=="CLEU" &
         .or.group_pep(i)%gtype=="CVAL".or.group_pep(i)%gtype=="CILE".or.group_pep(i)%gtype=="CMET" &
         .or.group_pep(i)%gtype=="CPHE".or.group_pep(i)%gtype=="CTYR".or.group_pep(i)%gtype=="CTRP") then
         cur_hyd(2)=cur_hyd(2)+1
         icpointer(cur_hyd(2),2)=i
      elseif(group_pep(i)%gtype=="ARG".or.group_pep(i)%gtype=="LYS".or.group_pep(i)%gtype=="HIP".or. &
         group_pep(i)%gtype=="NARG".or.group_pep(i)%gtype=="NLYS".or.group_pep(i)%gtype=="NHIP".or. & 
         group_pep(i)%gtype=="CARG".or.group_pep(i)%gtype=="CLYS".or.group_pep(i)%gtype=="CHIP") then
         cur_hyd(3)=cur_hyd(3)+1
         icpointer(cur_hyd(3),3)=i
      elseif(group_pep(i)%gtype=="SER".or.group_pep(i)%gtype=="THR".or.group_pep(i)%gtype=="ASN".or. &
         group_pep(i)%gtype=="GLN".or.group_pep(i)%gtype=="NSER".or.group_pep(i)%gtype=="NTHR".or. &
         group_pep(i)%gtype=="NASN".or.group_pep(i)%gtype=="NGLN".or.group_pep(i)%gtype=="CSER".or. & 
         group_pep(i)%gtype=="CTHR".or.group_pep(i)%gtype=="CASN".or.group_pep(i)%gtype=="CGLN") then
         cur_hyd(4)=cur_hyd(4)+1
         icpointer(cur_hyd(4),4)=i
      elseif(group_pep(i)%gtype=="PRO".or.group_pep(i)%gtype=="CYS".or.group_pep(i)%gtype=="ALA".or. & 
         group_pep(i)%gtype=="NPRO".or.group_pep(i)%gtype=="NCYS".or.group_pep(i)%gtype=="NALA" &
         .or.group_pep(i)%gtype=="CPRO".or.group_pep(i)%gtype=="CCYS".or.group_pep(i)%gtype=="CALA") then
         cur_hyd(5)=cur_hyd(5)+1
         icpointer(cur_hyd(5),5)=i
      elseif(group_pep(i)%gtype=="GLU".or.group_pep(i)%gtype=="ASP".or.group_pep(i)%gtype=="NGLU".or. & 
         group_pep(i)%gtype=="NASP".or.group_pep(i)%gtype=="CGLU".or.group_pep(i)%gtype=="CASP") then
         cur_hyd(6)=cur_hyd(6)+1
         icpointer(cur_hyd(6),6)=i
      endif
   enddo

elseif(ph_value.gt.6.0.and.ph_value.lt.8.3) then
   do i=1, pep_res
      if(group_pep(i)%gtype=="GLY".or.group_pep(i)%gtype=="NGLY".or.group_pep(i)%gtype=="CGLY") then
         cur_hyd(1)=cur_hyd(1)+1
         icpointer(cur_hyd(1),1)=i
      elseif(group_pep(i)%gtype=="LEU".or.group_pep(i)%gtype=="VAL".or.group_pep(i)%gtype=="ILE" & 
         .or.group_pep(i)%gtype=="MET".or.group_pep(i)%gtype=="PHE".or.group_pep(i)%gtype=="TYR" & 
         .or.group_pep(i)%gtype=="TRP".or.group_pep(i)%gtype=="NLEU".or.group_pep(i)%gtype=="NVAL" & 
         .or.group_pep(i)%gtype=="NILE".or.group_pep(i)%gtype=="NMET".or.group_pep(i)%gtype=="NPHE" & 
         .or.group_pep(i)%gtype=="NTYR".or.group_pep(i)%gtype=="NTRP".or.group_pep(i)%gtype=="CLEU" &
         .or.group_pep(i)%gtype=="CVAL".or.group_pep(i)%gtype=="CILE".or.group_pep(i)%gtype=="CMET" & 
         .or.group_pep(i)%gtype=="CPHE".or.group_pep(i)%gtype=="CTYR".or.group_pep(i)%gtype=="CTRP") then
         cur_hyd(2)=cur_hyd(2)+1
         icpointer(cur_hyd(2),2)=i
      elseif(group_pep(i)%gtype=="ARG".or.group_pep(i)%gtype=="LYS".or.group_pep(i)%gtype=="NARG".or. & 
         group_pep(i)%gtype=="NLYS".or.group_pep(i)%gtype=="CARG".or.group_pep(i)%gtype=="CLYS") then
         cur_hyd(3)=cur_hyd(3)+1
         icpointer(cur_hyd(3),3)=i
      elseif(group_pep(i)%gtype=="SER".or.group_pep(i)%gtype=="THR".or.group_pep(i)%gtype=="ASN".or. & 
         group_pep(i)%gtype=="GLN".or.group_pep(i)%gtype=="HIE".or.group_pep(i)%gtype=="NSER".or. & 
         group_pep(i)%gtype=="NTHR".or.group_pep(i)%gtype=="NASN".or.group_pep(i)%gtype=="NGLN".or. & 
         group_pep(i)%gtype=="NHIE".or.group_pep(i)%gtype=="CSER".or.group_pep(i)%gtype=="CTHR".or. & 
         group_pep(i)%gtype=="CASN".or.group_pep(i)%gtype=="CGLN".or.group_pep(i)%gtype=="CHIE") then
         cur_hyd(4)=cur_hyd(4)+1
         icpointer(cur_hyd(4),4)=i
      elseif(group_pep(i)%gtype=="PRO".or.group_pep(i)%gtype=="CYS".or.group_pep(i)%gtype=="ALA".or. & 
         group_pep(i)%gtype=="NPRO".or.group_pep(i)%gtype=="NCYS".or.group_pep(i)%gtype=="NALA" &
         .or.group_pep(i)%gtype=="CPRO".or.group_pep(i)%gtype=="CCYS".or.group_pep(i)%gtype=="CALA") then
         cur_hyd(5)=cur_hyd(5)+1
         icpointer(cur_hyd(5),5)=i
      elseif(group_pep(i)%gtype=="GLU".or.group_pep(i)%gtype=="ASP".or.group_pep(i)%gtype=="NGLU".or. & 
         group_pep(i)%gtype=="NASP".or.group_pep(i)%gtype=="CGLU".or.group_pep(i)%gtype=="CASP") then
         cur_hyd(6)=cur_hyd(6)+1
         icpointer(cur_hyd(6),6)=i
      endif
   enddo

elseif(ph_value.ge.8.3.and.ph_value.lt.10.1) then
   do i=1, pep_res
      if(group_pep(i)%gtype=="GLY".or.group_pep(i)%gtype=="NGLY".or.group_pep(i)%gtype=="CGLY") then
         cur_hyd(1)=cur_hyd(1)+1
         icpointer(cur_hyd(1),1)=i
      elseif(group_pep(i)%gtype=="LEU".or.group_pep(i)%gtype=="VAL".or.group_pep(i)%gtype=="ILE".or. & 
         group_pep(i)%gtype=="MET".or.group_pep(i)%gtype=="PHE".or.group_pep(i)%gtype=="TYR".or. & 
         group_pep(i)%gtype=="TRP".or.group_pep(i)%gtype=="NLEU".or.group_pep(i)%gtype=="NVAL".or. & 
         group_pep(i)%gtype=="NILE".or.group_pep(i)%gtype=="NMET".or.group_pep(i)%gtype=="NPHE".or. & 
         group_pep(i)%gtype=="NTYR".or.group_pep(i)%gtype=="NTRP".or.group_pep(i)%gtype=="CLEU" &
         .or.group_pep(i)%gtype=="CVAL".or.group_pep(i)%gtype=="CILE".or.group_pep(i)%gtype=="CMET" & 
         .or.group_pep(i)%gtype=="CPHE".or.group_pep(i)%gtype=="CTYR" &
         .or.group_pep(i)%gtype=="CTRP") then
         cur_hyd(2)=cur_hyd(2)+1
         icpointer(cur_hyd(2),2)=i
      elseif(group_pep(i)%gtype=="ARG".or.group_pep(i)%gtype=="LYS".or.group_pep(i)%gtype=="NARG".or. & 
         group_pep(i)%gtype=="NLYS".or.group_pep(i)%gtype=="CARG".or.group_pep(i)%gtype=="CLYS") then
         cur_hyd(3)=cur_hyd(3)+1
         icpointer(cur_hyd(3),3)=i
      elseif(group_pep(i)%gtype=="SER".or.group_pep(i)%gtype=="THR".or.group_pep(i)%gtype=="ASN".or. & 
         group_pep(i)%gtype=="GLN".or.group_pep(i)%gtype=="HIE".or.group_pep(i)%gtype=="NSER".or. & 
         group_pep(i)%gtype=="NTHR".or.group_pep(i)%gtype=="NASN".or.group_pep(i)%gtype=="NGLN".or. & 
         group_pep(i)%gtype=="NHIE".or.group_pep(i)%gtype=="CSER".or.group_pep(i)%gtype=="CTHR".or. & 
         group_pep(i)%gtype=="CASN".or.group_pep(i)%gtype=="CGLN".or.group_pep(i)%gtype=="CHIE") then
         cur_hyd(4)=cur_hyd(4)+1
         icpointer(cur_hyd(4),4)=i
      elseif(group_pep(i)%gtype=="PRO".or.group_pep(i)%gtype=="ALA".or.group_pep(i)%gtype=="NPRO".or. & 
         group_pep(i)%gtype=="NALA".or.group_pep(i)%gtype=="CPRO".or.group_pep(i)%gtype=="CALA") then
         cur_hyd(5)=cur_hyd(5)+1
         icpointer(cur_hyd(5),5)=i
      elseif(group_pep(i)%gtype=="GLU".or.group_pep(i)%gtype=="ASP".or.group_pep(i)%gtype=="CYT".or. & 
         group_pep(i)%gtype=="NGLU".or.group_pep(i)%gtype=="NASP".or.group_pep(i)%gtype=="NCYT" &
         .or.group_pep(i)%gtype=="CGLU".or.group_pep(i)%gtype=="CASP".or.group_pep(i)%gtype=="CCYT") then
         cur_hyd(6)=cur_hyd(6)+1
         icpointer(cur_hyd(6),6)=i
      endif
   enddo

elseif(ph_value.ge.10.1.and.ph_value.lt.10.5) then
   do i=1, pep_res
      if(group_pep(i)%gtype=="GLY".or.group_pep(i)%gtype=="NGLY".or.group_pep(i)%gtype=="CGLY") then
         cur_hyd(1)=cur_hyd(1)+1
         icpointer(cur_hyd(1),1)=i
      elseif(group_pep(i)%gtype=="LEU".or.group_pep(i)%gtype=="VAL".or.group_pep(i)%gtype=="ILE".or. & 
         group_pep(i)%gtype=="MET".or.group_pep(i)%gtype=="PHE".or.group_pep(i)%gtype=="TRP".or. & 
         group_pep(i)%gtype=="NLEU".or.group_pep(i)%gtype=="NVAL".or.group_pep(i)%gtype=="NILE".or. & 
         group_pep(i)%gtype=="NMET".or.group_pep(i)%gtype=="NPHE".or.group_pep(i)%gtype=="NTRP".or. & 
         group_pep(i)%gtype=="CLEU".or.group_pep(i)%gtype=="CVAL".or.group_pep(i)%gtype=="CILE" &
         .or.group_pep(i)%gtype=="CMET".or.group_pep(i)%gtype=="CPHE".or.group_pep(i)%gtype=="CTRP") then
         cur_hyd(2)=cur_hyd(2)+1
         icpointer(cur_hyd(2),2)=i
      elseif(group_pep(i)%gtype=="ARG".or.group_pep(i)%gtype=="LYS".or.group_pep(i)%gtype=="NARG".or. & 
         group_pep(i)%gtype=="NLYS".or.group_pep(i)%gtype=="CARG".or.group_pep(i)%gtype=="CLYS") then
         cur_hyd(3)=cur_hyd(3)+1
         icpointer(cur_hyd(3),3)=i
      elseif(group_pep(i)%gtype=="SER".or.group_pep(i)%gtype=="THR".or.group_pep(i)%gtype=="ASN".or. & 
         group_pep(i)%gtype=="GLN".or.group_pep(i)%gtype=="HIE".or.group_pep(i)%gtype=="NSER".or. & 
         group_pep(i)%gtype=="NTHR".or.group_pep(i)%gtype=="NASN".or.group_pep(i)%gtype=="NGLN".or. & 
         group_pep(i)%gtype=="NHIE".or.group_pep(i)%gtype=="CSER".or.group_pep(i)%gtype=="CTHR".or. & 
         group_pep(i)%gtype=="CASN".or.group_pep(i)%gtype=="CGLN".or.group_pep(i)%gtype=="CHIE") then
         cur_hyd(4)=cur_hyd(4)+1
         icpointer(cur_hyd(4),4)=i
      elseif(group_pep(i)%gtype=="PRO".or.group_pep(i)%gtype=="ALA".or.group_pep(i)%gtype=="NPRO".or. & 
         group_pep(i)%gtype=="NALA".or.group_pep(i)%gtype=="CPRO".or.group_pep(i)%gtype=="CALA") then
         cur_hyd(5)=cur_hyd(5)+1
         icpointer(cur_hyd(5),5)=i
      elseif(group_pep(i)%gtype=="GLU".or.group_pep(i)%gtype=="ASP".or.group_pep(i)%gtype=="CYT".or. & 
         group_pep(i)%gtype=="TYX".or.group_pep(i)%gtype=="NGLU".or.group_pep(i)%gtype=="NASP"  &
         .or.group_pep(i)%gtype=="NCYT".or.group_pep(i)%gtype=="NTYX".or.group_pep(i)%gtype=="CGLU".or. & 
         group_pep(i)%gtype=="CASP".or.group_pep(i)%gtype=="CCYT".or.group_pep(i)%gtype=="CTYX") then
         cur_hyd(6)=cur_hyd(6)+1
         icpointer(cur_hyd(6),6)=i
      endif
   enddo

elseif(ph_value.ge.10.5.and.ph_value.lt.12.5) then
   do i=1, pep_res
      if(group_pep(i)%gtype=="GLY".or.group_pep(i)%gtype=="NGLY".or.group_pep(i)%gtype=="CGLY") then
         cur_hyd(1)=cur_hyd(1)+1
         icpointer(cur_hyd(1),1)=i
      elseif(group_pep(i)%gtype=="LEU".or.group_pep(i)%gtype=="VAL".or.group_pep(i)%gtype=="ILE".or. & 
         group_pep(i)%gtype=="MET".or.group_pep(i)%gtype=="PHE".or.group_pep(i)%gtype=="TRP".or. & 
         group_pep(i)%gtype=="NLEU".or.group_pep(i)%gtype=="NVAL".or.group_pep(i)%gtype=="NILE".or. & 
         group_pep(i)%gtype=="NMET".or.group_pep(i)%gtype=="NPHE".or.group_pep(i)%gtype=="NTRP".or. & 
         group_pep(i)%gtype=="CLEU".or.group_pep(i)%gtype=="CVAL".or.group_pep(i)%gtype=="CILE" &
         .or.group_pep(i)%gtype=="CMET".or.group_pep(i)%gtype=="CPHE".or.group_pep(i)%gtype=="CTRP") then
         cur_hyd(2)=cur_hyd(2)+1
         icpointer(cur_hyd(2),2)=i
      elseif(group_pep(i)%gtype=="ARG".or.group_pep(i)%gtype=="NARG".or.group_pep(i)%gtype=="CARG") then
         cur_hyd(3)=cur_hyd(3)+1
         icpointer(cur_hyd(3),3)=i
      elseif(group_pep(i)%gtype=="SER".or.group_pep(i)%gtype=="THR".or.group_pep(i)%gtype=="ASN".or. & 
         group_pep(i)%gtype=="GLN".or.group_pep(i)%gtype=="HIE".or.group_pep(i)%gtype=="LYN" &
         .or.group_pep(i)%gtype=="NSER".or.group_pep(i)%gtype=="NTHR".or.group_pep(i)%gtype=="NASN".or. & 
         group_pep(i)%gtype=="NGLN".or.group_pep(i)%gtype=="NHIE".or.group_pep(i)%gtype=="NLYN" &
         .or.group_pep(i)%gtype=="CSER".or.group_pep(i)%gtype=="CTHR".or.group_pep(i)%gtype=="CASN".or. & 
         group_pep(i)%gtype=="CGLN".or.group_pep(i)%gtype=="CHIE".or.group_pep(i)%gtype=="CLYN") then
         cur_hyd(4)=cur_hyd(4)+1
         icpointer(cur_hyd(4),4)=i
      elseif(group_pep(i)%gtype=="PRO".or.group_pep(i)%gtype=="ALA".or.group_pep(i)%gtype=="NPRO".or. & 
         group_pep(i)%gtype=="NALA".or.group_pep(i)%gtype=="CPRO".or.group_pep(i)%gtype=="CALA") then
         cur_hyd(5)=cur_hyd(5)+1
         icpointer(cur_hyd(5),5)=i
      elseif(group_pep(i)%gtype=="GLU".or.group_pep(i)%gtype=="ASP".or.group_pep(i)%gtype=="CYT".or. & 
         group_pep(i)%gtype=="TYX".or.group_pep(i)%gtype=="NGLU".or.group_pep(i)%gtype=="NASP"  &
         .or.group_pep(i)%gtype=="NCYT".or.group_pep(i)%gtype=="NTYX".or.group_pep(i)%gtype=="CGLU".or. & 
         group_pep(i)%gtype=="CASP".or.group_pep(i)%gtype=="CCYT".or.group_pep(i)%gtype=="CTYX") then
         cur_hyd(6)=cur_hyd(6)+1
         icpointer(cur_hyd(6),6)=i
      endif
   enddo

elseif(ph_value.ge.12.5) then
   do i=1, pep_res
      if(group_pep(i)%gtype=="GLY".or.group_pep(i)%gtype=="NGLY".or.group_pep(i)%gtype=="CGLY") then
         cur_hyd(1)=cur_hyd(1)+1
         icpointer(cur_hyd(1),1)=i
      elseif(group_pep(i)%gtype=="LEU".or.group_pep(i)%gtype=="VAL".or.group_pep(i)%gtype=="ILE".or. & 
         group_pep(i)%gtype=="MET".or.group_pep(i)%gtype=="PHE".or.group_pep(i)%gtype=="TRP".or. & 
         group_pep(i)%gtype=="NLEU".or.group_pep(i)%gtype=="NVAL".or.group_pep(i)%gtype=="NILE".or. & 
         group_pep(i)%gtype=="NMET".or.group_pep(i)%gtype=="NPHE".or.group_pep(i)%gtype=="NTRP".or. & 
         group_pep(i)%gtype=="CLEU".or.group_pep(i)%gtype=="CVAL".or.group_pep(i)%gtype=="CILE" &
         .or.group_pep(i)%gtype=="CMET".or.group_pep(i)%gtype=="CPHE".or.group_pep(i)%gtype=="CTRP") then
         cur_hyd(2)=cur_hyd(2)+1
         icpointer(cur_hyd(2),2)=i
      elseif(group_pep(i)%gtype=="SER".or.group_pep(i)%gtype=="THR".or.group_pep(i)%gtype=="ASN".or. & 
         group_pep(i)%gtype=="GLN".or.group_pep(i)%gtype=="HIE".or.group_pep(i)%gtype=="LYN".or. & 
         group_pep(i)%gtype=="ARN".or.group_pep(i)%gtype=="NSER".or.group_pep(i)%gtype=="NTHR".or. & 
         group_pep(i)%gtype=="NASN".or.group_pep(i)%gtype=="NGLN".or.group_pep(i)%gtype=="NHIE".or. & 
         group_pep(i)%gtype=="NLYN".or.group_pep(i)%gtype=="NARN".or.group_pep(i)%gtype=="CSER".or. & 
         group_pep(i)%gtype=="CTHR".or.group_pep(i)%gtype=="CASN".or.group_pep(i)%gtype=="CGLN".or. & 
         group_pep(i)%gtype=="CHIE".or.group_pep(i)%gtype=="CLYN".or.group_pep(i)%gtype=="CARN") then
         cur_hyd(4)=cur_hyd(4)+1
         icpointer(cur_hyd(4),4)=i
      elseif(group_pep(i)%gtype=="PRO".or.group_pep(i)%gtype=="ALA".or.group_pep(i)%gtype=="NPRO".or. & 
         group_pep(i)%gtype=="NALA".or.group_pep(i)%gtype=="CPRO".or.group_pep(i)%gtype=="CALA") then
         cur_hyd(5)=cur_hyd(5)+1
         icpointer(cur_hyd(5),5)=i
      elseif(group_pep(i)%gtype=="GLU".or.group_pep(i)%gtype=="ASP".or.group_pep(i)%gtype=="CYT".or. & 
         group_pep(i)%gtype=="TYX".or.group_pep(i)%gtype=="NGLU".or.group_pep(i)%gtype=="NASP"  &
         .or.group_pep(i)%gtype=="NCYT".or.group_pep(i)%gtype=="NTYX".or.group_pep(i)%gtype=="CGLU".or. & 
         group_pep(i)%gtype=="CASP".or.group_pep(i)%gtype=="CCYT".or.group_pep(i)%gtype=="CTYX") then
         cur_hyd(6)=cur_hyd(6)+1
         icpointer(cur_hyd(6),6)=i
      endif
   enddo
endif

!for each hydration type, calculate how far current number is from minimum and maximum limit
!A positive number means hydration category is outside bound (either min or max)
!A negative number means the hydration category is within the bounds
below_hyd = min_hyd - cur_hyd
above_hyd = cur_hyd - max_hyd
positive_number=0
negative_number=0

!randomize order of amino acids in each group type (not efficient, but it works)
!do i=1, 6
!   if(above_hyd(i)>0) then
!
!      do j=1, (cur_hyd(i)-1)
!         call RANDOM_NUMBER(ran2)
!         ic1=int(ran2*cur_hyd(i))+1
!         
!         call RANDOM_NUMBER(ran2)
!         ic2=int(ran2*cur_hyd(i))+1
!
!         k=icpointer(ic1,i)
!         icpointer(ic1,i)=icpointer(i,ic2)
!         icpointer(ic2,i)=k
!      enddo
!      pos_type(positive_number+1:positive_number+above_hyd(i)) = i
!      positive_number=positive_number+above_hyd(i)
!
!   elseif(below_hyd(i)>0) then
!      neg_type(negative_number+1:negative_number+below_hyd(i)) = i
!      negative_number=negative_number+below_hyd(i)
!   endif
!enddo
!
!!randomize order of amino acids in pos_type (not efficient, but it works)
!do i=1, (positive_number-1)
!   call RANDOM_NUMBER(ran2)
!   ic1=int(ran2*positive_number)+1
!
!   call RANDOM_NUMBER(ran2)
!   ic2=int(ran2*positive_number)+1
!
!   k=pos_type(ic1)
!   pos_type(ic1)=pos_type(ic2)
!   pos_type(ic2)=k
!enddo

!randomize order of amino acids in pos_type (not efficient, but it works)
!do i=1, (negative_number-1)
!
!   call RANDOM_NUMBER(ran2)
!   ic1=int(ran2*negative_number)+1
!
!   call RANDOM_NUMBER(ran2)
!   ic2=int(ran2*negative_number)+1
!
!   k=neg_type(ic1)
!   neg_type(ic1)=neg_type(ic2)
!   neg_type(ic2)=k
!enddo

if (positive_number >= negative_number) then
   sub_cycle = positive_number
   pos_neg_flag = 1
else      
   sub_cycle = negative_number
   pos_neg_flag = 0
endif

!Replace amino acids that are outside the allowable hydration range
do i=1, sub_cycle
    replace_incr_count = 0
    replace_decr_count = 0
    if (pos_neg_flag.eq.1) then !there are more amino acids over hydration limit than under
        !select residue that is to be replaced 
        positive_number = positive_number - 1
        hyd_decr_types = pos_type(i)
        replace_decr_count = cur_hyd(pos_type(i))
        replace_locs(1:replace_decr_count) = icpointer(:,pos_type(i))
        !select residue that is to be added by...  
        if (negative_number.gt.0) then  !selecting from residue type that is below hydration limit (first priority)
            !rep_neg_group = neg_type(i)
            hyd_incr_types = neg_type(i)
            replace_incr_count = num_AAs_available(neg_type(i))
            call scmf_choose_aminoacid(neg_type(i), aminoacid_number, aminoacid_name)
            negative_number = negative_number - 1
        else !selecting random type that can be increased by 1 without being over hydration limit
            do j=1,6 !iterate through all hydration types and see which ones can be increased by 1 without exceeding hydration limits
                if (above_hyd(j) < 0) then
                aminoacid_name(replace_incr_count+1:replace_incr_count+num_AAs_available(j)) = AAs_available(1:num_AAs_available(j),j)
                hyd_incr_types(replace_incr_count+1:replace_incr_count+num_AAs_available(j)) = j
                replace_incr_count = replace_incr_count + num_AAs_available(j)
                endif
            enddo
        endif
    else !there are more amino acids over hydration limit than under
        !select residue that is to be added
        negative_number = negative_number - 1
        aminoacid_name(1:num_AAs_available(i)) = AAs_available(1:num_AAs_available(i),i)
        hyd_incr_types = neg_type(i)
        replace_incr_count = num_AAs_available(i)

        !select residue that is to be replaced by...  
        if (positive_number.gt.0) then !selecting from residue type that is over hydration limit (first priority)
            hyd_decr_types = pos_type(i)
            positive_number = positive_number - 1
            replace_decr_count = cur_hyd(pos_type(i))
            replace_locs(1:replace_decr_count) = icpointer(:,pos_type(i))
        else !selecting random type that can be decreased by 1 without being under hydration limit
            replace_decr_count = 0
            do j=1,6 !iterate through all hydration types and see which ones can be decreased by 1 without exceeding hydration limits
                if (below_hyd(j) < 0) then
                hyd_decr_types(replace_decr_count+1:replace_decr_count+cur_hyd(j)) = j
                replace_locs(replace_decr_count+1:replace_decr_count+cur_hyd(j)) = icpointer(1:cur_hyd(j), j)
                replace_decr_count = replace_decr_count + cur_hyd(j)
                endif
            enddo
        endif
    endif

    flag=0
    l=1

    do while(l<=replace_decr_count) !iterate through all members of amino acid type(s) that is to be replaced, stopping once one works
        para_backup = para_pep(replace_locs(l)) !store energy parameters of amino acid being replaced
        res_backup = group_pep(replace_locs(l))
        j=1
        do while(j<=replace_incr_count) !try replacing with all possible amino acids in the replacement group type(s), stopping once one works
            call groupinfo(group_pep(replace_locs(l))%gtype, group_name_1, flag1) !get name for amino acid that is going to be replaced
            call groupinfo(aminoacid_name(j), group_name_2, flag2)
            aminoacid_name_1=group_name_2(flag1)
            call findrotamer(replace_locs(l), group_pep, aminoacid_name_1, rotanum, aa_group_i)
            do k=1, rotanum !go through all rotamers of replacement amino acid to see if one can fit without steric clash with receptor or peptide
                call residue_replace_inplace(replace_locs(l), group_pep, backbone_backup, k, aa_group_i)
                if (k.eq.1) call energy_parameter_single(group_pep, para_pep, replace_locs(l)) !loads the energy parameters for replacement amino acid
                call check_overlap(0, replace_locs(l), 0, group_pep, para_pep, feedback)
                if(feedback==1) then !replacement was successful!
                    hyd_decr = hyd_decr_types(l)
                    hyd_incr = hyd_incr_types(j)
                    flag=1
                    goto 10
                endif
            enddo
            j=j+1
        enddo
        l=l+1
        para_pep(replace_locs(l)) = para_backup !replacement not successful, so restore energy parameters of old residue.
        group_pep(replace_locs(l)) = res_backup
    enddo

    10  continue

    if(flag==1) then
        !now fix icpointers so the amino acid types are correct
        icpointer(cur_hyd(hyd_incr)+1, hyd_incr) = replace_locs(l)

        do k=1, cur_hyd(hyd_decr)
            if (icpointer(k,hyd_decr).eq.replace_locs(l)) then
                icpointer(k,hyd_decr) = icpointer(cur_hyd(hyd_decr),hyd_decr)
            endif
        enddo

        !update hydration properties
        cur_hyd(hyd_decr)=cur_hyd(hyd_decr) - 1
        cur_hyd(hyd_incr)=cur_hyd(hyd_incr) + 1 

        above_hyd(hyd_incr) = above_hyd(hyd_incr) + 1
        below_hyd(hyd_incr) = below_hyd(hyd_incr) - 1

        above_hyd(hyd_decr) = above_hyd(hyd_decr) - 1
        below_hyd(hyd_decr) = below_hyd(hyd_decr) + 1
    else
        !open(20, file="error.txt", access="append")
            write(20,"(A)") "Unable to make substitution to satisfy hydration properties without a steric clash!" 
            write(20,"(A)") "Either adjust the hydration properties to permit this amino acid, or provide a different staring structure!"
        !close(20)
        stop
    endif
enddo

return
end subroutine fit_to_hydration_range

subroutine dihedral_gradient_sidechain(group_pep, para_pep, res_i, sc_coords_orig, &
    num_sc_angles, sc_groups, monitor, index, delta_chi, GB_flag, CG_flag, CA, energy_forward, energy_backward, gradient)
implicit none
integer                                         :: num_sc_angles, sc_groups(6), monitor(6)
integer                                         :: i, res_i, GB_flag, CG_flag
real                                            :: delta_chi
real                                            :: CA(3)
double precision                                :: energy_forward(:), energy_backward(:)
double precision                                :: gradient(:,:)
type(groupdetails)                              :: group_pep(pep_res)
type(index4sidechain)                           :: index(60)
type(conformer4sidechain)                       :: sc_coords(6), sc_coords_orig(6)
type(energyparameters)                          :: para_pep(:)

do i=1, num_sc_angles !for each dihedral angle (i.e. "num_sc_angles"), will calculate energy change due to increasing dihedral by this amount
    sc_coords=sc_coords_orig !reset coordinates of sidechain to original values
    call sidechain_rotation_4_gradient(num_sc_angles, i, sc_groups, delta_chi, sc_coords, CA, group_pep, monitor, index, res_i) !perform rotation of sidechain
    res_coords_flags(res_i) = 1; GB_area_calc_flags(res_i) = 1
    call bindingenergy_single(group_pep, para_pep, res_i, GB_flag, CG_flag, energy_forward(i)) !get energy at new sidechain position
enddo

!will now repeat the entire process, but with a negative angle change
do i=1, num_sc_angles
    sc_coords=sc_coords_orig !reset coordinates of sidechain to original values
    call sidechain_rotation_4_gradient(num_sc_angles, i, sc_groups, -delta_chi, sc_coords, CA, group_pep, monitor, index, res_i) !perform rotation of sidechain
    res_coords_flags(res_i) = 1; GB_area_calc_flags(res_i) = 1
    call bindingenergy_single(group_pep, para_pep, res_i, GB_flag, CG_flag, energy_backward(i))
enddo

!calculate numerical derivative for each dihedral angle
gradient(1:num_sc_angles,1) = (energy_forward(1:num_sc_angles)-energy_backward(1:num_sc_angles))/(2*delta_chi)

end subroutine dihedral_gradient_sidechain

subroutine dihedral_gradient_whole(group_pep, para_pep, res_i, sc_coords_orig, num_sc_angles, sc_groups, &
    dihedral, monitor, index, delta_chi, GB_flag, CG_flag, CG_flags, CA, energy_forward, energy_backward, gradient)
implicit none
real                                    :: delta_chi
integer                                 :: i, num_sc_angles, monitor(6), sc_groups(6), res_i, GB_flag, CG_flag, CG_flags(pep_res)
real                                    :: CA(3)
double precision                        :: energy_forward(:), energy_backward(:)
double precision                        :: gradient(:,:)
type(groupdetails)                      :: group_pep(pep_res)
type(energyparameters)                  :: para_pep(pep_res)
type(conformer4sidechain)               :: sc_coords(6),  sc_coords_orig(6)
type(index4sidechain)                   :: index(60)
type(dihedralparameters)                :: dihedral

!first get the forward difference
do i=1, num_sc_angles
    sc_coords=sc_coords_orig
    call sidechain_rotation_4_gradient(num_sc_angles, i, sc_groups, delta_chi, sc_coords, CA, group_pep, monitor, index, res_i) !perform rotation of sidechain
    res_coords_flags(res_i) = 1; GB_area_calc_flags(res_i) = 1
    call bindingenergy_sidechainoptimization(group_pep, para_pep, res_i, dihedral, num_sc_angles, GB_flag, CG_flag, CG_flags, energy_forward(i))
enddo

!now evaluate the backward value
do i=1, num_sc_angles
    sc_coords=sc_coords_orig
    call sidechain_rotation_4_gradient(num_sc_angles, i, sc_groups, delta_chi, sc_coords, CA, group_pep, monitor, index, res_i) !perform rotation of sidechain
    res_coords_flags(res_i) = 1; GB_area_calc_flags(res_i) = 1
    call bindingenergy_sidechainoptimization(group_pep, para_pep, res_i, dihedral, num_sc_angles, GB_flag, CG_flag, CG_flags, energy_backward(i))
enddo

!numerical gradient is estimated by finite difference between forward and backward value 
gradient(1:num_sc_angles,1) = (energy_forward(1:num_sc_angles)-energy_backward(1:num_sc_angles))/(2*delta_chi)

end subroutine dihedral_gradient_whole

subroutine concertedrotation_center(group_pep, original_pep, backbone_backup, ran_resi, phipsi_num, &
            group_candidates, backbone_backup_candidates)
implicit none
integer                     :: ran_resi(3), i, j, k, count, ic
integer                     :: phipsi_num, psi_num, H, feedback, flag
real                        :: C0(3), N1(3), CA1(3), C1(3), N2(3), CA2(3), C2(3) 
real                        :: N3(3), CA3(3), C3(3), N4(4), CA3_target(3)
real                        :: C10_temp(3), N1_temp(3), CA1_temp(3), C1_temp(3), N12_temp(3)
real                        :: C20_temp(3), N2_temp(3), CA2_temp(3), C2_temp(3), N22_temp(3)
real                        :: C30_temp(3), N3_temp(3), CA3_temp(3)
real                        :: c11, c22, c33, c44, dist, rmsd, ran2
real                        :: omg1, omg2, theta_omg0, theta_phi1, theta_psi1 
real                        :: theta_omg1, theta_phi2, theta_psi2, theta_omg2
real                        :: l_omg0, l_phi1, l_psi1, l_omg1, l_phi2, l_psi2, l_omg2, l_phi3
real                        :: U_phi(3), U_psi(3), U_omg(3), U_omg0(3), U_phi1(3), U_psi1(3)
real                        :: U_omg1(3), U_phi2(3), U_psi2(3), U_omg2(3), U_phi3(3)
real                        :: Lomg1(3), Lomg2(3), Q_omg1(3), Q_omg2(3), e(3), t(3) 
real                        :: t_prime(3), q2(3), q2_prime(3), m(3), n(3)
real                        :: Tlab(3,3), Tlab_inverse(3,3), Tphi1(3,3), Tpsi1(3,3), Tomg1(3,3)
real                        :: r_CA_N(3), r1(3), r1_old(3), r1_new(3), r_old_temp(3), r_new_temp(3)
real                        :: delta_theta, theta_old_temp, theta_new_temp
real                        :: delta_phi, phi1_old, phi1, phi2, phi(2), psi(2), psi1(2)
real                        :: phi1_temp, psi1_temp, phi2_temp, psi2_temp, phi3_temp, psi3_temp
integer                     :: phi_iter, phi_lim
type(groupdetails)          :: group_pep(pep_res), original_pep(pep_res), group_candidates(pep_res, 20), Tgroup(pep_res), Tgroup_new(pep_res)
type(databackup)            :: backbone_backup(pep_res), backbone_backup_candidates(pep_res, 20) 
type(databackup)            :: Tbackbone_backup(pep_res), Tbackbone_backup_new(pep_res)

Tgroup=group_pep(1:pep_res)
Tbackbone_backup=backbone_backup

phipsi_num=0

!get coordinates of the N, CA, C and atoms in the first pivot
do i=1, Tgroup(ran_resi(1))%cnum1
   if(Tgroup(ran_resi(1))%atype1(i)=="N") then
      N1(1:3)=Tgroup(ran_resi(1))%coo1(1:3,i)
   elseif(Tgroup(ran_resi(1))%atype1(i)=="CA") then
      CA1(1:3)=Tgroup(ran_resi(1))%coo1(1:3,i)
   endif
enddo
do i=1, Tgroup(ran_resi(1))%cnum3
   if(Tgroup(ran_resi(1))%atype3(i)=="C") then
      C1(1:3)=Tgroup(ran_resi(1))%coo3(1:3,i)
   endif
enddo

U_phi1(1:3)=CA1(1:3)-N1(1:3) !vector between CA1 and N1 (axis for phi angles)
U_psi1(1:3)=C1(1:3)-CA1(1:3) !vector between C1 and CA1 (axis for psi angles)
call normalize_vector(U_phi1)
call normalize_vector(U_psi1)

theta_phi1=acosd(dot_product(U_phi1, U_psi)) !get phi1 angle
Tlab(1,1:3)=U_phi1(1:3)

Tlab(3,1)=(U_phi(2)*U_psi(3)-U_phi(3)*U_psi(2))/sind(theta_phi1)
Tlab(3,2)=(U_phi(3)*U_psi(1)-U_phi(1)*U_psi(3))/sind(theta_phi1)
Tlab(3,3)=(U_phi(1)*U_psi(2)-U_phi(2)*U_psi(1))/sind(theta_phi1)

Tlab(2,1)=Tlab(3,2)*Tlab(1,3)-Tlab(3,3)*Tlab(1,2)
Tlab(2,2)=Tlab(3,3)*Tlab(1,1)-Tlab(3,1)*Tlab(1,3)
Tlab(2,3)=Tlab(3,1)*Tlab(1,2)-Tlab(3,2)*Tlab(1,1)

Tlab_inverse=transpose(Tlab)

!get angle between 60 and 74 (where delta_theta is between -4 and 4)
count=0
do while(count.le.20)
    call RANDOM_NUMBER(ran2)
    delta_theta=ran2*8.0-4.0
    if((theta_phi1+delta_theta).ge.60.0.and.(theta_phi1+delta_theta).le.74.0) exit
    count=count+1
    if(count.eq.21) goto 20
enddo

do ic=ran_resi(1), ran_resi(3)
    if(ic.ne.ran_resi(1)) then
        do i=1, Tgroup(ic)%cnum1
            r_old_temp(1:3)=Tgroup(ic)%coo1(1:3,i)-CA1(1:3)
            r1_old=matmul(Tlab, r_old_temp)
            theta_old_temp=acosd(r1_old(1)/sqrt(r1_old(1)*r1_old(1)+r1_old(2)*r1_old(2)))
            if(r1_old(2).lt.0.0) then
                theta_old_temp=360.00-theta_old_temp
            endif
            theta_new_temp=theta_old_temp+delta_theta
            r1_new(1)=sqrt(r1_old(1)*r1_old(1)+r1_old(2)*r1_old(2))*cosd(theta_new_temp)
            r1_new(2)=sqrt(r1_old(1)*r1_old(1)+r1_old(2)*r1_old(2))*sind(theta_new_temp)
            r1_new(3)=r1_old(3)
            r_new_temp=matmul(Tlab_inverse, r1_new)
            Tgroup(ic)%coo1(1:3,i)=r_new_temp(1)+CA1(1:3)
        enddo

        do i=1, Tgroup(ic)%cnum2
            r_old_temp(1:3)=Tgroup(ic)%coo2(1:3,i)-CA1(1:3)
            r1_old=matmul(Tlab, r_old_temp)
            theta_old_temp=acosd(r1_old(1)/sqrt(r1_old(1)*r1_old(1)+r1_old(2)*r1_old(2)))
            if(r1_old(2).lt.0.0) then
                theta_old_temp=360.00-theta_old_temp
            endif
            theta_new_temp=theta_old_temp+delta_theta
            r1_new(1)=sqrt(r1_old(1)*r1_old(1)+r1_old(2)*r1_old(2))*cosd(theta_new_temp)
            r1_new(2)=sqrt(r1_old(1)*r1_old(1)+r1_old(2)*r1_old(2))*sind(theta_new_temp)
            r1_new(3)=r1_old(3)
            r_new_temp=matmul(Tlab_inverse, r1_new)
            Tgroup(ic)%coo2(1:3,i)=r_new_temp(1:3)+CA1(1:3)
        enddo

        r_old_temp(1:3)=Tbackbone_backup(ic)%coo(1:3)-CA1(1:3)
        r1_old=matmul(Tlab, r_old_temp)
        theta_old_temp=acosd(r1_old(1)/sqrt(r1_old(1)*r1_old(1)+r1_old(2)*r1_old(2)))
        if(r1_old(2).lt.0.0) then
            theta_old_temp=360.00-theta_old_temp
        endif
        theta_new_temp=theta_old_temp+delta_theta
        r1_new(1)=sqrt(r1_old(1)*r1_old(1)+r1_old(2)*r1_old(2))*cosd(theta_new_temp)
        r1_new(2)=sqrt(r1_old(1)*r1_old(1)+r1_old(2)*r1_old(2))*sind(theta_new_temp)
        r1_new(3)=r1_old(3)
        r_new_temp=matmul(Tlab_inverse, r1_new)
        Tbackbone_backup(ic)%coo(1:3)=r_new_temp(1:3)+CA1(1:3)
    endif

    if(ic.ne.ran_resi(3)) then
        do i=1, Tgroup(ic)%cnum3
            r_old_temp(1:3)=Tgroup(ic)%coo3(1:3,i)-CA1(1:3)
            r1_old=matmul(Tlab, r_old_temp)
            theta_old_temp=acosd(r1_old(1)/sqrt(r1_old(1)*r1_old(1)+r1_old(2)*r1_old(2)))
            if(r1_old(2).lt.0.0) then
                theta_old_temp=360.00-theta_old_temp
            endif
            theta_new_temp=theta_old_temp+delta_theta
            r1_new(1)=sqrt(r1_old(1)*r1_old(1)+r1_old(2)*r1_old(2))*cosd(theta_new_temp)
            r1_new(2)=sqrt(r1_old(1)*r1_old(1)+r1_old(2)*r1_old(2))*sind(theta_new_temp)
            r1_new(3)=r1_old(3)
            r_new_temp=matmul(Tlab_inverse, r1_new)
            Tgroup(ic)%coo3(1:3,i)=r_new_temp(1:3)+CA1(1:3)
        enddo
    endif
enddo

!get position of C before N in first pivot
do i=1, Tgroup(ran_resi(1)-1)%cnum3
   if(Tgroup(ran_resi(1)-1)%atype3(i)=="C") then
      C0(1:3)=Tgroup(ran_resi(1)-1)%coo3(1:3,i)
      exit
   endif
enddo

!get position of N and alpha carbon in first pivot
do i=1, Tgroup(ran_resi(1))%cnum1
   if(Tgroup(ran_resi(1))%atype1(i)=="N") then
      N1(1:3)=Tgroup(ran_resi(1))%coo1(1:3,i)
   elseif(Tgroup(ran_resi(1))%atype1(i)=="CA") then
      CA1(1:3)=Tgroup(ran_resi(1))%coo1(1:3,i)
   endif
enddo

!get position of C in first pivot
do i=1, Tgroup(ran_resi(1))%cnum3
   if(Tgroup(ran_resi(1))%atype3(i)=="C") then
      C1(1:3)=Tgroup(ran_resi(1))%coo3(1:3,i)
      exit
   endif
enddo

!get position of N and alpha carbon in second pivot
do i=1, Tgroup(ran_resi(2))%cnum1
   if(Tgroup(ran_resi(2))%atype1(i)=="N") then
      N2(1:3)=Tgroup(ran_resi(2))%coo1(1:3,i)
   elseif(Tgroup(ran_resi(2))%atype1(i)=="CA") then
      CA2(1:3)=Tgroup(ran_resi(2))%coo1(1:3,i)
   endif
enddo

!get position of C in second pivot
do i=1, Tgroup(ran_resi(2))%cnum3
   if(Tgroup(ran_resi(2))%atype3(i)=="C") then
      C2(1:3)=Tgroup(ran_resi(2))%coo3(1:3,i)
      exit
   endif
enddo

!get position of N and alpha carbon in third pivot
do i=1, Tgroup(ran_resi(3))%cnum1
   if(Tgroup(ran_resi(3))%atype1(i)=="N") then
      N3(1:3)=Tgroup(ran_resi(3))%coo1(1:3,i)
   elseif(Tgroup(ran_resi(3))%atype1(i)=="CA") then
      CA3(1:3)=Tgroup(ran_resi(3))%coo1(1:3,i)
   endif
enddo

!get position of C in third pivot
do i=1, Tgroup(ran_resi(3))%cnum3
   if(Tgroup(ran_resi(3))%atype3(i)=="C") then
      C3(1:3)=Tgroup(ran_resi(3))%coo3(1:3,i)
      exit
   endif
enddo

!get position of N in residue after third pivot
do i=1, Tgroup(ran_resi(3)+1)%cnum1
   if(Tgroup(ran_resi(3)+1)%atype1(i)=="N") then
      N4(1:3)=Tgroup(ran_resi(3)+1)%coo1(1:3,i)
      exit
   endif
enddo

!get bond vectors and bond lengths for all rotations
U_omg0(1:3)=N1(1:3)-C0(1)
l_omg0=sqrt(U_omg0(1)**2+U_omg0(2)**2+U_omg0(3)**2)

!calculate bond vectors
U_phi1 = CA1 - N1
U_psi1 = C1  - CA1
U_omg1 = N2  - C1
U_phi2 = CA2 - N2
U_psi2 = C2  - CA2
U_omg2 = N3  - C2
U_phi3 = CA3 - N3

!get bond vector lengths and normalize vectors

call vec_norm(U_phi1, l_phi1)
call vec_norm(U_psi1, l_psi1)
call vec_norm(U_omg1, l_omg1)
call vec_norm(U_phi2, l_phi2)
call vec_norm(U_psi2, l_psi2)
call vec_norm(U_omg2, l_omg2)
call vec_norm(U_phi3, l_phi3)
U_phi1 = U_phi1/l_phi1
U_psi1 = U_psi1/l_psi1
U_omg1 = U_omg1/l_omg1
U_phi2 = U_phi2/l_phi2
U_psi2 = U_psi2/l_psi2
U_omg2 = U_omg2/l_omg2
U_phi3 = U_phi2/l_phi3

!get dihedrals for phi1, omega1, and omega2
call calc_dihedral(C0, N1, CA1, C1, phi1_old)
call calc_dihedral(CA1, C1, N2, CA2, omg1)
call calc_dihedral(CA2, C2, N3, CA3, omg2)

!get all other current dihedral angles
theta_omg0=acosd(dot_product(U_omg0, U_phi1))
theta_phi1=acosd(dot_product(U_phi1, U_psi1))
theta_psi1=acosd(dot_product(U_psi1, U_omg1))
theta_omg1=acosd(dot_product(U_omg1, U_phi2))
theta_phi2=acosd(dot_product(U_phi2, U_psi2))
theta_psi2=acosd(dot_product(U_psi2, U_omg2))
theta_omg2=acosd(dot_product(U_omg2, U_phi3))

!get matrix to express coordinates in "lab frame" for first pivot
Tlab(1,1:3)=U_phi1(1:3)

Tlab(3,1)=(U_phi(2)*U_omg(3)-U_phi(3)*U_omg(2))/sind(theta_omg0)
Tlab(3,2)=(U_phi(3)*U_omg(1)-U_phi(1)*U_omg(3))/sind(theta_omg0)
Tlab(3,3)=(U_phi(1)*U_omg(2)-U_phi(2)*U_omg(1))/sind(theta_omg0)
Tlab(2,1)=Tlab(3,2)*Tlab(1,3)-Tlab(3,3)*Tlab(1,2)
Tlab(2,2)=Tlab(3,3)*Tlab(1,1)-Tlab(3,1)*Tlab(1,3)
Tlab(2,3)=Tlab(3,1)*Tlab(1,2)-Tlab(3,2)*Tlab(1,1)

!store final target location for the alpha carbon of the third pivot point
do i=1, group_pep(ran_resi(3))%cnum1
   if(group_pep(ran_resi(3))%atype1(i)=="CA") then
      CA3_target(1:3)=group_pep(ran_resi(3))%coo1(1:3,i)
   endif
enddo

r_CA_N =CA3_target - N1 !vector point from N in first pivot to alpha carbon in third pivot
r1=matmul(Tlab,r_CA_N) !transform vector into "lab frame" coordinate space
r1 = r1 - l_phi1; ! r1(2)=r1(2); r1(3)=r1(3)

Lomg1(1)=l_omg1+l_phi2*cosd(theta_omg1)
Lomg1(2)=-l_phi2*sind(theta_omg1)*cosd(omg1)
Lomg1(3)=l_phi2*sind(theta_omg1)*sind(omg1)

Lomg2(1)=l_omg2+l_phi3*cosd(theta_omg2)
Lomg2(2)=-l_phi3*sind(theta_omg2)*cosd(omg2)
Lomg2(3)=l_phi3*sind(theta_omg2)*sind(omg2)

Q_omg1(1)=l_psi1*cosd(theta_psi1)+Lomg1(1)
Q_omg1(2)=l_psi1*sind(theta_psi1)+Lomg1(2)
Q_omg1(3)=Lomg1(3)

Q_omg2(1)=l_psi2*cosd(theta_psi2)+Lomg2(1)
Q_omg2(2)=l_psi2*sind(theta_psi2)+Lomg2(2)
Q_omg2(3)=Lomg2(3)

e(1)=1; e(2)=0; e(3)=0

!test different phi values to see if any lead to an acceptable backbone
phi_lim = int((delta_phi_max - delta_phi_min)/delta_phi_step)

!perhaps a useful place to add parallelization (if each loop iteration independent, not sure of that yet)
do phi_iter=1, phi_lim
    delta_phi = delta_phi_min + (phi_iter-1)*delta_phi_step
    if(abs(delta_phi).le.2.0) goto 10
    phi1=REAL(phi1_old+delta_phi)
    call transformatrix(theta_phi1, phi1, Tphi1)
    t=matmul(transpose(Tphi1), r1)

    t_prime(1)=(t(1)**2+t(2)**2+t(3)**2)+(Q_omg1(1)**2+Q_omg1(2)**2+Q_omg1(3)**2)-(Q_omg2(1)**2+Q_omg2(2)**2+ & 
                Q_omg2(3)**2)-2*t(1)*(cosd(theta_psi1)*Q_omg1(1)+sind(theta_psi1)*Q_omg1(2))

    t_prime(2)=2*t(2)*(cosd(theta_psi1)*Q_omg1(2)-sind(theta_psi1)*Q_omg1(1))+2*t(3)*Q_omg1(3)
    t_prime(3)=2*t(3)*(cosd(theta_psi1)*Q_omg1(2)-sind(theta_psi1)*Q_omg1(1))-2*t(2)*Q_omg1(3)

    c11=t_prime(1)/sqrt(t_prime(2)**2+t_prime(3)**2)
    if(abs(c11).le.1) then
        if(t_prime(3).gt.0) then
            H=0
        else
            H=180
        endif
        psi(1)=atand(t_prime(2)/t_prime(3))-asind(c11)+H
        psi(2)=atand(t_prime(2)/t_prime(3))+asind(c11)-180.0+H
        do j=1, 2
            do while (psi(j).le.(-179.5).or.psi(j).gt.(180.5))
                if(psi(j).le.(-179.5)) then
                psi(j)=psi(j)+360.0
                elseif(psi(j).gt.(180.5)) then
                psi(j)=psi(j)-360.0
                endif
            end do
        enddo

        psi_num=0
        do j=1, 2
            psi_num=psi_num+1
            psi1(psi_num)=psi(j)
        enddo
    else
        goto 10
    endif

    do i=1, psi_num
        call transformatrix(theta_psi1, psi1(i), Tpsi1)
        q2=matmul(transpose(Tpsi1), t)
        q2=q2-Q_omg1
        call transformatrix(theta_omg1, omg1, Tomg1)
        q2_prime=matmul(transpose(Tomg1), q2)

        m(1)=sind(theta_phi2)*(sind(theta_psi2)*Q_omg2(1)-cosd(theta_psi2)*Q_omg2(2))
        m(2)=sind(theta_phi2)*Q_omg2(3)
        m(3)=cosd(theta_phi2)*(cosd(theta_psi2)*Q_omg2(1)+sind(theta_psi2)*Q_omg2(2))-q2_prime(1)
        c22=m(3)/sqrt(m(1)**2+m(2)**2)
        if(abs(c22).le.1) then
            if (m(2).eq.0) then
                psi(1) = 90 - asind(c22)
                psi(2) = -90 + asind(c22)
            else
                if(m(2).gt.0) then
                H=0
                else
                H=180
                endif
                psi(1)=atand(m(1)/m(2))-asind(c22)+H
                psi(2)=atand(m(1)/m(2))+asind(c22)-180.0+H
            endif
            do j=1, 2
                do while (psi(j).le.(-179.5).or.psi(j).gt.(180.5))
                if(psi(j).le.(-179.5)) then
                    psi(j)=psi(j)+360.0
                elseif(psi(j).gt.(180.5)) then
                    psi(j)=psi(j)-360.0
                endif
                end do
            enddo

            do j=1, 2
                n(1)=(sind(theta_phi2)*cosd(theta_psi2)+cosd(theta_phi2)*sind(theta_psi2)*cosd(psi(j)))*Q_omg2(1)
                n(1)=n(1)+(sind(theta_phi2)*sind(theta_psi2)-cosd(theta_phi2)*cosd(theta_psi2)*cosd(psi(j)))*Q_omg2(2)
                n(1)=n(1)-cosd(theta_phi2)*sind(psi(j))*Q_omg2(3)
                n(2)=sind(theta_psi2)*sind(psi(j))*Q_omg2(1)-cosd(theta_psi2)*sind(psi(j))*Q_omg2(2)+cosd(psi(j))*Q_omg2(3)
                n(3)=-q2_prime(2)
                c33=n(3)/sqrt(n(1)**2+n(2)**2)
                c44=q2_prime(3)
                if(abs(c33).le.1) then
                if(n(2).gt.0) then
                    H=0
                else
                    H=180
                endif
                phi(1)=atand(n(1)/n(2))-asind(c33)+H
                phi(2)=atand(n(1)/n(2))+asind(c33)-180.0+H
                do k=1, 2
                    do while (phi(k).le.(-179.5).or.phi(k).gt.(180.5))
                        if(phi(k).lt.(-179.5)) then
                            phi(k)=phi(k)+360.0
                        elseif(phi(k).gt.(180.5)) then
                            phi(k)=phi(k)-360.0
                        endif
                    end do
                enddo

                flag=0
                do k=1, 2
                    if(abs((n(1)*sind(phi(k))+n(2)*cosd(phi(k)))-c44).le.0.005) then
                        phi2=phi(k)
                        flag=1
                    endif
                enddo

                if(flag==0) goto 40
                call backbone_rotation_center_multiangle(Tgroup, Tbackbone_backup, ran_resi, phi1, psi1(i), &
                        phi2, psi(j), Tgroup_new, Tbackbone_backup_new)

                do k=1, Tgroup_new(ran_resi(1)-1)%cnum3
                    if(Tgroup_new(ran_resi(1)-1)%atype3(k)=="C") then
                        C10_temp(1:3)=Tgroup_new(ran_resi(1)-1)%coo3(1:3,k)
                    endif
                enddo
                do k=1, Tgroup_new(ran_resi(1))%cnum1
                    if(Tgroup_new(ran_resi(1))%atype1(k)=="N") then
                        N1_temp(1:3)=Tgroup_new(ran_resi(1))%coo1(1:3,k)
                    elseif(Tgroup_new(ran_resi(1))%atype1(k)=="CA") then
                        CA1_temp(1:3)=Tgroup_new(ran_resi(1))%coo1(1:3,k)
                    endif
                enddo
                do k=1, Tgroup_new(ran_resi(1))%cnum3
                    if(Tgroup_new(ran_resi(1))%atype3(k)=="C") then
                        C1_temp(1:3)=Tgroup_new(ran_resi(1))%coo3(1:3,k)
                    endif
                enddo
                do k=1, Tgroup_new(ran_resi(1)+1)%cnum1
                    if(Tgroup_new(ran_resi(1)+1)%atype1(k)=="N") then
                        N12_temp(1:3)=Tgroup_new(ran_resi(1)+1)%coo1(1:3,k)
                    endif
                enddo

                do k=1, Tgroup_new(ran_resi(2)-1)%cnum3
                    if(Tgroup_new(ran_resi(2)-1)%atype3(k)=="C") then
                        C20_temp(1:3)=Tgroup_new(ran_resi(2)-1)%coo3(1:3,k)
                    endif
                enddo
                do k=1, Tgroup_new(ran_resi(2))%cnum1
                    if(Tgroup_new(ran_resi(2))%atype1(k)=="N") then
                        N2_temp(1:3)=Tgroup_new(ran_resi(2))%coo1(1:3,k)
                    elseif(Tgroup_new(ran_resi(2))%atype1(k)=="CA") then
                        CA2_temp(1:3)=Tgroup_new(ran_resi(2))%coo1(1:3,k)
                    endif
                enddo
                do k=1, Tgroup_new(ran_resi(2))%cnum3
                    if(Tgroup_new(ran_resi(2))%atype3(k)=="C") then
                        C2_temp(1:3)=Tgroup_new(ran_resi(2))%coo3(1:3,k)
                    endif
                enddo
                do k=1, Tgroup_new(ran_resi(2)+1)%cnum1
                    if(Tgroup_new(ran_resi(2)+1)%atype1(k)=="N") then
                        N22_temp(1:3)=Tgroup_new(ran_resi(2)+1)%coo1(1:3,k)
                    endif
                enddo

                do k=1, Tgroup_new(ran_resi(3)-1)%cnum3
                    if(Tgroup_new(ran_resi(3)-1)%atype3(k)=="C") then
                        C30_temp(1:3)=Tgroup_new(ran_resi(3)-1)%coo3(1:3,k)
                    endif
                enddo
                do k=1, Tgroup_new(ran_resi(3))%cnum1
                    if(Tgroup_new(ran_resi(3))%atype1(k)=="N") then
                        N3_temp(1:3)=Tgroup_new(ran_resi(3))%coo1(1:3,k)
                    elseif(Tgroup_new(ran_resi(3))%atype1(k)=="CA") then
                        CA3_temp(1:3)=Tgroup_new(ran_resi(3))%coo1(1:3,k)
                    endif
                enddo

                !check if the new backbone bone angles satisfy the ramachandran map
                call calc_dihedral(C10_temp, N1_temp, CA1_temp, C1_temp, phi1_temp)
                call calc_dihedral(N1_temp, CA1_temp, C1_temp, N12_temp, psi1_temp)
                if(rama_map(NINT(phi1_temp),NINT(psi1_temp))==1.or.rama_map(NINT(phi1_temp),NINT(psi1_temp))==2) then
                    call calc_dihedral(C20_temp, N2_temp, CA2_temp, C2_temp, phi2_temp)
                    call calc_dihedral(N2_temp, CA2_temp, C2_temp, N22_temp, psi2_temp)
                    if(rama_map(NINT(phi2_temp),NINT(psi2_temp))==1.or.rama_map(NINT(phi2_temp),NINT(psi2_temp))==2) then

                        !check the target positions are very close to the actual positions (i.e. residues that weren't supposed to move did not actually move)
                        dist=sqrt((CA3_temp(1)-CA3_target(1))**2+(CA3_temp(2)-CA3_target(2))**2+(CA3_temp(3)-CA3_target(3))**2)
                        if(dist.le.(0.05)) then
                            call calc_dihedral(C30_temp, N3_temp, CA3_temp, C3, phi3_temp)
                            call calc_dihedral(N3_temp, CA3_temp, C3, N4, psi3_temp)
                            if(rama_map(NINT(phi3_temp),NINT(psi3_temp))==1.or.rama_map(NINT(phi3_temp),NINT(psi3_temp))==2) then
                            call calc_rmsd(Tgroup_new, original_pep, rmsd)
                            call backbonemove_criterion(rmsd, feedback) !check RMSD is within allowable range
                            if(feedback==1) then
                                !all checks passed, so add to possible list!
                                phipsi_num=phipsi_num+1
                                if(phipsi_num.gt.20) then
                                    phipsi_num=20
                                    goto 20
                                endif
                                do k=1, pep_res
                                    group_candidates(k,phipsi_num)=Tgroup_new(k)
                                    backbone_backup_candidates(k,phipsi_num)%coo=Tbackbone_backup_new(k)%coo
                                enddo
                            endif
                            endif
                        endif
                    endif
                endif
                endif
    40                continue
            enddo
        endif
    enddo
10      continue

enddo
20   continue

return
end subroutine concertedrotation_center

subroutine concertedrotation_whole(group_pep, original_pep, backbone_backup, ran_resi, phipsi_num, group_candidates, backbone_backup_candidates)
implicit none
integer                         :: ran_resi(3), k
integer                         :: phipsi_num, flag1, flag2
type(groupdetails)              :: group_pep(pep_res), group_candidates(pep_res, 20)
type(databackup)                :: backbone_backup(pep_res), backbone_backup_candidates(pep_res,20)
type(groupdetails)              :: Tgroup1(pep_res), Tgroup2(pep_res), original_pep(pep_res)
type(databackup)                :: Tbackbone_backup1(pep_res), Tbackbone_backup2(pep_res)

phipsi_num=0
call concertedrotation_N2Cterm(group_pep, original_pep, backbone_backup, ran_resi, flag1, Tgroup1, Tbackbone_backup1)
if(flag1.ne.0) then
    call concertedrotation_C2Nterm(Tgroup1, original_pep, Tbackbone_backup1, ran_resi, flag2, Tgroup2, Tbackbone_backup2)
    if(flag2.ne.0) then
        phipsi_num=phipsi_num+1
        do k=1, pep_res
            group_candidates(k,phipsi_num)=Tgroup2(k)
            backbone_backup_candidates(k,phipsi_num)%coo=Tbackbone_backup2(k)%coo
        enddo
    endif
endif

call concertedrotation_C2Nterm(group_pep, original_pep, backbone_backup, ran_resi, flag1, Tgroup1, Tbackbone_backup1)
if(flag1.ne.0) then
    call concertedrotation_N2Cterm(Tgroup1, original_pep, Tbackbone_backup1, ran_resi, flag2, Tgroup2, Tbackbone_backup2)
    if(flag2.ne.0) then
        phipsi_num=phipsi_num+1
        do k=1, pep_res
            group_candidates(k,phipsi_num)=Tgroup2(k)
            backbone_backup_candidates(k,phipsi_num)%coo=Tbackbone_backup2(k)%coo
        enddo
    endif
endif

return
end subroutine concertedrotation_whole

subroutine concertedrotation_Nterm(group_pep, original_pep, backbone_backup, ran_resi, psiphi_num, group_candidates, backbone_backup_candidates)
implicit none
integer                             :: ran_resi(3), i, k, ic
integer                             :: psiphi_num, feedback, Tresinum
real                                :: N30(3), C3(3), CA3(3), N3(3), C32(3)
real                                :: N20(3), C2(3), CA2(3), N2(3), C22(3)
real                                :: rmsd, ran2, delta_psi3, delta_phi3, delta_psi2, delta_phi2
real                                :: psi3, phi3, psi2, phi2, Tpsi3, Tphi3, Tpsi2, Tphi2
type(groupdetails)                  :: group_pep(pep_res), group_candidates(pep_res, 20), Tgroup(pep_res), original_pep(pep_res)
type(databackup)                    :: backbone_backup(pep_res), backbone_backup_candidates(pep_res, 20), Tbackbone_backup(pep_res)

do k=1, group_pep(ran_resi(3)+1)%cnum1
    if(group_pep(ran_resi(3)+1)%atype1(k)=="N") then
        N30(1:3)=group_pep(ran_resi(3)+1)%coo1(1:3,k)
    endif
enddo
do k=1, group_pep(ran_resi(3))%cnum3
    if(group_pep(ran_resi(3))%atype3(k)=="C") then
        C3(1:3)=group_pep(ran_resi(3))%coo3(1:3,k)
    endif
enddo
do k=1, group_pep(ran_resi(3))%cnum1
    if(group_pep(ran_resi(3))%atype1(k)=="CA") then
        CA3(1:3)=group_pep(ran_resi(3))%coo1(1:3,k)
    elseif(group_pep(ran_resi(3))%atype1(k)=="N") then
        N3(1:3)=group_pep(ran_resi(3))%coo1(1:3,k)
    endif
enddo
do k=1, group_pep(ran_resi(3)-1)%cnum3
    if(group_pep(ran_resi(3)-1)%atype3(k)=="C") then
        C32(1:3)=group_pep(ran_resi(3)-1)%coo3(1:3,k)
    endif
enddo

do k=1, group_pep(ran_resi(2)+1)%cnum1
    if(group_pep(ran_resi(2)+1)%atype1(k)=="N") then
        N20(1:3)=group_pep(ran_resi(2)+1)%coo1(1:3,k)
    endif
enddo
do k=1, group_pep(ran_resi(2))%cnum3
    if(group_pep(ran_resi(2))%atype3(k)=="C") then
        C2(1:3)=group_pep(ran_resi(2))%coo3(1:3,k)
    endif
enddo
do k=1, group_pep(ran_resi(2))%cnum1
    if(group_pep(ran_resi(2))%atype1(k)=="CA") then
        CA2(1:3)=group_pep(ran_resi(2))%coo1(1:3,k)
    elseif(group_pep(ran_resi(2))%atype1(k)=="N") then
        N2(1:3)=group_pep(ran_resi(2))%coo1(1:3,k)
    endif
enddo
do k=1, group_pep(ran_resi(2)-1)%cnum3
    if(group_pep(ran_resi(2)-1)%atype3(k)=="C") then
        C22(1:3)=group_pep(ran_resi(2)-1)%coo3(1:3,k)
    endif
enddo

call calc_dihedral(N30, C3, CA3, N3, psi3)
call calc_dihedral(C3, CA3, N3, C32, phi3)

call calc_dihedral(N20, C2, CA2, N2, psi2)
call calc_dihedral(C2, CA2, N2, C22, phi2)

psiphi_num=0
Tresinum=ran_resi(3)-ran_resi(1)
do i=1, 20
    if(Tresinum.le.4) then
        call RANDOM_NUMBER(ran2)
        ic=int(ran2*flavoredregion_number)+1

        Tpsi3=flavored_region(ic,2)
        Tphi3=flavored_region(ic,1)

        delta_psi3=Tpsi3-psi3
        delta_phi3=Tphi3-phi3

        call RANDOM_NUMBER(ran2)
        ic=int(ran2*flavoredregion_number)+1

        Tpsi2=flavored_region(ic,2)
        Tphi2=flavored_region(ic,1)

        delta_psi2=Tpsi2-psi2
        delta_phi2=Tphi2-phi2
    else
        call RANDOM_NUMBER(ran2); delta_psi3=ran2*6.9-3.45
        call RANDOM_NUMBER(ran2); delta_phi3=ran2*6.9-3.45
        Tpsi3=psi3+delta_psi3; Tphi3=phi3+delta_phi3
        if(Tpsi3.ge.180.5) then
            Tpsi3=Tpsi3-360.0
            delta_psi3=Tpsi3-psi3
        elseif(Tpsi3.le.(-179.5)) then
            Tpsi3=Tpsi3+360.0
            delta_psi3=Tpsi3-psi3
        endif
        if(Tphi3.ge.180.5) then
            Tphi3=Tphi3-360.0
            delta_phi3=Tphi3-phi3
        elseif(Tphi3.le.(-179.5)) then
            Tphi3=Tphi3+360.0
            delta_phi3=Tphi3-phi3
        endif
        if(rama_map(NINT(Tphi3), NINT(Tpsi3))==0) goto 20

        call RANDOM_NUMBER(ran2); delta_psi2=ran2*10.9-5.45
        call RANDOM_NUMBER(ran2); delta_phi2=ran2*10.9-5.45
        Tpsi2=psi2+delta_psi2; Tphi2=phi2+delta_phi2
        if(Tpsi2.ge.180.5) then
            Tpsi2=Tpsi2-360.0
            delta_psi2=Tpsi2-psi2
        elseif(Tpsi2.le.(-179.5)) then
            Tpsi2=Tpsi2+360.0
            delta_psi2=Tpsi2-psi2
        endif
        if(Tphi2.ge.180.5) then
            Tphi2=Tphi2-360.0
            delta_phi2=Tphi2-phi2
        elseif(Tphi2.le.(-179.5)) then
            Tphi2=Tphi2+360.0
            delta_phi2=Tphi2-phi2
        endif
        if(rama_map(NINT(Tphi2), NINT(Tpsi2))==0) goto 20
    endif

    call backbone_rotation_Nterm_multiangle(group_pep, backbone_backup, ran_resi, delta_psi3, delta_phi3, &
    delta_psi2, delta_phi2, Tgroup, Tbackbone_backup)

    call calc_rmsd(Tgroup, original_pep, rmsd)
    call backbonemove_criterion(rmsd, feedback)
    if(feedback==1) then
        psiphi_num=psiphi_num+1
        do k=1, pep_res
            group_candidates(k,psiphi_num)=Tgroup(k)
            backbone_backup_candidates(k,psiphi_num)%coo=Tbackbone_backup(k)%coo
        enddo
        if(psiphi_num.eq.number4peptide_candidates) goto 10
    endif
20  continue
enddo
10 continue

return
end subroutine concertedrotation_Nterm

subroutine concertedrotation_C2Nterm(group_pep, original_pep,  backbone_backup, ran_resi, flag, Tgroup, Tbackbone_backup)
implicit none
integer                             :: ran_resi(3), i, k
integer                             :: flag, feedback
real                                :: N20(3), C2(3), CA2(3), N2(3), C22(3)
real                                :: rmsd, ran2, delta_psi3, delta_phi3, delta_psi2, delta_phi2
real                                :: psi2, phi2, Tpsi2, Tphi2
type(groupdetails)                  :: group_pep(pep_res), Tgroup(pep_res), original_pep(pep_res)
type(databackup)                    :: backbone_backup(pep_res), Tbackbone_backup(pep_res)

do k=1, group_pep(ran_resi(2)+1)%cnum1
    if(group_pep(ran_resi(2)+1)%atype1(k)=="N") then
        N20(1:3)=group_pep(ran_resi(2)+1)%coo1(1:3,k)
    endif
enddo
do k=1, group_pep(ran_resi(2))%cnum3
    if(group_pep(ran_resi(2))%atype3(k)=="C") then
        C2(1:3)=group_pep(ran_resi(2))%coo3(1:3,k)
    endif
enddo
do k=1, group_pep(ran_resi(2))%cnum1
    if(group_pep(ran_resi(2))%atype1(k)=="CA") then
        CA2(1:3)=group_pep(ran_resi(2))%coo1(1:3,k)
    elseif(group_pep(ran_resi(2))%atype1(k)=="N") then
        N2(1:3)=group_pep(ran_resi(2))%coo1(1:3,k)
    endif
enddo
do k=1, group_pep(ran_resi(2)-1)%cnum3
    if(group_pep(ran_resi(2)-1)%atype3(k)=="C") then
        C22(1:3)=group_pep(ran_resi(2)-1)%coo3(1:3,k)
    endif
enddo

call calc_dihedral(N20, C2, CA2, N2, psi2)
call calc_dihedral(C2, CA2, N2, C22, phi2)

flag=0
do i=1, 20
    delta_psi3=0.0
    call RANDOM_NUMBER(ran2); delta_phi3=ran2*4.9-2.45
    call RANDOM_NUMBER(ran2); delta_psi2=ran2*6.9-3.45
    call RANDOM_NUMBER(ran2); delta_phi2=ran2*6.9-3.45
    Tpsi2=psi2+delta_psi2; Tphi2=phi2+delta_phi2
    if(Tpsi2.ge.180.5) then
        Tpsi2=Tpsi2-360.0
        delta_psi2=Tpsi2-psi2
    elseif(Tpsi2.le.(-179.5)) then
        Tpsi2=Tpsi2+360.0
        delta_psi2=Tpsi2-psi2
    endif
    if(Tphi2.ge.180.5) then
        Tphi2=Tphi2-360.0
        delta_phi2=Tphi2-phi2
    elseif(Tphi2.le.(-179.5)) then
        Tphi2=Tphi2+360.0
        delta_phi2=Tphi2-phi2
    endif
    if(rama_map(NINT(Tphi2), NINT(Tpsi2))==0) goto 20

    call backbone_rotation_Nterm_multiangle(group_pep, backbone_backup, ran_resi, delta_psi3, delta_phi3, &
    delta_psi2, delta_phi2, Tgroup, Tbackbone_backup)

    call calc_rmsd(Tgroup, original_pep, rmsd)
    call backbonemove_criterion(rmsd, feedback)
    if(feedback==1) then
        flag=flag+1
        goto 10
    endif
20  continue
enddo
10 continue

return
end subroutine concertedrotation_C2Nterm

subroutine concertedrotation_Cterm(group_pep, original_pep, backbone_backup, ran_resi, phipsi_num, &
   group_candidates, backbone_backup_candidates)
implicit none
integer                             :: ran_resi(3), i, k, ic
integer                             :: phipsi_num, feedback, Tresinum
real                                :: C10(3), N1(3), CA1(3), C1(3), N12(3)
real                                :: C20(3), N2(3), CA2(3), C2(3), N22(3)
real                                :: rmsd, ran2, delta_phi1, delta_psi1, delta_phi2, delta_psi2
real                                :: phi1, psi1, phi2, psi2, Tphi1, Tpsi1, Tphi2, Tpsi2
type(groupdetails)                  :: group_pep(pep_res), group_candidates(pep_res, 20), Tgroup(pep_res), original_pep(pep_res)
type(databackup)                    :: backbone_backup(pep_res), backbone_backup_candidates(pep_res, 20), Tbackbone_backup(pep_res)

do k=1, group_pep(ran_resi(1)-1)%cnum3
    if(group_pep(ran_resi(1)-1)%atype3(k)=="C") then
        C10(1:3)=group_pep(ran_resi(1)-1)%coo3(1:3,k)
    endif
enddo
do k=1, group_pep(ran_resi(1))%cnum1
    if(group_pep(ran_resi(1))%atype1(k)=="N") then
        N1(1:3)=group_pep(ran_resi(1))%coo1(1:3,k)
    elseif(group_pep(ran_resi(1))%atype1(k)=="CA") then
        CA1(1:3)=group_pep(ran_resi(1))%coo1(1:3,k)
    endif
enddo
do k=1, group_pep(ran_resi(1))%cnum3
    if(group_pep(ran_resi(1))%atype3(k)=="C") then
        C1(1:3)=group_pep(ran_resi(1))%coo3(1:3,k)
    endif
enddo
do k=1, group_pep(ran_resi(1)+1)%cnum1
    if(group_pep(ran_resi(1)+1)%atype1(k)=="N") then
        N12(1:3)=group_pep(ran_resi(1)+1)%coo1(1:3,k)
    endif
enddo

do k=1, group_pep(ran_resi(2)-1)%cnum3
    if(group_pep(ran_resi(2)-1)%atype3(k)=="C") then
        C20(1:3)=group_pep(ran_resi(2)-1)%coo3(1:3,k)
    endif
enddo
do k=1, group_pep(ran_resi(2))%cnum1
    if(group_pep(ran_resi(2))%atype1(k)=="N") then
        N2(1:3)=group_pep(ran_resi(2))%coo1(1:3,k)
    elseif(group_pep(ran_resi(2))%atype1(k)=="CA") then
        CA2(1:3)=group_pep(ran_resi(2))%coo1(1:3,k)
    endif
enddo
do k=1, group_pep(ran_resi(2))%cnum3
    if(group_pep(ran_resi(2))%atype3(k)=="C") then
        C2(1:3)=group_pep(ran_resi(2))%coo3(1:3,k)
    endif
enddo
do k=1, group_pep(ran_resi(2)+1)%cnum1
    if(group_pep(ran_resi(2)+1)%atype1(k)=="N") then
        N22(1:3)=group_pep(ran_resi(2)+1)%coo1(1:3,k)
    endif
enddo

call calc_dihedral(C10, N1, CA1, C1, phi1)
call calc_dihedral(N1, CA1, C1, N12, psi1)
call calc_dihedral(C20, N2, CA2, C2, phi2)
call calc_dihedral(N2, CA2, C2, N22, psi2)

phipsi_num=0
Tresinum=ran_resi(3)-ran_resi(1)
do i=1, 20
    if(Tresinum.le.4) then
        call RANDOM_NUMBER(ran2)
        ic=int(ran2*flavoredregion_number)+1
        Tphi1=flavored_region(ic,1)
        Tpsi1=flavored_region(ic,2)

        delta_phi1=Tphi1-phi1
        delta_psi1=Tpsi1-psi1

        call RANDOM_NUMBER(ran2)
        ic=int(ran2*flavoredregion_number)+1
        Tphi2=flavored_region(ic,1)
        Tpsi2=flavored_region(ic,2)

        delta_phi2=Tphi2-phi2
        delta_psi2=Tpsi2-psi2
    else
        call RANDOM_NUMBER(ran2); delta_phi1=ran2*6.8-3.4
        call RANDOM_NUMBER(ran2); delta_psi1=ran2*6.8-3.4
        Tphi1=phi1+delta_phi1; Tpsi1=psi1+delta_psi1
        if(Tphi1.ge.180.5) then
            Tphi1=Tphi1-360.0
            delta_phi1=Tphi1-phi1
        elseif(Tphi1.le.(-179.5)) then
            Tphi1=Tphi1+360.0
            delta_phi1=Tphi1-phi1
        endif
        if(Tpsi1.ge.180.5) then
            Tpsi1=Tpsi1-360.0
            delta_psi1=Tpsi1-psi1
        elseif(Tpsi1.le.(-179.5)) then
            Tpsi1=Tpsi1+360.0
            delta_psi1=Tpsi1-psi1
        endif
        if(rama_map(NINT(Tphi1), NINT(Tpsi1))==0) goto 20

        call RANDOM_NUMBER(ran2); delta_phi2=ran2*10.8-5.4
        call RANDOM_NUMBER(ran2); delta_psi2=ran2*10.8-5.4
        Tphi2=phi2+delta_phi2; Tpsi2=psi2+delta_psi2
        if(Tphi2.ge.180.5) then
            Tphi2=Tphi2-360.0
            delta_phi2=Tphi2-phi2
        elseif(Tphi2.le.(-179.5)) then
            Tphi2=Tphi2+360.0
            delta_phi2=Tphi2-phi2
        endif
        if(Tpsi2.ge.180.5) then
            Tpsi2=Tpsi2-360.0
            delta_psi2=Tpsi2-psi2
        elseif(Tpsi2.le.(-179.5)) then
            Tpsi2=Tpsi2+360.0
            delta_psi2=Tpsi2-psi2
        endif
        if(rama_map(NINT(Tphi2), NINT(Tpsi2))==0) goto 20
    endif

    call backbone_rotation_Cterm_multiangle(group_pep, backbone_backup, ran_resi, delta_phi1, delta_psi1,  & 
    delta_phi2, delta_psi2, Tgroup, Tbackbone_backup)

    call calc_rmsd(Tgroup, original_pep, rmsd)
    call backbonemove_criterion(rmsd, feedback)
    if(feedback==1) then
        phipsi_num=phipsi_num+1
        do k=1, pep_res
            group_candidates(k,phipsi_num)=Tgroup(k)
            backbone_backup_candidates(k,phipsi_num)%coo=Tbackbone_backup(k)%coo
        enddo
        if(phipsi_num.eq.number4peptide_candidates) goto 10
    endif
20  continue
enddo
10    continue

return
end subroutine concertedrotation_Cterm

subroutine concertedrotation_N2Cterm(group_pep, original_pep, backbone_backup, ran_resi, flag, Tgroup, Tbackbone_backup)
implicit none
integer                         :: ran_resi(3), i, k
integer                         :: flag, feedback
real                            :: C20(3), N2(3), CA2(3), C2(3), N22(3)
real                            :: rmsd, ran2, delta_phi1, delta_psi1, delta_phi2, delta_psi2
real                            :: phi2, psi2, Tphi2, Tpsi2
type(groupdetails)              :: group_pep(pep_res), Tgroup(pep_res), original_pep(pep_res)
type(databackup)                :: backbone_backup(pep_res), Tbackbone_backup(pep_res)

do k=1, group_pep(ran_resi(2)-1)%cnum3
    if(group_pep(ran_resi(2)-1)%atype3(k)=="C") then
        C20(1:3)=group_pep(ran_resi(2)-1)%coo3(1:3,k)
    endif
enddo
do k=1, group_pep(ran_resi(2))%cnum1
    if(group_pep(ran_resi(2))%atype1(k)=="N") then
        N2(1:3)=group_pep(ran_resi(2))%coo1(1:3,k)
    elseif(group_pep(ran_resi(2))%atype1(k)=="CA") then
        CA2(1:3)=group_pep(ran_resi(2))%coo1(1:3,k)
    endif
enddo
do k=1, group_pep(ran_resi(2))%cnum3
    if(group_pep(ran_resi(2))%atype3(k)=="C") then
        C2(1:3)=group_pep(ran_resi(2))%coo3(1:3,k)
    endif
enddo
do k=1, group_pep(ran_resi(2)+1)%cnum1
    if(group_pep(ran_resi(2)+1)%atype1(k)=="N") then
        N22(1:3)=group_pep(ran_resi(2)+1)%coo1(1:3,k)
    endif
enddo

call calc_dihedral(C20, N2, CA2, C2, phi2)
call calc_dihedral(N2, CA2, C2, N22, psi2)

flag=0
do i=1, 20
    delta_phi1=0.0
    call RANDOM_NUMBER(ran2); delta_psi1=ran2*4.9-2.45
    call RANDOM_NUMBER(ran2); delta_phi2=ran2*6.9-3.45
    call RANDOM_NUMBER(ran2); delta_psi2=ran2*6.9-3.45
    Tphi2=phi2+delta_phi2; Tpsi2=psi2+delta_psi2
    if(Tphi2.ge.180.5) then
        Tphi2=Tphi2-360.0
        delta_phi2=Tphi2-phi2
    elseif(Tphi2.le.(-179.5)) then
        Tphi2=Tphi2+360.0
        delta_phi2=Tphi2-phi2
    endif
    if(Tpsi2.ge.180.5) then
        Tpsi2=Tpsi2-360.0
        delta_psi2=Tpsi2-psi2
    elseif(Tpsi2.le.(-179.5)) then
        Tpsi2=Tpsi2+360.0
        delta_psi2=Tpsi2-psi2
    endif
    if(rama_map(NINT(Tphi2), NINT(Tpsi2)).ne.0) then
        call backbone_rotation_Cterm_multiangle(group_pep, backbone_backup, ran_resi, delta_phi1, delta_psi1, &
        delta_phi2, delta_psi2, Tgroup, Tbackbone_backup)

        call calc_rmsd(Tgroup, original_pep, rmsd)
        call backbonemove_criterion(rmsd, feedback)
        if(feedback==1) then
            flag=1
            return
        endif
    endif
enddo

return
end subroutine concertedrotation_N2Cterm

end module advanced_function