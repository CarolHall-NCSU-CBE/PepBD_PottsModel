module input

use datatypes
use sys_vars

contains
subroutine inputfile
implicit none
integer*8             :: status, str_status, len, l_index, r_index
character*40          :: name, val
character*80          :: line


!read input file
open(10, file="input.txt", status="old")
    do while(.true.)
        read(10, '(A)', iostat=status) line
        if(status.ne.0) exit !end of file has been reached

        !split entry into the two parts (I do this instead of a list read because LIB_PATH variable may have a "/" character, which is not read properly in a list read)
        len = LEN_TRIM(line) !get length of input, without trailing spaces
        l_index = INDEX(line(1:len), ' ') !get location if first space
        if (l_index.eq.0) then !the values in the line not separated by spaces, so try tabs
            l_index = INDEX(line(1:len), achar(9)) !achar(9) is apparently the ASCII table entry for the tab character
            r_index = INDEX(line(1:len), achar(9), .true.) 
        else
            r_index = INDEX(line(1:len), ' ', .true.) !get location of last space, not considering trailing spaces
        endif

        !make sure the line was properly split
        if (l_index.eq.0 .or. r_index.eq.0) then
            write(*,*) 'Unable to parse data in input file!'
            write(*,*) 'Make sure each line has only two values, and there is either a space or tab between values'
            write(*,*) 'Terminating program'
        endif

        name = line(1:l_index-1) !get first entry in line
        val = line(r_index+1:len) !get second entry in line

        if(name=="PDBFILE") then
            system_pdb = val
        elseif(name=='LIB_PATH') then
            lib_path = val
        elseif(name=="PEP_RES") then
            read(val,*,iostat=str_status)  pep_res
        elseif(name=="RECEPTOR_NAME") then
            read(val,*,iostat=str_status) receptor_name
        elseif(name=="DEBUG_FLAG") then
            read(val,*,iostat=str_status) debug_flag
        elseif(name=="SOLVATION_FLAG") then
            read(val,*,iostat=str_status) solvation_flag
        elseif(name=="NP_FLAG") then
            read(val,*,iostat=str_status) NP_flag
        elseif(name=="WEIGHTING_FACTOR") then
            read(val,*,iostat=str_status) weighting_factor
        elseif(name=="INSERT_PEP_RIGOR") then
            read(val,*,iostat=str_status) insert_pep_rigor
        elseif(name=="EVAL_TWOBODY_FLAG") then
            read(val,*,iostat=str_status) eval_twobody_flag
        elseif(name=="CG_RBORN") then
            read(val,*,iostat=str_status) cg_rborn
        else
            open(20, file="error.txt", access="append")
            write(20,*) "Error in input file!"
            write(20,*) "Unrecognized parameter ", name, " in input file, terminating program."
            close(20)
            stop
        endif
end do



return
end subroutine inputfile
  
end module input