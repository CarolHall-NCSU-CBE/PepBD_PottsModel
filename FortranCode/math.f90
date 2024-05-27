module math
use datatypes
use sys_vars
use omp_lib
external SSYEVR
external SLAMCH

contains

subroutine DIST_SQUARE(ri, rj, dij)
implicit none
real            , intent(in)        :: ri(3), rj(3)
real                                :: dij

dij = (ri(1)-rj(1))**2 + (ri(2)-rj(2))**2 + (ri(3)-rj(3))**2
return

end subroutine DIST_SQUARE

subroutine vec_norm(v, norm)
implicit none
real        :: v(3), norm

norm = (v(1)*v(1) + v(2)*v(2) + v(3)*v(3))**0.5
return

end subroutine vec_norm

subroutine normalize_vector(v)
implicit none
real    :: v(3), mag

call vec_norm(v, mag)
v = v/mag

return

end subroutine normalize_vector

subroutine cross_product(u, v, w)
implicit none
real                ::u(3), v(3), w(3)

w(1) = u(2)*v(3) - u(3)*v(2)
w(2) = u(3)*v(1) - u(1)*v(3)
w(3) = u(1)*v(2) - u(2)*v(1)
return

end subroutine cross_product

subroutine normalvector(rsta, rmid, rend, r_nor)
implicit none
real                              :: rsta(3), rmid(3), rend(3), r_nor(3)
real                              :: a(3), b(3)

a=rsta-rmid
b=rend-rmid

call cross_product(a,b,r_nor)

return
end subroutine normalvector

subroutine vectorrotation(rsta, rend, m)
implicit none
real                              :: rsta(3), rend(3), m(3,3)
real                              :: r_cropro(3), a(3), a1(3,3), a2(3,3), a3(3,3)
real                              :: absrsta, absrend, r_dotpro, cos, sin

absrsta=SQRT(rsta(1)*rsta(1)+rsta(2)*rsta(2)+rsta(3)*rsta(3))
absrend=SQRT(rend(1)*rend(1)+rend(2)*rend(2)+rend(3)*rend(3))

r_dotpro=dot_product(rsta, rend)

r_cropro(1)=rsta(2)*rend(3)-rsta(3)*rend(2)
r_cropro(2)=rsta(3)*rend(1)-rsta(1)*rend(3)
r_cropro(3)=rsta(1)*rend(2)-rsta(2)*rend(1)

cos=r_dotpro/(absrsta*absrend)
sin=SQRT(r_cropro(1)*r_cropro(1)+r_cropro(2)*r_cropro(2)+r_cropro(3)*r_cropro(3))/(absrsta*absrend)

a(1)=r_cropro(1)/SQRT(r_cropro(1)*r_cropro(1)+r_cropro(2)*r_cropro(2)+r_cropro(3)*r_cropro(3))
a(2)=r_cropro(2)/SQRT(r_cropro(1)*r_cropro(1)+r_cropro(2)*r_cropro(2)+r_cropro(3)*r_cropro(3))
a(3)=r_cropro(3)/SQRT(r_cropro(1)*r_cropro(1)+r_cropro(2)*r_cropro(2)+r_cropro(3)*r_cropro(3))

a1(1,1)=a(1)*a(1)
a1(1,2)=a(1)*a(2)
a1(1,3)=a(1)*a(3)
a1(2,1)=a(2)*a(1)
a1(2,2)=a(2)*a(2)
a1(2,3)=a(2)*a(3)
a1(3,1)=a(3)*a(1)
a1(3,2)=a(3)*a(2)
a1(3,3)=a(3)*a(3)

a2=-a1
a2(1,1)=1+a2(1,1)
a2(2,2)=1+a2(2,2)
a2(3,3)=1+a2(3,3)
a2=cos*a2

a3(1,1)=0.0
a3(1,2)=-a(3)
a3(1,3)=a(2)
a3(2,1)=a(3)
a3(2,2)=0.0
a3(2,3)=-a(1)
a3(3,1)=-a(2)
a3(3,2)=a(1)
a3(3,3)=0
a3=sin*a3

m=a1+a2+a3

return
end subroutine vectorrotation

subroutine axisrotation(a, cos, sin, m)
implicit none
real                           :: cos, sin
real                           :: a(3), a1(3,3), a2(3,3), a3(3,3)
real                           :: m(3,3)

a1(1,1)=a(1)*a(1)
a1(1,2)=a(1)*a(2)
a1(1,3)=a(1)*a(3)
a1(2,1)=a(2)*a(1)
a1(2,2)=a(2)*a(2)
a1(2,3)=a(2)*a(3)
a1(3,1)=a(3)*a(1)
a1(3,2)=a(3)*a(2)
a1(3,3)=a(3)*a(3)

a2=-a1
a2(1,1)=1+a2(1,1)
a2(2,2)=1+a2(2,2)
a2(3,3)=1+a2(3,3)
a2=cos*a2

a3(1,1)=0.0
a3(1,2)=-a(3)
a3(1,3)=a(2)
a3(2,1)=a(3)
a3(2,2)=0.0
a3(2,3)=-a(1)
a3(3,1)=-a(2)
a3(3,2)=a(1)
a3(3,3)=0
a3=sin*a3

m=a1+a2+a3

return
end subroutine axisrotation

subroutine transformatrix(bond_angle, dihedral_angle, T)
implicit none
real                           :: bond_angle, dihedral_angle
real                           :: T(3,3)

T(1,1)=cosd(bond_angle)
T(1,2)=sind(bond_angle)
T(1,3)=0.0
T(2,1)=-sind(bond_angle)*cosd(dihedral_angle)
T(2,2)=cosd(bond_angle)*cosd(dihedral_angle)
T(2,3)=sind(dihedral_angle)
T(3,1)=sind(bond_angle)*sind(dihedral_angle)
T(3,2)=-cosd(bond_angle)*sind(dihedral_angle)
T(3,3)=cosd(dihedral_angle)

return
end subroutine transformatrix

subroutine calc_dihedral(p1, p2, p3, p4, angle)
implicit none
real                              :: p1(3), p2(3), p3(3), p4(3)
real                              :: angle, angle_T1, angle_T2
real                              :: rsta(3), rend(3), rmid(3)
real                              :: r_1(3), r_2(3)
real                              :: absrsta, absrend, r_dotpro, ratio

call normalvector(p1, p2, p3, rend)
call normalvector(p2, p3, p4, rsta)

absrsta=SQRT(rsta(1)*rsta(1)+rsta(2)*rsta(2)+rsta(3)*rsta(3))
absrend=SQRT(rend(1)*rend(1)+rend(2)*rend(2)+rend(3)*rend(3))
r_dotpro=dot_product(rsta, rend)

!added this logic since we could occasionally get the ratio to be greater than 1 due to floating point error, which would cause program to crash
ratio = r_dotpro/(absrsta*absrend)
if (abs(ratio).gt.1) then
    if (ratio.lt.0) then
        ratio = -1
    else
        ratio = 1
    endif
endif
angle_T1=acosd(ratio)

!round dihedral to 180 or 0 if it is very close to those numbers
if(abs(180.0-angle_T1).le.(0.5)) then
    angle=180.0
elseif(abs(angle_T1).le.(0.5)) then
    angle=0.0
else !will do some extra analysis

    !get vector normal to plane formed by vector rsta and vector rmid
    rmid=0.0
    r_2=p3-p2
    call normalvector(rsta, rmid, rend, r_1) 

    !get dot_
    absrsta=SQRT(r_1(1)*r_1(1)+r_1(2)*r_1(2)+r_1(3)*r_1(3))
    absrend=SQRT(r_2(1)*r_2(1)+r_2(2)*r_2(2)+r_2(3)*r_2(3))
    r_dotpro=dot_product(r_1, r_2)/(absrsta*absrend)

    if(abs(r_dotpro-1.0).le.(0.1)) then
        angle_T2=0.00
    elseif(abs(r_dotpro+1.0).le.(0.1)) then
        angle_T2=180
    else
        angle_T2=acosd(r_dotpro)
    endif

    if(angle_T2.gt.90) then
        angle=-angle_T1
    else
        angle=angle_T1
    endif
endif

return
end subroutine calc_dihedral

subroutine rotation_matrix(axis, theta, m)
implicit none

real            :: m(3,3), theta, s, c, zero=0.0, one=1.0
integer         :: axis

s = sin(theta); c = cos(theta)

!get rotation matrix
if (axis.eq.0) then         !rotating about x-axis
m = reshape((/  one,    zero,    zero, &
                zero,    c,       -s, &
                zero,    s,       c/), shape(m))

elseif(axis.eq.1) then      !rotating about y-axis
m = reshape((/  c,      zero,   s, &
                zero,  one,     zero, &
                 -s,   zero,    c/), shape(m))

else                        !rotating about z-axis
m = reshape((/  c ,     -s,     zero, &
                s,      c,      zero, &
                zero,   zero,   one/), shape(m))
endif

return
    
end subroutine rotation_matrix

subroutine Euler_Rotation_Matrix(angles, rot_m)
implicit none
real                        :: rot_m(3,3), angles(3)

!get matrix to perform Euler rotations
!Using matrix defined here: https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Variational_Principles_in_Classical_Mechanics_(Cline)/13%3A_Rigid-body_Rotation/13.13%3A_Euler_Angles

rot_m(1,1) =  COS(angles(1))*COS(angles(3)) - SIN(angles(1))*COS(angles(2))*SIN(angles(3))
rot_m(2,1) =  SIN(angles(1))*COS(angles(3)) + COS(angles(1))*COS(angles(2))*SIN(angles(3))
rot_m(3,1) =  SIN(angles(2))*SIN(angles(3))
rot_m(1,2) = -COS(angles(1))*SIN(angles(3)) - SIN(angles(1))*COS(angles(2))*COS(angles(3))
rot_m(2,2) = -SIN(angles(1))*SIN(angles(3)) + COS(angles(1))*COS(angles(2))*COS(angles(3))
rot_m(3,2) =  SIN(angles(2))*COS(angles(3))
rot_m(1,3) =  SIN(angles(1))*SIN(angles(2))
rot_m(2,3) = -COS(angles(1))*SIN(angles(2))
rot_m(3,3) =  COS(angles(2))

return 
end subroutine Euler_Rotation_Matrix

subroutine decompose_mom_tensor(group_pep, obj_axes, moments)
implicit none

real                        :: mom_tensor(3,3), obj_axes(3,3), moments(3), v(3), mag
integer                     :: ax, count
type(groupdetails)          :: group_pep(pep_res)

!debug variables
!real                        :: mom_tensor_transform(3,3), obj_axes_t(3,3), D(3,3), mom_tensor_orig(3,3)

!!!!!!!! LAPACK SSYEVR PARAMETERS !!!!!

!!     .. Parameters ..
integer            N, NSELECT
PARAMETER        ( N = 3, NSELECT = 3 )
integer            LDA, LDZ
PARAMETER        ( LDA = N, LDZ = N )
integer            LWMAX
PARAMETER        ( LWMAX = 1000 )

!!      .. Local Scalars ..
integer            INFO, LWORK, LIWORK, IL, IU, M
real   ABSTOL, VL, VU

!     .. Local Arrays ..
integer            ISUPPZ( N ), IWORK( LWMAX )
real               W( N ), Z( LDZ, NSELECT), WORK( LWMAX )

!     .. External Subroutines ..
EXTERNAL         SSYEVR
!     .. Intrinsic Functions ..
INTRINSIC        INT, ABS, SQRT

!!! STEP TWO: now get coordinates of peptide in terms of (x,y,z) and Euler angles
!!! NOTE: I use radius of gyration axis and x,y,z axis of system to define the Euler angles

!calculate peptide center of mass
call calc_com_backbone(group_pep, pep_com, pep_res)

!calculate peptide's moment of inertia tensor 
call calc_mom_tensor_backbone(group_pep, mom_tensor, pep_com) 

!Use SSYEVR to find eigenvalues/eigenvectors. This is copied and pasted from Intel page here https://www.intel.com/content/www/us/en/docs/onemkl/code-samples-lapack/2023-1/ssyevr-example-fortran.html
ABSTOL = -1.0 !Negative ABSTOL means using the default value
!Set IL, IU to compute NSELECT smallest eigenvalues
IL = 1
IU = NSELECT

LWORK = -1
LIWORK = -1
!query the workspace
CALL SSYEVR( 'Vectors', 'All', 'Upper', N, mom_tensor, LDA, VL, VU, IL, &
            IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, &
            LIWORK, INFO )
LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
LIWORK = MIN( LWMAX, IWORK( 1 ) )

! Solve eigenproblem. Just do it.
CALL SSYEVR( 'Vectors', 'Indices', 'Upper', N, mom_tensor, LDA, VL, VU, IL, &
            IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, &
            LIWORK, INFO )

!Check for convergence.
IF( INFO.GT.0 ) THEN
    WRITE(*,*)'The algorithm failed to compute eigenvalues.'
    STOP
END IF

!store moments of inertia, from largest to smallest
count = 0
do ax=1,3
    obj_axes(:,4-ax) = Z(:,ax)
    moments(4-ax) = W(ax)
enddo

!check the eignvectors were properly returned - not sure why, but SSYEVR stopped returning them properly 
do ax=1,3
    call vec_norm(obj_axes(:,ax), mag)
    if (abs(mag-1) .gt. 0.001) then
        write(*,*) "Fixing a missing eigenvector!!"
        call cross_product(obj_axes(:,MOD(ax+1,3)), obj_axes(:,MOD(ax+2,3)), obj_axes(:,ax))
        count = count + 1
        if (count .eq. 2) then
            write(*,*) "Multiple eigenvectors returned by SSYEVR incorrect! Stopping program"
            stop
        endif
    endif
enddo

!Since the direction of the eigenvectors can be flipped, need constraints to fully determine peptide orientation 
!I pick the below checks so the axes align with physical features of the peptide 

!Check 1: ensure obj_axes define a right handed system
call cross_product(obj_axes(:,1), obj_axes(:,2), v)
if (SUM(ABS(v - obj_axes(:,3))) > 0.1 ) then
    !(*,*) "Enforcing right-handedness of system!!!"
    obj_axes(:,3) = -1.*obj_axes(:,3)
endif

!Check 2: Always pick the first eigenvector direction so it has positive dot-product with the vector from peptide center of mass to N-terminus
mag = DOT_PRODUCT(obj_axes(:,1), group_pep(1)%coo1(:,1) - pep_com)
if (mag .lt. 0) then !need to change sign of two axes to keep system right-handed
    !write(*,*) "Forcing dot product of first peptide axis with N-terminus to be postive!!!"
    obj_axes(:,1) = -1*obj_axes(:,1)
    obj_axes(:,3) = -1*obj_axes(:,3)
endif

!Check 3: Always pick the second eigenvector direction so it has positive dot-product with the vector from peptide center of mass to C-terminus
mag = DOT_PRODUCT(obj_axes(:,2), group_pep(pep_res)%coo3(:,1) - pep_com)
if (mag .lt. 0) then 
    !write(*,*) "Forcing dot product of second peptide axis with the y-axis of system!!!"
    obj_axes(:,2) = -1*obj_axes(:,2)
    obj_axes(:,3) = -1*obj_axes(:,3)
endif

!if (debug_flag.eq.1) then
!    !verify the eigenvalues/eigenvectors actually diagonalizes the moment of inertia tensor
!    D = 0
!    obj_axes_t = TRANSPOSE(obj_axes)
!    D(1,1) = moments(1); D(2,2) = moments(2); D(3,3) = moments(3);
!    mom_tensor_transform = MATMUL(obj_axes_t, MATMUL(D, obj_axes))
!    mom_tensor_transform = MATMUL(obj_axes, MATMUL(D, obj_axes_t))
!    mom_tensor_transform = MATMUL(obj_axes_t, MATMUL(mom_tensor_orig, obj_axes))
!    mom_tensor_transform = MATMUL(obj_axes, MATMUL(mom_tensor_orig, obj_axes_t))
!endif

return
end subroutine decompose_mom_tensor

subroutine calc_Euler_angles(group_pep, angles, obj_axes, moments)
implicit none

real                        :: angles(3), N_Eul(3), moments(3), obj_axes(3,3)
real, parameter             :: PI = 3.141592653589793238462643383
real                        :: z_axis(3) = (/0.,0.,1./), x_axis(3) = (/1.,0.,0./), v(3)
type(groupdetails)          :: group_pep(pep_res)

!debug variables
!real                        :: foo

!get the peptide axes from its moment of inertia priniciple components
call decompose_mom_tensor(group_pep, obj_axes, moments)

!get the N-vector
call cross_product(z_axis,obj_axes(:,3), N_Eul)
call normalize_vector(N_Eul)

!calculate the Euler angles using the z-x-z order, following the steps here https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Variational_Principles_in_Classical_Mechanics_(Cline)/13%3A_Rigid-body_Rotation/13.13%3A_Euler_Angles
angles(1) = ACOS(DOT_PRODUCT(x_axis, N_Eul)) 
angles(2) = ACOS(DOT_PRODUCT(z_axis, obj_axes(:,3)))
angles(3) = ACOS(DOT_PRODUCT(N_Eul, obj_axes(:,1)))

!shift angles if they are actually in 3rd or 4th quadrant (ACOS always gives result in range of 0 to pi)
if (N_Eul(2) < 0) angles(1) = -angles(1) !shift if angle is in 3rd or 4th quadrant
call cross_product(z_axis, obj_axes(:,3), v)
if (v(1) * N_Eul(1) < 0) angles(2) = -angles(2) !shift if angle is in 3rd or 4th quadrant; Logic is checking if cross product is antiparallel, in which case sin theta is negative so theta is in third or 4th quadrant
call cross_product(N_Eul, obj_axes(:,1), v)
if (v(1) * obj_axes(1,3) < 0) angles(3) = -angles(3) !shift if angle is in 3rd or 4th quadrant; Logic is checking if cross product is antiparallel, in which case sin theta is negative so theta is in third or 4th quadrant

!make sure angles(3) is in range of -pi/2 to pi/2 to avoid redundancy
!angles(2) = angles(2) -PI*INT(angles(2)/(PI/2))

!and now the Euler angles are calculated, hooray!
return
end subroutine calc_Euler_angles

subroutine calc_com(group, com, num_res)
implicit none
real                            :: com(3), net_mass
integer                         :: num_res, i, j
type(groupdetails)              :: group(:)

com = 0.0; net_mass = 0.0
do i=1, num_res
    do j=1, group(i)%cnum1
        com = com + group(i)%masses1(j) * group(i)%coo1(:,j)
        net_mass = net_mass + group(i)%masses1(j)
    enddo

    do j=1, group(i)%cnum2
        com = com + group(i)%masses2(j) * group(i)%coo2(:,j)
        net_mass = net_mass + group(i)%masses2(j)
    enddo

    do j=1, group(i)%cnum3
        com = com + group(i)%masses3(j) * group(i)%coo3(:,j)
        net_mass = net_mass + group(i)%masses3(j)
    enddo
enddo

com = com/net_mass

return
end subroutine calc_com

subroutine calc_com_backbone(group, com, num_res)
implicit none
real                            :: com(3), net_mass
integer                         :: num_res, i, j
type(groupdetails)              :: group(:)

com = 0.0; net_mass = 0.0
do i=1, num_res
    do j=1, group(i)%cnum1
        com = com + group(i)%masses1(j) * group(i)%coo1(:,j)
        net_mass = net_mass + group(i)%masses1(j)
    enddo

    do j=1, group(i)%cnum3
        com = com + group(i)%masses3(j) * group(i)%coo3(:,j)
        net_mass = net_mass + group(i)%masses3(j)
    enddo
enddo

com = com/net_mass

return
end subroutine calc_com_backbone

subroutine calc_rmsd(group_pep, original_pep, rmsd)
implicit none
integer                         :: i, j
real                            :: rmsd, norm
real                            :: N_cur(3), N_orig(3), CA_cur(3), CA_orig(3), C_cur(3), C_orig(3), delta(3)
type(groupdetails)              :: group_pep(pep_res), original_pep(pep_res)

rmsd=0.0
do i=1, pep_res
    do j=1, group_pep(i)%cnum1
        if(group_pep(i)%atype1(j)=="N") then
            N_cur=group_pep(i)%coo1(:,j)
        elseif(group_pep(i)%atype1(j)=="CA") then
            CA_cur=group_pep(i)%coo1(:,j)
        endif
    enddo
    do j=1, group_pep(i)%cnum3
        if(group_pep(i)%atype3(j)=="C") then
            C_cur=group_pep(i)%coo3(:,j)
        endif
    enddo

    do j=1, original_pep(i)%cnum1
        if(original_pep(i)%atype1(j)=="N") then
            N_orig=original_pep(i)%coo1(:,j)
        elseif(original_pep(i)%atype1(j)=="CA") then
            CA_orig=original_pep(i)%coo1(:,j)
        endif
    enddo
    do j=1, original_pep(i)%cnum3
        if(original_pep(i)%atype3(j)=="C") then
            C_orig=original_pep(i)%coo3(:,j)
        endif
    enddo
    rmsd = 0 

    delta = N_cur-N_orig
    call vec_norm(delta, norm)
    rmsd = rmsd+norm

    delta = CA_cur-CA_orig
    call vec_norm(delta, norm)
    rmsd = rmsd+norm

    delta = C_cur-C_orig
    call vec_norm(delta, norm)
    rmsd = rmsd+norm
enddo

rmsd=sqrt(rmsd/(3*pep_res))

return
end subroutine calc_rmsd

subroutine calc_mom_tensor(group, M, com)
implicit none

integer                     :: res, atom
real                        :: M(3,3), com(3)
type(groupdetails)          :: group(pep_res)

do res=1,pep_res
    do atom=1,group(res)%cnum1
        M(1,1) = M(1,1) + group(res)%masses1(atom) * ((group(res)%coo1(2,atom)-com(2))**2 + (group(res)%coo1(3,atom)-com(3))**2)
        M(2,2) = M(2,2) + group(res)%masses1(atom) * ((group(res)%coo1(1,atom)-com(1))**2 + (group(res)%coo1(3,atom)-com(3))**2) 
        M(3,3) = M(3,3) + group(res)%masses1(atom) * ((group(res)%coo1(1,atom)-com(1))**2 + (group(res)%coo1(2,atom)-com(2))**2) 

        M(1,2) = M(1,2) - group(res)%masses1(atom) * (group(res)%coo1(1,atom)-com(1)) * (group(res)%coo1(2,atom)-com(2))    
        M(1,3) = M(1,3) - group(res)%masses1(atom) * (group(res)%coo1(1,atom)-com(1)) * (group(res)%coo1(3,atom)-com(3)) 
        M(2,3) = M(2,3) - group(res)%masses1(atom) * (group(res)%coo1(2,atom)-com(2)) * (group(res)%coo1(3,atom)-com(3))
    enddo

    do atom=1,group(res)%cnum2
        M(1,1) = M(1,1) + group(res)%masses2(atom) * ((group(res)%coo2(2,atom)-com(2))**2 + (group(res)%coo2(3,atom)-com(3))**2)
        M(2,2) = M(2,2) + group(res)%masses2(atom) * ((group(res)%coo2(1,atom)-com(1))**2 + (group(res)%coo2(3,atom)-com(3))**2) 
        M(3,3) = M(3,3) + group(res)%masses2(atom) * ((group(res)%coo2(1,atom)-com(1))**2 + (group(res)%coo2(2,atom)-com(2))**2) 
        M(1,2) = M(1,2) + group(res)%masses2(atom) * (group(res)%coo2(1,atom)-com(1)) * (group(res)%coo2(2,atom)-com(2))    
        M(1,3) = M(1,3) + group(res)%masses2(atom) * (group(res)%coo2(1,atom)-com(1)) * (group(res)%coo2(3,atom)-com(3)) 
        M(2,3) = M(2,3) + group(res)%masses2(atom) * (group(res)%coo2(2,atom)-com(2)) * (group(res)%coo2(3,atom)-com(3))  
    enddo

    do atom=1,group(res)%cnum3
        M(1,1) = M(1,1) + group(res)%masses3(atom) * ((group(res)%coo3(2,atom)-com(2))**2 + (group(res)%coo3(3,atom)-com(3))**2)
        M(2,2) = M(2,2) + group(res)%masses3(atom) * ((group(res)%coo3(1,atom)-com(1))**2 + (group(res)%coo3(3,atom)-com(3))**2)  
        M(3,3) = M(3,3) + group(res)%masses3(atom) * ((group(res)%coo3(1,atom)-com(1))**2 + (group(res)%coo3(2,atom)-com(2))**2)  

        M(1,2) = M(1,2) - group(res)%masses3(atom) * (group(res)%coo3(1,atom)-com(1)) * (group(res)%coo3(2,atom)-com(2)) 
        M(1,3) = M(1,3) - group(res)%masses3(atom) * (group(res)%coo3(1,atom)-com(1)) * (group(res)%coo3(3,atom)-com(3))
        M(2,3) = M(2,3) - group(res)%masses3(atom) * (group(res)%coo3(2,atom)-com(2)) * (group(res)%coo3(3,atom)-com(3))  
    enddo
enddo

M(2,1) = M(1,2)
M(3,1) = M(1,3)
M(3,2) = M(2,3) 

M = M/net_mass_backbone

end subroutine calc_mom_tensor

subroutine calc_mom_tensor_backbone(group, M, com)
implicit none

integer                     :: res, atom
real                        :: M(3,3), com(3)
type(groupdetails)          :: group(pep_res)

do res=1,pep_res
    do atom=1,group(res)%cnum1
        M(1,1) = M(1,1) + group(res)%masses1(atom) * ((group(res)%coo1(2,atom)-com(2))**2 + (group(res)%coo1(3,atom)-com(3))**2)
        M(2,2) = M(2,2) + group(res)%masses1(atom) * ((group(res)%coo1(1,atom)-com(1))**2 + (group(res)%coo1(3,atom)-com(3))**2) 
        M(3,3) = M(3,3) + group(res)%masses1(atom) * ((group(res)%coo1(1,atom)-com(1))**2 + (group(res)%coo1(2,atom)-com(2))**2) 

        M(1,2) = M(1,2) - group(res)%masses1(atom) * (group(res)%coo1(1,atom)-com(1)) * (group(res)%coo1(2,atom)-com(2))    
        M(1,3) = M(1,3) - group(res)%masses1(atom) * (group(res)%coo1(1,atom)-com(1)) * (group(res)%coo1(3,atom)-com(3)) 
        M(2,3) = M(2,3) - group(res)%masses1(atom) * (group(res)%coo1(2,atom)-com(2)) * (group(res)%coo1(3,atom)-com(3))
    enddo

    !don't include sidechain atoms in this calculation, as I only want backbone properties!
    !do atom=1,group(res)%cnum2
    !    M(1,1) = M(1,1) + group(res)%masses2(atom) * ((group(res)%coo2(2,atom)-com(2))**2 + (group(res)%coo2(3,atom)-com(3))**2)
    !    M(2,2) = M(2,2) + group(res)%masses2(atom) * ((group(res)%coo2(1,atom)-com(1))**2 + (group(res)%coo2(3,atom)-com(3))**2) 
    !    M(3,3) = M(3,3) + group(res)%masses2(atom) * ((group(res)%coo2(1,atom)-com(1))**2 + (group(res)%coo2(2,atom)-com(2))**2) 
    !    M(1,2) = M(1,2) + group(res)%masses2(atom) * (group(res)%coo2(1,atom)-com(1)) * (group(res)%coo2(2,atom)-com(2))    
    !    M(1,3) = M(1,3) + group(res)%masses2(atom) * (group(res)%coo2(1,atom)-com(1)) * (group(res)%coo2(3,atom)-com(3)) 
    !    M(2,3) = M(2,3) + group(res)%masses2(atom) * (group(res)%coo2(2,atom)-com(2)) * (group(res)%coo2(3,atom)-com(3))  
    !enddo

    do atom=1,group(res)%cnum3
        M(1,1) = M(1,1) + group(res)%masses3(atom) * ((group(res)%coo3(2,atom)-com(2))**2 + (group(res)%coo3(3,atom)-com(3))**2)
        M(2,2) = M(2,2) + group(res)%masses3(atom) * ((group(res)%coo3(1,atom)-com(1))**2 + (group(res)%coo3(3,atom)-com(3))**2)  
        M(3,3) = M(3,3) + group(res)%masses3(atom) * ((group(res)%coo3(1,atom)-com(1))**2 + (group(res)%coo3(2,atom)-com(2))**2)  

        M(1,2) = M(1,2) - group(res)%masses3(atom) * (group(res)%coo3(1,atom)-com(1)) * (group(res)%coo3(2,atom)-com(2)) 
        M(1,3) = M(1,3) - group(res)%masses3(atom) * (group(res)%coo3(1,atom)-com(1)) * (group(res)%coo3(3,atom)-com(3))
        M(2,3) = M(2,3) - group(res)%masses3(atom) * (group(res)%coo3(2,atom)-com(2)) * (group(res)%coo3(3,atom)-com(3))  
    enddo
enddo

M(2,1) = M(1,2)
M(3,1) = M(1,3)
M(3,2) = M(2,3) 

M = M/net_mass_backbone

end subroutine calc_mom_tensor_backbone

subroutine shuffle_array_int(array)
implicit none
integer                 :: array(:), index, i, temp
real                    :: rand

do i=size(array), 1, -1
    call RANDOM_NUMBER(rand)
    index= 1+INT(rand*i)
    temp = array(i)
    array(i) = array(index)
    array(index) = temp
enddo
return

end subroutine shuffle_array_int

end module math