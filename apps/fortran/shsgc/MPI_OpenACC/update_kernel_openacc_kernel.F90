!
! auto-generated by ops_fortran.py
!
MODULE UPDATE_KERNEL_MODULE
USE OPS_FORTRAN_DECLARATIONS
USE OPS_FORTRAN_RT_SUPPORT

USE OPS_CONSTANTS
USE ISO_C_BINDING

INTEGER(KIND=4) xdim1
#define OPS_ACC1(x) (x+1)
INTEGER(KIND=4) xdim2
#define OPS_ACC2(x) (x+1)
INTEGER(KIND=4) xdim3
#define OPS_ACC3(x) (x+1)

INTEGER(KIND=4) multi_d4
INTEGER(KIND=4) xdim4
#define OPS_ACC_MD4(d,x) ((x)*3+(d))

contains

!$ACC ROUTINE(update_kernel) SEQ
!user function
subroutine update_kernel(rho_new, rhou_new, rhoE_new, s)

  real (kind=8), DIMENSION(1) :: rho_new, rhou_new, rhoE_new
  real (kind=8), INTENT(in), DIMENSION(3) :: s

  rho_new(OPS_ACC1(0))  = rho_new(OPS_ACC1(0))  + s(OPS_ACC_MD4(1,0));
  rhou_new(OPS_ACC2(0)) = rhou_new(OPS_ACC2(0)) + s(OPS_ACC_MD4(2,0));
  rhoE_new(OPS_ACC3(0)) = rhoE_new(OPS_ACC3(0)) + s(OPS_ACC_MD4(3,0));

end subroutine

#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3

#undef OPS_ACC_MD4


subroutine update_kernel_wrap( &
& opsDat1Local, &
& opsDat2Local, &
& opsDat3Local, &
& opsDat4Local, &
& dat1_base, &
& dat2_base, &
& dat3_base, &
& dat4_base, &
& start, &
& end )
  IMPLICIT NONE
  real(8) :: opsDat1Local(*)
  real(8) :: opsDat2Local(*)
  real(8) :: opsDat3Local(*)
  real(8), INTENT(IN) :: opsDat4Local(*)
  integer :: dat1_base
  integer :: dat2_base
  integer :: dat3_base
  integer :: dat4_base
  integer(4) start(1)
  integer(4) end(1)
  integer n_x


  !$acc parallel deviceptr(opsDat1Local,opsDat2Local,opsDat3Local,opsDat4Local)  
  !$acc loop 
  DO n_x = 1, end(1)-start(1)+1
    call update_kernel( &
    & opsDat1Local(dat1_base+(n_x-1)*1), &
    & opsDat2Local(dat2_base+(n_x-1)*1), &
    & opsDat3Local(dat3_base+(n_x-1)*1), &
    & opsDat4Local(dat4_base+(n_x-1)*3) )
  END DO
  !$acc end parallel

end subroutine

!host subroutine
subroutine update_kernel_host( userSubroutine, block, dim, range, &
& opsArg1, &
& opsArg2, &
& opsArg3, &
& opsArg4)
  IMPLICIT NONE
  character(kind=c_char,len=*), INTENT(IN) :: userSubroutine
  type ( ops_block ), INTENT(IN) :: block
  integer(kind=4), INTENT(IN):: dim
  integer(kind=4)   , DIMENSION(dim), INTENT(IN) :: range
  real(kind=8) t1,t2,t3
  real(kind=4) transfer_total, transfer

  type ( ops_arg )  , INTENT(IN) :: opsArg1
  real(8), DIMENSION(:), POINTER :: opsDat1Local
  integer(kind=4) :: opsDat1Cardinality
  integer(kind=4), POINTER, DIMENSION(:)  :: dat1_size
  integer(kind=4) :: dat1_base

  type ( ops_arg )  , INTENT(IN) :: opsArg2
  real(8), DIMENSION(:), POINTER :: opsDat2Local
  integer(kind=4) :: opsDat2Cardinality
  integer(kind=4), POINTER, DIMENSION(:)  :: dat2_size
  integer(kind=4) :: dat2_base

  type ( ops_arg )  , INTENT(IN) :: opsArg3
  real(8), DIMENSION(:), POINTER :: opsDat3Local
  integer(kind=4) :: opsDat3Cardinality
  integer(kind=4), POINTER, DIMENSION(:)  :: dat3_size
  integer(kind=4) :: dat3_base

  type ( ops_arg )  , INTENT(IN) :: opsArg4
  real(8), DIMENSION(:), POINTER :: opsDat4Local
  integer(kind=4) :: opsDat4Cardinality
  integer(kind=4), POINTER, DIMENSION(:)  :: dat4_size
  integer(kind=4) :: dat4_base

  integer n_x
  integer start(1)
  integer end(1)
  integer(kind=4) :: n

  type ( ops_arg ) , DIMENSION(4) :: opsArgArray

  opsArgArray(1) = opsArg1
  opsArgArray(2) = opsArg2
  opsArgArray(3) = opsArg3
  opsArgArray(4) = opsArg4

  call setKernelTime(13,userSubroutine//char(0),0.0_8,0.0_8,0.0_4,0)
  call ops_timers_core(t1)

#ifdef OPS_MPI
  IF (getRange(block, start, end, range) < 0) THEN
    return
  ENDIF
#else
  DO n = 1, 1
    start(n) = range(2*n-1)
    end(n) = range(2*n);
  END DO
#endif

  call c_f_pointer(getDatSizeFromOpsArg(opsArg1),dat1_size,(/dim/))
  xdim1 = dat1_size(1)
  opsDat1Cardinality = opsArg1%dim * xdim1
  dat1_base = getDatBaseFromOpsArg1D(opsArg1,start,1)
  call c_f_pointer(opsArg1%data_d,opsDat1Local,(/opsDat1Cardinality/))

  call c_f_pointer(getDatSizeFromOpsArg(opsArg2),dat2_size,(/dim/))
  xdim2 = dat2_size(1)
  opsDat2Cardinality = opsArg2%dim * xdim2
  dat2_base = getDatBaseFromOpsArg1D(opsArg2,start,1)
  call c_f_pointer(opsArg2%data_d,opsDat2Local,(/opsDat2Cardinality/))

  call c_f_pointer(getDatSizeFromOpsArg(opsArg3),dat3_size,(/dim/))
  xdim3 = dat3_size(1)
  opsDat3Cardinality = opsArg3%dim * xdim3
  dat3_base = getDatBaseFromOpsArg1D(opsArg3,start,1)
  call c_f_pointer(opsArg3%data_d,opsDat3Local,(/opsDat3Cardinality/))

  call c_f_pointer(getDatSizeFromOpsArg(opsArg4),dat4_size,(/dim/))
  xdim4 = dat4_size(1)
  opsDat4Cardinality = opsArg4%dim * xdim4
  multi_d4 = getDatDimFromOpsArg(opsArg4) ! dimension of the dat
  dat4_base = getDatBaseFromOpsArg1D(opsArg4,start,multi_d4)
  call c_f_pointer(opsArg4%data_d,opsDat4Local,(/opsDat4Cardinality/))

  call ops_H_D_exchanges_device(opsArgArray,4)
  call ops_halo_exchanges(opsArgArray,4,range)
  call ops_H_D_exchanges_device(opsArgArray,4)

  call ops_timers_core(t2)

  call update_kernel_wrap( &
  & opsDat1Local, &
  & opsDat2Local, &
  & opsDat3Local, &
  & opsDat4Local, &
  & dat1_base, &
  & dat2_base, &
  & dat3_base, &
  & dat4_base, &
  & start, &
  & end )

  call ops_timers_core(t3)
  call ops_set_dirtybit_device(opsArgArray, 4)
  call ops_set_halo_dirtybit3(opsArg1,range)
  call ops_set_halo_dirtybit3(opsArg2,range)
  call ops_set_halo_dirtybit3(opsArg3,range)

  !Timing and data movement
  transfer_total = 0.0_4
  call ops_compute_transfer(1, start, end, opsArg1,transfer)
  transfer_total = transfer_total + transfer
  call ops_compute_transfer(1, start, end, opsArg2,transfer)
  transfer_total = transfer_total + transfer
  call ops_compute_transfer(1, start, end, opsArg3,transfer)
  transfer_total = transfer_total + transfer
  call ops_compute_transfer(1, start, end, opsArg4,transfer)
  transfer_total = transfer_total + transfer
  call setKernelTime(13,userSubroutine,t3-t2,t2-t1,transfer_total,1)
end subroutine
END MODULE
