program postproc_vort
! ======================================================================
! LINEARIZED IVP 
! ======================================================================
   USE OMP_LIB
   USE MPI
   USE MOD_MISC
   USE MOD_BANDMAT
   USE MOD_EIG
   USE MOD_FD
   USE MOD_LIN_LEGENDRE
   USE MOD_SCALAR3
   USE MOD_FFT 
   USE MOD_LAYOUT
   USE MOD_LEGOPS
   USE MOD_MARCH
   USE MOD_DIAGNOSTICS
   USE MOD_INIT
   USE MOD_EVP, only: SAVE_PERTURB
! -------------------------
implicit none
! -------------------------
integer:: iii,it,mm,kk,pp,II,JJ,FILESTATUS
real(p8),dimension(1,1):: status
real(P8):: time_start,time_end,time0,tmp
type(scalar):: psi,chi,psio,chio,dpsi,dchi
type(scalar):: rur,rup,uz
complex(P8),dimension(:),allocatable:: vec_1, vec_2, vec_bar
complex(P8),dimension(:,:,:),allocatable:: psi_glb, chi_glb
complex(P8),dimension(:,:,:),allocatable:: rur_glb, rup_glb, uz_glb
logical:: file_save = .TRUE.

! REFINEMENT
character(LEN=200):: FileName_R
real(p8), allocatable, dimension(:,:):: EMK, EMK0
integer:: converg_flag1, converg_flag2

CALL MPI_INIT_THREAD(MPI_THREAD_SERIALIZED,MPI_THREAD_MODE,IERR)
IF (MPI_THREAD_MODE.LT.MPI_THREAD_SERIALIZED) THEN
   WRITE(*,*) 'The threading support is lesser than that demanded.'
   CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
ENDIF
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_RANK, IERR)

CALL READCOM('NOECHO')
CALL READIN(5)
CALL LEGINIT()

if (MPI_RANK.eq.0) then
    WRITE(*,*) 'PROGRAM STARTED'
    CALL PRINT_REAL_TIME()
endif

allocate(EMK(NTCHOP,NXCHOPDIM))
! allocate(EMK0(NTCHOP,NXCHOPDIM))

files%n = 1
IF (MPI_RANK.EQ.0) THEN
    CALL SYSTEM('rm ./*.dat')
    CALL SYSTEM('rm ./output/*.output')
ENDIF
do while (file_save)

    call allocate(psi)
    call allocate(chi)

    ! load saved file
    call mload(TRIM(ADJUSTL(FILES%SAVEDIR))//files%psi(files%n),psi)
    call mload(TRIM(ADJUSTL(FILES%SAVEDIR))//files%chi(files%n),chi)

    ! stop the postprocessing once the last file is loaded
    if(files%n > files%ne) goto 999

    ! postprocessing
    ! CALL ALLOCATE(RUR)
    ! CALL ALLOCATE(RUP)
    ! CALL ALLOCATE(UZ)
    ! CALL PC2VOR(RUR,RUP,UZ) ! VORTICITY
    ! CALL TOFP(RUR)
    ! CALL TOFP(RUP)
    ! CALL TOFP(UZ)

    ! ALLOCATE(RUR_GLB(NDIMR,NDIMTH,NDIMX), RUP_GLB(NDIMR,NDIMTH,NDIMX), UZ_GLB(NDIMR,NDIMTH,NDIMX))
    ! CALL MASSEMBLE(RUR%E,RUR_GLB,2)
    ! CALL MASSEMBLE(RUP%E,RUP_GLB,2)
    ! CALL MASSEMBLE( UZ%E, UZ_GLB,2)

    ! CALL DEALLOCATE(RUR)
    ! CALL DEALLOCATE(RUP)
    ! CALL DEALLOCATE(UZ)

    ! IF (MPI_RANK.EQ.0) THEN
    !     RUR_GLB = RUR_GLB**2 + RUP_GLB**2
    !     DO II = 1,NTH

    !     ENDDO
    ! ENDIF
    ! DEALLOCATE(RUR_GLB,RUP_GLB,UZ_GLB)

    ! save energy
    EMK = ENERGY_SPEC_MODIFIED(psi,chi)
    ! IF (MPI_RANK.EQ.0) THEN   

    !     WRITE(*,*) 'FILE#', files%n,'/',files%ne
       
    !     ! save energy map
    !     ! OPEN(UNIT=667,FILE='especData_MK'//RTOA6(files%t(files%n))//'.dat',STATUS='UNKNOWN',&
    !     ! & IOSTAT=FILESTATUS, POSITION='APPEND')
    !     OPEN(UNIT=667,FILE=TRIM(ADJUSTL(FILES%SAVEDIR))//'especData_MK'//ITOA3(files%n)//'.output',STATUS='UNKNOWN',&
    !         & IOSTAT=FILESTATUS, POSITION='APPEND')
    !     ! WRITE(667, 788) tim%t, ((EMK(MM,KK),KK=1,NXCHOPDIM),MM=1,NTCHOP)
    !     DO MM = 1,NTCHOP
    !         WRITE(667, 788) (EMK(MM,KK),KK=1,NXCHOPDIM)
    !     ENDDO  
    !     CLOSE(667)     
    !     788 FORMAT(*(E23.15))

    ! endif
    IF (MPI_RANK.EQ.0) THEN   
        WRITE(*,*) 'FILE#', files%n,'/',files%ne
        ! save energy of 0,0
        OPEN(UNIT=667,FILE=TRIM(ADJUSTL(FILES%SAVEDIR))//'especData_00.output',STATUS='UNKNOWN',&
            & IOSTAT=FILESTATUS, POSITION='APPEND')
        WRITE(667, 788) real(files%n), EMK(1,1)
        CLOSE(667)     
        788 FORMAT(*(E23.15))
    endif

    ! ! save PFF space
    ! call allocate(psio); psio%e = psi%e; psio%ln = psi%ln; CALL RTRAN(psio,1)
    ! call allocate(chio); chio%e = chi%e; chio%ln = chi%ln; CALL RTRAN(chio,1)
    ! DO II = 1,SIZE(PSI%E,2)
    !     DO JJ = 1,SIZE(PSI%E,3)
    !         ! IF (((II+PSI%INTH .EQ. MONITOR_MK(2,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(2,2)+1)) &
    !         !     .OR. ((II+PSI%INTH .EQ. MONITOR_MK(3,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(3,2)+1))) THEN
    !         !     CALL save_mode(tim%t,psio%E(:,II,JJ),chio%E(:,II,JJ),II+PSI%INTH-1,JJ+PSI%INX-1)
    !         ! ENDIF
    !         IF ((II+PSI%INTH .EQ. MONITOR_MK(4,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(4,2)+1)) THEN
    !             CALL save_mode(real(files%n),psio%E(:,II,JJ),chio%E(:,II,JJ),II+PSI%INTH-1,JJ+PSI%INX-1)
    !         ENDIF
    !     ENDDO
    ! ENDDO
    ! CALL DEALLOCATE(psio); CALL DEALLOCATE(chio)

    ! Truncate modes
    DO II = 1,SIZE(PSI%E,2)
        DO JJ = 1,SIZE(PSI%E,3)
            IF ((II+PSI%INTH .EQ. MONITOR_MK(4,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(4,2)+1)) THEN
                ! psi%e(41:,II,JJ) = 0.d0
                ! chi%e(41:,II,JJ) = 0.d0
                CALL save_mode(real(files%n),psi%E(:,II,JJ),chi%E(:,II,JJ),II+PSI%INTH-1,JJ+PSI%INX-1)
            ENDIF
        ENDDO
    ENDDO
    EMK = ENERGY_SPEC_MODIFIED(psi,chi)
    if (MPI_RANK.EQ.0) then   
        ! save energy of 0,0
        OPEN(UNIT=667,FILE=TRIM(ADJUSTL(FILES%SAVEDIR))//'especData_00.output',STATUS='UNKNOWN',&
            & IOSTAT=FILESTATUS, POSITION='APPEND')
        WRITE(667, 788) real(files%n), EMK(1,1)
        CLOSE(667)     
    endif

    ! ! save PFF space
    ! call allocate(psio); psio%e = psi%e; psio%ln = psi%ln; CALL RTRAN(psio,1)
    ! call allocate(chio); chio%e = chi%e; chio%ln = chi%ln; CALL RTRAN(chio,1)
    ! DO II = 1,SIZE(PSI%E,2)
    !     DO JJ = 1,SIZE(PSI%E,3)
    !         ! IF (((II+PSI%INTH .EQ. MONITOR_MK(2,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(2,2)+1)) &
    !         !     .OR. ((II+PSI%INTH .EQ. MONITOR_MK(3,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(3,2)+1))) THEN
    !         !     CALL save_mode(tim%t,psio%E(:,II,JJ),chio%E(:,II,JJ),II+PSI%INTH-1,JJ+PSI%INX-1)
    !         ! ENDIF
    !         IF ((II+PSI%INTH .EQ. MONITOR_MK(4,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(4,2)+1)) THEN
    !             CALL save_mode(real(files%n),psio%E(:,II,JJ),chio%E(:,II,JJ),II+PSI%INTH-1,JJ+PSI%INX-1)
    !         ENDIF
    !     ENDDO
    ! ENDDO
    ! CALL DEALLOCATE(psio); CALL DEALLOCATE(chio)

    ! ! save FFF space
    ! DO II = 1,SIZE(PSI%E,2)
    !     DO JJ = 1,SIZE(PSI%E,3)
    !         ! IF (((II+PSI%INTH .EQ. MONITOR_MK(2,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(2,2)+1)) &
    !         !     .OR. ((II+PSI%INTH .EQ. MONITOR_MK(3,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(3,2)+1))) THEN
    !         !     CALL save_mode(tim%t,psio%E(:,II,JJ),chio%E(:,II,JJ),II+PSI%INTH-1,JJ+PSI%INX-1)
    !         ! ENDIF
    !         IF ((II+PSI%INTH .EQ. MONITOR_MK(4,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(4,2)+1)) THEN
    !             CALL save_mode(real(files%n),psi%E(:,II,JJ),chi%E(:,II,JJ),II+PSI%INTH-1,JJ+PSI%INX-1)
    !         ENDIF
    !     ENDDO
    ! ENDDO

    ! deallocate
    call deallocate(psi)
    call deallocate(chi)

    files%n = files%n + 1
    if(files%n > files%ne) file_save = .FALSE.

enddo

999 continue
time_end = mpi_wtime()

!> output
! !> 1. FFF space: perturbation file
! ALLOCATE(PSI_GLB(NRCHOPDIM,NTCHOPDIM,NXCHOPDIM),CHI_GLB(NRCHOPDIM,NTCHOPDIM,NXCHOPDIM))
! CALL MASSEMBLE(PSI%E,PSI_GLB,1); CALL MASSEMBLE(CHI%E,CHI_GLB,1)
! IF (MPI_RANK.EQ.0) THEN
!    ALLOCATE(vec_1(2*NRCHOPDIM),vec_2(2*NRCHOPDIM),vec_bar(2*NRCHOPDIM))
!    vec_1(1:NRCHOPDIM) = PSI_GLB(:,MONITOR_MK(2,1)+1,MONITOR_MK(2,2)+1)
!    vec_1(NRCHOPDIM+1:) = CHI_GLB(:,MONITOR_MK(2,1)+1,MONITOR_MK(2,2)+1)
!    vec_2(1:NRCHOPDIM) = PSI_GLB(:,MONITOR_MK(3,1)+1,MONITOR_MK(3,2)+1)
!    vec_2(NRCHOPDIM+1:) = CHI_GLB(:,MONITOR_MK(3,1)+1,MONITOR_MK(3,2)+1) 
!    vec_bar(1:NRCHOPDIM) = PSI_GLB(:,MONITOR_MK(1,1)+1,MONITOR_MK(1,2)+1)
!    vec_bar(NRCHOPDIM+1:) = CHI_GLB(:,MONITOR_MK(1,1)+1,MONITOR_MK(1,2)+1)
!    CALL SAVE_PERTURB('./converg/new_perturb_final.input', &
!       M(MONITOR_MK(2,1)), AK(MONITOR_MK(2,1),MONITOR_MK(2,2)),vec_1, CMPLX(eig1), &
!       M(MONITOR_MK(3,1)), AK(MONITOR_MK(3,1),MONITOR_MK(3,2)),vec_2, CMPLX(eig2), &
!       M(MONITOR_MK(1,1)), AK(MONITOR_MK(1,1),MONITOR_MK(1,2)),vec_bar, CMPLX(0.D0))
!    DEALLOCATE(vec_1,vec_2,vec_bar)
! ENDIF
! DEALLOCATE(PSI_GLB,CHI_GLB)

! !> 2. PFF space: PSI-CHI
! CALL RTRAN(psi,1); CALL RTRAN(chi,1)
! DO II = 1,SIZE(PSI%E,2)
!    DO JJ = 1,SIZE(PSI%E,3)
!       IF (((II+PSI%INTH .EQ. MONITOR_MK(2,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(2,2)+1)) &
!          .OR. ((II+PSI%INTH .EQ. MONITOR_MK(3,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(3,2)+1))) THEN
!          CALL save_mode(tim%t,psi%E(:,II,JJ),chi%E(:,II,JJ),II+PSI%INTH-1,JJ+PSI%INX-1)
!       ENDIF
!    ENDDO
! ENDDO
! CALL RTRAN(psi,-1); CALL RTRAN(chi,-1)

!> final printout
IF (MPI_RANK.eq.0) THEN
    print *,tim%n,' steps'
    WRITE(*,*) 'PROGRAM STARTED'
    CALL PRINT_REAL_TIME()
    WRITE(*,*) 'EXECUTION TIME: ',time_end-time_start,'seconds'
ENDIF
call MPI_BARRIER(MPI_COMM_IVP,IERR)
call MPI_FINALIZE(IERR)

! ======================================================================
contains
! ======================================================================

subroutine save_mode(t,psi_new,chi_new,m_i,k_i)
! ======================================================================
   complex(p8),DIMENSION(:):: psi_new,chi_new
   real:: t
   integer:: m_i,k_i
   integer:: nn

   open(UNIT=777,FILE='mode_track_'//ITOA3(m_i)//'_'//ITOA3(k_i)//'.dat',&
   &STATUS='UNKNOWN',ACTION='WRITE',ACCESS='APPEND')

   WRITE(777,320) t,psi_new(1:size(psi_new)),chi_new(1:size(psi_new))

   close(777)
320 FORMAT(F10.3,',',(S,E24.16E3,SP,E24.16E3,'i'),*(','S,E24.16E3,SP,E24.16E3,'i'))
! 320 FORMAT(F10.3,',',(S,E14.6E3,SP,E14.6E3,'i'),*(','S,E14.6E3,SP,E14.6E3,'i'))
end subroutine save_mode
! ======================================================================

end program postproc_vort
