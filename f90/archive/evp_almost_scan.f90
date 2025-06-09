program evp_scan
! 2 4.0 5 2 2.0

    implicit none
    INTEGER:: MREAD,EIG_BAR_IND,M_1
    REAL:: AKREAD,AK_1
    
    call system('make swipe')
    
    DO MREAD = 7,9
        DO AKREAD = 4.D0,10.D0,1.D0
            DO EIG_BAR_IND = 4,5
                DO M_1 = 1,9
                    DO AK_1 = 1.D0,10.D0,1.D0
                        IF ((AKREAD.EQ.AK_1).AND.(MREAD.EQ.M_1)) CYCLE

                        WRITE(*,*) MREAD,AKREAD,EIG_BAR_IND,M_1,AK_1

                        ! IF ((MREAD.LE.1).AND.(AKREAD.LT.7.D0)) CYCLE
                        ! IF ((MREAD.EQ.1).AND.(AKREAD.EQ.7.D0)) THEN
                        !     IF (EIG_BAR_IND.LT.4) CYCLE
                        !     IF ((EIG_BAR_IND.EQ.4).AND.(M_1.LT.2)) CYCLE
                        !     IF ((EIG_BAR_IND.EQ.4).AND.(M_1.EQ.2).AND.(AK_1.LT.3.D0)) CYCLE
                        ! ENDIF

                        call change_mk(MREAD,AKREAD,EIG_BAR_IND,M_1,AK_1)
                        call system('mpirun -np 4 ./bin/evp_almost_discrete_exec ')
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDDO
    
    contains
    subroutine change_mk(m_in,ak_in,eig_in,m_in2,ak_in2)
        INTEGER,INTENT(IN):: m_in,eig_in,m_in2
        REAL,INTENT(IN):: ak_in,ak_in2
        INTEGER:: fid
        ! CHARACTER*72:: STR1, STR2
    
        OPEN(UNIT=fid,FILE='almostdegen_read.input',STATUS='UNKNOWN')
            
        WRITE(fid,'(i4)') m_in
        WRITE(fid,'(f10.3)') ak_in
        WRITE(fid,'(i4)') eig_in
        WRITE(fid,'(i4)') m_in2
        WRITE(fid,'(f10.3)') ak_in2
        ! WRITE(STR1,'(i4)') MREAD
        ! WRITE(STR2,'(f10.3)') AKREAD
        ! STR1 = REPLACE(STR1,' ','', every=.TRUE.)
        ! WRITE(fid,*) trim(STR1)
        ! WRITE(fid,*) trim(STR2)
        WRITE(fid,*) 'End'
        
        CLOSE(fid)
        end subroutine
    
    end program evp_scan