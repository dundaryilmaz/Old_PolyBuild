module vec
    implicit none
      INTEGER  :: CMAX
      PARAMETER ( CMAX = 3) !Maximum components
      LOGICAL :: stpwrt
      contains

!C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!C  Gram-Schmidt
!=========================================================
    function random_unit_vec() result (vector)
    real*8 ,dimension (3) :: vector
     call random_number(vector)
        vector=vector-0.5d0
        vector=vector/sqrt(dot_product(vector,vector))
    end function random_unit_vec
!=========================================================
      SUBROUTINE gsch(v1,v2,v3,u1,u2,u3)

     ! USE vec

      IMPLICIT none

      REAL*8 :: inrm1,np2,inp2,inrm2,&
&               np13,inp13,nrm13,np2A3,nrm23,inp23,&
&               np23,nrm3,inrm3,&
&               nrm1,nrm2,&
&                NEG1
      DOUBLE PRECISION, DIMENSION(CMAX) :: v1,v2,v3,&
&               v2on1,nv2on1,v2i,&
&               v3on1,v3on2,v3i,v3ii,nv3on1,nv3on2,&
&               u1,u2,u3
      DOUBLE PRECISION, DIMENSION(CMAX,CMAX) :: proj1,&
&                proj2,proj3,proj12,proj123                 !projection oporators

!C Set needed varuables
      NEG1 = -1.0

!C Calculate u1

!C  Find the norm of the first vector
      CALL vnorm(v1,nrm1)
!C Normalize the first vector to get the unit vecto u1
      inrm1 = 1/nrm1
      CALL vsclmlt(inrm1,v1,u1)

!C Calculate u2

!C Calculate the projection of v2 on u1 for v2i
      CALL vouter(u1,u1,proj1)  ! |u1><u1|
      v2on1 = MATMUL(proj1,v2)  ! |u1><u1|v2>
      CALL vsclmlt(NEG1,v2on1,nv2on1)
      CALL vadd(v2,nv2on1,v2i)  ! |v2i> = |v2> - |u1><u1|v2>
!C Normalize v2on1 to get u2
      CALL vnorm(v2i,nrm2)
      inrm2 = 1/nrm2
      CALL vsclmlt(inrm2,v2i,u2) ! |u2> = |v2i> / <v2i|v2i>

!C Calculate u3
      CALL vouter(u2,u2,proj2)  ! |u2><u2|
      v3on1 = MATMUL(proj1,v3)  ! |u1><u1|v3>
      v3on2 = MATMUL(proj2,v3)  ! |u2><u2|v3>
      CALL vsclmlt(NEG1,v3on1,nv3on1)
      CALL vsclmlt(NEG1,v3on2,nv3on2)
      CALL vadd(v3,nv3on1,v3i)    ! ||v3> - |u1><u1|v3>
      CALL vadd(v3i,nv3on2,v3ii)  ! ||v3> - |u1><u1|v3> - |u2><u2|v3>
!C Normalize v3ii to get u3
      CALL vnorm(v3ii,nrm3)
      inrm3 = 1/nrm3
      CALL vsclmlt(inrm3,v3ii,u3) ! |u3> = |v3i> / <v3i|v3i>
!C Calculate projection opporater for u3
      CALL vouter(u3,u3,proj3)

!C Add projection opperaters
      CALL MATADD(proj1,proj2,proj12)
      CALL MATADD(proj12,proj3,proj123)


!C Write output
      IF (stpwrt) THEN
        WRITE(6,*) v1(:)
        WRITE(6,*) v2(:)
        WRITE(6,*) v3(:)
        WRITE(*,*)''
        WRITE(*,*)'u1 ='
        WRITE(6,*)u1
        WRITE(*,*)''
        WRITE(6,*) '|u1><u1| ='
        WRITE(6,*) proj1
        WRITE(6,*) ''
        WRITE(6,*) '|u1><u1|v2> ='
        WRITE(6,*) v2on1
        WRITE(6,*) ''
        WRITE(6,*) '-|u1><u1|v2> ='
        WRITE(6,*) nv2on1
        WRITE(6,*) ''
        WRITE(6,*) '|v2>-|u1><u1|v2> ='
        WRITE(6,*) v2i
        WRITE(6,*) ''
        WRITE(*,*)'sqrt<v2i|v2i> ='
        WRITE(6,*) nrm2
        WRITE(*,*)''
        WRITE(6,*) ''
        WRITE(*,*)'u2 ='
        WRITE(6,*)u2
        WRITE(*,*)''
        WRITE(6,*) '|u2><u2| ='
        WRITE(6,*) proj2
        WRITE(6,*) ''
        WRITE(6,*) '|u3><u3| ='
        WRITE(6,*) proj3
        WRITE(6,*) ''
        WRITE(6,*) '|u1><u1|+|u2><u2|+|u3><u3| ='
        WRITE(6,*) proj123
        WRITE(6,*) ''
      ENDIF



      END SUBROUTINE gsch

!C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!C  normalizes a given vector

      SUBROUTINE vnorm(vec1,n1)

      !USE vec, ONLY:CMAX

      IMPLICIT none
      REAL*8 :: vec1(CMAX),n1
      REAL*8 :: vsq

      CALL viner(vec1,vec1,vsq) ! <v1|v1>

      n1 = SQRT(vsq)            ! sqrt(<v1|v1>)

      ENDSUBROUTINE vnorm

!C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!C  transpose

      SUBROUTINE vtrns(kt1,br1)

     ! USE vec, ONLY:CMAX

      IMPLICIT none
      REAL*8 :: kt1(CMAX),br1(1,CMAX)
      INTEGER :: i

      DO i=1,CMAX
        br1(1,i) = kt1(i)
      ENDDO

      END SUBROUTINE vtrns


!C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!C adds two vectors

      SUBROUTINE vadd(vec1,vec2,vec3)

    !  USE vec, ONLY:CMAX

      IMPLICIT none
      REAL*8 :: vec1(CMAX),vec2(CMAX),vec3(CMAX)
      INTEGER :: i

      DO i=1,CMAX
        vec3(i) = vec1(i) + vec2(i)
      ENDDO

      END SUBROUTINE vadd

!C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!C  scaler multiplication of a vector with a real number

      SUBROUTINE vsclmlt(a,vec1,vec2)

   !   USE vec, ONLY:CMAX

      IMPLICIT none
      REAL*8 :: vec1(CMAX),vec2(CMAX),a
      INTEGER :: i

      DO i=1,CMAX
        vec2(i) = a*vec1(i)
      ENDDO


      END SUBROUTINE vsclmlt
!C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!C   inner product of two vectors

      SUBROUTINE viner(kt1,kt2,np12)

  !    USE vec, ONLY:CMAX

      IMPLICIT none
      REAL*8 :: kt1(CMAX),kt2(CMAX),np12,br1(1,CMAX)
      INTEGER :: i

      CALL vtrns(kt1,br1)

      np12 = 0.0 ! set inner product to zero

      DO i=1,CMAX
        np12 = np12 + br1(1,i)*kt2(i)
      ENDDO

      END SUBROUTINE viner

!C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!C  Outer product

      SUBROUTINE vouter(kt1,kt2,up12)


  !    USE vec, ONLY:CMAX

      IMPLICIT none
      REAL*8 :: kt1(CMAX),kt2(CMAX),up12(CMAX,CMAX),&
&                br2(1,CMAX)
      INTEGER :: i,j

      CALL vtrns(kt2,br2)

      DO i=1,CMAX
        DO j=1,CMAX
          up12(i,j) = kt1(i)*br2(1,j)
        ENDDO
      ENDDO

      END SUBROUTINE vouter

!C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!C  matrix addition

      SUBROUTINE MATADD(mat1,mat2,mat12)

  !    USE vec

      IMPLICIT none
      DOUBLE PRECISION, DIMENSION(CMAX,CMAX) :: mat1,mat2,mat12
      INTEGER :: i,j

      DO i=1,CMAX
        DO j=1,CMAX
          mat12(i,j) = mat1(i,j) + mat2(i,j)
        ENDDO
      ENDDO

      END SUBROUTINE MATADD

end module vec
