module global_variables
    use vec
	 !use ifport
    implicit none
    real*8,dimension(3) :: BOX,BOXI
    character(len=40) :: boxfile,merfile
    real*8 :: SCAL,avg_bead_distance,mer_width,radius_cont
    integer :: mers_per_chain,mpc_functional
    integer :: n_chains=0
    integer :: n_amorph_chains=0
    integer :: nchains_functional
    integer :: n_mc_steps,	mers_per_chain_cont
    integer :: n_cont_chains=0
    logical :: debug=.TRUE.
    logical :: rough_surface=.FALSE., fiber_mode=.false. , &
	          fiber_shell_mode=.false.
	real*8 , dimension(2) :: fiber_center
	real*8 :: fiber_radius,fiber_radius_in,fiber_radius_out
    !---------------------------------------------
    ! interaction penalties
    real*8 :: bead_dist_cutoff=1.0 ! 4 A for bead dist
contains

    function init_random_seed() result(iostat)
        integer :: iostat
        INTEGER :: i, clock,n
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed
		character(len=10) :: get_time
		character(len=8) :: get_date
		real*8 :: dble_time
		integer ,dimension(1) :: int_time
		call date_and_time(get_date,get_time)
		read(get_time,'(f10.3)') dble_time
		int_time=int(1000*dble_time)
        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))
        CALL SYSTEM_CLOCK(COUNT=clock)
        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        CALL RANDOM_SEED(PUT = seed)
        !CALL RANDOM_SEED(PUT = int_time(1:1))
        DEALLOCATE(seed)
        iostat=0
      !print*,' Seed ', seed(:)
    END function init_random_seed

    function init_random_seed_() result(iostat)
        implicit none
        integer :: iostat
        integer, allocatable :: seed(:)
        integer :: i, n, un, istat, dt(8), pid, t(2), s
        integer(8) :: count, tms

        call random_seed(size = n)
        allocate(seed(n))
        ! First try if the OS provides a random number generator
        open(unit=un, file="/dev/urandom", access="stream", &
&            form="unformatted", action="read", status="old", iostat=istat)
        if (istat == 0) then
            read(un) seed
            close(un)
        else
            ! Fallback to XOR:ing the current time and pid. The PID is
            ! useful in case one launches multiple instances of the same
            ! program in parallel.
            call system_clock(count)
            if (count /= 0) then
                t = transfer(count, t)
            else
                call date_and_time(values=dt)
                tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                    + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                    + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                    + dt(5) * 60 * 60 * 1000 &
                    + dt(6) * 60 * 1000 + dt(7) * 1000 &
                    + dt(8)
                t = transfer(tms, t)
            end if
            s = ieor(t(1), t(2))
           ! pid = getpid() + 1099279 ! Add a prime
            s = ieor(s, pid)
            if (n >= 3) then
                seed(1) = t(1) + 36269
                seed(2) = t(2) + 72551
                seed(3) = pid
                if (n > 3) then
                    seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
                end if
            else
                seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
            end if
        end if
        call random_seed(put=seed)
    end function init_random_seed_

    SUBROUTINE SVDCMP(A,M,N,MP,NP,W,V)
        integer :: NMAX=800
        integer :: np,M,N,MP
        real*8 , dimension(np,np) :: V
        real*8 , dimension(mp,np) :: A
        real*8 , dimension(800)  :: RV1
        real*8 , dimension(np)    :: W
        integer :: i ,L,NM,J,K,ITS
        real*8 :: X, Z,Y,S,SCALE,G,F,C,H,ANORM
        !      DIMENSION A(MP,NP),W(NP),V(NP,NP),RV1(NMAX)
        G=0.0
        SCALE=0.0
        ANORM=0.0
        DO 25 I=1,N
            L=I+1
            RV1(I)=SCALE*G
            G=0.0
            S=0.0
            SCALE=0.0
            IF (I.LE.M) THEN
                DO 11 K=I,M
                    SCALE=SCALE+ABS(A(K,I))
11              CONTINUE
                IF (SCALE.NE.0.0) THEN
                    DO 12 K=I,M
                        A(K,I)=A(K,I)/SCALE
                        S=S+A(K,I)*A(K,I)
12                  CONTINUE
                    F=A(I,I)
                    G=-SIGN(SQRT(S),F)
                    H=F*G-S
                    A(I,I)=F-G
                    IF (I.NE.N) THEN
                        DO 15 J=L,N
                            S=0.0
                            DO 13 K=I,M
                                S=S+A(K,I)*A(K,J)
13                          CONTINUE
                            F=S/H
                            DO 14 K=I,M
                                A(K,J)=A(K,J)+F*A(K,I)
14                          CONTINUE
15                      CONTINUE
                    ENDIF
                    DO 16 K= I,M
                        A(K,I)=SCALE*A(K,I)
16                  CONTINUE
                ENDIF
            ENDIF
            W(I)=SCALE *G
            G=0.0
            S=0.0
            SCALE=0.0
            IF ((I.LE.M).AND.(I.NE.N)) THEN
                DO 17 K=L,N
                    SCALE=SCALE+ABS(A(I,K))
17              CONTINUE
                IF (SCALE.NE.0.0) THEN
                    DO 18 K=L,N
                        A(I,K)=A(I,K)/SCALE
                        S=S+A(I,K)*A(I,K)
18                  CONTINUE
                    F=A(I,L)
                    G=-SIGN(SQRT(S),F)
                    H=F*G-S
                    A(I,L)=F-G
                    DO 19 K=L,N
                        RV1(K)=A(I,K)/H
19                  CONTINUE
                    IF (I.NE.M) THEN
                        DO 23 J=L,M
                            S=0.0
                            DO 21 K=L,N
                                S=S+A(J,K)*A(I,K)
21                          CONTINUE
                            DO 22 K=L,N
                                A(J,K)=A(J,K)+S*RV1(K)
22                          CONTINUE
23                      CONTINUE
                    ENDIF
                    DO 24 K=L,N
                        A(I,K)=SCALE*A(I,K)
24                  CONTINUE
                ENDIF
            ENDIF
            ANORM=MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
25      CONTINUE
        DO 32 I=N,1,-1
            IF (I.LT.N) THEN
                IF (G.NE.0.0) THEN
                    DO 26 J=L,N
                        V(J,I)=(A(I,J)/A(I,L))/G
26                  CONTINUE
                    DO 29 J=L,N
                        S=0.0
                        DO 27 K=L,N
                            S=S+A(I,K)*V(K,J)
27                      CONTINUE
                        DO 28 K=L,N
                            V(K,J)=V(K,J)+S*V(K,I)
28                      CONTINUE
29                  CONTINUE
                ENDIF
                DO 31 J=L,N
                    V(I,J)=0.0
                    V(J,I)=0.0
31              CONTINUE
            ENDIF
            V(I,I)=1.0
            G=RV1(I)
            L=I
32      CONTINUE
        DO 39 I=N,1,-1
            L=I+1
            G=W(I)
            IF (I.LT.N) THEN
                DO 33 J=L,N
                    A(I,J)=0.0
33              CONTINUE
            ENDIF
            IF (G.NE.0.0) THEN
                G=1.0/G
                IF (I.NE.N) THEN
                    DO 36 J=L,N
                        S=0.0
                        DO 34 K=L,M
                            S=S+A(K,I)*A(K,J)
34                      CONTINUE
                        F=(S/A(I,I))*G
                        DO 35 K=I,M
                            A(K,J)=A(K,J)+F*A(K,I)
35                      CONTINUE
36                  CONTINUE
                ENDIF
                DO 37 J=I,M
                    A(J,I)=A(J,I)*G
37              CONTINUE
            ELSE
                DO 38 J= I,M
                    A(J,I)=0.0
38              CONTINUE
            ENDIF
            A(I,I)=A(I,I)+1.0
39      CONTINUE
        DO 49 K=N,1,-1
            DO 48 ITS=1,30
                DO 41 L=K,1,-1
                    NM=L-1
                    IF ((ABS(RV1(L))+ANORM).EQ.ANORM)  GO TO 2
                    IF ((ABS(W(NM))+ANORM).EQ.ANORM)  GO TO 1
41              CONTINUE
1               C=0.0
                S=1.0
                DO 43 I=L,K
                    F=S*RV1(I)
                    IF ((ABS(F)+ANORM).NE.ANORM) THEN
                        G=W(I)
                        H=SQRT(F*F+G*G)
                        W(I)=H
                        H=1.0/H
                        C= (G*H)
                        S=-(F*H)
                        DO 42 J=1,M
                            Y=A(J,NM)
                            Z=A(J,I)
                            A(J,NM)=(Y*C)+(Z*S)
                            A(J,I)=-(Y*S)+(Z*C)
42                      CONTINUE
                    ENDIF
43              CONTINUE
2               Z=W(K)
                IF (L.EQ.K) THEN
                    IF (Z.LT.0.0) THEN
                        W(K)=-Z
                        DO 44 J=1,N
                            V(J,K)=-V(J,K)
44                      CONTINUE
                    ENDIF
                    GO TO 3
                ENDIF
                !          IF (ITS.EQ.30) PAUSE 'No convergence in 30 iterations'
                X=W(L)
                NM=K-1
                Y=W(NM)
                G=RV1(NM)
                H=RV1(K)
                F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y)
                G=SQRT(F*F+1.0)
                F=((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
                C=1.0
                S=1.0
                DO 47 J=L,NM
                    I=J+1
                    G=RV1(I)
                    Y=W(I)
                    H=S*G
                    G=C*G
                    Z=SQRT(F*F+H*H)
                    RV1(J)=Z
                    C=F/Z
                    S=H/Z
                    F= (X*C)+(G*S)
                    G=-(X*S)+(G*C)
                    H=Y*S
                    Y=Y*C
                    DO 45 NM=1,N
                        X=V(NM,J)
                        Z=V(NM,I)
                        V(NM,J)= (X*C)+(Z*S)
                        V(NM,I)=-(X*S)+(Z*C)
45                  CONTINUE
                    Z=SQRT(F*F+H*H)
                    W(J)=Z
                    IF (Z.NE.0.0) THEN
                        Z=1.0/Z
                        C=F*Z
                        S=H*Z
                    ENDIF
                    F= (C*G)+(S*Y)
                    X=-(S*G)+(C*Y)
                    DO 46 NM=1,M
                        Y=A(NM,J)
                        Z=A(NM,I)
                        A(NM,J)= (Y*C)+(Z*S)
                        A(NM,I)=-(Y*S)+(Z*C)
46                  CONTINUE
47              CONTINUE
                RV1(L)=0.0
                RV1(K)=F
                W(K)=X
48          CONTINUE
3       CONTINUE
49  CONTINUE
    RETURN
END subroutine


end module global_variables
