module polymer_type
    use global_variables
    use checkbox
    implicit none
    type monomer
        character(len=16) :: monomer_name
        integer,allocatable,dimension(:) :: deleted
        integer :: n_mer_atoms,headc,tailc
        integer ,allocatable, dimension(:) :: atype
        character(len=2) , allocatable, dimension(:) :: symbol
        real*8 ,allocatable, dimension(:,:) ::coord
        real*8 ,allocatable, dimension(:) ::mass
        real*8 , dimension(2,3) :: base_vecs
        real*8  :: mer_length
        real*8 , dimension(3) :: HEADH, TAILH
    contains 
        procedure , public :: create => read_monomer 
    end type
    type chain
        integer :: n_beads,chain_index,n_mers
        integer, allocatable , dimension(:) :: close_flag
        integer , allocatable , dimension(:) :: mer_flag
        logical  :: functional
        character*6 :: chaintype
        real*8 , allocatable , dimension(:,:) :: coord
        real*8 , dimension(3) :: head,tail
        type(monomer) , allocatable , dimension(:) :: mers
    contains
        procedure , public :: put_monomer
        procedure , public :: put_all_monomers
        procedure , public :: check_two_mers
        procedure , public :: print_beads_only
        procedure , public :: rearrange_mer
        procedure , public :: print_polymer_chain
        procedure , public :: map_atypes
    end type chain

    type(monomer) :: base_mer
    type(chain), dimension(:) , allocatable :: chains
contains



    function read_monomer(this,filename) result(iostat)
        class(monomer) :: this
        character(len=*) :: filename
        character(len=40) :: dump
        integer :: i,j,n_mer_atoms,iostat
        open(1234,file=filename)
        read(1234,*)this%monomer_name, n_mer_atoms
        this%n_mer_atoms=n_mer_atoms
        print*,n_mer_atoms
        allocate(this%atype(n_mer_atoms),&
&       this%symbol(n_mer_atoms), &
&       this%mass(n_mer_atoms), &
&       this%deleted(n_mer_atoms), &
&       this%coord(n_mer_atoms,3))
        this%deleted=0
       
        read(1234,*) dump,this%base_vecs(1,:)
        read(1234,*) dump,this%base_vecs(2,:)
        do i=1, n_mer_atoms
          read(1234,*) this%symbol(i),this%coord(i,:),this%atype(i)
        end do
        read(1234,*) dump, this%headc,this%headh(:)
        read(1234,*) dump, this%tailc,this%tailh(:)
        read(1234,*) dump, mer_width
        read(1234,*) dump, this%mer_length
        iostat=0
        close(1234)
    end function
    function map_atypes(this) result(iostat)
        class(chain) :: this
        integer :: iostat
        integer :: i, j
        do i=1,this%n_mers
            do j=1,this%mers(i)%n_mer_atoms
                if( this%mers(i)%symbol(j) .eq. 'C ') then
                                ! this%mers(i)%atype(j) = 1
                                 this%mers(i)%mass(j) = 12.0
				end if
                if( this%mers(i)%symbol(j) .eq. 'H ') then
                                ! this%mers(i)%atype(j) = 2
                                 this%mers(i)%mass(j) = 1.0
				end if
                if( this%mers(i)%symbol(j) .eq. 'O ') then
                                 !this%mers(i)%atype(j) = 3
                                 this%mers(i)%mass(j) = 16.0
				end if
                if( this%mers(i)%symbol(j) .eq. 'N ') then
                                 !this%mers(i)%atype(j) = 4
                                 this%mers(i)%mass(j) = 28.0
				end if
                if( this%mers(i)%symbol(j) .eq. 'S ') then
                                 !this%mers(i)%atype(j) = 5
                                 this%mers(i)%mass(j) = 22.0
				end if
                if( this%mers(i)%symbol(j) .eq. 'Br') then
                                 !this%mers(i)%atype(j) = 6
                                 this%mers(i)%mass(j) = 32.0
				end if
            end do
        end do
        iostat=0
    end function map_atypes



    function print_polymer_chain(this,filename,indx) result(iostat)
        class(chain) :: this
        character(len=*) filename
        integer :: i , j,cnt , k,indx,iostat
        real*8 ,dimension(3) :: vec
        real*8 :: dumy
        !    pause
        !   print*,this%n_mers
        cnt=indx
        if( indx == 0 ) then
            open( 133, file=filename)
            write(133,1001) 'CRYST1',box,90.0,90.0,90.0
        end if
        if( this%chaintype .eq. 'amorph' ) then
          vec=this%head(:)
          ! vec is the position of chain head cap hydrogen!
          cnt=cnt+1
                write(133,1002)'ATOM  ', cnt,'H '&
                    &         ,vec,'H '
        end if
        do i=1,this%n_mers
            !     print*,'i ',i,this%mers(i)%n_mer_atoms
            do j=1,this%mers(i)%n_mer_atoms
                !         print*,' j ',j
                cnt=cnt+1
                write(133,1002)'ATOM  ', cnt,this%mers(i)%symbol(j)&
                    &         ,this%mers(i)%coord(j,:),this%mers(i)%symbol(j)
            end do
        end do
         if( this%chaintype .eq. 'amorph' ) then

          vec=this%tail(:)
          ! vec is the position of chain tail cap hydrogen!
          cnt=cnt+1
                write(133,1002)'ATOM  ', cnt,'H '&
                    &         ,vec,'H '
        end if

        iostat=cnt

1001    format(A6,3(f9.3),x,3(f7.2,x))
1002    format(A6,i5,A4,15x,3(f8.3),23x,A2)
    end function print_polymer_chain

    function put_monomer(this,nth) result(iostat)
        use vec
		use checkbox
        class(chain) :: this
        integer:: nth,i,iostat,m,cnt
        real*8 , dimension(3) :: vecA,vecB,vecC,Ua,Ub,Uc,vectemp
        real*8 , dimension(3,3) :: VBS
		integer :: p
        !-----------------------------------------------
        ! puts ith monomer into chain
        !-----------------------------------------------
        if( nth > this%n_mers ) then
            iostat=1
            print*, ' Error at put monomer'
            return
        end if
        !-----------------------------------------------
        if( nth < this%n_mers ) then
        vecA=this%coord(nth+1,:) - this%coord(nth,:)
        else if( nth .eq. this%n_mers ) then
        vecA=this%coord(nth,:) - this%coord(nth-1,:)
        end if
		if( this%functional) then
		vecA=this%coord(2,:)-this%coord(1,:)
		!pause
		end if
		cnt=0
10      cnt=cnt+1
		vecB=random_unit_vec()
        vecC=random_unit_vec()
        call gsch(vecA,vecB,vecC,Ua,Ub,Uc)
        ! multiply with mer's base vector
        VBS(1,:)=Ua(:) !*base_mer%base_vecs(1,1)
        VBS(2,:)=Ub(:) !*base_mer%base_vecs(1,2)
        VBS(3,:)=Uc(:) ! *base_mer%base_vecs(1,3)
        ! print*,'U mag : ',vecmag(Ua),vecmag(Ub),vecmag(Uc)
		!write(*,'(A,x,9f6.2)') 'A : ', Ua,Ub,Uc
        this%mers(nth)=base_mer
		
        do i=1,base_mer%n_mer_atoms
            do m=1,3
                this%mers(nth)%coord(i,m)=base_mer%coord(i,1)*VBS(1,m)+&
                    &	   base_mer%coord(i,2)*VBS(2,m)+base_mer%coord(i,3)*VBS(3,m)
            end do
			vectemp=this%mers(nth)%coord(i,:)+this%coord(nth,:)
            p=insidevoids_buff(vectemp)
            if( p .eq. 0 .and. cnt < 100) then
			 ! print*,nth, ' atom inside box'
			  goto 10
			  end if
            this%mers(nth)%coord(i,:)=this%mers(nth)%coord(i,:)+this%coord(nth,:)
        end do
        if( this%chaintype .eq. 'amorph') then
            if( nth .eq. 1 ) then
            do m=1,3
                vectemp(m)=base_mer%headh(1)*VBS(1,m)+&
                    &      base_mer%headh(2)*VBS(2,m)+base_mer%headh(3)*VBS(3,m)
            end do
              !add head hydrogen to head c
              this%head(:)=this%mers(1)%coord(base_mer%headc,:)+vectemp
           !   print*,'Head H' , base_mer%headh(:)
           !   print*,'Head C ' ,sqrt(dot_product(vectemp,vectemp))
           end if
            if( nth .eq. this%n_mers ) then
            vectemp=0.d0
            !print*,'length tailh ', vecmag(base_mer%tailh)
            do m=1,3
                vectemp(m)=base_mer%tailh(1)*VBS(1,m)+&
                    &      base_mer%tailh(2)*VBS(2,m)+base_mer%tailh(3)*VBS(3,m)
            end do
              !add head hydrogen to head c
            !  print*,'Tail Vec : ', vectemp
            !  print*,' Distance to C ', sqrt(dot_product(vectemp,vectemp))
              this%tail(:)=this%coord(nth,:)+vectemp
              vectemp=this%tail(:)-this%mers(nth)%coord(base_mer%tailc,:)
            !  print*,'Second C : ', this%mers(nth)%coord(base_mer%tailc,:)
            !  print*,'Tail H ', this%tail(:)
            !  print*,'dist to second c ', vecmag(vectemp)
           end if


         end if


        iostat=0
        this%mer_flag(nth)=1
    end function

    !==================================================================
    function  check_two_mers(this,a,b) result( min_dist)
        class(chain) :: this
        integer :: a , b
        real*8 :: min_dist,dist
        integer :: i , j
        real*8 , dimension (3) :: vecA,vecB
        if(  (this%mer_flag(a) .ne. 1) .or. &
            &         (this%mer_flag(b) .ne. 1) ) then
            print*,' Error with check_two_mers '
            return
        end if
        min_dist=20.d0
        do i=1,this%mers(a)%n_mer_atoms
            vecA=this%mers(a)%coord(i,:)
            do j=1,this%mers(b)%n_mer_atoms
                vecB=this%mers(b)%coord(j,:)
                dist=pbc_dist(vecA,vecB)
                if( dist < min_dist) min_dist=dist
            !        print*,'i : ',i,vecA
            !        print*,' j ',j,vecB  !' dist ',dist
            !           write(*,'(A,i2,x,A,i2,x,f8.4)')' i : ', i ,&
            !&          ' j : ', j , dist
            end do
        end do
    end function

    !==================================================================
    function put_all_monomers(this) result(iostat)
        class(chain) :: this
        integer :: iostat
        integer :: i , j ,k,cnt,try
        logical :: flag=.FALSE.
        real*8  :: dist1,dist2,dist,dist3
        cnt=0
        try=0
9       try=try+1
        do i=1,this%n_mers
            ! print*,i
            iostat=this%put_monomer(i)
        end do
	!	pause 'put monomer'
10  continue
    flag=.FALSE.
    do i=1,this%n_mers
           !  print*,'----- ' , i
        dist=0
        !11    continue
        if ( i > 1) then
            dist1=this%check_two_mers(i,i-1)
        else
            dist1=10.d0
        end if
        if ( i < this%n_mers) then
            dist2=this%check_two_mers(i,i+1)
        else
            dist2=10.d0
        end if

        !dist3=pbc_dist(this%coord(2*(i-1)-1,:),this%coord(2*(i+1)-1,:))
        !print*,' i : ',i,dist1,dist2

        if ( min(dist1,dist2) < 1.2) then
            !       iostat=this%put_monomer(i-1)
            !       iostat=this%put_monomer(i)
            !       iostat=this%put_monomer(i+1)
            !       print*,'before ',dist1,dist2
            dist3=this%rearrange_mer(i,dist1,dist2)
            !  	 dist1=this%check_two_mers(i,i-1)
            !	 dist2=this%check_two_mers(i,i+1)

              !  print*,i,dist3
            if( dist3  < 1.5 ) flag=.TRUE.
              ! pause
        !       goto 11
        end if
    end do
    cnt=cnt+1
   ! print*,'Try : ',cnt
    !  pause
    if( cnt > 10 ) then
        cnt=0
         print*,' Return back ',try
        !      pause
        if( try  > 5 ) goto 13
        goto 9
    end if
    if( flag ) goto 10
13 continue
   !print*, 'Monomers put '


   end function put_all_monomers

   !==================================================================
   function rearrange_mer(this,nth_mer,dist1,dist2) result ( dist)
       class(chain) :: this
       integer :: nth_mer
       type(monomer) :: temp_monomer
       real*8 :: dist,dist1,dist2,min_dist
       integer ::i, iostat
       !    dist1=10.d0
       !    dist2=10.d0
       min_dist=min(dist1,dist2)
       temp_monomer=this%mers(nth_mer)
       do i=1,50
           iostat=this%put_monomer(nth_mer)
           if( nth_mer == 1 ) then
               dist2=this%check_two_mers(1,2)
               iostat=this%put_monomer(nth_mer+1)
           else if( nth_mer == this%n_mers) then
               dist1=this%check_two_mers(nth_mer-1,nth_mer)
               if( this%n_mers > 1) iostat=this%put_monomer(nth_mer-1)
           else
               iostat=this%put_monomer(nth_mer+1)
               dist1=this%check_two_mers(nth_mer-1,nth_mer)
               dist2=this%check_two_mers(nth_mer,nth_mer+1)
               !print*,'ra ',i,nth_mer,dist1,dist2
           end if
           if( min(dist1,dist2) >= min_dist )then
               !print*,'ra ',i,nth_mer,dist1,dist2
               min_dist=min(dist1,dist2)
               temp_monomer=this%mers(nth_mer)
           end if
		   !pause
       end do
       this%mers(nth_mer)=temp_monomer
       dist=min_dist
   end function rearrange_mer

   function debugger(de_chain) result(iostat)
       type(chain) :: de_chain
       integer :: i , j , k,cnt,cnt_atom,iostat
       real*8 :: dist,dist1,dist2
       logical :: flag=.FALSE.
       iostat=0
       write(*,*) 'Debugger '
       print*,'Number of Mers : ',de_chain%n_mers
       print*,'Bead Coordinates '
       do i=1,de_chain%n_beads
           print*,de_chain%coord(i,:)
       end do
       cnt_atom=0
       print*,'Monomer Coordinates'
       !	write(132,*)(de_chain%n_mers*2) !*base_mer%n_mer_atoms
       do j=1,de_chain%n_mers
           iostat=de_chain%put_monomer(j)
       end do
       cnt=0
       goto 15
10 continue
   flag=.FALSE.
   do i=2,de_chain%n_mers-1
       !    print*,'----- ' , i
       dist=0
11 continue
   dist1=de_chain%check_two_mers(i,i-1)
   dist2=de_chain%check_two_mers(i,i+1)
   print*,i,dist1,dist2

   if ( min(dist1,dist2) < 0.1) then
       iostat=de_chain%put_monomer(i-1)
       iostat=de_chain%put_monomer(i)
       iostat=de_chain%put_monomer(i+1)
       flag=.TRUE.
       goto 11
   end if
   end do
   cnt=cnt+1
   print*,'Try : ',cnt
   if( flag ) goto 10
15 continue
   !	print*,'Dist : ',i, dist
   !if( dist > 2.0 ) return
   !iostat=de_chain%print_beads_only('beads.pdb')
   do j=1,de_chain%n_mers
       do i=1,2 !base_mer%n_mer_atoms
           !print*,i
           cnt_atom=cnt_atom+1


       end do
   end do
   end function debugger

   function print_all_beads(chains,filename) result(iostat)
       type(chain),dimension(:) :: chains
       character(len=*) :: filename
       integer :: iostat
       integer :: i, j , k,indx
       iostat=chains(1)%print_beads_only(filename,'start')
       do i=2,n_chains
           indx=sum(chains(1:i-1)%n_beads)
           iostat=chains(i)%print_beads_only(filename,'mid',indx)
       end do
       write(132,*) 'TER'
   end function print_all_beads

   function print_beads_only(this,filename,chain_pos,indx) result(iostat)
       class(chain) :: this
       character(len=*) :: filename
       integer , optional :: indx
       integer :: cnt_atom,i,iostat
       character(len=*),optional :: chain_pos
       !-------------------------------------------------------------
       ! this function outputs only backbone atoms of the polymer
       ! in pdb format
       !------------------------------------------------------------
       cnt_atom=0
       if( present(indx) ) cnt_atom=indx

       if( present(chain_pos)) then
           if( chain_pos .eq. 'start') then
               open( 132, file=filename)
               write(132,1001) 'CRYST1',box,90.0,90.0,90.0
           end if
       end if
       do i=1,this%n_beads
           cnt_atom=cnt_atom+1
           write(132,1002)'ATOM  ', cnt_atom,'C ',this%coord(i,:),'C '
       end do
       if( present(chain_pos)) then
           if( chain_pos .eq. 'end') then
               write(132,*)'TER'
           end if
       end if
       iostat=0

       !   adjustl(de_chain%mers(j)%symbol(i)),&
       !	de_chain%mers(j)%coord(i,:),de_chain%mers(j)%symbol(i)


1001   format(A6,3(f9.3),x,3(f7.2,x))
1002   format(A6,i5,A4,15x,3(f8.3),23x,A2)
   end function print_beads_only

   end module polymer_type
