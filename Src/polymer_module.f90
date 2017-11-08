module polymer_module
    use global_variables
    use checkbox
    use polymer_type
    use omp_lib
    implicit none
	!----
	!---
!---------------------------------------------
contains
    !----------------------------------------
    function pre_analize() result(iostat)
        integer :: iostat,n_beads
        real*8 :: volume,density,vol_per_bead
        iostat=0
        volume=box(1)*box(2)*box(3)
        write(*,'(A,x,f12.2)')'Volume of the sim box : ' ,volume
        volume=volume-n_boxes*box_lengths(1)*box_lengths(2)*box_lengths(3)
		if( n_boxes > 0 ) then
         write(*,'(A26,f12.2)')'Volume for Polymer Chains:',volume
		end if
        n_beads=2*n_chains*mers_per_chain
		if( n_beads > 0 ) then
        vol_per_bead=volume/dble(n_beads)
        write(*,'(A,f12.4)')  'Volume per bead       : ',vol_per_bead
        avg_bead_distance=vol_per_bead**(1./3.)
		end if
        write(*,'(A,f8.2)')   'Dist                  : ',avg_bead_distance
		continue
    end function pre_analize
    function create_singlemerchain(head_position,direction,nth_chain) result(new_chain)
       type(chain) :: new_chain
	   integer :: nth_chain
	   real*8,dimension(3) :: head_position,direction
	    new_chain%n_beads=2
		new_chain%n_mers=1
		new_chain%chaintype='single'
		allocate(new_chain%mers(1),new_chain%coord(2,3),new_chain%close_flag(2),&
&                new_chain%mer_flag(1))
        new_chain%functional=.TRUE.
		new_chain%coord(1,:)=head_position
		new_chain%coord(2,:)=direction+head_position
	end function

    function create_chain(n_mers,nth_chain,surf_type,head_pos,success) result (new_chain)
        type(chain) :: new_chain
        integer :: nth_chain
        integer :: n_mers,i,p,n_beads
        integer ::cnt_interchain,cnt_second_bead,cnt_head
        integer , optional :: surf_type
        real*8, dimension(3) , optional :: head_pos
        real*8 :: rescale,dist_prev
        real*8 ,dimension(3) :: vec,new_pos,vectemp
        logical,  optional :: success
             !n_mers=mers_per_chain
        cnt_interchain=0
        cnt_head=0
        cnt_second_bead=0
        if(present(success) ) then
            success=.FALSE.
        end if
1000    format( i3,' th mer of chain ',i3)
             !-----------------------------------------------
             ! in previous version we were using two beads
             ! for each mer
             !-----------------------------------------------
        new_chain%n_beads=n_mers
        new_chain%n_mers=n_mers
        new_chain%chaintype='amorph'
        n_beads=n_mers
        new_chain%chain_index=nth_chain
        allocate(new_chain%mers(n_mers),new_chain%coord(n_beads,3),&
&       new_chain%close_flag(n_beads),new_chain%mer_flag(n_mers))
        new_chain%functional=.FALSE.
        if ( present(head_pos) ) then
            new_chain%functional=.TRUE.
        end if
!        if( .not. new_chain%functional)&
!&      	    print*,'Create Chain ',nth_chain
        new_chain%close_flag(:)=0
             !============================================
             ! Chain Head Creation
             !============================================
90      vec=random_pos_box()
        !vec=(/60.08,18.24,71.78/)
        p=insidevoids_buff(vec)
        !write(*,'(A,3f8.2,A,i2)')' Vec: ',vec,' Result : ',p
		!pause
        !write(*,'(A,3f8.2)')' Trying chain head : ',vec
        if(new_chain%functional ) then
            p=1
            vec=head_pos
        end if
        if( p .eq. 0 .and. new_chain%functional ) then
            print*,'p : ', p
            print*,' Trying to create functional group wrong pos'
            stop
        end if
        cnt_head=cnt_head+1
        new_chain%coord(1,1:3)=vec
        rescale=scaler(cnt_head)
        if( p .eq. 0 ) goto  90
        if( new_chain%functional) goto 91
        if( check_distance_with_allchains(new_chain,1,rescale)) goto 90
91  continue
  	    !print*,'Chain Head Started '
    	!write(*,1000)1,nth_chain
      ! write( *, *) 'Chain Head : ',vec
    if( debug ) then
      !pause
    end if
	if( n_mers .eq. 1) return
95  cnt_second_bead=cnt_second_bead+1
    vec=random_unit_vec()
    if( new_chain%functional) then
        if( surf_type .eq. 1 ) then ! upper xy plane
            vec=(/0.d0,0.d0,+1.d0/)
        else if( surf_type .eq. 2) then !lower xy plane
            vec=(/0.d0,0.d0,-1.d0/)
        end if
    end if
     ! scale with mer length
    vec=vec*base_mer%mer_length
    new_pos=vec+new_chain%coord(1,:)
    !check for box
    p=insidevoids_buff(new_pos)
    if( new_chain%functional ) p=1
    !new_pos=return_box(new_pos)
    new_chain%coord(2,:)=new_pos
!	write(*,'(A,3f12.2)')'Second Bead Position:',new_pos
    if( new_chain%functional) goto 96
    rescale=scaler(cnt_second_bead)
    if(  ( p .eq. 0 ) .or. &
&          (check_distance_with_allchains(new_chain,2,rescale))) then
        if( cnt_second_bead .gt. 10000) then
            cnt_second_bead=0
            cnt_head=0
            print*,'Trying again '
            goto 90
        end if
        goto 95
    end if
96 continue
   !write(*,1000)2,nth_chain
   if( debug) then
   ! write(*,*) ' Second bead : ', new_pos
    !print*, pbc_dist(new_pos,new_chain%coord(1,:))
   !pause
   end if
   do i=3,n_beads
       ! check wheater it is inside voids
       ! get a random unit vector
100    vec=random_unit_vec()
       cnt_interchain=cnt_interchain+1
       ! scale with mer length
       vec=vec*base_mer%mer_length
       new_pos=vec+new_chain%coord(i-1,:)
	   vectemp=new_pos
       !chain angle
       dist_prev= pbc_dist(new_pos,new_chain%coord(i-2,:)) 
       if(dist_prev < 1.53d0* base_mer%mer_length )then
           if( cnt_interchain > 10000) then
               cnt_second_bead=0
               cnt_interchain=0
               if( new_chain%functional ) then
                   success=.FALSE.
                   return
               end if
               goto 90
           end if
           goto 100
       end if
       !check for box
       p=insidevoids_buff(new_pos)
       new_chain%coord(i,:)=new_pos
       if( p .eq. 0 ) goto 100
       rescale=scaler(cnt_interchain)

       if( check_distance_within_chain(new_chain,i)) goto 100
       if( check_distance_with_allchains(new_chain,i,rescale)) goto 100
   end do

   if( present(success) ) success=.TRUE.

   end function create_chain
!------------------------------------------------------------------------
   function parallel_chain_create(s_c,n_chains,n_mers) result(iostat)
   use omp_lib
       integer :: n_chains,iostat,s_c,thread_id
       integer :: n_mers,i,n_threads,j,ch_id,id
       real*8  :: time_i,wtime
       real*8 , dimension(3) :: vec
       logical :: success
	   iostat=0
       print*,size(chains),'number of chains about to create : ', n_chains
       !n_threads=4

        do i=s_c,n_chains-1
          allocate(chains(i)%mers(n_mers),chains(i)%coord(n_mers,3))
        end do


  wtime = omp_get_wtime ( )

  do i=1,n_mers
!!$omp parallel shared(chains)
!!$omp do private(j)
 101 continue
     do j=s_c,s_c+n_chains
       vec=find_next_mer_pos(j,i,success)
       if(success) chains(j)%coord(i,:)=vec
     end do





!!$omp end do
!!$omp end parallel
end do
  wtime = omp_get_wtime ( )-wtime
  print*,'elapsed time ', wtime



   end function parallel_chain_create

   function find_next_mer_pos(ch_id,nth_mer,success) result(mer_pos)
     real*8 , dimension(3) :: mer_pos
     integer :: ch_id,id
     integer :: nth_mer,thread_id
     logical :: success
     real*8 , dimension(3) :: vec
101  continue !'try again'
     if( nth_mer .eq. 1 ) then
             mer_pos=random_pos_box()
         else
           vec=random_unit_vec()*base_mer%mer_length
           mer_pos=chains(ch_id)%coord(nth_mer-1,:)+vec
      end if
        if(check_distance_within_chain(chains(ch_id),nth_mer)) goto 101
      success=.TRUE.
   end function find_next_mer_pos
   function check_distance_with_allchains_parallel(ch_id,nth_mer,rescale) result(folded)
     integer :: ch_id, nth_mer
     real*8 :: rescale
     real*8 , dimension(3) :: vecA
     logical :: folded
     folded = .true.
   end function


   !==================================================================================
   function check_distance_with_allchains(checkedchain,nthbead,rescale) result(folded)
       type(chain) :: checkedchain
       integer :: nthbead,i,j
       logical :: folded
       integer :: nth_chain
       real*8 :: rescale
       real*8, dimension(3) :: veca,vecb
       real*8  :: distance
       !--------------------
       folded=.FALSE.

       nth_chain=checkedchain%chain_index
       !print*,'chain index ',nth_chain
       if( checkedchain%chain_index .eq. 1 ) then
           return
       end if
       veca=checkedchain%coord(nthbead,:)
       do i=1,nth_chain-1
           do j=1,chains(i)%n_beads
               vecb=chains(i)%coord(j,:)
               distance=pbc_dist(veca,vecb)

               !	 write(*,'(A,x,i4,A,x,i3,x,A,x,i2,x,A,x,i3,x,A,4f8.2)')&
               !&	 'Check chain : ',nth_chain,'Bead : ',nthbead,&
               !&     'To Chain ',i,' To bead ',j, ' vec : ',vecb,distance
!               if ( distance < rescale*1.5d0*base_mer%mer_length ) then
               if ( distance < rescale*mer_width ) then ! diger zincirlere olan uzakligi mer length
			   									   ! e baglamak dogru degil
                   !	  print*,'dist ',distance,rescale,rescale*2.d0*base_mer%mer_length
                   !	  pause
                   folded=.TRUE.
                   return
               end if
           end do
       end do

   end function check_distance_with_allchains
   !=============================================================

   function check_distance_within_chain(checkedchain,nthbead) result(folded)
       type(chain) :: checkedchain
       integer :: nthbead
       logical :: folded
       integer :: i
       real*8 , dimension(3) :: vec1,vec2
       real*8 :: dist
       !---------------------------------------------------------------
       ! checks the distances between nthbead and the previous beads
       ! if closer then cutoff folded=.TRUE.
       !---------------------------------------------------------------
       folded=.FALSE.
       if(nthbead .le. 2) return
       vec1=checkedchain%coord(nthbead,:)
       do  i = 1, nthbead-2
           vec2=checkedchain%coord(i,:)
           dist=pbc_dist(vec1,vec2)
           if( i == nthbead -2 ) then
               if( dist < 1.4*base_mer%mer_length ) then !cos(15)=1.93
                   folded=.TRUE.
                         !print*,' dist 1 ', dist
                   return
               end if
           else
               if( dist < mer_width ) then !cos(15)§
                   folded=.TRUE.
                     !print*,' dist 2 ' , i,dist
               end if
           end if
       end do
   end function check_distance_within_chain

   function scaler(count) result(scale)
       integer :: count
       real*8 :: scale
       if( count > 5000) then
           scale=0.5+exp((5000.d0-dble(count))/20000.d0)/2.d0
       else
           scale=1.d0
       end if
   end function scaler
   !===:========================================================================
   function readmer(merfile) result(iostat)
       character(len=*) :: merfile
       integer :: iostat
       real*8 , dimension(3) :: vec1,vec2
       integer :: n_atoms,i
       character(len=40) :: MTITLE,dumy

       write(*,*) 'READING BASE MONOMER DATA'
       OPEN(UNIT=11,FILE=merfile,STATUS='old')
       READ(11,*) MTITLE
       READ(11,*) vec1 !MBV(1,:)
       READ(11,*) vec2 !MBV(2,:)
       READ(11,*) ! rhogcm
       READ(11,*,err=99,end=99) SCAL ,rough_surface,n_mc_steps !,BLENG
       READ(11,*) BOX(:)
       READ(11,*) n_chains,mers_per_chain,nchains_functional,&
&                mpc_functional       !CAMORPH,MPC !/1000 scale up box
       READ(11,*)    !CTHERM,LTHERM,CONTBUFFZ
       READ(11,*) !ADDWALL
       READ(11,*)
       READ(11,*) n_atoms
       allocate(base_mer%atype(n_atoms),base_mer%coord(n_atoms,3),&
           &           base_mer%symbol(n_atoms) )
       do i=1,n_atoms
           read(11,*) base_mer%atype(i),base_mer%coord(i,:)
           !if( base_mer%atype(i) == 1 ) base_mer%symbol(i)='H '
           !if( base_mer%atype(i) == 6 ) base_mer%symbol(i)='C '
           !if( base_mer%atype(i) == 8 ) base_mer%symbol(i)='O '
       end do
      !   read(11,*) dumy,base_mer%HEADH(:)
      !   read(11,*) dumy,base_mer%TAILH(:)
       base_mer%n_mer_atoms=n_atoms
       base_mer%base_vecs(1,:)=vec1
       base_mer%base_vecs(2,:)=vec2
       base_mer%mer_length=SCAL
       !write(*,*) ' BASE MONOMER DATA READ'
       iostat=1
       return
99     stop 'Format change in the input file : SCAL, rough_surf_flag, n_mc_Steps'
   end function readmer
   function build_continuous_chain(n_mers,nth_chain,chainaxis) result(new_chain)
       type(chain) :: new_chain
       character*1 :: chainaxis
       character*2 :: walltype
       integer :: n_mers,i
       integer :: nth_chain,n_beads
       real*8  :: step_length,rescale
       real*8  :: corr_step_length,margin
       real*8  , dimension(3) :: vec,new_pos
       integer :: p
       logical :: success
       !     Here is the algorithm:
       !     Idea is to create a continious chain which exits the simulation box at +x plane and re-enters
       !     the box at -x plane. It will be continious along the x axis or any given axis.
       !     How to do this?
       !     Pick a position at -x wall. Divide length of the simulation box along x direction to number of
       !     mers in the chain. Lets call this step length. If step_length is less than mer length then quit.
       !     Normally it should be a fraction of the mer length.
       !     To find the second point pick two random numbers and normalize and add step length to the
       !     x axis. End of the day final mer will be at the +x wall.
       !     Find the vector to the head of the chain.
       !     Divide this vector to number of mers this will be correction step length : corr_step_length
       !     add this corr_step_length all mers step by step. End of the chain will move to the desired point.
       !
       integer :: cnt_head=0
       integer :: cnt_interchain=0
       integer :: cnt_second_bead=0
       !==========================
       new_chain%n_beads=n_mers
       new_chain%n_mers=n_mers
       new_chain%chaintype='conts'
       n_beads=n_mers
       new_chain%chain_index=nth_chain
       allocate(new_chain%mers(n_mers),new_chain%coord(n_beads,3),&
           &       new_chain%close_flag(n_beads),new_chain%mer_flag(n_mers))
       new_chain%functional=.FALSE.
       select case (chainaxis)
           case('x')
               step_length=(BOX(1)-base_mer%mer_length)/dble(n_mers)
               walltype='+x'
           case('y')
               step_length=(BOX(2)-base_mer%mer_length)/dble(n_mers)
               walltype='+y'
           case('z')
               !step_length=(BOX(3)-base_mer%mer_length)/dble(n_mers)
               step_length=(BOX(3))/dble(n_mers)
               walltype='+z'
       end select
       !------------------------------------------------------------
       if( step_length > base_mer%mer_length )then
           print*,step_length,base_mer%mer_length
           stop 'Insufficient number of mers to create cont. chain'
       end if
       !------------------------------------------------------------
       !print*,' Walltype : ', walltype, chainaxis,step_length
       !============================================
       ! Chain Head Creation
       !============================================
90     vec=random_pos_box()
       if(walltype .eq. '+x' ) vec(1)=0.d0
       if(walltype .eq. '+y' ) vec(2)=0.d0
       if(walltype .eq. '+z' ) vec(3)=0.d0
       !margin=20.d0
       p=insidevoids_buff(vec)

       !print*,' p  ',p
       !write(*,'(A,3f8.2)')' Trying chain head : ',vec
       cnt_head=cnt_head+1
       new_chain%coord(1,1:3)=vec
       rescale=scaler(cnt_head)
       if( p .eq. 0 ) goto  90
	   if (sqrt ( (vec(1)-Box(1)/2.d0)**2+(vec(2)-box(2)/2)**2.d0) .gt. radius_cont) goto 90
	   !if( vec(1)+margin .gt. BOX(1)) goto 90
	   !if( vec(2)+margin .gt. BOX(2)) goto 90
	   !if( vec(1)-margin .lt. 0) goto 90
	   !if( vec(2)-margin .lt. 0) goto 90

       if( check_distance_with_allchains(new_chain,1,rescale)) goto 90
!       print*,'Chain Head Started '
!       write(*,1000) nth_chain
!       write(*,'(A,3f8.4)') 'First bead : ', vec

95     cnt_second_bead=cnt_second_bead+1
       vec=random_unit_vec()
 !      if( chainaxis .eq. 'x')
       vec=step_adder(vec,base_mer%mer_length,step_length,chainaxis)

    ! scale with mer length
!   vec=vec*base_mer%mer_length
   new_pos=vec+new_chain%coord(1,:)
!   new_pos=return_box(new_pos)

   !check for box
   p=insidevoids_buff(new_pos)
   if( new_chain%functional ) p=1
   new_pos=return_box(new_pos)
   new_chain%coord(2,:)=new_pos
   if( new_chain%functional) goto 96
   rescale=scaler(cnt_second_bead)
   if(  ( p .eq. 0 ) .or. &
       &          (check_distance_with_allchains(new_chain,2,rescale))) then
       if( cnt_second_bead .gt. 100000) then
           cnt_second_bead=0
           cnt_head=0
           print*,'Trying again '
           goto 90
       end if
       goto 95
   end if
96 continue
   !print*, '2nd bead pos : ', new_pos
   !write(*,1000)2,nth_chain
   if( debug) then
   ! write(*,*) ' Second bead : ', new_pos
    !print*, pbc_dist(new_pos,new_chain%coord(1,:))
   !pause
   end if
   do i=3,n_beads
       ! check whether it is inside voids
       ! get a random unit vector
100    vec=random_unit_vec()
       vec=step_adder(vec,base_mer%mer_length,step_length,chainaxis)
       cnt_interchain=cnt_interchain+1
       ! scale with mer length
       !vec=vec*base_mer%mer_length
       new_pos=vec+new_chain%coord(i-1,:)

       !chain angle
       if( pbc_dist(new_pos,new_chain%coord(i-2,:)) < 1.73d0* base_mer%mer_length )then
!       if( pbc_dist(new_pos,new_chain%coord(i-2,:)) < 1.53d0* base_mer%mer_length )then
           if( cnt_interchain > 10000) then
               print*,'Trying beads'
               cnt_second_bead=0
               cnt_interchain=0
               if( new_chain%functional ) then
                   success=.FALSE.
                   return
               end if
               goto 90
           end if
           goto 100
       end if
       !check for box
       p=insidevoids_buff(new_pos)
!       new_pos=return_box(new_pos)
       new_chain%coord(i,:)=new_pos
       !print*,cnt_interchain
       if( p .eq. 0 ) goto 100
       !pause '6'
       !        if( cnt_interchain > 1000) then
       !          rescale=0.75+exp((1000.d0-dble(cnt_interchain))/10000.d0)/4.d0
       ! print*,cnt_interchain,rescale,nth_chain,i
       !stop
       !         end if
       rescale=scaler(cnt_interchain)

       if( check_distance_within_chain(new_chain,i)) goto 100
       if( check_distance_with_allchains(new_chain,i,rescale)) goto 100
       !write(*,1000)i,n_beads
       !pause 'end'
    !   write(*,'(A,3f8.2)') 'Next Bead Pos : ',new_pos
   end do
   ! vec=new_chain%coord(n_beads,:)-new_chain%coord(1,:)
   ! write(*,'(A,3f10.4)') 'Vector Head to Tail : ' , vec
   ! vec=pbc_x(vec)
   ! write(*,'(A,3f10.4)') 'Vector Head to Tail : ' , vec
    new_chain=join_head_to_tail(new_chain,chainaxis)
   do i=3,n_beads
     new_chain%coord(i,:)=return_box(new_chain%coord(i,:))
	 !write(*,'(3f9.2)') new_chain%coord(i,:)
   end do
   ! vec=new_chain%coord(n_beads,:)-new_chain%coord(1,:)
   ! write(*,'(A,3f10.4)') 'Vector Head to Tail : ' , vec
   ! vec=pbc_x(vec)
   ! write(*,'(A,3f10.4)') 'Vector Head to Tail : ' , vec


1000 format('Chain No ',i5, ' started.')

   end function build_continuous_chain
   function join_head_to_tail(new_chain,chainaxis) result(joined_chain)
     type(chain) :: new_chain
     type(chain) :: joined_chain
     character*1 :: chainaxis
     real*8 , dimension(3) :: rheadtail,delr
     integer :: i , j , k
     integer :: n_beads
       joined_chain=new_chain
       n_beads=new_chain%n_beads
       rheadtail=new_chain%coord(n_beads,:)-new_chain%coord(1,:)
       rheadtail=pbc_x(rheadtail)
       delr=rheadtail/dble(n_beads)
       if(chainaxis .eq. 'x') delr(1)=0.d0
       if(chainaxis .eq. 'y') delr(2)=0.d0
       if(chainaxis .eq. 'z') delr(3)=0.d0
       do i=2,n_beads
         joined_chain%coord(i,:)=joined_chain%coord(i,:)-dble(i)*delr
       end do

   end function join_head_to_tail


   function step_adder(unit_vec,seg_length,step_length,step_axis) result(stepped_vec)
     real*8 , dimension(3) :: unit_vec,stepped_vec
     real*8 :: seg_length,step_length
     integer :: step_dim
     character*1 :: step_axis
     real*8 , dimension(3) :: temp_vec
     real*8 , dimension(2) :: unit2d_vec

     if( step_axis .eq. 'x') then
       unit2d_vec(1)=unit_vec(2)/dsqrt(1.d0-unit_vec(1)**2)
       unit2d_vec(2)=unit_vec(3)/dsqrt(1.d0-unit_vec(1)**2)
!       print*,' unit 2d vec : ',unit2d_vec,dot_product(unit2d_vec,unit2d_vec)
       temp_vec(1)=step_length
       temp_vec(2)=unit2d_vec(1)*sqrt(seg_length**2-step_length**2)
       temp_vec(3)=unit2d_vec(2)*sqrt(seg_length**2-step_length**2)
     end if
     if( step_axis .eq. 'y') then
       unit2d_vec(1)=unit_vec(1)/dsqrt(1.d0-unit_vec(2)**2)
       unit2d_vec(2)=unit_vec(3)/dsqrt(1.d0-unit_vec(2)**2)
!       print*,' unit 2d vec : ',unit2d_vec,dot_product(unit2d_vec,unit2d_vec)
       temp_vec(2)=step_length
       temp_vec(1)=unit2d_vec(1)*sqrt(seg_length**2-step_length**2)
       temp_vec(3)=unit2d_vec(2)*sqrt(seg_length**2-step_length**2)
     end if
     if( step_axis .eq. 'z') then
       unit2d_vec(1)=unit_vec(1)/dsqrt(1.d0-unit_vec(3)**2)
       unit2d_vec(2)=unit_vec(2)/dsqrt(1.d0-unit_vec(3)**2)
!       print*,' unit 2d vec : ',unit2d_vec,dot_product(unit2d_vec,unit2d_vec)
       temp_vec(3)=step_length
       temp_vec(1)=unit2d_vec(1)*sqrt(seg_length**2-step_length**2)
       temp_vec(2)=unit2d_vec(2)*sqrt(seg_length**2-step_length**2)
     end if
       stepped_vec=temp_vec
     !  print*,'Step : ', temp_vec
     !  print*, ' Length of vec : ', sqrt(dot_product(temp_vec,temp_vec))


   end function step_adder

 function create_polymer_defects(chains,def_symbol,def_type,n_defects,def_treshold) result(iostat)
	    type(chain),dimension(:) :: chains
		character(len=2):: def_symbol
		integer :: def_type, n_defects
		integer :: del_atom,n_atoms,q
		integer :: iostat,i,j,cnt,del_chain,del_mer
		real :: rn,distance,def_cutoff
		real*8 , optional :: def_treshold
		logical :: near_flag=.false.
		  i=1
		  def_cutoff=2
		  if( present (def_treshold) ) def_cutoff=def_treshold
		  do while (i <= n_defects)
101		    call random_number(rn)
            ! pick the chain!
            del_chain=int(n_chains*rn)+1

            del_mer=int(chains(del_chain)%n_mers*rn)+1 ! picked monomer !
            call random_number(rn)
            del_atom=int(chains(del_chain)%mers(del_mer)%n_mer_atoms*rn)+1
			if( chains(del_chain)%mers(del_mer)%deleted(del_atom) .eq. 1) goto 101
			if( chains(del_chain)%mers(del_mer)%symbol(del_atom) .ne. def_symbol ) goto 101
			near_flag=.false.
            if( def_type .eq. 1) then  ! Creating mono vacancy defects
!		      do j=1, this%n_atoms
!			    if( j .ne. del_atom .and. this%deleted(j) .eq. 1) then
!			      distance=pbc_dist(this%coord(j,:),this%coord(del_atom,:))
!				  if (distance .lt. def_cutoff  ) near_flag=.true.
!				end if
 !             end do
!			  if ( .not. near_flag) then
			    i=i+1
				chains(del_chain)%mers(del_mer)%deleted(del_atom) = 1
!			  end if
            end if
			    cnt=0
			if ( def_type > 1) then
!                do q=1, this%n_atoms ! checking two monovacancy sites distance
!			        if( q .ne. del_atom .and. this%deleted(q) .eq. 1) then
!			            distance=pbc_dist(this%coord(q,:),this%coord(del_atom,:))
!				        if (distance .lt. 6*def_cutoff  ) near_flag=.true.
!                    end if
!                end do
!                if ( .not. near_flag) then
!                    this%deleted(del_atom)=1
!			        i=i+1
!			        do j=1, this%n_atoms
!				        if( this%symbol(j) .ne. def_symbol ) cycle
!			            distance=pbc_dist(this%coord(j,:),this%coord(del_atom,:))
!			            if( this%deleted(j) .eq. 0) then
!			                if( distance .lt. def_cutoff .and. cnt .lt. def_type-1) then
!			                    this%deleted(j)=1
!					            cnt=cnt+1
!			                end if
!			            end if
!			        end do
!			    end if
            end if
		  end do
		  write(*,'(A,x,i3,x,A,x,A)') 'Defects are created.' , n_defects, def_symbol, 'Atoms deleted'
		  iostat=0
	end function



   function omp_test() result (iostat)
use omp_lib
    implicit none

   integer :: iostat

  integer ( kind = 4 ) id
  real ( kind = 8 ) wtime
  iostat=0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HELLO_OPENMP'
  write ( *, '(a)' ) '  FORTRAN90/OpenMP version'

  wtime = omp_get_wtime ( )

  write ( *, '(a,i8)' ) &
    '  The number of processors available = ', omp_get_num_procs ( )
  write ( *, '(a,i8)' ) &
    '  The number of threads available    = ', omp_get_max_threads ( )
!
!  OUTSIDE THE PARALLEL REGION, have each thread say hello (there's only 1).
!
  id = omp_get_thread_num ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  OUTSIDE the parallel region.'
  write ( *, '(a)' ) ' '

  write ( *, '(a,i8,a,i8)' ) '  HELLO from process ', id

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Going INSIDE the parallel region:'
  write ( *, '(a)' ) ' '
!
!  INSIDE THE PARALLEL REGION, have each thread say hello.
!
!$omp parallel &
!$omp private ( id )
  id = omp_get_thread_num ( )

  write ( *, '(a,i8,a,i8)' ) '  HELLO from process ', id

!$omp end parallel
!
!  Finish up by measuring the elapsed time.
!
  wtime = omp_get_wtime ( ) - wtime

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Back OUTSIDE the parallel region.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Elapsed wall clock time = ', wtime
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HELLO_OPENMP'
  write ( *, '(a)' ) '  Normal end of execution.'
  end function omp_test

   end module polymer_module
