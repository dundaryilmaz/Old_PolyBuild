module ceramic_type
    use checkbox
    use IFPORT
    implicit none
    type ceramic
        integer :: n_atoms
        integer :: block_id
        integer , allocatable,dimension(:) :: atype
        integer , allocatable,dimension(:) :: deleted
        integer , allocatable,dimension(:) :: surf_type
        integer , allocatable,dimension(:) :: functional_site
        character(len=3),allocatable , dimension (:) :: symbol
        integer, allocatable, dimension(:) :: indx,chain
        real*8 , allocatable , dimension(:,:) :: coord
        real*8 , allocatable ,dimension(:) ::  mass
        real*8 , dimension(3,3) :: h_matrix
    contains
        procedure, public :: print_block
        procedure , public :: create_surface_roughness
        procedure , public :: find_functional_sites
		procedure, public :: create_defects
    end type
    type(ceramic ) , allocatable , dimension(:) :: blocks
 !   type(ceramic ) :: base_ceramic

contains
    function read_unit_cell(filename) result(new_ceramic)
        type(ceramic) :: new_ceramic
        character(len=*) :: filename
        character(len=256) :: cwd,fullname
        integer :: n_atoms,i
        real*8 ,dimension(3):: vec
        real*8 :: dumy
        i=getcwd(cwd)
        fullname=cwd(:LNBLNK(cwd))//"/"//filename(:LNBLNK(filename))
        print*,fullname
        open(13,file=filename)
        read(13,*) n_atoms
        new_ceramic%n_atoms=n_atoms
        allocate(&
            &  new_ceramic%symbol(n_atoms) , &
            &  new_ceramic%indx(n_atoms)   , &
            &  new_ceramic%mass(n_atoms)   , &
            &  new_ceramic%coord(n_atoms,3), &
            &  new_ceramic%atype(n_atoms) ,&
            &  new_ceramic%chain(n_atoms) &
            )
        read(13,*) !scale
        read(13,*) new_ceramic%h_matrix(1,1:3)
        read(13,*) new_ceramic%h_matrix(2,1:3)
        read(13,*) new_ceramic%h_matrix(3,1:3)
        read(13,*) ! Direct
		!print*,new_ceramic%h_matrix(1,1)
		!print*,new_ceramic%h_matrix(2,2)
		!print*,new_ceramic%h_matrix(3,3)
        do i=1,n_atoms
            read(13,*) new_ceramic%atype(i),vec,&
                &      new_ceramic%chain(i) , new_ceramic%symbol(i),new_ceramic%mass(i)
            new_ceramic%coord(i,1) =vec(1)*new_ceramic%h_matrix(1,1)
            new_ceramic%coord(i,2) =vec(2)*new_ceramic%h_matrix(2,2)
            new_ceramic%coord(i,3) =vec(3)*new_ceramic%h_matrix(3,3)
        end do
		close(13)
    end function
!=========================================================
	function fiber_shell(center,radius_in,radius_out,base_ceramic) result( new_shell)
	    type(ceramic) :: new_shell,base_ceramic
		real*8 ,dimension(2)  :: center
		real*8 , dimension(2,2):: bases
		real*8 :: radius_out,dist,radius_in
		integer  :: nx,ny,nz
		logical :: count_flag= .true.
		integer :: cnt=0,ix,iy,iz,q,qb
		real*8 , dimension(3) :: origin,vec
		  nz=int(BOX(3)/base_ceramic%h_matrix(3,3))
          nx=int(radius_out/base_ceramic%h_matrix(1,1))
          ny=int(radius_out/base_ceramic%h_matrix(2,2))
		  bases(1,1:2)=(/0.d0,0.d0/)
		  bases(2,1:2)=(/0.5d0,0.5d0/)
10		  continue !print*,'fiber core creating ' 
		  do iz=1,nz
		    do ix=-nx,nx
		      do iy=-ny,ny
                 origin=(/center(1),center(2),0.d0/)+&
              &         (/               & 
              &         dble(ix)*base_ceramic%h_matrix(1,1),&
              &         dble(iy)*base_ceramic%h_matrix(2,2),&
              &         dble(iz-1)*base_ceramic%h_matrix(3,3)/)
               ! base1
			    do qb=1,2
				  origin(1:2)=origin(1:2)+(/bases(qb,1)*base_ceramic%h_matrix(1,1),&
					&   bases(qb,2)*base_ceramic%h_matrix(2,2)/)
					dist=sqrt(dot_product(origin(1:2)-center,origin(1:2)-center))
					if ( dist >  radius_out .or. dist < radius_in) cycle
					!write(*,'(A,4f8.2)') 'Base : ',origin(1:2),dist
					do q=1,base_ceramic%n_atoms
			           cnt=cnt+1
			           if( count_flag ) cycle
                       new_shell%symbol(cnt)=base_ceramic%symbol(q)
                       new_shell%atype(cnt) =base_ceramic%atype(q)
                       new_shell%mass(cnt) =base_ceramic%mass(q)
                       new_shell%coord(cnt,:)=base_ceramic%coord(q,:)+origin
                    end do
		         end do
               end do
			end do
	    end do
		if( .not. count_flag) then
		 write(*,'(A,x,f6.2,x,A,x,f6.2)') 'Fiber shell created with inner radius ',radius_in, &
&        		 'outer radius',radius_out
		 return
		end if
		!print*,' Counted'
	    write(*,'(A,x,i8)')'Number of atoms at the shell region:', cnt 
		new_shell%n_atoms=cnt
		allocate(new_shell%symbol(cnt) ,&
		     &   new_shell%atype(cnt), &
		     &   new_shell%mass(cnt), &
            &  new_shell%deleted(cnt), &
            &  new_shell%indx(cnt), &
            &  new_shell%functional_site(cnt), &
		     &   new_shell%coord(cnt,3))
			 cnt=0
			 count_flag=.false.
             goto 10
		 
    end function 
    !===================================================
	function fiber_core(center,radius,base_ceramic) result( new_core)
	    type(ceramic) :: new_core,base_ceramic
		real*8 ,dimension(2)  :: center
		real*8 , dimension(2,2):: bases
		real*8 :: radius,dist
		integer  :: nx,ny,nz
		logical :: count_flag= .true.
		integer :: cnt=0,ix,iy,iz,q,qb
		real*8 , dimension(3) :: origin,vec
		  nz=int(BOX(3)/base_ceramic%h_matrix(3,3))
          nx=int(radius/base_ceramic%h_matrix(1,1))
          ny=int(radius/base_ceramic%h_matrix(2,2))
		  bases(1,1:2)=(/0.d0,0.d0/)
		  bases(2,1:2)=(/0.5d0,0.5d0/)
10		  continue !print*,'fiber core creating ' 
		  do iz=1,nz
		    do ix=-nx,nx
		      do iy=-ny,ny
                 origin=(/center(1),center(2),0.d0/)+&
              &         (/               & 
              &         dble(ix)*base_ceramic%h_matrix(1,1),&
              &         dble(iy)*base_ceramic%h_matrix(2,2),&
              &         dble(iz-1)*base_ceramic%h_matrix(3,3)/)
               ! base1
			    do qb=1,2
				  origin(1:2)=origin(1:2)+(/bases(qb,1)*base_ceramic%h_matrix(1,1),&
					&   bases(qb,2)*base_ceramic%h_matrix(2,2)/)
					dist=sqrt(dot_product(origin(1:2)-center,origin(1:2)-center))
					if ( dist >  radius) cycle
					!write(*,'(A,4f8.2)') 'Base : ',origin(1:2),dist
					do q=1,base_ceramic%n_atoms
			           cnt=cnt+1
			           if( count_flag ) cycle
                       new_core%symbol(cnt)=base_ceramic%symbol(q)
                       new_core%atype(cnt) =base_ceramic%atype(q)
                       new_core%mass(cnt) =base_ceramic%mass(q)
                       new_core%coord(cnt,:)=base_ceramic%coord(q,:)+origin
                    end do
		         end do
               end do
			end do
	    end do
		if( .not. count_flag) then
		 write(*,'(A)') 'Fiber core created'
		 return
		end if
		!print*,' Counted'
		write(*,'(A,i8)')'Total number of atoms  : ', cnt 
		new_core%n_atoms=cnt
		allocate(new_core%symbol(cnt) ,&
		     &   new_core%atype(cnt), &
		     &   new_core%mass(cnt), &
            &  new_core%deleted(cnt), &
            &  new_core%indx(cnt), &
            &  new_core%functional_site(cnt), &
		     &   new_core%coord(cnt,3))
			 cnt=0
			 count_flag=.false.
             goto 10
		 
    end function 

    function fill_blocks(box_id,base_ceramic) result(new_block)
        type(ceramic) :: new_block,base_ceramic
        integer :: box_id
        integer :: nx,ny,nz
        integer :: i , j ,k ,q
        integer :: n_atoms,n_ucell,cnt
        real*8  :: r_nx, r_ny, r_nz
        real*8  , dimension(3) :: pos
        real*8  , dimension(3) :: vec
        real*8  , dimension(3,3) :: h_matrix
        real*8  :: surf_depth=2.d0
        real*8  :: z_cor, z_limit

        nx=nint(box_dimensions(box_id,4)/base_ceramic%h_matrix(1,1))
        ny=nint(box_dimensions(box_id,5)/base_ceramic%h_matrix(2,2))
        nz=nint(box_dimensions(box_id,6)/base_ceramic%h_matrix(3,3))
        r_nx=box_dimensions(box_id,4)/base_ceramic%h_matrix(1,1)
        r_ny=box_dimensions(box_id,5)/base_ceramic%h_matrix(2,2)
        r_nz=box_dimensions(box_id,6)/base_ceramic%h_matrix(3,3)
        h_matrix=base_ceramic%h_matrix
         print*, box_id
		 print*,'X: ', box_dimensions(box_id,4)/base_ceramic%h_matrix(1,1)
		 print*,'Y :', box_dimensions(box_id,5)/base_ceramic%h_matrix(2,2)
		 print*, 'Z' , box_dimensions(box_id,6)/base_ceramic%h_matrix(3,3) 
                 write(*,*) nx,ny,nz
                 write(*,*) nx*ny*nz*base_ceramic%n_atoms
          
        if ( abs(r_nx-nx) .gt. 0.01) then
            write(*,'(A,x,f8.4,x,i4)') 'Lattice mismatch in X direction : ', r_nx,nx
               h_matrix(1,1)=box_dimensions(box_id,4)/dble(nx)
            write(*,'(A,x,f8.4,x,A,x,f8.4,x,A,x,f5.2)') &
&           'Old lattice constant in x:',  base_ceramic%h_matrix(1,1),&
&           'New Lattice Constant :',h_matrix(1,1),&
&           'Strain applied:(%)',(base_ceramic%h_matrix(1,1)-h_matrix(1,1))/h_matrix(1,1)*100 
        end if
        if ( abs(r_ny-ny) .gt. 0.01) then
            write(*,'(A,x,f8.4,x,i4)') 'Lattice mismatch in Y direction : ', r_ny,ny
               h_matrix(2,2)=box_dimensions(box_id,5)/dble(ny)
            write(*,'(A,x,f8.4,x,A,x,f8.4,x,A,x,f5.2)') &
&           'Old lattice constant in y:',  base_ceramic%h_matrix(2,2),&
&           'New Lattice Constant :',h_matrix(2,2),&
&           'Strain applied:(%)',(base_ceramic%h_matrix(2,2)-h_matrix(2,2))/h_matrix(2,2)*100 
        end if
        if ( abs(r_nz-nz) .gt. 0.01) then
            write(*,'(A,x,f8.4,x,i4)') 'Lattice mismatch in Z direction : ', r_nz,nz
               h_matrix(3,3)=box_dimensions(box_id,6)/dble(nz)
            write(*,'(A,x,f8.4,x,A,x,f8.4,x,A,x,f5.2)') &
&           'Old lattice constant in z:',  base_ceramic%h_matrix(3,3),&
&           'New Lattice Constant :',h_matrix(3,3),&
&           'Strain applied:(%)',(base_ceramic%h_matrix(3,3)-h_matrix(3,3))/h_matrix(3,3)*100 
        end if
!        nx=int(box_lengths(1)/base_ceramic%h_matrix(1,1))
!        ny=int(box_lengths(2)/base_ceramic%h_matrix(2,2))
!        nz=int(box_lengths(3)/base_ceramic%h_matrix(3,3))
        n_atoms=nx*ny*nz*base_ceramic%n_atoms
        n_ucell=base_ceramic%n_atoms
        new_block%n_atoms=n_atoms
        new_block%block_id=box_id
        write(*,'(A,x,i2)')'Filling Block : ', box_id
        !print*,' nx : ',nx,' ny : ',ny,' nz : ',nz
        !print*,(box_lengths(1)/base_ceramic%h_matrix(1,1))
        !print*,' Number of atoms in block : ',n_atoms
        !print*,'Size-x : ',dble(nx)*base_ceramic%h_matrix(1,1)
        !print*,'Size-y : ',dble(ny)*base_ceramic%h_matrix(2,2)
        !print*,'Size-z : ',dble(nz)*base_ceramic%h_matrix(3,3)
        allocate(&
            &  new_block%symbol(n_atoms) , &
            &  new_block%indx(n_atoms)   , &
            &  new_block%mass(n_atoms)   , &
            &  new_block%coord(n_atoms,3), &
            &  new_block%atype(n_atoms), &
            &  new_block%deleted(n_atoms), &
            &  new_block%surf_type(n_atoms), &
            &  new_block%functional_site(n_atoms) &
            )
        new_block%surf_type(:)=0
        new_block%functional_site(:)=0

        new_block%deleted(:)=0
        cnt=0
        do i=1,nx
            do j=1,ny
                do k=1,nz
                    !     print*,i,j,k
                    pos=box_dimensions(box_id,:)+&
                        (/&
                        &         dble(i-1)*h_matrix(1,1),&
                        &         dble(j-1)*h_matrix(2,2),&
                        &         dble(k-1)*h_matrix(3,3)/)
                    !          print*,i , j ,k,box_id
                    !          print*,box_dimensions(box_id,:)
                    !          print*,pos
                    do q=1,base_ceramic%n_atoms
             !           if( (i .eq. 1) .and. (base_ceramic%chain(q) .eq. 1) ) cycle
						cnt=cnt+1

                        new_block%symbol(cnt)=base_ceramic%symbol(q)
                        new_block%atype(cnt) =base_ceramic%atype(q)
                        new_block%mass(cnt) =base_ceramic%mass(q)
                        new_block%coord(cnt,:)=base_ceramic%coord(q,:)+pos

                    !        new_block%coord((cnt-1)*n_ucell+1:cnt*n_ucell,1)=base_ceramic%coord(:,1)+pos(1)
                    !        new_block%coord((cnt-1)*n_ucell+1:cnt*n_ucell,2)=base_ceramic%coord(:,2)+pos(2)
                    !        new_block%coord((cnt-1)*n_ucell+1:cnt*n_ucell,3)=base_ceramic%coord(:,3)+pos(3)
                    end do
                end do
            end do
        end do
		n_atoms=cnt
		!print*, 'N atoms : ',cnt
        !     pause
        !detect surface atoms
        !-----------------------------------------------
        ! surf_type = 1 upper xy plane
        ! surf_type = 2 lower xy plane
        ! surf_type = 0 bulk
        !-----------------------------------------------
        do i=1,n_atoms
            vec= new_block%coord(i,:)
            z_cor=new_block%coord(i,3)
            z_limit=box_dimensions(new_block%block_id,3)+box_lengths(3)-surf_depth
            if(z_cor > z_limit) new_block%surf_type(i)=1
            z_limit=box_dimensions(new_block%block_id,3)+surf_depth
            if(z_cor < z_limit) new_block%surf_type(i)=2
        end do
    end function
    function print_block(this,filename,indx) result(iostat)
        class(ceramic) :: this
        character(len=*) :: filename
        integer :: indx,cnt,i
        integer :: iostat
        if( indx == 0) then
            !---- print header ---- '
            open(134,file=filename)
            write(134,1001) 'CRYST1',box,90.0,90.0,90.0
        end if
        cnt=indx
		!print*,' Index : ', indx,this%n_atoms
        !  pause
        do i=1,this%n_atoms
            if( this%deleted(i) == 0 ) then
                cnt=cnt+1
                !print*,i,cnt

                write(134,1002) 'ATOM  ', cnt,this%symbol(i),this%surf_type(i),&
                    &                 this%coord(i,:),this%symbol(i)
            end if
        end do
        iostat=cnt


1001    format(A6,3(f9.3),x,3(f7.2,x))
1002    format(A6,i5,A4,3x,i2,10x,3(f8.3),23x,A2)
    end function print_block
    function find_functional_sites(this,surf_type,site_pos) result(iostat)
        class(ceramic) :: this
        integer :: surf_type
        real*8  :: surf_depth=2.d0
        real*8  :: rn
        real*8 , dimension(3) ,intent(out):: site_pos
        integer :: iostat
        integer :: pick_atom
        iostat=0
10      call random_number(rn)
        pick_atom=int(this%n_atoms*rn)+1

        if(this%surf_type(pick_atom) .ne. surf_type ) goto 10
        if(this%functional_site(pick_atom) .ne. 0 ) goto 10
        if(this%deleted(pick_atom) .ne. 0 ) goto 10
        this%functional_site(pick_atom) = 1
        !    print*,'Found a functional site : ',this%symbol(pick_atom),' ',&
        !&    	this%coord(pick_atom,:)
        if( surf_type .eq. 1 ) then ! upper xy plane
            site_pos=this%coord(pick_atom,:)+(/0.d0,0.d0,1.5d0/)
        end if
        if( surf_type .eq. 2 ) then ! lower xy plane
            site_pos=this%coord(pick_atom,:)+(/0.d0,0.d0,-1.5d0/)
        end if



    end function

    function create_defects(this,def_symbol,def_type,n_defects,def_treshold) result(iostat)
	    class(ceramic) :: this
		character(len=2):: def_symbol
		integer :: def_type, n_defects
		integer :: del_atom,n_atoms,q
		integer :: iostat,i,j,cnt
		real :: rn,distance,def_cutoff
		real*8 , optional :: def_treshold
		logical :: near_flag=.false.
          n_atoms=this%n_atoms
		  i=1
		  def_cutoff=2
		  if( present (def_treshold) ) def_cutoff=def_treshold
		  do while (i <= n_defects)
101		    call random_number(rn)
            del_atom=int(n_atoms*rn)+1 ! picked random atom
			if( this%deleted(del_atom) .eq. 1) goto 101
			if( this%symbol(del_atom) .ne. def_symbol ) goto 101
			near_flag=.false.
            if( def_type .eq. 1) then  ! Creating mono vacancy defects
		      do j=1, this%n_atoms
			    if( j .ne. del_atom .and. this%deleted(j) .eq. 1) then
			      distance=pbc_dist(this%coord(j,:),this%coord(del_atom,:))
				  if (distance .lt. def_cutoff  ) near_flag=.true.
				end if
              end do
			  if ( .not. near_flag) then
			    i=i+1
				this%deleted(del_atom) = 1
			  end if
            end if
			    cnt=0
			if ( def_type > 1) then
                do q=1, this%n_atoms ! checking two monovacancy sites distance
			        if( q .ne. del_atom .and. this%deleted(q) .eq. 1) then
			            distance=pbc_dist(this%coord(q,:),this%coord(del_atom,:))
				        if (distance .lt. 6*def_cutoff  ) near_flag=.true.
                    end if
                end do
                if ( .not. near_flag) then
                    this%deleted(del_atom)=1
			        i=i+1
			        do j=1, this%n_atoms
				        if( this%symbol(j) .ne. def_symbol ) cycle
			            distance=pbc_dist(this%coord(j,:),this%coord(del_atom,:))
			            if( this%deleted(j) .eq. 0) then
			                if( distance .lt. def_cutoff .and. cnt .lt. def_type-1) then
			                    this%deleted(j)=1
					            cnt=cnt+1
			                end if
			            end if
			        end do
			    end if
            end if
		  end do
		  write(*,'(A,x,i3,x,A,x,A)') 'Defects are created.' , n_defects, def_symbol, 'Atoms deleted'
		  iostat=0
	end function
    function create_surface_roughness(this, surf_type,rough_type ) result(iostat)
        class(ceramic) :: this
        integer :: surf_type
        integer :: rough_type
        integer :: iostat,cnt_al,cnt_o
        integer :: i , j , k,cnt
        integer , allocatable,dimension(:) :: surf_list
        real*8 :: surf_depth=2.d0
        real*8 , dimension(3) :: vec
        real*8 :: z_limit,z_cor,rn
        if( surf_type == 1 ) then ! upper surface
            z_limit=box_dimensions(this%block_id,3)+box_lengths(3)-surf_depth
        end if
        if( surf_type == 2 ) then ! lower surface
            z_limit=box_dimensions(this%block_id,3)+surf_depth
        end if
        cnt=0
        cnt_o=0
        cnt_al=0
        do i=1,this%n_atoms
            z_cor=this%coord(i,3)
            if( surf_type == 1 ) then
                if ( z_cor > z_limit ) then
                    call random_number(rn)
                    if( rn > 0.5d0 ) then
                        this%deleted(i)=1
                        cnt=cnt+1
                        if( this%symbol(i) .eq. 'Al') cnt_al=cnt_al+1
                        if( this%symbol(i) .eq. 'O ') cnt_o=cnt_o+1
                      ! print*,'Deleted ',cnt,this%symbol(i)
                    end if
                end if
            end if
            if( surf_type == 2 ) then
                if ( z_cor < z_limit ) then
                    call random_number(rn)
                    if( rn > 0.5d0 ) then
                        this%deleted(i)=1
                        cnt=cnt+1
                        if( this%symbol(i) .eq. 'Al') cnt_al=cnt_al+1
                        if( this%symbol(i) .eq. 'O ') cnt_o=cnt_o+1
                      !print*,'Deleted ',cnt,this%symbol(i)
                    end if
                end if
            end if
        end do
        print*,cnt_al, ' Al atoms ',cnt_o , 'Oxygen Atoms '
        iostat=0
    end function  create_surface_roughness

end module ceramic_type
