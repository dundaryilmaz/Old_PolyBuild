module print_files

    use polymer_type
    use ceramic_type
    use global_variables
    use polymer_module
	use nanoparticle_module
    implicit none
contains
    function print_all( chains, blocks , filename,blocknumber,particle) result (iostat)
        type(chain)   ,dimension(:) :: chains
        type(ceramic) ,dimension(:) :: blocks
		type(nanoparticle) , optional :: particle
		integer , optional :: blocknumber
		integer :: nblocks
        character(len=*) :: filename
		character(len=50) ::filexyz, filelammps,filecfg,filebgf
        integer :: iostat
        integer :: n_total_atoms,n_polymer_atoms,n_ceramic_atoms,n_nanopartatoms
        integer :: i , j , k,cnt,temp,maxtype,maxtypeo,q
		real*8 ,dimension(3) :: rvec
		nblocks=size(blocks)
		if( present(blocknumber) ) nblocks=blocknumber
!        print*,'-',nblocks
        filelammps=trim(filename)//'.lmp'
		filecfg=trim(filename)//'.cfg'
		filexyz=trim(filename)//'.xyz'
		filebgf=trim(filename)//'.bgf'
 !       print*,'--'
        temp=0
        n_ceramic_atoms=0
		if( nblocks .gt.0 ) then
        n_ceramic_atoms=sum(blocks(:)%n_atoms)
        do i=1,nblocks
            temp=temp+sum(blocks(i)%deleted(:))
			if( sum(blocks(i)%deleted(:)) > 0 ) &
            write(*,'(A,i3)')'Deleted from block ',i , sum(blocks(i)%deleted(:))
            write(*,'(A,x,i2,x,A,x,i8)')&
&			'# of atoms in Block',i,':',blocks(i)%n_atoms-sum(blocks(i)%deleted(:))

        maxtype=maxval(blocks(i)%atype(:))
        maxtypeo=max(maxtype,maxtypeo)
        print*,maxtype,maxtypeo
        end do
		end if
		n_nanopartatoms=0
		if( present(particle)) n_nanopartatoms=particle%n_atoms
        n_ceramic_atoms=n_ceramic_atoms-temp
        n_polymer_atoms=sum(chains(:)%n_mers)
        temp=0
        do i=1,n_chains
            do j=1,chains(i)%n_mers
                do q=1,chains(i)%mers(j)%n_mer_atoms
                    temp=temp+chains(i)%mers(j)%deleted(q)
           
                end do
            end do
             if ( chains(i)%chaintype .eq. 'amorph') temp=temp-2 ! add tail/head H atoms
        end do
        n_polymer_atoms=base_mer%n_mer_atoms*n_polymer_atoms-temp ! subtract deleted!
        n_total_atoms=n_polymer_atoms+n_ceramic_atoms+n_nanopartatoms
        maxtype=maxval(base_mer%atype(:))
        maxtypeo=max(maxtype,maxtypeo)
        !print*, maxtypeo
		if( n_polymer_atoms  > 0 ) &
		write(*,'(A,x,i6)')'Polymer Atoms:',n_polymer_atoms
		if( present(particle))&
&		write(*,'(A,x,i8)')'Nanoparticle Atoms :',n_nanopartatoms
        write(*,'(A,x,i8)')'Total Number of atoms : ',n_total_atoms
        open(144,file=filelammps)
        open(145,file=filecfg)
		open(146,file=filexyz)
		open(147,file=filebgf)
        write(144,1001)n_polymer_atoms,n_ceramic_atoms
        write(144,*)
        write(144,1002) n_total_atoms
        write(144,1003)
        write(144,1004) maxtypeo
        write(144,*)
        write(144,1005) 0.d0 , BOX(1)
        write(144,1006) 0.d0 , BOX(2)
        write(144,1007) 0.d0 , BOX(3)
        write(144,*)
        write(144,1008)
        write(144,*)
        write(145,'(A,x,i6)') 'Number of particles = ',n_total_atoms
		write(145,'(A)') 'A =1.0'
		write(145,'(A,x,f13.6)')'H0(1,1)=', BOX(1)
		write(145,'(A,x,f13.6)')'H0(1,2)=', 0.0 
		write(145,'(A,x,f13.6)')'H0(1,3)=', 0.0 
		write(145,'(A,x,f13.6)')'H0(2,1)=', 0.0 
		write(145,'(A,x,f13.6)')'H0(2,2)=', BOX(2)
		write(145,'(A,x,f13.6)')'H0(2,3)=', 0.0 
		write(145,'(A,x,f13.6)')'H0(3,1)=', 0.0 
		write(145,'(A,x,f13.6)')'H0(3,2)=', 0.0 
		write(145,'(A,x,f13.6)')'H0(3,3)=', BOX(3)
		write(145,'(A)') '.NO_VELOCITY.'
		write(145,'(A)')'entry_count = 4 '
		write(145,'(A)')'auxiliary[0] = q'
		write(146,*) n_total_atoms
		write(146,*)
		write(147,'(A)')'XTLGRF 200'
		write(147,'(A,x,A)')'DESCRP',filename
		write(147,'(A)')'REMARK BGF FILE CREATED WITH POLYBUILD'
		write(147,'(A,2x,6(f10.5,x))') 'CRYSTX',BOX(1:3),90.0,90.0,90.0
		write(147,'(A)')'FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)'
205  format(a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)		
        cnt=0
        do i=1,n_boxes
            do j=1,blocks(i)%n_atoms
                if( blocks(i)%deleted(j) .eq. 0) then
                    cnt=cnt+1
!                    write(144,1009) cnt,symbol_to_hell(blocks(i)%symbol(j)),0.d0, blocks(i)%coord(j,:)
                    write(144,1009) cnt,blocks(i)%atype(j),0.d0, blocks(i)%coord(j,:)
					write(145,'(f12.4)') blocks(i)%mass(j)
					write(145,'(A3)') blocks(i)%symbol(j)
					rvec(1)=blocks(i)%coord(j,1)/BOX(1)
					rvec(2)=blocks(i)%coord(j,2)/BOX(2)
					rvec(3)=blocks(i)%coord(j,3)/BOX(3)
					write(145,'(4(f13.6,x))') rvec,1.d0
					write(146,*) blocks(i)%symbol(j),blocks(i)%coord(j,:)
					write(147,205)'HETATM',cnt,blocks(i)%symbol(j),'','','',blocks(i)%coord(j,:),blocks(i)%symbol(j),0,0,0.0
                end if
            end do
        end do
        !write(144,*)'---------'
        do i=1,n_chains+nchains_functional
            do j=1,chains(i)%n_mers
                do k=1,chains(i)%mers(j)%n_mer_atoms
                    if(chains(i)%mers(j)%deleted(k) .eq. 0) then
                    cnt=cnt+1
                        write(144,1009) cnt,chains(i)%mers(j)%atype(k),&
                            &    0.d0, chains(i)%mers(j)%coord(k,:)
						    write(145,'(f12.4)') chains(i)%mers(j)%mass(k)
						    write(145,'(A3)') chains(i)%mers(j)%symbol(k)
					    rvec(1)=chains(i)%mers(j)%coord(k,1)/BOX(1)
					    rvec(2)=chains(i)%mers(j)%coord(k,2)/BOX(2)
					    rvec(3)=chains(i)%mers(j)%coord(k,3)/BOX(3)
					    write(145,'(4(f13.6,x))') rvec,1.d0
					    write(146,*) chains(i)%mers(j)%symbol(k),chains(i)%mers(j)%coord(k,:)
					    write(147,205)'HETATM',cnt,chains(i)%mers(j)%symbol(k),'','','',chains(i)%mers(j)%coord(k,:),&
&            					chains(i)%mers(j)%symbol(k),0,0,0.0
                    end if
                end do
            end do
            if ( chains(i)%chaintype .eq. 'amorph') then
                cnt=cnt+1
                write(144,1009) cnt,base_mer%htype, 0.d0, chains(i)%head(:)
                rvec=chains(i)%head/box
                write(145,'(f12.4)') 1.0
                write(145,'(A3)') 'H '
                write(145,'(4(f13.6,x))') rvec,1.d0
                write(146,*) 'H  ',chains(i)%head(:)
                write(147,205)'HETATM',cnt,'H  ','','','',chains(i)%head(:),'H  ',0,0,0.0

                cnt=cnt+1
                write(144,1009) cnt,base_mer%htype, 0.d0,chains(i)%tail(:)
                rvec=chains(i)%tail/box
                write(145,'(f12.4)') 1.0
                write(145,'(A3)') 'H '
                write(145,'(4(f13.6,x))') rvec,1.d0
                write(146,*) 'H  ',chains(i)%tail(:)
                write(147,205)'HETATM',cnt,'H  ','','','',chains(i)%tail(:),'H  ',0,0,0.0
             end if
                
        end do
		if( present(particle))then
		  cnt=cnt+1
		  do i=1,n_nanopartatoms
		  write(144,1009) cnt,particle%atype(i),0.d0 , particle%coord(i,:)
		  write(145,'(f12.4)') particle%mass(i)
		  write(145,'(A3)') particle%symbol(i)
		  rvec(1)=particle%coord(i,1)/BOX(1)
		  rvec(2)=particle%coord(i,2)/BOX(2)
		  rvec(3)=particle%coord(i,3)/BOX(3)
		  write(145,'(4(f13.6,x))') rvec,1.d0
		  write(146,*) particle%symbol(i), particle%coord(i,:)
		  write(147,205)'HETATM',cnt,particle%symbol(i),'','','',particle%coord(i,:),particle%symbol(i),0,0,0.0
		  end do
	    end if	  


        iostat=0
        !-----------------------------------------------------------
1001    format('Generated with NewCompositeBuilder ',i10,' Polymer Atoms ',i10,&
            &           ' Ceramic Atoms ')
1002    format(i10,' atoms ')
1003    format(5x,'0 impropers')
1004    format(i3,5x,' atom types')
1005    format(2(f16.4,2x) , 'xlo xhi')
1006    format(2(f16.4,2x) , 'ylo yhi')
1007    format(2(f16.4,2x) , 'zlo zhi')
1008    format('Atoms')
1009    format(i10,2x,i2,3x,f8.2,3f16.6)
    end function print_all
    !----------------------------
    function print_chains(chains,filename) result(iostat)
        type(chain) :: chains(:)
        character(len=*):: filename
        integer :: iostat
        integer :: i, j,k,cnt,n_polymer_atoms
        real*8 , dimension(3) :: vec
        n_polymer_atoms=sum(chains(:)%n_mers)
        n_polymer_atoms=base_mer%n_mer_atoms*n_polymer_atoms+n_amorph_chains*2
        cnt=0
        open(144,file=filename)

        write(144,1001)n_polymer_atoms
        write(144,*)
        write(144,1002) n_polymer_atoms
        write(144,*)
        write(144,1004)
        write(144,*)
        write(144,1005) 0.d0 , BOX(1)
        write(144,1006) 0.d0 , BOX(2)
        write(144,1007) 0.d0 , BOX(3)
        write(144,*)
        write(144,1008)
        write(144,*)
        do i=1,size(chains(:))
            !print*,' chain no : ', i
            if( chains(i)%chaintype .eq. 'amorph') then
                vec=chains(i)%coord(1,:)-chains(i)%coord(2,:)
                vec=vec/sqrt(dot_product(vec,vec))
                vec=1.15d0*vec
                vec=chains(i)%coord(1,:)+vec
                vec=chains(i)%head(:)
                 ! vec is the position of chain head cap hydrogen!
                cnt=cnt+1
                ! Hydrogen type for reax is 2 !
                write(144,1009)cnt,2 , 0.d0   ,vec
            end if

            do j=1,chains(i)%n_mers

                do k=1,chains(i)%mers(j)%n_mer_atoms
                    cnt=cnt+1
                    write(144,1009) cnt,chains(i)%mers(j)%atype(k),&
                        &    0.d0, chains(i)%mers(j)%coord(k,:)
                end do
            end do

            if( chains(i)%chaintype .eq. 'amorph') then
                vec=chains(i)%coord(chains(i)%n_mers,:)-chains(i)%coord(chains(i)%n_mers-1,:)
                vec=vec/sqrt(dot_product(vec,vec))
                vec=1.15d0*vec
                vec=chains(i)%coord(chains(i)%n_mers,:)+vec
                 ! vec is the position of chain head cap hydrogen!
                 vec=chains(i)%tail(:)
                cnt=cnt+1
                ! Hydrogen type for reax is 2 !
                write(144,1009)cnt,2 , 0.d0   ,chains(i)%tail(:)
            end if




        end do
        iostat=0
1001    format('Generated with NewCompositeBuilder ',i10,' Polymer Atoms ',i10,&
            &           ' Ceramic Atoms ')
1002    format(i10,' atoms ')
1004    format(5x,'4 atom types')
1005    format(2(f16.4,2x) , 'xlo xhi')
1006    format(2(f16.4,2x) , 'ylo yhi')
1007    format(2(f16.4,2x) , 'zlo zhi')
1008    format('Atoms')
1009    format(i10,2x,i2,3x,f8.2,3f16.6)



    end function print_chains


    function symbol_to_hell(symbol) result(hellid)
        character(len=*) :: symbol
        integer :: hellid

        hellid=0
        if( symbol .eq. 'O ') hellid = 1
        if( symbol .eq. 'C ') hellid = 2
        if( symbol .eq. 'H ') hellid = 3
        if( symbol .eq. 'Al') hellid = 4
    end function symbol_to_hell



end module print_files
