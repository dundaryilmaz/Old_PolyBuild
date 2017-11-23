module userinterface
use global_variables
use polymer_type
use polymer_module
use ceramic_type
use print_files
use mcmodule
use nanoparticle_module
implicit none
contains

function inputscript_reader(filename) result(iostat)
character(len=*) :: filename
character(len=120) ::cline
character(len=40) :: word,outputfilename
character(len=40) :: fiber_ucell
character(len=40) :: monomer_file
character(len=80),dimension(:),allocatable :: block_ucells
integer :: iostat,ilo,iho,indx,n_words,ihi,i
integer :: n_mer_atoms,j,n_mcsteps,q,qaccp,p,qq,k,nfunctional
integer :: n_def_groups
integer,allocatable ::n_defects(:),def_type(:),def_block(:)
character(len=2) ,allocatable :: def_symbol(:)
logical :: done
real*8 :: coord(3),save_penalty,total_penalty,inter_penalty,self_penalty,vec(3)
real*8 :: mctemp,rn
real*8, allocatable:: def_treshold(:)
character(len=2) :: symbol
character(len=40) ::npfilename, npsurfdata
logical :: nanoparticle_mode=.false.
type(chain) , dimension(:), allocatable :: save_chains
type(nanoparticle) :: particle
type(ceramic),allocatable :: base_unitcells(:)
open(11,file=filename,status='old')
open(12,file='log.txt', status='unknown')
nfunctional=0
100  read(11,'(A80)',end=101) cline
  if( cline(1:1) .eq. '#') goto 100
   n_words=0
   ilo=0
   ihi=0
   done=.true.
   call s_word_count(cline,n_words)
   call word_next_read(cline,word,done)
   if( word .eq. 'monomer') then
     call word_next_read(cline,word,done)
     read(word,*) monomer_file
     iostat=base_mer%create(monomer_file)
     write(*,'(A,x,A)') 'Monomer name :' , base_mer%monomer_name
     write(12,'(A,x,A)') 'Monomer name :' , base_mer%monomer_name
     write(*,'(A,x,i2)') 'Number of atoms in monomer : ',n_mer_atoms
     write(12,'(A,x,i2)') 'Number of atoms in monomer : ',n_mer_atoms
     write(12,'(A,x,3f8.4)')' Base Vector - 1 :',base_mer%base_vecs(1,:)
     write(12,'(A,x,3f8.4)')' Base Vector - 2 :',base_mer%base_vecs(2,:)
     write(12,'(A,3f8.4)')'Head ', base_mer%headh
     write(12,'(A,3f8.4)')'Tail ', base_mer%tailh
     write(*,'(A,f8.4)') 'Mer Width is :',mer_width
   end if
   if( word .eq. 'simulation_box') then
     call word_next_read(cline,word,done)
     read(word,*) BOX(1)
     call word_next_read(cline,word,done)
     read(word,*) BOX(2)
     call word_next_read(cline,word,done)
     read(word,*) BOX(3)
     write(*,'(A,x,3f8.2)') 'Simulation Box : ',BOX(1:3)
     write(12,'(A,x,3f8.2)') 'Simulation Box : ',BOX(1:3)
   end if
   !--------------------------------------------------------------
   !  Create Blocks which will be filled with crystal structures |
   !--------------------------------------------------------------
   
   if ( word .eq. 'blocks') then
     call word_next_read(cline,word,done)
     read(word,*) n_boxes
     allocate(base_unitcells(n_boxes),block_ucells(n_boxes))
     allocate(box_dimensions(n_boxes,6))
     if( n_boxes > 0 ) then
       do i=1,n_boxes
         read(11,'(A)') cline
         done=.true.
         call word_next_read(cline,word,done)
         read(word,*) box_dimensions(i,1)       
         call word_next_read(cline,word,done)
         read(word,*) box_dimensions(i,2)
         call word_next_read(cline,word,done)
         read(word,*) box_dimensions(i,3)
	 call word_next_read(cline,word,done)
         read(word,*) box_dimensions(i,4)
         call word_next_read(cline,word,done)
         read(word,*) box_dimensions(i,5)
         call word_next_read(cline,word,done)
         read(word,*) box_dimensions(i,6)
         call word_next_read(cline,word,done)             
         read(word,*) box_buff
         call word_next_read(cline,word,done)
         read(word,'(A)') block_ucells(i)
         write(*,'(A,x,3f8.2,x,A,x,3f8.2,x,A,x,A)')&
&              'box origins : ', box_dimensions(i,1:3),&
&              'box dimensions : ',box_dimensions(i,4:6),&
&              'Unit Cell File: ',trim(block_ucells(i))
       end do
     end if
   end if
   if( word .eq. 'defects') then
     call word_next_read(cline,word,done)
     read(word,*) n_def_groups
     allocate(n_defects(n_def_groups), def_symbol(n_def_groups),&
&             def_block(n_def_groups),  def_type(n_def_groups), &
&             def_treshold(n_def_groups))
     do i=1,n_def_groups
        read(11,*) def_block(i),def_symbol(i),def_type(i),n_defects(i),&
&                  def_treshold(i)
     end do
   end if
   if ( word .eq. 'amorph_chain') then
     call word_next_read(cline,word,done)
     read(word,*) n_amorph_chains
     call word_next_read(cline,word,done)
     read(word,*)mers_per_chain
   end if
   if( word .eq. 'fiber' ) then
     fiber_mode=.true.
     call word_next_read(cline,word,done)
     read(word,*) fiber_center(1)
     call word_next_read(cline,word,done)
     read(word,*) fiber_center(2)
     call word_next_read(cline,word,done)
     read(word,*) fiber_radius
     call word_next_read(cline,word,done)
     read(word,*) fiber_ucell
   end if
   if(word .eq. 'fiber_shell') then
     fiber_shell_mode=.true.
     call word_next_read(cline,word,done)
     read(word,*) fiber_center(1)
     call word_next_read(cline,word,done)
     read(word,*) fiber_center(2)
     call word_next_read(cline,word,done)
     read(word,*) fiber_radius_in
     call word_next_read(cline,word,done)
     read(word,*) fiber_radius_out
     call word_next_read(cline,word,done)
     read(word,*) fiber_ucell
   end if
   if ( word .eq. 'cont_chain') then
     call word_next_read(cline,word,done)
     read(word,*) n_cont_chains
     call word_next_read(cline,word,done)
     read(word,*)mers_per_chain_cont
     call word_next_read(cline,word,done)
     read(word,*)radius_cont
   end if
   if ( word .eq. 'mcsteps') then
     call word_next_read(cline,word,done)
     read(word,*) n_mc_steps
   end if
   if( word .eq. 'filename') then
     call word_next_read(cline,word,done)
     read(word,*) outputfilename
   end if
   if( word.eq. 'nanoparticle') then
     call word_next_read(cline,word,done)
     read(word,*) npfilename
     call word_next_read(cline,word,done)
     read(word,*) npsurfdata
     nanoparticle_mode=.TRUE.
   end if
   if(word .eq. 'functional') then
     call word_next_read(cline,word,done)
     read(word,*) nfunctional
  end if
  goto 100
101 continue
!------------------------------------------|
!           driver part                    |
!------------------------------------------|
!           create amorphous chains        |
!------------------------------------------|
iostat=init_random_seed()
if(nanoparticle_mode )&
&   call nanoparticle_driver(npfilename,npsurfdata,nfunctional,particle)
n_chains=n_cont_chains+n_amorph_chains+nfunctional
write(*,'(A,A)') '# of Amorphous chains |# of Continuous Chains |# of Functional Chains | Total # of Chains '
write(*,'(8x,i4,8x,8x,i4,8x,10x,i4,8x,10x,i4,8x)')n_amorph_chains,n_cont_chains,nfunctional,n_chains
allocate(chains(n_chains))
if( n_chains > 0 )iostat=pre_analize()
!------------------------------------------|
! Create Continious Chains                 |
!------------------------------------------|
if(n_cont_chains .ne. 0 ) then
  do i=1, n_cont_chains
    chains(i)=build_continuous_chain(mers_per_chain_cont,i,'z')
  end do
  write(*,'(A)') 'Continuous chains created '
end if
!------------------------------------------
! Create Amorphous Chains
!------------------------------------------
if( n_amorph_chains .gt. 0  ) then
  do i=n_cont_chains+1,n_cont_chains+n_amorph_chains
    chains(i)=create_chain(mers_per_chain,i)
    total_penalty=self_chain_interaction(chains(i))
  end do
  write(*,"(i3,x,'Amorphous Chains Created')") n_amorph_chains  
end if	
!------------------------------------------
! Create Nano-Particle
!------------------------------------------
if( nanoparticle_mode) then
  q=n_cont_chains+n_amorph_chains
  do i=1,nfunctional
    chains(q+i)=create_singlemerchain(particle%functional_sites(i,1:3),&
&               particle%functional_sites(i,4:6)	,q)
  end do
  write(*,'(A)') 'Functional chains created'
end if
do i=1,n_chains
  do j=1,chains(i)%n_beads
    vec=chains(i)%coord(j,:)
    p=insidevoids_buff(vec)
  end do
end do

allocate(save_chains(n_chains))
!n_mc_steps=20000
save_penalty=100
if( n_mc_steps .gt. 0) then
  qaccp=0
  mctemp=1000
  do qq=1,10
    mctemp=mctemp*0.5
    do q=1,n_mc_steps
      save_chains(1:n_amorph_chains) = chains(n_cont_chains+1:n_chains)
      do i=n_cont_chains+1,n_chains
        iostat=chain_move(chains(i))
      end do
      total_penalty=0
      self_penalty=0
      inter_penalty=0
      do i=n_cont_chains+1,n_chains-1
        do j=i+1,n_chains
          inter_penalty=inter_penalty+&
&                    two_chain_interaction(chains(i),chains(j))
        end do
      end do
      do i=n_cont_chains+1,n_chains
	self_penalty=self_penalty+self_chain_interaction(chains(i))
      end do
      total_penalty=self_penalty+inter_penalty
      if( q*qq .eq. 1) then
        if( total_penalty > save_penalty ) save_penalty = 2*total_penalty
          write(*,1002) total_penalty,save_penalty,q,n_mc_steps
	end if
	call random_number(rn)
        if( exp((-total_penalty+save_penalty)/mcTemp) >  rn ) then
          qaccp=qaccp+1
	  if( mod(qaccp,100)== 0) &
&           write(*,1002) total_penalty,save_penalty,qaccp,q,n_mc_steps
          save_penalty =  total_penalty
          else
          chains(n_cont_chains+1:n_chains)=save_chains
        end if
        if( mod(q,1000) == 0 )  print*,'---',q,n_mc_steps
      end do
    end do
  end if
1002 format('Current Penalty : ',e12.5,3x,'Previous Penalty',e12.5,3x,&
&                     'Accepted Step: ',i6,x,'Step :',i6,x,'of',x,i6,' mc steps')
if( n_chains .gt. 0) then
  iostat=print_all_beads(chains,'mbeads.pdb')
end if
box_buff=box_buff-mer_width/2
do i=1,n_chains
  if( chains(i)%functional) then
    iostat=chains(i)%put_monomer(1)
  else
    iostat=chains(i)%put_all_monomers()
  end if
  iostat=chains(i)%map_atypes()
end do
do i=1,n_chains
  do j=1,chains(i)%n_mers
    do q=1,n_mer_atoms
      vec=chains(i)%mers(j)%coord(q,:)
      p=insidevoids_buff(vec)
      if( p .eq. 0 ) &
&	write(*,'(A,i3,x,A,i3,A,x,3f8.2,x,A,x,3f8.2,x,A,i4)') &
&	'Chain: ',i,' Bead :',j,' Bead Pos : ',&
&       chains(i)%coord(j,:),' Atom Pos :', vec,' Inside : ',p
    end do
  end do
end do
indx=0
do i=1,n_chains
  indx=chains(i)%print_polymer_chain('polymer.pdb',indx)
end do
if( n_boxes == 0 ) &
&    iostat=print_chains(chains,'test.lmps')
iostat=1
do i=1,n_chains
  if(chains(i)%chaintype .eq. 'amorph') then
    vec=chains(i)%tail
    vec=return_box(vec)
    chains(i)%tail=vec
    vec=chains(i)%head
    vec=return_box(vec)
    chains(i)%head=vec
  end if
  do q=1,chains(i)%n_mers
    do k=1,chains(i)%mers(q)%n_mer_atoms
      vec=chains(i)%mers(q)%coord(k,:)
      vec=return_box(vec)
      chains(i)%mers(q)%coord(k,:)=vec
    end do
  end do
end do
if( allocated(n_defects)) then
  do i=1,n_def_groups
    if(def_block(i) .eq. 0) then
      iostat=create_polymer_defects(chains,def_symbol(i),def_type(i),n_defects(i),def_treshold(i))
    end if
  end do 
end if
!-----------------------------------------------------------
! ceramic building stage
!-----------------------------------------------------------
if( n_boxes > 0  ) then
  allocate(blocks(n_boxes))
  do i=1,n_boxes
    base_unitcells(i)=read_unit_cell(block_ucells(i))
    blocks(i)=fill_blocks(i,base_unitcells(i))
  end do
  if( allocated(n_defects)) then
    do i=1, n_def_groups
      iostat=blocks(def_block(i))%create_defects(def_symbol(i),def_type(i),n_defects(i),def_treshold(i))
    end do
  end if
  if( rough_surface) then
    write(*,*) 'Implementing Surface Roughness'
    iostat=blocks(1)%create_surface_roughness(1,1)
    iostat=blocks(2)%create_surface_roughness(2,1)
    print*,'Bottom Block: ',blocks(1)%n_atoms-sum(blocks(1)%deleted(:))
    print*,'Top    Block: ',blocks(2)%n_atoms-sum(blocks(2)%deleted(:))
  end if
  indx=0
  do i=1,n_boxes
    indx=blocks(i)%print_block('block.pdb',indx)
  end do
end if
if ( fiber_mode) then
  n_boxes=1
  allocate(base_unitcells(1))
  base_unitcells(1)=read_unit_cell(fiber_ucell)
  allocate(blocks(1))
  blocks(1)=fiber_core(fiber_center,fiber_radius,base_unitcells(1))
  if( allocated(n_defects)) &
&      iostat=blocks(def_block(i))%create_defects(def_symbol(i),&
&      def_type(i),n_defects(i),def_treshold(i))
end if
if ( fiber_shell_mode) then
  n_boxes=1
  allocate(base_unitcells(1))
  base_unitcells(1)=read_unit_cell(fiber_ucell)
  allocate(blocks(1))
  blocks(1)=fiber_shell(fiber_center,fiber_radius_in,fiber_radius_out,base_unitcells(1))
  if( allocated(n_defects)) then
    do i=1,n_def_groups
      if( def_block(i) .eq. 1) &
&       iostat=blocks(def_block(i))%create_defects(def_symbol(i),&
&              def_type(i),n_defects(i),def_treshold(i))
    end do
  end if
end if
if(allocated(blocks))then
  if(nanoparticle_mode ) then
    iostat=print_all(chains,blocks,outputfilename,0,particle)
  else
    iostat=print_all(chains,blocks,outputfilename)
  end if
  else
  if(nanoparticle_mode ) then
    iostat=print_all(chains,blocks,outputfilename,0,particle)
  else
    iostat=print_all(chains,blocks,outputfilename,0)
  end if
end if
end function inputscript_reader
subroutine nanoparticle_driver(npfilename,npsurfdata,nfunc,particle)
   character(len=*) :: npfilename,npsurfdata
   type(nanoparticle) :: particle
   integer  :: nfunc
   integer :: i, j , k,q,surindex,cnt,try
   real*8 :: rn,vec1(3),vec2(3),dist,normal_vec(3),mindistance
   integer , allocatable,dimension(:) :: temp_list
   real*8 , allocatable,dimension(:,:) :: vec_list
   real*8 ,allocatable, dimension(:,:) :: functional_sites
   logical :: closeto=.false.
    
     particle=read_nanoparticle(npfilename,npsurfdata)
	 allocate(particle%functional_sites(nfunc,6))
     allocate(temp_list(particle%n_surface),vec_list(particle%n_surface,6))	 
	 q=0
	 cnt=0
	   try=0
	   mindistance=4.0
	 do while(cnt .le. nfunc)
	   call random_number(rn)
	   !-----------------------
	   ! how to select functional sites
	   !1- get a coplanar atom
	   k=int(dble(particle%n_surface)*rn+1)
	   surindex=particle%surface_atoms(k)
	   q=particle%nearest_facet(k)
	   !print*,'Selected surface atom : ', k,' th surface atom'
	   !print*,'Index of selected atom : ' , surindex
	   !print*,'Nearest facet of ', k, 'th surface atom is ' ,q
	   normal_vec=particle%normalvectors(q,:)
	  
	   !print*,'Normal vector : ',vec1
       vec1=particle%coord(surindex,:)+normal_vec*1.5d0
	   closeto=.FALSE.
	   do j=1,cnt
		 dist=sqrt(dot_product(vec1-vec_list(j,:),(vec1-vec_list(j,:)))) 
	     ! print*,'j ' , j,'cnt',cnt,dist
		 if(dist .lt. mindistance) then

		! print*,'Close : ',closeto
		 closeto=.TRUE.
		 end if
		 try=try+1
		 if( try .eq. 100) then
		  mindistance=mindistance-0.1d0
		  !print*,try, cnt, mindistance
		  try=0
		  !pause
		  end if
	     !vec1=particle
       end do
	   if ( .not. closeto) then
	     if( mindistance .lt. 2.0 ) then
		  write(*,'(A)') 'WARNING: TOO CLOSE CHAINS'
		 end if
	     mindistance=3.5
	     vec_list(cnt,1:3)=vec1
		 vec_list(cnt,4:6)=normal_vec
		 cnt=cnt+1
	   !write(*,'(A,3f8.4)') 'Site vector: ',vec1
	   end if
	   !pause
	 end do
	   !allocate(functional_sites(cnt,6))
	   particle%functional_sites=vec_list(1:cnt,:)
	  ! print*,functional_sites
	  ! do i=1,nfunc

	
end subroutine



end module userinterface
