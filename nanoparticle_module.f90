module nanoparticle_module
	use global_variables
	type nanoparticle
	  integer :: n_atoms
	  integer, allocatable , dimension(:) :: atype
	  real*8 , allocatable, dimension(:,:) :: coord
	  real*8 , allocatable,dimension(:) :: mass
	  character(len=2) , allocatable,  dimension(:) :: symbol
	  integer :: n_vertices
	  integer :: n_surface
	  integer :: n_facets
	  integer,allocatable ,dimension (:) :: vertices, surface_atoms,nearest_facet,nearest_vertex
	  integer,allocatable,dimension(:,:) :: facets
	  real*8 ,allocatable ,dimension(:,:) :: normalvectors
	  real*8 , allocatable, dimension(:,:) :: functional_sites
	  integer ,allocatable,dimension(:):: surf_facet_list
	 end type nanoparticle
	contains
	!-----------------------
	function read_nanoparticle(filenanoparticle,surfacedata) result( particle)
      type(nanoparticle) :: particle
	  character(len=*) :: filenanoparticle,surfacedata
	  integer :: i, j , k,dd,ii,jj,kk
	  real*8 :: cm(3),shift(3)
	  open(11,file=filenanoparticle,status='old')
	  open(12,file=surfacedata,status='old')
      
	  Print*,'Box : ', BOX
	  read(11,*) dd
	  particle%n_atoms=dd
	  allocate(particle%coord(dd,3),particle%atype(dd),&
&   	       particle%symbol(dd),particle%mass(dd))
	  read(11,*) 
	  do i=1,dd
        read(11,*) particle%symbol(i),particle%coord(i,:)
	  end do
	  close(11)
	  cm(1)=sum(particle%coord(1,:))/dd
	  cm(2)=sum(particle%coord(2,:))/dd
	  cm(3)=sum(particle%coord(3,:))/dd
	  shift(1)=BOX(1)/2.d0-cm(1)
	  shift(2)=BOX(2)/2.d0-cm(2)
	  shift(3)=BOX(3)/2.d0-cm(3)
	  do i=1,dd
         particle%coord(i,:)=particle%coord(i,:)+shift
	  end do
	  read(12,*) dd
	  particle%n_vertices=dd
	  allocate(particle%vertices(dd))
	  read(12,*)
	  !This data come from Python, indexes start from 0
	  do i=1,dd
	    read(12,*) ii
		particle%vertices(i)=ii+1
	  end do
	  read(12,*) dd
	  particle%n_facets=dd
	  allocate(particle%facets(dd,3),particle%normalvectors(dd,3))
	  read(12,*)
	  do i=1,dd
    	 read(12,*) ii,jj,kk,particle%normalvectors(i,:)
	     particle%facets(i,1)=ii+1
	     particle%facets(i,2)=jj+1
	     particle%facets(i,3)=kk+1
	  end do

      read(12,*)
      read(12,*)
      read(12,*)
      read(12,*)dd
      read(12,*)
	  particle%n_surface=dd
	  allocate(particle%surface_atoms(dd),particle%nearest_facet(dd),particle%nearest_vertex(dd))
	  !allocate(particle%surf_facet_list(dd))
	  !surf_facet_list is the list of facets of surface atoms
	  ! particle%surf_facet_list(5)=7 means 5 th surface atom belongs to 7th facet!
	  do i=1,dd
        read(12,*) ii , jj,kk
	    particle%surface_atoms(i)=ii+1
	    particle%nearest_facet(i)=jj+1
	    particle%nearest_vertex(i)=kk+1
	  end do
    end function read_nanoparticle 
	  


end module nanoparticle_module
