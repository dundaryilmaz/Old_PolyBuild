program polybuild
!----------------------------------------------------------------------------------
!
!   COMPOSITE POLYMER/OXIDE BUILDER
!
!       DUNDAR YILMAZ
!       
!    
!   
!
!------------------------------------------------------------------------------------
   use global_variables
   use polymer_module
   use polymer_type
   use mcmodule
   use ceramic_type
   use print_files
   use userinterface
   use omp_lib
!------------------------------
   implicit none
   integer :: i, j ,q,indx,surf_type,ntotal_chains
   real*8 :: total_penalty,save_penalty,penalty
   type(chain),dimension(:),allocatable :: save_chains
   real*8, dimension(3)  :: Va, Vb, Vc, Ua,Ub,Uc
   real*8  :: pos(3),dist,mindist
   real*8, dimension(:,:),allocatable :: sites
   logical :: success=.FALSE.
!====================================================
! Initialize the run get command line arguments
   integer :: n_arguments,ic,len, status
   character*80 :: cline
   integer :: iostat
   character*60 :: inputfile
   real*8 , dimension(3) :: box_o,box_l,x
     write(*,'(A)') 'Polymer Builder V 1.0'
     n_arguments=iargc()
	 BOX=20.0
     if( n_arguments < 1 ) then
       write(*,*) 'usage : ./compbuild inputscript'
       stop
     end if
     ic=1
          call get_command_argument(ic,cline,len,status)
    read(cline,*) inputfile
    iostat=inputscript_reader(inputfile)

end program polybuild
