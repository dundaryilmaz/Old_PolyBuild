module checkbox
use global_variables
implicit none
    integer :: n_boxes
	real*8 , dimension(:,:) , allocatable :: box_dimensions
    real*8 , dimension(3) :: box_lengths
	real*8 :: box_buff

contains
    function random_pos_wall(walltype) result(vector)
      character*2 :: walltype
      real*8 , dimension(3) :: vector
      real*8 :: tmp_vec
      call random_number(vector)
        vector(1)=BOX(1)*vector(1)
        vector(2)=BOX(2)*vector(2)
        vector(3)=BOX(3)*vector(3)
        if( walltype .eq.'-x' ) vector(1)=0.d0
        if( walltype .eq.'-y' ) vector(2)=0.d0
        if( walltype .eq.'-z' ) vector(3)=0.d0
        if( walltype .eq.'+x' ) vector(1)=BOX(1)
        if( walltype .eq.'+y' ) vector(2)=BOX(2)
        if( walltype .eq.'+z' ) vector(3)=BOX(3)

    end function random_pos_wall
    function vecmag(vecA) result(magA)
    real*8 , dimension(:) :: vecA
    real*8 :: magA
      maga=sqrt(dot_product(vecA,vecA))
    end function vecmag
    function random_pos_box() result (vector)
    real*8 ,dimension (3) :: vector
    call random_number(vector)

        vector(1)=BOX(1)*vector(1)
        vector(2)=BOX(2)*vector(2)
        vector(3)=BOX(3)*vector(3)
    end function random_pos_box
    !=========================================================
    function set_boxes() result(iostat)
        !	character(len=*) :: boxfile
        integer:: i,iostat
        !      open(101,file=boxfile)
        !         read(101,*) n_boxes,box_lengths,box_buff
        ! no scaling on this version
        !          box_lengths=box_lengths/SCAL
        !          box_buff=box_buff/SCAL
        !           allocate(box_dimensions(n_boxes,3))
        do i=1,n_boxes
            read(101,*) box_dimensions(i,1:6)
        !         box_dimensions(i,:)=box_dimensions(i,:)-BOXI/2
        end do
        !          box_dimensions=box_dimensions/SCAL
        !          print*,box_lengths
        !          do i=1,n_boxes
        !           write(*,'(3f16.4)')box_dimensions(i,:)
        !          end do
        !  pause
        !          close(101)
        iostat=0
    end function set_boxes
    !=========================================================
    function pbc_x(X) result(Xp)
        ! returns periodic boundary condition
        ! applied Xp
        real*8 , dimension (3) :: X,Xp
        integer :: m
        Xp=X
11        do m=1,3
            if( X(M) .lt. -BOX(M)/2) Xp(M)=X(M)+BOX(M)
            if( X(M) .gt.  BOX(M)/2) Xp(M)=X(M)-BOX(M)
        end do
		 if( (Xp(1) .lt. -BOX(1)/2) .or. &
&            (Xp(2) .lt. -BOX(2)/2) .or. &
&            (Xp(3) .lt. -BOX(3)/2) .or. &
&            (Xp(1) .gt.  BOX(1)/2) .or. &
&            (Xp(2) .gt.  BOX(2)/2) .or. &
&            (Xp(3) .gt.  BOX(3)/2) ) then
!               write(*,'(A,3f8.2)')'Before: ',X
!               write(*,'(A,3f8.2)')'After: ',Xp
!			   pause
			   X=Xp
			   goto 11
		 end if


    end function pbc_x
    !========================================================
    function return_box(X) result( Xr)
        real*8 , dimension(3) :: X,Xr
        integer :: m
        Xr=X
12        do m=1,3
            if( X(M) .lt.  0.d0 ) Xr(M)=X(M)+BOX(M)
            if( X(M) .gt.  BOX(M)) Xr(M)=X(M)-BOX(M)
        end do
		   if( (Xr(1) .lt. 0) .or. &
&              (Xr(2) .lt. 0) .or. &		   
&              (Xr(3) .lt. 0) .or. &		   
&              (Xr(1) .gt. BOX(1)) .or. &		   
&              (Xr(2) .gt. BOX(2)) .or. &		   
&              (Xr(3) .gt. BOX(3)) ) then
!               write(*,'(A,3f8.2)')'Before: ',X
!               write(*,'(A,3f8.2)')'After: ',Xr
			   X=Xr
			   goto 12
		   end if

    end function return_box
    !=========================================================
    function pbc_dist(A,B) result(distance)
        !returns distance between A B in a periodic boundary box
        real*8 , dimension(3) :: A, B,AB,ABp
        real*8  :: distance

        AB=A-B
        ABp=pbc_x(AB)
        distance=sqrt(dot_product(ABp,ABp))

    end function pbc_dist
    !=========================================================
    function insidevoids_buff(X) result(iostat)
        !----------------------------------------------------------
        ! returns 0 if point X is inside one of the boxes defined
        ! here
        !----------------------------------------------------------
        integer:: iostat,i,ii
        real*8 , dimension (3) :: X,Xp,Xk
        real*8 , dimension (3) :: box1_c1,box1_c2,box1_c3,box1_c4
        real*8 , dimension (3) :: box1_c5,box1_c6,box1_c7,box1_c8
        real*8 , dimension (8,8,3) :: boxes
        real*8 , dimension (3) :: box_origin,high_corner,shift_vec,temp_x,l_boxl
		real*8 , dimension(3) :: box_dimen
        real*8  :: xlat,ylat,zlat,temp,dist
        logical :: box_contains_point_nd,inside
        !xlat=31.0 ; ylat=16.1 ; zlat=5.6;
        !xlat=31.0 ; ylat=16.5 ; zlat=6.2;
        iostat=0
		temp_x=X
		temp_x=return_box(temp_x)
        !temp_x=(X+BOX/2)
        !l_boxl=box_lengths*SCAL
        !do i=1,3
        !    if( temp_x(i) .lt. 0 ) temp_x(i)=temp_x(i)+box(i)
        !    if( temp_x(i) .gt. box(i)) temp_x(i)=temp_x(i)-box(i)
        !end do
            ! temp_x=temp_x*SCAL
           !  print*,'X :',temp_x*SCAL
		!-------------------------------
		!  fiber mode
		!-------------------------------
		if ( fiber_mode) then
		  dist=sqrt(dot_product(X(1:2)-fiber_center,X(1:2)-fiber_center))
		  if ( dist < fiber_radius + 1.5d0 ) then
		    iostat= 0
		    else
		    iostat=1
		  end if
		  return
		end if

		if ( fiber_shell_mode) then
		  dist=sqrt(dot_product(X(1:2)-fiber_center,X(1:2)-fiber_center))
		  if ( dist > fiber_radius_in - 1.5d0 ) then
		    iostat= 0
		    else
		    iostat=1
		  end if
		  return
		end if
        do i=1,n_boxes

            !print*,' box : ',i
            box_origin=box_dimensions(i,1:3)
            !inside=new_inside_voids_buff(temp_x,box_origin,box_lengths,box_buff)
            inside=new_inside_voids_buff(temp_x,box_origin,box_dimensions(i,4:6),box_buff)
			if ( inside ) return
	    end do
        iostat=1

           ! box_origin=box_origin+BOX/2
            !print*,'bo ',scal*box_origin
           ! shift_vec=BOX/2-(box_lengths/2+box_origin)
           ! box_origin=box_origin+shift_vec-box_buff
           ! l_boxl=box_lengths+2*box_buff
            !pause
            !print*,'sbo :',scal*(box_origin)
           ! high_corner=box_origin+l_boxl
           ! temp_x=X+BOX/2+shift_vec
           ! do ii=1,3
           !     if( temp_x(ii) .lt. 0 ) temp_x(ii)=temp_x(ii)+box(ii)
           !     if( temp_x(ii) .gt. box(ii)) temp_x(ii)=temp_x(ii)-box(ii)
           ! end do
           ! if( box_contains_point_nd(3,box_origin,high_corner,temp_x) )then
           !     iostat=0
                !               print*,'point is inside one of boxes',i
                !               pause
           !     return
           ! end if
        !end do
        !print*,X
                !print*,Xp
                !print*,'----'
      !         pause

    end function insidevoids_buff
   function new_inside_voids_buff(X,box_origin,box_length,buff) result(inside)
     real*8 ,dimension(3) :: X,box_origin,box_length
	 real*8 :: buff
	 integer :: i , j , k
	 logical :: inside
       i=1
       do while (i <= 3) 
	    inside=line_contains_point(X(i),box_origin(i),box_origin(i)+box_length(i),BOX(i),buff)
		if( .not. inside) return
		i=i+1
		
	   end do
	    
		 !if( box_origin(i)-buff > 0 ) th

  end function
   function line_contains_point( X,Xlo,Xhi,BOXL,buff) result (inside)
     real*8 :: X,Xlo,Xhi,BOXL,buff
     logical :: inside
	 real*8 :: Xlotemp,Xhitemp
	   inside=.FALSE.
    !    print*,'X = ', X
	!	print*,'Xlo= ', Xlo, ' Xhi = ',Xhi
	!	print*,'BOXL= ',BOXL, ' Buff = ',buff
	    if( (Xlo-buff) .ge. 0.d0 .and. (Xhi+buff) .le. BOXL) then
          if( (X .le. (Xlo-buff) ) .or. (X .ge. (Xhi+buff))) return		  
          elseif( ((Xlo-buff) .lt. 0 ).and. (Xhi+buff .le. BOXL)) then
		    if( (X .lt. (Xlo-buff+BOXL)) .and. (X .gt. (Xhi+buff))) return
          elseif( ((Xlo-buff) .gt. 0) .and.  (Xhi+buff) .gt. BOXL ) then
		    if( (X .lt. (Xlo-buff)) .and. (X .gt. (Xhi+buff-Boxl))) return
		end if
		inside=.TRUE.
  end function



end module checkbox

