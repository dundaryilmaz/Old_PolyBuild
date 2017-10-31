 module mcmodule
   !------------------------------------------------
   !monte carlo module for the new composite builder
   !first coded nov 30
   ! on the way from konya to istanbul
   !-----------------------------------------------
   use global_variables
   use checkbox
   use polymer_type
   implicit none
   contains
   function self_chain_interaction(chainA) result(penalty)
     real*8 :: penalty
	 type(chain) :: chainA
	 integer :: i , j , k
	 real*8 :: distance,mer_length
	 penalty =0.d0
	 mer_length=base_mer%mer_length
	   if(chainA%n_beads < 3 ) return
	   do i=3,chainA%n_beads
	     distance=pbc_dist(chainA%coord(i-2,:),chainA%coord(i,:))
!		 write(*,'(A,x,i4,x,A,x,i4,x,i4,7f8.2)')'Chain: ',chainA%chain_index,&
!& 		  'Distance beads : ',i-2,i,distance,chainA%coord(i-2,:),chainA%coord(i,:)
		 if( distance < mer_length*1.7d0 ) then
		   penalty=penalty+exp(-1.5*(distance-mer_length*1.7)/(mer_length*1.7))
		 end if
	  end do
!	  pause
   end function
	     
   function two_chain_interaction(chainA,chainB) result(penalty)
     real*8 :: penalty
     type(chain) :: chainA,chainB
     integer :: i , j , k
     real*8 :: distance
     penalty=0.d0
!	 print*,mer_width,'---'
     do i=1,chainA%n_beads
       do j=1,chainB%n_beads
         distance=pbc_dist(chainA%coord(i,:),chainB%coord(j,:))
         if( distance > mer_width ) then !avg_bead_distance) then
           penalty=penalty+exp(-1.5*(distance/mer_width))
!          write(*,'(A,x,i4,x,A,x,i4,x,f8.2,x,f8.2,x,f8.2)')&
!&   		  'i : ',i,' j  ',j,distance,exp(-2.0*(distance/mer_width)),penalty
          !     chainA%close_flag(i)=1
          !     chainB%close_flag(j)=1
		  else
!          write(*,'(A,x,i4,x,A,x,i4,x,f8.2,x,f16.2,x,f16.2)')&
!&   		  'i : ',i,' j  ',j,distance,exp(-2.0*(distance/mer_width)),penalty
		   penalty=penalty+exp(1.5*(mer_width/distance))
         end if
       end do
     end do
!  print*,penalty
!	 pause
   end function two_chain_interaction

   function chain_move(mc_chain) result(iostat)
     type(chain) :: mc_chain
     integer :: i , j , k
     integer :: iostat,cnt
     real*8 , dimension(3) :: vec,vecA,vecB,vecC,Ua,Ub,Uc
     real*8 :: scale_cf,length,comp1,comp2
     do i=1,mc_chain%n_beads
	   cnt=0
       if( .not. mc_chain%functional ) then
100		 continue
         vec=0.01d0*random_unit_vec()
         cnt=cnt+1
		 if( cnt == 10 ) goto 101
		 if( insidevoids_buff(mc_chain%coord(i,:)+vec) .eq. 0 ) goto 100
         mc_chain%coord(i,:)=mc_chain%coord(i,:)+vec
       end if
101	   continue	   
     end do
         ! rescaling
     do i=2,mc_chain%n_beads
       vecA=mc_chain%coord(i,:)-mc_chain%coord(i-1,:)
       vecB=pbc_x(vecA)
!	   write(*,'(2(i3,x,A,3f8.2),x,A,3f8.2)')&
!&	   i-1,'.th bead vec = > ',mc_chain%coord(i-1,:),&
!&      i,' th bead vec => ', mc_chain%coord(i,:),&
!&      'Distance vec : ', vecB
	   vec=vecB/sqrt(dot_product(vecB,vecB))
       length=pbc_dist(mc_chain%coord(i,:),mc_chain%coord(i-1,:)) !sqrt(dot_product(vec,vec))
!	   write(*,'(A,x,2f8.2)') 'Distance', length,sqrt(dot_product(vec,vec))
       scale_cf=base_mer%mer_length/length
	   !write(*,'(i4,3f8.3)'),i,length,base_mer%mer_length,scale_cf
       vecA=vec*base_mer%mer_length
	   vecB=mc_chain%coord(i-1,:)+vecA
	   vec=return_box(vecB)
       mc_chain%coord(i,:)=vec
!	   pause

     end do
     iostat=0
	 !pause
 end function chain_move
!  function chain_move(mc_chain) result(iostat)
!    type(chain) :: mc_chain
!    integer :: i , j , k
!    integer :: iostat,cnt
!    real*8 , dimension(3) :: vec,vecA,vecB,vecC,Ua,Ub,Uc
!    real*8 :: scale_cf,length,comp1,comp2
!    do i=1,mc_chain%n_beads
!      if( .not. mc_chain%functional .or. i > 2 ) then
!        vec=0.1d0*random_unit_vec()
!        if( i ==1 ) then
!          vecA=mc_chain%coord(2,:)-mc_chain%coord(1,:)
!          else
!          vecA=mc_chain%coord(i,:)-mc_chain%coord(i-1,:)
!        end if
		!cnt=0
!100		!continue
!        cnt=cnt+1
!        vecB=random_unit_vec()
!        vecC=random_unit_vec()
!        call gsch(VecA,VecB,VecC,Ua,Ub,Uc)
!        call random_number(comp1)
!        call random_number(comp2)
!        vec=comp1*Ub+comp2*Uc
!        vec=vec/sqrt(dot_product(vec,vec))
!        vec=vec*0.01d0
		!if( cnt == 100 ) cycle
		!if( insidevoids_buff(mc_chain%coord(i,:)+vec) .eq. 0 ) goto 100
!        mc_chain%coord(i,:)=mc_chain%coord(i,:)+vec
!      end if
!    end do
!        ! rescaling
!    do i=3,mc_chain%n_beads
!      vec=mc_chain%coord(i,:)-mc_chain%coord(i-1,:)
!      length=sqrt(dot_product(vec,vec))
!      scale_cf=base_mer%mer_length/length
!      vec=vec*scale_cf
!      mc_chain%coord(i,:)=mc_chain%coord(i-1,:)+vec
!    end do
!    iostat=0
!end function chain_move

end module mcmodule
