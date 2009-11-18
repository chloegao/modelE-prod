      module ent_make_struct

      use ent_types
      use ent_const
      use cohorts
      use patches
      use entcells

      public ent_readcsv

      contains

!************************************************************************
      subroutine skipstar (iu_entstruct)
      !* Skip a comment line in an Ent structure csv file.
      integer :: iu_entstruct
      !---
      character :: check

      read(iu_entstruct,*) check
      write(*,*) check
      do while (check.eq.'*') 
        read(iu_entstruct,*) check
      end do
      backspace(iu_entstruct)
      end subroutine skipstar

!************************************************************************


      subroutine read_entcell_struct ( ecp, iu_entstruct )
      type(entcelltype) :: ecp
      integer :: iu_entstruct
      !---
      !real*8 :: stext1,stext2,stext3,stext4,stext5 !soil_texture
      
      !read(iu_entstruct,*) stext1,stext2,stext3,stext4,stext5 !soil_texture
      !ecp%soil_texture(1) = stext1
      !ecp%soil_texture(2) = stext2
      !ecp%soil_texture(3) = stext3
      !ecp%soil_texture(4) = stext4
      !ecp%soil_texture(5) = stext5
      call skipstar(iu_entstruct)
      read(iu_entstruct,*) ecp%soil_texture 
      end subroutine read_entcell_struct

!************************************************************************
      

      subroutine read_patch_struct ( pp,iu_entstruct )
      type(patch),pointer :: pp
      integer :: iu_entstruct
      !---
      integer :: layer
!      real*8 :: age, area
!      real*8 :: tc1,tc2,tc3,tc4,tc5,tc6,tc7,tc8,tc9,tc10,tc11,tc12 !TpoolC
!      real*8 :: tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10,tn11,tn12 !TpoolN

      call skipstar(iu_entstruct)
      read(iu_entstruct,*) pp%age, pp%area, pp%soil_type
      do layer = 1,N_CASA_LAYERS
        call skipstar(iu_entstruct)
        read(iu_entstruct,*) pp%Tpool(CARBON,:,layer)
        read(iu_entstruct,*) pp%Tpool(NITROGEN,:,layer)
      end do
      end subroutine read_patch_struct

!************************************************************************

      subroutine read_cohort_struct ( cop, iu_entstruct )
      type(cohort) :: cop
      integer :: iu_entstruct
      !---
!      real*8 :: nm,Ntot,LAI,LMA,h,crown_dx,crown_dy,dbh,root_d,clump
!      real*8 :: fracroot !skip
!      real*8 :: C_fol,N_fol,C_sw,N_sw,C_hw,N_hw,C_lab,N_lab,
!     &     C_froot,N_froot,C_croot,N_croot

      call skipstar(iu_entstruct)
      read(iu_entstruct,*) cop%pft
      call skipstar(iu_entstruct)
      read(iu_entstruct,*) cop%h,cop%crown_dx,cop%crown_dy,
     &     cop%dbh,cop%root_d,cop%clump
      write(*,*) cop%h,cop%crown_dx,cop%crown_dy,
     &     cop%dbh,cop%root_d,cop%clump
      call skipstar(iu_entstruct)
      read(iu_entstruct,*) cop%C_fol,cop%N_fol,cop%C_sw,cop%N_sw,
     &     cop%C_hw,cop%N_hw,cop%C_lab,cop%N_lab,cop%
     &     C_froot,cop%N_froot,cop%C_croot,cop%N_croot

      end subroutine read_cohort_struct
!************************************************************************
      subroutine check_ij_bounds(i,j,I0,I1,J0,J1)
      use filemanager
      implicit none
      integer :: i,j,I0,I1,J0,J1

      if ((i.lt.I0).or.(i.gt.I1)) then
         call stop_model("ent_make_struct i out of bounds",255)
      elseif ((j.lt.J0).or.(j.gt.J1)) then
         call stop_model("ent_make_struct j out of bounds",255)
      endif

      end subroutine check_ij_bounds
!************************************************************************

      subroutine ent_struct_readcsv (IM,JM,I0,I1,J0,J1)
!@sum Read ent vegetation structure from ASCII CSV file
!@sum Order must be:
!@sum First line:  N_CASA_LAYERS
!@sum e, entcells in any order spanning dimension IM,JM
!@sum p, patches in order from oldest to youngest
!@sum c, cohorts in order from tallest to shortest
!@sum
!@sum This ordering is not critical, because the insert routines will
!@sum sort the patches and cohorts, given age and height.

      use FILEMANAGER
      implicit none
      integer :: IM,JM,I0,I1,J0,J1
      type(entcelltype) :: cells(I0:I1,J0:J1)
      type(entcelltype) :: cellsbuf(I0:I1,J0:J1)
      !----
      integer :: iu_entstruct
      character*80, parameter :: ENTSTRUCT="ENTSTRUCT"
      integer :: N_CASAlayers    !N_CASA_layers
      integer :: i,j
      type(entcelltype),pointer :: ecp
      type(patch),pointer :: pp
      type(cohort),pointer :: cop
      real*8 :: fracroot(N_DEPTH)
      logical :: MORE
      character :: next

      fracroot(1:N_DEPTH) = 0.d0
      MORE = .true.  !initialize

!      open(iu_entstruct,file="ent_struct_init.csv",status='unknown')
      call openunit(trim(ENTSTRUCT),iu_entstruct,.false.,.true.)
      read(iu_entstruct,*) N_CASAlayers
      write(*,*) 'N_CASAlayers',N_CASAlayers
      do while(MORE)
        read(iu_entstruct,*) next
        if (next.eq.'$') then !End of data
          MORE = .false.
        else if (next.eq.'*') then !skip comment
        else if (next.eq.'e') then !new entcell
          write(*,*) 'e'
          read(iu_entstruct,*) i,j
          write(*,*) i,j
          call check_ij_bounds(i,j,I0,I1,J0,J1)
          call entcell_construct(ecp)
          call zero_entcell( ecp )
          call read_entcell_struct( ecp,iu_entstruct )
          cells(i,j) = ecp
        else if (next.eq.'p') then !new patch
          write(*,*) 'p'
          call insert_patch(ecp,0.d0,0)
          pp => ecp%youngest
          call read_patch_struct(pp,iu_entstruct)
          write(*,*) 'patch area age soil',pp%area,pp%age,pp%soil_type
        else if (next.eq.'c') then !new cohort
          write(*,*) 'c'
          call insert_cohort(pp,0,
     &         0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
     &         fracroot,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
     &         0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
     &         0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)
          cop => pp%shortest
          call read_cohort_struct(cop,iu_entstruct)
          write(*,*) 'cohort pft height',cop%pft,cop%h
        end if
        cells(i,j) = ecp
      end do

!      write(*,*) 'Got here 0.'
!      do i = I0,I1
!        do j = J0,J1
!          call entcell_print(6,cells(i,j))
!        end do
!      end do

      write(*,*) 'Summarizing entcells...'
      do i = I0,I1
        do j = J0,J1
          call summarize_entcell(cells(i,j))
          write(*,*) '...after summarize_entcell'
          call entcell_print(6,cells(i,j))
          cellsbuf(i,j) = cells(i,j)
          write(*,*) '...after assign cellsbuf(i,j)'
        end do
      end do

      write(*,*) 'Printing cellsbuf...'
      do i = I0,I1
         do j = J0,J1
            call entcell_print(6, cellsbuf(i,j))
         end do
      end do 
      
      call closeunit(iu_entstruct)
      write(*,*) 'Done with creating Ent vegetation structure.'
      end subroutine ent_struct_readcsv


      subroutine ent_readcsv (ecp, iu_entstruct)
!@sum Read ent vegetation structure from ASCII CSV file
!@sum Order must be:
!@sum First line:  N_CASA_LAYERS
!@sum e, entcells in any order spanning dimension IM,JM
!@sum p, patches in order from oldest to youngest
!@sum c, cohorts in order from tallest to shortest
!@sum
!@sum This ordering is not critical, because the insert routines will
!@sum sort the patches and cohorts, given age and height.

      implicit none
      type(entcelltype),pointer :: ecp
      integer :: iu_entstruct      
      !----
      type(patch),pointer :: pp
      type(cohort),pointer :: cop
      real*8 :: fracroot(N_DEPTH)
      character :: next
      integer :: counter

      fracroot(1:N_DEPTH) = 0.d0

      call entcell_construct(ecp)
      call zero_entcell( ecp )
      call read_entcell_struct( ecp,iu_entstruct )
      
      counter = 0
      do
        counter = counter + 1
        if ( counter >= 1000) 
     &       call stop_model("ent_readcsv: infinite loop",255)
        read(iu_entstruct,*) next
        if (next.eq.'$') then
          exit
        else if (next.eq.'*') then !skip comment
        else if (next.eq.'p') then !new patch
          write(*,*) 'p'
          call insert_patch(ecp,0.d0,0)
          pp => ecp%youngest
          call read_patch_struct(pp,iu_entstruct)
          write(*,*) 'patch area age soil',pp%area,pp%age,pp%soil_type
        else if (next.eq.'c') then !new cohort
          write(*,*) 'c'
          call insert_cohort(pp,0,
     &         0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
     &         fracroot,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
     &         0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
     &         0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)
          cop => pp%shortest
          call read_cohort_struct(cop,iu_entstruct)
          write(*,*) 'cohort pft height',cop%pft,cop%h
        end if
      end do

      call summarize_entcell(ecp)

      end subroutine ent_readcsv


      end module ent_make_struct
