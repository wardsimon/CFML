!!----
!!----
!!----
SubModule (CFML_EoS) EoS_024
   Contains
   
   !!----
   !!---- WRITE_EOS_DATAFILE
   !!----   General routine to Write Data in a Lun iunit
   !!----
   !!---- 17/07/2015
   !!
   Module Subroutine Write_EoS_DataFile(dat,lun)
      !---- Arguments ----!
      type (eos_data_list_type), intent(in) :: dat  ! data structure
      integer,                   intent(in) :: lun  ! Unit to write the information

      !---- Variables ----!
      integer              :: ierr,i,j,k
      character(len=128)   :: text

      !> set up the labels in order
      integer, parameter           :: INI=4, IEND=21
      character(len=5),dimension(INI:IEND) :: lab=['T    ','sigT ','P    ','sigP ','V    ','sigV ',  &
                                                   'A    ','sigA ','B    ','sigB ','C    ','sigC ',  &
                                                   'ALPHA','sigAL','BETA ','sigBE','GAMMA','sigGA']

      !>
      !> assume that unit is connected and open.
      !>

      !> Write header info
      write(unit=lun,fmt='(a,a)',iostat=ierr) 'TITLE ',trim(dat%title)
      write(unit=lun,fmt='(a)',iostat=ierr)    '#'
      write(unit=lun,fmt='(a)',iostat=ierr)    '#  Data file written by CFML eos module'
      write(unit=lun,fmt='(a)',iostat=ierr)    '#'

      !> Crystal system
      if (len_trim(dat%system) > 0)then
         write(unit=lun,fmt='(a,a)',iostat=ierr)  'SYSTEM ',trim(dat%system)
         write(unit=lun,fmt='(a)',iostat=ierr)    '#'
      end if

      !> Original Tscale of data is not known, write out Tscale K if T data present
      if (dat%ic_dat(4) == 1)then
         write(unit=lun,fmt='(a)',iostat=ierr)  'TSCALE  K'
         write(unit=lun,fmt='(a)',iostat=ierr)    '#'
      end if

      !> Scales
      if (len_trim(dat%Pscale_name) > 0)then
         write(unit=lun,fmt='(a,a)',iostat=ierr)  'PSCALE ',trim(dat%Pscale_name)
         write(unit=lun,fmt='(a)',iostat=ierr)    '#'
      end if
      if (len_trim(dat%Vscale_name) > 0)then
         write(unit=lun,fmt='(a,a)',iostat=ierr)  'VSCALE ',trim(dat%Vscale_name)
         write(unit=lun,fmt='(a)',iostat=ierr)    '#'
      end if
      if (len_trim(dat%Lscale_name) > 0)then
         write(unit=lun,fmt='(a,a)',iostat=ierr)  'LSCALE ',trim(dat%Lscale_name)
         write(unit=lun,fmt='(a)',iostat=ierr)    '#'
      end if
      
      !> Datatype: we assume that all data are the same type: responsibility of calling program
      select case(dat%eosd(1)%xtype)
      case(1)
         write(unit=lun,fmt='(a)',iostat=ierr)  'DATATYPE MODULI ISOTHERMAL'
         write(unit=lun,fmt='(a)',iostat=ierr)    '#'          
      case(2)
         write(unit=lun,fmt='(a)',iostat=ierr)  'DATATYPE MODULI ADIABATIC'
         write(unit=lun,fmt='(a)',iostat=ierr)    '#'        
      end select
      
      
      !> build format line
      text='FORMAT 1'
      do i=ini,iend
         if (dat%ic_dat(i) ==1) text=trim(text)//' '//lab(i)
      end do
      write(unit=lun,fmt='(a)',iostat=ierr)    trim(text)

      !> write the data: all data written out
      do i=1,dat%n        ! loop over data points
         if (dat%ic_dat(4) == 1) write(unit=lun,fmt='(f12.6)',iostat=ierr,advance='no') dat%eosd(i)%t
         if (dat%ic_dat(5) == 1) write(unit=lun,fmt='(f12.6)',iostat=ierr,advance='no') dat%eosd(i)%sigt
         if (dat%ic_dat(6) == 1) write(unit=lun,fmt='(f12.6)',iostat=ierr,advance='no') dat%eosd(i)%p
         if (dat%ic_dat(7) == 1) write(unit=lun,fmt='(f12.6)',iostat=ierr,advance='no') dat%eosd(i)%sigp
         if (dat%ic_dat(8) == 1) write(unit=lun,fmt='(f12.6)',iostat=ierr,advance='no') dat%eosd(i)%v
         if (dat%ic_dat(9) == 1) write(unit=lun,fmt='(f12.6)',iostat=ierr,advance='no') dat%eosd(i)%sigv

         !> small loop over cell edges
         do j=1,3
            k=8+2*j
            if (dat%ic_dat(k) == 1)   write(unit=lun,fmt='(f12.6)',iostat=ierr,advance='no') dat%eosd(i)%cell(j)
            if (dat%ic_dat(k+1) == 1) write(unit=lun,fmt='(f12.6)',iostat=ierr,advance='no') dat%eosd(i)%sigc(j)
         end do

         !> small loop over cell angles
         do j=1,3
            k=14+2*j
            if (dat%ic_dat(k) == 1)   write(unit=lun,fmt='(f12.6)',iostat=ierr,advance='no') dat%eosd(i)%ang(j)
            if (dat%ic_dat(k+1) == 1) write(unit=lun,fmt='(f12.6)',iostat=ierr,advance='no') dat%eosd(i)%siga(j)
         end do
         write(unit=lun,fmt='(1x)',iostat=ierr)          ! forces new line
      end do
   End Subroutine Write_Eos_Datafile

   !!----
   !!---- WRITE_EOS_FILE
   !!----    General routine to Write EoS in a Lun iunit
   !!----
   !!---- 17/07/2015
   !!
   Module Subroutine Write_Eos_File(Eos,Lun)
      !---- Arguments ----!
      type (EoS_Type),intent(in)   :: Eos ! EoS object
      integer,intent(in)           :: lun ! Unit

      !---- Variables ----!
      character(len=12)            :: stext
      character(len=512)           :: text
      integer                      :: ierr,i,j
      real(kind=cp)                :: valp

      !>
      !> assume that unit is connected and open.
      !>

      !> Header info written to file from calling program
      write(unit=lun,fmt='(a)',iostat=ierr) ' '
      write(unit=lun,fmt='(a)',iostat=ierr) ' EosFit parameter file'
      write(unit=lun,fmt='(a)',iostat=ierr) ' '
      write(unit=lun,fmt='(a)',iostat=ierr) ' This is a fixed-format file. If you edit this file, make sure you do not move '// &
                                            'anything!'
      write(unit=lun,fmt='(a)',iostat=ierr) ' It is safer to change the parameters by loading the file to EosFit, and saving '// &
                                            'the file after making changes'
      write(unit=lun,fmt='(a)',iostat=ierr) ' _______________________________________________________________________________'// &
                                            '_____________________________'
      write(unit=lun,fmt='(a)',iostat=ierr) ' '


      !> title and comments
      write(unit=lun,fmt='(a,a)',iostat=ierr) 'Title =',trim(eos%title)
      write(unit=lun,fmt='(a,a)',iostat=ierr) 'Savedate =',trim(eos%savedate)

      do i=1,size(eos%doc)
         if (len_trim(eos%doc(i)) > 0)then
            write(unit=lun,fmt='(a)') 'Comment ='//trim(eos%doc(i))
         end if
      end do
      write(unit=lun,fmt='(a)',iostat=ierr) ' '

      !> eos type
      text=',  ('//trim(eos%model)//')'
      if (eos%imodel == 0) text=',  (none)'
      write(unit=lun,fmt='(a,i3,a)',iostat=ierr) 'Model =',eos%imodel,text
      write(unit=lun,fmt='(a,i3)',iostat=ierr) 'Order =',eos%iorder

      text=',  ('//trim(eos%tmodel)//')'
      if (eos%itherm == 0)text=',  (none)'
      write(unit=lun,fmt='(a,i3,a,a,a)',iostat=ierr) 'Thermal =',eos%itherm,text

      text=',  ('//trim(eos%cmodel)//')'
      if (eos%icross == 0)text=',  (none)'
      write(unit=lun,fmt='(a,i3,a,a,a)',iostat=ierr) 'Cross =',eos%icross,text

      text=',  ('//trim(eos%tranmodel)//')'
      if (eos%itran == 0)text=',  (none)'
      write(unit=lun,fmt='(a,i3,a,a,a)',iostat=ierr) 'Trans =',eos%itran,text

      text=',  ('//trim(eos%smodel)//')'
      if (eos%ishear == 0)text=',  (none)'
      write(unit=lun,fmt='(a,i3,a,a,a)',iostat=ierr) 'Shear =',eos%ishear,text

      if (eos%linear)then
         write(unit=lun,fmt='(a)',iostat=ierr) 'Type = Linear'
      else
         write(unit=lun,fmt='(a)',iostat=ierr) 'Type = Volume'
      end if
      write(unit=lun,fmt='(a,f10.5)',iostat=ierr) 'Pref =',eos%pref
      write(unit=lun,fmt='(a,f10.5)',iostat=ierr) 'Tref =',eos%tref
      write(unit=lun,fmt='(a,a)',iostat=ierr) 'Pscale =',trim(eos%pscale_name)
      write(unit=lun,fmt='(a,a)',iostat=ierr) 'Vscale =',trim(eos%vscale_name)
      write(unit=lun,fmt='(a,f10.5)',iostat=ierr) 'Stoich =',eos%stoich
      if (eos%density0 > tiny(0.)) then
         write(unit=lun,fmt='(a,f10.5)',iostat=ierr) 'Density0 =',eos%density0
      end if

      write(unit=lun,fmt='(a)',iostat=ierr) ' '

      !> Eos parameters
      do i=1,n_eospar
         valp=eos%params(i)*eos%factor(i)
         if(abs(valp) < 1.0E7_cp)then
            text=string_real(valp,precision(valp)+2)
         else
            write(text,'(''    Inf'')')
         end if

         if (eos%iuse(i) == 0) then
            write(unit=lun,fmt='(a,i2,a12,5a)')'Param =',i,text(1:12)
         else
            write(unit=lun,fmt='(a,i2,a12,5a)')'Param =',i,text(1:12),'     (',eos%parname(i),',  ',&
                 trim(eos%comment(i)),')'
         end if
      end do

      !> VCV: stored as scaled values for precision
      write(unit=lun,fmt='(a)',iostat=ierr) ' '
      write(unit=lun,fmt='(a)',iostat=ierr) '  Variance-Covariance matrix='
      do i=1,n_eospar
         text=''
         do j=1,n_eospar
            if (abs(eos%vcv(i,j)) > tiny(0.0))then
               write(stext,fmt='(e12.5)')eos%vcv(i,j )*eos%factor(i)*eos%factor(j)
               text=trim(text)//stext
            else              !   123456789012
               text=trim(text)//' 0.00000E+00'
            end if
         end do
         write(unit=lun,fmt='(a)',iostat=ierr)trim(text)
      end do
      write(unit=lun,fmt='(a)',iostat=ierr) ' '
      write(unit=lun,fmt='(a)',iostat=ierr) ' '
   End Subroutine Write_Eos_File

   !!----
   !!---- WRITE_EOSCAL
   !!----   Subroutine to write the calculated parameters of an eos to file at a series of PT points
   !!----   NO  header info is printed here.
   !!----
   !!----   Therefore the program header and write_info_eos have to be called first before calling this routine
   !!----   Then write_eoscal_header is called from here
   !!----
   !!----   Change: 06/10/2015 to make write_eoscal_header private, and change name from write_eoscal_file
   !!----   Change: 12/12/2017 created eoscal_text so that errors and values are printed when error state
   !!----   Change: 19/12/2018 added error flag to return to calling program, if warning or error on at least one calc
   !!---- 
   !!---- 17/07/2015
   !!
   Module Subroutine Write_Eoscal(Pmin,Pmax,Pstep,Tmin,Tmax,Tstep,Tscale_In,Eos,Lun,Nprint,eoscal_err)
      !---- Arguments ----!
      real(kind=cp),    intent(in)  ::  pmin, pmax, pstep   !P to calculate properties
      real(kind=cp),    intent(in)  ::  tmin,tmax,tstep     !T to calculate properties
      character(len=*), intent(in)  ::  tscale_in           ! Name of the Tscale for output, either C or K
                                                            ! If Pstep or Tstep  < tiny(0.0) then only Pmin (or Tmin) calculated
      type(EoS_Type),   intent(in)  ::  eos                 ! Eos
      integer,          intent(in)  :: lun                  ! logical unit for printing
      integer,          intent(out) :: nprint               ! Number of data printed
      logical,          intent(out) :: eoscal_err           ! error flag

      !---- Local variable ----!
      real(kind=cp)           :: p,t      ! The p and T of each calculation
      real(kind=cp)           :: pst,tst  ! local copy of tstep and pstep
      character(len=255)      :: text     ! local text variable
      character(len=1)        :: tscale   ! local name of tscale
      logical                 :: loop_p   ! loop indicator .true. for inner loop of calcs over P
      integer,dimension(19)   :: ip=[6,6,9,8,6,5,5,9,7,7,5,9,7,7,6,6,6,6,6] ! format for output
      integer                 :: i

      real(kind=cp),dimension(6) :: parvals(7)
      real(kind=cp),dimension(6) :: esd
      real(kind=cp),dimension(19):: parout,esdout
      real(kind=cp)              :: v0,fp,fs,agt

      !> init
      nprint=0    ! output counter
      eoscal_err=.false.

      !> Tscale for output: C or K
      if (len_trim(tscale_in) == 0)then
         tscale='K'
      else
         tscale=U_case(tscale_in)
         if (tscale /= 'K' .and. tscale /='C')tscale='K'
      end if

      !> Write file header
      call write_eoscal_header(eos,lun,tscale)

      !> copy Pstep/Tstep
      tst=tstep
      pst=pstep

      !> set up loop control variables
      if (abs(pst) > tiny(0.0))then       ! inner loop over P
         loop_p=.true.
         if (abs(tst) < tiny(0.0))then   ! no outerloop
            tst=10.*max((tmax-tmin),1.0)       ! set tstep big enough to stop loop
         end if
      else
         loop_p=.false.                  ! inner loop over T
         if (abs(pst) < tiny(0.0))then   ! no outerloop
            pst=10.*max((pmax-pmin),1.0)       ! set pstep big enough to stop loop
         end if
      end if

      !> Initialise loop variables
      p=pmin
      t=tmin

      !> Start of outer loop
      outer: do
         !> reset inner loop variable to start value
         if (loop_p)then
            p=pmin
         else
            t=tmin
         end if

         inner: do
             call clear_error()
             call physical_check(eos,Pin=p,Tin=T)
             if (err_CFML%IErr==1) then
                text=trim(string_real(p,6))//'  '//trim(string_real(t,6))//' :   '//trim(err_CFML%Msg)
                write(lun,'(a)')trim(text)
                eoscal_err=.true.
             else
                text=eoscal_text(p,t,Tscale_In,Eos)
                write(lun,'(a)')trim(text)      ! This way we get to see the calculated values even if error in calcs with valid eos
                if (err_CFML%IErr==1)then
                    write(lun,'(a)')'   *****WARNING:   '//trim(err_CFML%Msg)
                    eoscal_err=.true.
                endif
                
             endif
            nprint=nprint+1

            !> Now increment inner loop variable and test for completion
            if (loop_p) then
               p=p+pst                           ! inner loop over p
               if (p > pmax+0.99_cp*pst)exit inner
            else
               t=t+tst                           ! inner loop over t
               if (t > tmax+0.99_cp*tst)exit inner
            end if
         end do inner

         write(lun,'(/)')        ! blank line to help some plotting programs

         !> Now increment outer loop variable and test for completion
         if (loop_p)then
            t=t+tst                           ! outer loop over t
            if (t > tmax+0.99_cp*tst)exit outer
         else
            p=p+pst                           ! outer loop over p
            if (p > pmax+0.99_cp*pst)exit outer
         end if
      end do outer
   End Subroutine Write_Eoscal
   
   !!--++
   !!--++ WRITE_EOSCAL_HEADER
   !!--++
   !!--++ PUBLIC
   !!--++ Subroutine that print information on iout unit
   !!--++
   !!--++ 17/07/2015
   !!
   Module Subroutine Write_Eoscal_Header(Eos,Lun,Tscale_In)
      !---- Arguments ----!
      type(EoS_Type),intent(in)   :: eos         ! Eos information
      integer,       intent(in)   :: lun         ! logical unit for printing
      character(len=*),intent(in) :: tscale_in   ! Scale for Temp

      !---- Local Variables ----!
      character(len=1)     :: tscale
      character(len=255)   :: head     ! local text variable for column headers

      !> Warning for Tait or Murnaghan
      write(lun,'(//)')
      if (eos%itherm /= 0) then
         write(lun,'("  Note that values of alpha are multiplied by a factor of ",f5.1,"x10^5"//)')eos%alphafactor/1.0E5
      end if

      if (eos%imodel == 1 .or. eos%imodel == 5 .or. eos%imodel == 6) then
         write(lun,'("  Do not forget: Normalised Pressure and strain not defined for ",a," Eos")')trim(eos%model)
      else if(eos%itherm /= 0)then
         write(lun,'("  Normalised Pressure (NP) and finite strain (f) are defined relative to V at P=0 and same T")')
      else
         write(lun,'("  Normalised Pressure is NP and finite strain is f")')
      end if

      !> tscale for output: C or K
      if (len_trim(tscale_in) == 0) then
         tscale='K'
      else
         tscale=U_case(tscale_in)
         if (tscale /= 'K' .and. tscale /='C')tscale='K'
      end if

      !> create column header
      if (eos%linear) then
         write(head,'(" Press   Temp",a1,"     Length  esdL       L/L0T  esd(L/L0T)  M      esdM  Mprime ", &
             &"esdMp   Mpp  esdMpp  f         esdf       NP      esdNP   dM/dT    esddK   alpha  esda")' ) Tscale
      else
         write(head,'(" Press   Temp",a1,"     Volume  esdV       V/V0T  esd(V/V0T)  K      esdK  Kprime ", &
             &"esdKp   Kpp  esdKpp  f         esdf       NP      esdNP   dK/dT    esddK   alpha  esda")' ) Tscale
      end if
      if (eos%itran > 0)head=trim(head)//'   spstrain'

      if (eos%density0 > tiny(0.0)) head=trim(head)//'  density  esdden'

      if (eos%pthermaleos .and. eos%itran ==0)head=trim(head)//'  Ptherm'

      if (eos%itherm == 7)head=trim(head)//' MGDGamma DebyeT'

      if (eos%itherm > 0 .and. abs(eos%params(18)) > tiny(0.0)) then
         if (eos%linear)then
            head=trim(head)//'    Ms    Gamma'
         else
            head=trim(head)//'    Ks    Gamma'
         end if
      end if

      !> Write header
      write(lun,'(/a)')trim(head)
   End Subroutine Write_Eoscal_Header
   
   !!----
   !!---- WRITE_INFO_EOS
   !!----    Subroutine that print information on iout unit
   !!----
   !!---- 17/07/2015
   !!
   Module Subroutine Write_Info_Eos(Eospar,iout)
      !---- Arguments ----!
      type(Eos_Type),    intent(in) :: Eospar  ! EoS object
      integer, optional, intent(in) :: iout    ! Logical unit

      !---- Local Variables ----!
      character(len=30) :: line,string
      integer           :: i,lun

      !> Unit to print the information
      lun=6
      if (present(iout)) lun=iout

      !> Header / Title
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '  EOS Information'
      write(unit=lun,fmt='(a)') '-------------------'
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '    Title: '//trim(eospar%title)
      write(unit=lun,fmt='(a)') 'Eos Saved: '//trim(eospar%savedate)
      !> Doc Information
      do i=1,size(eospar%doc)
          if (len_trim(eospar%doc(i)) > 0) then
             write(unit=lun,fmt='(a)') '  Comment: '//trim(eospar%doc(i))
          end if
      end do
      write(unit=lun,fmt='(a)') ' '
      if (eospar%imodel /= 0) then
         if (len_trim(eospar%Pscale_name) > 0)write(unit=lun,fmt='(a)') '   Pscale: '//trim(eospar%Pscale_name)
         if (len_trim(eospar%Vscale_name) > 0)write(unit=lun,fmt='(a)') '   Vscale: '//trim(eospar%Vscale_name)
         write(unit=lun,fmt='(a,t27,f8.3)') '   Stoichiometry: ',eospar%stoich


         !> Reference Density
         if (eospar%density0 > tiny(0.0)) then
            write(unit=lun,fmt='(a,t27,f8.3)') '   Reference density: ',eospar%density0
         end if

         write(unit=lun,fmt='(a)') ' '
         write(unit=lun,fmt='(a)') '  Compressibility'
         write(unit=lun,fmt='(a)') '-------------------'
         write(unit=lun,fmt='(a)') ' '
         write(unit=lun,fmt='(a)') '   Model: '//trim(eospar%model)
         if (eospar%imodel > 0) then                ! no output if no p eos: all done in write_info_eos_thermal
            write(unit=lun,fmt='(a,i2)') '   Order: ',eospar%iorder
            if (eospar%linear) then
               write(unit=lun,fmt='(a)') '   Class: Linear'
            else
               write(unit=lun,fmt='(a)') '   Class: Volume'
            end if

            !> Pressure Parameters
            write(unit=lun,fmt='(a,t27,f8.3)') '   Pressure of reference: ',eospar%pref
            write(unit=lun,fmt='(a)') ' '

            do i=1,4
               if (eospar%iuse(i) /= 0) then
                  line=string_numstd(eospar%params(i)*eospar%factor(i),eospar%esd(i)*eospar%factor(i))
                  string=' '
                  select case(eospar%iuse(i))
                     case (2)
                        string=' [FIXED VALUE]'
                     case (3)
                        string=' [IMPLIED VALUE]'
                  end select
                  write(unit=lun,fmt='(3x,a5,": ",a,T30,":",a)') &
                       trim(eospar%parname(i)),trim(line),trim(eospar%comment(i))//trim(string)
               end if
            end do
            write(unit=lun,fmt='(a)') ' '
         end if
      end if

      !> Thermal EOS
      if (eospar%itherm > 0) call write_info_eos_thermal(eospar,lun)

      !> Cross terms
      if (eospar%icross > 0) call write_info_eos_cross(eospar,lun)

      !> Transition
      if (eospar%itran > 0) call write_info_eos_transition(eospar,lun)

      !> Shear
      if (eospar%ishear > 0) call write_info_eos_shear(eospar,lun)

      !> End
      write(unit=lun,fmt='(a)') ' '
   End Subroutine Write_Info_Eos

   !!--++
   !!--++ WRITE_INFO_EOS_CROSS
   !!--++
   !!--++ PRIVATE
   !!--++ Subroutine that print information on iout unit
   !!--++
   !!--++ Date: 11/07/2016
   !!
   Module Subroutine Write_Info_Eos_Cross(Eos,Iout)
      !---- Arguments ----!
      type(Eos_Type),    intent(in) :: Eos
      integer, optional, intent(in) :: iout

      !---- Local Variables ----!
      character(len=30) :: line,string
      integer           :: i,lun

      !> Check
      if (Eos%icross ==0) return

      !> Unit to print the information
      lun=6
      if (present(iout)) lun=iout

      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '      P-T cross-terms'
      write(unit=lun,fmt='(a)') '---------------------------'
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '   Model: '//trim(eos%cmodel)

      do i=5,6
         if (eos%iuse(i) /= 0) then
            line=string_numstd(eos%params(i)*eos%factor(i),eos%esd(i)*eos%factor(i))     ! include scaling
            string=' '
            select case(eos%iuse(i))
               case(2)
                  string=' [FIXED VALUE]'
               case(3)
                  string=' [IMPLIED VALUE]'
            end select
            write(unit=lun,fmt='(3x,a5,'': '',a,T30,'':'',a)') &
                 trim(eos%parname(i)),trim(line),trim(eos%comment(i))//trim(string)
         end if
      end do
      write(unit=lun,fmt='(a)') ' '
   End Subroutine Write_Info_Eos_Cross

   !!--++
   !!--++ WRITE_INFO_EOS_SHEAR
   !!--++
   !!--++ PRIVATE
   !!--++ Subroutine that print information on iout unit
   !!--++
   !!--++ 11/07/2016
   !!
   Module Subroutine Write_Info_Eos_Shear(Eos,Iout)
      !---- Arguments ----!
      type(Eos_Type),    intent(in) :: Eos
      integer, optional, intent(in) :: iout

      !---- Local Variables ----!
      character(len=30) :: line,string
      integer           :: i,lun

      !> Check
      if(Eos%ishear ==0) return

      !> Unit to print the information
      lun=6
      if (present(iout)) lun=iout

      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '  Shear modulus Information'
      write(unit=lun,fmt='(a)') '---------------------------'
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '   Model: '//trim(eos%smodel)

      do i=30,39
         if (eos%iuse(i) /= 0) then
            line=string_numstd(eos%params(i)*eos%factor(i),eos%esd(i)*eos%factor(i))     ! include scaling
            string=' '
            select case(eos%iuse(i))
                case(2)
                string=' [FIXED VALUE]'
                case(3)
                string=' [IMPLIED VALUE]'
            end select
            write(unit=lun,fmt='(3x,a5,'': '',a,T30,'':'',a)') &
                 trim(eos%parname(i)),trim(line),trim(eos%comment(i))//trim(string)
         end if
      end do
      write(unit=lun,fmt='(a)') ' '
   End Subroutine Write_Info_Eos_Shear

   !!--++
   !!--++ WRITE_INFO_EOS_THERMAL
   !!--++
   !!--++ PRIVATE
   !!--++ Subroutine that print information on iout unit
   !!--++
   !!--++ 17/07/2015
   !!
   Module Subroutine Write_Info_Eos_Thermal(Eospar,iout)
      !---- Arguments ----!
      type(Eos_Type),    intent(in) :: Eospar  ! EoS object
      integer, optional, intent(in) :: iout    ! Logical unit

      !---- Local Variables ----!
      character(len=30) :: line,string
      integer           :: i,lun,is

      !> Init
      lun=6
      if (present(iout)) lun=iout

      is=10                            ! If pmodel present, it was already reported
      if (eospar%imodel == 0) is=1     ! if no pmodel report all params here


      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '  Thermal Expansion'
      write(unit=lun,fmt='(a)') '-------------------'
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '   Model: '//trim(eospar%tmodel)
      write(unit=lun,fmt='(a)') ' '

      write(unit=lun,fmt='(a,f8.2,a)') '   Temperature of reference: ',eospar%tref,' K'
      write(unit=lun,fmt='(a)') ' '
      do i=is,19
         if (eospar%iuse(i) /= 0) then
            line=string_numstd(eospar%params(i)*eospar%factor(i),eospar%esd(i)*eospar%factor(i))     ! include scaling
            string=' '
            select case(eospar%iuse(i))
                case(2)
                string=' [FIXED VALUE]'
                case(3)
                string=' [IMPLIED VALUE]'
            end select
            write(unit=lun,fmt='(3x,a5,": ",a,T30,":",a)') &
                  trim(eospar%parname(i)),trim(line),trim(eospar%comment(i))//trim(string)
         end if
      end do
   End Subroutine Write_Info_Eos_Thermal

   !!--++
   !!--++ WRITE_INFO_EOS_TRANSITION
   !!--++
   !!--++ PRIVATE
   !!--++ Subroutine that print information on iout unit
   !!--++
   !!--++ 17/07/2015
   !!
   Module Subroutine Write_Info_Eos_Transition(Eospar,iout)
      !---- Arguments ----!
      type(Eos_Type),    intent(in) :: Eospar  ! EoS object
      integer, optional, intent(in) :: iout    ! Logical Unit

      !---- Local Variables ----!
      character(len=30) :: line,string
      integer           :: i,lun

      !> Check
      if(eospar%itran == 0)return

      !> Init
      lun=6
      if (present(iout)) lun=iout

      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '  Phase Transition '
      write(unit=lun,fmt='(a)') '-------------------'
      write(unit=lun,fmt='(a)') ' '
      write(unit=lun,fmt='(a)') '   Model: '//trim(eospar%tranmodel)
      write(unit=lun,fmt='(a)') ' '

      do i=20,29
         if (eospar%iuse(i) /= 0) then
            line=string_numstd(eospar%params(i)*eospar%factor(i),eospar%esd(i)*eospar%factor(i))     ! include scaling
            string=' '
            select case(eospar%iuse(i))
               case (2)
                  string=' [FIXED VALUE]'
               case (3)
                  string=' [IMPLIED VALUE]'
            end select
            write(unit=lun,fmt='(3x,a5,": ",a,T30,":",a)') &
                  trim(eospar%parname(i)),trim(line),trim(eospar%comment(i))//trim(string)
         end if
      end do

   End Subroutine Write_Info_Eos_Transition
   
   !!----
   !!---- EOSCAL_TEXT
   !!----   Subroutine to write the calculated parameters of an eos to file at one PT point
   !!----   NO  header info is printed here.
   !!----
   !!----   Normally called from Write_eoscal
   !!----   written 12/2017 by extracting code from write_eoscal
   !!----   in order to allow values and error messages to be written
   !!----
   !!---- 31/05/2019
   !!
   Module Function Eoscal_Text(P,T,Tscale_In,Eos) Result(text)
      !---- Arguments ----!
      real(kind=cp),    intent(in)  :: p                   !P to calculate properties
      real(kind=cp),    intent(in)  :: t                   !T to calculate properties
      character(len=*), intent(in)  :: tscale_in           ! Name of the Tscale for output, either C or K
      type(EoS_Type),   intent(in)  :: eos                 ! Eos
      character(len=:), allocatable :: text                ! character string with results

      !---- Local variable ----!
      character(len=1)        :: tscale                                     ! local name of tscale
      integer,dimension(19)   :: ip=[6,6,9,8,6,5,5,9,7,7,5,9,7,7,6,6,6,6,6] ! format for output
      integer                 :: i

      real(kind=cp),dimension(6) :: parvals(7)
      real(kind=cp),dimension(6) :: esd
      real(kind=cp),dimension(19):: parout,esdout
      real(kind=cp)              :: v0,fp,fs,agt

      !> Init
      text="  "
      
      !> Tscale for output: C or K
      if (len_trim(tscale_in) == 0)then
         tscale='K'
      else
         tscale=U_case(tscale_in)
         if (tscale /= 'K' .and. tscale /='C')tscale='K'
      end if

      call clear_error()
      esd=0.0_cp
      esdout=0.0_cp
      parout=0.0_cp

      !> Now do the calculations at P,T
      call EoS_Cal(P,T,eos,Parvals)    ! GET V,K ETC
      if (sum(eos%vcv) > tiny(0.0)) CALL eos_cal_esd(P,T,eos,esd)

      !> build ouput value array
      V0=Get_Volume(0.0,T,Eos)

      parout(1)=p
      parout(2)=t
      if (tscale =='C')parout(2)=parout(2)-273.16
      parout(3)=parvals(1)*eos%factor(1)      ! v
      esdout(3)=esd(1)*eos%factor(1)
      parout(4)=parvals(1)/v0                 ! v/V0 at this T
      esdout(4)=esdout(3)/v0

      !> convert  V,K,Kp,Kpp to output values
      do i=2,4
         parout(i+3)=parvals(i)*eos%factor(i)
         esdout(i+3)=esd(i)*eos%factor(i)
      end do

      !>deal with f-F:
      if (abs(p) < tiny(0.0) ) then
         call ffcal_eos(p,t,eos,fp,fs)      ! because F not defined numerically at P=0
         parout(9)=FP
         esdout(9)=0.0_cp
      else
         call ffcal_dat_esd(parvals(1),esd(1),V0,0.0_cp,P,0.0_cp,Eos, &          ! only esd input is esd(V) at this P
              parout(9),esdout(9),parout(8),esdout(8))
      end if

      !> dK/dT
      parout(10)=parvals(5)*eos%factor(5)
      esdout(10)=esd(5)*eos%factor(5)

      !> handle alpha
      parout(11)=parvals(6)*eos%alphafactor
      esdout(11)=esd(6)*eos%alphafactor

      !> spon strain
      if (eos%itran > 0) parout(12)=Get_Transition_Strain(P,T,Eos)

      !> density
      if (eos%density0 > tiny(0.0)) then
         parout(13)=eos%density0*eos%params(1)/parvals(1)
         parout(14)=parout(13)*esd(1)/parvals(1)
      end if

      !> Thermal pressure
      if (eos%pthermaleos .and. eos%itran ==0) parout(15)=p-get_pressure(parvals(1),eos%tref,eos)

      !>MGD EoS parameters
      if (eos%itherm == 7) then
          parout(16)=get_grun_v(parvals(1),eos)      ! Gruneisen gamma
          parout(17)=get_DebyeT(parvals(1),eos)      !Debye T
      end if

      !> Report adiabatic properties
      if (eos%itherm > 0 .and. abs(eos%params(18)) > tiny(0.0)) then
         parout(19)=Get_Grun_V(parvals(1),eos)             !Gruneisen for Kt--> Ks
         agt=parvals(6)*parout(19)*T        ! Get_Grun knows about linear
         if (eos%linear) agt=3.0_cp*agt
         parout(18)=(1.0_cp+agt)*parvals(2) ! Ks/Ms
      end if

      !> output this datum: dynamic formatting to text string
      !> pressure (no esd)
      text=trim(string_real(parout(1),ip(1)))

      !> T value (no esd)
      text=trim(text)//'  '//trim(string_real((parout(2)),ip(2)))

      !> other params
      do i=3,11
         text=trim(text)//'  '//trim(string_real(parout(i),ip(i)))//' '//&
              trim(string_real(esdout(i),ip(i)))
      end do

      if (eos%itran > 0) text=trim(text)//'  '//trim(string_real(parout(12),ip(12)))    ! spontaneous strain

      if (eos%density0 > tiny(0.0)) &
         text=trim(text)//'  '//trim(string_real(parout(13),ip(13)))//' '//&
              trim(string_real(parout(14),ip(14)))

      if (eos%pthermaleos .and. eos%itran ==0) text=trim(text)//'  '//&
              trim(string_real(parout(15),ip(15)))

      if (eos%itherm == 7) text=trim(text)//'  '//trim(string_real(parout(16),ip(16)))//'  '//&
              trim(string_real(parout(17),ip(17)))

      if (eos%itherm > 0 .and. abs(eos%params(18)) > tiny(0.)) &
              text=trim(text)//'  '//trim(string_real(parout(18),ip(18)))//'  '//&
              trim(string_real(parout(19),ip(19)))

   End Function Eoscal_text
   
End SubModule EoS_024   