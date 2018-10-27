! Template algorithm for constructing an arbitrary group characterized by matrices of
! whatever kind and dimensions.
       allocate(Op(maxnum_op))
       do i=1,maxnum_op
         call Allocate_Operator(d,Op(i))
       end do
       allocate(Mat(d,d))

       !Construct the list of the generators on top of Op. The identity is always the first operator
       do i=1,ngen
         call Get_Mat_From_Symb_Op(gen(i),Mat,invt)
         if(Err_group) return
         Op(i+1)%Mat=Mat
         Op(i+1)%time_inv=invt
       end do
       ngen=ngen+1
       !
       !Construct the raw Group
       call Get_Group_From_Generators(ngen,Op,multip)
       if(Err_group) return

       !Allocate provisionally to Multip the lattice translations and anti-Translations
       allocate(Lat_tr(d-1,multip), aLat_tr(d-1,multip))
       allocate(centre_coord(d-1))

       call Reorder_Operators(multip, Op, centred, centre_coord, Numops, num_lat, num_alat, Lat_tr, aLat_tr,mag_type)
       if(Err_group) return

       Grp%multip=multip
       Grp%d=d
       if(allocated(Grp%Op)) deallocate(Grp%Op)
       call Allocate_Operators(d,multip,Grp%Op)
       Grp%Op(1:multip)=Op(1:multip)

       !allocate(character(len=40) :: Grp%Symb_Op(multip)) !it doesn't work in gfortran
       if(allocated(Grp%Symb_Op)) Deallocate(Grp%Symb_Op)
       allocate(Grp%Symb_Op(multip))
       do i=1,multip
          Grp%Symb_Op(i)=trim(Symbol_Operator(Op(i)))
       end do

       if(num_lat > 0) then
         if(allocated(Grp%Lat_tr)) Deallocate(Grp%Lat_tr)
         allocate(Grp%Lat_tr(1:d-1,1:Num_Lat))
       end if
       if(num_alat > 0) then
         if(allocated(Grp%aLat_tr)) Deallocate(Grp%aLat_tr)
         allocate(Grp%aLat_tr(1:d-1,1:Num_aLat))
       end if
       if(allocated(Grp%centre_coord)) Deallocate(Grp%centre_coord)
       allocate(Grp%centre_coord(1:d-1))
       Grp%Numops      = Numops
       Grp%centred     = centred
       Grp%mag_type    = mag_type
       Grp%num_lat     = num_lat
       Grp%num_alat    = num_alat
       Grp%centre_coord= centre_coord
       if(num_lat  > 0)  Grp%Lat_tr = Lat_tr(1:d-1,1:Num_Lat)
       if(num_alat > 0)  Grp%aLat_tr=aLat_tr(1:d-1,1:Num_aLat)

