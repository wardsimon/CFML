!!----
!!---- Program: SPG_INFO
!!----          Example of simple program using CFML
!!----
!!---- Author: Juan Rodriguez-Carvajal
!!---- Revision: Nov-2008
!!
! => Enter a space group: 
! => Space Group (HM/Hall symbol or number): P n m a
! => Enter a a k-vector: 0.321 0 0.5


Program SSPG_Info
   !---- Use Modules ----!
   use CFML_Crystallographic_Symmetry, only: Space_Group_Type, set_SpaceGroup, &
                                             Write_SpaceGroup


   use CFML_Propagation_Vectors_TEST,       only: Group_k_Type, Set_Gk, K_Star, & 
                                             Write_Group_k, k_EQUIV  

   use CFML_GlobalDeps,                only: cp

   use CFML_Math_General, only: Determinant, Equal_Matrix, Equal_Vector, Modulo_Lat

   !---- Variables ----!
   implicit none
   integer                    :: i, j, m, n, sizegg, iter, nmax, k , t1, t2
   !real                       :: t

   !--- Variables of the original example:
   !
   !--- To prompt the user: "Enter a space group" 
   !
   character(len=20)      :: spgr
   type(Space_Group_type) :: grp_espacial


   !--- Variables for the new part pf the code:
   !
   !--- A vector to compute the K_Star() of entered space group.
   !    Another vector for its nint(), just in case.
   !
   real(kind=cp), dimension(3) :: vec_k
   real(kind=cp), dimension(3) :: vec_k_int

   real(kind=cp), dimension(3)   :: vec_k_prime
   real(kind=cp), dimension(3)   :: vec_k_condition1
   real(kind=cp), dimension(3)   :: vec_k_condition2
   real(kind=cp), dimension(3)   :: vec_H_test

   real, dimension(1) :: cond1 = [ 1 ]
   real, dimension(1) :: cond2 = [ 2 ]

   real, dimension(6) :: m_values = (/0,1,2,3,4,6/) ! 5 is excluded, hace otro n_values
   real, dimension(4) :: n_values = (/1,4,2,4/)! los indices de matrix de solo Gk del ejemplo PNMA (4 elementos), para probar; hacer subrutina recursiva de verdad porque estos son corectos per FIJOS
!see http://matrix.reshish.com/powCalculation.php
   real, dimension(6) :: t4_possible_values_1 ! todo el array de m_values dividido por indice de la matriz 1
   real, dimension(6) :: t4_possible_values_2 ! todo el array de m_values dividido por indice de la matriz 2

   !real, dimension(4) :: t4_values = (/                        ! {m_values/n_values} 
   !                                     0/1, 0/4 , 0/2 , 0/4 
   !                                     1/1, 1/4 , 1/2 , 1/4 
   !                                     2/1, 2/4 , 2/2 , 2/4 
   !                                     3/1, 3/4 , 3/2 , 3/4 
   !                                     4/1, 4/4 , 4/2 , 4/4 
   !                                     6/1, 6/4 , 6/2 , 6/4    /) 
   !
   !                                     0 0    0   0
   !                                     1 0.25 0.5 0.25
   !                                     2 0.5  1   0.5
   !                                     3 0.75 1.5 0.75
   !                                     4 1    2   1
   !                                     6 1.5  3   1.5
   !
   !                       eS DECIR: 0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, LUEGO:
   real, dimension(10) :: t4_values_totales = (/ 0.0, 0.25, 0.5, 0.75, 1., 1.5, 2., 3., 4., 6. /)
         

   real(kind=cp), dimension(4)     :: vec_H_R  ! a buffer vec, therea aer many, one for operator, allocat
   !real(kind=cp), dimension(4:4)   :: mat_R_S   ! a buffer matrix, therea aer many, one for operator, NOT USED

   real, dimension(:,:,:),   allocatable :: matrices_R_S   ! The {H(R)} que seran (4:4:NUMERO OPERADORES)
   real, dimension(:,:),     allocatable :: vectors_T      ! The translations
   real, dimension(:,:),     allocatable :: P ! a test matrix for multiplication

   real, dimension(4) :: vec_k_buffer_1 ! make allocatable, init to 0
   real, dimension(4) :: vec_k_buffer_2
   real, dimension(4) :: vec_k_buffer_1_sin 
   real, dimension(4) :: vec_k_buffer_2_sin
   
   real, dimension(4) :: vec_SSG_law ! 
   
   logical :: eq
   character(len=80) :: fmtt

   real(kind=cp), dimension(4)   :: vec_modulo_lat


   !
   !--- Declaring variables for Space Group and Group k
   !
   Type (Space_Group_Type) :: SPGk_in, SPGk_out
   Type (Group_k_Type)     :: Gk_out




   !---- Procedure ----!


   do

      WRITE(*,*) " "
      WRITE(*,*) " ============= STARTING CORRECT EXAMPLE: ========"
      WRITE(*,*) " "

      write(unit=*,fmt="(a)") " => Enter a space group: "
      write(unit=*,fmt="(a)",advance="no") " => Space Group (HM/Hall symbol or number): "
      read(unit=*,fmt="(a)") spgr
      if (len_trim(spgr)==0) exit
      write(unit=*,fmt="(a)",advance="no") " => Enter a a k-vector: "
      read(unit=*,fmt=*) vec_k
     
      !> Setting the Space Group Information and 
      !> writing the SpaceGroup Information.
      call set_spacegroup(spgr,grp_espacial)
      call Write_SpaceGroup(grp_espacial, full=.true.)

      WRITE(*,*) " "
      WRITE(*,*) " ============= END OF CORRECT EXAMPLE. ========"
      WRITE(*,*) " "

      WRITE(*,*) " "
      WRITE(*,*) " ============= HIC SUNT DRACONES: ============="
      WRITE(*,*) " "

      !> Computing the K-Star (i.e. the small group or Gk)
      !> we provide a fixed vector by now: K=(1,1,1)    for instance
      !vec_k     = (/ 1, 1, 1 /)
      !vec_k_int = nint(vec_k)


      !> Recast the SpaceGroup variable, to decouple this new part of the code
      !> from the original 'space_group_info' example.
      !> So, we just rename the space group input by the user as 'SPGk_in':
      SPGk_in=grp_espacial


      WRITE(*,*) " "
      WRITE(*,*) " =============== SMALL GROUP IS: =============="
      WRITE(*,*) " "


      !> Computing the small group. 
      !>The output will be loaded into Gk_out so we have to write it out:
      !call K_Star(vec_k_int, SPGk_in, Gk_out)
      !call Write_Group_k(Gk_out)

      WRITE(*,*) " "
      WRITE(*,*) " =============== CHECKING: ==================="
      WRITE(*,*) " "

      !> Checking: 
      !> convert the propagation vector group to a conventional space group
      !> and write it out. It has to be the initial (user-prompted) space group.

      WRITE(*,*) " "
      WRITE(*,*) " ================== END OF CHECKING. =========="
      WRITE(*,*) " "

      WRITE(*,*) " "
      WRITE(*,*) " ================== SMALL GROUP EXTENDED: ====="
      WRITE(*,*) " "
      
      !> Computing the small group. 
      !> The output will be loaded into Gk_out so we have to write it out:
      WRITE(*,*) " "
      WRITE(*,*) " +++++++++++++++++++++++++++++++1"
      WRITE(*,*) " "
      !call K_Star_EXT(vec_k_int, SPGk_in, Gk_out) !no sirve
      call K_Star(vec_k, SPGk_in, Gk_out,ext=.true.)  
      call Write_Group_k(Gk_out)
      call Set_Gk(Gk_out, SPGk_out,ext=.true.)
      call Write_SpaceGroup(SPGk_out, full=.true.)
     
      !donde supuestamente SPGk_out es Gk-k el grupo extendido

      !---> MOD(A,P) computes the remainder of the division of A by P.
      !---> vector v starts in 1: v(1), v(2), v(3)

      WRITE(*,*) " "
      WRITE(*,*) " +++++++++++++++++++++++++++++++2"
      WRITE(*,*) " now with the SPGk_out: STEP 3 find the H(R)"
      
      ! Pointer to the extended little group:  ",Gk%p(1:2*Gk%ngk)

      !do j=1, size(Gk_out%p(1:2*Gk_out%ngk)) 

      ! allocate here al the {R_S} matrices and set all to zero. Also the 3*3 T vectors   
      !se van a expandir las matrixec_R_S, hay que realocatar para expandir. Usamos un nmax gordo
      nmax=100
      !
      if(allocated(matrices_R_S)) deallocate(matrices_R_S)
         allocate(matrices_R_S(4,4,SPGk_out%multip+nmax)) ! 4,4,num operators; put a variabel
       matrices_R_S=0
       !write(*,*) "check-------------------", matrices_R_S 

     if(allocated(vectors_T)) deallocate(vectors_T)
         allocate(vectors_T(3, SPGk_out%multip)) ! 3,3,nomp
       vectors_T=0
       !write(*,*) "check-------------------", vectors_T 



      !set all to zero? no, we just set to zeros the buffer: mat_H_R
      ! bbuffer not used nymore
      

      do j=1, SPGk_out%multip
       
        
        WRITE(*,*) "  " !   limite=",SPGk_out%nk
        vec_k_prime=matmul(vec_k, SPGk_out%SymOp(j)%Rot) ! ese SPGk_out
        WRITE(*,*) " k_prime = ", vec_k_prime

        vec_k_condition1=vec_k_prime-vec_k ! k1=k'-k  <=>  k' = kR = +k + H(R)   ver si H es de red reciproca
        vec_k_condition2=vec_k_prime+vec_k ! k2=k'+k  <=>  k' = kR = -k + H(R)
  
        !para ver si es de la red reciproca usese k_Equiv(H,k) metiendo k=(000) fijo 

        WRITE(*,*) " k1 =      ", vec_k_condition1
        WRITE(*,*) " k2 =      ", vec_k_condition2
       
        vec_H_test=(/ 0, 0, 0 /)
        WRITE(*,*) "vec_H_test  =      ", vec_H_test


        if(   k_EQUIV(vec_k_condition1, vec_H_test, SPGk_out%SPG_lat)   ) then  ! si k1=H(R) es de la red reciproca
             WRITE(*,*) "     found  cond 1 . expanding and storing   " 
             vec_H_R(1:3)=vec_k_condition1
             vec_H_R(4:4)=1 ! for cond 1

             !join H(R(3*3)) and vec_H_R(1*4)  as a the j-esim matrices_R_S

             WRITE(*,*) " R 3*3 :"
               do  n=1,3
                 WRITE(*,"(a,3i4)") "      ",     SPGk_out%SymOp(j)%Rot(n,:) ! la matiz R
                 matrices_R_S(n,:,j)= SPGk_out%SymOp(j)%Rot(n,:)
               end do

            WRITE(*,*) "     H(R)=", vec_H_R 
            matrices_R_S(:,4,j) = vec_H_R  ! add thevector
            
            WRITE(*,*) "  "

            !show R_S j-esi matrix
            WRITE(*,*) " R_S 4*4 no ordenada    ",  matrices_R_S(:,:,j)  
            WRITE(*,*) " R_S 4*4 ordenada:"            
            do n=1,4
               write(*, '(1000F14.7)')( matrices_R_S(m,n, j) ,m=1,4)
            end do 

           
  
       else 
             WRITE(*,*) "     found cond 2  . (si no es la 1 es la 2) expanding and storing   "
             vec_H_R(1:3)=vec_k_condition2
             vec_H_R(4)=-1 ! for cond 2
  
             !join H(R(3*3)) and vec_H_R(1*4)  as a the j-esim matrices_R_S

             WRITE(*,*) " R 3*3 :"
               do  n=1,3
                 WRITE(*,"(a,3i4)") "      ",     SPGk_out%SymOp(j)%Rot(n,:) ! la matiz R
                 matrices_R_S(n,:,j)= SPGk_out%SymOp(j)%Rot(n,:)
               end do

            WRITE(*,*) "     H(R)=", vec_H_R 
            matrices_R_S(:,4,j) = vec_H_R  ! add thevector
             
            WRITE(*,*) "  "

            !show R_S j-esi matrix
            !show R_S j-esi matrix
            WRITE(*,*) " R_S 4*4 no ordenada    ",  matrices_R_S(:,:,j)  
            WRITE(*,*) " R_S 4*4 ordenada:"            
            do n=1,4
               write(*, '(1000F14.7)')( matrices_R_S(m,n, j) ,m=1,4)
            end do

       end if

 
 
        ! hacer lo del vector t fases    SPGk_out%SymOp(j)%Rot(n,:) pero con Trans
           
       
       do i=1, SPGk_out%multip
        !   !matrices_R_S = SPGk_out%Symop(i)%rot
            vectors_T(:,i)= SPGk_out%Symop(i)%tr
            WRITE(*,*) " Vectores de T (individuales):"  , SPGk_out%Symop(i)%tr
          
       end do
 
       WRITE(*,*) " Vectores de T (todos):"  , vectors_T    ! PONERLOS EN ORDEN
      

 
        ! hacer el tipo e grupo supersepacial

      end do ! del j



      WRITE(*,*) ""
      WRITE(*,*) ""

      ! consist grupo paso 4
      WRITE(*,*) " MATRICES fuera del bucle donde se consiguieron (exportada so t osys)" 

           do j=1,SPGk_out%multip
              WRITE(*,*) ""
           do n=1,4
              
               write(*, '(1000F14.7)')( matrices_R_S(m,n, j) ,m=1,4)
           end do
           end do




  !      
   !--- Now starts the expansion of the generator set to get a group
   ! 
   ! sizegg is the instantaneous size of the expanding set (a better name would be sizeEXPSET)
   ! it is initialized to 3 which is the initial number of generators
   !
   ! definir sizegg, iter, P, nmax

   if(allocated(P)) deallocate(P) ! matriz prueba
   allocate(P(4,4)) ! alocatar, ojo que sera variable, de momento 4
   P=0 ! initialization  to zero

   !se van a expandir las matrixec_R_S, hay que realocatar para expandir. Usamos un nmax gordo
   !nmax=100, hevHO ARRIBA DESDE EL PPRINCIPIO

   !define nmax:
   !write(unit=*,fmt="(a)") " => Enter the maximum dimension of the group (as sec. limit): "
   !read(unit=*,fmt=*) nmax
   !if (nmax==0) exit
   !
   !nmax=50


   !PLAN
   !antes de nada hay que sacar los indices de los operadores, bastan solo los R 3*3 de Gk
   ! buscar un ejemplo de funcion recursiva
   !P(:,:)= matmul(matrices_R_S(:,:,i),matrices_R_S(:,:,i))
   !   if(P isqual Identity)
   !   done
   !   else -no igual a Id   
   !    incrementar indice iniciado en 1
   !    volver a multiplicar P    


   write(*, *) ""
   write(*, *) " ---------------matmul"
   write(*, *) ""
   !basta el grupo de Gk, los 4 operadores!! en el ejemplo concreto Pnma 0.321 0 0.5, claro
   sizegg=SPGk_out%numops ! Empezamos solo con los Gk. Con G-k tenemos 8 operadores de sim en este ejemplo de Pnma 0.321 0 0.5
   write(*, *) "sizegg ",sizegg

   !allocate and initialize vector buffers ti=(t1, t2, t3, t4) and tj=(t1', t2', t3', t4'):
   ! (t1, t2, t3) are the T part, stored in vectors_T
   !
   

   iter=0
   fmtt=" "
   write(unit=fmtt,fmt="(a,i5,a)") "(a,i3,a,i3,a,",n*n,"i3,a,i4)"
   do_ext:do   
     iter=sizegg
     do_i:do i=1,sizegg 
       do_j: do j=1,sizegg ! multiply the matrices of the R_S set each other

          ! arrays de posibles valores para t4 de la primera matriz (i) y de la segunda (j):
          
          t4_possible_values_1 = m_values/n_values(i) !arrray m_values div por indice de la matriz i
          t4_possible_values_2 = m_values/n_values(j) !arrray m_values div por indice de la matriz i
 
WRITE(*,*) " "
WRITE(*,*) " "
WRITE(*,*) " m_values-------------", m_values
WRITE(*,*) " n_values-------------", n_values
WRITE(*,*) " t4_possible_values_1 = m_values/n_values(of matrix=",i,")-------------", t4_possible_values_1
!WRITE(*,*) " t4_possible_values_1-------------", t4_possible_values_1


          
          WRITE(*,*) "loop: multuplying ",i,"x",j
          WRITE(*,*) "    t4_possible_values of ",i," are: ",t4_possible_values_1
          WRITE(*,*) "    t4_possible_values of ",j," are: ",t4_possible_values_2

           !variar todos los t4
           do_t1:do t1=1, size(t4_possible_values_1) 
           do_t2:do t2=1, size(t4_possible_values_2) 

            !t_1
            vec_k_buffer_1(1:3)=vectors_T(:,i)
            vec_k_buffer_1(4)=t4_possible_values_1(t1) 
            WRITE(*,*) "  -------------------------- (t1,t2)=",t1," ",t2 
            WRITE(*,*) "  --------------------------vec_k_buffer_1",vec_k_buffer_1
            !t_2
            vec_k_buffer_2(1:3)=vectors_T(:,j)
            vec_k_buffer_2(4)=t4_possible_values_2(t2) 

            P(:,:)      = matmul(matrices_R_S(:,:,i), matrices_R_S(:,:,j)) 
            vec_SSG_law = matmul(matrices_R_S(:,:,i), vec_k_buffer_2) + vec_k_buffer_1
            WRITE(*,*) "        vec_k_buffer_1", vec_k_buffer_1
            WRITE(*,*) "        vec_k_buffer_2", vec_k_buffer_2


            eq=.false. ! set default guess: they are not equal
            do_k:do k=1,sizegg ! comprobar si esta en el grupo generador
                !WRITE(*,*) "test:", i, j, k

                if(Equal_Matrix(P(:,:), matrices_R_S(:,:,k), 4)) then ! 4=dim of matrix
                   eq=.true. ! in in the group, display values of t1, t2
                   write(*,*) "  P is IN the group. Values (t1,t2)=", vec_k_buffer_1(4), " ", vec_k_buffer_2(4)
                   exit !exit current loop of k-index
                end if  
             end do do_k
             if(eq) cycle do_j
             !else
             write(*,*) "  P is NOT the group"
        

           end do do_t2
           end do do_t1


           !t_1
           !vec_k_buffer_1(1:3)=vectors_T(:,i) 
           !vec_k_buffer_1(4)=t4_values(t) ! t4 array starts in 1  
           !vec_k_buffer_1_sin = MOdulo_Lat(vec_k_buffer_1)
           !write(*, *) "   t_i    =",vec_k_buffer_1
           !write(*, *) "   t_i_sin=",vec_k_buffer_1_sin

           !t_2
           !vec_k_buffer_2(1:3)=vectors_T(:,j)
           !vec_k_buffer_2(4)=m_values(t)/n_values(t) 
           !vec_k_buffer_2_sin = Modulo_Lat(vec_k_buffer_2)
           !write(*, *) "   t_j    =",vec_k_buffer_2
           !write(*, *) "   t_j_sin=",vec_k_buffer_2_sin

           !end do do_t

           !P(:,:)      = matmul(matrices_R_S(:,:,i), matrices_R_S(:,:,j))! assign
           !vec_SSG_law = matmul(matrices_R_S(:,:,i), vec_k_buffer_2) + vec_k_buffer_1

           !write(*, *) "   "
           !write(*, *) "   vector SSG law R_S_",i," x t_",j," + t_", i," = ", vec_SSG_law
           !write(*, *) "   matrix SSG law ",i,"x",j 
          ! do n=1,4 
           !write(*, '(1000F14.7)')( P(m,n), m=1,4)
          ! end do
           !write(*, *) "   "

           !end do do_t



           !define vectors with the t4 as a incongnit
           ! t_i=(vector de traslacion de i, t4_i) =  vec_k_buffer_1
           ! t_j=(vector de traslacion de j, t4_j) = vec_k_buffer_2
           ! v=R_j*t_j + t_i=v(mi/ni, mj/nj)
           ! v_sin=Modulo_Lat(v) 
           !
           ! if P esta en el grupo de operadores solo de Gk, ciclar
           ! if P no esta, cer que valores de t4_i, t4_j serian consistentes con que lo fuera  
           !
           ! recuerda que solo buscamos los {t4(R)}, seran en este ejemplo 4 t4's, 
           ! uno para cada operador de R simple rot de Gk          
 
           !show P matrix
           !WRITE(*,*) " "
           !WRITE(*,*) " multiplying: P=:", i,"x",j
           !do n=1,4             
           !    write(*, '(1000F14.7)')( P(m,n) ,m=1,4)
           !end do

           !eq=.false. ! set default guess: they are not equal
           !do_k:do k=1,sizegg ! comprobar si esta en el grupo generador
           !     !WRITE(*,*) "test:", i, j, k
!
           !     if(Equal_Matrix(P(:,:), matrices_R_S(:,:,k), 4)) then ! 4=dim of matrix
           !        eq=.true. ! if repeated it does not expand the set, keep multiplying
           !        write(*,*) "  P is IN the group, with indeex=", k 
           !        exit !exit current loop of k-index
           !     end if  
           !end do do_k  
           !if(eq) cycle do_j
           !
           !here is the else path (if NOT in the group, do the latticemod):
           !write(*,*) "  P is NOT the group: primero quitarle el modulolat, luego check if they are of the form R1t2+t1"
           !vector1=Modulo_Lat(vector2_operador_trans_correspondiente_a_esta_matriz)
           !vec_modulo_lat=0 ! no init here! just for test
           !vec_modulo_lat=Modulo_Lat(vectors_T(:,i)) ! i or j?? cambiar a un ejemplo en q salga
           
           !sizegg=sizegg+1   ! do we have to add it to group? NO
           !if(sizegg > nmax) exit do_ext 
          !if(sizegg > nmax) then 
          !  exit do_ext
          !  write(unit=*,fmt="(a)") "  Limit reached!!!"                   
          !end if   
          ! matrices_R_S(:,:,sizegg)=P  ! include it to the set 
          ! write(unit=*,fmt=fmtt) " =>  Including the matrix: R_S(",i,") x R_S(",j,")=", P(:,:)
          ! write(unit=*,fmt=fmtt) " =>  Current size of Set=",sizegg 
          ! cycle do_ext       

          ! see Get_T_SubGroups, in CFML_Symmetry.f90, line 6908, 6884

          
       !  
       end do do_j
     end do do_i
     if (iter == sizegg) exit
   end do do_ext !main loop Procedure




      WRITE(*,*) " "
      WRITE(*,*) " ================== END. ======================"
      WRITE(*,*) " "

      WRITE(*,*) " "
      WRITE(*,*) " "
      WRITE(*,*) " "


   end do




End Program SSPG_Info
