Program Testing02
    !---- Use modules ----!
    Use CFML_GlobalDeps
    Use CFML_Maths
    
    implicit none
    
    !---- Local Variables ----!
    integer, parameter              :: NITER=10000000
    character(len=132)              :: linea
    real(kind=cp)                   :: vini,vfin
    real(kind=cp),   dimension(3,3) :: b1,b2,b3 
    integer                         :: n1,n2,i,j
    integer,         dimension(3,3) :: a1
    complex(kind=cp),dimension(4,4) :: c1,c2,c3
    complex(kind=cp)                :: det 
    
    !> Header
    print*,"-------------------------------"
    print*,"----    CrysFML Testing    ----"
    print*,"---- JGP-JRC        Test02 ----"
    print*,"-------------------------------"
    
    !> GCD
    print*," GCD Calculations..."
    call cpu_time(vini)
    do i=1,NITER
       n1=gcd(180, 324)       ! Result=36
       if (i < 5) print*,n1
    end do   
    call cpu_time(vfin)
    print*,vini,vfin
    print '("OK!  Time = ",f6.3," seconds.")',vfin-vini
    pause
    
    !> LCM    
    print*," LCM Calculations..." 
    call cpu_time(vini) 
    do i=1,NITER
       n1=lcm(180, 324)       ! Result 1620
       if (i < 5) print*,n1
    end do   
    call cpu_time(vfin)
    print*,vini,vfin
    print '("OK!  Time = ",f6.3," seconds.")',vfin-vini
    pause
    
    !> Matrix operations
    a1(1,:)=[1,2,3]
    a1(2,:)=[0,9,2]
    a1(3,:)=[4,3,1]
    
    b1(1,:)=[1.0, -1.0, 0.0]
    b1(2,:)=[0.0 , 1.0, 0.0]
    b1(3,:)=[2.0 , 0.0, 1.0]
    
    b2(1,:)=[ 1.0,   1.0, 0.0]
    b2(2,:)=[ 0.0 ,  1.0, 0.0]
    b2(3,:)=[-2.0 , -2.0, 1.0]
    
    c1(1,1)=(47.0, -15.0)
    c1(1,2)=(62.0,   5.0)
    c1(1,3)=( 0.0, -72.0)
    c1(1,4)=(61.0,  20.0)
    
    c1(2,1)=(   6.0,  14.0)
    c1(2,2)=( -17.0,   3.0)
    c1(2,3)=(-102.0,  91.0)
    c1(2,4)=(   7.0, -12.0)
    
    c1(3,1)=(  13.0, -55.0)
    c1(3,2)=(  32.0,   8.0)
    c1(3,3)=(  41.0,   7.0)
    c1(3,4)=(  25.0,   1.0)

    c1(4,1)=( 111.0,  25.0)
    c1(4,2)=(  40.0,   0.0)
    c1(4,3)=(  12.0, -82.0)
    c1(4,4)=(  58.0, -30.0)
    
    c2(1,1)=(   0.0039,  -0.0080)
    c2(1,2)=(  -0.0084,  -0.0041)
    c2(1,3)=(  -0.0095,   0.0179)
    c2(1,4)=(   0.0002,  -0.0017)
    
    c2(2,1)=(   0.0368,   0.0218)
    c2(2,2)=(   0.0148,  -0.0271)
    c2(2,3)=(  -0.0338,  -0.0351)
    c2(2,4)=(  -0.0050,  -0.0162)
    
    c2(3,1)=(  -0.0042,  -0.0051)
    c2(3,2)=(  -0.0077,  -0.0018)
    c2(3,3)=(   0.0035,   0.0048)
    c2(3,4)=(   0.0007,   0.0036)

    c2(4,1)=(  -0.0197,  -0.0166)
    c2(4,2)=(  -0.0015,   0.0188)
    c2(4,3)=(   0.0337,   0.0154)
    c2(4,4)=(   0.0054,   0.0172)
    
    print*," Integer Determinat..." 
    call cpu_time(vini) 
    do i=1,NITER
       n1=determ(a1,3)
       if (i < 5) print*,n1    ! Result -89
    end do   
    call cpu_time(vfin)
    print*,vini,vfin
    print '("OK!  Time = ",f6.3," seconds.")',vfin-vini
    pause
    
    print*," Inverse of real matrix..."
    call cpu_time(vini)  
    do i=1,NITER
       b3=inverse_array(b1)                        !Result=b2
       if (mod(i,2)==0 .and. i <= 10) print*,b3
    end do   
    call cpu_time(vfin)
    print*, vini,vfin
    print '(" Time = ",f6.3," seconds.")',vfin-vini
    pause
    
    print*," Inverse of complex matrix..."
    call cpu_time(vini)  
    do i=1,NITER
       c3=inverse_array(c1)
       if (mod(i,2)==0 .and. i <= 10) print*,c3
    end do   
    call cpu_time(vfin)
    print*, vini, vfin
    print '(" Time = ",f6.3," seconds.")',vfin-vini
    pause
    
End Program Testing02


