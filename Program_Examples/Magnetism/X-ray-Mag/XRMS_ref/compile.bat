cls
gfortran -c XRMS_types.f90 -I%crysfml%/gfortran/libc 
gfortran -c global_data.f90 -I%crysfml%/gfortran/libc 
gfortran -c XRMS_procedures.f90 -cpp -I%crysfml%/gfortran/libc 
gfortran -c XRMS_ref.f90 -cpp -I%crysfml%/gfortran/libc 
gfortran *.o -o XRMS_ref -L%crysfml%/gfortran/libc -lcrysfml
del *.o *.mod > nul
