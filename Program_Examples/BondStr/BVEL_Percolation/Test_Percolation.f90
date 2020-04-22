Program Test_Percolation

  Use CFML_GlobalDeps
  Use CFML_Percolation

  Implicit None

  Integer i, Lun

  Real(kind=cp) Emin,Eini,Eend,dE,dE_ini
  Real(kind=cp), Dimension(3) :: E_percol,E_percol_aux
  Real(kind=cp), Dimension(:,:,:), Allocatable :: rho

  Character(len=256) pgrid_file, out_file
  Character(len=1), Dimension(3) :: axis

  Lun = 11
  pgrid_file  = "LiFePO4-D20_bvel.pgrid"

  ! Read BVEL map

  Call Read_BVEL(pgrid_file,Lun,rho)

  ! Set initial values for percolation analysis

  Emin        = Minval(rho)
  rho         = rho - Emin
  Emin        = 0.0
  Eini        = 0.0
  Eend        = 3.0
  dE          = 0.5
  E_percol(:) = -1.
  axis(1)     = "a"
  axis(2)     = "b"
  axis(3)     = "c"

  ! Percolation analysis

  Write(unit=*,fmt="(/,a)") "Computing first estimation of percolation energies (it can take some minutes) ...."

  Call Percol_Analysis(rho,Emin,Eini,Eend,dE,E_percol)

  Do i = 1 , 3
     If (E_percol(i) > 0.0) Then
        Write(unit=*,fmt="(tr4,a,a1,a,f6.2,a3)") "Percolation along ", axis(i), ": Yes, Percolation energy: ", E_percol(i), " eV"
     Else
        Write(unit=*,fmt="(tr4,a,a1,a)") "Percolation along ", axis(i), ": No"
     End If
  End Do

  Write(unit=*,fmt="(/,a)") "Refining energies...."

  dE_ini = dE
  Do i = 1 , 3
     dE = dE_ini
     If (E_percol(i) > 0.) Then
        Write(unit=*,fmt="(tr4,a,a1)") "axis ", axis(i)
        Eini = E_percol(i) - dE
        dE   = 0.1
        Eend = E_percol(i) + 0.01
        Write(unit=*,fmt="(tr6,a,f4.2,a,f4.2,a)") "Searching percolation between ",Eini," and ",Eend," eV"
        Call Percol_Analysis(rho,Emin,Eini,Eend,dE,E_percol_aux,axis=i)
        E_percol(i) = E_percol_aux(i)
        Write(unit=*,fmt="(tr8,a,f6.2,a3)") "Percolation energy: ", E_percol(i), " eV"
        Eini = E_percol(i) - dE
        dE   = 0.01
        Eend = E_percol(i) + 0.01
        Write(unit=*,fmt="(tr6,a,f4.2,a,f4.2,a)") "Searching percolation between ",Eini," and ",Eend," eV"
        Call Percol_Analysis(rho,Emin,Eini,Eend,dE,E_percol_aux,axis=i)
        E_percol(i) = E_percol_aux(i)
        Write(unit=*,fmt="(tr8,a,f6.2,a3,/)") "Percolation energy: ", E_percol(i), " eV"
     End If
  End Do

End Program Test_Percolation
