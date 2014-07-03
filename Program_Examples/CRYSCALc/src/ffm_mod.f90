! elements 3d:    f = A.exp(-ax**2) + B.exp(-bx**2) + C
!
!                 les parametres A, a, B, b et C sont extraits de la publication suivante:
!                   E.J. Lisher and J.B. Forsith
!                   Acta. Cryst. A27 (1971), 545
!
!              Sc  0, 1, 2, 3
!              Ti  0, 1, 2, 3
!              V   0, 1, 2, 3
!              Cr  0, 1, 2, 3
!              Mn  0, 1, 2, 3
!              Fe  0, 1, 2, 3, 4
!              Co  0, 1, 2, 3
!              Ni  0, 1, 2, 3


! scandium
   ffm%atom_name(1)    = 'scandium'
   ffm%atom_label(1,1) = 'Sc'
    ffm%coef(1,1,1:5) =  (/ .5498,  107.3 ,  .4493, 18.27, -.00253/)

   ffm%atom_label(1,2) = 'Sc1'
    ffm%coef(1,2,1:5) =  (/ .4633,   54.18,  .5418, 14.73, -.00633/)

   ffm%atom_label(1,3) = 'Sc2'
    ffm%coef(1,3,1:5) =  (/ .391,    36.16,  .6184, 12.82, -.00941/)

   ffm%atom_label(1,4) = 'Sc3'
    ffm%coef(1,4,1:5) =  (/ .655,    20.84,  .3533, 12.04, -.0097/)



! titane

   ffm%atom_name(2)    = 'titanium'
   ffm%atom_label(2,1) = 'Ti'
    ffm%coef(2,1,1:5) =  (/ .4809, 73.56, .52,   13.67,  -.0045/)

   ffm%atom_label(2,2) = 'Ti1'
    ffm%coef(2,2,1:5) =  (/ .4243, 42.69, .5823, 11.7,   -.00794/)

   ffm%atom_label(2,3) = 'Ti2'
    ffm%coef(2,3,1:5) =  (/ .3667, 30.66, .6443, 10.41,  -.0111/)

   ffm%atom_label(2,4) = 'Ti3'
    ffm%coef(2,4,1:5) =  (/ .3318, 22.84, .6821,  9.425, -.0137/)


! vanadium

   ffm%atom_name(3)    = 'Vanadium'
   ffm%atom_label(3,1) = 'V'
    ffm%coef(3,1,1:5) =  (/ .4438, 54.56, .5589, 10.86,  -.00607/)

   ffm%atom_label(2,2) = 'V1'
    ffm%coef(3,2,1:5) =  (/ .408,  34.87, .5998,  9.615, -.00913/)

   ffm%atom_label(2,3) = 'V2'
    ffm%coef(3,3,1:5) =  (/ .3595, 26.3,  .6524,  8.729, -.0122/)

   ffm%atom_label(2,4) = 'V3'
    ffm%coef(3,4,1:5) =  (/ .3542, 19.9,  .67,    7.86,  -.0151/)

! chrome

   ffm%atom_name(4)    = 'Chromium'
   ffm%atom_label(4,1) = 'Cr'
    ffm%coef(4,1,1:5) =  (/ .4326, 45.13, .571,  8.997, -.00681/)

   ffm%atom_label(4,2) = 'Cr1'
    ffm%coef(4,2,1:5) =  (/ .3957, 29.01, .6129, 8.082, -.0101/)

   ffm%atom_label(4,3) = 'Cr2'
    ffm%coef(4,3,1:5) =  (/ .3774, 21.97, .6352, 7.303, -.0132/)

   ffm%atom_label(4,4) = 'Cr3'
    ffm%coef(4,4,1:5) =  (/ .3574, 17.47, .6587, 6.704, -.016/)


chrome0:

chrome

chrome




