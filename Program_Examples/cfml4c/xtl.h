#ifndef __CBIND_H__
#define __CBIND_H__

typedef struct C_Crystal_Cell_Type_{
      float cell[3], ang[3];
      float cell_std[3], ang_std[3];
      float rcell[3], rang[3];
      float GD[9],GR[9];
      float Cr_Orth_cel[9];
      float Orth_Cr_cel[9];
      float BL_M[9];
      float BL_Minv[9];
      float CellVol;
      float RCellVol;
      char  CartType[1];
}     C_Crystal_Cell_Type;

typedef struct C_Atom_Type_{
   char    Lab[10] ;
   char    ChemSymb[2];
   char    SfacSymb[4];
   _Bool   Active;
   int     Z;
   int     Mult;
   float   X[3];
   float   X_Std[3];
   float   MX[3];
   int     LX[3];
   float   Occ;
   float   Occ_Std;
   float   MOcc;
   int     LOcc;
   float   Biso;
   float   Biso_std;
   float   MBiso;
   int     LBiso;
   char    Utype[4];
   char    ThType[5];
   float   U[6];
   float   U_std[6];
   float   Ueq;
   float   MU[6];
   int     LU[6];
   float   Charge;
   float   Moment;
   int     Ind[5] ;
   int     NVar;
   float   VarF[10];
   char    AtmInfo[40];
} C_Atom_Type;

typedef struct C_Atom_List_Type_{
    int         natoms;
    C_Atom_Type *atom;
}  C_Atom_List_Type;

typedef struct C_Wyck_Pos_Type_{ 
       int              multp    ;
       char             site[6]     ;
       int              norb     ;
       char             str_orig[40] ;
       char             *str_orbit[40];
}  C_Wyck_Pos_Type;

typedef struct C_Wyckoff_Type_{ 
       int                           num_orbit;
       C_Wyck_Pos_Type               orbit[26]    ;
}  C_Wyckoff_Type;

typedef struct C_Sym_Oper_Type_{ 
    int     Rot[9] ;
    float   Tr[3]  ;
}  C_Sym_Oper_Type;

typedef struct C_Space_Group_Type_{
      int                     NumSpg;
      char                    SPG_Symb[20];
      char                    Hall[16];
      char                    CrystalSys[12];
      char                    Laue[5];
      char                    PG[5];
      char                    Info[5];
      char                    SG_setting[80];
      _Bool                   Hexa;
      char                    SPG_lat[1];
      char                    SPG_latsy[2];
      int                     NumLat;
      float                   Latt_trans[36];
      char                    Bravais[51];
      char                    Centre[80];
      int                     Centred;
      float                   Centre_coord[3];
      int                     NumOps;
      int                     Multip;
      int                     Num_gen;
      C_Sym_Oper_Type         SymOp[192];
      char                    *SymopSymb[192]; 
      C_Wyckoff_Type          Wyckoff;
      float                   R_Asym_Unit[6];
} C_Space_Group_Type;


extern void C_Readn_Set_Xtal_Structure(char*,int,C_Crystal_Cell_Type*,C_Space_Group_Type*,C_Atom_List_Type*,float*);
//extern void xtl2cfml_(int*, char*, double*, double*);

#endif