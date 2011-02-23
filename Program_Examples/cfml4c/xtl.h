#ifndef __CBIND_H__
#define __CBIND_H__

typedef struct C_Crystal_Cell_Type_{
      int     var1;
      float   var2;
      int     var3[6];
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
} C_Crystal_Cell_Type;

extern void xtl2cfml_(int*, char*, double*, double*);

#endif