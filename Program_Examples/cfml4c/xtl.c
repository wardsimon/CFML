// Tutorial for CrysFML    
// ./make_xtl; ./xtl test

#include <stdio.h>
#include <string.h>
#include "xtl.h"

int main(int argc,char *argv[])
{  
  int   length;
  char  line[80];
  float  cell[3],ang[3],number; 
  C_Crystal_Cell_Type   ccell;
  C_Space_Group_Type    cspgr;
  C_Atom_List_Type      catom;
  
  if (argc>1) 
  {
      strcpy (line,argv[1]);
  }
  else 
  {
    printf("Give a filename, please! : ");
    scanf("%s",line); 
  };
  length=strlen(line);
  //xtl2cfml_(&length,line,cell,ang);
  //printf("%7.3lf %7.3lf %7.3lf %7.2lf %7.2lf %7.2lf\n", 
  //    cell[0],cell[1],cell[2],ang[0],ang[1],ang[2]);    
  ccell.cell[0]=1.0;
  number=1.0;
  printf("C program: %s %d \n ",line,length);   
  C_Readn_Set_Xtal_Structure(line,length,&ccell,&cspgr,&catom,&number);
  printf("C-Type in C  :  %7.3lf %7.3lf %7.3lf %7.2lf %7.2lf %7.2lf %7.2lf\n", 
        ccell.cell[0],ccell.cell[1],ccell.cell[2],ccell.ang[0],ccell.ang[1],ccell.ang[2],number);
  return 0;
} 