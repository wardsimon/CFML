// Tutorial for CrysFML    
// ./make_xtl; ./xtl test

#include <stdio.h>
#include <string.h>
#include "xtl.h"

int main(int argc,char *argv[])
{  
  int   length;
  char  line[80];
  double cell[3],ang[3];
  
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
  xtl2cfml_(&length,line,cell,ang);
  printf("%7.3lf %7.3lf %7.3lf %7.2lf %7.2lf %7.2lf\n", 
      cell[0],cell[1],cell[2],ang[0],ang[1],ang[2]); 
  return 0;
} 