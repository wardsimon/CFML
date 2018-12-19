#include <iostream>
#include "NexusToCFML.h"
#include <string.h>

using namespace std;


extern "C"
{

void print_C(char *string) /* equivalent: char string[]  */
       {
         cout << " print_C"  << string << endl;;
       }

void c_read_init_nxs(char *string ){
//	 cout << "cpp c_read_init_nxs  01 Nexus filename " <<  string << endl;
	 ArgumentHandler * m_AH = ArgumentHandler::getArgumentHandlerInstance();
	 m_AH->setNexusFileName(string);
//	 cout << "cpp c_read_init_nxs  02 "  << endl;
}

void c_read_header_nxs(){
	 ArgumentHandler * AH = ArgumentHandler::getArgumentHandlerInstance();
	 HeaderBloc* lHeaderBloc = HeaderBloc::getHeaderBlocInstance();
	 //cout << "cpp c_read_header_nxs  01 "  << endl;
	 lHeaderBloc->read();
	 //cout << "cpp c_read_header_nxs  02 "  << endl;
}

void c_get_header_numor_nxs(int *p1){
  HeaderBloc* lHeaderBloc = HeaderBloc::getHeaderBlocInstance();
  *p1 = lHeaderBloc->getHeaderNumor();
  //cout << "c_get_header_numor_nxs Numor  = " << p1 << endl;
}

void c_get_header_scantype_nxs(int *p1){
  HeaderBloc* lHeaderBloc = HeaderBloc::getHeaderBlocInstance();
  *p1 = lHeaderBloc->getHeaderScanType();
  //cout << "c_get_header_numor_nxs Numor  = " << p1 << endl;
}


void c_get_header_instr_name_nxs(char *lstring, int *size){
//	 cout << "c_get_header_instr_name_nxs " <<  lstring  << "size  " << *size  << endl;
	 HeaderBloc* lHeaderBloc = HeaderBloc::getHeaderBlocInstance();
	 memcpy(lstring, lHeaderBloc->getHeaderInstrName().c_str(), *size);
	 //cout << "c_get_header_instr_name_nxs  "<< lstring  << endl;
}

void c_get_header_user_nxs(char *lstring, int *size){
	 //cout << "c_get_header_instr_name_nxs " <<  lstring  << "size  " << *size  << endl;
	 HeaderBloc* lHeaderBloc = HeaderBloc::getHeaderBlocInstance();
	 memcpy(lstring, lHeaderBloc->getHeaderUser().c_str(), *size);
	 //cout << "c_get_header_instr_name_nxs  "<< lstring  << endl;
}
void c_get_header_local_contact_name_nxs(char *lstring, int *size){
	 //cout << "c_get_header_instr_name_nxs " <<  lstring  << "size  " << *size  << endl;
	 HeaderBloc* lHeaderBloc = HeaderBloc::getHeaderBlocInstance();
	 memcpy(lstring, lHeaderBloc->getHeaderLocalContact().c_str(), *size);
	 //cout << "c_get_header_instr_name_nxs  "<< lstring  << endl;
}
void c_get_header_date_nxs(char *lstring, int *size){
	 //cout << "c_get_header_instr_name_nxs " <<  lstring  << "size  " << *size  << endl;
	 HeaderBloc* lHeaderBloc = HeaderBloc::getHeaderBlocInstance();
	 memcpy(lstring, lHeaderBloc->getHeaderDate().c_str(), *size);
	 //cout << "c_get_header_instr_name_nxs  "<< lstring  << endl;
}

void c_get_header_subt_nxs(char *lstring, int *size){
	 //cout << "c_get_header_subt_nxs " <<  lstring << endl;
	 HeaderBloc* lHeaderBloc = HeaderBloc::getHeaderBlocInstance();
	 int size2 = sizeof(lHeaderBloc->getHeaderSubT().c_str());
 	 memset(lstring, ' ', *size );
	 memcpy(lstring, lHeaderBloc->getHeaderSubT().c_str(), size2);
}


void c_read_integer_bloc_nxs(){
	 ArgumentHandler * AH = ArgumentHandler::getArgumentHandlerInstance();
	 IntegerBloc* lIntegerBloc = IntegerBloc::getIntegerBlocInstance();
	 lIntegerBloc->read();
}


void c_get_integer_bloc_nxs(int p3[], int size){

//// For testing
//  for (int ki = 0 ; ki<size ; ki++){
//    cout << "cpp par fc_inttab = " << p3[ki] << endl;
//    p3[ki] = ki+100;
//  }
	IntegerBloc* lIntegerBloc = IntegerBloc::getIntegerBlocInstance();
	lIntegerBloc->getIntegerBloc(p3,size);
}



void c_read_float_bloc_nxs(){
	 ArgumentHandler * AH = ArgumentHandler::getArgumentHandlerInstance();
	 FloatBloc* lFloatBloc = FloatBloc::getFloatBlocInstance();
	 lFloatBloc->read();
}


void c_get_float_bloc_nxs(float p3[], int size){

//// For testing
//  for (int ki = 0 ; ki<size ; ki++){
//    cout << "cpp par fc_inttab = " << p3[ki] << endl;
//    p3[ki] = ki*100;
//  }
	FloatBloc* lFloatBloc = FloatBloc::getFloatBlocInstance();
	lFloatBloc->getFloatBloc(p3,size);
}



void c_read_data_bloc_param_nxs(){
	 ArgumentHandler * AH = ArgumentHandler::getArgumentHandlerInstance();
	 DataBloc* lDataBloc = DataBloc::getDataBlocInstance();
	 lDataBloc->read();
}


void c_get_data_bloc_time_nxs(float p3[], int size){

//// For testing
//  for (int ki = 0 ; ki<size ; ki++){
//    cout << "cpp par fc_inttab = " << p3[ki] << endl;
//    p3[ki] = ki*100;
//  }
	DataBloc* lDataBloc = DataBloc::getDataBlocInstance();
	lDataBloc->getDataBlocTime(p3,size);
}


void c_get_data_bloc_moni_nxs(float p3[], int size){

//// For testing
//  for (int ki = 0 ; ki<size ; ki++){
//    cout << "cpp par fc_inttab = " << p3[ki] << endl;
//    p3[ki] = ki*100;
//  }
	DataBloc* lDataBloc = DataBloc::getDataBlocInstance();
	lDataBloc->getDataBlocMoni(p3,size);
}


void c_get_data_bloc_total_count_nxs(float p3[], int size){

//// For testing
//  for (int ki = 0 ; ki<size ; ki++){
//    cout << "cpp par fc_inttab = " << p3[ki] << endl;
//    p3[ki] = ki*100;
//  }
	DataBloc* lDataBloc = DataBloc::getDataBlocInstance();
	lDataBloc->getDataBlocTotalCount(p3,size);
}


void c_get_data_bloc_angle1_nxs(float p3[], int size){

//// For testing
//  for (int ki = 0 ; ki<size ; ki++){
//    cout << "cpp par fc_inttab = " << p3[ki] << endl;
//    p3[ki] = ki*100;
//  }
	DataBloc* lDataBloc = DataBloc::getDataBlocInstance();
	lDataBloc->getDataBlocAngle1(p3,size);
}

void c_get_data_bloc_data_full_nxs(int p3[], int size){

	DataBloc* lDataBloc = DataBloc::getDataBlocInstance();
	lDataBloc->getDataBlocDataFull(p3,size);
}

}


