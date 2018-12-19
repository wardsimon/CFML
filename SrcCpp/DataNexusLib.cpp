/*
 * Nomad Instrument Control Software
 *
 * Copyright 2011 Institut Laue-Langevin
 *
 * Licensed under the EUPL, Version 1.1 only (the "License");
 * You may not use this work except in compliance with the Licence.
 * You may obtain a copy of the Licence at:
 *
 * http://joinup.ec.europa.eu/software/page/eupl
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the Licence is distributed on an "AS IS" basis,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the Licence for the specific language governing permissions and
 * limitations under the Licence.
 */


#include "DataNexusLib.h"

using namespace std;

int DataNexusLib::OpenPath(NXhandle myfile_id,string path1) {
	if (NXopenpath(myfile_id, path1.c_str()) != NX_OK) {
		cout << "Error open " <<path1  << endl;
		return -1;
	}
	return 1;
}

int DataNexusLib::OpenPath(NXhandle myfile_id,string path1, string path2) {
	if (NXopenpath(myfile_id, path1.c_str()) != NX_OK) {
			cout << "Error open " <<path1  << endl;
			return -1;
		}
	if (NXopenpath(myfile_id, path2.c_str()) != NX_OK) {
			cout << "Error open " <<path2  << endl;
			return -1;
		}
	return 1;
}

int DataNexusLib::OpenPath(NXhandle myfile_id,string path1, string path2, string path3) {

	if (NXopenpath(myfile_id, path1.c_str()) != NX_OK) {
		cout << "Error open " << path1 << endl;
		return -1;
	}
	if (NXopenpath(myfile_id, path2.c_str()) != NX_OK) {
		cout << "Error open " << path2 << endl;
		return -1;
	}
	if (NXopenpath(myfile_id, path3.c_str()) != NX_OK) {
		cout << "Error open " << path3 << endl;
		return -1;
	}

	return 1;
}

string DataNexusLib::getStringValue(NXhandle myfile_id,string champ){
	int NXtype, NXrank;
	int NXdims[32];
	void *data_buffer;
    string err = "777";
    string res;
	if (NXopendata (myfile_id, champ.c_str())!=NX_OK){ cout << "FRA03"<< endl; return err;}
	if (NXgetinfo (myfile_id, &NXrank, NXdims, &NXtype)!=NX_OK){ cout << "FRA04"<< endl; return err;}
	if (NXmalloc ((void **) &data_buffer, NXrank, NXdims, NXtype)!=NX_OK){ cout << "FRA05"<< endl; return err;}
	if (NXgetdata (myfile_id, data_buffer)!=NX_OK){ cout << "FRA06"<< endl; return err;}
	res=(string)((char *) data_buffer);

	// Debug cout << champ <<  "  " << res << endl;
	if (NXfree ((void **) &data_buffer) !=NX_OK) { cout << "FRA06"<< endl; return err;}
    return  res;

}

int DataNexusLib::getIntValue(NXhandle myfile_id,string champ){
	int NXtype, NXrank;
	int NXdims[32];
	void *data_buffer;
    int err = 777;
    int resI;
	if (NXopendata (myfile_id, champ.c_str())!=NX_OK){ cout << "FRA13"<< endl; return err;}
	if (NXgetinfo (myfile_id, &NXrank, NXdims, &NXtype)!=NX_OK){ cout << "FRA14"<< endl; return err;}
	if (NXmalloc ((void **) &data_buffer, NXrank, NXdims, NXtype)!=NX_OK){ cout << "FRA15"<< endl; return err;}
	if (NXgetdata (myfile_id, data_buffer)!=NX_OK){ cout << "FRA16"<< endl; return err;}
 	resI=*((int *)data_buffer);
 	// Debug cout << champ <<  " (Int)  =  " << resI << endl;
	if (NXfree ((void **) &data_buffer) !=NX_OK) { cout << "FRA17"<< endl; return err;}
    return  resI;

}


float DataNexusLib::getFloatValue(NXhandle myfile_id,string champ){
	int NXtype, NXrank;
	int NXdims[32];
	void *data_buffer;
    float err = 777;
    float resI;
	if (NXopendata (myfile_id, champ.c_str())!=NX_OK){ cout << "FRA23"<< endl; return err;}
	if (NXgetinfo (myfile_id, &NXrank, NXdims, &NXtype)!=NX_OK){ cout << "FRA24"<< endl; return err;}
	if (NXmalloc ((void **) &data_buffer, NXrank, NXdims, NXtype)!=NX_OK){ cout << "FRA25"<< endl; return err;}
	if (NXgetdata (myfile_id, data_buffer)!=NX_OK){ cout << "FRA26"<< endl; return err;}
	//resI=(int)((char *) data_buffer);
	resI=*((float *)data_buffer);
	// Debug cout << champ <<  " (float)  =  " << resI << endl;
	if (NXfree ((void **) &data_buffer) !=NX_OK) { cout << "FRA27"<< endl; return err;}
    return  resI;

}


 float DataNexusLib::getFloatIndexValue(NXhandle myfile_id,string champ,int index){

	int NXtype, NXrank;
	int NXdims[32];
	void *data_buffer;
    float err = 777;
    float resF;
	if (NXopendata (myfile_id, champ.c_str())!=NX_OK){ cout << "FRA33"<< endl; return err;}
	if (NXgetinfo (myfile_id, &NXrank, NXdims, &NXtype)!=NX_OK){ cout << "FRA34"<< endl; return err;}
	if (NXmalloc ((void **) &data_buffer, NXrank, NXdims, NXtype)!=NX_OK){ cout << "FRA35"<< endl; return err;}

	if (index > (NXdims[0]-1) ){
		cerr << "ERROR getFloatIndexValue bad index "<< index << endl;
		return err;
	}

	if (NXgetdata (myfile_id, data_buffer)!=NX_OK){ cout << "FRA37"<< endl; return err;}


	float * dataPtr = ((float *)data_buffer);
	resF =  dataPtr[index] ;
	//if (resF<0.0000000001) resF = 0.0; ?? Why
	// Debug cout << champ <<  " (float)  =  " <<  resF << endl;

	if (NXfree ((void **) &data_buffer) !=NX_OK) { cout << "FRA38"<< endl; return err;}
	return  resF;

}


