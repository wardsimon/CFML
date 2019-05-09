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

#include "NameServer.h"
#include "./DataBloc.h"
#include "./DataNexusLib.h"
#include "./blosc_filter.h"


#include <fstream>
#include <napi.h>
#include <boost/format.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "hdf5/serial/hdf5.h"
#include "blosc.h"


const std::string NAME_ENTRY0 			 = "/entry0";
const std::string NAME_DATA_SCAN 		 = "/entry0/data_scan";
const std::string NAME_SCAN_VAR 			 = "/entry0/data_scan/scanned_variables";
const std::string NAME_SCAN_VAR_DATA 	 = "/entry0/data_scan/scanned_variables/data";
const std::string NAME_SCAN_VAR_NAME 	 = "/entry0/data_scan/scanned_variables/variables_names";

const std::string NAME_SCAN_VAR_NAME_NAME 	 = "/entry0/data_scan/scanned_variables/variables_names/name";
const std::string NAME_SCAN_VAR_NAME_PROPERTY = "/entry0/data_scan/scanned_variables/variables_names/property";
const std::string NAME_SCAN_VAR_NAME_AXIS 	 = "/entry0/data_scan/scanned_variables/variables_names/axis";
const std::string NAME_SCAN_VAR_NAME_SCANNED  = "/entry0/data_scan/scanned_variables/variables_names/scanned";
const std::string NAME_SCAN_VAR_NAME_UNIT 	 = "/entry0/data_scan/scanned_variables/variables_names/unit";



const std::string NAME_DETECTOR_DATA 	 = "/entry0/data_scan/detector_data";
const std::string NAME_DETECTOR_DATA_DATA = "/entry0/data_scan/detector_data/data";


const std::string NAME_CAMERA_SCAN			= "/entry0/diff_camera";
const std::string NAME_CAMERA_SCAN_VAR		= "/entry0/diff_camera/scanned_variables";
const std::string NAME_CAMERA_SCAN_VAR_DATA	= "/entry0/diff_camera/scanned_variables/data";
const std::string NAME_CAMERA_SCAN_VAR_NAME 	= "/entry0/diff_camera/scanned_variables/variables_names";
const std::string NAME_CAMERA_SCAN_VAR_NAME_NAME = "/entry0/diff_camera/scanned_variables/variables_names/name";
const std::string NAME_CAMERA_SCAN_VAR_NAME_PROPERTY = "/entry0/diff_camera/scanned_variables/variables_names/property";
const std::string NAME_CAMERA_SCAN_VAR_NAME_AXIS = "/entry0/diff_camera/scanned_variables/variables_names/axis";
const std::string NAME_CAMERA_SCAN_VAR_NAME_SCANNED = "/entry0/diff_camera/scanned_variables/variables_names/scanned";
const std::string NAME_CAMERA_SCAN_VAR_NAME_UNIT = "/entry0/diff_camera/scanned_variables/variables_names/unit";

const std::string NAME_CAMERA_IMAGE_DATA	= "/entry0/diff_camera/Image_Data";
#include <boost/date_time/posix_time/posix_time.hpp>

using namespace boost::posix_time;


using namespace std;

#define DIM0            4

DataBloc* DataBloc::_singleton = 0;


DataBloc::DataBloc() {
	m_AH = ArgumentHandler::getArgumentHandlerInstance();

	m_index_time = -1;
	m_index_monitor = -1;
	m_index_totalCount = -1;
	m_index_angles1 = -1;
	m_index_angles2 = -1;
	m_index_angles3 = -1;
	m_index_angles4 = -1;

	m_dataRaw = false;
}


DataBloc::~DataBloc() {

}


DataBloc * DataBloc::getDataBlocInstance() {
	if (_singleton == 0) {
		_singleton = new DataBloc();
	}
	return _singleton;

}


void DataBloc::getDataBlocTime(float p4[], int size) {
	for (int kz = 0; kz < size; kz++) {
		p4 [kz] = m_time[kz];
	}

}
void DataBloc::getDataBlocMoni(float p4[], int size) {
	for (int kz = 0; kz < size; kz++) {
		p4 [kz] = m_monitor[kz];
	}

}
void DataBloc::getDataBlocTotalCount(float p4[], int size) {
	for (int kz = 0; kz < size; kz++) {
		p4 [kz] = m_totaCount[kz];
	}
}
void DataBloc::getDataBlocAngle1(float p4[], int size) {
	int indice_property = m_index_angles2 ; // Omega
	if (indice_property == m_index_angles1) {
		for (int kz = 0; kz < size; kz++) {
			p4[kz] = m_angles1[kz];
		}
	} else if (indice_property == m_index_angles2) {
		for (int kz = 0; kz < size; kz++) {
			p4[kz] = m_angles2[kz];
		}
	} else if (indice_property == m_index_angles3) {
		for (int kz = 0; kz < size; kz++) {
			p4[kz] = m_angles3[kz];
		}
	} else if (indice_property == m_index_angles4) {
		for (int kz = 0; kz < size; kz++) {
			p4[kz] = m_angles4[kz];
		}
	}


}


void DataBloc::getDataBlocDataFull(int p4[], int size) {
	boost::posix_time::ptime lStartedTime;
	lStartedTime = microsec_clock::local_time();
//	for (int kz = 0; kz < size; kz++) {
//		p4 [kz] = (int) m_dread_data[kz];
//	}

	read_data2(p4);
    //std::memcpy(p4, m_dread_data,sizeof(int)*size);


//	time_period timediff(lStartedTime, microsec_clock::local_time());
//    float timePrint = ((float)timediff.length().total_seconds()
//				+ (float)timediff.length().fractional_seconds() / 1000000.);
//    cout << " Time to  getDataBlocDataFull " << timePrint <<endl;


}

void DataBloc::setDataMode( bool dataRaw){
	m_dataRaw = dataRaw;
}


int DataBloc::read() {
	//cout << " DataBloc::read  01 " << endl;
	NXstatus myret;
	NXhandle myfile_id;
	m_nxFileName = m_AH->getNexusFileName();

	myret = NXopen(m_nxFileName.c_str(), NXACC_READ, &myfile_id);
	if (myret != NX_OK)
		return -1;

	string champ;
	string res;
	//cout << " DataBloc::read  02 " << endl;

	DataNexusLib::OpenPath(myfile_id,NAME_DATA_SCAN   );
	champ = "actual_step";
	m_actual_step = DataNexusLib::getIntValue(myfile_id, champ);
	champ = "total_steps";
	m_total_steps = DataNexusLib::getIntValue(myfile_id, champ);

	if ((m_total_steps == 1) && (m_actual_step == 1)) {
		m_is_scan = true;
	} else {
		m_is_scan = false;
	}
	//cout << " DataBloc::read  03 " << endl;

	DataNexusLib::OpenPath(myfile_id,NAME_ENTRY0 );
	champ = "run_number";
	m_RecordNumber = DataNexusLib::getIntValue(myfile_id, champ);

	if (NXclose(&myfile_id) != NX_OK)
		return -1;

	read_scan_var();
	//read_data();
}

// Fonction calculant l'écart-type des nombres contenus dans un vecteur
int DataBloc::save() {
	//cout << " DataBloc::save  " << endl;
	string dataFileName = m_AH->getAsciiFileName();

	ofstream dataFile(dataFileName.c_str(), std::ofstream::out | std::ofstream::app);

	for (int kstep = 1; kstep <= m_actual_step; kstep++) {

		dataFile << string(80, 'S') << "\n";
		int l1 = 0;
		int l2 = 1;
		ostringstream temp1;
		temp1 << std::right << std::setfill(' ') << setw(8) << kstep << setw(8) << m_total_steps - kstep << setw(8)
				<< m_total_steps << setw(8) << m_RecordNumber << setw(8) << l1 << setw(8) << l2;
		dataFile << temp1.str() << string(80 - temp1.str().length(), ' ') << "\n";

		dataFile << string(80, 'F') << "\n";
		int li_co = 4;
		int v = 1;
		dataFile << boost::format("%8i%8i") % li_co % v;
		dataFile << string(64, ' ') << "\n";
		dataFile << "           time         monitor       Total Cou     anglesx1000                 " << "\n";
		//dataFile << string(80, ' ') << "\n";

		ostringstream temp2, temp3;
		temp2 << std::right << std::setfill(' ') << scientific << setprecision(8) << setw(16) << m_time[kstep - 1] << setw(16)
				<< m_monitor[kstep - 1] << setw(16) << m_totaCount[kstep - 1] << setw(16) << m_angles2[kstep - 1] ;
		//<< setw(16) << m_angles2[kstep - 1];
		//cout << "HDS  "<< temp2.str() << endl;
		dataFile << temp2.str() << string(80 - temp2.str().length(), ' ') << "\n";

		//temp3 << std::right << std::setfill(' ') << scientific << setprecision(8) << setw(16) << m_angles3[kstep - 1]
		//		<< setw(16) << m_angles4[kstep - 1];
		//dataFile << temp3.str() << string(80 - temp3.str().length(), ' ') << "\n";

		// Save data
		if (m_dataRaw) {
			//cout << " write rawData" << endl;
		} else {
			//cout << " write CalibrateData" << endl;
		}

		dataFile << string(80, 'I') << "\n";
		l1 = 163840;
		ostringstream temp4;
		temp4 << std::right << std::setfill(' ') << setw(8) << l1;
		dataFile << temp4.str() << string(80 - temp4.str().length(), ' ') << "\n";

		int dataSize = 163840;
		char line[80 + 2];
		int nbvalue = 10;
		int valuesize = 8;

		int q = dataSize / nbvalue;
		int r = dataSize - q * nbvalue;

		int *data = &m_dread_data[dataSize * (kstep - 1)];

		string m_Format;
		m_Format = "%8i";
		char * linePtr = 0;

		for (int i = 0; i < q; ++i) {
			linePtr = line;

			for (int j = 0; j < nbvalue; ++j) {
				sprintf(linePtr, m_Format.c_str(), (int) (*data));
				++data;
				linePtr += valuesize;
			}

			*linePtr = '\n';
			++linePtr;
			*linePtr = '\0';

			dataFile << line;
		}

		linePtr = line;

		if (r != 0) {
			for (int j = 0; j < r; ++j) {
				sprintf(linePtr, m_Format.c_str(), (int) (*data));
				++data;
				linePtr += valuesize;
			}

			for (int j = r * valuesize; j < 80; ++j) {
				*linePtr = ' ';
				++linePtr;
			}

			*linePtr = '\n';
			++linePtr;
			*linePtr = '\0';
			dataFile << line;
		}

	}
	return 0;
}

void DataBloc::decode_var_name(){
	hid_t file_id, dataset_id, dataspace_id, memtype;
	int size1 = m_total_steps ;
	hid_t dataspace;
	herr_t status;
	int status_n, rank;
	hsize_t dims_out[2]; /* dataset dimensions */
	hid_t group_id1, group_id2, group_id3,group_id4;


	file_id =  H5Fopen(m_nxFileName.c_str() , H5F_ACC_RDWR, H5P_DEFAULT);

	group_id1 = H5Gopen1(file_id, "/entry0/");
	if (group_id1 < 0) {
		cout << " Error read  H5Gopen1 cannot open /entry" << endl;
		return;
	}

	group_id2 = H5Gopen1(group_id1, "/entry0/data_scan/");
	if (group_id2 < 0) {
		cout << " Error read  H5Gopen1 cannot open /entry0/data_scan/" << endl;
		return;
	}
	group_id3 = H5Gopen1(group_id2, "/entry0/data_scan/scanned_variables");
		if (group_id3 < 0) {
			cout << " Error read  H5Gopen1 cannot open /entry0/data_scan/scanned_variables/" << endl;
			return;
		}
	group_id4 = H5Gopen1(group_id3, "/entry0/data_scan/scanned_variables/variables_names");
				if (group_id3 < 0) {
					cout << " Error read  H5Gopen1 cannot open /entry0/data_scan/scanned_variables/variables_names" << endl;
					return;
				}
	dataset_id = H5Dopen2(group_id4, "/entry0/data_scan/scanned_variables/variables_names/property", H5P_DEFAULT);
	if (dataset_id < 0) {
		cout << " Error read  H5Dopen2 cannot open /entry0/data_scan/scanned_variables/variables_names/property" << endl;
		return;
	}




    hsize_t     dims[1] = {DIM0};
    char        **rdata;
    int         ndims;
	hid_t filetype;

    filetype = H5Dget_type (dataset_id);
    if (filetype < 0) {
    		cout << " Error read filetype " << endl;
    		return;
    }
    hid_t space = H5Dget_space (dataset_id);
    ndims = H5Sget_simple_extent_dims (space, dims, NULL);
    //cout << "FC22 ndims "  << dims[0]<< endl;
	m_nb_prop = dims[0];

    /*
     * Create the memory datatype.
     */
    memtype = H5Tcopy (H5T_C_S1);
    status = H5Tset_size (memtype, H5T_VARIABLE);

    rdata = (char **) malloc (dims[0] * sizeof (char *));
	/*
	 * Read the data.
	 */
	status = H5Dread(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata);
	for (int i = 0; i < dims[0]; i++) {
		string names = rdata[i];
		//cout << "FC22  i " << i << endl;
		//cout << " FC22 names " << names << endl;
		bool find_name = false;

		string str2 =  names;
		if (str2=="gamma"){
			//cout << "Indice gamma = " <<  i << endl;
			m_prop_index[i] = m_index_angles1;
		}else if (str2=="omega"){
			//cout << "Indice omega = " <<  i << endl;
			m_prop_index[i] =m_index_angles2 ;

		}else if(str2=="chi"){
			////cout << "Indice chi = " <<  i << endl;
			m_prop_index[i] =m_index_angles3 ;

		}else if(str2=="phi"){
			//cout << "Indice phi = " <<  i << endl;
			m_prop_index[i] = m_index_angles4;

		}else if(str2=="Time"){
			////cout << "Indice Time = " <<  i << endl;
			m_prop_index[i] = m_index_time;

		}else if(str2=="TotalCount"){
			//cout << "Indice TotalCount = " <<  i << endl;
			m_prop_index[i] =m_index_totalCount ;

		}else if(str2=="Monitor"){
			//cout << "Indice Monitor = " <<  i << endl;
			m_prop_index[i] = m_index_monitor;

		}else {
			//cout << " FC22 names " << names << "not stored" << endl;
		}
	}


    /*
	 * Close and release resources.  Note that H5Dvlen_reclaim works
	 * for variable-length strings as well as variable-length arrays.
	 * Also note that we must still free the array of pointers stored
	 * in rdata, as H5Tvlen_reclaim only frees the data these point to.
	 */
	status = H5Dvlen_reclaim(memtype, space, H5P_DEFAULT, rdata);
	free(rdata);

	status = H5Dclose(dataset_id);
	if (status < 0) {
		//cout << "error1 dataset_id" << endl;

	}
	status = H5Sclose (space);
	status = H5Tclose (memtype);
    status = H5Tclose (memtype);

	dataset_id = H5Dopen2(group_id4, "/entry0/data_scan/scanned_variables/variables_names/name", H5P_DEFAULT);
	if (dataset_id < 0) {
		//cout << " Error read  H5Dopen2 cannot open /entry0/data_scan/scanned_variables/variables_names/name" << endl;
		return;
	}


    filetype = H5Dget_type (dataset_id);
    if (filetype < 0) {
    		cout << " Error read filetype " << endl;
    		return;
    }
    space = H5Dget_space (dataset_id);
    ndims = H5Sget_simple_extent_dims (space, dims, NULL);
    //cout << "FC23 ndims "  << dims[0]<< endl;
	m_nb_prop = dims[0];

    /*
     * Create the memory datatype.
     */
    memtype = H5Tcopy (H5T_C_S1);
    status = H5Tset_size (memtype, H5T_VARIABLE);

    rdata = (char **) malloc (dims[0] * sizeof (char *));
	/*
	 * Read the data.
	 */
	status = H5Dread(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata);
	for (int i = 0; i < dims[0]; i++) {
		string names = rdata[i];
		//cout << "FC23  i " << i << endl;
		//cout << " FC23 names " << names << endl;
		bool find_name = false;

		string str2 =  names;
		if (str2=="gamma"){
			//cout << "Indice gamma = " <<  i << endl;
			m_prop_index[i] = m_index_angles1;
		}else if (str2=="omega"){
			//cout << "Indice omega = " <<  i << endl;
			m_prop_index[i] =m_index_angles2 ;

		}else if(str2=="chi"){
			////cout << "Indice chi = " <<  i << endl;
			m_prop_index[i] =m_index_angles3 ;

		}else if(str2=="phi"){
			//cout << "Indice phi = " <<  i << endl;
			m_prop_index[i] = m_index_angles4;

		}else if(str2=="Time"){
			////cout << "Indice Time = " <<  i << endl;
			m_prop_index[i] = m_index_time;

		}else if(str2=="TotalCount"){
			//cout << "Indice TotalCount = " <<  i << endl;
			m_prop_index[i] =m_index_totalCount ;

		}else if(str2=="Monitor"){
			//cout << "Indice Monitor = " <<  i << endl;
			m_prop_index[i] = m_index_monitor;

		}else {
			//cout << " FC23 names " << names << "not stored" << endl;
		}


	}


    /*
	 * Close and release resources.  Note that H5Dvlen_reclaim works
	 * for variable-length strings as well as variable-length arrays.
	 * Also note that we must still free the array of pointers stored
	 * in rdata, as H5Tvlen_reclaim only frees the data these point to.
	 */
	status = H5Dvlen_reclaim(memtype, space, H5P_DEFAULT, rdata);
	free(rdata);

    status = H5Dclose(dataset_id);
    	if (status < 0) {
    		cout << "error1 dataset_id" << endl;

    	}
    	status = H5Sclose (space);
    	status = H5Tclose (memtype);
        status = H5Tclose (memtype);

    status = H5Gclose(group_id4);
    	if (status < 0) {
    		cout << "error1 group_id4" << endl;
    	}
	status = H5Gclose(group_id3);
	if (status < 0) {
		cout << "error1 group_id3" << endl;
	}
	status = H5Gclose(group_id2);
	if (status < 0) {
		cout << "error1 group_id2" << endl;
	}
	status = H5Gclose(group_id1);
	if (status < 0) {
		cout << "error1 group_id1" << endl;
	}

	status = H5Fclose(file_id);
	if (status < 0) {
		cout << "error1 file_id" << endl;
	}

}

void DataBloc::read_scan_var(){

	m_time     = new float[m_total_steps];
	m_monitor  = new float[m_total_steps];
	m_totaCount= new float[m_total_steps];
	m_angles1  = new float[m_total_steps];
	m_angles2  = new float[m_total_steps];
	m_angles3  = new float[m_total_steps];
	m_angles4  = new float[m_total_steps];

	m_index_time = 0;
	m_index_monitor =1;
	m_index_totalCount = 2;
	m_index_angles1 = 3;
	m_index_angles2 = 4;
	m_index_angles3 =5;
	m_index_angles4 =6;

    // Default
	for (int ki=0; ki<7;ki++){
		m_prop_index[ki]  = -1;
	}
	m_nb_prop = 1;

	decode_var_name();


	for (int ki = 0; ki < m_total_steps; ki++) {
		m_time[ki] = 0.0;
		m_monitor[ki] = 0.0;
		m_totaCount[ki] = 0.0;
		m_angles1[ki] = 0.0;
		m_angles2[ki] = 0.0;
		m_angles3[ki] = 0.0;
		m_angles4[ki] = 0.0;

	}


	int indice_property = 0;
	int size1 = m_actual_step  ;
	hid_t dataspace;
	herr_t status;
	int status_n, rank;
	hsize_t dims_out[2]; /* dataset dimensions */
	hid_t group_id1, group_id2, group_id3;

	hid_t file_id, dataset_id, dataspace_id;
	hsize_t offset[2]; /* hyperslab offset in the file */

	file_id =  H5Fopen(m_nxFileName.c_str() , H5F_ACC_RDWR, H5P_DEFAULT);


	group_id1 = H5Gopen1(file_id, "/entry0/");
	if (group_id1 < 0) {
		cout << " Error read  H5Gopen1 cannot open /entry" << endl;
		return;
	}

	group_id2 = H5Gopen1(group_id1, "/entry0/data_scan");
	if (group_id2 < 0) {
		cout << " Error read  H5Gopen1 cannot open /entry0/data_scan" << endl;
		return;
	}
	group_id3 = H5Gopen1(group_id2, "/entry0/data_scan/scanned_variables");
		if (group_id3 < 0) {
			cout << " Error read  H5Gopen1 cannot open /entry0/data_scan/scanned_variables" << endl;
			return;
		}

	dataset_id = H5Dopen2(group_id2, "/entry0/data_scan/scanned_variables/data", H5P_DEFAULT);
	if (dataset_id < 0) {
		cout << " Error read  H5Dopen2 cannot open /entry0/data_scan/scanned_variables/data" << endl;
		return;
	}


	dataspace = H5Dget_space(dataset_id); /* dataspace handle */
	if (dataspace < 0) {
			cout << " Error H5Dget_space  " << endl;
			return;
		}
	rank = H5Sget_simple_extent_ndims(dataspace);
	if (rank < 0) {
			cout << " Error H5Sget_simple_extent_ndims" << endl;
			return;
		}
	status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
	//printf("\nRank: %d\nDimensions: %lu x %lu \n", rank, (unsigned long) (dims_out[0]), (unsigned long) (dims_out[1]));



	for (int kpr = 0; kpr < m_nb_prop ;kpr++ ) {
		float* dread_data = new float[size1];
		indice_property = kpr;
		/*
		 * Define hyperslab in the dataset.
		 */
		offset[0] = indice_property;
		offset[1] = 0;
		hsize_t count[2];
		count[0] = 1;
		count[1] = size1;
		status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
		if (status < 0) {
			cout << " Error H5Sselect_hyperslab" << endl;
			return;
		}
		/*
		 * Define the memory dataspace.
		 */
		hsize_t dimsm[2];
		dimsm[0] = 1;
		dimsm[1] = size1;
		hid_t memspace;

		memspace = H5Screate_simple(2, dimsm, NULL);
		if (memspace < 0) {
			cout << " Error H5Screate_simple" << endl;
			return;
		}

		/*
		 * Define memory hyperslab.
		 */
		hsize_t offset_out[2];
		offset_out[0] = 0;
		offset_out[1] = 0;
		hsize_t count_out[2];
		count_out[0] = 1;
		count_out[1] = size1;
		status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
		if (status < 0) {
			cout << " Error H5Sselect_hyperslab" << endl;
			return;
		}

		/*
		 * Read data from hyperslab in the file into the hyperslab in
		 * memory and display.
		 */
		status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, memspace, dataspace, H5P_DEFAULT, dread_data);
		if (status < 0) {
			cout << " Error H5Dread read_scan_var indice properties  " << kpr << endl;
			return;
		}

		//cout << "size1 = " << size1 << endl;
		if (m_prop_index[indice_property] == m_index_angles1) {

			for (int kz = 0; kz < size1; kz++) {
				//cout << "Angle1[" << kz << "]= " << dread_data[kz] << endl;
				m_angles1[kz] = dread_data[kz];
			}
		} else if (m_prop_index[indice_property] == m_index_angles2) {
			for (int kz = 0; kz < size1; kz++) {
				//cout << "Angle2[" << kz << "]= " << dread_data[kz] << endl;
				m_angles2[kz] = dread_data[kz];
			}

		} else if (m_prop_index[indice_property] == m_index_angles3) {

			for (int kz = 0; kz < size1; kz++) {
				//cout << "Angle3[" << kz << "]= " << dread_data[kz] << endl;
				m_angles3[kz] = dread_data[kz];
			}
		} else if (m_prop_index[indice_property] == m_index_angles4) {
			for (int kz = 0; kz < size1; kz++) {
				//cout << "Angle4[" << kz << "]= " << dread_data[kz] << endl;
				m_angles4[kz] = dread_data[kz];
			}
		} else if (m_prop_index[indice_property] == m_index_monitor) {
			for (int kz = 0; kz < size1; kz++) {
				//cout << "Monitor[" << kz << "]= " << dread_data[kz] << endl;
				m_monitor[kz] = dread_data[kz];
			}
		} else if (m_prop_index[indice_property] == m_index_time) {
			for (int kz = 0; kz < size1; kz++) {
				//cout << "Time[" << kz << "]= " << dread_data[kz] << endl;
				m_time[kz] = dread_data[kz];
			}
		} else if (m_prop_index[indice_property] == m_index_totalCount) {
			for (int kz = 0; kz < size1; kz++) {
				//cout << "totalCount[" << kz << "]= " << dread_data[kz] << endl;
				m_totaCount[kz] = dread_data[kz];
			}
		}

		delete[] dread_data;

	}

	status = H5Dclose(dataset_id);
	if (status < 0) {
		cout << "error1 dataset_id" << endl;
	}
	status = H5Gclose(group_id1);
	if (status < 0) {
		cout << "error1 group_id1" << endl;
	}
	status = H5Gclose(group_id2);
	if (status < 0) {
		cout << "error1 group_id2" << endl;
	}

	status = H5Fclose(file_id);
	if (status < 0) {
		cout << "error1 file_id" << endl;
	}

}

void  DataBloc::read_data(){

//	cout << "DataBloc::read_data "<< m_total_steps<< endl;
	int size1 = 640*256   *  m_total_steps;
	m_dread_data = new int[size1];

	herr_t status;

	char *version, *date;
	int r = register_blosc(&version, &date);

	hid_t group_id1, group_id2, group_id3;

	hid_t file_id, dataset_id, dataspace_id;
	hsize_t offset[2]; /* hyperslab offset in the file */

	file_id = H5Fopen(m_nxFileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

	group_id1 = H5Gopen1(file_id, "/entry0/");
	if (group_id1 < 0) {
		cout << " Error read  H5Gopen1 cannot open /entry" << endl;
		return;
	}

	group_id2 = H5Gopen1(group_id1, "/entry0/data_scan");
	if (group_id2 < 0) {
		cout << " Error read  H5Gopen1 cannot open /entry0/data_scan" << endl;
		return;
	}
	group_id3 = H5Gopen1(group_id2, "/entry0/data_scan/detector_data");
	if (group_id3 < 0) {
		cout << " Error read  H5Gopen1 cannot open /entry0/data_scan/detector_data" << endl;
		return;
	}

	hid_t HDF5_datatype;
	HDF5_datatype = H5T_NATIVE_UINT32;

	if (m_dataRaw ){
			cout << " read rawdata" << endl;
			dataset_id = H5Dopen2(group_id3, "/entry0/data_scan/detector_data/data_raw",H5P_DEFAULT);
	}else {
			cout << " read CalibrateData" << endl;
			dataset_id = H5Dopen2(group_id3, "/entry0/data_scan/detector_data/data",H5P_DEFAULT);
	}

	if (dataset_id < 0) {
			cout << "Error H5Dopen2   /entry0/data_scan/detector_data/data_raw " << endl;
	}

	status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			m_dread_data);

	if (status < 0){
			cout << "Error H5Dread data_raw" << endl;
	}

//  Debug
//	int V_x = 640;
//	int V_y = 256;
//	int V_z = m_total_steps;
//
//	int dres1_sum_raw[m_total_steps];
//
//	for (int kz = 0; kz < 3; kz++) {
//		dres1_sum_raw[kz] = 0;
//			for (int kl = 0; kl < (V_y * V_x); kl++) {
//				dres1_sum_raw[kz] = dres1_sum_raw[kz] + m_dread_data[kl + kz * (V_y * V_x)];
//			}
//			cout << " frame_raw[" << kz << "]sum " << dres1_sum_raw[kz] << "---";
//	cout <<  endl << " End Data raw"<< endl;
//
//	}
}




void  DataBloc::read_data2(int *m_dread_data){

	boost::posix_time::ptime lStartedTime;
	lStartedTime = microsec_clock::local_time();

//	cout << "DataBloc::read_data "<< m_total_steps<< endl;
	int size1 = 640*256   *  m_total_steps;
	//m_dread_data = new int[size1];

	herr_t status;

	char *version, *date;
	int r = register_blosc(&version, &date);

	hid_t group_id1, group_id2, group_id3;

	hid_t file_id, dataset_id, dataspace_id;
	hsize_t offset[2]; /* hyperslab offset in the file */

	file_id = H5Fopen(m_nxFileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

	group_id1 = H5Gopen1(file_id, "/entry0/");
	if (group_id1 < 0) {
		cout << " Error read  H5Gopen1 cannot open /entry" << endl;
		return;
	}

	group_id2 = H5Gopen1(group_id1, "/entry0/data_scan");
	if (group_id2 < 0) {
		cout << " Error read  H5Gopen1 cannot open /entry0/data_scan" << endl;
		return;
	}
	group_id3 = H5Gopen1(group_id2, "/entry0/data_scan/detector_data");
	if (group_id3 < 0) {
		cout << " Error read  H5Gopen1 cannot open /entry0/data_scan/detector_data" << endl;
		return;
	}

	hid_t HDF5_datatype;
	HDF5_datatype = H5T_NATIVE_UINT32;

	if (m_dataRaw ){
			cout << " read2 rawdata" << endl;
			dataset_id = H5Dopen2(group_id3, "/entry0/data_scan/detector_data/data_raw",H5P_DEFAULT);
	}else {
			cout << " read2 CalibrateData" << endl;
			dataset_id = H5Dopen2(group_id3, "/entry0/data_scan/detector_data/data",H5P_DEFAULT);
	}

	if (dataset_id < 0) {
			cout << "Error H5Dopen2   /entry0/data_scan/detector_data/data_raw " << endl;
	}

//	time_period timediff(lStartedTime, microsec_clock::local_time());
//	    float timePrint = ((float)timediff.length().total_seconds()
//					+ (float)timediff.length().fractional_seconds() / 1000000.);
//	    cout << " Time to  READ_DATA2 " << timePrint <<endl;

	status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			m_dread_data);

	if (status < 0){
			cout << "Error H5Dread data_raw" << endl;
	}

//	time_period timediff2(lStartedTime, microsec_clock::local_time());
//    float timePrint2 = ((float)timediff2.length().total_seconds()
//				+ (float)timediff2.length().fractional_seconds() / 1000000.);
//    cout << " Time to  READ_DATA2 " << timePrint2 <<endl;


//  Debug
//	int V_x = 640;
//	int V_y = 256;
//	int V_z = m_total_steps;
//
//	int dres1_sum_raw[m_total_steps];
//
//	for (int kz = 0; kz < 3; kz++) {
//		dres1_sum_raw[kz] = 0;
//			for (int kl = 0; kl < (V_y * V_x); kl++) {
//				dres1_sum_raw[kz] = dres1_sum_raw[kz] + m_dread_data[kl + kz * (V_y * V_x)];
//			}
//			cout << " frame_raw[" << kz << "]sum " << dres1_sum_raw[kz] << "---";
//	cout <<  endl << " End Data raw"<< endl;
//
//	}
}




