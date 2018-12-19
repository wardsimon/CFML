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
#include "./HeaderBloc.h"
#include "./DataNexusLib.h"

#include <fstream>
#include <boost/format.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <sstream>
#include <iostream>
#include <iomanip>

using namespace std;

HeaderBloc* HeaderBloc::_singleton = 0;

HeaderBloc::HeaderBloc() {
	m_AH = ArgumentHandler::getArgumentHandlerInstance();
}

HeaderBloc::~HeaderBloc() {

}

HeaderBloc * HeaderBloc::getHeaderBlocInstance() {
	if (_singleton == 0) {
		_singleton = new HeaderBloc();
	}
	return _singleton;

}

int HeaderBloc::read() {

	string nxFileName = m_AH->getNexusFileName();
	//cout << " HeaderBloc::read 01 " << nxFileName << endl;

	NXstatus myret;
	NXhandle myfile_id; //TODO: decommenter toute la partie ci dessous et commenter l'ancienne déclaration de mode
	int NXtype, NXrank;
	int NXdims[32];
	void *data_buffer;

	myret = NXopen(nxFileName.c_str(), NXACC_READ, &myfile_id);
	if (myret != NX_OK) {
		cout << " HeaderBloc::read   Failed to open : " << nxFileName << endl;
		return -1;
	}

	string nPath;
	string nPath2;
	string nPath3;

	string champ;
	string res;

	DataNexusLib::OpenPath(myfile_id, "/entry0/");
	DataNexusLib::OpenPath(myfile_id, "instrument/");
	champ = "name";
	m_instname = DataNexusLib::getStringValue(myfile_id, champ);
//	cout << " HeaderBloc::read  m_instname " << m_instname << endl;
	//////////////////
	DataNexusLib::OpenPath(myfile_id, "/entry0/");
	DataNexusLib::OpenPath(myfile_id, "user/");
	champ = "name";
	m_username = DataNexusLib::getStringValue(myfile_id, champ);
	champ = "namelocalcontact";
	m_localcontact = DataNexusLib::getStringValue(myfile_id, champ);
//	cout << " HeaderBloc::read  m_localcontact " << m_localcontact << endl;

	DataNexusLib::OpenPath(myfile_id, "/entry0/");
	champ = "end_time";
	m_endtime = DataNexusLib::getStringValue(myfile_id, champ);
//	cout << " HeaderBloc::read  m_endtime " << m_endtime << endl;

	champ = "experiment_identifier";
	m_subtitle = DataNexusLib::getStringValue(myfile_id, champ);

	champ = "run_number";
	m_RecordNumber = DataNexusLib::getIntValue(myfile_id, champ);

	m_ScanType = 1;

	if (NXclose(&myfile_id) != NX_OK)
		return -1;

}

int HeaderBloc::getHeaderNumor() {
	return m_RecordNumber;
}

int HeaderBloc::getHeaderScanType() {
	return m_ScanType;
}
const string& HeaderBloc::getHeaderInstrName() {
	//cout <<  " getHeaderInstrName  " << m_instname<< endl;
	return m_instname;
}
const string& HeaderBloc::getHeaderSubT() {
	//cout <<  " getHeaderSubT  " << m_instname<< endl;
	return m_subtitle;
}

const std::string& HeaderBloc::getHeaderUser(){
	return m_username;
}
const std::string& HeaderBloc::getHeaderLocalContact(){
	return m_localcontact;
}
const std::string& HeaderBloc::getHeaderDate(){
	return m_endtime;
}


int HeaderBloc::save() {
	cout << " HeaderBloc::save  " << endl;

	string dataFileName = m_AH->getAsciiFileName();

	ofstream dataFile(dataFileName.c_str());

	dataFile << string(80, 'R') << "\n";
	int li_co = 0;
	int version = 4;
	dataFile
			<< boost::format("  %06i%8i%8i") % (int) m_RecordNumber % li_co
					% version;
	dataFile << string(56, ' ') << "\n";

	dataFile << string(80, 'A') << "\n";
	int m_Size = 80;
	dataFile << boost::format("%8i") % m_Size << string(7, ' ') << "1"
			<< string(64, ' ') << "\n";
	dataFile
			<< "Inst User L.C. Date Time                                                        "
			<< "\n";

	ostringstream temp;
	temp << std::left << std::setfill(' ') << setw(4) << m_instname << setw(6)
			<< m_username.substr(0, 6) << setw(4) << m_localcontact.substr(0, 4)
			<< setw(18) << m_endtime.substr(0, 18);
	dataFile << temp.str() << string(80 - temp.str().length(), ' ') << "\n";

	dataFile << string(80, 'A') << "\n";
	m_Size = 80;
	dataFile << boost::format("%8i") % m_Size << string(7, ' ') << "1"
			<< string(64, ' ') << "\n";
	dataFile
			<< "Title                                                                   Scantype"
			<< "\n";
	dataFile << m_subtitle << string(72 - m_subtitle.size(), ' ') << "omega   "
			<< "\n";

	dataFile.close();

}

