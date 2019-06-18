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
#include "./IntegerBloc.h"
#include "./DataNexusLib.h"

#include <fstream>
#include <napi.h>
#include <boost/format.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <sstream>
#include <iostream>
#include <iomanip>

using namespace std;

const string NAME_DATA_SCAN 		 = "/entry0/data_scan";

IntegerBloc* IntegerBloc::_singleton = 0;


IntegerBloc::IntegerBloc() {
	m_AH = ArgumentHandler::getArgumentHandlerInstance();
}

IntegerBloc::~IntegerBloc() {

}

IntegerBloc * IntegerBloc::getIntegerBlocInstance() {
	if (_singleton == 0) {
		_singleton = new IntegerBloc();
	}
	return _singleton;

}

void IntegerBloc::getIntegerBloc(int p3[], int size) {

//	// For testing
//	  for (int ki = 0 ; ki<size ; ki++){
//	    cout << "cpp par fc_inttab = " << p3[ki] << endl;
//	    p3[ki] = ki+100;
//	  }
//      n%manip=ivalues(4)              ! 1: 2Theta, 2: Omega, 3:Chi, 4: Phi
//      n%nbang=ivalues(5)              ! Total number of angles moved during scan
//      n%nframes=ivalues(7)            ! Measured Frames. In general equal to those prescripted
//      n% icalc=ivalues(9)
//      n%nbdata=ivalues(24)            ! Number of Points
//      n%icdesc(1:7)=ivalues(25:31)
	p3[4-1] =m_manip;
	p3[5-1] =m_nbang;
	p3[7-1] =m_npdone;
	p3[9-1] =m_ipara;
	p3[24-1] =m_nbdata;

	p3[25-1] =m_icdesc1;
	p3[26-1] =m_icdesc2;
	p3[27-1] =m_icdesc3;
	p3[28-1] =m_icdesc4;
	p3[29-1] =m_icdesc5;
	p3[30-1] =m_icdesc6;
	p3[31-1] =m_icdesc7;


}

int IntegerBloc::read() {

	string nxFileName = m_AH->getNexusFileName();
//	cout << " IntegerBloc::read  " << nxFileName << endl;
	NXstatus myret;
	NXhandle myfile_id; //TODO: decommenter toute la partie ci dessous et commenter l'ancienne déclaration de mode
	myret = NXopen(nxFileName.c_str(), NXACC_READ, &myfile_id);
	if (myret != NX_OK)
		return -1;

	string nPath;
	string nPath2;
	string nPath3;

	string champ;
	string res;

	m_nvers = 4;
	m_ntype = 2;
	m_kctrl = 4;
	m_manip = 2;

	DataNexusLib::OpenPath(myfile_id, "/entry0/", "instrument/", "SingleCrystalSettings/");
	champ = "nbAngles";
	//m_nbang = DataNexusLib::getIntValue(myfile_id, champ);;
	m_nbang  = 1;


	m_nkmes = 1;
	m_npdone = 1;

	DataNexusLib::OpenPath(myfile_id,NAME_DATA_SCAN   );
	champ = "actual_step";
	m_npdone  = DataNexusLib::getIntValue(myfile_id, champ);
	champ = "total_steps";
	m_nkmes  = DataNexusLib::getIntValue(myfile_id, champ);


	DataNexusLib::OpenPath(myfile_id, "/entry0/");
	champ = "mode";
	m_jcode = DataNexusLib::getIntValue(myfile_id, champ);
	m_ipara = 1;
	m_ianal = 0;
	m_imode = 0;
	m_itgv = 0;
	DataNexusLib::OpenPath(myfile_id, "/entry0/", "sample/");
	champ = "use_temp";
	m_iregul = DataNexusLib::getIntValue(myfile_id, champ);
	m_ivolt = 0;
	DataNexusLib::OpenPath(myfile_id, "/entry0/", "instrument/", "SingleCrystalSettings/");
	champ = "naxe";
	m_naxe = DataNexusLib::getIntValue(myfile_id, champ);

	m_npstart = 0;
	m_ilast1 = 0;
	m_isa = 0;
	m_flgkif = 0;
	m_ih = 0;
	m_ik = 0;
	m_nbsqs = 0;
	m_nb_det = 1;

	DataNexusLib::OpenPath(myfile_id, "/entry0/", "instrument/", "Det1/");
	champ = "detsize";
	m_nbdata = DataNexusLib::getIntValue(myfile_id, champ);

//	DataNexusLib::OpenPath(myfile_id, "/entry0/", "instrument/", "SingleCrystalSettings/");
//	champ = "icdesc1";
//	m_icdesc1 = DataNexusLib::getIntValue(myfile_id, champ);
//	DataNexusLib::OpenPath(myfile_id, "/entry0/", "instrument/", "SingleCrystalSettings/");
//	champ = "icdesc2";
//	m_icdesc2 = DataNexusLib::getIntValue(myfile_id, champ);
//	DataNexusLib::OpenPath(myfile_id, "/entry0/", "instrument/", "SingleCrystalSettings/");
//	champ = "icdesc3";
//	m_icdesc3 = DataNexusLib::getIntValue(myfile_id, champ);
//	DataNexusLib::OpenPath(myfile_id, "/entry0/", "instrument/", "SingleCrystalSettings/");
//	champ = "icdesc4";
//	m_icdesc4 = DataNexusLib::getIntValue(myfile_id, champ);


	m_icdesc1 = 0;
	m_icdesc2 = 0;
	m_icdesc3 = 0;
	m_icdesc4 =0;
	m_icdesc5 = 0;
	m_icdesc6 = 0;
	m_icdesc7 = 0;

	if (NXclose(&myfile_id) != NX_OK)
		return -1;

}

int IntegerBloc::save() {
	cout << " IntegerBloc::save  " << endl;

	string dataFileName = m_AH->getAsciiFileName();

	ofstream dataFile(dataFileName.c_str(), std::ofstream::out | std::ofstream::app);

	dataFile << string(80, 'I') << "\n";
	int li_co = 31;
	int version = 4;
	dataFile << boost::format("%8i%8i") % li_co % version;
	dataFile << string(64, ' ') << "\n";

	dataFile << "   nvers   ntype   kctrl   manip   nbang   nkmes  npdone   jcode   ipara   ianal" << "\n";
	dataFile << "   imode    itgv  iregul   ivolt    naxe npstart  ilast1     isa  flgkif      ih" << "\n";
	dataFile << "      ik   nbsqs  nb_det  nbdata icdesc1 icdesc2 icdesc3 icdesc4 icdesc5 icdesc6" << "\n";
	dataFile << " icdesc7                                                                        " << "\n";

	ostringstream temp1, temp2, temp3, temp4;

	temp1 << std::right  << std::setfill(' ') << setw(8) << m_nvers << setw(8) << m_ntype << setw(8) << m_kctrl << setw(8)
			<< m_manip << setw(8) << m_nbang << setw(8) << m_nkmes << setw(8) << m_npdone << setw(8) << m_jcode << setw(8)
			<< m_ipara << setw(8) << m_ianal;
	dataFile << temp1.str() << string(80 - temp1.str().length(), ' ') << "\n";

	temp2 << std::right << std::setfill(' ') << setw(8) << m_imode << setw(8) << m_itgv << setw(8) << m_iregul << setw(8)
			<< m_ivolt << setw(8) << m_naxe << setw(8) << m_npstart << setw(8) << m_ilast1 << setw(8) << m_isa << setw(8)
			<< m_flgkif << setw(8) << m_ih;
	dataFile << temp2.str() << string(80 - temp2.str().length(), ' ') << "\n";

	temp3 << std::right << std::setfill(' ') << setw(8) << m_ik << setw(8) << m_nbsqs << setw(8) << m_nb_det << setw(8)
			<< m_nbdata << setw(8) << m_icdesc1 << setw(8) << m_icdesc2 << setw(8) << m_icdesc3 << setw(8) << m_icdesc4
			<< setw(8) << m_icdesc5 << setw(8) << m_icdesc6;
	dataFile << temp3.str() << string(80 - temp3.str().length(), ' ') << "\n";

	temp4 << std::right << std::setfill(' ') << setw(8) << m_icdesc7;
	dataFile << temp4.str() << string(80 - temp4.str().length(), ' ') << "\n";


	dataFile.close();
	return 0;

}

