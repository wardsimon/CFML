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
#include "./FloatBloc.h"
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

FloatBloc* FloatBloc::_singleton = 0;


FloatBloc::FloatBloc() {
	m_AH = ArgumentHandler::getArgumentHandlerInstance();
}

FloatBloc::~FloatBloc() {

}

FloatBloc * FloatBloc::getFloatBlocInstance() {
	if (_singleton == 0) {
		_singleton = new FloatBloc();
	}
	return _singleton;

}

void FloatBloc::getFloatBloc(float p4[], int size) {
//Testing
//	for (int ki = 0 ; ki<size ; ki++){
//	    cout << "cpp par fc_inttab = " << p4[ki] << endl;
//	    p4[ki] = 11.1*(ki+1);
//	  }
	p4[1-1]=m_qH;
	p4[2-1]=m_qK;
	p4[3-1]=m_qL;
	p4[4-1]= m_phi;
	p4[5-1]=m_chi;
	p4[6-1]= m_omega;
	p4[7-1]= m_2theta;
	p4[8-1]= m_psi;

	p4[9-1]= m_ub11 ; p4[10-1]= m_ub12 ; p4[11-1]= m_ub13;
	p4[12-1]= m_ub21 ; p4[13-1]= m_ub22 ; p4[14-1]= m_ub23;
	p4[15-1]= m_ub31 ; p4[16-1]= m_ub32 ; p4[17-1]= m_ub33;
//	cout << " m_ub1 " <<m_ub11 << "  " <<m_ub12 << "  " << m_ub13<< endl;
//	cout << " m_ub2 " <<m_ub21 << "  " <<m_ub22 << "  " << m_ub23<< endl;
//	cout << " m_ub3 " <<m_ub31 << "  " <<m_ub32 << "  " << m_ub33<< endl;
	p4[18-1] = m_wavelength;
	p4[22-1]= m_Hmax ; p4[23-1]= m_Kmax ; p4[24-1]= m_Lmax;
	p4[25-1]= m_DeltaH ;   p4[26-1]= m_DeltaK ;  p4[27-1]= m_DeltaL;
	p4[30-1]= m_Ddetector;
	p4[36-1]= m_scan_start;  p4[37-1]= m_scan_step;  p4[38-1]=m_scan_width;
	p4[39-1]=m_preset;
	p4[43-1]=m_couplingfactor;
	p4[46-1]= m_Tset;p4[47-1]= m_Treg;p4[48-1]= m_Tsample;p4[49-1]=m_Voltmeter ;p4[50-1]=m_Magfield ;


//    n%HMin=rvalues(1:3)          ! HKL min
//            n%angles=rvalues(4:8)        ! Phi, Chi, Omega, 2Theta, Psi
//            n%ub(1,:)=rvalues(9:11)      !
//            n%ub(2,:)=rvalues(12:14)     ! UB Matrix
//            n%ub(3,:)=rvalues(15:17)     !
//            n%wave=rvalues(18)           ! Wavelength
//            n%HMax=rvalues(22:24)        ! HKL max
//            n%dh=rvalues(25:27)          ! Delta HKL
//            n%dist=rvalues(30)           ! distance
//            n%scans=rvalues(36:38)       ! Scan start, Scan step, Scan width
//            n%preset=rvalues(39)         ! Preset
//            n%cpl_fact=rvalues(43)       ! Coupling factor
//            n%conditions=rvalues(46:50)  ! Temp-s, Temp-r, Temp-sample. Voltmeter, Mag.Field


}


int FloatBloc::read() {

	string nxFileName = m_AH->getNexusFileName();
//	cout << " FloatBloc::read  " << nxFileName << endl;
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
// L1
	DataNexusLib::OpenPath(myfile_id, "/entry0/", "instrument/", "SingleCrystalSettings/");
	champ = "actual_reflection";
	m_qH = DataNexusLib::getFloatIndexValue(myfile_id, champ, 0);
	m_qK = DataNexusLib::getFloatIndexValue(myfile_id, champ, 1);
	m_qL = DataNexusLib::getFloatIndexValue(myfile_id, champ, 2);
//	cout << " m_qL  "<< m_qL << endl;
	DataNexusLib::OpenPath(myfile_id, "/entry0/", "instrument/", "phi/");
	champ = "value";
	m_phi = DataNexusLib::getFloatValue(myfile_id, champ);
	DataNexusLib::OpenPath(myfile_id, "/entry0/", "instrument/", "chi/");
	champ = "value";
	m_chi = DataNexusLib::getFloatValue(myfile_id, champ);
// L2
	DataNexusLib::OpenPath(myfile_id, "/entry0/", "instrument/", "omega/");
	champ = "value";
	m_omega = DataNexusLib::getFloatValue(myfile_id, champ);
	DataNexusLib::OpenPath(myfile_id, "/entry0/", "instrument/", "gamma/");
	champ = "value";
	m_2theta = DataNexusLib::getFloatValue(myfile_id, champ);
	DataNexusLib::OpenPath(myfile_id, "/entry0/", "instrument/", "SingleCrystalSettings/");
	champ = "psi";
	m_psi = DataNexusLib::getFloatValue(myfile_id, champ);

	DataNexusLib::OpenPath(myfile_id, "/entry0/", "instrument/", "SingleCrystalSettings/");
	champ = "orientation_matrix";
	m_ub11 = DataNexusLib::getFloatIndexValue(myfile_id, champ, 0);
	m_ub12 = DataNexusLib::getFloatIndexValue(myfile_id, champ, 1);
	m_ub13 = DataNexusLib::getFloatIndexValue(myfile_id, champ, 2);
// L3
	m_ub21 = DataNexusLib::getFloatIndexValue(myfile_id, champ, 3);
	m_ub22 = DataNexusLib::getFloatIndexValue(myfile_id, champ, 4);
	m_ub23 = DataNexusLib::getFloatIndexValue(myfile_id, champ, 5);
//L4
	m_ub31 = DataNexusLib::getFloatIndexValue(myfile_id, champ, 6);
	m_ub32 = DataNexusLib::getFloatIndexValue(myfile_id, champ, 7);
	m_ub33 = DataNexusLib::getFloatIndexValue(myfile_id, champ, 8);
	DataNexusLib::OpenPath(myfile_id, "/entry0/");
	champ = "wavelength";
	m_wavelength = DataNexusLib::getFloatValue(myfile_id, champ);
	m_spare1 = 0.00;
	m_danalyser = 0.00;
//	DataNexusLib::OpenPath(myfile_id, "/entry0/", "instrument/", "/Monochromator");
//	champ = "dspacing";
//	m_danalyser= DataNexusLib::getFloatValue(myfile_id, champ);
//L5
	m_energy = 0.00;
	m_Hmax = 0.00;
	m_Kmax = 0.00;
	m_Lmax = 0.00;
	m_DeltaH = 0.00;
//L6
	m_DeltaK = 0.00;
	m_DeltaL = 0.00;
	m_Deltaenergy = 0.00;
	m_Ki = 0.00;
	m_Ddetector = 764.00;
//L7
	m_xoff = 0.00;
	m_zoff = 0.00;
	m_radius = 0.00;
	m_yoff = 0.00;
	m_attenuat = 0.00;
//L8
	bool is_scan = true;
	if (is_scan) {
		DataNexusLib::OpenPath(myfile_id, "/entry0/", "instrument/", "ScanInfo/");
		champ = "first_start_value";
		m_scan_start = DataNexusLib::getFloatValue(myfile_id, champ);
		champ = "first_step_value";
		m_scan_step = DataNexusLib::getFloatValue(myfile_id, champ);
		champ = "first_width_value";
		m_scan_width = DataNexusLib::getFloatValue(myfile_id, champ);
	} else {
		m_scan_start = 0.00;
		m_scan_step = 0.00;
		m_scan_width = 0.00;
	}
	DataNexusLib::OpenPath(myfile_id, "/entry0/");
	champ = "preset";
	m_preset = DataNexusLib::getFloatValue(myfile_id, champ);
	m_addbkgstep = 0.00;
//L9
	m_addbkgwidth = 0.00;
	m_addbkgpreset = 0.00;
	m_couplingfactor = 0.00;
	m_spare2 = 0.00;
	m_spare3 = 0.00;
//L10
	DataNexusLib::OpenPath(myfile_id, "/entry0/", "sample/");
	champ = "setpoint_temperature";
	m_Tset = DataNexusLib::getFloatValue(myfile_id, champ);
	champ = "regulation_temperature";
	m_Treg = DataNexusLib::getFloatValue(myfile_id, champ);
	champ = "temperature";
	m_Tsample = DataNexusLib::getFloatValue(myfile_id, champ);
	//champ = "preset";
	m_Voltmeter = 0.00;
//	champ = "field";
//	m_Magfield = DataNexusLib::getFloatValue(myfile_id, champ);
	m_Magfield = 0.0;
	if (NXclose(&myfile_id) != NX_OK)
		return -1;

}

int FloatBloc::save() {
	cout << " FloatBloc::save  " << endl;

	string dataFileName = m_AH->getAsciiFileName();

	ofstream dataFile(dataFileName.c_str(), std::ofstream::out | std::ofstream::app);

	dataFile << string(80, 'F') << "\n";
	int li_co = 50;
	int version = 10;
	dataFile << boost::format("%8i%8i") % li_co % version;
	dataFile << string(64, ' ') << "\n";

	dataFile << "              qH              qK              qL             phi             chi" << "\n";
	dataFile << "           omega  2theta (gamma)             psi         ub(1 1)         ub(1 2)" << "\n";
	dataFile << "         ub(1 3)         ub(2 1)         ub(2 2)         ub(2 3)         ub(3 1)" << "\n";
	dataFile << "         ub(3 2)         ub(3 3)      wavelength                       danalyser" << "\n";
	dataFile << "          energy            Hmax            Kmax            Lmax          DeltaH" << "\n";
	dataFile << "          DeltaK          DeltaL     Deltaenergy         Ki (Kf)       Ddetector" << "\n";
	dataFile << "            xoff            zoff          radius            yoff        attenuat" << "\n";
	dataFile << "      scan start       scan step      scan width          preset    add.bkg.step" << "\n";
	dataFile << "   add.bkg.width  add.bkg.preset  couplingfactor         (spare)         (spare)" << "\n";
	dataFile << "        Tset (k)        Treg (k)     Tsample (k)       Voltmeter       Mag.field" << "\n";



	ostringstream temp1 ;
	temp1 << std::right << std::setfill(' ') << scientific << setprecision(8) << setw(16)
		  << setw(16) <<m_qH  << setw(16) << m_qK  << setw(16) << m_qL << setw(16) << m_phi	<<  setw(16)	<< m_chi ;
	dataFile << temp1.str() << string(80 - temp1.str().length(), ' ') << "\n";
	temp1.str( std::string() );temp1.clear();
	temp1 << std::right << std::setfill(' ') << scientific << setprecision(8) << setw(16)
		  << setw(16) <<m_omega  << setw(16) <<  m_2theta << setw(16) << m_psi << setw(16)<< m_ub11	<<  setw(16)	<< m_ub12 ;
	dataFile << temp1.str() << string(80 - temp1.str().length(), ' ') << "\n";
	temp1.str( std::string() );temp1.clear();
	temp1 << std::right << std::setfill(' ') << scientific << setprecision(8) << setw(16)
		  << setw(16) <<m_ub13  << setw(16) << m_ub21  << setw(16) << m_ub22 << setw(16)<<m_ub23 	<<  setw(16)	<< m_ub31 ;
	dataFile << temp1.str() << string(80 - temp1.str().length(), ' ') << "\n";
	temp1.str( std::string() );temp1.clear();
	temp1 << std::right << std::setfill(' ') << scientific << setprecision(8) << setw(16)
		  << setw(16) << m_ub32 << setw(16) <<  m_ub33 << setw(16) << m_wavelength << setw(16)<< m_spare1	<<  setw(16)	<< m_danalyser ;
	dataFile << temp1.str() << string(80 - temp1.str().length(), ' ') << "\n";
	temp1.str( std::string() );temp1.clear();
	temp1 << std::right << std::setfill(' ') << scientific << setprecision(8) << setw(16)
		  << setw(16) << m_energy << setw(16) << m_Hmax  << setw(16) << m_Kmax << setw(16)<< m_Lmax <<  setw(16)	<< m_DeltaH ;
	dataFile << temp1.str() << string(80 - temp1.str().length(), ' ') << "\n";
	temp1.str( std::string() );temp1.clear();
	temp1 << std::right << std::setfill(' ') << scientific << setprecision(8) << setw(16)
		  << setw(16) << m_DeltaK << setw(16) << m_DeltaL  << setw(16) << m_Deltaenergy << setw(16)<< m_Ki	<<  setw(16)	<< m_Ddetector ;
	dataFile << temp1.str() << string(80 - temp1.str().length(), ' ') << "\n";
	temp1.str( std::string() );temp1.clear();
	temp1 << std::right << std::setfill(' ') << scientific << setprecision(8) << setw(16)
		  << setw(16) << m_xoff << setw(16) << m_zoff  << setw(16) << m_radius << setw(16) << m_yoff	<<  setw(16)	<<m_attenuat  ;
	dataFile << temp1.str() << string(80 - temp1.str().length(), ' ') << "\n";
	temp1.str( std::string() );temp1.clear();
	temp1 << std::right << std::setfill(' ') << scientific << setprecision(8) << setw(16)
		  << setw(16) << m_scan_start << setw(16) << m_scan_step  << setw(16) << m_scan_width << setw(16)<<	m_preset<<  setw(16)	<< m_addbkgstep ;
	dataFile << temp1.str() << string(80 - temp1.str().length(), ' ') << "\n";
	temp1.str( std::string() );temp1.clear();
	temp1 << std::right << std::setfill(' ') << scientific << setprecision(8) << setw(16)
		  << setw(16) << m_addbkgwidth << setw(16) <<  m_addbkgpreset << setw(16) << m_couplingfactor << setw(16)<< m_spare2	<<  setw(16)	<<m_spare3  ;
	dataFile << temp1.str() << string(80 - temp1.str().length(), ' ') << "\n";
	temp1.str( std::string() );temp1.clear();
	temp1 << std::right << std::setfill(' ') << scientific << setprecision(8) << setw(16)
		  << setw(16) << m_Tset << setw(16) << m_Treg  << setw(16) << m_Tsample << setw(16)<< m_Voltmeter	<<  setw(16)	<< m_Magfield ;
	dataFile << temp1.str() << string(80 - temp1.str().length(), ' ') << "\n";
	temp1.str( std::string() );temp1.clear();

	dataFile.close();
	return 0;
}




