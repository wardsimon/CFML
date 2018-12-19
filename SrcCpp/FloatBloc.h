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



#ifndef FLOATBLOC_H
#define FLOATBLOC_H

#include "ArgumentHandler.h"

class FloatBloc {

private:
	static FloatBloc *_singleton;
	/*!
	 * \Constructor
	 */
	FloatBloc();
	~FloatBloc();

public:
	static FloatBloc * getFloatBlocInstance();
	static void reset();
	void getFloatBloc(float p4[], int size) ;
	/*!
	 * \brief Getter of actual time of the counter is seconds
	 *
	 * \return actual time of the counter is seconds
	 */
	 int   save( ) ;
	 int   read( ) ;
//	 double mean_of();
//	 double variance_of();
//	 double standard_deviation_of();


private:

	 ArgumentHandler * m_AH;

    float m_qH ;
    float m_qK ;
    float m_qL ;
    float m_phi ;
    float m_chi;
	//
	float m_omega;
	float m_2theta ;
	float m_psi;
	float m_ub11;
	float m_ub12;
    //
    float m_ub13;
    float m_ub21;
    float m_ub22;
    float m_ub23;
    float m_ub31;
    //
	float m_ub32 ;
	float m_ub33  ;
	float m_wavelength  ;
	float m_spare1  ;
	float m_danalyser;
	//
	float m_energy  ;
	float m_Hmax    ;
	float m_Kmax  ;
	float m_Lmax  ;
	float m_DeltaH;
	//
	float m_DeltaK ;
	float m_DeltaL  ;
	float m_Deltaenergy ;
	float m_Ki ;
	float m_Ddetector;
	//
	float m_xoff ;
	float m_zoff ;
	float m_radius ;
	float m_yoff   ;
	float m_attenuat;
	//
	float m_scan_start;
	float m_scan_step ;
	float m_scan_width  ;
	float m_preset  ;
	float m_addbkgstep;
	//
	float m_addbkgwidth;
	float m_addbkgpreset ;
	float m_couplingfactor  ;
	float m_spare2  ;
	float m_spare3;
    //
	float m_Tset;
	float m_Treg;
	float m_Tsample;
	float m_Voltmeter  ;
	float m_Magfield;



};


#endif
