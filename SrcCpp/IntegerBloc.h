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



#ifndef INTEGERBLOC_H
#define INTEGERBLOC_H

#include "ArgumentHandler.h"

class  IntegerBloc {

private:
	static IntegerBloc *_singleton;
	/*!
	 * \Constructor
	 */
	IntegerBloc();
	~IntegerBloc();
public:
	static IntegerBloc * getIntegerBlocInstance();
	static void reset();
	void getIntegerBloc(int p3[], int size) ;

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

	 int m_nvers;
	 int m_ntype;
	 int m_kctrl;
	 int m_manip ;
	 int m_nbang;
	 int m_nkmes;
	 int m_npdone ;
	 int m_jcode ;
	 int m_ipara ;
	 int m_ianal;

	 int m_imode ;
	 int m_itgv ;
	 int m_iregul  ;
	 int m_ivolt ;
	 int m_naxe;
	 int m_npstart ;
	 int m_ilast1   ;
	 int m_isa ;
	 int m_flgkif  ;
	 int m_ih;

	 int m_ik  ;
	 int m_nbsqs ;
	 int m_nb_det  ;
	 int m_nbdata ;
	 int m_icdesc1 ;
	 int m_icdesc2 ;
	 int m_icdesc3 ;
	 int m_icdesc4  ;
	 int m_icdesc5;
	 int m_icdesc6;
	 int m_icdesc7;


};


#endif
