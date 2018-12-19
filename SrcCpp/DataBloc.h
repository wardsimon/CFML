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



#ifndef DATABLOC_H
#define DATABLOC_H

#include "ArgumentHandler.h"
#include <napi.h>


class  DataBloc {

private:
	static DataBloc *_singleton;



	/*!
	 * \Constructor
	 */
	DataBloc();
	~DataBloc();

public:

	static DataBloc * getDataBlocInstance();
	static void reset();
	void getDataBlocTime(float p4[], int size) ;
	void getDataBlocMoni(float p4[], int size) ;
	void getDataBlocTotalCount(float p4[], int size) ;
	void getDataBlocAngle1(float p4[], int size) ;
//	void getDataBlocDataFull(float p4[], int size) ;
	void getDataBlocDataFull(int p4[], int size) ;

	/*!
	 * \brief Getter of actual time of the counter is seconds
	 *
	 * \return actual time of the counter is seconds
	 */
	 int   save( ) ;
	 int   read( ) ;
	 void setDataMode( bool dataRaw);

//	 double mean_of();
//	 double variance_of();
//	 double standard_deviation_of();


private:

	 void read_scan_var();
	 void decode_var_name();
	 void   read_data();
	 void   read_data2(int *pp);

	 void fastMemcpy(void *pvDest, void *pvSrc, size_t nBytes) ;

	 ArgumentHandler * m_AH;
	 std::string m_nxFileName;

	 bool  m_is_scan;
	 int m_actual_step;
	 int m_total_steps;
	 int m_RecordNumber;

	 float* m_time;
	 float* m_monitor;
	 float* m_totaCount;
	 float* m_angles1;
	 float* m_angles2;
	 float* m_angles3;
	 float* m_angles4;

	 int m_index_time;
	 int m_index_monitor;
	 int m_index_totalCount;
	 int m_index_angles1;
	 int m_index_angles2;
	 int m_index_angles3;
	 int m_index_angles4;

	 int m_nb_prop;
	 int m_prop_index[7];

	 bool m_dataRaw;

	 int *m_dread_data;


};


#endif
