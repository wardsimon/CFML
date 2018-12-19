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



#ifndef DATANEXUSLIB_H
#define DATANEXUSLIB_H

#include "ArgumentHandler.h"
#include "./napi.h"


class  DataNexusLib {

public:

	static  int  OpenPath(NXhandle myfile_id,std::string path1);
	static  int  OpenPath(NXhandle myfile_id,std::string path1,std::string path2);
	static  int  OpenPath(NXhandle myfile_id,std::string path1,std::string path2,std::string path3);


	static std::string getStringValue(NXhandle myfile_id,std::string champ);
	static int getIntValue(NXhandle myfile_id,std::string champ);
	static float getFloatValue(NXhandle myfile_id,std::string champ);
	static float getFloatIndexValue(NXhandle myfile_id,std::string champ,int index);
private:

};


#endif
