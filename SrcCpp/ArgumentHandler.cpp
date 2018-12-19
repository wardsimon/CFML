/*
 * Nomad Instrument Control Software
 *
 * Copyright 2011 Institut Laue-Langevin
 *
 * Licensed under the EUPL, Version 1.1 only (the "License");
 * You may not use this work except in compliance with the Licence.
 * You may obtain a copy of the Licence at:
 *
 * http://www.osor.eu/eupl
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the Licence is distributed on an "AS IS" basis,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the Licence for the specific language governing permissions and
 * limitations under the Licence.
 */

#include "NameServer.h"
#include "ArgumentHandler.h"

using namespace std;

ArgumentHandler* ArgumentHandler::_singleton = 0;

ArgumentHandler::ArgumentHandler(){
	verboseMode=false;
	replaceMode=true;
	asciiFileName="";
	nexusFileName="";
	targetDirectory="";
}


ArgumentHandler::~ArgumentHandler(){

}

ArgumentHandler * ArgumentHandler::getArgumentHandlerInstance(){
	if (_singleton == 0){
		_singleton = new ArgumentHandler();
	}
	return _singleton;


}

void ArgumentHandler::reset()
{
  if (_singleton)
  {
    delete _singleton;
    _singleton = 0;
  }
}

void ArgumentHandler::verboseOn(){
	verboseMode = true;
}

void ArgumentHandler::replaceOff(){
	replaceMode = true;
}

void ArgumentHandler::setAsciiFileName(string name){
	asciiFileName = name;
}

void ArgumentHandler::setNexusFileName(string name){
	nexusFileName=name;
}

void ArgumentHandler::setTargetDirectory(string name){
	targetDirectory = name;
}
void ArgumentHandler::setDataFileConfig(string name){
	dataFileConfig=name;
}

bool ArgumentHandler::isVerbose(){
	return verboseMode;
}
bool ArgumentHandler::isReplace(){
	return replaceMode;
}
string ArgumentHandler::getAsciiFileName(){
	return asciiFileName;
}
string ArgumentHandler::getNexusFileName(){
	return nexusFileName;
}
string ArgumentHandler::getTargetDirectory(){
	return targetDirectory;
}
string ArgumentHandler::getDataFileConfig(){
	return dataFileConfig;
}


