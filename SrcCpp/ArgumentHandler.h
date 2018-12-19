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


#ifndef ARGUMENTHANDLER_H
#define ARGUMENTHANDLER_H

#include <iostream>
#include <string>

class ArgumentHandler{
private:

	static ArgumentHandler *_singleton;
	bool verboseMode;
	bool replaceMode;
	std::string asciiFileName;
	std::string nexusFileName;
	std::string targetDirectory;
	std::string dataFileConfig;
	ArgumentHandler();
	~ArgumentHandler();

public:
	static ArgumentHandler * getArgumentHandlerInstance();
	static void reset();
	void verboseOn();
	void replaceOff();
	void setAsciiFileName(std::string name);
	void setNexusFileName(std::string name);
	void setTargetDirectory(std::string name);
	void setDataFileConfig(std::string name);

	bool isVerbose();
	bool isReplace();
	std::string getAsciiFileName();
	std::string getNexusFileName();
	std::string getTargetDirectory();
	std::string getDataFileConfig();
};

#endif
