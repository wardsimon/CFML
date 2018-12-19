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

#ifndef NAMESERVER_H
#define NAMESERVER_H

#include <string>
#include <iostream>

// Absolute path for ICS configuration
static const std::string CONF_PATH = ".nomadserver/" ;
static const std::string USERFILES_CONF_PATH = CONF_PATH + "UserFiles/" ;

static const std::string DICO_NAME = "NewAscii2Nexus.txt";
static const std::string CHECK_FILE= CONF_PATH + "CheckFile.txt";
static const std::string DICO_PATH = CONF_PATH + DICO_NAME;

static bool VERBOSE_MODE;
static bool REPLACE_MODE;

// By default create Ics server, option argument -env
static std::string PrefixNameServer="NomadServer";
// Absolute path for ICS server configuration
static std::string SERVER_CONF_PATH = CONF_PATH + "Config/";
// Absolute path for ICS server configuration
static const std::string DATA_FILE_MANAGER_TYPE = SERVER_CONF_PATH + "DataFileManagerType.txt";
// Absolute path for ICS server configuration
static const std::string COMMAND_LINE_TYPE = SERVER_CONF_PATH + "CommandLineType.txt";

#endif
