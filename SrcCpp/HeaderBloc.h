
#ifndef HEADERBLOC_H
#define HEADERBLOC_H

#include "ArgumentHandler.h"

class  HeaderBloc {
private:
	static HeaderBloc *_singleton;
	/*!
	 * \Constructor
	 */
	HeaderBloc();
	~HeaderBloc();

public:
	static HeaderBloc * getHeaderBlocInstance();
	static void reset();
	/*!
	 * \brief Getter of actual time of the counter is seconds
	 *
	 * \return actual time of the counter is seconds
	 */
	 int   save( ) ;
	 int   read( ) ;

	 int getHeaderNumor();
	 int getHeaderScanType();
	 const std::string& getHeaderInstrName();
	 const std::string& getHeaderSubT();
	 const std::string& getHeaderUser();
	 const std::string& getHeaderLocalContact();
	 const std::string& getHeaderDate();

//	 double mean_of();
//	 double variance_of();
//	 double standard_deviation_of();

private:

	 ArgumentHandler * m_AH;

	 int m_RecordNumber;
	 int m_ScanType;

	 std::string m_instname;
	 std::string m_username;
	 std::string m_localcontact;
	 std::string m_endtime;

	 std::string m_subtitle;
};



#endif
