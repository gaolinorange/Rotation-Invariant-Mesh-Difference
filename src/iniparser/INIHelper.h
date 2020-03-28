// wrapper for the iniparser in third party libs

#pragma once

#include "iniparser.h"

namespace Base
{
	namespace Utility
	{
		// a simple ini parser based on third party lib ini parser
		// NOTE: the *.ini file should end with an extra \r\n, otherwise, the parser may fail to load it.
		class MyINIParser
		{
		public:
			MyINIParser();
			~MyINIParser();

			//////////////////////////////////////////////////////////////////////////
			// service functions
			bool LoadINIFile(const char *fname);    // load ini file
			int GetNumSections();   // get number of sections in the ini file
			char* GetSectionName(int id);   // get the name for section id
			char* GetString(const char* key, char* def);  // get char value
			int GetInt(const char* key, int notfound);   // get int value
			double GetDouble(char* key, double notfound);   // get double value
			bool GetBoolean(const char* key, int notfound);   // get bool value

			// add set and save functions, TODO

		private:
			dictionary* dict;   // data structure from ini parser
		};
	}
}
