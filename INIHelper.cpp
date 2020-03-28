//#include "stdafx.h"
#include "INIHelper.h"
#include <assert.h>

namespace Base
{
	namespace Utility
	{
		MyINIParser::MyINIParser()
		{
			dict = NULL;
		}

		MyINIParser::~MyINIParser()
		{
			if (dict!=NULL)
			{
				iniparser_freedict(dict);
			}
		}

		bool MyINIParser::LoadINIFile(const char *fname)
		{
			dict = iniparser_load(fname);
			return (dict!=NULL);
		}

		int MyINIParser::GetNumSections()
		{
			assert(dict!=NULL);
			return iniparser_getnsec(dict);
		}

		char* MyINIParser::GetSectionName(int id)
		{
			assert(dict!=NULL);
			return iniparser_getsecname (dict, id);
		}

		char* MyINIParser::GetString(const char* key, char* def)
		{
			assert(dict!=NULL);
			return iniparser_getstring(dict, key, def);
		}

		int MyINIParser::GetInt(const char* key, int notfound)
		{
			assert(dict!=NULL);
			return iniparser_getint(dict, key, notfound);
		}

		double MyINIParser::GetDouble(char* key, double notfound)
		{
			assert(dict!=NULL);
			return iniparser_getdouble(dict, key, notfound);
		}

		bool MyINIParser::GetBoolean(const char* key, int notfound)
		{
			assert(dict!=NULL);
			return (bool)iniparser_getboolean(dict, key, notfound);
		}
	}
}
