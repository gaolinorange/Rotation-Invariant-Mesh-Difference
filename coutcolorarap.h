#pragma once
#include <iostream>
#include <windows.h>
#include <string>

namespace coutcmd{
inline std::ostream& blue(std::ostream &s)
{
	HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE); 
	SetConsoleTextAttribute(hStdout, FOREGROUND_BLUE
		|FOREGROUND_GREEN|FOREGROUND_INTENSITY);
	return s;
}

inline std::ostream& red(std::ostream &s)
{
	HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
	SetConsoleTextAttribute(hStdout,FOREGROUND_RED|FOREGROUND_INTENSITY);
	return s;
}

inline std::ostream& green(std::ostream &s)
{
	HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE); 
	SetConsoleTextAttribute(hStdout, 
		FOREGROUND_GREEN|FOREGROUND_INTENSITY);
	return s;
}

inline std::ostream& yellow(std::ostream &s)
{
	HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE); 
	SetConsoleTextAttribute(hStdout, 
		FOREGROUND_GREEN|FOREGROUND_RED|FOREGROUND_INTENSITY);
	return s;
}

inline std::ostream& white(std::ostream &s)
{
	HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE); 
	SetConsoleTextAttribute(hStdout, 
		FOREGROUND_RED|FOREGROUND_GREEN|FOREGROUND_BLUE);
	return s;
}



struct color 
{
	color(WORD attribute):m_color(attribute){};
	WORD m_color;
};

template <class _Elem, class _Traits>
std::basic_ostream<_Elem,_Traits>& 
	operator<<(std::basic_ostream<_Elem,_Traits>& i, color& c)
{
	HANDLE hStdout=GetStdHandle(STD_OUTPUT_HANDLE); 
	SetConsoleTextAttribute(hStdout,c.m_color);
	return i;
}


inline bool multistringsplit(const std::string& input, std::vector<std::string>& stringsplited, char termin  = ' ')
{
	if (input.empty())
	{
		return false;
	}
	stringsplited.clear();
	int inputLength = input.length();
	std::string tmp;
	for (int i=0;i<inputLength;i++)
	{
		if (input[i]==termin )
		{
			if ( i+1< inputLength)
			{
				if (!tmp.empty() && input[i+1]!=termin )
				{
					stringsplited.push_back(tmp);
					tmp.clear();
				}
			}
		}
		else
		{
			tmp+=input[i];
		}
	}
	if (!tmp.empty())
	{
		stringsplited.push_back(tmp);
	}
	return true;
}

}