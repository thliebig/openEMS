/*
*	Copyright (C) 2010 Thorsten Liebig (Thorsten.Liebig@gmx.de)
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*!
\class ErrorMsg
\author Thorsten Liebig
\version $Revision: 1.2 $
\date $Date: 2006/01/25 11:47:07 $
*/

#ifndef _ERRORMSG_H_
#define _ERRORMSG_H_

class ErrorMsg
{
public:
	///Constructor defines number of error messages
	ErrorMsg(unsigned int NoMessage=0);
	///Deconstructor
	virtual ~ErrorMsg();
	///Methode for defining error messages
	/*! \param nr Number of defining error message \param *Message Set error message string  \sa Error */
	void SetMsg(unsigned int nr, const char *Message);
	///Call an error message. Will exit the program!
	/*! \param nr Number of called error message. default is 0 \sa SetMsg*/
	void Error(unsigned int nr=0,char *chAddMsg=0);

	void Error(unsigned int nr,int addNr);

protected:
	void ownError(void);
	unsigned int NoMsg;
	char **Msg;
};

#endif //_ERRORMSG_H_
