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
	void SetMsg(unsigned int nr, char *Message);
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
