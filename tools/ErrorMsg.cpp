#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ErrorMsg.h"

ErrorMsg::ErrorMsg(unsigned int NoMessage)
{
	NoMsg=NoMessage;
	if (NoMsg>0) Msg = new char*[NoMsg]; if (Msg==NULL) { fprintf(stderr,"Memory allocation failed!! exiting... \0"); exit(1); }
	for (unsigned int i=0;i<NoMsg;i++) Msg[i]=NULL;
}

ErrorMsg::~ErrorMsg()
{
	for (unsigned int i=0;i<NoMsg;i++) {delete[] Msg[i]; Msg[i]=NULL;};
	delete[] Msg; Msg=NULL;
}

void ErrorMsg::SetMsg(unsigned int nr, char *Message)
{
	if ((nr<1) || (nr>NoMsg) || (Message==NULL)) ownError();
	Msg[nr-1] = new char[strlen(Message)+1]; if (Msg[nr-1]==NULL) { fprintf(stderr,"Memory allocation failed!! exiting... \0"); exit(1); }
	Msg[nr-1]=strcpy(Msg[nr-1],Message);
}

void ErrorMsg::Error(unsigned int nr,char *chAddMsg)
{
	if ((nr>0) && (nr<=NoMsg))
	{
		if (Msg[nr-1]!=NULL) fprintf(stderr,"%s",Msg[nr-1]);
		else fprintf(stderr,"unkown error occured!! Error code: %d exiting...",nr);
		if (chAddMsg!=NULL) fprintf(stderr,"%s",chAddMsg);
		getchar();
		exit(nr);
	}
	else
	{
		fprintf(stderr,"unkown error occured!! Error code: %d exiting...",nr);
		getchar();
		exit(nr);
	}
}

void ErrorMsg::Error(unsigned int nr,int addNr)
{
	if ((nr>0) && (nr<=NoMsg))
	{
		if (Msg[nr-1]!=NULL) fprintf(stderr,"%s",Msg[nr-1]);
		else fprintf(stderr,"unkown error occured!! Error code: %d exiting...",nr);
		fprintf(stderr,"%d",addNr);
		getchar();
		exit(nr);
	}
	else
	{
		fprintf(stderr,"unkown error occured!! Error code: %d exiting...",nr);
		getchar();
		exit(nr);
	}
}

void ErrorMsg::ownError(void)
{
	fprintf(stdout," Error occured by using Error Message class!! ... exiting...");
	exit(-1);
}
