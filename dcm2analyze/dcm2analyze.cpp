#include <mlib3.h>
#include "../dcmlib/dcmlib.h"
using namespace dcmlib;
static char rcsid[]="$Id: " __FILE__ ", " __DATE__ ", " __TIME__;

//some global definitions
static char g_srcRoot[MAXL]="";
static char g_destRoot[MAXL]="";
static char g_SerInstUID[MAXL]="";
static char g_MsgBuf[1024]="";
static unsigned short g_imWid=0;
static unsigned short g_imHt=0;

static int g_x0=0;
static int g_x1=0;
static int g_y0=0;
static int g_y1=0;
static int g_z0=0;
static int g_z1=0;
static bool g_bSort=true;
static bool g_bDcmRT=false;

//static float g_minThreshold=4.0;
//static float g_maxThreshold=1e+20;
//static float g_minOpacity=0.3;
//static float g_maxOpacity=0.8;
//static float g_PatternContrast=0.2;

static const int _MSG_ERROR=0;
static const int _MSG_WARNING=1;
static const int _MSG_STATUS=2;

void Msg(char* text, int type=_MSG_STATUS)
{
	switch (type)
	{
		case _MSG_ERROR: fprintf(stderr,"\nError: %s",text); break;
		case _MSG_WARNING: fprintf(stderr,"\nWarning: %s",text); break;
		case _MSG_STATUS: fprintf(stderr,"\n%s",text);break;
	}
}
void PrintUsage(ConsoleConfig& cc)
{
	fprintf (stderr, "\nDICOM to raw volume preserving original DICOM orientation. Optionally, extract subvolume");
	fprintf (stderr, "\nUsage:\t%s [bounds] [options] -D <DICOM_folder> -A <analyze root>\n",cc.m_program);

	fprintf (stderr, "\n\bounds (in pixels):\n");
	fprintf (stderr, "\t-[x0|y0|z0]<int>\tsubvolume lower x-y-z bounds");
	fprintf (stderr, "\n\t-[x1|y1|z1]<int>\tsubvolume upper x-y-z bounds");

	fprintf (stderr, "\noptions:\n");
	fprintf (stderr, "\t-u <SerInstUID>\t input Series Instance UID\n");
//	fprintf (stderr, "\t-s\t\t sort DICOM input by slice location (sorting is off by default)\n");
}

bool TestFunction()
{
	return true;
}
bool GetArgs(int argc, char* argv[])
{
		char* ptr, str[MAXL],c;
		if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; else ptr++;

		int k,i;
		bool bNextDCMFile=false, bNextAnFile=false, bNextSerInstUID=false, bEnd=false;
		int nreqs=0;

		int nSet=0;
		for (k = 0, i = 1; i < argc; i++) 
		{
			if (*(argv[i]) == '-')
			{
				bEnd=false;
				strcpy(str,argv[i]); 
				ptr=str;
				while(c=*(ptr++))
				{
					switch(c)
					{
						case 'D': bNextDCMFile=true;bEnd=true;		break;
						case 'A': bNextAnFile=true;bEnd=true;		break;
						case 'u': bNextSerInstUID=true;bEnd=true;	break;
						case 'x':
							c = *ptr++;
							if( c == '0') {g_x0=atoi(argv[i+1]); nSet++;}
							else {g_x1=atoi(argv[i+1]);nSet++;}
							break;
						case 'y':
							c = *ptr++;
							if( c == '0') {g_y0=atoi(argv[i+1]);nSet++;}
							else {g_y1=atoi(argv[i+1]);nSet++;}
							break;
						case 'z':
							c = *ptr++;
							if( c == '0') {g_z0=atoi(argv[i+1]);nSet++;}
							else {g_z1=atoi(argv[i+1]);nSet++;}
							break;
						case 's':
							g_bSort=true;
							break;
						case 'r':
							g_bDcmRT=true;
							break;
					}
					if(bEnd) break;
				}
			}
			else if(bNextDCMFile){bNextDCMFile=false; strcpy(g_srcRoot,argv[i]);nreqs++;}
			else if(bNextAnFile){bNextAnFile=false; strcpy(g_destRoot,argv[i]);nreqs++;}
			else if(bNextSerInstUID){bNextSerInstUID=false; strcpy(g_SerInstUID, argv[i]);}
		}
		if(nreqs<2) {return false;}
		if(nSet!=0 && nSet!=5) return false;
		return true;
}
void ReleaseMemory(ML3Array<DcmFileFormat*>& f)
{
	for(int i=0; i<f.GetSize(); i++) delete f[i];
}
int main (int argc, char *argv[])
{
	if (!TestFunction()) return 0;
	ConsoleConfig cc(rcsid);
	if(!cc.ProcessCommandLine(argc, argv) || !GetArgs(argc,argv))
	{
		PrintUsage(cc);
		return 0;
	}
	FileList fl;
	if(!fl.MakeList(g_srcRoot))
	{
		Msg("No files found in the specified DICOM folder",_MSG_ERROR);
		return 0;
	}

	fl.Sort();

	ML3Array<DcmFileFormat*> fformat;

	ReadDcmDir(fl, fformat, g_srcRoot,  g_SerInstUID,false);	

	if(fformat.GetSize()<1)
	{
		sprintf(g_MsgBuf,"No DICOM files found in %s\n", g_srcRoot);
		Msg(g_MsgBuf,_MSG_ERROR);
		ReleaseMemory(fformat);
		return 0;
	}
	else 
	{
		sprintf(g_MsgBuf,"%d DICOM files found\n", fformat.GetSize());
		Msg(g_MsgBuf,_MSG_STATUS);
	}

	char serUIDmsg[MAXL]="";
	if(strlen(g_SerInstUID)>0) strcpy(serUIDmsg,"with specified Series Instance UID ");
	
	Volume* v=0;
	int nVols=1;

	enum iFormat { OneFramePerFile, OneVolumePerFile };
	std::vector<int> file_formats;
	bool is_mixed_format, is_single_frame,is_format_undefined;
	DetectDcmFrameFormat(fformat,file_formats,is_mixed_format,is_single_frame,is_format_undefined);

	if (is_mixed_format || is_format_undefined || ! is_single_frame)
	{
		cout << "dcm2analyze ERROR: input number of frames per file is either inconsistent, unsupported or undefined." << endl;
		return -1;
	}
	if(!Dcm2Vol(fformat,v,nVols) || nVols<1)		
	{
		ReleaseMemory(fformat);
		return 0;
	}
	ReleaseMemory(fformat);
	Msg("Writing raw file...",_MSG_STATUS);
	if(nVols!=1) v[0].WriteAnalyze(g_destRoot);
	else {
		Msg("Multi-volume DICOM is not supported in this version.");
	}		
	delete [] v;
	fprintf(stderr, "\n");
	return 0;
}