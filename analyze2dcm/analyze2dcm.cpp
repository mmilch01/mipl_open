//#ifndef _4DFP 
//#define _4DFP
//#endif
#ifndef _WIN32
#if !defined(HAVE_CONFIG_H)
#define HAVE_CONFIG_H
#endif
#endif

#include <mlib3.h>
#include <vector>
#include "../dcmlib/dcmlib.h"
#include "zlib.h"
#include "zconf.h"

#ifndef _LODEPNG_INCLUDED
#include "lodepng.h"
#define _LODEPNG_INCLUDED
#endif

static char rcsid[]="$Id: " __FILE__ ", " __DATE__ ", " __TIME__;

//some global definitions
static char g_DCMSrc[MAXL]="";
static char g_anSrc[MAXL]="";
static char g_OutRoot[MAXL]="";
static char g_PN[MAXL]="";
static char g_SerDescr[MAXL]="";
static char g_SerInstUID[MAXL]="";
static char g_SerNum[MAXL]="";
static char g_MsgBuf[1024]="";
static bool g_xFlip=false;
static bool g_yFlip=false;
static bool g_zFlip=false;
static int g_iPNG=0;
static bool g_bSC=false;
static int g_Max=-1;
using namespace dcmlib;

//static float g_minThreshold=4.0;
//static float g_maxThreshold=1e+20;
//static float g_minOpacity=0.3;
//static float g_maxOpacity=0.8;
//static float g_PatternContrast=0.2;
static int	 g_Origin[3]={0,0,0};

static const int _MSG_ERROR=0;
static const int _MSG_WARNING=1;
static const int _MSG_STATUS=2;
static const int _PNG_MODE_RGB=1;
static const int _PNG_MODE_MONO=2;

void Msg(char* text, int type=_MSG_STATUS)
{
	switch (type)
	{
		case _MSG_ERROR: fprintf(stderr,"\nError: %s",text); break;
		case _MSG_WARNING: fprintf(stderr,"\nWarning: %s",text); break;
		case _MSG_STATUS: fprintf(stderr,"\n%s",text);break;
	}
}

bool TestFunction()
{
	return true;
}
bool GetArgs(int argc, char* argv[])
{
    /* 3. SET THE USAGE/HELP   */
	AnyOption opt;
	opt.addUsage( "image -> DICOM" );
    opt.addUsage( "Usage: analyze2dcm [options] <DICOM_dir> <image root>");
	opt.addUsage( "Options:");
//	opt.addUsage( " \t\t\t-c <x,y,z>\t origin coordinates" );
	opt.addUsage( " \t\t-o <dir>\t output root" );
	opt.addUsage( " \t\t-u <SerInstUID>\t DICOM Series Instance UID" );
	opt.addUsage( " \t\t-p <PatName>\t set patient name" );
	opt.addUsage( " \t\t-d <Description> set series description" );
	opt.addUsage( " \t\t-n <int>\t set series number" );
	opt.addUsage( " \t\t-s <int>\t scale to maximum value" );
	opt.addUsage( " \t\t-g <c|m> \t assume PNG as input format, write (c)olor or (m)onochrome DICOM");
	opt.addUsage( " \t\t-c \t\t convert to SC [preserve original modality]");
	opt.addUsage( " \t\t-x \t\t flip x axis" );
	opt.addUsage( " \t\t-y \t\t flip y axis" );
	opt.addUsage( " \t\t-z \t\t flip z axis" );

//	opt.setOption('c');
	opt.setOption('o');
	opt.setOption('u');
	opt.setOption('s');
	opt.setOption('n');
	opt.setOption('p');
	opt.setOption('d');
	opt.setFlag('x');
	opt.setFlag('y');
	opt.setFlag('z');
	opt.setOption('g');
	opt.setFlag('c');

	opt.processCommandArgs(argc,argv);

	if(opt.getArgc()<2)// || (opt.getValue('p'))/*|| opt.getValue('f'))*/)
	{
		opt.printUsage();
		return false;
	}

	strcpy(g_anSrc,opt.getArgv(1));
	strcpy(g_DCMSrc,opt.getArgv(0));
	if(opt.getFlag('x')) g_xFlip=true;
	if(opt.getFlag('y')) g_yFlip=true;
	if(opt.getFlag('z')) g_zFlip=true;
	if(opt.getFlag('c')) g_bSC=true;

	char *s;
	if(s=opt.getValue('g'))
	{ 
		if (strcmp(s,"c")==0) 
			g_iPNG=_PNG_MODE_RGB; 	
		else if (strcmp(s,"m")==0)
			g_iPNG=_PNG_MODE_MONO;
		else
		{ opt.printUsage(); return false;}
	}
	if(s=opt.getValue('o')) strcpy(g_OutRoot,s);
	if(s=opt.getValue('u')) strcpy(g_SerInstUID,s);
	if(s=opt.getValue('p')) strcpy(g_PN,s);
	if(s=opt.getValue('d')) strcpy(g_SerDescr,s);
	if(s=opt.getValue('n')) strcpy(g_SerNum,s);
	if(s=opt.getValue('s')) g_Max=atoi(s);
	return true;
}
void RGB2Real(std::vector<unsigned char>& image, Matrix& slice, int w, int h)
{
	slice.ReSize(w,h);
	Real val,m=1.0/3.0;
	int ind;
	for(int j=0; j<w; j++)
	{
		ind=j*h;
		for(int i=0; i<h; i++)
		{
			val=(Real)(image[ind+i]);
			slice(j+1,i+1)=val;
		}
	}
}
void RGB2Real(std::vector<unsigned char>& image, Matrix& slice, int w, int h, int channel)
{
	//channel=0(R),1(G),2(B),3(alpha).
	slice.ReSize(w,h);
	Real val,m=1.0/3.0;
	int ind;//,ind1;
	for(int j=0; j<w; j++)
	{
		ind=j*h*4+channel;
		for(int i=0,ii=0; i<h; i++,ii+=4)
		{
			val=(Real)(image[ind+ii]);
			slice(j+1,i+1)=val;
		}
	}
}
bool ReadPNG_mono(FileList& fl, Volume& v)
{
	char* fil=fl.GetFirstFile();	
	std::vector<unsigned char> image;
	unsigned w, h;
	unsigned error;
	bool first=true;
	Matrix slice;
	int ind=0;
	std::string path;
	do
	{
		path=g_anSrc;
		path.append("/"); path.append(fil);
		image.clear();
		error=lodepng::decode(image,w,h,path,LCT_GREY,8);
		if(error) return false;
		if (first) v.InitMemory(w,h,fl.NFiles());
		RGB2Real(image,slice,w,h);
		v.SetSlice(slice,ind);
		first=false;
		ind++;
	} while(fil=fl.GetNextFile());
//	v.Flip(2);
	return true;
}

bool ReadPNG_rgb(FileList& fl, Volume* v)
{
	char* fil=fl.GetFirstFile();	
	std::vector<unsigned char> image;
	unsigned w, h;
	unsigned error;
	bool first=true;
	Matrix slice;
	int ind=0;
	std::string path;
	do
	{
		path=g_anSrc;
		path.append("/"); path.append(fil);
		image.clear();
        //cout << path << endl; 
		error=lodepng::decode(image,w,h,path);//,LCT_GREY,8);
		if(error) {
            cout << "decoding error " << error << endl;
            return false;
        }
		if (first) for (int i=0; i<3; i++) v[i].InitMemory(w,h,fl.NFiles());
		for (int i=0; i<3; i++) 
		{
			RGB2Real(image,slice,w,h,i);
			v[i].SetSlice(slice,ind);
		}
		first=false;
		ind++;
		fprintf(stderr,"\rReading PNG files: %d out of %d", ind, fl.NFiles());
	} while(fil=fl.GetNextFile());
//	v.Flip(2);
	return true;
}

// save single-frame RGB DICOM
bool SaveDCM_rgb(char* fname, DcmFileFormat* fformat, char* serInstUID, Volume* v, int frame)
{
	DcmDataset *dataset=fformat->getDataset();
	//verify that byte format is supported
	unsigned short par1=0,par2=0,dcmw=0,dcmh=0;

//	dataset->findAndGetUint16(DCM_BitsAllocated,par1);
//	dataset->findAndGetUint16(DCM_SamplesPerPixel,par2);

	dataset->findAndGetUint16(DCM_Rows,dcmh);
	dataset->findAndGetUint16(DCM_Columns,dcmw);
	long nFramesPerFile=1;
//	dataset->findAndGetLongInt(DCM_NumberOfFrames,nFramesPerFile);
	
	int start_frame=g_Origin[2],
		end_frame=start_frame+v[0].m_dims[2]-1;

	//modify voxels if this frame has non-empty intersection with subvolume
	if(frame>=start_frame && frame<=end_frame)
	{
//		unsigned char*& rbuf=buf;
		unsigned long count;
		unsigned short w=v[0].SX(),h=v[0].SY();

		unsigned char *buf=new unsigned char[3*v[0].SX()*v[0].SY()];
		int za=frame-g_Origin[2];
		int ind=0;
		for(int y=g_Origin[1],ya=0,yy=g_Origin[1]*3; y<g_Origin[1]+v[0].m_dims[1]; y++,ya++,yy+=3)
		{
			ind=w*yy;
			for(int x=g_Origin[0],xa=0,xx=g_Origin[0]*3; x<g_Origin[0]+v[0].m_dims[0]; x++,xa++,xx+=3)
			{
				buf[ind+xx]=(unsigned char)(v[0](xa,ya,za)+.5);
				buf[ind+xx+1]=(unsigned char)(v[1](xa,ya,za)+.5);
				buf[ind+xx+2]=(unsigned char)(v[2](xa,ya,za)+.5);
			}
		}
		dataset->putAndInsertUint8Array(DCM_PixelData, (Uint8*)buf,v[0].SX()*v[0].SY()*3);
		dataset->putAndInsertUint16(DCM_Rows,v[0].SY());
		dataset->putAndInsertUint16(DCM_Columns,v[0].SX());
		char s[20];
		sprintf(s,"%d",frame+1);
		dataset->putAndInsertOFStringArray(DCM_InstanceNumber,s);
		dataset->putAndInsertOFStringArray(DCM_AcquisitionNumber,s);
		delete[] buf;
	}
//	update DICOM tags
	if (g_bSC) DcmCodec::convertToSecondaryCapture(dataset);

	dataset->putAndInsertString(DCM_SeriesInstanceUID, serInstUID);
	dataset->putAndInsertUint16(DCM_NumberOfFrames,1);
	dataset->putAndInsertUint16(DCM_SamplesPerPixel,3);
	dataset->putAndInsertUint16(DCM_PlanarConfiguration,0);
	dataset->putAndInsertString(DCM_PhotometricInterpretation,"RGB");
	dataset->putAndInsertUint16(DCM_BitsAllocated,8);
	dataset->putAndInsertUint16(DCM_BitsStored,8);
	dataset->putAndInsertUint16(DCM_HighBit,7);
	dataset->putAndInsertUint16(DCM_PixelRepresentation,0);
	dataset->putAndInsertString(DCM_NumberOfFrames,"1");

	if(strlen(g_SerDescr)>0)
		dataset->putAndInsertString(DCM_SeriesDescription, g_SerDescr);
	if(strlen(g_PN)>0)
		dataset->putAndInsertString(DcmTagKey(0x0010, 0x0010),g_PN);
	if(strlen(g_SerNum)>0)
		dataset->putAndInsertString(DCM_SeriesNumber,g_SerNum);
	return SaveNewDcmFile(fname, fformat, serInstUID, frame, false);
}
bool SaveDCM(char* fname, DcmFileFormat* fformat, char* serInstUID, Volume& v, int frame, double slope, double intercept, bool bResize)
{
	DcmDataset *dataset=fformat->getDataset();
	//verify that byte format is supported
	unsigned short par1=0,par2=0,dcmw=0,dcmh=0;
	dataset->findAndGetUint16(DCM_Rows,dcmh);
	dataset->findAndGetUint16(DCM_Columns,dcmw);
	long nFramesPerFile=1;
	dataset->findAndGetLongInt(DCM_NumberOfFrames,nFramesPerFile);
	
	if(!Is16bpp(dataset))
	{
		Msg("Unsupported pixel format",_MSG_ERROR);
		return false;
	}
	int start_frame=g_Origin[2],
		end_frame=start_frame+v.m_dims[2]-1;

	//modify voxels if this frame has non-empty intersection with subvolume
	if(frame>=start_frame && frame<=end_frame)
	{
		unsigned char *buf;
		unsigned char*& rbuf=buf;
		unsigned long count;
		unsigned short w=v.SX(),h=v.SY();

		//check if dimensions of fformat and v match.
		if( bResize ) //pixel dimensions & number of frames will be resized
		{
			unsigned short *buf=new unsigned short[v.SX()*v.SY()];
			int za=frame-g_Origin[2];
			for(int y=g_Origin[1],ya=0; y<g_Origin[1]+v.m_dims[1]; y++,ya++)
			{
				for(int x=g_Origin[0],xa=0; x<g_Origin[0]+v.m_dims[0]; x++,xa++)
				{
					buf[w*y+x]=(unsigned short)(v(xa,ya,za)*slope+intercept+.5);
				}
			}
			dataset->putAndInsertUint8Array(DCM_PixelData, (Uint8*)buf,v.SX()*v.SY()*2);
			dataset->putAndInsertUint16(DCM_Rows,v.SY());
			dataset->putAndInsertUint16(DCM_Columns,v.SX());
			char s[20];
			sprintf(s,"%d",frame+1);
			dataset->putAndInsertOFStringArray(DCM_InstanceNumber,s);
			dataset->putAndInsertOFStringArray(DCM_AcquisitionNumber,s);
			delete[] buf;
		}
		else //pixel dimensions and buffer stay the same.
		{
			dataset->findAndGetUint8Array(DCM_PixelData,(const Uint8* &)rbuf, &count);
			unsigned short *ptr=(unsigned short*)buf;
		//	int za=v.m_dims[2]-(frame-g_Origin[2])-1;
			{
				if ( g_Origin[0]+v.m_dims[0]>dcmw || g_Origin[1]+v.m_dims[1]>dcmh )
				{
					cout << "Source and destination dimensions don't match, exiting";
					return false;
				}
				int za=frame-g_Origin[2];
				for(int y=g_Origin[1],ya=0; y<g_Origin[1]+v.m_dims[1]; y++,ya++)
				{
					for(int x=g_Origin[0],xa=0; x<g_Origin[0]+v.m_dims[0]; x++,xa++)
					{
						ptr[dcmw*y+x]=(unsigned short)(v(xa,ya,za)*slope+intercept+.5);
					}
				}
			}
		}
	}

//	update DICOM tags
	if(g_bSC) DcmCodec::convertToSecondaryCapture(dataset);

	OFCondition ofc=dataset->putAndInsertString(DCM_SeriesInstanceUID, serInstUID);
/*	char imgSOPInstUID[MAXL];
	GenerateUID(imgSOPInstUID);
	ofc=dataset->putAndInsertString(DCM_SOPInstanceUID, (const char*)imgSOPInstUID);	
	if(ofc.bad())
	{
		sprintf(g_MsgBuf,"Update of (0x0008, 0x0018) (value =\"%s\") tag failed. Code = %d, Module = %d, Message = \"%s\"",
			imgSOPInstUID,ofc.code(),ofc.module(),ofc.text());
		Msg(g_MsgBuf, _MSG_WARNING);
	}
*/
	if(strlen(g_SerDescr)>0)
		dataset->putAndInsertString(DCM_SeriesDescription, g_SerDescr);
	if(strlen(g_PN)>0)
		dataset->putAndInsertString(DcmTagKey(0x0010, 0x0010),g_PN);
	if(strlen(g_SerNum)>0)
		ofc=dataset->putAndInsertString(DCM_SeriesNumber,g_SerNum);
	return SaveNewDcmFile(fname, fformat, serInstUID, frame, false);
}
void ReleaseMemory(ML3Array<DcmFileFormat*>& f)
{
	for(int i=0; i<f.GetSize(); i++) delete f[i];
}

int main (int argc, char *argv[])
{
	if (!TestFunction()) return 0;
	ConsoleConfig cc(rcsid);

	if(!GetArgs(argc,argv)) return 0;

	cc.ProcessCommandLine(argc,argv);
	if(strlen(g_OutRoot)<1)
	{
		cc.GetRoot(g_anSrc,g_OutRoot);
		strcat(g_OutRoot,".dcm");
	}
	FileList fl;
	if(!fl.MakeList(g_DCMSrc))
	{
		cout << "No files found in specified DICOM folder";
		return -1;
	}
	
	ML3Array<DcmFileFormat*> fformat;
	ReadDcmDir(fl, fformat, g_DCMSrc,  g_SerInstUID,false);
	char serUIDmsg[MAXL]="";
	if(strlen(g_SerInstUID)>0) strcpy(serUIDmsg,"with specified Series Instance UID ");
	std::vector<ML3Array<DcmFileFormat*> > vols;
	long nFrames=0;
	int nVols;
	nVols=1;
	
	if(nVols<1) 
	{
		cout << "No valid volumes found in DICOM directory" << endl;
		return -1;
	}
	Volume* v=NULL;
	if (g_iPNG==_PNG_MODE_RGB) //read in PNG dir.
	{
		FileList fl1;
		nVols=3;
		v=new Volume[3];
		if(!fl1.MakeList(g_anSrc))
		{
			cout << "No files found in specified PNG folder";
			return -1;
		}
		cout << "Reading PNG files in RGB space..." << endl;
		if(!ReadPNG_rgb(fl1, v))
		{
			cout << "Cannot read a set of PNG files from directory " << g_anSrc << endl;
			return -1;
		}
		fl1.Sort();
	}
	else if (g_iPNG==_PNG_MODE_MONO)
	{
		FileList fl1;
		nVols=1;
		v=new Volume[1];
		if(!fl1.MakeList(g_anSrc))
		{
			cout << "No files found in specified PNG folder";
			return -1;
		}
		cout << "Reading PNG files in monochrome space..." << endl;
		if(!ReadPNG_mono(fl1, *v))
		{
			cout << "Cannot read a set of PNG files from directory " << g_anSrc << endl;
			return -1;
		}
		fl1.Sort();
	}
	else // read in Analyze volume.
	{
		v=new Volume[1];
		if(!cc.ReadVolume(v[0],g_anSrc))
		{
			cout << "Cannot read Analyze volume" << endl;
			delete[] v;
			return -1;
		}
		for (int i=0; i<nVols; i++)
			vols[i].Sort(CompareSliceLocations);
	}
//	fformat.Sort(CompareSliceLocations);

	Msg("Writing DICOM...",_MSG_STATUS);
	char dcmfile[MAXL], command[MAXL], serInstUID[MAXL];
	sprintf(command, "mkdir %s",g_OutRoot);
	system(command);
	GenerateUID(serInstUID);

//	for(int z=0; z<v.m_dims[2]; z++)
	fprintf(stderr,"\n");

	for(int i=0; i<nVols; i++)
	{
		if(g_xFlip) v[i].Flip(1);
		if(g_yFlip) v[i].Flip(2);
		if(g_zFlip) v[i].Flip(3);
	}
	Matrix stats;
	double a,b,mx,mn,diff;
	int inst;
	if (g_iPNG==_PNG_MODE_RGB)
	{
		for(int z=0; z<v[0].SZ(); z++)
		{
				fprintf(stderr,"\rWriting file %d out of %d", z+1, v[0].SZ());
				if(!SaveDCM_rgb(g_OutRoot,fformat[0],serInstUID,v,z))
				{
					cout << "Error writing " << dcmfile << endl;
					delete [] v;
					return -1;
				}
		}
	} //end if color PNG
	else //if not color PNG
	{
		for (int i=0; i<nVols; i++)
		{
			v[i].GetMaxMin(mx,mn);
			diff=mx-mn;
			if(diff==0) diff=1.0;
			a=g_Max/diff; b=-mn*g_Max/diff;
			if(g_Max<=0) {a=1; b=0;}
			cout << "Writing volume " << i+1 << " of " << nVols << endl;
			//multi-frame
			if (nFrames>1)
			{
				fprintf(stderr,"\rWriting multi-frame volume, file 1 of 1");
				if(!SaveDCM(g_OutRoot, fformat[0], serInstUID, v[i], 1, a,b,false))
				{
					cout << "Error writing " << dcmfile << endl;
					Msg(g_MsgBuf,_MSG_ERROR);
					ReleaseMemory(fformat);
					delete [] v;
					return -1;
				}			
			}
			else
			{
				for(int z=0; z<v[i].SZ(); z++)
				{
					if(g_iPNG==_PNG_MODE_MONO)
					{
						fprintf(stderr,"\rWriting file %d out of %d", z+1, v[i].SZ());
						if(!SaveDCM(g_OutRoot,fformat[0],serInstUID,v[i],z,a,b,true))
						{
							cout << "Error writing " << dcmfile << endl;
							delete [] v;
							return -1;
						}
					}
					else
					{
						inst=GetInstanceNumber(vols[i][z]);
			//			cout << "Writing instance " << inst << " out of total " << nVols*v[i].SZ() << ;
						fprintf(stderr,"\rWriting file %d out of %d, instance %d", z+1, v[i].SZ(),inst);
						if(!SaveDCM(g_OutRoot, vols[i][z], serInstUID, v[i], z,a,b,false))
						{
							cout << "Error writing " << dcmfile << endl;
							Msg(g_MsgBuf,_MSG_ERROR);
							ReleaseMemory(fformat);
							delete [] v;
							return -1;
						}
					}
				}
			}
			cout << endl;
		}
	} //if not PNG
	ReleaseMemory(fformat);
	cout << endl;
	delete [] v;
	return 0;
}
