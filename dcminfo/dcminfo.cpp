#include <mlib3.h>
#include "../dcmlib/dcmlib.h"
//test change -- to be removed later.

#define _MAXTAGS 100
char g_DcmDir[FILENAME_MAX];

unsigned short g_Group[_MAXTAGS], g_Elem[_MAXTAGS];
int g_nTags=0;
char g_Mode=0;
bool g_bQuiet=false;
using namespace dcmlib;

bool GetArgs(int argc, char* argv[])
{
	AnyOption opt;

    /* 3. SET THE USAGE/HELP   */
    opt.addUsage( "Usage:");
	opt.addUsage( "dcminfo [-q] -t <group> <elem> [<group> <elem> ... <group> <elem>] <DICOM dir>" );
	opt.addUsage( "Options:");
	opt.addUsage( "\t\tprint DICOM tags\r-t" );
	opt.addUsage( "\t\tquiet mode\r-q" );
	
	opt.setFlag('q');
	opt.setFlag('t');

	opt.processCommandArgs(argc,argv,_MAXTAGS*2+2);
	g_bQuiet=opt.getFlag('q');

	if(opt.getFlag('t'))
	{	
		for (int i=0,j=0; i<opt.getArgc()-1; i+=2, j++)
		{
				std::stringstream ss1,ss2;
				ss1<<std::hex<<opt.getArgv(i);
				ss1>>g_Group[j];
				ss2<<std::hex<<opt.getArgv(i+1);
				ss2>>g_Elem[j];
				g_nTags++;
		}
		g_Mode='t';
		strcpy(g_DcmDir,opt.getArgv(opt.getArgc()-1));
		return true;
	}
	else
	{
		opt.printUsage(); 
		return false;
	}
}
void PrintOrientation(_3DINFO& info)
{
	ColumnVector c[3];
	for(int i=0; i<3; i++) c[i].ReSize(3);
	c[0] << info.ex1 << info.ex2 << info.ex3;
	c[1] << info.ey1 << info.ey2 << info.ey3;
	c[2] << info.ez1 << info.ez2 << info.ez3;
	for(int i=0; i<3; i++) c[i]/=c[i].NormFrobenius();
//	cout << "Direction cosines: " << endl << c[0] << endl << c[1] << endl << c[2] << endl;
	ColumnVector xx[3];
//	Matrix 
	for(int i=0; i<3; i++) xx[i].ReSize(3);
	xx[0] << 1 << 0 << 0;
	xx[1] << 0 << 1 << 0;
	xx[2] << 0 << 0 << 1;

	//for now, we assume axial acquisition.
	double cos0=DotProduct(c[0],xx[0]), cos1=DotProduct(c[1],xx[1]), cos2=DotProduct(c[2],xx[2]);
	if(fabs(cos2)<.7 || fabs(cos1)<.7)
	{
//		cout << "Non-axial acquisition, not supported."<<endl;
//		return;
	}	
	//Cosines in RAI (DICOM) system.
	bool xflip=(DotProduct(c[0],xx[0])<0), yflip=(DotProduct(c[1],xx[1])<0), zflip=(DotProduct(c[2],xx[2])<0);
	cout << (xflip?"L":"R") << (yflip?"P":"A") << (zflip?"S":"I");
}

int main (int argc, char *argv[])
{
	static char rcsid[] = "$Id: " __FILE__ ", " __DATE__ ", " __TIME__;
	if(!GetArgs(argc, argv))
		return 0;
	if(!g_bQuiet) std::cout << rcsid << endl;
//	cout << "Print DICOM tag(s)" <<", ver.: "<< __DATE__ << " " << __TIME__ << endl;
	//first, read in the original file.
	FileList fl;
	if(!fl.MakeList(g_DcmDir))
	{
//		cout << "No files found in input directory" << endl;
		return 0;
	}
	if(g_Mode=='t')
	{
		DcmFileFormat dff;
		ReadDcmFile(fl, g_DcmDir, dff);
		DcmDataset *dd=dff.getDataset();
		OFCondition ofc;
		DcmElement* el;
		OFString str;
		unsigned long pos=0;
		for(int i=0; i<g_nTags; i++)
		{
			ofc=dd->findAndGetElement(DcmTagKey(g_Group[i],g_Elem[i]),el);
			if(ofc.good())	el->getOFStringArray(str,pos);
			else str="<empty>";
			std::cout << str.data() << " ";
		}
		std::cout << endl;
		return 0;
	}
}