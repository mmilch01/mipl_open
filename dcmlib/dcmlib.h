#ifndef _WIN32
#if !defined(HAVE_CONFIG_H)
#define HAVE_CONFIG_H
#endif
#endif
#include <vector>
#include <algorithm>
#include <iostream>

#include <dcmtk/config/osconfig.h>
//#ifndef HAVE_STRSTREAM
//#error "cfunix.h not included"
//#endif
#include <dcmtk/dcmdata/dctk.h>
#include <dcmtk/dcmdata/dccodec.h>
#include <dcmtk/ofstd/ofcond.h>
#include <dcmtk/ofstd/ofstdinc.h>

//#include <iostream>
//#ifdef WIN32
//#endif
#include "../mlib3/mlib3.h"
namespace dcmlib{

static const int COMPARE_INSTANCE=0;
static const int COMPARE_SLICE_LOCATION=1;

typedef struct 
{
	//width/height in image coordinates
	unsigned short w,h,d;
	//pixel spacing (voxel size in mm)
	double dx,dy,dz;
	//direction cosines
	double ex1,ex2,ex3,ey1,ey2,ey3,ez1,ez2,ez3;
	//image position (DICOM) for the middle slice;
	double xm,ym,zm;
	//volume center (in patient coordinates)	
	//double xc,yc,zc;
	//slice location for the 1st slice in dataset.
	double slx,sly,slz;
} _3DINFO;
/*
bool InitSlice(Volume& v,DcmDataset *dataset, int z);
bool SaveVolFrameDCM(const char* fname, DcmFileFormat* fformat, Volume& v, int frame);
bool ReadDcmImage(std::string& fname, DcmFileFormat& dff, Volume &v);
int	ExtractFrames(ML3Array<DcmFileFormat*>& ff, std::vector<ML3Array<DcmFileFormat*> >& frames);
int CompareSliceLocations(DcmFileFormat*& f1, DcmFileFormat*& f2);
int CompareInstanceNumbers(DcmFileFormat*& f1, DcmFileFormat*& f2);
int CompareEchos(DcmFileFormat*& f1, DcmFileFormat*& f2);
bool PreceedsSliceLocation(DcmFileFormat*& f1, DcmFileFormat*& f2);
bool PreceedsInstanceNumber(DcmFileFormat*& f1, DcmFileFormat*& f2);
bool PreceedsEchoNumber(DcmFileFormat*& f1, DcmFileFormat*& f2);
int find_vol(DcmFileFormat* dff, std::vector<int>& stind);

double GetSliceLocation(DcmFileFormat*& f);
int	GetEchoNumber(DcmFileFormat*& f);
bool ExtractFloatArray(DcmDataset* dataset, const DcmTagKey& key, double* arr, int sz);
//
int GetInstanceNumber(DcmFileFormat*& f);

bool Is16bpp(DcmDataset* dataset);
bool IsMultiFrame(ML3Array<DcmFileFormat*>& ff, long &nFrames);
bool SaveNewDcmFile(char* fname, DcmFileFormat* fformat, char* serInstUID, int frame, bool bMakeSC,vector< vector<std::string> >* modifyTags=NULL);
void GenerateUID(char* buf);
*/

int GetInstanceNumber(DcmFileFormat*& f);
void GenerateUID(char* buf);
int CompareSliceLocations(DcmFileFormat*& f1, DcmFileFormat*& f2);
double ExtractFloat(DcmDataset* dataset, const DcmTagKey& key);

OFCondition LoadDcmFile(char* file, DcmFileFormat& format);

//bool InitOrientationData(ML3Array<DcmFileFormat*>& dcmFiles, Volume& v, _3DINFO& S0, _3DINFO& S1, bool bRescale);
void ReadDcmDir(FileList& fl, ML3Array<DcmFileFormat*>& fformat, char* dir, char* serInstUID, bool bAll, bool bQuiet);
void ReadDcmDir(FileList& fl, ML3Array<DcmFileFormat*>& fformat, char* dir, char* serInstUID, bool bAll);

bool LooksLikeDICOMFile(char* filename);
void ReadDcmFile(FileList& fl, char* dir, DcmFileFormat& f);

bool SaveNewDcmFile(char* fname, DcmFileFormat* fformat, char* serInstUID, int frame, bool bMakeSC, vector< vector<std::string> >* modifyTags = NULL);
bool Is16bpp(DcmDataset* dataset);

unsigned short GetMaxIntensity(Volume& v, ML3Array<DcmFileFormat*>& fformat, bool bRobust = true);

}
