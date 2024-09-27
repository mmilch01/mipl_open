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

bool InitSlice(Volume& v, DcmDataset* dataset, int z);
double GetSliceLocation(DcmFileFormat*& f);
int CompareInstanceNumbers(DcmFileFormat*& f1, DcmFileFormat*& f2);
int CompareSliceLocations(DcmFileFormat*& f1, DcmFileFormat*& f2);
void DetectDcmFrameFormat(ML3Array<DcmFileFormat*>& fformat, std::vector<int>& file_formats, \
	bool& is_mixed_format, bool& is_single_frame, bool& is_format_undefined);
bool Dcm2Vol(ML3Array<DcmFileFormat*>& ff, Volume*& v, int& nVols, bool bQuiet = false);
int GetInstanceNumber(DcmFileFormat*& f);
void GenerateUID(char* buf);
int CompareSliceLocations(DcmFileFormat*& f1, DcmFileFormat*& f2);
double ExtractFloat(DcmDataset* dataset, const DcmTagKey& key);

OFCondition LoadDcmFile(char* file, DcmFileFormat& format);


void ReadDcmDir(FileList& fl, ML3Array<DcmFileFormat*>& fformat, char* dir, char* serInstUID, bool bAll, bool bQuiet);
void ReadDcmDir(FileList& fl, ML3Array<DcmFileFormat*>& fformat, char* dir, char* serInstUID, bool bAll);

bool LooksLikeDICOMFile(char* filename);
void ReadDcmFile(FileList& fl, char* dir, DcmFileFormat& f);

bool SaveNewDcmFile(char* fname, DcmFileFormat* fformat, char* serInstUID, int frame, bool bMakeSC, vector< vector<std::string> >* modifyTags = NULL);
bool Is16bpp(DcmDataset* dataset);

unsigned short GetMaxIntensity(Volume& v, ML3Array<DcmFileFormat*>& fformat, bool bRobust = true);

}
