#include "dcmlib.h"

namespace dcmlib{
//structure that contains 3D coordinate/slice transformation info.


OFCondition LoadDcmFile(char* file, DcmFileFormat& format)
{
	return format.loadFile(file);
}
void GenerateUID(char* buf)
{
	static int ind = 1;
	::dcmGenerateUniqueIdentifier(buf);
	int len = strlen(buf);
	sprintf(buf + len, "%d", ind++);
	if (ind > 9999) ind = 1;
}
unsigned short GetMaxIntensity(Volume& v, ML3Array<DcmFileFormat*>& fformat, bool bRobust/*=true*/)
{
	DcmDataset* dataset;
	unsigned short res = 0;
	unsigned char* buf;
	unsigned char*& rbuf = buf;
	Analyze an;
	an.m_dsr.dime.dim[0] = 4;
	an.m_dsr.dime.dim[1] = v.m_dims[0];
	an.m_dsr.dime.dim[2] = v.m_dims[1];
	an.m_dsr.dime.dim[3] = 1;
	an.m_dsr.dime.dim[4] = 1;
	an.m_dsr.dime.datatype = (short)Analyze::DT_SIGNED_SHORT;
	an.m_dsr.dime.bitpix = 16;

	Volume v0;
	v0 = v;
	for (int i = 0; i < fformat.GetSize(); i++)
	{
		dataset = fformat[i]->getDataset();
		dataset->findAndGetUint8Array(DCM_PixelData, (const Uint8*&)rbuf);
		v0.InitMemory(an, buf, i);
	}
	if (bRobust) {
		//99th percentile
		Matrix hist; 
		v0.Hist(hist);
		int ind = 1, last_ind = 1;
		double val, sum = 0, msum = hist.Column(1).Sum() * 0.99;
		while (sum < (long)msum)
		{
			val = hist(ind, 2);
			sum += hist(ind, 1);
			if (val > 0) last_ind = ind;
			ind++;
		}
		return hist(last_ind, 2);
	}
	else
		return v0.GetMax();
}
int CompareSliceLocations(DcmFileFormat*& f1, DcmFileFormat*& f2)
{
	DcmDataset* dd1 = f1->getDataset(), * dd2 = f2->getDataset();
	if (!dd1 || !dd2) return 0;
	OFString ofs;
	double num1, num2;
	OFCondition ofc;
	ofc = dd1->findAndGetOFString(DCM_SliceLocation, ofs);
	if (ofc.bad()) return 0;
	num1 = atof(ofs.data());
	ofc = dd2->findAndGetOFString(DCM_SliceLocation, ofs);
	if (ofc.bad()) return 0;
	num2 = atof(ofs.data());
	if (num1 < num2) return -1;
	else return (num1 == num2) ? 0 : 1;
}

int GetInstanceNumber(DcmFileFormat*& f)
{
	DcmDataset* dd = f->getDataset();
	if (!dd) return 0;
	OFString ofs;
	OFCondition ofc = dd->findAndGetOFString(DCM_InstanceNumber, ofs);
	if (ofc.bad()) return 0;
	return atoi(ofs.data());
}


double ExtractFloat(DcmDataset* dataset, const DcmTagKey& key)
{
	OFString ofs;
	OFCondition ofc = dataset->findAndGetOFString(key, ofs);
	if (ofc.bad()) return 0.0;
	return atof(ofs.data());
}

bool LooksLikeDICOMFile(char *filename)
{
	if(!filename)	 return false;

	// Read first 132 bytes
	FILE* fp = fopen(filename,"rb");
	if(!fp)	return false;
	unsigned char	top[132];
	int nread = fread(top,1,132,fp);
	fclose(fp);

	// Check 132 bytes for "DICOMness"
	if(nread<132)	return false;
	else			
	{
		if(top[128]=='D' && top[129]=='I' && 
				top[130]=='C' && top[131]=='M')	return true;
		else return false;
	}
}

void ReadDcmFile(FileList& fl, char* dir, DcmFileFormat& dff)
{
	ML3Array<DcmFileFormat*> fformat;
	OFCondition error;
	char fpath[MAXL];
	bool bFound = false;
	char* nextFile = fl.GetFirstFile();

	for (int i = 0; i < fl.NFiles(); i++)
	{
#ifdef _WIN32
		sprintf(fpath, "%s\\%s", dir, nextFile);
#else
		sprintf(fpath, "%s/%s", dir, nextFile);
#endif
		if (!LooksLikeDICOMFile(fpath)) continue;
		error = dff.loadFile(fpath);
		if (!error.bad())
		{
			bFound = true; break;
		}
	}
	if (!bFound)
	{
		//cout << "No DICOM files found" << endl;
		return;
	}
}

void ReadDcmDir(FileList& fl, ML3Array<DcmFileFormat*>& fformat, char* dir, char* serInstUID, bool bAll)
{
	ReadDcmDir(fl, fformat, dir, serInstUID, bAll, false);
};

void ReadDcmDir(FileList& fl, ML3Array<DcmFileFormat*>& fformat, char* dir, char* serInstUID, bool bAllSeries, bool bQuiet)
{
	int nDCMFiles = 0;
	OFCondition error;
	char fpath[MAXL];
	char rserUID[MAXL];
	OFString str;
	std::string* nextFile = fl.GetFirstFile_str();
	double echo_time = 0;
	if (!bQuiet)
		fprintf(stderr, "\nReading DICOM dir...\n");
	char buf[256] = "";
	if (serInstUID == NULL)
		serInstUID = buf;

	for (int i = 0; i < fl.NFiles(); i++)
	{
		if (!bQuiet) fprintf(stderr, "\rReading file %d out of %d", i + 1, fl.NFiles());
		if (i > 0) nextFile = fl.GetNextFile_str();
#ifdef _WIN32
		sprintf(fpath, "%s\\%s", dir, nextFile->c_str());
#else
		sprintf(fpath, "%s/%s", dir, nextFile->c_str());
#endif
		//		if(!LooksLikeDICOMFile(fpath)) continue;
		DcmFileFormat* f = new DcmFileFormat();
		error = f->loadFile(fpath);
		if (error.bad()) { delete f; continue; }
		//define series inst uid from the first valid DICOM file.
		OFString str;
		if (fformat.GetSize() < 1 && strlen(serInstUID) == 0)
		{
			f->getDataset()->findAndGetOFString(DCM_SeriesInstanceUID, str);
			strcpy(serInstUID, str.data());
			//			findAndGetString(DCM_SeriesInstanceUID,(const char *&) g_SerInstUID);
		}
		else if (!bAllSeries)
		{
			//check if series inst uid's match
			f->getDataset()->findAndGetOFString(DCM_SeriesInstanceUID, str);
			strcpy(rserUID, str.data());
			//this file belongs to a different series
			if (strcmp(rserUID, serInstUID)) { delete f; continue; }
			//for MR scans - check echo times.
//			if(ExtractFloat(f->getDataset(),DCM_EchoTime)!=echo_time) {delete f; continue;}
		}
		fformat.Push(f);
		if (fformat.GetSize() == 1) //for the first image, check echo time.
		{
			echo_time = ExtractFloat(f->getDataset(), DCM_EchoTime);
		}
	}
	if (!bQuiet) fprintf(stderr, "\n");
}

bool Is16bpp(DcmDataset* dataset)
{
	//verify that byte format is supported
	unsigned short par1 = 0, par2 = 0;
	dataset->findAndGetUint16(DCM_BitsAllocated, par1);
	dataset->findAndGetUint16(DCM_SamplesPerPixel, par2);
	return (par1 == 16 || par2 == 1);
}

//save (anonymized) DICOM file.
bool SaveNewDcmFile(char* fname, DcmFileFormat* fformat, char* serInstUID, int frame, bool bMakeSC, vector<vector<std::string> >* modifyTags)
{
	DcmDataset* dataset = fformat->getDataset();

	long nFramesPerFile = 1;
	dataset->findAndGetLongInt(DCM_NumberOfFrames, nFramesPerFile);
	bool bMultiFrame = (nFramesPerFile > 1);

	//	update DICOM tags
	OFCondition ofc;
	if (bMakeSC) DcmCodec::convertToSecondaryCapture(dataset);
	if (serInstUID) dataset->putAndInsertString(DCM_SeriesInstanceUID, serInstUID);
	char str[MAXL];
	GenerateUID(str);
	ofc = dataset->putAndInsertString(DCM_SOPInstanceUID, (const char*)str);
	if (ofc.bad()) return false;
	unsigned short g, e;
	if (modifyTags)
	{
		for (unsigned i = 0; i < modifyTags->size(); i++)
		{
			g = strtol(((*modifyTags)[i][0].c_str()), NULL, 16);
			e = strtol(((*modifyTags)[i][1].c_str()), NULL, 16);
			ofc = dataset->putAndInsertString(DcmTag(g, e), (*modifyTags)[i][2].c_str());
			if (ofc.bad()) return false;
		}
	}
	char lfname[1024];
#ifdef _WIN32
	sprintf(lfname, "%s\\%03d_%s.dcm", fname, frame, str);
#else 
	sprintf(lfname, "%s/%03d_%s.dcm", fname, frame, str);
#endif
	sprintf(str, "%d", frame + 1);
	//	ofc=dataset->putAndInsertString(DCM_InstanceNumber,str);
	if (ofc.bad()) return false;

	ofc = fformat->saveFile(lfname);
	return (!ofc.bad());
}

}//end of namespace dcmlib.
