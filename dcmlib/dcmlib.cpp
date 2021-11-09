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

/*
bool InitVolMultiFrame(Volume& v, DcmDataset *dataset, int w, int h, long nFrames)
{
	v.InitMemory(w,h,nFrames);
	unsigned char *buf;
	unsigned char*& rbuf=buf;
	unsigned long count;
	//extract voxel size
	bool bVoxSizeDefined=true;
	double dx,dy,dz;
	OFString ofs;
	//extract voxel size: pixel spacing
	OFCondition ofc=dataset->findAndGetOFStringArray(DCM_PixelSpacing, ofs);
	do
	{
		if(ofc.bad()) {bVoxSizeDefined=false; break;}
		size_t sep=ofs.find('\\');
		if(sep==string::npos) {bVoxSizeDefined=false;break;}
		char tmp[50];
		strncpy(tmp,ofs.data(),sep);
		tmp[sep]=0;
		dx=atof(tmp);
		strcpy(tmp,ofs.data()+sep+1);
		dy=atof(tmp);
		if(dx==0 || dy==0) bVoxSizeDefined=false;
	}while(false);
		
	//extract voxel size: slice thickness
	ofc=dataset->findAndGetOFString(DCM_SliceThickness,ofs);
	do
	{
		if(ofc.bad()) {bVoxSizeDefined=false; break;}
		dz=atof(ofs.data());
		if(dz==0) bVoxSizeDefined=false;
	}while(false);
	if(bVoxSizeDefined) v.SetVoxelDims(dx,dy,dz);
			//extract pixel array
	ofc=dataset->findAndGetUint8Array(DCM_PixelData,(const Uint8* &)rbuf, &count);
	if (ofc.bad()) return false;

	unsigned short *ptr=(unsigned short*)buf;

	long xy,xyz;
	for(int z=0; z<v.SZ(); z++)
	{
		xyz=v.SX()*v.SY()*z;
		for(int y=0; y<v.SY(); y++)
		{
			xy=xyz+v.SX()*y;
			for(int x=0; x<v.SX(); x++)
			{
				v(x,y,z)=(Real)ptr[xy+x];
			}
		}
	}
	return true;
}
bool InitSlice(Volume& v,DcmDataset *dataset, int z)
{
		unsigned char *buf;
		unsigned char*& rbuf=buf;
		unsigned long count;
		//extract voxel size
		bool bVoxSizeDefined=true;
		double dx,dy,dz;
		OFString ofs;
		//extract voxel size: pixel spacing
		OFCondition ofc=dataset->findAndGetOFStringArray(DCM_PixelSpacing, ofs);
		do
		{
			if(ofc.bad()) {bVoxSizeDefined=false; break;}
			size_t sep=ofs.find('\\');
			if(sep==string::npos) {bVoxSizeDefined=false;break;}
			char tmp[50];
			strncpy(tmp,ofs.data(),sep);
			tmp[sep]=0;
			dx=atof(tmp);
			strcpy(tmp,ofs.data()+sep+1);
			dy=atof(tmp);
			if(dx==0 || dy==0) bVoxSizeDefined=false;
		}while(false);
		
		//extract voxel size: slice thickness
		ofc=dataset->findAndGetOFString(DCM_SliceThickness,ofs);
		do
		{
			if(ofc.bad()) {bVoxSizeDefined=false; break;}
			dz=atof(ofs.data());
			if(dz==0) bVoxSizeDefined=false;
		}while(false);
		if(bVoxSizeDefined) v.SetVoxelDims(dx,dy,dz);

		//extract pixel array
		ofc=dataset->findAndGetUint8Array(DCM_PixelData,(const Uint8* &)rbuf, &count);
		if (ofc.bad()) return false;

		unsigned short *ptr=(unsigned short*)buf;

		for(int y=0; y<v.SY(); y++)
		{
			for(int x=0; x<v.SX(); x++)
			{
				v(x,y,z)=(Real)ptr[v.SX()*y+x];
			}
		}
		return true;
}
int	 ExtractFrames(ML3Array<DcmFileFormat*>& ff, std::vector<ML3Array<DcmFileFormat*> >& vols)
{
	int nVols=1,echo,echo0;
	DcmFileFormat* dff0;
	std::vector<int> echos;

	//check for multi-echo sequence.
//	ff.QuickSort(PreceedsEchoNumber);
	ff.Sort(dcmlib::CompareEchos);
	dff0=ff[0];
	echo0=GetEchoNumber(dff0);
	echos.push_back(echo0);
	for(int i=1; i<ff.GetSize(); i++)
	{
		echo=GetEchoNumber(ff[i]);
		if(echo0!=echo)
		{
			echos.push_back(echo);
			echo0=echo;
			nVols++;
		}
		dff0=ff[i];
	}
	if(nVols>1) //we have a multi-echo sequence, populate vols with echos.
	{
		dff0=ff[0];
		echo0=GetEchoNumber(dff0);
		vols.resize(nVols);
		for(int j=0; j<echos.size(); j++){if(echo0==echos[j]) {vols[j].Add(dff0); break;}}
		for(int i=1; i<ff.GetSize(); i++)
		{
			echo=GetEchoNumber(ff[i]);
			for(int j=0; j<echos.size(); j++){if(echo==echos[j]) {vols[j].Add(ff[i]); break;}}
		}
		return nVols;
	}

	//check for pseudo-multiframe
//	ff.QuickSort(PreceedsInstanceNumber);
	ff.Sort(CompareInstanceNumbers);
	int nSets=0, dir=1;
	unsigned short imHt, imWid;
	DcmDataset *dataset=ff[0]->getDataset();
	dataset->findAndGetUint16(DCM_Rows,imHt);
	dataset->findAndGetUint16(DCM_Columns,imWid);
	dir=(GetSliceLocation(ff[0])<GetSliceLocation(ff[1]))?1:-1;
	std::vector<int> stind;
	dff0=ff[0];
	stind.push_back(GetInstanceNumber(dff0));
	for(int i=1; i<ff.GetSize(); i++)
	{
//		cout<<"First slice instance: " << GetInstanceNumber(dff0) << "first slice location: " << GetSliceLocation(dff0) << endl;
//		cout<<"Instance: "<< GetInstanceNumber(ff[i]) << ", slice location: " << GetSliceLocation(ff[i]) << endl;
		if(CompareSliceLocations(dff0,ff[i])*dir>0)
		{
//			cout << "Added start instance number: " << GetInstanceNumber(ff[i]) << endl;
			stind.push_back(GetInstanceNumber(ff[i]));
			nVols++;
		}
		dff0=ff[i];
	}
//	if (nVols!=1) cout << "Warning: series of different size detected" << endl;
//	nVols=(int)(((double)ff.GetSize()/(double)nSets)+.5);
	vols.resize(nVols);
	//2nd pass

// USE THE CODE BELOW IF YOU WANT ONLY SLICES WITH EXACTLY SAME LOCATIONS TO GO INTO THE MULTI-FRAME VOLUME.
// BY DEFAULT, SLICES ARE SORTED BY INSTANCE NUMBER AND MONOTONOUS SLICE LOCATION DIFFERENTIAL.

	// THIS IS INSTANCE-BASED VOLUME SORTING.
	// comment all the following lines out and uncomment the previous block to turn on location-based slice sorting.
	dff0=ff[0];
	int vind=0;
	vols[0].Add(dff0);
	for(int i=1; i<ff.GetSize(); i++)
	{
//		cout<<"First slice instance: " << GetInstanceNumber(dff0) << "first slice location: " << GetSliceLocation(dff0) << endl;
//		cout<<"Instance: "<< GetInstanceNumber(ff[i]) << ", slice location: " << GetSliceLocation(ff[i]) << endl;
		if(CompareSliceLocations(dff0,ff[i])*dir>0)
		{
//			cout << "Added start instance number: " << GetInstanceNumber(ff[i]) << endl;
			stind.push_back(GetInstanceNumber(ff[i]));
			vind++;
		}
		dff0=ff[i];
		vols[vind].Add(ff[i]);
	}
	return nVols;
}
int find_vol(DcmFileFormat* dff, std::vector<int>& stind)
{
	if(stind.size()<2) return 0;
	int in=GetInstanceNumber(dff);
	for(int i=0; i<stind.size()-1; i++)
	{
		if(in<stind[i+1] && in>=stind[i]) return i;
	}
	return stind.size()-1;
}
bool IsMultiFrame(ML3Array<DcmFileFormat*>& ff, long &nFrames)
{
	nFrames=0;
	if (ff.GetSize()==1)
	{
		ff[0]->getDataset()->findAndGetLongInt(DCM_NumberOfFrames,nFrames);
		return (nFrames>1); //extract frames and exit.
	}
	else return false;
}
bool Dcm2Vol(ML3Array<DcmFileFormat*>& ff, Volume*& v, int& nVols, int mode, bool bQuiet)
{	
//	int num=GetInstanceNumber(ff[0]);
	unsigned short imHt, imWid;

	DcmDataset *dataset=ff[0]->getDataset();
	dataset->findAndGetUint16(DCM_Rows,imHt);
	dataset->findAndGetUint16(DCM_Columns,imWid);

	unsigned short par1=0,par2=0;
	dataset->findAndGetUint16(DCM_BitsAllocated,par1);
	dataset->findAndGetUint16(DCM_SamplesPerPixel,par2);
	if(par1!=16 || par2!=1) 
	{
		cout << "Current BitsAllocated: "<<par1<<" (supported 16) , SamplesPerPixel: "<<par2<< "supported 1" << endl;
		return false;
	}
	std::vector<ML3Array<DcmFileFormat*> > vols;

	// check if the image is multi frame
	long nFrames,instInVol;

	if (ff.GetSize()==1)
	{
		nVols=1;
		nFrames=1;
		ff[0]->getDataset()->findAndGetLongInt(DCM_NumberOfFrames,nFrames);
		if (nFrames>1) //extract frames and exit.
		{
			vols.resize(nVols);
			vols[0].Add(ff[0]);
			instInVol = vols[0].GetSize();
		}
	}
	else
	{
		nFrames=1;
		nVols=ExtractFrames(ff, vols);
		instInVol=vols[0].GetSize();
	}

	if (nFrames>1) //extract multi frame volume.
	{
		v=new Volume[1];
		if(!InitVolMultiFrame(v[0],ff[0]->getDataset(),imWid,imHt,nFrames)) 
		{
			delete [] v;
			return false;
		}
		return true;
	}

	for(int j=0; j<nVols; j++)
	{
		if (mode==dcmlib::COMPARE_INSTANCE)
			vols[j].Sort(CompareInstanceNumbers);
		else if (mode==dcmlib::COMPARE_SLICE_LOCATION)
			vols[j].Sort(CompareSliceLocations);
	}
	if (!bQuiet) cout << "Extracting pixels from DICOM..." << endl;
	v=new Volume[nVols];
	for (int i=0; i<nVols; i++)
		v[i].InitMemory(imWid,imHt,instInVol);

	for(int j=0; j<nVols; j++)
	{
		for(int z=0; z<vols[j].GetSize(); z++)
		{
			dataset=vols[j][z]->getDataset();
			if(!InitSlice(v[j],dataset,z))
			{
				delete [] v;
				return false;
			}
		}		
	}
	return true;
}

bool InitOrientationData(ML3Array<DcmFileFormat*>& dcmFiles, Volume& v, _3DINFO& S0, _3DINFO& S1, bool bRescale)
{
	//always calc the dimensions of DICOM dataset.
	S0.d=(unsigned short)dcmFiles.GetSize();
	DcmDataset *dataset1=dcmFiles[0]->getDataset(), *dataset2=dcmFiles[1]->getDataset();
	dataset1->findAndGetUint16(DCM_Rows,S0.h);
	dataset1->findAndGetUint16(DCM_Columns,S0.w);

	if(bRescale) v.Resample(S0.w,S0.h,S0.d);

	S1.w=v.m_dims[0]; S1.h=v.m_dims[1]; S1.d=v.m_dims[2];

	//the rest is only valid when voxel dimensions are specified for both reference and source volumes.
	if(!v.IsVoxDims()) return false;

	//slice thickness
	S0.dz=ExtractFloat(dataset1,DCM_SliceThickness);
	if(S0.dz==0) return false;

	//extract voxel size: pixel spacing
	double tmp[6];
	if(!ExtractFloatArray(dataset1,DCM_PixelSpacing,tmp,2)) return false;
	S0.dx=tmp[0]; S0.dy=tmp[1];
	
	//direction cosines.
	if(!ExtractFloatArray(dataset1,DCM_ImageOrientationPatient,tmp,6)) return false;
	S0.ex1=S1.ex1=tmp[0]; S0.ex2=S1.ex2=tmp[1]; S0.ex3=S1.ex3=tmp[2];
	S0.ey1=S1.ey1=tmp[3]; S0.ey2=S1.ey2=tmp[4]; S0.ey3=S1.ey3=tmp[5];

	//now calculate the direction cosine for z axis.
	ColumnVector v0(ArrayLengthSpecifier(3)), v1(ArrayLengthSpecifier(3));
	v0(1)=tmp[0]; v0(2)=tmp[1]; v0(3)=tmp[2];
	v1(1)=tmp[3]; v1(2)=tmp[4]; v1(3)=tmp[5];
	ColumnVector v2=crossproduct(v0,v1);
	S0.ez1=S1.ez1=v2(1); S0.ez2=S1.ez2=v2(2); S0.ez3=S1.ez3=v2(3);

	S1.dx=v.m_voxel_dims[0]; S1.dy=v.m_voxel_dims[1]; S1.dz=v.m_voxel_dims[2];
	
	//image position 
	if(!ExtractFloatArray(dataset1,DCM_ImagePositionPatient,tmp,3)) return false;
	S0.slx=tmp[0]; S0.sly=tmp[1]; S0.slz=tmp[2];
	if(!ExtractFloatArray(dataset2,DCM_ImagePositionPatient,tmp,3)) return false;
	S1.slx=tmp[0]; S1.sly=tmp[1]; S1.slz=tmp[2];	
	return true;
}

bool ExtractFloatArray(DcmDataset* dataset, const DcmTagKey& key, double* arr, int sz)
{
	OFString ofs;
	OFCondition ofc=dataset->findAndGetOFStringArray(key, ofs);
	if(ofc.bad()) return false;
	std::string s(ofs.data());
	vector<string> vs=ConsoleConfig::Tokenize(s,"\\");
	if(vs.size()<(unsigned int) sz) return false;
	int count=0;
	do
	{
		arr[count]=atof(vs[count].data());
		count++;
	} while(count<sz);
	return true;
}
bool PreceedsEchoNumber(DcmFileFormat*& f1, DcmFileFormat*& f2)
{
	return GetEchoNumber(f1)<GetEchoNumber(f2);
}
int CompareEchos(DcmFileFormat*& f1, DcmFileFormat*& f2)
{
	int num1=GetEchoNumber(f1), num2=GetEchoNumber(f2);
	if(num1<num2) return -1;
	else return (num1==num2)? 0 : 1;
}
bool PreceedsSliceLocation(DcmFileFormat*& f1, DcmFileFormat*& f2)
{
	OFString ofs;
	double num1, num2;
	f1->getDataset()->findAndGetOFString(DCM_SliceLocation,ofs);
	num1=atof(ofs.data());
	f2->getDataset()->findAndGetOFString(DCM_SliceLocation,ofs);
	num2=atof(ofs.data());
	return (num1<num2);
}
int	GetEchoNumber(DcmFileFormat*& f)
{
	DcmDataset *dd=f->getDataset();
	OFString ofs;
	OFCondition ofc=dd->findAndGetOFString(DCM_EchoNumbers,ofs);
	if(ofc.bad()) return 0;
	return atoi(ofs.data());
}
double GetSliceLocation(DcmFileFormat*& f)
{
	DcmDataset *dd=f->getDataset();
	OFString ofs;
	OFCondition ofc=dd->findAndGetOFString(DCM_SliceLocation,ofs);
	if(ofc.bad()) return 0;
	return atof(ofs.data());
}
bool PreceedsInstanceNumber(DcmFileFormat*& f1, DcmFileFormat*& f2)
{
	return GetInstanceNumber(f1)<GetInstanceNumber(f2);
}
int CompareInstanceNumbers(DcmFileFormat*& f1, DcmFileFormat*& f2)
{
	int num1=GetInstanceNumber(f1);
	int num2=GetInstanceNumber(f2);
	if(num1<num2) return -1;
	else return (num1==num2)?0:1;
}
*/
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

/*
unsigned short GetMaxIntensity(Volume& v, ML3Array<DcmFileFormat*>& fformat, bool bRobust )
{
	DcmDataset *dataset;
	unsigned short res=0;
	unsigned char* buf;
	unsigned char*& rbuf=buf;
	Analyze an;
	an.m_dsr.dime.dim[0]=4;
	an.m_dsr.dime.dim[1]=v.m_dims[0];
	an.m_dsr.dime.dim[2]=v.m_dims[1];
	an.m_dsr.dime.dim[3]=1;
	an.m_dsr.dime.dim[4]=1;
	an.m_dsr.dime.datatype=(short)Analyze::DT_SIGNED_SHORT;
	an.m_dsr.dime.bitpix=16;

	Volume v0;
	v0=v;
	for(int i=0; i<fformat.GetSize(); i++)
	{
		dataset=fformat[i]->getDataset();
		dataset->findAndGetUint8Array(DCM_PixelData,(const Uint8* &)rbuf);
		v0.InitMemory(an,buf,i);
	}
	if(bRobust)
		return (unsigned short)HistogramAnalysis::Percentile(v0,99.0);
	else
		return v0.GetMax();
}

bool ReadDcmImage(std::string& fname, DcmFileFormat& dff, Volume &v)
{
	OFCondition error = dff.loadFile(fname.c_str());
	if (error.bad()) return false;
	unsigned short imHt, imWid;

	DcmDataset *dataset = dff.getDataset();
	dataset->findAndGetUint16(DCM_Rows, imHt);
	dataset->findAndGetUint16(DCM_Columns, imWid);

	unsigned short par1 = 0, par2 = 0;
	dataset->findAndGetUint16(DCM_BitsAllocated, par1);
	dataset->findAndGetUint16(DCM_SamplesPerPixel, par2);

	if (par1 != 16 || par2 != 1)
	{
		cout << "Current BitsAllocated: " << par1 << " (supported 16) , SamplesPerPixel: " << par2 << "supported 1" << endl;
		return false;
	}
	v.InitMemory(imWid, imHt, 1);
	return InitSlice(v, dataset, 0);
}
//save specified frame to DICOM file (no modifications to DICOM tags)
bool SaveVolFrameDCM(const char* fname, DcmFileFormat* fformat, Volume& v, int frame)
{
	DcmDataset *dataset = fformat->getDataset();
	//verify that byte format is supported
	unsigned short par1 = 0, par2 = 0, dcmw = 0, dcmh = 0;
	dataset->findAndGetUint16(DCM_Rows, dcmh);
	dataset->findAndGetUint16(DCM_Columns, dcmw);
	long nFramesPerFile = 1;
	dataset->findAndGetLongInt(DCM_NumberOfFrames, nFramesPerFile);
	bool bMultiFrame = (nFramesPerFile>1);
	if (bMultiFrame) return false;
	if (!Is16bpp(dataset)) return false;

	int start_frame = 0,
		end_frame = start_frame + v.m_dims[2] - 1;

	unsigned char *buf;
	unsigned char*& rbuf = buf;
	unsigned long count;

	dataset->findAndGetUint8Array(DCM_PixelData, (const Uint8* &)rbuf, &count);
	unsigned short *ptr = (unsigned short*)buf;
	//	int za=v.m_dims[2]-(frame-g_Origin[2])-1;
	int za = frame;
	for (int y = 0; y<v.m_dims[1]; y++)
		for (int x = 0; x<v.m_dims[0]; x++)
			ptr[dcmw*y + x] = (unsigned short)(v(x, y, frame)+.5);
	cout << "saving " << fname << endl;
	OFCondition ofc=fformat->saveFile(fname);
	return !ofc.bad();
}
*/
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
