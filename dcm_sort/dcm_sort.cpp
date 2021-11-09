#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <set>
#include <map>
#include <sys/stat.h>
//#include <mlib3.h>
#include <mlib3.h>
#include "../dcmlib/dcmlib.h"


static char rcsid[]="$Id: " __FILE__ ", " __DATE__ ", " __TIME__;
//some global definitions
static char g_dcmdir[MAXL]="";
using namespace dcmlib;

bool GetArgs(int argc, char* argv[])
{
	AnyOption opt;

    /* 3. SET THE USAGE/HELP   */
    opt.addUsage( "Usage: dcm_sort <DICOM dir>" );
	opt.setOption("d");
	opt.processCommandArgs(argc,argv);
	if(opt.getArgc()<1)
	{
		opt.printUsage();
		return false;
	}
	strcpy(g_dcmdir,opt.getArgv(0));
	return true;
}

bool TestFunction()
{
	return true;
}

typedef struct sd
{
	char seqName[100];
	char serNum[20];
	char serDesc[100];
	char serUID[100];
	char serDir[100];
	int count;
} __SP;

void GetStr(DcmDataset* dataset, DcmTagKey key, char* str)
{
	str[0]=0;
	OFString ofs;
	dataset->findAndGetOFString(key,ofs);
	if(ofs.length()>0)
		strcpy(str,ofs.data());
}
void RemoveChar(char* str, char remove)
{
	int len=strlen(str);
	if(len<1) return;
	char* temp=new char[len+1];
	int pos1=0;
	for(int i=0; i<len; i++)
	{
		if(str[i]!=remove) temp[pos1++]=str[i];
	}
	temp[pos1]=0;
	strcpy(str,temp);
}

std::multimap<std::string,__SP*> sm;
bool CompareScans(std::string s1, std::string s2)
{
	std::multimap<std::string,__SP*>::iterator pos1, pos2;
	pos1=sm.find(s1);
	pos2=sm.find(s2);
	__SP *ser1=sm.find(s1)->second,
		*ser2=sm.find(s2)->second;
	int n1=atoi(ser1->serNum), n2=atoi(ser2->serNum);
	return n1<n2;
}

int main (int argc, char *argv[])
{
	static char rcsid[] = "$Id: " __FILE__ ", " __DATE__ ", " __TIME__;
	cout << rcsid << endl;
	if (!TestFunction()) return 0;
	if( !GetArgs(argc, argv)) return 0;
	FileList fl;
	if(!fl.MakeList(g_dcmdir))
	{
		cout << "No files found in specified DICOM folder";
		return 0;
	}
	DcmDataset *dataset;

	OFString str;
//	std::multimap<std::string,__SP*> sm;
	std::multimap<std::string,__SP*>::iterator pos;
	__SP* ser;
	std::vector<std::string> scans;
	std::string serUID;
	std::set<std::string> dirs;
	char command[FILENAME_MAX], stname[FILENAME_MAX];
	DcmFileFormat f;
	char fpath[FILENAME_MAX];
	char* nextFile=fl.GetFirstFile();
    OFCondition error;
	int cnt=0,len=0;
	char temp[20];
	typedef struct _sinfo {std::string suid; std::string dir;} SINFO;

	for(int i=0; i<fl.NFiles(); i++)
	{
		fprintf(stderr,"\rReading file %d out of %d", i+1, fl.NFiles());
		if(i>0) nextFile=fl.GetNextFile();
#ifdef _WIN32
		sprintf(fpath,"%s\\%s",g_dcmdir,nextFile);
#else
		sprintf(fpath,"%s/%s",g_dcmdir,nextFile);
#endif
		if(!LooksLikeDICOMFile(fpath)) continue;	
		f.clear();
		error=dcmlib::LoadDcmFile(fpath,f);
		if (error.bad())
		{
			cout << "Error " << error.text() << "reading DICOM from " << nextFile << endl;
			continue;
		}
		dataset=f.getDataset();
		//get series inst UID
		dataset->findAndGetOFString(DCM_SeriesInstanceUID,str);
		serUID=std::string(str.data());
		if(sm.find(serUID)==sm.end())
		{
			ser=new __SP;
			ser->count=1;
			GetStr(dataset,DCM_SequenceName,ser->seqName);
			GetStr(dataset,DCM_SeriesInstanceUID,ser->serUID);
			RemoveChar(ser->seqName,' ');
			RemoveChar(ser->seqName,'*');
			if(strlen(ser->seqName)==0) strcpy(ser->seqName,"none");
			GetStr(dataset,DCM_SeriesDescription,ser->serDesc);
			RemoveChar(ser->serDesc,' ');
			if(strlen(ser->serDesc)==0) strcpy(ser->serDesc,"none");
			GetStr(dataset,DCM_SeriesNumber,ser->serNum);
			RemoveChar(ser->serNum,' ');

			scans.push_back(serUID);			
			cnt=0;
			struct stat st;
			do
			{				
				if(cnt==0)
				{
					len=sprintf(stname, "study%s",ser->serNum);
					sprintf(temp,"%s",ser->serNum);
				}
				else
				{
					len=sprintf(stname, "study%s00%d",ser->serNum, cnt);
					sprintf(temp,"%s00%d",ser->serNum,cnt);
				}
				if (len<260) stname[len]=0;
				else {
					cout << "study name is too long, cannot conitnue" << endl;
					exit(-1);
				}
				cnt++;
			} while (dirs.find(std::string(stname))!=dirs.end());
			strcpy(ser->serDir, stname);
			strcpy(ser->serNum,temp);
			dirs.insert(std::string(stname));
			sm.insert(std::make_pair(serUID,ser));
			sprintf(command, "mkdir %s", stname);
//			cout << endl << command << endl;
			system(command);
		}
		else
		{
			pos=sm.find(serUID);
			ser=pos->second;
			(ser->count)++;
		}
#ifndef WIN32
		sprintf(command, "cp -P %s %s",fpath,ser->serDir);
//		cout << endl << command << endl;
		system(command);
#else
		sprintf(command, "copy %s %s",fpath,ser->serDir);
//		cout << endl << command << endl;
		system(command);
#endif
	}
	std::sort(scans.begin(),scans.end(),CompareScans);
	ofstream file("DICOM.studies.txt");
	if(!file)
	{
		cout << "Cannot write DICOM.studies.txt file" << endl;
	}
	// perform some file << "..." operations	
	cout << endl;
	for(unsigned int i=0; i<scans.size(); i++)
	{
		pos=sm.find(scans[i]);
		ser=pos->second;
		file << std::left << std::setw(10) << ser->serNum << " " <<std::setw(15)<< ser->seqName << " " <<
			std::setw(30) << ser->serDesc << " " << std::setw(5) << ser->count << endl;
		cout << std::left << std::setw(10) << ser->serNum << " " <<std::setw(15)<< ser->seqName << " " <<
			std::setw(30) << ser->serDesc << " " << std::setw(5) << ser->count << endl;
	}
	//free memory
	for(pos=sm.begin(); pos!=sm.end(); ++pos)
	{
		delete pos->second;
	}
	cout << endl;
	return 0;
}