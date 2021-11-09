#ifndef _CONSOLECONFIG_INCLUDED_
#define _CONSOLECONFIG_INCLUDED_
#include <vector>
#include <string>
#include <ostream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <map>

#include "mlib3array.h"

#include "mlib3volume.h"
//#include "../pugixml/pugixml.hpp"
#ifndef _LODEPNG_INCLUDED
#include "lodepng.h"
#define _LODEPNG_INCLUDED
#endif

#ifdef _WIN32
#else
#include <dirent.h>
#endif

#ifdef _4DFP
extern "C"
{
	#include <rec.h>
	#include <ifh.h>
}
#endif

#define MAXL 256 //max string length

struct MLIB3_JSW {
	std::stringstream ss;
	int nm = 0;

	void prep(char* name) {
		if (!nm) ss << "{"; else ss << ","; nm++;
		ss << "\"" << name << "\"" << ":";
	}

	template<class T>
	void add_val(char* name, T& val){
		prep(name);
		ss << val;
	};
	void add_vector(char* name, ColumnVector& v)
	{
		prep(name);
		ss << "[";
		for (int i = 1; i <= v.size(); i++)
			if (i == 1) ss << v(i); else ss << "," << v(i);
		ss << "]";
	};

	string get_str() {
		ss << "}";
		return ss.str();
	};

};
//JSON_WRITER;

class ConsoleConfig
{
public:
	bool m_bAHeader;
	char m_program[MAXL], m_root[MAXL], m_control, **m_argv, m_rcsid[MAXL];
	int m_argc;

	ConsoleConfig(char* rcsid)
	{
		m_control='\0';
		m_bAHeader=false;
		strcpy(m_rcsid,rcsid);
		cout << m_rcsid << endl;
	};
	~ConsoleConfig(){};
	//read a space separated array of numbers from file. Lines starting with '#' are ignored.
	static bool ReadRealMatrixFromFile(const char* file, Matrix& out)
	{
		std::vector<std::vector<Real> > data;
		std::ifstream ifs(file);
		if (!ifs.is_open()) return false;
		std::string ln; Real val;
		std::vector<Real> ldata;
		int w=0,h=0;
		while (std::getline(ifs,ln))
		{
			if (ln[0]=='#') continue;
			std::stringstream ls(ln);
			ldata.clear();
			while (ls >> val) ldata.push_back(val);
			if (ldata.size() < 1 ) continue;
			if (w==0) w=ldata.size();
			if (ldata.size() != w) return false;
			data.push_back(ldata);
		}
		h=data.size();
		out.resize(h,w);
		for (int j=0; j<h; j++)	for (int i=0; i<w; i++)	out(j+1,i+1)=(data[j])[i];
		return true;
	}
	static bool ReadMapFromFile(const char* file, std::map<Real,Real>& vmap)
	{
		Matrix m;
		if (!ReadRealMatrixFromFile(file,m)) return false;
		if (m.ncols()<1 || m.nrows()<2) return false;
		for (int i=1; i<=m.ncols(); i++)
		{
			//debug
			//cout << m(1,i) << " mapped to " << m(2,i) << endl; 
			vmap[m(1,i)]=m(2,i);
		}
		return true;
	}
	//return optimal image mosaic dimensions based on the number of images, 
	//maximum number of images (-1 if all images to be displayed), 
	//and aspect ratio computed as wid:ht
	//evenly distributed item indices are stored in includedItems.
	static bool ComputeMosaicDims(int nItems, vector<int>& includedItems, Real& wid, Real& ht, int maxItems=-1)
	{
		includedItems.clear();
		if (wid<=0 || ht <=0 || nItems<1 || maxItems<1) return false;
		if (maxItems==-1 || nItems<maxItems) {	for (int i=1; i<=nItems; i++) includedItems.push_back(i); }
		else
		{
			Real st=(Real)(nItems-1)/(Real)(maxItems-1); int n=0;
			int ind; Real k=1.0;
			for (int i=1; i<maxItems; i++, k+=st){
				ind=(int)round(k); includedItems.push_back(ind);
				n++;
			}
			includedItems.push_back(nItems); n++; nItems=n;
		}
		Real scale_factor=wid/ht; 
		Real h=sqrt(nItems/scale_factor); 
		Real w=h*scale_factor;
		int n1=(int)(ceil(h)*floor(w));
		int n2=(int)(floor(h)*ceil(w));
		int r1=n1-nItems, r2=n2-nItems;
		wid=ht=-1;
		if (r2 >= 0) {
			if ((r2 <= r1) || (r1<0)) { wid = ceil(w); ht = floor(h); }
		}
		if ((wid<0) && (r1>=0)){
			if ( (r1<=r2) || (r2<0) ) { wid = floor(w); ht = ceil(h); }
		}
		if (wid<0) { wid = ceil(w); ht = ceil(h); }
		return true;
	}
	static vector<string> Tokenize(const string& str, const string& delimiters)
	{
		vector<string> tokens;
	    	
		// skip delimiters at beginning.
    		string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	    	
		// find first "non-delimiter".
    		string::size_type pos = str.find_first_of(delimiters, lastPos);

    		while (string::npos != pos || string::npos != lastPos)
    		{
        		// found a token, add it to the vector.
        		tokens.push_back(str.substr(lastPos, pos - lastPos));
			
        		// skip delimiters.  Note the "not_of"
        		lastPos = str.find_first_not_of(delimiters, pos);
			
        		// find next "non-delimiter"
        		pos = str.find_first_of(delimiters, lastPos);
    		}

		return tokens;
	}
	static void Matr2UCVect(Matrix& m, std::vector<unsigned char>& vec, bool b16bit, Real scale=1.0)
	{
		unsigned char v;
		for (int i = 1; i <= m.nrows(); i++)
			for (int j = 1; j <= m.ncols(); j++)
			{
				//assuming big endian, according to lodepng documentation
				if (b16bit)
				{
					v = (unsigned char)((int)(m(i, j)*scale) / 256);
					vec.push_back(v);
					v = (unsigned char)(((int)(m(i, j)*scale)) % 256);
					vec.push_back(v);
				}
				else
				{
					vec.push_back((unsigned char)(m(i, j)*scale));
				}
			}
	}
	static void SaveMatrixAsImage(const char* fname, Matrix& matr)
	{
		Real mn=matr.minimum(),mx=matr.maximum();
		Real hm = 255.0 / (mx-mn);
		Matrix mnorm; mnorm=matr; mnorm -=mn; mnorm *= hm; 
		vector<unsigned char> vec;
		ConsoleConfig::Matr2UCVect(mnorm,vec,false);
		lodepng::encode(fname, vec, mnorm.ncols(), mnorm.nrows(), LCT_GREY, 8);
	}

	static void GetCoords(const char* param, double* buf, int cnt)
	{			
		string ln(param);
		string delim(",-");
		vector<string> vals=Tokenize(ln,delim);
//		int x0=0,y0=0,z0=0,x1=0,y1=0,z1=0;
		if(cnt>=3)
		{
			buf[0]=atof(vals.at(0).data());
			buf[1]=atof(vals.at(1).data());
			buf[2]=atof(vals.at(2).data());
			if(cnt>=6)
			{
				buf[3]=atof(vals.at(3).data());
				buf[4]=atof(vals.at(4).data());
				buf[5]=atof(vals.at(5).data());
			}
		}
	}

#ifdef _4DFP
	void WriteRec(char* suffix)
	{
#ifdef __GNUG__
		char name[MAXL], imgfile[MAXL];
		if (suffix) sprintf(name, "%s_%s.4dfp.rec", m_root, suffix);
		else sprintf(name, "%s.4dfp.rec", m_root);
		startrece(name, m_argc, m_argv, m_rcsid, IsBigEndian() ? 'b' : 'l');
		sprintf(imgfile, "%s.4dfp.img", m_root);
		catrec(imgfile);
		endrec();
#endif
	}
	static bool ReadVolume4dfp(Volume &v, char* name)
	{
		char root[MAXL],fname[MAXL];
		GetRoot(name,root);
		if(v.Read4dfp(root)) return true;
		sprintf(fname,"%s.4dfp",root);
		return v.Read4dfp(fname);
	};
#endif

	bool		ReadVolume(Volume &v, char* fname=NULL)
	{
		char name[MAXL];
		GetRoot((fname)?fname:m_root,name);
#ifdef _4DFP
		if(!m_bAHeader)
		{
			if(v.Read4dfp(name)) return true;
			sprintf(name,"%s.4dfp",name);
			if(v.Read4dfp(name)) return true;
		}
#endif
		GetRoot((fname)?fname:m_root,name);
		if(!v.ReadAnalyze(name))
		{
			sprintf(name,"%s.4dint",name);
			if(!v.ReadAnalyze(name))
			{
				fprintf(stderr, "\nError reading %s\n",m_root);
				return false;
			}
		}
		return true;
	}
	void		SaveVolume(Volume& v, char* suf, bool bWriteRec=false)
	{
		char nm[MAXL];
		GetRootShort(nm);
		printf("Writing %s_%s\n",nm,suf);
#ifdef _4DFP
		if(!m_bAHeader)
		{
			if(bWriteRec) WriteRec(suf);
			v.Write4dfp(m_root,m_argc,m_argv);
			return;
		}
#endif
		sprintf(nm+strlen(nm),"_%s",suf);
		v.WriteAnalyze(nm);
	}
	static void SaveMatr(Matrix& M, ofstream& ofs, int prec=6){
		int acc; Real val;
		int nRow = M.Nrows(), nCol = M.Ncols();
		vector<int> wid(nCol);

		for (int i = 0; i<nCol; i++) {
			wid[i] = 1;
			for (int j = 1; j <= nRow; j++) {
				val = M(j, i + 1); acc = 1;
				while (fabs(val) >= 10) { val *= .1; acc++; }
				if (val<0) acc++;
				wid[i] = _mx(wid[i], acc);
			}
		}
		for (int i = 1; i <= nRow; i++) {
			for (int j = 1; j <= nCol; j++)
				ofs << fixed << setprecision(prec) << setw(prec + 2 + wid[j - 1]) << M(i, j);
			ofs << endl;
		}
	}
	static void PrintMatr(Matrix& M, int prec=6, int nCol=-1, int nRow=-1){
		int acc; Real val;
		if (nCol==-1) nCol = M.Ncols();
		if (nRow==-1) nRow = M.Nrows();
//		int nRow=M.Nrows(), 
		int old_precision = (int)cout.precision(prec);
		vector<int> wid(nCol);

		for (int i = 0; i<nCol; i++) {
			wid[i] = 1;
			for (int j = 1; j <= nRow; j++) {
				val = M(j, i + 1); acc = 1;
				while (fabs(val) >= 10) { val *= .1; acc++; }
				if (val<0) acc++;
				wid[i] = _mx(wid[i], acc);
			}
		}
		for (int i = 1; i <= nRow; i++) {
			for (int j = 1; j <= nCol; j++)
				cout << fixed << setw(prec+2+wid[j - 1]) << M(i, j);
			cout << endl;
		}
		cout.precision(old_precision);
	}
	static void PrintMatr(string msg, Matrix& m, int prec=6){
		cout << endl << msg.c_str() << endl;
		PrintMatr(m,prec);
	};
	static void PrintCV(string msg, ColumnVector& cv, int prec){
		cout << msg.c_str() << " ";
		for(int i=0; i<cv.nrows(); i++)
			cout << fixed << setprecision(prec) << cv(i+1) << " ";
		cout << endl;
	}
	static void MPrint(double* data, int nCols, int n, int prec=6)
	{
		for(int i=0; i<n; i++)
		{
			cout << setprecision(prec) << data[i] << " ";
			if((i>0 || nCols==1) && !(i%nCols) && i!=(n-1)) cout << "..." << endl;
		}
		cout << endl;
	}
	void		GetRootShort(char* root)
	{
		char* ptr1=strrchr((char*)m_root,'/'),
			*ptr2=strrchr((char*)m_root,'\\');
		if(ptr1<ptr2) ptr1=ptr2;
		if(!ptr1) ptr1=(char*)m_root;
		ptr2=strstr(ptr1,".4dfp");
		if(!ptr2) ptr2=strstr(ptr1,".4dint");
		if(!ptr2) sprintf (root, "%s", m_root);
		else
		{
			int len=(int)(ptr2-m_root);
			strncpy(root,m_root,len);
			root[len]=0;
		}
	}
	static void GetRoot(const char* name, char* root)
	{
		char	*str;
		strcpy (root, name);
		if (str = strrchr (root, '.')) 
			*str='\0';
	};
	static std::string AddSuffix(char* fname, char* suff)
	{
		//append suffix to output file name if needed.
		std::stringstream ss;
		std::string outstr;
		unsigned int sl=strlen(suff);
		if (sl>0 && strlen(fname)>sl) {
				if (strncmp(fname + strlen(fname) - sl, suff, sl)) ss << fname << suff;
				else ss << fname;
		}
		else ss << fname;
		return ss.str();
	}

	static char* GetFile(char* name)
	{
		char* str;
		if(str=strrchr(name, '\\'))
			return ++str;
		else if (str=strrchr(name,'/'))
			return ++str;
		else 
			return name;
	}

	static void GetDir(const char* name, char* dir)
	{
		char* str;
		strcpy(dir,name);
		if(str=strrchr(dir, '\\'))
			*str='\0';
		else if (str=strrchr(dir,'/'))
			*str='\0';
		else
			dir[0]=0;
	}
	bool ProcessCommandLine(int argc, char* argv[])
	{
		m_argc=argc; m_argv=argv;
		char* ptr, str[MAXL],c;
		if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; else ptr++;
		strcpy (m_program, ptr);

		int k,i;
		for (k = 0, i = 1; i < argc; i++) 
		{
			if (*(argv[i]) == '-')
			{
				strcpy(str,argv[i]); 
				ptr=str;
				while(c=*(ptr++))
					switch(c)
					{
						case '@': m_control= *ptr++; *ptr='\0'; break;
						case 'a': m_bAHeader=true; break;
					}

			}
			else if(k==0) 
			{
				GetRoot (argv[i], m_root); k++;
			}
		}
		if(k<1) return false;

		return true;
	};

};

class FileList
{
public:
	FileList(){m_nFiles=0;};
	~FileList(){}//for(int i=0; i<m_nFiles; i++) delete[] m_Files[i];};

#ifndef _WIN32
	bool MakeList(char* folder)
	{
		DIR *dir; 
		char str[MAXL];
		struct dirent *ent;
		int len;
		m_nFiles=0;
		if (dir=opendir(folder)){
			while (ent=readdir(dir)){
				len=strlen(ent->d_name);
				strcpy(str, ent->d_name);
                //cout << "file: " << str << endl;
                if (!strcmp(str,".") || !strcmp(str,"..")) //ignore first two entries.
                {
                    //cout << " ignored" << endl;
                    continue;
                }
				m_Files.push_back(std::string(str));
				//m_Files.Add(str);
				m_nFiles++;
			}
			return true;
		}
		else return false;
	};
#else
	bool MakeList(char* folder)
	{
		
        cout << "win32 defined" << endl;
		static int maxlen=32768;
		char command[MAXL];
		FILE* fp;
		sprintf(command,"dir /A-D /B %s > _dirlist_.tmp",folder);
		system(command);
		fp=fopen("_dirlist_.tmp","r");
		if(!fp) return false;
		m_nFiles=0;
		int len;
		char str[MAXL];
		while(fgets(command, MAXL,fp))
		{
			if(!strcmp(command, "_dirlist_.tmp")) continue;
			if((len=(int)strlen(command))>0)
			{
				//str=new char[len];
				strncpy(str,command,len-1);
				str[len-1]=0;
				m_Files.push_back(std::string(str));
				//m_Files.Add(str);
				m_nFiles++;
			}
		}
		fclose(fp);
		system("del _dirlist_.tmp");
		return true;
	};
#endif
	std::string* GetFirstFile_str()
	{
		if(m_nFiles<1) return NULL;
		m_CurFile=1;
		return (std::string*)(&(m_Files[0]));		
	}
	
	char* GetFirstFile()
	{
		if(m_nFiles<1) return NULL;
		m_CurFile=1;
		return (char*)m_Files[0].c_str();
	};
	std::string* GetNextFile_str()
	{
		if(m_CurFile>=m_nFiles) return NULL;
		return (std::string*)(&(m_Files[m_CurFile++]));
	};	
	char* GetNextFile()
	{
		if(m_CurFile>=m_nFiles) return NULL;
		return (char*)m_Files[m_CurFile++].c_str();
	};
	int NFiles(){return m_nFiles;};
	static int cmp(char*& s1, char*& s2){int r=strcmp(s1,s2); if (r<0) return -1; return (r==0)?0:1;};
	static int cmp1(std::string& s1, std::string& s2)
	{
		return s1.compare(s2);
	}
	void Sort(){
		//m_Files.Sort(FileList::cmp);
		std::sort(m_Files.begin(),m_Files.end());
	};

private:
	int m_nFiles, m_CurFile;
	//ML3Array<char*> m_Files;
	std::vector<string> m_Files;
};
#endif //CONSOLECONFIG_INCLUDED
