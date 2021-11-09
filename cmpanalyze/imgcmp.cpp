#ifndef _WIN32
#if !defined(HAVE_CONFIG_H)
#define HAVE_CONFIG_H
#endif
#endif

#include <mlib3.h>

#define COEFF_NMI 1
#define COEFF_CC 0
#define COEFF_ABSDIFF 2
#define COEFF_SORENSEN 3
#define COEFF_ETA 4
#define COEFF_CR 5

static int g_x0=0;
static int g_x1=0;
static int g_y0=0;
static int g_y1=0;
static int g_z0=0;
static int g_z1=0;
static int g_NeibSz=0;
static Real g_thresh=0;
static int g_DiceVal=-1;
static char g_vol1[MAXL];
static char g_vol2[MAXL];
static char g_volOut[MAXL];
static string g_mask="";

bool g_bQuiet=false;
int  g_iCoeffType=0;
bool g_bPrintFlips=false;
bool g_bUncrop=false;
bool g_bDWICoef=false;
bool g_bVolOut=false;
bool g_bNeib=false;


bool GetArgs(int argc, char* argv[])
{
	AnyOption opt;

    /* 3. SET THE USAGE/HELP   */
	opt.addUsage( "Compute similarity metric between two Analyze or 4dfp volumes." );

    opt.addUsage( "Usage:");
	opt.addUsage( "imgcmp [options] <vol1> [<vol2>]" );
	opt.addUsage( "Options:");
	opt.addUsage( "\t\tquiet mode\r -q" );
	opt.addUsage( "\t\tprint needed flips\r -f" );
	opt.addUsage( "\t\tlower threshold\r -t <float>" );
	opt.addUsage( "\t\tPrint needed flips and resize parameters for vol2 cf. vol1.\r -u" );
	opt.addUsage( "\t\tuse NMI (CC is the defalut)\r -m" );
	opt.addUsage( "\t\tuse absolute difference (CC is the default)\r -a" );
	opt.addUsage( "\t\tcompute Sorensen-Dice coefficient\r -s" );
	opt.addUsage("\t\toutput to image\r -o <vol>");
	opt.addUsage( "\t\trestrict to specified mask\r --ma <vol>");
	
	opt.addUsage("\t\tisovalue for both volumes to compute Dice coefficient\r --dv <int>");
	opt.addUsage( "\t\tcompute similarity in vicinity of <int> size (-o required)\r --nb <int>" );
	opt.addUsage( "\t\tcompute eta\r --eta" );
	opt.addUsage("\t\tcompute correlation ratio (only with --nb) \r --cr");

	opt.noPOSIX();

	opt.setFlag("eta");
	opt.setFlag("q");
	opt.setFlag("f");
	opt.setFlag("p");
	opt.setFlag("m");
	opt.setFlag("a");
	opt.setFlag("u");	
	opt.setFlag("s");
	opt.setOption("t");
	opt.setOption("ma");
	opt.setOption("dv");
	opt.setOption("nb");
	opt.setFlag("cr");
	
	opt.setOption("o");

	opt.processCommandArgs(argc,argv);
	g_bDWICoef = false;

	if ((opt.getArgc() < 1 && g_bDWICoef) || (!g_bDWICoef && opt.getArgc()<2)) 
	{
		opt.printUsage();
		return false;
	}
	char *s;
	if (s=opt.getValue("o"))
	{
		g_bVolOut=true;
		strcpy(g_volOut,s);
	}
	if (s=opt.getValue("nb"))
	{
		if (!g_bVolOut)
		{
			cout << "imgcmp ERROR: -o option must be specified with --nb option" << endl;
			return false;
		}
		g_NeibSz=atoi(s);
		if (g_NeibSz<=0) {
			cout << "imgcmp ERROR: cannot read size for --nb option" << endl;
			return false;
		}
		g_bNeib=true;
	}

	if (s=opt.getValue("t"))
	{
		g_thresh=(Real)atof(s);
	}
	g_bQuiet=opt.getFlag("q");
	if(opt.getFlag("m")) 
		g_iCoeffType=COEFF_NMI;
	else if (opt.getFlag("a")) 
		g_iCoeffType=COEFF_ABSDIFF;
	else if (opt.getFlag("s")) 
		g_iCoeffType=COEFF_SORENSEN;
	else if (opt.getFlag("eta")) 
		g_iCoeffType=COEFF_ETA;
	else if (opt.getFlag("cr"))
		g_iCoeffType=COEFF_CR;
	g_bPrintFlips=opt.getFlag("f");
	g_bUncrop=opt.getFlag("u");
	if(g_bUncrop) g_bPrintFlips=true;
	strcpy(g_vol1,opt.getArgv(0));
	if(!g_bDWICoef)
	{
		strcpy(g_vol2,opt.getArgv(1));	
	}
	if (s=opt.getValue("ma"))
	{
		g_mask.assign(s);	
	}	
	if (s = opt.getValue("dv")) { g_DiceVal = atoi(s); if (g_DiceVal <= 0) { opt.printUsage(); return false; } }
	if ( g_iCoeffType == COEFF_CR && !g_bNeib )
	{
		cout << "imgcmp ERROR: -cr is only supported for -nb option" << endl;
		return false;
	}
	if (g_bNeib)
	{
		if ( g_iCoeffType != COEFF_NMI && g_iCoeffType != COEFF_ETA && g_iCoeffType != COEFF_CC && g_iCoeffType != COEFF_CR)
		{
			cout << "imgcmp ERROR: -nb option only works for nmi, eta, cc, cr" << endl;
			return false;
		}
	}
	return true;
}
int CompareCoef(Real c1, Real c2)
{
	switch(g_iCoeffType)
	{
		case COEFF_CC:
			return fabs(c1)>fabs(c2)?1:-1;
		default:
			return (c1>c2)? 1:-1;
	}
	return 0;
}

bool CalcCoefNeib(Volume& v1, Volume& v2, Volume &volCoef, int neib_sz)
{
	switch(g_iCoeffType)
	{
		case COEFF_NMI:	return MMath::MutualInfoNeib(v1,v2,volCoef,neib_sz); 
		case COEFF_CC:	return MMath::CorrCoeffNeib(v1,v2,volCoef,neib_sz);
		case COEFF_ETA:	return MMath::EtaNeib(v1,v2,volCoef,neib_sz);
		case COEFF_CR: return MMath::CorrRatioNeib(v1, v2, volCoef, neib_sz);
		default: return false;
	}
}

Real CalcCoef(Volume& v1, Volume& v2, Volume* mask)
{
	Volume res;
	int nEl=v2.m_nElements;
	switch(g_iCoeffType)
	{
		case COEFF_SORENSEN:
			if (g_DiceVal<0)
				return v1.Sorensen(v2,g_thresh,mask);
			else
				return v1.Sorensen1(v2,g_DiceVal);
		case COEFF_NMI:
			return v1.NMI(v2,100,mask);
		case COEFF_ABSDIFF:
			v1.Subtract(v2,res);
			if (mask) nEl=( (*mask) > 0 );
			return res.SumAbs(mask)/nEl;
		case COEFF_ETA:
			return MMath::Eta(v1,v2,mask);
		default:
			return v1.CorrCoeff(v2,mask);
	}
}

int main (int argc, char *argv[])
{
	if(!GetArgs(argc, argv))
		return 0;
	static char rcsid[] = "$Id: " __FILE__ ", " __DATE__ ", " __TIME__;
	if(!g_bQuiet) std::cout << rcsid << endl;
	else Volume::s_verbosity=Volume::VERB_NONE;

	Volume vol1,vol2;
	if(!vol1.Read(g_vol1))
	{
		std::cout << "none" << endl;
		return -1;
	}

	if(!vol2.Read(g_vol2))
	{
		std::cout << "none" << endl;
		return -1;
	}
	int nbins=50;
	Volume mask, *pmask;
	if (g_mask.length()<1) pmask=NULL;
	else if (!mask.Read((char*)g_mask.c_str()))
	{
		std::cout << "none" << endl;
		return -1;
	}
	else {
		pmask=&mask;		
	}
	if (g_bNeib)
	{
		Volume volcoef;
		if (!CalcCoefNeib(vol1,vol2,volcoef,g_NeibSz))
		{
			cout << "imgcmp ERROR computing local coefficient" << g_volOut << endl;
			exit (-1);
		}
		if (!volcoef.Write(g_volOut,vol1.m_coord_format,argc,argv)){
			cout << "imgcmp ERROR writing file " << g_volOut << endl;
			exit (-1);
		}
		exit (0);		
	}
	Real coef=CalcCoef(vol1,vol2,pmask);
	if( !g_bPrintFlips )
	{
		if (g_bQuiet)
			std::cout<<coef<<endl;
		else
			std::cout << "Coefficient: " << coef << endl;
		return 0;
	}
	if(g_bPrintFlips)
	{
		Volume res;
		Real cc, bestCC=0;
		int bestFl=0;
		Volume vol3;
		int cx,cy,cz;
		for(int fl=0; fl<8; fl++)
		{
			vol3=vol2;
			if( (fl & 1)>0 ) vol3.Flip(1);
			if( (fl & 2)>0 ) vol3.Flip(2);
			if( (fl & 4)>0 ) vol3.Flip(3);

			if(g_bUncrop)
			{
				if(vol1.Uncrop(vol3,cx,cy,cz))
				{
					if(g_bQuiet)
					{
						std::cout << (((fl & 1)>0)?1:0) << " " << (((fl & 2)>0)?1:0) << " " << (((fl & 4)>0)?1:0) << endl;
						std::cout << cx << " " << cy << " " << cz << " " << vol1.SX() <<" " << vol1.SY()<< " " << vol1.SZ() << endl;
					}
					else
					{
						std::cout << "Flips: " << (((fl & 1)>0)?1:0) << " " << (((fl & 2)>0)?1:0) << " " << (((fl & 4)>0)?1:0) << endl;
						std::cout << "Crop: " << cx << " " << cy << " " << cz << " " << vol1.SX() <<" " << vol1.SY()<< " " << vol1.SZ() << endl;
					}
//					getchar();
					return 0;
				}

			}
			cc=CalcCoef(vol1,vol3,pmask);
			if(CompareCoef(cc,bestCC)>0)
			{
				bestCC=cc;
				bestFl=fl;
			}
			if(!g_bQuiet)
				std::cout << "Flip: " << fl << ", c=" << cc << endl;
		}
		if(g_bQuiet) std::cout << (((bestFl & 1 )>0)?1:0) << " " << (((bestFl & 2)>0)?1:0) << " " << (((bestFl & 4)>0)?1:0) << " " << bestCC << endl;
		else	std::cout << "Flips: " << (((bestFl & 1 )>0)?1:0) << " " << (((bestFl & 2)>0)?1:0) << " " << (((bestFl & 4)>0)?1:0) << ", " << "coef = " << bestCC << endl;
	}
	return 0;
}