#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "Handy.h"

//g++ *.cpp -L/seq/mbrd/mbusby/Software/bamtools/lib -I/seq/mbrd/mbusby/Software/bamtools/include -I/seq/mbrd/mbusby/Software/bamtools/include/api

/*==========================================================================
This parses the bam file to generate values by transcript
Reports:
Reads
Unique reads
Pairs where both ends align 
Duplicates
Unique Duplicates


==========================================================================*/

unsigned int checkErrors();
void displayHelp();
void somethingsGoneWrong(unsigned int);
void processBamFile();
void getUnique();

//Variables taken from the inputs
string bamFileName="";
string outFileName="";
int nReads=0;
int seed=4503;
bool alignedOnly=true;
bool mateAligned=true;
bool ignoreDups=true;


map<string, double> randMap;

using namespace std;
using namespace BamTools;

/*
 *Useage: 
 * ./GetRandomByRead -bam /seq/mbrd/mbusby/analyseCoverage/RSEM_alignments_sorted/sample17.sorted.bam -out  /seq/mbrd/mbusby/analyseCoverage/chromosomeCounts/sample17cts.txt -nReads 1000000 -seed 5403 -aligned_only true
 *
 *
*/
int main(int argc, char* argv[]) 
{
		
	//Intialize library of useful things
	Handy h(0);
	
	int optind=1;

	while ((optind < argc) && (argv[optind][0]=='-')) 
	{	
        string sw = argv[optind];
				
		if (sw=="-h") 
		{	
            optind++;
			displayHelp();
			return 1;
        }
		
		else if(optind >=  argc-1)
		{
			cerr<<"Your final parameter, "<<sw<<" is missing a value."<<endl;
			return 1;
		}

		else if (h.cmpStringNoCase(sw, "-bam")==1)
		{	
            optind++;
			bamFileName = argv[optind];		
			optind++;
        }
		
		else if (h.cmpStringNoCase(sw, "-out")==1)
		{	
            optind++;
			outFileName = argv[optind];		
			optind++;
        }
		
		
		else if (h.cmpStringNoCase(sw, "-nReads")==1)		
		{	
            optind++;
			nReads=h.getIntFromString(argv[optind]);
			optind++;
        }
		
		
		else if (h.cmpStringNoCase(sw, "-seed")==1)
		{	
            optind++;
			seed=h.getIntFromString(argv[optind]);
			optind++;
        }
		
		else if (h.cmpStringNoCase(sw, "-aligned_only")==1)
		{	
            optind++;
			if(h.cmpStringNoCase(argv[optind], "false")==1)
			{
				alignedOnly=false;
				mateAligned=false;
			}
			optind++;
        }
		
		else if (h.cmpStringNoCase(sw, "-mate_aligned")==1)
		{	
            optind++;
			if(h.cmpStringNoCase(argv[optind], "false")==1)
			{
				mateAligned=false;
			}
			optind++;
        }
		
		else if (h.cmpStringNoCase(sw, "-ignore_dups")==1)
		{	
            optind++;
			if(h.cmpStringNoCase(argv[optind], "true")==1)
			{
				ignoreDups=true;
			}
			optind++;
        }
		
		else
		{
			cerr<<"Main: Unknown parameter:"<<sw<<endl;
			return 1;
		}
	}	
	
	checkErrors();
	
	cout<<"OutfileName after check errors:"<<outFileName<<endl;

	getUnique();
	cout<<"Unique Reads Found."<<endl;
	
	cout<<"OutfileName after getUnique:"<<outFileName<<endl;
	
	processBamFile();
	cout<<"Done\n";
		
}

void processBamFile()
{
	Handy h(0);
	BamReader reader;
	BamAlignment al;
	RefVector refVector;
	double cut=((double) nReads) / ((double) randMap.size() );
	
	cout<<"Cut is "<<cut;
	cout<<"Writing to: "<<outFileName<<endl;
	
	//Open output stream to write the findings for each transcript
	ofstream outputStream;
	outputStream.open(outFileName.c_str());
		
	if ( !reader.Open(bamFileName) ) {
		cerr << "Could not open input BAM file." << endl;
		return;
	}
	
	BamWriter::CompressionMode compressionMode = BamWriter::Compressed;
	
    // open BamWriter
    BamWriter writer;	
	
    // get BamReader metadata
    const string headerText = reader.GetHeaderText();
    const RefVector references = reader.GetReferenceData();
    if ( references.empty() ) {
        cout << "bamtools random ERROR: no reference data available...Trying to carry on." << endl;        
    }
	
    writer.SetCompressionMode(compressionMode);
    if ( !writer.Open(outFileName, headerText, references) ) {
        cerr << "bamtools random ERROR: could not open " << outFileName
             << " for writing... Aborting." << endl;
        reader.Close();
        return;
    }
	
	int i=0;
	
	if(alignedOnly == true || ignoreDups == true )
	{
			
		while(reader.GetNextAlignment(al))
		{	
			++i;	
				
			if(i%500000==0)
			{
				cout<<"Writing line "<<i<<" of alignment."<<endl;
				cout.flush();
			}

			if(randMap.find(al.Name)!=randMap.end() && randMap[al.Name]<cut )
			{
				writer.SaveAlignment(al);
			}
		}
	}
	else
	{
			
		while(reader.GetNextAlignment(al))
		{	
			++i;	
				
			if(i%500000==0)
			{
				cout<<"Reading line "<<i<<" of alignment."<<endl;
				cout.flush();
			}

			if(randMap[al.Name]<cut)
			{
				writer.SaveAlignment(al);
			}
		}
	}
	
	
	
	reader.Close();
	writer.Close();	
		
}


void getUnique()
{
	Handy h(0);
	BamReader reader;
	BamAlignment al;
	RefVector refVector;
	string readName;
	
	srand (seed);
	
	if ( !reader.Open(bamFileName) ) 
	{
		cerr << "Could not open input BAM file." << endl;
		return;
	}
	
	int i=0;
		
	if(alignedOnly==true && mateAligned == false)
	{
		while(reader.GetNextAlignment(al))
		{	
			++i;		
			if(i%10000000==0)
			{
				cout<<"Reading line "<<i<<" of alignment for coverage map."<<endl;
			}
			
			if( al.IsMapped() == true && al.IsMateMapped() == true   )
			{
				randMap.insert(map< string, double >::value_type(al.Name, (double)rand()  / ((double)RAND_MAX)  ));					
			}				
		
		}
	}
	else if(alignedOnly==true && mateAligned == true)
	{
		while(reader.GetNextAlignment(al))
		{	
			++i;		
			if(i%10000000==0)
			{
				cout<<"Reading line "<<i<<" of alignment for coverage map."<<endl;
			}
			
			if( al.IsMapped() == true && al.IsMateMapped() == true && randMap.find(al.Name)==randMap.end()  )
			{
				randMap.insert(map< string, double >::value_type(al.Name, (double)rand()  / ((double)RAND_MAX)  ));					
			}				
		
		}
	}
	else if(ignoreDups==true && alignedOnly==true && mateAligned == true)
	{
		while(reader.GetNextAlignment(al))
		{	
			++i;		
			if(i%10000000==0)
			{
				cout<<"Reading line "<<i<<" of alignment for coverage map."<<endl;
			}
			
			if( al.IsMapped() == true && al.IsMateMapped() == true &&  al.IsDuplicate() == false  && randMap.find(al.Name)==randMap.end()  )
			{
				randMap.insert(map< string, double >::value_type(al.Name, (double)rand()  / ((double)RAND_MAX)  ));					
			}				
		
		}
	}
	
	else 
	{
		while(reader.GetNextAlignment(al))
		{	
			++i;		
			if(i%10000000==0)
			{
				cout<<"Reading line "<<i<<" of alignment for coverage map."<<endl;
			}
			
			if(randMap.find(al.Name)==randMap.end())
			{
				randMap.insert(map< string, double >::value_type(al.Name, (double)rand()  / ((double)RAND_MAX)  ));					
			}				
		
		}
	}
	
	
	reader.Close();
}



/*===========================================================================================
Check that all of the necessary fields exist.
===========================================================================================*/

unsigned int checkErrors()
{
	//Errors 

	int err=0;
	Handy h(0);
	string problems="";
	
	if(bamFileName.length()==0)
	{
		problems.append("A bam file containing the reads is needed (-bam).\n");
		++err;
	}
	if(outFileName.length()==0)
	{
		if(bamFileName.length()>4)
		{	
			unsigned found = bamFileName.rfind("/");
			try
			{
				if (found!=std::string::npos)
				{		
					outFileName=bamFileName.substr(found+1, bamFileName.length()-5-(int) found)+".downsampled.bam";
					cout<<" first OutFileName: "<<outFileName<<endl;
					
				}
				
			}
			catch(int e)
			{
				outFileName='out.bam';
				cout<<"OutfileName after catch:"<<outFileName<<endl;
			}
			
			cout<<"OutfileName after try:"<<outFileName<<endl;
		}		
		else
		{
			outFileName='out.bam';
		}
		
		cout<<"OutfileName after inner loop:"<<outFileName<<endl;
	}
	
	if(nReads==0)
	{
		problems.append("No reads selected.  All done! or use -nReads option.");
		++err;
	}
		
		
	
	/*============================================================
	Check that all the relevant files can be read/written to
	==============================================================*/

	err=err+h.checkRead(bamFileName);	

	return err;
	
}


void displayHelp()
{

	cout<<"Required parameters:\n";	
	cout<<"-bam Name of the bam file\n";
	cout<<"-out Name of the output file"<<endl;
	cout<<"-nReads The number of reads you want."<<endl;
	cout<<"Optionall:\n";
	cout<<"-seed A random number seed."<<endl;
	cout<<"-aligned_only Include only aligned reads. Default is true.  Set to False to include all. "<<endl;
	cout<<"-mate_aligned Read and its mate are aligned. . "<<endl;
	cout<<"-ignore_dups Include only reads where both pairs are aligned an ignore marked duplicate reads and do not include them in output. Default is false. Assumes duplicates are marked."<<endl;


}

void somethingsGoneWrong(string whatsGoneWrong)
{

	cout<<"ERROR: Something has gone horribly wrong.\n";	
	cout<<whatsGoneWrong;
	cout<<"\n";	
	cerr<<"\nPlease try again.";
	
}
