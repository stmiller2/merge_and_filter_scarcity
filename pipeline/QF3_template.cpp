#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
using namespace std;

constexpr int q_floor = {{Q_FLOOR}};      // placeholder for q_floor, set by params.env
constexpr int q_cutoff = {{Q_CUTOFF}};    // placeholder for q_cutoff, set by params.env
constexpr float cutoff_pct = {{CUTOFF_PCT}}; // placeholder for cutoff_pct, set by params.env

int main()
{
	//Initialize variables and open file streams
	string inMPERs="combined.fastq";
	ifstream inRawReads(inMPERs.c_str());
	ofstream outGoodReads("good_reads.csv"), outPoorReads("poor_reads.csv");
	size_t split;				//position of the split between sequence & quality string
	double q, g, pct;			//Q-score, # bases > q_cutoff, % bases > q_cutoff
	int j;						//sequence position when looping through quality string
	string line, seq, qual;		//combined.fastq line, DNA sequence, quality string
	bool discarded;				//whether read contains bases < q_floor
	long readCount = 0;			//number of reads actually processed

	//Bail out clearly instead of silently writing empty output files if
	//combined.fastq is missing (e.g. an earlier merge/formatting step failed)
	if(!inRawReads.is_open())
	{
		cerr << "Error: could not open input file '" << inMPERs << "'. Was combined.fastq created by the merge/formatting step?" << '\n';
		return 1;
	}
	if(!outGoodReads.is_open() || !outPoorReads.is_open())
	{
		cerr << "Error: could not open 'good_reads.csv' / 'poor_reads.csv' for writing." << '\n';
		return 1;
	}

	//Print quality filtering variables to the console
	cout << "           Quality filtering settings:" << '\n' << '\t';
	cout << "           Proportion >=Q" << q_cutoff << " must be above " << cutoff_pct << '\n' << '\t';
	cout << "           All bases must have a Q-score >= Q" << q_floor << '\n';

	while(getline(inRawReads,line))
	{
		readCount++;
		split = line.find(" ");				//Get split position
		qual = line.substr(split + 1);		//Get quality sequence
		seq = line.substr(0, split);		//Get DNA sequence

		g = 0;
		discarded = false;
		
		for(j=0; j < qual.length(); j++)	//For each character in the quality string
		{
			q = double(qual[j]) - 33;		//Convert character to Q-score (ASCII - 33)
			if(q < q_floor)					//If base quality is < q_floor,
			{
				discarded = true;			//mark sequence as discarded
				break;						//and break the for loop (ignore further bases)
			} else
			{
				if(q >= q_cutoff)			//If base quality is >= q_cutoff
				{
					g = g + 1;				//Add one to g (good bases count)
				}
			}
		}
		if(discarded == false)				//If sequence isn't already marked discarded,
		{
			pct = g / qual.length();		//Calculate percent bases >= q_cutoff
			if(pct > cutoff_pct)			//If percent meets standards,
			{
				outGoodReads << seq << '\n';//Write seq to good_reads.csv
			} else							//If not,
			{								//Write seq to poor_reads.csv
				//Include "E2" classification and actual percent >=q_cutoff
				outPoorReads << seq << "," << "E2" << "," << pct << '\n';
			}
		} else								//If sequence is already marked discarded,
		{									//Write seq to poor_reads.csv
			//Include "E1" classification and position of first subthreshold base
			outPoorReads << seq << "," << "E1" << "," << j << '\n';
		}
	}
	inRawReads.close();
	outGoodReads.close();
	outPoorReads.close();

	//A merge/formatting failure often produces a present-but-empty
	//combined.fastq, which the is_open() checks above don't catch on their
	//own -- flag that case explicitly too.
	if(readCount == 0)
	{
		cerr << "Error: '" << inMPERs << "' contained no reads. Upstream merging/formatting may have failed." << '\n';
		return 1;
	}

	cout << "           High quality reads written to 'good_reads.csv'" << '\n';
	cout << "           Low quality reads written to 'poor_reads.csv'" << '\n' << '\t';
	cout << "           'E1' - 1+ bases < Q" << q_floor << " (first base index shown)" << '\n' << '\t';
	cout << "           'E2' - < " << cutoff_pct << " Q" << q_cutoff << " (actual percentage shown)" << '\n';

	return 0;
}