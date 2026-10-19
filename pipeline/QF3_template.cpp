#include <fstream>
#include <iostream>
#include <string>
using namespace std;

constexpr int q_floor = {{Q_FLOOR}};      // placeholder for q_floor, set by params.env
constexpr int q_cutoff = {{Q_CUTOFF}};    // placeholder for q_cutoff, set by params.env
constexpr float cutoff_pct = {{CUTOFF_PCT}}; // placeholder for cutoff_pct, set by params.env
constexpr int end_slop = {{END_SLOP}};    // placeholder for end_slop (FALSE -> 0), set by params.env
constexpr bool trim_slop = {{TRIM_SLOP}}; // placeholder for trim_slop, set by params.env

int main(int argc, char* argv[])
{
	//Which side(s) of each combined read hold a sequencing-read 5' end. Passed per
	//sample by merge_reads.sh, since it depends on the read mode:
	//  L  = left only  (single-end R1)
	//  R  = right only (single-end R2; reverse-complemented, so its 5' end is on the right)
	//  LR = both       (PEAR-merged reads, or R1 + revcomp(R2) concatenation)
	//end_slop is only ever applied at these 5' ends, never at a 3' end.
	string ends = (argc == 2) ? argv[1] : "";
	if(ends != "L" && ends != "R" && ends != "LR")
	{
		cerr << "Usage: QF3.out <L|R|LR>  (side(s) of the read that are sequencing 5' ends)" << '\n';
		return 1;
	}
	const int slop_left  = (ends == "L" || ends == "LR") ? end_slop : 0;
	const int slop_right = (ends == "R" || ends == "LR") ? end_slop : 0;

	//Input is streamed on stdin, one read per line as "SEQUENCE QUALITY"
	//(produced on the fly by merge_reads.sh; nothing is written to disk first)
	ios::sync_with_stdio(false);
	cin.tie(nullptr);
	istream& inRawReads = cin;
	ofstream outGoodReads("good_reads.csv"), outPoorReads("poor_reads.csv");
	size_t split;				//position of the split between sequence & quality string
	double q, g, pct;			//Q-score, # bases > q_cutoff, % bases > q_cutoff
	int j;						//sequence position when looping through quality string
	int len, start, stop;		//read length; quality-checked region is [start, stop)
	int offset;					//shift so reported positions index into the written sequence
	string line, seq, qual;		//input line, DNA sequence, quality string
	string outSeq;				//sequence as written to output (slop removed if trim_slop)
	bool discarded;				//whether read contains bases < q_floor
	long readCount = 0;			//number of reads actually processed

	if(!outGoodReads.is_open() || !outPoorReads.is_open())
	{
		cerr << "Error: could not open 'good_reads.csv' / 'poor_reads.csv' for writing." << '\n';
		return 1;
	}

	//Print quality filtering variables to the console
	cout << "           Quality filtering settings:" << '\n' << '\t';
	cout << "           Proportion >=Q" << q_cutoff << " must be above " << cutoff_pct << '\n' << '\t';
	cout << "           All bases must have a Q-score >= Q" << q_floor << '\n';
	if(end_slop > 0)
	{
		cout << '\t' << "           End slop: " << slop_left << " base(s) on left, " << slop_right
		     << " on right exempt from filtering (" << (trim_slop ? "trimmed from" : "kept in") << " output)" << '\n';
	}

	while(getline(inRawReads,line))
	{
		readCount++;
		split = line.find(" ");				//Get split position
		qual = line.substr(split + 1);		//Get quality sequence
		seq = line.substr(0, split);		//Get DNA sequence

		len = qual.length();
		start = slop_left;					//Skip 5' slop on the left
		stop = len - slop_right;			//Skip 5' slop on the right

		//Read too short to have anything left after removing slop
		if(stop <= start)
		{
			//Include "E3" classification and read length
			outPoorReads << seq << "," << "E3" << "," << len << '\n';
			continue;
		}

		outSeq = trim_slop ? seq.substr(start, stop - start) : seq;
		offset = trim_slop ? start : 0;

		g = 0;
		discarded = false;

		for(j=start; j < stop; j++)			//For each character in the checked region
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
			pct = g / (stop - start);		//Calculate percent of checked bases >= q_cutoff
			if(pct > cutoff_pct)			//If percent meets standards,
			{
				outGoodReads << outSeq << '\n';//Write seq to good_reads.csv
			} else							//If not,
			{								//Write seq to poor_reads.csv
				//Include "E2" classification and actual percent >=q_cutoff
				outPoorReads << outSeq << "," << "E2" << "," << pct << '\n';
			}
		} else								//If sequence is already marked discarded,
		{									//Write seq to poor_reads.csv
			//Include "E1" classification and position of first subthreshold base
			//(0-based index into the sequence as written)
			outPoorReads << outSeq << "," << "E1" << "," << j - offset << '\n';
		}
	}
	outGoodReads.close();
	outPoorReads.close();

	if(!outGoodReads || !outPoorReads)
	{
		cerr << "Error: failed while writing 'good_reads.csv' / 'poor_reads.csv' (disk full?)." << '\n';
		return 1;
	}
	if(readCount == 0)
	{
		cerr << "Error: no reads received on stdin. Upstream merging/formatting may have failed." << '\n';
		return 1;
	}

	cout << "           High quality reads written to 'good_reads.csv'" << '\n';
	cout << "           Low quality reads written to 'poor_reads.csv'" << '\n' << '\t';
	cout << "           'E1' - 1+ bases < Q" << q_floor << " (first base index shown)" << '\n' << '\t';
	cout << "           'E2' - < " << cutoff_pct << " Q" << q_cutoff << " (actual percentage shown)" << '\n';
	if(end_slop > 0)
	{
		cout << '\t' << "           'E3' - read too short to filter after removing end slop (read length shown)" << '\n';
	}

	return 0;
}