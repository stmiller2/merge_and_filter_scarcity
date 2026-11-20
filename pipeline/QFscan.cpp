// Used by preview_filters.sh. The per-read decision must stay identical to QF3_template.cpp.
// QF3 stores cutoff_pct as a float, so it is compared as (double)(float)cutoff_pct here too.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

const vector<int> FLOORS = {0, 5, 10, 15, 20};
const vector<int> CUTOFFS = {0, 5, 10, 15, 20, 25, 30};
const vector<double> PCTS = {0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95};

typedef vector<long> Key;
Key make_key(int f, int c, double p) { return Key{f, c, lround(p * 10000)}; }

struct Counts { double e1 = 0, e2 = 0, e3 = 0, good = 0; };

bool on_grid(int f, int c, double p)
{
	bool pf = find(FLOORS.begin(), FLOORS.end(), f) != FLOORS.end();
	bool pc = find(CUTOFFS.begin(), CUTOFFS.end(), c) != CUTOFFS.end();
	bool pp = false;
	for(double x : PCTS) if(fabs(x - p) < 1e-9) pp = true;
	return pf && pc && pp;
}

// -----------------------------------------------------------------------------
// scan
// -----------------------------------------------------------------------------
int scan(int argc, char* argv[])
{
	if(argc != 14)
	{
		cerr << "Usage: QFscan scan <sample> <L|R|LR> <end_slop> <cur_floor> <cur_cutoff> <cur_pct>"
		        " <reads_in_file> <reads_sampled> <pear_unassembled> <pear_discarded> <grid.csv> <samples.csv>\n";
		return 1;
	}
	string sample = argv[2], ends = argv[3];
	int end_slop = atoi(argv[4]);
	int cur_floor = atoi(argv[5]), cur_cutoff = atoi(argv[6]);
	double cur_pct = atof(argv[7]);
	long reads_in_file = atol(argv[8]), reads_sampled = atol(argv[9]);
	long pear_un = atol(argv[10]), pear_disc = atol(argv[11]);
	string grid_path = argv[12], samples_path = argv[13];

	if(ends != "L" && ends != "R" && ends != "LR") { cerr << "QFscan: bad ends '" << ends << "'\n"; return 1; }
	const int slop_left  = (ends == "L" || ends == "LR") ? end_slop : 0;
	const int slop_right = (ends == "R" || ends == "LR") ? end_slop : 0;

	vector<int> floors = FLOORS;   floors.push_back(cur_floor);
	vector<int> cutoffs = CUTOFFS; cutoffs.push_back(cur_cutoff);
	vector<double> thr;
	for(double p : PCTS) thr.push_back(double(float(p)));
	thr.push_back(double(float(cur_pct)));
	const size_t nf = floors.size(), nc = cutoffs.size(), np = thr.size();

	vector<double> good(nf * nc * np, 0), e2(nf * nc * np, 0), e1(nf, 0);
	double e3 = 0;
	long n_reads = 0;
	bool qseen[256] = {false};

	ios::sync_with_stdio(false);
	cin.tie(nullptr);
	string line, qual;
	vector<long> gc(nc);
	while(getline(cin, line))
	{
		n_reads++;
		size_t split = line.find(' ');
		qual = line.substr(split + 1);
		int len = qual.length();
		for(int j = 0; j < len; j++) qseen[(unsigned char)qual[j]] = true;
		int start = slop_left, stop = len - slop_right;
		if(stop <= start) { e3++; continue; }

		int minq = 1 << 30;
		fill(gc.begin(), gc.end(), 0);
		for(int j = start; j < stop; j++)
		{
			int q = int(qual[j]) - 33;
			if(q < minq) minq = q;
			for(size_t ci = 0; ci < nc; ci++) if(q >= cutoffs[ci]) gc[ci]++;
		}
		double region = stop - start;
		for(size_t fi = 0; fi < nf; fi++)
		{
			if(minq < floors[fi]) { e1[fi]++; continue; }
			for(size_t ci = 0; ci < nc; ci++)
			{
				double frac = gc[ci] / region;
				size_t base = (fi * nc + ci) * np;
				for(size_t pi = 0; pi < np; pi++)
				{
					if(frac > thr[pi]) good[base + pi]++;
					else e2[base + pi]++;
				}
			}
		}
	}

	vector<int> qvals;
	for(int k = 0; k < 256; k++) if(qseen[k]) qvals.push_back(k - 33);

	bool new_grid;
	{ ifstream t(grid_path); new_grid = !t.good() || t.peek() == ifstream::traits_type::eof(); }
	ofstream g(grid_path, ios::app);
	if(!g) { cerr << "QFscan: cannot write " << grid_path << '\n'; return 1; }
	if(new_grid) g << "sample,q_floor,q_cutoff,cutoff_pct,e1,e2,e3,good,qf_input_reads,total_input_reads,pct_good\n";
	auto write_row = [&](int f, int c, double p, size_t fi, size_t ci, size_t pi)
	{
		size_t k = (fi * nc + ci) * np + pi;
		char pbuf[32]; snprintf(pbuf, sizeof pbuf, "%.4g", p);
		char pct[32];  snprintf(pct, sizeof pct, "%.4f", reads_sampled > 0 ? 100.0 * good[k] / reads_sampled : 0.0);
		g << sample << ',' << f << ',' << c << ',' << pbuf << ',' << long(e1[fi]) << ',' << long(e2[k]) << ','
		  << long(e3) << ',' << long(good[k]) << ',' << n_reads << ',' << reads_sampled << ',' << pct << '\n';
	};
	for(size_t fi = 0; fi < FLOORS.size(); fi++)
		for(size_t ci = 0; ci < CUTOFFS.size(); ci++)
			for(size_t pi = 0; pi < PCTS.size(); pi++)
				write_row(FLOORS[fi], CUTOFFS[ci], PCTS[pi], fi, ci, pi);
	if(!on_grid(cur_floor, cur_cutoff, cur_pct))
		write_row(cur_floor, cur_cutoff, cur_pct, nf - 1, nc - 1, np - 1);
	g.close();

	bool new_samples;
	{ ifstream t(samples_path); new_samples = !t.good() || t.peek() == ifstream::traits_type::eof(); }
	ofstream s(samples_path, ios::app);
	if(!s) { cerr << "QFscan: cannot write " << samples_path << '\n'; return 1; }
	if(new_samples) s << "sample,reads_in_file,reads_sampled,pear_unassembled,pear_discarded,qf_input_reads,q_values\n";
	s << sample << ',' << reads_in_file << ',' << reads_sampled << ',' << pear_un << ',' << pear_disc << ',' << n_reads << ',';
	for(size_t k = 0; k < qvals.size(); k++) s << (k ? " " : "") << qvals[k];
	s << '\n';
	s.close();

	if(!g || !s) { cerr << "QFscan: error writing output (disk full?)\n"; return 1; }
	if(n_reads == 0) cerr << "QFscan: warning: no reads received for sample " << sample << '\n';
	return 0;
}

// -----------------------------------------------------------------------------
// render
// -----------------------------------------------------------------------------
struct Block
{
	string title;
	double total = 0, pear_un = 0, pear_disc = 0;
	map<Key, Counts> grid;
	set<int> qvals;
};

bool use_color = false, use_utf8 = false;

const int PAL[9] = {53, 55, 61, 25, 30, 36, 77, 149, 220};
const int PAL_FG[9] = {255, 255, 255, 255, 255, 16, 16, 16, 16};

string fmt_int(double v)
{
	string s = to_string(llround(v)), out;
	int n = s.size();
	for(int i = 0; i < n; i++)
	{
		out += s[i];
		if(s[i] != '-' && (n - i - 1) % 3 == 0 && i != n - 1) out += ',';
	}
	return out;
}

string strf(const char* fmt, double v) { char b[64]; snprintf(b, sizeof b, fmt, v); return b; }

string bold(const string& s) { return use_color ? "\033[1m" + s + "\033[0m" : s; }

string shade(const string& cell, int bin, bool cur)
{
	if(!use_color) return cell;
	return "\033[48;5;" + to_string(PAL[bin]) + "m\033[38;5;" + to_string(PAL_FG[bin]) + "m" +
	       (cur ? "\033[1m" : "") + cell + "\033[0m";
}

string rule(const string& title)
{
	string bar = use_utf8 ? "─" : "-";
	string s;
	for(int i = 0; i < 3; i++) s += bar;
	s += " " + title + " ";
	int pad = 80 - 5 - (int)title.size();
	for(int i = 0; i < max(pad, 3); i++) s += bar;
	return bold(s);
}

string bar(double frac, int width)
{
	if(frac <= 0) return "";
	if(!use_utf8) return string(max(1, (int)lround(frac * width)), '#');
	static const char* EIGHTHS[8] = {"▏", "▎", "▍", "▌", "▋", "▊", "▉", "█"};
	int eighths = max(1, (int)lround(frac * width * 8));
	string s;
	for(int i = 0; i < eighths / 8; i++) s += EIGHTHS[7];
	if(eighths % 8) s += EIGHTHS[eighths % 8 - 1];
	return s;
}

string qvals_str(const set<int>& q)
{
	if(q.empty()) return "none";
	string s;
	if(q.size() <= 10)
	{
		for(int v : q) s += (s.empty() ? "" : " ") + to_string(v);
		return s + "  (binned: thresholds between adjacent values behave identically)";
	}
	return to_string(*q.begin()) + "-" + to_string(*q.rbegin()) + " (" + to_string(q.size()) + " distinct values)";
}

void print_block(const Block& b, int cur_floor, int cur_cutoff, double cur_pct,
                 const string& end_slop_label, bool merge, bool slop_on)
{
	Key cur_key = make_key(cur_floor, cur_cutoff, cur_pct);
	auto pct = [&](double v) { return b.total > 0 ? 100.0 * v / b.total : 0.0; };

	cout << '\n' << rule(b.title) << '\n';
	cout << "Q values present: " << qvals_str(b.qvals) << '\n';

	auto it = b.grid.find(cur_key);
	if(it != b.grid.end())
	{
		const Counts& c = it->second;
		cout << "\nEstimated fates at current settings (q_floor=" << cur_floor << ", q_cutoff=" << cur_cutoff
		     << ", cutoff_pct=" << strf("%g", cur_pct) << ", end_slop=" << end_slop_label << "):\n";
		auto line = [&](const string& label, double v)
		{ cout << "  " << label << string(max(1, 44 - (int)label.size()), ' ') << strf("%6.2f%%", pct(v)) << '\n'; };
		line("Good (passed all filters)", c.good);
		if(merge)
		{
			line("PEAR: not assembled", b.pear_un);
			line("PEAR: discarded", b.pear_disc);
		}
		line("Poor E1: >=1 base below Q" + to_string(cur_floor), c.e1);
		line("Poor E2: <=" + strf("%g", cur_pct) + " of bases at >=Q" + to_string(cur_cutoff), c.e2);
		if(slop_on || c.e3 > 0) line("Poor E3: too short after removing end slop", c.e3);
	}

	cout << "\nE1 (>=1 base below q_floor), % of input reads. Same for every q_cutoff/cutoff_pct:\n";
	double e1max = 0;
	vector<double> e1s;
	for(int f : FLOORS)
	{
		double v = pct(b.grid.at(make_key(f, CUTOFFS[0], PCTS[0])).e1);
		e1s.push_back(v);
		e1max = max(e1max, v);
	}
	for(size_t i = 0; i < FLOORS.size(); i++)
		cout << "  q_floor " << strf("%2.0f", FLOORS[i]) << "  " << strf("%6.2f%%", e1s[i]) << "  "
		     << bar(e1max > 0 ? e1s[i] / e1max : 0, 50) << '\n';
	if(merge)
		cout << "  (PEAR losses, " << strf("%.2f%%", pct(b.pear_un + b.pear_disc))
		     << " of input reads, happen before filtering and are the same in every cell below)\n";

	double vmin = 1e300, vmax = -1e300;
	for(int f : FLOORS) for(int c : CUTOFFS) for(double p : PCTS)
	{
		double v = pct(b.grid.at(make_key(f, c, p)).good);
		vmin = min(vmin, v); vmax = max(vmax, v);
	}
	auto bin_of = [&](double v) { return vmax > vmin ? min(8, (int)((v - vmin) / (vmax - vmin) * 9)) : 8; };

	cout << "\nGOOD READS, % of input reads, by setting. [ ] = current settings.\n";
	if(use_color)
	{
		cout << "Colour scale: " << strf("%.1f%% ", vmin);
		for(int k = 0; k < 9; k++) cout << shade("   ", k, false);
		cout << strf(" %.1f%%", vmax) << '\n';
	}

	bool cur_on_grid = on_grid(cur_floor, cur_cutoff, cur_pct);
	for(int f : FLOORS)
	{
		cout << '\n' << bold("q_floor = Q" + to_string(f)) << (f == cur_floor ? "  (current)" : "") << '\n';
		cout << "q_cutoff  cutoff_pct:\n         ";
		for(double p : PCTS) cout << strf(" %5.2f ", p);
		cout << '\n';
		for(int c : CUTOFFS)
		{
			cout << (c == cur_cutoff ? bold(strf("%7.0f  ", c)) : strf("%7.0f  ", c));
			for(double p : PCTS)
			{
				double v = pct(b.grid.at(make_key(f, c, p)).good);
				bool cur = cur_on_grid && f == cur_floor && c == cur_cutoff && fabs(p - cur_pct) < 1e-9;
				string cell = cur ? strf("[%5.1f]", v) : strf(" %5.1f ", v);
				cout << shade(cell, bin_of(v), cur);
			}
			cout << '\n';
		}
	}
	if(!cur_on_grid)
		cout << "\n(Current settings are not on the grid; their estimate is in the breakdown above.)\n";
}

vector<string> split_csv(const string& line)
{
	vector<string> out;
	string field;
	stringstream ss(line);
	while(getline(ss, field, ',')) out.push_back(field);
	if(!line.empty() && line.back() == ',') out.push_back("");
	return out;
}

int render(int argc, char* argv[])
{
	if(argc != 12)
	{
		cerr << "Usage: QFscan render <grid.csv> <samples.csv> <cur_floor> <cur_cutoff> <cur_pct>"
		        " <end_slop_label> <merge TRUE|FALSE> <color 0|1> <utf8 0|1> <all 0|1>\n";
		return 1;
	}
	string grid_path = argv[2], samples_path = argv[3];
	int cur_floor = atoi(argv[4]), cur_cutoff = atoi(argv[5]);
	double cur_pct = atof(argv[6]);
	string end_slop_label = argv[7];
	bool merge = string(argv[8]) == "TRUE";
	use_color = string(argv[9]) == "1";
	use_utf8 = string(argv[10]) == "1";
	bool show_all = string(argv[11]) == "1";
	bool slop_on = end_slop_label != "FALSE" && atoi(end_slop_label.c_str()) > 0;

	vector<string> order;
	map<string, Block> blocks;
	map<string, double> weight;
	map<string, double> in_file;
	ifstream sf(samples_path);
	if(!sf) { cerr << "QFscan: cannot read " << samples_path << '\n'; return 1; }
	string line;
	getline(sf, line);
	while(getline(sf, line))
	{
		vector<string> v = split_csv(line);
		if(v.size() < 7) continue;
		Block& b = blocks[v[0]];
		order.push_back(v[0]);
		double file_n = atof(v[1].c_str()), sampled = atof(v[2].c_str());
		b.total = sampled;
		b.pear_un = atof(v[3].c_str());
		b.pear_disc = atof(v[4].c_str());
		stringstream qs(v[6]); int q;
		while(qs >> q) b.qvals.insert(q);
		weight[v[0]] = sampled > 0 ? file_n / sampled : 0;
		in_file[v[0]] = file_n;
	}
	if(order.empty()) { cerr << "QFscan: no samples to show\n"; return 1; }

	ifstream gf(grid_path);
	if(!gf) { cerr << "QFscan: cannot read " << grid_path << '\n'; return 1; }
	getline(gf, line);
	while(getline(gf, line))
	{
		vector<string> v = split_csv(line);
		if(v.size() < 8 || !blocks.count(v[0])) continue;
		Counts& c = blocks[v[0]].grid[make_key(atoi(v[1].c_str()), atoi(v[2].c_str()), atof(v[3].c_str()))];
		c.e1 = atof(v[4].c_str()); c.e2 = atof(v[5].c_str()); c.e3 = atof(v[6].c_str()); c.good = atof(v[7].c_str());
	}

	Block pooled;
	double total_in_file = 0;
	for(const string& s : order)
	{
		const Block& b = blocks[s];
		double w = weight[s];
		pooled.total += w * b.total;
		pooled.pear_un += w * b.pear_un;
		pooled.pear_disc += w * b.pear_disc;
		pooled.qvals.insert(b.qvals.begin(), b.qvals.end());
		total_in_file += in_file[s];
		for(const auto& kv : b.grid)
		{
			Counts& c = pooled.grid[kv.first];
			c.e1 += w * kv.second.e1; c.e2 += w * kv.second.e2;
			c.e3 += w * kv.second.e3; c.good += w * kv.second.good;
		}
	}
	pooled.title = "ALL " + to_string(order.size()) + " SAMPLES (weighted by read depth, ~" + fmt_int(total_in_file) + " reads)";

	for(const string& s : order)
		blocks[s].title = "SAMPLE " + s + " (" + fmt_int(blocks[s].total) + " of " + fmt_int(in_file[s]) + " reads sampled)";

	if(order.size() == 1)
	{
		print_block(blocks[order[0]], cur_floor, cur_cutoff, cur_pct, end_slop_label, merge, slop_on);
		return 0;
	}

	print_block(pooled, cur_floor, cur_cutoff, cur_pct, end_slop_label, merge, slop_on);

	Key cur_key = make_key(cur_floor, cur_cutoff, cur_pct);
	size_t w = 6;
	for(const string& s : order) w = max(w, s.size());
	bool any_e3 = slop_on;
	for(const string& s : order) if(blocks[s].grid.at(cur_key).e3 > 0) any_e3 = true;

	cout << '\n' << rule("PER SAMPLE, CURRENT SETTINGS (% of input reads)") << '\n';
	string hdr = "  " + string("Sample") + string(w - 6, ' ') + "     Reads    Good";
	if(merge) hdr += "  PEAR-NA  PEAR-D";
	hdr += "      E1      E2";
	if(any_e3) hdr += "      E3";
	cout << hdr << '\n';
	for(const string& s : order)
	{
		const Block& b = blocks[s];
		const Counts& c = b.grid.at(cur_key);
		auto p = [&](double v) { return b.total > 0 ? 100.0 * v / b.total : 0.0; };
		string fi = fmt_int(in_file[s]);
		cout << "  " << s << string(w - s.size(), ' ') << string(max(1, 10 - (int)fi.size()), ' ') << fi
		     << strf("%8.1f", p(c.good));
		if(merge) cout << strf("%9.1f", p(b.pear_un)) << strf("%8.1f", p(b.pear_disc));
		cout << strf("%8.1f", p(c.e1)) << strf("%8.1f", p(c.e2));
		if(any_e3) cout << strf("%8.1f", p(c.e3));
		cout << '\n';
	}

	if(show_all)
		for(const string& s : order)
			print_block(blocks[s], cur_floor, cur_cutoff, cur_pct, end_slop_label, merge, slop_on);
	return 0;
}

int main(int argc, char* argv[])
{
	string mode = argc > 1 ? argv[1] : "";
	if(mode == "scan") return scan(argc, argv);
	if(mode == "render") return render(argc, argv);
	cerr << "Usage: QFscan <scan|render> ...\n";
	return 1;
}
