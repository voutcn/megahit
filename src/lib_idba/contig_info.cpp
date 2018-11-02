#include "lib_idba/contig_info.h"

#include <fstream>
#include <istream>
#include <ostream>
#include <string>

#include "lib_idba/bit_edges.h"

using namespace std;


istream &operator >>(istream &is, ContigInfo &contig_info)
{
    is.read((char *)&contig_info.in_edges_, sizeof(BitEdges));
    is.read((char *)&contig_info.out_edges_, sizeof(BitEdges));
    is.read((char *)&contig_info.kmer_size_, sizeof(uint16_t));
    is.read((char *)&contig_info.kmer_count_, sizeof(uint32_t));

    int size = 0;
    if (!is.read((char *)&size, sizeof(int)))
        return is;

    contig_info.counts_.resize(size);
    for (int i = 0; i < size; ++i)
        is.read((char *)&contig_info.counts_[i], sizeof(SequenceCountUnitType));

    return is;
}

ostream &operator <<(ostream &os, const ContigInfo &contig_info)
{
    os.write((char *)&contig_info.in_edges_, sizeof(BitEdges));
    os.write((char *)&contig_info.out_edges_, sizeof(BitEdges));
    os.write((char *)&contig_info.kmer_size_, sizeof(uint16_t));
    os.write((char *)&contig_info.kmer_count_, sizeof(uint32_t));

    int size = contig_info.counts_.size();
    os.write((char *)&size, sizeof(int));
    for (int i = 0; i < size; ++i)
        os.write((char *)&contig_info.counts_[i], sizeof(SequenceCountUnitType));

    return os;
}

void ReadContigInfo(const string &filename, deque<ContigInfo> &contig_infos)
{
    contig_infos.clear();
    ifstream fin(filename.c_str(), ios_base::in | ios_base::binary);
    ContigInfo contig_info;
    while (fin >> contig_info)
        contig_infos.push_back(contig_info);
    
}

void WriteContigInfo(const string &filename, const deque<ContigInfo> &contig_infos)
{
    ofstream fout(filename.c_str(), ios_base::out | ios_base::binary);
    for (unsigned i = 0; i < contig_infos.size(); ++i)
        fout << contig_infos[i];
}

