#include <iostream>
#include <utility>
#include <fstream>

#include "bc.h"

using namespace std;

//int main()
int main(int argc, char* argv[])
{

    ifstream seqFile(argv[1]);
    ifstream barFile(argv[2]);
    //ifstream seqFile("../test.txt");
    //ifstream barFile("../bench_cells.txt");

    vector<string> sequences;
    vector<string> barcodes;
    vector<string> read_names;

    string read_name = "", seq = "", temp = "";

    while (!seqFile.eof())
    {
        seqFile >> read_name >> seq;
        read_names.push_back(read_name);
        sequences.push_back(seq);
    }
    while (getline(barFile, temp)) {
        barcodes.push_back(temp);
    }

    cout << "There are " << sequences.size() << " sequences and " << barcodes.size() << " barcodes" << endl;

    //barcodeMatch(sequences, barcodes, read_names);
    barcodeMatch(sequences, barcodes, read_names, atof(argv[3]), atof(argv[4]), atof(argv[5]),
        atoi(argv[6]), atoi(argv[7]), atoi(argv[8]), atof(argv[9]), 0.05, 6, argv[10]);

    return(0);
}
