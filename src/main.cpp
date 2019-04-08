#include <iostream>
#include <fstream>
#include "Solver.h"
#include "myUsage.h"

using namespace std;

int main(int argc, char *argv[]) {

    if (argc != 3) {
        cerr << "usage: ./solver <input_file> <output_file>" << endl;
        exit(-1);
    }
    MyUsage usage;
    usage.reset();
    ofstream file(argv[2]);
    Solver* solver = new Solver();
    solver->read(argv[1]);

    /*
    solver->debug_dump(file);
    return 0;
    */

    solver->solve();
    solver->dump(file);
    usage.report(1, 1);
    return 0;
}
