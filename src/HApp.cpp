#include "HDual.h"
#include "HTimer.h"
#include "HTester.h"

#include <set>
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

void solvePlain(const char *filename);
void solveTasks(const char *filename);
void solveMulti(const char *filename, const char *partitionfile = 0);

int main(int argc, char **argv) {
    if (argc == 2) {
        solvePlain(argv[1]);
    } else if (argc == 3) {
        string choice = argv[1];
        if (choice == "-sip") {
            solveTasks(argv[2]);
        }
        if (choice == "-pami") {
            solveMulti(argv[2]);
        }

    } else if (argc == 4) {
        string choice = argv[1];
        if (choice == "multi") {
            solveMulti(argv[2], argv[3]);
        }
        if (choice == "-pami") {
            HModel model;
            model.intOption[INTOPT_PRINT_FLAG] = 1;
            model.intOption[INTOPT_PERMUTE_FLAG] = 1;
            model.dblOption[DBLOPT_PAMI_CUTOFF] = atof(argv[3]);
            model.setup(argv[2]);

            HDual solver;
            solver.solve(&model, HDUAL_VARIANT_MULTI, 8);

            model.printResult();
        }
        if (choice == "-repeat") {
            HTester tester;
            tester.setup(argv[2]);
            tester.testUpdate(atoi(argv[3]));
        }
    } else if (argc == 5) {
        string choice1 = argv[1];
        string choice2 = argv[3];
        double cutoff = atof(argv[4]);
        if (choice1 == "multi" && choice2 == "cutoff") {
            HModel model;
            model.dblOption[DBLOPT_PAMI_CUTOFF] = cutoff;
            model.intOption[INTOPT_PRINT_FLAG] = 1;
            model.intOption[INTOPT_PERMUTE_FLAG] = 1;
            model.setup(argv[2]);

            HDual solver;
            solver.solve(&model, HDUAL_VARIANT_MULTI, 8);

            model.printResult();
        }
    }
    return 0;
}

void solvePlain(const char *filename) {
    HModel model;
    model.intOption[INTOPT_PRINT_FLAG] = 1;
//    model.intOption[INTOPT_PERMUTE_FLAG] = 1;
    model.setup(filename);

    HDual solver;
    solver.solve(&model);

    model.printResult();
    model.writePivots("plain");
}

void solveTasks(const char *filename) {
    HModel model;
    model.intOption[INTOPT_PRINT_FLAG] = 1;
    model.intOption[INTOPT_PERMUTE_FLAG] = 1;
    model.setup(filename);

    HDual solver;
    solver.solve(&model, HDUAL_VARIANT_TASKS, 8);

    model.printResult();
    model.writePivots("tasks");
}

void solveMulti(const char *filename, const char *partitionfile) {
    HModel model;
    model.intOption[INTOPT_PRINT_FLAG] = 1;
    model.intOption[INTOPT_PERMUTE_FLAG] = 1;
    if (partitionfile) {
        model.strOption[STROPT_PARTITION_FILE] = partitionfile;
    }
    model.setup(filename);

    HDual solver;
    solver.solve(&model, HDUAL_VARIANT_MULTI, 8);

    model.printResult();
    model.writePivots("multi");
}
