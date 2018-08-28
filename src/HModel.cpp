#include "HModel.h"
#include "HConst.h"
#include "HTimer.h"

#include <cctype>
#include <cmath>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <map>
#include <set>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
using namespace std;

HModel::HModel() {
    intOption[INTOPT_PRINT_FLAG] = 0; // no print
    intOption[INTOPT_TRANSPOSE_FLAG] = 0; // no transpose
    intOption[INTOPT_SCALE_FLAG] = 1; // do scale
    intOption[INTOPT_TIGHT_FLAG] = 1; // do tight
    intOption[INTOPT_PERMUTE_FLAG] = 0; // no permute
    intOption[INTOPT_PERTURB_FLAG] = 1; // do perturb

    dblOption[DBLOPT_PRIMAL_TOL] = 1e-7;
    dblOption[DBLOPT_DUAL_TOL] = 1e-7;
    dblOption[DBLOPT_PAMI_CUTOFF] = 0.95;

    strOption[STROPT_PARTITION_FILE] = "";
}

void HModel::setup(const char *filename) {
    // Setup timer
    timer.reset();
    totalTime = 0;
    problemStatus = -1;
    numberIteration = 0;
    modelName = filename;

    // Load the model
    setup_loadMPS(filename);

    // Check size
    numTot = numCol + numRow;
    if (numRow == 0)
        return;

    // Setup random buffers: shuffle the break
    intBreak.resize(numTot);
    for (int i = 0; i < numTot; i++)
        intBreak[i] = i;
    for (int i = numTot - 1; i >= 1; i--) {
        int j = random.intRandom() % (i + 1);
        swap(intBreak[i], intBreak[j]);
    }
    dblXpert.resize(numTot);
    for (int i = 0; i < numTot; i++)
        dblXpert[i] = random.dblRandom();

    // Setup other part
    setup_transposeLP();
    setup_scaleMatrix();
    //    setup_tightenBound();
    //    setup_shuffleColumn();
    setup_allocWorking();

    // Initialize the values
    initCost();
    initBound();
    initValue();

    // Save the input time
    totalTime += timer.getTime();

}

bool load_mpsLine(FILE *file, int lmax, char *line, char *flag, double *data) {
    int F1 = 1, F2 = 4, F3 = 14, F4 = 24, F5 = 39, F6 = 49;

    // check the buffer
    if (flag[1]) {
        flag[1] = 0;
        memcpy(&data[2], &line[F5], 8);
        data[0] = atof(&line[F6]);
        return true;
    }

    // try to read some to the line
    for (;;) {
        // Line input
        fgets(line, lmax, file);

        // Line trim   -- to delete tailing white spaces
        int lcnt = strlen(line) - 1;
        while (isspace(line[lcnt]) && lcnt >= 0)
            lcnt--;
        if (lcnt <= 0 || line[0] == '*')
            continue;

        // Line expand -- to get data easier
        lcnt++;
        while (lcnt < F4)
            line[lcnt++] = ' '; // For row and bound row name
        if (lcnt == F4)
            line[lcnt++] = '0'; // For bound value
        line[lcnt] = '\0';

        // Done with section symbol
        if (line[0] != ' ') {
            flag[0] = line[0];
            return false;
        }

        // Read major symbol & name
        flag[0] = line[F1 + 1] == ' ' ? line[F1] : line[F1 + 1];
        memcpy(&data[1], &line[F2], 8);

        // Read 1st minor name & value to output
        memcpy(&data[2], &line[F3], 8);
        data[0] = atof(&line[F4]);

        // Keep 2nd minor name & value for future
        if (lcnt > F5)
            flag[1] = 1;
        break;
    }

    return true;
}

void HModel::setup_loadMPS(const char *filename) {
    // MPS file buffer
    numRow = 0;
    numCol = 0;
    objOffset = 0;
    FILE *file = fopen(filename, "r");
    if (file == 0) {
        return;
    }

    // Input buffer
    const int lmax = 128;
    char line[lmax];
    char flag[2] = { 0, 0 };
    double data[3];

    // Load NAME and ROWS
    load_mpsLine(file, lmax, line, flag, data);
    load_mpsLine(file, lmax, line, flag, data);

    vector<char> rowType;
    map<double, int> rowIndex;
    double objName = 0;
    while (load_mpsLine(file, lmax, line, flag, data)) {
        if (flag[0] == 'N' && objName == 0) {
            objName = data[1];
        } else {
            rowType.push_back(flag[0]);
            rowIndex[data[1]] = numRow++;
        }
    }

    // Load COLUMNS
    map<double, int> colIndex;
    double lastName = 0;
    while (load_mpsLine(file, lmax, line, flag, data)) {
        if (lastName != data[1]) { // New column
            lastName = data[1];
            colIndex[data[1]] = numCol++;
            colCost.push_back(0);
            Astart.push_back(Aindex.size());
        }
        if (data[2] == objName) // Cost
            colCost.back() = data[0];
        else if (data[0] != 0) {
            Aindex.push_back(rowIndex[data[2]]);
            Avalue.push_back(data[0]);
        }
    }
    Astart.push_back(Aindex.size());

    // Load RHS
    vector<double> RHS(numRow, 0);
    while (load_mpsLine(file, lmax, line, flag, data)) {
        if (data[2] != objName) {
            int iRow = rowIndex[data[2]];
            RHS[iRow] = data[0];
        } else {
            objOffset = data[0]; // Objective offset
        }
    }

    // Load RANGES
    rowLower.resize(numRow);
    rowUpper.resize(numRow);
    if (flag[0] == 'R') {
        while (load_mpsLine(file, lmax, line, flag, data)) {
            int iRow = rowIndex[data[2]];
            if (rowType[iRow] == 'L' || (rowType[iRow] == 'E' && data[0] < 0)) {
                rowLower[iRow] = RHS[iRow] - fabs(data[0]);
                rowUpper[iRow] = RHS[iRow];
            } else {
                rowUpper[iRow] = RHS[iRow] + fabs(data[0]);
                rowLower[iRow] = RHS[iRow];
            }
            rowType[iRow] = 'X';
        }
    }

    // Setup bounds for row without 'RANGE'
    for (int iRow = 0; iRow < numRow; iRow++) {
        switch (rowType[iRow]) {
        case 'L':
            rowLower[iRow] = -HSOL_CONST_INF;
            rowUpper[iRow] = RHS[iRow];
            break;
        case 'G':
            rowLower[iRow] = RHS[iRow];
            rowUpper[iRow] = +HSOL_CONST_INF;
            break;
        case 'E':
            rowLower[iRow] = RHS[iRow];
            rowUpper[iRow] = RHS[iRow];
            break;
        case 'N':
            rowLower[iRow] = -HSOL_CONST_INF;
            rowUpper[iRow] = +HSOL_CONST_INF;
            break;
        case 'X':
            break;
        }
    }

    // Load BOUNDS
    colLower.assign(numCol, 0);
    colUpper.assign(numCol, HSOL_CONST_INF);
    if (flag[0] == 'B') {
        while (load_mpsLine(file, lmax, line, flag, data)) {
            int iCol = colIndex[data[2]];

            switch (flag[0]) {
            case 'O': /*LO*/
                colLower[iCol] = data[0];
                break;
            case 'I': /*MI*/
                colLower[iCol] = -HSOL_CONST_INF;
                break;
            case 'L': /*PL*/
                colUpper[iCol] = HSOL_CONST_INF;
                break;
            case 'X': /*FX*/
                colLower[iCol] = data[0];
                colUpper[iCol] = data[0];
                break;
            case 'R': /*FR*/
                colLower[iCol] = -HSOL_CONST_INF;
                colUpper[iCol] = HSOL_CONST_INF;
                break;
            case 'P': /*UP*/
                colUpper[iCol] = data[0];
                if (colLower[iCol] == 0 && data[0] < 0)
                    colLower[iCol] = -HSOL_CONST_INF;
                break;
            }
        }
    }

    // Load ENDATA and close file
    fclose(file);
}

void HModel::setup_transposeLP() {
    if (intOption[INTOPT_TRANSPOSE_FLAG] == 0)
        return;

    int transposeCancelled = 0;
    if (1.0 * numCol / numRow > 0.2) {
//        cout << "transpose-cancelled-by-ratio" << endl;
        transposeCancelled = 1;
        return;
    }

    // Convert primal cost to dual bound
    const double inf = HSOL_CONST_INF;
    vector<double> dualRowLower(numCol);
    vector<double> dualRowUpper(numCol);
    for (int j = 0; j < numCol; j++) {
        double lower = colLower[j];
        double upper = colUpper[j];

        /*
         * Primal      Dual
         * Free        row = c
         * x > 0       row < c
         * x < 0       row > c
         * x = 0       row free
         * other       cancel
         */

        if (lower == -inf && upper == inf) {
            dualRowLower[j] = colCost[j];
            dualRowUpper[j] = colCost[j];
        } else if (lower == 0 && upper == inf) {
            dualRowLower[j] = -inf;
            dualRowUpper[j] = colCost[j];
        } else if (lower == -inf && upper == 0) {
            dualRowLower[j] = colCost[j];
            dualRowUpper[j] = +inf;
        } else if (lower == 0 && upper == 0) {
            dualRowLower[j] = -inf;
            dualRowUpper[j] = +inf;
        } else {
            transposeCancelled = 1;
            break;
        }
    }

    // Check flag
    if (transposeCancelled == 1) {
//        cout << "transpose-cancelled-by-column" << endl;
        return;
    }

    // Convert primal row bound to dual variable cost
    vector<double> dualColLower(numRow);
    vector<double> dualColUpper(numRow);
    vector<double> dualCost(numRow);
    for (int i = 0; i < numRow; i++) {
        double lower = rowLower[i];
        double upper = rowUpper[i];

        /*
         * Primal      Dual
         * row = b     Free
         * row < b     y < 0
         * row > b     y > 0
         * row free    y = 0
         * other       cancel
         */

        if (lower == upper) {
            dualColLower[i] = -inf;
            dualColUpper[i] = +inf;
            dualCost[i] = -lower;
        } else if (lower == -inf && upper != inf) {
            dualColLower[i] = -inf;
            dualColUpper[i] = 0;
            dualCost[i] = -upper;
        } else if (lower != -inf && upper == inf) {
            dualColLower[i] = 0;
            dualColUpper[i] = +inf;
            dualCost[i] = -lower;
        } else if (lower == -inf && upper == inf) {
            dualColLower[i] = 0;
            dualColUpper[i] = 0;
            dualCost[i] = 0;
        } else {
            transposeCancelled = 1;
            break;
        }
    }

    // Check flag
    if (transposeCancelled == 1) {
//        cout << "transpose-cancelled-by-row" << endl;
        return;
    }

    // We can now really transpose things
    vector<int> iwork(numRow, 0);
    vector<int> ARstart(numRow + 1, 0);
    int AcountX = Aindex.size();
    vector<int> ARindex(AcountX);
    vector<double> ARvalue(AcountX);
    for (int k = 0; k < AcountX; k++)
        iwork[Aindex[k]]++;
    for (int i = 1; i <= numRow; i++)
        ARstart[i] = ARstart[i - 1] + iwork[i - 1];
    for (int i = 0; i < numRow; i++)
        iwork[i] = ARstart[i];
    for (int iCol = 0; iCol < numCol; iCol++) {
        for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
            int iRow = Aindex[k];
            int iPut = iwork[iRow]++;
            ARindex[iPut] = iCol;
            ARvalue[iPut] = Avalue[k];
        }
    }

    // Transpose the problem!
    swap(numRow, numCol);
    Astart.swap(ARstart);
    Aindex.swap(ARindex);
    Avalue.swap(ARvalue);
    colLower.swap(dualColLower);
    colUpper.swap(dualColUpper);
    rowLower.swap(dualRowLower);
    rowUpper.swap(dualRowUpper);
    colCost.swap(dualCost);
//    cout << "problem-transposed" << endl;
}

void HModel::setup_scaleMatrix() {
    if (intOption[INTOPT_SCALE_FLAG] == 0)
        return;

    // Reset all scaling to 1
    vector<double> colScale(numCol, 1);
    vector<double> rowScale(numRow, 1);

    // Find out min0 / max0, skip on if in [0.2, 5]
    const double inf = HSOL_CONST_INF;
    double min0 = inf, max0 = 0;
    for (int k = 0, AnX = Astart[numCol]; k < AnX; k++) {
        double value = fabs(Avalue[k]);
        min0 = min(min0, value);
        max0 = max(max0, value);
    }
    if (min0 >= 0.2 && max0 <= 5)
        return;

    // See if we want to include cost include if min-cost < 0.1
    double minc = inf;
    for (int i = 0; i < numCol; i++)
        if (colCost[i])
            minc = min(minc, fabs(colCost[i]));
    bool doCost = minc < 0.1;

    // Search up to 6 times
    vector<double> rowMin(numRow, inf);
    vector<double> rowMax(numRow, 1 / inf);
    for (int search_count = 0; search_count < 6; search_count++) {
        // Find column scale, prepare row data
        for (int iCol = 0; iCol < numCol; iCol++) {
            // For column scale (find)
            double colMin = inf;
            double colMax = 1 / inf;
            double myCost = fabs(colCost[iCol]);
            if (doCost && myCost != 0)
                colMin = min(colMin, myCost), colMax = max(colMax, myCost);
            for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
                double value = fabs(Avalue[k]) * rowScale[Aindex[k]];
                colMin = min(colMin, value), colMax = max(colMax, value);
            }
            colScale[iCol] = 1 / sqrt(colMin * colMax);

            // For row scale (only collect)
            for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
                int iRow = Aindex[k];
                double value = fabs(Avalue[k]) * colScale[iCol];
                rowMin[iRow] = min(rowMin[iRow], value);
                rowMax[iRow] = max(rowMax[iRow], value);
            }
        }

        // For row scale (find)
        for (int iRow = 0; iRow < numRow; iRow++)
            rowScale[iRow] = 1 / sqrt(rowMin[iRow] * rowMax[iRow]);
        rowMin.assign(numRow, inf);
        rowMax.assign(numRow, 1 / inf);
    }

    // Make it numerical better
    const double ln2 = log(2.0);
    for (int iCol = 0; iCol < numCol; iCol++)
        colScale[iCol] = pow(2.0, floor(log(colScale[iCol]) / ln2 + 0.5));
    for (int iRow = 0; iRow < numRow; iRow++)
        rowScale[iRow] = pow(2.0, floor(log(rowScale[iRow]) / ln2 + 0.5));

    // Apply scaling to matrix and bounds
    for (int iCol = 0; iCol < numCol; iCol++)
        for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++)
            Avalue[k] *= (colScale[iCol] * rowScale[Aindex[k]]);

    for (int iCol = 0; iCol < numCol; iCol++) {
        colLower[iCol] /= colLower[iCol] == -inf ? 1 : colScale[iCol];
        colUpper[iCol] /= colUpper[iCol] == +inf ? 1 : colScale[iCol];
        colCost[iCol] *= colScale[iCol];
    }
    for (int iRow = 0; iRow < numRow; iRow++) {
        rowLower[iRow] *= rowLower[iRow] == -inf ? 1 : rowScale[iRow];
        rowUpper[iRow] *= rowUpper[iRow] == +inf ? 1 : rowScale[iRow];
    }
}

void HModel::setup_tightenBound() {
    if (intOption[INTOPT_TIGHT_FLAG] == 0)
        return;

    // Make a AR copy
    vector<int> iwork(numRow, 0);
    vector<int> ARstart(numRow + 1, 0);
    int AcountX = Aindex.size();
    vector<int> ARindex(AcountX);
    vector<double> ARvalue(AcountX);
    for (int k = 0; k < AcountX; k++)
        iwork[Aindex[k]]++;
    for (int i = 1; i <= numRow; i++)
        ARstart[i] = ARstart[i - 1] + iwork[i - 1];
    for (int i = 0; i < numRow; i++)
        iwork[i] = ARstart[i];
    for (int iCol = 0; iCol < numCol; iCol++) {
        for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
            int iRow = Aindex[k];
            int iPut = iwork[iRow]++;
            ARindex[iPut] = iCol;
            ARvalue[iPut] = Avalue[k];
        }
    }

    // Save column bounds
    vector<double> colLower0 = colLower;
    vector<double> colUpper0 = colUpper;

    double big_B = 1e10;
    int iPass = 0;
    for (;;) {
        int numberChanged = 0;
        for (int iRow = 0; iRow < numRow; iRow++) {
            // SKIP free rows
            if (rowLower[iRow] < -big_B && rowUpper[iRow] > big_B)
                continue;

            // possible row
            int ninfU = 0;
            int ninfL = 0;
            double xmaxU = 0.0;
            double xminL = 0.0;
            int myStart = ARstart[iRow];
            int myEnd = ARstart[iRow + 1];
            // Compute possible lower and upper ranges

            for (int k = myStart; k < myEnd; ++k) {
                int iCol = ARindex[k];
                double value = ARvalue[k];
                double upper = value > 0 ? colUpper[iCol] : -colLower[iCol];
                double lower = value > 0 ? colLower[iCol] : -colUpper[iCol];
                value = fabs(value);
                if (upper < big_B)
                    xmaxU += upper * value;
                else
                    ++ninfU;
                if (lower > -big_B)
                    xminL += lower * value;
                else
                    ++ninfL;
            }

            // Build in a margin of error
            xmaxU += 1.0e-8 * fabs(xmaxU);
            xminL -= 1.0e-8 * fabs(xminL);

            double xminLmargin =
                    (fabs(xminL) > 1.0e8) ? 1e-12 * fabs(xminL) : 0;
            double xmaxUmargin =
                    (fabs(xmaxU) > 1.0e8) ? 1e-12 * fabs(xmaxU) : 0;

            // Skip redundant row : also need to consider U < L  case
            double comp_U = xmaxU + ninfU * 1.0e31;
            double comp_L = xminL - ninfL * 1.0e31;
            if (comp_U <= rowUpper[iRow] + 1e-7
                    && comp_L >= rowLower[iRow] - 1e-7)
                continue;

            double row_L = rowLower[iRow];
            double row_U = rowUpper[iRow];

            // Now see if we can tighten column bounds
            for (int k = myStart; k < myEnd; ++k) {
                double value = ARvalue[k];
                int iCol = ARindex[k];
                double col_L = colLower[iCol];
                double col_U = colUpper[iCol];
                double new_L = -HSOL_CONST_INF;
                double new_U = +HSOL_CONST_INF;

                if (value > 0.0) {
                    if (row_L > -big_B && ninfU <= 1
                            && (ninfU == 0 || col_U > +big_B))
                        new_L = (row_L - xmaxU) / value + (1 - ninfU) * col_U
                                - xmaxUmargin;
                    if (row_U < +big_B && ninfL <= 1
                            && (ninfL == 0 || col_L < -big_B))
                        new_U = (row_U - xminL) / value + (1 - ninfL) * col_L
                                + xminLmargin;
                } else {
                    if (row_L > -big_B && ninfU <= 1
                            && (ninfU == 0 || col_L < -big_B))
                        new_U = (row_L - xmaxU) / value + (1 - ninfU) * col_L
                                + xmaxUmargin;
                    if (row_U < +big_B && ninfL <= 1
                            && (ninfL == 0 || col_U > +big_B))
                        new_L = (row_U - xminL) / value + (1 - ninfL) * col_U
                                - xminLmargin;
                }

                if (new_U < col_U - 1.0e-12 && new_U < big_B) {
                    colUpper[iCol] = max(new_U, col_L);
                    numberChanged++;
                }
                if (new_L > col_L + 1.0e-12 && new_L > -big_B) {
                    colLower[iCol] = min(new_L, col_U);
                    numberChanged++;
                }
            }
        }

        if (numberChanged == 0)
            break;
        iPass++;
        if (iPass > 10)
            break;

    }

    double useTolerance = 1.0e-3;
    for (int iCol = 0; iCol < numCol; iCol++) {
        if (colUpper0[iCol] > colLower0[iCol] + useTolerance) {
            const double relax = 100.0 * useTolerance;
            if (colUpper[iCol] - colLower[iCol] < useTolerance + 1.0e-8) {
                colLower[iCol] = max(colLower0[iCol], colLower[iCol] - relax);
                colUpper[iCol] = min(colUpper0[iCol], colUpper[iCol] + relax);
            } else {
                if (colUpper[iCol] < colUpper0[iCol]) {
                    colUpper[iCol] = min(colUpper[iCol] + relax,
                            colUpper0[iCol]);
                }
                if (colLower[iCol] > colLower0[iCol]) {
                    colLower[iCol] = min(colLower[iCol] - relax,
                            colLower0[iCol]);
                }
            }
        }
    }
}

void HModel::setup_shuffleColumn() {
    if (intOption[INTOPT_PERMUTE_FLAG] == 0)
        return;

    // 1. Shuffle the column index
    HRandom localRandom;
    for (int i = 0; i < 10; i++)
        localRandom.intRandom();
    vector<int> iFrom(numCol);
    for (int i = 0; i < numCol; i++)
        iFrom[i] = i;
    for (int i = numCol - 1; i >= 1; i--) {
        int j = localRandom.intRandom() % (i + 1);
        swap(iFrom[i], iFrom[j]);
    }

    // 2. Save original copy
    vector<int> start = Astart;
    vector<int> index = Aindex;
    vector<double> value = Avalue;
    vector<double> lower = colLower;
    vector<double> upper = colUpper;
    vector<double> xcost = colCost;
    vector<int> ibreak = intBreak;
    vector<double> dxpert = dblXpert;

    // 3. Generate the permuted matrix
    int countX = 0;
    for (int i = 0; i < numCol; i++) {
        int ifrom = iFrom[i];
        Astart[i] = countX;
        for (int k = start[ifrom]; k < start[ifrom + 1]; k++) {
            Aindex[countX] = index[k];
            Avalue[countX] = value[k];
            countX++;
        }
        colLower[i] = lower[ifrom];
        colUpper[i] = upper[ifrom];
        colCost[i] = xcost[ifrom];
        intBreak[i] = ibreak[ifrom];
        dblXpert[i] = dxpert[ifrom];
    }
    assert(Astart[numCol] == countX);
}

void HModel::setup_allocWorking() {
    // Setup starting base
    basicIndex.resize(numRow);
    for (int iRow = 0; iRow < numRow; iRow++)
        basicIndex[iRow] = iRow + numCol;
    nonbasicFlag.assign(numTot, 0);
    nonbasicMove.resize(numTot);
    for (int i = 0; i < numCol; i++)
        nonbasicFlag[i] = 1;

    // Matrix, factor
    matrix.setup(numCol, numRow, &Astart[0], &Aindex[0], &Avalue[0]);
    factor.setup(numCol, numRow, &Astart[0], &Aindex[0], &Avalue[0],
            &basicIndex[0]);
    limitUpdate = 5000;

    // Setup other buffer
    buffer.setup(numRow);
    bufferLong.setup(numCol);

    // Setup bounds and solution spaces
    workCost.resize(numTot);
    workDual.resize(numTot);
    workShift.assign(numTot, 0);

    workLower.resize(numTot);
    workUpper.resize(numTot);
    workRange.resize(numTot);
    workValue.resize(numTot);

    baseLower.resize(numRow);
    baseUpper.resize(numRow);
    baseValue.resize(numRow);
}

void HModel::initCost(int perturb) {
    // Copy the cost
    for (int i = 0; i < numCol; i++)
        workCost[i] = colCost[i];
    for (int i = numCol; i < numTot; i++)
        workCost[i] = 0;
    workShift.assign(numTot, 0);

    // See if we want to skip perturbation
    problemPerturbed = 0;
    if (perturb == 0 || intOption[INTOPT_PERTURB_FLAG] == 0)
        return;
    problemPerturbed = 1;

    // Perturb the original costs, scale down if is too big
    double bigc = 0;
    for (int i = 0; i < numCol; i++)
        bigc = max(bigc, fabs(workCost[i]));
    if (bigc > 100)
        bigc = sqrt(sqrt(bigc));

    // If there's few boxed variables, we will just use Simple perturbation
    double boxedRate = 0;
    for (int i = 0; i < numTot; i++)
        boxedRate += (workRange[i] < 1e30);
    boxedRate /= numTot;
    if (boxedRate < 0.01)
        bigc = min(bigc, 1.0);
    if (bigc < 1) {
//        bigc = sqrt(bigc);
    }

    // Determine the perturbation base
    double base = 5e-7 * bigc;

    // Now do the perturbation
    for (int i = 0; i < numCol; i++) {
        double lower = colLower[i];
        double upper = colUpper[i];
        double xpert = (fabs(workCost[i]) + 1) * base * (1 + dblXpert[i]);
        if (lower == -HSOL_CONST_INF && upper == HSOL_CONST_INF) {
            // Free - no perturb
        } else if (upper == HSOL_CONST_INF) {  // Lower
            workCost[i] += xpert;
        } else if (lower == -HSOL_CONST_INF) { // Upper
            workCost[i] += -xpert;
        } else if (lower != upper) {  // Boxed
            workCost[i] += (workCost[i] >= 0) ? xpert : -xpert;
        } else {
            // Fixed - no perturb
        }
    }

    for (int i = numCol; i < numTot; i++) {
        workCost[i] += (0.5 - dblXpert[i]) * 1e-12;
    }
}

void HModel::initBound(int phase) {
    // Copy bounds
    for (int i = 0; i < numCol; i++) {
        workLower[i] = colLower[i];
        workUpper[i] = colUpper[i];
    }
    for (int i = 0, j = numCol; i < numRow; i++, j++) {
        workLower[j] = -rowUpper[i];
        workUpper[j] = -rowLower[i];
    }

    // Change to dual phase 1 bound
    if (phase == 1) {
        const double inf = HSOL_CONST_INF;
        for (int i = 0; i < numTot; i++) {
            if (workLower[i] == -inf && workUpper[i] == inf) {
                // Won't change for row variables: they should never
                // Become non basic
                if (i >= numCol)
                    continue;
                workLower[i] = -1000, workUpper[i] = 1000;  // FREE
            } else if (workLower[i] == -inf) {
                workLower[i] = -1, workUpper[i] = 0;        // UPPER
            } else if (workUpper[i] == inf) {
                workLower[i] = 0, workUpper[i] = 1;         // LOWER
            } else {
                workLower[i] = 0, workUpper[i] = 0;         // BOXED or FIXED
            }
        }
    }

    // Work out ranges
    for (int i = 0; i < numTot; i++)
        workRange[i] = workUpper[i] - workLower[i];
}

void HModel::initValue() {
    for (int i = 0; i < numTot; i++) {
        if (nonbasicFlag[i]) {
            if (workLower[i] == workUpper[i]) {
                // Fixed
                workValue[i] = workLower[i];
                nonbasicMove[i] = 0;
            } else if (workLower[i] != -HSOL_CONST_INF) {
                // Lower
                workValue[i] = workLower[i];
                nonbasicMove[i] = 1;
            } else if (workUpper[i] != +HSOL_CONST_INF) {
                // UPPER
                workValue[i] = workUpper[i];
                nonbasicMove[i] = -1;
            } else {
                // FREE
                workValue[i] = 0;
                nonbasicMove[i] = 0;
            }
        } else {
            nonbasicMove[i] = 0;
        }
    }
}

void HModel::computeFactor() {
    try {
        factor.build();
    } catch (runtime_error& error) {
        cout << error.what() << endl;
        problemStatus = 3;
        // TODO Change to return
        writePivots("failed");
        exit(0);
    }
    countUpdate = 0;
}

void HModel::computeDual() {
    buffer.clear();
    for (int iRow = 0; iRow < numRow; iRow++) {
        buffer.index[iRow] = iRow;
        buffer.array[iRow] = workCost[basicIndex[iRow]]
                + workShift[basicIndex[iRow]];
    }
    buffer.count = numRow;
    factor.btran(buffer, 1);

    bufferLong.clear();
    matrix.price_by_col(bufferLong, buffer);
    for (int i = 0; i < numCol; i++)
        workDual[i] = workCost[i] - bufferLong.array[i];
    for (int i = numCol; i < numTot; i++)
        workDual[i] = workCost[i] - buffer.array[i - numCol];
}

void HModel::computeDualInfeasInDual(int *dualInfeasCount) {
    int workCount = 0;
    const double inf = HSOL_CONST_INF;
    const double tau_d = dblOption[DBLOPT_DUAL_TOL];
    for (int i = 0; i < numTot; i++) {
        // Only for non basic variables
        if (!nonbasicFlag[i])
            continue;
        // Free
        if (workLower[i] == -inf && workUpper[i] == inf)
            workCount += (fabs(workDual[i]) >= tau_d);
        // In dual, assuming that boxed variables will be flipped
        if (workLower[i] == -inf || workUpper[i] == inf)
            workCount += (nonbasicMove[i] * workDual[i] <= -tau_d);
    }
    *dualInfeasCount = workCount;
}

void HModel::computeDualInfeasInPrimal(int *dualInfeasCount) {
    int workCount = 0;
    const double inf = HSOL_CONST_INF;
    const double tau_d = dblOption[DBLOPT_DUAL_TOL];
    for (int i = 0; i < numTot; i++) {
        // Only for non basic variables
        if (!nonbasicFlag[i])
            continue;
        // Free
        if (workLower[i] == -inf && workUpper[i] == inf)
            workCount += (fabs(workDual[i]) >= tau_d);
        // In primal don't assume flip
        workCount += (nonbasicMove[i] * workDual[i] <= -tau_d);
    }
    *dualInfeasCount = workCount;
}

void HModel::correctDual(int *freeInfeasCount) {
    const double tau_d = dblOption[DBLOPT_DUAL_TOL];
    const double inf = HSOL_CONST_INF;
    int workCount = 0;
    for (int i = 0; i < numTot; i++) {
        if (nonbasicFlag[i]) {
            if (workLower[i] == -inf && workUpper[i] == inf) {
                // FREE variable
                workCount += (fabs(workDual[i]) >= tau_d);
            } else if (nonbasicMove[i] * workDual[i] <= -tau_d) {
                if (workLower[i] != -inf && workUpper[i] != inf) {
                    // Boxed variable = flip
                    flipBound(i);
                } else {
                    // Other variable = shift
                    problemPerturbed = 1;
                    if (nonbasicMove[i] == 1) {
                        double dual = (1 + random.dblRandom()) * tau_d;
                        double shift = dual - workDual[i];
                        workDual[i] = dual;
                        workCost[i] = workCost[i] + shift;
                    } else {
                        double dual = -(1 + random.dblRandom()) * tau_d;
                        double shift = dual - workDual[i];
                        workDual[i] = dual;
                        workCost[i] = workCost[i] + shift;
                    }
                }
            }
        }
    }

    *freeInfeasCount = workCount;
}

void HModel::computePrimal() {
    buffer.clear();
    for (int i = 0; i < numTot; i++)
        if (nonbasicFlag[i] && workValue[i] != 0)
            matrix.collect_aj(buffer, i, workValue[i]);
    factor.ftran(buffer, 1);
    for (int i = 0; i < numRow; i++) {
        int iCol = basicIndex[i];
        baseValue[i] = -buffer.array[i];
        baseLower[i] = workLower[iCol];
        baseUpper[i] = workUpper[iCol];
    }
}

void HModel::computeObject(int phase) {
    objective = 0;
    for (int i = 0; i < numTot; i++)
        if (nonbasicFlag[i])
            objective += workValue[i] * workDual[i];
    if (phase != 1)
        objective -= objOffset;
}

void HModel::shiftCost(int iCol, double amount) {
    problemPerturbed = 1;
    assert(workShift[iCol] == 0);
    workShift[iCol] = amount;
}

void HModel::shiftBack(int iCol) {
    workDual[iCol] -= workShift[iCol];
    workShift[iCol] = 0;
}

void HModel::flipBound(int iCol) {
    const int move = nonbasicMove[iCol] = -nonbasicMove[iCol];
    workValue[iCol] = move == 1 ? workLower[iCol] : workUpper[iCol];
}

void HModel::updateFactor(HVector *column, HVector *row_ep, int *iRow,
        int *hint) {
    timer.recordStart(HTICK_UPDATE_FACTOR);
    factor.update(column, row_ep, iRow, hint);
    if (countUpdate >= limitUpdate)
        *hint = 1;
    timer.recordFinish(HTICK_UPDATE_FACTOR);
}

void HModel::updateMatrix(int columnIn, int columnOut) {
    timer.recordStart(HTICK_UPDATE_FACTOR);
    matrix.update(columnIn, columnOut);
    timer.recordFinish(HTICK_UPDATE_FACTOR);
}

void HModel::updatePivots(int columnIn, int rowOut, int sourceOut) {
    int columnOut = basicIndex[rowOut];

    // Incoming variable
    basicIndex[rowOut] = columnIn;
    nonbasicFlag[columnIn] = 0;
    nonbasicMove[columnIn] = 0;
    baseLower[rowOut] = workLower[columnIn];
    baseUpper[rowOut] = workUpper[columnIn];

    // Outgoing variable
    nonbasicFlag[columnOut] = 1;
    if (workLower[columnOut] == workUpper[columnOut]) {
        workValue[columnOut] = workLower[columnOut];
        nonbasicMove[columnOut] = 0;
    } else if (sourceOut == -1) {
        workValue[columnOut] = workLower[columnOut];
        nonbasicMove[columnOut] = 1;
    } else {
        workValue[columnOut] = workUpper[columnOut];
        nonbasicMove[columnOut] = -1;
    }

    countUpdate++;
}

void HModel::changeUpdate(int updateMethod) {
    factor.change(updateMethod);
}

void HModel::reportPivots(int columnIn, int columnOut, double alpha) {
    if (columnIn >= 0)
        numberIteration++;
    historyColumnIn.push_back(columnIn);
    historyColumnOut.push_back(columnOut);
    historyAlpha.push_back(alpha);
}

void HModel::reportStatus(int status) {
    problemStatus = status;
}

void HModel::printMessage(const char *message) {
    if (intOption[INTOPT_PRINT_FLAG])
        printf("%s\n", message);
}

void HModel::printObject() {
    if (intOption[INTOPT_PRINT_FLAG])
        printf("%10d  %20.10e\n", numberIteration, objective);
}

void HModel::printResult() {
    if (problemStatus == 0) {
        printf("HsolRpIt - OPTIMAL %16s %20.10e %10d %10.3f\n", modelName.c_str(),
                objective, numberIteration, totalTime);
    } else {
        printf("HsolRpIt - NOT-OPT status = %d\n", problemStatus);
    }
}

void HModel::printProgress() {
//    return;
//    static double nextReport = 0;
//    double currentTime = timer.getTime();
//    if (currentTime >= nextReport) {
//        computeObject();
//        printf("PROGRESS %16s %20.10e %10d %10.3f\n", modelName.c_str(),
//                objective, numberIteration, timer.getTime());
//        if (currentTime < 50) {
//            nextReport = ((int) (5 * currentTime + 1)) / 5.0 - 0.00001;
//        } else if (currentTime < 500) {
//            nextReport = ((int) (currentTime + 1)) - 0.00001;
//        } else {
//            nextReport = ((int) (0.2 * currentTime + 1)) / 0.2 - 0.00001;
//        }
//    }
}

void HModel::writePivots(const char *suffix) {
    string filename = "z-" + modelName + "-" + suffix;
    ofstream output(filename.c_str());
    int count = historyColumnIn.size();
    output << modelName << " " << count << "\t" << totalTime << endl;
    output << setprecision(12);
    for (int i = 0; i < count; i++) {
        output << historyColumnIn[i] << "\t";
        output << historyColumnOut[i] << "\t";
        output << historyAlpha[i] << endl;
    }
    output.close();
}

void HModel::writeMPS(const char *filename) {

}
