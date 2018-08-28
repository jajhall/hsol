#ifndef HDUAL_H_
#define HDUAL_H_

#include "HModel.h"
#include "HMatrix.h"
#include "HFactor.h"
#include "HVector.h"
#include "HDualRow.h"
#include "HDualRHS.h"

#include <set>
#include <string>
#include <vector>
using namespace std;

const int HSOL_THREAD_LIMIT = 32;
const int HSOL_SLICED_LIMIT = 100;

enum HDUAL_VARIANT {
    HDUAL_VARIANT_PLAIN = 0, HDUAL_VARIANT_TASKS, HDUAL_VARIANT_MULTI,
};

class HDual {
public:
    void solve(HModel *model, int variant = 0, int num_threads = 1);
public:
    void init(int num_threads);
    void init_slice(int init_sliced_num);

    void solve_phase1();
    void solve_phase2();

    void rebuild();
    void cleanup();

    void iterate();
    void iterate_tasks();
    void iterate_multi();

    void chooseRow();

    void chooseColumn(HVector *row_ep);
    void chooseColumn_slice(HVector *row_ep);

    void updateFtranBFRT();
    void updateFtran();
    void updateFtranDSE(HVector *dseVector);
    void updateVerify();
    void updateDual();
    void updatePrimal(HVector *dseVector);
    void updatePivots();

    void major_chooseRow();
    void major_chooseRowBtran();
    void minor_chooseRow();

    void minor_update();
    void minor_updateDual();
    void minor_updatePrimal();
    void minor_updatePivots();
    void minor_updateRows();

    void major_update();
    void major_updateFtranPrepare();
    void major_updateFtranParallel();
    void major_updateFtranFinal();
    void major_updatePrimal();
    void major_updateFactor();

    void major_rollback();

    // Variant choice
    int dual_variant;

    // Model
    HModel *model;
    double Tp; // Tolerance for primal
    double Td; // Tolerance for dual

    int numCol;
    int numRow;
    int numTot;
    const HMatrix *matrix;
    const HFactor *factor;

    const int *jMove;
    const double *workRange;
    const double *baseLower;
    const double *baseUpper;
    double *baseValue;
    double *workDual;

    int solvePhase;
    int invertHint;

    HVector row_ep;
    HVector row_ap;
    HVector column;
    HVector columnBFRT;
    HVector columnDSE;
    double row_epDensity;
    double rowdseDensity;
    double columnDensity;

    HDualRow dualRow;

    // Solving related buffers
    int dualInfeasCount;

    HDualRHS dualRHS;

    // Simplex pivotal information
    int rowOut;
    int columnOut;
    int sourceOut; // -1 from small to lower, +1 to upper
    int columnIn;
    double deltaPrimal;
    double thetaDual;
    double thetaPrimal;
    double alpha;
    double alphaRow;

    // Partitioned coefficient matrix
    int slice_num;
    int slice_PRICE;
    int slice_start[HSOL_SLICED_LIMIT + 1];
    HMatrix slice_matrix[HSOL_SLICED_LIMIT];
    HVector slice_row_ap[HSOL_SLICED_LIMIT];
    HDualRow slice_dualRow[HSOL_SLICED_LIMIT];

    // Multiple price data
    struct MChoice {
        int rowOut;
        double baseValue;
        double baseLower;
        double baseUpper;
        double infeasValue;
        double infeasDevex;
        double infeasLimit;
        HVector row_ep;
        HVector column;
        HVector columnBFRT;
    };

    struct MFinish {
        int moveIn;
        double shiftOut;
        vector<int> flipList;

        int rowOut;
        int columnOut;
        int columnIn;
        double alphaRow;
        double thetaPrimal;
        double basicBound;
        double basicValue;
        double devex;
        HVector_ptr row_ep;
        HVector_ptr column;
        HVector_ptr columnBFRT;
    };

    int multi_num;
    int multi_iChoice;
    int multi_nFinish;
    int multi_iteration;
    int multi_chooseAgain;
    MChoice multi_choice[HSOL_THREAD_LIMIT];
    MFinish multi_finish[HSOL_THREAD_LIMIT];

    double total_fake;
    double total_INVERT_TICK;
    double total_FT_inc_TICK;
};

#endif /* HDUAL_H_ */
