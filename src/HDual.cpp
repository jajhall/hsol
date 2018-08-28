#include "HDual.h"
#include "HConst.h"
#include "HTimer.h"
#include "HPrimal.h"

#include <iostream>
#include <cstdio>
#include <cmath>
#include <set>
#include <stdexcept>
using namespace std;

void HDual::solve(HModel *ptr_model, int variant, int num_threads) {
    dual_variant = variant;
    model = ptr_model;
    if (model->getNumRow() == 0)
        return;
    model->timer.reset();

    // Initialise working environment
    init(num_threads);

    model->initCost(1);
    model->computeFactor();
    model->computeDual();
    model->computeDualInfeasInDual(&dualInfeasCount);
    solvePhase = dualInfeasCount > 0 ? 1 : 2;

    // Find largest dual and adjust the dual tolerance accordingly
    double largeDual = 0;
    for (int i = 0; i < numTot; i++) {
        if (model->getNonbasicFlag()[i]) {
            double myDual = fabs(workDual[i] * jMove[i]);
            if (largeDual < myDual)
                largeDual = myDual;
        }
    }

    // The major solving loop
    while (solvePhase) {
        switch (solvePhase) {
        case 1:
            solve_phase1();
            break;
        case 2:
            solve_phase2();
            break;
        case 4:
            break;
        default:
            solvePhase = 0;
            break;
        }
        // Jump for primal
        if (solvePhase == 4)
            break;
    }

    // Report the ticks before primal
    if (dual_variant == HDUAL_VARIANT_PLAIN) {
        int reportList[] = { HTICK_INVERT, HTICK_CHUZR1, HTICK_BTRAN,
                HTICK_PRICE, HTICK_CHUZC1, HTICK_CHUZC2, HTICK_CHUZC3,
                HTICK_FTRAN, HTICK_FTRAN_MIX, HTICK_FTRAN_DSE,
                HTICK_UPDATE_DUAL, HTICK_UPDATE_PRIMAL, HTICK_UPDATE_WEIGHT,
                HTICK_UPDATE_FACTOR };
        int reportCount = sizeof(reportList) / sizeof(int);
        model->timer.report(reportCount, reportList);
    }

    if (dual_variant == HDUAL_VARIANT_TASKS) {
        int reportList[] = { HTICK_INVERT, HTICK_CHUZR1, HTICK_BTRAN,
                HTICK_PRICE, HTICK_CHUZC1, HTICK_CHUZC2, HTICK_CHUZC3,
                HTICK_FTRAN, HTICK_FTRAN_MIX, HTICK_FTRAN_DSE,
                HTICK_UPDATE_DUAL, HTICK_UPDATE_PRIMAL, HTICK_UPDATE_WEIGHT,
                HTICK_UPDATE_FACTOR, HTICK_GROUP1, HTICK_GROUP2 };
        int reportCount = sizeof(reportList) / sizeof(int);
        model->timer.report(reportCount, reportList);
    }

    if (dual_variant == HDUAL_VARIANT_MULTI) {
        int reportList[] = { HTICK_INVERT, HTICK_CHUZR1, HTICK_BTRAN,
                HTICK_PRICE, HTICK_CHUZC1, HTICK_CHUZC2, HTICK_CHUZC3,
                HTICK_FTRAN, HTICK_FTRAN_MIX, HTICK_FTRAN_DSE,
                HTICK_UPDATE_DUAL, HTICK_UPDATE_PRIMAL, HTICK_UPDATE_WEIGHT,
                HTICK_UPDATE_FACTOR, HTICK_UPDATE_ROW_EP };
        int reportCount = sizeof(reportList) / sizeof(int);
        model->timer.report(reportCount, reportList);
        printf("PAMI   %-20s    CUTOFF  %6g    PERSISTENSE  %6g\n",
                model->modelName.c_str(), model->dblOption[DBLOPT_PAMI_CUTOFF],
                model->numberIteration / (1.0 + multi_iteration));
    }

    // Use primal to clean up
    if (solvePhase == 4) {
        HPrimal hPrimal;
        hPrimal.solvePhase2(model);
    }

    // Save the solved results
    model->totalTime += model->timer.getTime();
}

void HDual::init(int num_threads) {
    // Copy size, matrix and factor
    numCol = model->getNumCol();
    numRow = model->getNumRow();
    numTot = model->getNumTot();
    matrix = model->getMatrix();
    factor = model->getFactor();

    // Copy pointers
    jMove = model->getNonbasicMove();
    workDual = model->getWorkDual();
    workRange = model->getWorkRange();
    baseLower = model->getBaseLower();
    baseUpper = model->getBaseUpper();
    baseValue = model->getBaseValue();

    // Copy tolerances
    Tp = model->dblOption[DBLOPT_PRIMAL_TOL];
    Td = model->dblOption[DBLOPT_DUAL_TOL];

    // Setup local vectors
    columnDSE.setup(numRow);
    columnBFRT.setup(numRow);
    column.setup(numRow);
    row_ep.setup(numRow);
    row_ap.setup(numCol);
    columnDensity = 0;
    row_epDensity = 0;
    rowdseDensity = 0;

    // Setup other buffers
    dualRow.setup(model);
    dualRHS.setup(model);

    // Initialize for tasks
    if (dual_variant == HDUAL_VARIANT_TASKS) {
        init_slice(num_threads - 2);
    }

    // Initialize for multi
    if (dual_variant == HDUAL_VARIANT_MULTI) {
        multi_num = num_threads;
        if (multi_num < 1)
            multi_num = 1;
        if (multi_num > HSOL_THREAD_LIMIT)
            multi_num = HSOL_THREAD_LIMIT;
        for (int i = 0; i < multi_num; i++) {
            multi_choice[i].row_ep.setup(numRow);
            multi_choice[i].column.setup(numRow);
            multi_choice[i].columnBFRT.setup(numRow);
        }
        init_slice(multi_num - 1);
    }
    multi_iteration = 0;
    string partitionFile = model->strOption[STROPT_PARTITION_FILE];
    if (partitionFile.size()) {
        dualRHS.setup_partition(partitionFile.c_str());
    }
}

void HDual::init_slice(int init_sliced_num) {
    // Number of slices
    slice_num = init_sliced_num;
    if (slice_num < 1)
        slice_num = 1;
    if (slice_num > HSOL_SLICED_LIMIT)
        slice_num = HSOL_SLICED_LIMIT;

    // Alias to the matrix
    const int *Astart = matrix->getAstart();
    const int *Aindex = matrix->getAindex();
    const double *Avalue = matrix->getAvalue();
    const int AcountX = Astart[numCol];

    // Figure out partition weight
    double sliced_countX = AcountX / slice_num;
    slice_start[0] = 0;
    for (int i = 0; i < slice_num - 1; i++) {
        int endColumn = slice_start[i] + 1; // At least one column
        int endX = Astart[endColumn];
        int stopX = (i + 1) * sliced_countX;
        while (endX < stopX) {
            endX = Astart[++endColumn];
        }
        slice_start[i + 1] = endColumn;
        if (endColumn >= numCol) {
            slice_num = i; // SHRINK
            break;
        }
    }
    slice_start[slice_num] = numCol;

    // Partition the matrix, row_ap and related packet
    vector<int> sliced_Astart;
    for (int i = 0; i < slice_num; i++) {
        // The matrix
        int mystart = slice_start[i];
        int mycount = slice_start[i + 1] - mystart;
        int mystartX = Astart[mystart];
        sliced_Astart.resize(mycount + 1);
        for (int k = 0; k <= mycount; k++)
            sliced_Astart[k] = Astart[k + mystart] - mystartX;
        slice_matrix[i].setup(mycount, numRow, &sliced_Astart[0],
                Aindex + mystartX, Avalue + mystartX);

        // The row_ap and its packages
        slice_row_ap[i].setup(mycount);
        slice_dualRow[i].setupSlice(model, mycount);
    }
}

void HDual::solve_phase1() {
    model->printMessage("dual-phase-1-start");
    // Switch to dual phase 1 bounds
    model->initBound(1);
    model->initValue();

    // Main solving structure
    for (;;) {
        rebuild();
        for (;;) {
            switch (dual_variant) {
            default:
            case HDUAL_VARIANT_PLAIN:
                iterate();
                break;
            case HDUAL_VARIANT_TASKS:
                iterate_tasks();
                break;
            case HDUAL_VARIANT_MULTI:
                iterate_multi();
                break;
            }
            if (invertHint)
                break;
        }
        if (model->countUpdate == 0)
            break;
    }

    if (rowOut == -1) {
        model->printMessage("dual-phase-1-optimal");
        // Go to phase 2
        if (model->objective == 0) {
            solvePhase = 2;
        } else {
            // We still have dual infeasible
            if (model->problemPerturbed) {
                // Clean up perturbation and go on
                cleanup();
                if (dualInfeasCount == 0)
                    solvePhase = 2;
            } else {
                // Report dual infeasible
                solvePhase = -1;
                model->printMessage("dual-infeasible");
                model->reportStatus(LP_Status_Unbounded);
            }
        }
    } else if (columnIn == -1) {
        // We got dual phase 1 unbounded - strange
        model->printMessage("dual-phase-1-unbounded");
        if (model->problemPerturbed) {
            // Clean up perturbation and go on
            cleanup();
            if (dualInfeasCount == 0)
                solvePhase = 2;
        } else {
            // Report strange issues
            solvePhase = -1;
            model->printMessage("dual-phase-1-not-solved");
            model->reportStatus(LP_Status_Failed);
        }
    }

    if (solvePhase == 2) {
        model->initBound();
        model->initValue();
    }

}

void HDual::solve_phase2() {
    model->printMessage("dual-phase-2-start");

    // Collect free variables
    dualRow.create_Freelist();

    // Main solving structure
    for (;;) {
        rebuild();
        if (dualInfeasCount > 0)
            break;
        for (;;) {
            model->printProgress();
            switch (dual_variant) {
            default:
            case HDUAL_VARIANT_PLAIN:
                iterate();
                break;
            case HDUAL_VARIANT_TASKS:
                iterate_tasks();
                break;
            case HDUAL_VARIANT_MULTI:
                iterate_multi();
                break;
            }

            if (invertHint)
                break;
        }
        if (model->countUpdate == 0)
            break;
    }

    if (dualInfeasCount > 0) {
        // We got free variables - need dual phase 1
        model->printMessage("dual-phase-2-found-free");
        solvePhase = 1;
    } else if (rowOut == -1) {
        // Tentative optimal
        model->printMessage("dual-phase-2-optimal");
        cleanup();
        if (dualInfeasCount > 0) {
            solvePhase = 4; // Do primal
        } else {
            solvePhase = 0;
            model->printMessage("problem-optimal");
            model->reportStatus(LP_Status_Optimal);
        }
    } else if (columnIn == -1) {
        model->printMessage("dual-phase-2-unbounded");
        if (model->problemPerturbed) {
            cleanup();
        } else {
            solvePhase = -1;
            model->printMessage("problem-infeasible");
            model->reportStatus(LP_Status_Infeasible);
        }
    }
}

void HDual::rebuild() {
    // Save history information
    model->reportPivots(-1, -1, 0); // Indicate REINVERT
    model->timer.recordStart(HTICK_INVERT);

    invertHint = 0;
    // Rebuild model->factor
    if (model->countUpdate > 0) {
        const int *baseIndex = model->getBaseIndex();
        for (int i = 0; i < numRow; i++)
            dualRHS.workDevexFull[baseIndex[i]] = dualRHS.workDevex[i];
        model->computeFactor();
        for (int i = 0; i < numRow; i++)
            dualRHS.workDevex[i] = dualRHS.workDevexFull[baseIndex[i]];
    }

    // Recompute dual solution
    model->computeDual();
    model->correctDual(&dualInfeasCount);
    model->computePrimal();

    // Collect primal infeasible as a list
    dualRHS.create_infeasArray();
    dualRHS.create_infeasList(columnDensity);

    // Compute the objective value
    model->computeObject(solvePhase);
    model->printObject();

    total_INVERT_TICK = factor->pseudoTick;
    total_FT_inc_TICK = 0;
    total_fake = 0;

    model->timer.recordFinish(HTICK_INVERT);
}

void HDual::cleanup() {
    // Remove perturbation and recompute the dual solution
    model->printMessage("dual-cleanup-shift");
    model->initCost();
    model->initBound();
    model->computeDual();
    model->computeObject(solvePhase);
    model->printObject();

    model->computeDualInfeasInPrimal(&dualInfeasCount);
}

void HDual::iterate() {
    chooseRow();
    chooseColumn(&row_ep);

    updateFtranBFRT();
    updateFtran();
    updateFtranDSE(&row_ep);

    updateVerify();
    updateDual();
    updatePrimal(&row_ep);
    updatePivots();
}

void HDual::iterate_tasks() {
    slice_PRICE = 1;

    // Group 1
    chooseRow();

    // Disable slice when too sparse
    if (1.0 * row_ep.count / numRow < 0.01)
        slice_PRICE = 0;

    model->timer.recordStart(HTICK_GROUP1);
#pragma omp parallel
#pragma omp single
    {
#pragma omp task
        {
            columnDSE.copy(&row_ep);
            updateFtranDSE(&columnDSE);
        }
#pragma omp task
        {
            if (slice_PRICE)
                chooseColumn_slice(&row_ep);
            else
                chooseColumn(&row_ep);
#pragma omp task
            updateFtranBFRT();
#pragma omp task
            updateFtran();
#pragma omp taskwait
        }
    }
    model->timer.recordFinish(HTICK_GROUP1);

    updateVerify();
    updateDual();
    updatePrimal(&columnDSE);
    updatePivots();
}

void HDual::chooseRow() {
    if (invertHint)
        return;
    for (;;) {
        // Choose row
        dualRHS.choose_normal(&rowOut);
        if (rowOut == -1) {
            invertHint = 1;
            return;
        }

        // Verify weight
        model->timer.recordStart(HTICK_BTRAN);
        row_ep.clear();
        row_ep.count = 1;
        row_ep.index[0] = rowOut;
        row_ep.array[rowOut] = 1;
        row_ep.packFlag = true;
        factor->btran(row_ep, row_epDensity);
        model->timer.recordFinish(HTICK_BTRAN);
        double u_weight = dualRHS.workDevex[rowOut];
        double c_weight = dualRHS.workDevex[rowOut] = row_ep.norm2();
        if (u_weight >= 0.25 * c_weight)
            break;
    }

    // Assign basic info
    columnOut = model->getBaseIndex()[rowOut];
    if (baseValue[rowOut] < baseLower[rowOut])
        deltaPrimal = baseValue[rowOut] - baseLower[rowOut];
    else
        deltaPrimal = baseValue[rowOut] - baseUpper[rowOut];
    sourceOut = deltaPrimal < 0 ? -1 : 1;
    row_epDensity *= 0.95;
    row_epDensity += 0.05 * row_ep.count / numRow;
}

void HDual::chooseColumn(HVector *row_ep) {
    if (invertHint)
        return;

    // Compute pivot row
    model->timer.recordStart(HTICK_PRICE);
    row_ap.clear();
    matrix->price_by_row(row_ap, *row_ep);
    model->timer.recordFinish(HTICK_PRICE);

    // Choose column - possible
    model->timer.recordStart(HTICK_CHUZC1);
    dualRow.clear();
    dualRow.workDelta = deltaPrimal;
    dualRow.create_Freemove(row_ep);
    dualRow.choose_makepack(&row_ap, 0);
    dualRow.choose_makepack(row_ep, numCol);
    dualRow.choose_possible();
    model->timer.recordFinish(HTICK_CHUZC1);

    // Choose column - check problem
    columnIn = -1;
    if (dualRow.workTheta <= 0 || dualRow.workCount == 0) {
        invertHint = 1;
        return;
    }

    // Choose column - final
    dualRow.choose_final();
    dualRow.delete_Freemove();

    columnIn = dualRow.workPivot;
    alphaRow = dualRow.workAlpha;
    thetaDual = dualRow.workTheta;
}

void HDual::chooseColumn_slice(HVector *row_ep) {
    if (invertHint)
        return;

    model->timer.recordStart(HTICK_CHUZC1);
    dualRow.clear();
    dualRow.workDelta = deltaPrimal;
    dualRow.create_Freemove(row_ep);

    // Row_ep:         PACK + CC1
#pragma omp task
    {

        dualRow.choose_makepack(row_ep, numCol);
        dualRow.choose_possible();
    }

    // Row_ap: PRICE + PACK + CC1
    for (int i = 0; i < slice_num; i++) {
#pragma omp task
        {
            slice_row_ap[i].clear();
            slice_matrix[i].price_by_row(slice_row_ap[i], *row_ep);

            slice_dualRow[i].clear();
            slice_dualRow[i].workDelta = deltaPrimal;
            slice_dualRow[i].choose_makepack(&slice_row_ap[i], slice_start[i]);
            slice_dualRow[i].choose_possible();
        }
    }
#pragma omp taskwait

    // Join CC1 results here
    for (int i = 0; i < slice_num; i++)
        dualRow.choose_joinpack(&slice_dualRow[i]);

    // Infeasible we created before
    columnIn = -1;
    if (dualRow.workTheta <= 0 || dualRow.workCount == 0) {
        invertHint = 1;
        return;
    }
    model->timer.recordFinish(HTICK_CHUZC1);

    // Choose column 2, This only happens if didn't go out
    dualRow.choose_final();
    dualRow.delete_Freemove();
    columnIn = dualRow.workPivot;
    alphaRow = dualRow.workAlpha;
    thetaDual = dualRow.workTheta;
}

void HDual::updateFtranBFRT() {
    if (invertHint)
        return;
    model->timer.recordStart(HTICK_FTRAN_MIX);
    dualRow.update_flip(&columnBFRT);
    if (columnBFRT.count)
        factor->ftran(columnBFRT, columnDensity);
    model->timer.recordFinish(HTICK_FTRAN_MIX);
}

void HDual::updateFtran() {
    if (invertHint)
        return;
    model->timer.recordStart(HTICK_FTRAN);
    column.clear();
    column.packFlag = true;
    matrix->collect_aj(column, columnIn, 1);
    factor->ftran(column, columnDensity);
    alpha = column.array[rowOut];
    model->timer.recordFinish(HTICK_FTRAN);
}

void HDual::updateFtranDSE(HVector *dseVector) {
    if (invertHint)
        return;
    model->timer.recordStart(HTICK_FTRAN_DSE);
    factor->ftran(*dseVector, rowdseDensity);
    model->timer.recordFinish(HTICK_FTRAN_DSE);
}

void HDual::updateVerify() {
    if (invertHint)
        return;

    // The alpha
    double aCol = fabs(alpha);
    double aRow = fabs(alphaRow);
    double aDiff = fabs(aCol - aRow);
    if (aDiff / min(aCol, aRow) > 1e-7 && model->countUpdate > 0) {
        invertHint = 1;
    }

    // We get this thing, but it is not actived by default.
//    // The dual reduced cost
//    double dualin_u = workDual[columnIn];
//    double dualin_c = model->getWorkCost()[columnIn];
//    for (int i = 0; i < column.count; i++) {
//        int iRow = column.index[i];
//        int iCol = model->getBaseIndex()[iRow];
//        double value = column.array[iRow];
//        double cost = model->getWorkCost()[iCol] + model->getWorkShift()[iCol];
//        dualin_c -= cost * value;
//    }
//    double dualin_diff = fabs(dualin_c - dualin_u);
////    if (dualin_diff > Td) {
////        cout << dualin_c << "\t" << dualin_u << "\t" << dualin_diff << endl;
////        invertHint = 1;
////    }
}

void HDual::updateDual() {
    if (invertHint)
        return;

    // Update - dual (shift and back)
    if (thetaDual == 0)
        model->shiftCost(columnIn, -workDual[columnIn]);
    else {
        dualRow.update_dual(thetaDual);
        if (dual_variant != HDUAL_VARIANT_PLAIN && slice_PRICE)
            for (int i = 0; i < slice_num; i++)
                slice_dualRow[i].update_dual(thetaDual);
    }
    workDual[columnIn] = 0;
    workDual[columnOut] = -thetaDual;
    model->shiftBack(columnOut);
}

void HDual::updatePrimal(HVector *dseVector) {
    if (invertHint)
        return;

    // Update - primal and weight
    dualRHS.update_primal(&columnBFRT, 1);
    dualRHS.update_infeasList(&columnBFRT);
    double x_out = baseValue[rowOut];
    double l_out = baseLower[rowOut];
    double u_out = baseUpper[rowOut];
    thetaPrimal = (x_out - (deltaPrimal < 0 ? l_out : u_out)) / alpha;
    dualRHS.update_primal(&column, thetaPrimal);
    double thisDevex = dualRHS.workDevex[rowOut] / (alpha * alpha);
    dualRHS.update_weight(&column, thisDevex, -2 / alpha, &dseVector->array[0]);
    dualRHS.workDevex[rowOut] = thisDevex;
    dualRHS.update_infeasList(&column);

    rowdseDensity = 0.95 * rowdseDensity + 0.05 * dseVector->count / numRow;
    columnDensity = 0.95 * columnDensity + 0.05 * column.count / numRow;

    total_fake += column.fakeTick;
    total_fake += dseVector->fakeTick;

    total_FT_inc_TICK += column.pseudoTick;
    total_FT_inc_TICK += dseVector->pseudoTick;
}

void HDual::updatePivots() {
    if (invertHint)
        return;

    // Update - pivots
    model->updatePivots(columnIn, rowOut, sourceOut);
    model->reportPivots(columnIn, columnOut, alpha);
    model->updateFactor(&column, &row_ep, &rowOut, &invertHint);
    model->updateMatrix(columnIn, columnOut);
    dualRow.delete_Freelist(columnIn);
    dualRHS.update_pivots(rowOut,
            model->getWorkValue()[columnIn] + thetaPrimal);

    if (total_fake >= factor->fakeTick && model->countUpdate >= 50) {
//        cout << total_fake << "\t" << factor->fakeTick << endl;
//        printf(
//                "%10d   %10.2e  %10.2e  %10.4f          %10.2e  %10.2e  %10.4f\n",
//                model->countUpdate, total_INVERT_TICK, total_FT_inc_TICK,
//                total_FT_inc_TICK / total_INVERT_TICK, factor->fakeTick,
//                total_fake, total_fake / factor->fakeTick);
        invertHint = 1;
    }
}

