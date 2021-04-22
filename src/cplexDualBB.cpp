// Jonathan Helgert (jhelgert@mail.uni-mannheim.de)

#include "cplexDualBB.hpp"

#include <blaze/Math.h>
#include <ilconcert/ilothread.h>
#include <ilcplex/ilocplex.h>

#include <cmath>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

#include "Helpers.hpp"

void BranchCallback::parseJson(std::string& filename) {
    // Read the json file
    std::ifstream i(filename);
    json j;
    i >> j;

    p  = j["p"];
    n  = j["n"];
    m1 = j["m1"];
    m2 = j["m2"];

    // Parse and initialize the Vectors and Matrices
    std::vector<std::vector<double>> Q0_tmp  = j["Q0"];
    std::vector<std::vector<double>> A1_tmp  = j["A1"];
    std::vector<std::vector<double>> A2_tmp2 = j["A2"];
    std::vector<double> c0_tmp               = j["c0"];
    std::vector<double> b1_tmp               = j["b1"];
    std::vector<double> b2_tmp2              = j["b2"];
    if (j.find("r0") != j.end()) {
        r0 = j["r0"];
    }

    // Initialize the blaze matrices and vectors
    Q0  = BMat(n, n, helpers::flatten(Q0_tmp).data());
    A1  = BMat(m1, n, helpers::flatten(A1_tmp).data());
    A1T = blaze::trans(A1);
    c0  = BVec(n, c0_tmp.data());
    b1  = BVec(m1, b1_tmp.data());
    rs  = blaze::zero<double>(p);

    // construct the matrix A2 and vector b2
    // Each thread will work on private copy of A2,b2 later.
    A2_tmp           = blaze::zero<double>(m2 + n, n);
    b2_tmp           = blaze::zero<double>(m2 + n);
    auto A2_tmp_view = blaze::submatrix(A2_tmp, 0UL, 0UL, m2, n);
    auto b2_tmp_view = blaze::subvector(b2_tmp, 0UL, m2);
    A2_tmp_view      = BMat(m2, n, helpers::flatten(A2_tmp2).data());
    b2_tmp_view      = BVec(m2, b2_tmp2.data());

    // Primal quadr. ineq. constraints
    for (size_t i = 0; i < p; ++i) {
        std::vector<std::vector<double>> Qi_tmp =
            j[std::string("Q") + std::to_string(i + 1)];
        std::vector<double> ci_tmp =
            j[std::string("c") + std::to_string(i + 1)];
        Qs.push_back(BMat(n, n, helpers::flatten(Qi_tmp).data()));
        cs.push_back(BVec(n, ci_tmp.data()));
        rs[i] = j[std::string("r") + std::to_string(i + 1)];
    }
}

void BranchCallback::setAlphas(std::vector<double>& _alphas,
                               bool use_nontrivial_alphas) {
    // Lagrange Multipliers
    for (const auto& val : _alphas) {
        // Create Vector alpha = (val, val, val, ...., val)
        // and add it to alphas
        alphas.push_back(blaze::uniform(p, val));
    }

    if (use_nontrivial_alphas) {
        // add the non-uniform alpha vectors, i.e.
        // alpha = (alpha_val 0 ... 0)
        // alpha = (0 ... 0 alpha_val 0 .. 0)
        for (const auto& val : _alphas) {
            for (size_t i = 0; i < p; ++i) {
                BVec tmp = blaze::zero<double>(p);
                tmp[i]   = val;
                alphas.push_back(tmp);
            }
        }
    }
}

void BranchCallback::preCalculateMatrices() {
    // Precalculate Qalpha, Qalphainv, calpha, ralpha for the given
    // alphas.
    for (const auto& alpha : alphas) {
        auto ralpha = r0 + blaze::dot(alpha, rs);
        auto calpha = c0;
        auto Qalpha = Q0;
        // calpha = c0 + alpha_1 * c_1 + .... + alpha_p * c_p
        // Qalpha = Q0 + alpha_1 * Q_1 + .... + alpha_p * Q_p
        for (size_t i = 0; i < p; ++i) {
            calpha += alpha[i] * cs[i];
            Qalpha += alpha[i] * Qs[i];
        }
        // Calculate the inverse Qalpha via Cholesky decomposition
        BMat Qalphainv(Qalpha);
        blaze::invert<blaze::byLLH>(Qalphainv);

        // Add to Qalphas, Qalphainvs, calphas, ralphas
        // Qalphas.push_back(Qalpha);
        Qalphainvs.push_back(Qalphainv);
        calphas.push_back(calpha);
        ralphas.push_back(ralpha);
    }
}

double BranchCallback::calcDualObjective(const BMat& Qalphainv,
                                         const BVec& calpha, double ralpha,
                                         const BMat& A2, const BMat& A2T,
                                         const BVec& b2) const noexcept {
    BVec tmp1 = A2 * Qalphainv * calpha + b2;
    // tmp2 = A2 * Qalpha(-1) * A2^T
    BMat tmp2 = A2 * Qalphainv * A2T;
    // Calculate the inverse of tmp2 via Cholesky decomposition
    // and calculate mu

    // DEBUGGING
    // std::cout << "A2 inv(Qalpha A2^T dims (" << tmp2.rows() << ","
    //           << tmp2.columns() << ")\n";

    // std::cout << tmp2 << "\n";
    // printMatrix(A2);
    // std::cout << A2 << "\n";

    // Invert via cholesky decomposition
    blaze::invert<blaze::byLLH>(tmp2);
    auto mu = -1.0 * tmp2 * tmp1;
    // ----------------------------------------------------------
    // tmp = calpha + A1^T * lambda + A2^T * mu
    BVec tmp = calpha;
    // tmp += A1T * lambda; // lambda = 0
    tmp += A2T * mu;
    // -0.5 * tmp^T * Qalpha^(-1) * tmp + ralpha - b1^T * lambda - b2^T * mu
    // We want to prevent unncessary temporary objects
    // (Unfortunately, BLAZE doesn't provide a suited expression template
    // for our kind of calculation).
    double objVal = 0.0;
    objVal += -0.5 * blaze::trans(tmp) * Qalphainv * tmp;
    objVal += ralpha;
    // objVal -= blaze::dot(b1, lambda); // lambda = 0
    objVal -= blaze::dot(b2, mu);
    // std::cout << "objVal = " << objVal << "\n";
    return objVal;
}

double BranchCallback::findBestDualBound(const BMat& A2,
                                         const BVec& b2) const noexcept {
    // DEBUGGING
    // std::cout << "A2 dims (" << A2.rows() << "," << A2.columns() << ")"
    //           << "\n";
    BMat A2T          = blaze::trans(A2);
    double bestObjVal = std::numeric_limits<double>::min();
    // printf("Before loop inside findBestDualBound()\n");
    for (size_t i = 0; i < alphas.size(); ++i) {
        // printf("inside loop findBestDualbOund()\n");

        // DEBUGGING
        // std::cout << "alpha[" << i << "]"
        //           << "\n";

        double objVal = calcDualObjective(Qalphainvs[i], calphas[i], ralphas[i],
                                          A2, A2T, b2);
        // printf("objVal = %lf\n", objVal);
        if (objVal > bestObjVal) {
            bestObjVal = objVal;
            // bestAlpha = alphas[i];
        }
    }
    return bestObjVal;
}

void BranchCallback::logIncumbents(IloCplex::Callback::Context const& context) {
    // Is a feasible solution known?
    if (context.getLongInfo(
            IloCplex::Callback::Context::Info::Infos::Feasible)) {
        auto obj = context.getIncumbentObjective();
        auto t   = context.getDoubleInfo(
            IloCplex::Callback::Context::Info::Infos::Time);
        if (std::abs(lastIncumbent - obj) > 0.0) {
            lck.lock();
            incumbentsInfo.push_back(std::make_tuple(obj, t - startTime));
            lastIncumbent = obj;
            lck.unlock();
        }
    }
}

BranchCallback::BranchCallback(IloNumVarArray _x, std::string filename,
                               std::vector<double>& _alphas_tmp,
                               double _startTime,
                               bool _only_single_branch       = false,
                               size_t _dual_branching_threads = 4)
    : x(_x),
      calls(0),
      branches(0),
      only_single_branch(_only_single_branch),
      dual_branching_threads(_dual_branching_threads) {
    // Parse the matrices from the json file
    parseJson(filename);
    // Read the lagrangian multipliers alpha
    setAlphas(_alphas_tmp);
    // Pre compute all inverses of Q(alpha) etc.
    preCalculateMatrices();
    // Set start timestamp for logging
    startTime = _startTime;
}

void BranchCallback::invoke(IloCplex::Callback::Context const& context) {
    // std::cout << "Hi from Callback! \n";
    if (context.inLocalProgress() || context.inCandidate()) {
        logIncumbents(context);
    }

    // If not in branching context, do nothing
    if (!context.inBranching()) {
        return;
    }

    // Only run the callback for half of the threads. The other half of threads
    // is used by cplex's branching decisions
    if (context.getLongInfo(IloCplex::Callback::Context::Info::ThreadId) >=
        dual_branching_threads) {
        return;
    }

    // Count the number of calls of the callback
    lck.lock();
    ++calls;
    lck.unlock();

    // Only branch if the current node relaxation could be solved to
    // optimality, otherwise let cplex decide how to branch
    auto status = context.getRelaxationStatus(0);
    if (status != IloCplex::Optimal && status != IloCplex::OptimalInfeas) {
        return;
    }

    // Get the thread id
    auto tid = context.getLongInfo(IloCplex::Callback::Context::Info::ThreadId);

    // Get the active bounds of the current node relaxation

    // Create the new child nodes

    // Create copy of A2, b2 for left and right branch:
    BMat A2_tmp_up(A2_tmp);
    BMat A2_tmp_down(A2_tmp);
    BVec b2_tmp_up(b2_tmp);
    BVec b2_tmp_down(b2_tmp);

    IloNumVar branchVar;
    bool can_branch = false;
    size_t branch_idx;

    // Contains the variable values for the current node relaxation.
    IloNumArray v(context.getEnv());

    lck.lock();  // Start lock for variable object x
    // Write the values of the variables x into v
    context.getRelaxationPoint(x, v);

    size_t k             = 0;
    constexpr double eps = 1.0e-5;
    for (size_t i = 0; i < x.getSize(); ++i) {
        if (x[i].getType() != IloNumVar::Float) {
            double rval = std::round(v[i]);
            if (std::abs(rval - v[i]) < eps) {
                // x[i] is integer and fixed, incorporate it into A2, b2
                A2_tmp_up(m2 + k, i)   = 1.0;
                b2_tmp_up[m2 + k]      = rval;
                A2_tmp_down(m2 + k, i) = 1.0;
                b2_tmp_down[m2 + k]    = rval;
                ++k;
            } else {
                // x[i] is integer and not fixed, so we can branch on it
                if (!can_branch) {
                    can_branch = true;
                    branch_idx = i;
                    branchVar  = x[i];
                }
            }
        }
    }
    lck.unlock();  // end lock for reading variable object x

    // No variable left we can branch on, i.e. stop
    // the callback
    if (!can_branch) {
        return;
    }

    // We branch on x[branch_indx], incorporate it into A2, b2
    double val_up                   = std::ceil(v[branch_idx]);
    double val_down                 = std::floor(v[branch_idx]);
    A2_tmp_up(m2 + k, branch_idx)   = 1.0;
    A2_tmp_down(m2 + k, branch_idx) = 1.0;
    b2_tmp_up[m2 + k]               = val_up;
    b2_tmp_down[m2 + k]             = val_down;

    // Transform A2_up, A2_down to full rank
    BMat A2_tmp_up_T   = blaze::trans(A2_tmp_up);
    BMat A2_tmp_down_T = blaze::trans(A2_tmp_down);
    auto jb_up         = helpers::rref(A2_tmp_up_T);
    auto jb_down       = helpers::rref(A2_tmp_down_T);

    //
    auto A2_up = BMat(blaze::submatrix(A2_tmp_up, 0UL, 0UL, jb_up.size(), n));
    auto A2_down =
        BMat(blaze::submatrix(A2_tmp_down, 0UL, 0UL, jb_down.size(), n));
    auto b2_up   = BVec(blaze::elements(b2_tmp_up, jb_up));
    auto b2_down = BVec(blaze::elements(b2_tmp_down, jb_down));

    // Calculate the dual bounds for each new child's relaxation
    auto dualBound_up   = findBestDualBound(A2_up, b2_up);
    auto dualBound_down = findBestDualBound(A2_down, b2_down);

    if (!only_single_branch) {
        // Branching, i.e. create two new child nodes
        context.makeBranch(branchVar, val_up, IloCplex::BranchUp, dualBound_up);
        lck.lock();
        ++branches;
        lck.unlock();

        context.makeBranch(branchVar, val_down, IloCplex::BranchDown,
                           dualBound_down);
        lck.lock();
        ++branches;
        lck.unlock();
    } else {
        // Create new child node only for the node with better dual bound:
        if (dualBound_up < dualBound_down) {
            context.makeBranch(branchVar, val_up, IloCplex::BranchUp,
                               dualBound_up);
            lck.lock();
            ++branches;
            lck.unlock();
        } else {
            context.makeBranch(branchVar, val_down, IloCplex::BranchDown,
                               dualBound_down);
            lck.lock();
            ++branches;
            lck.unlock();
        }
    }

    // Prune the current node if .... ?
    // context.pruneCurrentNode();
}

int BranchCallback::getCalls() const {
    return calls;
}

int BranchCallback::getBranches() const {
    return branches;
}

std::vector<std::tuple<double, double>>
BranchCallback::getIncumbentsAndTimings() const {
    return incumbentsInfo;
}