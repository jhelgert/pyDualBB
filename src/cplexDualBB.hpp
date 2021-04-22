// Jonathan Helgert (jhelgert@mail.uni-mannheim.de)

#include "Helpers.hpp"
#include <blaze/Math.h>
#include <cmath>
#include <ilconcert/ilothread.h>
#include <ilcplex/ilocplex.h>
#include <limits>
#include <nlohmann/json.hpp>
#include <string>
#include <tuple>
#include <vector>

// Since CPLEX already uses multiple threads,
// we should to enforce the serial execution of some BLAZE operations,
// i.e. the calculation of matrix inverses etc.
#define BLAZE_USE_SHARED_MEMORY_PARALLELIZATION 0

// Aliases
using BVec = blaze::DynamicVector<double>;
using BVecIdx = blaze::DynamicVector<size_t>;
using BMat = blaze::DynamicMatrix<double>;
using json = nlohmann::json;

/**
 * @brief Implementation of the dualBB algorithm via generic Callback
 *
 */
class BranchCallback : public IloCplex::Callback::Function {
  private:
    // Cplex data
    IloNumVarArray x;
    int calls;
    int branches;
    double startTime;
    double lastIncumbent;
    std::vector<std::tuple<double, double>> incumbentsInfo;
    IloFastMutex lck{};
    // ----
    // Relaxed QCQP data:
    size_t n{0};  //< number of primal variables
    size_t p{0};  //< number of quadratic inequality constraints
    size_t m1{0}; //< number of affin-linear inequality constraints
    size_t m2{0}; //< number of affin-linear equality constraints
    // Primal objective function: 0.5*x' * Q0 * x + c0' * x + r0
    BMat Q0;
    BVec c0;
    double r0{0.0};
    // Primal quadr. ineq. constraints:
    // g_i(x) = 0.5* x' * Qs[i] * x + cs[i]'*x + rs[i] <= 0, i = 1, ..., p
    std::vector<BMat> Qs;
    std::vector<BVec> cs;
    BVec rs;
    // Primal affin-linear ineq. constraint
    // h_1(x) = A1 * x - b1 <= 0
    BMat A1;
    BMat A1T; // transpose of A1
    BVec b1;
    // Primal affin-linear eq. constraint
    // h_2(x) = A2 * x - b2 == 0
    // Each thread has its private matrix A2 and vector b2
    BMat A2_tmp;
    BVec b2_tmp;

    // Lagrangian Multipliers
    std::vector<BVec> alphas;

    // dual problem objective
    std::vector<BMat> Qalphainvs;
    std::vector<BVec> calphas;
    std::vector<double> ralphas;

    void parseJson(std::string& filename);

    void setAlphas(std::vector<double>& _alphas,
                   bool use_nontrivial_alphas = false);

    void preCalculateMatrices();

    double calcDualObjective(const BMat& Qalphainv, const BVec& calpha,
                             double ralpha, const BMat& A2, const BMat& A2T,
                             const BVec& b2) const noexcept;

    double findBestDualBound(const BMat& A2, const BVec& b2) const noexcept;

    void logIncumbents(IloCplex::Callback::Context const& context);

  public:
    BranchCallback() = delete;
    BranchCallback(BranchCallback& other) = delete;
    BranchCallback(BranchCallback&& other) = delete;
    BranchCallback& operator=(BranchCallback&& other) = delete;
    ~BranchCallback() = default;

    BranchCallback(IloNumVarArray _x, std::string filename,
                   std::vector<double>& _alphas_tmp, double _startTime);

    void invoke(IloCplex::Callback::Context const& context) override;
    int getCalls() const;
    int getBranches() const;
    std::vector<std::tuple<double, double>> getIncumbentsAndTimings() const;
};