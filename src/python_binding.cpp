#include <ilconcert/ilothread.h>
#include <ilcplex/ilocplex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include "cplexDualBB.hpp"
#include "cplexLogging.hpp"

using ms = std::chrono::milliseconds;

namespace py = pybind11;

template <typename T>
inline py::array_t<T> toPyArray(std::vector<T>&& passthrough) {
    // Pass result back to python
    auto* transferToHeapGetRawPtr = new std::vector<T>(std::move(passthrough));

    const py::capsule freeWhenDone(transferToHeapGetRawPtr, [](void* toFree) {
        delete static_cast<std::vector<T>*>(toFree);
    });

    auto passthroughNumpy =
        py::array_t<T>({transferToHeapGetRawPtr->size()},  // shape
                       {sizeof(T)},                        // strides
                       transferToHeapGetRawPtr->data(),    // ptr
                       freeWhenDone);
    return passthroughNumpy;
}

py::tuple solveDualBB(std::string problem_file, std::string json_file,
                      py::array_t<double> _alpha_vals,
                      std::map<std::string, int>& cplex_options) {
    // Create STL vector from numpy array
    std::vector<double> alpha_vals(_alpha_vals.data(),
                                   _alpha_vals.data() + _alpha_vals.size());

    // create cplex model object
    IloEnv env;
    IloModel model(env);
    IloNumVarArray x(env);
    IloObjective obj(env);
    IloRangeArray cons(env);
    IloCplex cplex(model);

    // Read the problem file (.lp or .mps or .mps.gz)
    cplex.importModel(model, problem_file.c_str(), obj, x, cons);

    // Set the callback options:

    // Set the Cplex options:

    cplex.setParam(IloCplex::Param::MIP::Strategy::NodeSelect,
                   cplex_options["mip_strategy_nodeselect"]);
    cplex.setParam(IloCplex::Param::MIP::Strategy::Search,
                   cplex_options["mip_strategy_search"]);
    cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicEffort,
                   cplex_options["mip_strategy_heuristiceffort"]);
    cplex.setParam(IloCplex::Param::TimeLimit,
                   (double)cplex_options["timelimit"]);

    // Start timing for logging purposes
    double start_time = cplex.getCplexTime();

    // Register the dualBB callback function and read all the problem matrices
    // etc from the json_file
    BranchCallback cb(x, json_file, alpha_vals, start_time);

    // Set the different contexts where our callback gets invoked
    cplex.use(&cb, IloCplex::Callback::Context::Id::Branching |
                       IloCplex::Callback::Context::Id::LocalProgress |
                       IloCplex::Callback::Context::Id::Candidate);

    // solve the model
    cplex.solve();

    std::cout << "hi after solve()\n";

    // total runtime
    double total_time = cplex.getCplexTime() - start_time;

    std::cout << "hi before getIncumbentsAndTimings()\n";

    // Extract all the incumbent objective values and the runtimes
    // Extract the incumbents and the runtimes
    std::vector<double> incObjVals;
    std::vector<double> incRuntimes;
    for (const auto& [objVal, t] : cb.getIncumbentsAndTimings()) {
        incObjVals.push_back(objVal);
        incRuntimes.push_back(t);
    }

    std::stringstream buffer;
    buffer << cplex.getStatus();

    std::string status(buffer.str());

    std::cout << "hi before env.end()\n";

    // End cplex model object lifetime
    env.end();

    std::cout << "hi bevore py::make_tuple()\n";

    // Create python tuple
    return py::make_tuple(status, total_time,
                          toPyArray<double>(std::move(incObjVals)),
                          toPyArray<double>(std::move(incRuntimes)));
}

py::tuple solveCplex(std::string problem_file,
                     std::map<std::string, int>& cplex_options) {
    // create cplex model object
    IloEnv env;
    IloModel model(env);
    IloNumVarArray x(env);
    IloObjective obj(env);
    IloRangeArray cons(env);
    IloCplex cplex(model);

    // Read the problem file (.lp or .mps or .mps.gz)
    cplex.importModel(model, problem_file.c_str(), obj, x, cons);
    // Set the Cplex options:

    cplex.setParam(IloCplex::Param::MIP::Strategy::NodeSelect,
                   cplex_options["mip_strategy_nodeselect"]);
    cplex.setParam(IloCplex::Param::MIP::Strategy::Search,
                   cplex_options["mip_strategy_search"]);
    cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicEffort,
                   cplex_options["mip_strategy_heuristiceffort"]);
    cplex.setParam(IloCplex::Param::TimeLimit,
                   (double)cplex_options["timelimit"]);

    // Start timing for logging purposes
    double start_time = cplex.getCplexTime();

    // Register the dualBB callback function and read all the problem matrices
    // etc from the json_file
    LoggingCallback cb(x, start_time);

    // Set the different contexts where our callback gets invoked
    cplex.use(&cb, IloCplex::Callback::Context::Id::LocalProgress |
                       IloCplex::Callback::Context::Id::Candidate);

    // solve the model
    cplex.solve();

    std::cout << "hi after solve()\n";

    // total runtime
    double total_time = cplex.getCplexTime() - start_time;

    std::cout << "hi before getIncumbentsAndTimings()\n";

    // Extract all the incumbent objective values and the runtimes
    // Extract the incumbents and the runtimes
    std::vector<double> incObjVals;
    std::vector<double> incRuntimes;
    for (const auto& [objVal, t] : cb.getIncumbentsAndTimings()) {
        incObjVals.push_back(objVal);
        incRuntimes.push_back(t);
    }

    std::stringstream buffer;
    buffer << cplex.getStatus();

    std::string status(buffer.str());

    std::cout << "hi before env.end()\n";

    // End cplex model object lifetime
    env.end();

    std::cout << "hi bevore py::make_tuple()\n";

    // Create python tuple
    return py::make_tuple(status, total_time,
                          toPyArray<double>(std::move(incObjVals)),
                          toPyArray<double>(std::move(incRuntimes)));
}

void py_dict_example(std::map<std::string, int>& options) {
    for (const auto& [key, value] : options) {
        std::cout << key << " = " << value << "; ";
    }
    std::cout << "\n";
}

PYBIND11_MODULE(_cplex_dualbb_wrapper, m) {
    m.def("solveDualBB", &solveDualBB,
          "Solves a convex MIQCQP by the dualBB algorithm.");
    m.def("solveCplex", &solveCplex,
          "Solves a convex MIQCQP by Cplex and logs each found incumbent");
    m.def("py_dict_example", &py_dict_example, "bla bla");
}