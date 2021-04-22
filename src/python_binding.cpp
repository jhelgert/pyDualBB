#include "cplexDualBB.hpp"
#include <algorithm>
#include <fstream>
#include <ilconcert/ilothread.h>
#include <ilcplex/ilocplex.h>
#include <iostream>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <vector>

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
        py::array_t<T>({transferToHeapGetRawPtr->size()}, // shape
                       {sizeof(T)},                       // strides
                       transferToHeapGetRawPtr->data(),   // ptr
                       freeWhenDone);
    return passthroughNumpy;
}

py::tuple solveDualBB(std::string problem_file, std::string json_file,
                      py::array_t<double> _alpha_vals) {
    // Create STL vector from numpy array
    std::vector<double> alpha_vals(_alpha_vals.data(),
                                   _alpha_vals.data() + _alpha_vals.size());

    // create cplex model object
    IloEnv env;
    IloModel model(env);
    IloNumVarArray x(env);
    IloCplex cplex(model);

    // Read the problem file (.lp or .mps or .mps.gz)
    cplex.importModel(model, filename.c_str());

    // Set the callback options:

    // Use best-estimate search as node search strategx
    cplex.setParam(IloCplex::Param::MIP::Strategy::NodeSelect,
                   CPX_NODESEL_BESTEST);

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

    // total runtime
    double total_time = cplex.getCplexTime() - start_time;

    // Extract all the incumbent objective values and the runtimes
    // Extract the incumbents and the runtimes
    std::vector<double> incObjVals;
    std::vector<double> incRuntimes;
    for (const auto& [objVal, t] : cb.getIncumbentsAndTimings()) {
        incObjVals.push_back(objVal);
        incRuntimes.push_back(t);
    }

    // End cplex model object lifetime
    env.end();

    // Create python tuple
    return py::make_tuple(cplex.getStatus(), total_time,
                          toPyArray<double>(std::move(incObjVals)),
                          toPyArray<double>(std::move(incRuntimes)));
}

PYBIND11_MODULE(_cplexmiqcqp_wrapper, m) {
    m.def("solveDualBB", &solveDualBB,
          "Solves a convex MIQCQP by the dualBB algorithm.");
}