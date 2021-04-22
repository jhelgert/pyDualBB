// Jonathan Helgert (jhelgert@mail.uni-mannheim.de)

#include "cplexLogging.hpp"

#include <ilconcert/ilothread.h>
#include <ilcplex/ilocplex.h>

#include <cmath>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

void LoggingCallback::logIncumbents(
    IloCplex::Callback::Context const& context) {
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

LoggingCallback::LoggingCallback(IloNumVarArray _x, double _startTime)
    : x(_x), calls(0), branches(0) {
    // Set start timestamp for logging
    startTime = _startTime;
}

void LoggingCallback::invoke(IloCplex::Callback::Context const& context) {
    if (context.inLocalProgress() || context.inCandidate()) {
        logIncumbents(context);
    }
}

int LoggingCallback::getCalls() const {
    return calls;
}

int LoggingCallback::getBranches() const {
    return branches;
}

std::vector<std::tuple<double, double>>
LoggingCallback::getIncumbentsAndTimings() const {
    return incumbentsInfo;
}