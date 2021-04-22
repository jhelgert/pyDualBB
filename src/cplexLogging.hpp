#ifndef CPLEX_LOGGING_HPP
#define CPLEX_LOGGING_HPP

// Jonathan Helgert (jhelgert@mail.uni-mannheim.de)
#include <ilconcert/ilothread.h>
#include <ilcplex/ilocplex.h>

#include <cmath>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

class LoggingCallback : public IloCplex::Callback::Function {
   private:
    // Cplex data
    IloNumVarArray x;
    int calls;
    int branches;
    double startTime;
    double lastIncumbent;
    std::vector<std::tuple<double, double>> incumbentsInfo;
    IloFastMutex lck{};
    void logIncumbents(IloCplex::Callback::Context const& context);

   public:
    LoggingCallback()                        = delete;
    LoggingCallback(LoggingCallback& other)  = delete;
    LoggingCallback(LoggingCallback&& other) = delete;
    LoggingCallback& operator=(LoggingCallback&& other) = delete;
    ~LoggingCallback()                                  = default;

    LoggingCallback(IloNumVarArray _x, double _startTime);

    void invoke(IloCplex::Callback::Context const& context) override;
    int getCalls() const;
    int getBranches() const;
    std::vector<std::tuple<double, double>> getIncumbentsAndTimings() const;
};

#endif  // Header guard CPLEX_LOGGING_HPP