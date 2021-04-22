from ._cplex_dualbb_wrapper import solveDualBB, py_dict_example


def dualBB(problem_file, json_file, alpha_vals, options=None):
    """Calls the dualBB algorithm

    Args:
        problem_file (str): Path to the problem file (.lp / .mps / .mps.gz)
        json_file (str): Path to the json file containing all problem matrices
        alpha_vals (np.array): Fixed lagrangian multipliers alpha
        options (dict, optional): Parameters used for Cplex. Defaults to None. 
                                  In this case the cplex defaults will be used.

    Returns:
        tuple: cplex status, total runtime in sec, incumbent objvals, incumbent runtimes
    """
    if options is None:
        cplex_options = {
            "mip_strategy_nodeselect": 1,  # depth-first search
            "mip_strategy_search": 0,  # let cplex choose
            "mip_strategy_heuristiceffort": 1,  # default
            "timelimit": 2**31-1  # default
        }
    else:
        cplex_options = options
    # Call the C++ method
    return solveDualBB(problem_file, json_file, alpha_vals, cplex_options)
