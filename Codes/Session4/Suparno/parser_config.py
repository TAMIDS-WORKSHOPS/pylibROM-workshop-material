from mfem.common.arg_parser import ArgParser

def get_parser():
    # Initialize argument parser for configuring the wave equation example in MFEM.
    # This setup allows customization of the simulation through command-line arguments.
    
    parser = ArgParser(description="MFEM wave equation example.")
    
    # Mesh file: Defaults to 'star.mesh'. The choice of mesh impacts the geometry and accuracy of the simulation.
    parser.add_argument('-m', '--mesh',
                        default='star.mesh',
                        action='store', type=str,
                        help='Mesh file to use.')
    
    # Serial mesh refinement: Defaults to 2. More refinements increase resolution but also computational cost.
    parser.add_argument('-rs', '--refine-serial',
                        action='store', default=2, type=int,
                        help="Number of times to refine the mesh uniformly in serial")
    
    # Parallel mesh refinement: Defaults to 1. Parallel refinements distribute computational load across cores.
    parser.add_argument('-rp', '--refine-parallel',
                        action='store', default=1, type=int,
                        help="Number of times to refine the mesh uniformly in parallel")
    
    # Finite element order: Defaults to 2 (quadratic elements). Higher order increases accuracy but also complexity.
    parser.add_argument('-o', '--order',
                        action='store', default=2, type=int,
                        help="Finite element order (polynomial degree)")
    
    # ODE solver choice: Defaults to 3 (SDIRK3). Different solvers offer trade-offs in accuracy and stability.
    parser.add_argument('-s', '--ode-solver',
                        action='store', default=3, type=int,
                        help="ODE solver: 1 - Backward Euler, 2 - SDIRK2, 3 - SDIRK3, "
                             "11 - Forward Euler, 12 - RK2, 13 - RK3 SSP, 14 - RK4.")
    
    # Final simulation time: Defaults to 0.5. Determines the duration of the simulation.
    parser.add_argument('-t', '--t-final',
                        action='store', default=0.5, type=float,
                        help="Final time; start time is 0.")
    
    # Time step for the simulation: Defaults to 0.01. Smaller steps improve accuracy but increase computation time.
    parser.add_argument("-dt", "--time-step",
                        action='store', default=0.01, type=float,
                        help="Time step.")
    
    # Alpha coefficient: Defaults to 0.01. Influences the equation's terms, affecting the solution's behavior.
    parser.add_argument('-a', '--alpha',
                        action='store', default=0.01, type=float,
                        help='Alpha coefficient')
    
    # Kappa coefficient: Defaults to 0.5. Similar to alpha, it modulates the equation's terms.
    parser.add_argument('-k', '--kappa',
                        action='store', default=0.5, type=float,
                        help='Kappa coefficient')
    
    # Visualization: Enabled by default. Visual output helps in understanding the simulation's dynamics.
    parser.add_argument('-vis', '--visualization',
                        action='store_true', default=True,
                        help='Enable GLVis visualization')
    
    # Disabling visualization: Overrides the default if specified, useful for running simulations without GUI access.
    parser.add_argument('-no-vis', '--no-visualization',
                        action='store_false', dest='visualization',
                        help='Disable GLVis visualization')
    
    # Save data files for VisIt (visit.llnl.gov) visualization.
    parser.add_argument('-visit', '--visit-datafiles',
                        action='store_true', default=False,
                        help="Save data files for VisIt (visit.llnl.gov) visualization.")
    
    # Visualize every n-th timestep (defaults to every 5th time-step).
    parser.add_argument("-vs", "--visualization-steps",
                    action='store', default=5,  type=int,
                    help="Visualize every n-th timestep.")
    
    # Save data using adios2 streams.
    parser.add_argument("-adios2", "--adios2-streams",
                    action='store_true', default=False,
                    help="Save data using adios2 streams.")
    

    
    return parser
