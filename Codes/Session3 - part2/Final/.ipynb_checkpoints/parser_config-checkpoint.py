from mfem.common.arg_parser import ArgParser

def get_parser():

    parser = ArgParser(description="Projection ROM - MFEM Poisson equation example.")

    # -m, --mesh: Specifies the mesh file for the simulation.
    # Mesh files define the geometric and topological characteristics of the simulation domain. 
    # 'star.mesh' is the default mesh file, used when no other file is specified.
    parser.add_argument('-m', '--mesh',
                        default='star.mesh',
                        action='store', type=str,
                        help='Mesh file to use.')


    # -o, --order: Determines the polynomial degree for finite elements or selects isoparametric space with -1.
    # Finite element order affects solution accuracy and computational complexity. 
    # A higher order generally increases both accuracy and computational demand.
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree) or -1 for isoparametric space.")


    # -id: Sets a parametric identifier used in simulations.
    # This identifier can be used to distinguish among different parameter sets or simulation instances.
    parser.add_argument("-id", "--id",
                        action='store', default=0, type=int, help="Parametric id")


    # -ns, --nset: Specifies the number of parametric snapshot sets.
    # Snapshot sets are used in model order reduction to represent the solution space efficiently.
    parser.add_argument("-ns", "--nset",
                        action='store', default=0, type=int, help="Number of parametric snapshot sets")


    # -sc, --static-condensation: Enables static condensation if set.
    # Static condensation reduces the system size by eliminating internal degrees of freedom in the finite element assembly process.
    parser.add_argument("-sc", "--static-condensation",
                        action='store_true', default=False,
                        help="Enable static condensation.")


    # -pa, --partial-assembly: Activates Partial Assembly mode.
    # In Partial Assembly, the global system matrix is not fully assembled, reducing memory usage and possibly computation time. 
    parser.add_argument("-pa", "--partial-assembly",
                        action='store_true', default=False,
                        help="Enable Partial Assembly.")


    # -f, --frequency: Sets the frequency for the exact solution in simulations.
    # The frequency parameter can influence the behavior of wave-related simulations or other frequency-dependent analyses. 
    parser.add_argument("-f", "--frequency",
                        action='store', default=1.0, type=float,
                        help="Set the frequency for the exact solution.")


    # -cf, --coefficient: Assigns a coefficient value.
    parser.add_argument("-cf", "--coefficient",
                        action='store', default=1.0, type=float,
                        help="Coefficient.")


    # -d, --device: Configures the computational device (e.g., 'cpu' or 'gpu').
    # Device configuration can significantly affect performance by leveraging specialized hardware capabilities.
    parser.add_argument("-d", "--device",
                        action='store', default='cpu', type=str,
                        help="Device configuration string, see Device::Configure().")


    # -visit, --visit-datafiles: Toggles the saving of data files for VisIt visualization.
    # VisIt is a free interactive parallel visualization and graphical analysis tool for viewing scientific data.
    parser.add_argument("-visit", "--visit-datafiles",
                        action='store_true', default=False,
                        help="Save data files for VisIt (visit.llnl.gov) visualization.")


    # -vis, --visualization: Enables or disables GLVis visualization.
    # GLVis is a lightweight tool for accurate and flexible finite element visualization.
    parser.add_argument("-vis", "--visualization",
                        action='store_true', default=True,
                        help="Enable or disable GLVis visualization.")


    # -fom: Controls the Full Order Model (FOM) phase activation.
    # The FOM phase involves running the original high-fidelity model, often used as a benchmark for reduced models.
    parser.add_argument("-fom", "--fom",
                        action='store_true', default=False,
                        help="Enable or disable the fom phase.")


    # -offline: Enables or disables the offline phase of model order reduction.
    # The offline phase involves pre-computing basis functions or other reduction data, crucial for efficient online computation.
    parser.add_argument("-offline", "--offline",
                        action='store_true', default=False,
                        help="Enable or disable the offline phase.")


    # -online: Toggles the online phase, where the reduced model is actually used for simulation.
    # The online phase is computationally cheaper and faster, relying on data prepared during the offline phase.
    parser.add_argument("-online", "--online",
                        action='store_true', default=False,
                        help="Enable or disable the online phase.")


    # -merge: Enables or disables the merge phase in computational workflows.
    # The merge phase typically involves combining data from different simulation runs or models for comprehensive analysis.
    parser.add_argument("-merge", "--merge",
                        action='store_true', default=False,
                        help="Enable or disable the merge phase.")


    parser.add_argument("-paraview", "--paraview",
                        action='store_true', default=True,
                        help="Enable or disable the paraview visualization.")


    return parser