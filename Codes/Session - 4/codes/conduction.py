# Attempt to import the parallel version of mfem (PyMFEM). If unsuccessful, provide instructions for installation.
import mfem.par as mfem
import numpy as np

# Import specific functionalities from mfem after ensuring it is installed
from mfem.par import intArray

class ConductionOperator(mfem.PyTimeDependentOperator):
    """
    A class for simulating the heat conduction process in a given finite element space.

    This operator extends the PyTimeDependentOperator from MFEM to solve time-dependent heat conduction problems using finite element methods. It supports both explicit and implicit time integration schemes to compute the future state of the system based on its current state.

    Parameters:
    - fespace (mfem.ParFiniteElementSpace): The finite element space for the simulation.
    - alpha (float): A coefficient representing the thermal diffusivity of the material.
    - kappa (float): A coefficient related to the conductivity of the material.
    - u (mfem.Vector): The initial temperature distribution across the finite element space.

    Attributes:
    - alpha (float): Thermal diffusivity coefficient.
    - kappa (float): Conductivity coefficient.
    - T (mfem.HypreParMatrix or None): The total system matrix used in implicit solves. None until first use.
    - K (mfem.ParBilinearForm or None): Stiffness matrix bilinear form. None until `SetParameters` is called.
    - M (mfem.ParBilinearForm or None): Mass matrix bilinear form. None until initialized in `__init__`.
    - fespace (mfem.ParFiniteElementSpace): The finite element space for the simulation.
    - Mmat (mfem.HypreParMatrix): The Hypre parallel matrix representing the mass matrix.
    - Kmat (mfem.HypreParMatrix): The Hypre parallel matrix representing the stiffness matrix.
    - M_solver (mfem.CGSolver): Conjugate Gradient solver for the mass matrix.
    - M_prec (mfem.HypreSmoother): Preconditioner for the mass matrix solve.
    - T_solver (mfem.CGSolver): Conjugate Gradient solver used in the implicit solve step.
    - T_prec (mfem.HypreSmoother): Preconditioner for the total system matrix solve.
    - z (mfem.Vector): Temporary vector for intermediate calculations.

    Methods:
    - Mult(u, du_dt): Performs an explicit step of the time integration.
    - ImplicitSolve(dt, u, du_dt): Performs an implicit step of the time integration.
    - SetParameters(u): Updates the system matrices based on the current solution.

    The class utilizes the Conjugate Gradient solver with a Jacobi preconditioner for solving the linear system of equations. It is designed to handle both the construction of the system matrices and the solution of the equations for each time step, allowing for the simulation of heat conduction in complex geometries and materials.
    """

    def __init__(self, fespace, alpha, kappa, u):
        """
        Initializes the conduction operator with the given finite element space, material properties, and initial condition.
    
        """
        mfem.PyTimeDependentOperator.__init__(self, fespace.GetTrueVSize(), 0.0)        
        self.alpha = alpha
        self.kappa = kappa
        self.fespace = fespace
        self.T = None
        self.K = None
        self.M = None

        self.ess_tdof_list = intArray()
        self.Mmat = mfem.HypreParMatrix()
        self.Kmat = mfem.HypreParMatrix()
        self.M_solver = mfem.CGSolver(fespace.GetComm())
        self.M_prec = mfem.HypreSmoother()
        self.M_prec.SetType(mfem.HypreSmoother.Jacobi)
        self.T_solver = mfem.CGSolver(fespace.GetComm())
        self.T_prec = mfem.HypreSmoother()
        self.z = mfem.Vector(self.Height())

        # Initialize and assemble the mass matrix.
        self.M = mfem.ParBilinearForm(fespace)
        self.M.AddDomainIntegrator(mfem.MassIntegrator())
        self.M.Assemble()
        self.M.FormSystemMatrix(self.ess_tdof_list, self.Mmat)

        # Configure the solvers.
        self.configure_solver(self.M_solver, self.M_prec, 1e-8)
        self.M_solver.SetOperator(self.Mmat)

        self.configure_solver(self.T_solver, self.T_prec, 1e-8)

        # Set up the stiffness matrix based on the initial condition.
        self.SetParameters(u)

    def configure_solver(self, solver, preconditioner, rel_tol):
        """
        Configures a given solver with the specified relative tolerance, preconditioner, and other parameters.
        """
        solver.iterative_mode = False
        solver.SetRelTol(rel_tol)
        solver.SetAbsTol(0.0) # AbsTol is not functional only relTol matters.
        solver.SetMaxIter(100)
        solver.SetPrintLevel(0) # Do not print anything while solving
        solver.SetPreconditioner(preconditioner)

    def Mult(self, u, du_dt):
        """
        Performs an explicit step to compute the derivative of u with respect to time.
        """
        self.Kmat.Mult(u, self.z)  # Apply stiffness matrix to u.
        self.z.Neg()  # Negate z to prepare the right-hand side of M*du_dt = -K*u.
        self.M_solver.Mult(self.z, du_dt)  # Solve the equation.

    def ImplicitSolve(self, dt, u, du_dt):
        """
        Performs an implicit step to solve for du_dt, considering the current state u and time step dt.
        """
        if self.T is None:
            self.assemble_T_matrix(dt)
            current_dt = dt  # Update the time step.

        self.Kmat.Mult(u, self.z)  # Apply stiffness matrix to u.
        self.z.Neg()  # Prepare the right-hand side of the equation.
        self.T_solver.Mult(self.z, du_dt)  # Solve T*du_dt = -K*u.

    def assemble_T_matrix(self, dt):
        """
        Assembles the system matrix T for the implicit solve based on the current dt.
        """
        # Assemble the matrix T = M + dt*K. This matrix is used in the implicit solve step.
        self.T = mfem.Add(1.0, self.Mmat, dt, self.Kmat)
        self.T_solver.SetOperator(self.T)

    def SetParameters(self, u):
        """
        Updates the material and simulation parameters based on the current state u.
        
        This method updates the diffusion coefficient for the stiffness matrix based on the current solution u and
        the thermal properties alpha and kappa. The stiffness matrix K is then assembled with this updated coefficient.
        
        Parameters:
        - u: Current temperature distribution as a mfem.Vector.
        """
        # Initialize a grid function from the finite element space for updating diffusion coefficients.
        u_alpha_gf = mfem.ParGridFunction(self.fespace)
        u_alpha_gf.SetFromTrueDofs(u)

        # Update the grid function values based on kappa and alpha. In this example, we assume a simple model where
        # the diffusion coefficient is directly proportional to kappa and alpha. However, this can be modified to
        # include more complex relationships or spatial variations.
        
        for i in range(u_alpha_gf.Size()):
            u_alpha_gf[i] = self.kappa + self.alpha * u_alpha_gf[i]

        # Assemble the stiffness matrix K with the updated diffusion coefficient.
        self.K = mfem.ParBilinearForm(self.fespace)
        u_coeff = mfem.GridFunctionCoefficient(u_alpha_gf)
        self.K.AddDomainIntegrator(mfem.DiffusionIntegrator(u_coeff))
        self.K.Assemble(0)
        # self.K.Finalize()
        self.K.FormSystemMatrix(self.ess_tdof_list, self.Kmat)

        # Reset the total system matrix T to ensure it's updated in the next implicit solve.
        self.T = None

class InitialTemperature(mfem.PyCoefficient):
    """
    Defines an initial temperature distribution as a function of position.
    
    This class extends mfem.PyCoefficient and specifies an initial temperature distribution for a simulation. The temperature is determined based on the Euclidean norm of the position vector x. Inside a radius of 0.5 units from the origin, the temperature is set to 2.0. Outside this radius, it defaults to 1.0.
    """

    def __init__(self):
        """
        Initializes the InitialTemperature object.
        """
        super().__init__()

    def EvalValue(self, x):
        """
        Evaluates the temperature at a given position x.
        
        Parameters:
        - x (array): Position at which to evaluate the initial temperature, given as a NumPy array.
        
        Returns:
        - float: The evaluated temperature value at position x.
        """
        # Calculate the Euclidean norm (L2 norm) of the position vector x.
        norm2 = np.sqrt(np.sum(x**2))
        
        # Return a temperature of 2.0 for points within a radius of 0.5 units from the origin, and 1.0 otherwise.
        if norm2 < 0.5:
            return 2.0
        return 1.0