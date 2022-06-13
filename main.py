import numpy as np
import matplotlib.pyplot as plt
import pathlib


class AdvectionSolver:
    """Provides functions for a difference method solver for 1D wave problems

    Functions:
        __init__ - creates a new spatial-temporal grid
            nx - number of gridpoints (int)
            left_boundary - left grid boundary for plot (float)
            right_boundary - right grid boundary for plot (float)
            nstep - number of time steps in duration (int)
            duration - length of time for a c=1 wave to travel the entire grid (int)
            bc_type - not implemented, all boundaries are periodic (string)

        initial_value_square - sets up a square wave initial value problem
            left_piecewise - values are 0 until this point, and then 1 (float)
            right_piecewise - values stop being 1 at this point and become 0 (float)

        solve_advection - solves the problem setup at init for -1 < c < 1
            c = speed of the wave, where 1 is the numerical stability limit (float)
            gif = values at which to output a gif [optional] (int)
    """

    def __init__(self, nx: int, left_boundary :float, right_boundary :float,
                 nstep :int, duration :int, bc_type='periodic'):
        self.nx = nx
        self.duration = duration
        self.nstep = nstep
        self.left = left_boundary
        self.right = right_boundary
        self.dx = (abs(left_boundary) + abs(right_boundary)) / nx
        self.dt = duration / nstep
        self.grid = np.linspace(left_boundary, right_boundary, nx)
        self.bc_type = bc_type
        self.adjust_zero = abs(self.left) / self.dx
        self.courant_limit = self.dx / self.dt

    def initial_value_square(self, left_piecewise :float, right_piecewise :float):

        x = self.grid
        self.grid = np.piecewise(x, [x < left_piecewise,
                                    ((x > left_piecewise) & (x < right_piecewise)),
                                    x >= right_piecewise],
                                    [0, 1, 0])

        return

    def solve_advection(self, c: float, gif=0):

        def output_gif(t):
            plt.plot(np.linspace(self.left, self.right, self.nx), grid_array.grid)
            plt.xlim(self.left, self.right)
            cur_time = t * self.dt
            cur_time = "{:.2f}".format(cur_time)
            fn = f"file_{t}.png"
            cur_dir = pathlib.Path('.')
            save_dir = cur_dir / 'figures' / fn
            plt.title(f"Time: {cur_time}")
            plt.savefig(save_dir)
            plt.close()
            return

        def solve_step(courant):
            new_grid = self.grid.copy()
            old_grid = self.grid.copy()

            if c > 0:
                for i in range(self.nx):
                    if i == 0:
                        #catch left side boundary condition before beginning
                        new_grid[i] = old_grid[i] -\
                                      courant * (old_grid[i] - old_grid[self.nx - 1])
                        continue
                    new_grid[i] = old_grid[i] -\
                                  courant * (old_grid[i] - old_grid[i - 1])

            if c < 0:
                # catch right side boundary condition before beginning
                new_grid[nx - 1] = old_grid[nx - 1] +\
                                   courant * (old_grid[nx - 1] - old_grid[0])
                for i in range(self.nx - 1):
                    new_grid[i] = old_grid[i] +\
                                  courant * (old_grid[i] - old_grid[i + 1])

            self.grid = new_grid

        if abs(c) > 1:
            raise ValueError("c must be between -1 and +1.")
        c = c * self.courant_limit
        courant = c * self.dt / self.dx

        for t in range(self.nstep + 1):
            if (t % 100 == 0):
                print(f'{int(t / self.nstep * 100)}% complete.')
            solve_step(courant)
            if (gif > 0):
                if (t % gif == 0):
                    output_gif(t)
        return

    def __repr__(self):
        return self.grid


if __name__ == "__main__":
    nx = 5000               # gridpoints
    nstep = 5000            # timesteps
    duration = 10.0         # simulation run time
    left_boundary = -0.5
    right_boundary = 0.5

    grid_array = AdvectionSolver(nx,
                                 left_boundary,
                                 right_boundary,
                                 nstep,
                                 duration)
    grid_array.initial_value_square(-0.25, 0.25)

    grid_array.solve_advection(1,100)
