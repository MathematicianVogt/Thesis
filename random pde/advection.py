#!/usr/bin/env python
# encoding: utf-8
r"""
One-dimensional advection with variable velocity
================================================

Solve the conservative variable-coefficient advection equation:

.. math:: q_t + (u(x)q)_x = 0.

Here q is the density of some conserved quantity and u(x) is the velocity.
The velocity field used is

.. math:: u(x) = 2 + sin(2\pi x).

The boundary conditions are periodic.
The initial data get stretched and compressed as they move through the
fast and slow parts of the velocity field.
"""


from __future__ import absolute_import
import numpy as np
from clawpack.pyclaw import plot
import pylab as plt
import math
from pde import makeMoviespace
from pde import *



class advection_prob:
    def __init__(self,nx):
        global steps 
        global grid_list
        global time
        global sol_list
        global true_sol_list
        steps=nx
    global qinit
    def qinit(state):

        

        x =state.grid.x.centers
        
        state.q[0,:] = np.sin(x)

    global auxinit
    def auxinit(state):
        # Initilize petsc Structures for aux
        xc=state.grid.x.centers
        #state.aux[0,:] = 1.0
        state.aux[0,:] = np.random.uniform(-20.0*(1.0+xc-xc),7.0*(1.0+xc-xc))
        
    global custom_bc
    def custom_bc(state, dim, t, qbc, auxbc, num_ghost):
        for i in xrange(num_ghost):
            qbc[0,i] = 0.0
    global setup
    def setup(use_petsc=False,solver_type='classic',kernel_language='Python',outdir='./_output'):
        from clawpack import riemann

        if use_petsc:
            import clawpack.petclaw as pyclaw
        else:
            from clawpack import pyclaw

        if solver_type=='classic':
            if kernel_language == 'Fortran':
                solver = pyclaw.ClawSolver1D(riemann.vc_advection_1D)
            elif kernel_language=='Python': 
                solver = pyclaw.ClawSolver1D(riemann.vc_advection_1D_py.vc_advection_1D)
        elif solver_type=='sharpclaw':
            if kernel_language == 'Fortran':
                solver = pyclaw.SharpClawSolver1D(riemann.vc_advection_1D)
            elif kernel_language=='Python': 
                solver = pyclaw.SharpClawSolver1D(riemann.vc_advection_1D_py.vc_advection_1D)
            solver.weno_order=weno_order
        else: raise Exception('Unrecognized value of solver_type.')

        solver.kernel_language = kernel_language

        solver.limiters = pyclaw.limiters.tvd.MC
        solver.bc_lower[0] = pyclaw.BC.custom
        solver.user_bc_lower = custom_bc
        solver.bc_upper[0] = pyclaw.BC.custom
        solver.user_bc_upper = custom_bc
        solver.aux_bc_lower[0] = 2
        solver.aux_bc_upper[0] = 2

        xlower=-100.0; xupper=100.0; mx=steps
        x = pyclaw.Dimension(xlower,xupper,mx,name='x')
        domain = pyclaw.Domain(x)
        num_aux=1
        num_eqn = 1
        state = pyclaw.State(domain,num_eqn,num_aux)

        qinit(state)
        auxinit(state)

        claw = pyclaw.Controller()
        claw.outdir = outdir
        claw.solution = pyclaw.Solution(state,domain)
        claw.solver = solver

        claw.tfinal = 10.0
        claw.setplot = setplot
        claw.keep_copy = True
        
        
       # print claw.solution._get_solution_attribute
        return claw

    #--------------------------
    global setplot
    def setplot(plotdata):
    #--------------------------
        """ 
        Specify what is to be plotted at each frame.
        Input:  plotdata, an instance of visclaw.data.ClawPlotData.
        Output: a modified version of plotdata.
        """ 
        plotdata.clearfigures()  # clear any old figures,axes,items data

        # Figure for q[0]
        plotfigure = plotdata.new_plotfigure(name='q', figno=1)

        # Set up for axes in this figure:
        plotaxes = plotfigure.new_plotaxes()
        plotaxes.ylimits = [-.1,1.1]
        plotaxes.title = 'q'

        # Set up for item on these axes:
        plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
        plotitem.plot_var = 0
        plotitem.plotstyle = '-o'
        plotitem.color = 'b'
        plotitem.kwargs = {'linewidth':2,'markersize':5}
        
        return plotdata
    def get_sol(self):
        return (grid_list,sol_list,true_sol_list ,time)
    def run(self):
        x = setup()
        x.run()
        solutions = x.frames
        global grid_list
        grid_list=[]
        global sol_list
        sol_list=[]
        len(solutions)
        for sol in solutions:
            x = sol.states[0].grid.dimensions[0].centers
            q = sol.states[0].q[0, :]
            grid_list.append(x)
            sol_list.append(q)
            #plt.plot(x,q)
            #plt.show()

        a=1.0
        f = lambda x : math.exp(-x**2)
        true_sol = lambda t,x: f(x-a*t)
        global time
        time =np.linspace(0.0,10.0,11)
        print time
        print len(grid_list)
        global true_sol_list
        true_sol_list=[]
        counter=0


        for grid in grid_list:
            current_sol=[]
            for i in range(0,len(grid)):
                current_sol.append(true_sol(time[counter],grid[i]))
            counter+=1
            true_sol_list.append(current_sol)


        error_list = []
        counter=0
        for sol in true_sol_list:
            error=[]
            for i in range(0,len(sol)):
                error.append(math.fabs(sol[i]-sol_list[counter][i]))
                # print sol[i]
                # print sol_list[counter][i]
            counter+=1
            print max(error)
            error_list.append(max(error))


        return (grid_list,sol_list,true_sol_list ,time)
        #return max(error_list)


nx=8000
x = advection_prob(nx)
v= x.run()
#makeMoviespace(v[0],v[1],v[2],v[3],"Numerical Solution", "True Solution","x","u")
makeMoviead(v[0],v[1],v[3],"Numerical Solution","x", "u")
#-----IMPORTANT FOR CONVERGENCE RATE
# error_lists=[]
# for i in range(0,4):
#     nx = 1000*(2**i)
#     x = advection_prob(nx)
#     error_lists.append(x.run())


# orders =[]
# for i in range(0,len(error_lists)-1):
#     orders.append(error_lists[i]/error_lists[i+1])
# print orders

# order_list=[]

# a = np.zeros((4,4))

# for i in range(0,len(error_lists)):
#     for j in range(0, len(error_lists[i])-1):
#         try:
#             a[i,j] = error_lists[i][j]/error_lists[i][j+1]
#         except:
#             pass

# print a













#print sol_list[0]
#print true_sol_list[0]
















# if __name__=="__main__":
#     from clawpack.pyclaw.util import run_app_from_main
#     output = run_app_from_main(setup,setplot)
#     #plot.interactive_plot() 