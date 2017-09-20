#!/usr/bin/env python
# encoding: utf-8

r"""
Burgers' equation
=========================

Solve the inviscid Burgers' equation:

.. math:: 
    q_t + \frac{1}{2} (q^2)_x & = 0.

This is a nonlinear PDE often used as a very simple
model for fluid dynamics.

The initial condition is sinusoidal, but after a short time a shock forms
(due to the nonlinearity).
"""
from __future__ import absolute_import
import numpy as np
from clawpack import riemann
import pylab as plt
import math
import time as timeyo
from pde import makeMoviespace


class burger_prob:
    def __init__(self,nx):
        global steps 
        steps=nx
    global qinit
        
    global custom_bc
    def custom_bc(state,dim,t,qbc,num_ghost):
            
            qbc[0,0] = 0.0
            qbc[0,dim.num_cells-1]=0.0
    global setup
    def setup(use_petsc=0,kernel_language='Fortran',outdir='./_output',solver_type='classic'):

        if use_petsc:
            import clawpack.petclaw as pyclaw
        else:
            from clawpack import pyclaw

        if kernel_language == 'Python': 
            riemann_solver = riemann.burgers_1D_py.burgers_1D
        elif kernel_language == 'Fortran':
            riemann_solver = riemann.burgers_1D

        if solver_type=='sharpclaw':
            solver = pyclaw.SharpClawSolver1D(riemann_solver)
        else:
            solver = pyclaw.ClawSolver1D(riemann_solver)
            solver.limiters = pyclaw.limiters.tvd.vanleer

        solver.kernel_language = kernel_language
            
        solver.bc_lower[0] = pyclaw.BC.custom
        solver.user_bc_lower = custom_bc
        solver.bc_upper[0] = pyclaw.BC.custom
        solver.user_bc_upper = custom_bc

        x = pyclaw.Dimension(-100.0,100.0,steps,name='x')
        domain = pyclaw.Domain(x)
        num_eqn = 1
        state = pyclaw.State(domain,num_eqn)

        xc = state.grid.x.centers
        state.q[0,:] = np.exp(-xc**2)
        state.problem_data['efix']=True

        claw = pyclaw.Controller()
        claw.tfinal = 10.0
        claw.solution = pyclaw.Solution(state,domain)
        claw.solver = solver
        claw.outdir = outdir
        claw.setplot = setplot
        claw.keep_copy = True

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

    def run(self):
        x = setup()
        x.run()
        solutions = x.frames
        grid_list=[]
        sol_list=[]
        len(solutions)
        for sol in solutions:
            x = sol.states[0].grid.dimensions[0].centers
            q = sol.states[0].q[0, :]
            grid_list.append(x)
            sol_list.append(q)
            #plt.plot(x,q)
            #plt.show()

        a=100.0
        f = lambda x : math.exp(-x**2)
        true_sol = lambda t,x: f(x-a*t)
        time =np.linspace(0.0,1.0,len(grid_list))
        print time
        print len(grid_list)
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


        return (grid_list,sol_list,sol_list ,time)
        #return max(error_list)


nx=8000
x = burger_prob(nx)
v= x.run()
makeMoviespace(v[0],v[1],v[2],v[3],"Numerical Solution", "True Solution","x","u")

# error_lists=[]
# for i in range(0,4):
#     nx = 1000*(2**i)
#     x = burger_prob(nx)
#     error_lists.append(x.run())


# orders =[]
# for i in range(0,len(error_lists)-1):
#     orders.append(error_lists[i]/error_lists[i+1])
# print orders