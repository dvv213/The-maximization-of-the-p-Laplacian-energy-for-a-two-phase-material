#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 04:59:19 2021

@author: donato
"""

import mshr 
import numpy as np
import logging
logging.getLogger('FFC').setLevel(logging.ERROR)
from dolfin import *
import dijitso
cache_params = {"cache_dir":"*","lib_dir":"*","log_dir":"*","inc_dir":"*","src_dir":"*"} 
class pLaplacianEnergyMaxFeaDir:
    """
    Maximization of the p-Laplacian energy fot the mixture of two isotropic materials.
    
    The problem is solved using its relaxed formulation:
        
    $$min F(theta,u):=\frac{1}{p}\int_{\Omega}frac{|\nabla u|^{p}}{(1+ctheta)^{p-1}}dx-\int_{\Omega}fudx $$
    
    for $\theta \in L^{\inft}(\Omega;[0,1])$, $\int_{\Omega }\theta dx\leq k$, $u\in W_{0}^{1,p}(\Omega)$,
    using feasible directions method. 
        
    
    """
    def __init__(self,mesh,p,c,k,f,max_it=100,eps=10**-9,alpha=0.1,beta=0.8,max_it_init=100,tol=10**-1):
        """
        Initialize the space Functions.
        """
        self.mesh=mesh
        self.p=p
        self.f=f
        self.c=c
        self.k=k
        
        self.V = FunctionSpace(mesh, "Lagrange", 1)
        self.Theta= FunctionSpace(mesh, "Discontinuous Lagrange", 0)
        self.V_vec = VectorFunctionSpace(mesh, "Discontinuous Lagrange", 0)
        
        
        self.max_it=max_it
        self.eps=eps
        self.bc=DirichletBC(self.V, 0.0, lambda x,onb: onb)
        self.alpha=alpha
        self.beta=beta
        self.max_it_init=max_it_init
        self.tol=tol
    def __get_initial_condition(self):
        """
        Solve the p-Laplacian equation
        return: the solution of -div(|\nabla u|^{p-2}\nabla u)=f \in \Omega
        """
        u=TrialFunction(self.V)
        v=TestFunction(self.V)
        a= dot(grad(u),grad(v))*dx
        L=self.f*v*dx
        u=Function(self.V)
        solve(a==L,u,self.bc)
        print('Obtaining initial condition')
        for i in range(self.max_it_init):
            d=self.__get_max_des_dir(u,theta=0)
            t=self.__line_search(d,u,theta=0)
            u=project(u+t*d,self.V)
            dF=assemble(dot(grad(d),grad(d))*dx)

            if dF<= self.eps:
                return u
        return u
    def __objective(self,u,theta=0):
        """
        Objective Function
        """
        return assemble((dot(grad(u),grad(u))**(self.p/2)/(1+self.c*theta)**(self.p-1)/self.p-self.f*u)*dx)
    def __line_search(self,d,u,theta=0,dtheta=Constant(0)):
        """ 
        Armijo line search
        """
        t=1
        dF=assemble(dot(grad(d),grad(d))*dx)
        #if dtheta is not None:
        do=assemble(dot(grad(u),grad(u))**((self.p)/2)/(1+self.c*theta)**(self.p)*dtheta*dx)
        do=abs(do*self.c*(self.p-1)/self.p)

        dF+=do
        F=self.__objective(u,theta)
        F_new=self.__objective(u+t*d,theta+t*dtheta)
        s_err=F_new-F+self.alpha*dF*t
        while s_err>=0:
            
            t=t*self.beta
            F_new=self.__objective(u+t*d,theta+t*dtheta)
            s_err=F_new-F+self.alpha*dF*t
            #print('error',s_err,'t',t)
        return t
    def __get_max_des_dir(self,u,theta):
        
        
        d=TrialFunction(self.V)
        v=TestFunction(self.V)
        
        w=project(dot(grad(u),grad(u)),self.Theta)
        w_vector=w.vector()[:]
        w_vector=w_vector**((self.p-2)/2)
        w_vector[(w_vector>=np.infty)&(np.isnan(w_vector))]=0
        w.vector()[:]=w_vector
        #w=(dot(grad(u),grad(u))+self.eps)**((self.p-2)/2)
        a=dot(grad(d),grad(v))*dx
        L=(-w*dot(grad(u),grad(v))/(1+self.c*theta)**(self.p-1)+v*self.f)*dx
        d=Function(self.V)
        solve(a==L,d,self.bc)
        return d
    def __get_max_dir_theta(self,u,theta):
        dtheta=project((dot(grad(u),grad(u)))**(self.p/2)/(1+self.c*theta)**(self.p),
                       self.Theta)
        
        mu_sup=dtheta.vector()[:].max()

        mu_inf=0
        new_theta=Function(self.Theta)
        mu=mu_inf#(mu_sup+mu_inf)/2
        dtheta_vec=dtheta.vector()[:]
        new_theta_vec=new_theta.vector()[:]
        new_theta_vec[dtheta_vec>mu]=1
        new_theta_vec[dtheta_vec<=mu]=0
        new_theta.vector()[:]=new_theta_vec
        int_theta=assemble(new_theta*dx)
        mu_ant=mu
        i=0
        while abs(int_theta-self.k)>self.eps :
            i+=1
            #print(i,int_theta,self.k,self.max_it)
            if int_theta>self.k:
                mu_inf=mu
            else:
                mu_sup=mu
            mu=(mu_sup+mu_inf)/2
            if abs(mu_ant-mu)<self.eps*mu_ant:
                break
            new_theta_vec=new_theta.vector()[:]
            new_theta_vec[dtheta_vec>mu]=1
            new_theta_vec[dtheta_vec<=mu]=0
            new_theta.vector()[:]=new_theta_vec
            int_theta=assemble(new_theta*dx)
            mu_ant=mu
            
        return mu,project(new_theta-theta,self.Theta)
        
    def solve(self):
        """
        Solve the problem. 
        Return the functions evaluations, the length steps and the norma of the gradient.
        """
        ts=[]
        dFs=[]
        Fs=[]
        mus=[]
        u=self.__get_initial_condition()
        theta=Function(self.Theta)
        mu,theta=self.__get_max_dir_theta(u,theta)
        F_old=self.__objective(u,theta)
        for i in range(self.max_it):
            print('Iteration:',i)
            
            d=self.__get_max_des_dir(u,theta)
            mu,dtheta=self.__get_max_dir_theta(u,theta)

            t=self.__line_search(d,u,theta,dtheta)

            dF=assemble(dot(grad(d),grad(d))*dx)+\
abs(assemble(dot(grad(u),grad(u))**((self.p)/2)/(1+self.c*theta)**(self.p)*dtheta*dx)*self.c*(self.p-1)/self.p)

            dFs.append(dF)
            ts.append(t)
            mus.append(mu)
            u=project(u+t*d,self.V)

            theta=project(theta+t*dtheta,self.Theta)

            F=self.__objective(u,theta)

            Fs.append(F)
            print('Objective',F,'Grad_norm:',dF)
            dijitso.cache.clean_cache(cache_params, dryrun=False, 
                                      categories=(u'inc', u'src', u'lib', u'log'))
            #if abs(F_old-F)<self.tol*abs(F_old):
            if dF<self.tol:
                print('Toleranced reached in ',i,' iterations.')
                break
            #F_old=F
        self.u=u
        self.theta=theta
        return ts,dFs,Fs,mus