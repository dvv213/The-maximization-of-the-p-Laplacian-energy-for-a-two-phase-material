
import mshr 
import numpy as np
import logging
logging.getLogger('FFC').setLevel(logging.ERROR)
from dolfin import *
import dijitso
from time import process_time
set_log_level(1000)
cache_params = {"cache_dir":"*","lib_dir":"*","log_dir":"*","inc_dir":"*","src_dir":"*"} 

def create_mesh(resolution):
    N_vertices=30#verticex number for polygonal aproximation of the circle.
    circle =mshr.Circle(Point(0.,0.),1.0,N_vertices)
    mesh = mshr.generate_mesh(circle,resolution)
    return mesh
def create_spaces(mesh):
    V = FunctionSpace(mesh, "Lagrange", 1)
    Theta= FunctionSpace(mesh, "Discontinuous Lagrange", 0)
    return V,Theta

def solve_problem(Theta,V,
          p,c,k,f,method,
          max_it=100,eps=10**-9,
          alpha=0.01,beta=0.5,
          max_it_init=100,tol=10**-3,
          max_it_theta=100):
    """
    Solve the relaxed formulation of the maximization of the p-Lplacian energy problem.
    
    
    Theta : Finite element for the theta vatiable.
    V : Finite element for the u variable.
    p : constant.
    c : constant.
    k : volume (or area) restriction.
    f : source term. It has to be a finite element function or a Expression in dolinf.
    method : {'feasible directions','alternating optimization'} method to solve the problem.
    max_it : Maximun number of iterations. This does not include the number 
    for iterations fot the initial condition.
    eps : tolerance for get the optimal theta or the maimun descent direction in theta variable.
    alpha : alpha constant for the Armijo's rule.
    beta : beta constant fot the Armijo's rule.
    max_it_init : maximun number of iterations for obtain the initial condition.
    tol : the tolerance for the gradient. if (dF*(du,do))<=tol the algorithm ends.
    max_it_theta : maximun number of iterations for get the optimal theta 
    or the maximun descent direction in the theta varaible.

    Returns the solution and a resume.

    """
    
    print('Solving maximization of the p-Lplacian energy problem with method:'+method+' and p:'+str(p))
    
    bc=DirichletBC(V, 0.0, lambda x,onb: onb)
    def __get_initial_condition():
        """
        The initial conditions is obtained by solving the p-Laplacian equation
    
        $$ min \frac{1}{p}\int_{\Omega}|\nabla u|^{p}dx-\int_{\Omega}fudx$$
        
        return: u
        """
        u=TrialFunction(V)
        v=TestFunction(V)
        a= dot(grad(u),grad(v))*dx
        L=f*v*dx
        u=Function(V)
        solve(a==L,u,bc,solver_parameters={'linear_solver': 'gmres',
                         'preconditioner': 'ilu'})#, solver_parameters={'linear_solver': 'superlu_dist',
                        # 'preconditioner': 'ilu'})
        print('Obtaining initial condition')
        for i in range(max_it_init):
            d=__get_max_des_dir(u,theta=0)
            t=__line_search(d,u,theta=0)
            u.vector()[:]=u.vector()[:]+t*d.vector()[:]
            dF=assemble(dot(grad(d),grad(d))*dx)
            print('Iteration '+str(i)+' initial condition')
            print('Objective Function', __objective(u),'Grad norm initial condition:',dF,
                  'lengthstep:',t)
            if dF<= tol:
                return u
        return u
    def __objective(u,theta=0):
        """
        Objective Function
        """ 
        return assemble((dot(grad(u),grad(u))**(p/2)/(1+c*theta)**(p-1)/p-f*u)*dx)
    

    def __line_search(d,u,theta=0,dtheta=Constant(0),t=1):
        """ 
        Armijo line search
        """

        dF=assemble(dot(grad(d),grad(d))*dx)
        #if dtheta is not None:
        do=assemble(dot(grad(u),grad(u))**((p)/2)/(1+c*theta)**(p)*dtheta*dx)
        do=abs(do*c*(p-1)/p)

        dF+=do
        F=__objective(u,theta)
        F_new=__objective(u+t*d,theta+t*dtheta)
        s_err=F_new-F+alpha*dF*t
        while s_err>=0:
            
            t=t*beta
            F_new=__objective(u+t*d,theta+t*dtheta)
            s_err=F_new-F+alpha*dF*t
            #print('error',s_err,'t',t)
        return t
    def __get_max_des_dir(u,theta):
        """
        Return the Riez representation of dF/du, which correspond to the maximun descent direction
        with respect to the H_0^1 norm.
        
        The linear problem is solved with gmres method.
        """
        d=TrialFunction(V)
        v=TestFunction(V)
        
        w=project(dot(grad(u),grad(u)),Theta)
        w_vector=w.vector()[:]
        w_vector=w_vector**((p-2)/2)
        w_vector[(w_vector>=np.infty)&(np.isnan(w_vector))]=0
        w.vector()[:]=w_vector
        
        a=dot(grad(d),grad(v))*dx
        L=(-w*dot(grad(u),grad(v))/(1+c*theta)**(p-1)+v*f)*dx
        d=Function(V)

        solve(a==L,d,bc,solver_parameters={'linear_solver': 'gmres',
                         'preconditioner': 'ilu'})#, solver_parameters={'linear_solver': 'superlu_dist',
                       #  'preconditioner': 'ilu'})
        return d
    
    def __get_theta_opt(du,mu):
        """
        du: Dolfin Functon. For u, du corresponds to |\nabla u|.
        mu: Positive float.
        
        Returns theta=max(min(|\nabla u|/mu-1,1),0)/c
        """
        #du=project((dot(grad(u),grad(u)))**(0.5),self.Theta)
        theta=Function(Theta)
        theta_vec=(du.vector()[:]/mu-1)/c
        theta_vec[theta_vec<0]=0
        theta_vec[theta_vec>1]=1
        theta.vector()[:]=theta_vec
        return theta
    def __get_mu_theta(u):
        """
        u: Dolfin Function.
        
        Returns (mu,theta) , theta is optimal for u and mu is the corresponding multiplier. 

        """
        
        du=project((dot(grad(u),grad(u)))**(0.5),Theta)
        
        mu_up= assemble(dot(grad(u),grad(u))**(p/2)*dx)**(1/p)/k**(1/p)
        mu_inf=0
        mu=(mu_up+mu_inf)/2
        theta=__get_theta_opt(du,mu)
        int_theta=assemble(theta*dx)
        #while abs(int_theta-self.k)>self.eps:
        mu_prev=mu_inf
        i=0
        while abs(int_theta-k)>eps and max_it_theta>i:
            i+=1
            mu_prev=mu
            if int_theta>k:
                mu_inf=mu
            else:
                mu_up=mu
            mu=(mu_up+mu_inf)/2
            theta=__get_theta_opt(du,mu)
            int_theta=assemble(theta*dx)
        
        return mu,theta
    

    def __get_new_theta_max_dir(mu,dtheta_vec):
        new_theta=Function(Theta)
        new_theta_vec=new_theta.vector()[:]
        new_theta_vec[dtheta_vec>mu]=1
        new_theta_vec[dtheta_vec<=mu]=0
        new_theta.vector()[:]=new_theta_vec
        return new_theta

    def __get_max_dir_theta(u,theta):

        dtheta=project((dot(grad(u),grad(u)))**(p/2)/(1+c*theta)**(p),Theta)
        dtheta_vec=dtheta.vector()[:]

        unique_values=np.sort(np.unique(dtheta_vec))


        mu=0
        new_theta=__get_new_theta_max_dir(mu,dtheta_vec)
        int_theta=assemble(new_theta*dx)
        
        
        j_sup=unique_values.shape[0]
        j_inf=0
        i=0
        j=0

        while abs(int_theta-k)>eps and i<max_it_theta and j_sup>j_inf:
            i+=1
            
            if int_theta>k:
                j_inf=j+1
            else:
                j_sup=j
            j=int((j_sup+j_inf)/2)
            mu=unique_values[j]
            new_theta=__get_new_theta_max_dir(mu,dtheta_vec)
            int_theta=assemble(new_theta*dx)
            #print(i,'Area new:',int_theta,'j:',j,'j inf:',j_inf,'j sup:',j_sup)
        if int_theta<k:
            non_constant_theta=Function(Theta)
            aux=non_constant_theta.vector()[:]
            aux[dtheta_vec==mu]=1
            non_constant_theta.vector()[:]=aux
            area_non_constant_theta=assemble(non_constant_theta*dx)
            aux=new_theta.vector()[:]
            aux[dtheta_vec==mu]=(k-int_theta)/area_non_constant_theta
            new_theta.vector()[:]=aux
        return mu,new_theta
    
    
    ts=[]
    dFs=[]
    Fs=[]
    mus=[]
    u=__get_initial_condition()
    mu,theta=__get_mu_theta(u)

    intial_time=process_time()
    i=0
    while i<=max_it:
        
        print('Iteration:',i)
        
        
        d=__get_max_des_dir(u,theta)
        if method == 'alternating optimization':
            t=__line_search(d,u,theta)
            u.vector()[:]=u.vector()[:]+t*d.vector()[:]
            
            theta_ant=theta.copy()
            mu,theta=__get_mu_theta(u)
            
            dtheta=project(theta-theta_ant,Theta)
            do=assemble(dot(grad(u),grad(u))**((p)/2)/(1+c*theta)**(p)*dtheta*dx)*c*(p-1)/p
            if do<0:
                do=0
                theta=theta_ant.copy()
            
            
        elif method=='feasible directions':
            
            mu,new_theta=__get_max_dir_theta(u,theta)
            dtheta=project(new_theta-theta,Theta)
            do=assemble(dot(grad(u),grad(u))**((p)/2)/(1+c*theta)**(p)*dtheta*dx)*c*(p-1)/p
            are_new_theta=assemble(new_theta*dx)

            if do<0 :
                dtheta=Function(Theta)
                do=0
            
            t=__line_search(d,u,theta,dtheta)
            u.vector()[:]=u.vector()[:]+t*d.vector()[:]
            theta.vector()[:]=theta.vector()[:]+t*dtheta.vector()[:]
        dF=assemble(dot(grad(d),grad(d))*dx)+do

            

        area_theta=assemble(theta*dx)
        F=__objective(u,theta)
        dFs.append(dF)
        ts.append(t)
        if method=='feasible directions':
            mu=mu**(1/p)
        mus.append(mu)

        Fs.append(F)
        print('Method',method,'Objective',F,'Grad_norm:',dF,'Dtheta:',do,'length step:',t,'area theta:',area_theta)
        dijitso.cache.clean_cache(cache_params, dryrun=False, 
                                  categories=(u'inc', u'src', u'lib', u'log'))

        if dF<tol:
            print('Toleranced reached in ',i,' iterations.')
            break
        i+=1
    total_time=process_time()-intial_time
    if i >max_it:
        print('Maximun number of iterations reached.')
    print('Problem solved in ',total_time,'seconds')
    resp={}
    resp['lenthg step']=ts
    resp['grad norm']=dFs
    resp['objective']=Fs
    resp['lagrange multiplier']=mus
    resp['u vec']=u.vector()[:]
    resp['theta vec']=theta.vector()[:]
    resp['total time']=total_time
    return resp
 
    
